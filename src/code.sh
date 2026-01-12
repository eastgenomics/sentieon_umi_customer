#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -eo pipefail
apt-get update && apt-get install -y gnuplot
# Set up options
source /home/dnanexus/license_setup.sh

export LD_LIBRARY_PATH=/opt/dnanexus/assets/htslib/lib:$LD_LIBRARY_PATH

wait_pids()
{
    # wait for all pids in the array to finish with non-zero exit status
    for pid in "$@"; do
        echo "Waiting for $pid"
        wait "$pid"
    done
}

stream_input()
{
    INPUT=$1
    if [ "$INPUT" != "" ]; then
        dx cat "$INPUT"
    else
        echo ""
    fi
}

wait_uploads()
{
    # wait until all $FILE_UPLOADED files are present
    #TODO: add timeout?
    arr=("$@")
    FINISHED=1
    TIME=30
    until [ $FINISHED -eq 0 ]; do
        FINISHED=0
        for file in "${arr[@]}"; do
            if [ ! -e "$file" ]; then
                #check if error
                if [ -e "${file}.error" ]; then
                    #Do not fully error, as there are some files that are not required outputs; for required outputs, DNAnexus will take care of erroring out
                    echo -n "FILE UPLOAD ERROR: "; cat "${file}.error"
                    #after reporting error, stop waiting for this file
                    touch "$file"
                    #error_report < cat "${file}.error"
                else
                    FINISHED=$((FINISHED + 1)) #if at least 1 file missed, it will continue
                fi
            fi
        done
        #TODO: if $FINISHED > 0 check if dx upload is running, if not error out?
        if [ $FINISHED -ne 0 ]; then
            if [ $(ps -ef|grep "dx upload" |grep -v "grep dx upload" -c) -eq 0 ]; then
                echo "No dx upload Running, error out in $TIME tries"
                TIME=$((TIME - 1))
                if [ $TIME -eq 0 ]; then
                    #Do not fully error, as there are some files that are not required outputs; for required outputs, DNAnexus will take care of erroring out
                    echo "FILE UPLOAD ERROR: There has been some issue with the uploads, as some are still pending but there has been no dx upload running for some time"
                    FINISHED=0
                    #error_report "There has been some issue with the uploads, as some are still pending but there has been no dx upload running for some time"
                fi
             fi
        fi
        echo "Still waiting for $FINISHED uploads"
        sleep 2
    done
}

upload()
{
    HANDLE="$1"
    INPUT="$2"
    OUTPUT_FILENAME="$3"
    OUTPUT_UPLOAD_EXCEPT="$OUTPUT_UPLOAD_EXCEPT --except $HANDLE"
    wait_uploads_arr+=("uploaded_$HANDLE.file")
    handle_uploads_arr+=("$HANDLE")
    if [ ! -e "$INPUT" ]; then
        error_report "File $INPUT not found"
    else
        if [ "$output_md5sum" == "true" ] && [ "${HANDLE:0:12}" != "md5sumS_OUT_" ]; then
            #create and upload md5sum file and wait for it to be uploaded
            #do not exclude it from or add to handle_upload as the output does not exist alone but as part of the md5sum output folder
            wait_uploads_arr+=("uploaded_md5sumS_OUT_$HANDLE.file")
            #do not put md5sum in background, as subsequent commands may move things around
            md5sum "$INPUT" > "$INPUT".md5
            if [ $? -eq 0 ]; then
                streaming_upload md5sumS_OUT_"$HANDLE" "$INPUT".md5 "$OUTPUT_FILENAME".md5 &
            else
                echo "Error calculating Md5sum for $OUTPUT_FILENAME" > uploaded_md5sumS_OUT_$HANDLE.file.error
            fi
            #create_md5_and_upload md5sumS_OUT_$HANDLE $INPUT $OUTPUT_FILENAME &
        fi
        streaming_upload "$HANDLE" "$INPUT" "$OUTPUT_FILENAME" &
    fi
}

streaming_upload()
{
    HANDLE="$1"
    INPUT="$2"
    OUTPUT_FILENAME="$3"
    if [ "${HANDLE:0:12}" == "md5sumS_OUT_" ]; then
        OUTPUT_FOLDER="./out/md5sums";
    else
        OUTPUT_FOLDER="./out/$HANDLE"
    fi
    OUTPUT_FILE="$OUTPUT_FOLDER/$OUTPUT_FILENAME"
    FILE_ID_FILE="file_id_$HANDLE.file"
    FILE_UPLOADED="uploaded_$HANDLE.file" #used to trigger potential removal

    rm -f "${FILE_UPLOADED}".error
    rm -f "$FILE_UPLOADED"
    mkdir -p "$OUTPUT_FOLDER"
    mv "$INPUT" "$OUTPUT_FILE"
    ln -s "$OUTPUT_FILE $INPUT"
    dx upload "$OUTPUT_FILE" --brief > "$FILE_ID_FILE" && touch "$FILE_UPLOADED" || (echo "Error uploading output $OUTPUT_FILENAME" > uploaded_$HANDLE.file.error)
}

########################
### Pipeline modules ###
########################
module_download_inputs()
{
	# Download all inputs 
	EXCEPT_LIST="$EXCEPT_LIST --except genome_fastagz --except mappings_bam"
	dx-download-all-inputs --parallel $EXCEPT_LIST

    input_command="mv -f --"
    if [[ "$in_files" != "" ]]; then
        eval $input_command "${in_files_path[@]}" .
    fi

	# Stream and unpack FASTA genome while using samtools to index it
	wait_pids_arr=()
	mkdir genome
	genome_file="genome/${genome_fastagz_name%.gz}"
	mkfifo genome_temp.fa
	dx download "${genome_fastagz}" -o - | gunzip -c | tee genome_temp.fa > "${genome_file}" &
	gunzipGenomeFastagz_pid=$!
	wait_pids_arr+=(${gunzipGenomeFastagz_pid})
	samtools faidx genome_temp.fa &
	samtools_genome_faidx_pid=$!
	wait_pids_arr+=(${samtools_genome_faidx_pid})
	if [ "$genomebwaindex_targz" != "" ]; then
		tar --no-same-owner -zxvf "$genomebwaindex_targz_path" -C genome &
		untarGenomeIndexTargz_pid=$!
		wait_pids_arr+=(${untarGenomeIndexTargz_pid})
	fi
	
	# Wait until untar commands are done
	wait_pids "${wait_pids_arr[@]}"
	unset wait_pids_arr
	wait_pids_arr=()
	
	# Rename genome file
	genome_file_prefix="${genomebwaindex_targz_name%.bwa-index.tar.gz}"
	rename "s/genome\.fa/${genome_file_prefix}\.fa/" ./genome/genome* || true
	mv genome_temp.fa.fai "${genome_file}.fai"
	rm genome_temp.fa
	
	#check if the BWA indices and the FASTA match
	if [ "$genomebwaindex_targz" != "" ]; then
		found_right=false
		for bwt_file in genome/${genome_file_prefix}.fa*.bwt; do
			if [ "${bwt_file%.bwt}" == "${genome_file}" ]; then
				found_right=true
			fi
		done
		if [ "$found_right" == "false" ]; then
			error_report "The BWA reference genome index does not correspond to the Reference genome FASTA file. The app will exit now."
		fi
	fi
	
	# Set sample name if not given
	if [ "$sample" = "" ]; then
		sample="${reads_fastqgzs_prefix[0]}"
		sample="$(echo $sample|sed 's|_00[0-9]$||'|sed 's|_R1$||'|sed 's|_1$||')"
	else
		sample="${sample// /_}"
	fi
	if [ "$sample" = "" ]; then
		sample="sample"
	fi
}

module_umi_align()
#output is ${dedup_bam_name}.${dedup_bam_extension} and ${CASE_PREFIX}metrics_files
#expects reads_fastqgzs_name reads_fastqgzs_prefix reads_fastqgzs_path and reads2_fastqgzs_path reads_fastqgzs and reads2_fastqgzs arrays
{   	
	if [ "${CASE_NAME}" != "" ]; then
		CASE_PREFIX="${CASE_NAME}_"
		CASE_POSTFIX="_${CASE_NAME}"
	fi
	CASE_SAMPLE="$sample$CASE_POSTFIX"	
	# Run bwa mem
	# Determine correct bwt mem
	mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
	mem_kb=$((mem_kb<400000000 ? mem_kb : 400000000)) #, but limit to reduce loading time in very high memory machines?
	export bwt_max_mem="$((mem_kb / 1024 / 1024 - 6))g"
	# If read group information CSV file exists: (1) Load read group information into arrays (2) Check if all input filenames can map to a row of read group info
	declare -A read_group_id_array 
	declare -A read_group_library_array 
	declare -A read_group_pu_array
	if [[ "$rg_info_csv_path" != "" ]]; then
		tr -d '\r' < "${rg_info_csv_path}" > "rg_info_csv_wo_carriagereturns.csv"
		# Load the CSV file into dictionary arrays, using the filename of the first read mates as key
		i=0
		while read -r line; do
			i=$(( i + 1 ))
			IFS=',' read -r -a rg_info <<< "$line"
			if [[ ${#rg_info[@]} -ne 4 ]]; then
				error_report "${rg_info_csv_name} is misformatted. Row ${i} does not have 4 columns (filename, RGID, RGLB, RGPU). Please check the read group information file \"${rg_info_csv_name}\" before re-run."
				#dx-jobutil-report-error "${rg_info_csv_name} is misformatted. Row ${i} does not have 4 columns (filename, RGID, RGLB, RGPU). Please check the read group information file \"${rg_info_csv_name}\" before re-run."
			fi
			filename="$(echo -e "${rg_info[0]}" | sed -e 's/[[:space:]]*$//')";
			read_group_id_array[${filename}]="$(echo -e "${rg_info[1]}" | sed -e 's/[[:space:]]*$//')";
			read_group_library_array[${filename}]="$(echo -e "${rg_info[2]}" | sed -e 's/[[:space:]]*$//')";
			read_group_pu_array[${filename}]="$(echo -e "${rg_info[3]}" | sed -e 's/[[:space:]]*$//')";
		done < "rg_info_csv_wo_carriagereturns.csv"
	
		# Check if all input filanems have corresponding read group info in the CSV file
		for i in ${!reads_fastqgzs_name[@]}; do
			match=0
			filename=${reads_fastqgzs_name[$i]}
			for key in "${!read_group_id_array[@]}"; do 
				if [ "${filename}" = "${key}" ]; then
					match=1
					break
				fi
			done
	
			if [[ "$match" -eq 0 ]]; then 
				error_report "Cannot find read group info for FASTQ file \"${filename}\". Please check the read group information file \"${rg_info_csv_name}\" before re-run."
				#dx-jobutil-report-error "Cannot find read group info for FASTQ file \"${filename}\". Please check the read group information file \"${rg_info_csv_name}\" before re-run."
			fi 
		done
	fi
	#run UMI pipeline per pair of inputs
	number_fastq=${#reads_fastqgzs_name[@]}
	merged_bam_list=""
	if [ "$duplex_umi" = "true" ]; then
		duplex_arg="-d"
	fi

	sorted_bam_args=()
	for i in ${!reads_fastqgzs_name[@]}; do
		filename=${reads_fastqgzs_name[$i]}
		mock_lane=$(( i+1 ))
		read_group_id=${read_group_id_array[${filename}]}
		# If no read group info provided, generate read group info with the following logic
		if [ "$read_group_id" = "" ]; then
			read_group_id="${reads_fastqgzs_prefix[$i]}"
			read_group_id="${read_group_id/_R1/}"
			read_group_id="${read_group_id/_1/}"
			read_group_id="${read_group_id}_${mock_lane}"
		fi
		read_group_library=${read_group_library_array[${filename}]}
		# If no read group info provided, generate read group info with the following logic
		if [ "$read_group_library" = "" ]; then
			read_group_library="${CASE_SAMPLE}"
		fi
		read_group_pu=${read_group_pu_array[${filename}]}
		# If no read group info provided, generate read group info with the following logic
		if [ "$read_group_pu" = "" ]; then
			read_group_pu=${read_group_id}
		fi
        # Handle the optional index read
        read_idx_path=""
        if [ ${#readsidx_fastqgzs_path[@]} -ne 0 ]; then
            read_idx_path="${readsidx_fastqgzs_path[$i]}"
        fi		

		# Extract + aln with bwa
		if [ "$stream_inputs" = "true" ]; then
			{ (${SENTIEON_INSTALL_DIR}/bin/sentieon umi extract $duplex_arg ${read_template} $read_idx_path \
				<(stream_input "${reads_fastqgzs[$i]}") <(stream_input "${reads2_fastqgzs[$i]}") \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -t ${sentieon_procs} $extra_bwa_options \
				"${genome_file}" ${bwa_options} \
				-R "@RG\tID:${read_group_id}\tPL:${read_group_platform}\tPU:${read_group_pu}\tLB:${read_group_library}\tSM:${CASE_SAMPLE}" \
				-p -C - || echo -n 'error' ) 2>&3 \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon util sort -i - --sam2bam -o sorted_${i}.bam \
				-t ${sentieon_procs} ${util_sort_options} --bam_compression $bam_compression --block_size 2G ;}\
				3>&1 1>&2 | grep -v --line-buffered "^\[M::mem_pestat" \
				| grep -v --line-buffered "^\[M::process"| ([ "$bwa_nonverbose" != "true" ] \
				&& cat || grep -v --line-buffered "^\[M::mem_process" )
		else
			{ (${SENTIEON_INSTALL_DIR}/bin/sentieon umi extract $duplex_arg ${read_template} $read_idx_path \
				"${reads_fastqgzs_path[$i]}" "${reads2_fastqgzs_path[$i]}" \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -t ${sentieon_procs} $extra_bwa_options \
				"${genome_file}" ${bwa_options} \
				-R "@RG\\tID:${read_group_id}\\tPL:${read_group_platform}\\tPU:${read_group_pu}\\tLB:${read_group_library}\\tSM:${CASE_SAMPLE}" \
				-p -C - || echo -n 'error' ) 2>&3 \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon util sort -i - --sam2bam -o sorted_${i}.bam \
				-t ${sentieon_procs} ${util_sort_options} --bam_compression $bam_compression --block_size 2G ;}\
				3>&1 1>&2 | grep -v --line-buffered "^\[M::mem_pestat" \
				| grep -v --line-buffered "^\[M::process"| ([ "$bwa_nonverbose" != "true" ] \
				&& cat || grep -v --line-buffered "^\[M::mem_process" )
			rm ${reads_fastqgzs_path[$i]}
			if [ ${#reads2_fastqgzs_path[@]} -ne 0 ]; then
				rm ${reads2_fastqgzs_path[$i]}
			fi
		fi
		sorted_bam_args+=("-i")
		sorted_bam_args+=("sorted_${i}.bam")
		number_fastq=$((number_fastq - 1))
		echo "$number_fastq remaining FASTQ files to be mapped"
	done

	# LC with `--consensus` and `--umi_tag` + Dedup
	# Output a deduped BAM + metrics on the deduped file
	metrics_algo_args=()
	if [ "${output_metrics}" = "true" ]; then
		metrics_algo_args=(
			"--algo" "CoverageMetrics" "--omit_base_output" "${CASE_PREFIX}coverage"
			"--algo" "GCBias" "--summary" "${CASE_PREFIX}gc_summary.txt" "${CASE_PREFIX}gc_metric.txt"
			"--algo" "MeanQualityByCycle" "${CASE_PREFIX}mq_metric.txt"
			"--algo" "QualDistribution" "${CASE_PREFIX}qd_metric.txt"
			"--algo" "InsertSizeMetricAlgo" "${CASE_PREFIX}is_metric.txt"
			"--algo" "AlignmentStat" "${CASE_PREFIX}aln_metric.txt"
		)
	fi

	if [ "$output_format" = "CRAM" ]; then
		output_bam_extension="cram"
		output_bai_extension="crai"
	else
		output_bam_extension="bam"
		output_bai_extension="bai"
	fi
	dedup_bam_extension=${output_bam_extension}
    dedup_bam_name="${CASE_PREFIX}sorted"

	${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t ${sentieon_procs} -r "${genome_file}" \
		"${sorted_bam_args[@]}" --algo LocusCollector --consensus --umi_tag XR \
		--fun score_info score.txt.gz  "${metrics_algo_args[@]}"
	${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t ${sentieon_procs} -r "${genome_file}" \
		"${sorted_bam_args[@]}" --algo Dedup --score_info score.txt.gz --metrics \
		"${CASE_PREFIX}dedup_metrics.txt" "${dedup_bam_name}.${dedup_bam_extension}"

	rm sorted_*.bam || true
}


################################
###      data pipeline       ###
################################

set -x
# Error checking: Check that there are equal number of FQ1 and FQ2
if [ ${#reads2_fastqgzs_path[@]} -ne 0 ] && [ ${#reads_fastqgzs_path[@]} -ne ${#reads2_fastqgzs_path[@]} ]; then
	dx-jobutil-report-error "The number of paired FASTQ files are not the same. First read mate has ${#reads_fastqgzs_path[@]} FASTQ files while second read mate has ${#reads_fastqgzs_path[@]} FASTQ files. Please provide matched FASTQ files. The app will exit now."
fi
if [ ${#readsidx_fastqgzs_path[@]} -ne 0 ] && [ ${#reads_fastqgzs_path[@]} -ne ${#readsidx_fastqgzs_path[@]} ]; then
	dx-jobutil-report-error "The number of UMI index FASTQ files are not the same. First read mate has ${#reads_fastqgzs_path[@]} FASTQ files while the umi index read has ${#readsidx_fastqgzs_path[@]} FASTQ files. Please provide matched FASTQ files. The app will exit now."
fi

EXCEPT_LIST=""
if [ "$stream_inputs" = "true" ]; then
	EXCEPT_LIST="$EXCEPT_LIST --except reads_fastqgzs"
	EXCEPT_LIST="$EXCEPT_LIST --except reads2_fastqgzs"
fi

CASE_NAME="" #may use this to differentiate between tumor and normal
#output is ${dedup_bam_name}.${dedup_bam_extension} and ${CASE_PREFIX}metrics_files
#expects reads_fastqgzs_name reads_fastqgzs_prefix reads_fastqgzs and reads2_fastqgzs arrays
module_download_inputs
module_umi_align

# Upload outputs
wait_uploads_arr=()
handle_uploads_arr=()
if [ "${dedup_bam_extension}" != "${output_bam_extension}" ]; then
	error_report "Something went wrong and the output file does not have the right extension"
fi
upload mappings_bam "${dedup_bam_name}.${dedup_bam_extension}" "${sample}_${dedup_bam_name}.${output_bam_extension}"
upload mappings_bam_bai "${dedup_bam_name}.${dedup_bam_extension}.${output_bai_extension}" "${sample}_${dedup_bam_name}.${output_bam_extension}.${output_bai_extension}"

mkdir -p ~/out/metrics/${sample}_metrics
mv "${CASE_PREFIX}dedup_metrics.txt" ~/out/metrics/${sample}_metrics/"${sample}".dedup_metrics.txt
if [ "${output_metrics}" = "true" ]; then
        ${SENTIEON_INSTALL_DIR}/bin/sentieon plot metrics -o ${CASE_PREFIX}metrics.pdf gc=${CASE_PREFIX}gc_metric.txt mq=${CASE_PREFIX}mq_metric.txt qd=${CASE_PREFIX}qd_metric.txt isize=${CASE_PREFIX}is_metric.txt
        mv ${CASE_PREFIX}gc_metric.txt ~/out/metrics/${sample}_metrics/"$sample".GCBias_metrics.txt
        mv ${CASE_PREFIX}gc_summary.txt ~/out/metrics/${sample}_metrics/"$sample".GCBiasSummary_metrics.txt
        mv ${CASE_PREFIX}mq_metric.txt ~/out/metrics/${sample}_metrics/"$sample".MeanQualityByCycle_metrics.txt
        mv ${CASE_PREFIX}qd_metric.txt ~/out/metrics/${sample}_metrics/"$sample".QualDistribution_metrics.txt
        mv ${CASE_PREFIX}is_metric.txt ~/out/metrics/${sample}_metrics/"$sample".InsertSize_metrics.txt
        mv ${CASE_PREFIX}aln_metric.txt ~/out/metrics/${sample}_metrics/"$sample".AlignmentStat_metrics.txt
        mv ${CASE_PREFIX}metrics.pdf ~/out/metrics/${sample}_metrics/"$sample".metrics.pdf
		ff=(${CASE_PREFIX}coverage*)
	for f in ${ff[@]}; do 
		if [ -z "${CASE_PREFIX}" ]; then
			f1=${sample}.$f
		else
			f1=${sample}.${f##${CASE_PREFIX}}
		fi
		mv $f ~/out/metrics/${sample}_metrics/$f1
	done
fi
wait_uploads "${wait_uploads_arr[@]}"
for HANDLE in "${handle_uploads_arr[@]}"; do
        dx-jobutil-add-output "$HANDLE" "$(cat file_id_$HANDLE.file)" --class=file
done

# Upload remaining results
dx-upload-all-outputs --parallel $OUTPUT_UPLOAD_EXCEPT

