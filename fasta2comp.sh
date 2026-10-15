#!/bin/bash
comp_usage="$(dirname $0)/comp_usage.sh"
source ${comp_usage}
check_dependencies(){
    for dep in cusp geecee seqtk parallel
    do
        if ! command -v ${dep} &> /dev/null
        then
            echo "Error: ${dep} is not installed"
            exit 1
        fi
    done
    }
export -f check_dependencies
check_dependencies
while [ "$1" != "" ]
do
    case $1 in
        --fasta_file     )
            shift
            fasta_file=$(realpath $1)
            ;;
        --out_file )
            shift
            out_file=$1
            ;;
        --threads )
            shift
            threads=$1
            ;;
		--help          )
            comp_usage
            exit 0
            ;;
	esac
	shift
done

if [ ! -f "${fasta_file}" ]
then
    echo "CDS fasta file ${fasta_file} missing. Exiting"
    exit 2
fi

num_proc=$(nproc)
if [ -z "${threads}" ]
then
    threads=${num_proc}
else
    threads_test=$(echo -e "${threads}\t${num_proc}" | awk '{if((int($1)==$1) && ($1<=$2)){print $1}}')
    if [ -z "${threads_test}" ]
    then
        echo "Invalid threads value. Exiting"
        exit 3
    fi
fi

get_codon_stats(){
    gene_id=$1
    fasta_file=$2
    fasta_seq=$(seqtk subseq ${fasta_file} <(echo "${gene_id}") )
    cusp_out=$(echo "${fasta_seq}" | cusp -sequence /dev/stdin stdout 2> /dev/null | grep -v "#" | grep ^[A-Z] | grep -v "*" | awk 'BEGIN{OFS="\t"}{print $1,$2,$5}')
    aa_list="A C D E F G H I K L M N P Q R S T V W Y"
    frac_list="${gene_id}"
    for aa in ${aa_list}
    do
        codon_list=$(echo "${cusp_out}" | grep -w "${aa}" | cut -f1 | sort -V | uniq)
        aa_count=$(echo   "${cusp_out}" | awk -v aa="${aa}" 'BEGIN{FS="\t"}{if($2==aa){sum+=$3}}END{print sum}' )
        for codon in ${codon_list}
        do
            if [ "${aa_count}" -eq 0 ]
            then
                codon_frac="0"
            else
                codon_frac=$(echo "${cusp_out}" | grep -w ^"${codon}" | cut -f3 | awk -v aa_count="${aa_count}" '{print $1/aa_count}')
            fi
            frac_list=$(echo -e "${frac_list}\t${codon_frac}")
        done
        frac_list=$(echo "${frac_list}" | perl -pe 's/\t$//')
    done
    echo "${frac_list}"
    }
export -f get_codon_stats
uuid=$(uuidgen | cut -d\- -f5)
id_list=$(grep ">" "${fasta_file}" | perl -pe 's/\>//;s/ .*//' | sort -V | uniq) 
id_count=$(echo "${id_list}" | wc -l)
id_iters=$(echo ${id_count} | awk '{modulo=$1%1000}{if(modulo==0){print $1/1000}else{print (($1-modulo)/1000)+1 }}')
for iter in $(seq 1 ${id_iters})
do
    start_line=$(echo "${iter}" | awk '{print (($1-1)*1000)+1}')
    id_sub_list=$(echo "${id_list}" | tail -n+${start_line} | head -n1000)
    parallel -j ${threads} get_codon_stats ::: "${id_sub_list[@]}" ::: "${fasta_file}" >> ${uuid}cusp.tsv
done
exit_code=$?
if [ "${exit_code}" -ne 0 ]
then
    echo "Error running cusp in parallel"
    exit ${exit_code}
fi
echo "$(sort -V ${uuid}.cusp.tsv | uniq )" > ${uuid}.cusp.tsv
geecee ${fasta_file} stdout 2>/dev/null | grep -v "#" | awk 'BEGIN{OFS="\t"}{print $1,$2}' | sort -V | uniq > ${uuid}.geecee.tsv
exit_code=$?
if [ "${exit_code}" -ne 0 ]
then
    echo "Error running geecee"
    exit ${exit_code}
fi
echo -e "id\t%G+C\tGCA\tGCC\tGCG\tGCT\tTGC\tTGT\tGAC\tGAT\tGAA\tGAG\tTTC\tTTT\tGGA\tGGC\tGGG\tGGT\tCAC\tCAT\tATA\tATC\tATT\tAAA\tAAG\tCTA\tCTC\tCTG\tCTT\tTTA\tTTG\tATG\tAAC\tAAT\tCCA\tCCC\tCCG\tCCT\tCAA\tCAG\tAGA\tAGG\tCGA\tCGC\tCGG\tCGT\tAGC\tAGT\tTCA\tTCC\tTCG\tTCT\tACA\tACC\tACG\tACT\tGTA\tGTC\tGTG\tGTT\tTGG\tTAC\tTAT" > "${out_file}"
paste ${uuid}.geecee.tsv ${uuid}.cusp.tsv | awk 'BEGIN{FS="\t"}{if($1==$3){print $0}}' | cut -f3 --complement >> "${out_file}"
rm ${uuid}.cusp.tsv ${uuid}.geecee.tsv