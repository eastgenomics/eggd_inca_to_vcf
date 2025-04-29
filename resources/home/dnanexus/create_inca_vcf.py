import pandas as pd
import argparse
import subprocess

# TODO: move to separate file?
INFO_FIELDS = {
    'LATEST_GERMLINE': 'latest_germline',
    'LATEST_ONCOGENICITY': 'latest_oncogenicity',
    'LATEST_DATE': 'latest_date',
    'LATEST_SAMPLE_ID': 'latest_sample_id',
    'TOTAL_CLASSIFICATIONS': 'total_classifications',
    'AGGREGATED_HGVS': 'aggregated_hgvs'}

GRCh37_CONTIG = '''\
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>\
'''

GRCh38_CONTIG = '''\
##contig=<ID=1,length=248956422,assembly=hg38>
##contig=<ID=2,length=242193529,assembly=hg38>
##contig=<ID=3,length=198295559,assembly=hg38>
##contig=<ID=4,length=190214555,assembly=hg38>
##contig=<ID=5,length=181538259,assembly=hg38>
##contig=<ID=6,length=170805979,assembly=hg38>
##contig=<ID=7,length=159345973,assembly=hg38>
##contig=<ID=8,length=145138636,assembly=hg38>
##contig=<ID=9,length=138394717,assembly=hg38>
##contig=<ID=10,length=133797422,assembly=hg38>
##contig=<ID=11,length=135086622,assembly=hg38>
##contig=<ID=12,length=133275309,assembly=hg38>
##contig=<ID=13,length=114364328,assembly=hg38>
##contig=<ID=14,length=107043718,assembly=hg38>
##contig=<ID=15,length=101991189,assembly=hg38>
##contig=<ID=16,length=90338345,assembly=hg38>
##contig=<ID=17,length=83257441,assembly=hg38>
##contig=<ID=18,length=80373285,assembly=hg38>
##contig=<ID=19,length=58617616,assembly=hg38>
##contig=<ID=20,length=64444167,assembly=hg38>
##contig=<ID=21,length=46709983,assembly=hg38>
##contig=<ID=22,length=50818468,assembly=hg38>
##contig=<ID=X,length=156040895,assembly=hg38>
##contig=<ID=Y,length=57227415,assembly=hg38>
##contig=<ID=M,length=16569,assembly=hg38>\
'''

MINIMAL_VCF_HEADER = '''\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
'''


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(description="Create VCF from previous interpretations")

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        help=(
            "CSV file exported from the previous interpretations "
            "database"
        ),
    )

    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        help=(
            "Output VCF filename"
        ),
    )

    parser.add_argument(
        "-g",
        "--genome_build",
        type=str,
        choices=["GRCh37", "GRCh38"],
        help="Genome build the samples were run in",
    )

    # TODO: str best method for inputting multiple probesets?
    parser.add_argument(
        "-set",
        "--probeset",
        type=str,
        nargs='+',
        help=(
            "probset_id or allele_origin to filter"
        ),
    )

    args = parser.parse_args()

    return args


def clean_csv(input_file):
    '''
    Clean up the Inca database CSV by: 
    - Convert to tab separated instead of comma
    - Rename to CHROM, POS, REF, ALT
    - Move CHROM, POS, REF, ALT to be first 4 columns
    - Remove all new lines or tabs within cells

    Parameters
    ----------
    input_file : str
        Filepath to Inca CSV

    Returns
    -------
    pd.DataFrame
        Dataframe with cleaned up data
    '''
    # TODO: add required_columns?
    df = pd.read_csv(
        input_file,
        delimiter=",",
        parse_dates=['date_last_evaluated'],
        low_memory=False)
    df.rename(columns={"chromosome": "CHROM", 
                       "start": "POS", 
                       "reference_allele": "REF", 
                       "alternate_allele": "ALT"}, inplace=True)
    df = df[["CHROM", "POS", "REF", "ALT"] + [ col for col in df.columns if col not in ["CHROM", "POS", "REF", "ALT"]]]
    df = df.applymap(lambda x: x.replace("\n", " ").strip() if isinstance(x, str) else x)

    return df


def filter_probeset(cleaned_csv, probeset):
    '''
    Filter cleaned data to interpreted variants for specified germline/somatic probesets

    Parameters
    ----------
    cleaned_csv : pd.DataFrame
        Dataframe with cleaned up data
    probeset : str
        Germline or somatic choice

    Returns
    -------
    pd.DataFrame
        Dataframe filtered by probeset
    '''
    # TODO: ensure probeset choices are only the 4
    interpreted_df = cleaned_csv[cleaned_csv['interpreted'].str.lower() == "yes"]
    if (interpreted_df['germline_classification'].isnull() & interpreted_df['oncogenicity_classification'].isnull()).any():
        raise ValueError("Both germline and oncogenicity classification are null in at least one row.")

    all_dfs = []
    for type in probeset:
        if type in ["99347387", "96527893"]:
            column = "probeset_id"
        else:
            column = "allele_origin"
        filtered_df = interpreted_df.loc[interpreted_df[column] == type]
        all_dfs.append(filtered_df)
    
    probeset_df = pd.concat(all_dfs, ignore_index=True)
    probeset_df = probeset_df.drop_duplicates()

    return probeset_df


def get_latest_entry(sub_df):
    '''
    Get latest entry by date

    Parameters
    ----------
    sub_df : pd.DataFrame
        Dataframe per group of CHROM, POS, REF, ALT

    Returns
    -------
    pd.Series
        Latest entry per group by date
    '''
    latest_idx = sub_df['date_last_evaluated'].idxmax()
    latest_entry = sub_df.loc[latest_idx]
    return latest_entry


def aggregate_hgvs(hgvs_series):
    '''
    Aggregates all unique HGVS

    Parameters
    ----------
    hgvs_series : pd.Series
        HGVSc per variant

    Returns
    -------
    str
        All HGVSc per variant joined
    '''
    unique_hgvs = hgvs_series.dropna().unique()
    return "|".join(unique_hgvs)


def format_total_classifications(classifications):
    '''
    Counts all classifications, including the latest classification.
    Returns classifications in the format: classification(count)|classification(count)

    Parameters
    ----------
    classifications : pd.Series
        Germline or somatic classifications per variant

    Returns
    -------
    str
        All classifications per variant joined
    '''
    counts = classifications.value_counts()
    formatted_counts = [f"{classification}({count})" for classification, count in counts.items()]
    return "|".join(formatted_counts)


def sort_aggregated_data(aggregated_df):
    '''
    Sort aggregate data
    
    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data

    Returns
    -------
    pd.DataFrame
        Dataframe sorted by CHROM and POS
    '''
    # Define chromosome order: numeric first, then X and Y and sort
    chromosome_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    aggregated_df['CHROM'] = pd.Categorical(aggregated_df['CHROM'], categories=chromosome_order, ordered=True)
    aggregated_df = aggregated_df.sort_values(by=['CHROM', 'POS'])

    return aggregated_df


def aggregate_uniq_vars(probeset_df, probeset, aggregated_database):
    '''
    Aggregate data for each unique variant
    Similaritites to create_vcf_from_inca_csv.py by Raymond Miles
    
    Parameters
    ----------
    probeset_df : pd.DataFrame
        Dataframe filtered by probeset
    probeset : str
        Germline or somatic choice
    aggregated_database : str
        Output filename for aggregated data

    Returns
    -------
    pd.DataFrame
        Dataframe of aggregated data
    '''
    probeset_df.loc[:, 'germline_classification'] = probeset_df['germline_classification'].str.replace(' ', '_')
    probeset_df.loc[:, 'oncogenicity_classification'] = probeset_df['oncogenicity_classification'].str.replace(' ', '_')
    probeset_df.loc[:, 'CHROM'] = probeset_df['CHROM'].str.replace(' ', '')
    probeset_df = probeset_df.dropna(subset=['date_last_evaluated'])

    aggregated_data = []
    for type in probeset:
        if type in ["germline", "99347387"]:
            classification = "germline"
        else:
            classification = "oncogenicity"

        grouped = probeset_df.groupby(['CHROM', 'POS', 'REF', 'ALT', f'{classification}_classification'])

        for _, group in grouped:
            latest_entry = get_latest_entry(group)
            latest_germline = latest_entry['germline_classification']
            latest_oncogenicity = latest_entry['oncogenicity_classification']
            latest_date = latest_entry['date_last_evaluated']
            latest_sample_id = latest_entry['specimen_id']
            hgvs = aggregate_hgvs(group['hgvsc'])
            # TODO: do we need a column for each total classification?
            total_classifications = format_total_classifications(group[f'{classification}_classification'])

            aggregated_data.append({
                'CHROM': latest_entry['CHROM'],
                'POS': latest_entry['POS'],
                'REF': latest_entry['REF'],
                'ALT': latest_entry['ALT'],
                'latest_germline': latest_germline,
                'latest_oncogenicity': latest_oncogenicity,
                'latest_date': latest_date,
                'latest_sample_id': latest_sample_id,
                'total_classifications': total_classifications,
                'aggregated_hgvs': hgvs
        })

    aggregated_df = pd.DataFrame(aggregated_data)
    aggregated_df = sort_aggregated_data(aggregated_df)
    aggregated_df.to_csv(aggregated_database, sep="\t", index=False)

    return aggregated_df


def intialise_vcf(aggregated_df, minimal_vcf):
    ''' 
    Initialise minimal VCF with CHROM, POS, ID, REF, ALT with minimal header

    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data
    minimal_vcf : str
        Output filename for the minimal VCF
    '''
    vcf_lines = []
    for _, row in aggregated_df.iterrows():
        vcf_line = f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\t."
        vcf_lines.append(vcf_line)

    with open(minimal_vcf, "w") as vcf_file:
        vcf_file.write(MINIMAL_VCF_HEADER)
        vcf_file.write("\n".join(vcf_lines) + "\n")


def write_vcf_header(genome_build):
    '''
    Write VCF header by populating INFO fields and specifying contigs

    Parameters
    ----------
    genome_build : str
        Genome build to specify contigs in header
    '''
    with open("header.vcf", "w") as header_vcf:
        # TODO: dynamic type, change INFO_FIELDS dict above
        for id, description in INFO_FIELDS.items():
            info_line = f'##INFO=<ID={id},Number=1,Type=String,Description="{description}">\n'
            header_vcf.write(info_line)
        
        if genome_build == "GRCh37":
            header_vcf.write(GRCh37_CONTIG)
        else:
            header_vcf.write(GRCh38_CONTIG)


def index_annotations(aggregated_database):
    '''
    Index the file with aggregated data
    
    Parameters
    ----------
    aggregated_database : str
        Output filename of aggregated data
    '''
    # TODO: rename sorted_ file better
    commands = f"""
    tail -n +2 {aggregated_database} > sorted_{aggregated_database}
    bgzip sorted_{aggregated_database}
    tabix -s1 -b2 -e2 sorted_{aggregated_database}.gz
    """
    try:
        subprocess.run(commands,
                       shell=True,
                       executable="/bin/bash",
                       capture_output=True,
                       check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with a non-zero exit code: {e.returncode}")
        print("Error:", e.stderr)


def bcftools_annotate_vcf(aggregated_database, minimal_vcf, output_file):
    '''
    Run bcftools annotate to annotate the minimal VCF with the aggregated info

    Parameters
    ----------
    aggregated_database : str
        Output filename of aggregated database
    minimal_vcf : str
        Output filename for the minimal VCF
    output_file : str
        Output filename for annotated VCF
    '''
    # TODO: allow output file saved to diff directory
    info_fields = ",".join(INFO_FIELDS.keys())
    bcftools_cmd = f"bcftools annotate -a sorted_{aggregated_database}.gz -h header.vcf -c CHROM,POS,REF,ALT,{info_fields} {minimal_vcf} > {output_file}"
    print(bcftools_cmd)
    try:
        subprocess.run(bcftools_cmd,
                       shell=True,
                       capture_output=True,
                       check=False)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with a non-zero exit code: {e.returncode}")
        print("Error:", e.stderr)


def main():
    args = parse_args()

    minimal_vcf = "minimal_vcf.vcf"
    aggregated_database = f"{'_'.join(args.probeset)}_aggregated_database.tsv"

    cleaned_csv = clean_csv(args.input_file)
    probeset_df = filter_probeset(cleaned_csv, args.probeset)
    aggregated_df = aggregate_uniq_vars(probeset_df, args.probeset, aggregated_database)

    # TODO: remove temp files?
    # TODO: should these functions be split in another file
    intialise_vcf(aggregated_df, minimal_vcf)
    write_vcf_header(args.genome_build)
    index_annotations(aggregated_database)
    bcftools_annotate_vcf(aggregated_database, minimal_vcf, args.output_file)


if __name__ == "__main__":
    main()
