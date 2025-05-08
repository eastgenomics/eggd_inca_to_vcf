import pandas as pd
import argparse
import config
import pysam
import pysam.bcftools
import os.path


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
        required=True,
        help=(
            "CSV file exported from the previous interpretations "
            "database"
        ),
    )

    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        default="output.vcf",
        help=(
            "Output VCF filename"
        ),
    )

    parser.add_argument(
        "-g",
        "--genome_build",
        type=str,
        required=True,
        choices=["GRCh37", "GRCh38"],
        help="Genome build the samples were run in",
    )

    parser.add_argument(
        "-set",
        "--probeset",
        type=str,
        nargs='+',
        default="germline somatic",
        help=(
            "probset_id or allele_origin to filter"
        ),
    )

    parser.add_argument(
        "--keep_tmp",
        action='store_false',
        help=(
            "Bool to keep intermediate tsv annotation file"
        ),
    )

    parser.add_argument(
        '--output_dir', 
        default=os.getcwd(),
        help="path to output directory"
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
    interpreted_df = cleaned_csv[cleaned_csv['interpreted'].str.lower() == "yes"]
    if (interpreted_df['germline_classification'].isnull() & interpreted_df['oncogenicity_classification'].isnull()).any():
        raise ValueError("Both germline and oncogenicity classification are null in at least one row.")

    all_dfs = []
    for type in probeset:
        if type in ["99347387", "96527893"]:
            column = "probeset_id"
        elif type in ["germline", "somatic"]:
            column = "allele_origin"
        else:
            raise ValueError(f"Invalid argument: '{type}'. Expected one of: germline, somatic, 99347387, or 96527893.")
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
    aggregated_df.to_csv(aggregated_database, sep="\t", index=False, header=False)

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
        vcf_file.write(config.MINIMAL_VCF_HEADER)
        vcf_file.write("\n".join(vcf_lines) + "\n")


def write_vcf_header(genome_build, header_filename):
    '''
    Write VCF header by populating INFO fields and specifying contigs

    Parameters
    ----------
    genome_build : str
        Genome build to specify contigs in header
    '''
    with open(header_filename, "w") as header_vcf:
        for field_info in config.INFO_FIELDS.values():
            info_line = f'##INFO=<ID={field_info["id"]},Number={field_info["number"]},Type={field_info["type"]},Description="{field_info["description"]}">\n'
            header_vcf.write(info_line)
        
        if genome_build == "GRCh37":
            header_vcf.write(config.GRCh37_CONTIG)
        else:
            header_vcf.write(config.GRCh38_CONTIG)


def index_annotations(aggregated_database):
    '''
    Index the file with aggregated data
    
    Parameters
    ----------
    aggregated_database : str
        Output filename of aggregated data
    '''
    pysam.tabix_compress(f"{aggregated_database}", f"{aggregated_database}.gz")
    pysam.tabix_index(f"{aggregated_database}.gz", seq_col=0, start_col=1, end_col=1)


def bcftools_annotate_vcf(aggregated_database, minimal_vcf, output_file, output_dir):
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
    output_dir : str
        Output directory for annotated VCF
    '''
    # Set up output directory
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    # Run bcftools annotate with pysam
    info_fields = ",".join(item["id"] for item in config.INFO_FIELDS.values())
    annotate_output = pysam.bcftools.annotate("-a", f"{aggregated_database}.gz", "-h", "header.vcf", "-c", f"CHROM,POS,REF,ALT,{info_fields}", f"{minimal_vcf}")
    with open(output_path, 'w') as f:
        f.write(annotate_output)


def main():
    args = parse_args()

    minimal_vcf = "minimal_vcf.vcf"
    header_filename = "header.vcf"
    aggregated_database = f"{'_'.join(args.probeset)}_aggregated_database.tsv"

    cleaned_csv = clean_csv(args.input_file)
    probeset_df = filter_probeset(cleaned_csv, args.probeset)
    aggregated_df = aggregate_uniq_vars(probeset_df, args.probeset, aggregated_database)

    intialise_vcf(aggregated_df, minimal_vcf)
    write_vcf_header(args.genome_build, header_filename)
    index_annotations(aggregated_database)
    bcftools_annotate_vcf(aggregated_database, minimal_vcf, args.output_file, args.output_dir)

    # TODO: use keep_tmp arg in code.sh


if __name__ == "__main__":
    main()
