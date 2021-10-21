#!/usr/bin/env python3
"""
Create an interactive plotly visualisation of coverage values of defined 
regions from an augmented bedtools coverage output file.
"""

import argparse
import pandas as pd
import plotly.express as px
from pathlib import Path

def get_sample_name(f_path):
    '''Extract samplename from string e.g: /home/rada/Documents/TNscope/PVAL_65_S1/exon_cov/PVAL_65_S1.exon_cov.tsv'''
    return str(Path(f_path).stem).split(".")[0]

def create_dir_if_not_exist(dir_name):
    '''Create a directory based on input string'''
    if not Path(dir_name).is_dir():
        Path(dir_name).mkdir(parents=True, exist_ok=True)


def main(args):
    """Run the command line program."""
    sample_name = str(args.sample).replace("./","",1)
    in_f = str(args.infile).replace("./","",1)
    print(args.infile)
    df = pd.read_csv(args.infile, 
        sep="\t", 
        # Note this must have the same contents as in add_ids.awk's header_row variable
        names=["sample","name","chromosome","start","end","coverage","num_bases_at_depth","length","procent_of_coverage_in_region"],
        header=0)

    # Assign correct datatypes to each column
    df = df.astype({
        'sample':'str',
        'name':'str',
        'chromosome':'str',
        'start':'int',
        'end':'int',
        'coverage':'int',
        'num_bases_at_depth':'int',
        'length':'int',
        'procent_of_coverage_in_region':'float'
        })
    
    # Duplicate rows with several bases, this enables calculating averages in a more easier way
    df = df.loc[df.index.repeat(df.num_bases_at_depth)] # https://stackoverflow.com/a/57009491
    
    # Create temp directory for files for double checking things
    create_dir_if_not_exist("temp")

    # Print in order to double check that the previous step went OK
    df.to_csv("temp/chrs_" + sample_name + ".tumor_deduped_repeated_cov.csv")
    
    # Get a list of names so they can be looped through later on
    region = df['name'].unique()
    # Split the df into "name" subgroups
    grouped_names = df.groupby(df.name)
    # This will hold a list of df:s with depth metrics data such as 
    region_list = []
    print("Number of unique regions: " + str(len(region)))

    for r in region:
        c = grouped_names.get_group(r)
        # Scrape NCBI transcript ID
        name_id = c.iloc[0]['name']
        c.index.name = "row_no"
        # Extract metrics values for the depth column
        c = c.describe()['coverage'].to_frame(name_id).T
        region_list.append(c)

    # Join all metrics data into one df
    regions_metrics = (pd.concat(region_list, axis=0)
                            .rename(columns={'count': 'total_length_of_region'})
                            .astype({
                                    'total_length_of_region':'int',
                                    'max':'int',
                                    'min':'int'}
                                    ))
    regions_metrics.index.name = "ID"
    regions_metrics.to_csv("temp/chrs_" + sample_name + ".tumor_deduped_repeated_cov.metrics.csv")

    # Create bar plots
    bar_fig = px.bar(regions_metrics.reset_index(), 
        y='mean', 
        x='ID', 
        hover_data=["total_length_of_region", "mean", 'std', 'min', '25%', '50%', '75%', 'max',"ID"],
        title="Sample name: " + sample_name)

    # https://stackoverflow.com/a/59869358
    # Write all bar plots into one html page
    bar_fig.write_html(sample_name + ".html")

if __name__ == '__main__':
    """Parse arguments and run the main function"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', 
                        type=str,
                        help='Tsv input file containing coverages of regions of interest')
    parser.add_argument('sample',
                        type=str,
                        help='Name of the sample, this will be used in naming of the output files')
    args = parser.parse_args()
    main(args)
