from typing import Dict, List
import logging
from os.path import join

import pandas as pd
from snakemake.utils import validate

print('viral-ampliseq-assembly version 1.0.0')

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index)

def get_bam_file(wildcards):
    """Get BAM file path for a sample"""
    return samples.loc[wildcards.sample, ['bam_file']]

