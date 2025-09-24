from glob import glob
import pandas as pd

configfile: "config.yml"


DIR = "../data/"+config["dataset"]
xml = {xml_file: os.path.basename(xml_file) for xml_file in glob(DIR+"/*xml.gz")}

print(xml)

rule all:
    input:
        expand("../results/{dataset}/{dataset}_concat.tsv", dataset=config["dataset"]),
        expand("../results/{dataset}/{xml}_abstracts.tsv", dataset=config["dataset"], xml=xml.values())


rule parse_xml:
    input:
        "../data/{dataset}/{xml}"
    output: 
        "../results/{dataset}/{xml}_abstracts.tsv"
    shell : 
        "python xml_to_tsv.py -x {input} -o {output}"


rule aggregate:
    params:
        xml_files = list(xml.values()),
        dataset = config['dataset']
    input:
        expand("../results/{{dataset}}/{xml}_abstracts.tsv", xml=xml.values())
    output:
        "../results/{dataset}/{dataset}_concat.tsv"
    shell:
        "for el in {params.xml_files}; do cat ../results/{params.dataset}/$el\_abstracts.tsv | tail -n +2; done > {output}"


