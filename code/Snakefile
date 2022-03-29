import os

FASTA_DIR = ""  # ENTER HERE YOUR FASTA FILES DIRECTORY PATH
OUTPUT_DIR = "" # ENTER HERE PATH TO OUTPUT FOLDER (pompan directory will be created)

FASTA_DIR = os.path.abspath(FASTA_DIR)
OUTPUT_DIR = os.path.abspath(OUTPUT_DIR)

if not os.path.exists(OUTPUT_DIR):
    os.mkdirs(OUTPUT_DIR)

WORKING_DIR = os.path.join(OUTPUT_DIR, "pompan")
PROKKA_DIR = os.path.join(WORKING_DIR, "prokka")
PROTEINORTHO_DIR = os.path.join(WORKING_DIR, "protheinortho")
PROT_CLUSTERS_DIR = os.path.join(WORKING_DIR, "prot_clusters")
CDS_CLUSTERS_DIR = os.path.join(WORKING_DIR, "cds_clusters")

SAMPLES = [os.path.splitext(x)[0] for x in os.listdir(FASTA_DIR)]
PROKKA_EXT = ["gff", "faa"]

print(expand("{prokka_dir}/prokka_{sample}/{sample}.faa",sample=SAMPLES,prokka_dir=PROKKA_DIR))

rule all:
    input:
        expand(["{prokka_dir}/prokka_{sample}/{sample}.{ext}"],sample=SAMPLES,ext=PROKKA_EXT, prokka_dir=PROKKA_DIR),
        "%s/single_copy_orthologs.tsv" % PROTEINORTHO_DIR

rule anno:
    input:
        "%s/{sample}.fasta" % FASTA_DIR
    conda:
        "envs/prokka.yaml"
    threads: workflow.cores
    output:
        annotation_faa = "%s/prokka_{sample}/{sample}.faa" % PROKKA_DIR,
        annotation_gff = "%s/prokka_{sample}/{sample}.gff" % PROKKA_DIR
    params:
        dir = directory("prokka/prokka_{sample}/"),
    shell:
        """
        prokka \
            --force \
            --cpus {threads} \
            --outdir {params.dir} \
            --prefix {wildcards.sample} \
            --centre X --compliant {input}
        """

rule proteinortho:
    input:
        expand("{prokka_dir}/prokka_{sample}/{sample}.faa",sample=SAMPLES, prokka_dir=PROKKA_DIR)
    conda:
        "envs/proteinotho.yaml"
    threads:
        workflow.cores
    output:
        "%s/myproject.proteinortho.tsv" % PROTEINORTHO_DIR
    params:
        protortho_dir = f"{PROTEINORTHO_DIR}"
    shell:
        """
        cd {params.protortho_dir}
        
        proteinortho6.pl -cpus={threads} {input}
        """

rule extract_single_copy_orthologs:
    input:
        rules.proteinortho.output
    output:
        "%s/single_copy_orthologs.tsv" % PROTEINORTHO_DIR
    shell:
        """
        cd {PROTEINORTHO_DIR}
        
        awk -F'\t' '{{if ($2=='3') if ($1=='3') print $0}}' myproject.proteinortho.tsv > single_copy_orthologs.tsv
        """

i = 0
CLUSTERS_NUMBER = []
with open(rules.extract_single_copy_orthologs.output) as fh:
    for line in fh:
        if not line.startswith("#"):
            i += 1
            CLUSTERS_NUMBER.append(i)

rule clusters_extracter:
    input:
        fastas = expand("%s/{sample}.fasta" % FASTA_DIR, sample=SAMPLES),
        proteins = expand("%s/prokka_{sample}/{sample}.faa" % PROKKA_DIR, sample=SAMPLES),
        annotations = expand("%s/prokka_{sample}/{sample}.gff" % PROKKA_DIR, sample=SAMPLES),
        orthologs = rules.extract_single_copy_orthologs.output
    output:
        prot_clusters = expand("%s/prot_cluster_{cluster_number}.mfa" % PROT_CLUSTERS_DIR, clusters_number=CLUSTERS_NUMBER),
        cds_clusters = expand("%s/cds_cluster_{cluster_number}.mfa" % CDS_CLUSTERS_DIR, clusters_number=CLUSTERS_NUMBER)
    params:
        cds_dir = directory(CDS_CLUSTERS_DIR),
        prot_dir = directory(PROT_CLUSTERS_DIR)
    shell:
        """
        clusters_extracter.py -f {input.fastas} -p {input.proteins} -a {input.annotations} -o {input.orthologs} -w {WORKING_DIR}
        """
