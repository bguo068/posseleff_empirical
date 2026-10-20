"""
requires: bcftools (>1.2) and wget
"""

from subprocess import run
from pathlib import Path

url_prefix = "https://pf8-release.cog.sanger.ac.uk/vcf"
files = """
Pf3D7_01_v3.filt.vcf.gz
Pf3D7_02_v3.filt.vcf.gz
Pf3D7_03_v3.filt.vcf.gz
Pf3D7_04_v3.filt.vcf.gz
Pf3D7_05_v3.filt.vcf.gz
Pf3D7_06_v3.filt.vcf.gz
Pf3D7_07_v3.filt.vcf.gz
Pf3D7_08_v3.filt.vcf.gz
Pf3D7_09_v3.filt.vcf.gz
Pf3D7_10_v3.filt.vcf.gz
Pf3D7_11_v3.filt.vcf.gz
Pf3D7_12_v3.filt.vcf.gz
Pf3D7_13_v3.filt.vcf.gz
Pf3D7_14_v3.filt.vcf.gz
Pf3D7_API_v3.filt.vcf.gz
Pf3D7_MIT_v3.filt.vcf.gz
""".strip().split()

for vcf in files:
    chr = vcf.removesuffix(".filt.vcf.gz")
    bcf = f"{chr}.passonly.bcf"
    index = f"{chr}.passonly.bcf.csi"

    if Path(index).exists():
        continue

    cmd = f"""
    bcftools view -f PASS {url_prefix}/{vcf} -Ob --threads 10 -o {bcf}
    sleep 100
    bcftools index -f {bcf}
    sleep 100
    """
    print(cmd)
    ret = run(cmd, shell=True)
    assert ret.returncode == 0

meta_files = dict(
    annotation_file="https://pf8-release.cog.sanger.ac.uk/annotations/PlasmoDB-55_Pfalciparum3D7.gff.gz",
    ref_genome_file="https://pf8-release.cog.sanger.ac.uk/reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta",
    meta_table="https://pf8-release.cog.sanger.ac.uk/metadata/Pf8_samples.txt",
    fws_table="https://pf8-release.cog.sanger.ac.uk/Pf8_fws.tsv",
)

for _, url in meta_files.items():
    cmd = f"wget {url}"
    print(cmd)
    ret = run(cmd, shell=True)
    assert ret.returncode == 0
