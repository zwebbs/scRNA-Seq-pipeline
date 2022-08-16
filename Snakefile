# File Name: Snakefile
# Created By: HKim
# Created On: 2022-07-19
# Purpose: Runs the scRNA-Seq analysis Pipeline


# Module Imports
# ----------------------------------------------------------------------------
from gardnersnake.core import Configuration
from pathlib import Path


# Global Configuration
# ----------------------------------------------------------------------------
yaml_config_filepath = Path(config["yaml_config"])
cfg = Configuration(filepath=yaml_config_filepath)
cfg.load()

GLOBALS = cfg.global_params
WORKDIR = Path(GLOBALS.working_directory)
PLOT_DIR = Path(GLOBALS.misc.plot_directory)
DATA_DIR = Path(GLOBALS.misc.data_dir)
REF = Path(GLOBALS.misc.reference_directory)
METADATA = (WORKDIR / Path(GLOBALS.files.sample_metadata))

# set analysis workdir
workdir: WORKDIR

# build a map of the sequencing data
SEQ = GLOBALS.files.sequencing
RUN_IDS = [s["run_id"] for s in SEQ]  # get the run ids for each sequencing obj
SEQ_MAP = {run_id:idx for idx,run_id in enumerate(RUN_IDS)}  # get their index


# Function Definitions
# ----------------------------------------------------------------------------
get_fastq_dir = lambda run_id: SEQ[SEQ_MAP[run_id]]["fastq_dir"]


# Rule 0. Pipeline Global Returns
# ----------------------------------------------------------------------------
rule All:
    input:
        expand(WORKDIR / "data/align/cellranger_count.{run_id}.rc.out", run_id=RUN_IDS),
        multiext("data/plots/HCA_Integration_clust_", "pca.pdf", "tsne.pdf", "umap.pdf"),
        multiext("data/plots/HCA_Integration_clust_", "cardio.pdf", "fib.pdf","endo.pdf", "peri.pdf")
        expand(PLOT_DIR / "{run_id}_scrublet_hist.png", run_id=RUN_IDS),
        expand(PLOT_DIR / "{run_id}_scrublet_umap.png", run_id=RUN_IDS),
        expand(DATA_DIR / "scrublet/{run_id}_scrublet_doublets.csv", run_id=RUN_IDS)

# Rule 1. Align and Quantify from FASTQ.
# ---------------------------------------------------------------------------
cellranger_rp = cfg.get_rule_params(rulename="CellRanger_FASTQ_to_counts")
rule CellRanger_FASTQ_to_counts:
    input: transcriptome = (REF / GLOBALS.files.transcriptome),
        fastq_dir = lambda wildcards: (WORKDIR / get_fastq_dir(run_id=f"{wildcards.run_id}"))
    output: cellranger_count_rc = (WORKDIR / "data/align/cellranger_count.{run_id}.rc.out")
    params: **(cellranger_rp.parameters), sample = lambda wildcards: f"{wildcards.run_id}"
    resources: **(cellranger_rp.resources), job_id = lambda wildcards: f"{wildcards.run_id}"
    shell:
        "mkdir -p data/align/ && cd data/align/ && "
        "cellranger count --id={params.sample}"
        " --transcriptome={input.transcriptome}"
        " --fastqs={input.fastq_dir}"
        " --sample={params.sample}"
        " --expect-cells={params.ncells}"
        " {params.extra_args}"
        " && check_directory -o {output.cellranger_count_rc}"
        " {params.checkfiles} {params.sample}/outs/"


# Rule 4. Doublet detection using scrublet
# ----------------------------------------------------------------------------
scrublet_detection_rp = cfg.get_rule_params(rulename="Scrublet_Doublet_Detect")
rule Scrublet_Doublet_Detect:
    input: rules.CellRanger_FASTQ_to_counts.output.cellranger_count_rc
    params:
        **(scrublet_detection_rp.parameters),
        plotdir = PLOT_DIR,
        outdir = DATA_DIR / "scrublet/",
        sample_name = lambda wildcards: f"{wildcards.run_id}"
    resources:
        **(scrublet_detection_rp.resources),
        job_id = lambda wildcards: f"{wildcards.run_id}"
    output:
        PLOT_DIR / "{run_id}_scrublet_hist.png",
        PLOT_DIR / "{run_id}_scrublet_umap.png",
        DATA_DIR / "scrublet/{run_id}_scrublet_doublets.csv"
    shell:
        "mkdir -p {params.plotdir} {params.outdir} && "
        "scripts/doublet_detection.py -r {params.doublet_prior_rate}"
        " --datadir {params.datadir} --sample-name {params.sample_name}"
        " --outdir {params.outdir} --plotdir {params.plotdir}"


# Rule 3. Integrate data and plot clusters and markers using Seurat
# ----------------------------------------------------------------------------
seurat_integration_rp = cfg.get_rule_params(rulename="Seurat_Integration")
rule Seurat_Integration:
    input: metadata = METADATA
    params: analysis_name = "HCA_Integration",
        datadir = "data/align/", plotdir = "data/plots/",
        rdata_path = "data/rdata/integrated_seurat.rds"
    resources: **(seurat_integration_rp.resources), job_id = "glob"
    output:
        multiext("data/plots/HCA_Integration_clust_", "pca.pdf", "tsne.pdf", "umap.pdf"),
        multiext("data/plots/HCA_Integration_clust_", "cardio.pdf", "fib.pdf","endo.pdf", "peri.pdf")
    shell:
        "mkdir -p {params.plotdir} && mkdir -p data/rdata/ &&"
        "Rscript scripts/integrateAndPlot.R -f {input.metadata}"
        " -d {params.datadir} -o {params.plotdir}"
        " -n {params.analysis_name} -r {params.rdata_path}"
