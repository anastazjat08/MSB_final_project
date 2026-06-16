import glob
import os
import re

configfile: "config/config.yaml"

GENOMES_DIR = config["genomes_dir"]

FASTAS = glob.glob(os.path.join(GENOMES_DIR, "*.fasta"))
SAMPLES = [os.path.basename(f).replace(".fasta", "") for f in FASTAS]

LABELS = config["muttui"]["labels"]

SIG_LIST = config["signatures"]["list"]


# =====================
# MUST BE
# =====================
rule all:
    input:
        # AMR report and analysis
        config["amr_analysis"]["report"],
        config["amr_analysis"]["output"],
        

        "1_data/2_processed/parsnp/parsnp.xmfa",
        "1_data/2_processed/parsnp/parsnp.fasta",
        "1_data/2_processed/gubbins/gubbins.final_tree.tre",
        "1_data/2_processed/gubbins/gubbins.filtered_polymorphic_sites.fasta",
        "1_data/2_processed/gubbins/gubbins.summary_of_snp_distribution.vcf",
        "1_data/3_results/iqtree/iqtree_final.treefile",


        # MutTui outputs
        config["muttui"]["conversion_path"],
        config["muttui"]["single_spectrum_dir"],

        # MutTui labels
        expand(config["muttui"]["labels_dir"] + "/labels_{label}.csv", label=LABELS),
        expand(config["muttui"]["split_out_dir"] + "/{label}/done.txt", label=LABELS),

        # SUPER-TABLES (AMR, METADATA, SPECTRA)
        expand(
            config["super_tables"]["dir"] + "/{label}/{label}_{type}.csv",
            label=config["muttui"]["labels"],
            type=config["super_tables"]["types"]
        ),

        # SUMMARY OF SUPER TABLES FOR STUDY, ST
        config["sbs_stacked_plots_dir"] + "/sbs_study.png",
        config["sbs_stacked_plots_dir"] + "/sbs_st.png",


        # COSINE SIMILARITY
        config["cosine_similarity_dir"] + "/MDR_level/cosine_distances.csv",
        config["cosine_similarity_dir"] + "/MDR_level/cosine_similarity.csv",
        config["cosine_similarity_dir"] + "/MDR_level/SBS_PCA_coordinates.csv",
        config["cosine_similarity_dir"] + "/MDR_level/SBS_PCA.pdf",

        config["cosine_similarity_dir"] + "/amr_burden/cosine_distances.csv",
        config["cosine_similarity_dir"] + "/amr_burden/cosine_similarity.csv",
        config["cosine_similarity_dir"] + "/amr_burden/SBS_PCA_coordinates.csv",
        config["cosine_similarity_dir"] + "/amr_burden/SBS_PCA.pdf",


        # PCA on 96 channels spectra
        config["PCA_96"]["dir"] + "/pca_colored.png",
        config["PCA_96"]["dir"] + "/supertable_with_pca.csv",
        config["PCA_96"]["dir"] + "/pca_explained_variance.csv",

        # LINEAR MODELS
        config["PCA_96"]["dir"] + "/models/PC1_linear_models.txt",
        config["signatures"]["models_dir"] + "/signature_models.txt",


        expand(
            config["signatures"]["plots_dir"] + "/{sig}_vs_amr.png",
            sig=SIG_LIST),
        expand(
            config["signatures"]["plots_dir"] + "/{sig}_boxplot.png",
            sig=SIG_LIST)



# =====================
# ANTIMICROBIAL REPORT
# =====================
rule amrfinder:
    input:
        os.path.join(GENOMES_DIR, "{sample}.fasta")
    output:
        temp("1_data/2_processed/amrfinder/{sample}.tsv")
    params:
        organism = config["organism"]
    shell:
        """
        echo "Processing sample: {wildcards.sample}" >> logs/amrfinder.log

        amrfinder --nucleotide {input} \
                  --organism {params.organism} \
                  --plus \
                  --name {wildcards.sample} \
                  > {output} \
                    2>> logs/amrfinder.log
        """

rule merge_amrfinder:
    input:
        expand("1_data/2_processed/amrfinder/{sample}.tsv", sample=SAMPLES)
    output:
        "1_data/3_results/amrfinder_result.tsv"
    shell:
        """
        echo "Merging AMRFinderPlus results into final TSV"

        # With first file add headers
        head -n 1 {input[0]} > {output}

        # Append data from all files, skipping headers
        tail -q -n +2 {input} >> {output}
        """

# =====================
# AMR ANALYSIS
# =====================
rule amr_analysis:
    input:
        results = config["amr_analysis"]["report"],
        metadata = config["amr_analysis"]["metadata"]
    output:
        directory(config["amr_analysis"]["output"])
    shell:
        """
        python 2_scripts/amr_analysis.py \
            -r {input.results} \
            -m {input.metadata} \
            -o {output}
        """


# =====================
# EXTRACTING CORE GENOME ALIGNMENT
# =====================
rule parsnp:
    input:
        ref = config["reference"]
    output:
        aln = "1_data/2_processed/parsnp/parsnp.xmfa"
    params:
        outdir = "1_data/2_processed/parsnp_run"
    shell:
        """
        echo "Running Parsnp with reference {input.ref}"

        parsnp -c \
               -r {input.ref} \
               -d {GENOMES_DIR} \
               -o {params.outdir}

        cp {params.outdir}/parsnp.xmfa {output.aln}
        """

rule xmfa_to_fasta:
    input:
        "1_data/2_processed/parsnp/parsnp.xmfa"
    output:
        "1_data/2_processed/parsnp/parsnp.fasta"
    shell:
        """
        harvesttools -x {input} -M {output}
        """

# =====================
# REMOVING RECOMBINATIONS
# =====================
rule gubbins:
    input:
        aln = "1_data/2_processed/parsnp/parsnp.fasta"
    output:
        tree = "1_data/2_processed/gubbins/gubbins.final_tree.tre",
        aln = "1_data/2_processed/gubbins/gubbins.filtered_polymorphic_sites.fasta"
    params:
        outdir = "1_data/2_processed/gubbins"
    shell:
        """
        run_gubbins.py \
            -p {params.outdir}/gubbins \
            {input.aln}
        """

# =====================
# BUILDING PHYLOGENETIC TREE
# =====================
rule iqtree:
    input:
        aln = "1_data/2_processed/gubbins/gubbins.filtered_polymorphic_sites.fasta"
    output:
        tree = "1_data/3_results/iqtree/iqtree_final.treefile"
    params:
        outdir = "1_data/2_processed/iqtree"
    shell:
        """
        mkdir -p {params.outdir}

        iqtree2 \
            -s {input.aln} \
            --seqtype DNA \
            -m MFP \
            -B 1000 \
            -T AUTO \
            --prefix {params.outdir}/iqtree

        mv {params.outdir}/iqtree.treefile {output.tree}
        """

# =====================
# CALCULATING MUTATIONAL SPECTRA
# =====================

# Convert VCF positions to MutTui format
rule muttui_position_conversion:
    input:
        vcf = "1_data/2_processed/gubbins/gubbins.summary_of_snp_distribution.vcf"
    output:
        pc_file ="1_data/2_processed/muttui/conversion.txt"
    shell:
        """
        mkdir -p $(dirname {output.pc_file})

        echo "Converting VCF positions for MutTui"
        MutTui convert-vcf \
            -v {input.vcf} \
            -o {output.pc_file}
        
        echo "Position conversion completed. Output at {output.pc_file}"
        """

# Calculate one spectrum for the whole tree (all branches)
rule muttui_single_spectrum:
    input:
        aln = "1_data/2_processed/gubbins/gubbins.filtered_polymorphic_sites.fasta",
        tree = "1_data/3_results/iqtree/iqtree_final.treefile",
        ref = config["reference"],
        conv = "1_data/2_processed/muttui/conversion.txt"
    params:
        outdir = "1_data/3_results/muttui/muttui_single_spectrum"
    output:
        directory("1_data/3_results/muttui/muttui_single_spectrum")
    shell:
        """
        mkdir -p {params.outdir}

        MutTui run \
            -a {input.aln} \
            -t {input.tree} \
            -r {input.ref} \
            -c {input.conv} \
            -o {params.outdir}
        """
# =====================
# GENERATING LABELS AND SPLITTING SPECTRA
# =====================
rule generate_labels:
    input:
        data=config["amr_analysis"]["report_metadata_mdr_burden"]
    output:
        config["muttui"]["labels_dir"] + "/labels_{label}.csv"
    params:
        label="{label}",
        outdir=config["muttui"]["labels_dir"]
    shell:
        """
        python3 2_scripts/generate_labels.py \
            -i {input.data} \
            -c {params.label} \
            -o {params.outdir}
        """

rule muttui_postprocess:
    input:
        mutations = config["muttui"]["single_spectrum_dir"] + "/all_included_mutations.csv",
        labels = config["muttui"]["labels_dir"] + "/labels_{label}.csv"
    output:
        touch(config["muttui"]["split_out_dir"] + "/{label}/done.txt")
    params:
        outdir = config["muttui"]["split_out_dir"] + "/{label}"
    shell:
        """
        mkdir -p {params.outdir}

        python3 -m MutTui.post_process_branch_mutations \
            -m {input.mutations} \
            -l {input.labels} \
            -o {params.outdir}

        touch {output}
        """

# =====================
# COMBINING AMR RESULTS/ METADATA WITH SPLITTED SPECTRA
# =====================
rule build_super_tables:
    input:
        spectra_dir=config["muttui"]["split_out_dir"] + "/{label}",
        metadata=config["amr_analysis"]["report_metadata_mdr_burden"]
    output:
        config["super_tables"]["dir"] + "/{label}/{label}_{type}.csv"
    params:
        outdir = config["super_tables"]["dir"],
        build_profiles = lambda wc: "--build_profiles" if wc.label == "Sanger_ID" else "",
        collapse_flag = lambda wc: "--collapse_to_six" if wc.type == "6" else ""
    shell:
        """
        python 2_scripts/build_super_tables.py \
            --spectra_dir {input.spectra_dir} \
            --amr_metadata {input.metadata} \
            --output_dir {params.outdir} \
            {params.build_profiles} \
            {params.collapse_flag}
        """

# =====================
# AMR AND SPECTRA ANALYSIS
# =====================
rule plot_sbs_stacked:
    input:
        study=config["super_tables"]["dir"] + "/Collection_name/Collection_name_6.csv",
        st=config["super_tables"]["dir"] + "/ST/ST_6.csv"
    output:
        study_plot=config["sbs_stacked_plots_dir"] + "/sbs_study.png",
        st_plot=config["sbs_stacked_plots_dir"] + "/sbs_st.png"
    shell:
        """
        python 2_scripts/plot_sbs_stacked.py --input {input.study} --output {output.study_plot}
        python 2_scripts/plot_sbs_stacked.py --input {input.st} --output {output.st_plot}
        """


rule cosine_similarity_mdr:
    input:
        spectra_dir=config["muttui"]["split_out_dir"] + "/MDR_level"
    output:
        dist=config["cosine_similarity_dir"] + "/MDR_level/cosine_distances.csv",
        sim=config["cosine_similarity_dir"] + "/MDR_level/cosine_similarity.csv",
        pca=config["cosine_similarity_dir"] + "/MDR_level/SBS_PCA_coordinates.csv",
        pdf=config["cosine_similarity_dir"] + "/MDR_level/SBS_PCA.pdf"
    run:
        spectra_dir=input.spectra_dir
        # list all files created by MutTui split
        files = os.listdir(spectra_dir)

        # keep only numeric MDR levels
        spectra_files = [
            os.path.join(spectra_dir, f)
            for f in files
            if re.match(r"mutational_spectrum_label_\d+\.csv$", f)
        ]

        # extract MDR levels
        MDR_LEVELS = [
            re.search(r"mutational_spectrum_label_(\d+)\.csv", f).group(1)
            for f in files
            if re.match(r"mutational_spectrum_label_\d+\.csv$", f)
        ]


        outdir = config["cosine_similarity_dir"] + "/MDR_level"
        spectra_str = " ".join(spectra_files)
        labels_str = " ".join([f"MDR{l}" for l in MDR_LEVELS])

        shell(f"""
            MutTui cluster \
                -s {spectra_str} \
                -l {labels_str} \
                -o {outdir}
        """)

rule cosine_similarity_burden:
    input:
        spectra_dir=config["muttui"]["split_out_dir"] + "/amr_burden"
    output:
        dist=config["cosine_similarity_dir"] + "/amr_burden/cosine_distances.csv",
        sim=config["cosine_similarity_dir"] + "/amr_burden/cosine_similarity.csv",
        pca=config["cosine_similarity_dir"] + "/amr_burden/SBS_PCA_coordinates.csv",
        pdf=config["cosine_similarity_dir"] + "/amr_burden/SBS_PCA.pdf"
    run:
        spectra_dir=input.spectra_dir
        # list all files created by MutTui split
        files = os.listdir(spectra_dir)

        # keep only numeric MDR levels
        spectra_files = [
            os.path.join(spectra_dir, f)
            for f in files
            if re.match(r"mutational_spectrum_label_\d+\.csv$", f)
        ]

        # extract MDR levels
        AMR_BURDEN = [
            re.search(r"mutational_spectrum_label_(\d+)\.csv", f).group(1)
            for f in files
            if re.match(r"mutational_spectrum_label_\d+\.csv$", f)
        ]


        outdir = config["cosine_similarity_dir"] + "/amr_burden"
        spectra_str = " ".join(spectra_files)
        labels_str = " ".join([f"AMR_BURDEN{l}" for l in AMR_BURDEN])

        shell(f"""
            MutTui cluster \
                -s {spectra_str} \
                -l {labels_str} \
                -o {outdir}
        """)
# PCA on all spectra for each genome (separate)
rule pca_96_channels:
    input:
        supertable=config["super_tables"]["dir"] + "/Sanger_ID/Sanger_ID_96.csv"
    output:
        pca_plot=config["PCA_96"]["dir"] + "/pca_colored.png",
        pca_table=config["PCA_96"]["dir"] + "/supertable_with_pca.csv",
        pca_exp=config["PCA_96"]["dir"] + "/pca_explained_variance.csv"
    params:
        mut_num=10,
        outdir=config["PCA_96"]["dir"]
    shell:
        """
        python 2_scripts/spectra_pca_96_channels.py \
            --supertable {input.supertable} \
            --mut_num {params.mut_num} \
            --output_dir {params.outdir}
        """

rule pc1_linear_models:
    input:
        pca_table=config["PCA_96"]["dir"] + "/supertable_with_pca.csv",
    output:
        models=config["PCA_96"]["dir"] + "/models/PC1_linear_models.txt"
    params:
        outdir=config["PCA_96"]["dir"] + "/models"
    shell:
        """
        python 2_scripts/spectra_PC1_linear_models.py \
            --pca_table {input.pca_table} \
            --output_dir {params.outdir}
        """

rule signature_linear_models:
    input:
        table6=config["super_tables"]["dir"] + "/Sanger_ID/Sanger_ID_6.csv"
    output:
        models=config["signatures"]["models_dir"] + "/signature_models.txt",
        plots=expand(
            config["signatures"]["plots_dir"] + "/{sig}_vs_amr.png",
            sig=SIG_LIST
        ) + expand(
            config["signatures"]["plots_dir"] + "/{sig}_boxplot.png",
            sig=SIG_LIST)
    params:
        outdir=config["signatures"]["dir"],
        mut_num=10,
        siglist=",".join(SIG_LIST)
    shell:
        """
        python 2_scripts/spectra_sig_linear_models.py \
            --table6 {input.table6} \
            --mut_num {params.mut_num} \
            --signatures {params.siglist} \
            --output_dir {params.outdir}
        """




