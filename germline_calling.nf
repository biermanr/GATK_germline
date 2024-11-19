nextflow.enable.dsl=2

process HaplotypeCaller {
    cpus 1
    memory '10 GB'
    time '4h'

    input:
        tuple val(meta), path(cram), path(crai)

    output:
        tuple val(meta), path("output.g.vcf.gz")

    script:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
            --input $cram \
            --verbosity INFO \
            --reference $params.reference_fasta \
            --intervals chr22 \
            --emit-ref-confidence GVCF \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -G StandardAnnotation -G StandardHCAnnotation \
            --output output.g.vcf.gz
        """

    stub:
        """
        touch output.g.vcf.gz
        """
}

/*
process GenomicsDBImport {
    input:
        tuple val(meta), ???

    output:
        ???

    script:
        """
        """

    stub:
        """
        """
}


process GenotypeGVCFs {
    input:
        tuple val(meta), ???

    output:
        ???

    script:
        """
        """

    stub:
        """
        """
}
*/

//START
workflow {
    params.reference_fasta = file("/projects/AKEY/akey_vol2/GTExSomaticMutations/Broad160XWGS.nobackup/Homo_sapiens_assembly38.fasta")

    //Read in sample info and CRAMs
    //into a tuple of [metadata, cram, crai]
    samples = Channel
        .fromPath("/projects/AKEY/akey_vol2/GTExSomaticMutations/Broad160XWGS.nobackup/cram_samplesheet.csv")
        .splitCsv(header: true)
        .map(row -> tuple(
            [ 
                "sample_id":row["Sample ID"],
                "donor_id":row["Collaborator Participant ID"],
                "collaborator_sample_id":row["Collaborator Sample ID"],
                "tissue_site":row["Tissue Site"],
                "tissue_site_detail":row["Tissue Site Detail"],
                "sample_num":row["Sample Num"]
            ],
            file(row["cram_path"]), 
            file(row["crai_path"])
        ))


    gvcfs = samples.take(2) | HaplotypeCaller //NOTE TAKE!
}

