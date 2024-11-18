nextflow.enable.dsl=2

process HaplotypeCaller {
    input:
        tuple val(meta), path(cram), path(crai)

    output:
        tuple val(meta), path("GVCF")

    script:
        """
        gatk --java-options -XX:-UsePerfData \
            HaplotypeCaller \
            --verbosity INFO \
            --reference ??? \
            --input $cram \
            --output GVCF
        """

    stub:
        """
        du -Lhs $cram > GVCF
        du -Lhs $crai >> GVCF

        # gatk is accessible through the singularity image
        which gatk >> GVCF

        # this dumps the HaplotypeCaller help message to file
        gatk --java-options -XX:-UsePerfData \
            HaplotypeCaller --help 2>> GVCF
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


    gvcfs = samples.take(10) | HaplotypeCaller //NOTE TAKE!
}

