nextflow.enable.dsl=2

process HaplotypeCaller {
    cpus 1
    memory { 4.GB * task.attempt }
    time { 12.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(meta), val(stem), path(crams), path(crais), val(region)

    output:
        tuple val(meta), val(region), path("${stem}_${region}.g.vcf.gz")

    script:
        def avail_mem = task.memory.toGiga()-1 //give gatk 1GB less mem than process
        def inputs = crams.collect{"-I ${it}"}.join(" ") //add "-I" before each cram
        """
        gatk --java-options -Xmx${avail_mem}g HaplotypeCaller \
            --reference $params.reference_fasta \
            --emit-ref-confidence GVCF \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -G StandardAnnotation -G StandardHCAnnotation \
            --output ${stem}_${region}.g.vcf.gz \
            --intervals $region \
            $inputs
        """

    stub:
        """
        touch ${stem}_${region}.g.vcf.gz
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
        gatk --java-options "-Xmx" GenomicsDBImport \
            --verbosity INFO \
            --intervals chr
            --genomicsdb-workspace-path \
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

    //Define multiple regions to run in parallel
    regions = Channel.of(
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22"
    )


    //STEP 1:
    //Group CRAM files by "sample_id" to run HaplotypeCaller on _part1 _part2 etc
    //
    //This looks ugly, but it's creating a grouping key as the donor_id+sample_id
    //then performing the grouping with groupTuple, afterwards it cleans up
    //the output by taking the first metadata map and removing the sample_num key
    //
    //finally get all combinations against the regions channel
    grouped_donor_samples = samples.map { m, cram, crai -> 
        tuple(m["donor_id"]+"_"+m["sample_id"], m, cram, crai) 
    }.groupTuple().map { k, ms, crams, crais ->
        def m = ms.first()
        m.remove("sample_num")
        tuple(m, k, crams, crais)
    }.combine(regions)

    donor_sample_gvcfs = grouped_donor_samples | HaplotypeCaller

    //STEP 2:
    //Call GenomicsDBImport over all samples for each donor/region combination
    //
    //
}

