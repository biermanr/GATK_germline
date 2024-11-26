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
        tuple val(meta), val(region), path("${stem}_${region}.g.vcf.gz"), path("${stem}_${region}.g.vcf.gz.tbi")

    script:
        def avail_mem = task.memory.toGiga()-1 //give gatk 1GB less mem than process
        def inputs = crams.collect{"-I ${it}"}.join(" ") //add "-I" before each cram
        """
        gatk --java-options -Xmx${avail_mem}g HaplotypeCaller \
            --reference $params.reference_fasta \
            --emit-ref-confidence GVCF \
            --create-output-variant-index \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -G StandardAnnotation -G AS_StandardAnnotation \
            --output ${stem}_${region}.g.vcf.gz \
            --intervals $region \
            $inputs
        """

    stub:
        """
        touch ${stem}_${region}.g.vcf.gz ${stem}_${region}.g.vcf.gz.tbi
        """
}

process CombineGVCFs {
    cpus 1
    memory { 4.GB * task.attempt }
    time { 12.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(meta), val(donor_id), path(gvcfs), path(gvcf_indices)

    output:
        tuple val(meta), val(donor_id), path("${donor_id}.g.vcf.gz"), path("${donor_id}.g.vcf.gz.tbi")

    script:
        def avail_mem = task.memory.toGiga()-1 //give gatk 1GB less mem than process
        def inputs = gvcfs.collect{"-V ${it}"}.join(" ") //add "-V" before each gvcf
        """
        gatk --java-options -Xmx${avail_mem}g CombineGVCFs \
            --output ${donor_id}.g.vcf.gz \
            --reference $params.reference_fasta \
            --create-output-variant-index \
            -G StandardAnnotation -G AS_StandardAnnotation \
            $inputs
        """

    stub:
        """
        touch "${donor_id}.g.vcf.gz ${donor_id}.g.vcf.gz.tbi"
        """
}


process GenotypeGVCFs {
    cpus 1
    memory { 8.GB * task.attempt }
    time { 12.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(meta), val(donor_id), path(combined_gvcf), path(combined_gvcf_index)

    output:
        tuple val(meta), val(donor_id), path("called_${donor_id}.vcf.gz")

    script:
        def avail_mem = task.memory.toGiga()-1 //give gatk 1GB less mem than process
        """
        gatk --java-options -Xmx${avail_mem}g GenotypeGVCFs \
            --reference $params.reference_fasta \
            --output called_${donor_id}.vcf.gz \
            --all-sites \
            -V $combined_gvcf
        """

    stub:
        """
        touch called_${donor_id}.vcf.gz
        """
}

//START
workflow {
    params.reference_fasta = file("/projects/AKEY/akey_vol2/GTExSomaticMutations/Broad160XWGS.nobackup/Homo_sapiens_assembly38.fasta")

    //Read in sample info and CRAMs
    //into a tuple of [metadata, cram, crai]
    samples = Channel
        .fromPath("two_donor_samplesheet.csv")
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
        )).take(40) //NOTE TAKE!!!!

    //Define multiple regions to run in parallel
    regions = Channel.of(
        "chr22",
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
    //Call CombineGVCFs over all samples for each donor, over all regions
    //
    //Takes as input the donor/sample/region gvcfs and groupTuple's by donor
    //Deletes all other meta fields, to create a new channel `grouped_donors`
    //which are then piped into CombineGVCFs
    grouped_donors = donor_sample_gvcfs.map { m, region, gvcf, gvcf_index -> 
        tuple(m["donor_id"], m, gvcf, gvcf_index) 
    }.groupTuple().map { donor_id, ms, gvcfs, gvcf_indices ->
        def m = ms.first()
        m.remove("sample_id")
        m.remove("collaborator_sample_id")
        m.remove("tissue_site")
        m.remove("tissue_site_detail")
        tuple(m, donor_id, gvcfs, gvcf_indices)
    }

    combined_gvcfs = grouped_donors | CombineGVCFs

    //STEP 3:
    //Call GenotypeGVCFs on the combined_gvcfs one-to-one
    //
    called_gvcfs = combined_gvcfs | GenotypeGVCFs

    called_gvcfs.view()
}

