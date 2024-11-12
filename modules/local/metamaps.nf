process METAMAPS {
    tag "$meta.id"
    label 'process_high_memory'

    input:
    tuple val(meta), path(reads)
    path(ref)

    output:
    tuple val(meta), path('*.WIMP')                     , optional:true, emit: classification_results
    tuple val(meta), path('*species.EM.WIMP')           , optional:true, emit: species_results
    tuple val(meta), path('*reads2Taxon')           , optional:true, emit: reads_to_taxon
    path('*species.EM.WIMP')                            , optional:true, emit: species_results_mqc


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/metamaps/metamaps-1.0.sif metamaps mapDirectly --all -r ${ref}/DB.fa -q ${reads} -o ./"$prefix"_classification_results --maxmemory ${params.metamaps_mem} -t ${params.metamaps_threads}
    
    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/metamaps/metamaps-1.0.sif metamaps classify --mappings "$prefix"_classification_results --DB ${ref} -t ${params.metamaps_threads}

    header=\$(head -n 1 "$prefix"_classification_results.EM.WIMP)


    echo \$header > "$prefix"_classification_results_1_strain.EM.WIMP

    awk '\$1 == "definedGenomes"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_1_strain.EM.WIMP


    echo \$header > "$prefix"_classification_results_2_species.EM.WIMP

    awk '\$1 == "species"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_2_species.EM.WIMP
    """
    // echo \$header > "$prefix"_classification_results_3_genus.EM.WIMP

    // awk '\$1 == "genus"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_3_genus.EM.WIMP


    // echo \$header > "$prefix"_classification_results_4_family.EM.WIMP

    // awk '\$1 == "family"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_4_family.EM.WIMP


    // echo \$header > "$prefix"_classification_results_5_order.EM.WIMP

    // awk '\$1 == "order"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_5_order.EM.WIMP


    // echo \$header > "$prefix"_classification_results_6_class.EM.WIMP

    // awk '\$1 == "class"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_6_class.EM.WIMP


    // echo \$header > "$prefix"_classification_results_7_phylum.EM.WIMP

    // awk '\$1 == "phylum"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_7_phylum.EM.WIMP


    // echo \$header > "$prefix"_classification_results_8_kingdom.EM.WIMP

    // awk '\$1 == "kingdom"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_8_kingdom.EM.WIMP


    // echo \$header > "$prefix"_classification_results_9_domain.EM.WIMP

    // awk '\$1 == "domain"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_9_domain.EM.WIMP
    
}
