process mapQFilterReads {
    input:
        tuple val(id), path(breaksBED)
        val mapQThreshold

    output:
        tuple val(id), path(breaksBED), path("${id}.filtered.bed")

    script:
    """
        mapq_filter_reads.py "${breaksBED}" "${id}.filtered.bed" -q "${mapQThreshold}" -l "${id}.mapq_filter_reads.log"
    """

}

process bedIntersectResctrictionSites {
    conda 'bioconda::bedtools==2.31'

    input:
        tuple val(id), path(originalBED), path(filteredBED)
        path referenceBED

    output:
        tuple val(id), path(originalBED), path("${id}.filtered.intersected.bed")

    script:
    """
    bedtools intersect \
        -a "${filteredBED}" \
        -b "${referenceBED}" \
        > "${id}.filtered.intersected.bed"
    """
}

process countNormaliseBreaks {
    input:
        tuple val(id), path(originalBED), path(intersectedBED)
        val normalisation_value

    output:
        path('*.counts.tsv'), emit: countsTSV

    script:
    """
        normalise_break_counts.py \
            "${originalBED}" \
            "${intersectedBED}" \
            -o "${id}.counts.tsv" \
            -n "${normalisation_value}"
    """
}

process collectStats {
    input:
        path (countsTSVs)

    output:
        path('all_samples.counts.tsv'), emit: collectedStatsTSV

    script:
    """
    python3 <<EOF
    import sys
    import csv

    tsvs = "$countsTSVs".split(" ")

    with open('all_samples.counts.tsv', 'w') as out:
        writer = csv.writer(out, delimiter="\\t", lineterminator="\\n")
        writer.writerow(["sample_name", "total_breaks", "canonical_breaks", "normalised_breaks"])
        
        for tsv in tsvs:
            with open(tsv, 'r') as input:
                reader = csv.reader(input, delimiter="\\t")
                for row in reader:
                    writer.writerow(row)
    EOF
    """
}

process plotStats {
    conda "$projectDir/env/plot-stats.yml"

    publishDir 'results', mode: 'symlink'

    input:
        path collectedStatsTSV
        val experimentName

    output:
        file('*.qc_plots.png')
        file('*.assignments.tsv')
        file('*.summary.txt')

    script:
    """
        plot_stats.py "$collectedStatsTSV" -o "$experimentName"
    """
}

workflow {

    beds_ch = Channel.fromPath("${params.input_dir}/*.bed")
                     .map( bed -> tuple(bed.baseName, bed))
    
    filtered_ch = mapQFilterReads(beds_ch, params.mapQThreshold)

    intersected_ch = bedIntersectResctrictionSites(filtered_ch, params.referenceBED)

    normalised_ch = countNormaliseBreaks(intersected_ch, params.normalisationValue)

    collectStats(normalised_ch.collect())

    plotStats(collectStats.out.collectedStatsTSV, params.experimentName)

}
