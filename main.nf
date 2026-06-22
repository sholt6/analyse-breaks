process mapQFilterReads {
    input:
        tuple val(id), path(breaksBED)

    output:
        tuple val(id), path(breaksBED), path("${id}.filtered.bed")

    script:
    """
        mapq_filter_reads.py "${breaksBED}" "${id}.filtered.bed" -q "${params.mapQThreshold}" -l "${id}.mapq_filter_reads.log"
    """

}

process bedIntersectRestrictionSites {
    conda 'bioconda::bedtools==2.31'

    input:
        tuple val(id), path(originalBED), path(filteredBED)

    output:
        tuple val(id), path(originalBED), path("${id}.filtered.intersected.bed")

    script:
    """
    bedtools intersect \
        -a "${filteredBED}" \
        -b "${params.referenceBED}" \
        > "${id}.filtered.intersected.bed"
    """
}

process countNormaliseBreaks {
    input:
        tuple val(id), path(originalBED), path(intersectedBED)

    output:
        path('*.counts.tsv'), emit: countsTSV

    script:
    """
        normalise_break_counts.py \
            "${originalBED}" \
            "${intersectedBED}" \
            -o "${id}.counts.tsv" \
            -n "${params.normalisationValue}"
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

    beds = Channel.fromPath("${params.input_dir}/*.bed")
                     .map( bed -> tuple(bed.baseName, bed))

    processed = beds
        | mapQFilterReads
        | bedIntersectRestrictionSites
        | countNormaliseBreaks

    stats = collectStats(processed.collect())
    plotStats(stats.collectedStatsTSV, params.experimentName)
}
