#!/usr/bin/env python3

import argparse
import logging
import csv
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from sklearn.mixture import GaussianMixture

from scipy.stats import mannwhitneyu

from dataclasses import dataclass


# Set up argument parser
parser = argparse.ArgumentParser(
    prog="plot_stats.py",
    description="""Accepts a TSV of normalised breakpoint counts from
    INDUCE-Seq data and provides outputs to assist with interpretation.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

parser.add_argument('tsv',
                    help="""path to a TSV containing counts of normalised
                    breakpoints""")

parser.add_argument('-o', '--output-prefix', default="results",
                    help="""prefix to be used on all output files""")

@dataclass
class MWResult:
    statistic: float
    pvalue: float
    control_count: int
    treated_count: int
    ambiguous_count: int
    control_median: float
    treated_median: float
    median_difference: float

    def summary_lines(self):
        return [
            ("Control count", self.control_count),
            ("Treated count", self.treated_count),
            ("Ambiguous count", self.ambiguous_count),
            ("Control median", self.control_median),
            ("Treated median", self.treated_median),
            ("Median difference", self.median_difference),
            ("Mann-Whitney statistic", self.statistic),
            ("Mann-Whitney p-value", round(self.pvalue, 5))
        ]


def make_qc_plots(breaks) -> plt.figure:
    
    fig, axs = plt.subplots(1, 3, figsize=(10, 5))

    axs[0].boxplot(breaks)
    axs[0].set_title('Box Plot')
    axs[0].set_ylabel('Normalised breakpoint counts')
    axs[0].yaxis.grid(True)

    axs[1].hist(breaks)
    axs[1].set_title('Histogram')
    axs[1].set_xlabel('Normalised breakpoint counts')
    axs[1].set_ylabel('Observed samples')
    axs[1].yaxis.grid(True)

    axs[2].violinplot(breaks, showmedians=True)
    axs[2].set_title('Violin Plot')
    axs[2].set_xlabel('Sample density')
    axs[2].set_ylabel('Normalised breakpoint counts')
    axs[2].yaxis.grid(True)
    
    return fig


def cluster_samples(counts_tsv) -> pd.DataFrame:
    df = pd.read_csv(counts_tsv, sep="\t")

    breaks = df[['normalised_breaks']].to_numpy()

    gmm = GaussianMixture(
        n_components=2,
        reg_covar=0.05,
        random_state=42
    )
    gmm.fit(breaks)

    probs = gmm.predict_proba(breaks)
    df["cluster_prob_0"] = probs[:, 0]
    df["cluster_prob_1"] = probs[:, 1]

    means = gmm.means_.flatten()
    control_cluster = np.argmin(means)
    treated_cluster = np.argmax(means)

    df["p_control"] = probs[:, control_cluster]
    df["p_treated"] = probs[:, treated_cluster]

    df["max_prob"] = df[["p_control", "p_treated"]].max(axis=1)

    CONFIDENCE_THRESHOLD = 0.99

    df["group"] = np.where(
        df["max_prob"] < CONFIDENCE_THRESHOLD,
        "ambiguous",
        np.where(
            df["p_control"] > df["p_treated"],
            "control",
            "treated"
        )
    )

    return df


def run_mann_whitney(clustered_df) -> None:
    controls = clustered_df.loc[
        clustered_df["group"] == "control",
        "normalised_breaks"
    ]

    treated = clustered_df.loc[
        clustered_df["group"] == "treated",
        "normalised_breaks"
    ]

    ambiguous = clustered_df.loc[
        clustered_df["group"] == "ambiguous",
        "normalised_breaks"
    ]

    statistic, pvalue = mannwhitneyu(controls,
                                     treated,
                                     alternative="two-sided")

    control_median = controls.median()
    treated_median = treated.median()

    return MWResult(
        statistic = statistic,
        pvalue = pvalue,
        control_count = len(controls),
        treated_count = len(treated),
        ambiguous_count = len(ambiguous),
        control_median = control_median,
        treated_median = treated_median,
        median_difference = (treated_median - control_median)
    )


def report_results(qc_plots, clustering_results_df,
                   mann_whitney_result, output_prefix) -> None:
    
    qc_plots.savefig(f"{output_prefix}.qc_plots.png")

    clustering_results_df.to_csv(f"{output_prefix}.assignments.tsv", 
                                 sep="\t", float_format="%.6f")

    summary_lines = mann_whitney_result.summary_lines()

    with open(f"{output_prefix}.summary.txt", "w") as f:
        f.write("\n".join(f"{label}:\t{item}" for label, item in summary_lines) + "\n")


def main() -> None:
    args = parser.parse_args()

    counts_tsv = args.tsv
    output_prefix = args.output_prefix

    df = pd.read_csv(counts_tsv, sep="\t")
    breaks = df[['normalised_breaks']]

    qc_plots = make_qc_plots(breaks)

    clustering_results_df = cluster_samples(counts_tsv)

    mann_whitney_result = run_mann_whitney(clustering_results_df)

    report_results(qc_plots, clustering_results_df, mann_whitney_result, output_prefix)


if __name__ == "__main__":
    main()