% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crosshapviz_wrapper.R
\name{crosshap_viz}
\alias{crosshap_viz}
\title{Build 'crosshap' visualization}
\usage{
crosshap_viz(HapObject, plot_left = "allele")
}
\arguments{
\item{HapObject}{Haplotype object created by crosshap::run_haplotyping}
}
\value{
Returns a crosshap visualization
}
\description{
Builds five individual plots using various elements of a HapObject created by
crosshap::run_haplotyping. The central dotplot displays relationship between
clusters of linked SNPs (marker groups), and distinct haplotypes present
within the population. Vertical plots (top/bottom) visualize individuals and
populations, grouped by haplotype. Horizontal plots (left/right) visualize
SNP information, grouped by marker group cluster.
}