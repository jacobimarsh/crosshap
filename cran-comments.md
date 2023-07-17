## R CMD check results

0 errors | 0 warnings | 0 notes

--

Submitting updated version (1.0.0 > 1.2.2) 
with several essential bug fixes and 
a minor

* Improved algorithm to handle extra variable (`tagphenos.R`)
  * Changed labels in `mid_dotplot.R`
  * Fixed associated documentation (`tagphenos.R`, `run_haplotyping.R`, `mid_dotplot.R`) and vignettes (`Getting_started.R`)

* Removed unnecessary argument producing misleading results (`pseudohaps.R`, `run_haplotpying_wrapper.R`)
* Fixed bug with node colours on figure (`clustree_viz.R`)
  * Fixed associated documentation 

* Updated test snapshots to fit with new version

* Updated Readme to include CRAN install instructions
