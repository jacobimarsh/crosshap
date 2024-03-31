── R CMD check results ───────────────────────── crosshap 1.4.0 ────
Duration: 1m 12.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded

--

Submitting updated versions (1.2.2 > 1.4.0) with bug fixes and minor improvements.
 
Implementing revisions identified as necessary during peer review (1.2.2 > 1.3.0) and identified bugs/necessary warnings/error messages following release (1.3.0 > 1.4.0) 

This version is expected to be stable until the next major advancement.

* Changed bar plot borders and NA colour
  * Updated associated documentation, unit tests and vignette
  
* Removed duplicated text and minor grammatical inconsistencies in vignette

* Updated README with new email and link to documentation

* Fixed bugs for plots and warnings for data with large numbers of heterozygous sites

* Added specific errors when incorrect inputs provided or non-viable parameters chosen for data
