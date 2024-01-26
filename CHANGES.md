# Changelog
## 2.6 (2024-01-26)
* Support for the POD5 format.
* Handling of barcoding in BAM and FASTQ files.
* Support for Fast5 .tar files without compression.
* Fixed an issue with BAM tags extraction.
* Fixed the issue with the 'passes_filtering'-only sequencing summary.
* Fixed the bug related to numpy.bool with numpy version 1.24 or later.
* Improved compatibility with pandas 2.0.

## 2.5 (2023-11-03)
* Fixed error when no failed reads were found (Issue #20).
* Fixed error when unclassified barcodes were missing.
* FASTQ and BAM files can now be used instead of the sequencing summary file.
* Added the ability to specify a barcode range (e.g., `--barcodes barcode01:barcode48`).

## 2.4 (2023-04-26)
* In over time graphs (read length, PHRED score and translocation speed), now fill the gaps for the 75% and 25% to avoid filling glitch.
* Fix 2D density plot title style.
* Fix error when a summary file with barcode information was provided in addition of barcoding files (Issue #17), now the barcode files will be skipped with a warning message in case summary file with barcode information is provided.
* Add the selected speed and sample frequency of the run in the "Run statistics" table of the ToulligQC report

## 2.3 (2023-03-22)
* Numpy 1.24 is now supported (thanks to Sean Black).
* Scatter plot of read length vs PHRED score has been replaced by a 2D density plot.
* Add bases per barcode distribution graphs.

## 2.2.3 (2022-09-29)
* Fix error when no Fast5 file is found in a directory provided as argument. Now throw an understandable error message.

## 2.2.2 (2022-08-31)
* Fix when multiple sequencing summary barcode files were available, the type of the 'barcode_arrangement' column in the dataframe was not correct.

## 2.2.1 (2022-01-05)
* Generated images were not included in the main HTML report file.

## 2.2 (2022-01-03)
* Add some flexibility to barcode specification (Thanks to  Hunter Cameron).
* The sequencing_summary.txt and sequencing_telemetry.js files can now be read compressed in gzip or bzip2.
* Add some log information on stdout for duration of the sequencing summary extractors.
* Change the format of the duration in log.
* Now logs memory used by dataframes.
* Fix Docker image build issue when updating setuptools.
* QScores and durations are now stored in 32 bits floats instead of 64 bits to reduce memory consumption (≈25% for 1D data).
* Barcode arrangements are now stored as categories instead of strings.

## 2.1.1 (2021-08-18)
* Fix issue when barcode list argument contains non existing barcode(s) in input data or when all existing barcodes are used.

## 2.1 (2021-06-28)
* The channel occupancy of the flowcell graph code has been rewritten to use Plotly. Add all/pass/fail/fail ratio views. The flowcell graph can now also handle Flongle and PromethION flowcels in addition of standard MinION flowcell
* Add "Sequencing kit" and "Barcode kits" entries in the run statistics table in html report
* In the distribution of read lengths graph, add buttons to show base count distribution in linear/log modes.
* Fix scatterplot graph where the default max x-axis value was always the max value for fail reads
* Fix help links in demo report
* Fix the name of the "Device and software" and "Run statistics" sections

## 2.0.1 (2021-04-14)
* In setup.py, set the developement status for ToulligQC as Production/Stable instead of Beta
* Add MANIFEST.in file to add resources files in PyPi package
* Fix error with latest versions of NumPy by add a missing int casting

## 2.0 (2021-04-09)
* Fix duration computation
* For PHRED score distribution boxplots, remove unnecessary interpolation before creating boxplot
* Remove duplicated code for 1D and 1D2 in PHRED score distribution graphs
* Remove duplicated code for 1D and 1D2 in Correlation between read length and PHRED score graphs
* In read length distribution graphs, add buttons to switch between linear and log scale for xaxis
* Add minimal qscore threshold in the "Device and software information" report table
* In 1D/1D2 sequencing summary extractors, now replace NA values for barcode assignment by "unclassified". Print a warning message on console
* Update unit tests
* For read count histogram tables, replace "frequency" by "percent"
* Rename y-axis for "Distribution of read lengths" graphs to "Read count"
* Add new command line options to finely define output file paths 
* Update the sigma value for gaussian filters when smoothing plots
* In correlation scatterplots, now ponderate the number of pass/fail spots by the pass/fail ratio when using interpolation

## 2.0b3 (2021-03-22)
* New CSS for HTML report
* Add new plots (Read length and PHRED over time, translocation speed...)
* Enhancement of existing graphs
* Big refactoring code for sequencing summary file parsing
* Big refactoring code for creating plots
* Reduce memory usage and execution time with barcodes
* Fix Plotly dependency version requirement
* Add L50 computation
* Sequence lengths of reads was stored into np.int16 that cannot handle >=32kb reads. Now use np.uint32
* Add ToulligQC logo in HTML report
* A telemetry file or Fast5 file is no more required
* Add new fields in the two first tables of the report: Run ID, operating system and basecalling date
* Barcode distribution pie charts can now be visualised as histograms
* In table, float values have now comma separator for thousands
* Update the yield number format in run statistocs table
* Update colors in the graphs
* Add an information link in all the graph titles

## 2.0b2 (2020-11-20)
* Fix import bug
* Fix graph names partially hidden in HTML summary element
* Rewrite help and rename arguments for clarity
* Create required and optional argument groups
* Create default values for --report-name and --output command line arguments
* Update report.html example in Docs with the new version of ToulliQC
* Create new presentation image for README

## 2.0b1 (2020-11-17)
* Refactoring of the sequencing_summary_extractor
* Refactoring of the 1dsqr_sequencing_summary_extractor
* Many performance improvements (reducing memory usage)
* Graphs are now made with Plotly
* Removal of unused options (Albacore log, FASTQ files, configuration file and samplesheet file)
* N50 information added to report.html
* Removal of Albacore support
* Now handle PromethION data
* Update of required dependencies versions
* Add unit tests
* Add new plots (throughput sequencing time)
* Update graph colors

## 1.3 (2019-11-07)
* Add a --barcodes option that allow to avoid samplesheet file creation
* The size of the graphs are now set to 1000x600px
* Many small fixes in graph generation (remove titles, fix grids and layouts...)
* In the HTML report, replace the tooltip icon by an unicode character

## 1.2 (2019-07-25)
* MultiFast5 file can now be used to retrieve run information
* Reporting other barcodes in the "other" category
* Gathering information from Telemetry files in the HTML report

## 1.1 (2019-03-21)
* Add Guppy support for 1D and 1D2
* Telemetry files generated by Albacore or Guppy can now be used to retrieve run information instead of reading a FAST5 file and the pipeline.log file.
* Refactoring of the code of the extractors

## 1.0 (2018-10-23)
* Report.data log file reviewed

## 0.10 (2018-07-18)
* Add pipeline.log parsing option

## 0.9 (2018-03-21)
* Fix out of memory error when parsing big FASTQ files. The parsing of FASTQ files is now faster


## 0.8 (2018-03-14)
* Fix unexisting import in toulligqc.py
* Fix the not working "--quiet" option


## 0.7 (2018-03-14)
* Fix Dockerfile that used the Ubuntu 17.04 (Ubuntu 17.04 packages repository is no more available)


## 0.6 (2018-03-12)
* Update html.report for 1D and 1Dsquare data
* Fix issue when processing fast5 files directory
* Add pass/fail filter
* Add extractor and graphs for 1dsquare analysis


## 0.5 (2017-11-28)
* Fix exception when toulligqc was launched with no arguments
* Remove pypandoc dependency in setup.py
* Fix issue when checking if directory paths ends with a '/'
* Fix issue when checking missing arguments


## 0.4 (2017-11-27)
* Fix issue with the --version option of ToulligQC


## 0.3 (2017-11-27)
* Fix issue with setup.py and pip install


## 0.2 (2017-11-27)
* ToulligQC can now handle Albacore 2.0 output
* The run date is now extracted from a FAST5 file
* Update setup.cfg for PyPi package submission
* Update ToulligQC documentation


## 0.1 (2017-08-30)
* Initial version
