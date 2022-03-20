#

This code contains routines to easily extract acoustic features from audio files inside a specific directory.
I coded in R, with well-known packages, and created a wrapper to invoke the routines with Python.
It works with **mp3**, **wav**, and **flac** file formats.
Also, some of my researches used this R code and its versions to study natural sound patterns, such as [[1]](https://doi.org/10.1016/j.ecoinf.2020.101184), [[2]](https://doi.org/10.1007/s00521-021-06501-w), and [[3]](https://repositorio.usp.br/item/002997058).

## Features

Available features and their references are grouped by their respective R packages.

### [soundecology](https://ljvillanueva.github.io/soundecology/)

* Acoustic Complexity Index (ACI) [[ref]](https://doi.org/10.1016/j.ecolind.2010.11.005)

* Acoustic Diversity Index (ADI) [[ref]](https://doi.org/10.1007/s10980-012-9806-4)

* Acoustic Evenness Index (AEI) [[ref]](https://doi.org/10.1007/s10980-011-9636-9)

* Bioacoustic Index (BIO) [[ref]]( https://doi.org/10.1890/07-0004.1)

* Normalized Difference Soundscape Index (NDSI) [[ref]](https://doi.org/10.1016/j.ecoinf.2012.08.001)

### [seewave](https://rug.mnhn.fr/seewave/)

* Acoustic Entropy Index (H), Temporal Entropy (Ht), and Frequency Entropy (Hf) [[ref]](https://doi.org/10.1371/journal.pone.0004065)

* Acoustic Richness (AR) and Median of Amplitude Envelope (M) [[ref]](https://doi.org/10.1016/j.ecolind.2011.05.006)

* Number of Peaks (NP) [[ref]](https://doi.org/10.1371/journal.pone.0065311)

* Roughness [[ref]](https://doi.org/10.1146/annurev-statistics-041715-033624)

* Rugosity [[ref]](https://doi.org/10.5424%2Fsjar%2F2009074-1109)

* Zero-crossing rate (ZCR) [[ref]](http://recherche.ircam.fr/anasyn/peeters/ARTICLES/Peeters_2003_cuidadoaudiofeatures.pdf)

### [tuneR](https://cran.r-project.org/web/packages/tuneR/index.html)

* Mel-frequency Cepstrum Coefficients (MFCC) [[ref]](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.11.9216)

### others

* Background noise index (BGN) [[ref]](https://arxiv.org/abs/2201.02099)

* Power Spectral Density (PSD) [[ref]](https://doi.org/10.1109/TAU.1967.1161901)

* Signal-to-noise Ratio (SNR) [[ref]](https://doi.org/10.1016/j.ecolind.2016.12.018)

* Mean Square Pressure (MSP) and Sound Pressure Level (SPL) [[ref]](https://doi.org/10.1016/j.marpolbul.2016.02.055)

* Root Mean Square (RMS) [[ref]](https://doi.org/10.1016/j.dsr.2016.06.001)

The code also generate spectrogram images of the recordings.

## Prerequisites

Additionally to R and the previous packages, the user will need to install [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) to perform parallel processing and reduce time consuming.
Also, to use Python wrapper, you must install the package [rpy2](https://pypi.org/project/rpy2/).

## Main functions

These functions are available in R to extract features and generate spectrograms of recordings.

```R
#It takes a directory and process all recordings inside it.
process.dir <- function(
  source.path, #directory with recordings to be processed
  target.path, #directory where results will be saved (default = source.path)
  aquatic, #if the audio files are from underwater landscape (default = FALSE). Depending on the landscape, some routines has different parameters.
  generate, #define the output indices ("index") or spectrograms ("spec") (default = "index")
  slice.size, #how many slices a recording (in seconds) should be divided (default = 1). A slice must not be greater than 90s.
  batch.size, #save indices results in batchs of batch.size (default = 100)
  img.dim, #if generate = "spec", it defines the image size (default = c(1366, 758))
  palette #if generate = "spec", it defines the image color palette (default = spectro.colors)
)  
```  

```R
#It process a specific recording.
process.file <- function(
  path, #recording to be processed
  spec.path, #if generate = "spec", it defines where image will be saved (default = NULL, inside path)
  aquatic, #if the audio files are from underwater landscape (default = FALSE). Depending on the landscape, some routines has different parameters.
  generate, #define the output indices ("index") or spectrograms ("spec") (default = "index")
  slice.size, #how many slices a recording (in seconds) should be divided (default = 1). A slice must not be greater than 90s.
  img.dim, #if generate = "spec", it defines the image size (default = c(1366, 758))
  palette, #if generate = "spec", it defines the image color palette (default = spectro.colors)
  start.parallel #if TRUE, process file in parallel mode (default = FALSE)
)
```

Python code has the class **AcousticIndices** with the following methods that have the same meaning of previous R functions:

```Python
process_dir(
  source_path, 
  target_path = None, 
  aquatic = False, 
  slice_size = 1)
```

```Python
process_file(
  path, 
  aquatic = False, 
  slice_size = 1)
```

## Results

Functions generate a feature set with 35 values and one identifier for each defined *slice.size*.
Indices are represented by a single value and the functions also return NDSI components (anthrophony and biophony), the mean and standard deviation of the PSD, and 12 MFCCs.
If you have stereo recordings, functions return the mean of the two channels.

Functions return a matrix of features and generate a CSV file inside *target.path*.
Spectrograms are also generated inside a *target.path* subdir.

Function *process.dir* generates log files to track the process. If index calculation stops for any reason, just reinvoke that function and the calculation will be resumed due to log files.

## Examples

```R
source("indices.R")

#Process all recordings inside a directory and returns a matrix with the results
files = process.dir("/home/user/recordings")

#Process a specific recording and returns a matrix with the results
one_file = process.file("/home/user/recordings/file.wav")

#Divides recordings into 2-second clips before calculating indices 
process.dir("/home/user/recordings", slice.size = 2)

#Divides recording into 2-second clips before generating spectrograms
process.dir("/home/user/recordings", slice.size = 2, generate = "spec")
```  

```python
from indices import AcousticIndices

extractor = AcousticIndices()

files = extractor.process_dir("/home/user/recordings")

one_file = extractor.process_file("/home/user/recordings/file.wav")
```

## Contact

* [Fábio Felix Dias](https://scholar.google.com.br/citations?hl=pt-BR&user=uQ_qg2MAAAAJ) - e-mail: <f_diasfabio@usp.br>
