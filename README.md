# bin2fil-py
This is a python script I wrote to replace (or at least add an alternative) to the [old bin2fil](https://github.com/gio54321/bin2fil).

## Installation
To run this script you must have python (note this is only tested with python3) installed with `numpy` and `scipy` libraries. You can install these dependencied by
```
pip install numpy
pip install scipy
```
then you can clone the repository on to your PC
```
git clone https://github.com/gio54321/bin2fil.git
```

## Execution
As shown in the help output:

```
$python bin2fil.py -h
usage: bin2fil.py [-h] [-s] [-l] [-mjd] [-fh] [-fl] [-a] [-c] [-sr] [-f] [-n]
                  [-ra] [-de] [-e] [-o] [-cs]
                  input_file

description: bin2fil tool rewritten in python

positional arguments:
  input_file            input bin file

optional arguments:
  -h, --help            show this help message and exit
  -s , --start          offset in seconds to start the conversion
  -l , --length         length of the conversion in seconds (default is until
                        EOF)
  -mjd , --mjd-start    MJD start time (default is taken from creation time of
                        fhe input file)
  -fh , --filt-high     frequency of the highpass filter in Hz
  -fl , --filt-low      frequency of the lowpass filter in Hz
  -a , --ampli          amplification during conversion
  -c , --channels       number of channels
  -sr , --sample-rate   sample rate
  -f , --center-freq    center frequency
  -n , --source-name    source name
  -ra , --source-ra     source radial ascension
  -de , --source-de     source declination
  -e , --sr-error       error in ppm of the dongle
  -o , --output-file    output file name (default is replaced only the
                        extension)
  -cs , --chunk-size    size of the chunks that are elaborated
```

so the most basic command you can run is 
```
python bin2fil.py file.bin
```
Take in consideration that *if not specified the program will always take the default values*. See below if you want to change these default values.

## Customization
At the beginning of the file there are some global variables you have to modify if you want to run the script without always specify all the parameters. These variables are pretty self-explainatory and are correspondents to the CLI arguments.

```python
# -------------- GLOBAL DEFAULTS VATIABLES ---------------
conv_start = 0
conv_len = 7200
ampli = 3.0
sample_rate = 1000
f_h = 0.1
f_l = 60.0
channels = 25
center_freq = 415
sr_error = 3.0
mjd_time = -1 # -1 means take from modify date
source_name = 'B0329+54'
source_ra = 33259.37
source_de = 543443.57
chunk_size = 50000 # I think this is the best but feel free to experiment
# --------------------------------------------------------
```