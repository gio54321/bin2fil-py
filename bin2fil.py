import os
import sys
import struct
import argparse
import time
import math

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

# -------------- GLOBAL DEFAULTS VATIABLES ---------------
conv_start = 0
conv_len = 7200
ampli = 1.0
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
show_saturation = False
# --------------------------------------------------------

filters = []
infile = ''
outfile = ''
running = True

class MyParser(argparse.ArgumentParser):
	# custom class to override the error function
	# of argparse to simply show the help message on error
	def error(self, message):
		sys.stderr.write('ERROR:%s\n' % message)
		self.print_help()
		sys.exit(2)


def parse_arguments(debug=True):
	# function supposed to parse the CLI arguments

	# do update the global vars instead of the local ones
	global conv_start, conv_len, ampli, sample_rate, f_h, f_l, channels, center_freq, sr_error
	global mjd_time, source_name, source_ra, source_de, infile, outfile, chunk_size, show_saturation

	# define the argument parser
	parser = MyParser(description='description: bin2fil tool rewritten in python')

	# start time
	parser.add_argument('-s', '--start', metavar='', action='store', type=int, default=conv_start, dest='conv_start',
		help='offset in seconds to start the conversion')

	# length
	parser.add_argument('-l', '--length', metavar='', action='store', type=int, default=conv_len, dest='conv_len',
		help='length of the conversion in seconds (default is until EOF)')

	# mjd start time
	parser.add_argument('-mjd', '--mjd-start', metavar='', action='store', type=float, default=mjd_time, dest='mjd_time',
		help='MJD start time (default is taken from creation time of the input file)')

	# highpass filter frequency
	parser.add_argument('-fh', '--filt-high', metavar='', action='store', type=float, dest='filt_h', default=f_h,
		help='frequency of the highpass filter in Hz')

	# lowpass filter frequency
	parser.add_argument('-fl', '--filt-low', metavar='', action='store', type=float, dest='filt_l', default=f_l,
		help='frequency of the lowpass filter in Hz')

	# ampli
	parser.add_argument('-a', '--ampli', metavar='', action='store', type=float, dest='ampli', default=ampli,
		help='amplification during conversion')

	# number of channels
	parser.add_argument('-c', '--channels', metavar='', action='store', type=float, dest='channels', default=channels,
		help='number of channels')

	# sample rate
	parser.add_argument('-sr', '--sample-rate', metavar='', action='store', type=int, dest='sample_rate', default=sample_rate,
		help='sample rate')

	# center frequency
	parser.add_argument('-f', '--center-freq', metavar='', action='store', type=float, dest='center_freq', default=center_freq,
		help='center frequency')

	# source name
	parser.add_argument('-n', '--source-name', metavar='', action='store', type=str, dest='source_name', default=source_name,
		help='source name')

	# source ra
	parser.add_argument('-ra', '--source-ra', metavar='', action='store', type=float, dest='source_ra', default=source_ra,
		help='source radial ascension')

	# source de
	parser.add_argument('-de', '--source-de', metavar='', action='store', type=float, dest='source_de', default=source_de,
		help='source declination')

	# sr error
	parser.add_argument('-e', '--sr-error', metavar='', action='store', type=float, dest='sr_error', default=sr_error,
		help='error in ppm of the dongle')

	# output file
	parser.add_argument('-o', '--output-file', metavar='', action='store', type=str, dest='out_file', default='',
		help='output file name (default is replaced only the extension)')

	# chunk size
	parser.add_argument('-cs', '--chunk-size', metavar='', action='store', type=int, dest='chunk_size', default=chunk_size,
		help='size of the chunks that are elaborated')

	# show saturation percentage
	parser.add_argument('-sp', '--saturation-percentage', action='store_true', dest='show_saturation',
		help='show the percentage of samples that are saturated during the clipping and a graph of the channels')

	# input file
	parser.add_argument('input_file', action='store', help='input bin file')

	# parse the arguments
	res = parser.parse_args()

	# update the global variables
	conv_start =		res.conv_start
	conv_len =			res.conv_len
	ampli =				res.ampli
	sample_rate =		res.sample_rate
	f_h =				res.filt_h
	f_l =				res.filt_l
	channels =			res.channels
	center_freq =		res.center_freq
	sr_error =			res.sr_error
	mjd_time =			res.mjd_time
	source_name =		res.source_name
	source_ra =			res.source_ra
	source_de =			res.source_de
	infile =			res.input_file
	chunk_size =		res.chunk_size
	show_saturation =	res.show_saturation

	# calulate the output filename with the correct extension
	if res.out_file == '':
		tmp = res.input_file.replace('.bin', '.fil')
		if tmp == res.input_file:
			outfile = ''.join((tmp, '.fil'))
		else:
			outfile = tmp
	else:
		outfile = res.out_file

	# print the conversion parameters
	print('conversion parameters:')
	for arg in vars(res):
		if arg == 'out_file':
			print('{:15s}{}'.format(arg, outfile))
			continue
		print('{:15s}{}'.format(arg, getattr(res, arg)))


def setup_filters(debug=False):
	# There are [channels] double filters
	# len(filters) is [channels]
	# each filter is characterized by two coefficients b and a
	# there are 2 filters, one is highpass and the other is lowpass
	# so fo example filters[0] = {'ah': , 'bh': , 'al': , 'bl':}

	global filters

	order_h = 1
	order_l = 3

	Wn_h = (f_h / sample_rate) * 2
	Wn_l = (f_l / sample_rate) * 2

	filters = []
	for i in range(channels):
		# highpass
		bh, ah = sig.butter(order_h, Wn_h, btype='high', output='ba')

		# lowpass
		bl, al = sig.butter(order_l, Wn_l, btype='low', output='ba')

		filters.append({
			'bh': bh,
			'ah': ah,
			'bl': bl,
			'al': al
			})

	if debug:
		print(filters)

def write_header(debug=False):
	# write the header of the .fil file
	global mjd_time

	file_length = os.path.getsize(infile)
	recording_length_secs = file_length / (2 * channels * sample_rate)
	mjd_time = (os.path.getmtime(infile) / 86400) + 40587 - recording_length_secs/86400 + conv_start/86400

	channel_width = -(sample_rate / channels) * 1e-6
	fch_1 = center_freq + (sample_rate * 2e-6) + (0.5 * channel_width)
	tsamp = (1 / sample_rate) * (1 / (sr_error * 1e-6) + 1)

	with open(outfile, 'wb') as f:

		f.write(struct.pack('<I', 12))
		f.write(bytearray('HEADER_START', 'ascii'))

		f.write(struct.pack('<I', 9))
		f.write(bytearray('data_type', 'ascii'))
		f.write(struct.pack('<I', 1))

		f.write(struct.pack('<I', 4))
		f.write(bytearray('nifs', 'ascii'))
		f.write(struct.pack('<I', 1))

		f.write(struct.pack('<I', 12))
		f.write(bytearray('telescope_id', 'ascii'))
		f.write(struct.pack('<I', 0))

		f.write(struct.pack('<I', 5))
		f.write(bytearray('nbits', 'ascii'))
		f.write(struct.pack('<I', 8))

		f.write(struct.pack('<I', 4))
		f.write(bytearray('foff', 'ascii'))
		f.write(struct.pack('<d', channel_width))

		f.write(struct.pack('<I', 4))
		f.write(bytearray('fch1', 'ascii'))
		f.write(struct.pack('<d', fch_1))

		f.write(struct.pack('<I', 6))
		f.write(bytearray('nchans', 'ascii'))
		f.write(struct.pack('<I', channels))

		f.write(struct.pack('<I', 5))
		f.write(bytearray('tsamp', 'ascii'))
		f.write(struct.pack('<d', tsamp))

		f.write(struct.pack('<I', 6))
		f.write(bytearray('tstart', 'ascii'))
		f.write(struct.pack('<d', mjd_time))

		f.write(struct.pack('<I', 11))
		f.write(bytearray('source_name', 'ascii'))
		f.write(struct.pack('<I', len(source_name)))
		f.write(bytearray(source_name, 'ascii'))

		f.write(struct.pack('<I', 7))
		f.write(bytearray('src_raj', 'ascii'))
		f.write(struct.pack('<d', source_ra))

		f.write(struct.pack('<I', 7))
		f.write(bytearray('src_dej', 'ascii'))
		f.write(struct.pack('<d', source_de))

		f.write(struct.pack('<I', 10))
		f.write(bytearray('HEADER_END', 'ascii'))

# function to handle the matplotlib close_event
def handle_close(evt):
	global running
	running = False

def elaborate(chunksize, debug=True):
	# elaborate and create the final .fil file
	# if chunksize == -1 then the entire file is loaded in memory

	if show_saturation:
		# initialize matplotlib figure, axis and close event
		fig, ax = plt.subplots()
		fig.canvas.mpl_connect('close_event', handle_close)

	# define the int16LE dtype
	dt = np.dtype('<i2')

	# create a file handler for the input file in binary mode
	ifh = open(infile, 'rb')

	# calculate some basic information for dividing the file in chunks
	file_length = os.path.getsize(infile)
	total_samples = file_length / (2 * channels)
	total_chunks = math.ceil(total_samples / chunksize) # this includes the last chunk

	for chunk in range(total_chunks):
		if (not running):
			return

		# read the binary file
		d = np.fromfile(ifh, dtype=dt, count=chunksize*channels)

		# reshape the ndarray to have each channel's data on each row
		d = d.reshape(-1, channels).transpose()

		# convert the array type to float32 for the calculation
		d = d.astype(np.float32, casting='safe')

		# multiply all the values by the ampli value
		d = d * ampli

		print('\rprocessing chunk\t{}/{}\t({:3}%)'
				.format(chunk, total_chunks, int(((chunk+1) / (total_chunks+1) )*100)), end='')

		# filter and clip data on each channel
		for c in range(channels):
			# actual filtering
			d[c] = sig.lfilter(filters[c]['bl'], filters[c]['al'],
				sig.lfilter(filters[c]['bh'], filters[c]['ah'], d[c])
				)

		# first part of clipping
		d = np.add(d, 128)

		if show_saturation:
			# calclate the saturated percentage and print it out
			saturated_samples = ((d<0) | (d>255)).sum()
			saturation_rateo = saturated_samples / d.size
			print("\tSAT {:.2f}%".format(saturation_rateo*100), end='')

		# second part of clipping
		d = np.clip(d, 0, 255)

		if show_saturation:
			# draw the channel graph
			ax.cla()
			plt.ylim(0, 255)
			plt.xlim(0, channels-1)
			ax.plot(d[:,0])
			plt.pause(1e-7)

		# convert the matrix into uint8
		d = d.astype(np.uint8)
		# flip the channels and reshape to have only one row
		d = np.flip(d, 0).transpose().reshape(1, -1)

		# append to the output file
		of = open(outfile, 'ab')
		d.tofile(of)

def main():
	start_time = time.time()

	parse_arguments(debug=False)
	setup_filters(debug=False)
	write_header(debug=False)
	elaborate(chunk_size, debug=False)

	print('\nprocessing finished in {:.1f} seconds'.format(time.time()-start_time))

if __name__ == '__main__':
	main()
