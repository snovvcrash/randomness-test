#
# gen_bit_sequences.py
#
# Naive Bit Sequence Generator
# by snovvcrash
# 04.2017
#

import secrets # Python 3.6
import random
import os.path

subdir = "tests"
try:
	os.mkdir(subdir)
except Exception:
	pass

# /dev/urandom

rng = secrets.SystemRandom()

with open(os.path.join(subdir, "urnd.txt"), "w") as f:
	for i in range(32768):
		if (i % 64 == 0 and i):
			f.write('\n')
		f.write(str(rng.randint(0, 1)))

# Pseudo-random

random.seed()

with open(os.path.join(subdir, "pseudornd.txt"), "w") as f:
	for i in range(32768):
		if (i % 64 == 0 and i):
			f.write('\n')
		f.write(str(random.randint(0, 1)))

# "Bad" sequence

if (os.path.isfile(os.path.join(subdir, "bad.txt")) == 0):
	with open(os.path.join(subdir, "bad.txt"), "w") as f:
		for i in range(32768):
			if (i % 64 == 0 and i):
				f.write('\n')
			if (i % 2 == 0):
				f.write('0')
			else:
				f.write('1')
