# Program: wvtlister.py
# Purpose: read a .hex file containing a wang 2200 tape dump
#          and report what is in it
# Version: 1.1, 2016/01/14
# Author:  Jim Battle

########################################################################
# This is a very crude utility to accept the decoded tape images produced
# by the WangWav2Tape program and turn them into program listings.
# It does very little in the way of error checking.  If the input is badly
# formed, the program will likely throw an exception somewhere.
#
# The output of this program would need to be split up manually to break
# out each program.
#
# A typical use would be:
#    ./WangWav2Tape.exe -v 1 -2200 -o output.hex edge.wav > log.txt
# Make sure it decodes cleanly, then
#    python wvlister.py output.hex > listing.txt

# note: the 2200 format on tape appears to be like this:
#    (1) burst of 66 '1' bits
#    (2) ~30 bit gap
#    (3) 2048 edges of the first copy of the data
#    (4) ~20 bit gap
#    (5) 2048 edges of the second copy of the data
#    (6) ~400 bit gap

# TODO:
#   *) more error checking
#   *) option to spit out each program as a separate file
#   *) deal with "SAVE P" type programs
#   *) deal with data files
#   *) perhaps make it a utility along the lines of wvdutil.py

import sys
import os.path
import string
import re

from wvfilelist import *

verbose = False

#############################################################################
# the wang 2200 stores each block twice, with a marker before each pair.
# blocks contains all the decoded data blocks since the last marker.
# pick the one most likely to be correct.
#############################################################################

def list_winner(blocks):
    listing = []
    if len(blocks) == 0:
        return listing
    blk_num, block = pick_winner(blocks)
    if verbose:
        listing.append("\n# absolute block %d" % blk_num)
    listing.extend(listOneBlock(block))
    return listing

def pick_winner(blocks):
    # check for the best, and hopefully most frequent, case:
    #     two blocks, and they match
    if len(blocks) == 2:
        idx1, blk1 = blocks[0]
        idx2, blk2 = blocks[1]
        if blk1 == blk2:
            return idx1, blk1

    # if one is 256 bytes, use that
    for (blk_num, block) in blocks:
        if len(block) == 256:
            return blk_num, block

    # just pick the longest one, I guess
    long_len, long_idx = 0,0
    idx = 0
    for (abs_idx, block) in blocks:
        if len(block) > long_len:
            long_len, long_idx = len(block), idx
        idx += 1
    return blocks[long_idx]

#############################################################################
# main program
#############################################################################

basename = os.path.basename(sys.argv[0])

# if no arguments, print help message and quit
if len(sys.argv) != 2:
    print 'Usage: ' + basename + ' <tape.hex>'
    sys.exit()

filename = sys.argv[1]
try:
    fh = open(filename, 'r')
except:
    print "Couldn't open '", filename, "'"
    sys.exit()

#############################################################################
# scan for blocks of data and for block markers
# from these, try to collect the best logical sequence of blocks
#############################################################################

blocklist = []
block = bytearray()

for line in fh:
    line = line.rstrip()
    if (len(line) == 0) or (line[0] == '#'):
        if len(block) > 0:
            # terminated the block
            blocklist.append(block)
            block = bytearray()
    else:
        # append this line of hex to current block
        block.extend(bytearray.fromhex(line))

# flush any residue
if len(block) > 0:
    # terminated the block
    blocklist.append(block)

#############################################################################
# dump all the collected blocks
#############################################################################

print "# The tape has", len(blocklist), "blocks"

listing = []
blocks = []

for i in range(len(blocklist)):
    block = blocklist[i]

    # a burst of 66 zero bits or 66 one bits comes between block pairs
    blen = len(block)
    start_marker    = (7 <= blen <= 8) and (block == (chr(0x00) * blen))
    continue_marker = (7 <= blen <= 8) and (block == (chr(0xFF) * blen))
    if start_marker:
        listing.extend(list_winner(blocks))
        blocks = []
        listing.append("\n# ===================================================")
        continue
    if continue_marker:
        listing.extend(list_winner(blocks))
        blocks = []
        continue

    blocks.append([ i, block ])

# there is no marker after the lsat
listing.extend(list_winner(blocks))

for line in listing:
    print line

sys.exit()

# vim: et:ts=4:sw=4
