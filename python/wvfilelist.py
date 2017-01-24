# Program: wvfilelist.py, used by wvdutil.py
# Purpose: encapsulates some useful routines for manipulating .wvd files
# Version: 1.0, 2008/01/29
# Version: 1.1, 2016/01/14
# Author:  Jim Battle

import string

###############################################################################
# table of BASIC atoms
###############################################################################


def initBasicTokens():
    global token
    token = [None] * 256
    token[0x80] = "LIST "
    token[0x81] = "CLEAR "
    token[0x82] = "RUN "
    token[0x83] = "RENUMBER "
    token[0x84] = "CONTINUE "
    token[0x85] = "SAVE "
    token[0x86] = "LIMITS "
    token[0x87] = "COPY "
    token[0x88] = "KEYIN "
    token[0x89] = "DSKIP "
    token[0x8A] = "AND "
    token[0x8B] = "OR "
    token[0x8C] = "XOR "
    token[0x8D] = "TEMP"
    token[0x8E] = "DISK "
    token[0x8F] = "TAPE "

    token[0x90] = "TRACE "
    token[0x91] = "LET "
#   token[0x92] = "FIX("        # BASIC-2
    token[0x93] = "DIM "
    token[0x94] = "ON "
    token[0x95] = "STOP "
    token[0x96] = "END "
    token[0x97] = "DATA "
    token[0x98] = "READ "
    token[0x99] = "INPUT "
    token[0x9A] = "GOSUB "
    token[0x9B] = "RETURN "
    token[0x9C] = "GOTO "
    token[0x9D] = "NEXT "
    token[0x9E] = "FOR "
    token[0x9F] = "IF "

    token[0xA0] = "PRINT "
    token[0xA1] = "LOAD "
    token[0xA2] = "REM "
    token[0xA3] = "RESTORE "
    token[0xA4] = "PLOT "
    token[0xA5] = "SELECT "
    token[0xA6] = "COM "
    token[0xA7] = "PRINTUSING "
    token[0xA8] = "MAT "
    token[0xA9] = "REWIND "
    token[0xAA] = "SKIP "
    token[0xAB] = "BACKSPACE "
    token[0xAC] = "SCRATCH "
    token[0xAD] = "MOVE "
    token[0xAE] = "CONVERT "
    token[0xAF] = "PLOT "        # "PLOT" in a SELECT PLOT context

    token[0xB0] = "STEP "
    token[0xB1] = "THEN "
    token[0xB2] = "TO "
    token[0xB3] = "BEG "
    token[0xB4] = "OPEN "
    token[0xB5] = "CI "
    token[0xB6] = "R "
    token[0xB7] = "D "
    token[0xB8] = "CO "
#   token[0xB9] = "LGT("        # BASIC-2 only
    token[0xBA] = "OFF "
    token[0xBB] = "DBACKSPACE "
    token[0xBC] = "VERIFY "
    token[0xBD] = "DA "
    token[0xBE] = "BA "
    token[0xBF] = "DC "

    token[0xC0] = "FN"        # no trailing space
    token[0xC1] = "ABS("
    token[0xC2] = "SQR("
    token[0xC3] = "COS("
    token[0xC4] = "EXP("
    token[0xC5] = "INT("
    token[0xC6] = "LOG("
    token[0xC7] = "SIN("
    token[0xC8] = "SGN("
    token[0xC9] = "RND("
    token[0xCA] = "TAN("
    token[0xCB] = "ARC"
    token[0xCC] = "#PI"       # no trailing space
    token[0xCD] = "TAB("
    token[0xCE] = "DEFFN"     # no trailing space
    token[0xCF] = "ARCTAN("

    token[0xD0] = "ARCSIN("
    token[0xD1] = "ARCCOS("
    token[0xD2] = "HEX("
    token[0xD3] = "STR("
    token[0xD4] = "ATN("
    token[0xD5] = "LEN("
    token[0xD6] = "RE"         # no trailing space (used by REDIM?)
    token[0xD7] = "#"
    token[0xD8] = "%"          # no trailing space
    token[0xD9] = "P"          # no trailing space
    token[0xDA] = "BT"         # no trailing space
    token[0xDB] = "G"          # no trailing space
    token[0xDC] = "VAL("
    token[0xDD] = "NUM("
    token[0xDE] = "BIN("
    token[0xDF] = "POS("

    token[0xE0] = "LS="        # no trailing space
    token[0xE1] = "ALL"        # no trailing space
    token[0xE2] = "PACK"       # no trailing space
    token[0xE3] = "CLOSE"      # no trailing space
    token[0xE4] = "INIT"       # no trailing space
    token[0xE5] = "HEX"        # eg HEXPRINT, not HEX(, no trailing space
    token[0xE6] = "UNPACK"     # no trailing space
    token[0xE7] = "BOOL"       # no trailing space
    token[0xE8] = "ADD"        # no trailing space
    token[0xE9] = "ROTATE"     # no trailing space
    token[0xEA] = "$"          # no trailing space
    token[0xEB] = "ERROR"      # no trailing space
#   token[0xEC] = "ERR"        # BASIC-2 only, no trailing space
#   token[0xED] = "DAC "       # BASIC-2 only
#   token[0xEE] = "DSC "       # BASIC-2 only
#   token[0xEF] = "SUB"        # BASIC-2 only, no trailing space

#   token[0xF0] = "LINPUT "    # BASIC-2 only
#   token[0xF1] = "VER"        # BASIC-2 only
#   token[0xF2] = " ELSE "     # BASIC-2 only
#   token[0xF3] = "SPACE"      # BASIC-2 only, no trailing space
#   token[0xF4] = "ROUND("     # BASIC-2 only
#   token[0xF5] = "AT("        # BASIC-2 only
#   token[0xF6] = "HEXOF("     # BASIC-2 only
#   token[0xF7] = "MAX("       # BASIC-2 only
#   token[0xF8] = "MIN("       # BASIC-2 only
#   token[0xF9] = "MOD("       # BASIC-2 only

    # all others are undefined
    return

###############################################################################
# given a list of file blocks, return a program file listing
###############################################################################


def listOneBlock(block):

    blen = len(block)
    listing = []

    # it is just too short to be useful
    if blen <= 10:
        listing.append("# short block, len=%d" % blen)
        listing.extend(listHexBlock(block))
        return listing

    data_flag    = (block[0] & 0x80) == 0x80
    pgm_flag     = (block[0] & 0x80) == 0x00
    header_flag  = (block[0] & 0x40) == 0x40
    trailer_flag = (block[0] & 0x20) == 0x20
    # listing.append("# length=%d, block[0]=0x%02x --> data=%d, pgm=%d, header=%d, trailer=%d" %
    #                 (blen, block[0], data_flag, pgm_flag, header_flag, trailer_flag))

    if blen != 256:
        # listing.append("# wrong block size; just dumping it")
        listing.extend(listHexBlock(block))
        return listing

    if header_flag:
        listing.extend(listHeaderBlock(block))
        return listing

    if pgm_flag:
        listing.extend(listProgramBlock(block))
        return listing

    if data_flag:
        listing.extend(listDataBlock(block))
        return listing

###############################################################################
# given a header record, return the filename it contains
###############################################################################
# NB: disk headers pad names shorter than 8 characters with 0x20
#     tape headers pad names shorter than 8 characters with 0xFF


def listHeaderBlock(block):
    listing = []

    fd_byte = block[9]
    if fd_byte != 0xFD:
        listing.append("# file header block doesn't have expected format")
        return listing

    name = "".join(map(chr, block[1:8]))
    while len(name) > 0:
        if (name[-1] != chr(0xff)) and (name[-1] != ' '):
            break
        name = name[:-1]

    listing.append("# filename = '%s'" % name)

    pgm_flag       = ((block[0] & 0x80) == 0x00)
    protected_flag = ((block[0] & 0x10) == 0x10)
    if pgm_flag and protected_flag:
        listing.append("# protected file flag set")

    return listing

###############################################################################
# dump a block as hex bytes
###############################################################################
# TODO: add ASCII map to the right


def listHexBlock(block):
    line = ""
    listing = []
    blen = len(block)
    for n in range(blen):
        line += "%02x" % block[n]
        if ((n % 16) == 15) or (n == blen-1):
            listing.append(line)
            line = ""
    listing.append("")
    return listing

###############################################################################
# dump a block as tokenized BASIC
###############################################################################


def listProgramBlock(block):
    global token
    try:
        token
    except:
        initBasicTokens()

    # the remaining blocks contain the program
    limit = len(block)
    off = 1                 # skip control byte
    listing = []
    line = ""

    while (off < 255) and (off < limit) and (block[off] != 0xFD) and (block[off] != 0xFE):
        c = block[off]
        if c == 0x0D:
            off += 2  # there are two 00 bytes after line end
            listing.append(line)
            line = ""
        elif (c < 128):
            line += chr(c)
        elif c == 0xFF and (off+2 < limit):
            # line number
            linenum = 1000*(block[off+1] / 16) + \
                       100*(block[off+1] % 16) + \
                        10*(block[off+2] / 16) + \
                         1*(block[off+2] % 16)
            off += 2
            line += "%d" % linenum
        elif token[c] is None:
            line += "\\x%02x" % c
        else:
            line += token[c]
        off += 1

    # if len(line) > 0:
    #    listing.append(line)

    return listing

###############################################################################
# given a list of file blocks for a structured data file (that is,
# produced by DATASAVE), produce a human readable listing
###############################################################################


# TODO: convert the listDataFromBlocks function when I have an example
def listDataBlock(block):
    return listHexBlock(block)


def listDataFromBlocks(blocks):

    secOffset = -1
    logicalRecNum = 0
    listing = []

    for sec in (blocks[0:-1]):  # skip the last sector
        secOffset += 1
        c0, c1 = sec[0], sec[1]
        # print ">>> sector control bytes 0x%02X 0x%02X" % (c0, c1)
        if (c0 & 0xA0) == 0xA0:
            listing.append('### trailer record; done')
            return listing
        if (c0 & 0x80) == 0:
            listing.append('Sector %d of data file has bad header byte 0x%02X' % (secOffset, c0))
            return listing
        if c1 == 0x01:
            listing.append('### Logical Record %d' % logicalRecNum)
            logicalRecNum += 1

        i = 2  # skip control bytes
        while (i < 255) and (ord(sec[i]) < 0xFC):
            c = sec[i]
            # print "i=%d, SOV=0x%02X" % (i, c)
            if c == 0x08:  # numeric
                fp = sec[i+1:i+9]
                expSign  = (ord(fp[0]) & 0x80) != 0
                mantSign = (ord(fp[0]) & 0x10) != 0
                exp =  1 * ((ord(fp[0]) >> 0) & 0xF) \
                    + 10 * ((ord(fp[1]) >> 4) & 0xF)
                mant = ("%d" % ((ord(fp[1]) >> 0) & 0xF)) + "."
                for d in fp[2:8]:
                    mant += ("%02X" % ord(d))
                fpascii = (" ", "-")[mantSign] + mant + \
                          ("+", "-")[ expSign] + ("%02X" % exp)
                listing.append(fpascii)
                i += 9
            elif (c < 128):
                listing.append("Unexpected SOV value of 0x%02X" % c)
                return listing
            else:
                strlen = c-128
                if i+1+strlen > 255:
                    print "SOV indicates string length that doesn't fit in sector"
                    return
                str = sec[i+1 : i+1+strlen]
                # expand non-ascii characters
                foo = [None] * len(str)
                for n in xrange(len(str)):
                    if (ord(str[n]) < 32) or (ord(str[n]) >= 128):
                        foo[n] = "\\x%02X" % ord(str[n])
                    else:
                        foo[n] = str[n]
                str = ''.join(foo)
                listing.append('"' + str + '"')
                i += 1 + strlen

    return listing

# vim: et:ts=4:sw=4
