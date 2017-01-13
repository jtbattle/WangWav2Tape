Wang 2200/3300 .wav to tape decoder
===================================

WangWav2Tape.cpp is a program I wrote in order to recover the data off of
some Wang 2200 and Wang 3300 cassette tapes, as I no longer have a native
machine that can read them.

The main site for things related to the Wang 2200 family of computers is
[wang2200.org](http://www.wang2200.org)

The main site for things related to the Wang 2200 family of computers is
[wang3300.org](http://www.wang3300.org)


Advice on Capturing Tapes
-------------------------

Normal audio cassettes have two pairs of tracks, set out like this:

```
                       --------------------------
  r/w head 0 --->      >>>>>>>>>> left >>>>>>>>>>
  r/w head 1 --->      >>>>>>>>>> right >>>>>>>>>
                       <<<<<<<<<< right <<<<<<<<<
                       <<<<<<<<<< left <<<<<<<<<<
                       --------------------------
```

The read/write heads are positioned offset on one side of the tape such
that a matching pair of left and right tracks are picked up.  When the
tape is flipped over, that offset is right to pick up the tracks going
in the opposite direction.

Unfortunately, the Wang 3300 and Wang 2200 wrote a single pair of tracks,
somewhat between the others.

```
                       --------------------------
                       >>>>>>>>>> left >>>>>>>>>>
  r/w head 0 --->      >>>>>>>>>> left >>>>>>>>>>

  r/w head 1 --->      >>>>>>>>>> right >>>>>>>>>
                       <<<<<<<<<< right <<<<<<<<<
                       --------------------------
```

The solution (worked out by
[Steve Witham](http://www.tiac.net/~sw/)
) was to get a cheap four-track multitrack cassette deck, the kind starving
musicians use to record demo tracks.  These decks have four r/w heads, and the
tape is only used in a single direction:

```
                       --------------------------
  r/w head 1 --->      >>>>>>>>> track 1 >>>>>>>>
  r/w head 2 --->      >>>>>>>>> track 2 >>>>>>>>
  r/w head 3 --->      >>>>>>>>> track 3 >>>>>>>>
  r/w head 4 --->      >>>>>>>>> track 4 >>>>>>>>
                       --------------------------
```

By mixing tracks 1 and 2 together to the "left" output of the unit, and
mixing tracks 3 and 4 together to the "right" output, we can read the
original Wang track that was between those pairs.  Sometimes, based on
head/track alignment, it is better to just use tracks 1 and 3, or 2 and 4,
or 1 and 4, etc.  Always review the captured waveforms to make sure they
look clean and strong, and if they aren't, this may be one way to improve it.

If the deck has a doing Dolby noise reduction setting, disable it, as
precompensation wasn't used when Wang recorded the data to the tape.

A very important step is to clean the tape head frequently.  I do it before
every single tape I capture.  Invariably, these tapes are quite old, and
they can suffer from
[sticky-shed syndrome](https://en.wikipedia.org/wiki/Sticky-shed_syndrome).
What that means is that bits of the oxide and binder can come off the tape
surface and accumulate on the read head.  One sure sign you have this problem
is if the tape starts making a loud squealing sound on playback -- not the
audio, but the tape mechanically is squealing as it slides past the read
head.  Use a Q-Tip dipped in isopropyl alcohol, rub the tape head to scrub
off/dissolve any gunk, then use the dry end of the Q-Tip to wipe off any
excess alcohol.  Blow on it a bit to ensure it is all evaporated, rewind
the tape a bit, and capture some more.

Before capturing, play the tape a bit and set the recording levels so you
have a strong a strong signal with no clipping. Then rewind, and capture
the audio.

The output of the multi-track recorder is sent to a PC audio application
to capture the resulting stereo waves.  In my case, I used the excellent open
source application, Audacity.  I've had success capturing at 44.1 KHz and at
22050 Hz.  44.1 KHz is overkill as the signal is quite band-limited.
Capture the audio at 16b resolution.

When saving the final .wav to a disk, save it at whatever rate it was
captured at, 16-bit, uncompressed PCM format.

Preprocessing the Captured Waveform
-----------------------------------

First, take the captured wave and save it to "raw.wav", just in case you
mess up the subsequent filtering steps, or want to experiment with different
preprocessing steps.

Next, apply a high pass filter, with these parameters:

* 6 db/octave roll-off
* cutoff: 1500 Hz for 2200 tapes, 1000 Hz for 3300 tapes

These settings aren't very critical.

Next, use the Normalize effect
* set limit of 3dB or so
* check the box to remove D.C. bias
* check the box to normalize the left and right tracks independently

You are then left with a waveform that has a sharp peak corresponding to
each flux transition on the original tape.  Spend a few moments eyeballing
it to make sure there are no gross problems.

This is saved to a file; I use "edge.wav", but anything will do.

In this directory, raw.wav and edge.wav samples have been extracted from
an actual Wang 2200 tape and processed as described above.  Take a look
at them under Audacity or another waveform viewer to get an idea of what
things should look like.

Running the Program
-------------------

This program was originally developed in 2008 on a Windows XP machine,
using the then-current g++ compiler in a cygwin environment, but recently
I recompiled it using MS Visual Studio 2015 C++ compiler without any
changes.

To see the command line options, type

    ./WangWav2Tape.exe

Typical usage would be:

    ./WangWav2Tape.exe -v 1 -2200 -o output.hex edge.wav > log.txt

This says:

* produce log messages at verbosity level 1
* the tape is Wang 2200 format
* output the final data image to a file named output.hex
* process the file "edge.wav"
* save the log messages to "log.txt"

Certain important errors are logged to stderr.

After it finishes (about 10 seconds for a 140 MB input file), open the
output.hex and log.txt files in an editor.  Before each block of the
output.hex file, there is a small comment header indicating if that block
decoded correctly or not.  In "vi", one can use this command to print
just the header blocks to review it more easily:

```
    :g/^#
```

What if a Tape Doesn't Decode Cleanly?
--------------------------------------

Note that the comment will say `Raw tape block #<some number>`.
This number won't just increment from 1 to N.  That is, it isn't the
logical block number on the tape, but the detected block number on
the actual recording.

Both the 2200 and the 3300 record each data block twice in a row
for redundancy, to guard against tape errors.  If a problem
occurs in decoding, take that raw block number and look for that
same block number in "log.txt".  There you can find some more detailed
information about the decoding process.

You can also try running with the verbosity level higher -- up to
"-v 4".  Then look at the log file again.  Often one can find the
actual transition where things go wrong.

Open the raw.wav and edge.wav files in Audacity.  Say you identify
that the program found a peak at sample "@404138" in channel 0
(the left channel), but it got rejected for some reason.  In Audacity,
move the time cursor to sample 404138 and have a look at the edge
and raw representations at that time.

Perhaps there was a "pop" in the tape.  Use Audacity to fix it by hand.
Or maybe that block gets hit in the middle with "fading" -- lower
amplitude, perhaps due to dust or grime on that spot on the tape.
Using Audacity, try boosting the waves in just that section of the
waveform.

Or if you are really looking for a challenge, try and read my code,
play with some of the hard-coded parameters, or maybe even add some
heuristics to better deal with troublesome situations.

Future Changes
--------------

There are quite a few ways this code could be enhanced.  However, at the
moment it is working pretty well for the tapes that I have, so I don't
plan on adding those improvements until they are needed.

Some ideas:

    *) The current code looks for the gaps between blocks after
       the raw samples have already been converted into a list of
       peaks.  Statistics about the peak distribution is taken next.
       However, in the original pass of converting the wav into the
       peak list, some of that information would have been useful
       in determining what is a peak or not.

       A better approach would be to identify peaks, take stats,
       then reprocess the waves again using this information to help
       refine what are plausible peaks or not.

       Specifically, we know apriori approximately what the bit timing
       should look like.  The timing stats should ignore any gaps where
       it appears we've dropped a bit.  This should give a better estimate
       of the mean and variance for edge timing.  When reprocessing, we
       should sync on known peaks, and then use our knowledge of timing
       to know where to look for other peaks.  Although, truthfully, this
       might be overkill as the real hardware does nothing like this and
       simply detects each peak in isolation.

    *) Different blocks of tape could have been written at different
       times, even by different machines.  Currently, my procedure
       does dc and amplitude normalization across the entire tape.
       It would be pretty simple to renormalize the left and right
       peak amplitudes while doing block level processing.  It would
       be even better to combine it with the above idea and do that
       renormalization to the raw waves a block at a time.

    *) The list of peaks is generated for the left and right tracks
       independently.  The way the Wang encodes bits, there can be
       a transition on one or the other tracks, not both at the same
       time.  This information could be useful in discriminating
       peaks that are ambiguous.

    *) I could have spent a lot more time tweaking the hard-coded
       parameters that do appear in the program.  I could have exposed
       them and allowed them to be set on the command line.

    *) The more verbose settings spew a lot of information.  I should
       really hang all of that information into a list associated with
       each data block.  At the very end, if the block decoded cleanly,
       all of that diagnostic information could be suppressed, and
       printed only for those few blocks that have decoding problems.

    *) put the normalization, eq, and high pass filters into the code
       and optionally enable them from the command line such that
       in many cases, Audacity isn't even needed.

Credits
----------

This code was written by Jim Battle,
[wang2200.org](http://www.wang2200.org)
web master.

However, I owe
[Steve Witham](http://www.tiac.net/~sw/)
a lot of credit for leading the way with his own approach to decoding, and for
coming up with the idea of using a multitrack tape deck to get to the data
tracks.  He was also a great sounding board while I developed this code.
Thanks, Steve.

License
----------
This code is released under the MIT License.
