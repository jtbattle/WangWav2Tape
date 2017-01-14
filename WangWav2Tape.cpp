// Prerelease version, 1/13/16
// Author: Jim Battle

// This program reads a .wav file containing sampled data from a Wang 3300 or
// Wang 2200 cassette tape.
//
// The Wang tape format saves data in two tracks.  However, the track placement
// isn't the same as a normal audio cassette.  To recover the data, one must
// use a four track music recording deck.  One data track can be recovered by
// mixing together the channel 0 & 1 tracks, and the other data track can be
// recovered by mixing the channel 2 & 3 tracks.
//
// One channel contains a transition for each 0 bit, and the other channel
// transitions on each 1 bit.  There is no preamble pattern; a gap is just
// the absense of any transitions for an extended duration.
//
// There isn't a whole lot of high-powered thinking going on here.
// It is just a lot of guess work on my part as to what seems to work
// and what doesn't.  Start simple and add complexity as required.
//
// The procedure is this:
//
//   (1) capture the tape to a .wav file
//
//   (2) manually review the waves and fix things (audacity is good):
//          remove dc bias
//          normalize the amplitude
//          remove pops or obvious noise
//          apply some kind equalization (this requires some trial and error)
//              eg,   0 db from  dc    to 200 Hz
//                  +10 db from 200 Hz to inf
//          apply a high-pass filter to reveal the edges
//              eg, highpass with 3db/octave rolloff, q=0.7, cutoff: 1001 hz
//          normalize the amplitude again
//          save the modified wave to another file
//       this should result in a clearly defined spike for each rising or
//       falling edge; the peak corresponds to the part of the edge slope
//       with the maximum rate of change (typically the middle)
//
//   (3) feed the processed waveform as the input to this program
//
// The program does this:
//
//   (1) It scans the two channels, detecting each peak and building a list
//       logging where each of these peaks occurs in the waveform.  BTW,
//       rather than reading the full file into memory, a sliding window
//       is used.  This permits the program to think it is randomly accessing
//       the sample buffer, without the memory requirements.
//
//       At this point, there is no more use for the waves, and the file
//       can be closed.
//
//   (2) The separate lists of left and right peaks are merged into a time-
//       ordered list.
//
//   (3) The list of peaks are scanned for gaps, which demark data blocks.
//       For each isolated block,
//
//       (3a) Transitions may occur in the left and right tracks in any
//            order, but not both at the same time, we may have L->L, L->R,
//            R->R, and R->L transition sequences.  Statistics about the
//            average and variance of these transitions are kept for each
//            of the four categories.
//
//       (3b) It may be, for various reasons, there is a time shift between
//            the left and right tracks.  Based on the transition stats,
//            shift one track or the other by some amount of time such that
//            the average L->R and R->L transitions are identical.  We must
//            be prepared, though, for a block of data which doesn't have
//            all types of transitions (eg, all 0s).
//
//       (3c) Sweep through the transition lists again like (3a), measuring
//            the transition statistics, but with the time shift taken into
//            account.
//
//       (3d) Sweep through the transition lists again, and discard unlikely
//            peaks.  A peak is unlikely if it doesn't occur within about
//            three std deviations of the mean, yet the following edge does.
//            That means that the discarded edge was likely spurious.
//
//       (3e) Do one final sweep, turning the peaks into ones and zeros,
//            assembling each group of 8 into a byte.  The 3300 tape format
//            includes checksum bytes, which aids in verifying a successful
//            decode.  Both the 3300 and 2200 tape formats record each block
//            of data twice, which also helps verify the result.  While the
//            2200 doesn't have a checksum (why, oh why, not?), we know that
//            there should be exactly 256 bytes.  If we end up with the wrong
//            byte count, or if we have extra transition bits at the end,
//            it is a strong clue that something is amiss.
//
// TODO and further ideas:
//
//   *) sometimes a section of the tape will show fading -- localized loss
//      of amplitude.  It might last only a few dozen bits (probably due to
//      poor tape oxide thickness, or something stuck on the tape to move
//      the r/w head away from the tape surface), but it can cause edges to
//      be dropped.  one idea is to have a very localized normalization.
//      eg, have a rolling window of, about 10 bits wide.  based on the
//      average amplitude, or perhaps average peak amplitude, normalize
//      the signal level.  this is most easily done at the edge stage,
//      since I don't want to write out or store a modified copy of the
//      wav data.
//
//   *) the 2200 tape format records each block twice for redundancy,
//      yet neither has a checksum.  thus, if a bit got flipped in one,
//      there is no simple way to know which is correct.
//
//      it would be thinkable to parse the block contents to make sure it
//      has valid structure.  while not perfect, it is better than nothing.
//      BASIC and structured data files have a known format.

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>        // sqrt() requires this
#include <assert.h>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::ios;
using std::endl;
using std::hex;
using std::dec;
using std::fixed;
using std::setfill;
using std::setw;

// ========================================================================
// type definitions

typedef unsigned long  uint32;
typedef unsigned short uint16;
typedef unsigned char  uint8;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(a)   (((a)<(0))?(-a):(a))

// pack two and four byte items in little endian fashion
#define PACK2(a, b) (((uint8)a <<  0) + ((uint8)b <<  8))
#define PACK4(a, b, c, d) \
    (((uint8)a <<  0) + ((uint8)b <<  8) + ((uint8)c << 16) + ((uint8)d << 24))


// ========================================================================
// global variables

// command line options
int    opt_v;           // verbosity level
string opt_ofn;         // output filename
string opt_ifn;         // input filename
enum { FMT_3300, FMT_2200 } opt_fmt;

// any peak with amplitude less this is simply ignored
// FIXME: it is currently predicated on opt_fmt, but it might be useful to
//        expose it as a command line option.
float opt_min_thold;

// ugly globals to avoid passing them around as parameters
int g_bit_per;          // expected bit period
int g_gap_thold;        // threshold for gap detection

float g_devs_filter = 2.5f;  // # std devs before an runt interval is filtered out
float g_devs_warn   = 3.5f;  // how many std devs before an interval is flagged

// ========================================================================
// this class opens a .wav file given a filename.
// there are methods to return various properties of the file, as well
// as the left and right sample at offset "n".  The sample is returned
// as a pair of single precision floats, normalized in the range of 0.0
// to 1.0.  If the file is mono, the left and right samples are the same;
// sample depths of 8b and 16b are both supported.
//
// The assumption is that samples will be processed approximately sequentially.
// It will return the correct value, no matter what order, but performance
// may be terrible if the access pattern is poor.
//
// This code isn't very efficient, just simple.

struct samplePair_t {
    float channel[2];
};

class ReadWav
{
public:
    ReadWav(string filename);
    ~ReadWav();

    // samples/sec
    int SampleFrequency() const;

    // number of samples in file, as declared by header
    // (which can be wrong -- perhaps the file is truncated)
    long ExpectedFrames() const;
    // calculated actual sample count
    long NumFrames() const;

    // random access interface:
    // return the sample pair at specified sample offset
    samplePair_t Sample(long offset);

    // sequential access interface:
    // return the sample pair at specified sample offset
    void SetOffset(long off);
    samplePair_t NextFrame();

private:
    // methods to read items directly from the file stream
    // values interpret the stream in little endian format.
    uint8  read_raw_8b();
    uint16 read_raw_16b();
    uint32 read_raw_32b();

    // save it for error messages later
    string m_filename;

    // the wav file stream object
    ifstream m_istream;

    // cache attributes read in from the file header
    int  m_num_channels;
    int  m_bytes_per_sample;
    int  m_bytes_per_frame;     // == #channels * (bytes/channel)
    int  m_sample_rate;
    long m_expected_frames;
    long m_actual_frames;
    long m_data_start;
    long m_data_end;

    // we keep a direct cache of N sample buffers, each holding K samples.
    const static int N_buffers_log2 = 4;
    const static int K_samples_log2 = 10;
    const static int N_buffers = (1 << N_buffers_log2);
    const static int K_samples = (1 << K_samples_log2);
    samplePair_t m_cache[N_buffers][K_samples];
    bool m_cache_valid[N_buffers];
    long m_cache_tag[N_buffers];

    // sequential sample interface
    long m_seq_offset;
};


uint8
ReadWav::read_raw_8b()
{
    char buff[1];
    m_istream.read(buff, 1);
    if (!m_istream.good()) {
        cerr << "Error: failure to read 8b quantity from file stream" << endl;
        exit(-1);
    }

    return (uint8)buff[0];
}


uint16
ReadWav::read_raw_16b()
{
    char buff[2];
    m_istream.read(buff, 2);
    if (!m_istream.good()) {
        cerr << "Error: failure to read 16b quantity from file stream" << endl;
        exit(-1);
    }

    return PACK2( buff[0], buff[1] );
}


uint32
ReadWav::read_raw_32b()
{
    char buff[4];
    m_istream.read(buff, 4);
    if (!m_istream.good()) {
        cerr << "Error: failure to read 32b quantity from file stream" << endl;
        exit(-1);
    }

    return PACK4( buff[0], buff[1], buff[2], buff[3] );
}


// RIFF WAV file format
//  __________________________
// | RIFF WAVE Chunk          |
// |   groupID  = 'RIFF'      |
// |   riffType = 'WAVE'      |
// |    __________________    |
// |   | Format Chunk     |   |
// |   |   ckID = 'fmt '  |   |
// |   |__________________|   |
// |    __________________    |
// |   | Sound Data Chunk |   |
// |   |   ckID = 'data'  |   |
// |   |__________________|   |
// |__________________________|
//
// although it is legal to have more than one data chunk,
// this program assumes there is only one.

ReadWav::ReadWav(string filename)
{
    // save filename for error messages
    m_filename = filename;

    // invalidate the cache
    for (int i=0; i<N_buffers; i++)
        m_cache_valid[i] = false;

    if (opt_v >= 1) {
        cout << "File: '" << m_filename << "'\n";
    }

    m_istream.open(m_filename.c_str(), ios::in | ios::binary);
    if (!m_istream.is_open()) {
        cerr << "Error: couldn't open file '" << m_filename << "'" << endl;
        exit(-1);
    }

    int chunkNum = 0;
    bool done = false;

    while (!done) {
        // now read the next chunk ID to see what is up next
        chunkNum++;
        uint32 chunkId   = read_raw_32b();
        uint32 chunkSize = read_raw_32b();

        if (chunkId == PACK4('R','I','F','F')) {

            uint32 riffType = read_raw_32b();
            if (riffType != PACK4('W','A','V','E')) {
                cerr << "Error: input file not a WAV file" << endl;
                exit(-1);
            }

        } else if (chunkId == PACK4('f','m','t',' ')) {

            uint16 fmtTag = read_raw_16b();
            if (fmtTag != 1) {
                cerr << "Error: can't deal with compressed WAV files" << endl;
                exit(-1);
            }

            uint16 chan = read_raw_16b();
            if (chan < 1 || chan > 2) {
                cerr << "Error: can't handle more than two channels" << endl;
                exit(-1);
            }
            m_num_channels = chan;

            m_sample_rate = read_raw_32b();
            if (m_sample_rate < 11000) {
                cerr << "Warning: the sample rate is low -- it might hurt conversion" << endl;
            }

            uint32 avgpbs = read_raw_32b();
            uint16 blockAlign = read_raw_16b();

            uint16 sampbits = read_raw_16b();
            if (sampbits == 8) {
                m_bytes_per_sample = 1;
            } else if (sampbits == 16) {
                m_bytes_per_sample = 2;
            } else {
                cerr << "Error: samples must be either 8b or 16b" << endl;
                exit(-1);
            }
            m_bytes_per_frame = m_bytes_per_sample * m_num_channels;

        } else if (chunkId != PACK4('d','a','t','a')) {

            // the first chunk MUST be a RIFF type
            if (chunkNum == 1) {
                cerr << "Error: input file not a WAV file" << endl;
                exit(-1);
            }

            cerr << "Unknown chunk type: '" <<
                    (char)((chunkId >>  0) & 0xFF) <<
                    (char)((chunkId >>  8) & 0xFF) <<
                    (char)((chunkId >> 16) & 0xFF) <<
                    (char)((chunkId >> 24) & 0xFF) << "'" << endl;
            cerr << "    " << chunkSize << " bytes skipped\n";

            // skip rest of header
            m_istream.seekg(chunkSize, ios::cur);
            if (!m_istream.good()) {
                cerr << "Error: failure to seek to end of chunk" << endl;
                exit(-1);
            }

        } else {
            // got a data chunk

            // remember where data begins
            m_data_start = m_istream.tellg();

            // compute dependent parameters
            m_expected_frames = chunkSize / m_bytes_per_frame;

            m_istream.seekg(0, ios::end);
            m_data_end = m_istream.tellg();
            m_actual_frames = (m_data_end - m_data_start + 1) / m_bytes_per_frame;

            SetOffset(0);

            done = true;
        }

    } // while (!done)


    if (opt_v >= 1) {
        cout << "WAV format:\n";
        cout << "    uncompressed\n";
        cout << "    " << m_num_channels << " channels\n";
        cout << "    " << (8*m_bytes_per_sample) << " bits/sample\n";
        cout << "    " << m_sample_rate << " samples/sec\n";
        cout << "    " << m_expected_frames << " samples expected\n";
        cout << "    " << m_actual_frames << " samples in file\n";
        cout << "    " << (float)m_expected_frames/m_sample_rate << " seconds\n";
    }
}


ReadWav::~ReadWav()
{
    m_istream.close();
}


int
ReadWav::SampleFrequency() const
{
    return m_sample_rate;
}


long
ReadWav::ExpectedFrames() const
{
    return m_expected_frames;
}


long
ReadWav::NumFrames() const
{
    return m_actual_frames;
}


samplePair_t
ReadWav::Sample(long offset)
{
    // |<--- buf_tag --->|<--- buf --->|<--- buf_offset --->|
    const long clamped_offset = MAX(0, offset);
    const int buf_offset = (clamped_offset & (K_samples - 1));
    const int blk_offset = (clamped_offset >> K_samples_log2);
    const int buf        = blk_offset & (N_buffers - 1);
    const long buf_tag = (clamped_offset >> (N_buffers_log2 + K_samples_log2));

    if (!m_cache_valid[buf] || (m_cache_tag[buf] != buf_tag)) {

        // cache miss
        m_cache_valid[buf] = true;
        m_cache_tag[buf]   = buf_tag;
        long block_first   = m_data_start + (K_samples*m_bytes_per_frame)*blk_offset;
        long block_last    = block_first  + (K_samples*m_bytes_per_frame) - 1;

#if 0
        cout << "Cache miss: buf #" << hex << buf <<
                           ", tag=" << hex << buf_tag <<
                           ", off=" << hex << buf_offset <<
                           ", start=" << hex << block_first <<
                           ", end  =" << hex << block_last << endl;
#endif

        if (block_first > m_data_end) {

#if 0
            cout << "Filling end with 0.0\n";
#endif
            // none of the buffer exists
            for (int i=0; i<K_samples; i++) {
                m_cache[buf][i].channel[0] =
                m_cache[buf][i].channel[1] = 0.0f;
            }

        } else {

            // at least part of the buffer exists in the file
            m_istream.seekg(block_first);
            if (!m_istream.good()) {
                cerr << "Error: couldn't seek to valid file offset" << endl;
                exit(-1);
            }

            char raw[4*K_samples];
            m_istream.read( raw, K_samples*m_bytes_per_frame );
            m_istream.clear();

            if ((m_num_channels == 1) && (m_bytes_per_sample == 1)) {
                // mono 8b
                for (int n=0; n<K_samples; n++) {
                    int s0 = (uint8)raw[n] - 128;
                    m_cache[buf][n].channel[0] =
                    m_cache[buf][n].channel[1] = s0 / 128.0f;
                }
            } else if ((m_num_channels == 1) && (m_bytes_per_sample == 2)) {
                // mono 16b
                for (int n=0; n<K_samples; n++) {
                    int s0 = ((uint8)raw[2*n+1] << 8) + (uint8)raw[2*n+0];
                    int ss0 = (s0 & 0x7FFF) - (s0 & 0x8000);
                    m_cache[buf][n].channel[0] =
                    m_cache[buf][n].channel[1] = ss0 / 32768.0f;
                }
            } else if ((m_num_channels == 2) && (m_bytes_per_sample == 1)) {
                // 8b stereo
                for (int n=0; n<K_samples; n++) {
                    int s0 = (uint8)raw[2*n+0] - 128;
                    int s1 = (uint8)raw[2*n+1] - 128;
                    m_cache[buf][n].channel[0] = s0 / 128.0f;
                    m_cache[buf][n].channel[1] = s1 / 128.0f;
                }
            } else if ((m_num_channels == 2) && (m_bytes_per_sample == 2)) {
                // 16b stereo
                for (int n=0; n<K_samples; n++) {
                    int s0 = ((uint8)raw[4*n+1] << 8) + (uint8)raw[4*n+0];
                    int ss0 = (s0 & 0x7FFF) - (s0 & 0x8000);
                    int s1 = ((uint8)raw[4*n+3] << 8) + (uint8)raw[4*n+2];
                    int ss1 = (s1 & 0x7FFF) - (s1 & 0x8000);
                    m_cache[buf][n].channel[0] = ss0 / 32768.0f;
                    m_cache[buf][n].channel[1] = ss1 / 32768.0f;
                }
            }

            // check if the final block is a partial block
            if (block_last > m_data_end) {
#if 0
                cout << "Filling trailing end with 0.0\n";
#endif
                for (int n=m_data_end - block_first; n<K_samples; n++) {
                    m_cache[buf][n].channel[0] =
                    m_cache[buf][n].channel[1] = 0.0f;
                }
            }
        }
    }

    return m_cache[buf][buf_offset];
}


void
ReadWav::SetOffset(long off)
{
    assert(off >= 0);
    m_seq_offset = off;
}


samplePair_t
ReadWav::NextFrame()
{
    return Sample(m_seq_offset++);
}


// ========================================================================
// simple peak detector
// ========================================================================

// item in the peak list
struct peak_t {
    int   channel;      // 0=left, 1=right
    float maxima_orig;  // value at peak, including sign
    float maxima;       // value at peak, absolute value
    long  sample_orig;  // sample where local maxima ocurred in wav file
    long  sample;       // local maxima location after time shift
};

typedef vector<peak_t> peakvec_t;

class PeakDet
{
public:
    PeakDet(ReadWav &wavobj, int channel, float min_thold);
    ~PeakDet();

    // feed another value
    void Stuff(float v);

    // indicate no more values are coming
    void Flush();

    // return the collection of edges
    peakvec_t& PeakDet::GetList();

private:
    const int    m_channel;     // 0=left, 1=right
    const float  m_min_peak;    // peaks must be at least this tall

    float m_window[3];          // window of three most recent values
    long  m_samples;            // number of samples received
    static long m_prev_sample;  // shared by both channels

    peakvec_t peaks;            // list of detected peaks

    ReadWav& m_wavobj;          // associated wav file
};


// static member needs explicit initialization
long PeakDet::m_prev_sample = 0;

PeakDet::PeakDet(ReadWav &wavobj, int channel, float min_thold) :
    m_wavobj(wavobj),
    m_channel(channel),
    m_min_peak(min_thold),
    m_samples(0)
{
    // nothing
}


PeakDet::~PeakDet()
{
    // nothing
}


// as each new sample is read in, check to see if it is a local maxima,
// and log it if it eppears to be so.
void
PeakDet::Stuff(float v)
{
    // roll window
    m_window[0] = m_window[1];
    m_window[1] = m_window[2];
    m_window[2] = ABS(v);
    m_samples++;

    // can't do anything until the window is filled
    if (m_samples < 2) {
        return;
    }

    bool interesting = false && (m_samples-2 >= 438614 && m_samples-2 <= 438614 && m_channel == 1);
    if (interesting) cout << "@" << m_samples-2 << ", v=" << m_window[1] << endl;

    // fast test to discard samples that obviously aren't maxima.
    // we also discard any peaks that have too small of an amplitude.
    if ( (m_window[1] < m_window[0]) ||
         (m_window[1] < m_window[2]) ||
         (m_window[1] < m_min_peak) )
        return;

    int this_sample = m_samples - 2;    // -1 because we incremented m_samples,
                                        // and -1 because we are looking back one item in the window

    // do a more thorough job than the quick test above.
    // try and ensure we don't have a local bump in the middle of a larger trend.
    // search +/- 6/16ths of a bit period.
    int off_delta = (6*g_bit_per + 8) >> 4;
    int earlier = 0;
    for (int offset = -off_delta; offset <= off_delta; offset++) {
        float s = m_wavobj.Sample(this_sample + offset).channel[m_channel];
        s = ABS(s);
        if (m_window[1] < s) {
            if (opt_v >= 4) {
                cout << "    rejecting fake peak " << m_window[1] <<
                        " @" << this_sample << ", channel " << m_channel << endl;
	    }
            return;
        } else {
            // if two points have the same value, ignore the later one
            if ((offset < 0) && (m_window[1] == s)) {
                earlier = offset;
	    }
        }
    }
    if (earlier) {
        if (opt_v >= 4) {
            cout << "    rejecting redundant peak @" << this_sample <<
                    ", earlier @" << this_sample+earlier <<
                    ", channel " << m_channel << endl;
        }
        return;
    }

    // use a wider window and make sure that there is a significant dip
    // to both sides of the peak.  experience has shown that if there is
    // no quick flux reversal on the track, it can take a while for the
    // signal to droop.  thus we use the wider window, and we are less
    // stringent about the test on the samples after the peak.
    // search +/- 12/16ths of a bit period.
    int off_delt = (12*g_bit_per + 8) >> 4;
    bool higherL = false;
    bool higherR = false;
    for (int offset = -off_delt; offset <= off_delt; offset++) {
        float s = m_wavobj.Sample(this_sample + offset).channel[m_channel];
        s = ABS(s);
        if (interesting) {
            cout << "    @(" << this_sample << " + " << offset << "), v=" << s;
	}
        // make sure we aren't just stuck in a flat spot
        higherL |= (offset < 0) &&
                   (m_window[1] > 1.12f*s) && (m_window[1] > s + 0.07f);
        higherR |= (offset > 0) &&
                   (m_window[1] > 1.08f*s) && (m_window[1] > s + 0.05f);
        if (interesting) {
            cout << " hl=" << higherL << ", hr=" << higherR << endl;
	}
    }
    if (!higherL || !higherR) {
        if (opt_v >= 4) {
            cout << "    rejecting flat peak " << m_window[1] <<
                    " @" << this_sample << ", channel " << m_channel << endl;
	}
        return;
    }

    peak_t p;
    p.channel     = m_channel;
    p.sample_orig = this_sample;
    p.sample      = this_sample;        // no shift yet
    p.maxima_orig = m_wavobj.Sample(this_sample).channel[m_channel];
    p.maxima      = m_window[1];

    // sanity check: make sure the signs of peaks alternates
    if (!peaks.empty()) {
        bool prev_neg = (peaks.back().maxima_orig < 0.0f);
        bool curr_neg =            (p.maxima_orig < 0.0f);
        if ((opt_v >= 3) && (prev_neg == curr_neg)) {
            cout << "Adjacent peaks with like signs, channel " << m_channel <<
                    ", @" << peaks.back().sample << " and " << p.sample << endl;
	}
    }

    // we have a new peak to add
    peaks.push_back(p);

    if (opt_v >= 3) {
        char name = (m_channel) ? 'R' : 'L';
        cout << name << " peak: @" << dec << this_sample << ", delta=" << this_sample-m_prev_sample << ", " << fixed << p.maxima << endl;
    }

    // this information is shared between the two channels for reporting
    m_prev_sample = this_sample;
}


void
PeakDet::Flush()
{
    // nothing right now
    // let's assume that a peak doesn't occur on the last sample of the stream
}


// return ref to list; caller must make a copy or there will be hell to pay
peakvec_t&
PeakDet::GetList()
{
    return peaks;
}


// ========================================================================
// utility class to compute average and std deviation of a stream of #'s
// ========================================================================

class Stat
{
public:
    Stat();
    void Clear();               // wipe out any accumulated stats
    void Input(float v);        // feed a new term in the series
    float Mean() const;         // get arithmetic mean
    float StdDev() const;       // get standard deviation
    int SeqLength() const;      // return # of samples in sequence
private:
    long  count;        // number of items
    double sum;         // sum of all input values
    double sqsum;       // sum of (each input value squared)

    double M;
    double Q;
};

Stat::Stat()
{
    Clear();
}

void
Stat::Clear()
{
    count = 0;
    sum = 0;
    sqsum = 0;
}

void
Stat::Input(float v)
{
    count++;

    // simple but numerically unstable
    sum += v;
    sqsum += v*v;

    // a bit more complicated, but well behaved
    // see: http://www.cs.berkeley.edu/~mhoemmen/cs194/Tutorials/variance.pdf
    if (count == 1) {
        M = v;
        Q = 0;
    } else {
        double prev_M = M;
        double prev_Q = Q;
        M = prev_M + (v - prev_M) / count;
        Q = prev_Q + ((count-1)*(v-prev_M)*(v-prev_M)) / count;
    }
}


float
Stat::Mean() const
{
    if (count == 0) {
        return 0;
    }
    return float(sum / count);
}


float
Stat::StdDev() const
{
    if (count == 0) {
        return 0;
    }
#if 0
    return sqrt( sqrt((sqsum - sum*sum/count) / count) );
#else
    return sqrt( float(Q) / count );
#endif
}


int
Stat::SeqLength() const
{
    return count;
}


// =========================================================================
// interleave two peak vectors into a single time ordered vector
// =========================================================================

// zip together the two peaks in time order
peakvec_t
Merge(peakvec_t &p0, peakvec_t &p1)
{
    peakvec_t::iterator ptr0(p0.begin());
    peakvec_t::iterator ptr1(p1.begin());
    peakvec_t rslt;

    if (opt_v >= 1) {
        cout << "Merging..." << endl;
    }

    // merge until one list or the other has been fully scanned
    while ( (ptr0 < p0.end()) && (ptr1 < p1.end()) ) {
        if (ptr0->sample <= ptr1->sample) {
            rslt.push_back(*ptr0++);
        } else {
            rslt.push_back(*ptr1++);
        }
    }

    // just push the remainder from the non-empty list
    while (ptr0 < p0.end()) {
        rslt.push_back(*ptr0++);
    }
    while (ptr1 < p1.end()) {
        rslt.push_back(*ptr1++);
    }

    if (opt_v >= 1) {
        cout << "    merged list has " << rslt.size() << " edges" << endl;
    }

    return rslt;
}


// =========================================================================
// Filtering a list of peaks
// =========================================================================

// Filtering the block means
//     taking stats on the peaks
//     shifting the right track to make L->R and R->L similar
//     tossing out unlikely peaks
enum { LL=0, LR=1, RL=2, RR=3 }; // transition classification
void
Filter(const int blk_num, peakvec_t &peaks, Stat *stats)
{
    if (peaks.size() < 32) {
        return;
    }

    // run through the peaks, calculating the mean and standard deviation
    // of the amplitude of the peaks, and the L->L, L->R, R->L, and R->R
    // transition intervales.  It is possible that not all transitions
    // occur in a given block.
    Stat amp;
    for (peakvec_t::iterator i=peaks.begin(); i < peaks.end(); i++) {
        amp.Input(i->maxima);
        if (i > peaks.begin()) {      // delta is relative to previous sample
            int xition = 2*((i-1)->channel) + (i->channel);
            int tdelta = i->sample - (i-1)->sample;
            stats[xition].Input(float(tdelta));
        }
    }

    if (opt_v >= 1) {
        cout << "    " << amp.SeqLength() << " peaks" <<
                ", avg=" << amp.Mean() <<
                ", std dev=" << amp.StdDev() << endl;
        cout << "    " << stats[LL].SeqLength() << " LL" <<
                ", avg=" << stats[LL].Mean() <<
                ", std dev=" << stats[LL].StdDev() << endl;
        cout << "    " << stats[LR].SeqLength() << " LR" <<
                ", avg=" << stats[LR].Mean() <<
                ", std dev=" << stats[LR].StdDev() << endl;
        cout << "    " << stats[RL].SeqLength() << " RL" <<
                ", avg=" << stats[RL].Mean() <<
                ", std dev=" << stats[RL].StdDev() << endl;
        cout << "    " << stats[RR].SeqLength() << " RR" <<
                ", avg=" << stats[RR].Mean() <<
                ", std dev=" << stats[RR].StdDev() << endl;
    }

    // loop through the peaks again, shifting the tracks by
    // an amount that will equalize the LR and RL difference,
    // and toss out peaks that are not tall enough.
    int tshift = int(stats[RL].Mean() - stats[LR].Mean() + 0.5f) >> 1;
    if (opt_v >= 1) {
        cout << "    applying R channel shift of " << tshift << " samples" << endl;
    }

    peakvec_t shifted;
    for (peakvec_t::iterator i=peaks.begin(); i < peaks.end(); i++) {
        peak_t p = *i;
        if (p.channel == 1) {
            p.sample += tshift;
	}
        shifted.push_back(p);
    }

    // rerun the stats on the updated block
    amp.Clear();
    for (int i=0; i<4; i++)
        stats[i].Clear();

    for (peakvec_t::iterator i=shifted.begin(); i < shifted.end(); i++) {
        amp.Input(i->maxima);
        if (i > shifted.begin()) {      // delta is relative to previous sample
            int xition = 2*((i-1)->channel) + (i->channel);
            int tdelta = i->sample - (i-1)->sample;
            stats[xition].Input(float(tdelta));
        }
    }

    for (int i=0; i<4; i++) {
        if ( (g_bit_per < stats[i].Mean() - g_devs_warn*stats[i].StdDev()) ||
             (g_bit_per > stats[i].Mean() + g_devs_warn*stats[i].StdDev()) ) {
            if (stats[i].SeqLength() > 4) // guard against inadequate stats
                cerr << "Warning: block #" << blk_num <<
                        ", apriori samples/bit = " << g_bit_per <<
                        ", measured samples/bit = " << stats[i].Mean() << endl;
        }
    }

    // loop through the shifted peaks again, tossing out any implausible
    // edges.  an edge is implausible if it isn't within three std dev (or so)
    // of the normal transition time and the following edge is.

    int lower[4], upper[4];
    for (int i=0; i<4; i++) {
        lower[i] = int(stats[i].Mean() - g_devs_filter*stats[i].StdDev() + 0.5f);
        upper[i] = int(stats[i].Mean() + g_devs_filter*stats[i].StdDev() + 0.5f);
    }

    peakvec_t rslt;
    for (peakvec_t::iterator i=shifted.begin(); i < shifted.end(); i++) {
        if (!rslt.empty()) {
            peak_t prev = rslt.back();
            int xition = 2*((i-1)->channel) + (i->channel);
            if ( (i != shifted.end()-1) &&
                 ((i+0)->sample - prev.sample < lower[xition]) && // this one is too close
                 ((i+1)->sample - prev.sample > lower[xition]) && // but the next one is plausible
                 ((i+1)->sample - prev.sample < upper[xition]) ) {
                if (opt_v >= 3) {
                    cout << "    rejecting tweener peak @" << i->sample_orig <<
                            "; height=" << i->maxima <<
                            ", channel " << i->channel << endl;
                    cout << "        " << prev.sample_orig   << " -> " <<
                                          (i+0)->sample_orig << " -> " <<
                                          (i+1)->sample_orig << endl;
                }
                continue;       // skip this one
            }
        }

        // it seems plausible, so save it
        rslt.push_back(*i);
    }

    // return updated list, along with stats
    peaks = rslt;
}


// =========================================================================
// decode edge stream
// =========================================================================

// described a decoded block of bytes
struct datablk_t {
    int number;         // block number of tape
    int warnings;       // count of how many warnings occurred in decoding
    int dribble_bits;   // number of bits after last decoded byte
    int byte_count;     // number of decoded bytes
    uint8 *data;

    // these fields contain decoded data if the tape is 3300 format
    int record_length;
    int file_number;
    int record_number;
    int load_addr;
    bool cksum_ok;
    bool ctl_cksum_ok;
};

// string blocks together into a file or whole tape
typedef vector<datablk_t> blockvec_t;

// the decoder operates as a state machine
datablk_t
DecodeBlk(const int blk_num, peakvec_t& peaks, const Stat * const stats)
{
    int warnings = 0;   // number of fishy transitions
    int byte_cnt = 0;   // byte of block counter
    int bit_cnt  = 0;   // bit of byte counter
    uint8 byte;         // byte being assembled
    uint8 data[1024];   // blocks never get this big, but just in case...

    if (opt_v >= 1) {
        cout << "Starting block decode at frame " << peaks.begin()->sample_orig << "\n";
    }

    if (peaks.size() < 32) {
        cout << "    skipping: only " << dec << peaks.size() << " edges" << endl;
        datablk_t d;
        d.number = blk_num;
        d.warnings = 1;
        d.dribble_bits = 0;
        d.byte_count = 0;
        d.data = NULL;
        return d;
    }

    // set bounds for reasonable transition times
    int lower[4], upper[4];
    for (int i=0; i<4; i++) {
        lower[i] = int(stats[i].Mean() - g_devs_warn*stats[i].StdDev() + 0.5f);
        upper[i] = int(stats[i].Mean() + g_devs_warn*stats[i].StdDev() + 0.5f);
    }

    for (peakvec_t::iterator pk = peaks.begin(); pk < peaks.end(); pk++) {

        bool funny_interval = false;

        // handle first bit of byte
        if (bit_cnt == 0) {
            byte = (pk->channel) ? 0x01 : 0x00;
            bit_cnt = 1;
            continue;
        }

        // this is not the first bit of a byte
        int xition   = 2*((pk-1)->channel) + (pk->channel);
        int interval = (pk->sample - (pk-1)->sample);

        if (interval < lower[xition]) {
            // interval is too short
            if (opt_v >= 3) {
                cout << "    found runt interval " << dec << interval <<
                        " at frame " << pk->sample_orig << "\n";
                funny_interval = true;
            }
        }
        else if (interval > upper[xition]) {
            // interval is too long
            if (opt_v >= 3) {
                cout << "    found long interval " << dec << interval <<
                        " at frame " << pk->sample_orig << "\n";
                funny_interval = true;
            }
        }
        // even if it might be too short or too long, we decode it anyway

        if (opt_v >= 3) {
            cout << "    got bit " << ((pk->channel) ? 1 : 0) << "\n";
	}
        byte = (byte << 1) | ((pk->channel) ? 1 : 0);
        if (++bit_cnt == 8) {
            if (opt_v >= 2) {
                cout << hex << setfill('0') <<
                        "    byte " << setw(3) << byte_cnt <<
                        "= 0x" << setw(2) << (int)byte << "\n";
	    }
            if (byte_cnt < 1024) {
                data[byte_cnt] = byte;
            } else if ((byte_cnt == 768) && (opt_v >= 2)) {
                cout << "    Warning: run on block!\n";
	    }
            byte_cnt++;
            byte = 0x00;
            bit_cnt = 0;
        }

        // I don't really sign this getting triggered anymore ...
        // I think the earlier, more agressive detection of true maxima weeds
        // out the "too soon edges".  I can imagine that the missing edge
        // case could occur, though.
        if (funny_interval) {
            warnings++;
//  cerr << "Detected funny peak @" << pk->sample_orig << endl;
            if (opt_v >= 3) {
                for (int d=-3; d<3; d++) {
                    if ( (pk+d-1 >= peaks.begin()) && (pk+d < peaks.end()) ) {
                        int ts0 = (pk+d  )->sample_orig;
                        int tsp = (pk+d-1)->sample_orig;
                        cout << "        frame " << ts0 <<
                                " delta " << ts0 - tsp <<
                                " amplitude " << (pk+d)->maxima <<
                                " channel "   << (pk+d)->channel;
                        if (d == 0) {
                            cout << "  <<<\n";
                        } else {
                            cout << "\n";
                        }
                    }
                }
            }
        }

    } // pk iterator

    if (bit_cnt != 0) {
        if (opt_v >= 3) {
            cout << "    warning: " << bit_cnt << " residual bits" << endl;
	}
    }

    // return data
    datablk_t d;
    d.number = blk_num;
    d.warnings = warnings;
    d.byte_count = byte_cnt;
    d.dribble_bits = bit_cnt;
    d.data = new uint8[byte_cnt];
    assert(d.data != NULL);
    memmove( (void*)d.data, (void*)data, byte_cnt );
    return d;
}


void
DumpBlk(datablk_t &blk)
{
    cout << "## Block " << dec << blk.number <<
            ", length = " << blk.byte_count << " bytes, " <<
            blk.dribble_bits << " dribble bits, " <<
            blk.warnings << " warnings\n";

    if (blk.byte_count < 10) {
        cout << "## too small to decode\n";
        return;
    }

    if (opt_fmt == FMT_3300) {
        // decode the header
        int cksum         =                   blk.data[1];
        int record_length = 256*blk.data[2] + blk.data[3];
        int ctl_cksum     =                   blk.data[4];
        int file_number   =                   blk.data[5];
        int record_number = 256*blk.data[6] + blk.data[7];
        int load_addr     = 256*blk.data[8] + blk.data[9];

        // control check sum: sum of bytes 0, 2,3, 5..9
        int computed_ctl_cksum =
        computed_ctl_cksum = blk.data[0]
                           + blk.data[2] + blk.data[3]
                           + blk.data[5]
                           + blk.data[6] + blk.data[7]
                           + blk.data[8] + blk.data[9];
        computed_ctl_cksum &= 0xFF;

        // block check sum: sum of bytes 0, 2,3, 5..9, 10..10+length-1
        int computed_cksum = computed_ctl_cksum;
        for (int i=10; i < 10+record_length; i++)
            computed_cksum += blk.data[i];
        computed_cksum &= 0xFF;

        cout << "# block type: " << hex << setfill('0') << setw(2) << (int)blk.data[0] << "\n";

        cout << "# block cksum: ";
        if (cksum == computed_cksum) {
            cout << "OK\n";
        } else {
            cout << "BAD!  0x" << hex << cksum << " expected, 0x" << computed_cksum << " seen\n";
        }

        cout << "# record length = " << dec << record_length << " bytes\n";

        cout << "# control cksum: ";
        if (ctl_cksum == computed_ctl_cksum) {
            cout << "OK\n";
	} else {
            cout << "BAD!  0x" << hex << ctl_cksum << " expected, 0x" << computed_ctl_cksum << " seen\n";
	}

        cout << "# file number   = " << dec << file_number << "\n";
        cout << "# record number = " << dec << record_number << "\n";
        cout << "# load addr     = " << hex << setw(4) << load_addr << "\n";

        // tuck away the decoded header data
        blk.record_length = record_length;
        blk.file_number   = file_number;
        blk.record_number = record_number;
        blk.load_addr     = load_addr;
        blk.cksum_ok      = (cksum == computed_cksum);
        blk.ctl_cksum_ok  = (ctl_cksum == computed_ctl_cksum);
    }

    // report the data
    int start = 0;              // first byte to dump
    int end = blk.byte_count;   // one past end
    int cnt = 0;

    for (int i=start; i<end; i++, cnt++) {

        // print data byte
        cout << hex << setfill('0') << setw(2) << (int)blk.data[i];

        // if at the end of a line of data, print it in ASCII
        if ( (cnt%16 == 15) ||   // end of a line
             (i == end-1) ) {    // end of block
            if (cnt%16 < 15) {
                // pad to make sure ascii part lines up
                for (int k=(cnt+1)%16 ; k<16; k++)
                    cout << "  ";
            }
            cout << "  ";
            for (int j=i-(cnt%16); j <= i; j++) {
                if (blk.data[j] < 32 || blk.data[j] > 126) {
                    cout << '.';
                } else {
                    cout << (char)blk.data[j];
                }
            }
            cout << "\n";
        }
    }
}


// Decode the entire tape stream into a list of data blocks:
//     break the list into blocks, by detecting large gaps
//     clean up the edges for each block
//     send resulting block into block decoder
blockvec_t
DecodeBytes(peakvec_t& peaks)
{
    if (opt_v >= 1) {
        cout << "Decoding..." << endl;
    }

    blockvec_t rslt;    // assembled list of blocks

    int blk_num = 0;
    for (peakvec_t::iterator pk = peaks.begin(); pk < peaks.end(); blk_num++) {

        // list of peaks associated with one block
        peakvec_t blk_peaks;

        // record the start of the block
        peakvec_t::iterator sob = pk;

        // scan ahead looking for a gap, or end of tape.
        // after the while loop, pk points to either the end of list,
        // or the first sample after the gap.
        while (pk < peaks.end()) {
            if ( (pk != sob) &&
                 (pk->sample - (pk-1)->sample >= g_gap_thold) )
                break;
            blk_peaks.push_back(*pk);
            pk++;
        }

        if (opt_v >= 1) {
            peakvec_t::iterator eob = pk-1;
            float t_ms = (eob->sample - sob->sample) / 44.1f;
            cout << "Block #" << dec << blk_num <<
                    ": start @" << sob->sample <<
                    ", end @"   << eob->sample <<
                    " (" << t_ms << " ms)" << endl;
        }

        // clean up edges
        Stat stats[4];
        Filter(blk_num, blk_peaks, stats);

        // decode it the best we can into a block of bytes
        datablk_t data = DecodeBlk(blk_num, blk_peaks, stats);
        DumpBlk(data);          // report it
        rslt.push_back(data);   // save it for later processing
    }

    return rslt;
}


// ========================================================================
// produce output file
// ========================================================================

void
EmitBlk3300(ofstream &os, blockvec_t::iterator blk, string msg="")
{
    const int hdr_bytes = 10;   // number of header bytes in a 3300 block

    os << setfill('0'); // we only use width() when hex is in play

    os << "## Raw tape block #" << dec << blk->number << "\n";
    if (msg.length() > 0) {
        os << "#  " << msg << endl;
    }
    if (!blk->ctl_cksum_ok || !blk->cksum_ok) {
        os << "#  Bad checksum -- the entire block is untrustworthy\n";
    }
    os << "#  file number:   " << blk->file_number << "\n";
    os << "#  record number: " << blk->record_number << "\n";
    os << "#  record length: " << blk->record_length << "\n";
    os << "#  load address:  " << hex << setw(4) << blk->load_addr << "\n";
    if (blk->warnings > 0 || blk->dribble_bits > 0) {
        os << "#  data bytes:    " << dec << blk->byte_count - hdr_bytes << "\n";
    }

    // emit in intel hex format, kind of
    bool hdr = false;                           // haven't emitted header
    int end = blk->load_addr + blk->byte_count; // one past last byte
    uint8 sum;

    os << hex;
    for (int i=hdr_bytes; i< blk->byte_count; i++) {
        int a = blk->load_addr + i - hdr_bytes;

        if (!hdr) {
            int len = (end-a > 16) ? (16-a%16) : end-a;
            os << ":";
            os << setw(2) << len; // record length
            os << setw(4) << a;   // address
            os << "00";           // data record
            hdr = true;
            sum = len + (a >> 8) + (a & 0xFF) + 0x00;
        }

        // emit the next data byte
        os << setw(2) << (int)blk->data[i];
        sum += blk->data[a];

        // end of line?
        if ((a%16 == 15) || (a == end-1)) {
            hdr = false;
            sum = -sum;
            os << setw(2) << (int)sum << "\n";
        }
    }

    os << endl; // emit a blank line and flush
}


// given a list of blocks, try to sort out the redundant pairs,
// emit just the correct one, and note any problems on the bad ones.
// TODO: check tape file structure.
void
EmitTape3300(blockvec_t rawblocks)
{
    // open output file
    ofstream os(opt_ofn.c_str(), ios::out | ios::trunc);
    if (!os.is_open()) {
        cerr << "Error: couldn't create file '" << opt_ofn << "'" << endl;
        exit(-1);
    }

    // try and match up redundant pairs of data blocks, even with errors.
    // (a) if there is a single block, just report it.
    // (b) if the first block matches the second block,
    //     smile, report one of the blocks.
    // (c) if the first block and second block have good headers and they
    //     match, pick the one with the correct checksum.
    // (d) report the block containing more data

    for (blockvec_t::iterator blk = rawblocks.begin(); blk < rawblocks.end(); ) {
        blockvec_t::iterator nxt = blk+1;

        // case (a) -- single, unmatched block
        if (blk == rawblocks.end()-1) {
            string msg("Unmatched block");
            EmitBlk3300(os, blk, msg);
            blk++;
            continue;
        }

        // check if all the header records match
        bool hdr_match = (blk->record_length == nxt->record_length) &&
                         (blk->file_number   == nxt->file_number) &&
                         (blk->record_number == nxt->record_number) &&
                         (blk->load_addr     == nxt->load_addr);

        // cases (b) and (c) -- headers match
        if (blk->ctl_cksum_ok && nxt->ctl_cksum_ok && hdr_match) {
            if (blk->cksum_ok) {
                EmitBlk3300(os, blk);
            } else if (nxt->cksum_ok) {
                EmitBlk3300(os, nxt);
            } else {
                // both bad; report the one with more data
                // it is arguable that the criterion should be # warnings,
                // or # of dribble bits.
                if (blk->record_length >= nxt->record_length) {
                    EmitBlk3300(os, blk);
                } else {
                    EmitBlk3300(os, nxt);
                }
            }
            blk += 2;
            continue;
        }

        // cases (d) -- report the block as a singleton
        string msg("Unmatched block");
        EmitBlk3300(os, blk, msg);
        blk++;
        continue;

    } // for (blk)

    os.close();
}


void
EmitBlk2200(ofstream &os, blockvec_t::iterator blk, string msg="")
{
    os << setfill('0'); // we only use width() when hex is in play

    os << "## Raw tape block #" << dec << blk->number << "\n";
    if (msg.length() > 0) {
        os << "#  " << msg << endl;
    }
    if (blk->byte_count != 256 || blk->warnings != 0 || blk->dribble_bits != 0) {
        os << "#  Warning bad byte count, " <<
              blk->warnings << " warnings," <<
              blk->dribble_bits << " dribble bits\n";
    }
    os << "#  data bytes:    " << dec << blk->byte_count << "\n";

    // emit in intel hex format, kind of
    int end = blk->byte_count;   // one past last byte

    for (int i=0; i < blk->byte_count; i++) {
        os << hex << setw(2) << (int)blk->data[i];
        if ((i%16 == 15) || (i == blk->byte_count-1)) {
            os << "\n";
	}
    }

    os << endl; // emit a blank line and flush
}


// given a list of blocks, try to sort out the redundant pairs,
// emit just the correct one, and note any problems on the bad ones.
// TODO: check tape file structure.
void
EmitTape2200(blockvec_t rawblocks)
{
    // open output file
    ofstream os(opt_ofn.c_str(), ios::out | ios::trunc);
    if (!os.is_open()) {
        cerr << "Error: couldn't create file '" << opt_ofn << "'" << endl;
        exit(-1);
    }

    // try and match up redundant pairs of data blocks, even with errors.
    // (a) if the block is < 10 bytes, it is likely just a burst, and ignore it
    // (b) if there is a single block, just report it.
    // (c) if the first block matches the second block,
    //     smile, report one of the blocks.
    // (d) on the 2200, data blocks always have 256 bytes exactly.
    //     if one block has 256 bytes and no dribble bits, use that block.
    // (e) pick one based on some ad-hoc criteria and note the problem
    //
    // note: the 2200 format on tape appears to be  like this:
    //    (1) burst of 66 '1' bits
    //    (2) ~30 bit gap
    //    (3) 2048 edges of the first copy of the data
    //    (4) ~20 bit gap
    //    (5) 2048 edges of the second copy of the data
    //    (6) ~400 bit gap
    // thus we can use the short burst block to help group pairs of blocks.

    bool burst = false; // indicates that previous block was a burst

    for (blockvec_t::iterator blk = rawblocks.begin(); blk < rawblocks.end(); ) {
        blockvec_t::iterator nxt = blk+1;

        // case (a) -- burst of 1's
        if (blk->byte_count < 10) {
            burst = true;
            blk++;
            continue;
        }

        // case (b) -- single, unmatched block
        if (blk == rawblocks.end()-1) {
            string msg("Unmatched block");
            EmitBlk2200(os, blk, msg);
            burst = false;
            blk++;
            continue;
        }

        bool blk_ok = (blk->byte_count == 256) && (blk->dribble_bits == 0);
        bool nxt_ok = (nxt->byte_count == 256) && (nxt->dribble_bits == 0);
        bool data_match = blk_ok && nxt_ok
                       && (memcmp( blk->data, nxt->data, 256 ) == 0);

        // (c) the first block matches the second block
        if (blk_ok && nxt_ok && data_match) {
            EmitBlk2200(os, blk);
            blk += 2;
            burst = false;
            continue;
        }

        // both blocks appear valid at first blush, so just pick one
        if (blk_ok && nxt_ok && !data_match) {
            string msg("Warning: block pair didn't match!  See the log file.");
            EmitBlk2200(os, blk, msg);
            blk += 2;
            burst = false;
            continue;
        }

        // (d) if we are just after a burst, we have a pair,
        //     so pick the better one and skip both.
        if (burst) {
            if (blk_ok) {
                EmitBlk2200(os, blk);
            } else if (nxt_ok) {
                EmitBlk2200(os, nxt);
            } else if (blk->byte_count >= nxt->byte_count) {
                EmitBlk2200(os, blk);
            } else {
                EmitBlk2200(os, nxt);
	    }
            blk += 2;
            burst = false;
            continue;
        }

        // cases (e) -- report the block as a singleton
        string msg("Unmatched block");
        EmitBlk2200(os, blk, msg);
        blk++;
        continue;

    } // for (blk)

    os.close();
}


// =========================================================================
// main
// =========================================================================

void
usage(int code)
{
    cerr << "Usage: WangWav2Tape [<options> ...] <inname>\n";
    cerr << "  no -v         report only errors\n";
    cerr << "  -v            same as '-v 1'\n";
    cerr << "  -v 1          ... and warnings and block information\n";
    cerr << "  -v 2          ... and each byte decoded\n";
    cerr << "  -v 3          ... and each bit decoded\n";
    cerr << "  -2200         assume tape was produced by a Wang 2200 (default)\n";
    cerr << "  -3300         assume tape was produced by a Wang 3300\n";
    cerr << "  -o <outname>  specifies output filename;\n";
    cerr << "Version: January 24, 2008\n";
    exit(code);
}


// parse the command line arguments
void
ParseArgs(int argc, char **argv)
{
    // set command line defaults
    opt_v   = 0;
    opt_ofn = "";
    opt_fmt = FMT_2200;

#if 1
    if (argc < 2) {     // we need at least one parameter
        usage(-1);
    }
#endif

    for (int i=1; i<argc-1; i++) {

        if (strcmp(argv[i], "-v") == 0) {
            // verbose reporting
            opt_v = 1;  // at least
            if (i+1 < argc-1) {
                // check for optional number
                opt_v = atol(argv[i+1]);
                if (opt_v == 0) {
                    // if next arg wasn't a number, atol returns 0
                    // however, if the user did "-v 0", we mess up in this assumption
                    opt_v = 1;
                } else {
                    i++;        // skip numeric argument
                }
            }
        }

        else if (strcmp(argv[i], "-2200") == 0) {
            opt_fmt = FMT_2200;
        }

        else if (strcmp(argv[i], "-3300") == 0) {
            opt_fmt = FMT_3300;
        }

        else if (strcmp(argv[i], "-o") == 0) {
            if (i+1 < argc-1) {
                opt_ofn = argv[i+1];
                i++;
            } else {
                cerr << "Error: final argument is input filename, not part of -o specification\n";
                usage(-1);
            }
        }

        else {
            cerr << "Error: unrecognized option '" << argv[i] << "'\n";
            usage(-1);
        }
    }

    // one last parameter -- it must be the original wav file
    opt_ifn = argv[argc-1];

    // if -o wasn't specified, derive the output filename from the input filename
    if (opt_ofn.length() == 0) {
        // make sure the suffix is ".wav"
        int len = opt_ifn.length();
        if ((len > 4) && (opt_ifn.compare(len-4, 4, ".wav") ||
                          opt_ifn.compare(len-4, 4, ".WAV")) ) {
            // we just change the suffix of the input filename
            opt_ofn = opt_ifn.substr(0, len-4);
        } else {
            // don't understand suffix, so just tack on new suffix
            opt_ofn = opt_ifn;
        }
        opt_ofn = opt_ofn + ".hex";
    }
}


// ========================================================================
// main logic
// ========================================================================

int
main(int argc, char **argv)
{
    // parse command line arguments
    ParseArgs(argc, argv);

    // establish an interface to the wav file
    ReadWav wavobj(opt_ifn);

    // computed some dependent parameters
    // FIXME: allow command line override of g_bit_per
    switch (opt_fmt) {
        case FMT_3300:
            g_bit_per = (int)( 44           // measured as typical @ 44.1KHz
                               * wavobj.SampleFrequency()
                               / 44100.0f
                                + 0.5);     // rounding
            // the gap between blocks is pretty large ...
            // hundreds to thousands of bit times.
            g_gap_thold = 20*g_bit_per;
            // FIXME: empirical for the single tape I have
            opt_min_thold = 0.20f;
            break;
        case FMT_2200:
            g_bit_per = (int)( 30           // measured as typical @ 44.1KHz
                               * wavobj.SampleFrequency()
                               / 44100.0f
                                + 0.5);     // rounding
            // on the 2200, the gap between the redundant copies of
            // the block is really short -- under 20 bit times.
            // a gap is 20 bit periods of silence, or more
            g_gap_thold = 10*g_bit_per;
            // FIXME: empirical
            opt_min_thold = 0.12f;
            break;
        default:
            assert(0);
    }

    // scan the wav file and find all peaks
    PeakDet peakDet0(wavobj, 0, opt_min_thold);
    PeakDet peakDet1(wavobj, 1, opt_min_thold);

    if (opt_v >= 1) {
        cout << "Scanning wav file..." << endl;
    }

    const long len = wavobj.NumFrames();
    for (long c=0; c<len; c++) {
        samplePair_t s = wavobj.NextFrame();
        //cout << dec << "Sample " << c << ": L=" << fixed << s.channel[0] << endl;
        //cout << dec << "Sample " << c << ": R=" << fixed << s.channel[1] << endl;
        peakDet0.Stuff(s.channel[0]);
        peakDet1.Stuff(s.channel[1]);
    }

    peakDet0.Flush();
    peakDet1.Flush();

    peakvec_t peaks[2];
    peaks[0] = peakDet0.GetList();
    peaks[1] = peakDet1.GetList();

    if (opt_v >= 1) {
        cout << "    " << peaks[0].size() << " left peaks\n";
        cout << "    " << peaks[1].size() << " right peaks\n";
    }

    // zip the peaks of the two tracks together
    peakvec_t merged = Merge(peaks[0], peaks[1]);

    // decode the tape
    blockvec_t tape = DecodeBytes(merged);

    // generate final output
    if (opt_fmt == FMT_3300) {
        EmitTape3300(tape);
    } else {
        EmitTape2200(tape);
    }

    return 0;
}
