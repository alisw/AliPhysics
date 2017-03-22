// Streams.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Classes to implement reading from or writing to gzipped files.
// Adapted for Sherpa by Frank Siegert from:
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
// (http://www.cs.unc.edu/Research/compgeom/gzstream).
// Further adapted to PYTHIA by Stefan Prestel.

#ifndef Pythia8_Streams_H
#define Pythia8_Streams_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string.h>

namespace Pythia8 {

#ifdef GZIPSUPPORT
#include <zlib.h>
//==========================================================================

// Internal classes to implement gzstream. See below for user classes.

// -------------------------------------------------------------------------

class gzstreambuf : public std::streambuf {
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    gzFile           file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    char             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    gzstreambuf() : opened(0) {
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position
        // ASSERT: both input & output capabilities will not be used together
    }
    int is_open() { return opened; }
    gzstreambuf* open( const char* name, int open_mode);
    gzstreambuf* close();
    ~gzstreambuf() { close(); }

    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};

// -------------------------------------------------------------------------

class gzstreambase : virtual public std::ios {
protected:
    gzstreambuf buf;
public:
    gzstreambase() { init(&buf); }
    gzstreambase( const char* name, int open_mode);
    ~gzstreambase();
    void open( const char* name, int open_mode);
    void close();
    gzstreambuf* rdbuf() { return &buf; }
};

//==========================================================================

// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz*
// function interface of the zlib. Files are compatible with gzip compression.

// -------------------------------------------------------------------------

class igzstream : public gzstreambase, public std::istream {
public:
    igzstream() : std::istream( &buf) {}
    igzstream( const char* name, int mode = std::ios::in)
        : gzstreambase( name, mode), std::istream( &buf) {}
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int mode = std::ios::in) {
        gzstreambase::open( name, mode);
    }
};
// -------------------------------------------------------------------------

class ogzstream : public gzstreambase, public std::ostream {
public:
    ogzstream() : std::ostream( &buf) {}
    ogzstream( const char* name, int mode = std::ios::out)
        : gzstreambase( name, mode), std::ostream( &buf) {}
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int mode = std::ios::out) {
        gzstreambase::open( name, mode);
    }
};
#else
typedef std::ifstream igzstream;
typedef std::ofstream ogzstream;
#endif

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Streams_H
