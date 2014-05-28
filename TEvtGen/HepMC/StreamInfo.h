//--------------------------------------------------------------------------
#ifndef HEPMC_STREAM_INFO_H
#define HEPMC_STREAM_INFO_H

//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, March 2009
//
// This class contains the extra information needed when using streaming IO
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "HepMC/Units.h"

namespace HepMC {

/// The known_io enum is used to track which type of input is being read
enum known_io { gen=1, ascii, extascii, ascii_pdt, extascii_pdt };

//! StreamInfo contains extra information needed when using streaming IO.

///
/// \class  StreamInfo
/// This class contains the extra information needed when using streaming IO
/// to process HepMC GenEvents
///
class StreamInfo {
public:
    /// default constructor
    StreamInfo( );
    /// destructor
    ~StreamInfo() {}

    /// IO_GenEvent begin event block key
    std::string IO_GenEvent_Key()          const { return m_io_genevent_start; }
    /// IO_GenEvent end event block key
    std::string IO_GenEvent_End()          const { return m_io_genevent_end; }

    /// IO_Ascii begin event block key
    /// IO_Ascii has been removed, but we want to be able to read 
    /// existing files written by IO_Ascii
    std::string IO_Ascii_Key()             const { return m_io_ascii_start; }
    /// IO_Ascii end event block key
    std::string IO_Ascii_End()             const { return m_io_ascii_end; }
    /// IO_Ascii begin particle data block key
    std::string IO_Ascii_PDT_Key()             const { return m_io_ascii_pdt_start; }
    /// IO_Ascii end particle data block key
    std::string IO_Ascii_PDT_End()             const { return m_io_ascii_pdt_end; }

    /// IO_ExtendedAscii begin event block key
    /// IO_ExtendedAscii has been removed, but we want to be able to read 
    /// existing files written by IO_ExtendedAscii
    std::string IO_ExtendedAscii_Key()     const { return m_io_extendedascii_start; }
    /// IO_ExtendedAscii end event block key
    std::string IO_ExtendedAscii_End()     const { return m_io_extendedascii_end; }
    /// IO_ExtendedAscii begin particle data block key
    std::string IO_ExtendedAscii_PDT_Key()             const { return m_io_extendedascii_pdt_start; }
    /// IO_ExtendedAscii end particle data block key
    std::string IO_ExtendedAscii_PDT_End()             const { return m_io_extendedascii_pdt_end; }

    /// get IO type
    int io_type() const { return m_io_type; }
    /// set IO type
    void set_io_type( int );

    /// true if the stream has a file type key
    /// has_key is true by default
    bool has_key() const { return m_has_key; }
    /// set to false if the stream does not have a file type key
    void set_has_key( bool );
    
    /// get the I/O momentum units
    Units::MomentumUnit io_momentum_unit() const { return m_io_momentum_unit; }
    /// get the I/O length units
    Units::LengthUnit io_position_unit() const { return m_io_position_unit; }

    /// get the I/O stream id
    /// This is used for sanity checking.
    int stream_id() const { return m_stream_id; }
    
    /// Special information is processed the first time we use the IO
    bool finished_first_event() const { return m_finished_first_event_io; }
    /// Special information is processed the first time we use the IO
    void set_finished_first_event( bool b ) { m_finished_first_event_io = b; }

    /// needed when reading a file without units if those units are 
    /// different than the declared default units 
    /// (e.g., the default units are MeV, but the file was written with GeV)
    /// This method is not necessary if the units are written in the file
    void use_input_units( Units::MomentumUnit, Units::LengthUnit );
    
    /// reading_event_header will return true when streaming input is 
    /// processing the GenEvent header information
    bool reading_event_header();
    /// set the reading_event_header flag
    void set_reading_event_header(bool);

private: // data members
    bool        m_finished_first_event_io;
    // GenEvent I/O method keys
    std::string m_io_genevent_start;
    std::string m_io_ascii_start;
    std::string m_io_extendedascii_start;
    std::string m_io_genevent_end;
    std::string m_io_ascii_end;
    std::string m_io_extendedascii_end;
    // particle data I/O method keys
    std::string m_io_ascii_pdt_start;
    std::string m_io_extendedascii_pdt_start;
    std::string m_io_ascii_pdt_end;
    std::string m_io_extendedascii_pdt_end;
    // io information
    int         m_io_type;
    bool        m_has_key;
    // default io units - used only when reading a file with no units
    Units::MomentumUnit m_io_momentum_unit;
    Units::LengthUnit   m_io_position_unit;
    // used to keep identify the I/O stream
    unsigned int m_stream_id;
    static unsigned int m_stream_counter;
    // used to keep track when reading event
    bool m_reading_event_header;

};

} // HepMC

#endif  // HEPMC_STREAM_INFO_H
//--------------------------------------------------------------------------
