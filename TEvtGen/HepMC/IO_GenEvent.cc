//--------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, July 2006
// event input/output in ascii format for machine reading
// IO_GenEvent format contains HeavyIon and PdfInfo classes
//////////////////////////////////////////////////////////////////////////

#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_Exception.h"
#include "HepMC/GenEvent.h"
#include "HepMC/StreamHelpers.h"

namespace HepMC {

    IO_GenEvent::IO_GenEvent( const std::string& filename, std::ios::openmode mode ) 
    : m_mode(mode), 
      m_file(filename.c_str(), mode), 
      m_ostr(0),
      m_istr(0),
      m_iostr(0),
      m_have_file(false),
      m_error_type(IO_Exception::OK),
      m_error_message()
    {
	if ( (m_mode&std::ios::out && m_mode&std::ios::in) ||
	     (m_mode&std::ios::app && m_mode&std::ios::in) ) {
            m_error_type = IO_Exception::InputAndOutput;
	    m_error_message ="IO_GenEvent::IO_GenEvent Error, open of file requested of input AND output type. Not allowed. Closing file.";
	    std::cerr << m_error_message << std::endl;
	    m_file.close();
	    return;
	}
	// now we set the streams
	m_iostr = &m_file;
	if ( m_mode&std::ios::in ) {
	    m_istr = &m_file;
	    m_ostr = NULL;
	    detail::establish_input_stream_info(m_file);
	}
	if ( m_mode&std::ios::out ) {
	    m_ostr = &m_file;
	    m_istr = NULL;
	    detail::establish_output_stream_info(m_file);
	}
	m_have_file = true;
    }


    IO_GenEvent::IO_GenEvent( std::istream & istr ) 
    : m_ostr(0),
      m_istr(&istr),
      m_iostr(&istr),
      m_have_file(false),
      m_error_type(IO_Exception::OK),
      m_error_message()
    { 
        detail::establish_input_stream_info( istr );
    }

    IO_GenEvent::IO_GenEvent( std::ostream & ostr )
    : m_ostr(&ostr),
      m_istr(0),
      m_iostr(&ostr),
      m_have_file(false),
      m_error_type(IO_Exception::OK),
      m_error_message()
   {
        detail::establish_output_stream_info( ostr );
   }

    IO_GenEvent::~IO_GenEvent() {
    	if ( m_ostr != NULL ) {
	    write_HepMC_IO_block_end(*m_ostr);
	}
	if(m_have_file) m_file.close();
    }

    void IO_GenEvent::use_input_units( Units::MomentumUnit mom, 
                                       Units::LengthUnit len ) {
        if( m_istr != NULL ) {
            set_input_units( *m_istr, mom, len );
	}
    }

    void IO_GenEvent::print( std::ostream& ostr ) const { 
	ostr << "IO_GenEvent: unformated ascii file IO for machine reading.\n"; 
	if(m_have_file)    ostr  << "\tFile openmode: " << m_mode ;
	ostr << " stream state: " << m_ostr->rdstate()
	     << " bad:" << (m_ostr->rdstate()&std::ios::badbit)
	     << " eof:" << (m_ostr->rdstate()&std::ios::eofbit)
	     << " fail:" << (m_ostr->rdstate()&std::ios::failbit)
	     << " good:" << (m_ostr->rdstate()&std::ios::goodbit) << std::endl;
    }

    void IO_GenEvent::precision( int size )  { 
        if( size > 16 ) { 
	    std::cerr << "IO_GenEvent::precision Error, "
		      << "precision is greater than 16. "
		      << "Not allowed. Using default precision of 16."
		      << std::endl;
            size = 16;
	}
        if(m_ostr) {
	    m_ostr->precision(size);
	}
    }
	
    bool IO_GenEvent::fill_next_event( GenEvent* evt ){
	//
	// reset error type
        m_error_type = IO_Exception::OK;
	//
	// test that evt pointer is not null
	if ( !evt ) {
            m_error_type = IO_Exception::NullEvent;
	    m_error_message = "IO_GenEvent::fill_next_event error - passed null event.";
	    std::cerr << m_error_message << std::endl;
	    return false;
	}
	// make sure the stream is good, and that it is in input mode
	if ( !(*m_istr) ) return false;
	if ( !m_istr ) {
            m_error_type = IO_Exception::WrongFileType;
	    m_error_message = "HepMC::IO_GenEvent::fill_next_event attempt to read from output file.";
	    std::cerr << m_error_message << std::endl;
	    return false;
	}
	// use streaming input
        try {
	    *m_istr >> *evt;
	}
        catch (IO_Exception& e) {
            m_error_type = IO_Exception::InvalidData;
	    m_error_message = e.what();
	    evt->clear();
 	    return false;
        }
	if( evt->is_valid() ) return true;
	return false;
    }

    void IO_GenEvent::write_event( const GenEvent* evt ) {
	/// Writes evt to output stream. It does NOT delete the event after writing.
	//
	// make sure the state is good, and that it is in output mode
	if ( !evt  ) return;
	if ( m_ostr == NULL ) {
            m_error_type = IO_Exception::WrongFileType;
	    m_error_message = "HepMC::IO_GenEvent::write_event attempt to write to input file.";
	    std::cerr << m_error_message << std::endl;
	    return;
	}
	//
	// write event listing key before first event only.
	write_HepMC_IO_block_begin(*m_ostr);
	// explicit cast is necessary
	GenEvent e = *evt;
	*m_ostr << e ;
    }

    void IO_GenEvent::write_comment( const std::string comment ) {
	// make sure the stream is good, and that it is in output mode
	if ( !(*m_ostr) ) return;
	if ( m_ostr == NULL ) {
            m_error_type = IO_Exception::WrongFileType;
	    m_error_message = "HepMC::IO_GenEvent::write_event attempt to write to input file.";
	    std::cerr << m_error_message << std::endl;
	    return;
	}
	// write end of event listing key if events have already been written
	write_HepMC_IO_block_end(*m_ostr);
	// insert the comment key before the comment
	*m_ostr << "\n" << "HepMC::IO_GenEvent-COMMENT\n";
	*m_ostr << comment << std::endl;
    }
	
} // HepMC
