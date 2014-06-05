//--------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////
// Mikhail.Kirsanov@Cern.CH, 2006
// event input/output in ascii format for eye and machine reading
//
// for arguments mostly similar to IO_Ascii. Special value of
// argument filename in constructor: if it is "cout" the output is to std::cout
//////////////////////////////////////////////////////////////////////////////

#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Version.h"

namespace HepMC {

  IO_AsciiParticles::IO_AsciiParticles( const char* filename, std::ios::openmode mode ) 
  : m_precision(2),
    m_mode(mode), m_finished_first_event_io(0)
  {
    if(std::string(filename) == std::string("cout")) {
      m_outstream = &(std::cout);
      m_file = 0;
    } else {
      m_file = new std::fstream(filename, mode);
      m_outstream = m_file;
      if ( (m_mode&std::ios::out && m_mode&std::ios::in) ||
           (m_mode&std::ios::app && m_mode&std::ios::in) ) {
	    std::cerr << "IO_AsciiParticles::IO_AsciiParticles Error, open of file requested "
                  << "of input AND output type. Not allowed. Closing file."
                  << std::endl;
        m_file->close();
        delete m_file;
	return;
      }
    }
    // precision 16 (# digits following decimal point) is the minimum that
    // will capture the full information stored in a double
    // with precision <= 2 the width of output will be < 80 characters
    m_outstream->precision(m_precision);
    // we use decimal to store integers, because it is smaller than hex!
    m_outstream->setf(std::ios::dec,std::ios::basefield);
    m_outstream->setf(std::ios::scientific,std::ios::floatfield);
  }

  IO_AsciiParticles::~IO_AsciiParticles() {
    if(m_file) {
       m_file->close();
       delete m_file;
    }
  }

  void IO_AsciiParticles::print( std::ostream& ostr ) const { 
    ostr << "IO_AsciiParticles: formated ascii file IO for eye and machine reading.\n" 
         << "\tFile openmode: " << m_mode 
         << " file state: " << m_outstream->rdstate()
         << " bad:" << (m_outstream->rdstate()&std::ios::badbit)
         << " eof:" << (m_outstream->rdstate()&std::ios::eofbit)
         << " fail:" << (m_outstream->rdstate()&std::ios::failbit)
         << " good:" << (m_outstream->rdstate()&std::ios::goodbit) << std::endl;
  }

  void IO_AsciiParticles::write_event( const GenEvent* evt ) {
  // Writes evt to m_outstream. It does NOT delete the event after writing.
    //
	// check the state of m_outstream is good, and that it is in output mode
	if ( !evt || !m_outstream ) return;
	if ( !(m_mode&std::ios::out) ) {
	    std::cerr << "HepMC::IO_AsciiParticles::write_event "
		      << " attempt to write to input file." << std::endl;
	    return;
	}
	//
	// write event listing key before first event only.
	if ( !m_finished_first_event_io ) {
	    m_finished_first_event_io = 1;
        *m_outstream << "0 Run  HepMC::IO_AsciiParticles eye-readable events output"
                     << std::endl;
        *m_outstream << "#      HepMC::Version " << versionName() << std::endl;
        *m_outstream <<
    "  #  stat pdg  moth1   px        py         pz     energy    mass      eta"
                     << std::endl;
	}
	//
	// output the event data
	std::vector<long int> random_states = evt->random_states();
	*m_outstream << evt->event_number() << " Event" << std::endl;
#if 0
	*m_outstream << " " << evt->event_scale();
	output( evt->alphaQCD() );
	output( evt->alphaQED() );
	output( evt->signal_process_id() );
	output(   ( evt->signal_process_vertex() ?
		    evt->signal_process_vertex()->barcode() : 0 )   );
	output( evt->vertices_size() ); // total number of vertices.
	output( (int)random_states.size() );
	for ( std::vector<long int>::iterator rs = random_states.begin(); 
	      rs != random_states.end(); ++rs ) {
	    output( *rs );
	}
	output( (int)evt->weights().size() );
	for ( WeightContainer::const_iterator w = evt->weights().begin(); 
	      w != evt->weights().end(); ++w ) {
	    output( *w );
	}
	output('\n');
#endif
	//
    int nparticles=0, imoth=0, ip=0, istati;
    double xmassi, etai;
    *m_outstream << evt->particles_size() << " particles" << std::endl;
    GenVertex* orig;
    for(HepMC::GenEvent::particle_const_iterator part = evt->particles_begin();
        part != evt->particles_end(); ++part ) {
      //if( (*part)->status() != 1 ) continue;
      nparticles++;
      ip++;
      istati = (*part)->status();
      if( (*part)->end_vertex() && istati == 1) {
        std::cout << "final particle with end vertex!" << std::endl;
        istati = -100;
      }
      imoth=0;
      orig = (*part)->production_vertex();
      if(orig) {
        imoth = 0;
        bool ifound=false;
        for(HepMC::GenEvent::particle_const_iterator part1 =
                                                     evt->particles_begin();
                                                     part1 != part; part1++ ) {
          imoth++;
          if( (*part1)->end_vertex() == orig ) { ifound = true; break; }
        }
        if(!ifound) imoth = 0;
      }

      m_outstream->width(4);
      *m_outstream << ip << " ";

      m_outstream->width(3);
      *m_outstream << istati << " ";

      m_outstream->width(5);
      *m_outstream << (*part)->pdg_id() << " ";

      m_outstream->width(3);
      *m_outstream << imoth << "  ";

      if((*part)->momentum().px() >= 0.) *m_outstream << " ";
      *m_outstream << (*part)->momentum().px() << " ";
      if((*part)->momentum().py() >= 0.) *m_outstream << " ";
      *m_outstream << (*part)->momentum().py() << " ";
      if((*part)->momentum().pz() >= 0.) *m_outstream << " ";
      *m_outstream << (*part)->momentum().pz() << " "
             << (*part)->momentum().e() << " ";

      xmassi = (*part)->generatedMass();
      if(fabs(xmassi) < 0.0001) xmassi =0.;
      m_outstream->setf(std::ios::fixed);
      m_outstream->precision(3);
      m_outstream->width(8);
      *m_outstream << xmassi << " ";
      m_outstream->setf(std::ios::scientific,std::ios::floatfield);
      m_outstream->precision(m_precision);

      m_outstream->setf(std::ios::fixed);
      m_outstream->precision(3);
      m_outstream->width(6);
      etai = (*part)->momentum().eta();
      if(etai > 999.)etai = 999.;
      if(etai < -999.)etai = -999.;
      *m_outstream << etai << std::endl;
      m_outstream->setf(std::ios::scientific,std::ios::floatfield);
      m_outstream->precision(m_precision);

    }
  }

  bool IO_AsciiParticles::fill_next_event( GenEvent* evt ){
	//
	//
	// test that evt pointer is not null
	if ( !evt ) {
	    std::cerr 
		<< "IO_AsciiParticles::fill_next_event error - passed null event." 
		<< std::endl;
	    return false;
	}
	// check the state of m_outstream is good, and that it is in input mode
	if ( !m_file )
      std::cerr << "HepMC::IO_AsciiParticles::fill_next_event "
                << " no file for input" << std::endl;
	if ( !(m_mode&std::ios::in) ) {
	    std::cerr << "HepMC::IO_AsciiParticles::fill_next_event "
		      << " attempt to read from output file" << std::endl;
	    return false;
	}
    std::cerr << "IO_AsciiParticles input is not yet implemented" << std::endl;
    return false;
  }

  void IO_AsciiParticles::write_comment( const std::string comment ) {
	// check the state of *m_outstream is good, and that it is in output mode
	if ( !m_outstream ) return;
	if ( !(m_mode&std::ios::out) ) {
	    std::cerr << "HepMC::IO_AsciiParticles::write_particle_data_table "
		      << " attempt to write to input file." << std::endl;
	    return;
	}
	// write end of event listing key if events have already been written
	write_end_listing();
	// insert the comment key before the comment
	*m_outstream << "\n" << "HepMC::IO_AsciiParticles-COMMENT\n";
	*m_outstream << comment << std::endl;
  }

  bool IO_AsciiParticles::write_end_listing() {
	return false;
  }

} // HepMC

