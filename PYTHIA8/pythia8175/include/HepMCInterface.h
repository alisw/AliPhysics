// HepMCInterface.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch
// Header file for the I_Pythia8 class,
// which converts a PYTHIA event record to the standard HepMC format.

#ifndef Pythia8_HepMCInterface_H
#define Pythia8_HepMCInterface_H

#include <vector>
#include "HepMC/IO_BaseClass.h"
#include "Pythia.h"

namespace HepMC {

// Forward references to some classes.
class GenEvent;
class GenVertex;
class GenParticle;

//==========================================================================

// The I_Pythia8 class.

class I_Pythia8 : public IO_BaseClass {

public:

  // Constructor and destructor.
  I_Pythia8() : m_internal_event_number(0), 
    m_print_inconsistency(true), m_free_parton_warnings(true), 
    m_crash_on_problem(false),   m_convert_gluon_to_0(false),
    m_store_pdf(true), m_store_proc(true), m_store_xsec(true) {;}
  virtual ~I_Pythia8() {;}

  // The recommended method to convert Pythia events into HepMC ones.
  bool fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt, 
    int ievnum = -1 ) {return fill_next_event( pythia.event, evt, 
    ievnum, &pythia.info, &pythia.settings);}

  // Alternative method to convert Pythia events into HepMC ones.
  bool fill_next_event( Pythia8::Event& pyev, GenEvent* evt, 
    int ievnum = -1, Pythia8::Info* pyinfo = 0, 
    Pythia8::Settings* pyset = 0);

  // Read out values for some switches.
  bool print_inconsistency()  const {return m_print_inconsistency;}
  bool free_parton_warnings() const {return m_free_parton_warnings;}
  bool crash_on_problem()     const {return m_crash_on_problem;}
  bool convert_gluon_to_0()   const {return m_convert_gluon_to_0;}
  bool store_pdf()            const {return m_store_pdf;}
  bool store_proc()           const {return m_store_proc;}
  bool store_xsec()           const {return m_store_xsec;}  

  // Set values for some switches.
  void set_print_inconsistency(bool b = true)  {m_print_inconsistency = b;}
  void set_free_parton_warnings(bool b = true) {m_free_parton_warnings = b;}
  void set_crash_on_problem(bool b = false)    {m_crash_on_problem = b;}
  void set_convert_gluon_to_0(bool b = false)  {m_convert_gluon_to_0 = b;}
  void set_store_pdf(bool b = true)            {m_store_pdf = b;}
  void set_store_proc(bool b = true)           {m_store_proc = b;}
  void set_store_xsec(bool b = true)           {m_store_xsec = b;}  

private: 

  // Following methods are not implemented for this class.
  virtual bool fill_next_event( GenEvent*  ) { return 0; }
  virtual void write_event( const GenEvent* ) {;}

  // Use of copy constructor is not allowed.
  I_Pythia8( const I_Pythia8& ) : IO_BaseClass() {}

  // Data members.
  int  m_internal_event_number;
  bool m_print_inconsistency, m_free_parton_warnings,
       m_crash_on_problem, m_convert_gluon_to_0, 
       m_store_pdf, m_store_proc, m_store_xsec;

};

//==========================================================================

} // end namespace HepMC

#endif  // end Pythia8_HepMCInterface_H
