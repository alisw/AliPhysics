// HepMCInterface.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch
// Function definitions (not found in the header) for the I_Pythia8 class, 
// which converts a PYTHIA event record to the standard HepMC format.

#include "HepMCInterface.h"
#include "HepMC/GenEvent.h"

namespace HepMC {

//==========================================================================

// Main method for conversion from PYTHIA event to HepMC event.
// Read one event from Pythia8 and fill GenEvent, 
// and return T/F = success/failure.

bool I_Pythia8::fill_next_event( Pythia8::Event& pyev, GenEvent* evt, 
    int ievnum, Pythia8::Info* pyinfo, Pythia8::Settings* pyset) {

  // 1. Error if no event passed.
  if (!evt) {
    std::cerr << "I_Pythia8::fill_next_event error - passed null event." 
              << std::endl;
    return 0;
  }

  // Event number counter.
  if ( ievnum >= 0 ) {
    evt->set_event_number(ievnum);
    m_internal_event_number = ievnum;
  } else {
    evt->set_event_number(m_internal_event_number);
    ++m_internal_event_number;
  }

  // Conversion factors from Pythia units GeV and mm to HepMC ones.
  double momFac = HepMC::Units::conversion_factor(HepMC::Units::GEV, 
    evt->momentum_unit());
  double lenFac = HepMC::Units::conversion_factor(HepMC::Units::MM, 
    evt->length_unit());
    
  // 2. Create a particle instance for each entry and fill a map, and 
  // a vector which maps from the particle index to the GenParticle address.
  std::vector<GenParticle*> hepevt_particles( pyev.size() );
  for (int i = 1; i < pyev.size(); ++i) {

    // Fill the particle.
    hepevt_particles[i] = new GenParticle( 
      FourVector( momFac * pyev[i].px(), momFac * pyev[i].py(),
                  momFac * pyev[i].pz(), momFac * pyev[i].e()  ),
      pyev[i].id(), pyev.statusHepMC(i) );
    hepevt_particles[i]->suggest_barcode(i);
    hepevt_particles[i]->set_generated_mass( momFac * pyev[i].m() );

    // Colour flow uses index 1 and 2.
    int colType = pyev[i].colType();
    if (colType ==  1 || colType == 2)
      hepevt_particles[i]->set_flow(1, pyev[i].col());
    if (colType == -1 || colType == 2)
      hepevt_particles[i]->set_flow(2, pyev[i].acol());
  }

  // Here we assume that the first two particles in the list
  // are the incoming beam particles.
  evt->set_beam_particles( hepevt_particles[1], hepevt_particles[2] );
 
  // 3. Loop over particles AGAIN, this time creating vertices.
  // We build the production vertex for each entry in hepevt.
  // The HEPEVT pointers are bi-directional, so gives decay vertices as well.
  for (int i = 1; i < pyev.size(); ++i) {
    GenParticle *p = hepevt_particles[i];

    // 3a. Search to see if a production vertex already exists.
    std::vector<int> mothers = pyev.motherList(i);
    unsigned int imother = 0;
    int mother = -1; // note that in Pythia8 there is a particle number 0!
    if ( !mothers.empty() ) mother = mothers[imother];
    GenVertex* prod_vtx = p->production_vertex();
    while ( !prod_vtx && mother > 0 ) {
      prod_vtx = hepevt_particles[mother]->end_vertex();
      if ( prod_vtx ) prod_vtx->add_particle_out( p );
      mother = ( ++imother < mothers.size() ) ? mothers[imother] : -1;
    }

    // 3b. If no suitable production vertex exists - and the particle has
    // at least one mother or position information to store - make one.
    FourVector prod_pos( lenFac * pyev[i].xProd(), lenFac * pyev[i].yProd(),
                         lenFac * pyev[i].zProd(), lenFac * pyev[i].tProd() );
    if ( !prod_vtx && ( mothers.size() > 0 || prod_pos != FourVector() ) ) {
      prod_vtx = new GenVertex();
      prod_vtx->add_particle_out( p );
      evt->add_vertex( prod_vtx );
    }

    // 3c. If prod_vtx doesn't already have position specified, fill it.
    if ( prod_vtx && prod_vtx->position() == FourVector() )
      prod_vtx->set_position( prod_pos );

    // 3d. loop over mothers to make sure their end_vertices are consistent.
    imother = 0;
    mother = -1;
    if ( !mothers.empty() ) mother = mothers[imother];
    while ( prod_vtx && mother > 0 ) {

      // If end vertex of the mother isn't specified, do it now.
      if ( !hepevt_particles[mother]->end_vertex() ) {
        prod_vtx->add_particle_in( hepevt_particles[mother] );

      // Problem scenario: the mother already has a decay vertex which 
      // differs from the daughter's production vertex. This means there is
      // internal inconsistency in the HEPEVT event record. Print an error.
      // Note: we could provide a fix by joining the two vertices with a 
      // dummy particle if the problem arises often.
      } else if (hepevt_particles[mother]->end_vertex() != prod_vtx ) {
       if ( m_print_inconsistency ) std::cerr
          << "HepMC::I_Pythia8: inconsistent mother/daugher "
          << "information in Pythia8 event " << std::endl
          << "i = " << i << " mother = " << mother
          << "\n This warning can be turned off with the "
          << "I_Pythia8::print_inconsistency switch." << std::endl;
      }

      // End of vertex-setting loops.
      mother = ( ++imother < mothers.size() ) ? mothers[imother] : -1;
    }
  }

  // If hadronization switched on then no final coloured particles.
  bool doHadr = (pyset == 0) ? m_free_parton_warnings
    : pyset->flag("HadronLevel:all") && pyset->flag("HadronLevel:Hadronize");

  // 4. Check for particles which come from nowhere, i.e. are without 
  // mothers or daughters. These need to be attached to a vertex, or else 
  // they will never become part of the event. 
  for (int i = 1; i < pyev.size(); ++i) {
    if ( !hepevt_particles[i]->end_vertex() &&
         !hepevt_particles[i]->production_vertex() ) {
      std::cerr << "hanging particle " << i << std::endl;
      GenVertex* prod_vtx = new GenVertex();
      prod_vtx->add_particle_out( hepevt_particles[i] );
      evt->add_vertex( prod_vtx );
    }

    // Also check for free partons (= gluons and quarks; not diquarks?).
    if ( doHadr && m_free_parton_warnings ) {
      if ( hepevt_particles[i]->pdg_id() == 21 &&
        !hepevt_particles[i]->end_vertex() ) {
        std::cerr << "gluon without end vertex " << i << std::endl;
        if ( m_crash_on_problem ) exit(1);
      }
      if ( abs(hepevt_particles[i]->pdg_id()) <= 6 &&
        !hepevt_particles[i]->end_vertex()         ) {
        std::cerr << "quark without end vertex " << i << std::endl;
        if ( m_crash_on_problem ) exit(1);
      }
    }
  }

  // 5. Store PDF, weight, cross section and other event information.
  // Flavours of incoming partons.
  if (m_store_pdf && pyinfo != 0) {
    int id1pdf = pyinfo->id1pdf();
    int id2pdf = pyinfo->id2pdf();
    if ( m_convert_gluon_to_0 ) {
      if (id1pdf == 21) id1pdf = 0;
      if (id2pdf == 21) id2pdf = 0;
    }

    // Store PDF information. 
    evt->set_pdf_info( PdfInfo( id1pdf, id2pdf, pyinfo->x1pdf(),
      pyinfo->x2pdf(), pyinfo->QFac(), pyinfo->pdf1(), pyinfo->pdf2() ) );   
  }

  // Store process code, scale, alpha_em, alpha_s.
  if (m_store_proc && pyinfo != 0) {
    evt->set_signal_process_id( pyinfo->code() );
    evt->set_event_scale( pyinfo->QRen() );
    if (evt->alphaQED() <= 0) evt->set_alphaQED( pyinfo->alphaEM() );
    if (evt->alphaQCD() <= 0) evt->set_alphaQCD( pyinfo->alphaS() );
  }

  // Store cross-section information in pb and event weight. The latter is 
  // usually dimensionless, but in units of pb for Les Houches strategies +-4.
  if (m_store_xsec && pyinfo != 0) {
    HepMC::GenCrossSection xsec;
    xsec.set_cross_section( pyinfo->sigmaGen() * 1e9, pyinfo->sigmaErr() * 1e9);
    evt->set_cross_section(xsec);
    evt->weights().push_back( pyinfo->weight() );
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace HepMC
