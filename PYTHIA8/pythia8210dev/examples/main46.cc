// main46.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
//
// Author: S.Chekanov (ANL).
// An example of how to write ProMC files.
//
// The ProMC library is described at http://atlaswww.hep.anl.gov/asc/promc/
// A makefile can be found in the ProMC  package, in examples/pythia.

#include <map>

// ProMC file. Google does not like these warnings.
#pragma GCC diagnostic ignored "-pedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#include "ProMCBook.h"

// Pythia header and namespace.
#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

string getEnvVar( std::string const & key ) {
  char * val = getenv( key.c_str() );
  return (val == NULL) ? std::string("") : std::string(val);
}

//--------------------------------------------------------------------------

void readPDG( ProMCHeader * header  ) {

  string temp_string;
  istringstream curstring;
  string PdgTableFilename = getEnvVar("PROMC");
  if (PdgTableFilename.size() < 2) PdgTableFilename = string(PROMC);
  PdgTableFilename += "/data/particle.tbl";

  ifstream file_to_read(PdgTableFilename.c_str());
  if (!file_to_read.good()) {
    cout << "**        ERROR: PDG Table (" << PdgTableFilename
         <<  ") not found! exit.                        **" << endl;
    exit(1);
    return;
  }

  // First three lines of the file are useless.
  getline(file_to_read,temp_string);
  getline(file_to_read,temp_string);
  getline(file_to_read,temp_string);

  while (getline(file_to_read,temp_string)) {
    // Needed when using several times istringstream::str(string).
    curstring.clear();
    curstring.str(temp_string);
    long int ID; std::string name; int charge; float mass; float width;
    float lifetime;
    // ID name   chg       mass    total width   lifetime
    //  1 d      -1      0.33000     0.00000   0.00000E+00
    //  in the table, the charge is in units of e+/3
    //  the total width is in GeV
    //  the lifetime is ctau in mm
    curstring >> ID >> name >> charge >> mass >> width >> lifetime;
    ProMCHeader_ParticleData* pp= header->add_particledata();
    pp->set_id(ID);
    pp->set_mass(mass);
    pp->set_name(name);
    pp->set_width(width);
    pp->set_lifetime(lifetime);
    cout << ID << " " << name << " " << mass << endl;
  }

}

//--------------------------------------------------------------------------

int main() {

  int Ntot = 1000; // Total number of events.

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();

  // ****************  book ProMC file **********************
  ProMCBook* epbook = new ProMCBook("Pythia8.promc","w");
  epbook->setDescription(Ntot,"PYTHIA8");

  // Info on incoming beams and CM energy.
  ProMCHeader header;
  header.set_id1( pythia.info.idA() );
  header.set_id2( pythia.info.idB() );
  header.set_ecm( pythia.info.eCM() );
  header.set_s( pythia.info.s() );

 // Use the range 0.01 MeV to 20 TeV using varints (integers).
 // With particle in GeV, we multiply it by kEV, to get 0.01 MeV = 1 unit.
 const double kEV = 1000*100;
 // With particle in mm, we multiply it by kL to get 0.01 mm = 1 unit.
  const double kL = 100;

  // Set units.
  header.set_momentumunit( (int)kEV );
  header.set_lengthunit( (int)kL );

   // Store a map with PDG information (stored in the header).
  readPDG( &header );
  epbook->setHeader(header); // write header

  // Begin event loop. Generate event. Skip if error.
  for (int n = 0; n < Ntot; n++) {
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
        cout << " Aborted since reached end of Les Houches Event File\n";
        break;
      }
      continue;
    }

    //************  ProMC file ***************//
    ProMCEvent promc;

    // Fill event information.
    ProMCEvent_Event *eve = promc.mutable_event();
    eve->set_number(n);
    eve->set_process_id( pythia.info.code() );     // process ID
    eve->set_scale( pythia.info.pTHat( ));         // relevant for 2 -> 2 only
    eve->set_alpha_qed( pythia.info.alphaEM() );
    eve->set_alpha_qcd( pythia.info.alphaS() );
    eve->set_scale_pdf( pythia.info.QFac() );
    eve->set_x1( pythia.info.x1pdf() );
    eve->set_x2( pythia.info.x2pdf() );
    eve->set_id1( pythia.info.id1pdf() );
    eve->set_id2( pythia.info.id2pdf() );
    eve->set_pdf1( pythia.info.pdf1() );
    eve->set_pdf2( pythia.info.pdf2() );
    eve->set_weight( pythia.info.weight() );

    // Fill truth particle information, looping over all particles in event.
    ProMCEvent_Particles *pa= promc.mutable_particles();
    for (int i = 0; i < pythia.event.size(); i++) {

      // Fill information particle by particle.
      pa->add_id( i  );
      pa->add_pdg_id( pythia.event[i].id() );
      // Particle status in HEPMC style.
      // pa->add_status(  pythia.event.statusHepMC(i) );
      pa->add_status(  pythia.event[i].status() );
      pa->add_mother1( pythia.event[i].mother1() );
      pa->add_mother2( pythia.event[i].mother2() );
      pa->add_daughter1( pythia.event[i].daughter1() );
      pa->add_daughter2( pythia.event[i].daughter2() );
      // Only store three-momentum and mass, so need to calculate energy.
      pa->add_px( (int)(pythia.event[i].px()*kEV) );
      pa->add_py( (int)(pythia.event[i].py()*kEV) );
      pa->add_pz( (int)(pythia.event[i].pz()*kEV) );
      pa->add_mass( (int)(pythia.event[i].m()*kEV) );
      // Store production vertex; will often be the origin.
      pa->add_x( int(pythia.event[i].xProd()*kL) );
      pa->add_y( int(pythia.event[i].yProd()*kL) );
      pa->add_z( int(pythia.event[i].zProd()*kL) );
      pa->add_t( int(pythia.event[i].tProd()*kL) );

    } // end loop over particles in the event
    epbook->write(promc); // write event

  } // end loop over events

  // Print statistics.
  pythia.stat();
  double sigmapb = pythia.info.sigmaGen() * 1.0E9;
  cout << "== Cross section for this run = " <<  sigmapb << " pb" << endl;
  cout << "== Events for this run = " <<  Ntot << endl;
  double lumi = (Ntot/sigmapb)/1000;
  cout << "== Luminosity for this run = " <<  lumi  << " fb-1" << endl;
  cout << "\n\n";

  // Save post-generation statistics for ProMC.
  ProMCStat stat;
  stat.set_cross_section_accumulated( sigmapb ); // in pb
  stat.set_cross_section_error_accumulated( pythia.info.sigmaErr() * 1e9 );
  stat.set_luminosity_accumulated( Ntot/sigmapb );
  stat.set_ntried( pythia.info.nTried() );
  stat.set_nselected( pythia.info.nSelected() );
  stat.set_naccepted( pythia.info.nAccepted() );
  epbook->setStatistics( stat );

  // Close the ProMC file.
  epbook->close(); // close

  return 0;
}
