#ifndef _TauolaHEPEVTEvent_h_included_
#define _TauolaHEPEVTEvent_h_included_

/**
 * @class TauolaHEPEVTParticle
 *
 * @brief Single particle of HEPEVT event record
 *
 * This class implements the virtual methods of
 * TauolaEvent. In this way it provides an
 * interface between the generic TauolaEvent class
 * and information stored in HEPEVT event record.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include <iostream>
#include "TauolaEvent.h"
#include "TauolaParticle.h"
#include "TauolaHEPEVTParticle.h"

namespace Tauolapp
{

// Uncomment this line to use interface to common block HEPEVT
// But first be sure about suitable for you value of NMXHEP
// and whether phep, vhep should be declared float or double
//#define USE_HEPEVT_INTERFACE

#ifdef USE_HEPEVT_INTERFACE

// Change this value to match HEPEVT size
const int NMXHEP = 10000;

extern "C" struct {
  int    nevhep;            // serial number
  int    nhep;              // number of particles
  int    isthep[NMXHEP];    // status code
  int    idhep [NMXHEP];    // particle PDG ID
  int    jmohep[NMXHEP][2]; // parent particles
  int    jdahep[NMXHEP][2]; // childreen particles
  double phep  [NMXHEP][5]; // four-momentum, mass [GeV]
  double vhep  [NMXHEP][4]; // vertex [mm]
} hepevt_;

#endif

class TauolaHEPEVTParticle;

class TauolaHEPEVTEvent : public TauolaEvent {

 public:

  /** Default destructor */
  ~TauolaHEPEVTEvent();

  /** Default constructor */
  TauolaHEPEVTEvent();

  /** Add particle at the end of event record */
  void addParticle(TauolaHEPEVTParticle *p);

  /** Get particle at index 'i' */
  TauolaHEPEVTParticle *getParticle(int i);

  /** Get higher-most index of the particles in event (nhep) */
  int getParticleCount();

  /** Implementation of TauolaEvent virtual method.
      This returns a list of particles in the event with
      pdg id = "pdgID". */
  std::vector<TauolaParticle*> findParticles(int pdgID);

  /** Implementation of TauolaEven virtual method.
      This returns a list of particles in the event with
      pdg id = "pdgID" and stable status code. */
  std::vector<TauolaParticle*> findStableParticles(int pdgID);

  /** Print out list of particles in the event */
  void print();
  
  /** Remove all particles from the event */
  void clear();

#ifdef USE_HEPEVT_INTERFACE
  /** Fill TauolaHEPEVTEvent from HEPEVT common block */
  static void read_event_from_HEPEVT(TauolaHEPEVTEvent *evt);
  
  /** Write to HEPEVT common block content of TauolaHEPEVTEvent */
  static void write_event_to_HEPEVT(TauolaHEPEVTEvent *evt);
#endif

 private:

  /** List of all particles */
  std::vector<TauolaHEPEVTParticle*> particle_list;
};

} // namespace Tauolapp
#endif

