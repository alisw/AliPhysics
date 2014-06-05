#ifndef _TauolaHEPEVTParticle_h_included_
#define _TauolaHEPEVTParticle_h_included_

/**
 * @class TauolaHEPEVTParticle
 *
 * @brief Single particle of HEPEVT event record
 *
 * This class implements the virtual methods of
 * TauolaParticle. In this way it provides an
 * interface between the generic TauolaParticle class
 * and information stored in HEPEVT event record.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include <iostream>
#include <vector>

#include "f_Decay.h"
#include "TauolaParticle.h"
#include "TauolaHEPEVTEvent.h"

namespace Tauolapp
{

class TauolaHEPEVTEvent;

class TauolaHEPEVTParticle: public TauolaParticle {

 public:
  /** Default destructor */
  ~TauolaHEPEVTParticle();

  /** Default constructor */
  TauolaHEPEVTParticle(int pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de);

  /** Remove the decay branch from the event record and reset the particle status code to stable.
      WARNING: not implemented for HEPEVT. */
  void undecay();

  /** Set the mothers of this particle via a vector of TauolaParticle*/
  void setMothers(std::vector<TauolaParticle*> mothers);

  /** Set the daughters of this particle via a vector of TauolaParticle*/
  void setDaughters(std::vector<TauolaParticle*> daughters);

  /** Returns the mothers of this particle via a vector of TauolaParticle */
  std::vector<TauolaParticle*> getMothers();

  /** Returns the daughters of this particle via a vector of TauolaParticle */
  std::vector<TauolaParticle*> getDaughters();

  /** Check that the 4 momentum in conserved in the decay of this particle */
  void checkMomentumConservation();

  /** Creates a new particle of type TauolaHEPEVTParticle, with the given
      properties. The new particle bares no relations to this 
      particle, but `this particle' provides only a way of creating an instance of
      this derived class. eg. createNewParticle() is used inside
      filhep_() so that a TauolaHEPEVTParticle can be created without
      the method having explicit knowledge of the TauolaHEPEVTParticle
      class */
  TauolaHEPEVTParticle * createNewParticle(int pdg_id, int status, double mass,
                                          double px, double py,
                                          double pz, double e);

  /** Check if particle 'p' is daughter of this particle */
  bool isDaughterOf(TauolaHEPEVTParticle *p);

  /** Check if particle 'p' is mother of this particle */
  bool isMotherOf  (TauolaHEPEVTParticle *p);

  /** Print information on this particle into standard output */
  void print();

  /** Set the PDG ID code of this particle */
  void setPdgID(int pdg_id);

  /** Set the status of this particle */
  void setStatus(int statu);

  /** Set the mass of this particle */
  void setMass(double mass);

  /** Get the PDG ID code of this particle */
  int getPdgID();

  /** Get the status of this particle */
  int getStatus();

  /** Get the mass stored (i.e. not calculated from four vector) at generation step */
  double getMass();

  /** Returns the px component of the four vector*/
  double getPx();

  /** Returns the py component of the four vector */
  double getPy();

  /** Returns the pz component of the four vector */
  double getPz();

  /** Returns the energy component of the four vector */
  double getE();

  /** Set the px component of the four vector */
  void setPx( double px );

  /** Set the px component of the four vector */
  void setPy( double py );

  /** Set the pz component of the four vector */
  void setPz( double pz );

  /** Set the energy component of the four vector */
  void setE( double e );

  /** Get the barcode (position in list) of this particle */
  int getBarcode();

  /** Set barcode (position in  list) of this particle */
  void setBarcode(int barcode);

  /** Set event of this particle */
  void setEvent(TauolaHEPEVTEvent *event);
  
  /** Get index of first mother */
  int getFirstMotherIndex();
  
  /** Get index of second mother */
  int getSecondMotherIndex();

  /** Get index of first daughter */
  int getDaughterRangeStart();

  /** Get index of last daughter */
  int getDaughterRangeEnd();  

private:

  /** Event from which this particle is taken */
  TauolaHEPEVTEvent *m_event;

  /** Position in the event record */
  int m_barcode;

  /** Indexes of mothers (-1 if do not have mothers) */
  int m_first_mother, m_second_mother;

  /** Range of indexes of daughters (-1 if do not have daughters) */
  int m_daughter_start, m_daughter_end;

  /** PDG ID */
  int m_pdgid;

  /** Status (stable, decayed) */
  int m_status;

  /** Momentum */
  double m_px, m_py, m_pz, m_e;

  /** Mass saved at generation step */
  double m_generated_mass;

  /** List of created particles - if they are not in the event, they
      will be deleted when no longer needed */
  vector<TauolaHEPEVTParticle*> cache;
};

} // namespace Tauolapp
#endif

