#ifndef _PhotosHEPEVTParticle_h_included_
#define _PhotosHEPEVTParticle_h_included_

/**
 * @class PhotosHEPEVTParticle
 *
 * @brief Single particle of HEPEVT event record
 *
 * This class implements the virtual methods of
 * PhotosParticle. In this way it provides an
 * interface between the generic PhotosParticle class
 * and information stored in HEPEVT event record.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

#include "Photos.h"
#include "PhotosParticle.h"
#include "PhotosHEPEVTEvent.h"

namespace Photospp
{

class PhotosHEPEVTEvent;

class PhotosHEPEVTParticle: public PhotosParticle {

 public:
  /** Default destructor */
  ~PhotosHEPEVTParticle();

  /** Default constructor */
  PhotosHEPEVTParticle(int pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de);

  /** Add a new daughter to this particle */
  void addDaughter(PhotosParticle* daughter);

  /** Set the mothers of this particle via a vector of PhotosParticle*/
  void setMothers(std::vector<PhotosParticle*> mothers);

  /** Set the daughters of this particle via a vector of PhotosParticle*/
  void setDaughters(std::vector<PhotosParticle*> daughters);

  /** Returns the mothers of this particle via a vector of PhotosParticle */
  std::vector<PhotosParticle*> getMothers();

  /** Returns the daughters of this particle via a vector of PhotosParticle */
  std::vector<PhotosParticle*> getDaughters();

  /** Returns all particles in the decay tree of this particle
      via a vector of PhotosParticle */
  std::vector<PhotosParticle*> getAllDecayProducts();

  /** Check that the 4 momentum in conserved in the decay of this particle */
  bool checkMomentumConservation();

  /** Creates a new particle of type PhotosHEPEVTParticle, with the given
      properties. The new particle bares no relations to this 
      particle, but `this particle' provides only a way of creating an instance of
      this derived class. eg. createNewParticle() is used inside
      filhep_() so that a PhotosHEPEVTParticle can be created without
      the method having explicit knowledge of the PhotosHEPEVTParticle
      class */
  PhotosHEPEVTParticle * createNewParticle(int pdg_id, int status, double mass,
                                          double px, double py,
                                          double pz, double e);

  /** Creating history entries not implemented in HEPEVT */
  void createHistoryEntry();

  /** Create a self-decay vertex for this particle
      with 'out' being the outgoing particle in new vertex */
  void createSelfDecayVertex(PhotosParticle *out);

  /** Check if particle 'p' is daughter of this particle */
  bool isDaughterOf(PhotosHEPEVTParticle *p);

  /** Check if particle 'p' is mother of this particle */
  bool isMotherOf  (PhotosHEPEVTParticle *p);

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
  void setEvent(PhotosHEPEVTEvent *event);
  
  /** Get index of first mother */
  int getFirstMotherIndex();
  
  /** Get index of second mother */
  int getSecondMotherIndex();

  /** Get index of first daughter */
  int getDaughterRangeStart();

  /** Get index of last daughter */
  int getDaughterRangeEnd();  

private:

  /** Set index of first daughter */
  void setDaughterRangeStart(int i) { m_daughter_start=i; }
  
  /** Set index of last daughter */
  void setDaughterRangeEnd(int i)   { m_daughter_end  =i; }

  /** Event from which this particle is taken */
  PhotosHEPEVTEvent *m_event;

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
  vector<PhotosHEPEVTParticle*> cache;
};

} // namespace Photospp
#endif

