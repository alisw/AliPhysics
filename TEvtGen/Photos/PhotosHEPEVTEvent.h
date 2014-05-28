#ifndef _PhotosHEPEVTEvent_h_included_
#define _PhotosHEPEVTEvent_h_included_

/**
 * @class PhotosHEPEVTParticle
 *
 * @brief Single particle of HEPEVT event record
 *
 * This class implements the virtual methods of
 * PhotosEvent. In this way it provides an
 * interface between the generic PhotosEvent class
 * and information stored in HEPEVT event record.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include <iostream>
#include "f_Init.h"
#include "PhotosEvent.h"
#include "PhotosParticle.h"
#include "PhotosHEPEVTParticle.h"

namespace Photospp
{

class PhotosHEPEVTParticle;

class PhotosHEPEVTEvent : public PhotosEvent {

 public:

  /** Default destructor */
  ~PhotosHEPEVTEvent();

  /** Default constructor */
  PhotosHEPEVTEvent();

  /** Add particle at the end of event record */
  void addParticle(PhotosHEPEVTParticle *p);

  /** Get particle at index 'i' */
  PhotosHEPEVTParticle *getParticle(int i);

  /** Set particle at index 'i' */
  void setParticle(int i, PhotosHEPEVTParticle *p);

  /** Get higher-most index of the particles in event (nhep) */
  int getParticleCount();

  /** Get an unfiltered list of particles from the event record */
  std::vector<PhotosParticle*> getParticleList();

  /** Print out list of particles in the event */
  void print();
  
  /** Remove all particles from the event */
  void clear();

  /** Fill PhotosHEPEVTEvent from HEPEVT common block */
  static void read_event_from_HEPEVT(PhotosHEPEVTEvent *evt);
  
  /** Write to HEPEVT common block content of PhotosHEPEVTEvent */
  static void write_event_to_HEPEVT(PhotosHEPEVTEvent *evt);

 private:

  /** List of all particles */
  std::vector<PhotosHEPEVTParticle*> particle_list;
};

} // namespace Photospp
#endif

