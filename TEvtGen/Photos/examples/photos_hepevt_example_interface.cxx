/**
 * Interface for FORTRAN example of use of Photos++
 *
 * @author Tomasz Przedzinski
 * @date 21 August 2013
 */
#include "Photos/Photos.h"
#include "Photos/PhotosHEPEVTEvent.h"
using namespace Photospp;

extern "C" {

  void photos_init_()
  {
    Photos::initialize();
  }

  void photos_process_()
  {
    PhotosHEPEVTEvent *event = new PhotosHEPEVTEvent();

    PhotosHEPEVTEvent::read_event_from_HEPEVT(event);
    //event->print();

    event->process();
    //event->print();

    PhotosHEPEVTEvent::write_event_to_HEPEVT(event);

    delete event;
  }

  void photos_process_particle_(int *i)
  {
    PhotosHEPEVTEvent *event = new PhotosHEPEVTEvent();

    PhotosHEPEVTEvent::read_event_from_HEPEVT(event);
    //event->print();

    PhotosParticle *p = event->getParticle(*i - 1);

    Photos::processParticle(p);
    //event->print();

    PhotosHEPEVTEvent::write_event_to_HEPEVT(event);

    delete event;
  }

  void photos_process_branch_(int *i)
  {
    PhotosHEPEVTEvent *event = new PhotosHEPEVTEvent();

    PhotosHEPEVTEvent::read_event_from_HEPEVT(event);
    //event->print();

    PhotosParticle *p = event->getParticle(*i - 1);

    Photos::processBranch(p);
    //event->print();

    PhotosHEPEVTEvent::write_event_to_HEPEVT(event);

    delete event;
  }
}
