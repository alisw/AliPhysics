/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The ParticleCollection is the main component of the picoEvent
 *   It points to the particle objects in the picoEvent.
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.2  2000/02/01 00:33:32  laue
 * namespaces changed to run on the new Solaris Compiler CC5
 * since we can use member templates in franks1Histo.h we are doing it
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/
#ifndef AliFemtoParticleCollection_hh
#define AliFemtoParticleCollection_hh
#include "AliFemtoParticle.h"
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoParticle *, allocator<AliFemtoParticle *> >            AliFemtoParticleCollection;
typedef list<AliFemtoParticle *, allocator<AliFemtoParticle *> >::iterator  AliFemtoParticleIterator;
typedef list<AliFemtoParticle *, allocator<AliFemtoParticle *> >::const_iterator  AliFemtoParticleConstIterator;
#else
typedef list<AliFemtoParticle *>            AliFemtoParticleCollection;
typedef list<AliFemtoParticle *>::iterator  AliFemtoParticleIterator;
typedef list<AliFemtoParticle *>::const_iterator  AliFemtoParticleConstIterator;
#endif

#endif
