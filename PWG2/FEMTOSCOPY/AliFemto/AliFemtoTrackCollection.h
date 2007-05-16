/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The Collection of Tracks is the main component of the HbtEvent,
 *   which is essentially the transient microDST
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

#ifndef AliFemtoTrackCollection_hh
#define AliFemtoTrackCollection_hh
#include "AliFemtoTrack.h"
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoTrack*, allocator<AliFemtoTrack*> >            AliFemtoTrackCollection;
typedef list<AliFemtoTrack*, allocator<AliFemtoTrack*> >::iterator  AliFemtoTrackIterator;
#else
typedef list<AliFemtoTrack*>            AliFemtoTrackCollection;
typedef list<AliFemtoTrack*>::iterator  AliFemtoTrackIterator;
#endif

#endif
