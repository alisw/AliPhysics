/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
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
 * Revision 1.1  2000/07/16 21:44:11  laue
 * Collection and analysis for vertex dependent event mixing
 *
 *
 **************************************************************************/

#ifndef AliFemtoPicoEventCollectionVector_hh
#define AliFemtoPicoEventCollectionVector_hh
#include "AliFemtoPicoEventCollection.h"
#include <vector>
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::vector;
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef vector<AliFemtoPicoEventCollection*, allocator<AliFemtoPicoEventCollection*> >            AliFemtoPicoEventCollectionVector;  //!
typedef vector<AliFemtoPicoEventCollection*, allocator<AliFemtoPicoEventCollection*> >::iterator  AliFemtoPicoEventCollectionIterator;//!
#else
typedef vector<AliFemtoPicoEventCollection*>            AliFemtoPicoEventCollectionVector;//!
typedef vector<AliFemtoPicoEventCollection*>::iterator  AliFemtoPicoEventCollectionIterator;//!
#endif

#endif
