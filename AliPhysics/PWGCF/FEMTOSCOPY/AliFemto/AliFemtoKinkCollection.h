/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The Collection of Kinks is the a component of the HbtEvent,
 *   which is essentially the transient microDST
 *
 ****************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.1  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 *
 ***************************************************************************/


#ifndef AliFemtoKinkCollection_hh
#define AliFemtoKinkCollection_hh
#include "AliFemtoKink.h"
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoKink*, allocator<AliFemtoKink*> >            AliFemtoKinkCollection;
typedef list<AliFemtoKink*, allocator<AliFemtoKink*> >::iterator  AliFemtoKinkIterator;
#else
typedef list<AliFemtoKink*>            AliFemtoKinkCollection;
typedef list<AliFemtoKink*>::iterator  AliFemtoKinkIterator;
#endif

#endif

