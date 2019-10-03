/***************************************************************************
 *
 * $Id$
 *
 * Author: Frank Laue, BNL
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The Collection of v0s is the main component of the HbtEvent,
 *   which is essentially the transient microDST
 *
 ***************************************************************************/


#ifndef AliFemtoXiCollection_hh
#define AliFemtoXiCollection_hh
#include "AliFemtoXi.h"
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoXi*, allocator<AliFemtoXi*> >            AliFemtoXiCollection;
typedef list<AliFemtoXi*, allocator<AliFemtoXi*> >::iterator  AliFemtoXiIterator;
#else
typedef list<AliFemtoXi*>            AliFemtoXiCollection;
typedef list<AliFemtoXi*>::iterator  AliFemtoXiIterator;
#endif

#endif

