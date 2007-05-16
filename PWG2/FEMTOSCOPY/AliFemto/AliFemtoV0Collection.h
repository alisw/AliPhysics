/***************************************************************************
 *
 * $Id$
 *
 * Author: Tom Humanic, Ohio State, humanic@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The Collection of v0s is the main component of the HbtEvent,
 *   which is essentially the transient microDST
 *
 ***************************************************************************/


#ifndef AliFemtoV0Collection_hh
#define AliFemtoV0Collection_hh
#include "AliFemtoV0.h"
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoV0*, allocator<AliFemtoV0*> >            AliFemtoV0Collection;
typedef list<AliFemtoV0*, allocator<AliFemtoV0*> >::iterator  AliFemtoV0Iterator;
#else
typedef list<AliFemtoV0*>            AliFemtoV0Collection;
typedef list<AliFemtoV0*>::iterator  AliFemtoV0Iterator;
#endif

#endif

