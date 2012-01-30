/***************************************************************************
 *
 * $Id$
 *
 * Author: Frank Laue, Ohio State, laue@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *  The EventWriterCollection is pointed to by the Manager, and holds pointers
 *  to all EventWriter objects currently active
 *
 ***************************************************************************
 *
 **************************************************************************/

#ifndef AliFemtoEventWriterCollection_hh
#define AliFemtoEventWriterCollection_hh

#include "AliFemtoEventWriter.h"

#include <list>
#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoEventWriter*, allocator<AliFemtoEventWriter*> >            AliFemtoEventWriterCollection;
typedef list<AliFemtoEventWriter*, allocator<AliFemtoEventWriter*> >::iterator  AliFemtoEventWriterIterator;
#else
typedef list<AliFemtoEventWriter*>            AliFemtoEventWriterCollection;
typedef list<AliFemtoEventWriter*>::iterator  AliFemtoEventWriterIterator;
#endif

#endif
