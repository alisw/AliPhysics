///
/// \file AliFemtoCorrFctnCollection.h
/// \author Mike Lisa, Ohio State <lisa@mps.ohio-state.edu>
///
/// \brief Container of pointers to all correlation functions which are associated
///        with a particular Analysis object.
///

/***************************************************************************
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

#ifndef AliFemtoCorrFctnCollection_hh
#define AliFemtoCorrFctnCollection_hh


#include <list>
#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif
class AliFemtoCorrFctn;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoCorrFctn*, allocator<AliFemtoCorrFctn*> >            AliFemtoCorrFctnCollection;
typedef list<AliFemtoCorrFctn*, allocator<AliFemtoCorrFctn*> >::iterator  AliFemtoCorrFctnIterator;
#else
typedef list<AliFemtoCorrFctn*>            AliFemtoCorrFctnCollection;
typedef list<AliFemtoCorrFctn*>::iterator  AliFemtoCorrFctnIterator;
#endif

#endif
