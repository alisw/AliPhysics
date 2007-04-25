/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *  The AnalysisCollection is pointed to by the Manager, and holds pointers
 *  to all Analysis objects currently active
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.3  2000/03/17 17:23:05  laue
 * Roberts new three particle correlations implemented.
 *
 * Revision 1.2  2000/02/01 00:33:31  laue
 * namespaces changed to run on the new Solaris Compiler CC5
 * since we can use member templates in franks1Histo.h we are doing it
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef AliFemtoAnalysisCollection_hh
#define AliFemtoAnalysisCollection_hh


#include <list>
#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif
class AliFemtoBaseAnalysis;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoBaseAnalysis*, allocator<AliFemtoBaseAnalysis*> >            AliFemtoAnalysisCollection;
typedef list<AliFemtoBaseAnalysis*, allocator<AliFemtoBaseAnalysis*> >::iterator  AliFemtoAnalysisIterator;
#else
typedef list<AliFemtoBaseAnalysis*>            AliFemtoAnalysisCollection;
typedef list<AliFemtoBaseAnalysis*>::iterator  AliFemtoAnalysisIterator;
#endif

#endif
