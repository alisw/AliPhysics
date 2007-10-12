/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a do-nothing pair cut that simply says "true" to every pair           
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.2.2.1  2007/10/12 14:28:37  akisiel
 * New wave of cleanup and rule conformance
 *
 * Revision 1.2  2007/05/22 09:01:42  akisiel
 * Add the possibiloity to save cut settings in the ROOT file
 *
 * Revision 1.1  2007/05/16 10:22:11  akisiel
 * Making the directory structure of AliFemto flat. All files go into one common directory
 *
 * Revision 1.2  2007/05/03 09:41:06  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.5  2000/03/23 22:57:28  laue
 * Clone() function implemented
 *
 * Revision 1.4  2000/01/25 17:35:03  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.3  1999/10/15 01:57:05  lisa
 * Important enhancement of AliFemtoMaker - implement Franks CutMonitors
 * ----------------------------------------------------------
 * This means 3 new files in Infrastructure area (CutMonitor),
 * several specific CutMonitor classes in the Cut area
 * and a new base class in the Base area (AliFemtoCutMonitor).
 * This means also changing all Cut Base class header files from .h to .h
 * so we have access to CutMonitor methods from Cint command line.
 * This last means
 * 1) files which include these header files are slightly modified
 * 2) a side benefit: the TrackCuts and V0Cuts no longer need
 * a SetMass() implementation in each Cut class, which was stupid.
 * Also:
 * -----
 * Include Franks AliFemtoAssociationReader
 * ** None of these changes should affect any user **
 *
 * Revision 1.2  1999/07/06 22:33:21  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:56  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/


#ifndef ALIFEMTODUMMYPAIRCUT_H
#define ALIFEMTODUMMYPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoDummyPairCut : public AliFemtoPairCut{
public:
  AliFemtoDummyPairCut();
  AliFemtoDummyPairCut(const AliFemtoDummyPairCut&);
  //~AliFemtoDummyPairCut();

  virtual bool Pass(const AliFemtoPair*);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoDummyPairCut* Clone();

private:
  long fNPairsPassed;  // number of pairs analyzed by this cut that passed
  long fNPairsFailed;  // number of pairs analyzed by this cut that failed

#ifdef __ROOT__
  ClassDef(AliFemtoDummyPairCut, 1)
#endif
};

inline AliFemtoDummyPairCut::AliFemtoDummyPairCut(const AliFemtoDummyPairCut& c) : AliFemtoPairCut(c), fNPairsPassed(0), fNPairsFailed(0) { /* no-op */ }
inline AliFemtoDummyPairCut* AliFemtoDummyPairCut::Clone() { AliFemtoDummyPairCut* c = new AliFemtoDummyPairCut(*this); return c;}

#endif
