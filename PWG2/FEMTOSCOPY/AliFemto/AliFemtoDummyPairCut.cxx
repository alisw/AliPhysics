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
 * Revision 1.1.1.1  2007-03-07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.3  2000/01/25 17:35:02  laue
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
 * Revision 1.2  1999/07/06 22:33:21  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:56  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "AliFemtoDummyPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoDummyPairCut)
#endif

//__________________
AliFemtoDummyPairCut::AliFemtoDummyPairCut() :
  fNPairsPassed(0),
  fNPairsFailed(0)
{
  /* no-op */
}
//__________________
//AliFemtoDummyPairCut::~AliFemtoDummyPairCut(){
//  /* no-op */
//}
//__________________
bool AliFemtoDummyPairCut::Pass(const AliFemtoPair* pair){
  // Pass all pairs
  bool temp = true;
  temp ? fNPairsPassed++ : fNPairsFailed++;
  return true;
}
//__________________
AliFemtoString AliFemtoDummyPairCut::Report(){
  // prepare a report from the execution
  string stemp = "AliFemtoDummy Pair Cut - total dummy-- always returns true\n";
  char ctemp[100];
  sprintf(ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//__________________
TList *AliFemtoDummyPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  return tListSetttings;
}
