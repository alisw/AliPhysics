/***************************************************************************
 *
 * $Id$
 ***************************************************************************
 *          
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1  2007/05/16 10:25:06  akisiel
 * Making the directory structure of AliFemtoUser flat. All files go into one common directory
 *
 * Revision 1.4  2007/05/03 09:46:10  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.3  2007/04/27 07:25:59  akisiel
 * Make revisions needed for compilation from the main AliRoot tree
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 *
 **************************************************************************/

#include "AliFemtoQPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoQPairCut)
#endif
    
//__________________
AliFemtoQPairCut::AliFemtoQPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0)
{
  fNPairsPassed = fNPairsFailed = 0;
  fQlong[0]=-1.0; fQlong[1]=100.0;
  fQout[0]=-1.0;  fQout[1]=100.0;
  fQside[0]=-1.0; fQside[1]=100.0;
  fQinv[0]=-1.0;  fQinv[1]=100.0;
}
//__________________
AliFemtoQPairCut::~AliFemtoQPairCut()
{
//  /* no-op */
}
//__________________
bool AliFemtoQPairCut::Pass(const AliFemtoPair* pair)
{
  //bool temp = true;
  //temp ? fNPairsPassed++ : fNPairsFailed++;
  if ((fabs(pair->QLongCMS())<fQlong[0])||(fabs(pair->QLongCMS())>fQlong[1]))
  {
	fNPairsFailed++;
	return false;
  }
  if ((fabs(pair->QOutCMS())<fQout[0])||(fabs(pair->QOutCMS())>fQout[1]))
  {
	fNPairsFailed++;
	return false;
  }
  if ((fabs(pair->QSideCMS())<fQside[0])||(fabs(pair->QSideCMS())>fQside[1]))
  {
	fNPairsFailed++;
	return false;
  }
    if ((fabs(pair->KStar())<fQinv[0])||(fabs(pair->KStar())>fQinv[1]))
  {
	fNPairsFailed++;
	return false;
  }
  fNPairsPassed++;
  return true;
}
//__________________
AliFemtoString AliFemtoQPairCut::Report()
{
  string Stemp = "AliFemtoQ Pair Cut \n";
  char Ctemp[100];
  sprintf(Ctemp,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  Stemp += Ctemp;
  AliFemtoString returnThis = Stemp;
  return returnThis;
}
//__________________
TList *AliFemtoQPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoQPairCut.qout.maximum=%lf", fQout[0]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qout.minimum=%lf", fQout[1]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qside.maximum=%lf", fQside[0]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qside.minimum=%lf", fQside[1]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qlong.maximum=%lf", fQlong[0]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qlong.minimum=%lf", fQlong[1]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qinv.maximum=%lf", fQinv[0]);
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoQPairCut.qinv.minimum=%lf", fQinv[1]);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}
