/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutAntiGammaAlpha - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutAntiGammaAlpha.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
 *
 * Author: Adam Kisiel, Ohio State, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/

#include "AliFemtoPairCutAntiGammaAlpha.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutAntiGammaAlpha)
#endif

//__________________
AliFemtoPairCutAntiGammaAlpha::AliFemtoPairCutAntiGammaAlpha():
  AliFemtoShareQualityPairCut(),
  fMaxEEMinv(0.0),
  fMaxDAlpha(0.0),
  fDTPCMin(0),
  fDataType(kESD)
{
}
//__________________
AliFemtoPairCutAntiGammaAlpha::AliFemtoPairCutAntiGammaAlpha(const AliFemtoPairCutAntiGammaAlpha& c) : 
  AliFemtoShareQualityPairCut(c),
  fMaxEEMinv(0.0),
  fMaxDAlpha(0.0),
  fDTPCMin(0),
  fDataType(kESD)
{ 
  fMaxEEMinv = c.fMaxEEMinv;
  fMaxDAlpha = c.fMaxDAlpha;
  fDTPCMin = c.fDTPCMin;
  fDataType = c.fDataType;
}

AliFemtoPairCutAntiGammaAlpha& AliFemtoPairCutAntiGammaAlpha::operator=(const AliFemtoPairCutAntiGammaAlpha& c)
{
  if (this != &c) {
    fMaxEEMinv = c.fMaxEEMinv;
    fMaxDAlpha = c.fMaxDAlpha;
    fDTPCMin = c.fDTPCMin;
    fDataType = c.fDataType;
  }

  return *this;

}
//__________________
AliFemtoPairCutAntiGammaAlpha::~AliFemtoPairCutAntiGammaAlpha(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutAntiGammaAlpha::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC entrance separation and
  // quality and sharity
  bool temp = true;

  if(fDataType==kKine)
    return true;

  double me = 0.000511;

  if ((pair->Track1()->Track()->Charge() * pair->Track2()->Track()->Charge()) < 0.0) {
    double dalpha = TMath::Abs((pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x()
		   + pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y()
		   + pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()
				)/pair->Track1()->Track()->P().Mag()/pair->Track2()->Track()->P().Mag());

    
    double e1 = TMath::Sqrt(me*me + pair->Track1()->Track()->P().Mag2());
    double e2 = TMath::Sqrt(me*me + pair->Track2()->Track()->P().Mag2());
    
    double minv = (2*me*me + 2*(e1*e2 - 
			       pair->Track1()->Track()->P().x()*pair->Track2()->Track()->P().x() -
			       pair->Track1()->Track()->P().y()*pair->Track2()->Track()->P().y() -
					   pair->Track1()->Track()->P().z()*pair->Track2()->Track()->P().z()));

    double sminv = TMath::Sqrt(minv);
  
    if ((sminv < fMaxEEMinv) && (dalpha > fMaxDAlpha)) {
      temp = false;
    }
  }

  bool tempTPCEntrance = true;
  
  if(fDataType==kESD || fDataType==kAOD)
    {
      double distx = pair->Track1()->Track()->NominalTpcEntrancePoint().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
      double disty = pair->Track1()->Track()->NominalTpcEntrancePoint().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
      double distz = pair->Track1()->Track()->NominalTpcEntrancePoint().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
      double dist = sqrt(distx*distx + disty*disty + distz*distz);

      tempTPCEntrance = dist > fDTPCMin;
    }
 

  if (temp && tempTPCEntrance) {
    
    temp = AliFemtoShareQualityPairCut::Pass(pair);
    if (temp) {
      fNPairsPassed++;
    }
    else fNPairsFailed++;
    return temp;
  }
  else
    {
    fNPairsFailed++;
    return false;
    }



}
//__________________
AliFemtoString AliFemtoPairCutAntiGammaAlpha::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoPairCutAntiGammaAlpha Pair Cut - remove pairs possibly coming from Gamma conversions\n";  
  char ctemp[100];
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutAntiGammaAlpha::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutAntiGammaAlpha.maxeeminv=%f", fMaxEEMinv);
  snprintf(buf, 200, "AliFemtoPairCutAntiGammaAlpha.maxdalpha=%f", fMaxDAlpha);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutAntiGammaAlpha::SetMaxEEMinv(Double_t maxeeminv)
{
  fMaxEEMinv = maxeeminv;
}

 
void AliFemtoPairCutAntiGammaAlpha::SetMaxAlphaDiff(Double_t maxdalpha)
{
  fMaxDAlpha = maxdalpha;
}

void AliFemtoPairCutAntiGammaAlpha::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoPairCutAntiGammaAlpha::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}
