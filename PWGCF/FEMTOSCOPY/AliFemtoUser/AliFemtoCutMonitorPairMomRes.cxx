#include "AliFemtoCutMonitorPairMomRes.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include "AliFemtoPair.h"


AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes():
  fMomRes(0)
{
  // Default constructor
  fMomRes = new TH2D("PairMomRes","PairMomRes",100,0,2,100,0,2);

  
 
}

AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes(const char *aName, double qmin, double qmax, int nbins):
  fMomRes(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "PairMomRes%s", aName);
  fMomRes = new TH2D(name, "PairMomRes", nbins, qmin, qmax, nbins, qmin, qmax);

}

AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes(const AliFemtoCutMonitorPairMomRes &aCut):
  AliFemtoCutMonitor(aCut),
  fMomRes(0)
{
  // copy constructor
  fMomRes = new TH2D(*aCut.fMomRes);
}

AliFemtoCutMonitorPairMomRes::~AliFemtoCutMonitorPairMomRes()
{
  // Destructor
  delete fMomRes;
}

AliFemtoCutMonitorPairMomRes& AliFemtoCutMonitorPairMomRes::operator=(const AliFemtoCutMonitorPairMomRes& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitor::operator=(aCut);

  if (fMomRes) delete fMomRes;
  fMomRes = new TH2D(*aCut.fMomRes);

  return *this;
}

AliFemtoString AliFemtoCutMonitorPairMomRes::Report(){
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorMomRes report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoCutMonitorPairMomRes::Fill(const AliFemtoPair* aPair)
{
   double KStar_rec = aPair ->KStar();

  AliFemtoParticle *first = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *second = (AliFemtoParticle*)aPair->Track2();
  
  if(first!=NULL && second!=NULL) {

    
    AliFemtoModelHiddenInfo *tInfo1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo();
    AliFemtoModelHiddenInfo *tInfo2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

    if(tInfo1!=NULL && tInfo2!=NULL) {

      AliFemtoLorentzVector p1(::sqrt(tInfo1->GetTrueMomentum()->Mag2() + tInfo1->GetMass()*tInfo1->GetMass()), *tInfo1->GetTrueMomentum());
      AliFemtoLorentzVector p2(::sqrt(tInfo2->GetTrueMomentum()->Mag2() + tInfo2->GetMass()*tInfo2->GetMass()), *tInfo2->GetTrueMomentum());

      
      const double
	px1 = p1.x(),
	py1 = p1.y(),
	pz1 = p1.z(),
	pE1  = p1.e(),
	mass1_sqrd = std::max({0.0, p1.m2()}),

	px2 = p2.x(),
	py2 = p2.y(),
	pz2 = p2.z(),
	pE2  = p2.e(),
	mass2_sqrd = std::max({0.0, p2.m2()}),

	tPx = px1 + px2,
	tPy = py1 + py2,
	tPz = pz1 + pz2,
	tPE = pE1 + pE2;

      double tPtrans = tPx*tPx + tPy*tPy;
      double tMtrans = tPE*tPE - tPz*tPz;
      double tPinv = ::sqrt(tMtrans - tPtrans);
      tMtrans = ::sqrt(tMtrans);
      tPtrans = ::sqrt(tPtrans);

      double tQinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
	(py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);

      if( (tMtrans - tPtrans) > 0)
	{
      double tQ = (mass1_sqrd - mass2_sqrd)/tPinv;
  
      if((tQ*tQ - tQinvL)>=0.000000001)
	{
	  tQ = ::sqrt( tQ*tQ - tQinvL);

	  double KStar_true = tQ/2;

	  if(KStar_true!=0 || KStar_rec!=0)
	    fMomRes->Fill(KStar_rec,KStar_true);
	}
	}
    }
    
  }
  
   
}

void AliFemtoCutMonitorPairMomRes::Write()
{
  // Write out the relevant histograms
  fMomRes->Write();
}

TList *AliFemtoCutMonitorPairMomRes::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fMomRes);


  return tOutputList;
}
