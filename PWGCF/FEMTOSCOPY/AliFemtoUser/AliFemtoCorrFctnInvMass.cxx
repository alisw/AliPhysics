////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnInvMass.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnInvMass)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnInvMass::AliFemtoCorrFctnInvMass(const char* title, const int& aBins=3000, const double& aMin = 1000, const double& aMax = 5000, const double& aMass1 = 0, const double& aMass2 = 0):
  AliFemtoCorrFctn(),
  fNumInvMass(0),
  fDenInvMass(0),
  fBins(0),
  fMin(0),
  fMax(0),
  fMass1(0)
{

  fBins = aBins;
  fMin = aMin;
  fMax = aMax;

  fMass1 = aMass1;
  fMass2 = aMass2;
  

  // set up numerator
  char tTitNumD[101] = "NumInvMass";
  strncat(tTitNumD,title, 100);
  fNumInvMass = new TH1D(tTitNumD,title,aBins,aMin,aMax);
  // set up denominator
  char tTitDenD[101] = "DenInvMass";
  strncat(tTitDenD,title, 100);
  fDenInvMass = new TH1D(tTitDenD,title,aBins,aMin,aMax);


  // to enable error bar calculation...
  fNumInvMass->Sumw2();
  fDenInvMass->Sumw2();



}

//____________________________
AliFemtoCorrFctnInvMass::AliFemtoCorrFctnInvMass(const AliFemtoCorrFctnInvMass& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumInvMass(0),
  fDenInvMass(0),
  fBins(0),
  fMin(0),
  fMax(0),
  fMass1(0),
  fMass2(0)
{
  // copy constructor

  fBins = aCorrFctn.fBins;
  fMin = aCorrFctn.fMax;
  fMax = aCorrFctn.fMax;

  fMass1 = aCorrFctn.fMass1;
  fMass2 = aCorrFctn.fMass2;

  // copy constructor
  if (aCorrFctn.fNumInvMass)
    fNumInvMass = new TH1D(*aCorrFctn.fNumInvMass);
  else
    fNumInvMass= 0;
  if (aCorrFctn.fDenInvMass)
    fDenInvMass = new TH1D(*aCorrFctn.fDenInvMass);
  else
    fDenInvMass = 0;




}
//____________________________
AliFemtoCorrFctnInvMass::~AliFemtoCorrFctnInvMass(){
  // destructor

  delete fNumInvMass;
  delete fDenInvMass;

}

//_________________________
AliFemtoCorrFctnInvMass& AliFemtoCorrFctnInvMass::operator=(const AliFemtoCorrFctnInvMass& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  fBins = aCorrFctn.fBins;
  fMin = aCorrFctn.fMax;
  fMax = aCorrFctn.fMax;

  fMass1 = aCorrFctn.fMass1;
  fMass2 = aCorrFctn.fMass2;


  // copy constructor
  if (aCorrFctn.fNumInvMass)
    fNumInvMass = new TH1D(*aCorrFctn.fNumInvMass);
  else
    fNumInvMass= 0;
  if (aCorrFctn.fDenInvMass)
    fDenInvMass = new TH1D(*aCorrFctn.fDenInvMass);
  else
    fDenInvMass = 0;


  return *this;
}
//_________________________
void AliFemtoCorrFctnInvMass::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnInvMass::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fNumInvMass->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDenInvMass->GetEntries());
  stemp += ctemp;

  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnInvMass::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }



  double px1 = 999999;
  double py1 = 999999;
  double pz1 = 999999;
  if(pair->Track1()->Track())
  {
    px1 = pair->Track1()->Track()->P().x();
    py1 = pair->Track1()->Track()->P().y();
    pz1 = pair->Track1()->Track()->P().z();
  }
  else if(pair->Track1()->V0())
  {
    px1 = pair->Track1()->V0()->MomV0().x();
    py1 = pair->Track1()->V0()->MomV0().y();
    pz1 = pair->Track1()->V0()->MomV0().z();
  }
  double mass1 = fMass1;


  
  double E1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+mass1*mass1); 

  double px2 = 999999;
  double py2 = 999999;
  double pz2 = 999999;
  if(pair->Track2()->Track())
  {
    px2 = pair->Track2()->Track()->P().x();
    py2 = pair->Track2()->Track()->P().y();
    pz2 = pair->Track2()->Track()->P().z();
  }
  else if(pair->Track2()->V0())
  {
    px2 = pair->Track2()->V0()->MomV0().x();
    py2 = pair->Track2()->V0()->MomV0().y();
    pz2 = pair->Track2()->V0()->MomV0().z();
  }
  double mass2 = fMass2;
  
  double E2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+mass2*mass2);



  double minv = TMath::Sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));



  fNumInvMass->Fill(minv);





}
//____________________________
void AliFemtoCorrFctnInvMass::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }


  double px1 = 999999;
  double py1 = 999999;
  double pz1 = 999999;
  if(pair->Track1()->Track())
  {
    px1 = pair->Track1()->Track()->P().x();
    py1 = pair->Track1()->Track()->P().y();
    pz1 = pair->Track1()->Track()->P().z();
  }
  else if(pair->Track1()->V0())
  {
    px1 = pair->Track1()->V0()->MomV0().x();
    py1 = pair->Track1()->V0()->MomV0().y();
    pz1 = pair->Track1()->V0()->MomV0().z();
  }
  double mass1 = fMass1;
  
  double E1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+mass1*mass1);


  double px2 = 999999;
  double py2 = 999999;
  double pz2 = 999999;
  if(pair->Track2()->Track())
  {
    px2 = pair->Track2()->Track()->P().x();
    py2 = pair->Track2()->Track()->P().y();
    pz2 = pair->Track2()->Track()->P().z();
  }
  else if(pair->Track2()->V0())
  {
    px2 = pair->Track2()->V0()->MomV0().x();
    py2 = pair->Track2()->V0()->MomV0().y();
    pz2 = pair->Track2()->V0()->MomV0().z();
  }
  double mass2 = fMass2;
  
  double E2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+mass2*mass2);

  

  double minv = TMath::Sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));



  fDenInvMass->Fill(minv);

}


void AliFemtoCorrFctnInvMass::WriteHistos()
{
  // Write out result histograms
  fNumInvMass->Write();
  fDenInvMass->Write();
 

}

TList* AliFemtoCorrFctnInvMass::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumInvMass);
  tOutputList->Add(fDenInvMass);
  

  return tOutputList;

}






void AliFemtoCorrFctnInvMass::SetParticleMasses(double mass1, double mass2)
{
  fMass1=mass1;
  fMass2=mass2;
}


void AliFemtoCorrFctnInvMass::SetParticle1Mass(double mass)
{
  fMass1=mass;
}

void AliFemtoCorrFctnInvMass::SetParticle2Mass(double mass)
{
  fMass2=mass;
}
