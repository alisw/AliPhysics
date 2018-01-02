////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDEtaDPhiSimple.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnDEtaDPhiSimple)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnDEtaDPhiSimple::AliFemtoCorrFctnDEtaDPhiSimple(const char* title, const int& aPhiBins=20, const int& aEtaBins=20):
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaHiddenNumerator(0),
  fDPhiDEtaHiddenDenominator(0),
  fDPhiDEtaHiddenPrimaryNumerator(0),
  fDPhiDEtaHiddenPrimaryDenominator(0),
  fDPhiDEtaHiddenSecWeakNumerator(0),
  fDPhiDEtaHiddenSecWeakDenominator(0),
  fDPhiDEtaHiddenSecMatNumerator(0),
  fDPhiDEtaHiddenSecMatDenominator(0),
  fDPhiDEtaHiddenPrimaryNumeratorData(0),
  fDPhiDEtaHiddenPrimaryDenominatorData(0),
  fDPhiDEtaHiddenSecWeakNumeratorData(0),
  fDPhiDEtaHiddenSecWeakDenominatorData(0),
  fDPhiDEtaHiddenSecMatNumeratorData(0),
  fDPhiDEtaHiddenSecMatDenominatorData(0),
  fphiL(0),
  fphiT(0),
  fEtaBins(0),
  fPhiBins(0),
  fTitle(title),
  fReadHiddenInfo(false)
{

  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;

  fEtaBins = aEtaBins;
  fPhiBins = aPhiBins;

  // set up numerator
  char tTitNumD[101] = "NumDPhiDEta";
  strncat(tTitNumD,title, 100);
  fDPhiDEtaNumerator = new TH2D(tTitNumD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  char tTitDenD[101] = "DenDPhiDEta";
  strncat(tTitDenD,title, 100);
  fDPhiDEtaDenominator = new TH2D(tTitDenD,title,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);


  // to enable error bar calculation...
  fDPhiDEtaNumerator->Sumw2();
  fDPhiDEtaDenominator->Sumw2();



}

//____________________________
AliFemtoCorrFctnDEtaDPhiSimple::AliFemtoCorrFctnDEtaDPhiSimple(const AliFemtoCorrFctnDEtaDPhiSimple& aCorrFctn) :
  AliFemtoCorrFctn(),
  fDPhiDEtaNumerator(0),
  fDPhiDEtaDenominator(0),
  fDPhiDEtaHiddenNumerator(0),
  fDPhiDEtaHiddenDenominator(0),
  fDPhiDEtaHiddenPrimaryNumerator(0),
  fDPhiDEtaHiddenPrimaryDenominator(0),
  fDPhiDEtaHiddenSecWeakNumerator(0),
  fDPhiDEtaHiddenSecWeakDenominator(0),
  fDPhiDEtaHiddenSecMatNumerator(0),
  fDPhiDEtaHiddenSecMatDenominator(0),
  fDPhiDEtaHiddenPrimaryNumeratorData(0),
  fDPhiDEtaHiddenPrimaryDenominatorData(0),
  fDPhiDEtaHiddenSecWeakNumeratorData(0),
  fDPhiDEtaHiddenSecWeakDenominatorData(0),
  fDPhiDEtaHiddenSecMatNumeratorData(0),
  fDPhiDEtaHiddenSecMatDenominatorData(0),
  fphiL(0),
  fphiT(0),
  fEtaBins(0),
  fPhiBins(0),
  fTitle(aCorrFctn.fTitle),
  fReadHiddenInfo(false)
{
  // copy constructor

  fEtaBins = aCorrFctn.fEtaBins;
  fPhiBins = aCorrFctn.fPhiBins;

  // copy constructor
  if (aCorrFctn.fDPhiDEtaNumerator)
    fDPhiDEtaNumerator = new TH2D(*aCorrFctn.fDPhiDEtaNumerator);
  else
    fDPhiDEtaNumerator = 0;
  if (aCorrFctn.fDPhiDEtaDenominator)
    fDPhiDEtaDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
  else
    fDPhiDEtaDenominator = 0;

  fReadHiddenInfo = aCorrFctn.fReadHiddenInfo;
  if(fReadHiddenInfo)
    {

      if (aCorrFctn.fDPhiDEtaHiddenNumerator)
	fDPhiDEtaHiddenNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenNumerator);
      else
	fDPhiDEtaHiddenNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenDenominator)
	fDPhiDEtaHiddenDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenDenominator);
      else
	fDPhiDEtaHiddenDenominator = 0;


      if (aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator)
	fDPhiDEtaHiddenPrimaryNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator);
      else
	fDPhiDEtaHiddenPrimaryNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator)
	fDPhiDEtaHiddenPrimaryDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator);
      else
	fDPhiDEtaHiddenPrimaryDenominator = 0;


      if (aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator)
	fDPhiDEtaHiddenSecWeakNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator);
      else
	fDPhiDEtaHiddenSecWeakNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator)
	fDPhiDEtaHiddenSecWeakDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator);
      else
	fDPhiDEtaHiddenSecWeakDenominator = 0;

      if (aCorrFctn.fDPhiDEtaHiddenSecMatNumerator)
	fDPhiDEtaHiddenSecMatNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumerator);
      else
	fDPhiDEtaHiddenSecMatNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecMatDenominator)
	fDPhiDEtaHiddenSecMatDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominator);
      else
	fDPhiDEtaHiddenSecMatDenominator = 0;




      if (aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData)
	fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData);
      else
	fDPhiDEtaHiddenPrimaryNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData)
	fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData);
      else
	fDPhiDEtaHiddenPrimaryDenominatorData = 0;


      if (aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData)
	fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData);
      else
	fDPhiDEtaHiddenSecWeakNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData)
	fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData);
      else
	fDPhiDEtaHiddenSecWeakDenominatorData = 0;

      if (aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData)
	fDPhiDEtaHiddenSecMatNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData);
      else
	fDPhiDEtaHiddenSecMatNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData)
	fDPhiDEtaHiddenSecMatDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData);
      else
	fDPhiDEtaHiddenSecMatDenominatorData = 0;
    }



}
//____________________________
AliFemtoCorrFctnDEtaDPhiSimple::~AliFemtoCorrFctnDEtaDPhiSimple(){
  // destructor

  delete fDPhiDEtaNumerator;
  delete fDPhiDEtaDenominator;
  if(fReadHiddenInfo)
    {
      delete fDPhiDEtaHiddenNumerator;
      delete fDPhiDEtaHiddenDenominator;
      delete fDPhiDEtaHiddenPrimaryNumerator;
      delete fDPhiDEtaHiddenPrimaryDenominator;
      delete fDPhiDEtaHiddenSecWeakNumerator;
      delete fDPhiDEtaHiddenSecWeakDenominator;
      delete fDPhiDEtaHiddenSecMatNumerator;
      delete fDPhiDEtaHiddenSecMatDenominator;
      delete fDPhiDEtaHiddenPrimaryNumeratorData;
      delete fDPhiDEtaHiddenPrimaryDenominatorData;
      delete fDPhiDEtaHiddenSecWeakNumeratorData;
      delete fDPhiDEtaHiddenSecWeakDenominatorData;
      delete fDPhiDEtaHiddenSecMatNumeratorData;
      delete fDPhiDEtaHiddenSecMatDenominatorData;
    }





}

//_________________________
AliFemtoCorrFctnDEtaDPhiSimple& AliFemtoCorrFctnDEtaDPhiSimple::operator=(const AliFemtoCorrFctnDEtaDPhiSimple& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  fEtaBins = aCorrFctn.fEtaBins;
  fPhiBins = aCorrFctn.fPhiBins;

  fTitle = aCorrFctn.fTitle;

  if (aCorrFctn.fDPhiDEtaNumerator)
    fDPhiDEtaNumerator = new TH2D(*aCorrFctn.fDPhiDEtaNumerator);
  else
    fDPhiDEtaNumerator = 0;
  if (aCorrFctn.fDPhiDEtaDenominator)
    fDPhiDEtaDenominator = new TH2D(*aCorrFctn.fDPhiDEtaDenominator);
  else
    fDPhiDEtaDenominator = 0;

  fReadHiddenInfo = aCorrFctn.fReadHiddenInfo;
  if(fReadHiddenInfo)
    {

      if (aCorrFctn.fDPhiDEtaHiddenNumerator)
	fDPhiDEtaHiddenNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenNumerator);
      else
	fDPhiDEtaHiddenNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenDenominator)
	fDPhiDEtaHiddenDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenDenominator);
      else
	fDPhiDEtaHiddenDenominator = 0;



      if (aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator)
	fDPhiDEtaHiddenPrimaryNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator);
      else
	fDPhiDEtaHiddenPrimaryNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator)
	fDPhiDEtaHiddenPrimaryDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator);
      else
	fDPhiDEtaHiddenPrimaryDenominator = 0;


      if (aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator)
	fDPhiDEtaHiddenSecWeakNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator);
      else
	fDPhiDEtaHiddenSecWeakNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator)
	fDPhiDEtaHiddenSecWeakDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator);
      else
	fDPhiDEtaHiddenSecWeakDenominator = 0;

      if (aCorrFctn.fDPhiDEtaHiddenSecMatNumerator)
	fDPhiDEtaHiddenSecMatNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumerator);
      else
	fDPhiDEtaHiddenSecMatNumerator = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecMatDenominator)
	fDPhiDEtaHiddenSecMatDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominator);
      else
	fDPhiDEtaHiddenSecMatDenominator = 0;




      if (aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData)
	fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData);
      else
	fDPhiDEtaHiddenPrimaryNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData)
	fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData);
      else
	fDPhiDEtaHiddenPrimaryDenominatorData = 0;


      if (aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData)
	fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData);
      else
	fDPhiDEtaHiddenSecWeakNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData)
	fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData);
      else
	fDPhiDEtaHiddenSecWeakDenominatorData = 0;

      if (aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData)
	fDPhiDEtaHiddenSecMatNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData);
      else
	fDPhiDEtaHiddenSecMatNumeratorData = 0;
      if (aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData)
	fDPhiDEtaHiddenSecMatDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData);
      else
	fDPhiDEtaHiddenSecMatDenominatorData = 0;


    }


  return *this;
}
//_________________________
void AliFemtoCorrFctnDEtaDPhiSimple::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnDEtaDPhiSimple::Report(){
  // create report
  string stemp = "TPC Ncls Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fDPhiDEtaNumerator->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDPhiDEtaDenominator->GetEntries());
  stemp += ctemp;

  if(fReadHiddenInfo)
    {
      snprintf(ctemp , 100, "Number of entries in hidden numerator:\t%E\n",fDPhiDEtaHiddenNumerator->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden denominator:\t%E\n",fDPhiDEtaHiddenDenominator->GetEntries());
      stemp += ctemp;

      snprintf(ctemp , 100, "Number of entries in hidden primary numerator:\t%E\n",fDPhiDEtaHiddenPrimaryNumerator->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden primary denominator:\t%E\n",fDPhiDEtaHiddenPrimaryDenominator->GetEntries());
      stemp += ctemp;

      snprintf(ctemp , 100, "Number of entries in hidden second. weak numerator:\t%E\n",fDPhiDEtaHiddenSecWeakNumerator->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden second. weak denominator:\t%E\n",fDPhiDEtaHiddenSecWeakDenominator->GetEntries());
      stemp += ctemp;

      snprintf(ctemp , 100, "Number of entries in hidden second. material numerator:\t%E\n",fDPhiDEtaHiddenSecMatNumerator->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden second. material denominator:\t%E\n",fDPhiDEtaHiddenSecMatDenominator->GetEntries());
      stemp += ctemp;


      snprintf(ctemp , 100, "Number of entries in hidden primary numerator data:\t%E\n",fDPhiDEtaHiddenPrimaryNumeratorData->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden primary denominator data:\t%E\n",fDPhiDEtaHiddenPrimaryDenominatorData->GetEntries());
      stemp += ctemp;

      snprintf(ctemp , 100, "Number of entries in hidden second. weak numerator data:\t%E\n",fDPhiDEtaHiddenSecWeakNumeratorData->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden second. weak denominator data:\t%E\n",fDPhiDEtaHiddenSecWeakDenominatorData->GetEntries());
      stemp += ctemp;

      snprintf(ctemp , 100, "Number of entries in hidden second. material numerator data:\t%E\n",fDPhiDEtaHiddenSecMatNumeratorData->GetEntries());
      stemp += ctemp;
      snprintf(ctemp , 100, "Number of entries in hidden second. material denominator data:\t%E\n",fDPhiDEtaHiddenSecMatDenominatorData->GetEntries());
      stemp += ctemp;
    }

  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnDEtaDPhiSimple::AddRealPair( AliFemtoPair* pair){
  // add real (effect) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  /*double phi1 = pair->Track1()->Track()->P().Phi();
    double phi2 = pair->Track2()->Track()->P().Phi();
    double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    double eta2 = pair->Track2()->Track()->P().PseudoRapidity();*/

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double deta = eta1 - eta2;


  fDPhiDEtaNumerator->Fill(dphi, deta);

   if(fReadHiddenInfo)
    {
      AliFemtoModelHiddenInfo* hInfo1 = 0;
      AliFemtoModelHiddenInfo* hInfo2 = 0;
      if(pair->Track1()->Track())
	{
	    hInfo1 = (AliFemtoModelHiddenInfo*)pair->Track1()->Track()->GetHiddenInfo();
	}
      if(pair->Track1()->V0())
	{
	    hInfo1 = (AliFemtoModelHiddenInfo*)pair->Track1()->V0()->GetHiddenInfo();
	}

      if(pair->Track2()->Track())
	{
	    hInfo2 = (AliFemtoModelHiddenInfo*)pair->Track2()->Track()->GetHiddenInfo();
	}
      if(pair->Track2()->V0())
	{
	    hInfo2 = (AliFemtoModelHiddenInfo*)pair->Track2()->V0()->GetHiddenInfo();
	}

      if(hInfo1 && hInfo2)
	{
	  AliFemtoThreeVector *v1 = hInfo1->GetTrueMomentum();
	  AliFemtoThreeVector *v2 = hInfo2->GetTrueMomentum();

	  double hphi1 = v1->Phi();
	  double hphi2 = v2->Phi();
	  double heta1 = v1->PseudoRapidity();
	  double heta2 = v2->PseudoRapidity();


	  double dhphi = hphi1 - hphi2;
	  while (dhphi<fphiL) dhphi+=PIT;
	  while (dhphi>fphiT) dhphi-=PIT;

	  double dheta = heta1 - heta2;

	  fDPhiDEtaHiddenNumerator->Fill(dhphi, dheta);
	  if(hInfo1->GetOrigin()==0 && hInfo2->GetOrigin()==0)
	    {
	      fDPhiDEtaHiddenPrimaryNumerator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenPrimaryNumeratorData->Fill(dphi,deta);
	    }
	  else if(hInfo1->GetOrigin()==1 || hInfo2->GetOrigin()==1)
	    {
	      fDPhiDEtaHiddenSecWeakNumerator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenSecWeakNumeratorData->Fill(dphi,deta);
	    }
	  else if(hInfo1->GetOrigin()==2 || hInfo2->GetOrigin()==2)
	    {
	      fDPhiDEtaHiddenSecMatNumerator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenSecMatNumeratorData->Fill(dphi,deta);
	    }

	}
    }



}
//____________________________
void AliFemtoCorrFctnDEtaDPhiSimple::AddMixedPair( AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;

  /*double phi1 = pair->Track1()->Track()->P().Phi();
  double phi2 = pair->Track2()->Track()->P().Phi();
  double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = pair->Track2()->Track()->P().PseudoRapidity();*/

  double phi1 = pair->Track1()->FourMomentum().Phi();
  double phi2 = pair->Track2()->FourMomentum().Phi();
  double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

  double dphi = phi1 - phi2;
  while (dphi<fphiL) dphi+=PIT;
  while (dphi>fphiT) dphi-=PIT;

  double deta = eta1 - eta2;

  fDPhiDEtaDenominator->Fill(dphi, deta);


  if(fReadHiddenInfo)
    {
      AliFemtoModelHiddenInfo* hInfo1 = 0;
      AliFemtoModelHiddenInfo* hInfo2 = 0;
      if(pair->Track1()->Track())
	{
	    hInfo1 = (AliFemtoModelHiddenInfo*)pair->Track1()->Track()->GetHiddenInfo();
	}
      if(pair->Track1()->V0())
	{
	    hInfo1 = (AliFemtoModelHiddenInfo*)pair->Track1()->V0()->GetHiddenInfo();
	}

      if(pair->Track2()->Track())
	{
	    hInfo2 = (AliFemtoModelHiddenInfo*)pair->Track2()->Track()->GetHiddenInfo();
	}
      if(pair->Track2()->V0())
	{
	    hInfo2 = (AliFemtoModelHiddenInfo*)pair->Track2()->V0()->GetHiddenInfo();
	}


      if(hInfo1 && hInfo2)
	{
	  AliFemtoThreeVector *v1 = hInfo1->GetTrueMomentum();
	  AliFemtoThreeVector *v2 = hInfo2->GetTrueMomentum();


	  double hphi1 = v1->Phi();
	  double hphi2 = v2->Phi();
	  double heta1 = v1->PseudoRapidity();
	  double heta2 = v2->PseudoRapidity();


	  double dhphi = hphi1 - hphi2;
	  while (dhphi<fphiL) dhphi+=PIT;
	  while (dhphi>fphiT) dhphi-=PIT;

	  double dheta = heta1 - heta2;


	  fDPhiDEtaHiddenDenominator->Fill(dhphi, dheta);
	  if(hInfo1->GetOrigin()==0 && hInfo2->GetOrigin()==0)
	    {
	      fDPhiDEtaHiddenPrimaryDenominator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenPrimaryDenominatorData->Fill(dphi,deta);
	    }
	  else if(hInfo1->GetOrigin()==1 || hInfo2->GetOrigin()==1)
	    {
	      fDPhiDEtaHiddenSecWeakDenominator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenSecWeakDenominatorData->Fill(dphi,deta);
	    }
	  else if(hInfo1->GetOrigin()==2 || hInfo2->GetOrigin()==2)
	    {
	      fDPhiDEtaHiddenSecMatDenominator->Fill(dhphi,dheta);
	      fDPhiDEtaHiddenSecMatDenominator->Fill(dphi,deta);
	    }


	}
    }


}


void AliFemtoCorrFctnDEtaDPhiSimple::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
    if(fReadHiddenInfo)
    {
      fDPhiDEtaHiddenNumerator->Write();
      fDPhiDEtaHiddenDenominator->Write();

      fDPhiDEtaHiddenPrimaryNumerator->Write();
      fDPhiDEtaHiddenPrimaryDenominator->Write();

      fDPhiDEtaHiddenSecWeakNumerator->Write();
      fDPhiDEtaHiddenSecWeakDenominator->Write();

      fDPhiDEtaHiddenSecMatNumerator->Write();
      fDPhiDEtaHiddenSecMatDenominator->Write();

      fDPhiDEtaHiddenPrimaryNumeratorData->Write();
      fDPhiDEtaHiddenPrimaryDenominatorData->Write();

      fDPhiDEtaHiddenSecWeakNumeratorData->Write();
      fDPhiDEtaHiddenSecWeakDenominatorData->Write();

      fDPhiDEtaHiddenSecMatNumeratorData->Write();
      fDPhiDEtaHiddenSecMatDenominatorData->Write();
    }

}

TList* AliFemtoCorrFctnDEtaDPhiSimple::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
    if(fReadHiddenInfo)
    {
      tOutputList->Add(fDPhiDEtaHiddenNumerator);
      tOutputList->Add(fDPhiDEtaHiddenDenominator);

      tOutputList->Add(fDPhiDEtaHiddenPrimaryNumerator);
      tOutputList->Add(fDPhiDEtaHiddenPrimaryDenominator);

      tOutputList->Add(fDPhiDEtaHiddenSecWeakNumerator);
      tOutputList->Add(fDPhiDEtaHiddenSecWeakDenominator);

      tOutputList->Add(fDPhiDEtaHiddenSecMatNumerator);
      tOutputList->Add(fDPhiDEtaHiddenSecMatDenominator);

      tOutputList->Add(fDPhiDEtaHiddenPrimaryNumeratorData);
      tOutputList->Add(fDPhiDEtaHiddenPrimaryDenominatorData);

      tOutputList->Add(fDPhiDEtaHiddenSecWeakNumeratorData);
      tOutputList->Add(fDPhiDEtaHiddenSecWeakDenominatorData);

      tOutputList->Add(fDPhiDEtaHiddenSecMatNumeratorData);
      tOutputList->Add(fDPhiDEtaHiddenSecMatDenominatorData);
    }

  return tOutputList;

}






void AliFemtoCorrFctnDEtaDPhiSimple::SetParticleTypes(ParticleType partType1, ParticleType partType2)
{
  part1=partType1;
  part2=partType2;
}


void AliFemtoCorrFctnDEtaDPhiSimple::SetParticle1Type(ParticleType partType)
{
  part1=partType;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetParticle2Type(ParticleType partType)
{
  part2=partType;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetReadHiddenInfo(bool read)
{
  fReadHiddenInfo = read;

  int aEtaBins = fEtaBins;
  int aPhiBins = fPhiBins;

  // set up numerator
  TString tTitHNumD = "NumDPhiDEtaHidden" + fTitle;
  fDPhiDEtaHiddenNumerator = new TH2D(tTitHNumD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHDenD = "DenDPhiDEtaHidden" + fTitle;
  fDPhiDEtaHiddenDenominator = new TH2D(tTitHDenD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);

  // set up numerator
  TString tTitHPrimNumD = "NumDPhiDEtaHiddenPrimary" + fTitle;
  fDPhiDEtaHiddenPrimaryNumerator = new TH2D(tTitHPrimNumD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHPrimDenD = "DenDPhiDEtaHiddenPrimary" + fTitle;
  fDPhiDEtaHiddenPrimaryDenominator = new TH2D(tTitHPrimDenD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);

  // set up numerator
  TString tTitHSWNumD = "NumDPhiDEtaHiddenSecWeak" + fTitle;
  fDPhiDEtaHiddenSecWeakNumerator = new TH2D(tTitHSWNumD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHSWDenD = "DenDPhiDEtaHiddenSecWeak" + fTitle;
  fDPhiDEtaHiddenSecWeakDenominator = new TH2D(tTitHSWDenD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);


  // set up numerator
  TString tTitHSMNumD = "NumDPhiDEtaHiddenSecMat" + fTitle;
  fDPhiDEtaHiddenSecMatNumerator = new TH2D(tTitHSMNumD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHSMDenD = "DenDPhiDEtaHiddenSecMat" + fTitle;
  fDPhiDEtaHiddenSecMatDenominator = new TH2D(tTitHSMDenD,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);




  // set up numerator
  TString tTitHPrimNumDData = "NumDPhiDEtaHiddenPrimaryData" + fTitle;
  fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D(tTitHPrimNumDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHPrimDenDData = "DenDPhiDEtaHiddenPrimaryData" + fTitle;
  fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D(tTitHPrimDenDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);

  // set up numerator
  TString tTitHSWNumDData = "NumDPhiDEtaHiddenSecWeakData" + fTitle;
  fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D(tTitHSWNumDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHSWDenDData = "DenDPhiDEtaHiddenSecWeakData" + fTitle;
  fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D(tTitHSWDenDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);


  // set up numerator
  TString tTitHSMNumDData = "NumDPhiDEtaHiddenSecMatData" + fTitle;
  fDPhiDEtaHiddenSecMatNumeratorData = new TH2D(tTitHSMNumDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);
  // set up denominator
  TString tTitHSMDenDData = "DenDPhiDEtaHiddenSecMatData" + fTitle;
  fDPhiDEtaHiddenSecMatDenominatorData = new TH2D(tTitHSMDenDData,fTitle,aPhiBins,fphiL,fphiT,aEtaBins,-2.0,2.0);


  fDPhiDEtaHiddenNumerator->Sumw2();
  fDPhiDEtaHiddenDenominator->Sumw2();
  fDPhiDEtaHiddenPrimaryNumerator->Sumw2();
  fDPhiDEtaHiddenPrimaryDenominator->Sumw2();
  fDPhiDEtaHiddenSecWeakNumerator->Sumw2();
  fDPhiDEtaHiddenSecWeakDenominator->Sumw2();
  fDPhiDEtaHiddenSecMatNumerator->Sumw2();
  fDPhiDEtaHiddenSecMatDenominator->Sumw2();
  fDPhiDEtaHiddenPrimaryNumeratorData->Sumw2();
  fDPhiDEtaHiddenPrimaryDenominatorData->Sumw2();
  fDPhiDEtaHiddenSecWeakNumeratorData->Sumw2();
  fDPhiDEtaHiddenSecWeakDenominatorData->Sumw2();
  fDPhiDEtaHiddenSecMatNumeratorData->Sumw2();
  fDPhiDEtaHiddenSecMatDenominatorData->Sumw2();

}
