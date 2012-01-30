///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoModelCorrFctn3DLCMSSpherical: a class to calculate 3D correlation //
// for pairs of identical particles, binned in spherical coordinates.    //
// In analysis the function should be first created in a macro, then     //
// added to the analysis, and at the end of the macro the procedure to   //
// write out histograms should be called.                                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoModelCorrFctn3DLCMSSpherical.h"
#include "AliFemtoModelManager.h"
#include <TMath.h>
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoModelCorrFctn3DLCMSSpherical)
#endif

//____________________________
  AliFemtoModelCorrFctn3DLCMSSpherical::AliFemtoModelCorrFctn3DLCMSSpherical(char* title, const int& nqbins, const float& QLo, const float& QHi, const int& nphibins, const int& ncthetabins):
  fTrueNumeratorSph(0),
  fFakeNumeratorSph(0),
  fDenominatorSph(0),
  fPairCut(0x0)
{
  // set up numerator
  char tTitNum[101] = "NumTrue";
  strncat(tTitNum,title, 100);
  fTrueNumeratorSph = new TH3D(tTitNum,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);
  // set up numerator
  char tTitNumF[101] = "NumFake";
  strncat(tTitNumF,title, 100);
  fFakeNumeratorSph = new TH3D(tTitNumF,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);
  // set up denominator
  char tTitDen[101] = "Den";
  strncat(tTitDen,title, 100);
  fDenominatorSph = new TH3D(tTitDen,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);

  // to enable error bar calculation...
  fTrueNumeratorSph->Sumw2();
  fFakeNumeratorSph->Sumw2();
  fDenominatorSph->Sumw2();
}

AliFemtoModelCorrFctn3DLCMSSpherical::AliFemtoModelCorrFctn3DLCMSSpherical(const AliFemtoModelCorrFctn3DLCMSSpherical& aCorrFctn) :
  AliFemtoModelCorrFctn(),
  fTrueNumeratorSph(0),
  fFakeNumeratorSph(0),
  fDenominatorSph(0),
  fPairCut(0x0)
{
  // Copy constructor
  fTrueNumeratorSph = new TH3D(*aCorrFctn.fTrueNumeratorSph);
  fFakeNumeratorSph = new TH3D(*aCorrFctn.fFakeNumeratorSph);
  fDenominatorSph = new TH3D(*aCorrFctn.fDenominatorSph);
  fPairCut = aCorrFctn.fPairCut;
}
//____________________________
AliFemtoModelCorrFctn3DLCMSSpherical::~AliFemtoModelCorrFctn3DLCMSSpherical(){
  // Destructor
  delete fTrueNumeratorSph;
  delete fFakeNumeratorSph;
  delete fDenominatorSph;
}
//_________________________
AliFemtoModelCorrFctn3DLCMSSpherical& AliFemtoModelCorrFctn3DLCMSSpherical::operator=(const AliFemtoModelCorrFctn3DLCMSSpherical& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fTrueNumeratorSph) delete fTrueNumeratorSph;
  fTrueNumeratorSph = new TH3D(*aCorrFctn.fTrueNumeratorSph);
  if (fFakeNumeratorSph) delete fFakeNumeratorSph;
  fFakeNumeratorSph = new TH3D(*aCorrFctn.fFakeNumeratorSph);
  if (fDenominatorSph) delete fDenominatorSph;
  fDenominatorSph = new TH3D(*aCorrFctn.fDenominatorSph);
  
  fPairCut = aCorrFctn.fPairCut;
  
  return *this;
}

//_________________________
void AliFemtoModelCorrFctn3DLCMSSpherical::WriteOutHistos(){
  // Write out all histograms to file
  fTrueNumeratorSph->Write();
  fFakeNumeratorSph->Write();
  fDenominatorSph->Write();
}
//______________________________
TList* AliFemtoModelCorrFctn3DLCMSSpherical::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fTrueNumeratorSph); 
  tOutputList->Add(fFakeNumeratorSph); 
  tOutputList->Add(fDenominatorSph);  

  return tOutputList;
}

//_________________________
void AliFemtoModelCorrFctn3DLCMSSpherical::Finish(){
  // here is where we should normalize, fit, etc...
}

//____________________________
AliFemtoString AliFemtoModelCorrFctn3DLCMSSpherical::Report(){
  // Construct the report
  string stemp = "LCMS Frame Spherical 3D Model Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fTrueNumeratorSph->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDenominatorSph->GetEntries());
  stemp += ctemp;

  if (fPairCut){
    snprintf(ctemp , 100, "Here is the PairCut specific to this CorrFctn\n");
    stemp += ctemp;
    stemp += fPairCut->Report();
  }
  else{
    snprintf(ctemp , 100, "No PairCut specific to this CorrFctn\n");
    stemp += ctemp;
  }

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelCorrFctn3DLCMSSpherical::AddRealPair( AliFemtoPair* pair){
  // perform operations on real pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  Double_t weight = fManager->GetWeight(pair);


  double tKO = (pair->QOutCMS());
  double tKS = (pair->QSideCMS());
  double tKL = (pair->QLongCMS());

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fTrueNumeratorSph->Fill(tKR,tKP,tKC,weight);
}
//____________________________
void AliFemtoModelCorrFctn3DLCMSSpherical::AddMixedPair( AliFemtoPair* pair){
  // perform operations on mixed pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  Double_t weight = fManager->GetWeight(pair);

  double tKO = (pair->QOutCMS());
  double tKS = (pair->QSideCMS());
  double tKL = (pair->QLongCMS());

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fFakeNumeratorSph->Fill(tKR,tKP,tKC,weight);
  fDenominatorSph->Fill(tKR,tKP,tKC);
}


