///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DSpherical: a class to calculate 3D correlation      //
// for pairs of identical particles, binned in spherical coordinates.    //
// In analysis the function should be first created in a macro, then     //
// added to the analysis, and at the end of the macro the procedure to   //
// write out histograms should be called.                                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctn3DSpherical.h"
#include <TMath.h>
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctn3DSpherical)
#endif

//____________________________
  AliFemtoCorrFctn3DSpherical::AliFemtoCorrFctn3DSpherical(char* title, const int& nqbins, const float& QLo, const float& QHi, const int& nphibins, const int& ncthetabins):
  fNumerator(0),
  fDenominator(0),
  fPairCut(0x0)
{
  // set up numerator
  char tTitNum[100] = "Num";
  strcat(tTitNum,title);
  fNumerator = new TH3D(tTitNum,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);
  // set up denominator
  char tTitDen[100] = "Den";
  strcat(tTitDen,title);
  fDenominator = new TH3D(tTitDen,title,nqbins,QLo,QHi,nphibins,-TMath::Pi(),TMath::Pi(),ncthetabins,-1.0,1.0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}

AliFemtoCorrFctn3DSpherical::AliFemtoCorrFctn3DSpherical(const AliFemtoCorrFctn3DSpherical& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  fPairCut(0x0)
{
  // Copy constructor
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  fPairCut = aCorrFctn.fPairCut;
}
//____________________________
AliFemtoCorrFctn3DSpherical::~AliFemtoCorrFctn3DSpherical(){
  // Destructor
  delete fNumerator;
  delete fDenominator;
}
//_________________________
AliFemtoCorrFctn3DSpherical& AliFemtoCorrFctn3DSpherical::operator=(const AliFemtoCorrFctn3DSpherical& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  
  fPairCut = aCorrFctn.fPairCut;
  
  return *this;
}

//_________________________
void AliFemtoCorrFctn3DSpherical::WriteOutHistos(){
  // Write out all histograms to file
  fNumerator->Write();
  fDenominator->Write();
}
//______________________________
TList* AliFemtoCorrFctn3DSpherical::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator); 
  tOutputList->Add(fDenominator);  

  return tOutputList;
}

//_________________________
void AliFemtoCorrFctn3DSpherical::Finish(){
  // here is where we should normalize, fit, etc...
}

//____________________________
AliFemtoString AliFemtoCorrFctn3DSpherical::Report(){
  // Construct the report
  string stemp = "PRF Frame Spherical 3D Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;

  if (fPairCut){
    sprintf(ctemp,"Here is the PairCut specific to this CorrFctn\n");
    stemp += ctemp;
    stemp += fPairCut->Report();
  }
  else{
    sprintf(ctemp,"No PairCut specific to this CorrFctn\n");
    stemp += ctemp;
  }

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctn3DSpherical::AddRealPair( AliFemtoPair* pair){
  // perform operations on real pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  double tKO = pair->KOut();
  double tKS = pair->KSide();
  double tKL = pair->KLong();

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fNumerator->Fill(tKR,tKP,tKC);
}
//____________________________
void AliFemtoCorrFctn3DSpherical::AddMixedPair( AliFemtoPair* pair){
  // perform operations on mixed pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  double tKO = pair->KOut();
  double tKS = pair->KSide();
  double tKL = pair->KLong();

  double tKR = sqrt(tKO*tKO + tKS*tKS + tKL*tKL);
  double tKC;
  if ( fabs(tKR) < 1e-10 ) tKC = 0.0;
  else tKC=tKL/tKR;
  double tKP=atan2(tKS,tKO);

  fDenominator->Fill(tKR,tKP,tKC);
}


