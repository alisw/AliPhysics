/***************************************************************************
 *
 * $Id: AliFemtoQinvCorrFctnEMCIC.cxx  $
 *
 * Author: Nicolas Bock, Ohio State University, bock@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: Calculates of the Qinv Correlation Function, and also
 *              produces histograms to calculate EMCICs          
 *
 ***************************************************************************
 *
 **************************************************************************/


#include "AliFemtoQinvCorrFctnEMCIC.h"
//#include "AliFemtoHisto.h"
#include <cstdio>
#include <TVector2.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoQinvCorrFctnEMCIC)
#endif

//____________________________
AliFemtoQinvCorrFctnEMCIC::AliFemtoQinvCorrFctnEMCIC(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
AliFemtoQinvCorrFctn(title, nbins, QinvLo, QinvHi),
/*fESumReal(0),
  fEMultReal(0),
  fPtMultReal(0),
  fPzMultReal(0),*/
  fESumMix(0),
  fEMultMix(0),
  fPtMultMix(0),
  fPzMultMix(0)

{
  

  // set up emcic histograms
  /*char tTitESum[100] = "ESumReal";
  strncat(tTitESum,title, 100);
  fESumReal = new TH1D(tTitESum,title,nbins,QinvLo,QinvHi);
  char tTitEMult[100] = "EMultReal";
  strncat(tTitEMult,title, 100);
  fEMultReal = new TH1D(tTitEMult,title,nbins,QinvLo,QinvHi);
  char tTitPt[100] = "PtMultReal";
  strncat(tTitPt,title, 100);
  fPtMultReal = new TH1D(tTitPt,title,nbins,QinvLo,QinvHi);
  char tTitPz[100] = "PzMultReal";
  strncat(tTitPz,title, 100);
  fPzMultReal = new TH1D(tTitPz,title,nbins,QinvLo,QinvHi);*/
 
  char tTitESum2[101] = "ESumMix";
  strncat(tTitESum2,title, 100);
  fESumMix = new TH1D(tTitESum2,title,nbins,QinvLo,QinvHi);
  char tTitEMult2[101] = "EMultMix";
  strncat(tTitEMult2,title, 100);
  fEMultMix = new TH1D(tTitEMult2,title,nbins,QinvLo,QinvHi);
  char tTitPt2[101] = "PtMultMix";
  strncat(tTitPt2,title, 100);
  fPtMultMix = new TH1D(tTitPt2,title,nbins,QinvLo,QinvHi);
  char tTitPz2[101] = "PzMultMix";
  strncat(tTitPz2,title, 100);
  fPzMultMix = new TH1D(tTitPz2,title,nbins,QinvLo,QinvHi);



  // to enable error bar calculation...
  
  /*  fESumReal->Sumw2();
  fEMultReal->Sumw2();
  fPtMultReal->Sumw2();
  fPzMultReal->Sumw2();*/
  fESumMix->Sumw2();
  fEMultMix->Sumw2();
  fPtMultMix->Sumw2();
  fPzMultMix->Sumw2();
}

//____________________________
AliFemtoQinvCorrFctnEMCIC::AliFemtoQinvCorrFctnEMCIC(const AliFemtoQinvCorrFctnEMCIC& aCorrFctn) :
  AliFemtoQinvCorrFctn(aCorrFctn),
  /*fESumReal(0),
  fEMultReal(0),
  fPtMultReal(0),
  fPzMultReal(0),*/
  fESumMix(0),
  fEMultMix(0),
  fPtMultMix(0),
  fPzMultMix(0)
{
  // copy constructor
  
  /*fESumReal= new TH1D(*aCorrFctn.fESumReal);
  fEMultReal= new TH1D(*aCorrFctn.fEMultReal);
  fPtMultReal= new TH1D(*aCorrFctn.fPtMultReal);
  fPzMultReal= new TH1D(*aCorrFctn.fPzMultReal);*/
  fESumMix= new TH1D(*aCorrFctn.fESumMix);
  fEMultMix= new TH1D(*aCorrFctn.fEMultMix);
  fPtMultMix= new TH1D(*aCorrFctn.fPtMultMix);
  fPzMultMix= new TH1D(*aCorrFctn.fPzMultMix);

}
//____________________________
AliFemtoQinvCorrFctnEMCIC::~AliFemtoQinvCorrFctnEMCIC(){
  // destructor
  
  /*delete fESumReal;
  delete fEMultReal;
  delete fPtMultReal;
  delete fPzMultReal;*/
  delete fESumMix;
  delete fEMultMix;
  delete fPtMultMix;
  delete fPzMultMix;

}
//_________________________
AliFemtoQinvCorrFctnEMCIC& AliFemtoQinvCorrFctnEMCIC::operator=(const AliFemtoQinvCorrFctnEMCIC& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  /*if (fESumReal) delete fESumReal;
  fESumReal= new TH1D(*aCorrFctn.fESumReal);
  if (fEMultReal) delete fEMultReal;
  fEMultReal= new TH1D(*aCorrFctn.fEMultReal);
  if (fPtMultReal) delete fPtMultReal;
  fPtMultReal= new TH1D(*aCorrFctn.fPtMultReal);
  if (fPzMultReal) delete fPzMultReal;
  fPzMultReal= new TH1D(*aCorrFctn.fPzMultReal);
  if (fESumMix) delete fESumMix;*/
  
  fESumMix= new TH1D(*aCorrFctn.fESumMix);
  if (fEMultMix) delete fEMultMix;
  fEMultMix= new TH1D(*aCorrFctn.fEMultMix);
  if (fPtMultMix) delete fPtMultMix;
  fPtMultMix= new TH1D(*aCorrFctn.fPtMultMix);
  if (fPzMultMix) delete fPzMultMix;
  fPzMultMix= new TH1D(*aCorrFctn.fPzMultMix);

  return *this;
}

//____________________________
void AliFemtoQinvCorrFctnEMCIC::AddRealPair(AliFemtoPair* pair){
  // add true pair
  
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;
  AliFemtoQinvCorrFctn::AddRealPair(pair);
  
 
  //double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  
// The EMCICs are calculated here for real pairs
  /*AliFemtoLorentzVector tMom1 = pair->Track1()->FourMomentum();
  AliFemtoLorentzVector tMom2 = pair->Track2()->FourMomentum();
  double tE1 = tMom1.e();
  double tE2 = tMom2.e();
  double tPz1 = tMom1.pz();
  double tPz2 = tMom2.pz();
  
  TVector2 tPt1;  
  TVector2 tPt2; 
  tPt1.Set(tMom1.px(),tMom1.py());
  tPt2.Set(tMom2.px(),tMom2.py());
  double tPt1DotPt2 = tPt1*tPt2;
  
  fESumReal->Fill(tQinv,tE1+tE2);
  fEMultReal->Fill(tQinv,tE1*tE2);
  fPzMultReal->Fill(tQinv,tPz1*tPz2);
  fPtMultReal->Fill(tQinv,tPt1DotPt2);*/

}
//____________________________
void AliFemtoQinvCorrFctnEMCIC::AddMixedPair(AliFemtoPair* pair){
  // add mixed (background) pair
  if (fPairCut)
    if (!fPairCut->Pass(pair)) return;
  AliFemtoQinvCorrFctn::AddMixedPair(pair);
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...
  

  // The EMCICs are calculated here for mixed pairs
  AliFemtoLorentzVector tMom1 = pair->Track1()->FourMomentum();
  AliFemtoLorentzVector tMom2 = pair->Track2()->FourMomentum();
  double tE1 = tMom1.e();
  double tE2 = tMom2.e();
  double tPz1 = tMom1.pz();
  double tPz2 = tMom2.pz();
  
  TVector2 tPt1;  
  TVector2 tPt2; 
  tPt1.Set(tMom1.px(),tMom1.py());
  tPt2.Set(tMom2.px(),tMom2.py());
  double tPt1DotPt2 = tPt1*tPt2;
  
  fESumMix->Fill(tQinv,tE1+tE2);
  fEMultMix->Fill(tQinv,tE1*tE2);
  fPzMultMix->Fill(tQinv,tPz1*tPz2);
  fPtMultMix->Fill(tQinv,tPt1DotPt2);



}
//____________________________
void AliFemtoQinvCorrFctnEMCIC::Write(){
  // Write out neccessary objects
  AliFemtoQinvCorrFctn::Write();  //Write num and den
  /*fESumReal->Write();
  fEMultReal->Write();
  fPtMultReal->Write();
  fPzMultReal->Write(); */
  fESumMix->Write();
  fEMultMix->Write();
  fPtMultMix->Write();
  fPzMultMix->Write();

}
//______________________________
TList* AliFemtoQinvCorrFctnEMCIC::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList;
  tOutputList = (TList*)AliFemtoQinvCorrFctn::GetOutputList();
  cout << "Getting list from Qinv CF emicic" << endl;
  /*tOutputList->Add(fESumReal);
  tOutputList->Add(fEMultReal);
  tOutputList->Add(fPtMultReal);
  tOutputList->Add(fPzMultReal); */
  tOutputList->Add(fESumMix);
  tOutputList->Add(fEMultMix);
  tOutputList->Add(fPtMultMix);
  tOutputList->Add(fPzMultMix);
  return tOutputList;
}


