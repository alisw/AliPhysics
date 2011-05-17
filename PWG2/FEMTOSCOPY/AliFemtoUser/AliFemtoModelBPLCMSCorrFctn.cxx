////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoModelBPLCMSCorrFctn)
#endif

//____________________________
AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi)
  :
  AliFemtoModelCorrFctn(title, nbins, QLo, QHi),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0),
  fPairCut(0),
  fUseRPSelection(0)
{
  // set up true numerator
  char tTitNumT[101] = "Num3DTrue";
  strncat(tTitNumT,title, 100);
  fNumerator3DTrue = new TH3D(tTitNumT,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up fake numerator
  char tTitNumF[101] = "Num3DFake";
  strncat(tTitNumF,title, 100);
  fNumerator3DFake = new TH3D(tTitNumF,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char tTitDen[101] = "Den3D";
  strncat(tTitDen,title, 100);
  fDenominator3D = new TH3D(tTitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up ave qInv
  char tTitQinv[101] = "Qinv";
  strncat(tTitQinv,title, 100);
  fQinvHisto = new TH3D(tTitQinv,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  // to enable error bar calculation...
  fNumerator3DTrue->Sumw2();
  fNumerator3DFake->Sumw2();
  fDenominator3D->Sumw2();
}

AliFemtoModelBPLCMSCorrFctn::AliFemtoModelBPLCMSCorrFctn(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn),
  fNumerator3DTrue(0),
  fNumerator3DFake(0),
  fDenominator3D(0),
  fQinvHisto(0),
  fPairCut(0),
  fUseRPSelection(0)
{
  // Copy constructor
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  fDenominator3D   = new TH3D(*aCorrFctn.fDenominator3D);
  fQinvHisto       = new TH3D(*aCorrFctn.fQinvHisto);
  fPairCut         = aCorrFctn.fPairCut->Clone();
}
//____________________________
AliFemtoModelBPLCMSCorrFctn::~AliFemtoModelBPLCMSCorrFctn()
{
  // destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
  delete fNumerator3DTrue;
  delete fNumerator3DFake;
  delete fDenominator3D;
  delete fQinvHisto;
  if (fPairCut) delete fPairCut;
}
//_________________________
AliFemtoModelBPLCMSCorrFctn& AliFemtoModelBPLCMSCorrFctn::operator=(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;
  if (fNumerator3DTrue) delete fNumerator3DTrue;
  fNumerator3DTrue = new TH3D(*aCorrFctn.fNumerator3DTrue);
  if (fNumerator3DFake) delete fNumerator3DFake;
  fNumerator3DFake = new TH3D(*aCorrFctn.fNumerator3DFake);
  if (fDenominator3D) delete fDenominator3D;
  fDenominator3D = new TH3D(*aCorrFctn.fDenominator3D);
  if (fQinvHisto) delete fQinvHisto;
  fQinvHisto = new TH3D(*aCorrFctn.fQinvHisto);
  fPairCut = aCorrFctn.fPairCut->Clone();

  return *this;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::Write(){
  // Write out data histograms
  AliFemtoModelCorrFctn::Write();
  fNumerator3DTrue->Write();
  fNumerator3DFake->Write();
  fDenominator3D->Write();
  fQinvHisto->Write();
}
//________________________
TList* AliFemtoModelBPLCMSCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();

  tOutputList->Add(fNumerator3DTrue); 
  tOutputList->Add(fNumerator3DFake);  
  tOutputList->Add(fDenominator3D);  
  tOutputList->Add(fQinvHisto);  

  return tOutputList;
}

//_________________________
void AliFemtoModelBPLCMSCorrFctn::Finish(){
  fQinvHisto->Divide(fDenominator);
}

//____________________________
AliFemtoString AliFemtoModelBPLCMSCorrFctn::Report(){
  // Prepare a report from the execution
  string stemp = "LCMS Frame Bertsch-Pratt 3D Model Correlation Function Report:\n";
  char ctemp[100];
  snprintf(ctemp , 100, "Number of entries in numerator:\t%E\n",fNumeratorTrue->GetEntries());
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;
  /*  if (fCorrection)
      {
      float radius = fCorrection->GetRadius();
      snprintf(ctemp , 100, "Coulomb correction used radius of\t%E\n",radius);
      }
      else
      {
      snprintf(ctemp , 100, "No Coulomb Correction applied to this CorrFctn\n");
      }
      stemp += ctemp;
  */

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelBPLCMSCorrFctn::AddRealPair( AliFemtoPair* pair)
{
  // Store a real pair in numerator
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }
  
  Double_t weight = fManager->GetWeight(pair);

  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  fNumerator3DTrue->Fill(qOut, qSide, qLong, weight);
  fNumeratorTrue->Fill(pair->QInv(), weight);
}
//____________________________
void AliFemtoModelBPLCMSCorrFctn::AddMixedPair( AliFemtoPair* pair){
  // store mixed pair in denominator
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }

  Double_t weight = fManager->GetWeight(pair);

  double qOut = (pair->QOutCMS());
  double qSide = (pair->QSideCMS());
  double qLong = (pair->QLongCMS());

  fNumerator3DFake->Fill(qOut, qSide, qLong, weight);
  fDenominator3D->Fill(qOut, qSide, qLong, 1.0);
  fNumeratorFake->Fill(pair->QInv(), weight);
  fDenominator->Fill(pair->QInv(), 1.0);

}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelBPLCMSCorrFctn::Clone()
{
  // Clone the correlation function
  AliFemtoModelBPLCMSCorrFctn *tCopy = new AliFemtoModelBPLCMSCorrFctn(*this);
  
  return tCopy;
}

void AliFemtoModelBPLCMSCorrFctn::SetSpecificPairCut(AliFemtoPairCut* aCut)
{
  fPairCut = aCut;
}

void AliFemtoModelBPLCMSCorrFctn::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
