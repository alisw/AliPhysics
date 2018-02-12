////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// alifemtomodelcorrfctnsource - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctnSource, 1)
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
    
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(): 
  AliFemtoModelCorrFctn(),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0),
  fHistNumWS(0),
  fHistDenWS(0),
  fUseRPSelection(0)
{
  // default constructor
  char buf[100];
  char title[100] = "CFSource";
  snprintf(buf , 100,  "%sOut", title);
  fHistROut = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sSide", title);
  fHistRSide = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sLong", title);
  fHistRLong = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sInv", title);
  fHistRStar = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sdNdR", title);
  fHistdNdR = new TH1D(buf,buf,100,-50.0,50.0);

  snprintf(buf , 100,  "%sNWS", title);
  fHistNumWS = new TH2D(buf,buf,50,0.0,0.5,100,0.0,2.0);
  snprintf(buf , 100,  "%sDWS", title);
  fHistDenWS = new TH2D(buf,buf,50,0.0,0.5,100,0.0,2.0);

  fHistROut->Sumw2();
  fHistRSide->Sumw2();
  fHistRLong->Sumw2();
  fHistRStar->Sumw2();
  fHistdNdR->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0),
  fHistNumWS(0),
  fHistDenWS(0),
  fUseRPSelection(0)
{
  // basic constructor
  char buf[100];
  snprintf(buf , 100,  "%sOut", title);
  fHistROut = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sSide", title);
  fHistRSide = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sLong", title);
  fHistRLong = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sInv", title);
  fHistRStar = new TH1D(buf,buf,100,-50.0,50.0);
  snprintf(buf , 100,  "%sdNdR", title);
  fHistdNdR = new TH1D(buf,buf,100,-50.0,50.0);

  snprintf(buf , 100,  "%sNWS", title);
  fHistNumWS = new TH2D(buf,buf,50,0.0,0.5,100,0.0,2.0);
  snprintf(buf , 100,  "%sDWS", title);
  fHistDenWS = new TH2D(buf,buf,50,0.0,0.5,100,0.0,2.0);

  fHistROut->Sumw2();
  fHistRSide->Sumw2();
  fHistRLong->Sumw2();
  fHistRStar->Sumw2();
  fHistdNdR->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(const AliFemtoModelCorrFctnSource& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0),
  fHistNumWS(0),
  fHistDenWS(0),
  fUseRPSelection(0)
{
  // copy constructor
  fHistROut = new TH1D (*aCorrFctn.fHistROut);
  fHistRSide = new TH1D(*aCorrFctn.fHistRSide);
  fHistRLong = new TH1D(*aCorrFctn.fHistRLong);
  fHistRStar = new TH1D(*aCorrFctn.fHistRStar);
  fHistdNdR = new TH1D(*aCorrFctn.fHistdNdR);
  fHistNumWS = new TH2D(*aCorrFctn.fHistNumWS);
  fHistDenWS = new TH2D(*aCorrFctn.fHistDenWS);

  fUseRPSelection = aCorrFctn.fUseRPSelection;
}
//_______________________
AliFemtoModelCorrFctnSource::~AliFemtoModelCorrFctnSource()
{
  // destructor
  if (fHistROut) delete fHistROut;
  if (fHistRSide) delete fHistRSide;
  if (fHistRLong) delete fHistRLong;
  if (fHistRStar) delete fHistRStar;
  if (fHistdNdR) delete fHistdNdR;
  if (fHistNumWS) delete fHistNumWS;
  if (fHistDenWS) delete fHistDenWS;
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}

//_______________________
AliFemtoModelCorrFctnSource& AliFemtoModelCorrFctnSource::operator=(const AliFemtoModelCorrFctnSource& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) 
    return *this;
  if (aCorrFctn.fHistROut)
    fHistROut = new TH1D (*aCorrFctn.fHistROut);
  else fHistROut = 0;
  if (aCorrFctn.fHistRSide)
    fHistRSide = new TH1D(*aCorrFctn.fHistRSide);
  else fHistRSide = 0;
  if (aCorrFctn.fHistRLong)
    fHistRLong = new TH1D(*aCorrFctn.fHistRLong);
  else fHistRLong = 0;
  if (aCorrFctn.fHistRStar)
    fHistRStar = new TH1D(*aCorrFctn.fHistRStar);
  fHistRStar = 0;
  if (aCorrFctn.fHistdNdR)
    fHistdNdR = new TH1D(*aCorrFctn.fHistdNdR);
  else fHistdNdR = 0;
  if (aCorrFctn.fHistNumWS)
    fHistNumWS = new TH2D(*aCorrFctn.fHistNumWS);
  else fHistNumWS = 0;
  if (aCorrFctn.fHistDenWS)
    fHistDenWS = new TH2D(*aCorrFctn.fHistDenWS);
  else fHistDenWS = 0;

  fUseRPSelection = aCorrFctn.fUseRPSelection;

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnSource::Report()
{
  // construct report
  AliFemtoString tStr = "AliFemtoModelCorrFctnSource report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctnSource::AddRealPair(AliFemtoPair* aPair)
{
  // add real (effect) pair
//   if (fPairCut){
//     if (!(fPairCut->Pass(aPair))) return;
//   }
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(aPair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(aPair))) return;
	}
	else if (!(ktc->Pass(aPair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(aPair))) return;
  }
  
  AliFemtoModelCorrFctn::AddRealPair(aPair);

}
//_______________________
void AliFemtoModelCorrFctnSource::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
//   if (fPairCut){
//     if (!(fPairCut->Pass(aPair))) return;
//   }
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(aPair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(aPair))) return;
	}
	else if (!(ktc->Pass(aPair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(aPair))) return;
  }
  
  AliFemtoModelCorrFctn::AddMixedPair(aPair);
  // save the generated positions
  if (aPair->KStar() < 0.2) {
    fHistROut->Fill (fManager->GetWeightGenerator()->GetRStarOut());
    fHistRSide->Fill(fManager->GetWeightGenerator()->GetRStarSide());
    fHistRLong->Fill(fManager->GetWeightGenerator()->GetRStarLong());
    fHistRStar->Fill(fManager->GetWeightGenerator()->GetRStar());
    fHistdNdR->Fill (fManager->GetWeightGenerator()->GetRStar(),1.0/(fManager->GetWeightGenerator()->GetRStar()*fManager->GetWeightGenerator()->GetRStar()));
  }

  fHistDenWS->Fill(aPair->QInv(), 1.0);
  Double_t weight = fManager->GetWeight(aPair);
  fHistNumWS->Fill(aPair->QInv(), weight);
}
//_______________________
void AliFemtoModelCorrFctnSource::Write()
{
  // write out all the histograms
  fHistROut->Write();
  fHistRSide->Write();
  fHistRLong->Write();
  fHistRStar->Write();
  fHistdNdR->Write();
  fHistNumWS->Write();
  fHistDenWS->Write();

  AliFemtoModelCorrFctn::Write();
}
//________________________
TList* AliFemtoModelCorrFctnSource::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();

  tOutputList->Add(fHistROut); 
  tOutputList->Add(fHistRSide);  
  tOutputList->Add(fHistRLong);  
  tOutputList->Add(fHistRStar);  
  tOutputList->Add(fHistdNdR);  
  tOutputList->Add(fHistDenWS);
  tOutputList->Add(fHistNumWS);

  return tOutputList;
}
//_______________________
void AliFemtoModelCorrFctnSource::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
