////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnDirectYlm - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctnDirectYlm, 1)
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
    
//_______________________
AliFemtoModelCorrFctnDirectYlm::AliFemtoModelCorrFctnDirectYlm(): 
  AliFemtoModelCorrFctn(),
  fCYlmTrue(0),
  fCYlmFake(0),
  fUseLCMS(0)
{
  // default constructor

  fCYlmTrue = new AliFemtoCorrFctnDirectYlm();
  fCYlmFake = new AliFemtoCorrFctnDirectYlm();
  fCYlmTrue->SetUseLCMS(fUseLCMS);
  fCYlmFake->SetUseLCMS(fUseLCMS);
}
//_______________________
AliFemtoModelCorrFctnDirectYlm::AliFemtoModelCorrFctnDirectYlm(const char *title, Int_t aMaxL, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi, int aUseLCMS=0):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fCYlmTrue(0),
  fCYlmFake(0),
  fUseLCMS(aUseLCMS)
{
  // basic constructor
  char fname[1000];
  snprintf(fname, 1000, "%s%s", title, "True");
  fCYlmTrue = new AliFemtoCorrFctnDirectYlm(fname, aMaxL, aNbins, aQinvLo, aQinvHi, fUseLCMS);
  snprintf(fname, 1000, "%s%s", title, "Fake");
  fCYlmFake = new AliFemtoCorrFctnDirectYlm(fname, aMaxL, aNbins, aQinvLo, aQinvHi, fUseLCMS);
}
//_______________________
AliFemtoModelCorrFctnDirectYlm::AliFemtoModelCorrFctnDirectYlm(const AliFemtoModelCorrFctnDirectYlm& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fCYlmTrue(0),
  fCYlmFake(0),
  fUseLCMS(0)
{
  // copy constructor
  fUseLCMS = aCorrFctn.fUseLCMS;
  fCYlmTrue = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmTrue->Clone());
  fCYlmFake = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmFake->Clone());
}
//_______________________
AliFemtoModelCorrFctnDirectYlm::~AliFemtoModelCorrFctnDirectYlm()
{
  // destructor
  if (fCYlmTrue) delete fCYlmTrue;
  if (fCYlmFake) delete fCYlmFake;
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}

//_______________________
AliFemtoModelCorrFctnDirectYlm& AliFemtoModelCorrFctnDirectYlm::operator=(const AliFemtoModelCorrFctnDirectYlm& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) 
    return *this;

  fUseLCMS = aCorrFctn.fUseLCMS;

  if (aCorrFctn.fCYlmTrue)
    fCYlmTrue = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmTrue->Clone());
  else fCYlmTrue = 0;

  if (aCorrFctn.fCYlmFake)
    fCYlmFake = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmFake->Clone());
  else fCYlmFake = 0;

  if (aCorrFctn.fNumeratorTrue)
    fNumeratorTrue = new TH1D(*aCorrFctn.fNumeratorTrue);
  else
    fNumeratorTrue = 0;

  if (aCorrFctn.fNumeratorFake)
    fNumeratorFake = new TH1D(*aCorrFctn.fNumeratorFake);
  else
    fNumeratorFake = 0;

  if (aCorrFctn.fDenominator)
    fDenominator = new TH1D(*aCorrFctn.fDenominator);
  else
    fDenominator = 0;

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnDirectYlm::Report()
{
  // construct report
  AliFemtoString tStr = "AliFemtoModelCorrFctnDirectYlm report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctnDirectYlm::AddRealPair(AliFemtoPair* aPair)
{
  // add real (effect) pair
  if (fPairCut)
    if (!(fPairCut->Pass(aPair))) return;

  Double_t weight = fManager->GetWeight(aPair);
  
  if (fUseLCMS)
    fCYlmTrue->AddRealPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), weight);
  else
    fCYlmTrue->AddRealPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), weight);
}
//_______________________
void AliFemtoModelCorrFctnDirectYlm::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  if (fPairCut)
    if (!(fPairCut->Pass(aPair))) return;

  Double_t weight = fManager->GetWeight(aPair);

  if (fUseLCMS) {
    fCYlmTrue->AddMixedPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), 1.0);
    fCYlmFake->AddRealPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), weight);
    fCYlmFake->AddMixedPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), 1.0);
  }
  else {
    fCYlmTrue->AddMixedPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), 1.0);
    fCYlmFake->AddRealPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), weight);
    fCYlmFake->AddMixedPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), 1.0);
  }
}
//_______________________
void AliFemtoModelCorrFctnDirectYlm::Write()
{
  // write out all the histograms
  
  fCYlmTrue->Write();
  fCYlmFake->Write();
}
//_______________________
TList* AliFemtoModelCorrFctnDirectYlm::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();
  tOutputList->Clear();

  TList *tListCfTrue = fCYlmTrue->GetOutputList();
    
  TIter nextListCfTrue(tListCfTrue);
  while (TObject *obj = nextListCfTrue()) {
    tOutputList->Add(obj);
  }

  TList *tListCfFake = fCYlmFake->GetOutputList();
    
  TIter nextListCfFake(tListCfFake);
  while (TObject *obj = nextListCfFake()) {
    tOutputList->Add(obj);
  }
//   tOutputList->Add(fCYlmTrue->GetOutputList());
//   tOutputList->Add(fCYlmFake->GetOutputList());

  return tOutputList;
}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctnDirectYlm::Clone()
{
  // Clone the correlation function
  AliFemtoModelCorrFctnDirectYlm *tCopy = new AliFemtoModelCorrFctnDirectYlm(*this);
  
  return tCopy;
}
//_______________________
void AliFemtoModelCorrFctnDirectYlm::Finish()
{
  fCYlmTrue->Finish();
  fCYlmFake->Finish();
}
//_______________________
void AliFemtoModelCorrFctnDirectYlm::SetUseLCMS(int aUseLCMS)
{
  fUseLCMS = aUseLCMS;
  fCYlmTrue->SetUseLCMS(fUseLCMS);
  fCYlmFake->SetUseLCMS(fUseLCMS);
}
//_______________________
int  AliFemtoModelCorrFctnDirectYlm::GetUseLCMS()
{
  return fUseLCMS;
}

