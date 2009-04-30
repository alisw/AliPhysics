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
  fCYlmFake(0)
{
  // default constructor

  fCYlmTrue = new AliFemtoCorrFctnDirectYlm();
  fCYlmFake = new AliFemtoCorrFctnDirectYlm();
}
//_______________________
AliFemtoModelCorrFctnDirectYlm::AliFemtoModelCorrFctnDirectYlm(const char *title, Int_t aMaxL, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fCYlmTrue(0),
  fCYlmFake(0)
{
  // basic constructor
  char fname[1000];
  sprintf(fname, "%s%s", title, "True");
  fCYlmTrue = new AliFemtoCorrFctnDirectYlm(fname, aMaxL, aNbins, aQinvLo, aQinvHi);
  sprintf(fname, "%s%s", title, "Fake");
  fCYlmFake = new AliFemtoCorrFctnDirectYlm(fname, aMaxL, aNbins, aQinvLo, aQinvHi);
}
//_______________________
AliFemtoModelCorrFctnDirectYlm::AliFemtoModelCorrFctnDirectYlm(const AliFemtoModelCorrFctnDirectYlm& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fCYlmTrue(0),
  fCYlmFake(0)
{
  // copy constructor
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

  if (aCorrFctn.fCYlmTrue)
    fCYlmTrue = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmTrue->Clone());
  else fCYlmTrue = 0;

  if (aCorrFctn.fCYlmFake)
    fCYlmFake = dynamic_cast<AliFemtoCorrFctnDirectYlm*>(aCorrFctn.fCYlmFake->Clone());
  else fCYlmFake = 0;

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
  
  fCYlmTrue->AddRealPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), weight);
}
//_______________________
void AliFemtoModelCorrFctnDirectYlm::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  if (fPairCut)
    if (!(fPairCut->Pass(aPair))) return;

  Double_t weight = fManager->GetWeight(aPair);

  fCYlmTrue->AddMixedPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), 1.0);
  fCYlmFake->AddRealPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), weight);
  fCYlmFake->AddMixedPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), 1.0);
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

void AliFemtoModelCorrFctnDirectYlm::Finish()
{
  fCYlmTrue->Finish();
  fCYlmFake->Finish();
}
