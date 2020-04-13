////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnTrueQ - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctnTrueQ, 1)
#endif

#include <TH1D.h>
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnTrueQ.h"

//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ():
  AliFemtoModelCorrFctn("CF", 100, 0.0, 0.4)
{
}

//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(const char *title, Int_t aNbins, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, 0.0, aQinvHi)
{
}

//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fTrueNum(nullptr),
  fTrueDen(nullptr)
{
  // basic constructor
  TString num_name = Form("%sTrueQNum", title);
  TString den_name = Form("%sTrueQDen", title);

  fTrueNum = new TH1D(num_name, num_name + "; q_{inv}", aNbins,aQinvLo,aQinvHi);
  fTrueDen = new TH1D(den_name, den_name + "; q_{inv}", aNbins,aQinvLo,aQinvHi);

  fTrueNum->Sumw2();
  fTrueDen->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(const AliFemtoModelCorrFctnTrueQ& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fTrueNum(nullptr),
  fTrueDen(nullptr)
{
  // copy constructor
  fTrueNum = new TH1D(*aCorrFctn.fTrueNum);
  fTrueDen = new TH1D(*aCorrFctn.fTrueDen);
}
//_______________________
AliFemtoModelCorrFctnTrueQ::~AliFemtoModelCorrFctnTrueQ()
{
  // destructor
  delete fTrueNum;
  delete fTrueDen;
}

//_______________________
AliFemtoModelCorrFctnTrueQ& AliFemtoModelCorrFctnTrueQ::operator=(const AliFemtoModelCorrFctnTrueQ& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  *fTrueNum = *aCorrFctn.fTrueNum;
  *fTrueDen = *aCorrFctn.fTrueDen;

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnTrueQ::Report()
{
  // construct report
  AliFemtoString report = "AliFemtoModelCorrFctnTrueQ report:\n";
  report += Form("Number of entries in numerator: \t%E\n", fTrueNum->GetEntries());
  report += Form("Number of entries in denominator: \t%E\n", fTrueDen->GetEntries());
  report += AliFemtoModelCorrFctn::Report();

  return report;
}

//_______________________
void AliFemtoModelCorrFctnTrueQ::AddRealPair(AliFemtoPair* aPair)
{
  // add real (effect) pair
  AliFemtoModelCorrFctn::AddRealPair(aPair);
  fTrueNum->Fill(fManager->GetWeightGenerator()->GetKStar()*2);
}
//_______________________
void AliFemtoModelCorrFctnTrueQ::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  AliFemtoModelCorrFctn::AddMixedPair(aPair);
  // save the generated positions
  fTrueDen->Fill(fManager->GetWeightGenerator()->GetKStar()*2);
}
//_______________________
void AliFemtoModelCorrFctnTrueQ::Write()
{
  // write out all the histograms
  fTrueNum->Write();
  fTrueDen->Write();

  AliFemtoModelCorrFctn::Write();
}
//_______________________
TList* AliFemtoModelCorrFctnTrueQ::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();

  tOutputList->Add(fTrueNum);
  tOutputList->Add(fTrueDen);

  return tOutputList;
}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctnTrueQ::Clone() const
{
  // Clone the correlation function
  AliFemtoModelCorrFctnTrueQ *tCopy = new AliFemtoModelCorrFctnTrueQ(*this);

  return tCopy;
}
