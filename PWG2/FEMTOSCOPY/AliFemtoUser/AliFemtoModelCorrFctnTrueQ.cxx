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

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnTrueQ.h"
    
//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(): 
  AliFemtoModelCorrFctn(),
  fTrueNum(0),
  fTrueDen(0)
{
  // default constructor
  char buf[100];
  char title[100] = "CFTrueQ";
  sprintf(buf, "%sNum", title);
  fTrueNum = new TH1D(buf,buf,100,0.0,0.4);
  sprintf(buf, "%sDen", title);
  fTrueDen = new TH1D(buf,buf,100,0.0,0.4);

  fTrueNum->Sumw2();
  fTrueDen->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fTrueNum(0),
  fTrueDen(0)
{
  // basic constructor
  char buf[100];
  sprintf(buf, "%sTrueQNum", title);
  fTrueNum = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  sprintf(buf, "%sTrueQDen", title);
  fTrueDen = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);

  fTrueNum->Sumw2();
  fTrueDen->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnTrueQ::AliFemtoModelCorrFctnTrueQ(const AliFemtoModelCorrFctnTrueQ& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fTrueNum(0),
  fTrueDen(0)
{
  // copy constructor
  fTrueNum = new TH1D(*aCorrFctn.fTrueNum);
  fTrueDen = new TH1D(*aCorrFctn.fTrueDen);
}
//_______________________
AliFemtoModelCorrFctnTrueQ::~AliFemtoModelCorrFctnTrueQ()
{
  // destructor
  if (fTrueNum) delete fTrueNum;
  if (fTrueDen) delete fTrueDen;
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}

//_______________________
AliFemtoModelCorrFctnTrueQ& AliFemtoModelCorrFctnTrueQ::operator=(const AliFemtoModelCorrFctnTrueQ& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) 
    return *this;
  if (aCorrFctn.fTrueNum)
    fTrueNum = new TH1D (*aCorrFctn.fTrueNum);
  else fTrueNum = 0;
  if (aCorrFctn.fTrueDen)
    fTrueDen = new TH1D(*aCorrFctn.fTrueDen);
  else fTrueDen = 0;

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnTrueQ::Report()
{
  // construct report
  AliFemtoString tStr = "AliFemtoModelCorrFctnTrueQ report";

  return tStr;
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
AliFemtoModelCorrFctn* AliFemtoModelCorrFctnTrueQ::Clone()
{
  // Clone the correlation function
  AliFemtoModelCorrFctnTrueQ *tCopy = new AliFemtoModelCorrFctnTrueQ(*this);
  
  return tCopy;
}

