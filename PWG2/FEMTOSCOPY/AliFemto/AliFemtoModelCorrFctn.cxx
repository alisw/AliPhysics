////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctn - the base class for correlation function which    ///
/// uses the model framework and weight generation                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctn, 1)
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include <TH1D.h>
    
//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn(): 
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0)
{
  // Default constructor
  fNumeratorTrue = new TH1D("ModelNumTrue","ModelNumTrue",50,0.0,0.5);
  fNumeratorFake = new TH1D("ModelNumFake","ModelNumFake",50,0.0,0.5);
  fDenominator = new TH1D("ModelDen","ModelDen",50,0.0,0.5);

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();
}
//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0)
{
  // Normal constructor
  char buf[100];
  sprintf(buf, "NumTrue%s", title);
  fNumeratorTrue = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  sprintf(buf, "NumFake%s", title);
  fNumeratorFake = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  sprintf(buf, "Den%s", title);
  fDenominator = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();
}
//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn(const AliFemtoModelCorrFctn& aCorrFctn) :
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0)
{
  // Copy constructor
  if (aCorrFctn.fNumeratorTrue)
    fNumeratorTrue = new TH1D(*(aCorrFctn.fNumeratorTrue));
  if (aCorrFctn.fNumeratorFake)
    fNumeratorFake = new TH1D(*(aCorrFctn.fNumeratorFake));
  if (aCorrFctn.fDenominator)
    fDenominator = new TH1D(*(aCorrFctn.fDenominator));
  fManager = aCorrFctn.fManager;
}
//_______________________
AliFemtoModelCorrFctn::~AliFemtoModelCorrFctn()
{
  // Destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}
//_______________________
AliFemtoModelCorrFctn& AliFemtoModelCorrFctn::operator=(const AliFemtoModelCorrFctn& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;
  
  if (aCorrFctn.fNumeratorTrue)
    fNumeratorTrue = new TH1D(*(aCorrFctn.fNumeratorTrue));
  else
    fNumeratorTrue = 0;
  if (aCorrFctn.fNumeratorFake)
    fNumeratorFake = new TH1D(*(aCorrFctn.fNumeratorFake));
  else
    fNumeratorFake = 0;
  if (aCorrFctn.fDenominator)
    fDenominator = new TH1D(*(aCorrFctn.fDenominator));
  else
    fDenominator = 0;
  fManager = aCorrFctn.fManager;

  return *this;
}
//_______________________
void AliFemtoModelCorrFctn::ConnectToManager(AliFemtoModelManager *aManager)
{
  fManager = aManager;
}

//_______________________
AliFemtoString AliFemtoModelCorrFctn::Report()
{
  // Prepare report
  AliFemtoString tStr = "AliFemtoModelCorrFctn report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctn::AddRealPair(AliFemtoPair* aPair)
{
  Double_t weight = fManager->GetWeight(aPair);
  fNumeratorTrue->Fill(aPair->QInv(), weight);
}
//_______________________
void AliFemtoModelCorrFctn::AddMixedPair(AliFemtoPair* aPair)
{
  Double_t weight = fManager->GetWeight(aPair);
  fNumeratorFake->Fill(aPair->QInv(), weight);
  fDenominator->Fill(aPair->QInv(), 1.0);
}
//_______________________
void AliFemtoModelCorrFctn::EventBegin(const AliFemtoEvent* aEvent)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::EventEnd(const AliFemtoEvent* aEvent)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::Finish()
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::Write()
{
  // Write out data histos
  fNumeratorTrue->Write();
  fNumeratorFake->Write();
  fDenominator->Write();
}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctn::Clone()
{
  // Create clone
  AliFemtoModelCorrFctn *tCopy = new AliFemtoModelCorrFctn(*this);
  
  return tCopy;
}

