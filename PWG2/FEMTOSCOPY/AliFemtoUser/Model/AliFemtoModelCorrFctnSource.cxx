////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnSource - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctnSource, 1)
#endif

#include "Model/AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "Model/AliFemtoModelHiddenInfo.h"
#include "Model/AliFemtoModelCorrFctnSource.h"
    
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(): 
  AliFemtoModelCorrFctn(),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0)
{
}
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0)
{
  char buf[100];
  sprintf(buf, "%sOut", title);
  fHistROut = new TH1D(buf,buf,100,-50.0,50.0);
  sprintf(buf, "%sSide", title);
  fHistRSide = new TH1D(buf,buf,100,-50.0,50.0);
  sprintf(buf, "%sLong", title);
  fHistRLong = new TH1D(buf,buf,100,-50.0,50.0);
  sprintf(buf, "%sInv", title);
  fHistRStar = new TH1D(buf,buf,100,-50.0,50.0);
  sprintf(buf, "%sdNdR", title);
  fHistdNdR = new TH1D(buf,buf,100,-50.0,50.0);

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
  fHistdNdR(0)
{
  fHistROut = new TH1D (*aCorrFctn.fHistROut);
  fHistRSide = new TH1D(*aCorrFctn.fHistRSide);
  fHistRLong = new TH1D(*aCorrFctn.fHistRLong);
  fHistRStar = new TH1D(*aCorrFctn.fHistRStar);
  fHistdNdR = new TH1D(*aCorrFctn.fHistdNdR);
}
//_______________________
AliFemtoModelCorrFctnSource::~AliFemtoModelCorrFctnSource()
{
  if (fHistROut) delete fHistROut;
  if (fHistRSide) delete fHistRSide;
  if (fHistRLong) delete fHistRLong;
  if (fHistRStar) delete fHistRStar;
  if (fHistdNdR) delete fHistdNdR;
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}

//_______________________
AliFemtoModelCorrFctnSource& AliFemtoModelCorrFctnSource::operator=(const AliFemtoModelCorrFctnSource& aCorrFctn)
{
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

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnSource::Report()
{
  AliFemtoString tStr = "AliFemtoModelCorrFctnSource report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctnSource::AddRealPair(AliFemtoPair* aPair)
{
  AliFemtoModelCorrFctn::AddRealPair(aPair);
}
//_______________________
void AliFemtoModelCorrFctnSource::AddMixedPair(AliFemtoPair* aPair)
{
  AliFemtoModelCorrFctn::AddMixedPair(aPair);
  fHistROut->Fill (fManager->GetWeightGenerator()->GetRStarOut());
  fHistRSide->Fill(fManager->GetWeightGenerator()->GetRStarSide());
  fHistRLong->Fill(fManager->GetWeightGenerator()->GetRStarLong());
  fHistRStar->Fill(fManager->GetWeightGenerator()->GetRStar());
  fHistdNdR->Fill (fManager->GetWeightGenerator()->GetRStar(),1.0/(fManager->GetWeightGenerator()->GetRStar()*fManager->GetWeightGenerator()->GetRStar()));
}
//_______________________
void AliFemtoModelCorrFctnSource::Write()
{
  fHistROut->Write();
  fHistRSide->Write();
  fHistRLong->Write();
  fHistRStar->Write();
  fHistdNdR->Write();
  
  AliFemtoModelCorrFctn::Write();
}
//_______________________
AliFemtoModelCorrFctnSource* AliFemtoModelCorrFctnSource::Clone()
{
  AliFemtoModelCorrFctnSource *tCopy = new AliFemtoModelCorrFctnSource(*this);
  
  return tCopy;
}

