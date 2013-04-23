////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelGausLCMSFreezeOutGenerator, 1)
#endif

#include "math.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoLorentzVector.h"
#include <TMath.h>

//_______________________
AliFemtoModelGausLCMSFreezeOutGenerator::AliFemtoModelGausLCMSFreezeOutGenerator() :
  fSizeOut(0), fSizeSide(0), fSizeLong(0)
{
  // Default constructor
  fRandom = new TRandom2();
}

//_______________________
AliFemtoModelGausLCMSFreezeOutGenerator::AliFemtoModelGausLCMSFreezeOutGenerator(const AliFemtoModelGausLCMSFreezeOutGenerator &aModel):
  AliFemtoModelFreezeOutGenerator(aModel),
  fSizeOut(0), fSizeSide(0), fSizeLong(0)
{
  // Copy constructor
  fRandom = new TRandom2();
  SetSizeOut(aModel.GetSizeOut());
  SetSizeSide(aModel.GetSizeSide());
  SetSizeLong(aModel.GetSizeLong());
}
//_______________________
AliFemtoModelGausLCMSFreezeOutGenerator::~AliFemtoModelGausLCMSFreezeOutGenerator()
{
  if (fRandom) delete fRandom;
}
//_______________________
AliFemtoModelGausLCMSFreezeOutGenerator& AliFemtoModelGausLCMSFreezeOutGenerator::operator=(const AliFemtoModelGausLCMSFreezeOutGenerator &aModel)
{
  if (this != &aModel) {
    fRandom = new TRandom2();
    SetSizeOut(aModel.GetSizeOut());
    SetSizeSide(aModel.GetSizeSide());
    SetSizeLong(aModel.GetSizeLong());
  }

  return *this;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGenerator::GenerateFreezeOut(AliFemtoPair *aPair)
{
  // Generate two particle emission points with respect
  // to their pair momentum 
  // The source is the 3D Gaussian ellipsoid in the LCMS frame
  //AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  if ((!inf1) || (!inf2)) { cout << "Hidden info not created! "  << endl; exit(kFALSE); }

  Double_t tPx = inf1->GetTrueMomentum()->x() + inf2->GetTrueMomentum()->x();
  Double_t tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
  Double_t tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();
  Double_t tM1 = inf1->GetMass();
  Double_t tM2 = inf2->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + inf1->GetTrueMomentum()->Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + inf2->GetTrueMomentum()->Mag2());
  Double_t tEs = tE1 + tE2;

  Double_t tPt = sqrt(tPx*tPx + tPy*tPy);

  Double_t tRout = fRandom->Gaus(0.0, fSizeOut);
  Double_t tRside = fRandom->Gaus(0.0, fSizeSide);
  Double_t tRlong = fRandom->Gaus(0.0, fSizeLong);
  
  Double_t tXout = (tPx * tRout + tPy * tRside)/tPt;
  Double_t tXside = (tPy * tRout - tPx * tRside)/tPt;

  Double_t tBetaz = tPz/tEs;
  Double_t tGammaz = 1.0/TMath::Sqrt(1-tBetaz*tBetaz);
  
  Double_t tXlong = tGammaz * (tRlong + tBetaz * 0);
  Double_t tXtime = tGammaz * (0 + tBetaz * tRlong);
  
  if (!(inf1->GetEmissionPoint())) {
    AliFemtoLorentzVector *tPos = new AliFemtoLorentzVector(0,0,0,0);
    inf1->SetEmissionPoint(tPos);
    delete tPos;
  }
  else
    inf1->SetEmissionPoint(0,0,0,0);
  if (!(inf2->GetEmissionPoint())) {
    AliFemtoLorentzVector *tPos = new AliFemtoLorentzVector(tXout,tXside,tXlong,tXtime);
    inf2->SetEmissionPoint(tPos);
    delete tPos;
  }
  else
    inf2->SetEmissionPoint(tXout, tXside, tXlong, tXtime);
}

//_______________________
void AliFemtoModelGausLCMSFreezeOutGenerator::SetSizeOut(Double_t aSizeOut)
{
  fSizeOut = aSizeOut;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGenerator::SetSizeSide(Double_t aSizeSide)
{
  fSizeSide = aSizeSide;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGenerator::SetSizeLong(Double_t aSizeLong)
{
  fSizeLong = aSizeLong;
}

//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGenerator::GetSizeOut() const
{
  return fSizeOut;
}
//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGenerator::GetSizeSide() const
{
  return fSizeSide;
}
//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGenerator::GetSizeLong() const
{
  return fSizeLong;
}
//_______________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelGausLCMSFreezeOutGenerator::Clone() const
{ 
  return GetGenerator(); 
}
//_______________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelGausLCMSFreezeOutGenerator::GetGenerator() const 
{ 
  AliFemtoModelFreezeOutGenerator* tModel = new AliFemtoModelGausLCMSFreezeOutGenerator(*this); return tModel; 
}
