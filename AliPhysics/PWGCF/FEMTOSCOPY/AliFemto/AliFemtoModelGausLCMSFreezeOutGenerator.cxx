////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelGausLCMSFreezeOutGenerator, 1);
  /// \endcond
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
 // AliFemtoModelGlobalHiddenInfo

  //ml
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  if ((!inf1) || (!inf2)) { cout << "Hidden info not created! "  << endl; exit(kFALSE); }

  //std::cout<<" we are in Freeze-out Generator inf1 inf2  "<<inf1<<"  "<<inf2<<std::endl;
 //std::cout<<" inf1 GetMass "<<((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid()<<std::endl;
 //std::cout<<" true mom    " <<((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->x()<<std::endl;

  Double_t tPx = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->x()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->x();
  Double_t tPy = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->y()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->y();
  Double_t tPz = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->z()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->z();


 // Double_t tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
 // Double_t tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();

 //std::cout<<" tPx tPy tPz"<<tPx<<" "<<tPy<<" "<<tPz<<std::endl;
 if (!(tPx==0 && tPy==0 && tPz==0 )) {

  Double_t tM1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetMass();
  Double_t tM2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->Mag2());
  Double_t tEs = tE1 + tE2;

//std::cout<<" tM1 tM2 tE1 tE2"<<tM1<<" "<<tM2<<" "<<tE2<<std::endl;


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

//std::cout<<" tXout tXside before hidden infor "<<tXout<<" "<<tXside<<std::endl;

  inf1->SetEmissionPoint(0,0,0,0);
  inf2->SetEmissionPoint(tXout, tXside, tXlong, tXtime);

//std::cout<<" after all tXout tXside "<<tXout<<" "<<tXside<<std::endl;

 }
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
