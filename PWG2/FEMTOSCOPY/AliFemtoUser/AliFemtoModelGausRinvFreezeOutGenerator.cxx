////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausRinvFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelGausRinvFreezeOutGenerator, 1)
#endif

#include "math.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoLorentzVector.h"

//_______________________
AliFemtoModelGausRinvFreezeOutGenerator::AliFemtoModelGausRinvFreezeOutGenerator() :
  fSizeInv(0)
{
  // Default constructor
  fRandom = new TRandom2();
}

//_______________________
AliFemtoModelGausRinvFreezeOutGenerator::AliFemtoModelGausRinvFreezeOutGenerator(const AliFemtoModelGausRinvFreezeOutGenerator &aModel):
  AliFemtoModelFreezeOutGenerator(),
  fSizeInv(0)
{
  // Copy constructor
  fRandom = new TRandom2();
  SetSizeInv(aModel.GetSizeInv());
}
//_______________________
AliFemtoModelGausRinvFreezeOutGenerator::~AliFemtoModelGausRinvFreezeOutGenerator()
{
  if (fRandom) delete fRandom;
}
//_______________________
void AliFemtoModelGausRinvFreezeOutGenerator::GenerateFreezeOut(AliFemtoPair *aPair)
{
  // Generate two particle emission points with respect
  // to their pair momentum 
  // The source is the 3D Gaussian ellipsoid in the LCMS frame
  AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();

  if ((!inf1) || (!inf2)) { cout << "Hidden info not created! "  << endl; exit(kFALSE); }

  // Calculate sum momenta
  Double_t tPx = inf1->GetTrueMomentum()->x() + inf2->GetTrueMomentum()->x();
  Double_t tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
  Double_t tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();
  Double_t tM1 = inf1->GetMass();
  Double_t tM2 = inf2->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + inf1->GetTrueMomentum()->mag2());
  Double_t tE2 = sqrt(tM2*tM2 + inf2->GetTrueMomentum()->mag2());
  Double_t tEs = tE1 + tE2;

  Double_t tPt = sqrt(tPx*tPx + tPy*tPy);
  Double_t tMt = sqrt(tEs*tEs - tPz*tPz);

  // Generate positions in PRF from a Gaussian
  Double_t tROutS =  fRandom->Gaus(0,fSizeInv); // reuse of long
  Double_t tRSideS = fRandom->Gaus(0,fSizeInv);
  Double_t tRLongS = fRandom->Gaus(0,fSizeInv);
  Double_t tRTimeS = 0;
      
  Double_t tBetat  = tPt/tMt;
  Double_t tGammat = 1.0/sqrt(1.0-tBetat*tBetat);

  Double_t tBetaz  = tPz/tEs;
  Double_t tGammaz = 1.0/sqrt(1.0-tBetaz*tBetaz);

  Double_t tROut = tGammat * (tROutS + tBetat * tRTimeS);
  Double_t tDtL  = tGammat * (tRTimeS + tBetat * tROutS);
  Double_t tRSide = tRSideS;

  Double_t tRLong = tGammaz * (tRLongS + tBetaz * tDtL);
  Double_t tDt    = tGammaz * (tDtL + tBetaz * tRLongS);
	  
  tPx /= tPt;
  tPy /= tPt;
	  
  Double_t tXout  = tROut*tPx-tRSide*tPy;
  Double_t tXside = tROut*tPy+tRSide*tPx;
  Double_t tXlong = tRLong;
  Double_t tXtime = tDt;
  
  if (!(inf1->GetEmissionPoint())) {
    AliFemtoLorentzVector tPos(0,0,0,0);
    inf1->SetEmissionPoint(&tPos);
  }
  else
    inf1->SetEmissionPoint(0,0,0,0);
  if (!(inf2->GetEmissionPoint())) {
    AliFemtoLorentzVector tPos(tXout,tXside,tXlong,tXtime);
    inf2->SetEmissionPoint(&tPos);
  }
  else
    inf2->SetEmissionPoint(tXout, tXside, tXlong, tXtime);
}

//_______________________
void AliFemtoModelGausRinvFreezeOutGenerator::SetSizeInv(Double_t aSizeInv)
{
  fSizeInv = aSizeInv;
}
//_______________________
Double_t AliFemtoModelGausRinvFreezeOutGenerator::GetSizeInv() const
{
  return fSizeInv;
}
//_______________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelGausRinvFreezeOutGenerator::Clone() const
{ 
  return GetGenerator(); 
}
//_______________________
inline AliFemtoModelFreezeOutGenerator* AliFemtoModelGausRinvFreezeOutGenerator::GetGenerator() const 
{ 
  AliFemtoModelFreezeOutGenerator* tModel = new AliFemtoModelGausRinvFreezeOutGenerator(*this); 
  return tModel; 
}
