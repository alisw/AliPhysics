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
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoTrack.h"

//_______________________
AliFemtoModelGausRinvFreezeOutGenerator::AliFemtoModelGausRinvFreezeOutGenerator() :
  fSizeInv(0),
  fSelectPrimary(false)
{
  // Default constructor
  fRandom = new TRandom2();
}

//_______________________
AliFemtoModelGausRinvFreezeOutGenerator::AliFemtoModelGausRinvFreezeOutGenerator(const AliFemtoModelGausRinvFreezeOutGenerator &aModel):
  AliFemtoModelFreezeOutGenerator(),
  fSizeInv(0),
  fSelectPrimary(false)
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
AliFemtoModelGausRinvFreezeOutGenerator& AliFemtoModelGausRinvFreezeOutGenerator::operator=(const AliFemtoModelGausRinvFreezeOutGenerator &aModel)
{
  if (this != &aModel) {
    fRandom = new TRandom2();
    SetSizeInv(aModel.GetSizeInv());
  }

  return *this;

}
//_______________________
void AliFemtoModelGausRinvFreezeOutGenerator::GenerateFreezeOut(AliFemtoPair *aPair)
{
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  if ((!inf1) || (!inf2)) { cout << "Hidden info not created! "  << endl; exit(kFALSE); }

  if (fSelectPrimary) {
    const AliFemtoTrack *infg1 = dynamic_cast<const AliFemtoTrack *>(aPair->Track1()->Track());
    const AliFemtoTrack *infg2 = dynamic_cast<const AliFemtoTrack *>(aPair->Track2()->Track());
    
    if ((infg1) && (infg2)) {
      // assume the emission point is in [cm] and try to judge if
      // both particles are primary
      Double_t dist1 = infg1->GetGlobalEmissionPoint()->Perp();
      Double_t dist2 = infg2->GetGlobalEmissionPoint()->Perp();

      if ((dist1 > 0.05) && (dist2 > 0.05)) {
	// At least one particle is non primary
	if (!(inf1->GetEmissionPoint())) {
	  AliFemtoLorentzVector tPos(-1000,1000,-500,0);
	  inf1->SetEmissionPoint(&tPos);
	}
	else
	  inf1->SetEmissionPoint(-1000,1000,-500,0);
	if (!(inf2->GetEmissionPoint())) {
	  AliFemtoLorentzVector tPos(fRandom->Gaus(0,1000.0),fRandom->Gaus(0,1000),fRandom->Gaus(0,1000),0.0);
	  inf2->SetEmissionPoint(&tPos);
	}
	else
	  inf2->SetEmissionPoint(fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), 0.0);
	
	return;
      }
    }
  }

  // Generate two particle emission points with respect
  // to their pair momentum 
  // The source is the 3D Gaussian ellipsoid in the LCMS frame

  // Calculate sum momenta
  Double_t tPx = inf1->GetTrueMomentum()->x() + inf2->GetTrueMomentum()->x();
  Double_t tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
  Double_t tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();
  Double_t tM1 = inf1->GetMass();
  Double_t tM2 = inf2->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + inf1->GetTrueMomentum()->Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + inf2->GetTrueMomentum()->Mag2());
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
//_______________________
void AliFemtoModelGausRinvFreezeOutGenerator::SetSelectPrimaryFromHidden(bool aUse)
{
  fSelectPrimary = aUse;
}
Bool_t AliFemtoModelGausRinvFreezeOutGenerator::GetSelectPrimaryFromHidden()
{
  return fSelectPrimary;
}

