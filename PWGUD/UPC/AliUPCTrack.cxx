
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    contains parameters of the central track relevant for UPC analysis
//_____________________________________________________________________________

#include "TLorentzVector.h"

#include "AliUPCTrack.h"

ClassImp(AliUPCTrack)

//_____________________________________________________________________________
AliUPCTrack::AliUPCTrack()
 :TObject(),
  fCharge(0), fMaskMan(0), fFilterMap(0), fChi2perNDF(0),
  fTPCmomentum(0), fTPCsignal(0), fTPCncls(0),
  fTPCrows(0), fTPCnclsF(0), fTPCnclsS(0),
  fITSchi2perNDF(0), fITSClusterMap(0)
{

  // Default constructor

  for(Int_t i=0; i<2; i++) {fDZ[i] = 0.; fdzSPD[i] = 0.; fdzIP[i] = 0.;}
  for(Int_t i=0; i<3; i++) {fP[i] = 0.; fCov[i] = 0.; fCovSPD[i] = 0.; fCovIP[i] = 0.;}
}

//_____________________________________________________________________________
void AliUPCTrack::Clear(Option_t *)
{
  // clear track flags

  fMaskMan = 0;
  fFilterMap = 0;
}

//_____________________________________________________________________________
void AliUPCTrack::SetImpactParameters(Float_t *p, Float_t *cov)
{
  // set impact parameters in XY and Z to default primary vertex

  for(Int_t i=0; i<2; i++) {
    fDZ[i] = p[i];
    fCov[i] = cov[i];
  }
  fCov[2] = cov[2];
}

//_____________________________________________________________________________
void AliUPCTrack::SetImpactParametersSPD(Float_t *p, Float_t *cov)
{
  // set SPD impact parameters

  for(Int_t i=0; i<2; i++) {
    fdzSPD[i] = p[i];
    fCovSPD[i] = cov[i];
  }
  fCovSPD[2] = cov[2];
}

//_____________________________________________________________________________
void AliUPCTrack::SetImpactParametersIP(Float_t *p, Float_t *cov)
{
  // set impact parameters to nominal interaction point

  for(Int_t i=0; i<2; i++) {
    fdzIP[i] = p[i];
    fCovIP[i] = cov[i];
  }
  fCovIP[2] = cov[2];
}

//_____________________________________________________________________________
void AliUPCTrack::GetMomentum(TLorentzVector *v, Double_t mass) const
{
  // get track 4-momentum
  v->SetXYZM(fP[0],fP[1],fP[2],mass);
}

Int_t AliUPCTrack::GetITSNcls(void) const
{
  // number of points in ITS

  Int_t ncls = 0;

  for(UChar_t i=0; i<6; i++) if( fITSClusterMap & (1 << i) ) ncls++;

  return ncls;
}

//_____________________________________________________________________________
void AliUPCTrack::GetImpactParameters(Double_t *p, Double_t *cov) const
{
  // get impact parameters in XY and Z to default primary vertex

  for(Int_t i=0; i<2; i++) {
    p[i] = (Double_t) fDZ[i];
    cov[i] = (Double_t) fCov[i];
  }
  cov[2] = (Double_t) fCov[2];
}

//_____________________________________________________________________________
void AliUPCTrack::GetImpactParametersSPD(Double_t *p, Double_t *cov) const
{
  // get SPD impact parameters

  for(Int_t i=0; i<2; i++) {
    p[i] = (Double_t) fdzSPD[i];
    cov[i] = (Double_t) fCovSPD[i];
  }
  cov[2] = (Double_t) fCovSPD[2];
}

//_____________________________________________________________________________
void AliUPCTrack::GetImpactParametersIP(Double_t *p, Double_t *cov) const
{
  // get impact parameters to nominal interaction point

  for(Int_t i=0; i<2; i++) {
    p[i] = (Double_t) fdzIP[i];
    cov[i] = (Double_t) fCovIP[i];
  }
  cov[2] = (Double_t) fCovIP[2];
}







































