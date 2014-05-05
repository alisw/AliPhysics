//-------------------------------------------------------------------------
//               Implementation of the ITSU track class
//             based on the "cooked covariance" approach
//-------------------------------------------------------------------------
#include "AliITSUTrackCooked.h"
#include "AliCluster.h"
#include "AliESDtrack.h"

ClassImp(AliITSUTrackCooked)

AliITSUTrackCooked::AliITSUTrackCooked(): 
AliKalmanTrack()
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------
    for (Int_t i=0; i<2*AliITSUTrackerCooked::kNLayers; i++) {
        fIndex[i]=0;
    }
}

AliITSUTrackCooked::AliITSUTrackCooked(const AliESDtrack &t):
AliKalmanTrack()
{
    //--------------------------------------------------------------------
    // Constructior from an ESD track
    //--------------------------------------------------------------------
    Set(t.GetX(), t.GetAlpha(), t.GetParameter(), t.GetCovariance());
    SetLabel(t.GetITSLabel());
    SetChi2(t.GetITSchi2());
    SetNumberOfClusters(t.GetITSclusters(fIndex));
}

AliITSUTrackCooked::~AliITSUTrackCooked()
{
  //--------------------------------------------------------------------
  // Virtual destructor
  //--------------------------------------------------------------------
}

void AliITSUTrackCooked::SetClusterIndex(Int_t l, Int_t i)
{
    //--------------------------------------------------------------------
    // Set the cluster index
    //--------------------------------------------------------------------
    Int_t idx = (l<<28) + i;
    Int_t n=GetNumberOfClusters();
    fIndex[n]=idx;
    SetNumberOfClusters(n+1);
}

Double_t AliITSUTrackCooked::GetPredictedChi2(const AliCluster *c) const {
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};
  return AliExternalTrackParam::GetPredictedChi2(p,cov);
}

Bool_t AliITSUTrackCooked::PropagateTo(Double_t xk, Double_t t,Double_t x0rho) {
  //------------------------------------------------------------------
  // This function propagates a track
  // t is the material thicknes in units X/X0
  // x0rho is the material X0*density
  //------------------------------------------------------------------
  Double_t xOverX0,xTimesRho; 
  xOverX0 = t; xTimesRho = t*x0rho;
  if (!CorrectForMeanMaterial(xOverX0,xTimesRho,GetMass(),kTRUE)) return kFALSE;

  Double_t bz=GetBz();
  if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE;
  //Double_t b[3]; GetBxByBz(b);
  //if (!AliExternalTrackParam::PropagateToBxByBz(xk,b)) return kFALSE;

  return kTRUE;
}

Bool_t AliITSUTrackCooked::
Propagate(const AliCluster *c, Double_t t, Double_t x0rho) {
  //------------------------------------------------------------------
  // This function propagates a track to the plane the cluster belong to
  // t is the material thicknes in units X/X0
  // x0rho is the material X0*density
  //------------------------------------------------------------------
  Double_t xOverX0,xTimesRho; 
  xOverX0 = t; xTimesRho = t*x0rho;
  if (!CorrectForMeanMaterial(xOverX0,xTimesRho,GetMass(),kTRUE)) return kFALSE;

  Float_t xk, alpha;
  if (!c->GetXAlphaRefPlane(xk,alpha)) return kFALSE;

  Double_t bz=GetBz();
  if (!AliExternalTrackParam::Propagate(alpha, xk, bz)) return kFALSE;
  //Double_t b[3]; GetBxByBz(b);
  //if (!AliExternalTrackParam::PropagateBxByBz(alpha, xk, b)) return kFALSE;

  return kTRUE;
}

Bool_t AliITSUTrackCooked::Update(const AliCluster *c, Double_t chi2, Int_t idx)
{
  //--------------------------------------------------------------------
  // Update track params
  //--------------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), c->GetSigmaYZ(), c->GetSigmaZ2()};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  fIndex[n]=idx;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return kTRUE;
}
