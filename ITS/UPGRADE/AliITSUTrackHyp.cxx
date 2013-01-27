#include "AliITSUTrackHyp.h"
#include "AliESDtrack.h"
#include "AliCluster.h"
#include "AliITSUAux.h"

ClassImp(AliITSUTrackHyp)



//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(Int_t nlr) 
: fNLayers(nlr)
  ,fITSLabel(0)
  ,fESDTrack(0)
  ,fLayerSeeds(0)
{
  // def. c-tor
  if (fNLayers>0) fLayerSeeds = new TObjArray[fNLayers];
}

//__________________________________________________________________
AliITSUTrackHyp::~AliITSUTrackHyp() 
{
  // d-tor
  delete[] fLayerSeeds;
}

//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(const AliITSUTrackHyp &src)
  : AliKalmanTrack(src)
  , fNLayers(src.fNLayers)
  , fITSLabel(src.fITSLabel)
  , fESDTrack(src.fESDTrack)
  , fLayerSeeds(0)
{
  // copy c-tor
  if (fNLayers>0) {
    fLayerSeeds = new TObjArray[fNLayers];
    for (int ilr=fNLayers;ilr--;) {
      int ns = src.GetNSeeds(ilr);
      for (int isd=0;isd<ns;isd++) {
	AliITSUSeed* sd = src.GetSeed(ilr,isd);
	if (sd->IsKilled()) continue;
	AddSeed(sd,ilr);
      }      
    }
  }
  //
}

//__________________________________________________________________
AliITSUTrackHyp &AliITSUTrackHyp::operator=(const AliITSUTrackHyp &src)
{
  // copy 
  if (this == &src) return *this;
  this->~AliITSUTrackHyp();
  new(this) AliITSUTrackHyp(src);
  return *this;
  //
}

//__________________________________________________________________
void AliITSUTrackHyp::Print(Option_t* ) const
{
  printf("Track Hyp.#%4d. NSeeds:",GetUniqueID());
  for (int i=0;i<fNLayers;i++) printf(" (%d) %3d",i,GetNSeeds(i)); printf("\n");
}

//__________________________________________________________________
AliITSUSeed* AliITSUTrackHyp::GetWinner() const
{
  // Get best candidate
  return fLayerSeeds[0].GetEntriesFast()>0 ? GetSeed(0,0) : 0;
}

//__________________________________________________________________
void AliITSUTrackHyp::DefineWinner(int lr, int id)
{
  // assign best candidate
  AliITSUSeed* winner = GetSeed(lr,id);
  this->AliExternalTrackParam::operator=(*winner);
  SetChi2(winner->GetChi2GloNrm());
  SetNumberOfClusters(winner->GetNLayersHit());
}

//__________________________________________________________________
Double_t AliITSUTrackHyp::GetPredictedChi2(const AliCluster *cl) const
{
  // calculate chi2 to cluster
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  return AliExternalTrackParam::GetPredictedChi2(p,cov);
}

//__________________________________________________________________
Bool_t AliITSUTrackHyp::PropagateTo(Double_t /*xr*/, Double_t /*x0*/, Double_t /*rho*/)
{
  // NA
  AliFatal("Not to be used");
  return 0;
}

//__________________________________________________________________
Bool_t AliITSUTrackHyp::Update(const AliCluster* /*c*/, Double_t /*chi2*/, Int_t /*index*/)
{
  // NA
  AliFatal("Not to be used");
  return kFALSE;
}

//__________________________________________________________________
Bool_t AliITSUTrackHyp::Update(const AliCluster* cl)
{
  // update with cluster
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  double chi2 = AliExternalTrackParam::GetPredictedChi2(p,cov);
  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;
  SetChi2(GetChi2()+chi2);
  return kTRUE;
}

//__________________________________________________________________
Int_t AliITSUTrackHyp::FetchClusterInfo(Int_t *clIDarr) const
{
  // fill cl.id's in the array. The clusters of layer L will be set at slots
  // clID[2L] (and clID[2L+1] if there is an extra cluster).
  for (int i=fNLayers<<1;i--;) clIDarr[i]=-1;
  AliITSUSeed* seed = GetWinner();
  Int_t lr,ncl=0;
  while(seed) {
    int clID = seed->GetLrCluster(lr);
    if (clID>=0) {
      int slotLr = lr<<1;
      clIDarr[ clIDarr[slotLr]<0 ? slotLr : slotLr+1 ] = clID;
      ncl++;
    }
    seed = (AliITSUSeed*)seed->GetParent();
  }
  return ncl;
}
