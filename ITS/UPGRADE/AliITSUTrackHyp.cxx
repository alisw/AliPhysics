#include "AliITSUTrackHyp.h"
#include "AliESDtrack.h"
#include "AliCluster.h"
#include "AliITSUAux.h"
#include <TString.h>

ClassImp(AliITSUTrackHyp)



//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(Int_t nlr) 
: fNLayers(nlr)
  ,fITSLabel(0)
  ,fESDTrack(0)
  ,fWinner(0)
  ,fTPCSeed(0)
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
  delete fTPCSeed;
}

//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(const AliITSUTrackHyp &src)
  : AliKalmanTrack(src)
  , fNLayers(src.fNLayers)
  , fITSLabel(src.fITSLabel)
  , fESDTrack(src.fESDTrack)
  , fWinner(0)
  , fTPCSeed(src.fTPCSeed)
  , fLayerSeeds(0)
{
  // copy c-tor. Note: it is shallow
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
    fWinner = src.fWinner;
  }
  //
}
//__________________________________________________________________
void AliITSUTrackHyp::InitFrom(const AliITSUTrackHyp *src)
{
  // copy initial params
  fNLayers = src->fNLayers;
  fITSLabel = src->fITSLabel;
  fESDTrack = src->fESDTrack;
  fWinner = src->fWinner;
  fTPCSeed = src->fTPCSeed;
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
void AliITSUTrackHyp::Print(Option_t* opt) const
{
  printf("Track Hyp.#%4d. NSeeds:",GetUniqueID());
  TString opts = opt;
  opts.ToLower();
  Bool_t prSeeds = opts.Contains("l");
  for (int i=0;i<fNLayers;i++) {
    printf("Lr (%d) %3d ",i,GetNSeeds(i)); 
    if (prSeeds) {
      printf("\n");
      for (int isd=0;isd<GetNSeeds(i);isd++) ((AliITSUSeed*)GetSeed(i,isd))->Print(opt);
    }
  }
  if (!prSeeds) printf("\n");
}

//__________________________________________________________________
AliITSUSeed* AliITSUTrackHyp::GetWinner() const
{
  // Get best candidate. TODO
  return fWinner;
}

//__________________________________________________________________
AliITSUSeed* AliITSUTrackHyp::DefineWinner(int lr, int id)
{
  // assign best candidate
  if (GetNSeeds(lr)<=id) return 0;
  fWinner = GetSeed(lr,id);
  this->AliExternalTrackParam::operator=(*fWinner);
  SetChi2(fWinner->GetChi2GloNrm());
  SetNumberOfClusters(fWinner->GetNLayersHit());
  return fWinner;
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
Double_t AliITSUTrackHyp::Update(const AliCluster* cl)
{
  // update with cluster
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  double chi2 = AliExternalTrackParam::GetPredictedChi2(p,cov);
  if (!AliExternalTrackParam::Update(p,cov)) return -1;
  SetChi2(GetChi2()+chi2);
  return chi2;
}

//__________________________________________________________________
Int_t AliITSUTrackHyp::FetchClusterInfo(Int_t *clIDarr) const
{
  // fill cl.id's in the array. The clusters of layer L will be set at slots
  // clID[2L] (and clID[2L+1] if there is an extra cluster).
  for (int i=fNLayers<<1;i--;) clIDarr[i]=-1;
  return GetWinner() ? GetWinner()->FetchClusterInfo(clIDarr) : 0;
}

//__________________________________________________________________
Int_t AliITSUTrackHyp::GetNumberOfClusters() const
{
  // This is a temporary (slow) way of accessing number of clusters
  // TODO: add dedicated data members filled by winner
  AliITSUSeed* seed = GetWinner();
  int ncl = 0;
  if (!seed) {
    AliFatal("The winner is not set");
    return ncl;
  }
  return seed->GetNClusters();
  //
}

//__________________________________________________________________
Int_t AliITSUTrackHyp::GetClusterIndex(Int_t ind) const 
{
  // This is a temporary (slow) way of accessing cluster index
  // TODO: add dedicated data members filled by winner
  AliITSUSeed* seed = GetWinner();
  //  int ncl = 0;
  if (!seed) {
    AliFatal("The winner is not set");
    return -1;
  }
  return seed->GetClusterIndex(ind);
  //
}
