#include "AliSpherocityEstimator.h"
ClassImp(AliSpherocityEstimator)
AliSpherocityEstimator::AliSpherocityEstimator():
  TNamed(),
  fSphTrackCuts(0),
  fMinMulti(10),
  fTrackMulti(0),
  fMinPt(0.15)
{
  fSphTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fSphTrackCuts->SetRequireTPCRefit(kTRUE);
  fSphTrackCuts->SetRequireITSRefit(kTRUE);
  fSphTrackCuts->SetEtaRange(-0.8,0.8);
};
AliSpherocityEstimator::AliSpherocityEstimator(const AliSpherocityEstimator &source):
  TNamed(source.GetName(), source.GetTitle()),
  fSphTrackCuts(source.fSphTrackCuts),
  fMinMulti(source.fMinMulti),
  fTrackMulti(source.fTrackMulti),
  fMinPt(source.fMinPt)
{
};
AliSpherocityEstimator &AliSpherocityEstimator::operator=(const AliSpherocityEstimator &source) {
  if(&source == this) return *this;
  TNamed::operator=(source);
  *fSphTrackCuts=*source.fSphTrackCuts;
  fMinMulti=source.fMinMulti;
  fTrackMulti=source.fTrackMulti;
  fMinPt=source.fMinPt;
  return *this;
};
AliSpherocityEstimator::~AliSpherocityEstimator() {
  delete fSphTrackCuts;
};
void AliSpherocityEstimator::SetESDtrackCuts(AliESDtrackCuts *newcuts) {
  if(fSphTrackCuts)
    delete fSphTrackCuts;
  fSphTrackCuts = newcuts;
};
Double_t AliSpherocityEstimator::GetSpherocity(AliESDEvent *inevent) {
  if(!inevent) return -1; //If a problem with event, return -1
  if(inevent->GetNumberOfTracks()<fMinMulti) return -2;//if not enough tracks (1st check), return -2
  fTrackMulti=0;
  AliESDtrack *track; //Create pointer here, not for each track
  Double_t sumpt=0;
  std::vector<Double_t> nx;
  std::vector<Double_t> ny;
  std::vector<Double_t> px;
  std::vector<Double_t> py;
  for(Int_t itrk=0;itrk<inevent->GetNumberOfTracks();++itrk) {
    track = inevent->GetTrack(itrk);
    if(!track) continue;
    if(!fSphTrackCuts->AcceptTrack(track)) continue;
    Double_t pt = track->Pt();
    if(pt<fMinPt) continue;
    sumpt+=pt;
    Double_t phi = track->Phi();
    nx.push_back(TMath::Cos(phi));
    ny.push_back(TMath::Sin(phi));
    px.push_back(pt*nx[fTrackMulti]);
    py.push_back(pt*ny[fTrackMulti]);
    ++fTrackMulti;
  }; 
  if(fTrackMulti<fMinMulti) {
    px.clear();
    py.clear();
    nx.clear();
    ny.clear();
    return -2; //if not enought tracks 
  };
  //Calculating spherocity now
  Double_t retval=2; //Return value should be the minimal value in range 0..1, so 2 is a good starting point.
  for(Int_t i=0; i<fTrackMulti; i++) {
    Double_t num=0;
    for(Int_t j=0;j<fTrackMulti;j++)
      num+=TMath::Abs(ny[i]*px[j] -nx[i]*py[j]);
    Double_t sFull = TMath::Power((num/sumpt),2);
    if(sFull < retval)
      retval = sFull;
  };
  px.clear();
  py.clear();
  nx.clear();
  ny.clear();
  if(retval>1) return -3; //If sph>1, something went wrong
  retval=retval*TMath::Pi()*TMath::Pi()/4; //normalization
  return retval;
}; 

