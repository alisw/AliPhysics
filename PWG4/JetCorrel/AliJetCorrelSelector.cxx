/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

//__________________________________________
// Class for user selections. Object is created by ConfigJetCorrel.C macro
// and passed to AddTaskJetCorrel.C running macro.
//-- Author: Paul Constantin

#include "AliJetCorrelSelector.h"

using namespace std;
using namespace JetCorrelHD;

ClassImp(AliJetCorrelSelector)

AliJetCorrelSelector::AliJetCorrelSelector() : 
  fGenQA(kFALSE), fNumCorrel(0), nEvtTriggs(0), fPoolDepth(0), fCorrelType(NULL), fEvtTriggs(NULL),
  minTriggPt(0), maxTriggPt(0), bwTriggPt(0), minAssocPt(0), maxAssocPt(0), bwAssocPt(0),
  fITSRefit(kFALSE), fTPCRefit(kFALSE), fTRDRefit(kFALSE), fRejectKinkChild(kFALSE), fMaxEta(0),
  fMaxNsigmaVtx(0), fMaxTrkVtx(0), fMaxITSChi2(0), fMaxTPCChi2(0), fMinNClusITS(0), fMinNClusTPC(0), trkMinProx(0) {
  // (default) constructor
  fNumBins[centr] = 0; fBinning[centr] = NULL;
  fNumBins[zvert] = 0; fBinning[zvert] = NULL;
}

AliJetCorrelSelector::~AliJetCorrelSelector(){
  // destructor
  if(fCorrelType) delete [] fCorrelType;
  fNumCorrel = 0;
  if(fEvtTriggs) delete [] fEvtTriggs;
  nEvtTriggs = 0;
  if(fBinning[centr]) delete [] fBinning[centr];
  fNumBins[centr] = 0;
  if(fBinning[zvert]) delete [] fBinning[zvert];
  fNumBins[zvert] = 0;
}

void AliJetCorrelSelector::SetCorrelTypes(UInt_t s, UInt_t * const v){
  // fills the array of correlation types
  if(s<1){std::cerr<<"AliJetCorrelSelector::SetCorrelTypes - empty array"<<std::endl; exit(-1);}
  if(s>kMAXNUMCORREL){std::cerr<<"AliJetCorrelSelector: correlation array too big!"<<std::endl; exit(-1);}
  fNumCorrel = s;
  fCorrelType = new UInt_t[fNumCorrel];
  for(UInt_t k=0; k<fNumCorrel; k++){
    if(v[k]>99){
      std::cerr<<"AliJetCorrelSelector::SetCorrelTypes - read error? val["<<k<<"]="<<v[k]<<std::endl;
      exit(-1);
    }
    else fCorrelType[k] = v[k];
  } 
}

void AliJetCorrelSelector::SetTriggers(UInt_t s, TString * const v){
  // fills the array of event triggers
  if(s<1){std::cerr<<"AliJetCorrelSelector::SetTriggers - empty array"<<std::endl; exit(-1);}
  if(s>9){std::cerr<<"AliJetCorrelSelector: event trigger array too big!"<<std::endl; exit(-1);}
  nEvtTriggs = s;
  fEvtTriggs = new TString[nEvtTriggs];
  for(UInt_t k=0; k<nEvtTriggs; k++){
    if(!v[k].IsAscii()){
      std::cerr<<"AliJetCorrelSelector::SetTriggers - read error? val["<<k<<"]="<<v[k]<<std::endl;
      exit(-1);
    }
    else fEvtTriggs[k] = ToUpper(v[k]);
  } 
}


void AliJetCorrelSelector::SetBinningCentr(UInt_t s, Float_t * const v){
  // fills array of centrality bins
  if(s<1){std::cerr<<"AliJetCorrelSelector::SetBinningCentr - empty array"<<std::endl; exit(-1);}
  if(s>kMAXCENTBIN){std::cerr<<"AliJetCorrelSelector: centrality array too big!"<<std::endl; exit(-1);}
  fNumBins[centr] = s;
  fBinning[centr] = new Float_t[fNumBins[centr]]; 
  for(UInt_t k=0; k<fNumBins[centr]; k++){
    if(TMath::Abs(v[k])>999.){
      std::cerr<<"AliJetCorrelSelector::SetBinningCentr - read error? val["<<k<<"]="<<v[k]<<std::endl;
      exit(-1);
    }
    else fBinning[centr][k] = v[k];
  }
}

void AliJetCorrelSelector::SetBinningZvert(UInt_t s, Float_t * const v){
  // fills array of vertex bins
  if(s<1){std::cerr<<"AliJetCorrelSelector::SetBinningZvert - empty array"<<std::endl; exit(-1);}
  if(s>kMAXVERTBIN){std::cerr<<"AliJetCorrelSelector: vertex array too big!"<<std::endl; exit(-1);}
  fNumBins[zvert] = s;
  fBinning[zvert] = new Float_t[fNumBins[zvert]]; 
  for(UInt_t k=0; k<fNumBins[zvert]; k++){
    if(TMath::Abs(v[k])>999.){
      std::cerr<<"AliJetCorrelSelector::SetBinningZvert - read error? val["<<k<<"]="<<v[k]<<std::endl;
      exit(-1);
    }
    else fBinning[zvert][k] = v[k];
  }
}

void AliJetCorrelSelector::SetBinningTrigg(Float_t min, Float_t max, Float_t bw){
  // sets trigger Pt binning
  minTriggPt = min;
  maxTriggPt = max;
  bwTriggPt  = bw;
}

void AliJetCorrelSelector::SetBinningAssoc(Float_t min, Float_t max, Float_t bw){
  // sets associated Pt binning
  minAssocPt = min;
  maxAssocPt = max;
  bwAssocPt  = bw;
}

void AliJetCorrelSelector::Show(){
  // print out all user selections
  std::cout<<"Generic selections: "<<std::endl<<" PoolDepth="<<fPoolDepth;
  std::cout<<std::endl<<" Correlation Types: ";
  for(UInt_t k=0; k<fNumCorrel; k++) std::cout<<fCorrelType[k]<<" ";
  std::cout<<std::endl<<" Event Triggers: ";
  for(UInt_t k=0; k<nEvtTriggs; k++) std::cout<<fEvtTriggs[k]<<" ";
  std::cout<<std::endl<<" Centrality/Multiplicity binning: ";
  for(UInt_t k=0; k<fNumBins[centr]; k++) std::cout<<fBinning[centr][k]<<" ";
  std::cout<<std::endl<<" Vertex binning: ";
  for(UInt_t k=0; k<fNumBins[zvert]; k++) std::cout<<fBinning[zvert][k]<<" ";
  std::cout<<std::endl<<" Trigg binning:"<<minTriggPt<<"->"<<maxTriggPt<<"/"<<bwTriggPt;
  std::cout<<std::endl<<" Assoc binning:"<<minAssocPt<<"->"<<maxAssocPt<<"/"<<bwAssocPt;
  std::cout<<std::endl<<"Track selections: "<<std::endl
	   <<" MaxEta="<<fMaxEta<<" MaxTrkVtx="<<fMaxTrkVtx<<" MaxNsigmaVtx="<<fMaxNsigmaVtx<<std::endl
	   <<" MaxITSChi2="<<fMaxITSChi2<<" MaxTPCChi2="<<fMaxTPCChi2<<std::endl
	   <<" MinNClusITS="<<fMinNClusITS<<" MinNClusTPC="<<fMinNClusTPC<<std::endl
	   <<" ITSRefit="<<fITSRefit<<" TPCRefit="<<fTPCRefit<<" TRDRefit="<<fTRDRefit<<std::endl
	   <<" RejectKinkChild="<<fRejectKinkChild<<" minTrackPairTPCDist="<<trkMinProx<<std::endl;
}

Float_t AliJetCorrelSelector::BinBorder(BinType_t cType, UInt_t k){
  // returns bin margins
  if(k<=NoOfBins(cType)) return fBinning[cType][k];
  else {std::cerr<<"BinBorder Error: bin of type "<<cType<<" outside range "<<k<<std::endl; exit(0);}
}

Int_t AliJetCorrelSelector::GetBin(BinType_t cType, Float_t val){
  // returns bin number
  Int_t iBin=-1; UInt_t nBins=NoOfBins(cType);
  for(UInt_t i=0; i<nBins; i++)
    if(BinBorder(cType,i)<=val && val<BinBorder(cType,i+1)) iBin=i;
  return iBin;
}

//////////////////////////////////////////////////////////
// Cutting Methods
/////////////////////////////////////////////////////////

Bool_t AliJetCorrelSelector::SelectedEvtTrigger(AliVEvent *fEVT){
  // matches the event trigger classes with the user trigger classes
  if(fEVT->InheritsFrom("AliESDEvent")){
    const AliESDEvent *esd = (AliESDEvent*)fEVT;
    TString trigClass = esd->GetFiredTriggerClasses();
    if(nEvtTriggs==1 && fEvtTriggs[0].Contains("ALL")) return kTRUE;
    for(UInt_t k=0; k<nEvtTriggs; k++)
      if(trigClass.Contains(fEvtTriggs[k])) return kTRUE;
    return kFALSE;
  } else {std::cerr<<"AliJetCorrelSelector::SelectedEvtTrigger ERROR: not an ESD event!"<<std::endl; exit(0);}
}

Bool_t AliJetCorrelSelector::CloseTrackPair(Float_t dist){
  // applies two-track cut (dist at TPC entrance); it is possible that single-track cuts,
  // like fraction of shared TPC clusters, will avoid inclusion of split tracks...
  if(dist>trkMinProx) return kFALSE;
  return kTRUE;
}

Bool_t AliJetCorrelSelector::LowQualityTrack(AliESDtrack* track){
  // selects low quality tracks
  if(track->Eta()>fMaxEta) return kTRUE;
  UInt_t status = track->GetStatus();
  if(fITSRefit && !(status & AliESDtrack::kITSrefit)) return kTRUE;
  if(fTPCRefit && !(status & AliESDtrack::kTPCrefit)) return kTRUE;

//   UInt_t nClusITS = track->GetITSclusters(0);
//   if(nClusITS<fMinNClusITS) return kTRUE;
//   Float_t chi2ITS=-1.;
//   if(nClusITS!=0) chi2ITS = track->GetITSchi2()/Float_t(nClusITS);
//   if(chi2ITS<0 || chi2ITS>fMaxITSChi2) return kTRUE;

  UInt_t nClusTPC = track->GetTPCclusters(0); // or track->GetTPCNcls() ?
  if(nClusTPC<fMinNClusTPC) return kTRUE;
  Float_t chi2TPC=-1.;
  if(nClusTPC!=0) chi2TPC = track->GetTPCchi2()/Float_t(nClusTPC);
  if(chi2TPC<0 || chi2TPC>fMaxTPCChi2) return kTRUE;

  if(fRejectKinkChild && track->GetKinkIndex(0)>0) return kTRUE;
//  Float_t sigTrkVtx = GetSigmaToVertex(track);
//  if(sigTrkVtx<0 || sigTrkVtx>fMaxNsigmaVtx) return kTRUE;
  // instead of track-vertex DCA sigma cut, apply value-cut:
  Float_t b[2], bCov[3];
  track->GetImpactParameters(b,bCov);
  if((b[0]*b[0]+b[1]*b[1])>(fMaxTrkVtx*fMaxTrkVtx)) return kTRUE;

  return kFALSE;
}

Bool_t AliJetCorrelSelector::PassPID(AliESDtrack* track, PartType_t PartType){
  // checks if a track has the required ID
  Bool_t hasReqPID = kFALSE;
  Stat_t fPid;
  Stat_t fWeight;
  GetPID(track, fPid, fWeight);
  
  switch(PartType){
    case hadron:
      // if(fPID!=0) hasReqPID = kTRUE;
      hasReqPID = kTRUE;
      break;
    case electron:
      // if(fTRDRefit && !(status & AliESDtrack::kTRDrefit)) hasReqPID = kFALSE;
      // if(fPID!=0) hasReqPID = kFALSE;
      hasReqPID = kTRUE;
      break;
    case pion:
      // if(fPID!=2) hasReqPID = kFALSE;
      hasReqPID = kTRUE;
      break;
    case kaon:
      // if(fPID!=3) hasReqPID = kFALSE;
      hasReqPID = kTRUE;
      break;
    case proton:
      // if(fPID!=4) hasReqPID = kFALSE;
      hasReqPID = kTRUE;
      break;
    default:
      std::cerr<<"AliJetCorrelSelector::PassPID() - ERROR: wrong type!"<<std::endl;
      exit(-1);
  }
  return hasReqPID;
}

Float_t AliJetCorrelSelector::GetSigmaToVertex(AliESDtrack* track){
  // Calculates the number of sigma to the vertex; from ANALYSIS/AliESDtrackCuts
  Float_t b[2], bRes[2], bCov[3];
  track->GetImpactParameters(b,bCov);
  if(bCov[0]<=0 || bCov[2]<=0) return -1.;
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  if(TMath::Exp(-d*d/2)<1e-10) return 1000;

  Float_t nSigma = TMath::ErfInverse(1-TMath::Exp(-d*d/2))*TMath::Sqrt(2);
  return nSigma;
}

void AliJetCorrelSelector::GetPID(AliESDtrack* track, Stat_t& fpid, Stat_t& fweight){
  // Finds most probable particle: 0=Electron, 1=Muon, 2=Pion, 3=Kaon, 4=Proton
  fpid = -1;
  fweight = -1;
  
  Double_t wpart[5], wpartbayes[5];
  track->GetESDpid(wpart); // probability of the different particle types
  Double_t c[5]={1., 1., 1., 1., 1.}; // Tentative particle type "concentrations"
  //  Double_t c[5]={0.01, 0.01, 0.85, 0.10, 0.05};
  
  //Bayes formula
  Double_t rcc = 0.;
  for(Int_t i=0; i<5; i++) {rcc += c[i]*wpart[i];}
  if(TMath::Abs(rcc)<1e-10) return;
  for(Int_t i=0; i<5; i++) {wpartbayes[i] = c[i]*wpart[i]/rcc;}
  
  Int_t ipid=-1;
  Float_t max=0.;
  for(Int_t i=0; i<5; i++) {
    if(wpartbayes[i]>max) {
      ipid = i; 
      max = wpartbayes[i];
    }
  }

  fpid = ipid;
  fweight = max;
}
