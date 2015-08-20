/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Class for storing and handling D0 meson candidates properties         //
//  for estimating the feed-down fraction using several sets of cuts      //
//     Andrea Rossi <andrea.rossi@cern.ch>                                //
//     Felix Reidt  <felix.reidt@cern.ch>                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#include "AliHFsubtractBFDcuts.h"

#include <vector>

#include "TNamed.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TObjString.h"
#include "TMath.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"

ClassImp(AliHFsubtractBFDcuts);

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts()
  : TNamed()
  , fIsMC(kTRUE)
  , fCheckAcceptance(kTRUE)
  , fResolveResonances(kTRUE)
  , fPtMCGenStep(0x0)
  , fCutsData(0x0)
  , fCutsMC(0x0)
  , fMCarray(0x0)
  , fAODtracks(0x0)
  , fLabCand(-1)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fMotherPt(-1.)
  , fGenerateDecayList(kTRUE)
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAhists(0x0)
{
  // default constructor
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const char* name, const char* title)
  : TNamed(name,title)
  , fIsMC(kTRUE)
  , fCheckAcceptance(kTRUE)
  , fResolveResonances(kTRUE)
  , fPtMCGenStep(0x0)
  , fCutsData(0x0)
  , fCutsMC(0x0)
  , fMCarray(0x0)
  , fAODtracks(0x0)
  , fLabCand(-1)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fMotherPt(-1.)
  , fGenerateDecayList(kTRUE)
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAhists(0x0)
{
  // default constructor with name
  fDecayStrList = new TList();
  fDecayStrList->SetOwner();
  fQAhists = new TList();
  fQAhists->SetOwner();
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const AliHFsubtractBFDcuts& c)
  : TNamed(c.GetName(),c.GetTitle())
  , fIsMC(c.fIsMC)
  , fCheckAcceptance(c.fCheckAcceptance)
  , fResolveResonances(c.fResolveResonances)
  , fPtMCGenStep(c.fPtMCGenStep)
  , fCutsData(c.fCutsData)
  , fCutsMC(c.fCutsMC)
  , fMCarray(c.fMCarray)
  , fAODtracks(c.fAODtracks)
  , fLabCand(c.fLabCand)
  , fNprongs(c.fNprongs)
  , fNprongsInAcc(c.fNprongsInAcc)
  , fMotherPt(c.fMotherPt)
  , fGenerateDecayList(c.fGenerateDecayList)
  , fDecayProngs(c.fDecayProngs)
  , fDecayStrList(c.fDecayStrList) // FIXME: TList copy contructor not implemented
  , fQAhists(c.fQAhists)
{
  // copy constructor
}

AliHFsubtractBFDcuts AliHFsubtractBFDcuts::operator=(const AliHFsubtractBFDcuts& c)
{
  // assignment operator
  return c;
}

AliHFsubtractBFDcuts::~AliHFsubtractBFDcuts() {
  AliDebug(3, "TODO: implement me!");
}

void AliHFsubtractBFDcuts::InitHistos(){
  Int_t dimAxes[4]={500  ,24 ,30 ,400   };// mass, pt, normLXY, cosPointXY
  Double_t  min[4]={  1.7, 0., 0.,  0.96};
  Double_t  max[4]={  2.2,24.,30.,  1.  };
  fCutsData=new THnSparseF("fCutsDataFD","fCutsDataFD",4,dimAxes,min,max);
  fCutsData->GetAxis(0)->SetName("mass");
  fCutsData->GetAxis(0)->SetTitle("Mass (K,#pi) (GeV/#it{c^{2}})");
  fCutsData->GetAxis(1)->SetName("pt");
  fCutsData->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fCutsData->GetAxis(2)->SetName("NormDecLengthXY");
  fCutsData->GetAxis(2)->SetTitle("Normalized XY decay length");
  fCutsData->GetAxis(3)->SetName("CosPointXY");
  fCutsData->GetAxis(3)->SetTitle("Cos#theta_{point}^{XY}");

  Int_t dimAxesMC[5]={24 , 30,400   ,20 ,24 };// pt, normLXY, cosPointXY, #prongs, mother pt
  Double_t  minMC[5]={ 0., 0.,  0.96, 0., 0.};
  Double_t  maxMC[5]={24.,30.,  1.  ,20.,24.};
  fCutsMC=new THnSparseF("fCutsMCFD","fCutsMCFD",5,dimAxesMC,minMC,maxMC);
  fCutsMC->GetAxis(0)->SetName("pt");
  fCutsMC->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fCutsMC->GetAxis(1)->SetName("NormDecLengthXY");
  fCutsMC->GetAxis(1)->SetTitle("Normalized XY decay length");
  fCutsMC->GetAxis(2)->SetName("CosPointXY");
  fCutsMC->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  fCutsMC->GetAxis(3)->SetName("nProngs");
  fCutsMC->GetAxis(3)->SetTitle("Number of Decay Prongs");
  fCutsMC->GetAxis(4)->SetName("Mother_pt");
  fCutsMC->GetAxis(4)->SetTitle("Mother #it{p}_{T} (GeV/#it{c})");

  fPtMCGenStep=new TH3F("fPtMCGenStep","fPtMCGenStep;#it{p}_{T};Number of prongs;Mother #it{p}_{T}",24,0.,24.,20,0.,20.,24,0.,24.);

  fQAhists->Add(new TH1F("hRapidityDist", "All Particles;y;Counts (a.u.)",800,-4.0,4.0)); // 0
  fQAhists->Add(new TH1F("hRapidityDistStable", "All Stable Particles;y;Counts (a.u.)",800,-4.0,4.0)); // 1
  fQAhists->Add(new TH1F("hDecayLengthB", ";Decay Length B Meson (cm?);counts (a.u.)",200,0.,2.0)); // 2
  fQAhists->Add(new TH1F("hDecayLengthBxy", ";Decay Length B Meson in XY direction (cm?);counts (a.u.)",200,0.,2.0)); // 3
  fQAhists->Add(new TH2F("hDCABhypothesis", "All Particles;Distance of Closest Approach (cm?);Impact parameter of D0 x track (cm^{2});counts (a.u.)",1000,0.,1.0,10000,0.,0.01)); // 4
  fQAhists->Add(new TH2F("hDCABprongs", "Only Real B Decay Prongs;Distance of Closest Approach (cm?);Impact parameter of D0 x track (cm^{2});counts (a.u.)",1000,0.,1.0,10000,0.,0.01)); // 5
  return;
}

void AliHFsubtractBFDcuts::FillGenStep(AliAODMCParticle *dzeropart,Double_t pt/*=-1.*/,Double_t weight/*=1.*/,TClonesArray* mcArray/*=0x0*/){
  fMCarray=mcArray;
  if (pt<0) {
    pt=dzeropart->Pt();
  }
  if (fIsMC && fMCarray) {
    fNprongs=0;
    fNprongsInAcc=0;
    fDecayChain=kFALSE; // TODO: use this value
    fMotherPt=pt;
    fLabCand=dzeropart->GetLabel();
    if (!AnalyseDecay(fGenerateDecayList, kTRUE)) {
      AliDebug(3, "Error during the decay type determination!");
    }
    fPtMCGenStep->Fill(pt,(Double_t)fNprongsInAcc,fMotherPt,weight);
    // y distribution of all particles
    for (Int_t i=0; i<fMCarray->GetEntriesFast(); ++i) {
      Double_t y=((AliAODMCParticle*)fMCarray->UncheckedAt(i))->Y();
      ((TH1F*)fQAhists->At(0))->Fill(y); // all particles
      if (IsStable(i)) ((TH1F*)fQAhists->At(1))->Fill(y); // all stable particles
    }
  }
  else {
    fPtMCGenStep->Fill(pt,0.,weight);
  }
  fMCarray=0x0;
  return;
}

void AliHFsubtractBFDcuts::FillSparses(AliAODRecoDecayHF2Prong *dzerocand,Int_t isSelected,Double_t pt,Double_t massD0,Double_t massD0bar,Double_t weight,TClonesArray* mcArray,AliAODEvent* aodEvent){
  fMCarray=mcArray;
  if (aodEvent) {
    fAODtracks=aodEvent->GetTracks();
    fBkG      =aodEvent->GetMagneticField();
    fPriVtx   =aodEvent->GetPrimaryVertex();
  }
  if(isSelected<=0||isSelected>3){
    Printf("isSelected = %d", isSelected);
    return;
  }
  if(massD0<0){
    dzerocand->InvMassD0(massD0,massD0bar);
    Printf("mass D0 =%f, massD0bar = %f ", massD0, massD0bar);
  }
  if(pt<0){
    pt=dzerocand->Pt();
  }
  Double_t normalDecayLengXY=dzerocand->NormalizedDecayLengthXY();
  Double_t cptangXY=dzerocand->CosPointingAngleXY();
  Double_t pointData[4]={massD0,pt,normalDecayLengXY,cptangXY};

  if(isSelected==1||isSelected==3){
    fCutsData->Fill(pointData,weight);
  }
  if(isSelected==2||isSelected==3){
    pointData[0]=massD0bar;
    fCutsData->Fill(pointData,weight);
  }

  if(fIsMC && fMCarray) {
    fNprongs=0;
    fNprongsInAcc=0;
    fDecayChain=kFALSE; // TODO: use this value
    fMotherPt=pt;
    fD0Cand=dzerocand;
    if (fD0Cand) fD0CandParam=new AliNeutralTrackParam(fD0Cand);
    if (!GetCandidateLabel() || !AnalyseDecay(kFALSE, kFALSE)) {
      AliDebug(3, "Error during the decay type determination!");
    }
    Double_t pointMC[5]={pt,normalDecayLengXY,cptangXY,(Double_t)fNprongsInAcc,fMotherPt};
    fCutsMC->Fill(pointMC, weight);
  }
  fMCarray=0x0;
  fAODtracks=0x0;
  fBkG=0.;
  fD0Cand=0x0;
  return;
}

Bool_t AliHFsubtractBFDcuts::GetCandidateLabel(){
  if (fD0Cand) {
    Int_t labDau0=((AliAODTrack*)fD0Cand->GetDaughter(0))->GetLabel();
    AliAODMCParticle* firstDau=(AliAODMCParticle*)fMCarray->UncheckedAt(TMath::Abs(labDau0));
    fLabCand = firstDau->GetMother();
    return kTRUE;
  }
  else {
    AliDebug(3, "Could not obtain the label of the candidate");
    return kFALSE;
  }
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::AnalyseDecay(Bool_t generateString, Bool_t mcOnly) {
  AliAODMCParticle* cand=(AliAODMCParticle*)fMCarray->UncheckedAt(fLabCand);
  fLabMother = cand->GetMother();
  AliAODMCParticle* mother = (AliAODMCParticle*)fMCarray->UncheckedAt(fLabMother);
  Int_t pdgMother = mother->GetPdgCode();
  if (pdgMother<0) pdgMother*=-1; // treat particles and anti-particles the same way
  if (pdgMother==4 || pdgMother==2212 || pdgMother==2112) { // prompt production
    fNprongs=1;
    fNprongsInAcc=1;
    return kTRUE;
  }
  if ((pdgMother%1000)/100==4 || (pdgMother%10000)/1000==4) {
    // chained decay of charmed hadrons, using recursion to resolve it
    fDecayChain=kTRUE;
    fLabCand=fLabMother;
    return AnalyseDecay(generateString, mcOnly);
  }
  if ((pdgMother%1000)/100!=5 && (pdgMother%10000)/1000!=5) {
    AliDebug(3, "Found strange decay, expected the mother to be a beauty hadron!");
    fNprongs=0;
    fNprongsInAcc=0;
    return kFALSE;
  }
  CountProngs(fLabMother, fLabCand, generateString, mcOnly); // count the prongs

  if (generateString) {
    // Store the decay
    std::sort(fDecayProngs.begin(), fDecayProngs.end());
    TString decayStr = "";
    for (ULong64_t i=0; i<fDecayProngs.size(); ++i) {
      decayStr = (i==0) ? Form("%d", fDecayProngs[i]) : Form("%s_%d", decayStr.Data(), fDecayProngs[i]);
    }
    //decayStr = Form("%s__%d_%d_%d", decayStr.Data(), fNprongs, fNprongsInAcc, fDecayProngs.size());
    TObjString* str = new TObjString(decayStr);
    if (!fDecayStrList->FindObject(str)) fDecayStrList->Add(str); // only allow unique entries
    fDecayProngs.clear();
  }

  if (mcOnly) {
    Double_t decayLengthB;   // Decay length B meson
    Double_t decayLengthBxy; // Decay length B meson (xy-plane)
    Double_t originB[3] = {0.,0.,0.};
    if(!mother->XvYvZv(originB)) AliDebug(3, "Couldn't determine MC origin of the beauty hadron");
    Double_t originD[3] = {0.,0.,0.};
    if(!cand->XvYvZv(originD)) AliDebug(3, "Couldn't determine MC origin of the charmed hadron");
    decayLengthBxy = TMath::Sqrt((originB[0]-originD[0])*(originB[0]-originD[0])+
                                 (originB[1]-originD[1])*(originB[1]-originD[1]));
    decayLengthB   = TMath::Sqrt(decayLengthBxy*decayLengthBxy+(originB[2]-originD[2])*(originB[2]-originD[2]));
    ((TH1F*)fQAhists->At(2))->Fill(decayLengthB);
    ((TH1F*)fQAhists->At(3))->Fill(decayLengthBxy);
  }
  else {
    for (Int_t iAODtrack=0; iAODtrack<fAODtracks->GetEntriesFast(); ++iAODtrack) {
      CheckBhypothesis(iAODtrack, kFALSE);
    }
  }

  fMotherPt=mother->Pt();
  return kTRUE;
}

void AliHFsubtractBFDcuts::CountProngs(Int_t labCurrMother, Int_t labCurrExcl,
                                       Bool_t generateString, Bool_t mcOnly) {
  for (Int_t iMCParticle=0; iMCParticle<fMCarray->GetEntriesFast(); ++iMCParticle) {
    if (iMCParticle!=labCurrExcl) {
      if (((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetMother()==labCurrMother) {
        if (!fResolveResonances || IsStable(iMCParticle)) {
          if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetPdgCode());
          ++fNprongs;
          if (!mcOnly) {
            for (Int_t iAODtrack=0; iAODtrack<fAODtracks->GetEntriesFast(); ++iAODtrack) {
              AliAODTrack* aodTrack=(AliAODTrack*)fAODtracks->UncheckedAt(iAODtrack);
              if (aodTrack->GetLabel()==iMCParticle) CheckBhypothesis(iAODtrack, kTRUE);
            }
          }
          if (!fCheckAcceptance || IsInAcceptance(iMCParticle)) {
            ++fNprongsInAcc;
          }
        }
        else CountProngs(iMCParticle, -1, generateString, mcOnly);
      }
    }
    else {
      ++fNprongs; // candidate is only counted as a single prong
      ++fNprongsInAcc;
      if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetPdgCode());
    }
  }
}

Bool_t AliHFsubtractBFDcuts::IsStable(Int_t labProng) const {
  const Int_t stablePartPdgs[] = { 11, 13, 211, 321, 2212, 12, 14, 22, 111, 130 };
  const Int_t nStablePartPdgs  = sizeof(stablePartPdgs)/sizeof(Int_t);
  AliAODMCParticle* prong = (AliAODMCParticle*)fMCarray->UncheckedAt(labProng);
  Int_t pdgProng = prong->GetPdgCode();
  if (pdgProng<0) pdgProng*=-1; // treat particles and anti-particles the same way
  for (Int_t iPdg=0; iPdg<nStablePartPdgs; ++iPdg) {
    if (stablePartPdgs[iPdg] == pdgProng) return kTRUE;
  }
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::IsInAcceptance(Int_t labProng) const {
  AliDebug(1, "AliHFsubtractBFDcuts::IsInAcceptance(...) hasn't been implemented yet, prong");
  Double_t eta    = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Eta();
  Double_t pt     = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Pt();
  Short_t charge = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Charge();
  if ((pt>0.15) && (eta>-0.9) && (eta<0.9) && (charge!=0)) return kTRUE;
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::CheckBhypothesis(Int_t iAODtrack, Bool_t Bprong) {
  AliExternalTrackParam *t = new AliExternalTrackParam();
  t->CopyFromVTrack((AliVTrack*)fAODtracks->UncheckedAt(iAODtrack));

  TObjArray *tracks = new TObjArray(2);
  tracks->AddAt(t,           0);
  tracks->AddAt(fD0CandParam,1);

  AliAODVertex *bVtx = RecBvtx(tracks);
  if(!bVtx) {
    AliDebug(3, "Couldn't reconstruct B meson vertex!");
    delete t; t=0x0;
    return kFALSE;
  }

  const Double_t maxD = 1.;

  // Propagate candidates to secondary vertex
  Double_t dz[2],cov[3];
  t->PropagateToDCA(bVtx,fBkG,maxD,dz,cov);
  fD0CandParam->PropagateToDCA(bVtx,fBkG,maxD,dz,cov);

  // Impact parameters
  t->PropagateToDCA(fPriVtx,fBkG,maxD,dz,cov);
  Double_t d0Track=dz[0];
  Double_t d0TrackErr=TMath::Sqrt(cov[0]);
  fD0CandParam->PropagateToDCA(fPriVtx,fBkG,maxD,dz,cov);
  Double_t d0D0Cand=dz[0];
  Double_t d0D0CandErr=TMath::Sqrt(cov[0]);

  // distance of closest approach of the D0 and the track
  Double_t xDCAtrack, xDCAD0;
  Double_t dcaB=t->GetDCA(fD0CandParam,fBkG,xDCAtrack,xDCAD0);
  if (!Bprong) ((TH1F*)fQAhists->At(4))->Fill(dcaB,d0D0Cand*d0Track);
  else         ((TH1F*)fQAhists->At(5))->Fill(dcaB,d0D0Cand*d0Track);

  delete tracks; tracks=0x0;
  delete bVtx  ; bVtx  =0x0;
  delete t     ; t     =0x0;
  return kTRUE;
}

AliAODVertex* AliHFsubtractBFDcuts::RecBvtx(TObjArray *tracks) const {
  AliESDVertex* vtxESD=0x0;
  AliVertexerTracks* vertexer = new AliVertexerTracks(fBkG);
  vertexer->SetVtxStart((AliESDVertex*)fPriVtx);
  vtxESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);
  delete vertexer; vertexer=0x0;
  if(!vtxESD || vtxESD->GetNContributors()!=tracks->GetEntriesFast()) {
    AliDebug(3, "Couldn't reconstruct B Meson vertex");
    delete vtxESD; vtxESD=0x0;
    return 0x0;
  }
  Double_t rVtxSq=vtxESD->GetX()*vtxESD->GetX()+vtxESD->GetY()*vtxESD->GetY();
  if(rVtxSq>8.){
    // vertex outside beam pipe, reject candidate to avoid propagation through material
    delete vtxESD; vtxESD=0x0;
    return 0x0;
  }
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vtxESD->GetXYZ(pos);
  vtxESD->GetCovMatrix(cov);
  chi2perNDF=vtxESD->GetChi2toNDF();
  delete vtxESD; vtxESD=0x0;
  return new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,0);
}
