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
  , fLabCand(-1)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fMotherPt(-1.)
  , fGenerateDecayList(kTRUE)
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAHists(0x0)
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
  , fLabCand(-1)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fMotherPt(-1.)
  , fGenerateDecayList(kTRUE)
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAHists(0x0)
{
  // default constructor with name
  fDecayStrList = new TList();
  fDecayStrList->SetOwner();
  fQAHists = new TList();
  fQAHists->SetOwner();
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
  , fLabCand(c.fLabCand)
  , fNprongs(c.fNprongs)
  , fNprongsInAcc(c.fNprongsInAcc)
  , fMotherPt(c.fMotherPt)
  , fGenerateDecayList(c.fGenerateDecayList)
  , fDecayProngs(c.fDecayProngs)
  , fDecayStrList(c.fDecayStrList) // FIXME: TList copy contructor not implemented
  , fQAHists(c.fQAHists)
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

  fQAHists->Add(new TH1F("hRapidityDist", "All particles;y;Counts (a.u.)",800,-4.0,4.0)); // 0
  fQAHists->Add(new TH1F("hRapidityDistStable", "All stable particles;y;Counts (a.u.)",800,-4.0,4.0)); // 1
  fQAHists->Add(new TH1F("hDecayLengthB", ";Decay Length B Meson (cm?);counts (a.u.)",200,0.,2.0)); // 2
  fQAHists->Add(new TH1F("hDecayLengthBxy", ";Decay Length B Meson in XY direction (cm?);counts (a.u.)",200,0.,2.0)); // 3
  fQAHists->Add(new TH1F("hDecayLengthD", ";Decay Length D Meson (cm?);counts (a.u.)",200,0.,2.0)); // 4
  fQAHists->Add(new TH1F("hDecayLengthDxy", ";Decay Length D Meson in XY direction (cm?);counts (a.u.)",200,0.,2.0)); // 5
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
      Double_t y =((AliAODMCParticle*)fMCarray->UncheckedAt(i))->Y();
      ((TH1F*)fQAHists->At(0))->Fill(y); // all particles
      if (IsStable(i)) ((TH1F*)fQAHists->At(1))->Fill(y); // all stable particles
    }
  }
  else {
    fPtMCGenStep->Fill(pt,0.,weight);
  }
  fMCarray=0x0;
  return;
}

void AliHFsubtractBFDcuts::FillSparses(AliAODRecoDecayHF2Prong *dzerocand,Int_t isSelected,Double_t pt,Double_t massD0,Double_t massD0bar,Double_t weight,TClonesArray* mcArray){
  fMCarray=mcArray;
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
    GetCandidateLabel(dzerocand);
    if (!AnalyseDecay(kFALSE, kFALSE)) {
      AliDebug(3, "Error during the decay type determination!");
    }
    Double_t pointMC[5]={pt,normalDecayLengXY,cptangXY,(Double_t)fNprongsInAcc,fMotherPt};
    fCutsMC->Fill(pointMC, weight);
  }
  fMCarray=0x0;
  return;
}

void AliHFsubtractBFDcuts::GetCandidateLabel(AliAODRecoDecayHF2Prong *dzerocand) {
  Int_t labDau0=((AliAODTrack*)dzerocand->GetDaughter(0))->GetLabel();
  AliAODMCParticle* firstDau=(AliAODMCParticle*)fMCarray->UncheckedAt(TMath::Abs(labDau0));
  fLabCand = firstDau->GetMother();
}

Bool_t AliHFsubtractBFDcuts::AnalyseDecay(Bool_t generateString, Bool_t MConly) {
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
    return AnalyseDecay(generateString, MConly);
  }
  if ((pdgMother%1000)/100!=5 && (pdgMother%10000)/1000!=5) {
    AliDebug(3, "Found strange decay, expected the mother to be a beauty hadron!");
    fNprongs=0;
    fNprongsInAcc=0;
    return kFALSE;
  }
  CountProngs(fLabMother, fLabCand, generateString); // count the prongs

  if (generateString) {
    // Store the decay
    std::sort(fDecayProngs.begin(), fDecayProngs.end());
    TString decayStr = "";
    for (ULong64_t i=0; i<fDecayProngs.size(); ++i) {
      decayStr = (i==0) ? Form("%d", fDecayProngs[i]) : Form("%s_%d", decayStr.Data(), fDecayProngs[i]);
    }
    TObjString* str = new TObjString(decayStr);
    if (!fDecayStrList->FindObject(str)) fDecayStrList->Add(str); // only allow unique entries
    fDecayProngs.clear();
  }

  if (MConly) {
      Double_t decayLengthB;   // Decay length B meson
      Double_t decayLengthBxy; // Decay length B meson (xy-plane)
      Double_t decayLengthD;   // Decay length D meson
      Double_t decayLengthDxy; // Decay length D meson (xy-plane)
      Double_t originB[3] = {0.,0.,0.};
      if(!mother->XvYvZv(originB)) AliDebug(3, "Couldn't determine MC origin of the beauty hadron");
      Double_t originD[3] = {0.,0.,0.};
      if(!cand->XvYvZv(originD)) AliDebug(3, "Couldn't determine MC origin of the charmed hadron");
      decayLengthBxy = TMath::Sqrt((originB[0]-originD[0])*(originB[0]-originD[0])+
                                   (originB[1]-originD[1])*(originB[1]-originD[1]));
      decayLengthB   = TMath::Sqrt(decayLengthBxy*decayLengthBxy+(originB[2]-originD[2])*(originB[2]-originD[2]));
      ((TH1F*)fQAHists->At(2))->Fill(decayLengthB);
      ((TH1F*)fQAHists->At(3))->Fill(decayLengthBxy);
  }

  fMotherPt=mother->Pt();
  return kTRUE;
}

void AliHFsubtractBFDcuts::CountProngs(Int_t labCurrMother, Int_t labCurrExcl,
                                       Bool_t generateString) {
  for (Int_t i=0; i<fMCarray->GetEntriesFast(); ++i) {
    if (i!=labCurrExcl) {
      if (((AliAODMCParticle*)fMCarray->UncheckedAt(i))->GetMother()==labCurrMother) {
        if (!fResolveResonances || IsStable(i)) {
          if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(i))->GetPdgCode());
            ++fNprongs;
          if (!fCheckAcceptance || IsInAcceptance(i)) {
            ++fNprongsInAcc;
          }
        }
        else CountProngs(i, -1, generateString);
      }
    }
    else {
      ++fNprongs; // candidate is only counted as a single prong
      ++fNprongsInAcc;
      if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(i))->GetPdgCode());
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
