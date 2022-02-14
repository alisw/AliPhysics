/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <iostream>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TKDTree.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"

#include "AliEmcalJetTaggerTaskFast.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast)
/// \endcond

namespace PWGJE {

  namespace EMCALJetTasks{

  AliEmcalJetTaggerTaskFast::AliEmcalJetTaggerTaskFast() :
      AliAnalysisTaskEmcalJet("AliEmcalJetTaggerTaskFast", kTRUE),
      fJetTaggingType(kTag),
      fJetTaggingMethod(kGeo),
      fContainerBase(0),
      fContainerTag(1),
      fSpecPartContTag(-1),
      fMinFractionShared(0),
      fUseSumw2(0),
      fMatchingDone(0),
      fTypeAcc(kLimitBaseTagEtaPhi),
      fMaxDist(0.3),
      fInit(kFALSE),
      fh3PtJet1VsDeltaEtaDeltaPhi(nullptr),
      fh2PtJet1VsDeltaR(nullptr),
      fh2PtJet2VsFraction(nullptr),
      fh2PtJet1VsLeadPtAllSel(nullptr),
      fh2PtJet1VsLeadPtTagged(nullptr),
      fh2PtJet1VsPtJet2(nullptr),
      fh2PtJet2VsRelPt(nullptr),
      fh3PtJetDEtaDPhiConst(nullptr),
      fh3PtJetAreaDRConst(nullptr),
      fNAccJets(nullptr)
#ifdef JETTAGGERFAST_TEST
      , fIndexErrorRateBase(nullptr)
      , fIndexErrorRateTag(nullptr)
      , fContainerErrorRateBase(nullptr)
      , fContainerErrorRateTag(nullptr)
#endif
  {
    SetMakeGeneralHistograms(kTRUE);
  }

  AliEmcalJetTaggerTaskFast::AliEmcalJetTaggerTaskFast(const char *name) :
      AliAnalysisTaskEmcalJet(name, kTRUE),
      fJetTaggingType(kTag),
      fJetTaggingMethod(kGeo),
      fContainerBase(0),
      fContainerTag(1),
      fSpecPartContTag(-1),
      fMinFractionShared(0),
      fUseSumw2(0),
      fMatchingDone(0),
      fTypeAcc(kLimitBaseTagEtaPhi),
      fMaxDist(0.3),
      fInit(kFALSE),
      fh3PtJet1VsDeltaEtaDeltaPhi(nullptr),
      fh2PtJet1VsDeltaR(nullptr),
      fh2PtJet2VsFraction(nullptr),
      fh2PtJet1VsLeadPtAllSel(nullptr),
      fh2PtJet1VsLeadPtTagged(nullptr),
      fh2PtJet1VsPtJet2(nullptr),
      fh2PtJet2VsRelPt(nullptr),
      fh3PtJetDEtaDPhiConst(nullptr),
      fh3PtJetAreaDRConst(nullptr),
      fNAccJets(nullptr)
#ifdef JETTAGGERFAST_TEST
      , fIndexErrorRateBase(nullptr)
      , fIndexErrorRateTag(nullptr)
      , fContainerErrorRateBase(nullptr)
      , fContainerErrorRateTag(nullptr)
#endif
  {
    SetMakeGeneralHistograms(kTRUE);
  }

  void AliEmcalJetTaggerTaskFast::UserCreateOutputObjects() {
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    // Notify the user to be careful.
    AliErrorStream() << "This task isn't yet validated. Please use the standard AliAnalysisTaskEmcalJetTagger.\n";

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    const Int_t nBinsPt          = 40;
    const Int_t nBinsDPhi        = 72;
    const Int_t nBinsDEta        = 100;
    const Int_t nBinsDR          = 50;
    const Int_t nBinsFraction    = 101;

    const Double_t minPt       = -50.;
    const Double_t maxPt       = 150.;
    const Double_t minDPhi     = -0.5;
    const Double_t maxDPhi     =  0.5;
    const Double_t minDEta     = -0.5;
    const Double_t maxDEta     =  0.5;
    const Double_t minDR       =  0.;
    const Double_t maxDR       =  0.5;
    const Double_t minFraction =  -0.005;
    const Double_t maxFraction =  1.005;

    // Prepare histograms
    fh3PtJet1VsDeltaEtaDeltaPhi  = new TH3*[fNcentBins];
    fh2PtJet1VsDeltaR            = new TH2*[fNcentBins];
    fh2PtJet2VsFraction          = new TH2*[fNcentBins];
    fh2PtJet1VsLeadPtAllSel      = new TH2*[fNcentBins];
    fh2PtJet1VsLeadPtTagged      = new TH2*[fNcentBins];
    fh2PtJet1VsPtJet2            = new TH2*[fNcentBins];
    fh2PtJet2VsRelPt             = new TH2*[fNcentBins];

    for (Int_t i = 0; i < fNcentBins; i++) {
      fh3PtJet1VsDeltaEtaDeltaPhi[i] = 0;
      fh2PtJet1VsDeltaR[i]           = 0;
      fh2PtJet2VsFraction[i]         = 0;
      fh2PtJet1VsLeadPtAllSel[i]     = 0;
      fh2PtJet1VsLeadPtTagged[i]     = 0;
      fh2PtJet1VsPtJet2[i]           = 0;
      fh2PtJet2VsRelPt[i]            = 0;
    }

    TString histName = "";
    TString histTitle = "";
    for (Int_t i = 0; i < fNcentBins; i++) {

      histName = TString::Format("fh3PtJet1VsDeltaEtaDeltaPhi_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{#Delta#eta};#it{#Delta#varphi}",histName.Data());
      fh3PtJet1VsDeltaEtaDeltaPhi[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsDEta,minDEta,maxDEta,nBinsDPhi,minDPhi,maxDPhi);
      fOutput->Add(fh3PtJet1VsDeltaEtaDeltaPhi[i]);

      histName = TString::Format("fh2PtJet1VsDeltaR_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{#Delta R}",histName.Data());
      fh2PtJet1VsDeltaR[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsDR,minDR,maxDR);
      fOutput->Add(fh2PtJet1VsDeltaR[i]);

      histName = TString::Format("fh2PtJet2VsFraction_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet2};#it{f}_{shared}",histName.Data());
      fh2PtJet2VsFraction[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsFraction,minFraction,maxFraction);
      fOutput->Add(fh2PtJet2VsFraction[i]);

      histName = TString::Format("fh2PtJet1VsLeadPtAllSel_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
      fh2PtJet1VsLeadPtAllSel[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,20,0.,20.);
      fOutput->Add(fh2PtJet1VsLeadPtAllSel[i]);

      histName = TString::Format("fh2PtJet1VsLeadPtTagged_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
      fh2PtJet1VsLeadPtTagged[i] =  new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,20,0.,20.);
      fOutput->Add(fh2PtJet1VsLeadPtTagged[i]);

      histName = TString::Format("fh2PtJet1VsPtJet2_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{p}_{T,jet2}",histName.Data());
      fh2PtJet1VsPtJet2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
      fOutput->Add(fh2PtJet1VsPtJet2[i]);

      histName = TString::Format("fh2PtJet2VsRelPt_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet2};(#it{p}_{T,jet2}-#it{p}_{T,jet1})/#it{p}_{T,jet1}",histName.Data());
      fh2PtJet2VsRelPt[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,241,-2.41,2.41);
      fOutput->Add(fh2PtJet2VsRelPt[i]);
    }

    fh3PtJetDEtaDPhiConst = new TH3F("fh3PtJetDEtaDPhiConst","fh3PtJetDEtaDPhiConst;pT;#Delta #eta;#Delta #varphi",nBinsPt,minPt,maxPt,nBinsDEta,-1.,1.,nBinsDPhi,-1.,1.);
    fOutput->Add(fh3PtJetDEtaDPhiConst);

    fh3PtJetAreaDRConst = new TH3F("fh3PtJetAreaDRConst","fh3PtJetAreaDRConst;pT;A;#Delta R",nBinsPt,minPt,maxPt,50,0.,1.,50,0.,1.);
    fOutput->Add(fh3PtJetAreaDRConst);

    fNAccJets = new TH1F("fNAccJets","fNAccJets;N/ev",10,-0.5, 9.5);
    fOutput->Add(fNAccJets);

#ifdef JETTAGGERFAST_TEST
    fIndexErrorRateBase = new TH1F("indexErrorsBase", "Index errors nearest neighbor base jets", 1, 0.5, 1.5);
    fIndexErrorRateTag = new TH1F("indexErrorsTag", "Index errors nearest neighbors tag jets", 1, 0.5, 1.5);
    fContainerErrorRateBase = new TH1F("containerErrorsBase", "Matching errors container - kdtree base jets", 1, 0.5, 1.5);
    fContainerErrorRateTag = new TH1F("containerErrorsTag", "Matching errors container - kdtree tag jets", 1, 0.5, 1.5);
    fOutput->Add(fIndexErrorRateBase);
    fOutput->Add(fIndexErrorRateTag);
    fOutput->Add(fContainerErrorRateBase);
    fOutput->Add(fContainerErrorRateTag);
#endif

    if(fUseSumw2) {
      // =========== Switch on Sumw2 for all histos ===========
      for(auto it : *fOutput){
        TH1 *h1 = dynamic_cast<TH1*>(it);
        if (h1){
          h1->Sumw2();
          continue;
        }
        THnSparse *hn = dynamic_cast<THnSparse*>(it);
        if(hn)hn->Sumw2();
      }
    }

    TH1::AddDirectory(oldStatus);

    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  }

  void AliEmcalJetTaggerTaskFast::Init(){

    if(fInit) return;

    AliJetContainer *cont1 = GetJetContainer(fContainerBase);
    AliJetContainer *cont2 = GetJetContainer(fContainerTag);
    if(!cont1 || !cont2) AliError("Missing jet container");

    // when full azimuth, don't do anything
    Double_t phiMin1 = cont1->GetJetPhiMin(), phiMin2 = cont2->GetJetPhiMin();
    Bool_t isZeroTwoPi1 = kFALSE;
    //check only one side of phi, since the upper bound is not well defined
    if(phiMin1 > -1.e-6 && phiMin1 < 1.e-6) isZeroTwoPi1 = kTRUE;
    Bool_t isZeroTwoPi2 = kFALSE;
    if(phiMin2 > -1.e-6 && phiMin2 < 1.e-6) isZeroTwoPi2 = kTRUE;

    switch(fTypeAcc){
    case kNoLimit: break;
    case kLimitTagEta:
      cont2->SetJetEtaLimits(cont2->GetJetEtaMin()-0.1,cont2->GetJetEtaMax()+0.1);
      break;
    case kLimitTagEtaPhi:
      cont2->SetJetEtaLimits(cont2->GetJetEtaMin()-0.1,cont2->GetJetEtaMax()+0.1);
      if(!isZeroTwoPi2) cont2->SetJetPhiLimits(cont2->GetJetPhiMin()-0.1,cont2->GetJetPhiMax()+0.1);
      break;
    case kLimitBaseTagEtaPhi:
      cont1->SetJetEtaLimits(cont1->GetJetEtaMin()-0.1,cont1->GetJetEtaMax()+0.1);
      if(!isZeroTwoPi1) cont1->SetJetPhiLimits(cont1->GetJetPhiMin()-0.1,cont1->GetJetPhiMax()+0.1);
      cont2->SetJetEtaLimits(cont2->GetJetEtaMin()-0.1,cont2->GetJetEtaMax()+0.1);
      if(!isZeroTwoPi2) cont2->SetJetPhiLimits(cont2->GetJetPhiMin()-0.1,cont2->GetJetPhiMax()+0.1);
    };
    fInit = kTRUE;
    return;
  }

  Bool_t AliEmcalJetTaggerTaskFast::Run()
  {
    Init();
    AliJetContainer *contBase = GetJetContainer(fContainerBase),
                    *contTag = GetJetContainer(fContainerTag);

    ResetTagging(*contBase);
    ResetTagging(*contTag);

    fMatchingDone = MatchJetsGeo(*contBase, *contTag, fMaxDist);

    return kTRUE;
  }

  Bool_t AliEmcalJetTaggerTaskFast::FillHistograms()
  {
    // Fill histograms.

    AliEmcalJet *jet1 = NULL;
    AliJetContainer *jetCont = GetJetContainer(fContainerBase);
    if(!jetCont) return kFALSE;
    jetCont->ResetCurrentID();
    Int_t count = 0;
    while((jet1 = jetCont->GetNextAcceptJet())) {
      count++;
      Double_t ptJet1 =  jet1->Pt() - jetCont->GetRhoVal()*jet1->Area();
      fh2PtJet1VsLeadPtAllSel[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());

      //fill histo with angle between jet axis and constituents
      for(Int_t icc=0; icc<jet1->GetNumberOfTracks(); icc++) {
        AliVParticle *vp = static_cast<AliVParticle*>(jet1->Track(icc));
        if(!vp) continue;
        Double_t dEta = jet1->Eta()-vp->Eta();
        Double_t dPhi = jet1->Phi()-vp->Phi();
        if(dPhi<TMath::Pi()) dPhi+=TMath::TwoPi();
        if(dPhi>TMath::Pi()) dPhi-=TMath::TwoPi();
        fh3PtJetDEtaDPhiConst->Fill(ptJet1,dEta,dPhi);

        Double_t dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
        fh3PtJetAreaDRConst->Fill(ptJet1,jet1->Area(),dR);
      }

      if(jet1->GetTagStatus()<1 && fJetTaggingType==kTag)
        continue;

      AliEmcalJet *jet2 = NULL;
      if(fJetTaggingType==kTag)     jet2 = jet1->GetTaggedJet();
      if(fJetTaggingType==kClosest) jet2 = jet1->ClosestJet();
      if(!jet2) continue;

      Double_t ptJet2 =  jet2->Pt() - GetRhoVal(fContainerTag)*jet2->Area();

      Double_t fraction = -2;
      if(fSpecPartContTag > -1) fraction = jetCont->GetFractionSharedPt(jet1, GetParticleContainer(fSpecPartContTag));
      else fraction = jetCont->GetFractionSharedPt(jet1);

      fh2PtJet2VsFraction[fCentBin]->Fill(ptJet2,fraction);
      AliDebug(5, Form("Fraction = %f, minimum = %f", fraction, fMinFractionShared));
      //if(fJetTaggingType==kClosest) Printf("Fraction = %f, minimum = %f", fraction, fMinFractionShared);
      if(fraction<fMinFractionShared && fJetTaggingType==kClosest)
        continue;
      fh2PtJet1VsLeadPtTagged[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());
      fh2PtJet1VsPtJet2[fCentBin]->Fill(ptJet1,ptJet2);
      if(ptJet2>0.) fh2PtJet2VsRelPt[fCentBin]->Fill(ptJet2,(ptJet1-ptJet2)/ptJet2);

      Double_t dPhi = GetDeltaPhi(jet1->Phi(),jet2->Phi());
      if(dPhi>TMath::Pi())
        dPhi -= TMath::TwoPi();
      if(dPhi<(-1.*TMath::Pi()))
        dPhi += TMath::TwoPi();

      fh3PtJet1VsDeltaEtaDeltaPhi[fCentBin]->Fill(ptJet1,jet1->Eta()-jet2->Eta(),dPhi);
      fh2PtJet1VsDeltaR[fCentBin]->Fill(ptJet1,jet1->DeltaR(jet2));
    }
    fNAccJets->Fill(count);
    return kTRUE;
  }

  void AliEmcalJetTaggerTaskFast::ResetTagging(const AliJetContainer &c) const {
    for(auto j : c.all()){
      switch(fJetTaggingType){
      case kClosest:
        j->ResetMatching();
        break;
      case kTag:
        j->SetTaggedJet(nullptr);
        j->SetTagStatus(-1);
      };
    }
  }

  bool AliEmcalJetTaggerTaskFast::MatchJetsGeo(AliJetContainer &contBase, AliJetContainer &contTag, Float_t maxDist) const {
    const Int_t kNacceptedBase = contBase.GetNAcceptedJets(),
                kNacceptedTag = contTag.GetNAcceptedJets();
    if(!(kNacceptedBase && kNacceptedTag)) return false;

    // Build kd-trees
    TArrayD etaBase(kNacceptedBase), phiBase(kNacceptedBase),
            etaTag(kNacceptedTag), phiTag(kNacceptedTag);
    std::vector<AliEmcalJet *> jetsBase(kNacceptedBase), jetsTag(kNacceptedTag); // the storages are needed later for applying the tagging, in order to avoid multiple occurrence of jet selection

    int countBase(0), countTag(0);
    for(auto jb : contBase.accepted()) {
      etaBase[countBase] = jb->Eta();
      phiBase[countBase] = jb->Phi();
      jetsBase[countBase] = jb;
      countBase++;
    }
    for(auto jt : contTag.accepted()) {
      etaTag[countTag] = jt->Eta();
      phiTag[countTag] = jt->Phi();
      jetsTag[countTag] = jt;
      countTag++;
    }
    TKDTreeID treeBase(etaBase.GetSize(), 2, 1), treeTag(etaTag.GetSize(), 2, 1);
    treeBase.SetData(0, etaBase.GetArray());
    treeBase.SetData(1, phiBase.GetArray());
    treeBase.Build();
    treeTag.SetData(0, etaTag.GetArray());
    treeTag.SetData(1, phiTag.GetArray());
    treeTag.Build();

    TArrayI faMatchIndexTag(kNacceptedBase), faMatchIndexBase(kNacceptedTag);
    faMatchIndexBase.Reset(-1);
    faMatchIndexTag.Reset(-1);

    // find the closest distance to the full jet
    countBase = 0;
    for(auto j : contBase.accepted()) {
      Double_t point[2] = {j->Eta(), j->Phi()};
      Int_t index(-1); Double_t distance(-1);
      treeTag.FindNearestNeighbors(point, 1, &index, &distance);
      // test whether indices are matching:
      if(index >= 0 && distance < maxDist){
        AliDebugStream(1) << "Found closest tag jet for " << countBase << " with match index " << index << " and distance " << distance << std::endl;
        faMatchIndexTag[countBase]=index;
      } else {
        AliDebugStream(1) << "Not found closest tag jet for " << countBase << ", distance to closest " << distance << std::endl;
      }

#ifdef JETTAGGERFAST_TEST
      if(index>-1){
        Double_t distanceTest(-1);
        distanceTest = TMath::Sqrt(TMath::Power(etaTag[index] - j->Eta(), 2) +  TMath::Power(phiTag[index] - j->Phi(), 2));
        if(TMath::Abs(distanceTest - distance) > DBL_EPSILON){
          AliDebugStream(1) << "Mismatch in distance from tag jet with index from tree: " << distanceTest << ", distance from tree " << distance << std::endl;
          fIndexErrorRateBase->Fill(1);
        }
      }
#endif

      countBase++;
    }

    // other way around
    countTag = 0;
    for(auto j : contTag.accepted()){
      Double_t point[2] = {j->Eta(), j->Phi()};
      Int_t index(-1); Double_t distance(-1);
      treeBase.FindNearestNeighbors(point, 1, &index, &distance);
      if(index >= 0 && distance < maxDist){
        AliDebugStream(1) << "Found closest base jet for " << countBase << " with match index " << index << " and distance " << distance << std::endl;
        faMatchIndexBase[countTag]=index;
      } else {
        AliDebugStream(1) << "Not found closest tag jet for " << countBase << ", distance to closest " << distance << std::endl;
      }

#ifdef JETTAGGERFAST_TEST
      if(index>-1){
        Double_t distanceTest(-1);
        distanceTest = TMath::Sqrt(TMath::Power(etaBase[index] - j->Eta(), 2) +  TMath::Power(phiBase[index] - j->Phi(), 2));
        if(TMath::Abs(distanceTest - distance) > DBL_EPSILON){
          AliDebugStream(1) << "Mismatch in distance from base jet with index from tree: " << distanceTest << ", distance from tree " << distance << std::endl;
          fIndexErrorRateTag->Fill(1);
        }
      }
#endif

      countTag++;
    }

    // check for "true" correlations
    // these are pairs where the base jet is the closest to the tag jet and vice versa
    // As the lists are linear a loop over the outer base jet is sufficient.
    AliDebugStream(1) << "Starting true jet loop: nbase(" << kNacceptedBase << "), ntag(" << kNacceptedTag << ")\n";
    for(int ibase = 0; ibase < kNacceptedBase; ibase++) {
      AliDebugStream(2) << "base jet " << ibase << ": match index in tag jet container " << faMatchIndexTag[ibase] << "\n";
      if(faMatchIndexTag[ibase] > -1){
       AliDebugStream(2) << "tag jet " << faMatchIndexTag[ibase] << ": matched base jet " << faMatchIndexBase[faMatchIndexTag[ibase]] << "\n";
      }
      if(faMatchIndexTag[ibase] > -1 && faMatchIndexBase[faMatchIndexTag[ibase]] == ibase) {
        AliDebugStream(2) << "found a true match \n";
        AliEmcalJet *jetBase = jetsBase[ibase],
                    *jetTag = jetsTag[faMatchIndexTag[ibase]];
        if(jetBase && jetTag) {
#ifdef JETTAGGERFAST_TEST
          if(TMath::Abs(etaBase[ibase] - jetBase->Eta()) > DBL_EPSILON || TMath::Abs(phiBase[ibase] - jetBase->Phi()) > DBL_EPSILON){
            AliErrorStream() << "Selected incorrect base jet for tagging : eta test(" << jetBase->Eta() << ")/true(" << etaBase[ibase]
                             << "), phi test(" << jetBase->Phi() << ")/true(" << phiBase[ibase] << ")\n";
            fContainerErrorRateBase->Fill(1);
          }
          if(TMath::Abs(etaTag[faMatchIndexTag[ibase]] - jetTag->Eta()) > DBL_EPSILON || TMath::Abs(phiTag[faMatchIndexTag[ibase]] - jetTag->Phi()) > DBL_EPSILON){
            AliErrorStream() << "Selected incorrect tag jet for tagging : eta test(" << jetTag->Eta() << ")/true(" << etaTag[faMatchIndexTag[ibase]]
                             << "), phi test(" << jetTag->Phi() << ")/true(" << phiTag[faMatchIndexTag[ibase]] << ")\n";
            fContainerErrorRateTag->Fill(1);
          }
#endif
          // Test if the position of the jets correp
          Double_t dR = jetBase->DeltaR(jetTag);
          switch(fJetTaggingType){
          case kTag:
            jetBase->SetTaggedJet(jetTag);
            jetBase->SetTagStatus(1);

            jetTag->SetTaggedJet(jetBase);
            jetTag->SetTagStatus(1);
            break;
          case kClosest:
            jetBase->SetClosestJet(jetTag,dR);
            jetTag->SetClosestJet(jetBase,dR);
            break;
          };
        }
      }
    }
    return kTRUE;
  }

  Double_t AliEmcalJetTaggerTaskFast::GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2) {
    return GetDeltaPhi(jet1->Phi(),jet2->Phi());
  }

  Double_t AliEmcalJetTaggerTaskFast::GetDeltaPhi(Double_t phi1,Double_t phi2) {
    Double_t dPhi = phi1-phi2;
    if(dPhi <-0.5*TMath::Pi())  dPhi += TMath::TwoPi();
    if(dPhi > 1.5*TMath::Pi())  dPhi -= TMath::TwoPi();

    return dPhi;
  }

  AliEmcalJetTaggerTaskFast *AliEmcalJetTaggerTaskFast::AddTaskJetTaggerFast(const char * njetsBase,
      const char * njetsTag,
      const Double_t R,
      const char * nrhoBase,
      const char * nrhoTag,
      const char * ntracks,
      const char * nclusters,
      const char * type,
      const char * CentEst,
      Int_t        pSel,
      const char * trigClass){

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      std::cerr << "E-AddEmcalJetTaggerTaskFast: No analysis manager found.\n";
      return nullptr;
    }
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
      std::cerr << "E-AddEmcalJetTaggerTaskFast: This task requires an input event handler\n";
      return NULL;
    }

    TString wagonName = Form("JetTagger_%s_%s_TC%s",njetsBase,njetsTag,trigClass);

    //Configure jet tagger task
    AliEmcalJetTaggerTaskFast *task = new AliEmcalJetTaggerTaskFast(wagonName);

    task->SetNCentBins(4);
    //task->SetVzRange(-10.,10.);

    AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
    AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

    task->SetJetContainerBase(0);
    task->SetJetContainerTag(1);

    TString strType(type);
    AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetMaxTrackPt(10000.);
    }

    AliJetContainer *jetContTag = task->AddJetContainer(njetsTag,"TPC",R);
    if(jetContTag) {
      jetContTag->SetRhoName(nrhoTag);
      jetContTag->ConnectParticleContainer(trackCont);
      jetContTag->ConnectClusterContainer(clusterCont);
      jetContTag->SetMaxTrackPt(10000.);
    }
    for(Int_t i=0; i<2; i++) {
      task->SetPercAreaCut(0.6, i); //keep?
    }
    task->SetCentralityEstimator(CentEst);
    task->SelectCollisionCandidates(pSel);
    task->SetUseAliAnaUtils(kFALSE);

    mgr->AddTask(task);

    //Connnect input
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

    //Connect output
    TString contName(wagonName);
    TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,1,coutput1);

    return task;
  }


  }
}

