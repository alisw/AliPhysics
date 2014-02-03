//
// Jet tagger analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"

#include "AliAnalysisTaskEmcalJetTagger.h"

ClassImp(AliAnalysisTaskEmcalJetTagger)

//________________________________________________________________________
AliAnalysisTaskEmcalJetTagger::AliAnalysisTaskEmcalJetTagger() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetTagger", kTRUE),
  fJetTaggingType(kTag),
  fJetTaggingMethod(kGeo),
  fContainerBase(0),
  fContainerTag(1),
  fMatchingDone(0),
  fh3PtJet1VsDeltaEtaDeltaPhi(0),
  fh2PtJet1VsDeltaR(0),
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtJet1VsPtJet2(0),
  fh3PtJetDEtaDPhiConst(0),
  fh2PtJetDRConst(0),
  fh3PtJetAreaDRConst(0)
{
  // Default constructor.

  fh3PtJet1VsDeltaEtaDeltaPhi  = new TH3F*[fNcentBins];
  fh2PtJet1VsDeltaR            = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtAllSel      = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtTagged      = new TH2F*[fNcentBins];
  fh2PtJet1VsPtJet2            = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1VsDeltaEtaDeltaPhi[i] = 0;
    fh2PtJet1VsDeltaR[i]           = 0;
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtJet1VsPtJet2[i]           = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTagger::AliAnalysisTaskEmcalJetTagger(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fJetTaggingType(kTag),
  fJetTaggingMethod(kGeo),
  fContainerBase(0),
  fContainerTag(1),
  fMatchingDone(0),
  fh3PtJet1VsDeltaEtaDeltaPhi(0),
  fh2PtJet1VsDeltaR(0),
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtJet1VsPtJet2(0),
  fh3PtJetDEtaDPhiConst(0),
  fh2PtJetDRConst(0),
  fh3PtJetAreaDRConst(0)
{
  // Standard constructor.

  fh3PtJet1VsDeltaEtaDeltaPhi = new TH3F*[fNcentBins];
  fh2PtJet1VsDeltaR           = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtAllSel     = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtTagged     = new TH2F*[fNcentBins];
  fh2PtJet1VsPtJet2           = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1VsDeltaEtaDeltaPhi[i] = 0;
    fh2PtJet1VsDeltaR[i]           = 0;
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtJet1VsPtJet2[i]           = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTagger::~AliAnalysisTaskEmcalJetTagger()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTagger::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt          = 250;
  const Int_t nBinsDPhi        = 100;
  const Int_t nBinsDEta        = 100;
  const Int_t nBinsDR          = 50;

  const Double_t minPt = -50.;
  const Double_t maxPt = 200.;
  const Double_t minDPhi = -0.5;
  const Double_t maxDPhi =  0.5;
  const Double_t minDEta = -0.5;
  const Double_t maxDEta =  0.5;
  const Double_t minDR   =  0.;
  const Double_t maxDR   =  0.5;

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

  }

  fh3PtJetDEtaDPhiConst = new TH3F("fh3PtJetDEtaDPhiConst","fh3PtJetDEtaDPhiConst;pT;#Delta #eta;#Delta #varphi",nBinsPt,minPt,maxPt,nBinsDEta,-1.,1.,nBinsDPhi,-1.,1.);
  fOutput->Add(fh3PtJetDEtaDPhiConst);

  fh2PtJetDRConst = new TH2F("fh2PtJetDRConst","fh2PtJetDRConst;pT;#Delta R",nBinsPt,minPt,maxPt,100,0.,1.);
  fOutput->Add(fh2PtJetDRConst);

  fh3PtJetAreaDRConst = new TH3F("fh3PtJetAreaDRConst","fh3PtJetAreaDRConst;pT;A;#Delta R",nBinsPt,minPt,maxPt,100,0.,1.,100,0.,1.);
  fOutput->Add(fh3PtJetAreaDRConst);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTagger::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  if(fJetTaggingMethod==kGeo)
    MatchJetsGeo(fContainerBase,fContainerTag,0,0.3,2);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTagger::FillHistograms()
{
  // Fill histograms.

  for(int i = 0; i < GetNJets(fContainerBase);++i) {
    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, fContainerBase));
    if(!jet1) continue;

    Double_t ptJet1 =  jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
    fh2PtJet1VsLeadPtAllSel[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());

    //fill histo with angle between jet axis and constituents
    for(Int_t icc=0; icc<jet1->GetNumberOfTracks(); icc++) {
      AliVParticle *vp = static_cast<AliVParticle*>(jet1->TrackAt(icc, fTracks));
      if(!vp) continue;
      Double_t dEta = jet1->Eta()-vp->Eta();
      Double_t dPhi = jet1->Phi()-vp->Phi();
      if(dPhi<TMath::Pi()) dPhi+=TMath::TwoPi();
      if(dPhi>TMath::Pi()) dPhi-=TMath::TwoPi();
      fh3PtJetDEtaDPhiConst->Fill(ptJet1,dEta,dPhi);

      Double_t dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
      fh2PtJetDRConst->Fill(ptJet1,dR);
      fh3PtJetAreaDRConst->Fill(ptJet1,jet1->Area(),dR);
    }

    if(jet1->GetTagStatus()<1 && fJetTaggingType==kTag)
      continue;

    AliEmcalJet *jet2 = jet1->GetTaggedJet();
    if(!jet2) continue;
    Double_t ptJet2 =  jet2->Pt() - GetRhoVal(fContainerTag)*jet2->Area();
    fh2PtJet1VsLeadPtTagged[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());

    fh2PtJet1VsPtJet2[fCentBin]->Fill(ptJet1,ptJet2);

    Double_t dPhi = GetDeltaPhi(jet1->Phi(),jet2->Phi());
    if(dPhi>TMath::Pi())
      dPhi -= TMath::TwoPi();
    if(dPhi<(-1.*TMath::Pi()))
      dPhi += TMath::TwoPi();  
    
    fh3PtJet1VsDeltaEtaDeltaPhi[fCentBin]->Fill(ptJet1,jet1->Eta()-jet2->Eta(),dPhi);
    fh2PtJet1VsDeltaR[fCentBin]->Fill(ptJet1,GetDeltaR(jet1,jet2));
  }

  return kTRUE;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTagger::ResetTagging(const Int_t c) {

  //Reset tagging of container c

  for(int i = 0;i<GetNJets(c);i++){
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(GetJetFromArray(i, c));
    if(fJetTaggingType==kClosest)
      jet->ResetMatching();
    else if(fJetTaggingType==kTag) {
      jet->SetTaggedJet(0x0);
      jet->SetTagStatus(-1);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTagger::MatchJetsGeo(Int_t c1, Int_t c2,
						 Int_t iDebug, Float_t maxDist, Int_t type) {

  //
  // Match the full jets to the corresponding charged jets
  // Translation of AliAnalysisHelperJetTasks::GetClosestJets to AliEmcalJet objects
  // type: 0 = use acceptance cuts of container  1 = allow 0.1 one more for c2 in eta 2 = allow 0.1 more in eta and phi for c2

  if(c1<0) c1 = fContainerBase;
  if(c2<0) c2 = fContainerTag;

  const Int_t nJets1 = GetNJets(c1);
  const Int_t nJets2 = GetNJets(c2);

  if(nJets1==0 || nJets2==0) return;

  ResetTagging(c1);
  ResetTagging(c2);

  fMatchingDone = kFALSE;

  TArrayI faMatchIndex1;
  faMatchIndex1.Set(nJets2+1);
  faMatchIndex1.Reset(-1);

  TArrayI faMatchIndex2;
  faMatchIndex2.Set(nJets1+1);
  faMatchIndex2.Reset(-1);

  static TArrayS iFlag(nJets1*nJets2);
  if(iFlag.GetSize()<(nJets1*nJets2)){
    iFlag.Set(nJets1*nJets1+1);
  }
  iFlag.Reset(0);

  AliJetContainer *cont2 = GetJetContainer(c2);

  // find the closest distance to the full jet
  for(int i = 0;i<nJets1;i++){

    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, c1));
    if(!jet1) continue;

    Float_t dist = maxDist;
    
    for(int j = 0;j <nJets2; j++){
      AliEmcalJet *jet2 = 0x0;
      if(type==0)
	jet2 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(j, c2));
      else {
	jet2 = static_cast<AliEmcalJet*>(GetJetFromArray(j, c2));
	if(!jet2) continue;
	if(type>0) {
	  if(jet2->Eta()<(cont2->GetJetEtaMin()-0.1) || jet2->Eta()>(cont2->GetJetEtaMax()+0.1))
	    continue;
	  if(type==2) {
	    if(jet2->Phi()<(cont2->GetJetPhiMin()-0.1) || jet2->Phi()>(cont2->GetJetPhiMax()+0.1))
	      continue;
	  }
	}
      }
      if(!jet2)
	continue;

      Double_t dR = GetDeltaR(jet1,jet2);
      if(dR<dist && dR<maxDist){
	faMatchIndex2[i]=j;
	dist = dR;
      }
    }//j jet loop
    if(faMatchIndex2[i]>=0) {
      iFlag[i*nJets2+faMatchIndex2[i]]+=1;//j closest to i
      if(iDebug>10) Printf("Full Distance (%d)--(%d) %3.3f flag[%d] = %d",i,faMatchIndex2[i],dist,i*nJets2+faMatchIndex2[i],iFlag[i*nJets2+faMatchIndex2[i]]);
    }
  
  }//i jet loop


  // other way around
  for(int j = 0;j<nJets2;j++){
    AliEmcalJet *jet2 = 0x0;
    if(type==0)
      jet2 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(j, c2));
    else {
      jet2 = static_cast<AliEmcalJet*>(GetJetFromArray(j, c2));
      if(!jet2) continue;;
      if(type>0) {
	if(jet2->Eta()<(cont2->GetJetEtaMin()-0.1) || jet2->Eta()>(cont2->GetJetEtaMax()+0.1))
	  continue;
	if(type==2) {
	  if(jet2->Phi()<(cont2->GetJetPhiMin()-0.1) || jet2->Phi()>(cont2->GetJetPhiMax()+0.1))
	    continue;
	}
      }
    }
    if(!jet2)
      continue;

    Float_t dist = maxDist;
    for(int i = 0;i<nJets1;i++){
      AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, c1));
      if(!jet1)
	continue;

      Double_t dR = GetDeltaR(jet1,jet2); 
      if(dR<dist && dR<maxDist){
	faMatchIndex1[j]=i;
        dist = dR;
      }   
    }
    if(faMatchIndex1[j]>=0) {
      iFlag[faMatchIndex1[j]*nJets2+j]+=2;//i closest to j
      if(iDebug>10) Printf("Other way Distance (%d)--(%d) %3.3f flag[%d] = %d",faMatchIndex1[j],j,dist,faMatchIndex1[j]*nJets2+j,iFlag[faMatchIndex1[j]*nJets2+j]);
    }

  }
    
  // check for "true" correlations
  for(int i = 0;i<nJets1;i++){
    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetJetFromArray(i, c1));
    for(int j = 0;j<nJets2;j++){
      AliEmcalJet *jet2 = static_cast<AliEmcalJet*>(GetJetFromArray(j, c2));
      AliDebug(11,Form("%s: Flag[%d][%d] %d ",GetName(),i,j,iFlag[i*nJets2+j]));
      
      // we have a uniqe correlation
      if(iFlag[i*nJets2+j]==3){
	
	Double_t dR = GetDeltaR(jet1,jet2); 
	if(iDebug>1) Printf("closest jets %d  %d  dR =  %f",j,i,dR);
	
	if(fJetTaggingType==kTag) {
	  jet1->SetTaggedJet(jet2);
	  jet1->SetTagStatus(1);
	  
	  jet2->SetTaggedJet(jet1);
	  jet2->SetTagStatus(1);
	}
	else if(fJetTaggingType==kClosest) {
	  jet1->SetClosestJet(jet2,dR);
	  jet2->SetClosestJet(jet1,dR);
	}
      }
    }
  }

  fMatchingDone = kTRUE;
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTagger::GetDeltaR(const AliEmcalJet* jet1, const AliEmcalJet* jet2) const {
  //
  // Helper function to calculate the distance between two jets
  //

  Double_t dPhi = jet1->Phi() - jet2->Phi();
  if(dPhi>TMath::Pi())
    dPhi -= TMath::TwoPi();
  if(dPhi<(-1.*TMath::Pi()))
    dPhi += TMath::TwoPi();
  Double_t dEta = jet1->Eta() - jet2->Eta();
  Double_t dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
  return dR;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTagger::GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2) {
  //
  // Calculate azimuthal angle between the axises of the jets
  //

  return GetDeltaPhi(jet1->Phi(),jet2->Phi());

}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTagger::GetDeltaPhi(Double_t phi1,Double_t phi2) {
  //
  // Calculate azimuthal angle between the axises of the jets
  //

  Double_t dPhi = phi1-phi2;

  if(dPhi <-0.5*TMath::Pi())  dPhi += TMath::TwoPi();
  if(dPhi > 1.5*TMath::Pi())  dPhi -= TMath::TwoPi();

  return dPhi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTagger::GetFractionSharedPt(const AliEmcalJet *jet1, const AliEmcalJet *jet2) const
{
  //
  // Get fraction of shared pT between matched full and charged jet
  // Uses charged jet pT as baseline: fraction = \Sum_{const,full jet} pT,const,i / pT,jet,ch
  //

  Double_t fraction = 0.;
  Double_t jetPt2 = jet2->Pt();
 
  if(jetPt2>0) {

    AliDebug(2,Form("%s: nConstituents: %d, ch: %d  chne: %d ne: %d",GetName(),jet1->GetNumberOfConstituents(),jet2->GetNumberOfTracks(),jet1->GetNumberOfTracks(),jet1->GetNumberOfClusters()));
    
    Double_t sumPt = 0.;
    AliVParticle *vpf = 0x0;
    Int_t iFound = 0;
    for(Int_t icc=0; icc<jet2->GetNumberOfTracks(); icc++) {
      Int_t idx = (Int_t)jet2->TrackAt(icc);
      iFound = 0;
      for(Int_t icf=0; icf<jet1->GetNumberOfTracks(); icf++) {
	if(idx == jet1->TrackAt(icf) && iFound==0 ) {
	  iFound=1;
	  vpf = static_cast<AliVParticle*>(jet1->TrackAt(icf, fTracks));
	  if(!vpf) continue;
	  if(vpf->Charge()!=0)
	    sumPt += vpf->Pt();
	  continue;
	}
      }
    }
    
    fraction = sumPt/jetPt2;
  }

  return fraction;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTagger::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetTagger::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

