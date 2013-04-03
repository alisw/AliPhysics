// $Id$
//
// Emcal jet response matrix maker task.
//
// Author: S. Aiola

#include "AliJetResponseMaker.h"

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>
#include <TProfile.h>

#include "AliAnalysisManager.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliNamedArrayI.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
  fTracks2Name(""),
  fCalo2Name(""),
  fJets2Name(""),
  fRho2Name(""),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar1(0),
  fMatchingPar2(0),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fSelectPtHardBin(-999),
  fIsEmbedded(kFALSE),
  fIsPythia(kTRUE),
  fMCLabelShift(0),
  fUseCellsToMatch(kFALSE),
  fPythiaHeader(0),
  fPtHardBin(0),
  fNTrials(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fTracks2Map(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistEvents(0),
  fMCEnergy1vsEnergy2(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistLeadingJets1PtArea(0),
  fHistLeadingJets1CorrPtArea(0),
  fHistJets1NEFvsPt(0),
  fHistJets1CEFvsCEFPt(0),
  fHistJets1ZvsPt(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistLeadingJets2PtArea(0),
  fHistLeadingJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistLeadingJets2PtAreaAcceptance(0),
  fHistLeadingJets2CorrPtAreaAcceptance(0),
  fHistJets2NEFvsPt(0),
  fHistJets2CEFvsCEFPt(0),
  fHistJets2ZvsPt(0),
  fHistNonMatchedJets1PtArea(0),
  fHistNonMatchedJets2PtArea(0),
  fHistNonMatchedJets1CorrPtArea(0),
  fHistNonMatchedJets2CorrPtArea(0),
  fHistMissedJets2PtArea(0),
  fHistCommonEnergy1vsJet1Pt(0),
  fHistCommonEnergy2vsJet2Pt(0),
  fHistDistancevsJet1Pt(0),
  fHistDistancevsJet2Pt(0),
  fHistDistancevsCommonEnergy1(0),
  fHistDistancevsCommonEnergy2(0),
  fHistCommonEnergy1vsCommonEnergy2(0),
  fHistJet2PtOverJet1PtvsJet2Pt(0),
  fHistJet1PtOverJet2PtvsJet1Pt(0),
  fHistDeltaEtaPhi(0),
  fHistDeltaPtvsJet1Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtvsDistance(0),
  fHistDeltaPtvsCommonEnergy1(0),
  fHistDeltaPtvsCommonEnergy2(0),
  fHistDeltaPtvsArea1(0),
  fHistDeltaPtvsArea2(0),
  fHistDeltaPtvsDeltaArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistDeltaCorrPtvsJet1Pt(0),
  fHistDeltaCorrPtvsJet2Pt(0),
  fHistDeltaCorrPtvsDistance(0),
  fHistDeltaCorrPtvsCommonEnergy1(0),
  fHistDeltaCorrPtvsCommonEnergy2(0),
  fHistDeltaCorrPtvsArea1(0),
  fHistDeltaCorrPtvsArea2(0),
  fHistDeltaCorrPtvsDeltaArea(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistDeltaMCPtvsJet1Pt(0),
  fHistDeltaMCPtvsJet2Pt(0),
  fHistDeltaMCPtvsDistance(0),
  fHistDeltaMCPtvsCommonEnergy1(0),
  fHistDeltaMCPtvsCommonEnergy2(0),
  fHistDeltaMCPtvsArea1(0),
  fHistDeltaMCPtvsArea2(0),
  fHistDeltaMCPtvsDeltaArea(0),
  fHistJet1MCPtvsJet2Pt(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fTracks2Name("MCParticles"),
  fCalo2Name(""),
  fJets2Name("MCJets"),
  fRho2Name(""),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar1(0),
  fMatchingPar2(0),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fSelectPtHardBin(-999),
  fIsEmbedded(kFALSE),
  fIsPythia(kTRUE),
  fMCLabelShift(0),
  fUseCellsToMatch(kFALSE),
  fPythiaHeader(0),
  fPtHardBin(0),
  fNTrials(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fTracks2Map(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistEvents(0),
  fMCEnergy1vsEnergy2(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistLeadingJets1PtArea(0),
  fHistLeadingJets1CorrPtArea(0),
  fHistJets1NEFvsPt(0),
  fHistJets1CEFvsCEFPt(0),
  fHistJets1ZvsPt(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistLeadingJets2PtArea(0),
  fHistLeadingJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistLeadingJets2PtAreaAcceptance(0),
  fHistLeadingJets2CorrPtAreaAcceptance(0),
  fHistJets2NEFvsPt(0),
  fHistJets2CEFvsCEFPt(0),
  fHistJets2ZvsPt(0),
  fHistNonMatchedJets1PtArea(0),
  fHistNonMatchedJets2PtArea(0),
  fHistNonMatchedJets1CorrPtArea(0),
  fHistNonMatchedJets2CorrPtArea(0),
  fHistMissedJets2PtArea(0),
  fHistCommonEnergy1vsJet1Pt(0),
  fHistCommonEnergy2vsJet2Pt(0),
  fHistDistancevsJet1Pt(0),
  fHistDistancevsJet2Pt(0),
  fHistDistancevsCommonEnergy1(0),
  fHistDistancevsCommonEnergy2(0),
  fHistCommonEnergy1vsCommonEnergy2(0),
  fHistJet2PtOverJet1PtvsJet2Pt(0),
  fHistJet1PtOverJet2PtvsJet1Pt(0),
  fHistDeltaEtaPhi(0),
  fHistDeltaPtvsJet1Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtvsDistance(0),
  fHistDeltaPtvsCommonEnergy1(0),
  fHistDeltaPtvsCommonEnergy2(0),
  fHistDeltaPtvsArea1(0),
  fHistDeltaPtvsArea2(0),
  fHistDeltaPtvsDeltaArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistDeltaCorrPtvsJet1Pt(0),
  fHistDeltaCorrPtvsJet2Pt(0),
  fHistDeltaCorrPtvsDistance(0),
  fHistDeltaCorrPtvsCommonEnergy1(0),
  fHistDeltaCorrPtvsCommonEnergy2(0),
  fHistDeltaCorrPtvsArea1(0),
  fHistDeltaCorrPtvsArea2(0),
  fHistDeltaCorrPtvsDeltaArea(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistDeltaMCPtvsJet1Pt(0),
  fHistDeltaMCPtvsJet2Pt(0),
  fHistDeltaMCPtvsDistance(0),
  fHistDeltaMCPtvsCommonEnergy1(0),
  fHistDeltaMCPtvsCommonEnergy2(0),
  fHistDeltaMCPtvsArea1(0),
  fHistDeltaMCPtvsArea2(0),
  fHistDeltaMCPtvsDeltaArea(0),
  fHistJet1MCPtvsJet2Pt(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::~AliJetResponseMaker()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard)
{
  //
  // Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // Get the pt hard bin from the file path
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // (Partially copied from AliAnalysisHelperJetTasks)

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if(file.Contains(".zip#")){
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  Printf("%s",file.Data());

  // Get the pt hard bin
  TString strPthard(file);
  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec()) 
    pthard = strPthard.Atoi();
  else 
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));

  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec){
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec){
	// not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else{
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if(!key){
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list){
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::UserNotify()
{
  if (!fIsPythia)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection = 0;
  Float_t trials   = 0;
  Int_t   pthard   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain)
    tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthard);

  fHistTrials->Fill(pthard, trials);
  fHistXsection->Fill(pthard, xsection);
  fHistEvents->Fill(pthard, nevents);

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::UserCreateOutputObjects()
{
  // Create user objects.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // General histograms
  if (fIsPythia) {
    fHistTrialsAfterSel = new TH1F("fHistTrialsAfterSel", "fHistTrialsAfterSel", 11, 0, 11);
    fHistTrialsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrialsAfterSel->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsAfterSel);
    
    fHistEventsAfterSel = new TH1F("fHistEventsAfterSel", "fHistEventsAfterSel", 11, 0, 11);
    fHistEventsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEventsAfterSel->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEventsAfterSel);
    
    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
    fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
    fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
    fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

    const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
    const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
    
    for (Int_t i = 1; i < 12; i++) {
      fHistTrialsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEventsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      
      fHistTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistXsection->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    }
  }

  if (fIsEmbedded) {
    fMCEnergy1vsEnergy2 = new TH2F("fMCEnergy1vsEnergy2", "fMCEnergy1vsEnergy2", fNbins, fMinBinPt, fMaxBinPt*4, fNbins, fMinBinPt, fMaxBinPt*4);
    fMCEnergy1vsEnergy2->GetXaxis()->SetTitle("#Sigmap_{T,1}^{MC}");
    fMCEnergy1vsEnergy2->GetYaxis()->SetTitle("#Sigmap_{T,2}");
    fMCEnergy1vsEnergy2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fMCEnergy1vsEnergy2);
  }

  // Jets 1 histograms
  fHistLeadingJets1PtArea = new TH2F("fHistLeadingJets1PtArea", "fHistLeadingJets1PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets1PtArea->GetXaxis()->SetTitle("area");
  fHistLeadingJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistLeadingJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets1PtArea);

  fHistJets1PhiEta = new TH2F("fHistJets1PhiEta", "fHistJets1PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets1PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets1PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets1PhiEta);
  
  fHistJets1PtArea = new TH2F("fHistJets1PtArea", "fHistJets1PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1PtArea->GetXaxis()->SetTitle("area");
  fHistJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1PtArea);

  if (!fRhoName.IsNull()) {
    fHistLeadingJets1CorrPtArea = new TH2F("fHistLeadingJets1CorrPtArea", "fHistLeadingJets1CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistLeadingJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1}^{corr} (GeV/c)");
    fHistLeadingJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets1CorrPtArea);

    fHistJets1CorrPtArea = new TH2F("fHistJets1CorrPtArea", "fHistJets1CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1}^{corr} (GeV/c)");
    fHistJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets1CorrPtArea);
  }

  fHistJets1ZvsPt = new TH2F("fHistJets1ZvsPt", "fHistJets1ZvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1ZvsPt->GetXaxis()->SetTitle("Z");
  fHistJets1ZvsPt->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1ZvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1ZvsPt);
  
  fHistJets1NEFvsPt = new TH2F("fHistJets1NEFvsPt", "fHistJets1NEFvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1NEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistJets1NEFvsPt->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1NEFvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1NEFvsPt);
  
  fHistJets1CEFvsCEFPt = new TH2F("fHistJets1CEFvsCEFPt", "fHistJets1CEFvsCEFPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1CEFvsCEFPt->GetXaxis()->SetTitle("1-NEF");
  fHistJets1CEFvsCEFPt->GetYaxis()->SetTitle("(1-NEF)*p_{T,1} (GeV/c)");
  fHistJets1CEFvsCEFPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1CEFvsCEFPt);

  // Jets 2 histograms

  fHistLeadingJets2PtAreaAcceptance = new TH2F("fHistLeadingJets2PtAreaAcceptance", "fHistLeadingJets2PtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets2PtAreaAcceptance->GetXaxis()->SetTitle("area");
  fHistLeadingJets2PtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistLeadingJets2PtAreaAcceptance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets2PtAreaAcceptance);

  fHistLeadingJets2PtArea = new TH2F("fHistLeadingJets2PtArea", "fHistLeadingJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets2PtArea->GetXaxis()->SetTitle("area");
  fHistLeadingJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistLeadingJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets2PtArea);

  fHistJets2PhiEtaAcceptance = new TH2F("fHistJets2PhiEtaAcceptance", "fHistJets2PhiEtaAcceptance", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEtaAcceptance->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEtaAcceptance->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEtaAcceptance);

  fHistJets2PtAreaAcceptance = new TH2F("fHistJets2PtAreaAcceptance", "fHistJets2PtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtAreaAcceptance->GetXaxis()->SetTitle("area");
  fHistJets2PtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtAreaAcceptance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtAreaAcceptance);

  fHistJets2PhiEta = new TH2F("fHistJets2PhiEta", "fHistJets2PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEta);

  fHistJets2PtArea = new TH2F("fHistJets2PtArea", "fHistJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtArea->GetXaxis()->SetTitle("area");
  fHistJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtArea);

  if (!fRho2Name.IsNull()) {
    fHistJets2CorrPtAreaAcceptance = new TH2F("fHistJets2CorrPtAreaAcceptance", "fHistJets2CorrPtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtAreaAcceptance->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtAreaAcceptance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtAreaAcceptance);

    fHistLeadingJets2CorrPtAreaAcceptance = new TH2F("fHistLeadingJets2CorrPtAreaAcceptance", "fHistLeadingJets2CorrPtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets2CorrPtAreaAcceptance->GetXaxis()->SetTitle("area");
    fHistLeadingJets2CorrPtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistLeadingJets2CorrPtAreaAcceptance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets2CorrPtAreaAcceptance);

    fHistLeadingJets2CorrPtArea = new TH2F("fHistLeadingJets2CorrPtArea", "fHistLeadingJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistLeadingJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistLeadingJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets2CorrPtArea);

    fHistJets2CorrPtArea = new TH2F("fHistJets2CorrPtArea", "fHistJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtArea);
  }

  fHistJets2ZvsPt = new TH2F("fHistJets2ZvsPt", "fHistJets2ZvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2ZvsPt->GetXaxis()->SetTitle("Z");
  fHistJets2ZvsPt->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2ZvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2ZvsPt);
  
  fHistJets2NEFvsPt = new TH2F("fHistJets2NEFvsPt", "fHistJets2NEFvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2NEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistJets2NEFvsPt->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2NEFvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2NEFvsPt);
  
  fHistJets2CEFvsCEFPt = new TH2F("fHistJets2CEFvsCEFPt", "fHistJets2CEFvsCEFPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2CEFvsCEFPt->GetXaxis()->SetTitle("1-NEF");
  fHistJets2CEFvsCEFPt->GetYaxis()->SetTitle("(1-NEF)*p_{T,2} (GeV/c)");
  fHistJets2CEFvsCEFPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2CEFvsCEFPt);

  // Matching histograms

  fHistNonMatchedJets1PtArea = new TH2F("fHistNonMatchedJets1PtArea", "fHistNonMatchedJets1PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJets1PtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistNonMatchedJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJets1PtArea);

  fHistNonMatchedJets2PtArea = new TH2F("fHistNonMatchedJets2PtArea", "fHistNonMatchedJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJets2PtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistNonMatchedJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJets2PtArea);

  if (!fRhoName.IsNull()) {  
    fHistNonMatchedJets1CorrPtArea = new TH2F("fHistNonMatchedJets1CorrPtArea", "fHistNonMatchedJets1CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistNonMatchedJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistNonMatchedJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
    fHistNonMatchedJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistNonMatchedJets1CorrPtArea);
  }

  if (!fRho2Name.IsNull()) {  
    fHistNonMatchedJets2CorrPtArea = new TH2F("fHistNonMatchedJets2CorrPtArea", "fHistNonMatchedJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistNonMatchedJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistNonMatchedJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
    fHistNonMatchedJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistNonMatchedJets2CorrPtArea);
  }

  fHistMissedJets2PtArea = new TH2F("fHistMissedJets2PtArea", "fHistMissedJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistMissedJets2PtArea->GetXaxis()->SetTitle("area");  
  fHistMissedJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistMissedJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMissedJets2PtArea);

  fHistCommonEnergy1vsJet1Pt = new TH2F("fHistCommonEnergy1vsJet1Pt", "fHistCommonEnergy1vsJet1Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistCommonEnergy1vsJet1Pt->GetXaxis()->SetTitle("Common energy 1 (%)");
  fHistCommonEnergy1vsJet1Pt->GetYaxis()->SetTitle("p_{T,1}");  
  fHistCommonEnergy1vsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy1vsJet1Pt);

  fHistCommonEnergy2vsJet2Pt = new TH2F("fHistCommonEnergy2vsJet2Pt", "fHistCommonEnergy2vsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistCommonEnergy2vsJet2Pt->GetXaxis()->SetTitle("Common energy 2 (%)");
  fHistCommonEnergy2vsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistCommonEnergy2vsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy2vsJet2Pt);

  fHistDistancevsJet1Pt = new TH2F("fHistDistancevsJet1Pt", "fHistDistancevsJet1Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistDistancevsJet1Pt->GetXaxis()->SetTitle("Distance");
  fHistDistancevsJet1Pt->GetYaxis()->SetTitle("p_{T,1}");  
  fHistDistancevsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsJet1Pt);

  fHistDistancevsJet2Pt = new TH2F("fHistDistancevsJet2Pt", "fHistDistancevsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistDistancevsJet2Pt->GetXaxis()->SetTitle("Distance");
  fHistDistancevsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistDistancevsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsJet2Pt);

  fHistDistancevsCommonEnergy1 = new TH2F("fHistDistancevsCommonEnergy1", "fHistDistancevsCommonEnergy1", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistDistancevsCommonEnergy1->GetXaxis()->SetTitle("Distance");
  fHistDistancevsCommonEnergy1->GetYaxis()->SetTitle("Common energy 1 (%)");  
  fHistDistancevsCommonEnergy1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsCommonEnergy1);

  fHistDistancevsCommonEnergy2 = new TH2F("fHistDistancevsCommonEnergy2", "fHistDistancevsCommonEnergy2", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistDistancevsCommonEnergy2->GetXaxis()->SetTitle("Distance");
  fHistDistancevsCommonEnergy2->GetYaxis()->SetTitle("Common energy 2 (%)");  
  fHistDistancevsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsCommonEnergy2);

  fHistCommonEnergy1vsCommonEnergy2 = new TH2F("fHistCommonEnergy1vsCommonEnergy2", "fHistCommonEnergy1vsCommonEnergy2", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistCommonEnergy1vsCommonEnergy2->GetXaxis()->SetTitle("Common energy 1 (%)");
  fHistCommonEnergy1vsCommonEnergy2->GetYaxis()->SetTitle("Common energy 2 (%)");  
  fHistCommonEnergy1vsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy1vsCommonEnergy2);

  fHistJet2PtOverJet1PtvsJet2Pt = new TH2F("fHistJet2PtOverJet1PtvsJet2Pt", "fHistJet2PtOverJet1PtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, 300, 0, 1.5);
  fHistJet2PtOverJet1PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistJet2PtOverJet1PtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2} / p_{T,1}");
  fHistJet2PtOverJet1PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet2PtOverJet1PtvsJet2Pt);

  fHistJet1PtOverJet2PtvsJet1Pt = new TH2F("fHistJet1PtOverJet2PtvsJet1Pt", "fHistJet1PtOverJet2PtvsJet1Pt", fNbins, fMinBinPt, fMaxBinPt, 300, 0, 1.5);
  fHistJet1PtOverJet2PtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
  fHistJet1PtOverJet2PtvsJet1Pt->GetYaxis()->SetTitle("p_{T,1} / p_{T,2}");
  fHistJet1PtOverJet2PtvsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet1PtOverJet2PtvsJet1Pt);

  fHistDeltaEtaPhi = new TH2F("fHistDeltaEtaPhi", "fHistDeltaEtaPhi", 200, -1, 1, 250, -1.6, 4.8);
  fHistDeltaEtaPhi->GetXaxis()->SetTitle("#delta#eta");
  fHistDeltaEtaPhi->GetYaxis()->SetTitle("#delta#phi");
  fHistDeltaEtaPhi->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaEtaPhi);

  fHistDeltaPtvsJet1Pt = new TH2F("fHistDeltaPtvsJet1Pt", "fHistDeltaPtvsJet1Pt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
  fHistDeltaPtvsJet1Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsJet1Pt);

  fHistDeltaPtvsJet2Pt = new TH2F("fHistDeltaPtvsJet2Pt", "fHistDeltaPtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistDeltaPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsJet2Pt);

  fHistDeltaPtvsDistance = new TH2F("fHistDeltaPtvsDistance", "fHistDeltaPtvsDistance", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsDistance->GetXaxis()->SetTitle("Distance");  
  fHistDeltaPtvsDistance->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsDistance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsDistance);

  fHistDeltaPtvsCommonEnergy1 = new TH2F("fHistDeltaPtvsCommonEnergy1", "fHistDeltaPtvsCommonEnergy1", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
  fHistDeltaPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsCommonEnergy1);

  fHistDeltaPtvsCommonEnergy2 = new TH2F("fHistDeltaPtvsCommonEnergy2", "fHistDeltaPtvsCommonEnergy2", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
  fHistDeltaPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsCommonEnergy2);

  fHistDeltaPtvsArea1 = new TH2F("fHistDeltaPtvsArea1", "fHistDeltaPtvsArea1", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
  fHistDeltaPtvsArea1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsArea1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsArea1);

  fHistDeltaPtvsArea2 = new TH2F("fHistDeltaPtvsArea2", "fHistDeltaPtvsArea2", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
  fHistDeltaPtvsArea2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsArea2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsArea2);

  fHistDeltaPtvsDeltaArea = new TH2F("fHistDeltaPtvsDeltaArea", "fHistDeltaPtvsDeltaArea", 
				     80, -fJetRadius * fJetRadius * TMath::Pi() * 3, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
  fHistDeltaPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsDeltaArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsDeltaArea);

  fHistJet1PtvsJet2Pt = new TH2F("fHistJet1PtvsJet2Pt", "fHistJet1PtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistJet1PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,1}");
  fHistJet1PtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");
  fHistJet1PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet1PtvsJet2Pt);

  if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {  
    fHistDeltaCorrPtvsJet1Pt = new TH2F("fHistDeltaCorrPtvsJet1Pt", "fHistDeltaCorrPtvsJet1Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
    fHistDeltaCorrPtvsJet1Pt->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet1Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet1Pt);

    fHistDeltaCorrPtvsJet2Pt = new TH2F("fHistDeltaCorrPtvsJet2Pt", "fHistDeltaCorrPtvsJet2Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistDeltaCorrPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet2Pt);

    fHistDeltaCorrPtvsDistance = new TH2F("fHistDeltaCorrPtvsDistance", "fHistDeltaCorrPtvsDistance", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsDistance->GetXaxis()->SetTitle("Distance");  
    fHistDeltaCorrPtvsDistance->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsDistance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsDistance);

    fHistDeltaCorrPtvsCommonEnergy1 = new TH2F("fHistDeltaCorrPtvsCommonEnergy1", "fHistDeltaCorrPtvsCommonEnergy1", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
    fHistDeltaCorrPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsCommonEnergy1);

    fHistDeltaCorrPtvsCommonEnergy2 = new TH2F("fHistDeltaCorrPtvsCommonEnergy2", "fHistDeltaCorrPtvsCommonEnergy2", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
    fHistDeltaCorrPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsCommonEnergy2);

    fHistDeltaCorrPtvsArea1 = new TH2F("fHistDeltaCorrPtvsArea1", "fHistDeltaCorrPtvsArea1", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
    fHistDeltaCorrPtvsArea1->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsArea1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsArea1);
    
    fHistDeltaCorrPtvsArea2 = new TH2F("fHistDeltaCorrPtvsArea2", "fHistDeltaCorrPtvsArea2", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
    fHistDeltaCorrPtvsArea2->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsArea2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsArea2);
    
    fHistDeltaCorrPtvsDeltaArea = new TH2F("fHistDeltaCorrPtvsDeltaArea", "fHistDeltaCorrPtvsDeltaArea", 
					   80, -fJetRadius * fJetRadius * TMath::Pi() * 3, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
    fHistDeltaCorrPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsDeltaArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsDeltaArea);
    
    if (fRhoName.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else if (fRho2Name.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 2*fNbins, -fMaxBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    else
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 2*fNbins, -fMaxBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJet1CorrPtvsJet2CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetYaxis()->SetTitle("p_{T,2}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJet1CorrPtvsJet2CorrPt);
  }

  if (fIsEmbedded) {
    fHistDeltaMCPtvsJet1Pt = new TH2F("fHistDeltaMCPtvsJet1Pt", "fHistDeltaMCPtvsJet1Pt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
    fHistDeltaMCPtvsJet1Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsJet1Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsJet1Pt);

    fHistDeltaMCPtvsJet2Pt = new TH2F("fHistDeltaMCPtvsJet2Pt", "fHistDeltaMCPtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistDeltaMCPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsJet2Pt);

    fHistDeltaMCPtvsDistance = new TH2F("fHistDeltaMCPtvsDistance", "fHistDeltaMCPtvsDistance", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsDistance->GetXaxis()->SetTitle("Distance");  
    fHistDeltaMCPtvsDistance->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsDistance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsDistance);

    fHistDeltaMCPtvsCommonEnergy1 = new TH2F("fHistDeltaMCPtvsCommonEnergy1", "fHistDeltaMCPtvsCommonEnergy1", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
    fHistDeltaMCPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsCommonEnergy1);

    fHistDeltaMCPtvsCommonEnergy2 = new TH2F("fHistDeltaMCPtvsCommonEnergy2", "fHistDeltaMCPtvsCommonEnergy2", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
    fHistDeltaMCPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsCommonEnergy2);

    fHistDeltaMCPtvsArea1 = new TH2F("fHistDeltaMCPtvsArea1", "fHistDeltaMCPtvsArea1", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
    fHistDeltaMCPtvsArea1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsArea1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsArea1);

    fHistDeltaMCPtvsArea2 = new TH2F("fHistDeltaMCPtvsArea2", "fHistDeltaMCPtvsArea2", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
    fHistDeltaMCPtvsArea2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsArea2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsArea2);

    fHistDeltaMCPtvsDeltaArea = new TH2F("fHistDeltaMCPtvsDeltaArea", "fHistDeltaMCPtvsDeltaArea", 
					 80, -fJetRadius * fJetRadius * TMath::Pi() * 3, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
    fHistDeltaMCPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsDeltaArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsDeltaArea);

    fHistJet1MCPtvsJet2Pt = new TH2F("fHistJet1MCPtvsJet2Pt", "fHistJet1MCPtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistJet1MCPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,1}^{MC}");
    fHistJet1MCPtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");
    fHistJet1MCPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJet1MCPtvsJet2Pt);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut)
    return kFALSE;
  if (jet->Area() <= fJetAreaCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::AcceptBiasJet2(AliEmcalJet *jet) const
{ 
  // Accept jet with a bias.

  if (fLeadingHadronType == 0) {
    if (jet->MaxTrackPt() < fPtBiasJet2Track) return kFALSE;
  }
  else if (fLeadingHadronType == 1) {
    if (jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }
  else {
    if (jet->MaxTrackPt() < fPtBiasJet2Track && jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Execute once.

  if (!fJets2Name.IsNull() && !fJets2) {
    fJets2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJets2Name));
    if (!fJets2) {
      AliError(Form("%s: Could not retrieve jets2 %s!", GetName(), fJets2Name.Data()));
      return;
    }
    else if (!fJets2->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJets2Name.Data())); 
      fJets2 = 0;
      return;
    }
  }

  if (!fTracks2Name.IsNull() && !fTracks2) {
    fTracks2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracks2Name));
    if (!fTracks2) {
      AliError(Form("%s: Could not retrieve tracks2 %s!", GetName(), fTracks2Name.Data())); 
      return;
    }
    else {
      TClass *cl = fTracks2->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fTracks2Name.Data())); 
	fTracks2 = 0;
	return;
      }
    }

    if (fAreCollections2MC) {
      fTracks2Map = dynamic_cast<AliNamedArrayI*>(InputEvent()->FindListObject(fTracks2Name + "_Map"));
      // this is needed to map the MC labels with the indexes of the MC particle collection
      // if teh map is not given, the MC labels are assumed to be consistent with the indexes (which is not the case if AliEmcalMCTrackSelector is used)
      if (!fTracks2Map) {
	AliWarning(Form("%s: Could not retrieve map for tracks2 %s! Will assume MC labels consistent with indexes...", GetName(), fTracks2Name.Data())); 
	fTracks2Map = new AliNamedArrayI("tracksMap",9999);
	for (Int_t i = 0; i < 9999; i++) {
	  fTracks2Map->AddAt(i,i);
	}
      }
    }
  }

  if (!fCalo2Name.IsNull() && !fCaloClusters2) {
    fCaloClusters2 =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCalo2Name));
    if (!fCaloClusters2) {
      AliError(Form("%s: Could not retrieve clusters %s!", GetName(), fCalo2Name.Data())); 
      return;
    } else {
      TClass *cl = fCaloClusters2->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fCalo2Name.Data())); 
	fCaloClusters2 = 0;
	return;
      }
    }
  }

  if (!fRho2Name.IsNull() && !fRho2) {
    fRho2 = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRho2Name));
    if (!fRho2) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRho2Name.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (fJet2MinEta == -999)
    fJet2MinEta = fJetMinEta - fJetRadius;
  if (fJet2MaxEta == -999)
    fJet2MaxEta = fJetMaxEta + fJetRadius;
  if (fJet2MinPhi == -999)
    fJet2MinPhi = fJetMinPhi - fJetRadius;
  if (fJet2MaxPhi == -999)
    fJet2MaxPhi = fJetMaxPhi + fJetRadius;

  AliAnalysisTaskEmcalJet::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::IsEventSelected()
{
  // Check if event is selected

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin) 
    return kFALSE;

  return AliAnalysisTaskEmcalJet::IsEventSelected();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  if (fRho2)
    fRho2Val = fRho2->GetVal();

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
  
  if (MCEvent()) {
    fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
    if (!fPythiaHeader) {
      // Check if AOD
      AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

      if (aodMCH) {
        for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
          fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
          if (fPythiaHeader) break;
        }
      }
    }
  }

  if (fPythiaHeader) {
    Double_t pthard = fPythiaHeader->GetPtHard();
    
    for (fPtHardBin = 0; fPtHardBin < 11; fPtHardBin++) {
      if (pthard >= ptHardLo[fPtHardBin] && pthard < ptHardHi[fPtHardBin])
	break;
    }
    
    fNTrials = fPythiaHeader->Trials();
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::Run()
{
  // Find the closest jets

  if (fMatching == kNoMatching) 
    return kTRUE;
  else
    return DoJetMatching();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::DoJetMatching()
{
  DoJetLoop(kFALSE);

  const Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;

    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (jet1->ClosestJet() && jet1->ClosestJet()->ClosestJet() == jet1 && 
        jet1->ClosestJetDistance() < fMatchingPar1 && jet1->ClosestJet()->ClosestJetDistance() < fMatchingPar2) {    // Matched jet found
      jet1->SetMatchedToClosest(fMatching);
      jet1->ClosestJet()->SetMatchedToClosest(fMatching);
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(Bool_t order)
{
  // Do the jet loop.

  TClonesArray *jets1 = 0;
  TClonesArray *jets2 = 0;

  if (order) {
    jets1 = fJets2;
    jets2 = fJets;
  }
  else {
    jets1 = fJets;
    jets2 = fJets2;
  }

  Int_t nJets1 = jets1->GetEntriesFast();
  Int_t nJets2 = jets2->GetEntriesFast();

  for (Int_t j = 0; j < nJets2; j++) {
      
    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
    if (!jet2) {
      AliError(Form("Could not receive jet %d", j));
      continue;
    }

    jet2->ResetMatching();


    if (!AcceptJet(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;
  }
    
  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(jets1->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }

    jet1->ResetMatching();

    if (!AcceptJet(jet1))
      continue;

    if (jet1->MCPt() < 1)
      continue;

    if (order) {
     if (jet1->Eta() < fJet2MinEta || jet1->Eta() > fJet2MaxEta || jet1->Phi() < fJet2MinPhi || jet1->Phi() > fJet2MaxPhi)
	continue;
    }
    else {
      if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
	continue;
    }

    for (Int_t j = 0; j < nJets2; j++) {
      
      AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
      if (!jet2) {
	AliError(Form("Could not receive jet %d", j));
	continue;
      }
      
      if (!AcceptJet(jet2))
	continue;

      if (order) {
	if (jet2->Eta() < fJetMinEta || jet2->Eta() > fJetMaxEta || jet2->Phi() < fJetMinPhi || jet2->Phi() > fJetMaxPhi)
	  continue;
      }
      else {
	if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
	  continue;
      }

      SetMatchingLevel(jet1, jet2, fMatching);
    } // jet2 loop

  } // jet1 loop
}

//________________________________________________________________________
void AliJetResponseMaker::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
  Double_t deta = jet2->Eta() - jet1->Eta();
  Double_t dphi = jet2->Phi() - jet1->Phi();
  d = TMath::Sqrt(deta * deta + dphi * dphi);
}

//________________________________________________________________________
void AliJetResponseMaker::GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();
  Double_t totalPt1 = d1; // the total pt of the reconstructed jet will be cleaned from the background

  for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
    Bool_t track2Found = kFALSE;
    Int_t index2 = jet2->TrackAt(iTrack2);
    for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
      AliVParticle *track = jet1->TrackAt(iTrack,fTracks);
      if (!track) {
	AliWarning(Form("Could not find track %d!", iTrack));
	continue;
      }
      Int_t MClabel = TMath::Abs(track->GetLabel());
      if (MClabel > fMCLabelShift)
	MClabel -= fMCLabelShift;
      Int_t index = -1;
	  
      if (MClabel == 0) {// this is not a MC particle; remove it completely
	AliDebug(3,Form("Track %d (pT = %f) is not a MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	totalPt1 -= track->Pt();
	d1 -= track->Pt();
	continue;
      }
      else if (MClabel < fTracks2Map->GetSize()) {
	index = fTracks2Map->At(MClabel);
      }
	  
      if (index < 0) {
	AliDebug(2,Form("Track %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	continue;
      }

      if (index2 == index) { // found common particle
	track2Found = kTRUE;
	d1 -= track->Pt();
	AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	AliDebug(3,Form("Track %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			iTrack,track->Pt(),track->Eta(),track->Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
	d2 -= MCpart->Pt();
	break;
      }
    }
    for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
      AliVCluster *clus = jet1->ClusterAt(iClus,fCaloClusters);
      if (!clus) {
	AliWarning(Form("Could not find cluster %d!", iClus));
	continue;
      }
      TLorentzVector part;
      clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	  
      if (fUseCellsToMatch && fCaloCells) { // if the cell colection is available, look for cells with a matched MC particle
	for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
	  Int_t cellId = clus->GetCellAbsId(iCell);
	  Double_t cellFrac = clus->GetCellAmplitudeFraction(iCell);

	  Int_t MClabel = TMath::Abs(fCaloCells->GetCellMCLabel(cellId));
	  if (MClabel > fMCLabelShift)
	    MClabel -= fMCLabelShift;
	  Int_t index = -1;
	  
	  if (MClabel == 0) {// this is not a MC particle; remove it completely
	    AliDebug(3,Form("Cell %d (frac = %f) is not a MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    totalPt1 -= part.Pt() * cellFrac;
	    d1 -= part.Pt() * cellFrac;
	    continue;
	  }
	  else if (MClabel < fTracks2Map->GetSize()) {
	    index = fTracks2Map->At(MClabel);
	  }

	  if (index < 0) {
	    AliDebug(3,Form("Cell %d (frac = %f) does not have an associated MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    continue;
	  }
	  if (index2 == index) { // found common particle
	    d1 -= part.Pt() * cellFrac;
		
	    if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	      AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	      AliDebug(3,Form("Cell %d belonging to cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			      iCell,iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));		  
	      d2 -= MCpart->Pt() * cellFrac;
	    }
	    break;
	  }
	}
      }
      else { //otherwise look for the first contributor to the cluster, and if matched to a MC label remove it
	Int_t MClabel = TMath::Abs(clus->GetLabel());
	if (MClabel > fMCLabelShift)
	  MClabel -= fMCLabelShift;
	Int_t index = -1;
	    
	if (MClabel == 0) {// this is not a MC particle; remove it completely
	  AliDebug(3,Form("Cluster %d (pT = %f) is not a MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  totalPt1 -= part.Pt();
	  d1 -= part.Pt();
	  continue;
	}
	else if (MClabel < fTracks2Map->GetSize()) {
	  index = fTracks2Map->At(MClabel);
	}
	 
	if (index < 0) {
	  AliDebug(3,Form("Cluster %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  continue;
	}
	if (index2 == index) { // found common particle
	  d1 -= part.Pt();

	  if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	    AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	    AliDebug(3,Form("Cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			    iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
		
	    d2 -= MCpart->Pt();
	  }
	  break;
	}
      }
    }
  }

  if (d1 < 0)
    d1 = 0;

  if (d2 < 0)
    d2 = 0;

  if (totalPt1 < 1)
    d1 = -1;
  else
    d1 /= totalPt1;

  if (jet2->Pt() < 1)
    d2 = -1;
  else
    d2 /= jet2->Pt();
}

//________________________________________________________________________
void AliJetResponseMaker::GetSameCollectionsMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();

  if (fTracks && fTracks2) {

    for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
      Int_t index2 = jet2->TrackAt(iTrack2);
      for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
	Int_t index = jet1->TrackAt(iTrack);
	if (index2 == index) { // found common particle
	  AliVParticle *part = static_cast<AliVParticle*>(fTracks->At(index));
	  if (!part) {
	    AliWarning(Form("Could not find track %d!", index));
	    continue;
	  }
	  AliVParticle *part2 = static_cast<AliVParticle*>(fTracks2->At(index2));
	  if (!part2) {
	    AliWarning(Form("Could not find track %d!", index2));
	    continue;
	  }

	  d1 -= part->Pt();
	  d2 -= part2->Pt();
	  break;
	}
      }
    }

  }

  if (fCaloClusters && fCaloClusters2) {

    if (fUseCellsToMatch) {
      const Int_t nClus1 = jet1->GetNumberOfClusters();

      Int_t ncells1[nClus1];
      UShort_t *cellsId1[nClus1];
      Double_t *cellsFrac1[nClus1];
      Int_t *sortedIndexes1[nClus1];
      Double_t ptClus1[nClus1];
      for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
	Int_t index1 = jet1->ClusterAt(iClus1);
	AliVCluster *clus1 =  static_cast<AliVCluster*>(fCaloClusters->At(index1));
	if (!clus1) {
	  AliWarning(Form("Could not find cluster %d!", index1));
	  continue;
	}
	ncells1[iClus1] = clus1->GetNCells();
	cellsId1[iClus1] = clus1->GetCellsAbsId();
	cellsFrac1[iClus1] = clus1->GetCellsAmplitudeFraction();
	sortedIndexes1[iClus1] = new Int_t[ncells1[iClus1]];

	TLorentzVector part1;
	clus1->GetMomentum(part1, const_cast<Double_t*>(fVertex));
	ptClus1[iClus1] = part1.Pt();

	TMath::Sort(ncells1[iClus1], cellsId1[iClus1], sortedIndexes1[iClus1], kFALSE);
      }
      
      const Int_t nClus2 = jet2->GetNumberOfClusters();

      const Int_t maxNcells2 = 11520;
      Int_t sortedIndexes2[maxNcells2];
      for (Int_t iClus2 = 0; iClus2 < nClus2; iClus2++) {
	Int_t index2 = jet2->ClusterAt(iClus2);
	AliVCluster *clus2 =  static_cast<AliVCluster*>(fCaloClusters2->At(index2));
	if (!clus2) {
	  AliWarning(Form("Could not find cluster %d!", index2));
	  continue;
	}
	Int_t ncells2 = clus2->GetNCells();
	if (ncells2 > maxNcells2) {
	  AliError(Form("Number of cells in the cluster %d > %d",ncells2,maxNcells2));
	  continue;
	}
	UShort_t *cellsId2 = clus2->GetCellsAbsId();
	Double_t *cellsFrac2 = clus2->GetCellsAmplitudeFraction();

	TLorentzVector part2;
	clus2->GetMomentum(part2, const_cast<Double_t*>(fVertex));
	Double_t ptClus2 = part2.Pt();

	TMath::Sort(ncells2, cellsId2, sortedIndexes2, kFALSE);

	for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
	  
	  Int_t iCell1 = 0, iCell2 = 0;
	  Bool_t common=kFALSE;
	  while (iCell1 < ncells1[iClus1] && iCell2 < ncells2) {
	    if (cellsId1[iClus1][sortedIndexes1[iClus1][iCell1]] == cellsId2[sortedIndexes2[iCell2]]) { // found a common cell
	      d1 -= cellsFrac1[iClus1][sortedIndexes1[iClus1][iCell1]] * ptClus1[iClus1];
	      d2 -= cellsFrac2[sortedIndexes2[iCell2]] * ptClus2;
	      iCell1++;
	      iCell2++;
	      common = kTRUE;
	    }
	    else if (cellsId1[iClus1][sortedIndexes1[iClus1][iCell1]] > cellsId2[sortedIndexes2[iCell2]]) { 
	      iCell2++;
	    }
	    else {
	      iCell1++;
	    }
	  }
	}
      }
    }
    else {
      for (Int_t iClus2 = 0; iClus2 < jet2->GetNumberOfClusters(); iClus2++) {
	Int_t index2 = jet2->ClusterAt(iClus2);
	for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
	  Int_t index = jet1->ClusterAt(iClus);
	  if (index2 == index) { // found common particle
	    AliVCluster *clus =  static_cast<AliVCluster*>(fCaloClusters->At(index));
	    if (!clus) {
	      AliWarning(Form("Could not find cluster %d!", index));
	      continue;
	    }
	    AliVCluster *clus2 =  static_cast<AliVCluster*>(fCaloClusters2->At(index2));
	    if (!clus2) {
	      AliWarning(Form("Could not find cluster %d!", index2));
	      continue;
	    }
	    TLorentzVector part, part2;
	    clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	    clus2->GetMomentum(part2, const_cast<Double_t*>(fVertex));
	    
	    d1 -= part.Pt();
	    d2 -= part2.Pt();
	    break;
	  }
	}
      }
    }
  }

  if (d1 < 0)
    d1 = 0;

  if (d2 < 0)
    d2 = 0;

  if (jet1->Pt() > 0)
    d1 /= jet1->Pt();
  else
    d1 = -1;

  if (jet2->Pt() > 0)
    d2 /= jet2->Pt();
  else
    d2 = -1;
}

//________________________________________________________________________
void AliJetResponseMaker::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching) 
{
  Double_t d1 = -1;
  Double_t d2 = -1;

  switch (matching) {
  case kGeometrical:
    GetGeometricalMatchingLevel(jet1,jet2,d1);
    d2 = d1;
    break;
  case kMCLabel: // jet1 = detector level and jet2 = particle level!
    GetMCLabelMatchingLevel(jet1,jet2,d1,d2);
    break;
  case kSameCollections:
    GetSameCollectionsMatchingLevel(jet1,jet2,d1,d2);
    break;
  default:
    ;
  }

  if (d1 >= 0) {

    if (d1 < jet1->ClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
      jet1->SetClosestJet(jet2, d1);
    }
    else if (d1 < jet1->SecondClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet2, d1);
    }
  }
  
  if (d2 >= 0) {
    
    if (d2 < jet2->ClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
      jet2->SetClosestJet(jet1, d2);
    }
    else if (d2 < jet2->SecondClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet1, d2);
    }
  }
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  static Int_t indexes[9999] = {-1};

  if (fHistEventsAfterSel)
    fHistEventsAfterSel->SetBinContent(fPtHardBin + 1, fHistEventsAfterSel->GetBinContent(fPtHardBin + 1) + 1);
  if (fHistTrialsAfterSel)
    fHistTrialsAfterSel->SetBinContent(fPtHardBin + 1, fHistTrialsAfterSel->GetBinContent(fPtHardBin + 1) + fNTrials);

  GetSortedArray(indexes, fJets2, fRho2Val);

  const Int_t nJets2 = fJets2->GetEntriesFast();

  Int_t naccJets2 = 0;
  Int_t naccJets2Acceptance = 0;

  Double_t totalPt2 = 0;

  for (Int_t i = 0; i < nJets2; i++) {

    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(fJets2->At(indexes[i]));

    if (!jet2) {
      AliError(Form("Could not receive jet2 %d", i));
      continue;
    }

    if (!AcceptJet(jet2))
      continue;

    if (AcceptBiasJet(jet2) &&
	(jet2->Eta() > fJetMinEta && jet2->Eta() < fJetMaxEta && jet2->Phi() > fJetMinPhi && jet2->Phi() < fJetMaxPhi)) {
      
      fHistJets2PtAreaAcceptance->Fill(jet2->Area(), jet2->Pt());
      fHistJets2PhiEtaAcceptance->Fill(jet2->Eta(), jet2->Phi());

      totalPt2 += jet2->Pt();
      
      if (naccJets2Acceptance < fNLeadingJets)
	fHistLeadingJets2PtAreaAcceptance->Fill(jet2->Area(), jet2->Pt());
      
      if (!fRho2Name.IsNull()) {
	fHistJets2CorrPtAreaAcceptance->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
	if (naccJets2Acceptance < fNLeadingJets)
	  fHistLeadingJets2CorrPtAreaAcceptance->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
      }

      if (fTracks2) {
	for (Int_t it = 0; it < jet2->GetNumberOfTracks(); it++) {
	  AliVParticle *track2 = jet2->TrackAt(it, fTracks2);
	  if (track2) 
	    fHistJets2ZvsPt->Fill(track2->Pt() / jet2->Pt(), jet2->Pt());
	}
      }

      if (fCaloClusters2) {
	for (Int_t ic = 0; ic < jet2->GetNumberOfClusters(); ic++) {
	  AliVCluster *cluster2 = jet2->ClusterAt(ic, fCaloClusters2);
	  
	  if (cluster2) {
	    TLorentzVector nPart2;
	    cluster2->GetMomentum(nPart2, fVertex);
	    fHistJets2ZvsPt->Fill(nPart2.Et() / jet2->Pt(), jet2->Pt());
	  }
	}
      }

      fHistJets2NEFvsPt->Fill(jet2->NEF(), jet2->Pt());
      fHistJets2CEFvsCEFPt->Fill(1-jet2->NEF(), (1-jet2->NEF())*jet2->Pt());
      
      naccJets2Acceptance++;
    }

    if (!AcceptBiasJet2(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;
    
    fHistJets2PtArea->Fill(jet2->Area(), jet2->Pt());
    fHistJets2PhiEta->Fill(jet2->Eta(), jet2->Phi());
    
    if (naccJets2 < fNLeadingJets)
      fHistLeadingJets2PtArea->Fill(jet2->Area(), jet2->Pt());

    if (!fRho2Name.IsNull()) {
      fHistJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
      if (naccJets2 < fNLeadingJets)
	fHistLeadingJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
    }

    naccJets2++;

    if (jet2->MatchedJet()) {

      if (!AcceptJet(jet2->MatchedJet()) ||
	  jet2->MatchedJet()->Eta() < fJetMinEta || jet2->MatchedJet()->Eta() > fJetMaxEta || 
	  jet2->MatchedJet()->Phi() < fJetMinPhi || jet2->MatchedJet()->Phi() > fJetMaxPhi ||
	  !AcceptBiasJet(jet2->MatchedJet()) || 
	  jet2->MatchedJet()->MaxTrackPt() > fMaxTrackPt || jet2->MatchedJet()->MaxClusterPt() > fMaxClusterPt) {
	fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      }
      else {
	if (jet2->MatchedJet()->Pt() > fMaxBinPt)
	  fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());

	Double_t d=-1, ce1=-1, ce2=-1;
	if (jet2->GetMatchingType() == kGeometrical) {
	  if (fAreCollections2MC && !fAreCollections1MC) // the other way around is not supported
	    GetMCLabelMatchingLevel(jet2->MatchedJet(), jet2, ce1, ce2);
	  else if ((fAreCollections1MC && fAreCollections2MC) || (!fAreCollections1MC && !fAreCollections2MC))
	    GetSameCollectionsMatchingLevel(jet2->MatchedJet(), jet2, ce1, ce2);

	  d = jet2->ClosestJetDistance();
	}
	else if (jet2->GetMatchingType() == kMCLabel || jet2->GetMatchingType() == kSameCollections) {
	  GetGeometricalMatchingLevel(jet2->MatchedJet(), jet2, d);

	  ce1 = jet2->MatchedJet()->ClosestJetDistance();
	  ce2 = jet2->ClosestJetDistance();
	}

	fHistCommonEnergy1vsJet1Pt->Fill(ce1, jet2->MatchedJet()->Pt());
	fHistCommonEnergy2vsJet2Pt->Fill(ce2, jet2->Pt());

	fHistDistancevsJet1Pt->Fill(d, jet2->MatchedJet()->Pt());
	fHistDistancevsJet2Pt->Fill(d, jet2->Pt());

	fHistDistancevsCommonEnergy1->Fill(d, ce1);
	fHistDistancevsCommonEnergy2->Fill(d, ce2);
	fHistCommonEnergy1vsCommonEnergy2->Fill(ce1, ce2);

	fHistJet2PtOverJet1PtvsJet2Pt->Fill(jet2->Pt(), jet2->Pt() / jet2->MatchedJet()->Pt());
	fHistJet1PtOverJet2PtvsJet1Pt->Fill(jet2->MatchedJet()->Pt(), jet2->MatchedJet()->Pt() / jet2->Pt());

	Double_t deta = jet2->MatchedJet()->Eta() - jet2->Eta();
	Double_t dphi = jet2->MatchedJet()->Phi() - jet2->Phi();
	fHistDeltaEtaPhi->Fill(deta, dphi, jet2->Pt());

	Double_t dpt = jet2->MatchedJet()->Pt() - jet2->Pt();
	fHistDeltaPtvsJet1Pt->Fill(jet2->MatchedJet()->Pt(), dpt);
	fHistDeltaPtvsJet2Pt->Fill(jet2->Pt(), dpt);

	fHistDeltaPtvsDistance->Fill(d, dpt);
	fHistDeltaPtvsCommonEnergy1->Fill(ce1, dpt);
	fHistDeltaPtvsCommonEnergy2->Fill(ce2, dpt);

	fHistDeltaPtvsArea1->Fill(jet2->MatchedJet()->Area(), dpt);
	fHistDeltaPtvsArea2->Fill(jet2->Area(), dpt);

	Double_t darea = jet2->MatchedJet()->Area() - jet2->Area();
	fHistDeltaPtvsDeltaArea->Fill(darea, dpt);

	fHistJet1PtvsJet2Pt->Fill(jet2->MatchedJet()->Pt(), jet2->Pt());

	if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
	  Double_t corrpt1 = jet2->MatchedJet()->Pt() - fRhoVal * jet2->MatchedJet()->Area();
	  Double_t corrpt2 = jet2->Pt() - fRho2Val * jet2->Area();
	  Double_t dcorrpt = corrpt1 - corrpt2;
	  fHistDeltaCorrPtvsJet1Pt->Fill(jet2->MatchedJet()->Pt(), dcorrpt);
	  fHistDeltaCorrPtvsJet2Pt->Fill(jet2->Pt(), dcorrpt);
	  fHistDeltaCorrPtvsDistance->Fill(d, dcorrpt);
	  fHistDeltaCorrPtvsCommonEnergy1->Fill(ce1, dcorrpt);
	  fHistDeltaCorrPtvsCommonEnergy2->Fill(ce2, dcorrpt);
	  fHistDeltaCorrPtvsArea1->Fill(jet2->MatchedJet()->Area(), dcorrpt);
	  fHistDeltaCorrPtvsArea2->Fill(jet2->Area(), dcorrpt);
	  fHistDeltaCorrPtvsDeltaArea->Fill(darea, dcorrpt);
	  fHistJet1CorrPtvsJet2CorrPt->Fill(corrpt1, corrpt2);
	}

	if (fIsEmbedded) {
	  Double_t dmcpt = jet2->MatchedJet()->MCPt() - jet2->Pt();
	  fHistDeltaMCPtvsJet1Pt->Fill(jet2->MatchedJet()->MCPt(), dmcpt);
	  fHistDeltaMCPtvsJet2Pt->Fill(jet2->Pt(), dmcpt);
	  fHistDeltaMCPtvsDistance->Fill(d, dmcpt);
	  fHistDeltaMCPtvsCommonEnergy1->Fill(ce1, dmcpt);
	  fHistDeltaMCPtvsCommonEnergy2->Fill(ce2, dmcpt);
	  fHistDeltaMCPtvsArea1->Fill(jet2->MatchedJet()->Area(), dmcpt);
	  fHistDeltaMCPtvsArea2->Fill(jet2->Area(), dmcpt);
	  fHistDeltaMCPtvsDeltaArea->Fill(darea, dmcpt);
	  fHistJet1MCPtvsJet2Pt->Fill(jet2->MatchedJet()->MCPt(), jet2->Pt());
	}
      }
    }
    else {
      fHistNonMatchedJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());

      if (!fRho2Name.IsNull())
	fHistNonMatchedJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRhoVal * jet2->Area());
    }
  }

  GetSortedArray(indexes, fJets, fRhoVal);

  const Int_t nJets1 = fJets->GetEntriesFast();
  Int_t naccJets1 = 0;
  Double_t totalMCPt1 = 0;

  for (Int_t i = 0; i < nJets1; i++) {

    AliDebug(2,Form("Processing jet %d", i));
    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(indexes[i]));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  
    
    if (!AcceptJet(jet1))
      continue;

    if (!AcceptBiasJet(jet1))
      continue;

    if (jet1->MaxTrackPt() > fMaxTrackPt || jet1->MaxClusterPt() > fMaxClusterPt)
      continue;

    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (!jet1->MatchedJet()) {
      fHistNonMatchedJets1PtArea->Fill(jet1->Area(), jet1->Pt());
      if (!fRhoName.IsNull())
	fHistNonMatchedJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());
    }

    fHistJets1PtArea->Fill(jet1->Area(), jet1->Pt());
    fHistJets1PhiEta->Fill(jet1->Eta(), jet1->Phi());

    totalMCPt1 += jet1->MCPt();

    if (naccJets1 < fNLeadingJets)
      fHistLeadingJets1PtArea->Fill(jet1->Area(), jet1->Pt());

    if (!fRhoName.IsNull()) {
      fHistJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());

      if (naccJets1 < fNLeadingJets)
	fHistLeadingJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());
    }

    if (fTracks) {
      for (Int_t it = 0; it < jet1->GetNumberOfTracks(); it++) {
	AliVParticle *track1 = jet1->TrackAt(it, fTracks2);
	if (track1) 
	  fHistJets1ZvsPt->Fill(track1->Pt() / jet1->Pt(), jet1->Pt());
      }
    }
    
    if (fCaloClusters) {
      for (Int_t ic = 0; ic < jet1->GetNumberOfClusters(); ic++) {
	AliVCluster *cluster1 = jet1->ClusterAt(ic, fCaloClusters);
	
	if (cluster1) {
	  TLorentzVector nPart1;
	  cluster1->GetMomentum(nPart1, fVertex);
	  fHistJets2ZvsPt->Fill(nPart1.Et() / jet1->Pt(), jet1->Pt());
	}
      }
    }
    
    fHistJets1NEFvsPt->Fill(jet1->NEF(), jet1->Pt());
    fHistJets1CEFvsCEFPt->Fill(1-jet1->NEF(), (1-jet1->NEF())*jet1->Pt());

    naccJets1++;
  }

  if (fIsEmbedded) 
    fMCEnergy1vsEnergy2->Fill(totalMCPt1, totalPt2);

  return kTRUE;
}
