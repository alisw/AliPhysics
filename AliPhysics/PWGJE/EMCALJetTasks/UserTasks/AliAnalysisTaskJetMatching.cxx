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

// 
// General task to match two sets of jets
//
// This task takes two TClonesArray's as input: 
// [1] fSourceJets - e.g. pythia jets
// [2] fTargetJets - e.g. a samle containing pythia jets embedded in a pbpb event
// the task will try to match jets from the source array to the target array, and
// save the found TARGET jets in a new array 'fMatchedJets'
// the purpose of this task is twofold
// [1] matching for embedding studies
// [2] matching to create a detector response function
//
// matching is done in three steps
// [1] geometric matching, where jets are matched by R = sqrt(dphi*dphi-deta*deta) or directly via dphi and deta
// [2] optional injection / bijection check 
//     in this check, it is made sure that fSourceJets -> fMatchedJets is either an injective non-surjective 
//     or bijective map, depending on the matching resolution which is chosen for step [1]
//     so that each source jet has a unique target and vice-versa.
//     if R (or dphi, deta) are proportional to, or larger than, the jet radius, matching tends to be 
//     bijective (each source has a target), if R is chosen to be very small, source jets might be lost 
//     as no target can be found.
//     the mapping is done in such a way that each target is matched to its closest source and each source
//     is mapped to its closest target
// [3] constituent matching
//     - how many constituents of the source jet are present in the matched jet? 
//       a cut on the constituent fraction can be performed (not recommended)
//     - how much of the original pt is recovered in the matched jet? 
//       a cut on the fraction of recovered / original pt can be performed (recommended)
//
// detector response
//     to obtain a detector respose function, supply
// [1] fSourceJets - particle level jets
// [2] fTargetJets - detector level jets
// 
// Author: Redmer Alexander Bertens, Utrecht Univeristy, Utrecht, Netherlands
// rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 

// root includes
#include <TClonesArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TSystem.h>
// aliroot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
// emcal jet framework includes
#include <AliEmcalJet.h>
#include <AliAnalysisTaskJetMatching.h>
#include <AliLocalRhoParameter.h>

class AliAnalysisTaskJetMatching;
using namespace std;

ClassImp(AliAnalysisTaskJetMatching)

AliAnalysisTaskJetMatching::AliAnalysisTaskJetMatching() : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetMatching", kTRUE), 
    fSourceJets(0), fSourceJetsName(0), fTargetJets(0), fTargetJetsName(0), fMatchedJets(0), fMatchedJetsName(GetName()), fSourceRho(0), fSourceRhoName(0), fTargetRho(0), fTargetRhoName(0), fUseScaledRho(0), fSourceRadius(0.3), fTargetRadius(0.3), fMatchingScheme(kGeoEtaPhi), fUseEmcalBaseJetCuts(kFALSE), fSourceBKG(kNoSourceBKG), fTargetBKG(kNoTargetBKG), fOutputList(0), fHistUnsortedCorrelation(0), fHistMatchedCorrelation(0), fHistSourceJetPt(0), fHistMatchedSourceJetPt(0), fHistTargetJetPt(0), fHistMatchedJetPt(0), fHistSourceMatchedJetPt(0), fHistDetectorResponseProb(0), fHistNoConstSourceJet(0), fHistNoConstTargetJet(0), fHistNoConstMatchJet(0), fProfFracPtMatched(0), fProfFracPtJets(0), fProfFracNoMatched(0), fProfFracNoJets(0), fHistDiJet(0), fHistDiJetLeadingJet(0), fHistDiJetDPhi(0), fHistDiJetDPt(0), fHistAnalysisSummary(0), fProfQAMatched(0), fProfQA(0), fNoMatchedJets(200), fMatchEta(.3), fMatchPhi(.3), fMatchR(.08), fDoDetectorResponse(kFALSE), fMatchConstituents(kTRUE), fMinFracRecoveredConstituents(.5), fMinFracRecoveredConstituentPt(0.5), fGetBijection(kTRUE), fh1Trials(0x0), fh1AvgTrials(0x0), fh1Xsec(0x0), fAvgTrials(0) {
    // default constructor
    ClearMatchedJetsCache();
}
//_____________________________________________________________________________
AliAnalysisTaskJetMatching::AliAnalysisTaskJetMatching(const char* name) : AliAnalysisTaskEmcalJet(name, kTRUE),
    fSourceJets(0), fSourceJetsName(0), fTargetJets(0), fTargetJetsName(0), fMatchedJets(0), fMatchedJetsName(GetName()), fSourceRho(0), fSourceRhoName(0), fTargetRho(0), fTargetRhoName(0), fUseScaledRho(0), fSourceRadius(0.3), fTargetRadius(0.3), fMatchingScheme(kGeoEtaPhi), fUseEmcalBaseJetCuts(kFALSE), fSourceBKG(kNoSourceBKG), fTargetBKG(kNoTargetBKG), fOutputList(0), fHistUnsortedCorrelation(0), fHistMatchedCorrelation(0), fHistSourceJetPt(0), fHistMatchedSourceJetPt(0), fHistTargetJetPt(0), fHistMatchedJetPt(0), fHistSourceMatchedJetPt(0), fHistDetectorResponseProb(0), fHistNoConstSourceJet(0), fHistNoConstTargetJet(0), fHistNoConstMatchJet(0), fProfFracPtMatched(0), fProfFracPtJets(0), fProfFracNoMatched(0), fProfFracNoJets(0), fHistDiJet(0), fHistDiJetLeadingJet(0), fHistDiJetDPhi(0), fHistDiJetDPt(0), fHistAnalysisSummary(0), fProfQAMatched(0), fProfQA(0), fNoMatchedJets(200), fMatchEta(.3), fMatchPhi(.3), fMatchR(.08), fDoDetectorResponse(kFALSE), fMatchConstituents(kTRUE), fMinFracRecoveredConstituents(0.5), fMinFracRecoveredConstituentPt(0.5), fGetBijection(kTRUE), fh1Trials(0x0), fh1AvgTrials(0x0), fh1Xsec(0x0), fAvgTrials(0) {
    // constructor
    ClearMatchedJetsCache();
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskJetMatching::~AliAnalysisTaskJetMatching()
{
    // destructor
    if(fOutputList)             delete fOutputList;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::ExecOnce() 
{
    // initialize the anaysis
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    // get the stand alone jets from the input event
    fSourceJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fSourceJetsName.Data()));
    if(!fSourceJets) AliFatal(Form("%s: Container with name %s not found. Aborting", GetName(), fSourceJetsName.Data()));
    fTargetJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTargetJetsName.Data()));
    if(!fTargetJets) AliFatal(Form("%s: Container with name %s not found. Aborting", GetName(), fSourceJetsName.Data()));
    // append the list of matched jets to the event
    fMatchedJets->Delete();
    if (!(InputEvent()->FindListObject(fMatchedJetsName))) InputEvent()->AddObject(fMatchedJets);
    else AliFatal(Form("%s: Object with name %s already in event! Aborting", GetName(), fMatchedJetsName.Data()));
    FillAnalysisSummaryHistogram();
    switch (fSourceBKG) {
        case kSourceLocalRho : {
            fSourceRho = dynamic_cast<AliLocalRhoParameter*>(InputEvent()->FindListObject(fSourceRhoName));
            if(!fSourceRho) AliFatal(Form("%s: Object with name %s requested but not found! Aborting", GetName(), fSourceRhoName.Data()));
        } break;
        default : break;
    }
    switch (fTargetBKG) {
        case kTargetLocalRho : {
            fTargetRho = dynamic_cast<AliLocalRhoParameter*>(InputEvent()->FindListObject(fTargetRhoName));
            if(!fTargetRho) AliFatal(Form("%s: Object with name %s requested but not found! Aborting", GetName(), fTargetRhoName.Data()));
        } break;
        default : break;
    }
    AliAnalysisTaskEmcalJet::ExecOnce(); // init base class
    if(fDoDetectorResponse) fMatchConstituents = kFALSE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::UserCreateOutputObjects()
{
    // create output objects
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    // add the matched jets array to the event
    fMatchedJets = new TClonesArray("AliEmcalJet");
    fMatchedJets->SetName(fMatchedJetsName);
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    // add analysis histograms
    fHistUnsortedCorrelation = BookTH1F("fHistUnsortedCorrelation", "#Delta #varphi unsorted", 50, 0, TMath::Pi());
    fHistMatchedCorrelation = BookTH1F("fHistMatchedCorrelation", "#Delta #varphi matched", 50, 0, TMath::Pi());
    fHistSourceJetPt = (fDoDetectorResponse) ? BookTH1F("fHistParticleLevelJetPt", "p_{t}^{gen} [GeV/c]", 150, -50, 250) :  BookTH1F("fHistSourceJetPt", "p_{t} [GeV/c]", 150, -50, 250);      
    fHistMatchedSourceJetPt = (fDoDetectorResponse) ? BookTH1F("fHistMatchedParticleLevelJetPt", "p_{t}^{gen} [GeV/c]", 150, -50, 250) : BookTH1F("fHistMatchedSourceJetPt", "p_{t} [GeV/c]", 150, -50, 250);
    fHistTargetJetPt = BookTH1F("fHistTargetJetPt", "p_{t} [GeV/c]", 150, -50, 250);
    fHistMatchedJetPt = (fDoDetectorResponse) ? BookTH1F("fHistDetectorLevelJet", "p_{t}^{rec}", 150, -50, 250) : BookTH1F("fHistMatchedJetPt", "p_{t} [GeV/c]", 150, -50, 250);
    fHistSourceMatchedJetPt = (fDoDetectorResponse) ? BookTH2F("fHistDetectorResponse", "particle level jet p_{t}^{gen} [GeV/c]", "detector level jet p_{t}^{rec} [GeV/c]", 300, -50, 250, 300, -50, 250) : BookTH2F("fHistSourceMatchedJetPt", "source jet p_{t} [GeV/c]", "matched jet p_{t} [GeV/c]", 300, -50, 250, 300, -50, 250);
    if(fDoDetectorResponse) {
        fHistDetectorResponseProb = BookTH2F("fHistDetectorResponseProb", "(p_{t}^{det} - p_{t}^{part})/p_{t}^{part}", "p_{t}^{part}", 100, -1.5, 1., 20, 0, 200);
        fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
        fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
        fOutputList->Add(fh1Xsec);
        fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
        fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
        fOutputList->Add(fh1Trials);
        fh1AvgTrials = new TH1F("fh1AvgTrials","avg trials root file",1,0,1);
        fh1AvgTrials->GetXaxis()->SetBinLabel(1,"#sum{avg ntrials}");
        fOutputList->Add(fh1AvgTrials);
    }
    fHistNoConstSourceJet = BookTH1F("fHistNoConstSourceJet", "number of constituents source jets", 50, 0, 50);
    fHistNoConstTargetJet = BookTH1F("fHistNoConstTargetJet", "number of constituents target jets", 50, 0, 50);
    fHistNoConstMatchJet = BookTH1F("fHistNoConstMatchJet", "number of constituents matched jets", 50, 0, 50);
    if(!fDoDetectorResponse) { // these observables cannot be measured in current detector response mode
        fProfFracPtMatched = new TProfile("fProfFracPtMatched", "recovered target p_{T} / source p_{T}", 15, -50, 250);
        fOutputList->Add(fProfFracPtMatched);
        fProfFracNoMatched = new TProfile("fProfFracNoMatched", "recovered target constituents / source constituents", 15, -50, 250);
        fOutputList->Add(fProfFracNoMatched);
    }
    // the analysis summary histo which stores all the analysis flags is always written to file
    fHistAnalysisSummary = BookTH1F("fHistAnalysisSummary", "flag", 51, -0.5, 15.5);
    fProfQAMatched = new TProfile("fProfQAMatched", "fProfQAMatched", 3, 0, 3);
    fProfQAMatched->GetXaxis()->SetBinLabel(1, "<#delta p_{t} >");
    fProfQAMatched->GetXaxis()->SetBinLabel(2, "<#delta #eta>");
    fProfQAMatched->GetXaxis()->SetBinLabel(3, "<#delta #varphi>");
    fOutputList->Add(fProfQAMatched);
    fProfQA = new TProfile("fProfQA", "fProfQA", 3, 0, 3);
    fProfQA->GetXaxis()->SetBinLabel(1, "<#delta p_{t} >");
    fProfQA->GetXaxis()->SetBinLabel(2, "<#delta #eta>");
    fProfQA->GetXaxis()->SetBinLabel(3, "<#delta #varphi>");
    fOutputList->Add(fProfQA);
    fProfFracPtJets = new TProfile("fProfFracPtJets", "source p_{T} / target p_{T}", 15, -50, 250);
    fOutputList->Add(fProfFracPtJets);
    fProfFracNoJets = new TProfile("fProfFracNoJets", "source constituents / target constituents", 15, -50, 250);
    fOutputList->Add(fProfFracNoJets);
    switch (fMatchingScheme) {
        case kDiJet : {
            fHistDiJet = BookTH3F("fHistDiJet", "matched di-jet #varphi", "matched di-jet #eta", "leading jet p_{T} (GeV/c)", 100, 0., TMath::TwoPi(), 100, -5., 5., 100, 0, 200);
            fHistDiJetLeadingJet = BookTH3F("fHistDiJetLeadingJet", "leading jet #varphi", "leadingd jet #eta", "leading jet p_{T} (GeV/c)", 100, 0., TMath::TwoPi(), 100, -5., 5., 100, 0, 200);
            fHistDiJetDPhi = BookTH2F("fHistDiJetDPhi", "leading jet #varphi - (matched jet #varphi - #pi)", "leading jet p_{T} (GeV/c)", 100, -1.*TMath::Pi(), TMath::Pi(), 100, 0, 200);
            fHistDiJetDPt = BookTH2F("fHistDiJetDPt", "leading jet p_{T} - sub leading jet p_{T} (GeV/c)", "leading jet p_{T} (GeV/c)", 100, -25, 25, 100, 0, 200);
        } break;
        default : break;
    }
    fOutputList->Sort();
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetMatching::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Bool_t append)
{
    // book a TH1F and connect it to the output container
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    if(append && !fOutputList) return 0x0;
    TString title(name);
    title += Form(";%s;[counts]", x);
    TH1F* histogram = new TH1F(name, title.Data(), bins, min, max);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskJetMatching::BookTH2F(const char* name, const char* x, const char*y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Bool_t append)
{
    // book a TH2F and connect it to the output container
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    if(append && !fOutputList) return 0x0;
    TString title(name);
    title += Form(";%s;%s", x, y);
    TH2F* histogram = new TH2F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
TH3F* AliAnalysisTaskJetMatching::BookTH3F(const char* name, const char* x, const char* y, const char* z, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t binsz, Double_t minz, Double_t maxz, Bool_t append)
{
    // book a TH2F and connect it to the output container
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    if(append && !fOutputList) return 0x0;
    TString title(name);
    title += Form(";%s;%s;%s", x, y, z);
    TH3F* histogram = new TH3F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy, binsz, minz, maxz);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetMatching::Notify()
{
    // for each file get the number of trials and pythia cross section
    // see  $ALICE_PHYSICS/PWGJE/AliAnalysisTaskJetProperties.cxx
    //      $ALICE_PHYSICS/PWG/Tools/AliAnalysisHelperJetTasks.cxx
    // this function is implenented here temporarily to avoid introducing a dependency
    // later on this could just be a call to a static helper function
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    Float_t xsection(0), ftrials(1);
    fAvgTrials = ftrials;
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    if(tree) {
        TFile *curfile = tree->GetCurrentFile();
        if (!curfile) return kFALSE;
        TString file(curfile->GetName()); 
        if(file.Contains("root_archive.zip#")) file.Replace(file.Index("#",1,file.Index("root_archive",12,0,TString::kExact),TString::kExact)+1,file.Index(".root",5,TString::kExact)-file.Index("root_archive",12,0,TString::kExact),"");
        else file.ReplaceAll(gSystem->BaseName(file.Data()),"");
        TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root"));
        if(!fxsec) {
            fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
            if(fxsec) {
                TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
                if(key) {
                      TList *list = dynamic_cast<TList*>(key->ReadObj());
                      if(list) {
                          xsection = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
                          ftrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
                      }
                }
                fxsec->Close();
            }
        } else {
            TTree *xtree = (TTree*)fxsec->Get("Xsection");
            if(xtree) {
                UInt_t _ftrials  = 0;
                Double_t _xsection  = 0;
                xtree->SetBranchAddress("xsection",&_xsection);
                xtree->SetBranchAddress("ntrials",&_ftrials);
                xtree->GetEntry(0);
                ftrials = _ftrials;
                xsection = _xsection;
            }
            fxsec->Close();
        }
        fh1Xsec->Fill("<#sigma>",xsection);
        Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
        if(ftrials >= nEntries && nEntries > 0.) fAvgTrials = ftrials/nEntries;
        fh1Trials->Fill("#sum{ntrials}",ftrials); 
    }  
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetMatching::Run()
{
    // execute once for each event
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    if(!(InputEvent() && fSourceJets && fTargetJets && IsEventSelected())) return kFALSE;
    if(fh1AvgTrials) fh1AvgTrials->Fill("#sum{avg ntrials}", fAvgTrials);
    // step one: do geometric matching 
    switch (fMatchingScheme) {
        // FIXME having separate dedicated functions is not necessary and historical
        // cluttered code - should be cleaned up and merged into one
        case kGeoEtaPhi : {
            DoGeometricMatchingEtaPhi();
            // break if no jet was matched
            if(!fMatchedJetContainer[1][0]) return kTRUE;
        } break;
        case kGeoR : {
            DoGeometricMatchingR();
            // break if no jet was matched
            if(!fMatchedJetContainer[1][0]) return kTRUE;
        } break;
        case kDiJet : {
            DoDiJetMatching();
        } break;
        default : break;
    }
    // optional step two: get a bijection (avoid duplicate matches)
    if(fGetBijection)           GetBijection();
    // optional step three: match constituents within matched jets
    if(fMatchConstituents)      DoConstituentMatching();
    // stream data to output
    PostMatchedJets();
    FillMatchedJetHistograms();
    #ifdef DEBUGTASK
        PrintInfo();
    #endif
    PostData(1, fOutputList);
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoGeometricMatchingEtaPhi()
{
    // do geometric matching based on eta phi distance between jet centers
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    fNoMatchedJets = 0; // reset the matched jet counter
    Int_t iSource(fSourceJets->GetEntriesFast()), iTarget(fTargetJets->GetEntriesFast());
    for(Int_t i(0); i < iSource; i++) {
        AliEmcalJet* sourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(i)));
        if(!PassesCuts(sourceJet, 0)) continue;
        for(Int_t j(0); j < iTarget; j++) {
            AliEmcalJet* targetJet(static_cast<AliEmcalJet*>(fTargetJets->At(j)));
            if(!PassesCuts(targetJet, 1)) continue;
            if (fUseEmcalBaseJetCuts && !AcceptJet(targetJet, 1)) continue;
            if((TMath::Abs(sourceJet->Eta() - targetJet->Eta()) < fMatchEta )) {
                Double_t sourcePhi(sourceJet->Phi()), targetPhi(targetJet->Phi());
                if(TMath::Abs(sourcePhi-targetPhi) > TMath::Abs(sourcePhi-targetPhi+TMath::TwoPi())) sourcePhi+=TMath::TwoPi();
                if(TMath::Abs(sourcePhi-targetPhi) > TMath::Abs(sourcePhi-targetPhi-TMath::TwoPi())) sourcePhi-=TMath::TwoPi();
                if(TMath::Abs(sourcePhi-targetPhi) < fMatchPhi) {       // accept the jets as matching 
                Bool_t isBestMatch(kTRUE);
                    if(fGetBijection) { // match has been found, for bijection only store it if there's no better match
    #ifdef DEBUGTASK
        printf(" > Entering first bijection test \n");
    #endif
                        for(Int_t k(i); k < iSource; k++) {
                            AliEmcalJet* candidateSourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(k)));
                            if(PassesCuts(candidateSourceJet, 0)) continue;
    #ifdef DEBUGTASK
        printf("source distance %.2f \t candidate distance %.2f \n", GetR(sourceJet, targetJet),GetR(candidateSourceJet, targetJet));
    #endif
                            if(GetR(sourceJet, targetJet) > GetR(candidateSourceJet, targetJet)) {
                                isBestMatch = kFALSE;
                                break;
                            }
                        }
    #ifdef DEBUGTASK
        (isBestMatch) ? printf(" kept source \n ") : printf(" we can do better (rejected source) \n");
    #endif
                    }
                    if(isBestMatch) {
                        fMatchedJetContainer[fNoMatchedJets][0] = sourceJet;
                        fMatchedJetContainer[fNoMatchedJets][1] = targetJet;
                        fNoMatchedJets++;
                    }
                    if(fNoMatchedJets > 199) {   // maximum to keep the cache at reasonable size
                        AliError(Form("%s: Found too many matched jets (> 100). Adjust matching criteria !", GetName()));
                        return;
                    }
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoGeometricMatchingR()
{
    // do geometric matching based on shortest path between jet centers
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    fNoMatchedJets = 0; // reset the matched jet counter
    Int_t iSource(fSourceJets->GetEntriesFast()), iTarget(fTargetJets->GetEntriesFast());
    for(Int_t i(0); i < iSource; i++) {
        AliEmcalJet* sourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(i)));
        if(!PassesCuts(sourceJet, 0)) continue;
        for(Int_t j(0); j < iTarget; j++) {
            AliEmcalJet* targetJet(static_cast<AliEmcalJet*>(fTargetJets->At(j)));
            if(!PassesCuts(targetJet, 1)) continue;
            else if (fUseEmcalBaseJetCuts && !AcceptJet(targetJet, 1)) continue;
            if(GetR(sourceJet, targetJet) <= fMatchR) {
                Bool_t isBestMatch(kTRUE);
                if(fGetBijection) { // match has been found, for bijection only store it if there's no better match
    #ifdef DEBUGTASK
        printf(" > Entering first bijection test \n");
    #endif
                    for(Int_t k(i); k < iSource; k++) {
                        AliEmcalJet* candidateSourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(k)));
                        if(!PassesCuts(candidateSourceJet, 0)) continue;
    #ifdef DEBUGTASK
        printf("source distance %.2f \t candidate distance %.2f \n", GetR(sourceJet, targetJet),GetR(candidateSourceJet, targetJet));
    #endif
                        if(GetR(sourceJet, targetJet) > GetR(candidateSourceJet, targetJet)) {
                            isBestMatch = kFALSE;
                            break;
                        }
                    }
    #ifdef DEBUGTASK
        (isBestMatch) ? printf(" kept source \n ") : printf(" we can do better (rejected source) \n");
    #endif
                }
                if(isBestMatch) {
                    fMatchedJetContainer[fNoMatchedJets][0] = sourceJet;
                    fMatchedJetContainer[fNoMatchedJets][1] = targetJet;
                    fNoMatchedJets++;
                }
                if(fNoMatchedJets > 99) {
                    AliError(Form("%s: Found too many matched jets (> 100). Adjust matching criteria !", GetName()));
                    return;
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoDiJetMatching()
{
    // match dijets. this is in a sense a 'special' mode of the task as both jet supplied jet arrays 
    // (target and source) are the same jet collection, matching will be performed on distribution in
    // azimuth
    // no ouptut array is produced, dedicated dijet histo's are filled in this function instead
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    // get the leading jet in a given acceptance (TPC is default)
    Int_t leadingJetIndex(-1), subLeadingJetIndex(-1);
    // retrieve the leading jet, leadingJetIndex points to the leading jet
    AliEmcalJet* leadingJet(GetLeadingJet(fSourceJets, leadingJetIndex));
    if(!leadingJet) return;
    // fill phi and eta of leading jet (should always be in selected acceptance)
    fHistDiJetLeadingJet->Fill(leadingJet->Phi(), leadingJet->Eta(), leadingJet->Pt());
    Double_t sourcePhi(leadingJet->Phi()), targetPhi(-1);
    // get the sub-leading jet - faster when leading jet is also provided
    AliEmcalJet* subLeadingJet(GetSubLeadingJet(fSourceJets, leadingJetIndex, subLeadingJetIndex));
    if(!subLeadingJet) return;
    else { // check if the sub-leading jet is actually the away-side jet
        targetPhi = subLeadingJet->Phi() + TMath::Pi();
        // rotate jets to common phase
        if(TMath::Abs(sourcePhi) > TMath::Abs(sourcePhi+TMath::TwoPi())) sourcePhi+=TMath::TwoPi();
        if(TMath::Abs(sourcePhi) > TMath::Abs(sourcePhi-TMath::TwoPi())) sourcePhi-=TMath::TwoPi();
        if(TMath::Abs(targetPhi) > TMath::Abs(targetPhi+TMath::TwoPi())) targetPhi+=TMath::TwoPi();
        if(TMath::Abs(targetPhi) > TMath::Abs(targetPhi-TMath::TwoPi())) targetPhi-=TMath::TwoPi();
        if(TMath::Abs(sourcePhi-targetPhi) < fMatchPhi) {
            fHistDiJet->Fill(subLeadingJet->Phi(), subLeadingJet->Eta(), leadingJet->Pt());
            fHistDiJetDPhi->Fill(sourcePhi-targetPhi, leadingJet->Pt());
            fHistDiJetDPt->Fill(leadingJet->Pt() - subLeadingJet->Pt(), leadingJet->Pt());
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoConstituentMatching()
{
    // match constituents within matched jets 
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    if(!fTracks) {
        AliFatal(Form("%s Fatal error! To do deep matching, supply jet constituents ! \n", GetName()));
        return; // coverity ...
    }
    for(Int_t i(0); i < fNoMatchedJets; i++) {
        AliEmcalJet* sourceJet = fMatchedJetContainer[i][0];
        AliEmcalJet* targetJet = fMatchedJetContainer[i][1];
        if(sourceJet && targetJet) {    // duplicate check: slot migth be NULL
            Double_t targetPt(0);
            Int_t iSJ(sourceJet->GetNumberOfTracks());
            Int_t iTJ(targetJet->GetNumberOfTracks());
            Int_t overlap(0), alreadyFound(0);
            for(Int_t j(0); j < iSJ; j++) {
                alreadyFound = 0;
                Int_t idSource((Int_t)sourceJet->TrackAt(j));
                for(Int_t k(0); k < iTJ; k++) { // compare all target tracks to the source track
                    if(idSource == targetJet->TrackAt(k) && alreadyFound == 0) {
                        overlap++;
                        alreadyFound++; // avoid possible duplicate matching
                        AliVParticle* vp(static_cast<AliVParticle*>(targetJet->TrackAt(k, fTracks)));
                        if(vp) targetPt += vp->Pt();
                        continue;
                    }
                }
            }
            if((float)overlap/(float)iSJ < fMinFracRecoveredConstituents || targetPt/sourceJet->Pt() < fMinFracRecoveredConstituentPt) {
    #ifdef DEBUGTASK
        printf("  \n > Purging jet, recovered constituents ratio  %i / %i = %.2f <  or pt ratio %.2f / %.2f = %.2f < %.2f", overlap, iSJ, (float)overlap/(float)iSJ, targetPt, sourceJet->Pt(), targetPt/sourceJet->Pt(), fMinFracRecoveredConstituentPt);
    #endif
                fMatchedJetContainer[i][0] = 0x0;
                fMatchedJetContainer[i][1] = 0x0;
                continue;
            }
            if(sourceJet->Pt() > 0) {
                Double_t sourceRho(0), targetRho(0);
                if(fSourceRho) sourceRho = fSourceRho->GetLocalVal(sourceJet->Phi(), fSourceRadius)*sourceJet->Area();
                if(fTargetRho) targetRho = fTargetRho->GetLocalVal(targetJet->Phi(), fTargetRadius)*targetJet->Area();
                fProfFracPtMatched->Fill(sourceJet->Pt()-sourceRho, (targetPt-targetRho) / (sourceJet->Pt()-sourceRho));
                fProfFracPtJets->Fill(sourceJet->Pt()-sourceRho, (targetJet->Pt()-targetRho) / (sourceJet->Pt()-sourceRho));
                fProfFracNoMatched->Fill(sourceJet->Pt()-sourceRho, (double)overlap / (double)sourceJet->GetNumberOfTracks());
                fProfFracNoJets->Fill(sourceJet->Pt()-sourceRho, (double)targetJet->GetNumberOfTracks() / (double)sourceJet->GetNumberOfTracks());
            }
    #ifdef DEBUGTASK
        if(fDebug > 0) {
              printf("\n > Jet A: %i const\t", iSJ);
              printf(" > Jet B %i const\t", iTJ);
              printf(" -> OVERLAP: %i tracks <- \n", overlap);
        }
    #endif
       }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::GetBijection()
{
    // bijection of source and matched jets, based on closest distance between jets
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    for(Int_t i(0); i < fNoMatchedJets; i++) {
        for(Int_t j(i+1); j < fNoMatchedJets; j++) {
            if(fMatchedJetContainer[i][0] == fMatchedJetContainer[j][0]) {
                // found source with two targets, now see which target is closer to the source
                if(!(fMatchedJetContainer[i][0] && fMatchedJetContainer[i][1] && fMatchedJetContainer[j][0] && fMatchedJetContainer[j][1] )) continue;
                Double_t rA(GetR(fMatchedJetContainer[i][0], fMatchedJetContainer[i][1]));      // distance between connection A = i
                Double_t rB(GetR(fMatchedJetContainer[j][0], fMatchedJetContainer[j][1]));      // distance between connection B = j
                if (rB > rA) {  // jet two is far away, clear it from both target and source list
                    fMatchedJetContainer[j][0] = 0x0;
                    fMatchedJetContainer[j][1] = 0x0;
                } else {                // jet one is far away, clear it from both target and source list
                    fMatchedJetContainer[i][0] = fMatchedJetContainer[j][0];    // switch to avoid breaking loop
                    fMatchedJetContainer[i][1] = fMatchedJetContainer[j][1];    
                    fMatchedJetContainer[j][0] = 0x0;                           // clear duplicate jet from cache
                    fMatchedJetContainer[j][1] = 0x0;
                }
    #ifdef DEBUGTASK
         printf(" found duplicate jet, chose %.2f over %.2f \n" , (rB > rA) ? rA : rB, (rB > rA) ? rB : rA);
    #endif
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::PostMatchedJets()
{
    // post matched jets
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    fMatchedJets->Delete();     // should be a NULL operation, but added just in case
    for(Int_t i(0), p(0); i < fNoMatchedJets; i++) {
        if(fMatchedJetContainer[i][1]) {        // duplicate jets will have NULL value here and are skipped
            new((*fMatchedJets)[p]) AliEmcalJet(*fMatchedJetContainer[i][1]);
            p++; 
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::FillAnalysisSummaryHistogram() const
{
    // fill the analysis summary histrogram, saves all relevant analysis settigns
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(1, "fUseScaledRho");
    fHistAnalysisSummary->SetBinContent(1, (int)fUseScaledRho);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fMatchingScheme");
    fHistAnalysisSummary->SetBinContent(2, (int)fMatchingScheme);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fSourceBKG");
    fHistAnalysisSummary->SetBinContent(3, (int)fSourceBKG);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "fTargetBKG");
    fHistAnalysisSummary->SetBinContent(4, (int)fTargetBKG);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(5, "fMatchEta");
    fHistAnalysisSummary->SetBinContent(5, fMatchEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(6, "fMatchPhi");
    fHistAnalysisSummary->SetBinContent(6, fMatchPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(7, "fMatchR");
    fHistAnalysisSummary->SetBinContent(7, fMatchR);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(8, "iter");
    fHistAnalysisSummary->SetBinContent(8, 1.);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::FillMatchedJetHistograms() 
{
    // fill matched jet histograms
    #ifdef DEBUGTASK
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    for(Int_t i(0); i < fSourceJets->GetEntriesFast(); i++) {
        AliEmcalJet* source = static_cast<AliEmcalJet*>(fSourceJets->At(i));
        if(!source) continue;
        else if(fUseEmcalBaseJetCuts && !AcceptJet(source, 0)) continue;
        Double_t sourceRho(0), targetRho(0);
        if(fSourceRho) sourceRho = fSourceRho->GetLocalVal(source->Phi(), fSourceRadius)*source->Area();
        fHistSourceJetPt->Fill(source->Pt()-sourceRho);
        fHistNoConstSourceJet->Fill(source->GetNumberOfConstituents());
        for(Int_t j(0); j < fTargetJets->GetEntriesFast(); j++) {
            AliEmcalJet* target = static_cast<AliEmcalJet*>(fTargetJets->At(j));
            if(target) {
                if(fUseEmcalBaseJetCuts && !AcceptJet(target, 1)) continue;
                if(fTargetRho) targetRho = fTargetRho->GetLocalVal(target->Phi(), fTargetRadius)*target->Area();
                fProfQA->Fill(0.5, TMath::Abs((source->Pt()-sourceRho)-(target->Pt()-targetRho)));  
                fProfQA->Fill(1.5, TMath::Abs(source->Eta()-target->Eta()));
                fProfQA->Fill(2.5, TMath::Abs(source->Phi()-target->Phi()));
                fHistUnsortedCorrelation->Fill(PhaseShift(source->Phi()-target->Phi(), 2));
                if(j==0) {
                    fHistTargetJetPt->Fill(target->Pt()-targetRho);
                    fHistNoConstTargetJet->Fill(target->GetNumberOfConstituents());
                }
            }
        }
    }
    for(Int_t i(0); i < fNoMatchedJets; i++) {
        if(fMatchedJetContainer[i][0] && fMatchedJetContainer[i][1]) {
            Double_t sourceRho(0), targetRho(0);
            if(fSourceRho) sourceRho = fSourceRho->GetLocalVal(fMatchedJetContainer[i][0]->Phi(), fSourceRadius)*fMatchedJetContainer[i][0]->Area();
            if(fTargetRho) targetRho = fTargetRho->GetLocalVal(fMatchedJetContainer[i][1]->Phi(), fTargetRadius)*fMatchedJetContainer[i][1]->Area();
            fHistMatchedCorrelation->Fill(PhaseShift(fMatchedJetContainer[i][0]->Phi()-fMatchedJetContainer[i][1]->Phi(), 2));
            fHistMatchedSourceJetPt->Fill(fMatchedJetContainer[i][0]->Pt()-sourceRho);
            fHistMatchedJetPt->Fill(fMatchedJetContainer[i][1]->Pt()-targetRho);
            fHistNoConstMatchJet->Fill(fMatchedJetContainer[i][1]->Pt()-targetRho);
            fProfQAMatched->Fill(0.5, TMath::Abs((fMatchedJetContainer[i][0]->Pt()-sourceRho)-(fMatchedJetContainer[i][1]->Pt()-targetRho)));
            fProfQAMatched->Fill(1.5, TMath::Abs(fMatchedJetContainer[i][0]->Eta()-fMatchedJetContainer[i][1]->Eta()));
            fProfQAMatched->Fill(2.5, TMath::Abs(fMatchedJetContainer[i][0]->Phi()-fMatchedJetContainer[i][1]->Phi()));
            
            fHistSourceMatchedJetPt->Fill(fMatchedJetContainer[i][0]->Pt()-sourceRho, fMatchedJetContainer[i][1]->Pt()-targetRho);
            if(fDoDetectorResponse) {
                fProfFracPtJets->Fill(fMatchedJetContainer[i][0]->Pt()-sourceRho, (fMatchedJetContainer[i][1]->Pt()-targetRho) / (fMatchedJetContainer[i][0]->Pt()-sourceRho));
                fProfFracNoJets->Fill(fMatchedJetContainer[i][0]->Pt()-sourceRho, (double)fMatchedJetContainer[i][1]->GetNumberOfTracks() / (double)fMatchedJetContainer[i][0]->GetNumberOfTracks());
                fHistDetectorResponseProb->Fill((fMatchedJetContainer[i][1]->Pt()-fMatchedJetContainer[i][0]->Pt())/fMatchedJetContainer[i][0]->Pt(), fMatchedJetContainer[i][0]->Pt());
            }
        }
    }
}
//_____________________________________________________________________________
AliEmcalJet* AliAnalysisTaskJetMatching::GetLeadingJet(TClonesArray* source, Int_t &leadingJetIndex, Double_t etaMin, Double_t etaMax) 
{
    // return the leading jet within giiven acceptance
    Int_t iJets(source->GetEntriesFast());
    Double_t pt(0);
    AliEmcalJet* leadingJet(0x0);
    for(Int_t i(0); i < iJets; i++) {
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(source->At(i));
        if(jet->Eta() < etaMin || jet->Eta() > etaMax) continue;
        if(jet->Pt() > pt) {
           leadingJet = jet;
           pt = leadingJet->Pt();
           leadingJetIndex = i;
        }
    }
    return leadingJet;
}
//_____________________________________________________________________________
AliEmcalJet* AliAnalysisTaskJetMatching::GetSubLeadingJet(TClonesArray* source, Int_t leadingJetIndex, Int_t &subLeadingJetIndex) 
{
    // return the sub-leading jet within given acceptance
    // same as GetLeadingJet() but skips the leading jet (so returned jet is 
    // sub-leading by design)
    // there is no eta requirement on the location of the sub-leading jet
    Int_t iJets(source->GetEntriesFast());
    // if the leading jet isn't given, retrieve it
    if(leadingJetIndex < 0) GetLeadingJet(source, leadingJetIndex);
    AliEmcalJet *leadingJet(0x0), *prevLeadingJet(static_cast<AliEmcalJet*>(source->At(leadingJetIndex)));
    if(!prevLeadingJet) return 0x0;
    Double_t pt(0), leadingPt(prevLeadingJet->Pt());
    for(Int_t i(0); i < iJets; i++) {
        if(i == leadingJetIndex) continue;
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(source->At(i));
        // check if jet is actually sub-leading
        if(jet->Pt() > pt && jet->Pt() <= leadingPt) {
           leadingJet = jet;
           pt = leadingJet->Pt();
           subLeadingJetIndex = i;
        }
    }
    return leadingJet;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::PrintInfo() const
{
    // print some info 
    printf("\n > No. of source jets from %s \n \t %i \n ", fSourceJetsName.Data(), fSourceJets->GetEntriesFast());
    printf(" > No. of target jets from %s \n \t %i \n ", fTargetJetsName.Data(), fTargetJets->GetEntriesFast());
    printf(" > No. of matched jets from %s \n \t %i \n ", fMatchedJetsName.Data(), fMatchedJets->GetEntriesFast());
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________________
