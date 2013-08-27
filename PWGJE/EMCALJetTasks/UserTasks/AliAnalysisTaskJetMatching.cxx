// 
// General task to match two sets of jets
//
// This task takes two TClonesArray's as input: 
// [1] fSourceJets - e.g. pythia jets
// [2] fTargetJets - e.g. a samle containing pythia jets embedded in a pbpb event
// the task will try to match jets from the source array to the target array, and
// save the found TARGET jets in a new array 'fMatchedJets'
//
// Author: Redmer Alexander Bertens, Utrecht Univeristy, Utrecht, Netherlands
// rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 

// root includes
#include <TClonesArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
// aliroot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
// emcal jet framework includes
#include <AliEmcalJet.h>
#include <AliAnalysisTaskJetMatching.h>

class AliAnalysisTaskJetMatching;
using namespace std;

ClassImp(AliAnalysisTaskJetMatching)

AliAnalysisTaskJetMatching::AliAnalysisTaskJetMatching() : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetMatching", kTRUE), 
    fDebug(0), fSourceJets(0), fSourceJetsName(0), fTargetJets(0), fTargetJetsName(0), fMatchedJets(0), fMatchedJetsName(GetName()), fUseScaledRho(0), fMatchingScheme(kGeoEtaPhi), fDuplicateJetRecoveryMode(kDoNothing), fUseEmcalBaseJetCuts(kFALSE), fSourceBKG(kNoSourceBKG), fTargetBKG(kNoTargetBKG), fLocalJetMinEta(-10), fLocalJetMaxEta(-10), fLocalJetMinPhi(-10), fLocalJetMaxPhi(-10), fOutputList(0), fHistUnsortedCorrelation(0), fHistMatchedCorrelation(0), fHistSourceJetPt(0), fHistTargetJetPt(0), fHistMatchedJetPt(0), fHistNoConstSourceJet(0), fHistNoConstTargetJet(0), fHistNoConstMatchJet(0), fProfFracPtMatched(0), fProfFracPtJets(0), fProfFracNoMatched(0), fProfFracNoJets(0), fHistAnalysisSummary(0), fProfQAMatched(0), fProfQA(0), fNoMatchedJets(100), fMatchEta(.03), fMatchPhi(.03), fMatchR(.03), fMatchArea(0), fMaxRelEnergyDiff(.1), fMaxAbsEnergyDiff(5) {
    // default constructor
    ClearMatchedJetsCache();
}
//_____________________________________________________________________________
AliAnalysisTaskJetMatching::AliAnalysisTaskJetMatching(const char* name) : AliAnalysisTaskEmcalJet(name, kTRUE),
    fDebug(0), fSourceJets(0), fSourceJetsName(0), fTargetJets(0), fTargetJetsName(0), fMatchedJets(0), fMatchedJetsName(GetName()), fUseScaledRho(0), fMatchingScheme(kGeoEtaPhi), fDuplicateJetRecoveryMode(kDoNothing), fUseEmcalBaseJetCuts(kFALSE), fSourceBKG(kNoSourceBKG), fTargetBKG(kNoTargetBKG), fLocalJetMinEta(-10), fLocalJetMaxEta(-10), fLocalJetMinPhi(-10), fLocalJetMaxPhi(-10), fOutputList(0), fHistUnsortedCorrelation(0), fHistMatchedCorrelation(0), fHistSourceJetPt(0), fHistTargetJetPt(0), fHistMatchedJetPt(0), fHistNoConstSourceJet(0), fHistNoConstTargetJet(0), fHistNoConstMatchJet(0), fProfFracPtMatched(0), fProfFracPtJets(0), fProfFracNoMatched(0), fProfFracNoJets(0), fHistAnalysisSummary(0), fProfQAMatched(0), fProfQA(0), fNoMatchedJets(100), fMatchEta(.03), fMatchPhi(.03), fMatchR(.03), fMatchArea(0), fMaxRelEnergyDiff(.1), fMaxAbsEnergyDiff(5) {
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
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fLocalJetMinEta > -10 && fLocalJetMaxEta > -10) SetJetEtaLimits(fLocalJetMinEta, fLocalJetMaxEta);
    if(fLocalJetMinPhi > -10 && fLocalJetMaxPhi > -10) SetJetPhiLimits(fLocalJetMinPhi, fLocalJetMaxPhi);
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
    AliAnalysisTaskEmcalJet::ExecOnce(); // init base class
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::UserCreateOutputObjects()
{
    // create output objects
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
     // add the matched jets array to the event
    fMatchedJets = new TClonesArray("AliEmcalJet");
    fMatchedJets->SetName(fMatchedJetsName);
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    // add analysis histograms
    fHistUnsortedCorrelation = BookTH1F("fHistUnsortedCorrelation", "#Delta #varphi unsorted", 50, 0, TMath::Pi());
    fHistMatchedCorrelation = BookTH1F("fHistMatchedCorrelation", "#Delta #varphi matched", 50, 0, TMath::Pi());
    fHistSourceJetPt = BookTH1F("fHistSourceJetPt", "p_{t} [GeV/c]", 50, 0, 150);
    fHistTargetJetPt = BookTH1F("fHistTargetJetPt", "p_{t} [GeV/c]", 50, 0, 150);
    fHistMatchedJetPt = BookTH1F("fHistMatchedJetPt", "p_{t} [GeV/c]", 50, 0, 150);
    fHistNoConstSourceJet = BookTH1F("fHistNoConstSourceJet", "number of constituents", 50, 0, 50);
    fHistNoConstTargetJet = BookTH1F("fHistNoConstTargetJet", "number of constituents", 50, 0, 50);
    fHistNoConstMatchJet = BookTH1F("fHistNoConstMatchJet", "number of constituents", 50, 0, 50);
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
    fProfFracPtMatched = new TProfile("fProfFracPtMatched", "fProfFracPtMatched", 10, 0, 200);
    fOutputList->Add(fProfFracPtMatched);
    fProfFracPtJets = new TProfile("fProfFracPtJets", "fProfFracPtJets", 10, 0, 200);
    fOutputList->Add(fProfFracPtJets);
    fProfFracNoMatched = new TProfile("fProfFracNoMatched", "fProfFracNoMatched", 10, 0, 200);
    fOutputList->Add(fProfFracNoMatched);
    fProfFracNoJets = new TProfile("fProfFracNoJets", "fProfFracNoJets", 10, 0, 200);
    fOutputList->Add(fProfFracNoJets);
    fOutputList->Sort();
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetMatching::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Bool_t append)
{
    // book a TH1F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fOutputList) return 0x0;
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
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fOutputList) return 0x0;
    TString title(name);
    title += Form(";%s;%s", x, y);
    TH2F* histogram = new TH2F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetMatching::Run()
{
    // execute once for each event
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!(InputEvent() && fSourceJetsName && fTargetJets && IsEventSelected())) return kFALSE;
    // single event loop
    switch (fMatchingScheme) {
        case kGeoEtaPhi : {
            DoGeometricMatchingEtaPhi();
        } break;
        case kGeoR : {
            DoGeometricMatchingR();
        } break;
        case kGeoEtaPhiArea : {
            DoGeometricMatchingEtaPhi(kTRUE);
            } break;
        case kGeoRArea : {
            DoGeometricMatchingR(kTRUE);
            } break;
        case kDeepMatching : {
            DoGeometricMatchingEtaPhi();
            } break;
       default : break;
    }
    if(fMatchedJetContainer[1][0]) {       // if matched jets are found, fill some more histograms
        switch (fDuplicateJetRecoveryMode) {
            case kDoNothing : break;
            default : {
                DuplicateJetRecovery();
                break;
            }
        }
    }
    // if required do deep matching, i.e. match constituents in source and target jets 
    switch (fMatchingScheme) {
        case kDeepMatching : {
            DoDeepMatching();
            break; }
        default : break;
    }
    // stream data to output
    PostMatchedJets();
    FillMatchedJetHistograms();
    if(fDebug > 0) PrintInfo();
    PostData(1, fOutputList);
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoGeometricMatchingEtaPhi(Bool_t pairCuts)
{
    // do geometric matching
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fNoMatchedJets = 0; // reset the matched jet counter
    Int_t iSource(fSourceJets->GetEntriesFast()), iTarget(fTargetJets->GetEntriesFast());
    for(Int_t i(0); i < iSource; i++) {
        AliEmcalJet* sourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(i)));
        if(!PassesCuts(sourceJet)) continue;
        for(Int_t j(0); j < iTarget; j++) {
            AliEmcalJet* targetJet(static_cast<AliEmcalJet*>(fTargetJets->At(j)));
            if(!PassesCuts(targetJet)) continue;
            if((TMath::Abs(sourceJet->Eta() - targetJet->Eta()) < fMatchEta ) && (TMath::Abs(sourceJet->Phi()-targetJet->Phi()) < fMatchPhi)) {
                if(pairCuts && !PassesCuts(sourceJet, targetJet)) continue;
                fMatchedJetContainer[fNoMatchedJets][0] = sourceJet;
                fMatchedJetContainer[fNoMatchedJets][1] = targetJet;
                fNoMatchedJets++;
                if(fNoMatchedJets > 99) {
                    AliError(Form("%s: Found too many matched jets (> 100). Adjust matching criteria !", GetName()));
                    return;
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoGeometricMatchingR(Bool_t pairCuts)
{
    // do geometric matching based on shortest path between jet centers
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fNoMatchedJets = 0; // reset the matched jet counter
    Int_t iSource(fSourceJets->GetEntriesFast()), iTarget(fTargetJets->GetEntriesFast());
    for(Int_t i(0); i < iSource; i++) {
        AliEmcalJet* sourceJet(static_cast<AliEmcalJet*>(fSourceJets->At(i)));
        if(!PassesCuts(sourceJet)) continue;
        for(Int_t j(0); j < iTarget; j++) {
            AliEmcalJet* targetJet(static_cast<AliEmcalJet*>(fTargetJets->At(j)));
            if(!PassesCuts(targetJet)) continue;
            Double_t etaS(sourceJet->Eta()), etaT(targetJet->Eta());
            Double_t phiS(sourceJet->Phi()), phiT(targetJet->Phi());
            // if necessary change phase
            if(TMath::Abs(phiS-phiT) > TMath::Abs(phiS-phiT + TMath::TwoPi())) phiS+=TMath::TwoPi();
            if(TMath::Abs(phiS-phiT) > TMath::Abs(phiS-phiT - TMath::TwoPi())) phiS-=TMath::TwoPi();
            if(TMath::Sqrt(TMath::Abs((etaS-etaT)*(etaS-etaT)+(phiS-phiT)*(phiS-phiT)) <= fMatchR)) {
                if(pairCuts && !PassesCuts(sourceJet, targetJet)) continue;
                fMatchedJetContainer[fNoMatchedJets][0] = sourceJet;
                fMatchedJetContainer[fNoMatchedJets][1] = targetJet;
                fNoMatchedJets++;
                if(fNoMatchedJets > 99) {
                    AliError(Form("%s: Found too many matched jets (> 100). Adjust matching criteria !", GetName()));
                    return;
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DoDeepMatching()
{
    // match constituents, can be VERY slow 
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
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
            if(sourceJet->Pt() > 0) {
                fProfFracPtMatched->Fill(sourceJet->Pt(), targetPt / sourceJet->Pt());
                fProfFracPtJets->Fill(sourceJet->Pt(), targetJet->Pt() / sourceJet->Pt());
                fProfFracNoMatched->Fill(sourceJet->Pt(), (double)overlap / (double)sourceJet->GetNumberOfTracks());
                fProfFracNoJets->Fill(sourceJet->Pt(), (double)targetJet->GetNumberOfTracks() / (double)sourceJet->GetNumberOfTracks());
            }
            if(fDebug > 0) {
                printf("\n\n > Jet a has %i constituents \n", iSJ);
                printf(" > Jet b has %i constituents \n", iTJ);
                printf("  -OVERLAP %i tracks-\n\n", overlap);
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::DuplicateJetRecovery()
{
    // find target jets that have been matched to a source jet more than once - uses nested loops!
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t iDuplicateJets(0);            // counter for duplicate jets
    for(Int_t i(0); i < fNoMatchedJets; i++) {
        for(Int_t j(i+1); j < fNoMatchedJets; j++) {
            if(fMatchedJetContainer[i][1] == fMatchedJetContainer[j][1]) {
                iDuplicateJets++;
                switch (fDuplicateJetRecoveryMode) {
                    case kTraceDuplicates : { 
                        printf(" > found duplicate jet <\n");
                        break; }
                    case kRemoveDuplicates : { 
                         fMatchedJetContainer[j][1] = NULL;
                         break; }
                    default : break;
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::PostMatchedJets()
{
    // post matched jets
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
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
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(1, "fUseScaledRho");
    fHistAnalysisSummary->SetBinContent(1, (int)fUseScaledRho);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fMatchingScheme");
    fHistAnalysisSummary->SetBinContent(2, (int)fMatchingScheme);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fSourceBKG");
    fHistAnalysisSummary->SetBinContent(3, (int)fSourceBKG);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "fTargetBKG");
    fHistAnalysisSummary->SetBinContent(4, (int)fTargetBKG);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(5, "fLocalJetMinEta");
    fHistAnalysisSummary->SetBinContent(5, fLocalJetMinEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(6, "fLocalJetMaxEta");
    fHistAnalysisSummary->SetBinContent(6, fLocalJetMaxEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(7, "fLocalJetMinPhi");
    fHistAnalysisSummary->SetBinContent(7, fLocalJetMinPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(8, "fLocalJetMaxPhi");
    fHistAnalysisSummary->SetBinContent(8, fLocalJetMaxPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(9, "fMatchEta");
    fHistAnalysisSummary->SetBinContent(9, fMatchEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(10, "fMatchPhi");
    fHistAnalysisSummary->SetBinContent(10, fMatchPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(11, "fMatchR");
    fHistAnalysisSummary->SetBinContent(11, fMatchR);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(12, "fMatchArea");
    fHistAnalysisSummary->SetBinContent(12, fMatchArea);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(13, "fMaxRelEnergyDiff");
    fHistAnalysisSummary->SetBinContent(13, fMaxRelEnergyDiff);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(13, "fMaxAbsEnergyDiff");
    fHistAnalysisSummary->SetBinContent(13, fMaxAbsEnergyDiff);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(14, "iter");
    fHistAnalysisSummary->SetBinContent(13, 1.);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::FillMatchedJetHistograms() const
{
    // fill matched jet histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    for(Int_t i(0); i < fSourceJets->GetEntriesFast(); i++) {
        AliEmcalJet* source = static_cast<AliEmcalJet*>(fSourceJets->At(i));
        if(!source) continue;
        fHistSourceJetPt->Fill(source->Pt());
        fHistNoConstSourceJet->Fill(source->GetNumberOfConstituents());
        for(Int_t j(0); j < fTargetJets->GetEntriesFast(); j++) {
            AliEmcalJet* target = static_cast<AliEmcalJet*>(fTargetJets->At(j));
            if(target) {
            fProfQA->Fill(0.5, TMath::Abs(source->Pt()-target->Pt()));
            fProfQA->Fill(1.5, TMath::Abs(source->Eta()-target->Eta()));
            fProfQA->Fill(2.5, TMath::Abs(source->Phi()-target->Phi()));
 
                fHistUnsortedCorrelation->Fill(PhaseShift(source->Phi()-target->Phi(), 2));
                if(j==0) {
                    fHistTargetJetPt->Fill(target->Pt());
                    fHistNoConstTargetJet->Fill(target->GetNumberOfConstituents());
                }
            }
        }
    }
    for(Int_t i(0); i < fNoMatchedJets; i++) {
        if(fMatchedJetContainer[i][0] && fMatchedJetContainer[i][1]) {
            fHistMatchedCorrelation->Fill(PhaseShift(fMatchedJetContainer[i][0]->Phi()-fMatchedJetContainer[i][1]->Phi(), 2));
            fHistMatchedJetPt->Fill(fMatchedJetContainer[i][1]->Pt());
            fHistNoConstMatchJet->Fill(fMatchedJetContainer[i][1]->Pt());
            fProfQAMatched->Fill(0.5, TMath::Abs(fMatchedJetContainer[i][0]->Pt()-fMatchedJetContainer[i][1]->Pt()));
            fProfQAMatched->Fill(1.5, TMath::Abs(fMatchedJetContainer[i][0]->Eta()-fMatchedJetContainer[i][1]->Eta()));
            fProfQAMatched->Fill(2.5, TMath::Abs(fMatchedJetContainer[i][0]->Phi()-fMatchedJetContainer[i][1]->Phi()));
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::PrintInfo() const
{
    // print some info 
    printf(" > No. of source jets from %s \n \t %i \n ", fSourceJetsName.Data(), fSourceJets->GetEntriesFast());
    printf(" > No. of target jets from %s \n \t %i \n ", fTargetJetsName.Data(), fTargetJets->GetEntriesFast());
    printf(" > No. of matched jets from %s \n \t %i \n ", fMatchedJetsName.Data(), fMatchedJets->GetEntriesFast());
    if(fDebug > 3) InputEvent()->GetList()->ls();
}
//_____________________________________________________________________________
void AliAnalysisTaskJetMatching::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________________
