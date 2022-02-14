// /*
//  * AliAnalysisTaskPOmegaPenne.cxx
//  *
//  *  Created on: 11 Dec 2019
//  *      Author: Boris Bajtl
//  */



#include "AliAnalysisTaskPOmegaPenne.h"
#include <string.h>
#include "AliNanoAODTrack.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
// #include <chrono>

using std::cout;
using std::endl;

// carefull: not MC and second set of cuts possible!!
// #define RUN_SECOND_SET_OF_CUTS
ClassImp(AliAnalysisTaskPOmegaPenne)

    AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne() :  AliAnalysisTaskSE(),
                                                                fmixBeforePC(0),
                                                                fmixAfterPC(0),
                                                                fisInvMassPairClean(0),
                                                                fmultTrigger(0),
                                                                ffullBlastQA(0),
                                                                VEvent(0),
                                                                VTrack(0),
                                                                fEvent(0),
                                                                fTrack(0),
                                                                fConfig(0),
                                                                fGTI(0),
                                                                fTrackBufferSize(5000),
                                                                fEventCuts(0),
                                                                fv0(0),
                                                                fLambdaV0Cuts(0),
                                                                fAntiLambdaV0Cuts(0),
                                                                fPairCleaner(0),
                                                                fPartColl(0),
                                                                vLambda(0),
                                                                vAntiLambda(0),
                                                                fPairCleaner2(0),
                                                                fPartColl2(0),
                                                                fPartColl3(0),
                                                                tlEventCuts(0),
                                                                tlLambdaList(0),
                                                                tlAntiLambdaList(0),
                                                                tlResults(0),
                                                                tlResults2(0),
                                                                tlResults3(0),
                                                                tlResultsQA(0),
                                                                tlLambdaMC(0),
                                                                tlAntiLambdaMC(0),
                                                                tlRecombination_before(0),
                                                                tlRecombination_after(0),
                                                                hInvMassLambda_sanityCheck_before(0),
                                                                hInvMassLambda_total_before(0),
                                                                hInvMassLambda_shared_pion_before(0),
                                                                hInvMassLambda_shared_proton_before(0),
                                                                hInvMassLambda_shared_lambda_before(0),
                                                                hInvMassXi_sanityCheck_before(0),
                                                                hInvMassXi_total_before(0),
                                                                hInvMassXi_shared_bach_before(0),
                                                                hInvMassXi_shared_pi_daugh_before(0),
                                                                hInvMassXi_shared_prot_daugh_before(0),
                                                                hInvMassXi_shared_Lambda_before(0),
                                                                hInvMassXi_shared_pion_bach_prot_daugh_before(0),
                                                                hInvMassXi_nothing_shared(0),
                                                                hInvMassAntiLambda_sanityCheck_before(0),
                                                                hInvMassAntiLambda_total_before(0),
                                                                hInvMassAntiLambda_shared_pion_before(0),
                                                                hInvMassAntiLambda_shared_proton_before(0),
                                                                hInvMassAntiLambda_shared_lambda_before(0),
                                                                hInvMassAntiXi_sanityCheck_before(0),
                                                                hInvMassAntiXi_total_before(0),
                                                                hInvMassAntiXi_shared_bach_before(0),
                                                                hInvMassAntiXi_shared_pi_daugh_before(0),
                                                                hInvMassAntiXi_shared_prot_daugh_before(0),
                                                                hInvMassAntiXi_shared_Lambda_before(0),
                                                                hInvMassAntiXi_shared_pion_bach_prot_daugh_before(0),
                                                                hInvMassAntiXi_nothing_shared(0),
                                                                fEvtCounterBefore(0),
                                                                tlLambdaRecombination_after(0),
                                                                tlAntiLambdaRecombination_after(0),
                                                                tlXiRecombination_after(0),
                                                                tlAntiXiRecombination_after(0),
                                                                hInvMassLambda_sanityCheck_after(0),
                                                                hInvMassLambda_pi_bach_Xi_after(0),
                                                                hInvMassLambda_pi_daugh_Xi_after(0),
                                                                hInvMassLambda_prot_Xi_after(0),
                                                                hInvMassLambda_full_lambda_from_Xi_after(0),
                                                                hInvMassXi_sanityCheck_after(0),
                                                                hInvMassXi_Lamda_pi_daugh_after(0),
                                                                hInvMassXi_Lamda_prot_daugh_after(0),
                                                                hInvMassXi_Lamda_pi_bach_after(0),
                                                                hInvMassXi_Lamda_full_after(0),
                                                                hInvMassXi_Lamda_pi_no_correctLambdaMass(0),
                                                                hInvMassXi_Lamda_prot_no_correctLambdaMass(0),
                                                                hInvMassAntiLambda_sanityCheck_after(0),
                                                                hInvMassAntiLambda_pi_bach_Xi_after(0),
                                                                hInvMassAntiLambda_pi_daugh_Xi_after(0),
                                                                hInvMassAntiLambda_prot_Xi_after(0),
                                                                hInvMassAntiLambda_full_lambda_from_Xi_after(0),
                                                                hInvMassAntiXi_sanityCheck_after(0),
                                                                hInvMassAntiXi_AntiLamda_antipi_daugh_after(0),
                                                                hInvMassAntiXi_AntiLamda_antiprot_daugh_after(0),
                                                                hInvMassAntiXi_AntiLamda_antipi_bach_after(0),
                                                                hInvMassAntiXi_AntiLamda_full_after(0),
                                                                hInvMassAntiXi_AntiLamda_antipi_no_correctAntiLambdaMass(0),
                                                                hInvMassAntiXi_AntiLamda_antiprot_no_correctAntiLambdaMass(0),
                                                                fEvtCounterAfter(0),
                                                                tlInvMassPairClean(0),
                                                                tlCleanDecay(0),
                                                                hLambdaCleanedPartMassDiffToPDG_Decay(0),
                                                                hAntiLambdaCleanedPartMassDiffToPDG_Decay(0),
                                                                hXiCleanedPartMassDiffToPDG_Decay(0),
                                                                hAntiXiCleanedPartMassDiffToPDG_Decay(0),
                                                                hLambdaCleanedPartMass_Decay(0),
                                                                hAntiLambdaCleanedPartMass_Decay(0),
                                                                hXiCleanedPartMass_Decay(0),
                                                                hAntiXiCleanedPartMass_Decay(0),
                                                                tlCleanDecayAndDecay(0),
                                                                hLambdaCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                hXiCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                hAntiXiCleanedPartMassDiffToPDG_DecayDecay(0),    
                                                                hLambdaCleanedPartMass_DecayDecay(0),
                                                                hAntiLambdaCleanedPartMass_DecayDecay(0),
                                                                hXiCleanedPartMass_DecayDecay(0),
                                                                hAntiXiCleanedPartMass_DecayDecay(0),
                                                                tlCPA_PairClean_stats(0),
                                                                tlLambda_CPA_stats(0),
                                                                tlAntiLambda_CPA_stats(0),
                                                                tlXi_CPA_stats(0),
                                                                tlAntiXi_CPA_stats(0),
                                                                tlLambdaCPA_MC(0),
                                                                tlAntiLambdaCPA_MC(0),
                                                                tlXiCPA_MC(0),
                                                                tlAntiXiCPA_MC(0),
                                                                tlCPA_pT_Pairclean_CPA(0),
                                                                tlCPA_pT_Pairclean_InvMass(0),
                                                                h2_CPA_pt(0),
                                                                CPAPtBinningPrim_lambda(0),
                                                                CPAPtBinningMat_lambda(0),
                                                                CPAPtBinningSec_lambda(0),
                                                                CPAPtBinningCont_lambda(0),
                                                                CPAPtBinningPrim_lambda_dump(0),
                                                                CPAPtBinningMat_lambda_dump(0),
                                                                CPAPtBinningSec_lambda_dump(0),
                                                                CPAPtBinningCont_lambda_dump(0),
                                                                CPAPtBinningPrim_antilambda(0),
                                                                CPAPtBinningMat_antilambda(0),
                                                                CPAPtBinningSec_antilambda(0),
                                                                CPAPtBinningCont_antilambda(0),
                                                                CPAPtBinningPrim_xi(0),
                                                                CPAPtBinningMat_xi(0),
                                                                CPAPtBinningSec_xi(0),
                                                                CPAPtBinningCont_xi(0),
                                                                CPAPtBinningPrim_xi_dump(0),
                                                                CPAPtBinningMat_xi_dump(0),
                                                                CPAPtBinningSec_xi_dump(0),
                                                                CPAPtBinningCont_xi_dump(0),
                                                                CPAPtBinningPrim_antixi(0),
                                                                CPAPtBinningMat_antixi(0),
                                                                CPAPtBinningSec_antixi(0),
                                                                CPAPtBinningCont_antixi(0),
                                                                h2_pt_invMass(0)
                                                                // kStarXiLambda_unchanged(0),
                                                                // kStarXiLambda_changed(0),
                                                                // kStarAntiXiAntiLambda_unchanged(0),
                                                                // kStarAntiXiAntiLambda_changed(0)
{
 
}
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const char *name, bool isMC) : AliAnalysisTaskSE(name),
                                                                                        fmixBeforePC(0),
                                                                                        fmixAfterPC(0),
                                                                                        fisInvMassPairClean(0),
                                                                                        fmultTrigger(0),
                                                                                        ffullBlastQA(0),
                                                                                        VEvent(0),
                                                                                        VTrack(0),
                                                                                        fEvent(0),
                                                                                        fTrack(0),
                                                                                        fConfig(0),
                                                                                        fGTI(0),
                                                                                        fTrackBufferSize(5000),
                                                                                        fEventCuts(0),
                                                                                        fv0(0),
                                                                                        fLambdaV0Cuts(0),
                                                                                        fAntiLambdaV0Cuts(0),
                                                                                        fPairCleaner(0),
                                                                                        fPartColl(0),
                                                                                        vLambda(0),
                                                                                        vAntiLambda(0),
                                                                                        fPairCleaner2(0),
                                                                                        fPartColl2(0),
                                                                                        fPartColl3(0),
                                                                                        tlEventCuts(0),
                                                                                        tlLambdaList(0),
                                                                                        tlAntiLambdaList(0),
                                                                                        tlResults(0),
                                                                                        tlResults2(0),
                                                                                        tlResults3(0),
                                                                                        tlResultsQA(0),
                                                                                        tlLambdaMC(0),
                                                                                        tlAntiLambdaMC(0),
                                                                                        tlRecombination_before(0),
                                                                                        tlRecombination_after(0),
                                                                                        hInvMassLambda_sanityCheck_before(0),
                                                                                        hInvMassLambda_total_before(0),
                                                                                        hInvMassLambda_shared_pion_before(0),
                                                                                        hInvMassLambda_shared_proton_before(0),
                                                                                        hInvMassLambda_shared_lambda_before(0),
                                                                                        hInvMassXi_sanityCheck_before(0),
                                                                                        hInvMassXi_total_before(0),
                                                                                        hInvMassXi_shared_bach_before(0),
                                                                                        hInvMassXi_shared_pi_daugh_before(0),
                                                                                        hInvMassXi_shared_prot_daugh_before(0),
                                                                                        hInvMassXi_shared_Lambda_before(0),
                                                                                        hInvMassXi_shared_pion_bach_prot_daugh_before(0),
                                                                                        hInvMassXi_nothing_shared(0),
                                                                                        hInvMassAntiLambda_sanityCheck_before(0),
                                                                                        hInvMassAntiLambda_total_before(0),
                                                                                        hInvMassAntiLambda_shared_pion_before(0),
                                                                                        hInvMassAntiLambda_shared_proton_before(0),
                                                                                        hInvMassAntiLambda_shared_lambda_before(0),
                                                                                        hInvMassAntiXi_sanityCheck_before(0),
                                                                                        hInvMassAntiXi_total_before(0),
                                                                                        hInvMassAntiXi_shared_bach_before(0),
                                                                                        hInvMassAntiXi_shared_pi_daugh_before(0),
                                                                                        hInvMassAntiXi_shared_prot_daugh_before(0),
                                                                                        hInvMassAntiXi_shared_Lambda_before(0),
                                                                                        hInvMassAntiXi_shared_pion_bach_prot_daugh_before(0),
                                                                                        hInvMassAntiXi_nothing_shared(0),
                                                                                        fEvtCounterBefore(0),
                                                                                        tlLambdaRecombination_after(0),
                                                                                        tlAntiLambdaRecombination_after(0),
                                                                                        tlXiRecombination_after(0),
                                                                                        tlAntiXiRecombination_after(0),
                                                                                        hInvMassLambda_sanityCheck_after(0),
                                                                                        hInvMassLambda_pi_bach_Xi_after(0),
                                                                                        hInvMassLambda_pi_daugh_Xi_after(0),
                                                                                        hInvMassLambda_prot_Xi_after(0),
                                                                                        hInvMassLambda_full_lambda_from_Xi_after(0),
                                                                                        hInvMassXi_sanityCheck_after(0),
                                                                                        hInvMassXi_Lamda_pi_daugh_after(0),
                                                                                        hInvMassXi_Lamda_prot_daugh_after(0),
                                                                                        hInvMassXi_Lamda_pi_bach_after(0),
                                                                                        hInvMassXi_Lamda_full_after(0),
                                                                                        hInvMassXi_Lamda_pi_no_correctLambdaMass(0),
                                                                                        hInvMassXi_Lamda_prot_no_correctLambdaMass(0),
                                                                                        hInvMassAntiLambda_sanityCheck_after(0),
                                                                                        hInvMassAntiLambda_pi_bach_Xi_after(0),
                                                                                        hInvMassAntiLambda_pi_daugh_Xi_after(0),
                                                                                        hInvMassAntiLambda_prot_Xi_after(0),
                                                                                        hInvMassAntiLambda_full_lambda_from_Xi_after(0),
                                                                                        hInvMassAntiXi_sanityCheck_after(0),
                                                                                        hInvMassAntiXi_AntiLamda_antipi_daugh_after(0),
                                                                                        hInvMassAntiXi_AntiLamda_antiprot_daugh_after(0),
                                                                                        hInvMassAntiXi_AntiLamda_antipi_bach_after(0),
                                                                                        hInvMassAntiXi_AntiLamda_full_after(0),
                                                                                        hInvMassAntiXi_AntiLamda_antipi_no_correctAntiLambdaMass(0),
                                                                                        hInvMassAntiXi_AntiLamda_antiprot_no_correctAntiLambdaMass(0),
                                                                                        fEvtCounterAfter(0),
                                                                                        tlInvMassPairClean(0),
                                                                                        tlCleanDecay(0),
                                                                                        hLambdaCleanedPartMassDiffToPDG_Decay(0),
                                                                                        hAntiLambdaCleanedPartMassDiffToPDG_Decay(0),
                                                                                        hXiCleanedPartMassDiffToPDG_Decay(0),
                                                                                        hAntiXiCleanedPartMassDiffToPDG_Decay(0),
                                                                                        hLambdaCleanedPartMass_Decay(0),
                                                                                        hAntiLambdaCleanedPartMass_Decay(0),
                                                                                        hXiCleanedPartMass_Decay(0),
                                                                                        hAntiXiCleanedPartMass_Decay(0),
                                                                                        tlCleanDecayAndDecay(0),
                                                                                        hLambdaCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                                        hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                                        hXiCleanedPartMassDiffToPDG_DecayDecay(0),
                                                                                        hAntiXiCleanedPartMassDiffToPDG_DecayDecay(0),    
                                                                                        hLambdaCleanedPartMass_DecayDecay(0),
                                                                                        hAntiLambdaCleanedPartMass_DecayDecay(0),
                                                                                        hXiCleanedPartMass_DecayDecay(0),
                                                                                        hAntiXiCleanedPartMass_DecayDecay(0),
                                                                                        tlCPA_PairClean_stats(0),
                                                                                        tlLambda_CPA_stats(0),
                                                                                        tlAntiLambda_CPA_stats(0),
                                                                                        tlXi_CPA_stats(0),
                                                                                        tlAntiXi_CPA_stats(0),
                                                                                        tlLambdaCPA_MC(0),
                                                                                        tlAntiLambdaCPA_MC(0),
                                                                                        tlXiCPA_MC(0),
                                                                                        tlAntiXiCPA_MC(0),
                                                                                        tlCPA_pT_Pairclean_CPA(0),
                                                                                        tlCPA_pT_Pairclean_InvMass(0),
                                                                                        h2_CPA_pt(0),
                                                                                        CPAPtBinningPrim_lambda(0),
                                                                                        CPAPtBinningMat_lambda(0),
                                                                                        CPAPtBinningSec_lambda(0),
                                                                                        CPAPtBinningCont_lambda(0),
                                                                                        CPAPtBinningPrim_lambda_dump(0),
                                                                                        CPAPtBinningMat_lambda_dump(0),
                                                                                        CPAPtBinningSec_lambda_dump(0),
                                                                                        CPAPtBinningCont_lambda_dump(0),
                                                                                        CPAPtBinningPrim_antilambda(0),
                                                                                        CPAPtBinningMat_antilambda(0),
                                                                                        CPAPtBinningSec_antilambda(0),
                                                                                        CPAPtBinningCont_antilambda(0),
                                                                                        CPAPtBinningPrim_xi(0),
                                                                                        CPAPtBinningMat_xi(0),
                                                                                        CPAPtBinningSec_xi(0),
                                                                                        CPAPtBinningCont_xi(0),
                                                                                        CPAPtBinningPrim_xi_dump(0),
                                                                                        CPAPtBinningMat_xi_dump(0),
                                                                                        CPAPtBinningSec_xi_dump(0),
                                                                                        CPAPtBinningCont_xi_dump(0),
                                                                                        CPAPtBinningPrim_antixi(0),
                                                                                        CPAPtBinningMat_antixi(0),
                                                                                        CPAPtBinningSec_antixi(0),
                                                                                        CPAPtBinningCont_antixi(0),
                                                                                        h2_pt_invMass(0)
{
    DefineOutput(1, TList::Class());    // Event Cuts
    DefineOutput(2, TList::Class());    // Lambda Track Cuts
    DefineOutput(3, TList::Class());    // Anti Lambda Track Cuts
    DefineOutput(4, TList::Class());    // Results - PairCleaner
    DefineOutput(5, TList::Class());    // QA Results
    DefineOutput(6, TList::Class());    // reconstruction from daugthers histograms AFTER PairCleaner
    DefineOutput(7, TList::Class());    // Results2 
    DefineOutput(8, TList::Class());    // Results3

    if (isMC)
    {
        DefineOutput(9, TList::Class());    // MC V0 - Lambda
        DefineOutput(10, TList::Class());    // MC V0 - AntiLambda
    }
}

AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()       // Destructor
{
// alle Objecte die zu einer TList hinzugefügt wurden mit TList::Add() werden automatisch vom TList Destructor aufgelöst
//
// das hier überall nach 'if' gefragt wird liegt daran, weil die Objecte erst in 'UserCreateOutputObjects' erzeugt werden
// und die nicht zwingend aufgerufen werden müssen, wenn was schiefgeht.
    if(fEvent)                  delete fEvent;
    if(fGTI)                    delete fGTI;
    if(fv0)                     delete fv0;
    if(fPairCleaner)            delete fPairCleaner;
    if(fPartColl)               delete fPartColl;
    if(fPairCleaner2)           delete fPairCleaner2;
    if(fPartColl2)              delete fPartColl2;
    if(fPartColl3)              delete fPartColl3;
    if(tlLambdaList)            delete tlLambdaList;
    if(tlAntiLambdaList)        delete tlAntiLambdaList;
    if(tlRecombination_before)  delete tlRecombination_before;
    if(tlRecombination_after)   delete tlRecombination_after;
    if(tlResultsQA)             delete tlResultsQA;
    if(tlResults)               delete tlResults;
    if(tlResults2)              delete tlResults2;
    if(tlResults3)              delete tlResults3;
    if(h2_CPA_pt)               delete[] h2_CPA_pt;
    if(h2_pt_invMass)       delete[] h2_pt_invMass;
}

void AliAnalysisTaskPOmegaPenne::UserCreateOutputObjects()
{   
    fEvent = new AliFemtoDreamEvent(true, ffullBlastQA, GetCollisionCandidates(), true);
    fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

    fGTI = new AliVTrack *[fTrackBufferSize];
    
    fEventCuts->InitQA();
    
 
    // Lambda Cutys    ###########
    if (!fLambdaV0Cuts){AliFatal("Track Cuts for Particle Lambda not set!");}
    fLambdaV0Cuts->Init();
    fLambdaV0Cuts->SetName("Lambda");
    // ##

    // AntiLambda Cutys    ###########
    if (!fAntiLambdaV0Cuts){AliFatal("Track Cuts for Particle AntiLambda not set!");}
    fAntiLambdaV0Cuts->Init();
    fAntiLambdaV0Cuts->SetName("AntiLambda");
    // ##

    // V0 Candidates
    fv0 = new AliFemtoDreamv0();
    fv0->                  SetUseMCInfo(fLambdaV0Cuts->GetIsMonteCarlo() || fAntiLambdaV0Cuts->GetIsMonteCarlo());
    fv0->GetPosDaughter()->SetUseMCInfo(fLambdaV0Cuts->GetIsMonteCarlo() || fAntiLambdaV0Cuts->GetIsMonteCarlo()); 
    fv0->GetNegDaughter()->SetUseMCInfo(fLambdaV0Cuts->GetIsMonteCarlo() || fAntiLambdaV0Cuts->GetIsMonteCarlo()); 
    fv0->SetPDGCode(3122);
    fv0->SetPDGDaughterPos(2212);
    fv0->SetPDGDaughterNeg(211);
    // ##

    // Xi Cuts    ###########    
    
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 6, false);
    fPartColl = new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());
    fPartColl2 = new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());
    fPartColl3 = new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());
    // ##

    /////////////////////////////
    // BEFORE Paircleaning histos
    /////////////////////////////
    tlRecombination_before = new TList();        // Lambda and Xi recombination statistic histogramms for interchanged daughters
    tlRecombination_before->SetName("Recombination_before_pairclean");
    tlRecombination_before->SetOwner(kTRUE);
    
    if (fmixBeforePC)
    {
    }
    // $$$ END - BEFORE Paircleaning $$$

    ////////////////////////////
    // AFTER Paircleaning histos / Lists
    ////////////////////////////
    tlRecombination_after = new TList();
    tlRecombination_after->SetName("Recombination_after_pairclean");
    tlRecombination_after->SetOwner(kTRUE);

    //////////////////////
    // Inv Mass PC   /////
    //////////////////////
    tlInvMassPairClean = new TList();
    tlInvMassPairClean->SetName("PairCleaner_Stats");
    tlInvMassPairClean->SetOwner(kTRUE);

    tlCleanDecay = new TList();
    tlCleanDecay->SetName("CleanDecay");
    tlCleanDecay->SetOwner(kTRUE);

    tlCleanDecayAndDecay = new TList();
    tlCleanDecayAndDecay->SetName("CleanDecayAndDecay");
    tlCleanDecayAndDecay->SetOwner(kTRUE);

    tlCPA_PairClean_stats = new TList();
    tlCPA_PairClean_stats->SetName("extendedPCstats");
    tlCPA_PairClean_stats->SetOwner(kTRUE);


    if (fmixAfterPC)
    {
    }
    ///////////////////////////////////////////////////
    // CPA Distributions pT binned
    /////////////////////////////////////////////////
    ////
    // folders
    tlLambda_CPA_stats = new TList();
    tlLambda_CPA_stats->SetName("CPAstats_beforeAndAfterCPA_Cleaning");
    tlLambda_CPA_stats->SetOwner(kTRUE);

    tlAntiLambda_CPA_stats = new TList();
    tlAntiLambda_CPA_stats->SetName("CPAstat_InvMassAndRandom_Cleaning");
    tlAntiLambda_CPA_stats->SetOwner(kTRUE);

    tlXi_CPA_stats = new TList();
    tlXi_CPA_stats->SetName("MC_Origin_lambda_invMassCleaner");
    tlXi_CPA_stats->SetOwner(kTRUE);

    tlAntiXi_CPA_stats = new TList();
    tlAntiXi_CPA_stats->SetName("MC_Origin_ANTIlambda_invMassCleaner");
    tlAntiXi_CPA_stats->SetOwner(kTRUE);

    tlCPA_PairClean_stats->Add(tlLambda_CPA_stats);
    tlCPA_PairClean_stats->Add(tlAntiLambda_CPA_stats);
    tlCPA_PairClean_stats->Add(tlXi_CPA_stats);
    tlCPA_PairClean_stats->Add(tlAntiXi_CPA_stats);
    // histos
    h2_CPA_pt = new TH2F*[20];
    // 
    h2_CPA_pt[0] = new TH2F("CPAPtBinning_lambda_before_Paircleaning","CPAPtBinning_lambda_before_Paircleaning", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[1] = new TH2F("CPAPtBinning_anti_lambda_before_Paircleaning","CPAPtBinning_lambda_after_DecayCleaning_GOOD_paritcles", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[2] = new TH2F("CPAPtBinning_lambda_after_DecayCleaning","CPAPtBinning_lambda_after_DecayCleaning", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[3] = new TH2F("CPAPtBinning_anti_lambda_after_DecayCleaning","CPAPtBinning_anti_lambda_after_DecayCleaning", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[4] = new TH2F("nothing_anymore","nothing_anymore", 8, 0.3, 4.3, 800, 0.987, 1.0);
    // 
    h2_CPA_pt[5] = new TH2F("CPAPtBinning_lambda_after_InvMassCleaning","CPAPtBinning_lambda_after_InvMassCleaning", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[6] = new TH2F("CPAPtBinning_anti_lambda_after_InvMassCleaning","CPAPtBinning_anti_lambda_after_InvMassCleaning", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[7] = new TH2F("CPAPtBinning_lambda_after_Cleaning_At_Random","CPAPtBinning_lambda_after_Cleaning_At_Random", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[8] = new TH2F("CPAPtBinning_anti_lambda_after_Cleaning_At_Random","CPAPtBinning_anti_lambda_after_Cleaning_At_Random", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[9] = new TH2F("nothing_anymore","nothing_anymore", 8, 0.3, 4.3, 800, 0.987, 1.0);
    // 
    h2_CPA_pt[10] = new TH2F("CPAPtBinning_Lambda_Primary_invMassCleaner","CPAPtBinning_Lambda_Primary_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[11] = new TH2F("CPAPtBinning_Lambda_Material_invMassCleaner","CPAPtBinning_Lambda_Material_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[12] = new TH2F("CPAPtBinning_Lambda_Secondary_invMassCleaner","CPAPtBinning_Lambda_Secondary_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[13] = new TH2F("CPAPtBinning_Lambda_Contamination_invMassCleaner","CPAPtBinning_Lambda_Contamination_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[14] = new TH2F("nix","nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    // 
    h2_CPA_pt[15] = new TH2F("CPAPtBinning_Anti_Lambda_Primary_invMassCleaner","CPAPtBinning_Anti_Lambda_Primary_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[16] = new TH2F("CPAPtBinning_Anti_Lambda_Material_invMassCleaner","CPAPtBinning_Anti_Lambda_Material_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[17] = new TH2F("CPAPtBinning_Anti_Lambda_Secondary_invMassCleaner","CPAPtBinning_Anti_Lambda_Secondary_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[18] = new TH2F("CPAPtBinning_Anti_Lambda_Contamination_invMassCleaner","CPAPtBinning_Anti_Lambda_Contamination_invMassCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_CPA_pt[19] = new TH2F("nix","nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    // Axis Label
    for (size_t i = 0; i < 20; i++)
    {
        h2_CPA_pt[i]->GetXaxis()->SetTitle("P_{T}");
        h2_CPA_pt[i]->GetYaxis()->SetTitle("CPA");
    }
    
    // Connect
    for (size_t i = 0; i < 5; i++)
    {
        tlLambda_CPA_stats->Add(h2_CPA_pt[i]);
    }
    for (size_t i = 5; i < 10; i++)
    {
        tlAntiLambda_CPA_stats->Add(h2_CPA_pt[i]);
    }
    for (size_t i = 10; i < 15; i++)
    {
        tlXi_CPA_stats->Add(h2_CPA_pt[i]);
    }
    for (size_t i = 15; i < 20; i++)
    {
        tlAntiXi_CPA_stats->Add(h2_CPA_pt[i]);
    }
    

        // Decay Diff To PDG Mass
    hLambdaCleanedPartMassDiffToPDG_Decay = new TH1F("LambdaCleanedParticleDifferenceToPDGMass", "Cleaned Lambda Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hAntiLambdaCleanedPartMassDiffToPDG_Decay = new TH1F("AntiLambdaCleanedParticleDifferenceToPDGMass", "Cleaned Anti Lambda Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hXiCleanedPartMassDiffToPDG_Decay = new TH1F("XiCleanedParticleDifferenceToPDGMass", "Cleaned Xi Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hAntiXiCleanedPartMassDiffToPDG_Decay = new TH1F("AntiXiCleanedParticleDifferenceToPDGMass", "Cleaned Anti Xi Difference To PDG Mass", 300, -3.0, 3.0);

        // Decay Mass
    hLambdaCleanedPartMass_Decay = new TH1F("LambdaCleanedParticleDifferenceToPDGMass", "Lambda Cleaned Particle Mass", 400, 1.00, 1.20);
    hAntiLambdaCleanedPartMass_Decay = new TH1F("AntiLambdaCleanedParticleDifferenceToPDGMass", "Anti Lambda Cleaned Mass", 400, 1.00, 1.20);
    hXiCleanedPartMass_Decay = new TH1F("CleanedXiMass", "Cleaned Xi Particle Mass", 200, 1.321-0.06, 1.321+0.06);
    hAntiXiCleanedPartMass_Decay = new TH1F("CleanedAntiXiMass", "Cleaned Anti Xi Particle Mass", 200, 1.321-0.06, 1.321+0.06);

        // DecayAndDecay Diff To PDG Mass
    hLambdaCleanedPartMassDiffToPDG_DecayDecay = new TH1F("LambdaCleanedParticleDifferenceToPDGMass", "Cleaned Lambda Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay = new TH1F("AntiLambdaCleanedParticleDifferenceToPDGMass", "Cleaned Anti Lambda Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hXiCleanedPartMassDiffToPDG_DecayDecay = new TH1F("XiCleanedParticleDifferenceToPDGMass", "Cleaned Xi Mass Difference To PDG Mass", 300, -3.0, 3.0);
    hAntiXiCleanedPartMassDiffToPDG_DecayDecay = new TH1F("AntiXiCleanedParticleDifferenceToPDGMass", "Cleaned Anti Xi Difference To PDG Mass", 300, -3.0, 3.0);

        // DecayAndDecay Mass                             
    hLambdaCleanedPartMass_DecayDecay = new TH1F("LambdaCleanedParticleMass", "Lambda Cleaned Particle Mass", 800, 1.00, 1.20);
    hAntiLambdaCleanedPartMass_DecayDecay = new TH1F("AntiLambdaCleanedParticleMass", "Anti Lambda Cleaned Particle Mass", 800, 1.00, 1.20);
    hXiCleanedPartMass_DecayDecay = new TH1F("XiCleanedParticleMass", "Xi Cleaned Particle Mass", 200, 1.321-0.06, 1.321+0.06);
    hAntiXiCleanedPartMass_DecayDecay = new TH1F("AntiXiCleanedParticleMass", "Anti Cleaned Particle Mass", 200, 1.321-0.06, 1.321+0.06);

    /////////////////////////////////////////
    /////// MC CPA pt Binning AFTER paircleaing  -- checking becoz Femtodream Histos are filled before Paircleaing
    ////////////////////////////////////////
    tlCPA_pT_Pairclean_CPA = new TList();
    tlCPA_pT_Pairclean_CPA->SetName("CPA_pT_Pairclean_CPA");
    tlCPA_pT_Pairclean_CPA->SetOwner(kTRUE);
    // Lambda
    tlLambdaCPA_MC = new TList();
    tlLambdaCPA_MC->SetName("MC_LambdaCPA");
    tlLambdaCPA_MC->SetOwner(kTRUE);

    CPAPtBinningPrim_lambda = new TH2F("CPAPtBinningPrim_lambda_atRandomCleaner", "CPAPtBinningPrim_lambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningPrim_lambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_lambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_lambda = new TH2F("CPAPtBinningMat_lambda_atRandomCleaner", "CPAPtBinningMat_lambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningMat_lambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_lambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_lambda = new TH2F("CPAPtBinningSec_lambda_atRandomCleaner", "CPAPtBinningSec_lambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningSec_lambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_lambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_lambda = new TH2F("CPAPtBinningCont_lambda_atRandomCleaner", "CPAPtBinningCont_lambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningCont_lambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_lambda->GetYaxis()->SetTitle("CPA");

    // Lambda > dumps - TList tlLambdaCPA_MC
    CPAPtBinningPrim_lambda_dump = new TH2F("CPAPtBinningPrim_lambda_dump", "nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningPrim_lambda_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_lambda_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_lambda_dump = new TH2F("CPAPtBinningMat_lambda_dump", "nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningMat_lambda_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_lambda_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_lambda_dump = new TH2F("CPAPtBinningSec_lambda_dump", "nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningSec_lambda_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_lambda_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_lambda_dump = new TH2F("CPAPtBinningCont_lambda_dump", "nix", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningCont_lambda_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_lambda_dump->GetYaxis()->SetTitle("CPA");

    // Anti-Lambda
    tlAntiLambdaCPA_MC = new TList();
    tlAntiLambdaCPA_MC->SetName("AntiLambdaCPA_MC");
    tlAntiLambdaCPA_MC->SetOwner(kTRUE);

    CPAPtBinningPrim_antilambda = new TH2F("CPAPtBinningPrim_antilambda_atRandomCleaner", "CPAPtBinningPrim_antilambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningPrim_antilambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_antilambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_antilambda = new TH2F("CPAPtBinningMat_antilambda_atRandomCleaner", "CPAPtBinningMat_antilambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningMat_antilambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_antilambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_antilambda = new TH2F("CPAPtBinningSec_antilambda_atRandomCleaner", "CPAPtBinningSec_antilambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningSec_antilambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_antilambda->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_antilambda = new TH2F("CPAPtBinningCont_antilambda_atRandomCleaner", "CPAPtBinningCont_antilambda_atRandomCleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    CPAPtBinningCont_antilambda->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_antilambda->GetYaxis()->SetTitle("CPA");

    // Xi
    tlXiCPA_MC = new TList();
    tlXiCPA_MC->SetName("XiCPA_MC");
    tlXiCPA_MC->SetOwner(kTRUE);

    CPAPtBinningPrim_xi = new TH2F("CPAPtBinningPrim_xi", "CPAPtBinningPrim_xi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningPrim_xi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_xi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_xi = new TH2F("CPAPtBinningMat_xi", "CPAPtBinningMat_xi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningMat_xi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_xi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_xi = new TH2F("CPAPtBinningSec_xi", "CPAPtBinningSec_xi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningSec_xi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_xi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_xi = new TH2F("CPAPtBinningCont_xi", "CPAPtBinningCont_xi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningCont_xi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_xi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningPrim_xi_dump = new TH2F("CPAPtBinningPrim_xi_dump", "CPAPtBinningPrim_xi_dump", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningPrim_xi_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_xi_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_xi_dump = new TH2F("CPAPtBinningMat_xi_dump", "CPAPtBinningMat_xi_dump", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningMat_xi_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_xi_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_xi_dump = new TH2F("CPAPtBinningSec_xi_dump", "CPAPtBinningSec_xi_dump", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningSec_xi_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_xi_dump->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_xi_dump = new TH2F("CPAPtBinningCont_xi_dump", "CPAPtBinningCont_xi_dump", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningCont_xi_dump->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_xi_dump->GetYaxis()->SetTitle("CPA");

    // Anti-Xi
    tlAntiXiCPA_MC = new TList();
    tlAntiXiCPA_MC->SetName("AntiXiCPA_MC");
    tlAntiXiCPA_MC->SetOwner(kTRUE);

    CPAPtBinningPrim_antixi = new TH2F("CPAPtBinningPrim_antixi", "CPAPtBinningPrim_antixi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningPrim_antixi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningPrim_antixi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningMat_antixi = new TH2F("CPAPtBinningMat_antixi", "CPAPtBinningMat_antixi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningMat_antixi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningMat_antixi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningSec_antixi = new TH2F("CPAPtBinningSec_antixi", "CPAPtBinningSec_antixi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningSec_antixi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningSec_antixi->GetYaxis()->SetTitle("CPA");

    CPAPtBinningCont_antixi = new TH2F("CPAPtBinningCont_antixi", "CPAPtBinningCont_antixi", 8, 0.3, 4.3, 800, 1.00, 1.2);
    CPAPtBinningCont_antixi->GetXaxis()->SetTitle("P_{T}");
    CPAPtBinningCont_antixi->GetYaxis()->SetTitle("CPA");

    
    h2_pt_invMass = new TH2F*[64];

// [0-15] for MC    -   pT vs CPA stats     = before and CPA cleaner
    h2_pt_invMass[0] = new TH2F("CPAPtBinningPrim_lambda_before", "CPAPtBinningPrim_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[1] = new TH2F("CPAPtBinningMat_lambda_before", "CPAPtBinningMat_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[2] = new TH2F("CPAPtBinningSec_lambda_before", "CPAPtBinningSec_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[3] = new TH2F("CPAPtBinningCont_lambda_before", "CPAPtBinningCont_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[4] = new TH2F("CPAPtBinningPrim_Anti_lambda_before", "CPAPtBinningPrim_Anti_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[5] = new TH2F("CPAPtBinningMat_Anti_lambda_before", "CPAPtBinningMat_Anti_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[6] = new TH2F("CPAPtBinningSec_Anti_lambda_before", "CPAPtBinningSec_Anti_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[7] = new TH2F("CPAPtBinningCont_Anti_lambda_before", "CPAPtBinningCont_Anti_lambda_before", 8, 0.3, 4.3, 800, 0.987, 1.0);

    h2_pt_invMass[8] = new TH2F("CPAPtBinningPrim_Lambda_CPAcleaner", "CPAPtBinningPrim_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[9] = new TH2F("CPAPtBinningMat_Lambda_CPAcleaner", "CPAPtBinningMat_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[10] = new TH2F("CPAPtBinningSec_Lambda_CPAcleaner", "CPAPtBinningSec_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[11] = new TH2F("CPAPtBinningCont_Lambda_CPAcleaner", "CPAPtBinningCont_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[12] = new TH2F("CPAPtBinningPrim_Anti_Lambda_CPAcleaner", "CPAPtBinningPrim_Anti_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[13] = new TH2F("CPAPtBinningMat_Anti_Lambda_CPAcleaner", "CPAPtBinningMat_Anti_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[14] = new TH2F("CPAPtBinningSec_Anti_Lambda_CPAcleaner", "CPAPtBinningSec_Anti_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);
    h2_pt_invMass[15] = new TH2F("CPAPtBinningCont_Anti_Lambda_CPAcleaner", "CPAPtBinningCont_Anti_Lambda_CPAcleaner", 8, 0.3, 4.3, 800, 0.987, 1.0);

// [16-31] Lambda pT vs invariant Mass stats
    h2_pt_invMass[16] = new TH2F("Lambda_Invariant_Mass_pT_before", "Lambda_Invariant_Mass_pT_before", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[17] = new TH2F("Lambda_Invariant_Mass_pT_Clean_Decay_CPA", "Lambda_Invariant_Mass_pT_Clean_Decay_CPA", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[18] = new TH2F("Lambda_Invariant_Mass_pT_Clean_Decay_Inv_Mass", "Lambda_Invariant_Mass_pT_Clean_Decay_Inv_Mass", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[19] = new TH2F("Lambda_Invariant_Mass_pT_Clean_Decay_At_Random", "Lambda_Invariant_Mass_pT_Clean_Decay_At_Random", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[20] = new TH2F("Anti_Lambda_Invariant_Mass_pT_before", "Anti_Lambda_Invariant_Mass_pT_before", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[21] = new TH2F("Anti_Lambda_Invariant_Mass_pT_Clean_Decay_CPA", "Anti_Lambda_Invariant_Mass_pT_Clean_Decay_CPA", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[22] = new TH2F("Anti_Lambda_Invariant_Mass_pT_Clean_Decay_Inv_Mass", "Anti_Lambda_Invariant_Mass_pT_Clean_Decay_Inv_Mass", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[23] = new TH2F("Anti_Lambda_Invariant_Mass_pT_Clean_Decay_At_Random", "Anti_Lambda_Invariant_Mass_pT_Clean_Decay_At_Random", 8, 0.3, 4.3, 800, 1.00, 1.2);

    h2_pt_invMass[24] = new TH2F("Lambda_Invariant_Mass_kSTAR_before", "Lambda_Invariant_Mass_kSTAR_before", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[25] = new TH2F("Lambda_Invariant_Mass_kSTAR_Clean_Decay_CPA", "Lambda_Invariant_Mass_kSTAR_Clean_Decay_CPA", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[26] = new TH2F("Lambda_Invariant_Mass_kSTAR_Clean_Decay_Inv_Mass", "Lambda_Invariant_Mass_kSTAR_Clean_Decay_Inv_Mass", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[27] = new TH2F("Lambda_Invariant_Mass_kSTAR_Clean_Decay_At_Random", "Lambda_Invariant_Mass_kSTAR_Clean_Decay_At_Random", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[28] = new TH2F("Anti_Lambda_Invariant_Mass_kSTAR_before", "Anti_Lambda_Invariant_Mass_kSTAR_before", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[29] = new TH2F("Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_CPA", "Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_CPA", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[30] = new TH2F("Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_Inv_Mass", "Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_Inv_Mass", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);
    h2_pt_invMass[31] = new TH2F("Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_At_Random", "Anti_Lambda_Invariant_Mass_kSTAR_Clean_Decay_At_Random", 512, 1.07568, 1.15568, 1000, 0.00, 1.00);

// [32-64] for MC  - pT vs invMass - Origin Binning
    // before
    h2_pt_invMass[32] = new TH2F("MC_Lambda_Invariant_Mass_PRIM_beforePC", "MC_Lambda_Invariant_Mass_PRIM_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[33] = new TH2F("MC_Lambda_Invariant_Mass_MAT_beforePC", "MC_Lambda_Invariant_Mass_MAT_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[34] = new TH2F("MC_Lambda_Invariant_Mass_SEC_beforePC", "MC_Lambda_Invariant_Mass_SEC_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[35] = new TH2F("MC_Lambda_Invariant_Mass_CONT_beforePC", "MC_Lambda_Invariant_Mass_CONT_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[36] = new TH2F("MC_ANTILambda_Invariant_Mass_PRIM_beforePC", "MC_ANTILambda_Invariant_Mass_PRIM_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[37] = new TH2F("MC_ANTILambda_Invariant_Mass_MAT_beforePC", "MC_ANTILambda_Invariant_Mass_MAT_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[38] = new TH2F("MC_ANTILambda_Invariant_Mass_SEC_beforePC", "MC_ANTILambda_Invariant_Mass_SEC_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[39] = new TH2F("MC_ANTILambda_Invariant_Mass_CONT_beforePC", "MC_ANTILambda_Invariant_Mass_CONT_beforePC", 8, 0.3, 4.3, 800, 1.00, 1.2);

    // CPA cleaning
    h2_pt_invMass[40] = new TH2F("MC_Lambda_Invariant_Mass_PRIM_CPACleaning", "MC_Lambda_Invariant_Mass_PRIM_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[41] = new TH2F("MC_Lambda_Invariant_Mass_MAT_CPACleaning", "MC_Lambda_Invariant_Mass_MAT_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[42] = new TH2F("MC_Lambda_Invariant_Mass_SEC_CPACleaning", "MC_Lambda_Invariant_Mass_SEC_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[43] = new TH2F("MC_Lambda_Invariant_Mass_CONT_CPACleaning", "MC_Lambda_Invariant_Mass_CONT_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[44] = new TH2F("MC_ANTILambda_Invariant_Mass_PRIM_CPACleaning", "MC_ANTILambda_Invariant_Mass_PRIM_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[45] = new TH2F("MC_ANTILambda_Invariant_Mass_MAT_CPACleaning", "MC_ANTILambda_Invariant_Mass_MAT_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[46] = new TH2F("MC_ANTILambda_Invariant_Mass_SEC_CPACleaning", "MC_ANTILambda_Invariant_Mass_SEC_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[47] = new TH2F("MC_ANTILambda_Invariant_Mass_CONT_CPACleaning", "MC_ANTILambda_Invariant_Mass_CONT_CPACleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);

    // inv mass cleaning
    h2_pt_invMass[48] = new TH2F("MC_Lambda_Invariant_Mass_PRIM_invMassCleaning", "MC_Lambda_Invariant_Mass_PRIM_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[49] = new TH2F("MC_Lambda_Invariant_Mass_MAT_invMassCleaning", "MC_Lambda_Invariant_Mass_MAT_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[50] = new TH2F("MC_Lambda_Invariant_Mass_SEC_invMassCleaning", "MC_Lambda_Invariant_Mass_SEC_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[51] = new TH2F("MC_Lambda_Invariant_Mass_CONT_invMassCleaning", "MC_Lambda_Invariant_Mass_CONT_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[52] = new TH2F("MC_ANTILambda_Invariant_Mass_PRIM_invMassCleaning", "MC_ANTILambda_Invariant_Mass_PRIM_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[53] = new TH2F("MC_ANTILambda_Invariant_Mass_MAT_invMassCleaning", "MC_ANTILambda_Invariant_Mass_MAT_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[54] = new TH2F("MC_ANTILambda_Invariant_Mass_SEC_invMassCleaning", "MC_ANTILambda_Invariant_Mass_SEC_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[55] = new TH2F("MC_ANTILambda_Invariant_Mass_CONT_invMassCleaning", "MC_ANTILambda_Invariant_Mass_CONT_invMassCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);

    // random cleaning 
    h2_pt_invMass[56] = new TH2F("MC_Lambda_Invariant_Mass_PRIM_atRandomCleaning", "MC_Lambda_Invariant_Mass_PRIM_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[57] = new TH2F("MC_Lambda_Invariant_Mass_MAT_atRandomCleaning", "MC_Lambda_Invariant_Mass_MAT_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[58] = new TH2F("MC_Lambda_Invariant_Mass_SEC_atRandomCleaning", "MC_Lambda_Invariant_Mass_SEC_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[59] = new TH2F("MC_Lambda_Invariant_Mass_CONT_atRandomCleaning", "MC_Lambda_Invariant_Mass_CONT_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[60] = new TH2F("MC_ANTILambda_Invariant_Mass_PRIM_atRandomCleaning", "MC_ANTILambda_Invariant_Mass_PRIM_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[61] = new TH2F("MC_ANTILambda_Invariant_Mass_MAT_atRandomCleaning", "MC_ANTILambda_Invariant_Mass_MAT_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[62] = new TH2F("MC_ANTILambda_Invariant_Mass_SEC_atRandomCleaning", "MC_ANTILambda_Invariant_Mass_SEC_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);
    h2_pt_invMass[63] = new TH2F("MC_ANTILambda_Invariant_Mass_CONT_atRandomCleaning", "MC_ANTILambda_Invariant_Mass_CONT_atRandomCleaning", 8, 0.3, 4.3, 800, 1.00, 1.2);


    for(int i = 16; i < 24; i++)
    {
        h2_pt_invMass[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        h2_pt_invMass[i]->GetYaxis()->SetTitle("M_{p#pi} (GeV/c^{2})");
    }
    for(int i = 24; i < 32; i++)
    {
        h2_pt_invMass[i]->GetYaxis()->SetTitle("k* (GeV/c)");
        h2_pt_invMass[i]->GetXaxis()->SetTitle("M_{p#pi} (GeV/c^{2})");
    }
        //
        // Connect Histogramms to Lists
        //
    tlCPA_pT_Pairclean_InvMass = new TList();
    tlCPA_pT_Pairclean_InvMass->SetName("InvMassPtAndkStar_AllThreePaircleaner");
    tlCPA_pT_Pairclean_InvMass->SetOwner(kTRUE);
    
    // Lists
    tlInvMassPairClean->Add(tlCleanDecay);
    tlInvMassPairClean->Add(tlCleanDecayAndDecay);
    tlInvMassPairClean->Add(tlCPA_PairClean_stats);

    tlCPA_pT_Pairclean_CPA->Add(tlLambdaCPA_MC);
    tlCPA_pT_Pairclean_CPA->Add(tlAntiLambdaCPA_MC);
    tlCPA_pT_Pairclean_CPA->Add(tlXiCPA_MC);
    tlCPA_pT_Pairclean_CPA->Add(tlAntiXiCPA_MC);
    
    if (fLambdaV0Cuts->GetIsMonteCarlo() || fAntiLambdaV0Cuts->GetIsMonteCarlo())
    {
        for (int i = 0; i < 64; i++)
        {
            tlCPA_pT_Pairclean_InvMass->Add(h2_pt_invMass[i]);
        }
    }
    else
    {
        for (int i = 0; i < 32; i++)
        {
            tlCPA_pT_Pairclean_InvMass->Add(h2_pt_invMass[i]);
        }
    }
    
    tlCPA_PairClean_stats->Add(tlCPA_pT_Pairclean_CPA);
    tlCPA_PairClean_stats->Add(tlCPA_pT_Pairclean_InvMass);

    // Histos

    // decay Cleaning
    tlCleanDecay->Add(hLambdaCleanedPartMassDiffToPDG_Decay);
    tlCleanDecay->Add(hAntiLambdaCleanedPartMassDiffToPDG_Decay);
    tlCleanDecay->Add(hXiCleanedPartMassDiffToPDG_Decay);
    tlCleanDecay->Add(hAntiXiCleanedPartMassDiffToPDG_Decay);
    
    tlCleanDecay->Add(hLambdaCleanedPartMass_Decay);
    tlCleanDecay->Add(hAntiLambdaCleanedPartMass_Decay);
    tlCleanDecay->Add(hXiCleanedPartMass_Decay);
    tlCleanDecay->Add(hAntiXiCleanedPartMass_Decay);
    
    // decayAndDecay Cleaning
    tlCleanDecayAndDecay->Add(hLambdaCleanedPartMassDiffToPDG_DecayDecay);
    tlCleanDecayAndDecay->Add(hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay);
    tlCleanDecayAndDecay->Add(hXiCleanedPartMassDiffToPDG_DecayDecay);
    tlCleanDecayAndDecay->Add(hAntiXiCleanedPartMassDiffToPDG_DecayDecay);
    
    tlCleanDecayAndDecay->Add(hLambdaCleanedPartMass_DecayDecay);
    tlCleanDecayAndDecay->Add(hAntiLambdaCleanedPartMass_DecayDecay);
    tlCleanDecayAndDecay->Add(hXiCleanedPartMass_DecayDecay);
    tlCleanDecayAndDecay->Add(hAntiXiCleanedPartMass_DecayDecay);

    // CPA MC Binning - after PC
    tlLambdaCPA_MC->Add(CPAPtBinningPrim_lambda);
    tlLambdaCPA_MC->Add(CPAPtBinningMat_lambda);
    tlLambdaCPA_MC->Add(CPAPtBinningSec_lambda);
    tlLambdaCPA_MC->Add(CPAPtBinningCont_lambda);

    tlLambdaCPA_MC->Add(CPAPtBinningPrim_lambda_dump);
    tlLambdaCPA_MC->Add(CPAPtBinningMat_lambda_dump);
    tlLambdaCPA_MC->Add(CPAPtBinningSec_lambda_dump);
    tlLambdaCPA_MC->Add(CPAPtBinningCont_lambda_dump);

    tlAntiLambdaCPA_MC->Add(CPAPtBinningPrim_antilambda);
    tlAntiLambdaCPA_MC->Add(CPAPtBinningMat_antilambda);
    tlAntiLambdaCPA_MC->Add(CPAPtBinningSec_antilambda);
    tlAntiLambdaCPA_MC->Add(CPAPtBinningCont_antilambda);

    tlXiCPA_MC->Add(CPAPtBinningPrim_xi);
    tlXiCPA_MC->Add(CPAPtBinningMat_xi);
    tlXiCPA_MC->Add(CPAPtBinningSec_xi);
    tlXiCPA_MC->Add(CPAPtBinningCont_xi);

    tlXiCPA_MC->Add(CPAPtBinningPrim_xi_dump);
    tlXiCPA_MC->Add(CPAPtBinningMat_xi_dump);
    tlXiCPA_MC->Add(CPAPtBinningSec_xi_dump);
    tlXiCPA_MC->Add(CPAPtBinningCont_xi_dump);

    tlAntiXiCPA_MC->Add(CPAPtBinningPrim_antixi);
    tlAntiXiCPA_MC->Add(CPAPtBinningMat_antixi);
    tlAntiXiCPA_MC->Add(CPAPtBinningSec_antixi);
    tlAntiXiCPA_MC->Add(CPAPtBinningCont_antixi);

        // connect to output List tlRecombination_after
    tlRecombination_after->Add(tlInvMassPairClean);

    // weird stuff
    // kStarXiLambda_unchanged = new TH1F("kStarXiLambda_unchanged", "kStarXiLambda_unchanged", 1200, 0.00, 1.200);
    // kStarXiLambda_changed = new TH1F("kStarXiLambda_changed", "kStarXiLambda_changed", 1200, 0.00, 1.20);
    // tlRecombination_after->Add(kStarXiLambda_unchanged);
    // tlRecombination_after->Add(kStarXiLambda_changed);
    // kStarAntiXiAntiLambda_unchanged = new TH1F("kStarAntiXiAntiLambda_unchanged", "kStarAntiXiAntiLambda_unchanged", 1200, 0.00, 1.20);
    // kStarAntiXiAntiLambda_changed = new TH1F("kStarAntiXiAntiLambda_changed", "kStarAntiXiAntiLambda_changed", 1200, 0.00, 1.20);
    // tlRecombination_after->Add(kStarAntiXiAntiLambda_unchanged);
    // tlRecombination_after->Add(kStarAntiXiAntiLambda_changed);

    ///////////////////////////////////////
    // Connect Cuts to OutputContainers ///
    ///////////////////////////////////////

    tlResultsQA = new TList();
    tlResultsQA->SetName("ResultsQA");
    tlResultsQA->SetOwner();

    // tlResultsQA2 = new TList();
    // tlResultsQA2->SetName("ResultsQA2");
    // tlResultsQA2->SetOwner();

    
    if(!fEventCuts->GetMinimalBooking())
    {
        tlEventCuts             = fEventCuts->GetHistList();
    } else
    {
        tlEventCuts = new TList();
        tlEventCuts->SetName("EventCuts");
        tlEventCuts->SetOwner();
    }
    
    tlLambdaList            = fLambdaV0Cuts->GetQAHists();
    tlAntiLambdaList        = fAntiLambdaV0Cuts->GetQAHists();

// initialize and connect RESULTS
    // RESULTS 1 /////////////////////////////////////////////
    if (fConfig->GetUseEventMixing())
    {
        tlResults = fPartColl->GetHistList();
        tlResults->SetName("ResultsCPACleaner");
        if (!fConfig->GetMinimalBookingME())
        {
            tlResultsQA->Add(fPartColl->GetQAList());
            // tlResultsQA->Add(fPartColl2->GetQAList());
            // tlResultsQA->Add(fPartColl3->GetQAList());
            tlResultsQA->Add(fPairCleaner->GetHistList());
        }
    }
    else
    {
        tlResults = new TList();
        tlResults->SetOwner();
        tlResults->SetName("Results");
    }
    // RESULTS 2 ////////////////////////////////////////////
    if (fConfig->GetUseEventMixing())
    {
        tlResults2 = fPartColl2->GetHistList();
        tlResults2->SetName("ResultsInvMassCleaner");
        if (!fConfig->GetMinimalBookingME())
        {
            tlResults2->Add(fPartColl2->GetQAList());
        }
    }
    else
    {
        tlResults2 = new TList();
        tlResults2->SetOwner();
        tlResults2->SetName("Results2");
    }
    // RESULTS 3 ///////////////////////////////////////////
    if (fConfig->GetUseEventMixing())
    {
        tlResults3 = fPartColl3->GetHistList();
        tlResults3->SetName("ResultsRandomCleaner");
        if (!fConfig->GetMinimalBookingME())
        {
            tlResults3->Add(fPartColl3->GetQAList());
        }
    }
    else
    {
        tlResults3 = new TList();
        tlResults3->SetOwner();
        tlResults3->SetName("Results3");
    }

    PostData(1, tlEventCuts);           // cuts keeping Lambda
    PostData(2, tlLambdaList);
    PostData(3, tlAntiLambdaList);
    PostData(4, tlResults);             // only Decay Cleaning
    PostData(5, tlResultsQA);
    PostData(6, tlRecombination_after);
    PostData(7, tlResults2);     
    PostData(8, tlResults3);     
    
    
    /////////////////
    // Monte Carlo //
    /////////////////
    if (fLambdaV0Cuts->GetIsMonteCarlo())
    {
        if (!fLambdaV0Cuts->GetMinimalBooking())
        {
                tlLambdaMC = fLambdaV0Cuts->GetMCQAHists();
        }
        else
        {
            tlLambdaMC = new TList();
            tlLambdaMC->SetName("v0CutsMC");
            tlLambdaMC->SetOwner();
        }
        PostData(9, tlLambdaMC);
    }
    if (fAntiLambdaV0Cuts->GetIsMonteCarlo())
    {
        if (!fAntiLambdaV0Cuts->GetMinimalBooking())
        {
            tlAntiLambdaMC = fAntiLambdaV0Cuts->GetMCQAHists();
        }
        else{
            tlAntiLambdaMC = new TList();
            tlAntiLambdaMC->SetName("Antiv0CutsMC");
            tlAntiLambdaMC->SetOwner();
        }
        PostData(10, tlAntiLambdaMC);
    }
}         

//  #######################################################################
//  #######################################################################
//  #######################################################################
//  ##
//  ##                     USER EXEC
//  ##
//  ##
//  #######################################################################
//  #######################################################################
//  #######################################################################
//  #######################################################################

// static int genericCounter = 1;
// static int multsOfHundred = 0;


void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{

    VEvent = fInputEvent;
    
    if (!fInputEvent)
    {
        AliWarning("No Input VEvent");
        return;
    }

    fEvent->SetEvent(fInputEvent);
    if (fEventCuts->isSelected(fEvent))
    {
        ResetGlobalTrackReference();
        for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
        {
            VTrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
            if (!VTrack)
            {
                AliFatal("No Standard AOD");
                return;
            }
            StoreGlobalTrackReference(VTrack);
        }
        // timer_event_selection_end = std::chrono::high_resolution_clock::now();
        vLambda.clear();
        vAntiLambda.clear();
    
        // irgendwie benötigt um GetV0s() und GetCascade() zu holen
        AliAODEvent *aodEvent = dynamic_cast<AliAODEvent *>(fInputEvent); // caste input event auf ein AODEvent

        //###########################################
        //#
        //# Particle Selections
        //#
        //###########################################
        // ## Lambda Selection ## keep Lambdas
        fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);


        for (int iv0 = 0; iv0 < dynamic_cast<TClonesArray *>(aodEvent->GetV0s())->GetEntriesFast(); ++iv0)
        {
            AliAODv0 *v0 = aodEvent->GetV0(iv0);
            
            fv0->Setv0(fInputEvent, v0);

            if (fLambdaV0Cuts->isSelected(fv0))
            {
                vLambda.push_back(*fv0);
            }
            if (fAntiLambdaV0Cuts->isSelected(fv0))
            {
                vAntiLambda.push_back(*fv0);
            }
        }

        //  ######################################################################
        //  ##
        //  ##
        //  ##
        //  ##                     PAIRCLEANING
        //  ##
        //  ##
        //  ##
        //  ##
        //  ##
        //  #######################################################################
        fPairCleaner ->ResetArray();
        
// MC
        if(fLambdaV0Cuts->GetIsMonteCarlo())
        {
            for(auto it : vLambda)
            {
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                    h2_pt_invMass[0]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[32]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                    h2_pt_invMass[1]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[33]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                    h2_pt_invMass[2]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[34]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                    h2_pt_invMass[3]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[35]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
            }
            for(auto it : vAntiLambda)
            {
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                    h2_pt_invMass[4]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[36]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                    h2_pt_invMass[5]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[37]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                    h2_pt_invMass[6]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[38]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                    h2_pt_invMass[7]->Fill(it.GetPt(), it.GetCPA()); 
                    h2_pt_invMass[39]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
            }
        }

// non-MC
        for(size_t i = 0; i < vLambda.size(); i++)
        {
            h2_CPA_pt[0]->Fill(vLambda[i].GetPt(), vLambda[i].GetCPA());
            h2_pt_invMass[16]->Fill(vLambda[i].GetPt(), CalculateInvMassLambda(vLambda[i], false));
                // h2_pt_invMass[24]->Fill(CalculateInvMassLambda(vLambda[i], false), RelativePairMomentum(vLambda[i], 3122, vLambda[j], 3122));
            vLambda[i].SetUse(true);
        }

        for(size_t i = 0; i < vAntiLambda.size(); i++)
        {
            h2_CPA_pt[1]->Fill(vAntiLambda[i].GetPt(), vAntiLambda[i].GetCPA());
            h2_pt_invMass[20]->Fill(vAntiLambda[i].GetPt(), CalculateInvMassLambda(vAntiLambda[i], true));
                // h2_pt_invMass[28]->Fill(CalculateInvMassLambda(vAntiLambda[i], true), RelativePairMomentum(vAntiLambda[i], -3122, vAntiLambda[j], -3122));
            vAntiLambda[i].SetUse(true);
        }
        
        

        // Part Collection #1   -   Clean Decay CPA
        fPairCleaner->CleanDecay(&vLambda, 0);
        fPairCleaner->CleanDecay(&vAntiLambda, 1);
        fPairCleaner->StoreParticle(vLambda);
        fPairCleaner->StoreParticle(vAntiLambda);
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); 


        if(fLambdaV0Cuts->GetIsMonteCarlo())
        {
// Origin binning for FemtoDreamPaircleaning
            for(auto it : vLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        h2_pt_invMass[8]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[40]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        h2_pt_invMass[9]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[41]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        h2_pt_invMass[10]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[42]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        h2_pt_invMass[11]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[43]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                }
            }
            for(auto it : vAntiLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        h2_pt_invMass[12]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[44]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                        }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        h2_pt_invMass[13]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[45]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                        }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        h2_pt_invMass[14]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[46]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                        }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        h2_pt_invMass[15]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[47]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                        }
                }
            }
        }


        for(size_t i = 0; i < vLambda.size(); i++)
        {
            if (vLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[2]->Fill(vLambda[i].GetPt(), vLambda[i].GetCPA());
                h2_pt_invMass[17]->Fill(vLambda[i].GetPt(), CalculateInvMassLambda(vLambda[i], false));
                    // h2_pt_invMass[25]->Fill(CalculateInvMassLambda(vLambda[i], false), RelativePairMomentum(vLambda[i], 3122, vLambda[j], 3122));
            }
            vLambda[i].SetUse(true);
        }
        for(size_t i = 0; i < vAntiLambda.size(); i++)
        {
            if (vAntiLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[3]->Fill(vAntiLambda[i].GetPt(), vAntiLambda[i].GetCPA());
                h2_pt_invMass[21]->Fill(vAntiLambda[i].GetPt(), CalculateInvMassLambda(vAntiLambda[i], true));
                    // h2_pt_invMass[29]->Fill(CalculateInvMassLambda(vAntiLambda[i], true), RelativePairMomentum(vAntiLambda[i], 3122, vAntiLambda[j], 3122));
            }
            vAntiLambda[i].SetUse(true);
        }

        
        // Part Collection #2   -   Clean Decay Inv Mass
        fPairCleaner->ResetArray();
        fPairCleaner->CleanDecayInvMass(&vLambda, 3122, 2);
        fPairCleaner->CleanDecayInvMass(&vAntiLambda, 3122, 3);
        fPairCleaner->StoreParticle(vLambda);
        fPairCleaner->StoreParticle(vAntiLambda);
        fPartColl2->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); 
        
        if(fLambdaV0Cuts->GetIsMonteCarlo())
        {
            for(auto it : vLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        h2_CPA_pt[10]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[48]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        h2_CPA_pt[11]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[49]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        h2_CPA_pt[12]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[50]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        h2_CPA_pt[13]->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[51]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                }
            }
            for(auto it : vAntiLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        h2_CPA_pt[15]->Fill(it.GetPt(), it.GetCPA());
                    h2_pt_invMass[52]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        h2_CPA_pt[16]->Fill(it.GetPt(), it.GetCPA());
                    h2_pt_invMass[53]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        h2_CPA_pt[17]->Fill(it.GetPt(), it.GetCPA());
                    h2_pt_invMass[54]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        h2_CPA_pt[18]->Fill(it.GetPt(), it.GetCPA());
                    h2_pt_invMass[55]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                }
            }
        }


        for(size_t i = 0; i < vLambda.size(); i++)
        {
            if (vLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[5]->Fill(vLambda[i].GetPt(), vLambda[i].GetCPA());
                h2_pt_invMass[18]->Fill(vLambda[i].GetPt(), CalculateInvMassLambda(vLambda[i], false));
                    // h2_pt_invMass[26]->Fill(CalculateInvMassLambda(vLambda[i], false), RelativePairMomentum(vLambda[i], 3122, vLambda[j], 3122));
            }
            vLambda[i].SetUse(true);
        }
        for(size_t i = 0; i < vAntiLambda.size(); i++)
        {
            if (vAntiLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[6]->Fill(vAntiLambda[i].GetPt(), vAntiLambda[i].GetCPA());
                h2_pt_invMass[22]->Fill(vAntiLambda[i].GetPt(), CalculateInvMassLambda(vAntiLambda[i], true));
                    // h2_pt_invMass[30]->Fill(CalculateInvMassLambda(vAntiLambda[i], true), CalculateInvMassLambda(vAntiLambda[i], true));
            }
            vAntiLambda[i].SetUse(true);
        }

        
        // Part Collection #3   -   Clean Decay At Random
        fPairCleaner->ResetArray();
        fPairCleaner->CleanDecayAtRandom(&vLambda, 4);
        fPairCleaner->CleanDecayAtRandom(&vAntiLambda, 5);
        fPairCleaner->StoreParticle(vLambda);
        fPairCleaner->StoreParticle(vAntiLambda);
        fPartColl3->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); 


        if(fLambdaV0Cuts->GetIsMonteCarlo())
        {
            for(auto it : vLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        CPAPtBinningPrim_lambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[56]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        CPAPtBinningMat_lambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[57]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        CPAPtBinningSec_lambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[58]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        CPAPtBinningCont_lambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[59]->Fill(it.GetPt(), CalculateInvMassLambda(it, false));
                    }
                }
            }
            for(auto it : vAntiLambda)
            {
                if(it.UseParticle())
                {
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { 
                        CPAPtBinningPrim_antilambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[60]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { 
                        CPAPtBinningMat_antilambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[61]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { 
                        CPAPtBinningSec_antilambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[62]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                    if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { 
                        CPAPtBinningCont_antilambda->Fill(it.GetPt(), it.GetCPA()); 
                        h2_pt_invMass[63]->Fill(it.GetPt(), CalculateInvMassLambda(it, true));
                    }
                }
            }
        }

        for(size_t i = 0; i < vLambda.size(); i++)
        {
            if (vLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[7]->Fill(vLambda[i].GetPt(), vLambda[i].GetCPA());
                h2_pt_invMass[19]->Fill(vLambda[i].GetPt(), CalculateInvMassLambda(vLambda[i], false));
                    // h2_pt_invMass[27]->Fill(CalculateInvMassLambda(vLambda[i], false), RelativePairMomentum(vLambda[i], 3122, vLambda[j], 3122));
            }
            vLambda[i].SetUse(true);
        }
        for(size_t i = 0; i < vAntiLambda.size(); i++)
        {
            if (vAntiLambda[i].UseParticle() == true)
            {
                h2_CPA_pt[8]->Fill(vAntiLambda[i].GetPt(), vAntiLambda[i].GetCPA());
                h2_pt_invMass[23]->Fill(vAntiLambda[i].GetPt(), CalculateInvMassLambda(vAntiLambda[i], true));
                    // h2_pt_invMass[31]->Fill(vAntiLambda[i].GetPt(), CalculateInvMassLambda(vAntiLambda[i], true));
            }
            vAntiLambda[i].SetUse(true);
        }
        
        

        // std::cout << "dimensions of cleanParticles: " << fPairCleaner->GetCleanParticles().size() << "x" << fPairCleaner->GetCleanParticles()[0].size() << std::endl;        ### makes particleTypeXparticleNumber
        // if(fPairCleaner->GetCleanParticles()[0].size() || fPairCleaner->GetCleanParticles()[1].size() || fPairCleaner->GetCleanParticles()[2].size() || fPairCleaner->GetCleanParticles()[3].size())
        // {
        //     std::cout << "### New EVENT ###" << std::endl;
        //     if(!fPairCleaner->GetCleanParticles()[0].size()==0) std::cout << "Lambdas: " << fPairCleaner->GetCleanParticles()[0].size() << std::endl;
        //     if(!fPairCleaner->GetCleanParticles()[1].size()==0)std::cout << "AntiLambdas: " << fPairCleaner->GetCleanParticles()[1].size() << std::endl;
        //     if(!fPairCleaner->GetCleanParticles()[2].size()==0)std::cout << "Xi: " << fPairCleaner->GetCleanParticles()[2].size() << std::endl;
        //     if(!fPairCleaner->GetCleanParticles()[3].size()==0)std::cout << "AntiXi: " << fPairCleaner->GetCleanParticles()[3].size() << std::endl;
        // }
        // if(vLambda.size() > 0 && vXi.size() > 0)
        // {
        //     std::cout << "after paircleaner vLambda und vXi größe: " << vLambda.size() << " und " << vXi.size() << std::endl;
        // }

        
        //  ######################################################################
        //  ##
        //  ##
        //  ##
        //  ##                     DAUGHTER COMBINATION MIXING
        //  ##                          after Paircleaning
        //  ##
        //  ##
        //  ##
        //  ##
        //  #######################################################################

//         if(fLambdaV0Cuts->GetIsMonteCarlo())
//         {
// // Origin binning for FemtoDreamPaircleaning
//             for(auto it : vLambda)
//             {
//                 if(it.UseParticle())
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { CPAPtBinningPrim_lambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { CPAPtBinningMat_lambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { CPAPtBinningSec_lambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { CPAPtBinningCont_lambda->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//                 else
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { CPAPtBinningPrim_lambda_dump->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { CPAPtBinningMat_lambda_dump->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { CPAPtBinningSec_lambda_dump->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { CPAPtBinningCont_lambda_dump->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//                 it.SetUse(true);
//             }
//             for(auto it : vAntiLambda)
//             {
//                 if(it.UseParticle())
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { CPAPtBinningPrim_antilambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { CPAPtBinningMat_antilambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { CPAPtBinningSec_antilambda->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { CPAPtBinningCont_antilambda->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//                 it.SetUse(true);
//             }


// // origin Binning after InvMass Paircleaning /////////////////////////////////////
//             CleanDecay(&vLambda, "Lambda");
//             CleanDecay(&vAntiLambda, "AntiLambda");
//             for(auto it : vLambda)
//             {
//                 if(it.UseParticle())
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { h2_pt_invMass[0]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { h2_pt_invMass[1]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { h2_pt_invMass[2]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { h2_pt_invMass[3]->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//                 else
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { h2_pt_invMass[4]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { h2_pt_invMass[5]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { h2_pt_invMass[6]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { h2_pt_invMass[7]->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//             }
//             for(auto it : vAntiLambda)
//             {
//                 if(it.UseParticle())
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { h2_pt_invMass[8]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { h2_pt_invMass[9]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { h2_pt_invMass[10]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { h2_pt_invMass[11]->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//                 else
//                 {
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kPhysPrimary)   { h2_pt_invMass[12]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kMaterial)      { h2_pt_invMass[13]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kWeak)          { h2_pt_invMass[14]->Fill(it.GetPt(), it.GetCPA()); }
//                     if(it.GetParticleOrigin() == AliFemtoDreamBasePart::kFake)          { h2_pt_invMass[15]->Fill(it.GetPt(), it.GetCPA()); }
//                 }
//             }
//         }
        //###########################################
        //###########################################
        //
        //
        //           Postdata
        //
        //
        //###########################################
        //###########################################
        PostData(1, tlEventCuts);           // cuts keeping Lambda
        PostData(2, tlLambdaList);
        PostData(3, tlAntiLambdaList);
        PostData(4, tlResults);             // only Decay Cleaning
        PostData(5, tlResultsQA);
        PostData(6, tlRecombination_after);
        PostData(7, tlResults2);     
        PostData(8, tlResults3);     
    
        if (fLambdaV0Cuts->GetIsMonteCarlo())
        {
            PostData(9, tlLambdaMC);
        }
        if (fAntiLambdaV0Cuts->GetIsMonteCarlo())
        {
            PostData(10, tlAntiLambdaMC);
        }
    } // ende Event
}



//  #######################################################################
//  #######################################################################
//  #######################################################################
//  ##
//  ##                     CUSTOM FUNCTIONS
//  ##                     
//  ##
//  #######################################################################
//  #######################################################################
//  #######################################################################

void AliAnalysisTaskPOmegaPenne::ResetGlobalTrackReference()
{
    //This method was inherited form H. Beck analysis
    for (UShort_t i = 0; i < fTrackBufferSize; i++)
    {
        fGTI[i] = nullptr;
        // std::fill(fGTI.begin(),fGTI.end(), nullptr);
    }
}

//  Stores TrackID in Global Track Reference Array 'fGTI' if ID > 0
//
void AliAnalysisTaskPOmegaPenne::StoreGlobalTrackReference(AliVTrack *vTrack)
{
    //This method was inherited form H. Beck analysis
    AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(vTrack);
    const int trackID = vTrack->GetID();
    if (trackID < 0)
    {
        return;
    }
    if (trackID >= fTrackBufferSize)
    {
        printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID, fTrackBufferSize);
        return;
    }

    if (fGTI[trackID])
    {
        if ((!nanoTrack->GetFilterMap()) && (!vTrack->GetTPCNcls()))
        {
            return;
        }
        if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
        {
            printf("WARNING! global track info already there!");
            printf("    ###     TPCNcls track1 %u Track2 %u", (fGTI[trackID])->GetTPCNcls(), vTrack->GetTPCNcls());
            printf("   ###     FilterMap Track1 %u track2 %u\n", dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(), nanoTrack->GetFilterMap());
        }
    }
    fGTI[trackID] = vTrack;

}

// Always Negative Daughter First - Second Argument is the Positive Daughter
// -> in Baseparts GetMomentum(1), GetMomentum(2)
float AliAnalysisTaskPOmegaPenne::CalculateInvMassLambda(TVector3 momNegDaughter, int PDGnegDaughter, TVector3 momPosDaughter, int PDGposDaughter)
{
    float invMass = 0;
    
    float massDP = TDatabasePDG::Instance()->GetParticle(PDGposDaughter)->Mass(); // Proton
    float massDN = TDatabasePDG::Instance()->GetParticle(PDGnegDaughter)->Mass();  // Pion
    float EDaugP = TMath::Sqrt(
        massDP * massDP + 
        momPosDaughter.X() * momPosDaughter.X() + 
        momPosDaughter.Y() * momPosDaughter.Y() + 
        momPosDaughter.Z() * momPosDaughter.Z()
    );
    float EDaugN = TMath::Sqrt(
        massDN * massDN + 
        momNegDaughter.X() * momNegDaughter.X() + 
        momNegDaughter.Y() * momNegDaughter.Y() + 
        momNegDaughter.Z() * momNegDaughter.Z()
    );
    float energysum = EDaugP + EDaugN;
    float pSum2 = 
        ( momNegDaughter.X() + momPosDaughter.X() ) * 
        ( momNegDaughter.X() + momPosDaughter.X() ) 
        +
        ( momNegDaughter.Y() + momPosDaughter.Y() ) * 
        ( momNegDaughter.Y() + momPosDaughter.Y() ) 
        + 
        ( momNegDaughter.Z() + momPosDaughter.Z() ) * 
        ( momNegDaughter.Z() + momPosDaughter.Z() )
    ;
    invMass = TMath::Sqrt(energysum * energysum - pSum2);
    return invMass;
}
float AliAnalysisTaskPOmegaPenne::CalculateInvMassLambda(AliFemtoDreamBasePart *lambdaParticle, bool isAntiParticle)
{
    if(!isAntiParticle)
    {
        return CalculateInvMassLambda(lambdaParticle->GetMomentum(1), 211,
                                      lambdaParticle->GetMomentum(2), 2212);
    }
    else
    {
        return CalculateInvMassLambda(lambdaParticle->GetMomentum(1), 2212,
                                      lambdaParticle->GetMomentum(2), 211);
    }
}
float AliAnalysisTaskPOmegaPenne::CalculateInvMassLambda(AliFemtoDreamBasePart lambdaParticle, bool isAntiParticle)
{
    if(!isAntiParticle)
    {
        return CalculateInvMassLambda(lambdaParticle.GetMomentum(1), 211,
                                      lambdaParticle.GetMomentum(2), 2212);
    }
    else
    {
        return CalculateInvMassLambda(lambdaParticle.GetMomentum(1), 2212,
                                      lambdaParticle.GetMomentum(2), 211);
    }
}
// Parameter = Bachelor , Positive Daughter , Negative Dautgher
// -> in BaseParts GetMomentum(3) , GetMomentum(1), GetMomentum(2)
float AliAnalysisTaskPOmegaPenne::CalculateInvMassXi(TVector3 momBach, int PGGbach, TVector3 momPosDaughter, int PDGposDaughter, TVector3 momNegDaughter, int PDGnegDaughter)
{
    // float massPosDaugh = TDatabasePDG::Instance()->GetParticle(PDGposDaughter)->Mass();  // Proton 2212 or antiPion 211
    // float massNegDaugh = TDatabasePDG::Instance()->GetParticle(PDGnegDaughter)->Mass();   // Pion 211 or antiProton 2212
    float massBach = TDatabasePDG::Instance()->GetParticle(PGGbach)->Mass();            // Pion 211 or antiPion 211
    // float massV0 = CalculateInvMassLambda(momNegDaughter, PDGnegDaughter, momPosDaughter, PDGposDaughter);                 // Lambda
    float massV0 = TDatabasePDG::Instance()->GetParticle(3122)->Mass();     // lambda 3122
    
    TVector3 PtotV0 = (momPosDaughter + momNegDaughter);
    float Ev0 = TMath::Sqrt(massV0 * massV0 + PtotV0.Mag2());

    float EBach = TMath::Sqrt(TMath::Power(massBach, 2) + momBach.Mag2());

    float Ptot2Casc = (PtotV0 + momBach).Mag2();

    return TMath::Sqrt(TMath::Power(Ev0 + EBach,2) - Ptot2Casc);
}

float AliAnalysisTaskPOmegaPenne::CalculateInvMassXi(AliFemtoDreamBasePart *xiParticle, bool isAntiParticle)
{
    if(!isAntiParticle)
    {
        return CalculateInvMassXi(xiParticle->GetMomentum(3), 211, 
                                  xiParticle->GetMomentum(2), 2212,
                                  xiParticle->GetMomentum(1), 211);
    }
    else
    {
        return CalculateInvMassXi(xiParticle->GetMomentum(3), 211, 
                                  xiParticle->GetMomentum(2), 211,
                                  xiParticle->GetMomentum(1), 2212);
    }
}

void AliAnalysisTaskPOmegaPenne::CleanDecay(std::vector<AliFemtoDreamBasePart> *Decay, string particleSteering)
{
    float fPDGMassPart = 1.0;
    float fWeightPart1 = 1.0;
    float fWeightPart2 = 1.0;
    float fMassPart1 = 0.0;
    float fMassPart2 = 0.0;
    float fMassToPDG1 = 0.0;
    float fMassToPDG2 = 0.0;
    std::vector<int> IDDaug1;
    std::vector<int> IDDaug2;

    if(particleSteering == "Lambda" || particleSteering == "AntiLambda")
    {
        fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    }
    else if(particleSteering == "Xi" || particleSteering == "AntiXi")
    {
        fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    }
    else
    {
        std::cout << std::endl;
        std::cout << "######################################################" << std::endl;       
        std::cout << "Teilchensorte nicht korrekt ausgewählt: " << particleSteering << std::endl;
        std::cout << "kenne nur (Anti-)Lambda und (Anti-)Xi" << std::endl;       
        std::cout << "######################################################" << std::endl;       
        std::cout << std::endl;
        return;
    }
    
    for (std::vector<AliFemtoDreamBasePart>::iterator itDecay1 = Decay->begin();
         itDecay1 != Decay->end(); ++itDecay1)
    {
        if (itDecay1->UseParticle() == true)
        {
            for (auto itDecay2 = itDecay1 + 1; itDecay2 != Decay->end(); ++itDecay2)
            {
                if (itDecay1->UseParticle() == false) // break if particle 1 has lost the selection and been set to false
                {   //statistics on how much this has happened and how many other particles would have cleaned by an already cleaned particle 1 could be interesting
                    break;
                }
                if (itDecay2->UseParticle() == false)
                {
                    continue;
                }
                IDDaug1 = itDecay1->GetIDTracks();
                IDDaug2 = itDecay2->GetIDTracks();
                for (auto itID1s = IDDaug1.begin(); itID1s != IDDaug1.end(); ++itID1s)
                {
                    for (auto itID2s = IDDaug2.begin(); itID2s != IDDaug2.end(); ++itID2s)
                    {
                        if (*itID1s == *itID2s)
                        {
                            if (particleSteering == "Lambda")
                            {
                                fMassPart1 = CalculateInvMassLambda(itDecay1->GetMomentum(1), 211, itDecay1->GetMomentum(2), 2212);
                                fMassPart2 = CalculateInvMassLambda(itDecay2->GetMomentum(1), 211, itDecay2->GetMomentum(2), 2212);
                                fWeightPart1 = WeightLambda(itDecay1->GetPt());
                                fWeightPart2 = WeightLambda(itDecay2->GetPt());
                                // PDG - 3122 - Lambda
                                // fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0) / fWeightPart1;
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0) / fWeightPart2;
                                if (::abs(fMassToPDG1) >= ::abs(fMassToPDG2))
                                {
                                    itDecay1->SetUse(false);
                                    hLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    // h2_CPA_pt[2]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    // h2_CPA_pt[2]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            else if (particleSteering == "AntiLambda")
                            {
                                fMassPart1 = CalculateInvMassLambda(itDecay1->GetMomentum(1), 2212, itDecay1->GetMomentum(2), 211);
                                fMassPart2 = CalculateInvMassLambda(itDecay2->GetMomentum(1), 2212, itDecay2->GetMomentum(2), 211);
                                fWeightPart1 = WeightAntiLambda(itDecay1->GetPt());
                                fWeightPart2 = WeightAntiLambda(itDecay2->GetPt());
                                // PDG - 3122 - Lambda
                                // fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0) / fWeightPart1;
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0) / fWeightPart2;
                                if (::abs(fMassToPDG1) >= ::abs(fMassToPDG2))
                                {
                                    itDecay1->SetUse(false);
                                    hAntiLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    // h2_CPA_pt[7]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hAntiLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    // h2_CPA_pt[7]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            else if (particleSteering == "Xi")
                            {
                                fMassPart1 = CalculateInvMassXi(itDecay1->GetMomentum(3), 211, itDecay1->GetMomentum(2), 2212, itDecay1->GetMomentum(1), 211);
                                fMassPart2 = CalculateInvMassXi(itDecay2->GetMomentum(3), 211, itDecay2->GetMomentum(2), 2212, itDecay2->GetMomentum(1), 211);
                                fWeightPart1 = WeightXi(itDecay1->GetPt());
                                fWeightPart2 = WeightXi(itDecay2->GetPt());
                                // PDG - 3312 - Xi

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0) / fWeightPart1;
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0) / fWeightPart2;
                                if (::abs(fMassToPDG1) >= ::abs(fMassToPDG2))
                                {
                                    itDecay1->SetUse(false);
                                    hXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    // h2_CPA_pt[12]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    // h2_CPA_pt[12]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            else if (particleSteering == "AntiXi")
                            {
                                fMassPart1 = CalculateInvMassXi(itDecay1->GetMomentum(3), 211, itDecay1->GetMomentum(2), 211, itDecay1->GetMomentum(1), 2212);
                                fMassPart2 = CalculateInvMassXi(itDecay2->GetMomentum(3), 211, itDecay2->GetMomentum(2), 211, itDecay2->GetMomentum(1), 2212);
                                fWeightPart1 = WeightAntiXi(itDecay1->GetPt());
                                fWeightPart2 = WeightAntiXi(itDecay2->GetPt());
                                // PDG - 3312 - Xi
                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0) / fWeightPart1;
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0) / fWeightPart2;
                                if (::abs(fMassToPDG1) >= ::abs(fMassToPDG2))
                                {
                                    itDecay1->SetUse(false);
                                    hAntiXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    // h2_CPA_pt[17]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hAntiXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    // h2_CPA_pt[17]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            // std::cout << "######################################################" << std::endl;
                            // std::cout << "*************** CleanDecay ***************" << std::endl;
                            // if(particleSteering == "Lambda") std::cout << "Lambda" << std::endl;
                            // if(particleSteering == "AntiLambda") std::cout << "AntiLambda" << std::endl;
                            // if(particleSteering == "Xi") std::cout << "Xi" << std::endl;
                            // if(particleSteering == "AntiXi") std::cout << "AntiXi" << std::endl;
                            // std::cout << "fWeightPart1: " << fWeightPart1 << std::endl;
                            // std::cout << "fWeightPart2: " << fWeightPart2 << std::endl;
                            // std::cout << "itDecay1->Pt: " << itDecay1->GetPt() << std::endl;
                            // std::cout << "itDecay2->Pt: " << itDecay2->GetPt() << std::endl;
                            // std::cout << "fMassPart1: " << fMassPart1 << std::endl;
                            // std::cout << "fMassPart2: " << fMassPart2 << std::endl;
                            // std::cout << "fMassToPDG1: " << fMassToPDG1 << std::endl;
                            // std::cout << "fMassToPDG2: " << fMassToPDG2 << std::endl;
                            // std::cout << "######################################################" << std::endl;
                        }
                    }
                }
            }
        }
        else
            continue;
    }
}

void AliAnalysisTaskPOmegaPenne::CleanDecayAtRandom(std::vector<AliFemtoDreamBasePart> *Decay, string particleSteering)
{
    float fPDGMassPart = 1.0;
    // float fWeightPart1 = 1.0;
    // float fWeightPart2 = 1.0;
    float fMassPart1 = 0.0;
    float fMassPart2 = 0.0;
    float fMassToPDG1 = 0.0;
    float fMassToPDG2 = 0.0;
    std::vector<int> IDDaug1;
    std::vector<int> IDDaug2;

    if(particleSteering == "Lambda" || particleSteering == "AntiLambda")
    {
        fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    }
    else if(particleSteering == "Xi" || particleSteering == "AntiXi")
    {
        fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    }
    else
    {
        std::cout << std::endl;
        std::cout << "######################################################" << std::endl;       
        std::cout << "Teilchensorte nicht korrekt ausgewählt: " << particleSteering << std::endl;
        std::cout << "kenne nur (Anti-)Lambda und (Anti-)Xi" << std::endl;       
        std::cout << "######################################################" << std::endl;       
        std::cout << std::endl;
        return;
    }
    
    for (std::vector<AliFemtoDreamBasePart>::iterator itDecay1 = Decay->begin();
         itDecay1 != Decay->end(); ++itDecay1)
    {
        if (itDecay1->UseParticle() == true)
        {
            for (auto itDecay2 = itDecay1 + 1; itDecay2 != Decay->end(); ++itDecay2)
            {
                if (itDecay1->UseParticle() == false) // break if particle 1 has lost the selection and been set to false
                {   //statistics on how much this has happened and how many other particles would have cleaned by an already cleaned particle 1 could be interesting
                    break;
                }
                if (itDecay2->UseParticle() == false)
                {
                    continue;
                }
                IDDaug1 = itDecay1->GetIDTracks();
                IDDaug2 = itDecay2->GetIDTracks();
                for (auto itID1s = IDDaug1.begin(); itID1s != IDDaug1.end(); ++itID1s)
                {
                    for (auto itID2s = IDDaug2.begin(); itID2s != IDDaug2.end(); ++itID2s)
                    {
                        if (*itID1s == *itID2s)
                        {
                            if (particleSteering == "Lambda")
                            {
                                fMassPart1 = CalculateInvMassLambda(itDecay1->GetMomentum(1), 211, itDecay1->GetMomentum(2), 2212);
                                fMassPart2 = CalculateInvMassLambda(itDecay2->GetMomentum(1), 211, itDecay2->GetMomentum(2), 2212);
                                // PDG - 3122 - Lambda
                                // fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0);
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0);
                                if (gRandom->Uniform(0., 1.) < 0.5)
                                {
                                    itDecay1->SetUse(false);
                                    hLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    h2_CPA_pt[10]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    h2_CPA_pt[10]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                            }
                            else if (particleSteering == "AntiLambda")
                            {
                                fMassPart1 = CalculateInvMassLambda(itDecay1->GetMomentum(1), 2212, itDecay1->GetMomentum(2), 211);
                                fMassPart2 = CalculateInvMassLambda(itDecay2->GetMomentum(1), 2212, itDecay2->GetMomentum(2), 211);
                                // PDG - 3122 - Lambda
                                // fPDGMassPart = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0);
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0);
                                if (gRandom->Uniform(0., 1.) < 0.5)
                                {
                                    itDecay1->SetUse(false);
                                    hAntiLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    h2_CPA_pt[11]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hAntiLambdaCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    h2_CPA_pt[11]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            else if (particleSteering == "Xi")
                            {
                                fMassPart1 = CalculateInvMassXi(itDecay1->GetMomentum(3), 211, itDecay1->GetMomentum(2), 2212, itDecay1->GetMomentum(1), 211);
                                fMassPart2 = CalculateInvMassXi(itDecay2->GetMomentum(3), 211, itDecay2->GetMomentum(2), 2212, itDecay2->GetMomentum(1), 211);
                                // PDG - 3312 - Xi

                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0);
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0);
                                if (gRandom->Uniform(0., 1.) < 0.5)
                                {
                                    itDecay1->SetUse(false);
                                    hXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    h2_CPA_pt[12]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    h2_CPA_pt[12]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                            else if (particleSteering == "AntiXi")
                            {
                                fMassPart1 = CalculateInvMassXi(itDecay1->GetMomentum(3), 211, itDecay1->GetMomentum(2), 211, itDecay1->GetMomentum(1), 2212);
                                fMassPart2 = CalculateInvMassXi(itDecay2->GetMomentum(3), 211, itDecay2->GetMomentum(2), 211, itDecay2->GetMomentum(1), 2212);
                                // PDG - 3312 - Xi
                                fMassToPDG1 = ((fMassPart1 - fPDGMassPart) * 1000.0);
                                fMassToPDG2 = ((fMassPart2 - fPDGMassPart) * 1000.0);
                                if (gRandom->Uniform(0., 1.) < 0.5)
                                {
                                    itDecay1->SetUse(false);
                                    hAntiXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG1);
                                    h2_CPA_pt[17]->Fill(itDecay1->GetPt(), itDecay1->GetCPA());
                                }
                                else
                                {
                                    itDecay2->SetUse(false);
                                    hAntiXiCleanedPartMassDiffToPDG_Decay->Fill(fMassToPDG2);
                                    h2_CPA_pt[17]->Fill(itDecay2->GetPt(), itDecay2->GetCPA());
                                }
                            }
                        }
                    }
                }
            }
        }
        else
            continue;
    }
}

void AliAnalysisTaskPOmegaPenne::CleanDecayAndDecay(std::vector<AliFemtoDreamBasePart> *vecLambda,
                                                    std::vector<AliFemtoDreamBasePart> *vecXi,
                                                    bool isAntiParticle)
{
    // int counter = 0;
    float fPDGMassLambda = 1.0;
    float fPDGMassXi = 1.0;
    float fWeightLambda = 1.0;
    float fWeightXi = 1.0;
    float fMassLambda = 0.0;
    float fMassXi = 0.0;
    float fMassToPDGLambda = 0.0;
    float fMassToPDGXi = 0.0;

    fPDGMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    fPDGMassXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();

    std::vector<int> IDDaug1;
    std::vector<int> IDDaug2;

    for (auto itLambdaPart = vecLambda->begin(); itLambdaPart != vecLambda->end(); ++itLambdaPart)
    {
        if (itLambdaPart->UseParticle() == true)
        {
            for (auto itXiPart = vecXi->begin(); itXiPart != vecXi->end(); ++itXiPart)
            {
                if (itXiPart->UseParticle() == false)
                {
                    continue;
                }
                IDDaug1 = itLambdaPart->GetIDTracks();
                IDDaug2 = itXiPart->GetIDTracks();
                for (auto itID1s = IDDaug1.begin(); itID1s != IDDaug1.end(); ++itID1s)
                {
                    for (auto itID2s = IDDaug2.begin(); itID2s != IDDaug2.end(); ++itID2s)
                    {
                        if (*itID1s == *itID2s)
                        {
                            if (isAntiParticle == false)
                            {
                                fMassLambda = CalculateInvMassLambda(itLambdaPart->GetMomentum(1), 211, itLambdaPart->GetMomentum(2), 2212);
                                fMassXi = CalculateInvMassXi(itXiPart->GetMomentum(3), 211, itXiPart->GetMomentum(2), 2212, itXiPart->GetMomentum(1), 211);

                                fWeightLambda = WeightLambda(itLambdaPart->GetPt());
                                fWeightXi = WeightXi(itXiPart->GetPt());

                                fMassToPDGLambda = ((fMassLambda - fPDGMassLambda) * 1000.0) / fWeightLambda;
                                fMassToPDGXi = ((fMassXi - fPDGMassXi) * 1000.0) / fWeightXi;
                                
                                if (TMath::Abs(fMassToPDGLambda) < TMath::Abs(fMassToPDGXi))
                                {
                                    itXiPart->SetUse(false);
                                    hXiCleanedPartMassDiffToPDG_DecayDecay->Fill(fMassToPDGXi);
                                    hXiCleanedPartMass_DecayDecay->Fill(fMassXi);
                                    h2_CPA_pt[14]->Fill(itXiPart->GetPt(), itXiPart->GetCPA());
                                }
                                else
                                {
                                    itLambdaPart->SetUse(false);
                                    hLambdaCleanedPartMassDiffToPDG_DecayDecay->Fill(fMassToPDGLambda);
                                    hLambdaCleanedPartMass_DecayDecay->Fill(fMassLambda);
                                    h2_CPA_pt[4]->Fill(itLambdaPart->GetPt(), itLambdaPart->GetCPA());
                                }
                            }
                            else if (isAntiParticle == true)
                            {
                                fMassLambda = CalculateInvMassLambda(itLambdaPart->GetMomentum(1), 2212, itLambdaPart->GetMomentum(2), 211);
                                fMassXi = CalculateInvMassXi(itXiPart->GetMomentum(3), 211, itXiPart->GetMomentum(2), 211, itXiPart->GetMomentum(1), 2212);

                                fWeightLambda = WeightAntiLambda(itLambdaPart->GetPt());
                                fWeightXi = WeightAntiXi(itXiPart->GetPt());

                                fMassToPDGLambda = ((fMassLambda - fPDGMassLambda) * 1000.0) / fWeightLambda;
                                fMassToPDGXi = ((fMassXi - fPDGMassXi) * 1000.0) / fWeightXi;

                                if (TMath::Abs(fMassToPDGLambda) < TMath::Abs(fMassToPDGXi))
                                {
                                    itXiPart->SetUse(false);
                                    hAntiXiCleanedPartMassDiffToPDG_DecayDecay->Fill(fMassToPDGXi);
                                    hAntiXiCleanedPartMass_DecayDecay->Fill(fMassXi);
                                    h2_CPA_pt[19]->Fill(itXiPart->GetPt(), itXiPart->GetCPA());
                                }
                                else
                                {
                                    itLambdaPart->SetUse(false);
                                    hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay->Fill(fMassToPDGLambda);
                                    hAntiLambdaCleanedPartMass_DecayDecay->Fill(fMassLambda);
                                    h2_CPA_pt[9]->Fill(itLambdaPart->GetPt(), itLambdaPart->GetCPA());
                                }
                            }
                            // std::cout << "######################################################" << std::endl;
                            // std::cout << "*************** CleanDecayAndDecay ***************" << std::endl;
                            // if(isAntiParticle == true ) std::cout << "*** ANTI Teilchen ***" << std::endl;
                            // if(isAntiParticle == false) std::cout << "*** Teilchen ***" << std::endl;
                            // std::cout << "fWeightLambda: " << fWeightLambda << std::endl;
                            // std::cout << "fWeightXi: " << fWeightXi << std::endl;
                            // std::cout << "itLambdaPart->Pt: " << itLambdaPart->GetPt() << std::endl;
                            // std::cout << "itXiPart->Pt: " << itXiPart->GetPt() << std::endl;
                            // std::cout << "fMassLambda: " << fMassLambda << std::endl;
                            // std::cout << "fMassXi: " << fMassXi << std::endl;
                            // std::cout << "fMassToPDGLambda: " << fMassToPDGLambda << std::endl;
                            // std::cout << "fMassToPDGXi: " << fMassToPDGXi << std::endl;
                            // std::cout << "######################################################" << std::endl;
                        }
                    }
                }
                if (itLambdaPart->UseParticle() == false)
                {
                    break;
                }
            }
        }
        else
            continue;
    }
}

//                                            //
// weights from nanoAOD run 503_20200611-1233 //
//          12.06.2020                        //
//                                            //
float AliAnalysisTaskPOmegaPenne::WeightLambda(float pT)
{
    if (pT > 4.05)
    {
        return 1.8;
    }
    else if (pT < 0.55)
    {
        return 1.7;
    }
    else
    {    
        return ( 
             - 0.05079408f  * ::pow(pT,7)
             + 0.7466736f   * ::pow(pT,6)
             - 4.390268f    * ::pow(pT,5) 
             + 13.22549f    * ::pow(pT,4) 
             - 21.79303f    * ::pow(pT,3) 
             + 19.59184f    * ::pow(pT,2) 
             - 9.1271373f   *       pT 
             + 3.4100528f
        );
    }
}
float AliAnalysisTaskPOmegaPenne::WeightAntiLambda(float pT)
{
    if (pT > 4.05)
    {
        return 1.8;
    }
    else if (pT < 0.55)
    {
        return 1.7;
    }
    else
    {    
        return (  
                  0.002539743f  * ::pow(pT,7)
                - 0.07200098f   * ::pow(pT,6)
                + 0.7325843f    * ::pow(pT,5) 
                - 3.624922f     * ::pow(pT,4) 
                + 9.481334f     * ::pow(pT,3) 
                - 12.90322f     * ::pow(pT,2) 
                + 8.182929f     *       pT 
                - 0.1780557f
        );
    }
}
float AliAnalysisTaskPOmegaPenne::WeightXi(float pT)
{
    if (pT > 6.3)
    {
        return 2.8;
    }
    else if (pT < 0.55)
    {
        return 2.3;
    }
    else
    {    
        return (
            - 0.00005897109f * ::pow(pT,7)
            + 0.002190119f   * ::pow(pT,6)
            - 0.02780911f    * ::pow(pT,5) 
            + 0.1742865f     * ::pow(pT,4) 
            - 0.6289920f     * ::pow(pT,3) 
            + 1.407963f      * ::pow(pT,2) 
            - 1.710865f      *       pT 
            + 2.903313f
        );
    }
}
float AliAnalysisTaskPOmegaPenne::WeightAntiXi(float pT)
{
    if (pT > 6.3)
    {
        return 2.8;
    }
    else if (pT < 0.55)
    {
        return 2.3;
    }
    else
    {    
        return (  
              0.001766952f  * ::pow(pT,9)
            - 0.05378248f   * ::pow(pT,8)
            + 0.6983620f    * ::pow(pT,7)
            - 5.045875f     * ::pow(pT,6)
            + 22.18263f     * ::pow(pT,5) 
            - 60.92189f     * ::pow(pT,4) 
            + 103.1763f     * ::pow(pT,3) 
            - 102.1255f     * ::pow(pT,2) 
            + 52.47001f     *       pT 
            - 8.243706f
        );
    }
}

float AliAnalysisTaskPOmegaPenne::RelativePairMomentum(AliFemtoDreamBasePart *part1, const int pdg1, AliFemtoDreamBasePart *part2, const int pdg2) 
{
  TLorentzVector PartOne, PartTwo;

  PartOne.SetXYZM(part1->GetMomentum().X(), part1->GetMomentum().Y(), part1->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2->GetMomentum().X(), part2->GetMomentum().Y(), part2->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());

  TLorentzVector trackSum = PartOne + PartTwo;
  
  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}

float AliAnalysisTaskPOmegaPenne::RelativePairMomentum(AliFemtoDreamBasePart part1, const int pdg1, AliFemtoDreamBasePart part2, const int pdg2) 
{
  TLorentzVector PartOne, PartTwo;

  PartOne.SetXYZM(part1.GetMomentum().X(), part1.GetMomentum().Y(), part1.GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2.GetMomentum().X(), part2.GetMomentum().Y(), part2.GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());

  TLorentzVector trackSum = PartOne + PartTwo;
  
  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}