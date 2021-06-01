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

/* AliAnaysisTaskHypertritonKFTree
 *
 */

#include <iostream>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliTimeRangeCut.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "TChain.h"

#include "AliAnalysisTaskHypertritonKFTree.h"

/// your analysis class
class AliAnalysisTaskHypertritonKFTree;
/// classimp: necessary for root
ClassImp(AliAnalysisTaskHypertritonKFTree)
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTree::AliAnalysisTaskHypertritonKFTree() : AliAnalysisTaskSE(), 
fEventCuts(),
fTimeRangeCut(),
fPIDResponse(nullptr),
PrimVertex(),
He3CandTrackIdPos(10),
He3CandTrackIdNeg(10),
PionCandTrackIdPos(1500),
PionCandTrackIdNeg(1500),
DeuteronCandTrackIdPos(10),
DeuteronCandTrackIdNeg(10),
ProtonCandTrackIdPos(300),
ProtonCandTrackIdNeg(300),
HyperTritonCandidates(50),
DeuteronDaughter(50),
ProtonDaughter(50),
kRun2BodyDecay(false),
kRun3BodyDecay(false),
kDoQA(false),
kIsMC(false),
fESDtrackCuts(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
fMCEvent(nullptr),
hNumberOfEvents(nullptr),
histoEventCentrality(nullptr),
CandidateTree(nullptr),
CentralityPercentile(-999),
mass(-999),
ErrorMass(-999),
px(-999),
py(-999),
pz(-999),
Rapidity(-999),
Charge(0),
Chi2PerNDF(-999),
CosPointingAngle(-999),
DistanceToPV(-999),
DeviationFromPV(-999),
DistanceToPVXY(-999),
DeviationFromPVXY(-999),
massTopo(-999),
ErrorMassTopo(-999),
pxTopo(-999),
pyTopo(-999),
pzTopo(-999),
RapidityTopo(-999),
Chi2PerNDFTopo(-999),
CosPointingAngleTopo(-999),
DistanceToPVTopo(-999),
DeviationFromPVTopo(-999),
DistanceToPVXYTopo(-999),
DeviationFromPVXYTopo(-999),
DecayLength(-999),
ErrorDecayLength(-999),
DecayLengthXY(-999),
ErrorDecayLengthXY(-999),
DistanceOfDaughters(-999),
DeviationOfDaughters(-999),
DistanceOfDaughtersXY(-999),
DeviationOfDaughtersXY(-999),
pxPion(-999),
pyPion(-999),
pzPion(-999),
ChargePion(-999),
DCAPion(-999),
//DistanceToPVPion(-999),
DeviationFromPVPion(-999),
DCAPionXY(-999),
//DistanceToPVPionXY(-999),
DeviationFromPVPionXY(-999),
NCrossedRowsTPCPion(-999),
NPIDClusterTPCPion (-999),
TPCMomPion(-999),
TPCnSigmaPion(-999),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
PIDForTrackingPion(-999),
DistanceToSecVertPion(-999),
DeviationToSecVertPion(-999),
pxHe(-999),
pyHe(-999),
pzHe(-999),
ChargeHe(-999),
DCA3He(-999),
//DistanceToPV3He(-999),
DeviationFromPV3He(-999),
DCA3HeXY(-999),
//DistanceToPV3HeXY(-999),
DeviationFromPV3HeXY(-999),
NCrossedRowsTPC3He(-999),
NPIDClusterTPC3He(-999),
TPCMom3He(-999),
TPCnSigma3He(-999),
TPCnSigma3H(-999),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
PIDForTrackingHe3(-999),
DistanceToSecVertHe(-999),
DeviationToSecVertHe(-999),
///three body decay
CandidateTree_3Body(nullptr),
mass_Deuteron_Proton(-999),
mass_Proton_Pion(-999),
CosPointingAngle_Deuteron_Proton(-999),
DistanceOfDaughters_Deuteron_Proton(-999),
DeviationOfDaughters_Deuteron_Proton(-999),
DistanceOfDaughtersXY_Deuteron_Proton(-999),
DeviationOfDaughtersXY_Deuteron_Proton(-999),
DistanceOfDaughters_Deuteron_Pion(-999),
DeviationOfDaughters_Deuteron_Pion(-999),
DistanceOfDaughtersXY_Deuteron_Pion(-999),
DeviationOfDaughtersXY_Deuteron_Pion(-999),
DistanceOfDaughters_Proton_Pion(-999),
DeviationOfDaughters_Proton_Pion(-999),
DistanceOfDaughtersXY_Proton_Pion(-999),
DeviationOfDaughtersXY_Proton_Pion(-999),
OpeningAngle_Pion_Proton(-999),
OpeningAngle_Pion_Deuteron(-999),
OpeningAngle_Proton_Deuteron(-999),
pxDeuteron(-999),
pyDeuteron(-999),
pzDeuteron(-999),
ChargeDeuteron(-999),
DCADeuteron(-999),
//DistanceToPVDeuteron(-999),
DeviationFromPVDeuteron(-999),
DCADeuteronXY(-999),
//DistanceToPVDeuteronXY(-999),
DeviationFromPVDeuteronXY(-999),
NCrossedRowsTPCDeuteron(-999),
NPIDClusterTPCDeuteron (-999),
TPCMomDeuteron(-999),
TPCnSigmaDeuteron(-999),
TOFnSigmaDeuteron(-999),
HasPointOnITSLayer0Deuteron(false),
HasPointOnITSLayer1Deuteron(false),
HasPointOnITSLayer2Deuteron(false),
HasPointOnITSLayer3Deuteron(false),
HasPointOnITSLayer4Deuteron(false),
HasPointOnITSLayer5Deuteron(false),
PIDForTrackingDeuteron(-999),
DistanceToSecVertDeuteron(-999),
DeviationToSecVertDeuteron(-999),
pxProton(-999),
pyProton(-999),
pzProton(-999),
ChargeProton(-999),
DCAProton(-999),
//DistanceToPVProton(-999),
DeviationFromPVProton(-999),
DCAProtonXY(-999),
//DistanceToPVProtonXY(-999),
DeviationFromPVProtonXY(-999),
NCrossedRowsTPCProton(-999),
NPIDClusterTPCProton (-999),
TPCMomProton(-999),
TPCnSigmaProton(-999),
TOFnSigmaProton(-999),
HasPointOnITSLayer0Proton(false),
HasPointOnITSLayer1Proton(false),
HasPointOnITSLayer2Proton(false),
HasPointOnITSLayer3Proton(false),
HasPointOnITSLayer4Proton(false),
HasPointOnITSLayer5Proton(false),
PIDForTrackingProton(-999),
DistanceToSecVertProton(-999),
DeviationToSecVertProton(-999),
///Monte Carlo
pxMC(-999),
pyMC(-999),
pzMC(-999),
pxHeMC(-999),
pyHeMC(-999),
pzHeMC(-999),
pxPionMC(-999),
pyPionMC(-999),
pzPionMC(-999),
pxProtonMC(-999),
pyProtonMC(-999),
pzProtonMC(-999),
pxDeuteronMC(-999),
pyDeuteronMC(-999),
pzDeuteronMC(-999),
pxVariance(-999),
pyVariance(-999),
pzVariance(-999),
pxHeVariance(-999),
pyHeVariance(-999),
pzHeVariance(-999),
pxPionVariance(-999),
pyPionVariance(-999),
pzPionVariance(-999),
pxProtonVariance(-999),
pyProtonVariance(-999),
pzProtonVariance(-999),
pxDeuteronVariance(-999),
pyDeuteronVariance(-999),
pzDeuteronVariance(-999),
pMC(-999),
pTMC(-999),
RapidityMC(-999),
DecayLengthMC(-999),
DecayLengthXYMC(-999),
xSecVertexMC(-999),
ySecVertexMC(-999),
zSecVertexMC(-999),
xSecVertex(-999),
ySecVertex(-999),
zSecVertex(-999),
xSecVertexVariance(-999),
ySecVertexVariance(-999),
zSecVertexVariance(-999),
NumberOfDeuteronProtonCandidates(-1),
GeneratedTreeMC(nullptr),
GeneratedTreeMC_3Body(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistNsigmaTPCvsPDeuteron(nullptr),
fHistNsigmaTPCvsPProton(nullptr),
fHistNsigmaTOFvsPDeuteron(nullptr),
fHistNsigmaTOFvsPProton(nullptr),
fHistpxTrueRecHe3(nullptr),
fHistpyTrueRecHe3(nullptr),
fHistpzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr),
fHistMomDeuteron(nullptr),
fHistMomProton(nullptr),
fHistMomPion_3Body(nullptr),
fHistNDeuteronsPerEvent(nullptr),
fHistNDeuteronProtonCandidatesPerEvent(nullptr)
{
  /// default constructor, don't allocate memory here!
  /// this is used by root for IO purposes, it needs to remain empty
}
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTree::AliAnalysisTaskHypertritonKFTree(const char* name) : AliAnalysisTaskSE(name),
fEventCuts(),
fTimeRangeCut(),
fPIDResponse(nullptr),
PrimVertex(),
He3CandTrackIdPos(10),
He3CandTrackIdNeg(10),
PionCandTrackIdPos(1500),
PionCandTrackIdNeg(1500),
DeuteronCandTrackIdPos(10),
DeuteronCandTrackIdNeg(10),
ProtonCandTrackIdPos(300),
ProtonCandTrackIdNeg(300),
HyperTritonCandidates(50),
DeuteronDaughter(50),
ProtonDaughter(50),
kRun2BodyDecay(false),
kRun3BodyDecay(false),
kDoQA(false),
kIsMC(false),
fESDtrackCuts(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
fMCEvent(nullptr),
hNumberOfEvents(nullptr),
histoEventCentrality(nullptr),
CandidateTree(nullptr),
CentralityPercentile(-999),
mass(-999),
ErrorMass(-999),
px(-999),
py(-999),
pz(-999),
Rapidity(-999),
Charge(0),
Chi2PerNDF(-999),
CosPointingAngle(-999),
DistanceToPV(-999),
DeviationFromPV(-999),
DistanceToPVXY(-999),
DeviationFromPVXY(-999),
massTopo(-999),
ErrorMassTopo(-999),
pxTopo(-999),
pyTopo(-999),
pzTopo(-999),
RapidityTopo(-999),
Chi2PerNDFTopo(-999),
CosPointingAngleTopo(-999),
DistanceToPVTopo(-999),
DeviationFromPVTopo(-999),
DistanceToPVXYTopo(-999),
DeviationFromPVXYTopo(-999),
DecayLength(-999),
ErrorDecayLength(-999),
DecayLengthXY(-999),
ErrorDecayLengthXY(-999),
DistanceOfDaughters(-999),
DeviationOfDaughters(-999),
DistanceOfDaughtersXY(-999),
DeviationOfDaughtersXY(-999),
pxPion(-999),
pyPion(-999),
pzPion(-999),
ChargePion(-999),
DCAPion(-999),
//DistanceToPVPion(-999),
DeviationFromPVPion(-999),
DCAPionXY(-999),
//DistanceToPVPionXY(-999),
DeviationFromPVPionXY(-999),
NCrossedRowsTPCPion(-999),
NPIDClusterTPCPion (-999),
TPCMomPion(-999),
TPCnSigmaPion(-999),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
PIDForTrackingPion(-999),
DistanceToSecVertPion(-999),
DeviationToSecVertPion(-999),
pxHe(-999),
pyHe(-999),
pzHe(-999),
ChargeHe(-999),
DCA3He(-999),
//DistanceToPV3He(-999),
DeviationFromPV3He(-999),
DCA3HeXY(-999),
//DistanceToPV3HeXY(-999),
DeviationFromPV3HeXY(-999),
NCrossedRowsTPC3He(-999),
NPIDClusterTPC3He(-999),
TPCMom3He(-999),
TPCnSigma3He(-999),
TPCnSigma3H(-999),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
PIDForTrackingHe3(-999),
DistanceToSecVertHe(-999),
DeviationToSecVertHe(-999),
///three body decay
CandidateTree_3Body(nullptr),
mass_Deuteron_Proton(-999),
mass_Proton_Pion(-999),
CosPointingAngle_Deuteron_Proton(-999),
DistanceOfDaughters_Deuteron_Proton(-999),
DeviationOfDaughters_Deuteron_Proton(-999),
DistanceOfDaughtersXY_Deuteron_Proton(-999),
DeviationOfDaughtersXY_Deuteron_Proton(-999),
DistanceOfDaughters_Deuteron_Pion(-999),
DeviationOfDaughters_Deuteron_Pion(-999),
DistanceOfDaughtersXY_Deuteron_Pion(-999),
DeviationOfDaughtersXY_Deuteron_Pion(-999),
DistanceOfDaughters_Proton_Pion(-999),
DeviationOfDaughters_Proton_Pion(-999),
DistanceOfDaughtersXY_Proton_Pion(-999),
DeviationOfDaughtersXY_Proton_Pion(-999),
OpeningAngle_Pion_Proton(-999),
OpeningAngle_Pion_Deuteron(-999),
OpeningAngle_Proton_Deuteron(-999),
pxDeuteron(-999),
pyDeuteron(-999),
pzDeuteron(-999),
ChargeDeuteron(-999),
DCADeuteron(-999),
//DistanceToPVDeuteron(-999),
DeviationFromPVDeuteron(-999),
DCADeuteronXY(-999),
//DistanceToPVDeuteronXY(-999),
DeviationFromPVDeuteronXY(-999),
NCrossedRowsTPCDeuteron(-999),
NPIDClusterTPCDeuteron (-999),
TPCMomDeuteron(-999),
TPCnSigmaDeuteron(-999),
TOFnSigmaDeuteron(-999),
HasPointOnITSLayer0Deuteron(false),
HasPointOnITSLayer1Deuteron(false),
HasPointOnITSLayer2Deuteron(false),
HasPointOnITSLayer3Deuteron(false),
HasPointOnITSLayer4Deuteron(false),
HasPointOnITSLayer5Deuteron(false),
PIDForTrackingDeuteron(-999),
DistanceToSecVertDeuteron(-999),
DeviationToSecVertDeuteron(-999),
pxProton(-999),
pyProton(-999),
pzProton(-999),
ChargeProton(-999),
DCAProton(-999),
//DistanceToPVProton(-999),
DeviationFromPVProton(-999),
DCAProtonXY(-999),
//DistanceToPVProtonXY(-999),
DeviationFromPVProtonXY(-999),
NCrossedRowsTPCProton(-999),
NPIDClusterTPCProton (-999),
TPCMomProton(-999),
TPCnSigmaProton(-999),
TOFnSigmaProton(-999),
HasPointOnITSLayer0Proton(false),
HasPointOnITSLayer1Proton(false),
HasPointOnITSLayer2Proton(false),
HasPointOnITSLayer3Proton(false),
HasPointOnITSLayer4Proton(false),
HasPointOnITSLayer5Proton(false),
PIDForTrackingProton(-999),
DistanceToSecVertProton(-999),
DeviationToSecVertProton(-999),
///Monte Carlo
pxMC(-999),
pyMC(-999),
pzMC(-999),
pxHeMC(-999),
pyHeMC(-999),
pzHeMC(-999),
pxPionMC(-999),
pyPionMC(-999),
pzPionMC(-999),
pxProtonMC(-999),
pyProtonMC(-999),
pzProtonMC(-999),
pxDeuteronMC(-999),
pyDeuteronMC(-999),
pzDeuteronMC(-999),
pxVariance(-999),
pyVariance(-999),
pzVariance(-999),
pxHeVariance(-999),
pyHeVariance(-999),
pzHeVariance(-999),
pxPionVariance(-999),
pyPionVariance(-999),
pzPionVariance(-999),
pxProtonVariance(-999),
pyProtonVariance(-999),
pzProtonVariance(-999),
pxDeuteronVariance(-999),
pyDeuteronVariance(-999),
pzDeuteronVariance(-999),
pMC(-999),
pTMC(-999),
RapidityMC(-999),
DecayLengthMC(-999),
DecayLengthXYMC(-999),
xSecVertexMC(-999),
ySecVertexMC(-999),
zSecVertexMC(-999),
xSecVertex(-999),
ySecVertex(-999),
zSecVertex(-999),
xSecVertexVariance(-999),
ySecVertexVariance(-999),
zSecVertexVariance(-999),
NumberOfDeuteronProtonCandidates(-1),
GeneratedTreeMC(nullptr),
GeneratedTreeMC_3Body(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistNsigmaTPCvsPDeuteron(nullptr),
fHistNsigmaTPCvsPProton(nullptr),
fHistNsigmaTOFvsPDeuteron(nullptr),
fHistNsigmaTOFvsPProton(nullptr),
fHistpxTrueRecHe3(nullptr),
fHistpyTrueRecHe3(nullptr),
fHistpzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr),
fHistMomDeuteron(nullptr),
fHistMomProton(nullptr),
fHistMomPion_3Body(nullptr),
fHistNDeuteronsPerEvent(nullptr),
fHistNDeuteronProtonCandidatesPerEvent(nullptr)
{
  /// constructor
  /// define the input of the analysis: in this case we take a 'chain' of events
  /// this chain is created by the analysis manager, so no need to worry about it,
  /// it does its work automatically
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  /// trees
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  /// if run 2-body and 3-body in MC at the same time
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
}
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTree::~AliAnalysisTaskHypertritonKFTree()
{
  /// destructor
  /// at the end of your task, it is deleted from memory by calling this function
  if(fESDtrackCuts) delete fESDtrackCuts;
  if(fQAList) delete fQAList;
  if(fOutputList) delete fOutputList;
  if(CandidateTree) delete CandidateTree;
  if(GeneratedTreeMC) delete GeneratedTreeMC;
  
  if(CandidateTree_3Body) delete CandidateTree_3Body;
  if(GeneratedTreeMC_3Body) delete GeneratedTreeMC_3Body;
}

///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::UserCreateOutputObjects()
{
  /// create output objects
  ///
  /// this function is called ONCE at the start of your analysis (RUNTIME)
  /// here you ceate the histograms that you want to use
  ///
  /// the histograms are in this case added to a tlist, this list is in the end saved
  /// to an output file
  ///
  
  /// Runs which are bad for TPC PID due to a problem in the TPC gain in one or two sectors towards the end of the run
  //  fEventCuts.UseTimeRangeCut(); /// at the moment implemented by hand
  
  /// Create object for track selection
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(50);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5.);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  
  fQAList = new TList();          /// list which will contain all of QA histograms
  fQAList->SetOwner(kTRUE);
  
  fOutputList = new TList();     /// list which will contain all of your histograms at the end of the analysis, the contents of this list are written to the output file
  fOutputList->SetOwner(kTRUE);
  
  hNumberOfEvents = new TH1F("hNumberOfEvents","Events at different selection steps ",7,0,7);
  hNumberOfEvents->GetXaxis()->SetTitle("Selection step");
  hNumberOfEvents->Sumw2();
  fOutputList->Add(hNumberOfEvents);
  
  histoEventCentrality = new TH1I("histoEventCentrality","Events vs centrality percentile (V0M)",100,0,100);
  histoEventCentrality->GetXaxis()->SetTitle("Centrality percentile");
  histoEventCentrality->Sumw2();
  fOutputList->Add(histoEventCentrality);
  
  /// QA histograms
  if (kDoQA) {
    /// Add event selection QA plots
    fEventCuts.AddQAplotsToList(fQAList);
    /// Add own QA plots
    fHistNsigmaTPCvsPPion = new TH2F("fHistNsigmaTPCvsPPion","#it{N}#sigma^{TPC}_{#pi} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{#pi}", 100,0,10,100,-5.,5.);
    fHistNsigmaTPCvsPPion->Sumw2();
    fQAList->Add(fHistNsigmaTPCvsPPion);
    
    if (kRun2BodyDecay) {
      fHistNsigmaTPCvsP3He = new TH2F("fHistNsigmaTPCvsP3He","#it{N}#sigma^{TPC}_{^{3}He} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{^{3}He}", 100,0,10,100,-5.,5.);
      fHistNsigmaTPCvsP3He->Sumw2();
      fQAList->Add(fHistNsigmaTPCvsP3He);
    }
    
    if (kRun3BodyDecay) {
      fHistNsigmaTPCvsPDeuteron = new TH2F("fHistNsigmaTPCvsPDeuteron","#it{N}#sigma^{TPC}_{d} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{d}", 100,0,10,100,-5.,5.);
      fHistNsigmaTPCvsPDeuteron->Sumw2();
      fQAList->Add(fHistNsigmaTPCvsPDeuteron);
      
      fHistNsigmaTPCvsPProton = new TH2F("fHistNsigmaTPCvsPProton","#it{N}#sigma^{TPC}_{p} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{p}", 100,0,10,100,-5.,5.);
      fHistNsigmaTPCvsPProton->Sumw2();
      fQAList->Add(fHistNsigmaTPCvsPProton);
      
      fHistNsigmaTOFvsPDeuteron = new TH2F("fHistNsigmaTOFvsPDeuteron","#it{N}#sigma^{TOF}_{d} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TOF}_{d}", 100,0,10,100,-5.,5.);
      fHistNsigmaTOFvsPDeuteron->Sumw2();
      fQAList->Add(fHistNsigmaTOFvsPDeuteron);
      
      fHistNsigmaTOFvsPProton = new TH2F("fHistNsigmaTOFvsPProton","#it{N}#sigma^{TOF}_{p} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TOF}_{p}", 100,0,10,100,-5.,5.);
      fHistNsigmaTOFvsPProton->Sumw2();
      fQAList->Add(fHistNsigmaTOFvsPProton);
      
      fHistNDeuteronsPerEvent = new TH1I("fHistNDeuteronsPerEvent","Number of deuterons per event",11,-0.5,10.5);
      fHistNDeuteronsPerEvent->Sumw2();
      fQAList->Add(fHistNDeuteronsPerEvent);
      
      fHistNDeuteronProtonCandidatesPerEvent = new TH1I("fHistNDeuteronProtonCandidatesPerEvent","Number of deuteron+proton candidates per event",101,-0.5,100.5);
      fHistNDeuteronProtonCandidatesPerEvent->Sumw2();
      fQAList->Add(fHistNDeuteronProtonCandidatesPerEvent);
      
    }
  }
  
  /// Data
  if (kRun2BodyDecay) {
    ///  Tree with hypertriton candidates
    CandidateTree = new TTree("CandidateTree","CandidateTree"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
    CandidateTree->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    CandidateTree->Branch("mass",&mass,"mass/F");
    CandidateTree->Branch("ErrorMass",&ErrorMass,"ErrorMass/F");
    CandidateTree->Branch("px",&px,"px/F");
    CandidateTree->Branch("py",&py,"py/F");
    CandidateTree->Branch("pz",&pz,"pz/F");
    CandidateTree->Branch("Rapidity",&Rapidity,"Rapidity/F");
    CandidateTree->Branch("Charge",&Charge,"Charge/I");
    CandidateTree->Branch("Chi2PerNDF",&Chi2PerNDF,"Chi2PerNDF/F");
    CandidateTree->Branch("CosPointingAngle",&CosPointingAngle,"CosPointingAngle/F");
    CandidateTree->Branch("DistanceToPV",&DistanceToPV,"DistanceToPV/F");
    CandidateTree->Branch("DistanceToPVXY",&DistanceToPVXY,"DistanceToPVXY/F");
    CandidateTree->Branch("DeviationFromPV",&DeviationFromPV,"DeviationFromPV/F");
    CandidateTree->Branch("DeviationFromPVXY",&DeviationFromPVXY,"DeviationFromPVXY/F");
    
    CandidateTree->Branch("massTopo",&massTopo,"massTopo/F");
    CandidateTree->Branch("ErrorMassTopo",&ErrorMassTopo,"ErrorMassTopo/F");
    CandidateTree->Branch("pxTopo",&pxTopo,"pxTopo/F");
    CandidateTree->Branch("pyTopo",&pyTopo,"pyTopo/F");
    CandidateTree->Branch("pzTopo",&pzTopo,"pzTopo/F");
    CandidateTree->Branch("RapidityTopo",&RapidityTopo,"RapidityTopo/F");
    CandidateTree->Branch("Chi2PerNDFTopo",&Chi2PerNDFTopo,"Chi2PerNDFTopo/F");
    CandidateTree->Branch("CosPointingAngleTopo",&CosPointingAngleTopo,"CosPointingAngleTopo/F");
    CandidateTree->Branch("DistanceToPVTopo",&DistanceToPVTopo,"DistanceToPVTopo/F");
    CandidateTree->Branch("DistanceToPVXYTopo",&DistanceToPVXYTopo,"DistanceToPVXYTopo/F");
    CandidateTree->Branch("DeviationFromPVTopo",&DeviationFromPVTopo,"DeviationFromPVTopo/F");
    CandidateTree->Branch("DeviationFromPVXYTopo",&DeviationFromPVXYTopo,"DeviationFromPVXYTopo/F");
    
    CandidateTree->Branch("DecayLength",&DecayLength,"DecayLength/F");
    CandidateTree->Branch("ErrorDecayLength",&ErrorDecayLength,"ErrorDecayLength/F");
    CandidateTree->Branch("DecayLengthXY",&DecayLengthXY,"DecayLengthXY/F");
    CandidateTree->Branch("ErrorDecayLengthXY",&ErrorDecayLengthXY,"ErrorDecayLengthXY/F");
    
    /// Daughter variables
    CandidateTree->Branch("DistanceOfDaughters",&DistanceOfDaughters,"DistanceOfDaughters/F");
    CandidateTree->Branch("DistanceOfDaughtersXY",&DistanceOfDaughtersXY,"DistanceOfDaughtersXY/F");
    
    CandidateTree->Branch("DeviationOfDaughters",&DeviationOfDaughters,"DeviationOfDaughters/F");
    CandidateTree->Branch("DeviationOfDaughtersXY",&DeviationOfDaughtersXY,"DeviationOfDaughtersXY/F");
    
    CandidateTree->Branch("pxPion",&pxPion,"pxPion/F");
    CandidateTree->Branch("pyPion",&pyPion,"pyPion/F");
    CandidateTree->Branch("pzPion",&pzPion,"pzPion/F");
    CandidateTree->Branch("ChargePion",&ChargePion,"ChargePion/I");
    CandidateTree->Branch("DCAPion",&DCAPion,"DCAPion/F");
    //    CandidateTree->Branch("DistanceToPVPion",&DistanceToPVPion,"DistanceToPVPion/F");
    CandidateTree->Branch("DeviationFromPVPion",&DeviationFromPVPion,"DeviationFromPVPion/F");
    CandidateTree->Branch("DCAPionXY",&DCAPionXY,"DCAPionXY/F");
    //    CandidateTree->Branch("DistanceToPVPionXY",&DistanceToPVPionXY,"DistanceToPVPionXY/F");
    CandidateTree->Branch("DeviationFromPVPionXY",&DeviationFromPVPionXY,"DeviationFromPVPionXY/F");
    CandidateTree->Branch("NCrossedRowsTPCPion",&NCrossedRowsTPCPion,"NCrossedRowsTPCPion/I");
    CandidateTree->Branch("NPIDClusterTPCPion",&NPIDClusterTPCPion,"NPIDClusterTPCPion/I");
    CandidateTree->Branch("TPCMomPion",&TPCMomPion,"TPCMomPion/F");
    CandidateTree->Branch("TPCnSigmaPion",&TPCnSigmaPion,"TPCnSigmaPion/F");
    CandidateTree->Branch("HasPointOnITSLayer0Pion",&HasPointOnITSLayer0Pion,"HasPointOnITSLayer0Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer1Pion",&HasPointOnITSLayer1Pion,"HasPointOnITSLayer1Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer2Pion",&HasPointOnITSLayer2Pion,"HasPointOnITSLayer2Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer3Pion",&HasPointOnITSLayer3Pion,"HasPointOnITSLayer3Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer4Pion",&HasPointOnITSLayer4Pion,"HasPointOnITSLayer4Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer5Pion",&HasPointOnITSLayer5Pion,"HasPointOnITSLayer5Pion/O");
    CandidateTree->Branch("PIDForTrackingPion",&PIDForTrackingPion,"PIDForTrackingPion/I");
    CandidateTree->Branch("DistanceToSecVertPion",&DistanceToSecVertPion,"DistanceToSecVertPion/F");
    CandidateTree->Branch("DeviationToSecVertPion",&DeviationToSecVertPion,"DeviationToSecVertPion/F");
    
    CandidateTree->Branch("pxHe",&pxHe,"pxHe/F");
    CandidateTree->Branch("pyHe",&pyHe,"pyHe/F");
    CandidateTree->Branch("pzHe",&pzHe,"pzHe/F");
    CandidateTree->Branch("ChargeHe",&ChargeHe,"ChargeHe/I");
    CandidateTree->Branch("DCA3He",&DCA3He,"DCA3He/F");
    //    CandidateTree->Branch("DistanceToPV3He",&DistanceToPV3He,"DistanceToPV3He/F");
    CandidateTree->Branch("DeviationFromPV3He",&DeviationFromPV3He,"DeviationFromPV3He/F");
    CandidateTree->Branch("DCA3HeXY",&DCA3HeXY,"DCA3HeXY/F");
    //    CandidateTree->Branch("DistanceToPV3HeXY",&DistanceToPV3HeXY,"DistanceToPV3HeXY/F");
    CandidateTree->Branch("DeviationFromPV3HeXY",&DeviationFromPV3HeXY,"DeviationFromPV3HeXY/F");
    CandidateTree->Branch("NCrossedRowsTPC3He",&NCrossedRowsTPC3He,"NCrossedRowsTPC3He/I");
    CandidateTree->Branch("NPIDClusterTPC3He",&NPIDClusterTPC3He,"NPIDClusterTPC3He/I");
    CandidateTree->Branch("TPCMom3He",&TPCMom3He,"TPCMom3He/F");
    CandidateTree->Branch("TPCnSigma3He",&TPCnSigma3He,"TPCnSigma3He/F");
    CandidateTree->Branch("TPCnSigma3H",&TPCnSigma3H,"TPCnSigma3H/F");
    CandidateTree->Branch("HasPointOnITSLayer0He3",&HasPointOnITSLayer0He3,"HasPointOnITSLayer0He3/O");
    CandidateTree->Branch("HasPointOnITSLayer1He3",&HasPointOnITSLayer1He3,"HasPointOnITSLayer1He3/O");
    CandidateTree->Branch("HasPointOnITSLayer2He3",&HasPointOnITSLayer2He3,"HasPointOnITSLayer2He3/O");
    CandidateTree->Branch("HasPointOnITSLayer3He3",&HasPointOnITSLayer3He3,"HasPointOnITSLayer3He3/O");
    CandidateTree->Branch("HasPointOnITSLayer4He3",&HasPointOnITSLayer4He3,"HasPointOnITSLayer4He3/O");
    CandidateTree->Branch("HasPointOnITSLayer5He3",&HasPointOnITSLayer5He3,"HasPointOnITSLayer5He3/O");
    CandidateTree->Branch("PIDForTrackingHe3",&PIDForTrackingHe3,"PIDForTrackingHe3/I");
    CandidateTree->Branch("DistanceToSecVertHe",&DistanceToSecVertHe,"DistanceToSecVertHe/F");
    CandidateTree->Branch("DeviationToSecVertHe",&DeviationToSecVertHe,"DeviationToSecVertHe/F");
    
    CandidateTree->Branch("xSecVertex",&xSecVertex,"xSecVertex/F");
    CandidateTree->Branch("ySecVertex",&ySecVertex,"ySecVertex/F");
    CandidateTree->Branch("zSecVertex",&zSecVertex,"zSecVertex/F");
    
    CandidateTree->Branch("xSecVertexVariance",&xSecVertexVariance,"xSecVertexVariance/F");
    CandidateTree->Branch("ySecVertexVariance",&ySecVertexVariance,"ySecVertexVariance/F");
    CandidateTree->Branch("zSecVertexVariance",&zSecVertexVariance,"zSecVertexVariance/F");
    
    if (kIsMC){
      ///MC truth information
      /// true momentun of the hypertriton candidate
      CandidateTree->Branch("pxMC",&pxMC,"pxMC/F");
      CandidateTree->Branch("pyMC",&pyMC,"pyMC/F");
      CandidateTree->Branch("pzMC",&pzMC,"pzMC/F");
      /// true momentum of teh dauther tracks
      CandidateTree->Branch("pxHeMC",&pxHeMC,"pxHeMC/F");
      CandidateTree->Branch("pyHeMC",&pyHeMC,"pyHeMC/F");
      CandidateTree->Branch("pzHeMC",&pzHeMC,"pzHeMC/F");
      CandidateTree->Branch("pxPionMC",&pxPionMC,"pxPionMC/F");
      CandidateTree->Branch("pyPionMC",&pyPionMC,"pyPionMC/F");
      CandidateTree->Branch("pzPionMC",&pzPionMC,"pzPionMC/F");
      
      /// Variance of the momentun of the hypertriton candidate
      CandidateTree->Branch("pxVariance",&pxVariance,"pxVariance/F");
      CandidateTree->Branch("pyVariance",&pyVariance,"pyVariance/F");
      CandidateTree->Branch("pzVariance",&pzVariance,"pzVariance/F");
      /// Variance of the momentum of the dauther tracks
      CandidateTree->Branch("pxHeVariance",&pxHeVariance,"pxHeVariance/F");
      CandidateTree->Branch("pyHeVariance",&pyHeVariance,"pyHeVariance/F");
      CandidateTree->Branch("pzHeVariance",&pzHeVariance,"pzHeVariance/F");
      CandidateTree->Branch("pxPionVariance",&pxPionVariance,"pxPionVariance/F");
      CandidateTree->Branch("pyPionVariance",&pyPionVariance,"pyPionVariance/F");
      CandidateTree->Branch("pzPionVariance",&pzPionVariance,"pzPionVariance/F");
      
      CandidateTree->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      CandidateTree->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
      CandidateTree->Branch("DecayLengthXYMC",&DecayLengthXYMC,"DecayLengthXYMC/F");
      
      CandidateTree->Branch("xSecVertexMC",&xSecVertexMC,"xSecVertexMC/F");
      CandidateTree->Branch("ySecVertexMC",&ySecVertexMC,"ySecVertexMC/F");
      CandidateTree->Branch("zSecVertexMC",&zSecVertexMC,"zSecVertexMC/F");
    }
  }
  
  if (kRun3BodyDecay) {
    ///  Tree with hypertriton candidates
    CandidateTree_3Body = new TTree("CandidateTree_3Body","CandidateTree_3Body"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
    CandidateTree_3Body->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    CandidateTree_3Body->Branch("mass",&mass,"mass/F");
    CandidateTree_3Body->Branch("ErrorMass",&ErrorMass,"ErrorMass/F");
    CandidateTree_3Body->Branch("px",&px,"px/F");
    CandidateTree_3Body->Branch("py",&py,"py/F");
    CandidateTree_3Body->Branch("pz",&pz,"pz/F");
    CandidateTree_3Body->Branch("Rapidity",&Rapidity,"Rapidity/F");
    CandidateTree_3Body->Branch("Charge",&Charge,"Charge/I");
    CandidateTree_3Body->Branch("Chi2PerNDF",&Chi2PerNDF,"Chi2PerNDF/F");
    CandidateTree_3Body->Branch("CosPointingAngle",&CosPointingAngle,"CosPointingAngle/F");
    CandidateTree_3Body->Branch("DistanceToPV",&DistanceToPV,"DistanceToPV/F");
    CandidateTree_3Body->Branch("DistanceToPVXY",&DistanceToPVXY,"DistanceToPVXY/F");
    CandidateTree_3Body->Branch("DeviationFromPV",&DeviationFromPV,"DeviationFromPV/F");
    CandidateTree_3Body->Branch("DeviationFromPVXY",&DeviationFromPVXY,"DeviationFromPVXY/F");
    
    CandidateTree_3Body->Branch("massTopo",&massTopo,"massTopo/F");
    CandidateTree_3Body->Branch("ErrorMassTopo",&ErrorMassTopo,"ErrorMassTopo/F");
    CandidateTree_3Body->Branch("pxTopo",&pxTopo,"pxTopo/F");
    CandidateTree_3Body->Branch("pyTopo",&pyTopo,"pyTopo/F");
    CandidateTree_3Body->Branch("pzTopo",&pzTopo,"pzTopo/F");
    CandidateTree_3Body->Branch("RapidityTopo",&RapidityTopo,"RapidityTopo/F");
    CandidateTree_3Body->Branch("Chi2PerNDFTopo",&Chi2PerNDFTopo,"Chi2PerNDFTopo/F");
    CandidateTree_3Body->Branch("CosPointingAngleTopo",&CosPointingAngleTopo,"CosPointingAngleTopo/F");
    CandidateTree_3Body->Branch("DistanceToPVTopo",&DistanceToPVTopo,"DistanceToPVTopo/F");
    CandidateTree_3Body->Branch("DistanceToPVXYTopo",&DistanceToPVXYTopo,"DistanceToPVXYTopo/F");
    CandidateTree_3Body->Branch("DeviationFromPVTopo",&DeviationFromPVTopo,"DeviationFromPVTopo/F");
    CandidateTree_3Body->Branch("DeviationFromPVXYTopo",&DeviationFromPVXYTopo,"DeviationFromPVXYTopo/F");
    
    CandidateTree_3Body->Branch("DecayLength",&DecayLength,"DecayLength/F");
    CandidateTree_3Body->Branch("ErrorDecayLength",&ErrorDecayLength,"ErrorDecayLength/F");
    CandidateTree_3Body->Branch("DecayLengthXY",&DecayLengthXY,"DecayLengthXY/F");
    CandidateTree_3Body->Branch("ErrorDecayLengthXY",&ErrorDecayLengthXY,"ErrorDecayLengthXY/F");
    
    CandidateTree_3Body->Branch("mass_Deuteron_Proton",&mass_Deuteron_Proton,"mass_Deuteron_Proton/F");
    CandidateTree_3Body->Branch("mass_Proton_Pion",&mass_Proton_Pion,"mass_Proton_Pion/F");
    CandidateTree_3Body->Branch("CosPointingAngle_Deuteron_Proton",&CosPointingAngle_Deuteron_Proton,"CosPointingAngle_Deuteron_Proton/F");
    
    /// Daughter variables
    CandidateTree_3Body->Branch("DistanceOfDaughters_Deuteron_Proton",&DistanceOfDaughters_Deuteron_Proton,"DistanceOfDaughters_Deuteron_Proton/F");
    CandidateTree_3Body->Branch("DistanceOfDaughtersXY_Deuteron_Proton",&DistanceOfDaughtersXY_Deuteron_Proton,"DistanceOfDaughtersXY_Deuteron_Proton/F");
    CandidateTree_3Body->Branch("DeviationOfDaughters_Deuteron_Proton",&DeviationOfDaughters_Deuteron_Proton,"DeviationOfDaughters_Deuteron_Proton/F");
    CandidateTree_3Body->Branch("DeviationOfDaughtersXY_Deuteron_Proton",&DeviationOfDaughtersXY_Deuteron_Proton,"DeviationOfDaughtersXY_Deuteron_Proton/F");
    
    CandidateTree_3Body->Branch("DistanceOfDaughters_Deuteron_Pion",&DistanceOfDaughters_Deuteron_Pion,"DistanceOfDaughters_Deuteron_Pion/F");
    CandidateTree_3Body->Branch("DistanceOfDaughtersXY_Deuteron_Pion",&DistanceOfDaughtersXY_Deuteron_Pion,"DistanceOfDaughtersXY_Deuteron_Pion/F");
    CandidateTree_3Body->Branch("DeviationOfDaughters_Deuteron_Pion",&DeviationOfDaughters_Deuteron_Pion,"DeviationOfDaughters_Deuteron_Pion/F");
    CandidateTree_3Body->Branch("DeviationOfDaughtersXY_Deuteron_Pion",&DeviationOfDaughtersXY_Deuteron_Pion,"DeviationOfDaughtersXY_Deuteron_Pion/F");
    
    CandidateTree_3Body->Branch("DistanceOfDaughters_Proton_Pion",&DistanceOfDaughters_Proton_Pion,"DistanceOfDaughters_Proton_Pion/F");
    CandidateTree_3Body->Branch("DistanceOfDaughtersXY_Proton_Pion",&DistanceOfDaughtersXY_Proton_Pion,"DistanceOfDaughtersXY_Proton_Pion/F");
    CandidateTree_3Body->Branch("DeviationOfDaughters_Proton_Pion",&DeviationOfDaughters_Proton_Pion,"DeviationOfDaughters_Proton_Pion/F");
    CandidateTree_3Body->Branch("DeviationOfDaughtersXY_Proton_Pion",&DeviationOfDaughtersXY_Proton_Pion,"DeviationOfDaughtersXY_Proton_Pion/F");
    
    CandidateTree_3Body->Branch("OpeningAngle_Pion_Proton",&OpeningAngle_Pion_Proton,"OpeningAngle_Pion_Proton/F");
    CandidateTree_3Body->Branch("OpeningAngle_Pion_Deuteron",&OpeningAngle_Pion_Deuteron,"OpeningAngle_Pion_Deuteron/F");
    CandidateTree_3Body->Branch("OpeningAngle_Proton_Deuteron",&OpeningAngle_Proton_Deuteron,"OpeningAngle_Proton_Deuteron/F");
    
    CandidateTree_3Body->Branch("pxDeuteron",&pxDeuteron,"pxDeuteron/F");
    CandidateTree_3Body->Branch("pyDeuteron",&pyDeuteron,"pyDeuteron/F");
    CandidateTree_3Body->Branch("pzDeuteron",&pzDeuteron,"pzDeuteron/F");
    CandidateTree_3Body->Branch("ChargeDeuteron",&ChargeDeuteron,"ChargeDeuteron/I");
    CandidateTree_3Body->Branch("DCADeuteron",&DCADeuteron,"DCADeuteron/F");
    //    CandidateTree_3Body->Branch("DistanceToPVDeuteron",&DistanceToPVDeuteron,"DistanceToPVDeuteron/F");
    CandidateTree_3Body->Branch("DeviationFromPVDeuteron",&DeviationFromPVDeuteron,"DeviationFromPVDeuteron/F");
    CandidateTree_3Body->Branch("DCADeuteronXY",&DCADeuteronXY,"DCADeuteronXY/F");
    //    CandidateTree_3Body->Branch("DistanceToPVDeuteronXY",&DistanceToPVDeuteronXY,"DistanceToPVDeuteronXY/F");
    CandidateTree_3Body->Branch("DeviationFromPVDeuteronXY",&DeviationFromPVDeuteronXY,"DeviationFromPVDeuteronXY/F");
    CandidateTree_3Body->Branch("NCrossedRowsTPCDeuteron",&NCrossedRowsTPCDeuteron,"NCrossedRowsTPCDeuteron/I");
    CandidateTree_3Body->Branch("NPIDClusterTPCDeuteron",&NPIDClusterTPCDeuteron,"NPIDClusterTPCDeuteron/I");
    CandidateTree_3Body->Branch("TPCMomDeuteron",&TPCMomDeuteron,"TPCMomDeuteron/F");
    CandidateTree_3Body->Branch("TPCnSigmaDeuteron",&TPCnSigmaDeuteron,"TPCnSigmaDeuteron/F");
    CandidateTree_3Body->Branch("TOFnSigmaDeuteron",&TOFnSigmaDeuteron,"TOFnSigmaDeuteron/F");
    CandidateTree_3Body->Branch("HasPointOnITSLayer0Deuteron",&HasPointOnITSLayer0Deuteron,"HasPointOnITSLayer0Deuteron/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer1Deuteron",&HasPointOnITSLayer1Deuteron,"HasPointOnITSLayer1Deuteron/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer2Deuteron",&HasPointOnITSLayer2Deuteron,"HasPointOnITSLayer2Deuteron/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer3Deuteron",&HasPointOnITSLayer3Deuteron,"HasPointOnITSLayer3Deuteron/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer4Deuteron",&HasPointOnITSLayer4Deuteron,"HasPointOnITSLayer4Deuteron/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer5Deuteron",&HasPointOnITSLayer5Deuteron,"HasPointOnITSLayer5Deuteron/O");
    CandidateTree_3Body->Branch("PIDForTrackingDeuteron",&PIDForTrackingDeuteron,"PIDForTrackingDeuteron/I");
    CandidateTree_3Body->Branch("DistanceToSecVertDeuteron",&DistanceToSecVertDeuteron,"DistanceToSecVertDeuteron/F");
    CandidateTree_3Body->Branch("DeviationToSecVertDeuteron",&DeviationToSecVertDeuteron,"DeviationToSecVertDeuteron/F");
    
    CandidateTree_3Body->Branch("pxProton",&pxProton,"pxProton/F");
    CandidateTree_3Body->Branch("pyProton",&pyProton,"pyProton/F");
    CandidateTree_3Body->Branch("pzProton",&pzProton,"pzProton/F");
    CandidateTree_3Body->Branch("ChargeProton",&ChargeProton,"ChargeProton/I");
    CandidateTree_3Body->Branch("DCAProton",&DCAProton,"DCAProton/F");
    //    CandidateTree_3Body->Branch("DistanceToPVProton",&DistanceToPVProton,"DistanceToPVProton/F");
    CandidateTree_3Body->Branch("DeviationFromPVProton",&DeviationFromPVProton,"DeviationFromPVProton/F");
    CandidateTree_3Body->Branch("DCAProtonXY",&DCAProtonXY,"DCAProtonXY/F");
    //    CandidateTree_3Body->Branch("DistanceToPVProtonXY",&DistanceToPVProtonXY,"DistanceToPVProtonXY/F");
    CandidateTree_3Body->Branch("DeviationFromPVProtonXY",&DeviationFromPVProtonXY,"DeviationFromPVProtonXY/F");
    CandidateTree_3Body->Branch("NCrossedRowsTPCProton",&NCrossedRowsTPCProton,"NCrossedRowsTPCProton/I");
    CandidateTree_3Body->Branch("NPIDClusterTPCProton",&NPIDClusterTPCProton,"NPIDClusterTPCProton/I");
    CandidateTree_3Body->Branch("TPCMomProton",&TPCMomProton,"TPCMomProton/F");
    CandidateTree_3Body->Branch("TPCnSigmaProton",&TPCnSigmaProton,"TPCnSigmaProton/F");
    CandidateTree_3Body->Branch("TOFnSigmaProton",&TOFnSigmaProton,"TOFnSigmaProton/F");
    CandidateTree_3Body->Branch("HasPointOnITSLayer0Proton",&HasPointOnITSLayer0Proton,"HasPointOnITSLayer0Proton/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer1Proton",&HasPointOnITSLayer1Proton,"HasPointOnITSLayer1Proton/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer2Proton",&HasPointOnITSLayer2Proton,"HasPointOnITSLayer2Proton/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer3Proton",&HasPointOnITSLayer3Proton,"HasPointOnITSLayer3Proton/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer4Proton",&HasPointOnITSLayer4Proton,"HasPointOnITSLayer4Proton/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer5Proton",&HasPointOnITSLayer5Proton,"HasPointOnITSLayer5Proton/O");
    CandidateTree_3Body->Branch("PIDForTrackingProton",&PIDForTrackingProton,"PIDForTrackingProton/I");
    CandidateTree_3Body->Branch("DistanceToSecVertProton",&DistanceToSecVertProton,"DistanceToSecVertProton/F");
    CandidateTree_3Body->Branch("DeviationToSecVertProton",&DeviationToSecVertProton,"DeviationToSecVertProton/F");
    
    CandidateTree_3Body->Branch("pxPion",&pxPion,"pxPion/F");
    CandidateTree_3Body->Branch("pyPion",&pyPion,"pyPion/F");
    CandidateTree_3Body->Branch("pzPion",&pzPion,"pzPion/F");
    CandidateTree_3Body->Branch("ChargePion",&ChargePion,"ChargePion/I");
    CandidateTree_3Body->Branch("DCAPion",&DCAPion,"DCAPion/F");
    //    CandidateTree_3Body->Branch("DistanceToPVPion",&DistanceToPVPion,"DistanceToPVPion/F");
    CandidateTree_3Body->Branch("DeviationFromPVPion",&DeviationFromPVPion,"DeviationFromPVPion/F");
    CandidateTree_3Body->Branch("DCAPionXY",&DCAPionXY,"DCAPionXY/F");
    //    CandidateTree_3Body->Branch("DistanceToPVPionXY",&DistanceToPVPionXY,"DistanceToPVPionXY/F");
    CandidateTree_3Body->Branch("DeviationFromPVPionXY",&DeviationFromPVPionXY,"DeviationFromPVPionXY/F");
    CandidateTree_3Body->Branch("NCrossedRowsTPCPion",&NCrossedRowsTPCPion,"NCrossedRowsTPCPion/I");
    CandidateTree_3Body->Branch("NPIDClusterTPCPion",&NPIDClusterTPCPion,"NPIDClusterTPCPion/I");
    CandidateTree_3Body->Branch("TPCMomPion",&TPCMomPion,"TPCMomPion/F");
    CandidateTree_3Body->Branch("TPCnSigmaPion",&TPCnSigmaPion,"TPCnSigmaPion/F");
    CandidateTree_3Body->Branch("HasPointOnITSLayer0Pion",&HasPointOnITSLayer0Pion,"HasPointOnITSLayer0Pion/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer1Pion",&HasPointOnITSLayer1Pion,"HasPointOnITSLayer1Pion/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer2Pion",&HasPointOnITSLayer2Pion,"HasPointOnITSLayer2Pion/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer3Pion",&HasPointOnITSLayer3Pion,"HasPointOnITSLayer3Pion/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer4Pion",&HasPointOnITSLayer4Pion,"HasPointOnITSLayer4Pion/O");
    CandidateTree_3Body->Branch("HasPointOnITSLayer5Pion",&HasPointOnITSLayer5Pion,"HasPointOnITSLayer5Pion/O");
    CandidateTree_3Body->Branch("PIDForTrackingPion",&PIDForTrackingPion,"PIDForTrackingPion/I");
    CandidateTree_3Body->Branch("DistanceToSecVertPion",&DistanceToSecVertPion,"DistanceToSecVertPion/F");
    CandidateTree_3Body->Branch("DeviationToSecVertPion",&DeviationToSecVertPion,"DeviationToSecVertPion/F");
    
    CandidateTree_3Body->Branch("xSecVertex",&xSecVertex,"xSecVertex/F");
    CandidateTree_3Body->Branch("ySecVertex",&ySecVertex,"ySecVertex/F");
    CandidateTree_3Body->Branch("zSecVertex",&zSecVertex,"zSecVertex/F");
    
    CandidateTree_3Body->Branch("xSecVertexVariance",&xSecVertexVariance,"xSecVertexVariance/F");
    CandidateTree_3Body->Branch("ySecVertexVariance",&ySecVertexVariance,"ySecVertexVariance/F");
    CandidateTree_3Body->Branch("zSecVertexVariance",&zSecVertexVariance,"zSecVertexVariance/F");
    
    if (kIsMC){
      ///MC truth information
      /// true momentun of the hypertriton candidate
      ///Mc truth information
      /// true momentun of the hypertriton candidate
      CandidateTree_3Body->Branch("pxMC",&pxMC,"pxMC/F");
      CandidateTree_3Body->Branch("pyMC",&pyMC,"pyMC/F");
      CandidateTree_3Body->Branch("pzMC",&pzMC,"pzMC/F");
      /// true momentum of the dauther tracks
      CandidateTree_3Body->Branch("pxDeuteronMC",&pxDeuteronMC,"pxDeuteronMC/F");
      CandidateTree_3Body->Branch("pyDeuteronMC",&pyDeuteronMC,"pyDeuteronMC/F");
      CandidateTree_3Body->Branch("pzDeuteronMC",&pzDeuteronMC,"pzDeuteronMC/F");
      CandidateTree_3Body->Branch("pxProtonMC",&pxProtonMC,"pxProtonMC/F");
      CandidateTree_3Body->Branch("pyProtonMC",&pyProtonMC,"pyProtonMC/F");
      CandidateTree_3Body->Branch("pzProtonMC",&pzProtonMC,"pzProtonMC/F");
      CandidateTree_3Body->Branch("pxPionMC",&pxPionMC,"pxPionMC/F");
      CandidateTree_3Body->Branch("pyPionMC",&pyPionMC,"pyPionMC/F");
      CandidateTree_3Body->Branch("pzPionMC",&pzPionMC,"pzPionMC/F");
      
      /// Variance of the momentun of the hypertriton candidate
      CandidateTree_3Body->Branch("pxVariance",&pxVariance,"pxVariance/F");
      CandidateTree_3Body->Branch("pyVariance",&pyVariance,"pyVariance/F");
      CandidateTree_3Body->Branch("pzVariance",&pzVariance,"pzVariance/F");
      /// Variance of the momentum of the dauther tracks
      CandidateTree_3Body->Branch("pxDeuteronVariance",&pxDeuteronVariance,"pxDeuteronVariance/F");
      CandidateTree_3Body->Branch("pyDeuteronVariance",&pyDeuteronVariance,"pyDeuteronVariance/F");
      CandidateTree_3Body->Branch("pzDeuteronVariance",&pzDeuteronVariance,"pzDeuteronVariance/F");
      CandidateTree_3Body->Branch("pxProtonVariance",&pxProtonVariance,"pxProtonVariance/F");
      CandidateTree_3Body->Branch("pyProtonVariance",&pyProtonVariance,"pyProtonVariance/F");
      CandidateTree_3Body->Branch("pzProtonVariance",&pzProtonVariance,"pzProtonVariance/F");
      
      CandidateTree_3Body->Branch("pxPionVariance",&pxPionVariance,"pxPionVariance/F");
      CandidateTree_3Body->Branch("pyPionVariance",&pyPionVariance,"pyPionVariance/F");
      CandidateTree_3Body->Branch("pzPionVariance",&pzPionVariance,"pzPionVariance/F");
      
      CandidateTree_3Body->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      CandidateTree_3Body->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
      CandidateTree_3Body->Branch("DecayLengthXYMC",&DecayLengthXYMC,"DecayLengthXYMC/F");
      
      CandidateTree_3Body->Branch("xSecVertexMC",&xSecVertexMC,"xSecVertexMC/F");
      CandidateTree_3Body->Branch("ySecVertexMC",&ySecVertexMC,"ySecVertexMC/F");
      CandidateTree_3Body->Branch("zSecVertexMC",&zSecVertexMC,"zSecVertexMC/F");
    }
  }
  
  if (kIsMC){
    /// MC output
    if (kRun2BodyDecay) {
      /// 2-body decay
      /// generated MC
      GeneratedTreeMC = new TTree("GeneratedTreeMC","GeneratedTreeMC"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
      GeneratedTreeMC->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
      GeneratedTreeMC->Branch("pMC",&pMC,"pMC/F");
      GeneratedTreeMC->Branch("pTMC",&pTMC,"pTMC/F");
      GeneratedTreeMC->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      GeneratedTreeMC->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
      //      GeneratedTreeMC->Branch("DecayLengthXYMC",&DecayLengthXYMC,"DecayLengthXYMC/F");
    }
    
    if (kRun3BodyDecay) {
      /// 3-body decay
      /// generated MC
      GeneratedTreeMC_3Body = new TTree("GeneratedTreeMC_3Body","GeneratedTreeMC_3Body"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
      GeneratedTreeMC_3Body->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
      GeneratedTreeMC_3Body->Branch("pMC",&pMC,"pMC/F");
      GeneratedTreeMC_3Body->Branch("pTMC",&pTMC,"pTMC/F");
      GeneratedTreeMC_3Body->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      GeneratedTreeMC_3Body->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
      GeneratedTreeMC_3Body->Branch("DecayLengthXYMC",&DecayLengthXYMC,"DecayLengthXYMC/F");
    }
    
    if (kDoQA) {
      /// QA and corrections
      if (kRun2BodyDecay) {
        fHistpxTrueRecHe3 = new TH2F("fHistpxTrueRecHe3", "(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) vs #it{p}_{x}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpxTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{x}^{rec} (GeV/#it{c})");
        fHistpxTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) (GeV/#it{c})");
        fHistpxTrueRecHe3->Sumw2();
        fQAList->Add(fHistpxTrueRecHe3);
        
        fHistpyTrueRecHe3 = new TH2F("fHistpyTrueRecHe3", "(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) vs #it{p}_{y}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpyTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{y}^{rec} (GeV/#it{c})");
        fHistpyTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) (GeV/#it{c})");
        fHistpyTrueRecHe3->Sumw2();
        fQAList->Add(fHistpyTrueRecHe3);
        
        fHistpzTrueRecHe3 = new TH2F("fHistpzTrueRecHe3", "(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) vs #it{p}_{z}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpzTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{z}^{rec} (GeV/#it{c})");
        fHistpzTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) (GeV/#it{c})");
        fHistpzTrueRecHe3->Sumw2();
        fQAList->Add(fHistpzTrueRecHe3);
        
        fHistMomPion = new TH1F("fHistMomPion","Momentum distribution of #pi from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomPion->Sumw2();
        fQAList->Add(fHistMomPion);
        
        fHistMomHe3 = new TH1F("fHistMomHe3","Momentum distribution of ^{3}He from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomHe3->Sumw2();
        fQAList->Add(fHistMomHe3);
      }
      
      if (kRun3BodyDecay){
        fHistMomDeuteron = new TH1F("fHistMomDeuteron","Momentum distribution of deuteron from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomDeuteron->Sumw2();
        fQAList->Add(fHistMomDeuteron);
        
        fHistMomProton = new TH1F("fHistMomProton","Momentum distribution of proton from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomProton->Sumw2();
        fQAList->Add(fHistMomProton);
        
        fHistMomPion_3Body = new TH1F("fHistMomPion_3Body","Momentum distribution of #pi from ^{3}_{#Lambda}H (3-body);#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomPion_3Body->Sumw2();
        fQAList->Add(fHistMomPion_3Body);
      }
    }
  }
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  
  int NumberOfContainer = 3;
  
  if (kRun2BodyDecay) {
    PostData(NumberOfContainer++,CandidateTree);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC);
  }
  if (kRun3BodyDecay) {
    PostData(NumberOfContainer++,CandidateTree_3Body);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC_3Body);
  }
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::UserExec(Option_t *)
{
  /// user exec
  /// this function is called once for each event
  /// the manager will take care of reading the events from file, and with the static function InputEvent() you
  /// have access to the current event.
  /// once you return from the UserExec function, the manager will retrieve the next event from the chain
  
  if ( !PassedEventSelection() ) {
    PostData(1,fQAList);
    PostData(2,fOutputList);
    return;
  }
  
  /// set Magnetic field for KF
  KFParticle::SetField(fInputEvent->GetMagneticField());
  
  /// Create Primary Vertex with KF
  PrimVertex = CreateKFVertex(fInputEvent->GetPrimaryVertex());
  
  if (kIsMC){
    ProcessMC();
  } else {
    ProcessESD();
  }
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  
  int NumberOfContainer = 3;
  if (kRun2BodyDecay) {
    PostData(NumberOfContainer++,CandidateTree);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC);
  }
  if (kRun3BodyDecay) {
    PostData(NumberOfContainer++,CandidateTree_3Body);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC_3Body);
  }
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::Terminate(Option_t *)
{
  /// terminate
  /// called at the END of the analysis (when all events are processed)
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ProcessESD()
{
  /// container for found 3He candidates
  He3CandTrackIdPos.clear();
  He3CandTrackIdNeg.clear();
  
  /// container for found pion candidates
  PionCandTrackIdPos.clear();
  PionCandTrackIdNeg.clear();
  
  /// container for found deuteron candidates
  DeuteronCandTrackIdPos.clear();
  DeuteronCandTrackIdNeg.clear();
  /// container for found proton candidates
  ProtonCandTrackIdPos.clear();
  ProtonCandTrackIdNeg.clear();
  
  /// check for helium-3 candidates and store their track number
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(i));         /// get a track (type AliESDtrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;              /// skip track is it fails basic track selection
    
    if (kDoQA) {
      /// Fill QA histogram
      Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
      if(dEdxSigmaPion < 5) fHistNsigmaTPCvsPPion->Fill(track->P(),dEdxSigmaPion);
      
      if (kRun2BodyDecay){
        Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
        if(dEdxSigmaHe3 < 5) fHistNsigmaTPCvsP3He->Fill(track->P(), dEdxSigmaHe3);
      }
      
      if(kRun3BodyDecay){
        Double_t dEdxSigmaDeuteron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
        if(dEdxSigmaDeuteron < 5) fHistNsigmaTPCvsPDeuteron->Fill(track->P(),dEdxSigmaDeuteron);
        
        Double_t dEdxSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        if(dEdxSigmaProton < 5) fHistNsigmaTPCvsPProton->Fill(track->P(),dEdxSigmaProton);
        
        Double_t TOFSigmaDeuteron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
        if(TOFSigmaDeuteron < 5) fHistNsigmaTOFvsPDeuteron->Fill(track->P(),TOFSigmaDeuteron);
        
        Double_t TOFSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        if(TOFSigmaProton < 5) fHistNsigmaTOFvsPProton->Fill(track->P(),TOFSigmaProton);
      }
    }
    
    Int_t ChargeTrack = (Int_t) track->Charge();
    
    /// Needed for both decay channels
    if(PionSelection(track)){
      if (ChargeTrack > 0) {
        PionCandTrackIdPos.push_back(i);
      } else{
        PionCandTrackIdNeg.push_back(i);
      }
    }
    
    if (kRun2BodyDecay){
      /// 2-body decay
      if(Helium3Selection(track)){
        if (ChargeTrack > 0) {
          He3CandTrackIdPos.push_back(i);
        } else{
          He3CandTrackIdNeg.push_back(i);
        }
      }
    }
    
    if(kRun3BodyDecay){
      /// 3-body decay
      if(DeuteronSelection(track)){
        if (ChargeTrack > 0) {
          DeuteronCandTrackIdPos.push_back(i);
        } else{
          DeuteronCandTrackIdNeg.push_back(i);
        }
      }
      if(ProtonSelection(track)){
        if (ChargeTrack > 0) {
          ProtonCandTrackIdPos.push_back(i);
        } else{
          ProtonCandTrackIdNeg.push_back(i);
        }
      }
    }
  }
  
  if(kRun2BodyDecay){
    /// Loop only over correct charge combinations -> unlike sign pairs
    /// Automatically removes like sign pairs
    TwoBodyDecay(He3CandTrackIdPos,PionCandTrackIdNeg);
    TwoBodyDecay(He3CandTrackIdNeg,PionCandTrackIdPos);
  }
  if(kRun3BodyDecay){
    if (kDoQA) {
      fHistNDeuteronsPerEvent->Fill(DeuteronCandTrackIdPos.size()+DeuteronCandTrackIdNeg.size());
    }
    if(kDoQA) NumberOfDeuteronProtonCandidates = 0;
    /// Loop only over correct charge combinations -> Like sign proton and deuteron + unlike sign pion
    ThreeBodyDecay(DeuteronCandTrackIdPos, ProtonCandTrackIdPos, PionCandTrackIdNeg);
    ThreeBodyDecay(DeuteronCandTrackIdNeg, ProtonCandTrackIdNeg, PionCandTrackIdPos);
    
    if (kDoQA) {
      fHistNDeuteronProtonCandidatesPerEvent->Fill(NumberOfDeuteronProtonCandidates);
    }
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::TwoBodyDecay(const vector<Int_t>& He3CandTrackId, const vector<Int_t>& PionCandTrackId){
  /// next event if there are no 3He or pion candidates
  if ( (Int_t) He3CandTrackId.size() <= 0) return;
  if ( (Int_t) PionCandTrackId.size() <= 0) return;
  
  /// combine pion candidates with 3He candidates
  for(unsigned int iPion=0; iPion < (UInt_t) PionCandTrackId.size(); iPion++) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(PionCandTrackId[iPion]));         /// get a track (type AliESDtrack) from the event
    
    Int_t ChargePionTrack = (Int_t) track->Charge();
    /// create KFParticle for pion
    Int_t pdgPion = 211;
    if (ChargePionTrack < 0) {
      pdgPion = -211;
    }
    KFParticle kfpDaughter2 = CreateKFTrack(track, pdgPion);
    
    /// Fill pion variables
    FillPionVariables(track, kfpDaughter2);
    
    for (unsigned int jHe3=0; jHe3 < (UInt_t) He3CandTrackId.size(); jHe3++) {
      /// Make sure not to use a single track twice
      if ( PionCandTrackId[iPion]==He3CandTrackId[jHe3] ) continue;
      
      AliESDtrack* ESDtrackHe3 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(He3CandTrackId[jHe3]));
      Int_t Charge3He = (Int_t) ESDtrackHe3->Charge();
      /// create KFParticle for 3He
      Int_t pdg3He = 1000020030;
      if (Charge3He < 0) {
        pdg3He = -1000020030;
      }
      
      KFParticle kfpDaughter1 = CreateKFTrack(ESDtrackHe3, pdg3He);
      /// Fill 3He variables
      FillHe3Variables(ESDtrackHe3, kfpDaughter1);
      /// Fill daughter varibales
      FillDaughterVariables(kfpDaughter1,kfpDaughter2);
      //  if ( !DaughterSelection() ) continue; /// currently no selection
      
      /// create the kfp mother and histogram pt and mass
      KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
      
      Charge = kfpMother.GetQ();
      /// remove reconstructed candidates with wrong charge
      if ( abs(Charge)!= 1 ) continue;
      
      /// Selection on the hypertriton properties
      if (!HypertritonCandidateSelection(kfpMother)) continue;
      
      FillDistanceToSececondaryVertex(kfpDaughter1,kfpDaughter2,kfpMother);
      
      CandidateTree->Fill();
    }
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ThreeBodyDecay(const vector<Int_t>& DeuteronCandTrackId, const vector<Int_t>& ProtonCandTrackId, const vector<Int_t>& PionCandTrackId){
  
  /// next event if there are no 3He or pion candidates
  if ( (Int_t) DeuteronCandTrackId.size() <= 0) return;
  if ( (Int_t) ProtonCandTrackId.size() <= 0) return;
  if ( (Int_t) PionCandTrackId.size() <= 0) return;
  
  /// container for hypertriton candidates and daughters to combine with the pions
  HyperTritonCandidates.clear();
  DeuteronDaughter.clear();
  ProtonDaughter.clear();
  
  /// loop over proton candidates and deutreron candidates storing the combined mother and the tracks if mass is okay
  for (unsigned int jProton=0; jProton < (UInt_t) ProtonCandTrackId.size(); jProton++) {
    Int_t ProtonID = ProtonCandTrackId[jProton];
    AliESDtrack* ESDtrackProton = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ProtonCandTrackId[jProton]));
    
    Int_t ChargeProton = (Int_t) ESDtrackProton->Charge();
    /// create KFParticle for Proton
    Int_t pdgProton = 2212;
    if (ChargeProton < 0) {
      pdgProton = -2212;
    }
    
    KFParticle kfpDaughter2 = CreateKFTrack(ESDtrackProton, pdgProton);
    
    for (unsigned int lDeuteron=0; lDeuteron < (UInt_t) DeuteronCandTrackId.size(); lDeuteron++) {
      /// Do not fill tree variables in here. They cannot be linked to the correct branch
      Int_t DeuteronID = DeuteronCandTrackId[lDeuteron];
      /// Make sure not to use a single track twice
      if ( ProtonID==DeuteronID ) continue;
      /// Get deuteron track
      AliESDtrack* ESDtrackDeuteron = static_cast<AliESDtrack*>(fInputEvent->GetTrack(DeuteronID));
      Int_t ChargeDeuteron = (Int_t) ESDtrackDeuteron->Charge();
      /// implemented via split containers:  if ( ChargeProton*ChargeDeuteron <= 0) continue; /// select only same charge pairs
      /// create KFParticle for Proton
      Int_t pdgDeuteron = 1000010020;
      if (ChargeDeuteron < 0) {
        pdgDeuteron = -1000010020;
      }
      
      KFParticle kfpDaughter1 = CreateKFTrack(ESDtrackDeuteron, pdgDeuteron);
      
      /// preselection based on the opening angle between deuteron and proton
      Float_t OpeningAngle = kfpDaughter2.GetAngle(kfpDaughter1);
      if ( OpeningAngle > 0.5 ) continue;
      
      /// create the kfp mother from deuteron and proton
      Float_t ErrorMass = -999;
      Float_t Mass = -999;;
      KFParticle kfpMother_Deuteron_Proton(kfpDaughter1,kfpDaughter2);
      Bool_t MassCalculated = kfpMother_Deuteron_Proton.GetMass(Mass, ErrorMass);
      if (MassCalculated !=0) Mass = -999.;
      /// preselection based on the mass (to be implemented based on MC) !!!
      if (Mass < 2.75) continue; /// according to MC it is larger than 2.75 GeV/c^2
      if (Mass >= 2.9) continue; /// according to MC it is most of the times smaller than 2.9 GeV/c^2
      
      /// store mother and daughters to combine it later with the pions
      HyperTritonCandidates.push_back(kfpMother_Deuteron_Proton);
      DeuteronDaughter.push_back(DeuteronID);
      ProtonDaughter.push_back(ProtonID);
    }
  }
  
  /// No pion loop if no hypertriton candidate was found
  if ( (Int_t) HyperTritonCandidates.size() <= 0) return;
  if(kDoQA) NumberOfDeuteronProtonCandidates +=  HyperTritonCandidates.size();
  
  /// combine pion candidates with reconstructed deuteron-proton mothers
  Int_t PionID = -1;
  for(unsigned int iPion=0; iPion < (UInt_t) PionCandTrackId.size(); iPion++) {
    PionID = PionCandTrackId[iPion];
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(PionID));         /// get a track (type AliESDtrack) from the event
    
    Int_t ChargePionTrack = (Int_t) track->Charge();
    /// create KFParticle for pion
    Int_t pdgPion = 211;
    if (ChargePionTrack < 0) {
      pdgPion = -211;
    }
    KFParticle kfpDaughter3 = CreateKFTrack(track, pdgPion);
    
    /// Fill pion variables
    FillPionVariables(track, kfpDaughter3);
    
    for (unsigned int iHypertTritonCandidate=0; iHypertTritonCandidate < (UInt_t) HyperTritonCandidates.size(); iHypertTritonCandidate++) {
      
      /// Make sure not to use a single track twice
      if ( DeuteronDaughter[iHypertTritonCandidate]==PionID ) continue;
      if ( ProtonDaughter[iHypertTritonCandidate]==PionID ) continue;
      
      KFParticle kfpMother = (KFParticle) HyperTritonCandidates[iHypertTritonCandidate];
      
      Int_t ChargeCandidate = kfpMother.GetQ();
      /// implemented via split containers:  if ( ChargeCandidate*ChargePionTrack >= 0) continue; /// select only opposite charge pairs
      
      /// Fill variables from first step of reconstruction(deuteron + Proton)
      CosPointingAngle_Deuteron_Proton = CalculatePointingAngle(kfpMother,PrimVertex);
      Float_t ErrorMass = -999;
      Bool_t MassCalculated = kfpMother.GetMass(mass_Deuteron_Proton, ErrorMass);
      if (MassCalculated !=0) mass_Deuteron_Proton = -999.;
      
      /// Add pion
      kfpMother += kfpDaughter3;
      
      Charge = kfpMother.GetQ();
      /// remove LS information to reduce the size of the output
      if ( abs(Charge)!= 1 ) continue; /// should be ensured already by charge selection on daughter pairs
      
      /// Selection on the hypertriton properties
      if (!HypertritonCandidateSelection(kfpMother, false)) continue;
      
      AliESDtrack* ESDtrackDeuteron = static_cast<AliESDtrack*>(fInputEvent->GetTrack(DeuteronDaughter[iHypertTritonCandidate]));
      AliESDtrack* ESDtrackProton = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ProtonDaughter[iHypertTritonCandidate]));
      
      Int_t ChargeDeuteron = (Int_t) ESDtrackDeuteron->Charge();
      /// create KFParticle for Proton
      Int_t pdgDeuteron = 1000010020;
      if (ChargeDeuteron < 0) {
        pdgDeuteron = -1000010020;
      }
      
      Int_t ChargeProton = (Int_t) ESDtrackProton->Charge();
      Int_t pdgProton = 2212;
      if (ChargeProton < 0) {
        pdgProton = -2212;
      }
      
      KFParticle kfpDaughter1 = CreateKFTrack(ESDtrackDeuteron, pdgDeuteron);
      KFParticle kfpDaughter2 = CreateKFTrack(ESDtrackProton, pdgProton);
      
      /// Fill daughter varibales
      FillDaughterVariables(kfpDaughter1,kfpDaughter2, kfpDaughter3);
      
      if( OpeningAngle_Pion_Proton > 0.7 ) continue;
      if( OpeningAngle_Pion_Deuteron > 0.7 ) continue;
      
      /// Reject Lambda candidates
      /// create the kfp mother from proton and pion
      Float_t ErrorMass_Proton_Pion;
      KFParticle kfpMother_Proton_Pion(kfpDaughter2,kfpDaughter3);
      Bool_t MassCalculated2 = kfpMother_Proton_Pion.GetMass(mass_Proton_Pion, ErrorMass_Proton_Pion);
      //      if (MassCalculated2 == 0 && mass_Proton_Pion > 1.11 && mass_Proton_Pion < 1.12) continue; /// Mass window has to be adjusted, mass:Lambda = 1.116 GeV/c^2
      
      /// Fill deuteron variables
      FillDeuteronVariables(ESDtrackDeuteron, kfpDaughter1);
      /// Fill proton variables
      FillProtonVariables(ESDtrackProton, kfpDaughter2);
      
      FillDistanceToSececondaryVertex(kfpDaughter1, kfpDaughter2, kfpDaughter3, kfpMother);
      
      CandidateTree_3Body->Fill();
    }
  }
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ProcessMC()
{
  /// Monte Carlo only part of the task
  
  /// Protect against missing MC trees
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!mcH){
    AliError("No MC Event Handler available");
    return;
  }
  //  if(!mcH->InitOk()) return;
  //  if(!mcH->TreeK()) return;
  //  if(!mcH->TreeTR()) return;
  
  fMCEvent = mcH->MCEvent();
  
  if(!fMCEvent){
    AliError("No MC Event, but MC Data required");
    return;
  }
  
  ///----------------------------------------------------------------------------------------------
  /// create a translation table: ESDTrackId(mcTrackId) = ESDTrackIndex, or = -1 if there is no ESDTrack
  ///--------------------------------------------------------------------------------------------------
  vector<Int_t> ESDTrackId(fMCEvent->GetNumberOfTracks(), -1);
  
  /// loop over all ESD tracks and fill ESDTrackId
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(i)); /// get a track (type AliESDtrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;  /// skip track if it fails basic track selection
    
    ESDTrackId[abs(track->GetLabel())] = i;
  }
  
  /// loop over MC all particles
  for (Int_t igen = 0; igen < fMCEvent->GetNumberOfTracks(); igen++){
    AliVParticle *mcParticle = fMCEvent->GetTrack(igen);
    if(!mcParticle) continue;
    
    /// check if we have a hypertriton
    if( abs(mcParticle->PdgCode()) != 1010010030 ) continue;
    /// preselection on rapidity
    RapidityMC = mcParticle->Y();
    if (abs(RapidityMC) > 1.) continue;
    
    pMC = mcParticle->P();
    pTMC = mcParticle->Pt();
    
    //For MC candidate tree
    pxMC = mcParticle->Px();
    pyMC = mcParticle->Py();
    pzMC = mcParticle->Pz();
    
    Int_t Id3He = -1;
    Int_t IdDeuteron = -1;
    Int_t IdProton = -1;
    Int_t IdPion = -1;
    
    /// due to delta electrons a loop over daughters is needed
    AliVParticle* mcDaughter;
    Int_t pdgDaughter = -1;
    Int_t n3He=0;
    Int_t nPion=0;
    Int_t nDeutron=0;
    Int_t nProton=0;
    
    for (int idaughter = mcParticle->GetDaughterFirst(); idaughter <= mcParticle->GetDaughterLast(); idaughter++) {
      mcDaughter = fMCEvent->GetTrack(idaughter);
      if(!mcDaughter) continue;
      pdgDaughter = mcDaughter->PdgCode();
      
      if (abs(pdgDaughter) == 1000020030) { /// 3He
        n3He++;
        Id3He = idaughter;
      } else if (abs(pdgDaughter) == 1000010020){ /// deuteron
        nDeutron++;
        IdDeuteron = idaughter;
      } else if (abs(pdgDaughter) == 211){ /// pion
        nPion++;
        IdPion = idaughter;
      } else if (abs(pdgDaughter) == 2212){ /// proton
        nProton++;
        IdProton = idaughter;
      }
    }
    if (n3He > 1 || nDeutron > 1 || nProton > 1 || nPion > 1) continue;
    
    Bool_t IsTwoBody = false;
    Bool_t IsThreeBody = false;
    
    if (Id3He!=-1 && IdPion!=-1 && IdDeuteron==-1 && IdProton==-1) IsTwoBody = true;
    if (IdDeuteron!=-1 && IdProton!=-1 && IdPion!=-1 && Id3He==-1) IsThreeBody = true;
    
    if ( kRun2BodyDecay && IsTwoBody ) TwoBodyDecayMC(ESDTrackId, Id3He, IdPion);
    if ( kRun3BodyDecay && IsThreeBody ) ThreeBodyDecayMC(ESDTrackId, IdDeuteron, IdProton, IdPion);
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::TwoBodyDecayMC(const vector<Int_t>& ESDTrackId, Int_t Id3He, Int_t IdPion){
  
  AliVParticle* mcDaughter1 = fMCEvent->GetTrack(Id3He); /// 3He
  AliVParticle* mcDaughter2 = fMCEvent->GetTrack(IdPion); /// pion
  if (!mcDaughter1) return;
  if (!mcDaughter2) return;
  
  Int_t pdgDaughter1 = mcDaughter1->PdgCode();
  if (abs(pdgDaughter1)!= 1000020030) return;
  Int_t pdgDaughter2 = mcDaughter2->PdgCode();
  if (abs(pdgDaughter2)!= 211) return;
  
  /// fill variables for all generated MC hypertritons
  /// calculated decay length
  Double_t prodVertex[3];
  mcDaughter1->XvYvZv(prodVertex);
  
  xSecVertexMC = prodVertex[0];
  ySecVertexMC = prodVertex[1];
  zSecVertexMC = prodVertex[2];
  
  Double_t primVertex[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(primVertex);
  
  Double_t Difference[3];
  for (Int_t iCord=0; iCord < 3; iCord++) {
    Difference[iCord] = primVertex[iCord]-prodVertex[iCord];
  }
  
  DecayLengthMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]) + (Difference[2]*Difference[2])); /// approximation of the helix curve
  DecayLengthXYMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]));
  
  GeneratedTreeMC->Fill();
  
  /// look for the ESD trackId
  Int_t ESDId1 = ESDTrackId[Id3He];
  Int_t ESDId2 = ESDTrackId[IdPion];
  /// check that both daughters produced an ESD track (with basic track cut)
  if( ESDId1 == -1 || ESDId2 == -1 ) return;
  
  AliESDtrack* ESDTrack1 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId1));
  AliESDtrack* ESDTrack2 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId2));
  if (!ESDTrack1) return;
  if (!ESDTrack2) return;
  
  pxHeMC = mcDaughter1->Px();
  pyHeMC = mcDaughter1->Py();
  pzHeMC = mcDaughter1->Pz();
  
  pxPionMC =  mcDaughter2->Px();
  pyPionMC =  mcDaughter2->Py();
  pzPionMC =  mcDaughter2->Pz();
  
  if (kDoQA) {
    /// QA histograms
    /// Fill histograms for momentum range
    fHistMomHe3->Fill(mcDaughter1->P());
    
    fHistpxTrueRecHe3->Fill(2*ESDTrack1->Px(), pxHeMC-2*ESDTrack1->Px());
    fHistpyTrueRecHe3->Fill(2*ESDTrack1->Py(), pyHeMC-2*ESDTrack1->Py());
    fHistpzTrueRecHe3->Fill(2*ESDTrack1->Pz(), pzHeMC-2*ESDTrack1->Pz());
    
    fHistMomPion->Fill(mcDaughter2->P());
  }
  
  ///    Check identity of daughters
  if ( !Helium3Selection(ESDTrack1) ) return;
  if ( !PionSelection(ESDTrack2) ) return;
  
  KFParticle kfpDaughter1 = CreateKFTrack(ESDTrack1, pdgDaughter1);
  KFParticle kfpDaughter2 = CreateKFTrack(ESDTrack2, pdgDaughter2);
  
  FillHe3Variables(ESDTrack1, kfpDaughter1);
  FillPionVariables(ESDTrack2, kfpDaughter2);
  FillDaughterVariables(kfpDaughter1,kfpDaughter2);
  
  //  if ( !DaughterSelection() ) continue; /// currently no selection
  
  /// create the kfp mother and histogram pt and mass
  KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
  if (!HypertritonCandidateSelection(kfpMother)) return;
  
  FillDistanceToSececondaryVertex(kfpDaughter1,kfpDaughter2,kfpMother);
  
  Charge = kfpMother.GetQ();
  CandidateTree->Fill();
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ThreeBodyDecayMC(const vector<Int_t>& ESDTrackId, Int_t IdDeuteron, Int_t IdProton, Int_t IdPion){
  
  AliVParticle* mcDaughter1 = fMCEvent->GetTrack(IdDeuteron); /// Deuteron
  AliVParticle* mcDaughter2 = fMCEvent->GetTrack(IdProton); /// Proton
  AliVParticle* mcDaughter3 = fMCEvent->GetTrack(IdPion); /// Pion
  if (!mcDaughter1) return;
  if (!mcDaughter2) return;
  if (!mcDaughter3) return;
  
  Int_t pdgDaughter1 = mcDaughter1->PdgCode();
  if (abs(pdgDaughter1)!= 1000010020) return;
  Int_t pdgDaughter2 = mcDaughter2->PdgCode();
  if (abs(pdgDaughter2)!= 2212) return;
  Int_t pdgDaughter3 = mcDaughter3->PdgCode();
  if (abs(pdgDaughter3)!= 211) return;
  
  pxDeuteronMC = mcDaughter1->Px();
  pyDeuteronMC = mcDaughter1->Py();
  pzDeuteronMC = mcDaughter1->Pz();
  
  pxProtonMC =  mcDaughter2->Px();
  pyProtonMC =  mcDaughter2->Py();
  pzProtonMC =  mcDaughter2->Pz();
  
  pxPionMC =  mcDaughter3->Px();
  pyPionMC =  mcDaughter3->Py();
  pzPionMC =  mcDaughter3->Pz();
  
  /// fill variables for all generated MC hypertritons
  /// calculated decay length => not correct, different result depending on the daughter used
  /// best would be a decay channel independent method
  Double_t prodVertex[3];
  mcDaughter1->XvYvZv(prodVertex);
  
  xSecVertexMC = prodVertex[0];
  ySecVertexMC = prodVertex[1];
  zSecVertexMC = prodVertex[2];
  
  Double_t primVertex[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(primVertex);
  
  Double_t Difference[3];
  for (Int_t iCord=0; iCord < 3; iCord++) {
    Difference[iCord] = primVertex[iCord]-prodVertex[iCord];
  }
  DecayLengthMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]) + (Difference[2]*Difference[2])); /// approximation of the helix curve
  DecayLengthXYMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]));
  
  GeneratedTreeMC_3Body->Fill();
  
  /// look for the ESD trackId
  Int_t ESDId1 = ESDTrackId[IdDeuteron];
  Int_t ESDId2 = ESDTrackId[IdProton];
  Int_t ESDId3 = ESDTrackId[IdPion];
  /// check that daughters produced an ESD track (with basic track cut)
  if( ESDId1 == -1 || ESDId2 == -1 || ESDId3 == -1 ) return;
  
  AliESDtrack* ESDTrack1 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId1));
  AliESDtrack* ESDTrack2 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId2));
  AliESDtrack* ESDTrack3 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId3));
  if (!ESDTrack1) return;
  if (!ESDTrack2) return;
  if (!ESDTrack3) return;
  
  if (kDoQA) {
    /// QA histograms
    /// Fill histograms for momentum range
    fHistMomDeuteron->Fill(mcDaughter1->P());
    fHistMomProton->Fill(mcDaughter2->P());
    fHistMomPion_3Body->Fill(mcDaughter3->P());
  }
  
  /// Check identity of daughters
  if ( !DeuteronSelection(ESDTrack1) ) return;
  if ( !ProtonSelection(ESDTrack2) ) return;
  if ( !PionSelection(ESDTrack3) ) return;
  
  KFParticle kfpDaughter1 = CreateKFTrack(ESDTrack1, pdgDaughter1);
  KFParticle kfpDaughter2 = CreateKFTrack(ESDTrack2, pdgDaughter2);
  KFParticle kfpDaughter3 = CreateKFTrack(ESDTrack3, pdgDaughter3);
  
  FillDeuteronVariables(ESDTrack1, kfpDaughter1);
  FillProtonVariables(ESDTrack2, kfpDaughter2);
  FillPionVariables(ESDTrack3, kfpDaughter3);
  
  FillDaughterVariables(kfpDaughter1,kfpDaughter2, kfpDaughter3);
  
  /// create the kfp mother and histogram pt and mass
  KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
  /// evaluate mass of the pion and proton combination
  Float_t ErrorMass_Deuteron_Proton;
  Bool_t MassCalculated = kfpMother.GetMass(mass_Deuteron_Proton, ErrorMass_Deuteron_Proton);
  if (MassCalculated !=0) mass_Deuteron_Proton = -999.;
  
  CosPointingAngle_Deuteron_Proton = CalculatePointingAngle(kfpMother,PrimVertex);
  
  /// create the kfp mother from proton and pion to check if Lambda rejection is reasonable
  Float_t ErrorMass_Proton_Pion;
  KFParticle kfpMother_Proton_Pion(kfpDaughter2,kfpDaughter3);
  Bool_t MassCalculated2 = kfpMother_Proton_Pion.GetMass(mass_Proton_Pion, ErrorMass_Proton_Pion);
  //      if (MassCalculated2 == 0 && mass_Proton_Pion > 1.11 && mass_Proton_Pion < 1.12) continue; /// Mass window has to be adjusted, mass:Lambda = 1.116 GeV/c^2
  
  
  /// add deuteron
  kfpMother += kfpDaughter3;
  if (!HypertritonCandidateSelection(kfpMother, false)) return;
  
  Charge = kfpMother.GetQ();
  
  FillDistanceToSececondaryVertex(kfpDaughter1, kfpDaughter2, kfpDaughter3, kfpMother);
  
  CandidateTree_3Body->Fill();
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::PassedEventSelection() {
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return false;
  }
  hNumberOfEvents->Fill(0.5);
  
  /// Event selection using the AliEventCuts class
  if (!fEventCuts.AcceptEvent(fInputEvent)) {
    return false;
  }
  hNumberOfEvents->Fill(1.5);
  
  ///Vertex Contributors (olny for first checks
  if ( fInputEvent->GetPrimaryVertex()->GetNContributors() < 2 ) return false;
  hNumberOfEvents->Fill(2.5);
  
  CentralityPercentile = 300;
  AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
  if(!MultSelection) {
    /// If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    PostData(1,fQAList);
    return false;
  }else{
    CentralityPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  hNumberOfEvents->Fill(3.5);
  /// centrality selection 0--90% (maybe drop later)
  if (CentralityPercentile < 0.0) return false;
  if (CentralityPercentile >= 90.0) return false;
  hNumberOfEvents->Fill(4.5);
  
  ///Time-Range Selection (for LHC18r) (for nor standalone but could be integrated into AliEventCuts: fEventCuts.UseTimeRangeCut(); )
  fTimeRangeCut.InitFromEvent(InputEvent());
  const Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent());
  if (cutThisEvent) return false;
  hNumberOfEvents->Fill(5.5);
  
  /// Get PID response object
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) {
    AliError("No PID Response found");
    return false;
  }
  hNumberOfEvents->Fill(6.5);
  
  histoEventCentrality->Fill(CentralityPercentile);
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::PassedBasicTrackQualityCuts (AliESDtrack* track) {
  //Basic Track selection for the analysis
  
  Bool_t Passed = fESDtrackCuts->AcceptTrack(track);
  
  if (GetDCA(track,"3D") > 15.) Passed = false; /// upper limit for the DCA (?) to avoid KF floating point exception errors
  
  /// TPC nPIDcluster cut
  if ( (Int_t) track->GetTPCsignalN() < 50 ) Passed = false;
  
  return Passed;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::Helium3Selection (AliESDtrack* track)  {
  
  /// To reject He3 from spallation
  if ( track->P() < 1./2. ) return false; /// track information is rigidity p/q
  /// To ensure a clean selection
  if ( track->P() > 10. ) return false; /// track information is rigidity p/q
  
  /// PID selection
  Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
  if ( abs(dEdxSigmaHe3) > 4 ) return false;
  Double_t dEdxSigmaH3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
  if ( dEdxSigmaH3 < 2 ) return false;
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::PionSelection (AliESDtrack* track)  {
  
  /// pions from hypertriton are usally below p = 1.2 GeV/c
  if ( track->P() < 0.1 ) return false;
  if ( track->P() > 1.2 ) return false;
  
  /// Apply stronger preselection on number of crossed rows in the TPC for pions
  if ( (Int_t) track->GetTPCCrossedRows() < 70) return false;
  
  Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  /// Select pions with its specific energy loss in the TPC
  if (abs(dEdxSigmaPion) > 3.) return false;
  
  /// reject primary pions
  if (!kIsMC){ if (GetDCA(track,"3D") < 0.2 ) return false;}
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::DeuteronSelection (AliESDtrack* track)  {
  
  if ( track->P() < 1. ) return false; /// has to be adjusted (?)
  if ( track->P() > 10. ) return false; /// has to be adjusted (?)
  
  /// PID selection
  Double_t dEdxSigmaDeuteron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
  if ( abs(dEdxSigmaDeuteron) > 3. ) return false;
  
  if (track->P() > 1.5) {
    Double_t TOFSigmaDeuteron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
    if ( abs(TOFSigmaDeuteron) > 4. ) return false;
  }
  
  //  if (!kIsMC){ if (GetDCA(track,"3D") < 0.1 ) return false;}
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::ProtonSelection(AliESDtrack* track)  {
  
  /// proton from hypertriton are usally below p = 5.5 GeV/c
  if ( track->P() < 0.4 ) return false;
  if ( track->P() > 5.5 ) return false;
  
  /// PID selection
  Double_t dEdxSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  if ( abs(dEdxSigmaProton) > 3. ) return false;
  
  if (track->P() > 1.) {
    Double_t TOFSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    if ( abs(TOFSigmaProton) > 4. ) return false;
  }
  
  //  if (!kIsMC){ if (GetDCA(track,"3D") < 0.1 ) return false;}
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::DaughterSelection(){
  /// Put a preselection based on daughte rvariables here (2-body)
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::HypertritonCandidateSelection (KFParticle kfpMother, Bool_t TwoBody ){
  
  /// Select hypertriton candidates and fill variables
  if (kfpMother.GetNDF()<0) return false;
  if (kfpMother.GetChi2()<0) return false;
  if (kfpMother.GetChi2()>10000) return false; /// protection against infinite
  
  Bool_t MassCalculated = kfpMother.GetMass(mass, ErrorMass);
  if (MassCalculated !=0) return false;
  /// preselection
  if ( mass < 2.94 || mass > 3.05 ) return false;
  
  CosPointingAngle = CalculatePointingAngle(kfpMother,PrimVertex);
  /// preselection
  if (!kIsMC){
    /// different preselection for 2- & 3-body
    if (TwoBody && CosPointingAngle < 0.95) return false;
    if (!TwoBody && CosPointingAngle < 0.999) return false;
  }
  
  /// rapidity selection
  if ( ((kfpMother.E() - kfpMother.Pz()) > 0) && (kfpMother.E() + kfpMother.Pz()) >= 0 ) Rapidity = kfpMother.GetRapidity(); /// Missing protection for GetRapidity
  /// preselection
  if (!kIsMC){ if ( abs(Rapidity) > 1. ) return false;}
  Chi2PerNDF = kfpMother.GetChi2()/kfpMother.GetNDF();
  /// preselection
  if (!kIsMC){ if ( Chi2PerNDF > 20) return false;}
  
  /// fill variables
  DistanceToPV =  kfpMother.GetDistanceFromVertex(PrimVertex);
  DistanceToPVXY =  kfpMother.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPV = kfpMother.GetDeviationFromVertex(PrimVertex);
  DeviationFromPVXY = kfpMother.GetDeviationFromVertexXY(PrimVertex);
  
  ///preselection on DCA of hypertriton to the primary vertex
  if (DistanceToPV > 5) return false;
  
  px = kfpMother.GetPx();
  py = kfpMother.GetPy();
  pz = kfpMother.GetPz();
  
  if (kIsMC) {
    pxVariance = kfpMother.Covariance(9);
    pyVariance = kfpMother.Covariance(14);
    pzVariance = kfpMother.Covariance(20);
  }
  
  KFParticle kfpMotherTopo;
  kfpMotherTopo = kfpMother;
  
  kfpMotherTopo.SetProductionVertex(PrimVertex); /// if primary vertex is modified above use the modified version
  
  Bool_t MassCalculatedTopo = kfpMotherTopo.GetMass(massTopo, ErrorMassTopo);
  if (MassCalculatedTopo!=0) {
    massTopo = -999;
    ErrorMassTopo = -999;
  }
  
  Chi2PerNDFTopo = kfpMotherTopo.GetChi2()/kfpMotherTopo.GetNDF();
  
  CosPointingAngleTopo = CalculatePointingAngle(kfpMotherTopo,PrimVertex);
  
  DistanceToPVTopo =  kfpMotherTopo.GetDistanceFromVertex(PrimVertex);
  DistanceToPVXYTopo =  kfpMotherTopo.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVTopo = kfpMotherTopo.GetDeviationFromVertex(PrimVertex);
  DeviationFromPVXYTopo = kfpMotherTopo.GetDeviationFromVertexXY(PrimVertex);
  
  
  pxTopo = kfpMotherTopo.GetPx();
  pyTopo = kfpMotherTopo.GetPy();
  pzTopo = kfpMotherTopo.GetPz();
  if ( ((kfpMotherTopo.E() - kfpMotherTopo.Pz()) > 0) && (kfpMotherTopo.E() + kfpMotherTopo.Pz()) >= 0 ) RapidityTopo = kfpMotherTopo.GetRapidity();
  
  bool DecayLengthCalculated = kfpMotherTopo.GetDecayLength(DecayLength, ErrorDecayLength); /// returns 0 is correctly calculated
  if (DecayLengthCalculated != 0) {
    DecayLength = -999;
    ErrorDecayLength = -999;
  }
  bool DecayLengthXYCalculated = kfpMotherTopo.GetDecayLengthXY(DecayLengthXY, ErrorDecayLengthXY); /// returns 0 is correctly calculated
  if (DecayLengthXYCalculated != 0){
    DecayLengthXY = -999;
    ErrorDecayLengthXY = -999;
  }
  
  return true;
  
}
///____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTree::GetDCA (AliESDtrack *track , TString type)  {
  
  /// Calculate DCA to primary vertex
  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),fInputEvent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -9999;
  
  if(!strncmp(type,"xy",2)) {
    return impactParameter[0];
  } else if (!strncmp(type,"z",1)){
    return impactParameter[1];
  } else if (!strncmp(type,"3D",2)){
    return sqrt(impactParameter[0]*impactParameter[0] + impactParameter[1]*impactParameter[1]);
  }
  
  AliError("### Error ####: Wrong DCA type given \n");
  return -999;
}
///____________________________________________________________
KFParticle AliAnalysisTaskHypertritonKFTree::CreateKFTrack(AliESDtrack *track, int pdgCode){
  /// Input variables for KF should be float instead of double because not all functions are compatible with both
  
  /// GetTrack parameters
  Double_t trackParameter[6];
  Double_t covMatrix[21];
  
  track->GetXYZ(trackParameter);
  track->GetPxPyPz(&trackParameter[3]);
  track->GetCovarianceXYZPxPyPz(covMatrix);
  
  Int_t Charge = (Int_t) track->Charge();
  if (abs(pdgCode) == 1000020030){
    Charge = Charge*2; /// The exact value seems to be not relevant
    for (int i=3; i<6; i++) {
      trackParameter[i] = trackParameter[i]*2;
    }
    for (int i=6; i<21; i++) {
      covMatrix[i] = covMatrix[i]*2;  /// scale mom space entries of cov matrix by 2
      if (i==9 || i==13 || i==14 || i==18 || i==19 || i==20 ) {
        covMatrix[i] = covMatrix[i]*2;  /// scale mom mom entries of cov matrix by 4
      }
    }
  }
  
  /// Interface to KFParticle
  KFPTrack kfpTrk;
  /// Set the values
  kfpTrk.SetParameters((Float_t) trackParameter[0],(Float_t) trackParameter[1],(Float_t) trackParameter[2],(Float_t) trackParameter[3],(Float_t) trackParameter[4],(Float_t) trackParameter[5]);
  Float_t covF[21];
  for (Int_t i = 0; i<21;i++) { covF[i] = (Float_t) covMatrix[i]; }
  kfpTrk.SetCovarianceMatrix(covF);
  kfpTrk.SetCharge(Charge);
  kfpTrk.SetNDF(1); /// where is this coming from why 1  /// track should be 2?
  
  /// Get Chi2perNDF like for AOD
  Float_t TrackChi2perNDF = 999;
  Int_t  nClustersTPC = track->GetTPCNcls();
  if ( nClustersTPC > 5) {
    TrackChi2perNDF = (Float_t) track->GetTPCchi2()/Float_t(nClustersTPC - 5);
  }
  kfpTrk.SetChi2(TrackChi2perNDF);
  
  /// Build KFParticle
  KFParticle KFTrk(kfpTrk, pdgCode);
  
  return KFTrk;
}
///____________________________________________________________
KFVertex AliAnalysisTaskHypertritonKFTree::CreateKFVertex(const AliVVertex* vertex){
  
  /// GetTrack parameters
  Double_t param[6];
  Double_t cov[6];
  
  vertex->GetXYZ(param);
  vertex->GetCovarianceMatrix(cov);
  
  KFPVertex kfpVtx;
  /// Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
    (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}
///____________________________________________________________
Float_t AliAnalysisTaskHypertritonKFTree::CalculatePointingAngle(KFParticle KFPart, KFVertex KFVtx){
  
  KFPart.TransportToDecayVertex(); /// After SetProductionVertex the particle is stored at its production vertex but the information at the decay vertex is needed
  
  /// Store position of secondary vertex
  xSecVertex = KFPart.GetX();
  ySecVertex = KFPart.GetY();
  zSecVertex = KFPart.GetZ();
  
  xSecVertexVariance = KFPart.Covariance(0);
  ySecVertexVariance = KFPart.Covariance(2);
  zSecVertexVariance = KFPart.Covariance(5);
  
  Double_t v[3];
  v[0] = KFPart.GetX() - KFVtx.GetX();
  v[1] = KFPart.GetY() - KFVtx.GetY();
  v[2] = KFPart.GetZ() - KFVtx.GetZ();
  
  Double_t p[3];
  p[0] = KFPart.GetPx();
  p[1] = KFPart.GetPy();
  p[2] = KFPart.GetPz();
  
  Float_t vnorm3 = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  Float_t pnorm3 = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  
  Double_t pointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3); /// cos(pointing angle)
  
  return (Float_t) pointAngle;
  
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillHe3Variables (AliESDtrack* track, KFParticle KFPart){
  
  pxHe = track->Px();
  pyHe = track->Py();
  pzHe = track->Pz();
  
  ChargeHe = 2*track->Charge();
  
  TPCMom3He = track->GetTPCmomentum();
  TPCnSigma3He = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
  TPCnSigma3H = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
  DCA3He = GetDCA(track,"3D");
  //  DistanceToPV3He =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPV3He = KFPart.GetDeviationFromVertex(PrimVertex);
  DCA3HeXY = GetDCA(track,"xy");
  //  DistanceToPV3HeXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPV3HeXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0He3 = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1He3 = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2He3 = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3He3 = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4He3 = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5He3 = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPC3He = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPC3He = (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingHe3 = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxHeVariance = covMatrix[9];
    pyHeVariance = covMatrix[14];
    pzHeVariance = covMatrix[20];
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillPionVariables (AliESDtrack* track, KFParticle KFPart){
  
  pxPion = track->Px();
  pyPion = track->Py();
  pzPion = track->Pz();
  
  ChargePion = track->Charge();
  
  TPCMomPion = track->GetTPCmomentum();
  TPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  
  DCAPion = GetDCA(track,"3D");
  //  DistanceToPVPion =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPVPion = KFPart.GetDeviationFromVertex(PrimVertex);
  DCAPionXY = GetDCA(track,"xy");
  //  DistanceToPVPionXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVPionXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0Pion = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1Pion = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2Pion = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3Pion = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4Pion = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5Pion = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPCPion = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPCPion =  (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingPion = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxPionVariance = covMatrix[9];
    pyPionVariance = covMatrix[14];
    pzPionVariance = covMatrix[20];
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDeuteronVariables (AliESDtrack* track, KFParticle KFPart){
  
  pxDeuteron = track->Px();
  pyDeuteron = track->Py();
  pzDeuteron = track->Pz();
  
  ChargeDeuteron = track->Charge();
  
  TPCMomDeuteron = track->GetTPCmomentum();
  TPCnSigmaDeuteron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
  TOFnSigmaDeuteron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
  
  DCADeuteron = GetDCA(track,"3D");
  //  DistanceToPVDeuteron =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPVDeuteron = KFPart.GetDeviationFromVertex(PrimVertex);
  DCADeuteronXY = GetDCA(track,"xy");
  //  DistanceToPVDeuteronXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVDeuteronXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0Deuteron = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1Deuteron = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2Deuteron = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3Deuteron = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4Deuteron = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5Deuteron = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPCDeuteron = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPCDeuteron =  (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingDeuteron = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxDeuteronVariance = covMatrix[9];
    pyDeuteronVariance = covMatrix[14];
    pzDeuteronVariance = covMatrix[20];
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillProtonVariables (AliESDtrack* track, KFParticle KFPart){
  
  pxProton = track->Px();
  pyProton = track->Py();
  pzProton = track->Pz();
  
  ChargeProton = track->Charge();
  
  TPCMomProton = track->GetTPCmomentum();
  TPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  TOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  
  DCAProton = GetDCA(track,"3D");
  //  DistanceToPVProton =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPVProton = KFPart.GetDeviationFromVertex(PrimVertex);
  DCAProtonXY = GetDCA(track,"xy");
  //  DistanceToPVProtonXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVProtonXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0Proton = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1Proton = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2Proton = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3Proton = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4Proton = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5Proton = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPCProton = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPCProton =  (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingProton = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxProtonVariance = covMatrix[9];
    pyProtonVariance = covMatrix[14];
    pzProtonVariance = covMatrix[20];
  }
}

///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDaughterVariables (KFParticle kfpDaughter1, KFParticle kfpDaughter2){
  
  /// two body decay (3He, pion)
  DistanceOfDaughters = kfpDaughter1.GetDistanceFromParticle(kfpDaughter2);
  DeviationOfDaughters = kfpDaughter1.GetDeviationFromParticle(kfpDaughter2);
  
  DistanceOfDaughtersXY = kfpDaughter1.GetDistanceFromParticleXY(kfpDaughter2);
  DeviationOfDaughtersXY = kfpDaughter1.GetDeviationFromParticleXY(kfpDaughter2);
  
  //  Float_t DistanceOfDaughtersXY2 = kfpDaughter2.GetDistanceFromParticleXY(kfpDaughter1);
  //  if (DistanceOfDaughtersXY != DistanceOfDaughtersXY2) {
  //    if (abs(DistanceOfDaughtersXY2-DistanceOfDaughtersXY) > 0.01) {
  //      printf("########## %.2f vs %.2f; chi^2/ndf: %.2f \n", DistanceOfDaughtersXY, DistanceOfDaughtersXY2, Chi2PerNDF);
  //    }
  //  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDaughterVariables (KFParticle kfpDeuteron, KFParticle kfpProton, KFParticle kfpPion){
  /// three body decay
  /// expects:
  /// kfpDaughter1 == deuteron
  /// kfpDaughter2 == proton
  /// kfpDaughter3 == pion
  
  
  DistanceOfDaughters_Deuteron_Proton = kfpDeuteron.GetDistanceFromParticle(kfpProton);
  DeviationOfDaughters_Deuteron_Proton = kfpDeuteron.GetDeviationFromParticle(kfpProton);
  
  DistanceOfDaughtersXY_Deuteron_Proton = kfpDeuteron.GetDistanceFromParticleXY(kfpProton);
  DeviationOfDaughtersXY_Deuteron_Proton = kfpDeuteron.GetDeviationFromParticleXY(kfpProton);
  
  DistanceOfDaughters_Deuteron_Pion = kfpDeuteron.GetDistanceFromParticle(kfpPion);
  DeviationOfDaughters_Deuteron_Pion = kfpDeuteron.GetDeviationFromParticle(kfpPion);
  
  DistanceOfDaughtersXY_Deuteron_Pion = kfpDeuteron.GetDistanceFromParticleXY(kfpPion);
  DeviationOfDaughtersXY_Deuteron_Pion = kfpDeuteron.GetDeviationFromParticleXY(kfpPion);
  
  DistanceOfDaughters_Proton_Pion = kfpProton.GetDistanceFromParticle(kfpPion);
  DeviationOfDaughters_Proton_Pion = kfpProton.GetDeviationFromParticle(kfpPion);
  
  DistanceOfDaughtersXY_Proton_Pion = kfpProton.GetDistanceFromParticleXY(kfpPion);
  DeviationOfDaughtersXY_Proton_Pion = kfpProton.GetDeviationFromParticleXY(kfpPion);
  
  OpeningAngle_Pion_Proton = kfpPion.GetAngle(kfpProton);
  OpeningAngle_Pion_Deuteron = kfpPion.GetAngle(kfpDeuteron);
  OpeningAngle_Proton_Deuteron = kfpProton.GetAngle(kfpDeuteron);
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDistanceToSececondaryVertex(KFParticle kfpHelium, KFParticle kfpPion, KFParticle kfpMother){
  
  kfpMother.TransportToDecayVertex(); /// After SetProductionVertex the particle is stored at its production vertex but the information at the decay vertex is needed
  KFVertex SecondaryVertex(kfpMother);
  
  DistanceToSecVertHe = kfpHelium.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertHe = kfpHelium.GetDeviationFromVertex(SecondaryVertex);
  
  DistanceToSecVertPion = kfpPion.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertPion = kfpPion.GetDeviationFromVertex(SecondaryVertex);
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDistanceToSececondaryVertex(KFParticle kfpDeuteron, KFParticle kfpProton, KFParticle kfpPion, KFParticle kfpMother){
  
  kfpMother.TransportToDecayVertex(); /// After SetProductionVertex the particle is stored at its production vertex but the information at the decay vertex is needed
  KFVertex SecondaryVertex(kfpMother);
  
  DistanceToSecVertDeuteron = kfpDeuteron.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertDeuteron = kfpDeuteron.GetDeviationFromVertex(SecondaryVertex);
  
  DistanceToSecVertProton = kfpProton.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertProton = kfpProton.GetDeviationFromVertex(SecondaryVertex);
  
  DistanceToSecVertPion = kfpPion.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertPion = kfpPion.GetDeviationFromVertex(SecondaryVertex);
}
