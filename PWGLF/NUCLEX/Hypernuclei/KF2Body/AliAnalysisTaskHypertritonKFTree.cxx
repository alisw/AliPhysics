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

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TTree.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"

#include <vector>
#include<iostream>

#include "AliAnalysisTaskHypertritonKFTree.h"

class AliAnalysisTaskHypertritonKFTree;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskHypertritonKFTree) // classimp: necessary for root

AliAnalysisTaskHypertritonKFTree::AliAnalysisTaskHypertritonKFTree() : AliAnalysisTaskSE(), 
fEventCuts(),
fPIDResponse(nullptr),
PrimVertex(),
kRunAsData(false),
kDoQA(false),
kIsMC(false),
kDoMCQA(false),
fMCEvent(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
histoEventCentrality(nullptr),
fCandidateTree(nullptr),
CentralityPercentile(-1),
mass(-1),
ErrorMass(-1),
p(-1),
pT(-1),
Rapidity(-1),
massTopo(-1),
ErrorMassTopo(-1),
pTopo(-1),
pTTopo(-1),
RapidityTopo(-1),
SignOfPair(0),
CosPointingAngle(-1),
Chi2PerNDF(-1),
Chi2PerNDFTopo(-1),
//Chi2PerNDFMass(-1),
DecayLength(-1),
ErrorDecayLength(-1),
DecayLengthXY(-99),
ErrorDecayLengthXY(-99),
DistanceToPV(-1),
DeviationFromPV(-1),
DistanceToPVXY(-99),
DeviationFromPVXY(-99),
DistanceOfDaughters(-1),
DeviationOfDaughters(-1),
DistanceOfDaughtersXY(-99),
DeviationOfDaughtersXY(-99),
pPion(-1),
pTPion(-1),
DCAPion(-1),
DistanceToPVPion(-1),
DeviationFromPVPion(-1),
DCAPionXY(-99),
DistanceToPVPionXY(-99),
DeviationFromPVPionXY(-99),
NClusterTPCPion(-1),
NPIDClusterTPCPion (-1),
TPCMomPion(-1),
TPCnSigmaPion(-1),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
p3He(-1),
pT3He(-1),
DCA3He(-1),
DistanceToPV3He(-1),
DeviationFromPV3He(-1),
DCA3HeXY(-99),
DistanceToPV3HeXY(-99),
DeviationFromPV3HeXY(-99),
NClusterTPC3He(-1),
NPIDClusterTPC3He(-1),
TPCMom3He(-1),
TPCnSigma3He(-1),
TPCnSigma3H(-1),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
fGeneratedTreeMC(nullptr),
fCandidateTreeMC(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistPxTrueRecHe3(nullptr),
fHistPyTrueRecHe3(nullptr),
fHistPzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTree::AliAnalysisTaskHypertritonKFTree(const char* name) : AliAnalysisTaskSE(name),
fEventCuts(),
fPIDResponse(nullptr),
PrimVertex(),
kRunAsData(false),
kDoQA(false),
kIsMC(false),
kDoMCQA(false),
fMCEvent(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
histoEventCentrality(nullptr),
fCandidateTree(nullptr),
CentralityPercentile(-1),
mass(-1),
ErrorMass(-1),
p(-1),
pT(-1),
Rapidity(-1),
massTopo(-1),
ErrorMassTopo(-1),
pTopo(-1),
pTTopo(-1),
RapidityTopo(-1),
SignOfPair(0),
CosPointingAngle(-1),
Chi2PerNDF(-1),
Chi2PerNDFTopo(-1),
//Chi2PerNDFMass(-1),
DecayLength(-1),
ErrorDecayLength(-1),
DecayLengthXY(-99),
ErrorDecayLengthXY(-99),
DistanceToPV(-1),
DeviationFromPV(-1),
DistanceToPVXY(-99),
DeviationFromPVXY(-99),
DistanceOfDaughters(-1),
DeviationOfDaughters(-1),
DistanceOfDaughtersXY(-1),
DeviationOfDaughtersXY(-1),
pPion(-1),
pTPion(-1),
DCAPion(-1),
DistanceToPVPion(-1),
DeviationFromPVPion(-1),
DCAPionXY(-1),
DistanceToPVPionXY(-1),
DeviationFromPVPionXY(-1),
NClusterTPCPion(-1),
NPIDClusterTPCPion (-1),
TPCMomPion(-1),
TPCnSigmaPion(-1),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
p3He(-1),
pT3He(-1),
DCA3He(-1),
DistanceToPV3He(-1),
DeviationFromPV3He(-1),
DCA3HeXY(-1),
DistanceToPV3HeXY(-1),
DeviationFromPV3HeXY(-1),
NClusterTPC3He(-1),
NPIDClusterTPC3He(-1),
TPCMom3He(-1),
TPCnSigma3He(-1),
TPCnSigma3H(-1),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
fGeneratedTreeMC(nullptr),
fCandidateTreeMC(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistPxTrueRecHe3(nullptr),
fHistPyTrueRecHe3(nullptr),
fHistPzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr)
{
  // constructor
  // define the input of the analysis: in this case we take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it,
  // it does its work automatically
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTree::~AliAnalysisTaskHypertritonKFTree()
{
  // destructor
  // at the end of your task, it is deleted from memory by calling this function
  if(fQAList) delete fQAList;
  if(fOutputList) delete fOutputList;
  if(fCandidateTree) delete fCandidateTree;
  if(fCandidateTreeMC) delete fCandidateTreeMC;
  if(fGeneratedTreeMC) delete fGeneratedTreeMC;
}

//_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //
  
//  // Runs which are bad for TPC PID due to a problem in the TPC gain in one or two sectors towards the end of the run
//  fEventCuts.UseTimeRangeCut()
  
  fQAList = new TList();          // list which will contain all of QA histograms
  fQAList->SetOwner(kTRUE);
  
  fOutputList = new TList();          // list which will contain all of your histograms at the end of the analysis, the contents of this list are written to the output file
  fOutputList->SetOwner(kTRUE);
  
  histoEventCentrality = new TH1I("histoEventCentrality","Events vs centrality percentile (V0M)",100,0,100);
  histoEventCentrality->GetXaxis()->SetTitle("Centrality percentile");
  histoEventCentrality->Sumw2();
  fOutputList->Add(histoEventCentrality);
  
  // QA histograms
  fEventCuts.AddQAplotsToList(fQAList); /// Add event selection QA plots
  if (kRunAsData) {
    if (kDoQA) {
      fHistNsigmaTPCvsP3He = new TH2F("fHistNsigmaTPCvsP3He","#it{N}#sigma^{TPC}_{^{3}He} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{^{3}He}", 100,0,10,100,-5.,5.);
      fHistNsigmaTPCvsP3He->Sumw2();
      fQAList->Add(fHistNsigmaTPCvsP3He);
      
      fHistNsigmaTPCvsPPion = new TH2F("fHistNsigmaTPCvsPPion","#it{N}#sigma^{TPC}_{#pi} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{#pi}", 100,0,10,100,-5.,5.);
      fHistNsigmaTPCvsPPion->Sumw2();
      fQAList->Add(fHistNsigmaTPCvsPPion);
    }
    
    //  Tree with hypertriton candidates
    fCandidateTree = new TTree("fCandidateTree","fCandidateTree");
    fCandidateTree->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    fCandidateTree->Branch("mass",&mass,"mass/F");
    fCandidateTree->Branch("ErrorMass",&ErrorMass,"ErrorMass/F");
    fCandidateTree->Branch("p",&p,"p/F");
    fCandidateTree->Branch("pT",&pT,"pT/F");
    fCandidateTree->Branch("Rapidity",&Rapidity,"Rapidity/F");
    fCandidateTree->Branch("massTopo",&massTopo,"massTopo/F");
    fCandidateTree->Branch("ErrorMassTopo",&ErrorMassTopo,"ErrorMassTopo/F");
    fCandidateTree->Branch("pTopo",&pTopo,"pTopo/F");
    fCandidateTree->Branch("pTTopo",&pTTopo,"pTTopo/F");
    fCandidateTree->Branch("RapidityTopo",&RapidityTopo,"RapidityTopo/F");
    fCandidateTree->Branch("SignOfPair",&SignOfPair,"SignOfPair/I");
    fCandidateTree->Branch("CosPointingAngle",&CosPointingAngle,"CosPointingAngle/F");
    fCandidateTree->Branch("Chi2PerNDF",&Chi2PerNDF,"Chi2PerNDF/F");
    fCandidateTree->Branch("Chi2PerNDFTopo",&Chi2PerNDFTopo,"Chi2PerNDFTopo/F");
//    fCandidateTree->Branch("Chi2PerNDFMass",&Chi2PerNDFMass,"Chi2PerNDFMass/F");
    fCandidateTree->Branch("DecayLength",&DecayLength,"DecayLength/F");
    fCandidateTree->Branch("ErrorDecayLength",&ErrorDecayLength,"ErrorDecayLength/F");
    fCandidateTree->Branch("DecayLengthXY",&DecayLengthXY,"DecayLengthXY/F");
    fCandidateTree->Branch("ErrorDecayLengthXY",&ErrorDecayLengthXY,"ErrorDecayLengthXY/F");
    fCandidateTree->Branch("DistanceToPV",&DistanceToPV,"DistanceToPV/F");
    fCandidateTree->Branch("DistanceToPVXY",&DistanceToPVXY,"DistanceToPVXY/F");
    fCandidateTree->Branch("DeviationFromPV",&DeviationFromPV,"DeviationFromPV/F");
    fCandidateTree->Branch("DeviationFromPVXY",&DeviationFromPVXY,"DeviationFromPVXY/F");
    
    // Daughter variables
    fCandidateTree->Branch("DistanceOfDaughters",&DistanceOfDaughters,"DistanceOfDaughters/F");
    fCandidateTree->Branch("DeviationOfDaughters",&DeviationOfDaughters,"DeviationOfDaughters/F");
    fCandidateTree->Branch("DistanceOfDaughtersXY",&DistanceOfDaughtersXY,"DistanceOfDaughtersXY/F");
    fCandidateTree->Branch("DeviationOfDaughtersXY",&DeviationOfDaughtersXY,"DeviationOfDaughtersXY/F");
    
    fCandidateTree->Branch("pPion",&pPion,"pPion/F");
    fCandidateTree->Branch("pTPion",&pTPion,"pTPion/F");
    fCandidateTree->Branch("DCAPion",&DCAPion,"DCAPion/F");
    fCandidateTree->Branch("DistanceToPVPion",&DistanceToPVPion,"DistanceToPVPion/F");
    fCandidateTree->Branch("DeviationFromPVPion",&DeviationFromPVPion,"DeviationFromPVPion/F");
    fCandidateTree->Branch("DCAPionXY",&DCAPionXY,"DCAPionXY/F");
    fCandidateTree->Branch("DistanceToPVPionXY",&DistanceToPVPionXY,"DistanceToPVPionXY/F");
    fCandidateTree->Branch("DeviationFromPVPionXY",&DeviationFromPVPionXY,"DeviationFromPVPionXY/F");
    fCandidateTree->Branch("NClusterTPCPion",&NClusterTPCPion,"NClusterTPCPion/I");
    fCandidateTree->Branch("NPIDClusterTPCPion",&NPIDClusterTPCPion,"NPIDClusterTPCPion/I");
    fCandidateTree->Branch("TPCMomPion",&TPCMomPion,"TPCMomPion/F");
    fCandidateTree->Branch("TPCnSigmaPion",&TPCnSigmaPion,"TPCnSigmaPion/F");
    fCandidateTree->Branch("HasPointOnITSLayer0Pion",&HasPointOnITSLayer0Pion,"HasPointOnITSLayer0Pion/O");
    fCandidateTree->Branch("HasPointOnITSLayer1Pion",&HasPointOnITSLayer1Pion,"HasPointOnITSLayer1Pion/O");
    fCandidateTree->Branch("HasPointOnITSLayer2Pion",&HasPointOnITSLayer2Pion,"HasPointOnITSLayer2Pion/O");
    fCandidateTree->Branch("HasPointOnITSLayer3Pion",&HasPointOnITSLayer3Pion,"HasPointOnITSLayer3Pion/O");
    fCandidateTree->Branch("HasPointOnITSLayer4Pion",&HasPointOnITSLayer4Pion,"HasPointOnITSLayer4Pion/O");
    fCandidateTree->Branch("HasPointOnITSLayer5Pion",&HasPointOnITSLayer5Pion,"HasPointOnITSLayer5Pion/O");
    
    fCandidateTree->Branch("p3He",&p3He,"p3He/F");
    fCandidateTree->Branch("pT3He",&pT3He,"pT3He/F");
    fCandidateTree->Branch("DCA3He",&DCA3He,"DCA3He/F");
    fCandidateTree->Branch("DistanceToPV3He",&DistanceToPV3He,"DistanceToPV3He/F");
    fCandidateTree->Branch("DeviationFromPV3He",&DeviationFromPV3He,"DeviationFromPV3He/F");
    fCandidateTree->Branch("DCA3HeXY",&DCA3HeXY,"DCA3HeXY/F");
    fCandidateTree->Branch("DistanceToPV3HeXY",&DistanceToPV3HeXY,"DistanceToPV3HeXY/F");
    fCandidateTree->Branch("DeviationFromPV3HeXY",&DeviationFromPV3HeXY,"DeviationFromPV3HeXY/F");
    fCandidateTree->Branch("NClusterTPC3He",&NClusterTPC3He,"NClusterTPC3He/I");
    fCandidateTree->Branch("NPIDClusterTPC3He",&NPIDClusterTPC3He,"NPIDClusterTPC3He/I");
    fCandidateTree->Branch("TPCMom3He",&TPCMom3He,"TPCMom3He/F");
    fCandidateTree->Branch("TPCnSigma3He",&TPCnSigma3He,"TPCnSigma3He/F");
    fCandidateTree->Branch("TPCnSigma3H",&TPCnSigma3H,"TPCnSigma3H/F");
    fCandidateTree->Branch("HasPointOnITSLayer0He3",&HasPointOnITSLayer0He3,"HasPointOnITSLayer0He3/O");
    fCandidateTree->Branch("HasPointOnITSLayer1He3",&HasPointOnITSLayer1He3,"HasPointOnITSLayer1He3/O");
    fCandidateTree->Branch("HasPointOnITSLayer2He3",&HasPointOnITSLayer2He3,"HasPointOnITSLayer2He3/O");
    fCandidateTree->Branch("HasPointOnITSLayer3He3",&HasPointOnITSLayer3He3,"HasPointOnITSLayer3He3/O");
    fCandidateTree->Branch("HasPointOnITSLayer4He3",&HasPointOnITSLayer4He3,"HasPointOnITSLayer4He3/O");
    fCandidateTree->Branch("HasPointOnITSLayer5He3",&HasPointOnITSLayer5He3,"HasPointOnITSLayer5He3/O");
    
  }
  
  if (kIsMC){
    // MC output
    // generated MC
    
    fGeneratedTreeMC = new TTree("fGeneratedTreeMC","fGeneratedTreeMC");
    fGeneratedTreeMC->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    fGeneratedTreeMC->Branch("p",&p,"p/F");
    fGeneratedTreeMC->Branch("pT",&pT,"pT/F");
    fGeneratedTreeMC->Branch("Rapidity",&Rapidity,"Rapidity/F");
    fGeneratedTreeMC->Branch("DecayLength",&DecayLength,"DecayLength/F");
    
    //  Tree with hypertriton candidates
    fCandidateTreeMC = new TTree("fCandidateTreeMC","fCandidateTreeMC");
    fCandidateTreeMC->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    fCandidateTreeMC->Branch("mass",&mass,"mass/F");
    fCandidateTreeMC->Branch("ErrorMass",&ErrorMass,"ErrorMass/F");
    fCandidateTreeMC->Branch("p",&p,"p/F");
    fCandidateTreeMC->Branch("pT",&pT,"pT/F");
    fCandidateTreeMC->Branch("Rapidity",&Rapidity,"Rapidity/F");
    fCandidateTreeMC->Branch("massTopo",&massTopo,"massTopo/F");
    fCandidateTreeMC->Branch("ErrorMassTopo",&ErrorMassTopo,"ErrorMassTopo/F");
    fCandidateTreeMC->Branch("pTopo",&pTopo,"pTopo/F");
    fCandidateTreeMC->Branch("pTTopo",&pTTopo,"pTTopo/F");
    fCandidateTreeMC->Branch("RapidityTopo",&RapidityTopo,"RapidityTopo/F");
    fCandidateTreeMC->Branch("SignOfPair",&SignOfPair,"SignOfPair/I");
    fCandidateTreeMC->Branch("CosPointingAngle",&CosPointingAngle,"CosPointingAngle/F");
    fCandidateTreeMC->Branch("Chi2PerNDF",&Chi2PerNDF,"Chi2PerNDF/F");
    fCandidateTreeMC->Branch("Chi2PerNDFTopo",&Chi2PerNDFTopo,"Chi2PerNDFTopo/F");
//    fCandidateTreeMC->Branch("Chi2PerNDFMass",&Chi2PerNDFMass,"Chi2PerNDFMass/F");
    fCandidateTreeMC->Branch("DecayLength",&DecayLength,"DecayLength/F");
    fCandidateTreeMC->Branch("ErrorDecayLength",&ErrorDecayLength,"ErrorDecayLength/F");
    fCandidateTreeMC->Branch("DecayLengthXY",&DecayLengthXY,"DecayLengthXY/F");
    fCandidateTreeMC->Branch("ErrorDecayLengthXY",&ErrorDecayLengthXY,"ErrorDecayLengthXY/F");
    fCandidateTreeMC->Branch("DistanceToPV",&DistanceToPV,"DistanceToPV/F");
    fCandidateTreeMC->Branch("DistanceToPVXY",&DistanceToPVXY,"DistanceToPVXY/F");
    fCandidateTreeMC->Branch("DeviationFromPV",&DeviationFromPV,"DeviationFromPV/F");
    fCandidateTreeMC->Branch("DeviationFromPVXY",&DeviationFromPVXY,"DeviationFromPVXY/F");

    // Daughter variables
    fCandidateTreeMC->Branch("DistanceOfDaughters",&DistanceOfDaughters,"DistanceOfDaughters/F");
    fCandidateTreeMC->Branch("DeviationOfDaughters",&DeviationOfDaughters,"DeviationOfDaughters/F");
    fCandidateTreeMC->Branch("DistanceOfDaughtersXY",&DistanceOfDaughtersXY,"DistanceOfDaughtersXY/F");
    fCandidateTreeMC->Branch("DeviationOfDaughtersXY",&DeviationOfDaughtersXY,"DeviationOfDaughtersXY/F");
    
    fCandidateTreeMC->Branch("pPion",&pPion,"pPion/F");
    fCandidateTreeMC->Branch("pTPion",&pTPion,"pTPion/F");
    fCandidateTreeMC->Branch("DCAPion",&DCAPion,"DCAPion/F");
    fCandidateTreeMC->Branch("DistanceToPVPion",&DistanceToPVPion,"DistanceToPVPion/F");
    fCandidateTreeMC->Branch("DeviationFromPVPion",&DeviationFromPVPion,"DeviationFromPVPion/F");
    fCandidateTreeMC->Branch("DCAPionXY",&DCAPionXY,"DCAPionXY/F");
    fCandidateTreeMC->Branch("DistanceToPVPionXY",&DistanceToPVPionXY,"DistanceToPVPionXY/F");
    fCandidateTreeMC->Branch("DeviationFromPVPionXY",&DeviationFromPVPionXY,"DeviationFromPVPionXY/F");
    fCandidateTreeMC->Branch("NClusterTPCPion",&NClusterTPCPion,"NClusterTPCPion/I");
    fCandidateTreeMC->Branch("NPIDClusterTPCPion",&NPIDClusterTPCPion,"NPIDClusterTPCPion/I");
    fCandidateTreeMC->Branch("TPCMomPion",&TPCMomPion,"TPCMomPion/F");
    fCandidateTreeMC->Branch("TPCnSigmaPion",&TPCnSigmaPion,"TPCnSigmaPion/F");
    fCandidateTreeMC->Branch("HasPointOnITSLayer0Pion",&HasPointOnITSLayer0Pion,"HasPointOnITSLayer0Pion/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer1Pion",&HasPointOnITSLayer1Pion,"HasPointOnITSLayer1Pion/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer2Pion",&HasPointOnITSLayer2Pion,"HasPointOnITSLayer2Pion/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer3Pion",&HasPointOnITSLayer3Pion,"HasPointOnITSLayer3Pion/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer4Pion",&HasPointOnITSLayer4Pion,"HasPointOnITSLayer4Pion/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer5Pion",&HasPointOnITSLayer5Pion,"HasPointOnITSLayer5Pion/O");
    
    fCandidateTreeMC->Branch("p3He",&p3He,"p3He/F");
    fCandidateTreeMC->Branch("pT3He",&pT3He,"pT3He/F");
    fCandidateTreeMC->Branch("DCA3He",&DCA3He,"DCA3He/F");
    fCandidateTreeMC->Branch("DistanceToPV3He",&DistanceToPV3He,"DistanceToPV3He/F");
    fCandidateTreeMC->Branch("DeviationFromPV3He",&DeviationFromPV3He,"DeviationFromPV3He/F");
    fCandidateTreeMC->Branch("DCA3HeXY",&DCA3HeXY,"DCA3HeXY/F");
    fCandidateTreeMC->Branch("DistanceToPV3HeXY",&DistanceToPV3HeXY,"DistanceToPV3HeXY/F");
    fCandidateTreeMC->Branch("DeviationFromPV3HeXY",&DeviationFromPV3HeXY,"DeviationFromPV3HeXY/F");
    fCandidateTreeMC->Branch("NClusterTPC3He",&NClusterTPC3He,"NClusterTPC3He/I");
    fCandidateTreeMC->Branch("NPIDClusterTPC3He",&NPIDClusterTPC3He,"NPIDClusterTPC3He/I");
    fCandidateTreeMC->Branch("TPCMom3He",&TPCMom3He,"TPCMom3He/F");
    fCandidateTreeMC->Branch("TPCnSigma3He",&TPCnSigma3He,"TPCnSigma3He/F");
    fCandidateTreeMC->Branch("TPCnSigma3H",&TPCnSigma3H,"TPCnSigma3H/F");
    fCandidateTreeMC->Branch("HasPointOnITSLayer0He3",&HasPointOnITSLayer0He3,"HasPointOnITSLayer0He3/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer1He3",&HasPointOnITSLayer1He3,"HasPointOnITSLayer1He3/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer2He3",&HasPointOnITSLayer2He3,"HasPointOnITSLayer2He3/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer3He3",&HasPointOnITSLayer3He3,"HasPointOnITSLayer3He3/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer4He3",&HasPointOnITSLayer4He3,"HasPointOnITSLayer4He3/O");
    fCandidateTreeMC->Branch("HasPointOnITSLayer5He3",&HasPointOnITSLayer5He3,"HasPointOnITSLayer5He3/O");
    
    if (kDoMCQA) {
      // QA and corrections
      fHistPxTrueRecHe3 = new TH2F("fHistPxTrueRecHe3", "(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) vs #it{p}_{x}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
      fHistPxTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{x}^{rec} (GeV/#it{c})");
      fHistPxTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) (GeV/#it{c})");
      fHistPxTrueRecHe3->Sumw2();
      fQAList->Add(fHistPxTrueRecHe3);
      
      fHistPyTrueRecHe3 = new TH2F("fHistPyTrueRecHe3", "(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) vs #it{p}_{y}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
      fHistPyTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{y}^{rec} (GeV/#it{c})");
      fHistPyTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) (GeV/#it{c})");
      fHistPyTrueRecHe3->Sumw2();
      fQAList->Add(fHistPyTrueRecHe3);
      
      fHistPzTrueRecHe3 = new TH2F("fHistPzTrueRecHe3", "(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) vs #it{p}_{z}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
      fHistPzTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{z}^{rec} (GeV/#it{c})");
      fHistPzTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) (GeV/#it{c})");
      fHistPzTrueRecHe3->Sumw2();
      fQAList->Add(fHistPzTrueRecHe3);
      
      fHistMomPion = new TH1F("fHistMomPion","Momentum distribution of #pi from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
      fHistMomPion->Sumw2();
      fQAList->Add(fHistMomPion);
      fHistMomHe3 = new TH1F("fHistMomHe3","Momentum distribution of ^{3}He from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
      fHistMomHe3->Sumw2();
      fQAList->Add(fHistMomHe3);
    }
  }
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  if (kRunAsData) PostData(3,fCandidateTree);
  if (kIsMC) PostData(4,fCandidateTreeMC);
  if (kIsMC) PostData(5,fGeneratedTreeMC);
}
//_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you
  // have access to the current event.
  // once you return from the UserExec function, the manager will retrieve the next event from the chain
  
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }
  
  // Event selection using the AliEventCuts class
  if (!fEventCuts.AcceptEvent(fInputEvent)) {
    PostData(1,fQAList);
    return;
  }
  
  CentralityPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection");
  if(!MultSelection) {
     //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
     AliWarning("AliMultSelection object not found!");
    return;
  }else{
     CentralityPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  
  // Get PID response object
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) {
    AliError("No PID Response found");
    return;
  }
  
  histoEventCentrality->Fill(CentralityPercentile);
  
  // set Magnetic field for KF
  KFParticle::SetField(fInputEvent->GetMagneticField());
  
  // Create Primary Vertex with KF
  PrimVertex = CreateKFVertex(fInputEvent->GetPrimaryVertex());
  
  if (kRunAsData) ProcessAOD();
  if (kIsMC) ProcessMC();
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  if (kRunAsData) PostData(3,fCandidateTree);
  if (kIsMC) PostData(4,fCandidateTreeMC);
  if (kIsMC) PostData(5,fGeneratedTreeMC);
}
//_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ProcessAOD()
{
  // container for found 3He candidates
  vector<Int_t> He3CandTrackId;
  
  // check for helium-3 candidates and store their track number
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;              // skip track is it fails basic track selection
    
    if (kDoQA) {
      // Fill QA histogram
      Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
      if(dEdxSigmaHe3 < 5) fHistNsigmaTPCvsP3He->Fill(track->P(), dEdxSigmaHe3);
    }
    
    if(!Helium3Selection(track) ) continue;

    // create KFParticle for pion
    Int_t pdg3He = 1000020030;
    if ((Int_t) track->Charge() < 0) {
      pdg3He = -1000020030;
    }

    KFParticle KFPart = CreateKFTrack(track, pdg3He);
    // pre-selection
    DeviationFromPV3He = KFPart.GetDeviationFromVertex(PrimVertex);
    if( DeviationFromPV3He < 3 )  continue;
    
    He3CandTrackId.push_back(i);
  }
  // next event if there are no 3He candidates
  if ( (Int_t) He3CandTrackId.size() <= 0) return;
  
  // combine pion candidates with 3He candidates
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue; // skip track is it fails basic track selection
    
    if (kDoQA) {
      Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
      if(dEdxSigmaPion < 5) fHistNsigmaTPCvsPPion->Fill(track->P(),dEdxSigmaPion);
    }
    
    if(!PionSelection(track)) continue;
    
    Int_t ChargePion = (Int_t) track->Charge();
    // create KFParticle for pion
    Int_t pdgPion = 211;
    if (ChargePion > 0) {
      pdgPion = -211;
    }
    KFParticle kfpDaughter1 = CreateKFTrack(track, pdgPion);
    
    // pre-selection
    DistanceToPVPion =  kfpDaughter1.GetDistanceFromVertex(PrimVertex);
    if( DistanceToPVPion < 0.2 )  continue;
    
    for (unsigned int jHe3=0; jHe3 < He3CandTrackId.size(); jHe3++) {
      AliAODTrack* AODtrackHe3 = static_cast<AliAODTrack*>(fInputEvent->GetTrack(He3CandTrackId[jHe3]));
      
      Int_t Charge3He = (Int_t) AODtrackHe3->Charge();
      // create KFParticle for pion
      Int_t pdg3He = 1000020030;
      if (Charge3He < 0) {
        pdg3He = -1000020030;
      }
      
      KFParticle kfpDaughter2 = CreateKFTrack(AODtrackHe3, pdg3He);

//      if ( !DaughterSelection(kfpDaughter1,kfpDaughter2) ) continue; // Currently no selection
      // create the kfp mother and histogram pt and mass
      KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
      // Selection on the hypertriton properties
      if (!HypertritonCandidateSelection(kfpMother, kfpDaughter1, kfpDaughter2)) continue;
      
     // Calculate variables only for selected candidates
      FillPionVariables(track, kfpDaughter1);
      FillHe3Variables(AODtrackHe3, kfpDaughter2);
      FillDaughterVariables(kfpDaughter1,kfpDaughter2);
      
      SignOfPair = -1;
      if (ChargePion*Charge3He > 0) {
        SignOfPair = 1;
      }
      
      fCandidateTree->Fill();
      
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTree::ProcessMC()
{
  // Monte Carlo only part of the task
  
  // get the corresponding MC event fMCEvent
  fMCEvent = MCEvent();
  if(!fMCEvent) return;
  //----------------------------------------------------------------------------------------------
  // create a translation table: aodTrackId(mcTrackId) = aodTrackIndex, or = -1 if there is no aodTrack
  //--------------------------------------------------------------------------------------------------
  vector<Int_t> aodTrackId(fMCEvent->GetNumberOfTracks(), -1);
  
  // loop over all aod tracks and fill aodTrackId
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;              // skip track if it fails basic track selection
    aodTrackId[abs(track->GetLabel())]=i;
  }
  
  // next event if there are no good tracks
  if ( (Int_t) aodTrackId.size() <= 0) return;

  
  TClonesArray* mcParticles = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if ( mcParticles == NULL ) return;
  
  Double_t primVertex[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(primVertex);
  
  // loop over all particles
  for (Int_t i=0; i < mcParticles->GetEntriesFast(); i++) {
    AliAODMCParticle* mcParticle = static_cast<AliAODMCParticle*>(mcParticles->At(i));
    if(!mcParticle) continue;
    // check if we have a hyper-triton with two daughters
    if( (abs(mcParticle->GetPdgCode()) == 1010010030)  && (mcParticle->GetNDaughters() == 2) ) {
      
      Int_t iDaughter1 = mcParticle->GetDaughterFirst();
      Int_t iDaughter2 = mcParticle->GetDaughterLast();
      
      AliAODMCParticle* mcDaughter1 = static_cast<AliAODMCParticle*>(mcParticles->At(iDaughter1));
      AliAODMCParticle* mcDaughter2 = static_cast<AliAODMCParticle*>(mcParticles->At(iDaughter2));
      if (!mcDaughter1) continue;
      if (!mcDaughter2) continue;
      
      Int_t pdgDaughter1 = mcDaughter1->GetPdgCode();
      Int_t pdgDaughter2 = mcDaughter2->GetPdgCode();
      
      // check if one of the daugthers is a charged pion
      if( abs(pdgDaughter1) == 1000020030 || abs(pdgDaughter2) == 1000020030 ) {
        // rapidity cut
        Rapidity = mcParticle->Y();
        if (abs(Rapidity) > 0.5) continue;
        // fill variables for all generated MC hypertritons
        //calculated decay length
        Double_t prodVertex[3];
        mcDaughter1->XvYvZv(prodVertex);
        Double_t Difference[3];
        for (Int_t iCord=0; iCord < 3; iCord++) {
          Difference[iCord] = primVertex[iCord]-prodVertex[iCord];
        }

        DecayLength = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]) + (Difference[2]*Difference[2]));
        p = mcParticle->P();
        pT = mcParticle->Pt();
        
        fGeneratedTreeMC->Fill();
        
        // look for the aod trackId
        Int_t aodId1 = aodTrackId[iDaughter1];
        Int_t aodId2 = aodTrackId[iDaughter2];
        // check that both daughters produced an AOD track (with basic track cut)
        if( aodId1 == -1 || aodId2 == -1 ) continue;
        
        AliAODTrack* aodTrack1 = static_cast<AliAODTrack*>(fInputEvent->GetTrack(aodId1));
        AliAODTrack* aodTrack2 = static_cast<AliAODTrack*>(fInputEvent->GetTrack(aodId2));
        if (!aodTrack1) continue;
        if (!aodTrack2) continue;
        
        if (kDoMCQA) {
          // QA histograms
          // Fill histograms for momentum range
          if ( abs(pdgDaughter1) == 1000020030) {
            fHistMomHe3->Fill(mcDaughter1->P());
          }
          if ( abs(pdgDaughter2) == 1000020030) {
            fHistMomHe3->Fill(mcDaughter2->P());
          }
          
          if (abs(pdgDaughter1) == 211) {
            fHistMomPion->Fill(mcDaughter1->P());
          }
          if (abs(pdgDaughter2) == 211) {
            fHistMomPion->Fill(mcDaughter2->P());
          }
        }
        
        //    Check identity of daughters
        if (!(Helium3Selection(aodTrack1) || PionSelection(aodTrack1))) continue;
        if (!(Helium3Selection(aodTrack2) || PionSelection(aodTrack2))) continue;
        
        // select true He3 and pions
        if ( !((abs(pdgDaughter1)==1000020030 && abs(pdgDaughter2)==211) || (abs(pdgDaughter1)==211 && abs(pdgDaughter2)==1000020030)) ) continue;
        
        if (kDoMCQA) {
          // QA histograms
          // Fill histograms for momentum correction
          if (abs(pdgDaughter1) == 1000020030) {
            fHistPxTrueRecHe3->Fill(2*aodTrack1->Px(), 2*aodTrack1->Px()-mcDaughter1->Px());
            fHistPyTrueRecHe3->Fill(2*aodTrack1->Py(), 2*aodTrack1->Py()-mcDaughter1->Py());
            fHistPzTrueRecHe3->Fill(2*aodTrack1->Pz(), 2*aodTrack1->Pz()-mcDaughter1->Pz());
          }
          if (abs(pdgDaughter2) == 1000020030) {
            fHistPxTrueRecHe3->Fill(2*aodTrack2->Px(), 2*aodTrack2->Px()-mcDaughter2->Px());
            fHistPyTrueRecHe3->Fill(2*aodTrack2->Py(), 2*aodTrack2->Py()-mcDaughter2->Py());
            fHistPzTrueRecHe3->Fill(2*aodTrack2->Pz(), 2*aodTrack2->Pz()-mcDaughter2->Pz());
          }
        }
        
        KFParticle kfpDaughter1 = CreateKFTrack(aodTrack1, mcDaughter1->GetPdgCode());
        KFParticle kfpDaughter2 = CreateKFTrack(aodTrack2, mcDaughter2->GetPdgCode());
        
        if ( abs(pdgDaughter1) == 1000020030) {
          FillHe3Variables(aodTrack1, kfpDaughter1);
        } else {
          FillPionVariables(aodTrack1, kfpDaughter1);
        }
        
        if ( abs(pdgDaughter2) == 1000020030) {
          FillHe3Variables(aodTrack2, kfpDaughter2);
        } else {
          FillPionVariables(aodTrack2, kfpDaughter2);
        }
        
        // pre-selection
        if( DistanceToPVPion < 0.2 )  continue;
        if( DeviationFromPV3He < 3 )  continue;

//        if ( !DaughterSelection(kfpDaughter1, kfpDaughter2) ) continue; // currently no selection
        
        // create the kfp mother and histogram pt and mass
        KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
        if (!HypertritonCandidateSelection(kfpMother, kfpDaughter1, kfpDaughter2)) continue;
        
        FillDaughterVariables(kfpDaughter1,kfpDaughter2);
        
        SignOfPair = -1;
        if (aodTrack1->Charge()*aodTrack2->Charge() > 0) {
          SignOfPair = 1;
        }
                
        fCandidateTreeMC->Fill();
      }
    }
  }
}
//____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::PassedBasicTrackQualityCuts (AliAODTrack* track)  {
  //Basic Track selection for the analysis (filterbit, pseudorapidity range, ITS and TPC clusters)
  
  //Filterbit
  if(!(track->TestFilterMask(AliAODTrack::kTrkTPCOnly))) return false; // AliAODTrack::kTrkGlobalNoDCA // AliAODTrack::kTrkTPCOnly
  
  //  additional track selection
  if(abs(track->Eta()) > 0.9) return false; // Geometrical acceptance
//  if(track->GetTPCNcls() < 70) return false; // Additional constraint on the number of clusters in the TPC; 25.11.2019 (not used by default anymore -> Issue wit MC matching)
  //  if(track->GetITSNcls() < 2) return false; // Additional constraint on the number of hits in the ITS
  
  return true;
}
//____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::Helium3Selection (AliAODTrack* track)  {
  
  // To reject He3 from spallation
  if ( track->P() < 1./2. ) return false; // track information is rigidity p/q
  // To ensure a clean selection (needs to be adjusted)
  if ( track->P() > 8. ) return false; // track information is rigidity p/q

  // store TPC nsigma vs momentum to adjust selection
  Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
  if ( abs(dEdxSigmaHe3) > 5 ) return false;
  
//    Double_t dEdxSigmaH3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
//    if ( abs(dEdxSigmaH3) < 0 ) return false;
    
  return true;
}
//____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::PionSelection (AliAODTrack* track)  {
  
  //  // pions from hypertriton are usally below p = 2 GeV/c
  if ( track->P() > 1.2 ) return false;
  if ( track->P() < 0.1 ) return false;

  Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  // Select pions with its specific energy loss in the TPC
  if (abs(dEdxSigmaPion) > 4) return false;
  
  // reject primary pions
  //  if (GetDCA(track,"3D") < 0.1 ) return false;
  
  return true;
}
//____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::DaughterSelection (KFParticle kfpDaughter1, KFParticle kfpDaughter2){
    
  return true;
}
//____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTree::HypertritonCandidateSelection (KFParticle kfpMother, KFParticle kfpDaughter1, KFParticle kfpDaughter2){
  
  // Select hypertriton candidates and fill variables
  
  if (kfpMother.GetNDF()<0) return false;
  if (kfpMother.GetChi2()<0) return false;
  if (kfpMother.GetChi2()>10000) return false; // protection against infinite
  // rapidity selection
  if ( (kfpMother.E() - kfpMother.Pz()) <= 0 || (kfpMother.E() + kfpMother.Pz()) < 0 ) return false; // Missing protection for GetRapidity
  Rapidity = kfpMother.GetRapidity();
  Bool_t MassCalculated = kfpMother.GetMass(mass, ErrorMass);
  if (MassCalculated !=0) return false;
//
  // pre-selection
  if( mass < 2.94 || mass > 3.2 ) return false;
  if ( abs(Rapidity) > 0.5 ) return false;
  Chi2PerNDF = kfpMother.GetChi2()/kfpMother.GetNDF();
  if ( Chi2PerNDF > 100) return false;
  CosPointingAngle = CalculatePointingAngle(kfpMother,PrimVertex);
  if ( CosPointingAngle < 0.9) return false;
    
  KFParticle kfpMotherTopo;
  kfpMotherTopo = kfpMother;
    
//  KFVertex PrimVertexModified = CreateKFVertex(fInputEvent->GetPrimaryVertex());
//  PrimVertexModified -= kfpDaughter1;
//  PrimVertexModified -= kfpDaughter2;
//  PrimVertexModified += kfpMotherTopo; // Causes the decay length to be zero (?)
  
  kfpMotherTopo.SetProductionVertex(PrimVertex); // if primary vertex is modified above use the modified version
  Chi2PerNDFTopo = kfpMotherTopo.GetChi2()/kfpMotherTopo.GetNDF();
  if ( Chi2PerNDFTopo < 0) return false; // use for analysis as a function of ct or which use topological constraint
  
  // fill variables
  DistanceToPV =  kfpMother.GetDistanceFromVertex(PrimVertex);
  DistanceToPVXY =  kfpMother.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPV = kfpMother.GetDeviationFromVertex(PrimVertex);
  DeviationFromPVXY = kfpMother.GetDeviationFromVertexXY(PrimVertex);
  
  p = kfpMother.GetP();
  pT = kfpMother.GetPt();

  Bool_t MassCalculatedTopo = kfpMotherTopo.GetMass(massTopo, ErrorMassTopo);
  if (MassCalculatedTopo!=0) {
    massTopo = -1;
    ErrorMassTopo = -1;
  }
  
  pTopo = kfpMotherTopo.GetP();
  pTTopo = kfpMotherTopo.GetPt();
  if ( (kfpMotherTopo.E() - kfpMotherTopo.Pz()) <= 0 || (kfpMotherTopo.E() + kfpMotherTopo.Pz()) < 0 ) return false; // Missing protection for GetRapidity
  RapidityTopo = kfpMotherTopo.GetRapidity();
  
  bool DecayLengthCalculated = kfpMotherTopo.GetDecayLength(DecayLength, ErrorDecayLength); // returns 0 is correctly calculated
  if (DecayLengthCalculated != 0) {
    DecayLength = -1;
    ErrorDecayLength = -1;
  }
  bool DecayLengthXYCalculated = kfpMotherTopo.GetDecayLengthXY(DecayLengthXY, ErrorDecayLengthXY); // returns 0 is correctly calculated
  if (DecayLengthXYCalculated != 0){
    DecayLengthXY = -99;
    ErrorDecayLengthXY = -99;
  }

//  KFParticle kfpMotherMass;
//  kfpMotherMass = kfpMother;
//  // Nonlinear mass constraint with floating point exception protection
//  const float px = kfpMotherMass.GetPx();
//  const float py = kfpMotherMass.GetPy();
//  const float pz = kfpMotherMass.GetPz();
//  const float energy  = kfpMotherMass.GetE();
//  const float residual = (energy*energy - px*px - py*py - pz*pz) - mass*mass;
//
//  const float dm2 = float(4.f) * ( px*px*kfpMotherMass.GetCovariance(9) + py*py*kfpMotherMass.GetCovariance(14) + pz*pz*kfpMotherMass.GetCovariance(20) + energy*energy*kfpMotherMass.GetCovariance(27) + float(2.f) * ( px*py*kfpMotherMass.GetCovariance(13) + pz*(px*kfpMotherMass.GetCovariance(18)+py*kfpMotherMass.GetCovariance(19)) - energy*(px*kfpMotherMass.GetCovariance(24)+py*kfpMotherMass.GetCovariance(25)+pz*kfpMotherMass.GetCovariance(26)) ) ); // version matching V1.1
//
//  if (dm2 != 0) {
//    kfpMotherMass.SetNonlinearMassConstraint(2.99131);
//    Chi2PerNDFMass = kfpMotherMass.GetChi2()/kfpMotherMass.GetNDF();
//  } else{
//    Chi2PerNDFMass = -1;
//  }
    
  return true;
  
}
//____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTree::GetDCA (AliAODTrack *track , TString type)  {
  
  /// Calculate DCA to primary vertex
  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),fInputEvent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
  
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
//____________________________________________________________
KFParticle AliAnalysisTaskHypertritonKFTree::CreateKFTrack(AliAODTrack *track, int pdgCode){
  // Input variables for KF should be float instead of double because not all functions are compatible with both
  
  // GetTrack parameters
  Double_t trackParameter[6];
  Double_t covMatrix[21];
  
  track->GetXYZ(trackParameter);
  track->GetPxPyPz(&trackParameter[3]);
  track->GetCovarianceXYZPxPyPz(covMatrix);
  
  Int_t Charge = (Int_t) track->Charge();
  if (abs(pdgCode) == 1000020030){
    
    Charge = Charge*2; // The exact value seems to be not relevant
    
    for (int i=3; i<6; i++) {
      trackParameter[i] = trackParameter[i]*2;
    }
    
    for (int i=6; i<21; i++) {
      covMatrix[i] = covMatrix[i]*2;  // scale mom space entries of cov matrix by 2
      
      if (i==9 || i==13 || i==14 || i==18 || i==19 || i==20 ) {
        covMatrix[i] = covMatrix[i]*2;  // scale mom mom entries of cov matrix by 4
      }
    }
  }
  
  // Interface to KFParticle
  KFPTrack kfpTrk;
  // Set the values
  kfpTrk.SetParameters((Float_t) trackParameter[0],(Float_t) trackParameter[1],(Float_t) trackParameter[2],(Float_t) trackParameter[3],(Float_t) trackParameter[4],(Float_t) trackParameter[5]);
  Float_t covF[21];
  for (Int_t i = 0; i<21;i++) { covF[i] = (Float_t) covMatrix[i]; }
  kfpTrk.SetCovarianceMatrix(covF);
  kfpTrk.SetCharge(Charge);
  kfpTrk.SetNDF(1); // where is this coming from why 1  // track should be 2?
  kfpTrk.SetChi2((Float_t) track->Chi2perNDF()); // have to be scaled by ndf?
  
  // Build KFParticle
  KFParticle KFTrk(kfpTrk, pdgCode);
  return KFTrk;
}
//____________________________________________________________
KFVertex AliAnalysisTaskHypertritonKFTree::CreateKFVertex(const AliVVertex* vertex){
  
  // GetTrack parameters
  Double_t param[6];
  Double_t cov[21];
  
  vertex->GetXYZ(param);
  vertex->GetCovarianceMatrix(cov);
  
  KFPVertex kfpVtx;
  // Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
    (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}
//____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTree::CalculatePointingAngle(KFParticle KFPart, KFVertex KFVtx){
  
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
  
  Double_t pointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3); // cos(pointing angle)
  
  return pointAngle;
  
}
//____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillHe3Variables (AliAODTrack* track, KFParticle KFPart){
  
  p3He = track->P();
  pT3He = track->Pt();
  TPCMom3He = track->GetTPCmomentum();
  TPCnSigma3He = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
  TPCnSigma3H = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
  DCA3He = GetDCA(track,"3D");
  DistanceToPV3He =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPV3He = KFPart.GetDeviationFromVertex(PrimVertex);
  DCA3HeXY = GetDCA(track,"xy");
  DistanceToPV3HeXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPV3HeXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0He3 = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1He3 = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2He3 = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3He3 = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4He3 = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5He3 = track->HasPointOnITSLayer(5);
  
  NClusterTPC3He = (Int_t) track->GetTPCNcls();
  NPIDClusterTPC3He = (Int_t) track->GetTPCsignalN();

}

//____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillPionVariables (AliAODTrack* track, KFParticle KFPart){
  
  pPion = track->P();
  pTPion = track->Pt();
  TPCMomPion = track->GetTPCmomentum();
  TPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  
  DCAPion = GetDCA(track,"3D");
  DistanceToPVPion =  KFPart.GetDistanceFromVertex(PrimVertex);
  DeviationFromPVPion = KFPart.GetDeviationFromVertex(PrimVertex);
  DCAPionXY = GetDCA(track,"xy");
  DistanceToPVPionXY =  KFPart.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVPionXY = KFPart.GetDeviationFromVertexXY(PrimVertex);
  
  HasPointOnITSLayer0Pion = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1Pion = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2Pion = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3Pion = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4Pion = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5Pion = track->HasPointOnITSLayer(5);
  
  NClusterTPCPion = (Int_t) track->GetTPCNcls();
  NPIDClusterTPCPion =  (Int_t) track->GetTPCsignalN();

}
//____________________________________________________________
void AliAnalysisTaskHypertritonKFTree::FillDaughterVariables (KFParticle kfpDaughter1, KFParticle kfpDaughter2){
  
  DistanceOfDaughters = kfpDaughter1.GetDistanceFromParticle(kfpDaughter2);
  DeviationOfDaughters = kfpDaughter1.GetDeviationFromParticle(kfpDaughter2);
  
  DistanceOfDaughtersXY = kfpDaughter1.GetDistanceFromParticleXY(kfpDaughter2);
  DeviationOfDaughtersXY = kfpDaughter1.GetDeviationFromParticleXY(kfpDaughter2);
}
