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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to acquire MC-level predictions in general
// for several LF-related particle species. First deployed to deal with
// the Pb-Pb 5 TeV strangeness analysis. Adapted to several other use cases
// afterwards, including prompt/non-prompt HF
//
// Please report any bugs, complaints, suggestions to:
// --- david.dobrigkeit.chinellato@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskMCPredictions.h"
////
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions)

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions()
: AliAnalysisTaskSE(),
fListHist(0),
fHistEventCounter(0),
fHistChargedEta(0),
fSmallMultRange(1000),
fLargeMultRange(2000),
fRebinFactor(1),
fkNBBins(1),
fkNNpartBins(1),
fkNEtaBins(800),
fkSelectINELgtZERO(kTRUE),
fkALICE3SiliconMode(kTRUE),
fkWideRapiditySpeciesStudy(kFALSE),
fkDoImpactParameterStudy(kFALSE),
fkDoNpartStudy(kFALSE),
fkDoNMPIStudy(kTRUE),
fkDoRapidityStudy(kFALSE),
fkMinimumMultiplicity(-1),
fHistV0MMult(0),
fHistSPDMult(0),
fHistNchVsV0MMult(0),
fHistNchVsSPDMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0),
fHistNMPI(0),
fHistNchVsNMPI(0),
fkDo2pc(kTRUE),
fMinPtTrigger(2.0),
fMaxPtTrigger(4.0),
fHistPtTriggerD0(0),
fHistPtTriggerXiC(0),
fHistPtTriggerXiB(0),
fHist3d2pcD0Proton(0),
fHist3d2pcD0AntiProton(0),
fHist3d2pcD0D0(0),
fHist3d2pcD0D0bar(0),
fHist3d2pcD0KMinus(0),
fHist3d2pcD0KPlus(0),
fHist3d2pcXiCProton(0),
fHist3d2pcXiCAntiProton(0),
fHist3d2pcXiCD0(0),
fHist3d2pcXiCD0bar(0),
fHist3d2pcXiCKMinus(0),
fHist3d2pcXiCKPlus(0),
fHist3d2pcXiBProton(0),
fHist3d2pcXiBAntiProton(0),
fHist3d2pcXiBBMinus(0),
fHist3d2pcXiBBPlus(0),
fHist3d2pcXiBKMinus(0),
fHist3d2pcXiBKPlus(0),
fEMBufferFullD0(kFALSE),
fEMBufferCycleD0(0),
fEMBufferFullXiC(kFALSE),
fEMBufferCycleXiC(0),
fEMBufferFullXiB(kFALSE),
fEMBufferCycleXiB(0),
fHistMixed3d2pcD0Proton(0),
fHistMixed3d2pcD0AntiProton(0),
fHistMixed3d2pcD0D0(0),
fHistMixed3d2pcD0D0bar(0),
fHistMixed3d2pcD0KMinus(0),
fHistMixed3d2pcD0KPlus(0),
fHistMixed3d2pcXiCProton(0),
fHistMixed3d2pcXiCAntiProton(0),
fHistMixed3d2pcXiCD0(0),
fHistMixed3d2pcXiCD0bar(0),
fHistMixed3d2pcXiCKMinus(0),
fHistMixed3d2pcXiCKPlus(0),
fHistMixed3d2pcXiBProton(0),
fHistMixed3d2pcXiBAntiProton(0),
fHistMixed3d2pcXiBBMinus(0),
fHistMixed3d2pcXiBBPlus(0),
fHistMixed3d2pcXiBKMinus(0),
fHistMixed3d2pcXiBKPlus(0)
{
  for(Int_t ii=0; ii<10; ii++){
    fEMBufferEtaD0[ii]=0;
    fEMBufferPhiD0[ii]=0;
    fEMBufferEtaXiC[ii]=0;
    fEMBufferPhiXiC[ii]=0;
    fEMBufferEtaXiB[ii]=0;
    fEMBufferPhiXiB[ii]=0;
  }
  for(Int_t ih=0; ih<72; ih++){
    fHistPt[ih]          = 0x0;
    fHistEta[ih]         = 0x0;
    fHistEtaTriggeredMeson[ih]= 0x0;
    fHistEtaTriggeredCharm[ih]= 0x0;
    fHistEtaTriggeredBeauty[ih]= 0x0;
    fHistPtVsV0MMult[ih] = 0x0;
    fHistPtVsSPDMult[ih] = 0x0;
    fHistEtaVsSPDMult[ih]= 0x0;
    fHistYVsSPDMult[ih]  = 0x0;
    fHistPtVsNpart[ih]   = 0x0;
    fHistPtVsB[ih]       = 0x0;
    fHistPtVsNMPI[ih]   = 0x0;
    //    fHist3d2pcSE[ih]     = 0x0;
    //    fHist3d2pcXiSE[ih]   = 0x0;
    //    fHist3d2pcPhiSE[ih]  = 0x0;
  }
}

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions(const char *name, Int_t lNSmallBinning, Int_t lNLargeBinning, Int_t lRebinFactor, Int_t lNBBins, Int_t lNNpartBins, Int_t lNEtaBins)
: AliAnalysisTaskSE(name),
fListHist(0),
fHistEventCounter(0),
fHistChargedEta(0),
fSmallMultRange(lNSmallBinning),
fLargeMultRange(lNLargeBinning),
fRebinFactor(lRebinFactor),
fkNBBins(lNBBins),
fkNNpartBins(lNNpartBins),
fkNEtaBins(lNEtaBins),
fkSelectINELgtZERO(kTRUE),
fkALICE3SiliconMode(kTRUE),
fkWideRapiditySpeciesStudy(kFALSE),
fkDoImpactParameterStudy(kFALSE),
fkDoNpartStudy(kFALSE),
fkDoNMPIStudy(kTRUE),
fkDoRapidityStudy(kFALSE),
fkMinimumMultiplicity(-1),
fHistV0MMult(0),
fHistSPDMult(0),
fHistNchVsV0MMult(0),
fHistNchVsSPDMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0),
fHistNMPI(0),
fHistNchVsNMPI(0),
fkDo2pc(kTRUE),
fMinPtTrigger(2.0),
fMaxPtTrigger(4.0),
fHistPtTriggerD0(0),
fHistPtTriggerXiC(0),
fHistPtTriggerXiB(0),
fHist3d2pcD0Proton(0),
fHist3d2pcD0AntiProton(0),
fHist3d2pcD0D0(0),
fHist3d2pcD0D0bar(0),
fHist3d2pcD0KMinus(0),
fHist3d2pcD0KPlus(0),
fHist3d2pcXiCProton(0),
fHist3d2pcXiCAntiProton(0),
fHist3d2pcXiCD0(0),
fHist3d2pcXiCD0bar(0),
fHist3d2pcXiCKMinus(0),
fHist3d2pcXiCKPlus(0),
fHist3d2pcXiBProton(0),
fHist3d2pcXiBAntiProton(0),
fHist3d2pcXiBBMinus(0),
fHist3d2pcXiBBPlus(0),
fHist3d2pcXiBKMinus(0),
fHist3d2pcXiBKPlus(0),
fEMBufferFullXiC(kFALSE),
fEMBufferCycleXiC(0),
fEMBufferFullXiB(kFALSE),
fEMBufferCycleXiB(0),
fHistMixed3d2pcD0Proton(0),
fHistMixed3d2pcD0AntiProton(0),
fHistMixed3d2pcD0D0(0),
fHistMixed3d2pcD0D0bar(0),
fHistMixed3d2pcD0KMinus(0),
fHistMixed3d2pcD0KPlus(0),
fHistMixed3d2pcXiCProton(0),
fHistMixed3d2pcXiCAntiProton(0),
fHistMixed3d2pcXiCD0(0),
fHistMixed3d2pcXiCD0bar(0),
fHistMixed3d2pcXiCKMinus(0),
fHistMixed3d2pcXiCKPlus(0),
fHistMixed3d2pcXiBProton(0),
fHistMixed3d2pcXiBAntiProton(0),
fHistMixed3d2pcXiBBMinus(0),
fHistMixed3d2pcXiBBPlus(0),
fHistMixed3d2pcXiBKMinus(0),
fHistMixed3d2pcXiBKPlus(0)
{
  for(Int_t ii=0; ii<10; ii++){
    fEMBufferEtaD0[ii]=0;
    fEMBufferPhiD0[ii]=0;
    fEMBufferEtaXiC[ii]=0;
    fEMBufferPhiXiC[ii]=0;
    fEMBufferEtaXiB[ii]=0;
    fEMBufferPhiXiB[ii]=0;
  }
  for(Int_t ih=0; ih<72; ih++){
    fHistPt[ih]          = 0x0;
    fHistEta[ih]         = 0x0;
    fHistEtaTriggeredMeson[ih]= 0x0;
    fHistEtaTriggeredCharm[ih]= 0x0;
    fHistEtaTriggeredBeauty[ih]= 0x0;
    fHistPtVsV0MMult[ih] = 0x0;
    fHistPtVsSPDMult[ih] = 0x0;
    fHistEtaVsSPDMult[ih]= 0x0;
    fHistYVsSPDMult[ih]  = 0x0;
    fHistPtVsNpart[ih]   = 0x0;
    fHistPtVsB[ih]       = 0x0;
    fHistPtVsNMPI[ih]   = 0x0;
    //    fHist3d2pcSE[ih]     = 0x0;
    //    fHist3d2pcXiSE[ih]   = 0x0;
    //    fHist3d2pcPhiSE[ih]  = 0x0;
  }
  DefineOutput(1, TList::Class()); // Event Counter Histo
}


AliAnalysisTaskMCPredictions::~AliAnalysisTaskMCPredictions()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
  
  if (fListHist) {
    delete fListHist;
    fListHist = 0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions::UserCreateOutputObjects()
{
  //------------------------------------------------
  // Histograms: Basic Analysis Output
  //------------------------------------------------
  // Create histograms
  fListHist = new TList();
  fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
  
  //Settings for transverse momentum
  Int_t lNPtBins = 200;
  Double_t lMaxPt = 20.0;
  
  Int_t lNEtaBins = fkNEtaBins;
  Double_t lMaxAbsEta = 4;
  
  //Settings for charged particle counters (integers!)
  Int_t lNNchBins = fSmallMultRange/fRebinFactor;
  Double_t lLowNchBound  = -0.5;
  Double_t lHighNchBound = -0.5 + ((double)(fSmallMultRange));
  
  Int_t lNNchBinsV0M = fLargeMultRange/fRebinFactor;
  Double_t lLowNchBoundV0M  = -0.5;
  Double_t lHighNchBoundV0M = -0.5 + ((double)(fLargeMultRange));
  
  if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",1,0,1);
    //Keeps track of some basics
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fListHist->Add(fHistEventCounter);
  }
  if(! fHistChargedEta && fkDoRapidityStudy ) {
    //Histogram Output: Event-by-Event
    fHistChargedEta = new TH1D( "fHistChargedEta", ";#eta;Count",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
    fListHist->Add(fHistChargedEta);
  }
  //___________________________________________________
  if(! fHistV0MMult ) {
    //Histogram Output: Event-by-Event
    fHistV0MMult = new TH1D( "fHistV0MMult", ";V0M Mult;Count",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistV0MMult);
  }
  if(! fHistSPDMult ) {
    //Histogram Output: Event-by-Event
    fHistSPDMult = new TH1D( "fHistSPDMult", ";SPD Mult;Count",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistSPDMult);
  }
  if(! fHistNchVsV0MMult ) {
    //Histogram Output: Event-by-Event
    fHistNchVsV0MMult = new TH2D( "fHistNchVsV0MMult", ";V0M Mult;Count",
                                 lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,
                                 lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsV0MMult);
  }
  if(! fHistNchVsSPDMult ) {
    //Histogram Output: Event-by-Event
    fHistNchVsSPDMult = new TH2D( "fHistNchVsSPDMult", ";SPD Mult;Count",
                                 lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,
                                 lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsSPDMult);
  }
  //___________________________________________________
  if(! fHistNpart && fkDoNpartStudy) {
    //Histogram Output: Event-by-Event
    fHistNpart = new TH1D( "fHistNpart", ";N_{part};Count",fkNNpartBins,-0.5,((double)(fkNNpartBins))-0.5);
    //Keeps track of some basics
    fListHist->Add(fHistNpart);
  }
  if(! fHistNchVsNpart && fkDoNpartStudy ) {
    //Histogram Output: Event-by-Event
    fHistNchVsNpart = new TH2D( "fHistNchVsNpart", ";N_{part};Count",fkNNpartBins,-0.5,((double)(fkNNpartBins))-0.5,lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsNpart);
  }
  //___________________________________________________
  if(! fHistB && fkDoImpactParameterStudy ) {
    //Histogram Output: Event-by-Event
    fHistB = new TH1D( "fHistB", ";b;Count",fkNBBins,0,20);
    //Keeps track of some basics
    fListHist->Add(fHistB);
  }
  if(! fHistNchVsB && fkDoImpactParameterStudy ) {
    //Histogram Output: Event-by-Event
    fHistNchVsB = new TH2D( "fHistNchVsB", ";b;Count",fkNBBins,0,20,lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsB);
  }
  //___________________________________________________
  if(! fHistNMPI && fkDoNMPIStudy ) {
    //Histogram Output: Event-by-Event
    fHistNMPI = new TH1D( "fHistNMPI", ";N_{MPI};Count",50,-0.5,49.5);
    //Keeps track of some basics
    fListHist->Add(fHistNMPI);
  }
  if(! fHistNchVsNMPI && fkDoNMPIStudy ) {
    //Histogram Output: Event-by-Event
    fHistNchVsNMPI = new TH2D( "fHistNchVsNMPI", ";N_{part};Count",50,-0.5,49.5,lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsNMPI);
  }
  //___________________________________________________
  
  //Identified Particles
  TString lPartNames[72] = {
    "PiPlus", "PiMinus", "KPlus", "KMinus", "Proton", "AntiProton",
    "K0Short", "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
    "Phi", "KStar", "AntiKStar",
    "D0", "AntiD0", "DPlus", "DMinus", "D0s", "AntiD0s", "DStarPlus", "DStarMinus",
    "Lambdac", "AntiLambdac", "JPsi",
    "Pi0", "AntiPi0", "Eta", "AntiEta",
    "EtaPrime", "OmegaMeson",
    "Omegac", "Omegacbar", "Xic", "XicBar",
    "Xicc", "Xiccbar", "Omegacc", "Omegaccbar",
    "Omegaccc", "Omegacccbar",
    "Xib", "Xibbar",
    "Omegab", "Omegabbar",
    //Prompt measurements
    "PromptD0", "PromptAntiD0", "PromptDPlus", "PromptDMinus", "PromptD0s", "PromptAntiD0s", "PromptDStarPlus", "PromptDStarMinus",
    "PromptLambdac", "PromptAntiLambdac", "PromptJPsi",
    "PromptOmegac", "PromptOmegacbar", "PromptXic", "PromptXicBar",
    "PromptXicc", "PromptXiccbar", "PromptOmegacc", "PromptOmegaccbar",
    "PromptOmegaccc", "PromptOmegacccbar",
    "PromptXib", "PromptXibbar",
    "PromptOmegab", "PromptOmegabbar"
  };
  
  //Main Output: Histograms
  
  //Event counter histogram: Multiplicity, Npart, b (if available)
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPt[ih] ) {
      fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPt[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistEta[ih] ) {
      fHistEta[ih] = new TH1D(Form("fHistEta_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
      fListHist->Add(fHistEta[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistEtaTriggeredMeson[ih] ) {
      fHistEtaTriggeredMeson[ih] = new TH2D(Form("fHistEtaTriggeredMeson_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta, 10,0,10);
      fListHist->Add(fHistEtaTriggeredMeson[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistEtaTriggeredCharm[ih] ) {
      fHistEtaTriggeredCharm[ih] = new TH2D(Form("fHistEtaTriggeredCharm_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta, 10,0,10);
      fListHist->Add(fHistEtaTriggeredCharm[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistEtaTriggeredBeauty[ih] ) {
      fHistEtaTriggeredBeauty[ih] = new TH2D(Form("fHistEtaTriggeredBeauty_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta, 10,0,10);
      fListHist->Add(fHistEtaTriggeredBeauty[ih]);
    }
  }
  
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPtVsV0MMult[ih] ) {
      fHistPtVsV0MMult[ih] = new TH2D(Form("fHistPtVsV0MMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsV0MMult[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPtVsSPDMult[ih] ) {
      fHistPtVsSPDMult[ih] = new TH2D(Form("fHistPtVsSPDMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistEtaVsSPDMult[ih] && fkDoRapidityStudy ) {
      fHistEtaVsSPDMult[ih] = new TH2D(Form("fHistEtaVsSPDMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,1,-10,10);
      fListHist->Add(fHistEtaVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistYVsSPDMult[ih] && fkDoRapidityStudy ) {
      fHistYVsSPDMult[ih] = new TH2D(Form("fHistYVsSPDMult%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,1,-10,10);
      fListHist->Add(fHistYVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPtVsNpart[ih] && fkDoNpartStudy ) {
      fHistPtVsNpart[ih] = new TH2D(Form("fHistPtVsNpart_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",fkNNpartBins,-0.5,fkNNpartBins-0.5,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsNpart[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPtVsB[ih] && fkDoImpactParameterStudy ) {
      fHistPtVsB[ih] = new TH2D(Form("fHistPtVsB_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",fkNBBins,0,20,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsB[ih]);
    }
  }
  for(Int_t ih=0; ih<72; ih++){
    if(! fHistPtVsNMPI[ih] && fkDoNMPIStudy ) {
      fHistPtVsNMPI[ih] = new TH2D(Form("fHistPtVsNMPI_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",50,-0.5,49.5,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsNMPI[ih]);
    }
  }
  
  if(! fHistPtTriggerD0 ) {
    //Histogram Output: Event-by-Event
    fHistPtTriggerD0 = new TH1D( "fHistPtTriggerD0", ";#eta;Count",200,0,20);
    fListHist->Add(fHistPtTriggerD0);
  }
  if(! fHistPtTriggerXiC ) {
    //Histogram Output: Event-by-Event
    fHistPtTriggerXiC = new TH1D( "fHistPtTriggerXiC", ";#eta;Count",200,0,20);
    fListHist->Add(fHistPtTriggerXiC);
  }
  if(! fHistPtTriggerXiB ) {
    //Histogram Output: Event-by-Event
    fHistPtTriggerXiB = new TH1D( "fHistPtTriggerXiB", ";#eta;Count",200,0,20);
    fListHist->Add(fHistPtTriggerXiB);
  }
  
  if(! fHist3d2pcD0Proton ) {
    fHist3d2pcD0Proton = new TH3D("fHist3d2pcD0Proton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0Proton);
  }
  if(! fHist3d2pcD0AntiProton ) {
    fHist3d2pcD0AntiProton = new TH3D("fHist3d2pcD0AntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0AntiProton);
  }
  if(! fHist3d2pcD0D0 ) {
    fHist3d2pcD0D0 = new TH3D("fHist3d2pcD0D0","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0D0);
  }
  if(! fHist3d2pcD0D0bar ) {
    fHist3d2pcD0D0bar = new TH3D("fHist3d2pcD0D0bar","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0D0bar);
  }
  if(! fHist3d2pcD0KMinus ) {
    fHist3d2pcD0KMinus = new TH3D("fHist3d2pcD0KMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0KMinus);
  }
  if(! fHist3d2pcD0KPlus ) {
    fHist3d2pcD0KPlus = new TH3D("fHist3d2pcD0KPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcD0KPlus);
  }
  
  if(! fHist3d2pcXiCProton ) {
    fHist3d2pcXiCProton = new TH3D("fHist3d2pcXiCProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCProton);
  }
  if(! fHist3d2pcXiCAntiProton ) {
    fHist3d2pcXiCAntiProton = new TH3D("fHist3d2pcXiCAntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCAntiProton);
  }
  if(! fHist3d2pcXiCD0 ) {
    fHist3d2pcXiCD0 = new TH3D("fHist3d2pcXiCD0","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCD0);
  }
  if(! fHist3d2pcXiCD0bar ) {
    fHist3d2pcXiCD0bar = new TH3D("fHist3d2pcXiCD0bar","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCD0bar);
  }
  if(! fHist3d2pcXiCKMinus ) {
    fHist3d2pcXiCKMinus = new TH3D("fHist3d2pcXiCKMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCKMinus);
  }
  if(! fHist3d2pcXiCKPlus ) {
    fHist3d2pcXiCKPlus = new TH3D("fHist3d2pcXiCKPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiCKPlus);
  }
  
  if(! fHist3d2pcXiBProton ) {
    fHist3d2pcXiBProton = new TH3D("fHist3d2pcXiBProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBProton);
  }
  if(! fHist3d2pcXiBAntiProton ) {
    fHist3d2pcXiBAntiProton = new TH3D("fHist3d2pcXiBAntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBAntiProton);
  }
  if(! fHist3d2pcXiBBMinus ) {
    fHist3d2pcXiBBMinus = new TH3D("fHist3d2pcXiBBMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBBMinus);
  }
  if(! fHist3d2pcXiBBPlus ) {
    fHist3d2pcXiBBPlus = new TH3D("fHist3d2pcXiBBPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBBPlus);
  }
  if(! fHist3d2pcXiBKMinus ) {
    fHist3d2pcXiBKMinus = new TH3D("fHist3d2pcXiBKMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBKMinus);
  }
  if(! fHist3d2pcXiBKPlus ) {
    fHist3d2pcXiBKPlus = new TH3D("fHist3d2pcXiBKPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHist3d2pcXiBKPlus);
  }
  
  //event mixing
  if(! fHistMixed3d2pcD0Proton ) {
    fHistMixed3d2pcD0Proton = new TH3D("fHistMixed3d2pcD0Proton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0Proton);
  }
  if(! fHistMixed3d2pcD0AntiProton ) {
    fHistMixed3d2pcD0AntiProton = new TH3D("fHistMixed3d2pcD0AntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0AntiProton);
  }
  if(! fHistMixed3d2pcD0D0 ) {
    fHistMixed3d2pcD0D0 = new TH3D("fHistMixed3d2pcD0D0","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0D0);
  }
  if(! fHistMixed3d2pcD0D0bar ) {
    fHistMixed3d2pcD0D0bar = new TH3D("fHistMixed3d2pcD0D0bar","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0D0bar);
  }
  if(! fHistMixed3d2pcD0KMinus ) {
    fHistMixed3d2pcD0KMinus = new TH3D("fHistMixed3d2pcD0KMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0KMinus);
  }
  if(! fHistMixed3d2pcD0KPlus ) {
    fHistMixed3d2pcD0KPlus = new TH3D("fHistMixed3d2pcD0KPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcD0KPlus);
  }
  
  if(! fHistMixed3d2pcXiCProton ) {
    fHistMixed3d2pcXiCProton = new TH3D("fHistMixed3d2pcXiCProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCProton);
  }
  if(! fHistMixed3d2pcXiCAntiProton ) {
    fHistMixed3d2pcXiCAntiProton = new TH3D("fHistMixed3d2pcXiCAntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCAntiProton);
  }
  if(! fHistMixed3d2pcXiCD0 ) {
    fHistMixed3d2pcXiCD0 = new TH3D("fHistMixed3d2pcXiCD0","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCD0);
  }
  if(! fHistMixed3d2pcXiCD0bar ) {
    fHistMixed3d2pcXiCD0bar = new TH3D("fHistMixed3d2pcXiCD0bar","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCD0bar);
  }
  if(! fHistMixed3d2pcXiCKMinus ) {
    fHistMixed3d2pcXiCKMinus = new TH3D("fHistMixed3d2pcXiCKMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCKMinus);
  }
  if(! fHistMixed3d2pcXiCKPlus ) {
    fHistMixed3d2pcXiCKPlus = new TH3D("fHistMixed3d2pcXiCKPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiCKPlus);
  }
  
  if(! fHistMixed3d2pcXiBProton ) {
    fHistMixed3d2pcXiBProton = new TH3D("fHistMixed3d2pcXiBProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBProton);
  }
  if(! fHistMixed3d2pcXiBAntiProton ) {
    fHistMixed3d2pcXiBAntiProton = new TH3D("fHistMixed3d2pcXiBAntiProton","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBAntiProton);
  }
  if(! fHistMixed3d2pcXiBBMinus ) {
    fHistMixed3d2pcXiBBMinus = new TH3D("fHistMixed3d2pcXiBBMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBBMinus);
  }
  if(! fHistMixed3d2pcXiBBPlus ) {
    fHistMixed3d2pcXiBBPlus = new TH3D("fHistMixed3d2pcXiBBPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBBPlus);
  }
  if(! fHistMixed3d2pcXiBKMinus ) {
    fHistMixed3d2pcXiBKMinus = new TH3D("fHistMixed3d2pcXiBKMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBKMinus);
  }
  if(! fHistMixed3d2pcXiBKPlus ) {
    fHistMixed3d2pcXiBKPlus = new TH3D("fHistMixed3d2pcXiBKPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),10,0,10);
    fListHist->Add(fHistMixed3d2pcXiBKPlus);
  }
  
  //List of Histograms: Normal
  PostData(1, fListHist);
  
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictions::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  
  
  // Connect to the InputEvent
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.
  
  // Appropriate for ESD analysis!
  
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  //------------------------------------------------
  // Multiplicity Information Acquistion
  //------------------------------------------------
  
  //Monte Carlo Level information !
  //--------- GENERATED NUMBER OF CHARGED PARTICLES
  // ---> Variable Definition
  
  Long_t lNchEta5   = 0;
  Long_t lNchEta8   = 0;
  Long_t lNchEta8to15   = 0;
  Long_t lNchEta10  = 0;
  Long_t lNchEtaWide  = 0;
  Long_t lNchVZEROA = 0;
  Long_t lNchVZEROC = 0;
  Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
  
  Double_t lWideEta = 1.4;
  if( fkALICE3SiliconMode ) lWideEta = 4.0; //ALICE 3 mode: |eta|<4 => "SPD"
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {   // This is the begining of the loop on tracks
    TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
    if(!particleOne) continue;
    if(!particleOne->GetPDG()) continue;
    Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
    if(TMath::Abs(lThisCharge)<0.001) continue;
    if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
    
    Double_t gpt = particleOne -> Pt();
    Double_t geta = particleOne -> Eta();
    
    if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
    if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
    if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta8to15++;
    if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
    if( TMath::Abs(geta) < lWideEta ) lNchEtaWide++;
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
    if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
    if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------
  
  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
  
  if( lNchEtaWide < fkMinimumMultiplicity ) return;
  
  //------------------------------------------------
  // Acquire information on Npart, Ncoll, b
  //------------------------------------------------
  
  //Npart and Ncoll information
  AliGenHijingEventHeader* hHijing=0;
  AliGenDPMjetEventHeader* dpmHeader=0;
  AliGenEventHeader* mcGenH = lMCevent->GenEventHeader();
  
  Int_t fMC_NPart = -1;
  Int_t fMC_NColl = -1;
  Float_t fMC_b = -1;
  Int_t fMC_NMPI = -1;
  
  if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
    AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
    if(fMcPythiaHeader){
      fMC_NMPI = fMcPythiaHeader->GetNMPI();
    }
  }
  
  //DPMJet/HIJING info if available
  if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class()))
    hHijing = (AliGenHijingEventHeader*)mcGenH;
  else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
    TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
    hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing"));
    if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing pPb_0"));
    if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing_0"));
  }
  else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
    dpmHeader = (AliGenDPMjetEventHeader*)mcGenH;
  }
  if(hHijing)   {
    fMC_NPart = hHijing->ProjectileParticipants()+hHijing->TargetParticipants();
    fMC_NColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
  }
  if(dpmHeader) {
    fMC_NPart =dpmHeader->ProjectileParticipants()+dpmHeader->TargetParticipants();
    fMC_NColl =dpmHeader->NN()+dpmHeader->NNw()+dpmHeader->NwN()+dpmHeader->NwNw();
  }
  
  //check EPOS info, if available
  if ( IsEPOSLHC() ){
    AliGenHepMCEventHeader *lHepMCHeader = 0x0;
    if (mcGenH->InheritsFrom(AliGenHepMCEventHeader::Class()))
      lHepMCHeader = (AliGenHepMCEventHeader*)mcGenH;
    
    if (lHepMCHeader ){
      fMC_NPart = lHepMCHeader->Npart_proj()+lHepMCHeader->Npart_targ();
      fMC_NColl = lHepMCHeader->N_Nwounded_collisions() +
      lHepMCHeader->Nwounded_N_collisions() +
      lHepMCHeader->Nwounded_Nwounded_collisions();
      
      fMC_b = lHepMCHeader->impact_parameter();
    }
  }
  
  //------------------------------------------------
  // Fill Event Counters
  //------------------------------------------------
  
  //Basics: All Processed
  fHistEventCounter->Fill(0.5);
  
  if(fHistV0MMult)      fHistV0MMult        -> Fill ( lNchVZEROA+lNchVZEROC );
  if(fHistSPDMult)      fHistSPDMult        -> Fill ( lNchEtaWide );
  if(fHistNchVsV0MMult) fHistNchVsV0MMult   -> Fill ( lNchVZEROA+lNchVZEROC, lNchEta5  );
  if(fHistNchVsSPDMult) fHistNchVsSPDMult   -> Fill ( lNchEtaWide, lNchEta5  );
  if(fHistNpart)        fHistNpart          -> Fill ( fMC_NPart );
  if(fHistNchVsNpart)   fHistNchVsNpart     -> Fill ( fMC_NPart, lNchEta5  );
  if(fHistB)            fHistB              -> Fill ( fMC_b );
  if(fHistNchVsB)       fHistNchVsB         -> Fill ( fMC_b, lNchEta5  );
  if(fHistNMPI)         fHistNMPI           -> Fill ( fMC_NMPI );
  if(fHistNchVsNMPI)    fHistNchVsNMPI      -> Fill ( fMC_NMPI, lNchEta5  );
  
  //------------------------------------------------
  // Fill Spectra as Needed
  //------------------------------------------------
  
  //~All relevant PWG-LF Identified Particle Information (for looping)
  Int_t lPDGCodes[72] = {
    211, -211, 321, -321, 2212, -2212,
    310, 3122, -3122,
    3312, -3312, 3334, -3334,
    333, 313, -313,
    421, -421, 411, -411, 431, -431, 413, -413,
    4122, -4122, 443,
    111,-111,221,-221,
    331, 223,
    4332, -4332, 4232, -4232,
    4422, -4422, 4432, -4432,
    4444, -4444,
    5132, -5132,
    5332, -5332,
    421, -421, 411, -411, 431, -431, 413, -413,
    4122, -4122, 443,
    4332, -4332, 4232, -4232,
    4422, -4422, 4432, -4432,
    4444, -4444,
    5132, -5132,
    5332, -5332
  };
  TString lPartNames[72] = {
    "PiPlus", "PiMinus", "KaPlus", "KaMinus", "Proton", "AntiProton",
    "K0Short", "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
    "Phi", "KStar", "AntiKStar",
    "D0", "AntiD0", "DPlus", "DMinus", "D0s", "AntiD0s", "DStarPlus", "DStarMinus",
    "Lambdac", "AntiLambdac", "JPsi",
    "Pi0", "AntiPi0", "Eta", "AntiEta",
    "EtaPrime", "OmegaMeson",
    "Omegac", "Omegacbar", "Xic", "Xicbar",
    "Xicc", "Xiccbar", "Omegacc", "Omegaccbar",
    "Omegaccc", "Omegacccbar",
    "Xib", "Xibbar",
    "Omegab", "Omegabbar",
    //Prompt measurements
    "PromptD0", "PromptAntiD0", "PromptDPlus", "PromptDMinus", "PromptD0s", "PromptAntiD0s", "PromptDStarPlus", "PromptDStarMinus",
    "PromptLambdac", "PromptAntiLambdac", "PromptJPsi",
    "PromptOmegac", "PromptOmegacbar", "PromptXic", "PromptXicbar",
    "PromptXicc", "PromptXiccbar", "PromptOmegacc", "PromptOmegaccbar",
    "PromptOmegaccc", "PromptOmegacccbar",
    "PromptXib", "PromptXibbar",
    "PromptOmegab", "PromptOmegabbar"
  };
  Bool_t lCheckIsPhysicalPrimary[72] = {
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE
  };
  Bool_t lCheckHFFeeddown[72] = {
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE
  };
  
  Int_t lThisPDG  = 0;
  Double_t lThisRap  = 0;
  Double_t lThisPt   = 0;
  Bool_t lIsPhysicalPrimary = kFALSE;
  
  //===== Start 2pc =================
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // step 1: find particles to be correlated
  //         looking for XiC, antiprotons, D0bar, K+
  Long_t lNXiC=0, lNXiB=0;
  Long_t lNProtons=0, lNAntiProtons=0;
  Long_t lND0=0, lND0bar=0;
  Long_t lND0trigger=0, lND0bartrigger=0;
  Long_t lNBMinus=0, lNBPlus=0;
  Long_t lNKMinus=0, lNKPlus=0;
  
  TArrayI lXiB( lMCstack->GetNtrack() );
  TArrayI lXiC( lMCstack->GetNtrack() );
  TArrayI lProtons( lMCstack->GetNtrack() );
  TArrayI lAntiProtons( lMCstack->GetNtrack() );
  TArrayI lD0( lMCstack->GetNtrack() );
  TArrayI lD0bar( lMCstack->GetNtrack() );
  TArrayI lD0trigger( lMCstack->GetNtrack() );
  TArrayI lD0bartrigger( lMCstack->GetNtrack() );
  TArrayI lBMinus( lMCstack->GetNtrack() );
  TArrayI lBPlus( lMCstack->GetNtrack() );
  TArrayI lKMinus( lMCstack->GetNtrack() );
  TArrayI lKPlus( lMCstack->GetNtrack() );
  
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {
    // Determine if within acceptance, otherwise fully reject from list
    // done such that this check is done O(N) and not O(N^2)
    TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
    AliMCParticle* lMCPart = (AliMCParticle*)lMCevent->GetTrack(iCurrentLabelStack);
    if(!lThisParticle) continue;
    lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(iCurrentLabelStack);
    Double_t geta = lThisParticle -> Eta();
    Double_t gpt = lThisParticle -> Pt();
    
    //gotta reject any and all offspring of decays
    //simplest implementation: reject based on decay position
    Double_t lDistanceFromZero = TMath::Sqrt(
                                             TMath::Power( lThisParticle->Vx() , 2) +
                                             TMath::Power( lThisParticle->Vy() , 2) +
                                             TMath::Power( lThisParticle->Vz() , 2)
                                             );
    
    if(lDistanceFromZero>1e-12) continue; //remove everything outside of zero, should remove decay daus
    
    if( TMath::Abs(geta)<4.0 ){
      if(lThisParticle->GetPdgCode()==421) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        if( gpt < fMinPtTrigger || fMaxPtTrigger < gpt ) continue;
        lD0trigger[lND0trigger++] = iCurrentLabelStack;
      }
      if(lThisParticle->GetPdgCode()==4232) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        if( gpt < fMinPtTrigger || fMaxPtTrigger < gpt ) continue;
        lXiC[lNXiC++] = iCurrentLabelStack;
      }
      if(lThisParticle->GetPdgCode()==5132) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        if( gpt < fMinPtTrigger || fMaxPtTrigger < gpt ) continue;
        lXiB[lNXiB++] = iCurrentLabelStack;
      }
      if(lThisParticle->GetPdgCode()== 2212 && lIsPhysicalPrimary ) lProtons[lNProtons++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()==-2212 && lIsPhysicalPrimary ) lAntiProtons[lNAntiProtons++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()==  421 && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)==4 ) lD0[lND0++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()== -421 && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)==4 ) lD0bar[lND0bar++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()== +521 && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)==4 ) lBPlus[lNBPlus++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()== -521 && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)==4 ) lBMinus[lNBMinus++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()== +321 && lIsPhysicalPrimary ) lKPlus[lNKPlus++] = iCurrentLabelStack;
      if(lThisParticle->GetPdgCode()== -321 && lIsPhysicalPrimary ) lKMinus[lNKMinus++] = iCurrentLabelStack;
    }
  }
  //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  //----- Loop on Stack Starts Here ---------------
  for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
  {   // This is the begining of the loop on tracks
    
    TParticle* lPart = 0x0;
    lPart = lMCstack->Particle( ilab );
    //For CheckOrigin
    AliMCParticle* lMCPart = (AliMCParticle*)lMCevent->GetTrack(ilab);
    if(!lPart) {
      Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
      continue;
    }
    
    lThisPDG = lPart->GetPdgCode();
    //Continue if this is not a particle of the right PDG Code (avoids y-calculation problems)
    Bool_t lContinue = kTRUE;
    for(Int_t ih=0; ih<72; ih++) if( lThisPDG == lPDGCodes[ih] ) lContinue = kFALSE;
    if ( lContinue ) continue;
    
    lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
    //lThisRap   = lPart->Y();
    lThisPt    = lPart->Pt();
    
    //Use Physical Primaries only for filling These Histos
    //if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;
    lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(ilab);
    
    Double_t lDistanceFromZero = TMath::Sqrt(
                                             TMath::Power( lPart->Vx() , 2) +
                                             TMath::Power( lPart->Vy() , 2) +
                                             TMath::Power( lPart->Vz() , 2)
                                             );
    
    for(Int_t ih=0; ih<72; ih++){
      if( lThisPDG == lPDGCodes[ih] ) {
        //Check if primary (if needed) and if not don't use this particle
        if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
        if( lCheckHFFeeddown[ih] == kTRUE && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        //Fill Histograms
        if(fHistEta[ih] && lDistanceFromZero > 1e-12 ) fHistEta[ih] -> Fill ( lPart -> Eta() );
        if(fHistEtaTriggeredMeson[ih] && lNXiC > 0 && lDistanceFromZero > 1e-12 ) fHistEtaTriggeredMeson[ih] -> Fill ( lPart -> Eta(), lThisPt );
        if(fHistEtaTriggeredCharm[ih] && lNXiC > 0 && lDistanceFromZero > 1e-12 ) fHistEtaTriggeredCharm[ih] -> Fill ( lPart -> Eta(), lThisPt );
        if(fHistEtaTriggeredBeauty[ih] && lNXiB > 0 && lDistanceFromZero > 1e-12 ) fHistEtaTriggeredBeauty[ih] -> Fill ( lPart -> Eta(), lThisPt );
        if(fHistEtaVsSPDMult[ih]) fHistEtaVsSPDMult[ih] -> Fill( lNchEtaWide, lPart -> Eta() );
        if(fHistYVsSPDMult[ih]) fHistYVsSPDMult[ih] -> Fill( lNchEtaWide, lThisRap );
        if( TMath::Abs(lThisRap) < 0.5 && !fkWideRapiditySpeciesStudy ) {
          if( fHistPt[ih] ) fHistPt[ih]->Fill(lThisPt);
          if( fHistPtVsV0MMult[ih] ) fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
          if( fHistPtVsSPDMult[ih] ) fHistPtVsSPDMult[ih]->Fill(lNchEtaWide,lThisPt);
          if( fHistPtVsNpart[ih] ) fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
          if( fHistPtVsB[ih] ) fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
          if( fHistPtVsNMPI[ih] ) fHistPtVsNMPI[ih]->Fill(fMC_NMPI,lThisPt);
        }
        if( TMath::Abs(lThisRap) < 4.0 && fkWideRapiditySpeciesStudy ) {
          if( fHistPt[ih] ) fHistPt[ih]->Fill(lThisPt);
          if( fHistPtVsV0MMult[ih] ) fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
          if( fHistPtVsSPDMult[ih] ) fHistPtVsSPDMult[ih]->Fill(lNchEtaWide,lThisPt);
          if( fHistPtVsNpart[ih] ) fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
          if( fHistPtVsB[ih] ) fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
          if( fHistPtVsNMPI[ih] ) fHistPtVsNMPI[ih]->Fill(fMC_NMPI,lThisPt);
        }
        
      }
    }
  }//End of loop on tracks
  //----- End Loop on Stack ----------------------
  
  //  Actually correlate stuff with stuff
  for (Int_t iTrigger = 0;  iTrigger < lND0trigger; iTrigger++){   // trigger loop
    TParticle* lTriggerParticle = lMCstack->Particle(lD0trigger[iTrigger]);
    
    Double_t geta = lTriggerParticle -> Eta();
    Double_t gphi = lTriggerParticle -> Phi();
    fHistPtTriggerD0->Fill( lTriggerParticle -> Pt() );
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0Proton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0AntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lND0; iassoc++){   // associated loop
      if (lD0[iassoc] == lD0trigger[iTrigger]) continue;
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lD0[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0D0->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lND0bar; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lD0bar[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0D0bar->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0KMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcD0KPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  }
  
  //Event mixing for D0
  if( fEMBufferFullD0 && lND0trigger>0 ){ //require also that a trigger exists
    for (Int_t iTrigger = 0;  iTrigger < 10; iTrigger++){   // trigger loop
      Double_t geta = fEMBufferEtaD0[iTrigger]; //from previous events
      Double_t gphi = fEMBufferPhiD0[iTrigger]; //from previous events
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0Proton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0AntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lND0; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lD0[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0D0->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lND0bar; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lD0bar[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0D0bar->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0KMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcD0KPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    }
  }
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iTrigger = 0;  iTrigger < lNXiC; iTrigger++){   // trigger loop
    TParticle* lTriggerParticle = lMCstack->Particle(lXiC[iTrigger]);
    
    Double_t geta = lTriggerParticle -> Eta();
    Double_t gphi = lTriggerParticle -> Phi();
    fHistPtTriggerXiC->Fill( lTriggerParticle -> Pt() );
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCAntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lND0; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lD0[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCD0->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lND0bar; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lD0bar[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCD0bar->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCKMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiCKPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  }
  
  //Event mixing for XiC
  if( fEMBufferFullXiC && lNXiC > 0 ){ //demand event triggered
    for (Int_t iTrigger = 0;  iTrigger < 10; iTrigger++){   // trigger loop
      Double_t geta = fEMBufferEtaXiC[iTrigger]; //from previous events
      Double_t gphi = fEMBufferPhiXiC[iTrigger]; //from previous events
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCAntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lND0; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lD0[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCD0->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lND0bar; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lD0bar[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCD0bar->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCKMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiCKPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    }
  }//end event mixing loop for XiC
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iTrigger = 0;  iTrigger < lNXiB; iTrigger++){   // trigger loop
    TParticle* lTriggerParticle = lMCstack->Particle(lXiB[iTrigger]);
    
    Double_t geta = lTriggerParticle -> Eta();
    Double_t gphi = lTriggerParticle -> Phi();
    fHistPtTriggerXiB->Fill( lTriggerParticle -> Pt() );
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBAntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNBPlus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lBPlus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBBPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNBMinus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lBMinus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBBMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBKMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
      TParticle* lAssociatedParticle = 0x0;
      lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
      if(!lAssociatedParticle) {
        Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
        continue;
      }
      
      Double_t geta2 = lAssociatedParticle -> Eta();
      Double_t gphi2 = lAssociatedParticle -> Phi();
      lThisPt    = lAssociatedParticle->Pt();
      fHist3d2pcXiBKPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
    }
    //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  }
  
  //Event mixing loop for XiB
  if( fEMBufferFullXiB && lNXiB > 0 ){ //demand event triggered
    for (Int_t iTrigger = 0;  iTrigger < 10; iTrigger++){   // trigger loop
      Double_t geta = fEMBufferEtaXiB[iTrigger]; //from previous events
      Double_t gphi = fEMBufferPhiXiB[iTrigger]; //from previous events
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNAntiProtons; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lAntiProtons[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBAntiProton->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNBPlus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lBPlus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBBPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNBMinus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lBMinus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBBMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKMinus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKMinus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBKMinus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      for (Int_t iassoc = 0;  iassoc < lNKPlus; iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( lKPlus[iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        lThisPt    = lAssociatedParticle->Pt();
        fHistMixed3d2pcXiBKPlus->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
      //  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    }
  }//end event mixing
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // fill EM buffer
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {
    // Determine if within acceptance, otherwise fully reject from list
    // done such that this check is done O(N) and not O(N^2)
    TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
    AliMCParticle* lMCPart = (AliMCParticle*)lMCevent->GetTrack(iCurrentLabelStack);
    if(!lThisParticle) continue;
    lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(iCurrentLabelStack);
    Double_t geta = lThisParticle -> Eta();
    Double_t gpt = lThisParticle -> Pt();
    
    //gotta reject any and all offspring of decays
    //simplest implementation: reject based on decay position
    Double_t lDistanceFromZero = TMath::Sqrt(
                                             TMath::Power( lThisParticle->Vx() , 2) +
                                             TMath::Power( lThisParticle->Vy() , 2) +
                                             TMath::Power( lThisParticle->Vz() , 2)
                                             );
    
    if(lDistanceFromZero>1e-12) continue; //remove everything outside of zero, should remove decay daus
    
    if( gpt < fMinPtTrigger || fMaxPtTrigger < gpt ) continue;
    
    if( TMath::Abs(geta)<4.0 ){
      if(lThisParticle->GetPdgCode()==421) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        //Add to buffer
        fEMBufferEtaD0[ fEMBufferCycleD0 ] = lThisParticle->Eta();
        fEMBufferPhiD0[ fEMBufferCycleD0 ] = lThisParticle->Phi();
        fEMBufferCycleD0++;
        if(fEMBufferCycleD0>=10) fEMBufferFullD0 = kTRUE;
        fEMBufferCycleD0 = fEMBufferCycleD0%10;
      }
      if(lThisParticle->GetPdgCode()==4232) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        //Add to buffer
        fEMBufferEtaXiC[ fEMBufferCycleXiC ] = lThisParticle->Eta();
        fEMBufferPhiXiC[ fEMBufferCycleXiC ] = lThisParticle->Phi();
        fEMBufferCycleXiC++;
        if(fEMBufferCycleXiC>=10) fEMBufferFullXiC = kTRUE;
        fEMBufferCycleXiC = fEMBufferCycleXiC%10;
      }
      if(lThisParticle->GetPdgCode()==5132) {
        if( AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, kTRUE)!=4 ) continue;
        fEMBufferEtaXiB[ fEMBufferCycleXiB ] = lThisParticle->Eta();
        fEMBufferPhiXiB[ fEMBufferCycleXiB ] = lThisParticle->Phi();
        fEMBufferCycleXiB++;
        if(fEMBufferCycleXiB>=10) fEMBufferFullXiB = kTRUE;
        fEMBufferCycleXiB = fEMBufferCycleXiB%10;
      }
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  //===== Start 2pc nested loops =================
  //  if( fkDo2pc ) {
  //    //Apply the eta cut first or go home
  //    Long_t lNValidParticles = 0;
  //    TArrayI lValidParticles(lMCstack->GetNtrack());
  //    Long_t lNValidPhi = 0;
  //    TArrayI lValidPhi(lMCstack->GetNtrack());
  //    Long_t lNValidXi = 0;
  //    TArrayI lValidXi(lMCstack->GetNtrack());
  //    //----- Determine valid triggers ----------------------------------------------------------------
  //    for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  //    {
  //      // Determine if within acceptance, otherwise fully reject from list
  //      // done such that this check is done O(N) and not O(N^2)
  //      TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
  //      if(!lThisParticle) continue;
  //      Double_t geta = lThisParticle -> Eta();
  //      if( TMath::Abs(geta)<0.8 ) lValidParticles[lNValidParticles++]=iCurrentLabelStack;
  //    }
  //    //----- Loop on Stack ----------------------------------------------------------------
  //    for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < lNValidParticles; iCurrentLabelStack++)
  //    {   // This is the begining of the loop on tracks
  //      TParticle* lTriggerParticle = lMCstack->Particle(lValidParticles[iCurrentLabelStack]);
  //      if(!lTriggerParticle) continue;
  //      if(!lTriggerParticle->GetPDG()) continue;
  //      Double_t lThisCharge = lTriggerParticle->GetPDG()->Charge()/3.;
  //      //if(TMath::Abs(lThisCharge)<0.001) continue;
  //      //if(! (lMCstack->IsPhysicalPrimary(lValidParticles[iCurrentLabelStack])) ) continue;
  //
  //      Bool_t lTrigIsCharged = kTRUE;
  //      if( TMath::Abs(lThisCharge)<0.001 ) lTrigIsCharged = kFALSE;
  //      Bool_t lTrigIsPrimary = kTRUE;
  //      if ( !lMCstack->IsPhysicalPrimary(lValidParticles[iCurrentLabelStack]) ) lTrigIsPrimary = kFALSE;
  //      Bool_t lTrigIsPhi = kTRUE;
  //      if (lTriggerParticle->GetPdgCode()!=333) lTrigIsPhi = kFALSE;
  //
  //      if( ((!lTrigIsCharged)||(!lTrigIsPrimary)) && !lTrigIsPhi ) continue;
  //
  //      Double_t geta = lTriggerParticle -> Eta();
  //      Double_t gphi = lTriggerParticle -> Phi();
  //
  //      if( lTriggerParticle -> Pt() > fMinPtTriggerCharged && lTrigIsCharged && lTrigIsPrimary )
  //        fEtaTriggerCharged -> Fill( geta );
  //      if( lTriggerParticle -> Pt() > fMinPtTriggerXi && lTrigIsPrimary && TMath::Abs(lTriggerParticle->GetPdgCode())==3312 )
  //        fEtaTriggerXi      -> Fill( geta );
  //      if( lTriggerParticle -> Pt() > fMinPtTriggerPhi && TMath::Abs(lTriggerParticle->GetPdgCode())==333 )
  //        fEtaTriggerPhi     -> Fill( geta );
  //
  //      for (Int_t ilab = 0;  ilab < lNValidParticles; ilab++)
  //      {   // This is the begining of the loop on tracks
  //
  //        if(ilab == iCurrentLabelStack) continue; //remove auto-correlations
  //        TParticle* lAssociatedParticle = 0x0;
  //        lAssociatedParticle = lMCstack->Particle( lValidParticles[ilab] );
  //        if(!lAssociatedParticle) {
  //          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
  //          continue;
  //        }
  //
  //        lThisPDG = lAssociatedParticle->GetPdgCode();
  //
  //        //Continue if this is not a particle of the right PDG Code (avoids y-calculation problems)
  //        Bool_t lContinue = kTRUE;
  //        for(Int_t ih=0; ih<52; ih++) if( lThisPDG == lPDGCodes[ih] ) lContinue = kFALSE;
  //        if ( lContinue ) continue;
  //
  //        Double_t geta2 = lAssociatedParticle -> Eta();
  //        Double_t gphi2 = lAssociatedParticle -> Phi();
  //
  //        lThisPt    = lAssociatedParticle->Pt();
  //
  //        lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(lValidParticles[ilab]);
  //
  //        if( lTrigIsCharged && lTrigIsPrimary ){
  //          for(Int_t ih=0; ih<52; ih++){
  //            if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
  //              //Check if primary (if needed) and if not don't use this particle
  //              if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
  //              //Fill 2pc same-event histograms, please
  //              fHist3d2pcSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
  //            }
  //          }
  //        }
  //        if( lTriggerParticle->GetPdgCode() == 3312 ){
  //          for(Int_t ih=0; ih<52; ih++){
  //            if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
  //              //Check if primary (if needed) and if not don't use this particle
  //              if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
  //              //Fill 2pc same-event histograms, please
  //              fHist3d2pcXiSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
  //            }
  //          }
  //        }
  //        if( TMath::Abs( lTriggerParticle->GetPdgCode() ) == 333 ){
  //          for(Int_t ih=0; ih<52; ih++){
  //            if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
  //              //Check if primary (if needed) and if not don't use this particle
  //              if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
  //              //Fill 2pc same-event histograms, please
  //              fHist3d2pcPhiSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
  //            }
  //          }
  //        }
  //      }//End of loop on tracks
  //    }//End of loop on tracks
  //    //----- End Loop on Stack ------------------------------------------------------------
  //  }
  //===== End 2pc nested loops ===================
  
  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
    Printf("ERROR - AliAnalysisTaskMCPredictions : ouput data container list not available\n");
    return;
  }
  
  fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
  if (!fHistEventCounter) {
    Printf("ERROR - AliAnalysisTaskMCPredictions : fHistEventCounter not available");
    return;
  }
  
  TCanvas *canCheck = new TCanvas("AliAnalysisTaskMCPredictions","Event Multiplicity",10,10,510,510);
  canCheck->cd(1)->SetLogy();
  
  fHistEventCounter->SetMarkerStyle(22);
  fHistEventCounter->DrawCopy("E");
}

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
    ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
  }
  return ReturnValue;
}


//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsHijing() const {
  //Function to check if this is Hijing MC
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
      //Option 1: Just Hijing
      lReturnValue = kTRUE;
    } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
      //Option 2: cocktail involving Hijing
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      TIter next(headers);
      while (const TObject *obj=next()){
        //Look for an object inheriting from the hijing header class
        if ( obj->InheritsFrom(AliGenHijingEventHeader::Class()) ){ lReturnValue = kTRUE; }
      }
    }
  }
  return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsDPMJet() const {
  //Function to check if this is DPMJet
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
      //DPMJet Header is there!
      lReturnValue = kTRUE;
    }
  }
  return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsEPOSLHC() const {
  //Function to check if this is DPMJet
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    //A bit uncivilized, but hey, if it works...
    TString lHeaderTitle = mcGenH->GetName();
    if (lHeaderTitle.Contains("EPOSLHC")) {
      //This header has "EPOS" in its title!
      lReturnValue = kTRUE;
    }
  }
  return lReturnValue;
}

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions::ComputeDeltaPhi( Double_t phi1, Double_t phi2) const {
  //To be completely sure, use inner products
  Double_t x1, y1, x2, y2;
  x1 = TMath::Cos( phi1 );
  y1 = TMath::Sin( phi1 );
  x2 = TMath::Cos( phi2 );
  y2 = TMath::Sin( phi2 );
  Double_t lInnerProd  = x1*x2 + y1*y2;
  Double_t lVectorProd = x1*y2 - x2*y1;
  
  Double_t lReturnVal = 0;
  if( lVectorProd > 1e-8 ){
    lReturnVal = TMath::ACos(lInnerProd);
  }
  if( lVectorProd < -1e-8 ){
    lReturnVal = -TMath::ACos(lInnerProd);
  }
  if( lReturnVal < -TMath::Pi()/2 ) lReturnVal += 2*TMath::Pi();
  return lReturnVal;
}
