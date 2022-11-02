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
//#include "AliPythia8.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions)

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions()
: AliAnalysisTaskSE(),
fListHist(nullptr), fTree(nullptr), fkSaveTree(kFALSE),
//Variables for fTree
fPt(0),
fPID(0),
fStatus(0),
fNParents(0),
fV0MMultiplicity(0),
fHistEventCounter(nullptr),
fHistChargedEta(nullptr),
fSmallMultRange(1000),
fLargeMultRange(2000),
fRebinFactor(1),
fkNBBins(1),
fkNNpartBins(1),
fkNEtaBins(1),
fkNSpecies(78),
fkNIntervals(1),
fkSelectINELgtZERO(kTRUE),
fkWideRapiditySpeciesStudy(kFALSE),
fkDoImpactParameterStudy(kFALSE),
fkDoNpartStudy(kFALSE),
fkDoNMPIStudy(kTRUE),
fkDoRapidityStudy(kFALSE),
fkMinimumMultiplicity(-1),
fCheckOriginThirdArgument(kTRUE),
fHistV0MMult(nullptr),
fHistSPDMult(nullptr),
fHistNchVsV0MMult(nullptr),
fHistNchVsSPDMult(nullptr),
fHistNpart(nullptr),
fHistNchVsNpart(nullptr),
fHistB(nullptr),
fHistNchVsB(nullptr),
fHistNMPI(nullptr),
fHistNchVsNMPI(nullptr),
fkDDRebin(1),
fkMaxMultDDV0M(fLargeMultRange),
fkMaxMultDDSPD(fLargeMultRange),
fHistV0MvsSPD(nullptr),
fHistDDNch(nullptr),
fHistDDNMPI(nullptr),
fHistDDQ2(nullptr)
{
  for(Int_t ii=0; ii<10; ii++) fPIDMother[ii]=-1;
  for(Int_t ii=0; ii<10; ii++) fStatusMother[ii]=-1;
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii] = 0;
    fkIntervalMaxEta[ii] = 0;
    fkIntervalWeight[ii] = 1.0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;
  for(Int_t ih=0; ih<78; ih++){
    fkSpeciesSwitch[ih] = kTRUE;
    fHistPt[ih]          = nullptr;
    fHistEta[ih]         = nullptr;
    fHistPtVsV0MMult[ih] = nullptr;
    fHistPtVsSPDMult[ih] = nullptr;
    fHistEtaVsSPDMult[ih]= nullptr;
    fHistYVsSPDMult[ih]  = nullptr;
    fHistPtVsNpart[ih]   = nullptr;
    fHistPtVsB[ih]       = nullptr;
    fHistPtVsNMPI[ih]   = nullptr;
    fHistDDYield[ih]   = nullptr;
    fHistDDPt[ih]   = nullptr;
  }
}

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions(const char *name, Int_t lNSmallBinning, Int_t lNLargeBinning, Int_t lRebinFactor, Int_t lNBBins, Int_t lNNpartBins, Int_t lNEtaBins)
: AliAnalysisTaskSE(name),
fListHist(nullptr), fTree(nullptr), fkSaveTree(kFALSE),
fPt(0),
fPID(0),
fStatus(0),
fNParents(0),
fV0MMultiplicity(0),
fHistEventCounter(nullptr),
fHistChargedEta(nullptr),
fSmallMultRange(lNSmallBinning),
fLargeMultRange(lNLargeBinning),
fRebinFactor(lRebinFactor),
fkNBBins(lNBBins),
fkNNpartBins(lNNpartBins),
fkNEtaBins(lNEtaBins),
fkNSpecies(78),
fkNIntervals(1),
fkSelectINELgtZERO(kTRUE),
fkWideRapiditySpeciesStudy(kFALSE),
fkDoImpactParameterStudy(kFALSE),
fkDoNpartStudy(kFALSE),
fkDoNMPIStudy(kTRUE),
fkDoRapidityStudy(kFALSE),
fkMinimumMultiplicity(-1),
fCheckOriginThirdArgument(kTRUE),
fHistV0MMult(nullptr),
fHistSPDMult(nullptr),
fHistNchVsV0MMult(nullptr),
fHistNchVsSPDMult(nullptr),
fHistNpart(nullptr),
fHistNchVsNpart(nullptr),
fHistB(nullptr),
fHistNchVsB(nullptr),
fHistNMPI(nullptr),
fHistNchVsNMPI(nullptr),
fkDDRebin(lRebinFactor),
fkMaxMultDDV0M(fLargeMultRange),
fkMaxMultDDSPD(fLargeMultRange),
fHistV0MvsSPD(nullptr),
fHistDDNch(nullptr),
fHistDDNMPI(nullptr),
fHistDDQ2(nullptr)
{
  for(Int_t ii=0; ii<10; ii++) fPIDMother[ii]=-1;
  for(Int_t ii=0; ii<10; ii++) fStatusMother[ii]=-1;
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii] = 0;
    fkIntervalMaxEta[ii] = 0;
    fkIntervalWeight[ii] = 1.0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;
  for(Int_t ih=0; ih<78; ih++){
    fkSpeciesSwitch[ih] = kTRUE;
    fHistPt[ih]          = nullptr;
    fHistEta[ih]         = nullptr;
    fHistPtVsV0MMult[ih] = nullptr;
    fHistPtVsSPDMult[ih] = nullptr;
    fHistEtaVsSPDMult[ih]= nullptr;
    fHistYVsSPDMult[ih]  = nullptr;
    fHistPtVsNpart[ih]   = nullptr;
    fHistPtVsB[ih]       = nullptr;
    fHistPtVsNMPI[ih]   = nullptr;
    fHistDDYield[ih]   = nullptr;
    fHistDDPt[ih]   = nullptr;
  }
  DefineOutput(1, TList::Class()); // Event Counter Histo
  DefineOutput(2, TTree::Class()); // Event Counter Histo
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
  Int_t lNPtBins = 250;
  Double_t lMaxPt = 25.0;
  
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
  TString lPartNames[78] = {
    "PiPlus", "PiMinus", "KPlus", "KMinus", "Proton", "AntiProton",
    "K0Short", "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
    "Phi", "KStar", "AntiKStar", "LambdaStar", "AntiLambdaStar",
    "D0", "AntiD0", "DPlus", "DMinus", "D0s", "AntiD0s", "DStarPlus", "DStarMinus",
    "Lambdac", "AntiLambdac", "JPsi",
    "Pi0", "AntiPi0", "Eta", "AntiEta",
    "EtaPrime", "OmegaMeson",
    "Omegac", "Omegacbar", "Xic", "XicBar",
    "Xicc", "Xiccbar", "Omegacc", "Omegaccbar",
    "Omegaccc", "Omegacccbar",
    "Xib", "Xibbar",
    "Omegab", "Omegabbar",
    "Lambdab", "Lambdabbar",
    //Prompt measurements
    "PromptD0", "PromptAntiD0", "PromptDPlus", "PromptDMinus", "PromptD0s", "PromptAntiD0s", "PromptDStarPlus", "PromptDStarMinus",
    "PromptLambdac", "PromptAntiLambdac", "PromptJPsi",
    "PromptOmegac", "PromptOmegacbar", "PromptXic", "PromptXicBar",
    "PromptXicc", "PromptXiccbar", "PromptOmegacc", "PromptOmegaccbar",
    "PromptOmegaccc", "PromptOmegacccbar",
    "PromptXib", "PromptXibbar",
    "PromptOmegab", "PromptOmegabbar",
    "PromptLambdab", "PromptLambdabbar"
  };
  
  //Main Output: Histograms
  
  //Event counter histogram: Multiplicity, Npart, b (if available)
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPt[ih] && fkSpeciesSwitch[ih] ) {
      fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPt[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistEta[ih] && fkSpeciesSwitch[ih] ) {
      fHistEta[ih] = new TH1D(Form("fHistEta_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
      fListHist->Add(fHistEta[ih]);
    }
  }
  
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPtVsV0MMult[ih] && fkSpeciesSwitch[ih] ) {
      fHistPtVsV0MMult[ih] = new TH2D(Form("fHistPtVsV0MMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsV0MMult[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPtVsSPDMult[ih] && fkSpeciesSwitch[ih] ) {
      fHistPtVsSPDMult[ih] = new TH2D(Form("fHistPtVsSPDMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistEtaVsSPDMult[ih] && fkDoRapidityStudy && fkSpeciesSwitch[ih] ) {
      fHistEtaVsSPDMult[ih] = new TH2D(Form("fHistEtaVsSPDMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,1,-10,10);
      fListHist->Add(fHistEtaVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistYVsSPDMult[ih] && fkDoRapidityStudy && fkSpeciesSwitch[ih] ) {
      fHistYVsSPDMult[ih] = new TH2D(Form("fHistYVsSPDMult%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,1,-10,10);
      fListHist->Add(fHistYVsSPDMult[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPtVsNpart[ih] && fkDoNpartStudy && fkSpeciesSwitch[ih] ) {
      fHistPtVsNpart[ih] = new TH2D(Form("fHistPtVsNpart_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",fkNNpartBins,-0.5,fkNNpartBins-0.5,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsNpart[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPtVsB[ih] && fkDoImpactParameterStudy && fkSpeciesSwitch[ih] ) {
      fHistPtVsB[ih] = new TH2D(Form("fHistPtVsB_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",fkNBBins,0,20,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsB[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistPtVsNMPI[ih] && fkDoNMPIStudy && fkSpeciesSwitch[ih] ) {
      fHistPtVsNMPI[ih] = new TH2D(Form("fHistPtVsNMPI_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",50,-0.5,49.5,lNPtBins,0,lMaxPt);
      fListHist->Add(fHistPtVsNMPI[ih]);
    }
  }
  
  //Double-differential study: particle yields
  if(! fHistV0MvsSPD ) {
    fHistV0MvsSPD = new TH2D( "fHistV0MvsSPD", "",
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M,
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistV0MvsSPD);
  }
  if(! fHistDDNch ) {
    fHistDDNch = new TH2D( "fHistDDNch", "",
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M,
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistDDNch);
  }
  if(! fHistDDNMPI ) {
    fHistDDNMPI = new TH2D( "fHistDDNMPI", "",
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M,
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistDDNMPI);
  }
  if(! fHistDDQ2 ) {
    fHistDDQ2 = new TH2D( "fHistDDQ2", "",
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M,
                             lNNchBinsV0M/fkDDRebin,lLowNchBoundV0M,lHighNchBoundV0M);
    //Keeps track of some basics
    fListHist->Add(fHistDDQ2);
  }
  
  Int_t lNNchBinsDDV0M = fkMaxMultDDV0M/fkDDRebin;
  Double_t lLowNchBoundDDV0M  = -0.5;
  Double_t lHighNchBoundDDV0M = -0.5 + ((double)(fkMaxMultDDV0M));
  Int_t lNNchBinsDDSPD = fkMaxMultDDSPD/fkDDRebin;
  Double_t lLowNchBoundDDSPD  = -0.5;
  Double_t lHighNchBoundDDSPD = -0.5 + ((double)(fkMaxMultDDSPD));
  
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistDDYield[ih] && fkSpeciesSwitch[ih]) {
      fHistDDYield[ih] = new TH2D(Form("fHistDDYield_%s",lPartNames[ih].Data()), "",
                                  lNNchBinsDDSPD,lLowNchBoundDDSPD,lHighNchBoundDDSPD,
                                  lNNchBinsDDV0M,lLowNchBoundDDV0M,lHighNchBoundDDV0M);
      fListHist->Add(fHistDDYield[ih]);
    }
  }
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(! fHistDDPt[ih] && fkSpeciesSwitch[ih]) {
      fHistDDPt[ih] = new TH2D(Form("fHistDDPt_%s",lPartNames[ih].Data()), "",
                               lNNchBinsDDSPD,lLowNchBoundDDSPD,lHighNchBoundDDSPD,
                               lNNchBinsDDV0M,lLowNchBoundDDV0M,lHighNchBoundDDV0M);
      fListHist->Add(fHistDDPt[ih]);
    }
  }
  
  if(!fTree){
    fTree = new TTree("fTree","particles of interest");
    fTree->Branch("fPt", &fPt, "fPt/F");
    fTree->Branch("fNParents", &fNParents, "fNParents/I");
    fTree->Branch("fStatus", &fStatus, "fStatus/I");
    fTree->Branch("fStatusMother", &fStatusMother, "fStatusMother[fNParents]/I");
    fTree->Branch("fPID", &fPID, "fPID/I");
    fTree->Branch("fPIDMother", &fPIDMother, "fPIDMother[fNParents]/I");
    fTree->Branch("fV0MMultiplicity", &fV0MMultiplicity, "fV0MMultiplicity/I");
  }
  
  //List of Histograms: Normal
  PostData(1, fListHist);
  PostData(2, fTree);
  
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
  
  //this keeps multiplicity over a wide range
  //(multiple intervals as configured by the fkInterval... vars) 
  Double_t lNchEtaWide  = 0;
  
  Long_t lNchVZEROA = 0;
  Long_t lNchVZEROC = 0;
  Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
  
  Double_t lWideEta = 1.4;
  
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
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
    if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
    if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    
    //Special treatment: multiple intervals 
    for(Int_t ii = 0; ii<fkNIntervals; ii++ )
      if( fkIntervalMinEta[ii] < geta && geta < fkIntervalMaxEta[ii] ) lNchEtaWide+=fkIntervalWeight[ii];
  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------
  
  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
  
  if( lNchEtaWide < fkMinimumMultiplicity+1e-10 ) return;
  
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
  Float_t fMC_AvQ2 = -1;
  
  if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
    AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
    if(fMcPythiaHeader){
      fMC_NMPI = fMcPythiaHeader->GetNMPI();
      
//      AliPythia8 *gAliPythiaObject = AliPythia8::Instance();
//      if(gAliPythiaObject){
//        Double_t lAverageQ2=0;
//        for(Long_t iMPI=0; iMPI<gAliPythiaObject->Pythia8()->info.nMPI(); iMPI++){
//          lAverageQ2 += gAliPythiaObject->Pythia8()->info.pTMPI(iMPI);
//        }
//        if( gAliPythiaObject->Pythia8()->info.nMPI() > 0 )
//          lAverageQ2 /= gAliPythiaObject->Pythia8()->info.nMPI();
//      }
      
      //alternative while AliPythia8 not viable
      fMC_AvQ2 = fMcPythiaHeader->GetPtHard();
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
  if( !fHistEventCounter ) {
    Printf("ERROR: Could not retrieve fHistEventCounter! This will crash!\n");
  }
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
  if(fHistV0MvsSPD)     fHistV0MvsSPD       -> Fill ( lNchEtaWide, lNchVZEROA+lNchVZEROC);
  if(fHistDDNch)        fHistDDNch          -> Fill ( lNchEtaWide, lNchVZEROA+lNchVZEROC, lNchEta5);
  if(fHistDDNMPI)       fHistDDNMPI         -> Fill ( lNchEtaWide, lNchVZEROA+lNchVZEROC, fMC_NMPI);
  if(fHistDDQ2)         fHistDDQ2           -> Fill ( lNchEtaWide, lNchVZEROA+lNchVZEROC, fMC_AvQ2);
  
  //Save V0M multiplicity, please
  fV0MMultiplicity = lNchVZEROA+lNchVZEROC;
  
  //------------------------------------------------
  // Fill Spectra as Needed
  //------------------------------------------------
  
  //~All relevant PWG-LF Identified Particle Information (for looping)
  Int_t lPDGCodes[78] = {
    211, -211, 321, -321, 2212, -2212,
    310, 3122, -3122,
    3312, -3312, 3334, -3334,
    333, 313, -313, 3124, -3124,
    421, -421, 411, -411, 431, -431, 413, -413,
    4122, -4122, 443,
    111,-111,221,-221,
    331, 223,
    4332, -4332, 4232, -4232,
    4422, -4422, 4432, -4432,
    4444, -4444,
    5132, -5132,
    5332, -5332,
    5122, -5122,
    421, -421, 411, -411, 431, -431, 413, -413,
    4122, -4122, 443,
    4332, -4332, 4232, -4232,
    4422, -4422, 4432, -4432,
    4444, -4444,
    5132, -5132,
    5332, -5332,
    5122, -5122
  };
  TString lPartNames[78] = {
    "PiPlus", "PiMinus", "KaPlus", "KaMinus", "Proton", "AntiProton",
    "K0Short", "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
    "Phi", "KStar", "AntiKStar", "LambdaStar", "AntiLambdaStar",
    "D0", "AntiD0", "DPlus", "DMinus", "D0s", "AntiD0s", "DStarPlus", "DStarMinus",
    "Lambdac", "AntiLambdac", "JPsi",
    "Pi0", "AntiPi0", "Eta", "AntiEta",
    "EtaPrime", "OmegaMeson",
    "Omegac", "Omegacbar", "Xic", "Xicbar",
    "Xicc", "Xiccbar", "Omegacc", "Omegaccbar",
    "Omegaccc", "Omegacccbar",
    "Xib", "Xibbar",
    "Omegab", "Omegabbar",
    "Lambdab", "Lambdabbar",
    //Prompt measurements
    "PromptD0", "PromptAntiD0", "PromptDPlus", "PromptDMinus", "PromptD0s", "PromptAntiD0s", "PromptDStarPlus", "PromptDStarMinus",
    "PromptLambdac", "PromptAntiLambdac", "PromptJPsi",
    "PromptOmegac", "PromptOmegacbar", "PromptXic", "PromptXicbar",
    "PromptXicc", "PromptXiccbar", "PromptOmegacc", "PromptOmegaccbar",
    "PromptOmegaccc", "PromptOmegacccbar",
    "PromptXib", "PromptXibbar",
    "PromptOmegab", "PromptOmegabbar",
    "PromptLambdab", "PromptLambdabbar"
  };
  Bool_t lCheckIsPhysicalPrimary[78] = {
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE
  };
  Bool_t lCheckHFFeeddown[78] = {
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE
  };
  
  Int_t lThisPDG  = 0;
  Double_t lThisRap  = 0;
  Double_t lThisPt   = 0;
  Bool_t lIsPhysicalPrimary = kFALSE;
  
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
    for(Int_t ih=0; ih<78; ih++) if( lThisPDG == lPDGCodes[ih] ) lContinue = kFALSE;
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
    
    for(Int_t ih=0; ih<78; ih++){
      if( lThisPDG == lPDGCodes[ih] ) {
        //Check if primary (if needed) and if not don't use this particle
        if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
        if( lCheckHFFeeddown[ih] == kTRUE && AliVertexingHFUtils::CheckOrigin(lMCevent, lMCPart, fCheckOriginThirdArgument)!=4 ) continue;
        //Fill Histograms
        if(fHistEta[ih] && lDistanceFromZero > 1e-12 ) fHistEta[ih] -> Fill ( lPart -> Eta() );

        if(fHistEtaVsSPDMult[ih]) fHistEtaVsSPDMult[ih] -> Fill( lNchEtaWide, lPart -> Eta() );
        if(fHistYVsSPDMult[ih]) fHistYVsSPDMult[ih] -> Fill( lNchEtaWide, lThisRap );
        if( TMath::Abs(lThisRap) < 0.5 && !fkWideRapiditySpeciesStudy ) {
          if( fHistPt[ih] ) fHistPt[ih]->Fill(lThisPt);
          if( fHistPtVsV0MMult[ih] ) fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
          if( fHistPtVsSPDMult[ih] ) fHistPtVsSPDMult[ih]->Fill(lNchEtaWide,lThisPt);
          if( fHistPtVsNpart[ih] ) fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
          if( fHistPtVsB[ih] ) fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
          if( fHistPtVsNMPI[ih] ) fHistPtVsNMPI[ih]->Fill(fMC_NMPI,lThisPt);
          if( fHistDDYield[ih] ) fHistDDYield[ih]->Fill(lNchEtaWide,lNchVZEROA+lNchVZEROC);
          if( fHistDDPt[ih] ) fHistDDPt[ih]->Fill(lNchEtaWide,lNchVZEROA+lNchVZEROC, lThisPt);
        }
        if( TMath::Abs(lThisRap) < 4.0 && fkWideRapiditySpeciesStudy ) {
          if( fHistPt[ih] ) fHistPt[ih]->Fill(lThisPt);
          if( fHistPtVsV0MMult[ih] ) fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
          if( fHistPtVsSPDMult[ih] ) fHistPtVsSPDMult[ih]->Fill(lNchEtaWide,lThisPt);
          if( fHistPtVsNpart[ih] ) fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
          if( fHistPtVsB[ih] ) fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
          if( fHistPtVsNMPI[ih] ) fHistPtVsNMPI[ih]->Fill(fMC_NMPI,lThisPt);
          if( fHistDDYield[ih] ) fHistDDYield[ih]->Fill(lNchEtaWide,lNchVZEROA+lNchVZEROC);
          if( fHistDDPt[ih] ) fHistDDPt[ih]->Fill(lNchEtaWide,lNchVZEROA+lNchVZEROC, lThisPt);
        }
      }
    }
    
    if(fkSaveTree){
      //Store D0, D+, D- separately, please; PDG = 421, -421, 411, -411
      fNParents = 0;
      for(Int_t ii=0; ii<10; ii++) fPIDMother[ii]=-1;
      if(TMath::Abs(lThisRap)<0.5 && (TMath::Abs(lThisPDG)==421 || TMath::Abs(lThisPDG)==411) ){
        fPt = lThisPt;
        fPID = lThisPDG;
        fStatus = lMCPart->MCStatusCode();
        Int_t lblMother = lMCPart->GetMother();
        while(lblMother>=0){
          AliMCParticle* lParentParticle = (AliMCParticle*)lMCevent->GetTrack( lblMother );
          if(!lParentParticle) break;
          fPIDMother[fNParents] = lParentParticle->PdgCode();
          fStatusMother[fNParents] = lParentParticle->MCStatusCode();
          fNParents++;
          if(fNParents>=10) break;
          lblMother = lParentParticle->GetMother();
        }
        fTree->Fill();
      }
    }
  }//End of loop on tracks
  //----- End Loop on Stack ----------------------
  
  // Post output data.
  PostData(1, fListHist);
  PostData(2, fTree);
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

//______________________________________________________________________
void AliAnalysisTaskMCPredictions::PrintEtaIntervals(){
  for(Int_t ii=0; ii<fkNIntervals; ii++)
    cout<<"Interval #"<<ii<<"\tmin = "<<fkIntervalMinEta[ii]<<"\tmax = "<<fkIntervalMaxEta[ii]<<"\tweight = "<<fkIntervalWeight[ii]<<endl;
}

//______________________________________________________________________
void AliAnalysisTaskMCPredictions::SetStandardSpecies(){
  Bool_t lSaveThisSpecies[78] = {
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE
  };
  for(Int_t ii=0; ii<78; ii++){
    fkSpeciesSwitch[ii] = lSaveThisSpecies[ii]; 
  }
}
