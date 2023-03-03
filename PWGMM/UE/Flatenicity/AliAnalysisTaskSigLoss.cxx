/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
using std::cout;
using std::endl;


#include "AliAnalysisTaskSigLoss.h"

TString fgSpcName[AliAnalysisTaskSigLoss::kNspc] = {"Pi", "Ka", "Pr", "K0", "Ks", "Phi", "Lambda", "Xi", "Om"};

using namespace std;

ClassImp(AliAnalysisTaskSigLoss) // classimp: necessary for root

//___________________________________________________________________________________
AliAnalysisTaskSigLoss::AliAnalysisTaskSigLoss():
  AliAnalysisTaskSE(),
  fESD(0x0),
  fOutput(0x0),
  fEventCuts(0x0),
  fMCstack(0), 
  fMC(0), 
  fHistNEvents(0x0),
  fHistMCEvents(0x0),
  fHistVtxZ(0x0),
  fHistPrimMCGenVtxZallCh(0x0),
  fHistPrimMCGenVtxZcutCh(0x0),
  fMultEstimator("V0M"),
  fDoMultSel(kFALSE),
  fHistMultBeforeEvtSel(0x0),
  fHistMultAfterEvtSel(0x0),
  fLowMult(-1),
  fHighMult(-1),
  fEvtMult(-999),
  fTriggerSel(AliVEvent::kINT7),
  fMaxVtxZCut(10),
  fMinRapCut(-.5),
  fMaxRapCut( .5),
  fCMSRapFct(.0),
  fPlpType(kNoPileup),
  fMinPlpContribSPD(5),
  fMinPlpZdistSPD(.8)
{
  for (Int_t iSpc = 0; iSpc < kNspc; ++iSpc) {
    fHistPrimMCGenVtxZall[iSpc] = 0x0;
    fHistPrimMCGenVtxZcut[iSpc] = 0x0;
  }
}

//___________________________________________________________________________________
AliAnalysisTaskSigLoss::AliAnalysisTaskSigLoss(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fOutput(0x0),
  fEventCuts(0x0),
  fMCstack(0), 
  fMC(0), 
  fHistNEvents(0x0),
  fHistMCEvents(0x0),
  fHistVtxZ(0x0),
  fHistPrimMCGenVtxZallCh(0x0),
  fHistPrimMCGenVtxZcutCh(0x0),
  fMultEstimator("V0M"),
  fDoMultSel(kFALSE),
  fHistMultBeforeEvtSel(0x0),
  fHistMultAfterEvtSel(0x0),
  fLowMult(-1),
  fHighMult(-1),
  fEvtMult(-999),
  fTriggerSel(AliVEvent::kINT7),
  fMaxVtxZCut(10),
  fMinRapCut(-.5),
  fMaxRapCut( .5),
  fCMSRapFct(.0),
  fPlpType(kNoPileup),
  fMinPlpContribSPD(5),
  fMinPlpZdistSPD(.8)
{
  for (Int_t iSpc = 0; iSpc < kNspc; ++iSpc) {
    fHistPrimMCGenVtxZall[iSpc] = 0x0;
    fHistPrimMCGenVtxZcut[iSpc] = 0x0;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//___________________________________________________________________________________
AliAnalysisTaskSigLoss::~AliAnalysisTaskSigLoss()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0x0;
  }
}

//___________________________________________________________________________________
void AliAnalysisTaskSigLoss::UserCreateOutputObjects()
{
  const  double fPtBinLimits[kNbin+1] = {
    0.0, 0.1, 0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.6, 0.7, 0.8,
    0.9, 1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 30.0, 40.0, 50.0};

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("SigLoss");

  fEventCuts.AddQAplotsToList(fOutput); //Histograms for event Selection

  TString plpName[2] = {"Sel", "SPD"};
  fHistNEvents = new TH1I("fHistNEvents", "Number of processed events;Ev. Sel. Step;Counts", kNEvtCuts, .5, (kNEvtCuts + .5));
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(kIsReadable, "Readable");
  fHistNEvents->GetXaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistNEvents->GetXaxis()->SetBinLabel(kPassMultSel, "PassMultSel");
  fHistNEvents->GetXaxis()->SetBinLabel(kIsNotPileup, Form("IsNotPileup_%s", plpName[fPlpType].Data()));
  fHistNEvents->GetXaxis()->SetBinLabel(kHasRecVtx, "HasVertex");
  fHistNEvents->GetXaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fOutput->Add(fHistNEvents);

  fHistMCEvents = new TH1I("fHistMCEvents", "Number of processed events;Ev. Sel. Step;Counts", kNEvtCuts, .5, (kNEvtCuts + .5));
  fHistMCEvents->Sumw2();
  fHistMCEvents->SetMinimum(0);
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsReadable, "Readable");
  fHistMCEvents->GetXaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistMCEvents->GetXaxis()->SetBinLabel(kPassMultSel, "PassMultSel");
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsNotPileup, Form("IsNotPileup_%s", plpName[fPlpType].Data()));
  fHistMCEvents->GetXaxis()->SetBinLabel(kHasRecVtx, "HasVertex");
  fHistMCEvents->GetXaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fOutput->Add(fHistMCEvents);

  Int_t kNmultBins = 115;
  Float_t lMultBinLimit[kNmultBins];
  for (Int_t ibin = 0; ibin < kNmultBins; ++ibin)
    lMultBinLimit[ibin] = (ibin < 100) ? (0. + (1.) * ibin) : (100. + 100 * (ibin - 100));

  fHistMultBeforeEvtSel = new TH1F("fHistMultBeforeEvtSel", "Event Multiplicity before event selection", kNmultBins - 1, lMultBinLimit);
  fHistMultBeforeEvtSel->Sumw2();
  fOutput->Add(fHistMultBeforeEvtSel);

  fHistMultAfterEvtSel = new TH1F("fHistMultAfterEvtSel", "Event Multiplicity after event selection", kNmultBins - 1, lMultBinLimit);
  fHistMultAfterEvtSel->Sumw2();
  fOutput->Add(fHistMultAfterEvtSel);

  fHistVtxZ = new TH1F("fHistVtxZ", "Vtx Z distribution", 400, -20., 20.);
  fHistVtxZ->Sumw2();
  fHistVtxZ->SetMinimum(0);
  fOutput->Add(fHistVtxZ);

  TString histName("");
  for (Int_t ispc = 0; ispc < kNspc; ++ispc) {
    //Histograms MC part Gen before and after all selection Good Vertex Gen
    histName = Form("fHistPrimMCGenVtxZall%s", fgSpcName[ispc].Data());
    fHistPrimMCGenVtxZall[ispc] = new TH2F(histName, histName, kNEvtCuts, .5, (kNEvtCuts + .5), kNbin, fPtBinLimits);

    histName = Form("fHistPrimMCGenVtxZcut%s", fgSpcName[ispc].Data());
    fHistPrimMCGenVtxZcut[ispc] = new TH2F(histName, histName, kNEvtCuts, .5, (kNEvtCuts + .5), kNbin, fPtBinLimits);

    fOutput->Add(fHistPrimMCGenVtxZall[ispc]);
    fOutput->Add(fHistPrimMCGenVtxZcut[ispc]);
  }

  fHistPrimMCGenVtxZallCh = new TH2F("fHistPrimMCGenVtxZallCh", "fHistPrimMCGenVtxZallCh", kNEvtCuts, .5, (kNEvtCuts + .5), kNbin, fPtBinLimits);
  fOutput->Add(fHistPrimMCGenVtxZallCh);
  fHistPrimMCGenVtxZcutCh = new TH2F("fHistPrimMCGenVtxZcutCh", "fHistPrimMCGenVtxZcutCh", kNEvtCuts, .5, (kNEvtCuts + .5), kNbin, fPtBinLimits);
  fOutput->Add(fHistPrimMCGenVtxZcutCh);
  
  PostData(1, fOutput);
}

//___________________________________________________________________________________
void AliAnalysisTaskSigLoss::UserExec(Option_t*)
{
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    AliDebug(3, "fESD not available");
    PostData(1, fOutput);
    return;
  }

  Bool_t lHasGoodVtxGen = kFALSE;
  AliMCEventHandler* lMCevtHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!lMCevtHandler) {
    AliDebug(3, "Could not retrieve MC event handler");
    PostData(1, fOutput);
    return;
  }

  fMC = lMCevtHandler->MCEvent();
  if (!fMC) {
    AliDebug(3, "Could not retrieve MC event");
    PostData(1, fOutput);
    return;
  }

  fMCstack = fMC->Stack();
  if (!fMCstack) {
    AliDebug(3, "MC stack not available");
    PostData(1, fOutput);
    return;
  }

  const AliVVertex* lGenVtx = fMC->GetPrimaryVertex();
  if (!lGenVtx) {
    AliDebug(3, "MC Vtx not available");
    PostData(1, fOutput);
    return;
  }

  // Vertex z-cut
  if (!(TMath::Abs(lGenVtx->GetZ()) > fMaxVtxZCut)) {
    lHasGoodVtxGen = kTRUE;
  }

  //Event selection
  EEvtCut_Type lastEvtCutPassed = kIsReadable;
  Bool_t lIsEventSelected = IsEventAccepted(lastEvtCutPassed);

  DoMCPart(fMCstack, lastEvtCutPassed, lHasGoodVtxGen);
  
  for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep) {
    fHistNEvents->Fill(istep);
    if (lHasGoodVtxGen) 
        fHistMCEvents->Fill(istep);
  }

  if (!lIsEventSelected)  {
    AliDebugF(3,"Event rejected by event selection after %d step",(int)lastEvtCutPassed);
    PostData(1, fOutput);
    return;
  }
  
  if (fDoMultSel)
    fHistMultAfterEvtSel->Fill(fEvtMult);

  PostData(1, fOutput);
  return;
}

//___________________________________________________________________________________
void AliAnalysisTaskSigLoss::Terminate(Option_t*)
{
  return;
}

//___________________________________________________________________________________
Bool_t AliAnalysisTaskSigLoss::IsEventAccepted(EEvtCut_Type& evtSel)
{
  evtSel = kIsReadable;

  // Sel trigger
  //________________________________________
  UInt_t maskPhysSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  maskPhysSel &= fTriggerSel;
  if (!maskPhysSel) {
    AliDebugF(3, "No event passing physics evt. sel. for trigger %d", fTriggerSel);
    PostData(1, fOutput);
    return kFALSE;
  }
  evtSel = kPassTrig;
  
  // Sel mult
  //________________________________________
  if (fDoMultSel && IsMultSelected()) {
    if ((fLowMult > -1 && fEvtMult < fLowMult) || (fHighMult > -1 && fEvtMult > fHighMult)) {
      AliDebugF(3, "Event with multiplicity = %.2f outside range [%.2f,%.2f]", fEvtMult, fLowMult, fHighMult);
      PostData(1, fOutput);
      return kFALSE;
    }
    fHistMultBeforeEvtSel->Fill(fEvtMult);
  }
  evtSel = kPassMultSel;

  // Sel good events
  //________________________________________
  if (fEventCuts.AcceptEvent(fESD)) {
      evtSel = kHasGoodVtxZ;
      return kTRUE;
  }
  AliDebug(3, "Event not accepted by AliEventCuts");

  // Sel rec vertex
  //________________________________________
  if (!HasRecVertex()) {
      AliDebug(3, "No event vtx quality sel.");
      PostData(1, fOutput);
      return kFALSE;
  }
  evtSel = kHasRecVtx;

  // Sel good vertex Z
  //________________________________________
  if (!IsGoodVtxZ()) {
      AliDebugF(3, "Vertex with Z>%f cm", fMaxVtxZCut);
      PostData(1, fOutput);
      return kFALSE;
  }
  evtSel = kHasGoodVtxZ;

  // Sel PileUp
  //________________________________________
  if (IsPileup()) {
      AliDebug(3, "Event with PileUp");
      PostData(1, fOutput);
      return kFALSE;
  }
  evtSel = kIsNotPileup;
  
  return kTRUE;
}

//___________________________________________________________________________________
Bool_t AliAnalysisTaskSigLoss::IsMultSelected()
{
  fEvtMult = -999.;
  AliMultSelection* fMultSel = (AliMultSelection*) fESD->FindListObject("MultSelection");
  if (!fMultSel) {
      AliWarning("AliMultSelection object not found!");
      return kFALSE;
  } else {
      //Event selection is embedded in the Multiplicity estimator so that
      //the Multiplicity percentiles are well defined and refer to the same sample
      fEvtMult = fMultSel->GetMultiplicityPercentile(fMultEstimator.Data());
      return kTRUE;
  }
}

//___________________________________________________________________________________
Bool_t AliAnalysisTaskSigLoss::IsPileup()
{
  AliAnalysisUtils utils;

  switch (fPlpType) {
    case kPileupSPD :
      utils.SetMinPlpContribSPD(fMinPlpContribSPD);
      utils.SetMinPlpZdistSPD(fMinPlpZdistSPD);
      utils.SetUseMVPlpSelection(kFALSE);
      break;
    default :
      return kFALSE;
      break;
  }

  return utils.IsPileUpEvent(fESD);
}

//___________________________________________________________________________________
Bool_t AliAnalysisTaskSigLoss::HasRecVertex()
{
  float fMaxDeltaSpdTrackAbsolute = 0.5f;
  float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  float fMaxResolutionSPDvertex = 0.25f;
  float fMaxDispersionSPDvertex = 1.e14f;

  Bool_t fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag = BIT(AliEventCuts::kNoCuts);

  const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
  Bool_t isTrackV = true;
  if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
    isTrackV = false;
  const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

  if (vtSPD->GetNContributors() > 0)
    fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV)
    fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
      (fFlag & BIT(AliEventCuts::kVertexSPD)))
    fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex *&vtx =
      Bool_t(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
  if (!fPrimaryVertex)
    return kFALSE;

  /// Vertex quality cuts
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = Bool_t(fFlag & AliEventCuts::kVertexSPD) &&
                      Bool_t(fFlag & AliEventCuts::kVertexTracks)
                  ? vtTrc->GetZ() - vtSPD->GetZ()
                  : 0.; /// If one of the two vertices is not available this cut
                        /// is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      Bool_t(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <=
           fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() ||
       TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() ||
       vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
                                                   /// run1, only for ESD
      ) // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);

  Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
                  (TESTBIT(fFlag, AliEventCuts::kVertexQuality));

  return hasVtx;
}

//___________________________________________________________________________________
Bool_t AliAnalysisTaskSigLoss::IsGoodVtxZ()
{
  const AliESDVertex* vtx = fESD->GetPrimaryVertex();
  fHistVtxZ->Fill(vtx->GetZ());
  if (TMath::Abs(vtx->GetZ()) > fMaxVtxZCut) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSigLoss::DoMCPart(AliStack* fMCstack, EEvtCut_Type lastEvtCutPassed, Bool_t lHasGoodVtxGen)
{
  int nTrackMC = (fMCstack) ? fMCstack->GetNtrack() : 0;
  for (int i_mcTrk = 0; i_mcTrk < nTrackMC; ++i_mcTrk) {
    TParticle* mcTrk = fMCstack->Particle(i_mcTrk);
    if (!mcTrk) {
        printf("----ERROR: mcTrk not available------------------\n");
        continue;
    }
    
    int mcPDG = mcTrk->GetPdgCode();

    Int_t iPart = -1;
    if      (TMath::Abs(mcPDG) == 211)  iPart = 0; // pi+-
    else if (TMath::Abs(mcPDG) == 321)  iPart = 1; // K+-
    else if (TMath::Abs(mcPDG) == 2212) iPart = 2; // p/pbar
    else if (TMath::Abs(mcPDG) == 310)  iPart = 3; // K0
    else if (TMath::Abs(mcPDG) == 313)  iPart = 4; // Kstar
    else if (TMath::Abs(mcPDG) == 333)  iPart = 5; // Phi
    else if (TMath::Abs(mcPDG) == 3122) iPart = 6; // Lambda
    else if (TMath::Abs(mcPDG) == 3312) iPart = 7; // Xi
    else if (TMath::Abs(mcPDG) == 3334) iPart = 8; // Omega
    else  continue;

    Double_t mcPt   = mcTrk->Pt();
    Double_t mcEta  = mcTrk->Eta();
    Double_t mcMass = mcTrk->GetMass();
//     printf("\t ===> \t McMass = %.4f\n",mcMass);
// //     Double_t mcRap = EtaToY(mcPt, mcMass, mcEta) + fCMSRapFct;
//     printf("\t ===> \t mcRap = %.4f\n",mcRap);
// //     if (!IsRapIn(mcRap)) 
// //         continue;
    
    if (TMath::Abs(mcEta) > 1.)
        continue;    

    Bool_t lIsPhysPrimary = (iPart == 4 || iPart == 5) ? kTRUE : fMCstack->IsPhysicalPrimary(i_mcTrk);
    if (!lIsPhysPrimary) 
        continue;

    for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep) {
      fHistPrimMCGenVtxZallCh->Fill(istep, TMath::Abs(mcPt));
      fHistPrimMCGenVtxZall[iPart]->Fill(istep, TMath::Abs(mcPt));
      
      if (lHasGoodVtxGen){
        fHistPrimMCGenVtxZcutCh->Fill(istep, TMath::Abs(mcPt));
        fHistPrimMCGenVtxZcut[iPart]->Fill(istep, TMath::Abs(mcPt));
      }
    }
  }
}

//___________________________________________________________________________________
void AliAnalysisTaskSigLoss::SetupStandardEventCutsForRun2()
{
  fDoMultSel        = kTRUE;
  fPlpType          = kPileupSPD;
  fMinPlpContribSPD = 3;
  fMinPlpZdistSPD   = .8;
  fTriggerSel       = AliVEvent::kINT7;
  fMaxVtxZCut       = 10.;
}

