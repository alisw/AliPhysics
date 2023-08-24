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
// This task is meant to acquire MC-level predictions for strangeness studies vs multiplicity
// and ZDC Energy.
//
// Modified version of AliAnalysisTaskMCPredictions.cxx from David Dobrigkeit Chinellato
//
// --- Francesca Ercolessi
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
#include "TProfile2D.h"
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
#include "AliAnalysisTaskMCPredictionsStrgVsMultVsZDC.h"
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
#include "AliVVertex.h"
//#include "AliPythia8.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictionsStrgVsMultVsZDC)

AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::AliAnalysisTaskMCPredictionsStrgVsMultVsZDC()
: AliAnalysisTaskSE(),
fkSelectINELgtZERO(kTRUE),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistV0AMult(0),
fHistV0CMult(0),
fHistMult05(0),
fHistMult08(0),
fHistMult10(0),
fHistMult08to15(0),
fHistSPDClusters(0),
fHistSPDCl0(0),
fHistSPDCl1(0),
fHistNMPI(0),
fHistQ2(0),
fHistb(0),
fHistLeadingE(0),
fHistEffEnergy(0),
fHistRxy(0),
f2DHistINELgt0SPDV0M(0),
f2DHistLeadingESPDV0M(0),
f2DHistLeadingERecoPercSPDV0M(0),
f2DHistEffEnergySPDV0M(0),
f2DHistNchSPDV0M(0),
f2DHistNchRecoPercSPDV0M(0),
f2DHistINELgt0RecoPercSPDV0M(0),
f2DHistNMPISPDV0M(0),
f2DHistQ2SPDV0M(0),
f2DHistbSPDV0M(0),
f2DHistINELgt0Nch0815V0M(0),
f2DHistLeadingENch0815V0M(0),
f2DHistEffEnergyNch0815V0M(0),
f2DHistNchNch0815V0M(0),
f2DHistNMPINch0815V0M(0),
f2DHistINELgt0Nch0815ZDC(0),
f2DHistLeadingENch0815ZDC(0),
f2DHistEffEnergyNch0815ZDC(0),
f2DHistNchNch0815ZDC(0),
f2DHistNMPINch0815ZDC(0),
f2dHistZDCVsLE(0),
f2dHistZDCVsEE(0),
f2dHistZDCVsLEA(0),
f2dHistZDCVsLEC(0),
f2dHistZPVsLP(0),
f2dHistZNVsLN(0),
f2dHistZDCVsLEnoacc(0),
f2dHistZPVsLPnoacc(0),
f2dHistZNVsLNnoacc(0),
f2dHistSPDClRecoVsTrue(0),
f2dHistV0MRecoVsTrue(0),
f2dHistTrueVsRecoSPDCl(0),
f2dHistTrueVsRecoSPDCl0(0),
f2dHistTrueVsRecoSPDCl1(0),
f3dHistPi0SPDMultSPDCl(0),
p2dHistPi0SPDMultSPDCl(0),
f2dHistTrueVsRecoNch0815(0)
{
  for(Int_t ih=0; ih<22; ih++){
    fHistPt[ih] = 0x0;
    fHistDecayVtxPos[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
    f2DHistAvPtSPDV0M[ih] = 0x0;
    f2DHistPartNch0815V0M[ih] = 0x0;
    f2DHistPartRecoPercSPDV0M[ih] = 0x0;
    f2DHistPartNch0815ZDC[ih] = 0x0;
  }
}

AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const char *name, Float_t lCenterOfMassEnergy, Bool_t kDoPythia, Bool_t kDoEPOS)
: AliAnalysisTaskSE(name),
fCenterOfMassEnergy(lCenterOfMassEnergy),
fkSelectINELgtZERO(kTRUE),
fkDoPythia(kDoPythia),
fkDoEPOS(kDoEPOS),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistV0AMult(0),
fHistV0CMult(0),
fHistMult05(0),
fHistMult08(0),
fHistMult10(0),
fHistMult08to15(0),
fHistSPDClusters(0),
fHistSPDCl0(0),
fHistSPDCl1(0),
fHistNMPI(0),
fHistQ2(0),
fHistb(0),
fHistLeadingE(0),
fHistEffEnergy(0),
fHistRxy(0),
f2DHistINELgt0SPDV0M(0),
f2DHistLeadingESPDV0M(0),
f2DHistLeadingERecoPercSPDV0M(0),
f2DHistEffEnergySPDV0M(0),
f2DHistNchSPDV0M(0),
f2DHistNchRecoPercSPDV0M(0),
f2DHistINELgt0RecoPercSPDV0M(0),
f2DHistNMPISPDV0M(0),
f2DHistQ2SPDV0M(0),
f2DHistbSPDV0M(0),
f2DHistINELgt0Nch0815V0M(0),
f2DHistLeadingENch0815V0M(0),
f2DHistEffEnergyNch0815V0M(0),
f2DHistNchNch0815V0M(0),
f2DHistNMPINch0815V0M(0),
f2DHistINELgt0Nch0815ZDC(0),
f2DHistLeadingENch0815ZDC(0),
f2DHistEffEnergyNch0815ZDC(0),
f2DHistNchNch0815ZDC(0),
f2DHistNMPINch0815ZDC(0),
f2dHistZDCVsLE(0),
f2dHistZDCVsEE(0),
f2dHistZDCVsLEA(0),
f2dHistZDCVsLEC(0),
f2dHistZPVsLP(0),
f2dHistZNVsLN(0),
f2dHistZDCVsLEnoacc(0),
f2dHistZPVsLPnoacc(0),
f2dHistZNVsLNnoacc(0),
f2dHistSPDClRecoVsTrue(0),
f2dHistV0MRecoVsTrue(0),
f2dHistTrueVsRecoSPDCl(0),
f2dHistTrueVsRecoSPDCl0(0),
f2dHistTrueVsRecoSPDCl1(0),
f3dHistPi0SPDMultSPDCl(0),
p2dHistPi0SPDMultSPDCl(0),
f2dHistTrueVsRecoNch0815(0)
{
  for(Int_t ih=0; ih<22; ih++){
    fHistPt[ih] = 0x0;
    fHistDecayVtxPos[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
    f2DHistAvPtSPDV0M[ih] = 0x0;
    f2DHistPartNch0815V0M[ih] = 0x0;
    f2DHistPartRecoPercSPDV0M[ih] = 0x0;
    f2DHistPartNch0815ZDC[ih] = 0x0;
  }
  //
  DefineOutput(1, TList::Class());
}


AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::~AliAnalysisTaskMCPredictionsStrgVsMultVsZDC()
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
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::UserCreateOutputObjects()
{
  //------------------------------------------------
  // Histograms: Basic Analysis Output
  //------------------------------------------------
  // Create histograms
  fListHist = new TList();
  fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

  TString lPartNames[22] = {
    "PiPlus", "PiMinus",
    "KaPlus", "KaMinus",
    "Proton", "AntiProton",
    "K0Short",
    "Lambda", "AntiLambda",
    "XiMinus", "XiPlus",
    "OmegaMinus", "OmegaPlus",
    "Phi",
    "D0", "AntiD0",
    "DPlus", "DMinus",
    "Lambdac", "AntiLambdac",
    "JPsi",
    "Pi0"
  };

  //-----------------------------------------------------------------------------
  if(! fHistEventCounter ) {
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "INEL>0");
    fListHist->Add(fHistEventCounter);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0MMult ) {
    fHistV0MMult = new TH1D( "fHistV0MMult", ";V0M Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0MMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0AMult ) {
    fHistV0AMult = new TH1D( "fHistV0AMult", ";V0A Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0AMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0CMult ) {
    fHistV0CMult = new TH1D( "fHistV0CMult", ";V0C Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0CMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult05 ) {
    fHistMult05 = new TH1D( "fHistMult05", ";Nch in |#eta|<0.5 ;Count",1000,0,1000);
    fListHist->Add(fHistMult05);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult08 ) {
    fHistMult08 = new TH1D( "fHistMult08", ";Nch in |#eta|<0.8 ;Count",1000,0,1000);
    fListHist->Add(fHistMult08);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult10 ) {
    fHistMult10 = new TH1D( "fHistMult10", ";Nch in |#eta|<1.0 ;Count",1000,0,1000);
    fListHist->Add(fHistMult10);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult08to15 ) {
    fHistMult08to15 = new TH1D( "fHistMult08to15", ";Nch in 0.8<|#eta|<1.5 ;Count",1000,0,1000);
    fListHist->Add(fHistMult08to15);
  }

  //-----------------------------------------------------------------------------
  if(! fHistSPDClusters ) {
    fHistSPDClusters = new TH1D( "fHistSPDClusters", ";SPD Clusters ;Count",1000,0,1000);
    fListHist->Add(fHistSPDClusters);
  }

  //-----------------------------------------------------------------------------
  if(! fHistSPDCl0 ) {
    fHistSPDCl0 = new TH1D( "fHistSPDCl0", ";SPD Clusters ;Count",1000,0,1000);
    fListHist->Add(fHistSPDCl0);
  }

  //-----------------------------------------------------------------------------
  if(! fHistSPDCl1 ) {
    fHistSPDCl1 = new TH1D( "fHistSPDCl1", ";SPD Clusters ;Count",1000,0,1000);
    fListHist->Add(fHistSPDCl1);
  }

  //-----------------------------------------------------------------------------
  if(!fHistNMPI && fkDoPythia) {
    fHistNMPI = new TH1D( "fHistNMPI", ";NMPI;Count",50,0,50);
    fListHist->Add(fHistNMPI);
  }

  //-----------------------------------------------------------------------------
  if(!fHistQ2 && fkDoPythia) {
    fHistQ2 = new TH1D( "fHistQ2", "; Q^{2} ;Count",1000,0,1000);
    fListHist->Add(fHistQ2);
  }

  //-----------------------------------------------------------------------------
  if(!fHistb && fkDoEPOS) {
    fHistb = new TH1D( "fHistb", ";b;Count",20,0,20);
    fListHist->Add(fHistb);
  }

  //-----------------------------------------------------------------------------
  if(!fHistLeadingE ) {
    fHistLeadingE = new TH1D( "fHistLeadingE", ";Leading energy (GeV);Count",1300,0,fCenterOfMassEnergy);
    fListHist->Add(fHistLeadingE);
  }

  //-----------------------------------------------------------------------------
  if(!fHistEffEnergy ) {
    fHistEffEnergy = new TH1D( "fHistEffEnergy", ";Effective energy (GeV);Count",1300,0,fCenterOfMassEnergy);
    fListHist->Add(fHistEffEnergy);
  }

  //-----------------------------------------------------------------------------
  if(!fHistRxy ) {
    fHistRxy = new TH1D( "fHistRxy", ";Rxy (?);Count",200,0,50);
    fListHist->Add(fHistRxy);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistINELgt0SPDV0M) {
    f2DHistINELgt0SPDV0M = new TH2D( "f2DHistINELgt0SPDV0M", "INEL>0 events;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistINELgt0SPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistINELgt0RecoPercSPDV0M) {
    f2DHistINELgt0RecoPercSPDV0M = new TH2D( "f2DHistINELgt0RecoPercSPDV0M", "INEL>0 events;SPD percentile; V0M percentile", 100, 0, 100., 100, 0, 100.);
    fListHist->Add(f2DHistINELgt0RecoPercSPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistLeadingESPDV0M) {
    f2DHistLeadingESPDV0M = new TH2D( "f2DHistLeadingESPDV0M", "Leading energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistLeadingESPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistLeadingERecoPercSPDV0M) {
    f2DHistLeadingERecoPercSPDV0M = new TH2D( "f2DHistLeadingERecoPercSPDV0M", "Leading energy (GeV);SPD percentile; V0M percentile", 100, 0, 100., 100, 0, 100.);
    fListHist->Add(f2DHistLeadingERecoPercSPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistEffEnergySPDV0M) {
    f2DHistEffEnergySPDV0M = new TH2D( "f2DHistEffEnergySPDV0M", "Effective energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistEffEnergySPDV0M);
  }


  //-----------------------------------------------------------------------------
  if(!f2DHistNchSPDV0M) {
    f2DHistNchSPDV0M = new TH2D( "f2DHistNchSPDV0M", "Nch (|#eta|<0.5);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistNchSPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistNchRecoPercSPDV0M) {
    f2DHistNchRecoPercSPDV0M = new TH2D( "f2DHistNchRecoPercSPDV0M", "Nch (|#eta|<0.5);SPD percentile; V0M percentile", 100, 0, 100., 100, 0, 100.);
    fListHist->Add(f2DHistNchRecoPercSPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistNMPISPDV0M && fkDoPythia) {
    f2DHistNMPISPDV0M = new TH2D( "f2DHistNMPISPDV0M", "NMPI;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistNMPISPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistQ2SPDV0M && fkDoPythia) {
    f2DHistQ2SPDV0M = new TH2D( "f2DHistQ2SPDV0M", "Q2;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistQ2SPDV0M);
  }


  //-----------------------------------------------------------------------------
  if(!f2DHistbSPDV0M && fkDoEPOS) {
    f2DHistbSPDV0M = new TH2D( "f2DHistbSPDV0M", "b;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistbSPDV0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistINELgt0Nch0815V0M) {
    f2DHistINELgt0Nch0815V0M = new TH2D( "f2DHistINELgt0Nch0815V0M", "INEL>0 events;Nch 0.8 < |#eta| < 1.5; V0M multiplicity", 200, 0, 200., 500, 0, 500.);
    fListHist->Add(f2DHistINELgt0Nch0815V0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistLeadingENch0815V0M) {
    f2DHistLeadingENch0815V0M = new TH2D( "f2DHistLeadingENch0815V0M", "Leading energy (GeV);Nch 0.8 < |#eta| < 1.5; V0M multiplicity", 200, 0, 200., 500, 0, 500.);
    fListHist->Add(f2DHistLeadingENch0815V0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistEffEnergyNch0815V0M) {
    f2DHistEffEnergyNch0815V0M = new TH2D( "f2DHistEffEnergyNch0815V0M", "Effective energy (GeV);Nch 0.8 < |#eta| < 1.5; V0M multiplicity", 200, 0, 200., 500, 0, 500.);
    fListHist->Add(f2DHistEffEnergyNch0815V0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistNchNch0815V0M) {
    f2DHistNchNch0815V0M = new TH2D( "f2DHistNchNch0815V0M", "Nch (|#eta|<0.5);Nch 0.8 < |#eta| < 1.5; V0M multiplicity", 200, 0, 200., 500, 0, 500.);
    fListHist->Add(f2DHistNchNch0815V0M);
  }

  //-----------------------------------------------------------------------------
  if(!f2DHistNMPINch0815V0M && fkDoPythia) {
    f2DHistNMPINch0815V0M = new TH2D( "f2DHistNMPINch0815V0M", "NMPI;Nch 0.8 < |#eta| < 1.5; V0M multiplicity", 200, 0, 200., 500, 0, 500.);
    fListHist->Add(f2DHistNMPINch0815V0M);
  }

  //-----------------------------------------------------------------------------
  if (!f2DHistINELgt0Nch0815ZDC)
  {
    f2DHistINELgt0Nch0815ZDC = new TH2D("f2DHistINELgt0Nch0815ZDC", "INEL>0 events;Nch 0.8 < |#eta| < 1.5; Leading energy (within ZP and ZN acc)", 200, 0, 200., 1300, 0, fCenterOfMassEnergy);
    fListHist->Add(f2DHistINELgt0Nch0815ZDC);
  }

  //-----------------------------------------------------------------------------
  if (!f2DHistLeadingENch0815ZDC)
  {
    f2DHistLeadingENch0815ZDC = new TH2D("f2DHistLeadingENch0815ZDC", "Leading energy (GeV);Nch 0.8 < |#eta| < 1.5; Leading energy (within ZP and ZN acc)", 200, 0, 200., 1300, 0, fCenterOfMassEnergy);
    fListHist->Add(f2DHistLeadingENch0815ZDC);
  }

  //-----------------------------------------------------------------------------
  if (!f2DHistEffEnergyNch0815ZDC)
  {
    f2DHistEffEnergyNch0815ZDC = new TH2D("f2DHistEffEnergyNch0815ZDC", "Effective energy (GeV);Nch 0.8 < |#eta| < 1.5; Leading energy (within ZP and ZN acc)", 200, 0, 200., 1300, 0, fCenterOfMassEnergy);
    fListHist->Add(f2DHistEffEnergyNch0815ZDC);
  }

  //-----------------------------------------------------------------------------
  if (!f2DHistNchNch0815ZDC)
  {
    f2DHistNchNch0815ZDC = new TH2D("f2DHistNchNch0815ZDC", "Nch (|#eta|<0.5);Nch 0.8 < |#eta| < 1.5; Leading energy (within ZP and ZN acc)", 200, 0, 200., 1300, 0, fCenterOfMassEnergy);
    fListHist->Add(f2DHistNchNch0815ZDC);
  }

  //-----------------------------------------------------------------------------
  if (!f2DHistNMPINch0815ZDC && fkDoPythia)
  {
    f2DHistNMPINch0815ZDC = new TH2D("f2DHistNMPINch0815ZDC", "NMPI;Nch 0.8 < |#eta| < 1.5; Leading energy (within ZP and ZN acc)", 200, 0, 200., 1300, 0, fCenterOfMassEnergy);
    fListHist->Add(f2DHistNMPINch0815ZDC);
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!fHistPt[ih]) {
      fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()), Form("Generated %s;p_{T} (GeV/c)",lPartNames[ih].Data()), 250, 0, 25.);
      fListHist->Add(fHistPt[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!fHistDecayVtxPos[ih]) {
      fHistDecayVtxPos[ih] = new TH1D(Form("fHistDecayVtxPos_%s", lPartNames[ih].Data()), Form("Generated %s;Log_{10}(Decay Vtx XY)", lPartNames[ih].Data()), 40, -20, 20.);
      fListHist->Add(fHistDecayVtxPos[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!f2DHistPartSPDV0M[ih]) {
      f2DHistPartSPDV0M[ih] = new TH2D(Form("f2DHistPartSPDV0M_%s",lPartNames[ih].Data()), Form("Generated %s;SPD Clusters; V0M multiplicity",lPartNames[ih].Data()), 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistPartSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!f2DHistPartNch0815V0M[ih]) {
      f2DHistPartNch0815V0M[ih] = new TH2D(Form("f2DHistPartNch0815V0M_%s",lPartNames[ih].Data()), Form("Generated %s;Nch 0.8 < |#eta| < 1.5; V0M multiplicity",lPartNames[ih].Data()), 200, 0, 200., 500, 0, 500.);
      fListHist->Add(f2DHistPartNch0815V0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!f2DHistPartNch0815ZDC[ih]) {
      f2DHistPartNch0815ZDC[ih] = new TH2D(Form("f2DHistPartNch0815ZDC_%s",lPartNames[ih].Data()), Form("Generated %s;Nch 0.8 < |#eta| < 1.5; Leading energy (within ZN and ZP acc)",lPartNames[ih].Data()), 200, 0, 200., 1300, 0., fCenterOfMassEnergy);
      fListHist->Add(f2DHistPartNch0815ZDC[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!f2DHistAvPtSPDV0M[ih]) {
      f2DHistAvPtSPDV0M[ih] = new TH2D(Form("f2DHistAvPtSPDV0M_%s",lPartNames[ih].Data()), "#LT Pt #GT (GeV/c);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistAvPtSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<22; ih++){
    if(!f2DHistPartRecoPercSPDV0M[ih]) {
      f2DHistPartRecoPercSPDV0M[ih] = new TH2D(Form("f2DHistPartRecoPercSPDV0M_%s",lPartNames[ih].Data()), Form("Generated %s;SPD Clusters percentile (reco); V0M percentile (reco)",lPartNames[ih].Data()), 100, 0, 100., 100, 0, 0.);
      fListHist->Add(f2DHistPartRecoPercSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLE) {
    f2dHistZDCVsLE = new TH2D("f2dHistZDCVsLE", ";ZDC Energy Sum (a.u.);Leading energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLE);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsEE) {
    f2dHistZDCVsEE = new TH2D("f2dHistZDCVsEE", ";ZDC Energy Sum (a.u.); Effective energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsEE);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLEA) {
    f2dHistZDCVsLEA = new TH2D("f2dHistZDCVsLEA", ";ZDC-A Energy Sum (a.u.);Leading energy A-side (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLEA);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLEC) {
    f2dHistZDCVsLEC = new TH2D("f2dHistZDCVsLEC", ";ZDC-C Energy Sum (a.u.);Leading energy C-side (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLEC);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZPVsLP) {
    f2dHistZPVsLP = new TH2D("f2dHistZPVsLP", ";ZP Energy Sum (a.u.);Leading proton energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZPVsLP);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZNVsLN) {
    f2dHistZNVsLN = new TH2D("f2dHistZNVsLN", ";ZN Energy Sum (a.u.);Leading neutron (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZNVsLN);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLEnoacc) {
    f2dHistZDCVsLEnoacc = new TH2D("f2dHistZDCVsLEnoacc", ";ZDC Energy Sum (a.u.);Leading energy (|#eta|>8) (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLEnoacc);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZPVsLPnoacc) {
    f2dHistZPVsLPnoacc = new TH2D("f2dHistZPVsLPnoacc", ";ZP Energy Sum (a.u.);Leading proton energy (|#eta|>8) (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZPVsLPnoacc);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZNVsLNnoacc) {
    f2dHistZNVsLNnoacc = new TH2D("f2dHistZNVsLNnoacc", ";ZN Energy Sum (a.u.);Leading neutron (|#eta|>8) (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZNVsLNnoacc);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistSPDClRecoVsTrue) {
    f2dHistSPDClRecoVsTrue = new TH2D("f2dHistSPDClRecoVsTrue", "; SPDClusters centrality (%); SPD Clusters (true)", 100,0.,100., 800, 0, 800.);
    fListHist->Add(f2dHistSPDClRecoVsTrue);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistV0MRecoVsTrue) {
    f2dHistV0MRecoVsTrue = new TH2D("f2dHistV0MRecoVsTrue", ";V0M centrality (%); V0M Multiplicity (true)",100,0.,100., 500, 0, 500.);
    fListHist->Add(f2dHistV0MRecoVsTrue);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistTrueVsRecoSPDCl) {
    f2dHistTrueVsRecoSPDCl = new TH2D("f2dHistTrueVsRecoSPDCl", ";SPD Clusters (reco); SPD Clusters (true)",800,0.,800., 800, 0, 800.);
    fListHist->Add(f2dHistTrueVsRecoSPDCl);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistTrueVsRecoSPDCl0) {
    f2dHistTrueVsRecoSPDCl0 = new TH2D("f2dHistTrueVsRecoSPDCl0", ";SPD Cluster 0 (reco); SPD Cluster 0 (true)",800,0.,800., 800, 0, 800.);
    fListHist->Add(f2dHistTrueVsRecoSPDCl0);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistTrueVsRecoSPDCl1) {
    f2dHistTrueVsRecoSPDCl1 = new TH2D("f2dHistTrueVsRecoSPDCl1", ";SPD Cluster 1 (reco); SPD Cluster 1 (true)",800,0.,800., 800, 0, 800.);
    fListHist->Add(f2dHistTrueVsRecoSPDCl1);
  }

  //-----------------------------------------------------------------------------
  if(!f3dHistPi0SPDMultSPDCl) {
    f3dHistPi0SPDMultSPDCl = new TH3D("f3dHistPi0SPDMultSPDCl", ";SPD Multiplicity |#eta|<1 (true); SPD Clusters (reco);Pi0 Counts",500,0.,500., 500, 0, 500.,200,0,200.);
    fListHist->Add(f3dHistPi0SPDMultSPDCl);
  }

  //-----------------------------------------------------------------------------
  if(!p2dHistPi0SPDMultSPDCl) {
    p2dHistPi0SPDMultSPDCl = new TProfile2D("p2dHistPi0SPDMultSPDCl", "SPD Clusters (reco);SPD Multiplicity |#eta|<1 (true); N Pi0",200,0.,200.,500,0.,500.,0.,500.);
    fListHist->Add(p2dHistPi0SPDMultSPDCl);
  }

  //-----------------------------------------------------------------------------
  if (!f2dHistTrueVsRecoNch0815) {
    f2dHistTrueVsRecoNch0815 = new TH2D("f2dHistTrueVsRecoNch0815", ";Nch 0.8 < |#eta| < 1.5 (reco); Nch 0.8 < |#eta| < 1.5 (true)", 200, 0., 200., 200, 0., 200.);
    fListHist->Add(f2dHistTrueVsRecoNch0815);
  }

  //List of Histograms: Normal
  PostData(1, fListHist);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::UserExec(Option_t *)
{
  // Main loop --> called for each event
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  AliESDEvent *lESDevent = 0x0;
  //
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  //
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  //
  AliGenEventHeader* mcGenH = lMCevent->GenEventHeader();
  //
  if (lMCstack->GetNprimary() != lMCstack->GetNtrack()) lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );

  //Get primary vertex position
  Double_t mcVx = lMCevent->GetPrimaryVertex()->GetX();
  Double_t mcVy = lMCevent->GetPrimaryVertex()->GetY();
  Double_t mcVz = lMCevent->GetPrimaryVertex()->GetZ();

  //Events Processed
  fHistEventCounter->Fill(0.5);

  //Multiplicity
  Long_t lNchEta05         = 0;
  Long_t lNchEta10         = 0;
  Long_t lNchEta08         = 0;
  Long_t lNchEta08to15     = 0;
  Long_t lSPDCl0           = 0;
  Long_t lSPDCl1           = 0;
  Long_t lNchVZEROA        = 0;
  Long_t lNchVZEROC        = 0;
  Bool_t lEvSel_INELgtZERO = kFALSE;

  //Eff energy
  Float_t fLeadingE = 0.;
  Float_t fLeadingE_Aside = 0.;
  Float_t fLeadingE_Cside = 0.;
  Float_t fLeadingP = 0.;
  Float_t fLeadingN = 0.;
  Float_t fEffEnergy = fCenterOfMassEnergy;
  Float_t fLeadingEnoacc = 0.;
  Float_t fLeadingPnoacc = 0.;
  Float_t fLeadingNnoacc = 0.;

  //Utility
  const Float_t c = 299792458; //m/s
  Bool_t lIsPhysicalPrimary = kFALSE;

  //Particle info
  Int_t lPartCounter[22];
  for(Int_t i = 0; i<22; i++){
    lPartCounter[i] = 0;
  }

  TString lPartNames[22] = {
    "PiPlus", "PiMinus",
    "KaPlus", "KaMinus",
    "Proton", "AntiProton",
    "K0Short",
    "Lambda", "AntiLambda",
    "XiMinus", "XiPlus",
    "OmegaMinus", "OmegaPlus",
    "Phi",
    "D0", "AntiD0",
    "DPlus", "DMinus",
    "Lambdac", "AntiLambdac",
    "JPsi",
    "Pi0"
  };
  Int_t lPDGCodes[22] = {
    211, -211,
    321, -321,
    2212, -2212,
    310,
    3122, -3122,
    3312, -3312,
    3334, -3334,
    333,
    421, -421,
    411, -411,
    4122, -4122,
    443,
    111
  };

  Bool_t lCheckIsPhysicalPrimary[22] = {
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kTRUE, kTRUE,
    kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE, kFALSE,
    kFALSE,
    kFALSE
  };

  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  { // This is the begining of the loop on tracks

    TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
    if(!particleOne) continue;
    if(!particleOne->GetPDG()) continue;

    Double_t geta = particleOne -> Eta();
    Int_t charge = particleOne->GetPDG()->Charge();
    Double_t partcharge = particleOne->GetPDG()->Charge()/3.;

    //Vtx position of the particle
    Double_t vtxpart = particleOne->R();

    //Vtx position of daughters
    Int_t idau = particleOne->GetFirstDaughter();
    TParticle* firstdau;
    Double_t vtxdau = 9999;
    if(idau>=0 && idau<lMCstack->GetNtrack()) {
      firstdau = lMCstack->Particle(idau);
      vtxdau = firstdau->R();
    }
    Float_t AbsVz = TMath::Abs(particleOne->Vz() - mcVz);

    if( (TMath::Abs(partcharge)>0.001) ) {
      if(vtxpart <= 7. && vtxdau > 7.) {
        if( (TMath::Abs(geta) < 1.5928278) && AbsVz < 16.5 ) lSPDCl1++;
      }
      if(vtxpart <= 4. && vtxdau > 4.) {
        if( (TMath::Abs(geta) < 2.1245920) && AbsVz < 16.5 ) lSPDCl0++;
      }
    }

    //if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) )continue;
    //if(particleOne->GetFirstDaughter()>0) continue;
    if(CheckIsNotPrimary(!lESDevent,lMCstack, iCurrentLabelStack, particleOne)) continue;

    if( TMath::Abs(geta)>8. ) {
      fLeadingEnoacc += particleOne -> Energy();
      if (charge>0) fLeadingPnoacc += particleOne -> Energy();
      if (charge==0) fLeadingNnoacc += particleOne -> Energy();
    }

    if( (charge>0 && TMath::Abs(geta)>7. && TMath::Abs(geta)<8.7) ||
        (charge==0 && TMath::Abs(geta)>8.8)
      ) {
        fEffEnergy -= particleOne -> Energy();
        fLeadingE += particleOne -> Energy();

        if( (charge>0 && TMath::Abs(geta)>7. && TMath::Abs(geta)<8.7) ) fLeadingP += particleOne -> Energy();
        if (charge==0 && TMath::Abs(geta)>8.8) fLeadingN += particleOne -> Energy();

        if((charge>0 && geta>7. && geta<8.7) || (charge==0 && geta>8.8)) {
          fLeadingE_Aside += particleOne -> Energy();
        } else if((charge>0 && geta>-8.7 && geta<-7) || (charge==0 && geta<-8.8)) {
          fLeadingE_Cside += particleOne -> Energy();
        }
      }

    if(TMath::Abs(partcharge)<0.001) continue; //now only charged primaries

    //
    if( TMath::Abs(geta) < 0.5 ) lNchEta05++;
    if( TMath::Abs(geta) < 0.8 ) lNchEta08++;
    if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
    if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta08to15++;
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZERO = kTRUE;
    if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
    if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;

  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------

  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZERO && fkSelectINELgtZERO ) return;

  //
  Long_t lSPDClusters = lSPDCl0 + lSPDCl1;

  //Events INEL>0
  fHistEventCounter->Fill(1.5);

  //
  Int_t fMC_NMPI = -1;
  Float_t fMC_Q2 = -1;
  Float_t fMC_b = -1;
  //Int_t fMC_22 = -1;
  //Int_t fMC_NColl = -1;

  //PYTHIA
  if ( fkDoPythia ){
    if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
      AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
      if(fMcPythiaHeader){
        fMC_NMPI = fMcPythiaHeader->GetNMPI();
        fMC_Q2 = fMcPythiaHeader->GetPtHard();
      }
    }
  }

  //EPOS
  if ( fkDoEPOS ){
    if (mcGenH->InheritsFrom(AliGenHepMCEventHeader::Class())) {
      AliGenHepMCEventHeader * lHepMCHeader = dynamic_cast <AliGenHepMCEventHeader*> (mcGenH);
      if (lHepMCHeader ){
        //fMC_22 = lHepMCHeader->22_proj()+lHepMCHeader->22_targ();
        //fMC_NColl = lHepMCHeader->N_Nwounded_collisions() + lHepMCHeader->Nwounded_N_collisions() + lHepMCHeader->Nwounded_Nwounded_collisions();
        fMC_b = lHepMCHeader->impact_parameter();
      }
    }
  }

  //1D histograms
  if(fHistV0MMult)      fHistV0MMult        -> Fill ( lNchVZEROA+lNchVZEROC );
  if(fHistV0AMult)      fHistV0AMult        -> Fill ( lNchVZEROA );
  if(fHistV0CMult)      fHistV0CMult        -> Fill ( lNchVZEROC );
  if(fHistMult05)       fHistMult05         -> Fill ( lNchEta05 );
  if(fHistMult08)       fHistMult08         -> Fill ( lNchEta08 );
  if(fHistMult10)       fHistMult10         -> Fill ( lNchEta10 );
  if(fHistMult08to15)   fHistMult08to15     -> Fill ( lNchEta08to15 );
  if(fHistSPDClusters)  fHistSPDClusters    -> Fill ( lSPDClusters );
  if(fHistSPDCl0)       fHistSPDCl0         -> Fill ( lSPDCl0 );
  if(fHistSPDCl1)       fHistSPDCl1         -> Fill ( lSPDCl1 );
  if(fHistNMPI)         fHistNMPI           -> Fill ( fMC_NMPI );
  if(fHistQ2)           fHistQ2             -> Fill ( fMC_Q2 );
  if(fHistb)            fHistb              -> Fill ( fMC_b );
  if(fHistLeadingE)     fHistLeadingE       -> Fill ( fLeadingE );
  if(fHistEffEnergy)    fHistEffEnergy      -> Fill ( fEffEnergy );

  //----- Loop on Stack ----------------------------------------------------------------
  Int_t nPrimaries = lMCstack->GetNprimary();
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  { // This is the begining of the loop on tracks

    TParticle* part = lMCstack->Particle(iCurrentLabelStack);
    if(!part) continue;
    if(!part->GetPDG()) continue;
    lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(iCurrentLabelStack);
    Bool_t IsPrimary = part->IsPrimary();

    Double_t charge = part -> GetPDG()->Charge();
    Float_t px      = part -> Px();
    Float_t py      = part -> Py();
    Float_t pz      = part -> Pz();
    Float_t pt      = part -> Pt();
    Float_t energy  = part -> Energy();
    Float_t eta     = part -> Eta();
    Int_t pdg       = (Int_t) part -> GetPdgCode();
    Float_t y       = Rapidity(energy, pz);
    Float_t Vx      = part->Vx();
    Float_t Vy      = part->Vy();
    Double_t DeltaR = TMath::Sqrt( (Vx - mcVx)*(Vx - mcVx) + (Vy - mcVy)*(Vy - mcVy) );

    // Vtx position of daughters
    Int_t idau = part->GetFirstDaughter();
    TParticle *firstdau;
    Double_t vtxdau = 9999;
    if (idau >= 0 && idau < lMCstack->GetNtrack())
    {
      firstdau = lMCstack->Particle(idau);
      vtxdau = firstdau->R();
    }


    if(TMath::Abs(y)<0.5){
      for(Int_t ih=0; ih<22; ih++){ //loop over pdg codes
        if( pdg == lPDGCodes[ih] ) {

          Bool_t IsPrimary = kFALSE;
          if (!lESDevent){
            // if Fast generator
            IsPrimary = lIsPhysicalPrimary;
          } else {
            // if Full MC
            if (part->GetFirstDaughter() >= nPrimaries) IsPrimary = (iCurrentLabelStack < nPrimaries); // drop if the particle has a daughter among the primaries
            if (part->GetStatusCode() != 1) IsPrimary = kFALSE; // drop non final state particles
            if (DeltaR > 1E-6) IsPrimary = kFALSE; // drop if the particle is not produced in the primary vertex
          }

          //Check if Phyisical Primary if needed
          if( lCheckIsPhysicalPrimary[ih] == kTRUE && IsPrimary == kFALSE) continue;
          if( lCheckIsPhysicalPrimary[ih] == kFALSE && DeltaR>1) continue;

          //Fill histos
          if(fHistPt[ih]) fHistPt[ih] -> Fill( pt );
          if(fHistDecayVtxPos[ih]) fHistDecayVtxPos[ih] -> Fill( TMath::Log10(vtxdau+1E-20) );
          if(f2DHistAvPtSPDV0M[ih]) f2DHistAvPtSPDV0M[ih] -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, pt );
          if(f2DHistPartSPDV0M[ih]) f2DHistPartSPDV0M[ih] -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC );
          if(f2DHistPartNch0815V0M[ih]) f2DHistPartNch0815V0M[ih] -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC );
          if(f2DHistPartNch0815ZDC[ih]) f2DHistPartNch0815ZDC[ih] -> Fill( lNchEta08to15, fLeadingE );

          //Counter for specific particles in this event
          lPartCounter[ih]++;
        }
      }//end loop over pdg codes
    }
  } //end loop on tracks

  if(f2DHistINELgt0SPDV0M)   f2DHistINELgt0SPDV0M   -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC );
  if(f2DHistLeadingESPDV0M)  f2DHistLeadingESPDV0M  -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fLeadingE );
  if(f2DHistEffEnergySPDV0M) f2DHistEffEnergySPDV0M -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fEffEnergy );
  if(f2DHistNchSPDV0M)       f2DHistNchSPDV0M       -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, lNchEta05 );
  if(f2DHistNMPISPDV0M)      f2DHistNMPISPDV0M      -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_NMPI );
  if(f2DHistQ2SPDV0M)        f2DHistQ2SPDV0M        -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_Q2 );
  if(f2DHistbSPDV0M)         f2DHistbSPDV0M         -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_b );
  if(f2DHistINELgt0Nch0815V0M)   f2DHistINELgt0Nch0815V0M   -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC );
  if(f2DHistLeadingENch0815V0M)  f2DHistLeadingENch0815V0M  -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC, fLeadingE );
  if(f2DHistEffEnergyNch0815V0M) f2DHistEffEnergyNch0815V0M -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC, fEffEnergy );
  if(f2DHistNchNch0815V0M)       f2DHistNchNch0815V0M       -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC, lNchEta05 );
  if(f2DHistNMPINch0815V0M)      f2DHistNMPINch0815V0M      -> Fill( lNchEta08to15, lNchVZEROA+lNchVZEROC, fMC_NMPI );
  if(f2DHistINELgt0Nch0815ZDC)   f2DHistINELgt0Nch0815ZDC   -> Fill( lNchEta08to15, fLeadingE );
  if(f2DHistLeadingENch0815ZDC)  f2DHistLeadingENch0815ZDC  -> Fill( lNchEta08to15, fLeadingE, fLeadingE );
  if(f2DHistEffEnergyNch0815ZDC) f2DHistEffEnergyNch0815ZDC -> Fill( lNchEta08to15, fLeadingE, fEffEnergy );
  if(f2DHistNchNch0815ZDC)       f2DHistNchNch0815ZDC       -> Fill( lNchEta08to15, fLeadingE, lNchEta05 );
  if(f2DHistNMPINch0815ZDC)      f2DHistNMPINch0815ZDC      -> Fill( lNchEta08to15, fLeadingE, fMC_NMPI );


  //Reco information
  if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
  }

  Float_t lV0MPercentile = 300;
  Float_t lSPDClusterspercentile = 300;

  AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
  if( !MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    return;
  } else {
      lV0MPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      lSPDClusterspercentile = MultSelection->GetMultiplicityPercentile("SPDClusters");
  }

  fCentrality_V0M = lV0MPercentile;
  fCentrality_SPDClusters = lSPDClusterspercentile;

  AliMultiplicity* multiplicity =  lESDevent->GetMultiplicity();
  if( !multiplicity) {
    AliWarning("AliMultiplicity object not found!");
    return;
  }
  Double_t lSPDCl = multiplicity->GetNumberOfITSClusters(0) + multiplicity->GetNumberOfITSClusters(1);
  Int_t lSPDtracklets0815 = 0;
	for (auto it = 0; it<multiplicity->GetNumberOfTracklets(); it++) {
        Double_t eta = multiplicity->GetEta(it);
        if ( TMath::Abs(eta)<1.5 && TMath::Abs(eta)>0.8 ){
            lSPDtracklets0815++;
        }
	}

  // ZDC info ==========================================================
  const Double_t *aZDCN1 = lESDevent->GetESDZDC()->GetZNCTowerEnergy();
  fZNCpp = aZDCN1[0];
  const Double_t *aZDCN2 = lESDevent->GetESDZDC()->GetZNATowerEnergy();
  fZNApp = aZDCN2[0];
  const Double_t *aZDCP1 = lESDevent->GetESDZDC()->GetZPCTowerEnergy();
  fZPCpp = aZDCP1[0];
  const Double_t *aZDCP2 = lESDevent->GetESDZDC()->GetZPATowerEnergy();
  fZPApp = aZDCP2[0];

  if(f2dHistZDCVsLE)    f2dHistZDCVsLE         -> Fill ( fZNCpp+fZNApp+fZPCpp+fZPApp , fLeadingE );
  if(f2dHistZDCVsEE)    f2dHistZDCVsEE         -> Fill ( fZNCpp+fZNApp+fZPCpp+fZPApp , fEffEnergy );
  if(f2dHistZDCVsLEA)   f2dHistZDCVsLEA        -> Fill ( fZNApp+fZPApp , fLeadingE_Aside );
  if(f2dHistZDCVsLEC)   f2dHistZDCVsLEC        -> Fill ( fZNCpp+fZPCpp , fLeadingE_Cside );
  if(f2dHistZPVsLP)     f2dHistZPVsLP          -> Fill ( fZPCpp+fZPApp , fLeadingP );
  if(f2dHistZNVsLN)     f2dHistZNVsLN          -> Fill ( fZNCpp+fZNApp , fLeadingN );
  if(f2dHistZPVsLPnoacc)     f2dHistZPVsLPnoacc          -> Fill ( fZPCpp+fZPApp , fLeadingPnoacc );
  if(f2dHistZNVsLNnoacc)     f2dHistZNVsLNnoacc          -> Fill ( fZNCpp+fZNApp , fLeadingNnoacc );
  if(f2dHistZDCVsLEnoacc)    f2dHistZDCVsLEnoacc         -> Fill ( fZNCpp+fZNApp+fZPCpp+fZPApp , fLeadingEnoacc );
  if(f2dHistSPDClRecoVsTrue)    f2dHistSPDClRecoVsTrue    -> Fill ( fCentrality_SPDClusters , lSPDClusters );
  if(f2dHistV0MRecoVsTrue)      f2dHistV0MRecoVsTrue      -> Fill ( fCentrality_V0M , lNchVZEROA+lNchVZEROC );
  if(f2dHistTrueVsRecoSPDCl)    f2dHistTrueVsRecoSPDCl    -> Fill ( lSPDCl , lSPDClusters );
  if(f2dHistTrueVsRecoSPDCl0)    f2dHistTrueVsRecoSPDCl0  -> Fill ( multiplicity->GetNumberOfITSClusters(0), lSPDCl0 );
  if(f2dHistTrueVsRecoSPDCl1)    f2dHistTrueVsRecoSPDCl1  -> Fill ( multiplicity->GetNumberOfITSClusters(1), lSPDCl1 );
  if(f3dHistPi0SPDMultSPDCl)    f3dHistPi0SPDMultSPDCl    -> Fill ( lNchEta10 , lPartCounter[21], lSPDCl );
  if(p2dHistPi0SPDMultSPDCl)    p2dHistPi0SPDMultSPDCl    -> Fill ( lNchEta10 , lPartCounter[21], lSPDCl );
  if(f2DHistNchRecoPercSPDV0M)   f2DHistNchRecoPercSPDV0M  -> Fill( fCentrality_SPDClusters, fCentrality_V0M, lNchEta05 );
  if(f2DHistINELgt0RecoPercSPDV0M)   f2DHistINELgt0RecoPercSPDV0M  -> Fill( fCentrality_SPDClusters, fCentrality_V0M );
  if(f2dHistTrueVsRecoNch0815)   f2dHistTrueVsRecoNch0815  -> Fill( lSPDtracklets0815, lNchEta08to15 );
  for(Int_t ih=0; ih<22; ih++){ //loop over pdg codes
    if(f2DHistPartRecoPercSPDV0M[ih])   f2DHistPartRecoPercSPDV0M[ih]  -> Fill( fCentrality_SPDClusters, fCentrality_V0M, lPartCounter[ih]);
  }

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
    Printf("ERROR - AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : ouput data container list not available\n");
    return;
  }

  fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
  if (!fHistEventCounter) {
    Printf("ERROR - AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : fHistEventCounter not available");
    return;
  }

  TCanvas *canCheck = new TCanvas("AliAnalysisTaskMCPredictionsStrgVsMultVsZDC","Event Multiplicity",10,10,510,510);
  canCheck->cd(1)->SetLogy();

  fHistEventCounter->SetMarkerStyle(22);
  fHistEventCounter->DrawCopy("E");
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::IsEPOSLHC() const {
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
Double_t AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::Rapidity(Double_t E, Double_t Pz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (E - Pz + 1.e-13) != 0 && (E + Pz) != 0 ) {
    ReturnValue =  0.5*TMath::Log((E + Pz)/( E - Pz + 1.e-13));
  }
  return ReturnValue;
}

//_______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::CheckIsNotPrimary(Bool_t IsFastGenerator, AliStack* MCstack, Int_t iLabelStack, TParticle* part) const
{
  //Returns true is the particle is NOT primary, returns false otherwise

  Bool_t isnotprimary = kFALSE;
  Int_t nPrimaries = MCstack->GetNprimary();

  if (IsFastGenerator){
    if(! (MCstack->IsPhysicalPrimary(iLabelStack)) ) isnotprimary = kTRUE;
    if(part->GetFirstDaughter()>0) isnotprimary = kTRUE; // needed to exclude long-lived particles which are enabled using https://github.com/alisw/AliPhysics/blob/master/PWGLF/STRANGENESS/Cascades/Run2/macros/AddCustomMCGenPythia8.C#L253
    // See https://github.com/alisw/AliRoot/blob/48577f11e3d3c23f3ec117b11c0ff5fbd2ea9227/PYTHIA8/AliPythia8/AliPythia8.h#L69
  } else {
    if(part->GetPDG()->Charge() == 0) {
      if (part->GetFirstDaughter() == -1 || (part->GetFirstDaughter() >= nPrimaries)) {
        isnotprimary = (iLabelStack >= nPrimaries);
      } else {
        isnotprimary = kTRUE;
      }
      if (part->GetPdgCode() == 21) isnotprimary = kTRUE;
      if (part->GetStatusCode() != 1) isnotprimary = kTRUE;
    } else{
         isnotprimary = (!AliPWG0Helper::IsPrimaryCharged(part, nPrimaries)); // official definition of charged primary
    }
  }

  return isnotprimary;
}
