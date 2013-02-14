/************************************************************************* 
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
 *                                                                        * 
 * Author: X. Sanchez Castro                                              * 
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

#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TROOT.h>

#include "AliOADBContainer.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliCFContainer.h"
#include "AliCentrality.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskLambdaOverK0sJets.h"

extern TROOT *gROOT;


ClassImp(AliAnalysisTaskLambdaOverK0sJets)
ClassImp(AliMiniParticle)

// Global variables:
static Int_t    nbins = 100;                 // Number of bins for l, pt, mass for V0
static Int_t    nbinsPhi = 120;              // Number of bins for Phi
static Int_t    nbinsdPhi = 42;              // Number of bins for dPhi
static Int_t    nbinsdEta = 40;              // Number of bins for dEta
static Int_t    nbinPtLP = 200;
static Int_t    nbinsVtx = 20;

static Float_t pMin = 0.0;                  // Lower cut for transverse momentum
static Float_t pMax = 10.;                  // Max cut for transverse momentum for V0
static Float_t ptMaxLP = 50.;               // Max cut for transverse momentum LP

static Float_t lMin = 0.0;                  // Limits in the histo for fidutial volume
static Float_t lMax = 100.;                 // Limits in the fidutial volume

static Int_t   nMaxEvMix = 250;

//
//  
//

AliAnalysisTaskLambdaOverK0sJets::AliAnalysisTaskLambdaOverK0sJets(const char *name) :
  AliAnalysisTaskSE(name),

  fAOD(0), fIsMC(kFALSE), fData("PbPb2010"), fUsePID(kFALSE), fCentMin(0.), fCentMax(90.), fDoQA(kFALSE), fDoMixEvt(kFALSE), fTrigPtMin(8.), fTrigPtMax(20.), fTrigEtaMax(0.8), fCheckIDTrig(kFALSE), fSeparateInjPart(kTRUE), fEndOfHijingEvent(-1),  fPIDResponse(0),

  fMinPtDaughter(0.160), fMaxEtaDaughter(0.8), fMaxDCADaughter(1.0), fYMax(0.5), fDCAToPrimVtx(0.1), fMinCPA(0.998), fNSigma(3.0), fMinCtau(0.), fMaxCtau(3.), fIdTrigger(-1), fIsV0LP(0), fPtV0LP(0.), fIsSndCheck(0),

  fOutput(0), fOutputQA(0), fOutputME(0), fMEList(0x0), fTriggerParticles(0x0), fAssocParticles(0x0), fEvents(0), fCentrality(0), fPrimaryVertexX(0), fPrimaryVertexY(0), fPrimaryVertexZ(0), fCentMult(0), fdEdx(0), fdEdxPid(0), 

  fTriggerMCPtCent(0), fTriggerPtCent(0), fNTrigPerEvt(0), fTriggerWiSPDHit(0), fTriggerEtaPhi(0), fCheckTriggerFromV0Daug(0), fTriggerComingFromDaug(0), fTriggerIsV0(0), fCheckIDTrigPtK0s(0), fCheckIDTrigPhiK0s(0),  fCheckIDTrigEtaK0s(0), fCheckIDTrigPtLambda(0), fCheckIDTrigPhiLambda(0), fCheckIDTrigEtaLambda(0), fCheckIDTrigPtAntiLambda(0), fCheckIDTrigPhiAntiLambda(0),fCheckIDTrigEtaAntiLambda(0), fInjectedParticles(0),

  fK0sMCPt(0), fK0sMCPtRap(0), fK0sMCPtPhiEta(0), fK0sAssocPt(0), fK0sAssocPtArm(0), fK0sAssocPtRap(0), fK0sAssocPtPhiEta(0), fK0sMCResPhi(0), fLambdaMCPt(0), fLambdaMCPtRap(0), fLambdaMCPtPhiEta(0), fLambdaAssocPt(0), fLambdaAssocPtArm(0), fLambdaAssocPtRap(0), fLambdaAssocPtPhiEta(0), fLambdaMCResPhi(0), fAntiLambdaMCPt(0), fAntiLambdaMCPtRap(0), fAntiLambdaMCPtPhiEta(0), fAntiLambdaAssocPt(0), fAntiLambdaAssocPtArm(0), fAntiLambdaAssocPtRap(0), fAntiLambdaAssocPtPhiEta(0), fAntiLambdaMCResPhi(0),

  fHistArmenterosPodolanski(0), fHistArmPodBckg(0),
  
  fK0sMass(0), fK0sPtvsEta(0), fK0sPtvsRap(0), fK0sEtaPhi(0), fK0sMassPtPhi(0), fK0sSiPtL(0), fK0sDaughtersPt(0), fK0sDCADaugToPrimVtx(0), fK0sSpatialRes(0), fK0sBckgDecLength(0), fK0sBckgDCADaugToPrimVtx(0), fK0sdEdxPosDaug(0), fK0sdEdxNegDaug(0), fK0sBckgEtaPhi(0), fK0sBckgPhiRadio(0), fK0sBckgDCANegDaugToPrimVtx(0), fK0sBckgDCAPosDaugToPrimVtx(0), fV0MassCascade(0),
  
  fLambdaMass(0), fLambdaPtvsEta(0), fLambdaPtvsRap(0), fLambdaEtaPhi(0), fLambdaMassPtPhi(0), fLambdaSiPtL(0), fLambdaDaughtersPt(0), fLambdaDCADaugToPrimVtx(0), fLambdaSpatialRes(0), fLambdaBckgDecLength(0), fLambdaBckgDCADaugToPrimVtx(0), fLambdadEdxPosDaug(0), fLambdadEdxNegDaug(0), fLambdaBckgEtaPhi(0), fLambdaBckgPhiRadio(0), fLambdaBckgDCANegDaugToPrimVtx(0), fLambdaBckgDCAPosDaugToPrimVtx(0), 

  fAntiLambdaMass(0), fAntiLambdaPtvsEta(0), fAntiLambdaPtvsRap(0), fAntiLambdaEtaPhi(0), fAntiLambdaMassPtPhi(0), fAntiLambdaSiPtL(0), fAntiLambdaDaughtersPt(0), fAntiLambdaDCADaugToPrimVtx(0), fAntiLambdaSpatialRes(0), fAntiLambdaBckgDecLength(0), fAntiLambdaBckgDCADaugToPrimVtx(0), fAntiLambdadEdxPosDaug(0), fAntiLambdadEdxNegDaug(0), fAntiLambdaBckgEtaPhi(0), fAntiLambdaBckgPhiRadio(0), fAntiLambdaBckgDCANegDaugToPrimVtx(0), fAntiLambdaBckgDCAPosDaugToPrimVtx(0), 

  fK0sPIDPosDaug(0), fK0sPIDNegDaug(0), fK0sBckgPIDPosDaug(0), fK0sBckgPIDNegDaug(0), fK0sPhiEtaPosDaug(0), fK0sPhiEtaNegDaug(0), fK0sBckgPhiEtaPosDaug(0), fK0sBckgPhiEtaNegDaug(0), fK0sDCAPosDaug(0), fK0sDCANegDaug(0), fK0sBckgDCAPosDaug(0), fK0sBckgDCANegDaug(0), fK0sDifPtPosDaug(0), fK0sDifPtNegDaug(0), fK0sBckgDifPtPosDaug(0), fK0sBckgDifPtNegDaug(0), fK0sDecayPos(0), fK0sBckgDecayPos(0), fK0sDecayVertex(0), fK0sBckgDecayVertex(0), fK0sDecayVertexZoom(0), fK0sBckgDecayVertexZoom(0), fK0sCPA(0), fK0sBckgCPA(0), fK0sDCAV0Daug(0), fK0sBckgDCAV0Daug(0), fK0sNClustersTPC(0), fK0sBckgNClustersTPC(0), fK0sNClustersITSPos(0), fK0sNClustersITSNeg(0), fK0sBckgNClustersITSPos(0),  fK0sBckgNClustersITSNeg(0),   

  fLambdaPIDPosDaug(0), fLambdaPIDNegDaug(0), fLambdaBckgPIDPosDaug(0), fLambdaBckgPIDNegDaug(0), fLambdaPhiEtaPosDaug(0),fLambdaPhiEtaNegDaug(0), fLambdaBckgPhiEtaPosDaug(0),fLambdaBckgPhiEtaNegDaug(0), fLambdaDCAPosDaug(0),fLambdaDCANegDaug(0), fLambdaBckgDCAPosDaug(0), fLambdaBckgDCANegDaug(0), fLambdaDifPtPosDaug(0), fLambdaDifPtNegDaug(0), fLambdaBckgDifPtPosDaug(0), fLambdaBckgDifPtNegDaug(0), fLambdaDecayPos(0), fLambdaBckgDecayPos(0), fLambdaDecayVertex(0), fLambdaBckgDecayVertex(0), fLambdaDecayVertexZoom(0), fLambdaBckgDecayVertexZoom(0), fLambdaCPA(0), fLambdaBckgCPA(0), fLambdaDCAV0Daug(0), fLambdaBckgDCAV0Daug(0), fLambdaNClustersTPC(0), fLambdaBckgNClustersTPC(0), fLambdaNClustersITSPos(0), fLambdaNClustersITSNeg(0), fLambdaBckgNClustersITSPos(0),  fLambdaBckgNClustersITSNeg(0),

  fAntiLambdaPIDPosDaug(0), fAntiLambdaPIDNegDaug(0), fAntiLambdaBckgPIDPosDaug(0), fAntiLambdaBckgPIDNegDaug(0), fAntiLambdaPhiEtaPosDaug(0),fAntiLambdaPhiEtaNegDaug(0), fAntiLambdaBckgPhiEtaPosDaug(0),fAntiLambdaBckgPhiEtaNegDaug(0), fAntiLambdaDCAPosDaug(0),fAntiLambdaDCANegDaug(0), fAntiLambdaBckgDCAPosDaug(0), fAntiLambdaBckgDCANegDaug(0), fAntiLambdaDifPtPosDaug(0), fAntiLambdaDifPtNegDaug(0), fAntiLambdaBckgDifPtPosDaug(0), fAntiLambdaBckgDifPtNegDaug(0), fAntiLambdaDecayPos(0), fAntiLambdaBckgDecayPos(0), fAntiLambdaDecayVertex(0), fAntiLambdaBckgDecayVertex(0), fAntiLambdaDecayVertexZoom(0), fAntiLambdaBckgDecayVertexZoom(0), fAntiLambdaCPA(0), fAntiLambdaBckgCPA(0), fAntiLambdaDCAV0Daug(0), fAntiLambdaBckgDCAV0Daug(0), fAntiLambdaNClustersTPC(0), fAntiLambdaBckgNClustersTPC(0), fAntiLambdaNClustersITSPos(0), fAntiLambdaNClustersITSNeg(0), fAntiLambdaBckgNClustersITSPos(0),  fAntiLambdaBckgNClustersITSNeg(0)
  
{
  // Dummy Constructor

  for (Int_t i=0; i<kNCent*kN1; i++){ 
    
    // K0s
    fK0sdPhidEtaMC[i] = 0;
    fK0sdPhidEtaPtL[i] = 0;
    fK0sdPhidEtaPtLBckg[i] = 0;
    
    // Lambda
    fLambdadPhidEtaMC[i] = 0;
    fLambdadPhidEtaPtL[i] = 0;
    fLambdadPhidEtaPtLBckg[i] = 0;
    
    // AntiLambda
    fAntiLambdadPhidEtaMC[i] = 0;
    fAntiLambdadPhidEtaPtL[i] = 0;
    fAntiLambdadPhidEtaPtLBckg[i] = 0;  
  }

  // Gamma Conversion
  for (Int_t i=0; i<kNCent; i++)
    fGammaConversiondPhidEta[i] = 0;

  // Mixed events distributions
  for (Int_t i=0; i<(kN1*kNVtxZ*kNCent); i++){ 
    fK0sdPhidEtaME[i] = 0;
    fLambdadPhidEtaME[i] = 0;
    fAntiLambdadPhidEtaME[i] = 0;
  }
 
  // Constructor. Initialization of pointers
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
 

}

//___________________________________________________________________________________________

AliAnalysisTaskLambdaOverK0sJets::~AliAnalysisTaskLambdaOverK0sJets() 
{

  // Destructor
  if(fMEList){
    
    for(Int_t icent=0; icent<kNCent; icent++){
      for(Int_t iz=0; iz<kNVtxZ; iz++){
	fMEList[icent*kNVtxZ+iz]->Delete();  delete fMEList[icent*kNVtxZ+iz];
      }
    }
    delete[] fMEList; fMEList=0x0;
  }
  
  if(fTriggerParticles) {
    delete fTriggerParticles;
    fTriggerParticles=0x0;
  }

  if(fAssocParticles) {
    delete fAssocParticles;
    fAssocParticles=0x0;
  }
    
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::UserCreateOutputObjects()
{ 
  // Creating the histograms that are needed for the output 
  
  fOutput = new TList(); 
  fOutput->SetOwner();

  fOutputQA = new TList(); 
  fOutputQA->SetOwner();

  fOutputME = new TList(); 
  fOutputME->SetOwner();


  fMEList = new TList*[kNCent*kNVtxZ];
  for(Int_t icent=0; icent<kNCent; icent++){
    for(Int_t iz=0; iz<kNVtxZ; iz++){
      fMEList[icent*kNVtxZ+iz] = new TList();
      fMEList[icent*kNVtxZ+iz]->SetOwner(kFALSE);
    }
  }


  // ====== General characteristics of the event and tracks ====== //

  // Counter for the number of events in each step:
  fEvents=new TH1F("fEvents","Number of events",14,-0.5,13.5);
  fEvents->GetXaxis()->SetBinLabel(1,"calls to UserExec()");
  fEvents->GetXaxis()->SetBinLabel(2,"AOD available");
  fEvents->GetXaxis()->SetBinLabel(3,"CINT1B");
  fEvents->GetXaxis()->SetBinLabel(4,"V0M Cent");
  fEvents->GetXaxis()->SetBinLabel(5,"Vtx > 3 part");
  fEvents->GetXaxis()->SetBinLabel(6,"|VtxZ| < 10 cm");
  fEvents->GetXaxis()->SetBinLabel(7,"Mult && Cent");
  fEvents->GetXaxis()->SetBinLabel(8,"Bad ID Trigger");
  fEvents->GetXaxis()->SetBinLabel(9,"V0 is LP");
  fEvents->GetXaxis()->SetBinLabel(10,"Trigger is V0 daug");
  fEvents->GetXaxis()->SetBinLabel(11,"Trigger is V0 daug && 2nd check");
  fEvents->GetXaxis()->SetBinLabel(12,"Triggered");
  fEvents->GetXaxis()->SetBinLabel(13,"NOT Triggered");
  fEvents->GetXaxis()->SetBinLabel(14,"V0 is LP in MC");
  fEvents->GetYaxis()->SetTitle("Counts"); 
  fOutput->Add(fEvents);

  // Centrality:
  fCentrality = new TH1F("fCentrality","Centrality;Centrality (%);Events",100,0.,100.);
  fOutput->Add(fCentrality);

  // Primary Vertex:
  fPrimaryVertexX = new TH1F("fPrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fOutput->Add(fPrimaryVertexX);
  
  fPrimaryVertexY = new TH1F("fPrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fOutput->Add(fPrimaryVertexY);
  
  fPrimaryVertexZ = new TH1F("fPrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fOutput->Add(fPrimaryVertexZ);


  // ====== Trigger Particle characteristics ====== //

  // MC: Pt Trigger particle vs centrality:
  if(fIsMC){
    fTriggerMCPtCent = new TH3F("fTriggerMCPtCent","Trigger particle MC;p_{T} (GeV/c);#eta;centrality (%)",2*nbinPtLP,pMin,2*ptMaxLP,100,-1.,1.,100,0.,100.);
    fOutput->Add(fTriggerMCPtCent);
  }

  // Pt Trigger particle vs centrality:
  fTriggerPtCent = new TH3F("fTriggerPtCent","Trigger particle;p_{T} (GeV/c);centrality (%);Vertex Z (cm)",nbinPtLP,pMin,ptMaxLP,100,0.,100.,nbinsVtx,-10.,10.);
  fOutput->Add(fTriggerPtCent);

  fNTrigPerEvt = new TH1F("fNTrigPerEvt","Number of Trigger Particles Per Event",21,-0.5,20.5);
  fOutput->Add(fNTrigPerEvt);

  fTriggerWiSPDHit = new TH1F("fTriggerWiSPDHit","Number of Trigger Particles wi SPD Hits",3,0.,3.);
  fOutput->Add(fTriggerWiSPDHit);

  // Phi vs pseudorapidity:
  fTriggerEtaPhi = new TH2F("fTriggerEtaPhi","Trigger particle;#phi (rad);#eta",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fOutput->Add(fTriggerEtaPhi);
  
  // Check if Trigger particle comes from a V0 daughter:
  fCheckTriggerFromV0Daug = 
    new TH1F("fCheckTriggerFromV0Daug","Trigger particle from a V0 daughter;;Counts",4,-0.5,3.5);
  fCheckTriggerFromV0Daug->GetXaxis()->SetTitle("Flag"); 
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(1,"NOT V0 daug");
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(2,"V0 daug");
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(3,"V0 daug & V0 LP");
  fOutput->Add(fCheckTriggerFromV0Daug);
  
  fTriggerComingFromDaug = new TH1F("fTriggerComingFromDaug","Trigger particle from a V0 daughter;p_{T} (GeV/c);Counts",240, 0, 12);
  fOutput->Add(fTriggerComingFromDaug);

  fTriggerIsV0 = new TH1F("fTriggerIsV0","V0 candidate is a LP;p_{T} (GeV/c);Counts",nbinPtLP,pMin,ptMaxLP);
  fOutput->Add(fTriggerIsV0);

  // ------------------- > Comaring properties of this trigger with the daughters
  //   K0s
  fCheckIDTrigPtK0s = new TH3F("fCheckIDTrigPtK0s","K^{0}_{S}",200,-1.,1.,6,-0.5,5.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPtK0s);

  fCheckIDTrigPhiK0s = new TH3F("fCheckIDTrigPhiK0s","K^{0}_{S}",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPhiK0s);

  fCheckIDTrigEtaK0s = new TH3F("fCheckIDTrigEtaK0s","K^{0}_{S}",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigEtaK0s);

  //   Lambda
  fCheckIDTrigPtLambda = new TH3F("fCheckIDTrigPtLambda","#Lambda",200,-1.,1.,6,-0.5,5.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPtLambda);

  fCheckIDTrigPhiLambda  = new TH3F("fCheckIDTrigPhiLambda","#Lambda",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPhiLambda);

  fCheckIDTrigEtaLambda  = new TH3F("fCheckIDTrigEtaLambda","#Lambda",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigEtaLambda);

  //   AntiLambda
  fCheckIDTrigPtAntiLambda = new TH3F("fCheckIDTrigPtAntiLambda","#bar{#Lambda}",200,-1.,1.,6,-0.5,5.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPtAntiLambda);

  fCheckIDTrigPhiAntiLambda  = new TH3F("fCheckIDTrigPhiAntiLambda","#bar{#Lambda}",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigPhiAntiLambda);

  fCheckIDTrigEtaAntiLambda  = new TH3F("fCheckIDTrigEtaAntiLambda","#bar{#Lambda}",200,-1,1.,3,-0.5,2.5,100,1.,11.);
  fOutput->Add(fCheckIDTrigEtaAntiLambda);

  // ====== MC-true and  MC-Association information ====== //
  if(fIsMC){

    fInjectedParticles = new TH1F("fInjectedParticles","Injected particles;;Counts",2,0.,2.);
    fInjectedParticles->GetXaxis()->SetBinLabel(1,"Injected");
    fInjectedParticles->GetXaxis()->SetBinLabel(2,"Natural");
    fOutput->Add(fInjectedParticles);
    
    // K0s MC-true:
    fK0sMCPt       = new TH1F("fK0sMCPt", "K^{0}_{S} MC;p_{T} (GeV/c);Counts",nbins,pMin,pMax);
    fOutput->Add(fK0sMCPt);

    fK0sMCPtRap    = new TH3F("fK0sMCPtRap", "K^{0}_{S} MC;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fK0sMCPtRap);
  
    fK0sMCPtPhiEta    = new TH3F("fK0sMCPtPhiEta", "K^{0}_{S} MC;#phi (rad);#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fK0sMCPtPhiEta);
  
    // K0s MC-Association:
    fK0sAssocPt = 
      new TH1F("fK0sAssocPt","K^{0}_{S} Assoc: L_{T} vs p_{T};p_{T} (GeV/c);Counts",nbins,pMin,pMax);
    fOutput->Add(fK0sAssocPt);

    fK0sAssocPtArm = 
      new TH3F("fK0sAssocPtArm","K^{0}_{S} Assoc: p_{T} vs y vs centrality;p_{T} (GeV/c);y,centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fK0sAssocPtArm);

    fK0sAssocPtRap    = new TH3F("fK0sAssocPtRap","K^{0}_{S} Assoc;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fK0sAssocPtRap);
  
    fK0sAssocPtPhiEta    = new TH3F("fK0sAssocPtPhiEta","K^{0}_{S} Assoc;#phi;#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fK0sAssocPtPhiEta);

    fK0sMCResPhi     = new TH3F("fK0sMCResPhi","K^{0}_{S} Assoc: #phi resolution; #phi_{MC}-#phi_{Rec};p_{T} (GeV/c); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sMCResPhi);

    // Lambda MC-true: 
    fLambdaMCPt = new TH1F("fLambdaMCPt","#Lambda MC;p_{T} (GeV/c);Counts",nbins,pMin,pMax);
    fOutput->Add(fLambdaMCPt);

    fLambdaMCPtRap = new TH3F("fLambdaMCPtRap","#Lambda MC;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fLambdaMCPtRap);
  
    fLambdaMCPtPhiEta = new TH3F("fLambdaMCPtPhiEta","#Lambda MC;#phi (rad);#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fLambdaMCPtPhiEta);
  

    // Lambda MC-Association:
    fLambdaAssocPt = 
      new TH1F("fLambdaAssocPt","#Lambda Assoc: L_{T} vs p_{T};p_{T} (GeV/c);Counts",nbins,pMin,pMax);
    fOutput->Add(fLambdaAssocPt);

    fLambdaAssocPtRap = new TH3F("fLambdaAssocPtRap", "#Lambda Assoc;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fLambdaAssocPtRap);
  
    fLambdaAssocPtPhiEta = new TH3F("fLambdaAssocPtPhiEta", "#Lambda Assoc;#phi (rad);#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fLambdaAssocPtPhiEta);

    fLambdaMCResPhi     = new TH3F("fLambdaMCResPhi","#Lambda Assoc: #phi resolution;#phi_{MC}-#phi_{Rec};p_{T} (GeV/c); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaMCResPhi);
  
    // AntiLambda MC-true: 
    fAntiLambdaMCPt = new TH1F("fAntiLambdaMCPt","#bar{#Lambda} MC;p_{T} (GeV/c);Counts",nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaMCPt);

    fAntiLambdaMCPtRap = new TH3F("fAntiLambdaMCPtRap","#bar{#Lambda} MC;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fAntiLambdaMCPtRap);
  
    fAntiLambdaMCPtPhiEta = new TH3F("fAntiLambdaMCPtPhiEta","#bar{#Lambda} MC;#phi (rad);#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaMCPtPhiEta);

    // AntiLambda MC-Association:
    fAntiLambdaAssocPt = 
      new TH1F("fAntiLambdaAssocPt","#bar{#Lambda} Assoc: L_{T} vs p_{T};p_{T} (GeV/c)",nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaAssocPt);

    fAntiLambdaAssocPtRap = new TH3F("fAntiLambdaAssocPtRap", "#bar{#Lambda} Assoc;p_{T} (GeV/c);y;centrality",nbins,pMin,pMax,30,-1.5,1.5,100,0.,100.);
    fOutput->Add(fAntiLambdaAssocPtRap);
  
    fAntiLambdaAssocPtPhiEta = new TH3F("fAntiLambdaAssocPtPhiEta", "#bar{#Lambda} Assoc;#phi (rad);#eta;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),30,-1.5,1.5,nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaAssocPtPhiEta);

    fAntiLambdaMCResPhi     = new TH3F("fAntiLambdaMCResPhi","#bar{#Lambda} Assoc: #phi resolution;#phi_{MC}-#phi_{Rec};p_{T} (GeV/c); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaMCResPhi);

  }  

  // ====== Reconstruction information in AOD ====== //
  fHistArmenterosPodolanski  =
    new TH3F("fHistArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",
	     100,-1.0,1.0,50,0,0.5,7,-0.5,6.5);
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(1,"K^{0}_{S} Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(2,"K^{0}_{S} Bckg");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(3,"#Lambda Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(4,"#Lambda Bckg");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(5,"#bar{#Lambda} Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(6,"#bar{#Lambda} Bckg");
  fOutput->Add(fHistArmenterosPodolanski);

  fHistArmPodBckg =
    new TH3F("fHistArmPodBckg","Background: Armenteros-Podolanski phase space;#alpha;p_{t} arm",
	     100,-1.0,1.0,50,0,0.5,4,-0.5,3.5);
  fHistArmPodBckg->GetZaxis()->SetBinLabel(1,"K^{0}_{S}: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(2,"#Lambda: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(3,"#bar{#Lambda}: Trig events");
  fOutput->Add(fHistArmPodBckg);
  
  // ****** K0s ******
  fK0sMass = 
    new TH3F("fK0sMass", "K^{0}_{s}: mass vs p_{T}",nbins/2,0.448,0.548,nbins,pMin,pMax,6,-0.5,5.5);
  fK0sMass->GetXaxis()->SetTitle("Mass (GeV/c^2)"); 
  fK0sMass->GetYaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sMass->GetZaxis()->SetBinLabel(1,"All events");
  fK0sMass->GetZaxis()->SetBinLabel(2,"Central events");
  fK0sMass->GetZaxis()->SetBinLabel(3,"Arm-Pod Cut");
  fK0sMass->GetZaxis()->SetBinLabel(4,"Triggered events");
  fK0sMass->GetZaxis()->SetBinLabel(5,"Triggered central events");
  fOutput->Add(fK0sMass);

  fK0sPtvsEta =
    new TH3F("fK0sPtvsEta","K^{0}_{s}: p_{T} vs #eta",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fK0sPtvsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fK0sPtvsEta->GetYaxis()->SetTitle("#eta");
  fK0sPtvsEta->GetZaxis()->SetBinLabel(1,"All events");
  fK0sPtvsEta->GetZaxis()->SetBinLabel(2,"Triggered events");
  fK0sPtvsEta->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fK0sPtvsEta->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fK0sPtvsEta);

  fK0sPtvsRap =
    new TH3F("fK0sPtvsRap","K^{0}_{s}: p_{T} vs y",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fK0sPtvsRap->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fK0sPtvsRap->GetYaxis()->SetTitle("y");
  fK0sPtvsRap->GetZaxis()->SetBinLabel(1,"All events");
  fK0sPtvsRap->GetZaxis()->SetBinLabel(2,"Triggered events");
  fK0sPtvsRap->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fK0sPtvsRap->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fK0sPtvsRap);

  fK0sEtaPhi = 
    new TH2F("fK0sEtaPhi","K^{0}_{s}: #phi vs #eta;#phi (rad);#eta",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fOutput->Add(fK0sEtaPhi);

  fK0sMassPtPhi  = 
    new TH3F("fK0sMassPtPhi","K^{0}_{s}: mass vs pt vs #phi;Mass (GeV/c^2);p_{T} (GeV/c);#phi (rad)",
	     nbins/2,0.448,0.548,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fOutput->Add(fK0sMassPtPhi);


  // Correlations:
  fK0sDCADaugToPrimVtx  
    = new TH3F("fK0sDCADaugToPrimVtx","K^{0}_{S} Bckg: dca daughter vs. p_{T,l};DCA Pos daug (cm);DCA Neg daug (cm);p_{T,l} (GeV/c)",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fOutput->Add(fK0sDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fK0sSpatialRes = new TH3F("fK0sSpatialRes","K^{0}_{S}: Spatial resolution;#Delta#phi (rad);trig-assoc. resolution (cm);dec. length (cm)",
			    20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fK0sSpatialRes);

  char hNameHist[100];
  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fK0sdPhidEtaMC_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fK0sdPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"K^{0}_{S} MC: #Delta#phi vs #Delta#eta vs p_{T,l}",
					    nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					    nbinsdEta,-1.5,1.5,
					    nbinsVtx,-10.,10.);
	fK0sdPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
	fK0sdPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fK0sdPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fK0sdPhidEtaMC[jj*kN1+k]);
      }
  
      // Reconstruction level:
      snprintf(hNameHist,100, "fK0sdPhidEtaPtL_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fK0sdPhidEtaPtL[jj*kN1+k] = new TH3F(hNameHist,"K^{0}_{S}: #Delta#phi vs #Delta#eta vs p_{T,l}",
					   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					   nbinsdEta,-1.5,1.5,
					   nbinsVtx,-10.,10.);
      fK0sdPhidEtaPtL[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fK0sdPhidEtaPtL[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fK0sdPhidEtaPtL[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
      fOutput->Add(fK0sdPhidEtaPtL[jj*kN1+k]);
  
      snprintf(hNameHist,100, "fK0sdPhidEtaPtL_Bckg_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fK0sdPhidEtaPtLBckg[jj*kN1+k] = new TH3F(hNameHist,"K^{0}_{S}: #Delta#phi vs #Delta#eta vs p_{T,l}",
					       nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					       nbinsdEta,-1.5,1.5,
					       nbinsVtx,-10.,10.);
      fK0sdPhidEtaPtLBckg[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fK0sdPhidEtaPtLBckg[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fK0sdPhidEtaPtLBckg[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)");
      fOutput->Add(fK0sdPhidEtaPtLBckg[jj*kN1+k]);

    }
  }

  // Correlations (side-band):
  fK0sBckgDecLength
    = new TH2F("fK0sBckgDecLength","K^{0}_{S} Bckg: c#tau vs. p_{T,l}",
	       100,0.,15.,nbinPtLP,pMin,ptMaxLP);
  fK0sBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fK0sBckgDecLength->GetYaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fK0sBckgDecLength);

  fK0sBckgDCADaugToPrimVtx  
    = new TH3F("fK0sBckgDCADaugToPrimVtx","K^{0}_{S} Bckg: dca daughter vs. p_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fK0sBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fK0sBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fK0sBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fK0sBckgDCADaugToPrimVtx);
  
  fK0sBckgEtaPhi = 
    new TH2F("fK0sBckgEtaPhi","K^{0}_{s} Bckg: #phi vs #eta",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fK0sBckgEtaPhi->GetXaxis()->SetTitle("#phi (rad)"); 
  fK0sBckgEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fK0sBckgEtaPhi);

  fK0sBckgPhiRadio
    = new TH2F("fK0sBckgPhiRadio","K^{0}_{S} Bckg: #phi vs l_{T}",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fK0sBckgPhiRadio->GetXaxis()->SetTitle("#phi (rad)"); 
  fK0sBckgPhiRadio->GetYaxis()->SetTitle("l_{T} (cm)"); 
  fOutput->Add(fK0sBckgPhiRadio);
 
  fK0sBckgDCANegDaugToPrimVtx  
    = new TH2F("fK0sBckgDCANegDaugToPrimVtx","K^{0}_{S} Bckg: dca NegDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fK0sBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fK0sBckgDCANegDaugToPrimVtx);

  fK0sBckgDCAPosDaugToPrimVtx  
    = new TH2F("fK0sBckgDCAPosDaugToPrimVtx","K^{0}_{S} Bckg: dca PosDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fK0sBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fK0sBckgDCAPosDaugToPrimVtx);
        
  fV0MassCascade
    = new TH2F("fV0MassCascade","Cascade Reconstruction wi V0's candiates;Invariant Mass (GeV/c^{2});Cascade type",650, 1.2, 2.5,12,0.5,12.5);
  fOutput->Add(fV0MassCascade);


  // ****** Lambda ******
  fLambdaMass = 
    new TH3F("fLambdaMass","Mass vs p_{T} for \\Lambda",nbins,1.065,1.165,nbins,pMin,pMax,6,-0.5,5.5);
  fLambdaMass->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fLambdaMass->GetYaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaMass->GetZaxis()->SetBinLabel(1,"All events");
  fLambdaMass->GetZaxis()->SetBinLabel(2,"Central events");
  fLambdaMass->GetZaxis()->SetBinLabel(3,"Arm-Pod cut");
  fLambdaMass->GetZaxis()->SetBinLabel(4,"Triggered events");
  fLambdaMass->GetZaxis()->SetBinLabel(5,"Triggered central events");
  fOutput->Add(fLambdaMass);
  
  fLambdaPtvsEta =
    new TH3F("fLambdaPtvsEta","\\Lambda: p_{T} vs #eta",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fLambdaPtvsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fLambdaPtvsEta->GetYaxis()->SetTitle("#eta");
  fK0sPtvsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fLambdaPtvsEta->GetYaxis()->SetTitle("#eta");
  fLambdaPtvsEta->GetZaxis()->SetBinLabel(1,"All events");
  fLambdaPtvsEta->GetZaxis()->SetBinLabel(2,"Triggered events");
  fLambdaPtvsEta->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fLambdaPtvsEta->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fLambdaPtvsEta);

  fLambdaPtvsRap =
    new TH3F("fLambdaPtvsRap","\\Lambda: p_{T} vs y",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fLambdaPtvsRap->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fLambdaPtvsRap->GetYaxis()->SetTitle("y");
  fLambdaPtvsRap->GetYaxis()->SetTitle("#eta");
  fLambdaPtvsRap->GetZaxis()->SetBinLabel(1,"All events");
  fLambdaPtvsRap->GetZaxis()->SetBinLabel(2,"Triggered events");
  fLambdaPtvsRap->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fLambdaPtvsRap->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fLambdaPtvsRap);

  fLambdaEtaPhi =
    new TH2F("fLambdaEtaPhi","#phi vs #eta for #Lambda",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fLambdaEtaPhi->GetXaxis()->SetTitle("#phi (rad)"); 
  fLambdaEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fLambdaEtaPhi);


  fLambdaMassPtPhi  = 
    new TH3F("fLambdaMassPtPhi","#Lambda: mass vs pt vs #phi",
	     nbins,1.065,1.165,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fLambdaMassPtPhi->GetXaxis()->SetTitle("Mass (GeV/c^2)"); 
  fLambdaMassPtPhi->GetYaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaMassPtPhi->GetZaxis()->SetTitle("#phi (rad)"); 
  fOutput->Add(fLambdaMassPtPhi);

  // Correlations:
  fLambdaDCADaugToPrimVtx  
    = new TH3F("fLambdaDCADaugToPrimVtx","#Lambda Bckg: dca daughter vs. p_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fLambdaDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fLambdaDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fLambdaDCADaugToPrimVtx->GetZaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fLambdaDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fLambdaSpatialRes = new TH3F("fLambdaSpatialRes","#Lambda: Spatial resolution;#Delta#phi (rad);trig-assoc. resolution (cm);dec. length (cm)",
			       20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fLambdaSpatialRes);


  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fLambdadPhidEtaMC_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fLambdadPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"#Lambda MC: #Delta#phi vs #Delta#eta vs p_{T,l}",
					       nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					       nbinsdEta,-1.5,1.5,
					       nbinsVtx,-10.,10.);
	fLambdadPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
	fLambdadPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fLambdadPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fLambdadPhidEtaMC[jj*kN1+k]);
      }

      // Reconstruction level:
      snprintf(hNameHist,100, "fLambdadPhidEtaPtL_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fLambdadPhidEtaPtL[jj*kN1+k] = new TH3F(hNameHist,"#Lambda: #Delta#phi vs #Delta#eta vs p_{T,l}",
					      nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					      nbinsdEta,-1.5,1.5,
					      nbinsVtx,-10.,10.);
      fLambdadPhidEtaPtL[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fLambdadPhidEtaPtL[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fLambdadPhidEtaPtL[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)");
      fOutput->Add(fLambdadPhidEtaPtL[jj*kN1+k]);
  
      snprintf(hNameHist,100, "fLambdadPhidEtaPtL_Bckg_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fLambdadPhidEtaPtLBckg[jj*kN1+k] = new TH3F(hNameHist,"#Lambda: #Delta#phi  vs #Delta#eta vs p_{T,l}",
						  nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
						  nbinsdEta,-1.5,1.5,
						  nbinsVtx,-10.,10.);
      fLambdadPhidEtaPtLBckg[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fLambdadPhidEtaPtLBckg[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fLambdadPhidEtaPtLBckg[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)");
      fOutput->Add(fLambdadPhidEtaPtLBckg[jj*kN1+k]);
    }
  }

  // Correlations (side-band):
  fLambdaBckgDecLength
    = new TH2F("fLambdaBckgDecLength","#Lambda Bckg: c#tau vs. p_{T,l}",
	       100,0.,25.,nbinPtLP,pMin,ptMaxLP);
  fLambdaBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fLambdaBckgDecLength->GetYaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fLambdaBckgDecLength);
  
  fLambdaBckgDCADaugToPrimVtx  
    = new TH3F("fLambdaBckgDCADaugToPrimVtx","#Lambda Bckg: dca daughter vs. p_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fLambdaBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fLambdaBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fLambdaBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fLambdaBckgDCADaugToPrimVtx);
  
  fLambdaBckgEtaPhi = 
    new TH2F("fLambdaBckgEtaPhi","#Lambda Bckg: #phi vs #eta",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fLambdaBckgEtaPhi->GetXaxis()->SetTitle("#phi (rad)"); 
  fLambdaBckgEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fLambdaBckgEtaPhi);
    
  fLambdaBckgPhiRadio
    = new TH2F("fLambdaBckgPhiRadio","#Lambda Bckg: #phi vs l_{T}",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fLambdaBckgPhiRadio->GetXaxis()->SetTitle("#phi (rad)"); 
  fLambdaBckgPhiRadio->GetYaxis()->SetTitle("l_{T} (cm)"); 
  fOutput->Add(fLambdaBckgPhiRadio);


  fLambdaBckgDCANegDaugToPrimVtx  
    = new TH2F("fLambdaBckgDCANegDaugToPrimVtx","#Lambda Bckg: dca NegDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fLambdaBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fLambdaBckgDCANegDaugToPrimVtx);


  fLambdaBckgDCAPosDaugToPrimVtx  
    = new TH2F("fLambdaBckgDCAPosDaugToPrimVtx","#Lambda Bckg: dca PosDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fLambdaBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fLambdaBckgDCAPosDaugToPrimVtx);


  // ****** AntiLambda ******
  fAntiLambdaMass = 
    new TH3F("fAntiLambdaMass","Mass vs p_{T} for #bar{#Lambda}",nbins,1.065,1.165,nbins,pMin,pMax,6,-0.5,5.5);
  fAntiLambdaMass->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fAntiLambdaMass->GetYaxis()->SetTitle("p_{T} (GeV/c)"); 
  fAntiLambdaMass->GetZaxis()->SetBinLabel(1,"All events");
  fAntiLambdaMass->GetZaxis()->SetBinLabel(2,"Central events");
  fAntiLambdaMass->GetZaxis()->SetBinLabel(3,"Arm-Pod cut");
  fAntiLambdaMass->GetZaxis()->SetBinLabel(4,"Triggered events");
  fAntiLambdaMass->GetZaxis()->SetBinLabel(5,"Triggered central events");
  fOutput->Add(fAntiLambdaMass);
  
  fAntiLambdaPtvsEta =
    new TH3F("fAntiLambdaPtvsEta","#bar{#Lambda}: p_{T} vs #eta",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fAntiLambdaPtvsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fAntiLambdaPtvsEta->GetYaxis()->SetTitle("#eta");
  fK0sPtvsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fAntiLambdaPtvsEta->GetYaxis()->SetTitle("#eta");
  fAntiLambdaPtvsEta->GetZaxis()->SetBinLabel(1,"All events");
  fAntiLambdaPtvsEta->GetZaxis()->SetBinLabel(2,"Triggered events");
  fAntiLambdaPtvsEta->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fAntiLambdaPtvsEta->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fAntiLambdaPtvsEta);

  fAntiLambdaPtvsRap =
    new TH3F("fAntiLambdaPtvsRap","#bar{#Lambda}: p_{T} vs y",nbins,pMin,pMax,30,-1.5,1.5,4,-0.5,3.5);
  fAntiLambdaPtvsRap->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fAntiLambdaPtvsRap->GetYaxis()->SetTitle("y");
  fAntiLambdaPtvsRap->GetYaxis()->SetTitle("#eta");
  fAntiLambdaPtvsRap->GetZaxis()->SetBinLabel(1,"All events");
  fAntiLambdaPtvsRap->GetZaxis()->SetBinLabel(2,"Triggered events");
  fAntiLambdaPtvsRap->GetZaxis()->SetBinLabel(3,"All ev wi Inv Mass cut");
  fAntiLambdaPtvsRap->GetZaxis()->SetBinLabel(4,"Trig ev wi Inv Mass cut");
  fOutput->Add(fAntiLambdaPtvsRap);

  fAntiLambdaEtaPhi =
    new TH2F("fAntiLambdaEtaPhi","#phi vs #eta for #bar{#Lambda}",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fAntiLambdaEtaPhi->GetXaxis()->SetTitle("#phi (rad)"); 
  fAntiLambdaEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fAntiLambdaEtaPhi);


  fAntiLambdaMassPtPhi  = 
    new TH3F("fAntiLambdaMassPtPhi","#bar{#Lambda}: mass vs pt vs #phi",
	     nbins,1.065,1.165,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fAntiLambdaMassPtPhi->GetXaxis()->SetTitle("Mass (GeV/c^2)"); 
  fAntiLambdaMassPtPhi->GetYaxis()->SetTitle("p_{T} (GeV/c)"); 
  fAntiLambdaMassPtPhi->GetZaxis()->SetTitle("#phi (rad)"); 
  fOutput->Add(fAntiLambdaMassPtPhi);


  // Correlations:
  fAntiLambdaDCADaugToPrimVtx  
    = new TH3F("fAntiLambdaDCADaugToPrimVtx","#bar{#Lambda} Bckg: dca daughter vs. p_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fAntiLambdaDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fAntiLambdaDCADaugToPrimVtx->GetZaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fAntiLambdaDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fAntiLambdaSpatialRes = new TH3F("fAntiLambdaSpatialRes","#bar{#Lambda}: Spatial resolution;#Delta#phi (rad);trig-assoc. resolution (cm);dec. length (cm)",
				   20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fAntiLambdaSpatialRes);

  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fAntiLambdadPhidEtaMC_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fAntiLambdadPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"#bar{#Lambda} MC: #Delta#phi vs #Delta#eta vs p_{T,l}",
						   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
						   nbinsdEta,-1.5,1.5,
						   nbinsVtx,-10.,10.);
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fAntiLambdadPhidEtaMC[jj*kN1+k]);
      }

      // Reconstruction level:
      snprintf(hNameHist,100, "fAntiLambdadPhidEtaPtL_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fAntiLambdadPhidEtaPtL[jj*kN1+k] = new TH3F(hNameHist,"#bar{#Lambda}: #Delta#phi vs #Delta#eta vs p_{T,l}",
						  nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
						  nbinsdEta,-1.5,1.5,
						  nbinsVtx,-10.,10.);
      fAntiLambdadPhidEtaPtL[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fAntiLambdadPhidEtaPtL[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fAntiLambdadPhidEtaPtL[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)");
      fOutput->Add(fAntiLambdadPhidEtaPtL[jj*kN1+k]);
  
      snprintf(hNameHist,100, "fAntiLambdadPhidEtaPtL_Bckg_%.1f_%.1f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
      fAntiLambdadPhidEtaPtLBckg[jj*kN1+k] = new TH3F(hNameHist,"#bar{#Lambda}: #Delta#phi  vs #Delta#eta vs p_{T,l}",
						      nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
						      nbinsdEta,-1.5,1.5,
						      nbinsVtx,-10.,10.);
      fAntiLambdadPhidEtaPtLBckg[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#phi (rad)"); 
      fAntiLambdadPhidEtaPtLBckg[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
      fAntiLambdadPhidEtaPtLBckg[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)");
      fOutput->Add(fAntiLambdadPhidEtaPtLBckg[jj*kN1+k]);

    }
  }

  // Correlations (side-band):
  fAntiLambdaBckgDecLength
    = new TH2F("fAntiLambdaBckgDecLength","#bar{#Lambda} Bckg: c#tau vs. p_{T,l}",
	       100,0.,25.,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fAntiLambdaBckgDecLength->GetYaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fAntiLambdaBckgDecLength);
  
  fAntiLambdaBckgDCADaugToPrimVtx  
    = new TH3F("fAntiLambdaBckgDCADaugToPrimVtx","#bar{#Lambda} Bckg: dca daughter vs. p_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fAntiLambdaBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fAntiLambdaBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("p_{T,l} (GeV/c)"); 
  fOutput->Add(fAntiLambdaBckgDCADaugToPrimVtx);
  
  fAntiLambdaBckgEtaPhi = 
    new TH2F("fAntiLambdaBckgEtaPhi","#bar{#Lambda} Bckg: #phi vs #eta;#phi (rad);l_{T} (cm)",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fOutput->Add(fAntiLambdaBckgEtaPhi);
    
  fAntiLambdaBckgPhiRadio
    = new TH2F("fAntiLambdaBckgPhiRadio","#bar{#Lambda} Bckg: #phi vs l_{T};#phi (rad);l_{T} (cm)",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fOutput->Add(fAntiLambdaBckgPhiRadio);


  fAntiLambdaBckgDCANegDaugToPrimVtx  
    = new TH2F("fAntiLambdaBckgDCANegDaugToPrimVtx","#bar{#Lambda} Bckg: dca NegDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fAntiLambdaBckgDCANegDaugToPrimVtx);


  fAntiLambdaBckgDCAPosDaugToPrimVtx  
    = new TH2F("fAntiLambdaBckgDCAPosDaugToPrimVtx","#bar{#Lambda} Bckg: dca PosDaughter",
	       6,-0.5,5.5,90,0.,3.3);
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"Weak Decay");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Gamma conv.");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Unidentified mother");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Other");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fAntiLambdaBckgDCAPosDaugToPrimVtx);


  // Gamma conversion
  for(Int_t jj=0;jj<kNCent;jj++){
    snprintf(hNameHist,100, "fGammaConversiondPhidEta_Cent_%.0f_%.0f",kBinCent[jj],kBinCent[jj+1]); 
    fGammaConversiondPhidEta[jj] = new TH3F(hNameHist,"Gamma Conversion: #Delta#phi vs #Delta#eta;#Delta#phi (rad);#Delta#eta;Vertex Z (cm)",
					    2*nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					    nbinsdEta,-1.5,1.5,
					    nbinsVtx,-10.,10.);
    fOutput->Add(fGammaConversiondPhidEta[jj]);
  }

  // ============================================================= //

  // K0s in ME:  
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){
      
	snprintf(hNameHist,100,"fK0sdPhidEtaME_%.1f_%.1f_%.0f_%.0f_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);                  
	fK0sdPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH2F(hNameHist,"K^{0}_{S}: #Delta#phi vs #Delta#eta in ME;#Delta#phi (rad);#Delta#eta",
								nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								nbinsdEta,-1.5,1.5);
	fOutputME->Add(fK0sdPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  // Lambda in ME:  
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){

	snprintf(hNameHist,100,"fLambdadPhidEtaME_%.1f_%.1f_%.0lf_%.0lf_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);
	fLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH2F(hNameHist,"#Lambda: #Delta#phi vs #Delta#eta in ME;#Delta#phi (rad);#Delta#eta",
								   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								   nbinsdEta,-1.5,1.5);
	fOutputME->Add(fLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  // AntiLambda in ME:
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){

	snprintf(hNameHist,100,"fAntiLambdadPhidEtaME_%.1f_%.1f_%.0lf_%.0lf_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);
	fAntiLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH2F(hNameHist,"#bar{#Lambda}: #Delta#phi vs #Delta#eta in ME;#Delta#phi (rad);#Delta#eta",
								       nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								       nbinsdEta,-1.5,1.5);
	fOutputME->Add(fAntiLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  
  // ============================================================= //

  if(fDoQA){

    // Quality Assurance K0s:
    // Track PID :
    //    --- signal ---
    fK0sPIDPosDaug  = 
      new TH3F("fK0sPIDPosDaug","K^{0}_{S}: dE/dx Pos. Daug.;TPC Momentum (GeV/c);a.u.;p_{T} V0",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPIDPosDaug);

    fK0sPIDNegDaug  = 
      new TH3F("fK0sPIDNegDaug","K^{0}_{S}: dE/dx Neg. Daug.;TPC Momentum (GeV/c);a.u.;p_{T} V0",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPIDNegDaug);

    //    --- background ---
    fK0sBckgPIDPosDaug  = 
      new TH3F("fK0sBckgPIDPosDaug","K^{0}_{S} Bckg: dE/dx Pos. Daug.;TPC Momentum (GeV/c);a.u.;p_{T} V0",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPIDPosDaug);

    fK0sBckgPIDNegDaug  = 
      new TH3F("fK0sBckgPIDNegDaug","K^{0}_{S} Bckg: dE/dx Neg. Daug.;TPC Momentum (GeV/c);a.u.;p_{T} V0",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPIDNegDaug);

    // Phi Eta
    //     --- signal ---
    fK0sPhiEtaPosDaug = 
      new TH3F("fK0sPhiEtaPosDaug","K^{0}_{S}: #phi vs #eta Pos. Daug.;#phi;#eta;p_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPhiEtaPosDaug);

    fK0sPhiEtaNegDaug  = 
      new TH3F("fK0sPhiEtaNegDaug","K^{0}_{S}: #phi vs #eta Neg. Daug.;#phi;#eta;p_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPhiEtaNegDaug);

    //     --- background ---
    fK0sBckgPhiEtaPosDaug = 
      new TH3F("fK0sBckgPhiEtaPosDaug","K^{0}_{S} Bckg: #phi vs #eta Pos. Daug.;#phi;#eta;p_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPhiEtaPosDaug);

    fK0sBckgPhiEtaNegDaug  = 
      new TH3F("fK0sBckgPhiEtaNegDaug","K^{0}_{S} Bckg: #phi vs #eta Neg. Daug.;#phi;#eta;p_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPhiEtaNegDaug);

    // Distance of closest approach:
    //     --- signal ---
    fK0sDCAPosDaug = 
      new TH2F("fK0sDCAPosDaug","K^{0}_{S}: dca Pos;dca;p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDCAPosDaug);

    fK0sDCANegDaug =  
      new TH2F("fK0sDCANegDaug","K^{0}_{S}: dca Neg;dca;p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDCANegDaug);
    
    //     --- background ---
    fK0sBckgDCAPosDaug = 
      new TH2F("fK0sBckgDCAPosDaug","K^{0}_{S} Bckg: dca Pos;dca;p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDCAPosDaug);

    fK0sBckgDCANegDaug =  
      new TH2F("fK0sBckgDCANegDaug","K^{0}_{S} Bckg: dca Neg;dca;p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDCANegDaug);

    // Difference in momentum:
    //     --- signal ---
    fK0sDifPtPosDaug =  
      new TH2F("fK0sDifPtPosDaug","K^{0}_{S}: dif. p_{T};#Delta p_{T};p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDifPtPosDaug);

    fK0sDifPtNegDaug =  
      new TH2F("fK0sDifPtNegDaug","K^{0}_{S}: dif. p_{T};#Delta p_{T};p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDifPtNegDaug);

    //     --- background ---
    fK0sBckgDifPtPosDaug =  
      new TH2F("fK0sBckgDifPtPosDaug","K^{0}_{S} Bckg: dif. p_{T};#Delta p_{T};p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDifPtPosDaug);

    fK0sBckgDifPtNegDaug =  
      new TH2F("fK0sBckgDifPtNegDaug","K^{0}_{S} Bckg: dif. p_{T};#Delta p_{T};p_{T} V0",90,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDifPtNegDaug);


    // Decay vertex reconstruction:
    //     --- signal ---
    fK0sDecayPos  =  
      new TH3F("fK0sDecayPos","K^{0}_{S}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fK0sDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fK0sDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fK0sDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sDecayPos);

    fK0sDecayVertex  =  
      new TH2F("fK0sDecayVertex","K^{0}_{S}: decay length",100,0.,100.,nbins,pMin,pMax);
    fK0sDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fK0sDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sDecayVertex);

    fK0sDecayVertexZoom  =  
      new TH2F("fK0sDecayVertexZoom","K^{0}_{S}: decay length",20,0.,1.,nbins,pMin,pMax);
    fK0sDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fK0sDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sDecayVertexZoom);

    //     --- background ---
    fK0sBckgDecayPos  =  
      new TH3F("fK0sBckgDecayPos","K^{0}_{S}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fK0sBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fK0sBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fK0sBckgDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sBckgDecayPos);

    fK0sBckgDecayVertex  =  
      new TH2F("fK0sBckgDecayVertex","K^{0}_{S} Bckg: decay vertex",100,0.,100.,nbins,pMin,pMax);
    fK0sBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fK0sBckgDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sBckgDecayVertex);

    fK0sBckgDecayVertexZoom  =  
      new TH2F("fK0sBckgDecayVertexZoom","K^{0}_{S} Bckg: decay length",20,0.,1.,nbins,pMin,pMax);
    fK0sBckgDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fK0sBckgDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sBckgDecayVertexZoom);

    // Cosine of the Pointing Angle:
    //     --- signal ---
    fK0sCPA  =  
      new TH2F("fK0sCPA","K^{0}_{S}: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fK0sCPA->GetXaxis()->SetTitle("cpa"); 
    fK0sCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sCPA);
    //     --- background ---
    fK0sBckgCPA  =  
      new TH2F("fK0sBckgCPA","K^{0}_{S} Bckg: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fK0sBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fK0sBckgCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sBckgCPA);

    // DCA between daughters:
    //     --- signal ---
    fK0sDCAV0Daug  =  
      new TH2F("fK0sDCAV0Daug","K^{0}_{S}: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fK0sDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fK0sDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sDCAV0Daug);
    //     --- background ---
    fK0sBckgDCAV0Daug  =  
      new TH2F("fK0sBckgDCAV0Daug","K^{0}_{S} Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fK0sBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fK0sBckgDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fK0sBckgDCAV0Daug);

    // Number of TPC clusters:
    //     --- signal ---
    fK0sNClustersTPC = 
      new TH3F("fK0sNClustersTPCPos","K^{0}_{S};#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fK0sNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fK0sNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fK0sNClustersTPC);
    //     --- background ---
    fK0sBckgNClustersTPC = 
      new TH3F("fK0sBckgNClustersTPC","K^{0}_{S} Bckg;#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fK0sBckgNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fK0sBckgNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fK0sBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fK0sNClustersITSPos = 
      new TH3F("fK0sNClustersITSPos","K^{0}_{S}: Pos. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sNClustersITSPos);

    fK0sNClustersITSNeg = 
      new TH3F("fK0sNClustersITSNeg","K^{0}_{S}: Neg. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sNClustersITSNeg);
    //     --- background ---
    fK0sBckgNClustersITSPos = 
      new TH3F("fK0sBckgNClustersITSPos","K^{0}_{S} Bckg: Pos. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sBckgNClustersITSPos);

    fK0sBckgNClustersITSNeg = 
      new TH3F("fK0sBckgNClustersITSNeg","K^{0}_{S} Bckg: Neg. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sBckgNClustersITSNeg);
  

    // Quality Assurance Lambda:
    // Track PID :
    //    --- signal ---
    fLambdaPIDPosDaug  = 
      new TH3F("fLambdaPIDPosDaug","#Lambda: dE/dx Pos. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fLambdaPIDPosDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fLambdaPIDPosDaug->GetYaxis()->SetTitle("a.u."); 
    fLambdaPIDPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaPIDPosDaug);

    fLambdaPIDNegDaug  = 
      new TH3F("fLambdaPIDNegDaug","#Lambda: dE/dx Neg. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fLambdaPIDNegDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fLambdaPIDNegDaug->GetYaxis()->SetTitle("a.u."); 
    fLambdaPIDNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaPIDNegDaug);

    //    --- background ---
    fLambdaBckgPIDPosDaug  = 
      new TH3F("fLambdaBckgPIDPosDaug","#Lambda: dE/dx Pos. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fLambdaBckgPIDPosDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fLambdaBckgPIDPosDaug->GetYaxis()->SetTitle("a.u."); 
    fLambdaBckgPIDPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPIDPosDaug);

    fLambdaBckgPIDNegDaug  = 
      new TH3F("fLambdaBckgPIDNegDaug","#Lambda: dE/dx Neg. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fLambdaBckgPIDNegDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fLambdaBckgPIDNegDaug->GetYaxis()->SetTitle("a.u."); 
    fLambdaBckgPIDNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPIDNegDaug);

    // Phi Eta
    //     --- signal ---
    fLambdaPhiEtaPosDaug = 
      new TH3F("fLambdaPhiEtaPosDaug","#Lambda: #phi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaPhiEtaPosDaug->GetXaxis()->SetTitle("#phi"); 
    fLambdaPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaPhiEtaPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaPhiEtaPosDaug);

    fLambdaPhiEtaNegDaug  = 
      new TH3F("fLambdaPhiEtaNegDaug","#Lambda: #phi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaPhiEtaNegDaug->GetXaxis()->SetTitle("#phi"); 
    fLambdaPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaPhiEtaNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaPhiEtaNegDaug);

    //     --- background ---
    fLambdaBckgPhiEtaPosDaug = 
      new TH3F("fLambdaBckgPhiEtaPosDaug","#Lambda: #phi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaBckgPhiEtaPosDaug->GetXaxis()->SetTitle("#phi"); 
    fLambdaBckgPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaBckgPhiEtaPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPhiEtaPosDaug);

    fLambdaBckgPhiEtaNegDaug  = 
      new TH3F("fLambdaBckgPhiEtaNegDaug","#Lambda: #phi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaBckgPhiEtaNegDaug->GetXaxis()->SetTitle("#phi"); 
    fLambdaBckgPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaBckgPhiEtaNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPhiEtaNegDaug);

    // Distance of closest approach
    //     --- signal ---
    fLambdaDCAPosDaug = 
      new TH2F("fLambdaDCAPosDaug","#Lambda: dca Pos",90,0.,3.3,nbins,pMin,pMax);
    fLambdaDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaDCAPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDCAPosDaug);

    fLambdaDCANegDaug =  
      new TH2F("fLambdaDCANegDaug","#Lambda: dca Neg",90,0.,3.3,nbins,pMin,pMax);
    fLambdaDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaDCANegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDCANegDaug);
    
    //     --- background ---
    fLambdaBckgDCAPosDaug = 
      new TH2F("fLambdaBckgDCAPosDaug","#Lambda Bckg: dca Pos",90,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaBckgDCAPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCAPosDaug);

    fLambdaBckgDCANegDaug =  
      new TH2F("fLambdaBckgDCANegDaug","#Lambda Bckg: dca Neg",90,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaBckgDCANegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCANegDaug);

    // Difference in momentum
    //     --- signal ---
    fLambdaDifPtPosDaug =  
      new TH2F("fLambdaDifPtPosDaug","#Lambda: dif. p_{T} Pos",90,0.,3.3,nbins,pMin,pMax);
    fLambdaDifPtPosDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fLambdaDifPtPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDifPtPosDaug);

    fLambdaDifPtNegDaug =  
      new TH2F("fLambdaDifPtNegDaug","#Lambda: dif. p_{T} Neg",90,0.,3.3,nbins,pMin,pMax);
    fLambdaDifPtNegDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fLambdaDifPtNegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDifPtNegDaug);

    //     --- background ---
    fLambdaBckgDifPtPosDaug =  
      new TH2F("fLambdaBckgDifPtPosDaug","#Lambda Bckg: dif. p_{T} Pos.",90,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDifPtPosDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fLambdaBckgDifPtPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDifPtPosDaug);

    fLambdaBckgDifPtNegDaug =  
      new TH2F("fLambdaBckgDifPtNegDaug","#Lambda Bckg: dif. p_{t} Neg",90,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDifPtNegDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fLambdaBckgDifPtNegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDifPtNegDaug);


    // Decay vertex reconstruction
    //     --- signal ---
    fLambdaDecayPos  =  
      new TH3F("fLambdaDecayPos","#Lambda: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fLambdaDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fLambdaDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fLambdaDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDecayPos);

    fLambdaDecayVertex  =  
      new TH2F("fLambdaDecayVertex","#Lambda: decay length",100,0.,100.,nbins,pMin,pMax);
    fLambdaDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDecayVertex);

    fLambdaDecayVertexZoom  =  
      new TH2F("fLambdaDecayVertexZoom","#Lambda: decay length",20,0.,1.,nbins,pMin,pMax);
    fLambdaDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDecayVertexZoom);

    //     --- background ---
    fLambdaBckgDecayPos  =  
      new TH3F("fLambdaBckgDecayPos","#Lambda Bckg: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fLambdaBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fLambdaBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fLambdaBckgDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDecayPos);

    fLambdaBckgDecayVertex  =  
      new TH2F("fLambdaBckgDecayVertex","#Lambda Bckg: decay length",100,0.,100.,nbins,pMin,pMax);
    fLambdaBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaBckgDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDecayVertex);

    fLambdaBckgDecayVertexZoom  =  
      new TH2F("fLambdaBckgDecayVertexZoom","#Lambda Bckg: decay length",20,0.,1.,nbins,pMin,pMax);
    fLambdaBckgDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaBckgDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDecayVertexZoom);

    // Cosine of the Pointing Angle
    //     --- signal ---
    fLambdaCPA  =  
      new TH2F("fLambdaCPA","#Lambda: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fLambdaCPA->GetXaxis()->SetTitle("cpa"); 
    fLambdaCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaCPA);
    //     --- background ---
    fLambdaBckgCPA  =  
      new TH2F("fLambdaBckgCPA","#Lambda Bckg: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fLambdaBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fLambdaBckgCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgCPA);

    // DCA between daughters
    //     --- signal ---
    fLambdaDCAV0Daug  =  
      new TH2F("fLambdaDCAV0Daug","#Lambda: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fLambdaDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fLambdaDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaDCAV0Daug);
    //     --- background ---
    fLambdaBckgDCAV0Daug  =  
      new TH2F("fLambdaBckgDCAV0Daug","#Lambda Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fLambdaBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fLambdaBckgDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCAV0Daug);
  
    // Number of TPC clusters:
    //     --- signal ---
    fLambdaNClustersTPC = 
      new TH3F("fLambdaNClustersTPCPos","#Lambda;#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fLambdaNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fLambdaNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fLambdaNClustersTPC);
    //     --- background ---
    fLambdaBckgNClustersTPC = 
      new TH3F("fLambdaBckgNClustersTPC","#Lambda Bckg;#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fLambdaBckgNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fLambdaBckgNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fLambdaBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fLambdaNClustersITSPos = 
      new TH3F("fLambdaNClustersITSPos","#Lambda: Pos. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaNClustersITSPos);

    fLambdaNClustersITSNeg = 
      new TH3F("fLambdaNClustersITSNeg","#Lambda: Neg. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaNClustersITSNeg);
    //     --- background ---
    fLambdaBckgNClustersITSPos = 
      new TH3F("fLambdaBckgNClustersITSPos","#Lambda Bckg: Pos. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaBckgNClustersITSPos);

    fLambdaBckgNClustersITSNeg = 
      new TH3F("fLambdaBckgNClustersITSNeg","#Lambda Bckg: Neg. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaBckgNClustersITSNeg);



    // Quality Assurance AntiLambda:
    // Track PID :
    //    --- signal ---
    fAntiLambdaPIDPosDaug  = 
      new TH3F("fAntiLambdaPIDPosDaug","#bar{#Lambda}: dE/dx Pos. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fAntiLambdaPIDPosDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fAntiLambdaPIDPosDaug->GetYaxis()->SetTitle("a.u."); 
    fAntiLambdaPIDPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPIDPosDaug);

    fAntiLambdaPIDNegDaug  = 
      new TH3F("fAntiLambdaPIDNegDaug","#bar{#Lambda}: dE/dx Neg. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fAntiLambdaPIDNegDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fAntiLambdaPIDNegDaug->GetYaxis()->SetTitle("a.u."); 
    fAntiLambdaPIDNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPIDNegDaug);

    //    --- background ---
    fAntiLambdaBckgPIDPosDaug  = 
      new TH3F("fAntiLambdaBckgPIDPosDaug","#bar{#Lambda}: dE/dx Pos. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fAntiLambdaBckgPIDPosDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fAntiLambdaBckgPIDPosDaug->GetYaxis()->SetTitle("a.u."); 
    fAntiLambdaBckgPIDPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPIDPosDaug);

    fAntiLambdaBckgPIDNegDaug  = 
      new TH3F("fAntiLambdaBckgPIDNegDaug","#bar{#Lambda}: dE/dx Neg. Daug.",50,0.2,3,50,0.,6.,nbins,pMin,pMax);
    fAntiLambdaBckgPIDNegDaug->GetXaxis()->SetTitle("TPC Momentum (GeV/c)"); 
    fAntiLambdaBckgPIDNegDaug->GetYaxis()->SetTitle("a.u."); 
    fAntiLambdaBckgPIDNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPIDNegDaug);

    // Phi Eta
    //     --- signal ---
    fAntiLambdaPhiEtaPosDaug = 
      new TH3F("fAntiLambdaPhiEtaPosDaug","#bar{#Lambda}: #phi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaPhiEtaPosDaug->GetXaxis()->SetTitle("#phi"); 
    fAntiLambdaPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaPhiEtaPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPhiEtaPosDaug);

    fAntiLambdaPhiEtaNegDaug  = 
      new TH3F("fAntiLambdaPhiEtaNegDaug","#bar{#Lambda}: #phi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaPhiEtaNegDaug->GetXaxis()->SetTitle("#phi"); 
    fAntiLambdaPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaPhiEtaNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPhiEtaNegDaug);

    //     --- background ---
    fAntiLambdaBckgPhiEtaPosDaug = 
      new TH3F("fAntiLambdaBckgPhiEtaPosDaug","#bar{#Lambda}: #phi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaBckgPhiEtaPosDaug->GetXaxis()->SetTitle("#phi"); 
    fAntiLambdaBckgPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaBckgPhiEtaPosDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPhiEtaPosDaug);

    fAntiLambdaBckgPhiEtaNegDaug  = 
      new TH3F("fAntiLambdaBckgPhiEtaNegDaug","#bar{#Lambda}: #phi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaBckgPhiEtaNegDaug->GetXaxis()->SetTitle("#phi"); 
    fAntiLambdaBckgPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaBckgPhiEtaNegDaug->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPhiEtaNegDaug);

    // Distance of closest approach
    //     --- signal ---
    fAntiLambdaDCAPosDaug = 
      new TH2F("fAntiLambdaDCAPosDaug","#bar{#Lambda}: dca Pos",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaDCAPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCAPosDaug);

    fAntiLambdaDCANegDaug =  
      new TH2F("fAntiLambdaDCANegDaug","#bar{#Lambda}: dca Neg",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaDCANegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCANegDaug);
    
    //     --- background ---
    fAntiLambdaBckgDCAPosDaug = 
      new TH2F("fAntiLambdaBckgDCAPosDaug","#bar{#Lambda} Bckg: dca Pos",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaBckgDCAPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCAPosDaug);

    fAntiLambdaBckgDCANegDaug =  
      new TH2F("fAntiLambdaBckgDCANegDaug","#bar{#Lambda} Bckg: dca Neg",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaBckgDCANegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCANegDaug);

    // Difference in momentum
    //     --- signal ---
    fAntiLambdaDifPtPosDaug =  
      new TH2F("fAntiLambdaDifPtPosDaug","#bar{#Lambda}: dif. p_{T} Pos",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDifPtPosDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fAntiLambdaDifPtPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDifPtPosDaug);

    fAntiLambdaDifPtNegDaug =  
      new TH2F("fAntiLambdaDifPtNegDaug","#bar{#Lambda}: dif. p_{T} Neg",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDifPtNegDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fAntiLambdaDifPtNegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDifPtNegDaug);

    //     --- background ---
    fAntiLambdaBckgDifPtPosDaug =  
      new TH2F("fAntiLambdaBckgDifPtPosDaug","#bar{#Lambda} Bckg: dif. p_{T} Pos.",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDifPtPosDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fAntiLambdaBckgDifPtPosDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDifPtPosDaug);

    fAntiLambdaBckgDifPtNegDaug =  
      new TH2F("fAntiLambdaBckgDifPtNegDaug","#bar{#Lambda} Bckg: dif. p_{t} Neg",90,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDifPtNegDaug->GetXaxis()->SetTitle("#Delta p_{T}"); 
    fAntiLambdaBckgDifPtNegDaug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDifPtNegDaug);


    // Decay vertex reconstruction
    //     --- signal ---
    fAntiLambdaDecayPos  =  
      new TH3F("fAntiLambdaDecayPos","#bar{#Lambda}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fAntiLambdaDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fAntiLambdaDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fAntiLambdaDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDecayPos);

    fAntiLambdaDecayVertex  =  
      new TH2F("fAntiLambdaDecayVertex","#bar{#Lambda}: decay length",100,0.,100.,nbins,pMin,pMax);
    fAntiLambdaDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDecayVertex);

    fAntiLambdaDecayVertexZoom  =  
      new TH2F("fAntiLambdaDecayVertexZoom","#bar{#Lambda}: decay length",20,0.,1.,nbins,pMin,pMax);
    fAntiLambdaDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDecayVertexZoom);

    //     --- background ---
    fAntiLambdaBckgDecayPos  =  
      new TH3F("fAntiLambdaBckgDecayPos","#bar{#Lambda} Bckg: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fAntiLambdaBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fAntiLambdaBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fAntiLambdaBckgDecayPos->GetZaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDecayPos);

    fAntiLambdaBckgDecayVertex  =  
      new TH2F("fAntiLambdaBckgDecayVertex","#bar{#Lambda} Bckg: decay length",100,0.,100.,nbins,pMin,pMax);
    fAntiLambdaBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaBckgDecayVertex->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDecayVertex);

    fAntiLambdaBckgDecayVertexZoom  =  
      new TH2F("fAntiLambdaBckgDecayVertexZoom","#bar{#Lambda} Bckg: decay length",20,0.,1.,nbins,pMin,pMax);
    fAntiLambdaBckgDecayVertexZoom->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaBckgDecayVertexZoom->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDecayVertexZoom);

    // Cosine of the Pointing Angle
    //     --- signal ---
    fAntiLambdaCPA  =  
      new TH2F("fAntiLambdaCPA","#bar{#Lambda}: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fAntiLambdaCPA->GetXaxis()->SetTitle("cpa"); 
    fAntiLambdaCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaCPA);
    //     --- background ---
    fAntiLambdaBckgCPA  =  
      new TH2F("fAntiLambdaBckgCPA","#bar{#Lambda} Bckg: cosine of the pointing angle",100,0.9,1.,nbins,pMin,pMax);
    fAntiLambdaBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fAntiLambdaBckgCPA->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgCPA);

    // DCA between daughters
    //     --- signal ---
    fAntiLambdaDCAV0Daug  =  
      new TH2F("fAntiLambdaDCAV0Daug","#bar{#Lambda}: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fAntiLambdaDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fAntiLambdaDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCAV0Daug);
    //     --- background ---
    fAntiLambdaBckgDCAV0Daug  =  
      new TH2F("fAntiLambdaBckgDCAV0Daug","#bar{#Lambda} Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fAntiLambdaBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fAntiLambdaBckgDCAV0Daug->GetYaxis()->SetTitle("p_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCAV0Daug);

    // Number of TPC clusters:
    //     --- signal ---
    fAntiLambdaNClustersTPC = 
      new TH3F("fAntiLambdaNClustersTPCPos","#bar{#Lambda};#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fAntiLambdaNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fAntiLambdaNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fAntiLambdaNClustersTPC);
    //     --- background ---
    fAntiLambdaBckgNClustersTPC = 
      new TH3F("fAntiLambdaBckgNClustersTPC","#bar{#Lambda} Bckg;#phi;Num. TPC Clusters",nbinsPhi,0.,2.*TMath::Pi(),201,-0.5,200.5,3,-0.5,2.5); 
    fAntiLambdaBckgNClustersTPC->GetZaxis()->SetBinLabel(1,"Pos. Daug.");
    fAntiLambdaBckgNClustersTPC->GetZaxis()->SetBinLabel(3,"Neg. Daug.");
    fOutputQA->Add(fAntiLambdaBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fAntiLambdaNClustersITSPos = 
      new TH3F("fAntiLambdaNClustersITSPos","#bar{#Lambda}: Pos. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaNClustersITSPos);

    fAntiLambdaNClustersITSNeg = 
      new TH3F("fAntiLambdaNClustersITSNeg","#bar{#Lambda}: Neg. Daug;#phi;Num. ITS Clusters;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaNClustersITSNeg);
    //     --- background ---
    fAntiLambdaBckgNClustersITSPos = 
      new TH3F("fAntiLambdaBckgNClustersITSPos","#bar{#Lambda} Bckg: Pos. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaBckgNClustersITSPos);

    fAntiLambdaBckgNClustersITSNeg = 
      new TH3F("fAntiLambdaBckgNClustersITSNeg","#bar{#Lambda} Bckg: Neg. Daug;#phi;Num. ITS Clusters;;p_{T} (GeV/c)",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaBckgNClustersITSNeg);

  }

  // ============================================================= //
  
  PostData(1, fOutput);
  PostData(2, fOutputME);
  PostData(3, fOutputQA);
  
}

//___________________________________________________________________________________________

static Int_t VtxBin(Double_t vtx)
{
  // Bin in vertez position Z
  Int_t bin = -1;
  for(Int_t i=0;i<kNVtxZ;i++)
    if ((vtx>=kBinVtxZ[i]) && (vtx<kBinVtxZ[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

static Int_t PtBin(Double_t pt)
{
  // Bin in pt
  Int_t bin = -1;
  for(Int_t i=0;i<kNVtxZ;i++)
    if ((pt>=kPtBinV0[i]) && (pt<kPtBinV0[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

static Int_t CentBin(Double_t cent)
{
  // Bin in pt
  Int_t bin = -1;
  for(Int_t i=0;i<kNCent;i++)
    if ((cent>=kBinCent[i]) && (cent<kBinCent[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

Bool_t AliAnalysisTaskLambdaOverK0sJets::AcceptTrack(AliAODTrack *t) 
{
  // Track criteria for primaries particles 
  if (TMath::Abs(t->Eta())>0.8 )  return kFALSE; 
  if (!(t->TestFilterMask(1<<7))) return kFALSE; 

  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC < 70) return kFALSE;
  
  // Point in the SPD
  Int_t SPDHits = t->HasPointOnITSLayer(0) + t->HasPointOnITSLayer(1);

  // Propagate the global track to the DCA.
  /*
  Double_t PosAtDCA[2] = {-999,-999};
  Double_t covar[3] = {-999,-999,-999};
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  t->PropagateToDCA(vtx,fAOD->GetMagneticField(),100.,PosAtDCA,covar);
  */

  // 5) DCA cut (See R_AA paper).
  //Double_t DCAcutvalue[2];
  //DCAcutvalue[0] = 0.018 + 0.035*TMath::Power(t->Pt(),-1.01);
  //DCAcutvalue[1] = 2.; 
    
  //if( SPDHits && (TMath::Abs(PosAtDCA[0])>DCAcutvalue[0] || TMath::Abs(PosAtDCA[1])>DCAcutvalue[1])  ){ 
  if( SPDHits )
    fTriggerWiSPDHit->Fill(1.5);
  
  return kTRUE;   
}

//___________________________________________________________________________________________

static Bool_t AcceptTrackV0(const AliAODTrack *t) 
{ 
  // Track criteria for daughter particles of V0 candidate 
  if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC < 70) return kFALSE;

  return kTRUE;   
}

//___________________________________________________________________________________________

Bool_t AliAnalysisTaskLambdaOverK0sJets::AcceptV0(AliAODVertex *vtx, const AliAODv0 *v1) 
{ 
  // Selection for accepting V0 candidates 

  if (v1->GetOnFlyStatus()) return kFALSE;
  
  //if (v1->Pt() < pMin) return kFALSE; ***
  
  const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughter(1);
  const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughter(0);
    
  if( !ntrack1 || !ptrack1 ) return kFALSE;
  if( !AcceptTrackV0(ntrack1) ) return kFALSE;
  if( !AcceptTrackV0(ptrack1) ) return kFALSE;
  
  if( ntrack1->Charge() == ptrack1->Charge()) 
    return kFALSE;

  // Daughters: pseudo-rapidity cut
  if ( TMath::Abs(ntrack1->Eta()) > fMaxEtaDaughter  ||
       TMath::Abs(ptrack1->Eta()) > fMaxEtaDaughter  )
    return kFALSE;

  // Daughters: transverse momentum cut
  if ( ( ntrack1->Pt() < fMinPtDaughter ) || 
       ( ptrack1->Pt() < fMinPtDaughter )  ) 
    return kFALSE;
  
  // Daughters: Impact parameter of daughter to prim vtx
  Float_t xy = v1->DcaNegToPrimVertex();
  if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE;
  xy = v1->DcaPosToPrimVertex();
  if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE;

  // Daughters: DCA
  Float_t dca = v1->DcaV0Daughters();
  if (dca>fMaxDCADaughter) return kFALSE;

  // V0: Cosine of the pointing angle
  Float_t cpa=v1->CosPointingAngle(vtx);
  if (cpa<fMinCPA) return kFALSE;

  // V0: Fiducial volume
  Double_t xyz[3]; v1->GetSecondaryVtx(xyz);
  Float_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
  if (r2<5.*5.) return kFALSE;
  if (r2>lMax*lMax) return kFALSE;

  return kTRUE;
}

//___________________________________________________________________________________________

static Float_t dPHI(Float_t phi1, Float_t phi2) 
{ 
  // Calculate the phi difference between two tracks  
  Float_t deltaPhi = phi1 - phi2;
  
  if (deltaPhi<-TMath::PiOver2())    deltaPhi = deltaPhi + 2*(TMath::Pi());
  if (deltaPhi>(3*TMath::PiOver2()))  deltaPhi = deltaPhi - 2*(TMath::Pi());
  return deltaPhi;
}

//___________________________________________________________________________________________

static Float_t MyRapidity(Float_t rE, Float_t rPz)
{ 
  // Local method for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 

//___________________________________________________________________________________________

static Int_t SameTrack(AliAODTrack *trk, const AliAODTrack *daugTrk)
{ 
  // Local method to compaire the momentum between two tracks

  //double const kEpsilon = 0.01;
  Int_t    isSamePt = 0;

  /*
    Float_t p[3];     trk->GetPxPyPz(p);
    Float_t pNegTrk[3]; nTrk->GetPxPyPz(pNegTrk);
    Float_t pPosTrk[3]; pTrk->GetPxPyPz(pPosTrk);
  
     if( (  fabs(p[0]-pNegTrk[0])<kEpsilon && 
    fabs(p[1]-pNegTrk[1])<kEpsilon && 
    fabs(p[2]-pNegTrk[2])<kEpsilon ) 
    isSamePt = 1;
  */
    
  if(  (TMath::Abs(daugTrk->GetID())+1)==(TMath::Abs(trk->GetID()))  )
    isSamePt = 1;
  
  /*
    if(  (TMath::Abs(nTrk->GetID()))==(TMath::Abs(trk->GetID()))  ||
    (TMath::Abs(pTrk->GetID()))==(TMath::Abs(trk->GetID())) )  isSamePt = 1;
  */

  return isSamePt;

}

//___________________________________________________________________________________________

static Float_t SpatialResolution(Float_t p1x,Float_t p1y,Float_t p2x,Float_t p2y,Float_t dist)
{
  // Obtains the spacial resolution between trigger and V0
  // within a distance in (deltaPhi,deltaEta) < 0.1

  Float_t res = -100.;

  res = TMath::Sqrt( p1x*p1x + p1y*p1y )*TMath::Sqrt( p2x*p2x + p2y*p2y );
  res = (p1x*p2x + p1y*p2y)/res;

  res = TMath::ACos(res);
  
  return res = TMath::Sin(res)*dist;
 
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::RecCascade(AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo)
{
  // Local method to reconstruct cascades candidates from the combinations of three tracks
  // The input tracks correspond to the trigger particle and the daughter tracks of the V0 candidate (correlation step)
  // The trigger particle track will be always consider as a possible daughter of the V0 which coming from the Cascade decay.
  // The daughters of the V0 candidates are switched to be the bachelor track for the Cascade reconstruction.

  Float_t lMassBach=0., lPtot2Bach=0., lEBach=0.;
  Float_t lMassLambda=0., lPtot2Lambda=0., lELambda = 0.; 
  Float_t pLambda[3] = {0.,0.,0.};
  Float_t pCascade[3] = {0.,0.,0.};
  Float_t lMassCascade = 0., lPtot2Cascade=0.;

  // Two loops are done to consider the posibility to reconstruct a Xi or an Omega
  for(Int_t i=0;i<2;i++){

    // 0. Check the charge for both tracks: trk1 & trk2. 
    //    Usefull in the Lambda step.
    if( trk1->Charge() == trk2->Charge() ) 
      continue;
   
    // 1. Bachelor: Allocation for the track
    if(i==0) // Xi 
      lMassBach = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
    else if(i==1) //Omega
      lMassBach = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();

    lPtot2Bach = TMath::Power(trkBch->P(),2);

    lEBach = TMath::Sqrt(lPtot2Bach + lMassBach*lMassBach);

    // 2. Lambda: Kinematical properties
    lMassLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
      
    pLambda[0] = trk1->Px() + trk2->Px();
    pLambda[1] = trk1->Py() + trk2->Py();
    pLambda[2] = trk1->Pz() + trk2->Pz();

    lPtot2Lambda = pLambda[0]*pLambda[0] +  pLambda[1]*pLambda[1] +  pLambda[2]*pLambda[2];

    lELambda = TMath::Sqrt(lPtot2Lambda + lMassLambda*lMassLambda);

    // 3. Cascade: Reconstruction
    pCascade[0] = pLambda[0] + trkBch->Px();
    pCascade[1] = pLambda[1] + trkBch->Py();
    pCascade[2] = pLambda[2] + trkBch->Pz();

    lPtot2Cascade = pCascade[0]*pCascade[0] + pCascade[1]*pCascade[1] + pCascade[2]*pCascade[2];

    lMassCascade = TMath::Sqrt( TMath::Power(lEBach+lELambda,2) - lPtot2Cascade );
   
    // 4. Filling histograms
    if( histo.Contains("K0s") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,1);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,3);
    }
    else if( histo.Contains("AntiLambda") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,9);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,11);
    }
    else if( histo.Contains("Lambda") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,5);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,7);
    }

  }
  
}

//___________________________________________________________________________________________
 
void AliAnalysisTaskLambdaOverK0sJets::V0Loop(V0LoopStep_t step, Bool_t isTriggered, Int_t iArray, Int_t idTrig) 
{ 
  // Three options for the 'step' variable:
  // 1) TriggerCheck
  // 2) Reconstruction

  AliAODTrack *trkTrig;
  Float_t  ptTrig  = -100.;
  Float_t  phiTrig = -100.;
  Float_t  etaTrig = -100.; 
  Double_t pTrig[3]; 

  if( (step==kTriggerCheck || isTriggered) && idTrig>=0 ){
    trkTrig = (AliAODTrack*)fAOD->GetTrack(idTrig); 
    ptTrig  = trkTrig->Pt();
    phiTrig = trkTrig->Phi();
    etaTrig = trkTrig->Eta();
    trkTrig->GetPxPyPz(pTrig); 
  }
  
  // *************************************************
  // Centrality selection
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");

  // *************************************************
  // MC Event
  TClonesArray *stackMC = 0x0;
  Float_t mcXv=0., mcYv=0., mcZv=0.;
   
  if(fIsMC){
    TList *lst = fAOD->GetList();
    stackMC = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stackMC) {
      Printf("ERROR: stack not available");
    }

    AliAODMCHeader *mcHdr = 
      (AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
    
    mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
  }
  
  // *************************************************
  // V0 loop - AOD
  const AliAODVertex *vtx=fAOD->GetPrimaryVertex();
  Float_t xv=vtx->GetX(), yv=vtx->GetY(), zv=vtx->GetZ();
  Int_t nV0sTot = fAOD->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nV0sTot; iV0++) {
    
    AliAODv0 *v0=fAOD->GetV0(iV0);
    if (!v0) continue;
    if (!AcceptV0(fAOD->GetPrimaryVertex(),v0)) continue;
    
    const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
    const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);

    Float_t etaNeg = ntrack->Eta();
    Float_t etaPos = ptrack->Eta();
    Float_t lPtNeg = ntrack->Pt();
    Float_t lPtPos = ptrack->Pt();

    // Decay vertex
    Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
    Float_t dx=xyz[0]-xv, dy=xyz[1]-yv;//, dz=xyz[2]-zv;
    // Decay length: 2D & 3D 
    Float_t lt=TMath::Sqrt(dx*dx + dy*dy); 
    //Float_t dl=TMath::Sqrt(dx*dx + dy*dy + dz*dz);  
    // Momentum: 2D & 3D
    Float_t pt=TMath::Sqrt(v0->Pt2V0());
    //Float_t p=v0->P();
    
    Float_t dlK = 0.4977*lt/pt;
    Float_t dlL = 1.1157*lt/pt; 
    /*
      Float_t dlK  = v0->MassK0Short()*dl/p;
      Float_t dlL  = v0->MassLambda()*dl/p;
      Float_t dlAL = v0->MassAntiLambda()*dl/p;
    */

    // ctau
    Bool_t ctK=kTRUE;  if (dlK > fMaxCtau*2.68 || dlK < fMinCtau*2.68) ctK=kFALSE; 
    Bool_t ctL=kTRUE;  if (dlL > fMaxCtau*7.89 || dlL < fMinCtau*7.89) ctL=kFALSE; 
    Bool_t ctAL=kTRUE; if (dlL > fMaxCtau*7.89 || dlL < fMinCtau*7.89) ctAL=kFALSE;    

    // Armenteros variables:
    Float_t lAlphaV0      =  v0->AlphaV0();
    Float_t lPtArmV0      =  v0->PtArmV0();
    
    // MC Association:
    Bool_t lComeFromSigma     = kFALSE; 
    Bool_t lCheckMcK0Short    = kFALSE;
    Bool_t lCheckMcLambda     = kFALSE;
    Bool_t lCheckMcAntiLambda = kFALSE;
	
    Int_t lMCAssocNegDaug = 0;
    Int_t lMCAssocPosDaug = 0;  
    
    // ********* MC - Association *********
    // In case of injected-MC, the correlations might be done with only natural particles 
    Bool_t isNaturalPart = kTRUE;
    if(step==kReconstruction){
      
      if(fIsMC){        
      	if(!stackMC) goto noas;

	Int_t ntrkMC=stackMC->GetEntriesFast();
	
	Int_t nlab = TMath::Abs(ntrack->GetLabel());//** UInt_t
	Int_t plab = TMath::Abs(ptrack->GetLabel());
  
	// To avoid futher problems 
	if ( (nlab<0 || plab<0) ||
	     (nlab>=ntrkMC || plab>=ntrkMC) )
	  goto noas;      

	AliAODMCParticle *nPart=(AliAODMCParticle*)stackMC->UncheckedAt(nlab);
	AliAODMCParticle *pPart=(AliAODMCParticle*)stackMC->UncheckedAt(plab);
	// MC origin of daughters: Primaries?
	if( nPart->IsPhysicalPrimary() ) lMCAssocNegDaug = 1;
	if( pPart->IsPhysicalPrimary() ) lMCAssocPosDaug = 1;

	if(!nPart || !pPart)   goto noas;
	
	if ( TMath::Abs(nPart->Eta()) > fMaxEtaDaughter ||
	     TMath::Abs(pPart->Eta()) > fMaxEtaDaughter )
	  goto noas;
	
	// Daughter momentum cut
	if ( ( nPart->Pt() < fMinPtDaughter ) || 
	     ( pPart->Pt() < fMinPtDaughter )  ) 
	  goto noas;
	
	// ----------------------------------------
	
	Int_t lPDGCodeNegDaughter = nPart->GetPdgCode();
	Int_t lPDGCodePosDaughter = pPart->GetPdgCode();
	
	Int_t ipMother = pPart->GetMother();
	Int_t inMother = nPart->GetMother();
	
	if(inMother<0 || inMother>=ntrkMC) lMCAssocNegDaug = 5;
	if(ipMother<0 || ipMother>=ntrkMC) lMCAssocPosDaug = 5;

	if(inMother<0 || inMother>=ntrkMC) {  goto noas;}
	if(inMother != ipMother) { // did the negative daughter decay ?
	  AliAODMCParticle *negMotherOfMotherPart = (AliAODMCParticle*)stackMC->UncheckedAt(inMother);
	  if (negMotherOfMotherPart->GetMother() != ipMother) 
	    goto noas;
	}
	
	if (ipMother<0 || ipMother>=ntrkMC)
	  goto noas;     
	
	AliAODMCParticle *p0=(AliAODMCParticle*)stackMC->UncheckedAt(ipMother);
	if(!p0) 
	  goto noas; 

	// ----------------------------------------
	
	if ( (ipMother>=fEndOfHijingEvent) && 
	     (fEndOfHijingEvent!=-1)     && 
	     (p0->GetMother()<0) ) 
	  isNaturalPart = kFALSE; 
	else  isNaturalPart = kTRUE; 

	// ----------------------------------------
	
	if(fSeparateInjPart && !isNaturalPart) goto noas;     
	
	Int_t lPDGCodeV0 = p0->GetPdgCode();
	
	// MC origin of daughters:
	//Decay from Weak Decay?
	if( (TMath::Abs(lPDGCodeV0) == kK0Short) || (TMath::Abs(lPDGCodeV0) == kLambda0) || 
	    (TMath::Abs(lPDGCodeV0) == kSigmaMinus) || (TMath::Abs(lPDGCodeV0) == kSigmaPlus) ||
	    (TMath::Abs(lPDGCodeV0) == kSigma0) || (TMath::Abs(lPDGCodeV0) == kXiMinus) ||
	    (TMath::Abs(lPDGCodeV0) == kOmegaMinus) )
	  { lMCAssocNegDaug = 2; 	  lMCAssocPosDaug = 2; }
	// Gamma conversion
	else if( TMath::Abs(lPDGCodeV0) == kGamma )
	  { lMCAssocNegDaug = 3; 	  lMCAssocPosDaug = 3; }
	// Unidentied mother:
	else 
	  { lMCAssocNegDaug = 4; 	  lMCAssocPosDaug = 4; }


	Int_t lIndexMotherOfMother   = p0->GetMother();
	Int_t lPdgcodeMotherOfMother = 0;
	if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
	else {
	  AliAODMCParticle *lMCAODMotherOfMother=(AliAODMCParticle*)stackMC->UncheckedAt(lIndexMotherOfMother);
	  if (!lMCAODMotherOfMother) {lPdgcodeMotherOfMother=0;}
	  lPdgcodeMotherOfMother = lMCAODMotherOfMother->GetPdgCode();
	}
	
	// Daughter momentum cut: ! FIX it in case of AOD ! //MC or REc
	if ( (nPart->Pt()  < fMinPtDaughter ) ||
	     (pPart->Pt()  < fMinPtDaughter ) )
	  goto noas;
	
	if( (lPDGCodeV0 != kK0Short) &&
	    (lPDGCodeV0 != kLambda0) &&
	    (lPDGCodeV0 != kLambda0Bar) ) 
	  goto noas;
	

	lComeFromSigma     = kFALSE; 
	lCheckMcK0Short    = kFALSE;
	lCheckMcLambda     = kFALSE;
	lCheckMcAntiLambda = kFALSE;
	     
	// ----------------------------------------
      
	// K0s
	if( (lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211) &&
	    (inMother==ipMother) && (lPDGCodeV0==310) ) {
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary()  )
	    lCheckMcK0Short  = kTRUE;
	  
	}
	// Lambda
	else if( (lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211)  &&
		 (inMother==ipMother) && (lPDGCodeV0==3122)  ){
	  
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	       ) lComeFromSigma = kTRUE;
	  else lComeFromSigma = kFALSE; 
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() || 
	       ( !(((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() ) 
		 && (lComeFromSigma==kTRUE) )
	       ) lCheckMcLambda  = kTRUE; 
	  
	}
	// AntiLambda
	else if( (lPDGCodePosDaughter==211) && (lPDGCodeNegDaughter==-2212) &&
		 (inMother==ipMother) && (lPDGCodeV0==-3122) ) {
	  
	  
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	       ) lComeFromSigma = kTRUE;
	  else lComeFromSigma = kFALSE;  
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() || 
	       ( (!((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary()) 
		 && (lComeFromSigma==kTRUE) )
	       ) lCheckMcAntiLambda  = kTRUE;
	  
	}
	
	//  ----------------------------------------
	
	if ((p0->Pt())<pMin) goto noas;
	if (TMath::Abs(p0->Y())>fYMax ) goto noas;
	
	Float_t dxAs = mcXv - p0->Xv(),  dyAs = mcYv - p0->Yv(),  dzAs = mcZv - p0->Zv();
	Float_t l = TMath::Sqrt(dxAs*dxAs + dyAs*dyAs + dzAs*dzAs);
	
	dxAs = mcXv - pPart->Xv(); dyAs = mcYv - pPart->Yv();
	//Float_t ltAs = TMath::Sqrt(dxAs*dxAs + dyAs*dyAs);
	Float_t ptAs = p0->Pt();
	Float_t rapAs = p0->Y();
	Float_t etaAs = p0->Eta();
	// phi resolution for V0-reconstruction
	Float_t resPhi = p0->Phi() - v0->Phi();	

	if (l < 0.01) { // Primary V0
	  
	  // K0s:
	  if(ctK && lCheckMcK0Short){ 
	    fK0sAssocPt->Fill(ptAs);
	    fK0sAssocPtRap->Fill(ptAs,rapAs,centrality);
	    if(centrality<10)
	      fK0sAssocPtPhiEta->Fill(p0->Phi(),etaAs,ptAs);
	    
	    if(lPtArmV0 > TMath::Abs(0.2*lAlphaV0) )
	      fK0sAssocPtArm->Fill(ptAs,rapAs,centrality);

	    fK0sMCResPhi->Fill(resPhi,pt,centrality);
	  } // End K0s selection
	  // Lambda:
	  else if(ctL && lCheckMcLambda) {    
	    fLambdaAssocPt->Fill(ptAs);
	    fLambdaAssocPtRap->Fill(ptAs,rapAs,centrality);
	    if(centrality<10)
	      fLambdaAssocPtPhiEta->Fill(p0->Phi(),etaAs,ptAs);

	    fLambdaMCResPhi->Fill(resPhi,pt,centrality);
	  }// End Lambda:
	  // AntiLambda:
	  else if (ctAL && lCheckMcAntiLambda){
	    fAntiLambdaAssocPt->Fill(ptAs);
	    fAntiLambdaAssocPtRap->Fill(ptAs,rapAs,centrality);
	    if(centrality<10)
	      fAntiLambdaAssocPtPhiEta->Fill(p0->Phi(),etaAs,ptAs);

	    fAntiLambdaMCResPhi->Fill(resPhi,pt,centrality);
	  } // End AntiLambda:
	  
	} // End Primary V0 selection
	
	// After the kinematical selection of K0s and Lambdas
	// it might be that the daugthers are not identified through MC Association
	if(lMCAssocNegDaug==0)
	  lMCAssocNegDaug = 5;
	if(lMCAssocPosDaug==0)
	  lMCAssocPosDaug = 5;
	
	
      } // End MC-Association 
      
      // ************************************
      
    }// End Correlation Step
    
  noas:

    Float_t pPos = -100.;
    Float_t pNeg = -100.;
    Float_t dedxPos = -1000.;
    Float_t dedxNeg = -1000.;
    
    //Float_t nsigPosPion   = 0.;
    Float_t nsigPosProton = 0.;
    //Float_t nsigNegPion   = 0.;
    Float_t nsigNegProton = 0.;

    /*  
	if(fUsePID && !fIsMC) {     
	const AliAODPid *pidNeg = ntrack->GetDetPid();
	const AliAODPid *pidPos = ptrack->GetDetPid();
      
	if (pidNeg && pidPos) {
	pPos = pidPos->GetTPCmomentum();
	pNeg = pidNeg->GetTPCmomentum();
	dedxPos = pidPos->GetTPCsignal()/47.; 
	dedxNeg = pidNeg->GetTPCsignal()/47.; 
  

	if(pPos<1.){
	nsigPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kPion));
	nsigPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton));
	}
	if(pNeg<1.){
	nsigNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kPion));
	nsigNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kProton));
	}

	}

	}
    */

    Double_t pNegDaug[3];  ntrack->GetPxPyPz(pNegDaug);                  
    Double_t pPosDaug[3];  ptrack->GetPxPyPz(pPosDaug);

    Float_t dcaNeg = v0->DcaNegToPrimVertex();
    Float_t dcaPos = v0->DcaPosToPrimVertex();

    Float_t phiNeg = ntrack->Phi();
    Float_t phiPos = ptrack->Phi();
    Float_t difPtNeg = pt - lPtNeg;
    Float_t difPtPos = pt - lPtPos;

    Float_t dca   = v0->DcaV0Daughters();
    Float_t cpa   = v0->CosPointingAngle(fAOD->GetPrimaryVertex());
    Float_t lEta  = v0->PseudoRapV0();
    Float_t lPhi  = v0->Phi();
    //lPhi  = ( (lPhi < 0) ? lPhi + 2*TMath::Pi() : lPhi );    

    /*
    // Good regions
    lPhi  = v0->Phi();
    if(lPhi>0. && lPhi<1.8) continue;
    if(lPhi>2.2 && lPhi<3.8) continue;
    if(lPhi>4.2 && lPhi<5.8) continue;
    */

    /*
    // Bad regions
    lPhi  = v0->Phi();
    if(lPhi>1.8 && lPhi<2.2) continue;
    if(lPhi>3.8 && lPhi<4.2) continue;
    if(lPhi>5.8 && lPhi<6.2) continue;
    */

    // Comparing the pt of the trigger particle wrt the v0-candidate's daughter:
    // It is used as well for the side-band subtraction
    Int_t isSameTrkPosDaug = -1;
    Int_t isSameTrkNegDaug = -1;
    if( step==kTriggerCheck ){
      isSameTrkPosDaug = SameTrack(trkTrig,ptrack);
      isSameTrkNegDaug = SameTrack(trkTrig,ntrack);
    }

    // **********************************
    // Disentangle the V0 candidate
    Float_t massK0s = 0., mK0s = 0., sK0s = 0.;
    Float_t massLambda = 0., mLambda = 0., sLambda = 0.;
    Float_t massAntiLambda = 0.;

    Bool_t isCandidate2K0s = kFALSE;
    massK0s = v0->MassK0Short();
    mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    sK0s = 0.0044 + (0.008-0.0044)/(10-1)*(pt - 1.);
    if ( TMath::Abs(mK0s-massK0s) < 3*sK0s )
      isCandidate2K0s = kTRUE;     
    
    Bool_t isCandidate2Lambda = kFALSE;
    massLambda = v0->MassLambda();
    mLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    sLambda=0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
    if (TMath::Abs(mLambda-massLambda) < 3*sLambda)
      isCandidate2Lambda = kTRUE;  
    
    Bool_t isCandidate2LambdaBar = kFALSE;
    massAntiLambda = v0->MassAntiLambda();
    if (TMath::Abs(mLambda-massAntiLambda) < 3*sLambda)
      isCandidate2LambdaBar = kTRUE; 

    // *******************
    //   Gamma conversion
    // *******************

    if(step==kReconstruction)
      if( (TMath::Sqrt(lPtArmV0*lPtArmV0 + lAlphaV0*lAlphaV0) < 0.2)  && isNaturalPart ){
	fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, 2) );
      }

    // *******************
    //   K0s selection
    // *******************
    if (ctK && (TMath::Abs(v0->RapK0Short())<fYMax) &&
        ( lPtArmV0 > TMath::Abs(0.2*lAlphaV0) ) ) {
      
      switch(step) {
      case kTriggerCheck: 

	if (isCandidate2K0s){
	  if(pt>ptTrig){
	    fIsV0LP = 1; 
	    fPtV0LP = pt;
	  }       
	  
	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1 ){
	    Printf("  The LP has the same momentum in X and Y as one of the K0s daughters *** iV0 %d",iV0); 

	    if(fCheckIDTrig){  // Compare properties of daughters nad 
              Float_t difNegP[3];
	      difNegP[0] = pTrig[0]-pNegDaug[0];  difNegP[1] = pTrig[1]-pNegDaug[1]; difNegP[2] = pTrig[2]-pNegDaug[2];
              Float_t difPosP[3]; 
	      difPosP[0] = pTrig[0]-pPosDaug[0];  difPosP[1] = pTrig[1]-pPosDaug[1]; difPosP[2] = pTrig[2]-pPosDaug[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;
	    
	      Printf("  The LP has the same momentum in X and Y as one of the K0s daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		     iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		     TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		     TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		     TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		     TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		     );
	    
	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtK0s->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiK0s->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaK0s->Fill(posDeltaEta,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtK0s->Fill(difPosP[i],i+3,pt); 
		fCheckIDTrigPhiK0s->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaK0s->Fill(negDeltaEta,2.,pt);
	      }
	      
	    } // End check ID


	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0), iArray);
	    
	  } // Close isTrigFromV0daug

	}// End K0s Mass cut

	break; // End K0s selection for TriggerCheck
      case kReconstruction:
	
	fK0sMass->Fill(massK0s,pt,0);
	fK0sPtvsEta->Fill(pt,lEta,0);
	fK0sPtvsRap->Fill(pt,v0->RapK0Short(),0);
	if(centrality<10) fK0sMass->Fill(massK0s,pt,1);
	
	fK0sMassPtPhi->Fill(massK0s,pt,lPhi);


	if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	  fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, 3) );

	// Only for triggered events and in case of MC K0s is not an embeded particle
	if( isTriggered && isNaturalPart ){
	  
	  fK0sMass->Fill(massK0s,pt,3);
	  fK0sPtvsEta->Fill(pt,lEta,1);
	  fK0sPtvsRap->Fill(pt,v0->RapK0Short(),1);
	  
	  if(centrality<10) fK0sMass->Fill(massK0s,pt,4);
	  
	}

	// Invariant Mass cut
	if (TMath::Abs(mK0s-massK0s) < 3*sK0s) {

	  fK0sEtaPhi->Fill(lPhi,lEta);
	  fK0sPtvsEta->Fill(pt,lEta,2);
	  fK0sPtvsRap->Fill(pt,v0->RapK0Short(),2);

	  if(fDoQA){ // Quality Assurance
	    fK0sPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fK0sPIDNegDaug->Fill(pNeg,dedxNeg,pt);
	    
	    fK0sPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fK0sPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	    
	    fK0sDCAPosDaug->Fill(dcaPos,pt);
	    fK0sDCANegDaug->Fill(dcaNeg,pt);
	    
	    fK0sDifPtPosDaug->Fill(difPtPos,pt);
	    fK0sDifPtNegDaug->Fill(difPtNeg,pt);
	    
	    fK0sDecayPos->Fill(dx,dy,pt);
	    fK0sDecayVertex->Fill(lt,pt);
	    
	    if(lt<1.0)
	      fK0sDecayVertexZoom->Fill(lt,pt); //**
	    fK0sCPA->Fill(cpa,pt); //**
	    fK0sDCAV0Daug->Fill(dca,pt); //**

	  }// End QA 
	  
	} // End selection in mass
	if( TMath::Abs(mK0s-massK0s + 6.5*sK0s) < 1.5*sK0s ||
	    TMath::Abs(mK0s-massK0s - 6.5*sK0s) < 1.5*sK0s  ) {

	  fK0sEtaPhi->Fill(lPhi,lEta,-1);
	  
	  if(fDoQA){ // Quality Assurance
	    fK0sBckgPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fK0sBckgPIDNegDaug->Fill(pNeg,dedxNeg,pt);
            
	    fK0sBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fK0sBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	      
	    fK0sBckgDCAPosDaug->Fill(dcaPos,pt);
	    fK0sBckgDCANegDaug->Fill(dcaNeg,pt);
	      
	    fK0sBckgDifPtPosDaug->Fill(difPtPos,pt);
	    fK0sBckgDifPtNegDaug->Fill(difPtNeg,pt);
	      
	    fK0sBckgDecayPos->Fill(dx,dy,pt);
	    fK0sBckgDecayVertex->Fill(lt,pt);
	      
	    if(lt<1.0)
	      fK0sBckgDecayVertexZoom->Fill(lt,pt); //**
	    fK0sBckgCPA->Fill(cpa,pt); //**
	    fK0sBckgDCAV0Daug->Fill(dca,pt); //**
	      
	  } // End QA

	}// End selection in outside the mass cut
	
	break; // End K0s selection for Corrleation
      default:
	Printf( " Selection of 'step' is not set properly");
	break;
	
      }// End switch

    } // End K0s selection

    // *******************
    // Lambda selection
    // *******************
    if ( ctL && (TMath::Abs(v0->RapLambda())<fYMax)  &&
	 (nsigPosProton < fNSigma) ){

      switch(step) {
      case kTriggerCheck: 
	
	if (isCandidate2Lambda){
	  if(pt>ptTrig) {
	    fIsV0LP = 1;
	    fPtV0LP = pt;
	  }

	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1 ){
	    Printf("  The LP has the same momentum in X and Y as one of the Lambda daughters *** iV0 %d",iV0); 

	    if(fCheckIDTrig){  // Compare properties of daughters nad 
	      Float_t difNegP[3];
	      difNegP[0] = pTrig[0]-pNegDaug[0];  difNegP[1] = pTrig[1]-pNegDaug[1]; difNegP[2] = pTrig[2]-pNegDaug[2];
	      Float_t difPosP[3]; 
	      difPosP[0] = pTrig[0]-pPosDaug[0];  difPosP[1] = pTrig[1]-pPosDaug[1]; difPosP[2] = pTrig[2]-pPosDaug[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;
    
	      Printf("  The LP has the same momentum in X and Y as one of the Lambda daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		     iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		     TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		     TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		     TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		     TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		     );
	      
	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtLambda->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiLambda->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaLambda->Fill(posDeltaEta,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtLambda->Fill(difPosP[i],i+3,pt); 
		fCheckIDTrigPhiLambda->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaLambda->Fill(negDeltaEta,2.,pt);
	      }

	    } // End check ID

	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0), iArray);

	  } // Close isTrigFromV0daug

	} // End Lambda Mass cut	
	break; // End Lambda selection for TriggerCheck
      case kReconstruction:
	
	fLambdaMass->Fill(massLambda,pt,0);
	fLambdaPtvsEta->Fill(pt,lEta,0);
	fLambdaPtvsRap->Fill(pt,v0->RapLambda(),0);
	if(centrality<10) fLambdaMass->Fill(massLambda,pt,1);
	
	fLambdaMassPtPhi->Fill(massLambda,pt,lPhi);

	if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	    fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, 4) );

	// Only for triggered events and in case of MC Lambda is not a embeded particle
	if( isTriggered && isNaturalPart ){
	  
	  fLambdaMass->Fill(massLambda,pt,3);
	  fLambdaPtvsEta->Fill(pt,lEta,1);
	  fLambdaPtvsRap->Fill(pt,v0->RapLambda(),1);
	  
	  if(centrality<10)  fLambdaMass->Fill(massLambda,pt,4);
	  
	} 
	// Invariant Mass cut
	if (TMath::Abs(mLambda-massLambda) < 3*sLambda) {

	  fLambdaEtaPhi->Fill(lPhi,lEta);
	  fLambdaPtvsEta->Fill(pt,lEta,2);
	  fLambdaPtvsRap->Fill(pt,v0->RapLambda(),2);
	  
	  if(fDoQA){ // Quality Assurance
	    fLambdaPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fLambdaPIDNegDaug->Fill(pNeg,dedxNeg,pt);
          
	    fLambdaPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fLambdaPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);

	    fLambdaDCAPosDaug->Fill(dcaPos,pt);
	    fLambdaDCANegDaug->Fill(dcaNeg,pt);

	    fLambdaDifPtPosDaug->Fill(difPtPos,pt);
	    fLambdaDifPtNegDaug->Fill(difPtNeg,pt);

	    fLambdaDecayPos->Fill(dx,dy,pt);
	    fLambdaDecayVertex->Fill(lt,pt);

	    if(lt<1.0)
	      fLambdaDecayVertexZoom->Fill(lt,pt); //**
	    fLambdaCPA->Fill(cpa,pt); //**
	    fLambdaDCAV0Daug->Fill(dca,pt); //**

	  } //End QA
	 	  
	} // End selection in mass
	
	if( (TMath::Abs(mLambda-massLambda + 6.5*sLambda) < 1.5*sLambda) ||
	    (TMath::Abs(mLambda-massLambda - 6.5*sLambda) < 1.5*sLambda) ){

	  fLambdaEtaPhi->Fill(lPhi,lEta,-1);
	  
	  if(fDoQA){ // Quality Assurance
	    fLambdaBckgPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fLambdaBckgPIDNegDaug->Fill(pNeg,dedxNeg,pt);
          
	    fLambdaBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fLambdaBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);

	    fLambdaBckgDCAPosDaug->Fill(dcaPos,pt);
	    fLambdaBckgDCANegDaug->Fill(dcaNeg,pt);

	    fLambdaBckgDifPtPosDaug->Fill(difPtPos,pt);
	    fLambdaBckgDifPtNegDaug->Fill(difPtNeg,pt);

	    fLambdaBckgDecayPos->Fill(dx,dy,pt);
	    fLambdaBckgDecayVertex->Fill(lt,pt);

	    if(lt<1.0)
	      fLambdaBckgDecayVertexZoom->Fill(lt,pt); //**
	    fLambdaBckgCPA->Fill(cpa,pt); //**
	    fLambdaBckgDCAV0Daug->Fill(dca,pt); //**

	  }

	}// End selection in outside the mass cut

	break; // End Lambda selection for Correlation
      default:
	Printf(" Selection of 'step' is not set properly");
	break;
	
      }// End switch
      
    } // End Lambda selection

    // *******************
    // AntiLambda selection
    // *******************
    if ( ctAL && (TMath::Abs(v0->RapLambda())<fYMax)  &&
	 (nsigNegProton < fNSigma)  ) {
      
      switch(step) {
      case kTriggerCheck: 
	
	if (isCandidate2LambdaBar){
	  if(pt>ptTrig) {
	    fIsV0LP = 1;
	    fPtV0LP = pt;
	  }
	  
	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1 ){
	    Printf("  The LP has the same momentum in X and Y as one of the AntiLambda daughters *** iV0 %d",iV0); 

	    if(fCheckIDTrig){  // Compare properties of daughters nad 
	      Float_t difNegP[3];
	      difNegP[0] = pTrig[0]-pNegDaug[0];  difNegP[1] = pTrig[1]-pNegDaug[1]; difNegP[2] = pTrig[2]-pNegDaug[2];
	      Float_t difPosP[3]; 
	      difPosP[0] = pTrig[0]-pPosDaug[0];  difPosP[1] = pTrig[1]-pPosDaug[1]; difPosP[2] = pTrig[2]-pPosDaug[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;
    
	      Printf("  The LP has the same momentum in X and Y as one of the AntiLambda daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		     iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		     TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		     TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		     TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		     TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		     );
	      
	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtAntiLambda->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiAntiLambda->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaAntiLambda->Fill(posDeltaEta,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtAntiLambda->Fill(difPosP[i],i+3,pt); 
		fCheckIDTrigPhiAntiLambda->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaAntiLambda->Fill(negDeltaEta,2.,pt);
	      }

	    } // End check ID  

	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0), iArray);

	  }// Close isTrigFromV0daug
	  
	}// End AntiLambda Mass cut
	break; // End AntiLambda selection for CheckTrigger
      case kReconstruction: 
	
	fAntiLambdaMass->Fill(massAntiLambda,pt,0);
	fAntiLambdaPtvsEta->Fill(pt,lEta,0);
	fAntiLambdaPtvsRap->Fill(pt,v0->RapLambda(),0);
	if(centrality<10) fAntiLambdaMass->Fill(massAntiLambda,pt,1);

	fAntiLambdaMassPtPhi->Fill(massAntiLambda,pt,lPhi);
	  
	if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	  fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, 5) );

	// Only for triggered events and in case of MC AntiLambda is not a embeded particle
	if( isTriggered && isNaturalPart ){
	  
	  fAntiLambdaMass->Fill(massAntiLambda,pt,3);
	  fAntiLambdaPtvsEta->Fill(pt,lEta,1);
	  fAntiLambdaPtvsRap->Fill(pt,v0->RapLambda(),1);
	  
	  if(centrality<10)  fAntiLambdaMass->Fill(massAntiLambda,pt,4);
	  
	} 
	// Invariant Mass cut
	if (TMath::Abs(mLambda-massAntiLambda) < 3*sLambda) {
	  fAntiLambdaEtaPhi->Fill(lPhi,lEta);
	  fAntiLambdaPtvsEta->Fill(pt,lEta,2);
	  fAntiLambdaPtvsRap->Fill(pt,v0->RapLambda(),2);

	  if(fDoQA){ // Quality Assurance
	    fAntiLambdaPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fAntiLambdaPIDNegDaug->Fill(pNeg,dedxNeg,pt);
          
	    fAntiLambdaPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fAntiLambdaPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);

	    fAntiLambdaDCAPosDaug->Fill(dcaPos,pt);
	    fAntiLambdaDCANegDaug->Fill(dcaNeg,pt);

	    fAntiLambdaDifPtPosDaug->Fill(difPtPos,pt);
	    fAntiLambdaDifPtNegDaug->Fill(difPtNeg,pt);

	    fAntiLambdaDecayPos->Fill(dx,dy,pt);
	    fAntiLambdaDecayVertex->Fill(lt,pt);

	    if(lt<1.0)
	      fAntiLambdaDecayVertexZoom->Fill(lt,pt); //**
	    fAntiLambdaCPA->Fill(cpa,pt); //**
	    fAntiLambdaDCAV0Daug->Fill(dca,pt); //**

	  } //End QA

	} // End selection in mass
	
	if( (TMath::Abs(mLambda-massAntiLambda + 6.5*sLambda) < 1.5*sLambda) ||
	    (TMath::Abs(mLambda-massAntiLambda - 6.5*sLambda) < 1.5*sLambda) ){

	  fAntiLambdaEtaPhi->Fill(lPhi,lEta,-1);
	  
	  if(fDoQA){ // Quality Assurance
	    fAntiLambdaBckgPIDPosDaug->Fill(pPos,dedxPos,pt);
	    fAntiLambdaBckgPIDNegDaug->Fill(pNeg,dedxNeg,pt);
          
	    fAntiLambdaBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fAntiLambdaBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);

	    fAntiLambdaBckgDCAPosDaug->Fill(dcaPos,pt);
	    fAntiLambdaBckgDCANegDaug->Fill(dcaNeg,pt);

	    fAntiLambdaBckgDifPtPosDaug->Fill(difPtPos,pt);
	    fAntiLambdaBckgDifPtNegDaug->Fill(difPtNeg,pt);

	    fAntiLambdaBckgDecayPos->Fill(dx,dy,pt);
	    fAntiLambdaBckgDecayVertex->Fill(lt,pt);

	    if(lt<1.0)
	      fAntiLambdaBckgDecayVertexZoom->Fill(lt,pt); //**
	    fAntiLambdaBckgCPA->Fill(cpa,pt); //**
	    fAntiLambdaBckgDCAV0Daug->Fill(dca,pt); //**

	  } // End QA
		
	}// End selection in outside the mass cut
	break;
      default:
	Printf( " Selection of 'step' is not set properly");
	break;
      }// End switch
      
    } // End AntiLambda selection
      
  } // End V0 loop
  
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::TriggerParticle() 
{ 
  // Obtain the trigger particles of the event to perform the correlations in phi and eta

  // ----------------------------
  // 1. Trigger particles 
  Int_t nTrk= fAOD->GetNumberOfTracks();
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  Float_t zv=vtx->GetZ();

  for (Int_t i=0; i<nTrk; i++) {
    AliAODTrack *t = fAOD->GetTrack(i);
    if(!AcceptTrack(t)) continue;
    Double_t pt=t->Pt();
    Double_t eta=t->Eta();
  
    if( (pt>fTrigPtMin)  && (pt<fTrigPtMax) &&  (TMath::Abs(eta)<fTrigEtaMax) ) 
      fTriggerParticles->Add( new AliMiniParticle(centrality, zv, i, pt, t->Phi(), eta, 1) );    
  }

  // ----------------------------
  // 2. Checking if the trigger particle 
  // might be a daughter from the V0-candidate
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);
    Int_t id = trig->ID();
    V0Loop(kTriggerCheck,kFALSE,i,id);
  }
    
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::UserExec(Option_t *)
{
  // Main loop for the Analysis

  // Initializing global variables for the correlation studies (mandatory for each event).
  // ---- 1) Trigger Particle: id track
  fIdTrigger  = -1;
  // ---- 2) TriggerCheck: Variables used to crosscheck if trigger particle is a V0 daughter ---- //
  fIsV0LP     = 0;
  fPtV0LP     = -10.;
  fIsSndCheck = 0;

  // Getting AOD Event
  fAOD = (AliAODEvent *)InputEvent();
  fEvents->Fill(0); //event counter  

  if (!fAOD) {
    Printf("ERROR: aod not available");
    return;
  }
  fEvents->Fill(1);
  
  // Physics selection
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *hdr=(AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t maskIsSelected = hdr->IsEventSelected();
  Bool_t isSelected = kFALSE;

  if( fData.Contains("PbPb2010") )
    isSelected = (maskIsSelected & AliVEvent::kMB);
  else if( fData.Contains("PbPb2011") ){
      Bool_t isSelectedCentral   =  maskIsSelected  &  AliVEvent::kCentral;
      Bool_t isSelectedSemiCentral   =  maskIsSelected  &  AliVEvent::kSemiCentral;
      Bool_t isSelectedMB   =  maskIsSelected   &  AliVEvent::kMB;
      if( isSelectedCentral || isSelectedSemiCentral || isSelectedMB ) 	isSelected = kTRUE;
  }

  if(!isSelected) return;
  fEvents->Fill(2);

  // Centrality selection
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");
  fCentrality->Fill(centrality);

  if (!cent->IsEventInCentralityClass(fCentMin,fCentMax,"V0M")) return;
  fEvents->Fill(3);

  // Primary vertex
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  if (vtx->GetNContributors()<3) return;
  fEvents->Fill(4);

  Float_t xv=vtx->GetX(), yv=vtx->GetY(), zv=vtx->GetZ();
  
  fPrimaryVertexX->Fill(xv);
  fPrimaryVertexY->Fill(yv);
  fPrimaryVertexZ->Fill(zv);

  if (TMath::Abs(zv) > 10.) return;   
  fEvents->Fill(5);

  // Getting PID Response
  fPIDResponse = hdr->GetPIDResponse();

  Int_t curVtxBin = VtxBin(zv);
  Int_t curCentBin = CentBin(centrality);
  
  // **********************************************
  // Triggered Particle -  Trigger Particle
  fTriggerParticles = new TObjArray();
  fTriggerParticles->SetOwner(kTRUE);
  TriggerParticle(); 
 
  // V0-candidate is the highest particle in the event:
  if(fIsV0LP) { fEvents->Fill(8); fTriggerIsV0->Fill(fPtV0LP);}

  Int_t NtrigPerEvnt = 0;
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);

    if(trig->WhichCandidate()==0){
      fEvents->Fill(10);
      fTriggerComingFromDaug->Fill(trig->Pt());
      fCheckTriggerFromV0Daug->Fill(1);
      if(fIsV0LP)  fCheckTriggerFromV0Daug->Fill(2);
    }
    else if( trig->WhichCandidate()==1){
      fEvents->Fill(12); 
      fTriggerEtaPhi->Fill(trig->Phi(),trig->Eta());
      fTriggerPtCent->Fill(trig->Pt(),centrality,zv);
      fCheckTriggerFromV0Daug->Fill(0);
      NtrigPerEvnt++;
    }

  }

  fNTrigPerEvt->Fill(NtrigPerEvnt);
  
  // ******************************************
  // Start loop over MC particles
  
  fEndOfHijingEvent = -1;
  TClonesArray *stack = 0x0;
  Float_t mcXv=0., mcYv=0., mcZv=0.;
  
  if(fIsMC) {

    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
    
    AliAODMCHeader *mcHdr = 
      (AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
  
    mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
  
    Int_t nTrkMC = stack->GetEntriesFast();
    // -----------------------------------------
    // --------- Trigger particle --------------
    // -----------------------------------------

    Float_t triggerMCPt   = -1000.;
    Float_t triggerMCPhi  = -1000.;
    Float_t triggerMCEta  = -1000.;
    Bool_t  isTriggeredMC = kFALSE;

    for (Int_t iTrkMC = 0; iTrkMC < nTrkMC; iTrkMC++){
      
      AliAODMCParticle *p0 = (AliAODMCParticle*)stack->At(iTrkMC);
      if(!p0) continue;
      if(TMath::Abs(p0->Eta())>0.8) continue;
      if(p0->Pt()<0.15) continue;

      // ----------------------------------------
      
      // For injected MC: determine where HIJING event ends 
      if (fEndOfHijingEvent==-1) { 
        if ( ( p0->GetStatus() == 21 ) ||
	     ( (p0->GetPdgCode() == 443) &&
	       (p0->GetMother() == -1)   &&
	       (p0->GetDaughter(0) ==  (iTrkMC+1))) ) {
	  fEndOfHijingEvent = iTrkMC; 
        }
      }

      // ----------------------------------------
      
      Int_t isNaturalPart = 1;
      if ( (iTrkMC>=fEndOfHijingEvent) && 
	   (fEndOfHijingEvent!=-1)     && 
	   (p0->GetMother()<0) ) 
	isNaturalPart = 0; 
     
      // ----------------------------------------
      
      Int_t lPdgcodeCurrentPart = p0->GetPdgCode();       
      if ( (lPdgcodeCurrentPart == kK0Short) ||
	   (lPdgcodeCurrentPart == kLambda0) ||
	   (lPdgcodeCurrentPart == kLambda0Bar) ) continue;
      
      // ----------------------------------------
	
      if(isNaturalPart == 0) continue;
      if( !p0->IsPhysicalPrimary() ) continue;
      
      Float_t ptPrim = p0->Pt();
      if(ptPrim>triggerMCPt){

	triggerMCPt  = p0->Pt();
	triggerMCPhi = p0->Phi();
	triggerMCEta = p0->Eta();

      } // End trigger selection

    } // End loop over charged particles


    // ----------------------------------------------------------------
    // Check if Primary particle has larger pt than strange particles 
    // ----------------------------------------------------------------
    fEndOfHijingEvent = -1;
    Bool_t isV0LPMC = kFALSE;
    if(triggerMCPt>1.)
      for (Int_t iTrkMC = 0; iTrkMC < nTrkMC; iTrkMC++){
      
	AliAODMCParticle *p0 = (AliAODMCParticle*)stack->At(iTrkMC);
	if(!p0) continue;
	
	// ----------------------------------------

	Int_t lPdgcodeCurrentPart = p0->GetPdgCode();       
	if ( (lPdgcodeCurrentPart != kK0Short) &&
	     (lPdgcodeCurrentPart != kLambda0) &&
	     (lPdgcodeCurrentPart != kLambda0Bar) ) continue;
	
	// ----------------------------------------
	
	// For injected MC: determine where HIJING event ends 
	if (fEndOfHijingEvent==-1) { 
	  if ( ( p0->GetStatus() == 21 ) ||
	       ( (p0->GetPdgCode() == 443) &&
		 (p0->GetMother() == -1)   &&
		 (p0->GetDaughter(0) ==  (iTrkMC+1))) ) {
	    fEndOfHijingEvent = iTrkMC; 
	  }
	}

	// ----------------------------------------
	
	Int_t isNaturalPart = 1;
	if ( (iTrkMC>=fEndOfHijingEvent) && 
	     (fEndOfHijingEvent!=-1)     && 
	     (p0->GetMother()<0) ) 
	  isNaturalPart = 0; 
	
	if(fSeparateInjPart)
	  if(isNaturalPart == 0) continue;

	// ----------------------------------------
	
	Float_t lRapCurrentPart = MyRapidity(p0->E(),p0->Pz());      
	Float_t lPtCurrentPart  = p0->Pt();
	
	Int_t iCurrentMother = p0->GetMother();       
	AliAODMCParticle *pCurrentMother = (AliAODMCParticle *)stack->At(iCurrentMother);
	Int_t lPdgCurrentMother = 0;    
	if (iCurrentMother == -1) { lPdgCurrentMother = 0;}
	else { lPdgCurrentMother = pCurrentMother->GetPdgCode(); }
	
	Int_t id0  = p0->GetDaughter(0);
	Int_t id1  = p0->GetDaughter(1);
	
	//if ( id0 ==  id1 ) continue;
	if ( (id0 < 0 || id1 < 0) ||
	     (id0 >=nTrkMC  || id1 >= nTrkMC) ) continue;
	
	AliAODMCParticle *pDaughter0 = (AliAODMCParticle *)stack->UncheckedAt(id0);
	AliAODMCParticle *pDaughter1 = (AliAODMCParticle *)stack->UncheckedAt(id1);
	if (!pDaughter0 || !pDaughter1) continue; 
	
	if ( TMath::Abs(pDaughter0->Eta()) > fMaxEtaDaughter ||
	     TMath::Abs(pDaughter1->Eta()) > fMaxEtaDaughter )
	  continue;	

	// Daughter momentum cut: ! FIX it in case of AOD !
	if ( ( pDaughter0->Pt() < fMinPtDaughter ) || 
	     ( pDaughter1->Pt() < fMinPtDaughter )  ) 
	  continue;
	
	if ((p0->Pt())<pMin) continue;  
	if (TMath::Abs(lRapCurrentPart) > fYMax)  continue;
	
	Float_t dx = mcXv-p0->Xv(),  dy = mcYv-p0->Yv(),  dz = mcZv-p0->Zv();
	Float_t l = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	
	//Cut in the 3D-distance of the secondary vertex to primary vertex
	if (l > 0.01) continue; // secondary V0 
	
	// is V0 the LP?
	if(lPtCurrentPart>triggerMCPt)
	  isV0LPMC = kTRUE;	  

      }
    
    // If V0 is not the particle with highest pt
    // then the event is triggered
    if(!isV0LPMC)
      fTriggerMCPtCent->Fill(triggerMCPt,triggerMCEta,centrality);
    else
      fEvents->Fill(13);
    
    if( (triggerMCPt>=fTrigPtMin)  && (triggerMCPt<=fTrigPtMax) &&
	(TMath::Abs(triggerMCEta)<fTrigEtaMax) &&
	(isV0LPMC==kFALSE) ) 
      isTriggeredMC=kTRUE;


    // -----------------------------------------
    // ---------- Strange particles ------------
    // -----------------------------------------
    Float_t dPhiMC = -100.;
    Float_t dEtaMC = -100.;
    fEndOfHijingEvent = -1;
    for (Int_t iTrkMC = 0; iTrkMC < nTrkMC; iTrkMC++){
      
      AliAODMCParticle *p0 = (AliAODMCParticle*)stack->At(iTrkMC);
      if(!p0) continue;
      
      // ----------------------------------------

      Int_t lPdgcodeCurrentPart = p0->GetPdgCode();       
      if ( (lPdgcodeCurrentPart != kK0Short) &&
	   (lPdgcodeCurrentPart != kLambda0) &&
	   (lPdgcodeCurrentPart != kLambda0Bar) ) continue;
      
      // ----------------------------------------

      Int_t isNaturalPart = 1;
      if ( (iTrkMC>=fEndOfHijingEvent) && 
	   (fEndOfHijingEvent!=-1)     && 
	   (p0->GetMother()<0) ) 
	isNaturalPart = 0; 
      
      fInjectedParticles->Fill(isNaturalPart);

      if(fSeparateInjPart && !isNaturalPart) continue;

      // ----------------------------------------

      Float_t lRapCurrentPart = MyRapidity(p0->E(),p0->Pz());      
      Float_t lEtaCurrentPart = p0->Eta();
      Float_t lPhiCurrentPart = p0->Phi();
      Float_t lPtCurrentPart  = p0->Pt();

      Int_t iCurrentMother = p0->GetMother();       
      AliAODMCParticle *pCurrentMother = (AliAODMCParticle *)stack->At(iCurrentMother);
      Int_t lPdgCurrentMother = 0;    
      if (iCurrentMother == -1) { lPdgCurrentMother = 0;}
      else { lPdgCurrentMother = pCurrentMother->GetPdgCode(); }

      Int_t id0  = p0->GetDaughter(0);
      Int_t id1  = p0->GetDaughter(1);
    
      //if ( id0 ==  id1 ) continue;
      if ( (id0 < 0 || id1 < 0) ||
	   (id0 >=nTrkMC  || id1 >= nTrkMC) ) continue;

      AliAODMCParticle *pDaughter0 = (AliAODMCParticle *)stack->UncheckedAt(id0);
      AliAODMCParticle *pDaughter1 = (AliAODMCParticle *)stack->UncheckedAt(id1);
      if (!pDaughter0 || !pDaughter1) continue; 
   
      if ( TMath::Abs(pDaughter0->Eta()) > fMaxEtaDaughter ||
	   TMath::Abs(pDaughter1->Eta()) > fMaxEtaDaughter )
	continue;	

      // Daughter momentum cut: ! FIX it in case of AOD !
      if ( ( pDaughter0->Pt() < fMinPtDaughter ) || 
	   ( pDaughter1->Pt() < fMinPtDaughter )  ) 
	continue;
      
      if ((p0->Pt())<pMin) continue;  
      if (TMath::Abs(lRapCurrentPart) > fYMax)  continue;
    
      Float_t dx = mcXv-p0->Xv(),  dy = mcYv-p0->Yv(),  dz = mcZv-p0->Zv();
      Float_t l = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      
      //Cut in the 3D-distance of the secondary vertex to primary vertex
      if (l > 0.01) continue; // secondary V0 
     
      //Transverse distance to vertex
      dx = mcXv-pDaughter0->Xv(); dy = mcYv-pDaughter0->Yv();
      //Float_t lt=TMath::Sqrt(dx*dx + dy*dy);

      // K0s
      if (lPdgcodeCurrentPart == kK0Short) {

	fK0sMCPt->Fill(lPtCurrentPart);
	fK0sMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality);

	if(centrality<10)
	  fK0sMCPtPhiEta->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);

	// Triggered Event
	if(isTriggeredMC){
	  dPhiMC = dPHI(triggerMCPhi,lPhiCurrentPart);
	  dEtaMC = triggerMCEta - lEtaCurrentPart;
	  
	  // Pt bin
	  for(Int_t k=0;k<kN1;k++)
	    if( (lPtCurrentPart>kPtBinV0[k]) && (lPtCurrentPart<kPtBinV0[k+1]) ){	      
	      fK0sdPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);    
	    } // End pt bin
	  
	} // End trigger event

      } // End K0s selection
      
      // Lambda
      if (lPdgcodeCurrentPart == kLambda0) {

	fLambdaMCPt->Fill(lPtCurrentPart);
	fLambdaMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality);

	if(centrality<10)
	  fLambdaMCPtPhiEta->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);

	// Triggered Event
	if(isTriggeredMC){

	  dPhiMC = dPHI(triggerMCPhi,lPhiCurrentPart);
	  dEtaMC = triggerMCEta - lEtaCurrentPart;
	  
	  // Pt bin
	  for(Int_t k=0;k<kN1;k++)
	    if( (lPtCurrentPart>kPtBinV0[k]) && (lPtCurrentPart<kPtBinV0[k+1]) ){
     	      fLambdadPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);
	    } //End pt bin
	  
	} // End trigger event

      } // End Lambda
      // AntiLambda
      if (lPdgcodeCurrentPart == kLambda0Bar) {

	fAntiLambdaMCPt->Fill(lPtCurrentPart);
	fAntiLambdaMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality);

	if(centrality<10)
	  fAntiLambdaMCPtPhiEta->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);

	// Triggered Event
	if(isTriggeredMC){

	  dPhiMC = dPHI(triggerMCPhi,lPhiCurrentPart);
	  dEtaMC = triggerMCEta - lEtaCurrentPart;
	  
	  // Pt bin
	  for(Int_t k=0;k<kN1;k++)
	    if( (lPtCurrentPart>kPtBinV0[k]) && (lPtCurrentPart<kPtBinV0[k+1]) ){
     	      fAntiLambdadPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);
	    } //End pt bin
	  
	} // End trigger event

      } // End AntiLambda
     
    } // End loop over MC
    
  } // End MC condition


  // *************************************************
  // V0 loop - AOD
  fAssocParticles = new TObjArray(); 
  fAssocParticles->SetOwner(kTRUE);
  if(NtrigPerEvnt>0)
    V0Loop(kReconstruction,kTRUE,-1,-1);
  else 
    V0Loop(kReconstruction,kFALSE,-1,-1);
  

  //-------------------------------------------------------------
  // Correlations
  //-------------------------------------------------------------
  Float_t ptTrig=0., pxTrig=0., pyTrig=0.;
  Float_t massK0s=0., mK0s=0., sK0s=0.;
  Float_t massL=0.,   mL=0.,   sL=0.;
  Float_t massAL=0.; //,  mAL=0.,  sAL=0.;
  Float_t pt=-100., pxAssoc=-1000., pyAssoc=-1000.;
  Float_t lPhi=0., lEta=0.;
  Float_t lAlphaV0=0., lPtArmV0=0, dcaPos=0., dcaNeg=0.;
  Float_t dx=-100., dy=-100., dz=-100., lt=-100., res=-100.;
  Float_t dlK=-100., dlL=-100.;
  Float_t dPhi=-100., dEta=-100., radio=-100.;

  
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);
    if( trig->WhichCandidate() == 0 ) continue;

    AliAODTrack *tTrig = (AliAODTrack*)fAOD->GetTrack(trig->ID());
    ptTrig = tTrig->Pt();  pxTrig = tTrig->Px();  pyTrig = tTrig->Py(); 

    for(Int_t j=0; j<fAssocParticles->GetEntriesFast(); j++){
      AliMiniParticle* trackAssocME = (AliMiniParticle*) (fAssocParticles->At(j));
      AliAODv0 *tAssoc=fAOD->GetV0(trackAssocME->ID());
      const AliAODTrack *ntrack=(AliAODTrack *)tAssoc->GetDaughter(1);
      const AliAODTrack *ptrack=(AliAODTrack *)tAssoc->GetDaughter(0);

      if( SameTrack(tTrig,ntrack) || SameTrack(tTrig,ptrack) )
	continue;

      lPhi = trackAssocME->Phi();
      lEta = trackAssocME->Eta();

      // Correlation in deltaPhi & deltaEta
      dPhi = dPHI(trig->Phi(),lPhi);
      dEta = trig->Eta() - lEta;
      radio    = TMath::Sqrt(dPhi*dPhi + dEta*dEta);
     
      // Armenteros variables: 
      lAlphaV0      =  tAssoc->AlphaV0();
      lPtArmV0      =  tAssoc->PtArmV0();

      // 2D momentum
      pt = trackAssocME->Pt(); pxAssoc = tAssoc->Px(); pyAssoc = tAssoc->Py(); 
      // Decay vertex
      Double_t xyz[3]; tAssoc->GetSecondaryVtx(xyz);
      dx=xyz[0]-xv; dy=xyz[1]-yv; dz=xyz[2]-zv;
      // Decay length: 2D 
      lt=TMath::Sqrt(dx*dx + dy*dy); 
      // Spatial resolution trigger-V0 point decay
      res = SpatialResolution(pxTrig,pyTrig,pxAssoc,pyAssoc,lt);
      // Ctau
      dlK = 0.4977*lt/pt;
      dlL = 1.1157*lt/pt; 

      Int_t binPtv0 = PtBin( pt );
      if(binPtv0==-1) continue;

      Int_t lMCAssocNegDaug=0, lMCAssocPosDaug=0;

      // *******************
      //   Gamma conversion
      // *******************
      if( trackAssocME->WhichCandidate() == 2 )
	fGammaConversiondPhidEta[curCentBin]->Fill(dPhi,dEta,zv);

      // *******************
      //   K0s selection
      // *******************
      if( trackAssocME->WhichCandidate() == 3 ){

	massK0s = tAssoc->MassK0Short();
	mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
	sK0s = 0.0044 + (0.008-0.0044)/(10-1)*(pt - 1.);
	
	// ==== Correlations K0s invariant mass peak ==== //
	if (TMath::Abs(mK0s-massK0s) < 3*sK0s) {
	  
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,0);
	    fK0sPtvsEta->Fill(pt,lEta,3);
	    fK0sPtvsRap->Fill(pt,tAssoc->RapK0Short(),3);
	
	    // Pt bin & centrality
	    fK0sdPhidEtaPtL[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);

	    if(radio<0.1)
	      fK0sSpatialRes->Fill(dPhi,res,lt);
	    if(radio < 0.4){
	      fK0sDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);	    
	      RecCascade(tTrig,ntrack,ptrack,"K0s");
	      RecCascade(tTrig,ptrack,ntrack,"K0s");	
	    }
	
		
	}
	// ==== Correlations K0s background ==== //
	if( TMath::Abs(mK0s-massK0s + 6.5*sK0s) < 1.5*sK0s ||
	    TMath::Abs(mK0s-massK0s - 6.5*sK0s) < 1.5*sK0s  ) {
	  
	  fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,1);
	  
	  // Pt bin & centrality
	  fK0sdPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	
	  if(radio < 0.4){ // Under the correlation peak
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,0);
	    fK0sBckgDecLength->Fill(dlK,ptTrig);
	    fK0sBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    fK0sBckgEtaPhi->Fill(lPhi,lEta);
	    fK0sBckgPhiRadio->Fill(lPhi,lt);

	    // MC Association of daughter particles 
	    for(Int_t ii=0;ii<5;ii++){
	      if(lMCAssocNegDaug == ii){
		fK0sBckgDCANegDaugToPrimVtx->Fill(ii,dcaNeg);
		//if eta selection for the trigger particle
	      }
	      if(lMCAssocPosDaug == ii){
		fK0sBckgDCAPosDaugToPrimVtx->Fill(ii,dcaPos);
		//if eta selection for the trigger particle
	      }

	    }

	    //RecCascade(trkTrig,ntrack,ptrack,"K0s");
	    //RecCascade(trkTrig,ptrack,ntrack,"K0s");

	  }// End selection in the correlation peak
		
	} // End background selection
	
      } // End K0s selection

      // *******************
      // Lambda selection
      // *******************
      if( trackAssocME->WhichCandidate() == 4 ){
	massL = tAssoc->MassLambda();
	mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	sL = 0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
	
	// ==== Correlations Lambda invariant mass peak ==== //
	if (TMath::Abs(mL-massL) < 3*sL) {

	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,2);
	    fLambdaPtvsEta->Fill(pt,lEta,3);
	    fLambdaPtvsRap->Fill(pt,tAssoc->RapLambda(),3);
		
	    // Pt bin & centrality
	    fLambdadPhidEtaPtL[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
		
	    if(radio<0.1)
	      fLambdaSpatialRes->Fill(dPhi,res,lt);
	    if(radio < 0.4){
	      fLambdaDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	      RecCascade(tTrig,ntrack,ptrack,"Lambda");
	      RecCascade(tTrig,ptrack,ntrack,"Lambda");
	    }
	    
	} // End mass peak selection
	// ==== Correlations Lambda background ==== //
	if( TMath::Abs(mL-massL + 6.5*sL) < 1.5*sL ||
	    TMath::Abs(mL-massL - 6.5*sL) < 1.5*sL ) {
	    
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,3);
	  
	    // Pt bin & centrality
	    fLambdadPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	
	    if(radio < 0.4){ // Under the peak
	      fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,1);
	      fLambdaBckgDecLength->Fill(dlL,ptTrig);
	      fLambdaBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	      fLambdaBckgEtaPhi->Fill(lPhi,lEta);
	      fLambdaBckgPhiRadio->Fill(lPhi,lt);
		  
	      // MC Association of daughter particles 
	      for(Int_t ii=0;ii<5;ii++){
		if(lMCAssocNegDaug == ii)
		  fLambdaBckgDCANegDaugToPrimVtx->Fill(ii,dcaNeg);
		    
		if(lMCAssocPosDaug == ii)
		  fLambdaBckgDCAPosDaugToPrimVtx->Fill(ii,dcaPos);
	      }

	      //RecCascade(trkTrig,ntrack,ptrack,"Lambda");
	      //RecCascade(trkTrig,ptrack,ntrack,"Lambda");

	    }// End selection in the correlation peak
		
	} // End bacground selection
	
      }// End Lambda selection
       // *******************
      // AntiLambda selection
      // *******************
      if( trackAssocME->WhichCandidate() == 5 ){
	massAL = tAssoc->MassAntiLambda();
	mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	sL = 0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
	
	// ==== Correlations AntiLambda invariant mass peak ==== //
	if (TMath::Abs(mL-massAL) < 3*sL) {
	
	  fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,4);
	  fAntiLambdaPtvsEta->Fill(pt,lEta,3);
	  fAntiLambdaPtvsRap->Fill(pt,tAssoc->RapLambda(),3);
	  
	  // Pt bin & centrality
	  fAntiLambdadPhidEtaPtL[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);

	  if(radio<0.1)
	    fAntiLambdaSpatialRes->Fill(dPhi,res,lt);	      
	  if(radio < 0.4){
	    fAntiLambdaDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    RecCascade(tTrig,ntrack,ptrack,"AntiLambda");
	    RecCascade(tTrig,ptrack,ntrack,"AntiLambda");
	  }
	      
	} // End AntiLambda mass peak
	// ==== Correlations AntiLambda background ==== //
	if( (TMath::Abs(mL-massAL + 6.5*sL) < 1.5*sL) ||
	    (TMath::Abs(mL-massAL - 6.5*sL) < 1.5*sL) ){
	    
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,5);
	    
	    // Pt bin & centrality
	    fAntiLambdadPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	
	    if(radio < 0.4){ // Under the peak
	      fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,2);
	      fAntiLambdaBckgDecLength->Fill(dlL,ptTrig);
	      fAntiLambdaBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	      fAntiLambdaBckgEtaPhi->Fill(lPhi,lEta);
	      fAntiLambdaBckgPhiRadio->Fill(lPhi,lt);
		  
	      // MC Association of daughter particles 
	      for(Int_t ii=0;ii<5;ii++){
		if(lMCAssocNegDaug == ii)
		  fAntiLambdaBckgDCANegDaugToPrimVtx->Fill(ii,dcaNeg);
		    
		if(lMCAssocPosDaug == ii)
		  fAntiLambdaBckgDCAPosDaugToPrimVtx->Fill(ii,dcaPos);
	      }

	      //RecCascade(trkTrig,ntrack,ptrack,"AntiLambda");
	      //RecCascade(trkTrig,ptrack,ntrack,"AntiLambda");

	    }// End selection in the correlation peak
		
	}// End AntiLambda background

      } // End AntiLambda selection

    } // End loop over associated particles
    
  } // End loop over trigger particles
  
  //-------------------------------------------------------------
  // Mixing
  //-------------------------------------------------------------
  
  TList *evMixList = fMEList[curCentBin*kNVtxZ+curVtxBin];
  Int_t nMixed = evMixList->GetSize(); 
 
  if( nMixed>0 && fAssocParticles->GetEntriesFast() >= 0 ){
    
    for(Int_t ii=0; ii<nMixed; ii++){     
      
      AliMiniParticle* trackTriggerME = (AliMiniParticle*) (evMixList->At(ii));
      Double_t phiTrigME = trackTriggerME->Phi();
      Double_t etaTrigME = trackTriggerME->Eta();
           
      for(Int_t j=0; j<fAssocParticles->GetEntriesFast(); j++){
          
	AliMiniParticle* trackAssocME = (AliMiniParticle*) (fAssocParticles->At(j));
	if( CentBin(trackTriggerME->Centrality()) != CentBin(trackAssocME->Centrality()) ) continue;
	if( VtxBin(trackTriggerME->VtxZ()) != VtxBin(trackAssocME->VtxZ()) ) continue;
	if( trackAssocME->WhichCandidate() ==  2 ) continue;

	AliAODv0 *tAssoc=fAOD->GetV0(trackAssocME->ID());
	pt = tAssoc->Pt();

	Bool_t IsSelected = kFALSE;
	// K0s
	if( trackAssocME->WhichCandidate() == 3 ){
	  massK0s = tAssoc->MassK0Short();
	  mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
	  sK0s = 0.0044 + (0.008-0.0044)/(10-1)*(pt - 1.);
	  if (TMath::Abs(mK0s-massK0s) < 3*sK0s) IsSelected = kTRUE;
	}
	// Lambda
	if( trackAssocME->WhichCandidate() == 4 ){
	  massL = tAssoc->MassLambda();
	  mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	  sL = 0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
	  if (TMath::Abs(mL-massL) < 3*sL) IsSelected = kTRUE;
	}
	// AntiLambda
	if( trackAssocME->WhichCandidate() == 5 ){
	  massAL = tAssoc->MassAntiLambda();
	  mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	  sL = 0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
	  if (TMath::Abs(mL-massAL) < 3*sL) IsSelected = kTRUE;
	}

	if(!IsSelected) continue;

	Double_t phiAssocME = trackAssocME->Phi();
	Double_t etaAssocME = trackAssocME->Eta();
	 
	Double_t deltaPhi = dPHI(phiTrigME,phiAssocME);
	Double_t deltaEta = etaTrigME - etaAssocME;

	Int_t binPtv0 = PtBin( trackAssocME->Pt() );
	if(binPtv0==-1) continue;
    
	if( trackAssocME->WhichCandidate() == 3 ) {
	fK0sdPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta);}
	else if( trackAssocME->WhichCandidate() == 4 )
	  fLambdadPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta);
	else if( trackAssocME->WhichCandidate() == 5 )
	  fAntiLambdadPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta);
	             
      }
      
    }
    
  }
 
  //--------------------------------------------------------
  //Add the current event to the list of events for mixing
  //--------------------------------------------------------  
  
  //Add current  event to buffer and Remove redundant events 
  if(fTriggerParticles->GetEntriesFast()>=0){
    
    for(Int_t ii=0; ii<(fTriggerParticles->GetEntriesFast()); ii++){
      AliMiniParticle* trkTrig = (AliMiniParticle*) fTriggerParticles->At(ii);
      //cout << trkTrig->Pt() << "          " << ii << endl;
    
      evMixList->AddFirst(trkTrig);
      if(evMixList->GetSize() >= nMaxEvMix) {
	AliMiniParticle *tmp = (AliMiniParticle*) (evMixList->Last()) ;
	evMixList->RemoveLast();
	delete tmp;
      }
      
    }// End loop over fTriggerParticles

  }// End adding trigger particles to buffers
  
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  fOutput=(TList*)GetOutputData(1);
  fOutputME=(TList*)GetOutputData(2);
  fOutputQA=(TList*)GetOutputData(3);

  if (fOutput || fOutputME || fOutputQA) {

    if(fOutput)
      Printf("\n\t *** DONE: fOutput available *** \n");
    if(fOutputME)
      Printf("\n\t *** DONE: fOutputME available *** \n");
    if (fOutputQA)
      Printf("\n\t *** DONE: fOutputQA available *** \n");
  }
  if (!fOutput || !fOutputME || !fOutputQA) {

    if(!fOutput)
      Printf("\n\t *** ERROR: fOutput not available *** \n");
    if(!fOutputME) 
      Printf("\n\t *** ERROR: fOutputME available *** \n");
    if(!fOutputQA)
      Printf("\n\t *** ERROR: fOutputQA not available  *** \n");  
  }

  
  return;

}
