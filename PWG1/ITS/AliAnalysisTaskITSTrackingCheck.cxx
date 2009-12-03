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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTask to extract from ESD tracks the information
// on ITS tracking efficiency and resolutions.
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>  
#include <TCanvas.h>
#include <TParticle.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliVertexerTracks.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliTrackPointArray.h"
#include "../ITS/AliITSRecPoint.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"

#include "AliGenEventHeader.h" 
#include "AliAnalysisTaskITSTrackingCheck.h"


ClassImp(AliAnalysisTaskITSTrackingCheck)

//________________________________________________________________________
AliAnalysisTaskITSTrackingCheck::AliAnalysisTaskITSTrackingCheck(const char *name) : 
AliAnalysisTask(name, "ITSTrackingCheckTask"), 
fReadMC(kFALSE),
fReadRPLabels(kFALSE),
fFillNtuples(kFALSE),
fESD(0), 
fESDfriend(0),
fOutput(0), 
fHistNtracks(0),
fHistNclsITSMI(0),
fHistNclsITSSA(0),
fHistNclsITSSAInAcc(0),
fHistClusterMapITSMI(0),
fHistClusterMapITSMIok(0),
fHistClusterMapITSMIbad(0),
fHistClusterMapITSMIskipped(0),
fHistClusterMapITSMIoutinz(0),
fHistClusterMapITSMInorefit(0),
fHistClusterMapITSMInocls(0),
fHistClusterMapITSMIokoutinzbad(0),
fHistClusterMapITSSA(0),
fHistClusterMapITSSAok(0),
fHistClusterMapITSSAbad(0),
fHistClusterMapITSSAskipped(0),
fHistClusterMapITSSAoutinz(0),
fHistClusterMapITSSAnorefit(0),
fHistClusterMapITSSAnocls(0),
fHistClusterMapITSSAokoutinzbad(0),
fHistClusterMapITSSAInAcc(0),
fHistClusterMapITSSAokInAcc(0),
fHistClusterMapITSSAbadInAcc(0),
fHistClusterMapITSSAskippedInAcc(0),
fHistClusterMapITSSAoutinzInAcc(0),
fHistClusterMapITSSAnorefitInAcc(0),
fHistClusterMapITSSAnoclsInAcc(0),
fHistClusterMapITSSAokoutinzbadInAcc(0),
fHistClusterMapModuleITSSAokInAcc(0),
fHistClusterMapModuleITSSAbadInAcc(0),
fHistClusterMapModuleITSSAnoclsInAcc(0),
fHistPhiTPC(0),
fHistPtTPC(0),
fHistPtITSMI2(0),
fHistPtITSMI3(0),
fHistPtITSMI4(0),
fHistPtITSMI5(0),
fHistPtITSMI6(0),
fHistPtITSMISPD(0),
fNtupleESDTracks(0),
fNtupleITSAlignExtra(0),
fNtupleITSAlignSPDTracklets(0)
{
  // Constructor

  for(Int_t i=0; i<10; i++) fCountsPerPtBin[i]=0;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  //My private output
}
//________________________________________________________________________
AliAnalysisTaskITSTrackingCheck::~AliAnalysisTaskITSTrackingCheck()
{
  // Destructor

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read

    tree->SetBranchStatus("ESDfriend*", 1);
    tree->SetBranchAddress("ESDfriend.",&fESDfriend);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else {
      fESD = esdH->GetEvent();

    }
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::CreateOutputObjects()
{
  // Create histograms
  // Called once

  for(Int_t i=0; i<10; i++) fCountsPerPtBin[i]=0;

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();

  fHistNtracks = new TH1F("fHistNtracks", "N ESD tracks; N tracks; Events",5000, -0.5, 4999.5);
  fHistNtracks->Sumw2();
  fHistNtracks->SetMinimum(0);
  fOutput->Add(fHistNtracks);

  fHistNclsITSMI = new TH1F("fHistNclsITSMI", "N ITS clusters per track (MI); N clusters; Counts",7, -0.5, 6.5);
  fHistNclsITSMI->Sumw2();
  fHistNclsITSMI->SetMinimum(0);
  fOutput->Add(fHistNclsITSMI);

  fHistNclsITSSAInAcc = new TH1F("fHistNclsITSSAInAcc", "N ITS clusters per track (SA); N clusters; Counts",7, -0.5, 6.5);
  fHistNclsITSSAInAcc->Sumw2();
  fHistNclsITSSAInAcc->SetMinimum(0);
  fOutput->Add(fHistNclsITSSAInAcc);  

  fHistNclsITSSA = new TH1F("fHistNclsITSSA", "N ITS clusters per track (SA); N clusters; Counts",7, -0.5, 6.5);
  fHistNclsITSSA->Sumw2();
  fHistNclsITSSA->SetMinimum(0);
  fOutput->Add(fHistNclsITSSA);  

  fHistClusterMapITSMI = new TH1F("fHistClusterMapITSMI", "N tracks with point on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMI->Sumw2();
  fHistClusterMapITSMI->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMI);
  
  fHistClusterMapITSSA = new TH1F("fHistClusterMapITSSA", "N tracks with point on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSA->Sumw2();
  fHistClusterMapITSSA->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSA);

  fHistClusterMapITSSAInAcc = new TH1F("fHistClusterMapITSSAInAcc", "N tracks with point on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAInAcc->Sumw2();
  fHistClusterMapITSSAInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAInAcc);

  fHistClusterMapITSMIok = new TH1F("fHistClusterMapITSMIok", "N tracks with ok on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMIok->Sumw2();
  fHistClusterMapITSMIok->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIok);
  
  fHistClusterMapITSSAokInAcc = new TH1F("fHistClusterMapITSSAokInAcc", "N tracks with ok on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAokInAcc->Sumw2();
  fHistClusterMapITSSAokInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokInAcc);

  fHistClusterMapModuleITSSAokInAcc = new TH1F("fHistClusterMapModuleITSSAokInAcc", "N tracks with ok on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAokInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAokInAcc);

  fHistClusterMapITSSAok = new TH1F("fHistClusterMapITSSAok", "N tracks with ok on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAok->Sumw2();
  fHistClusterMapITSSAok->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAok);

  fHistClusterMapITSMIbad = new TH1F("fHistClusterMapITSMIbad", "N tracks with bad on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMIbad->Sumw2();
  fHistClusterMapITSMIbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIbad);
  
  fHistClusterMapITSSAbadInAcc = new TH1F("fHistClusterMapITSSAbadInAcc", "N tracks with bad on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAbadInAcc->Sumw2();
  fHistClusterMapITSSAbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAbadInAcc);

  fHistClusterMapModuleITSSAbadInAcc = new TH1F("fHistClusterMapModuleITSSAbadInAcc", "N tracks with bad on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAbadInAcc);

  fHistClusterMapITSSAbad = new TH1F("fHistClusterMapITSSAbad", "N tracks with bad on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAbad->Sumw2();
  fHistClusterMapITSSAbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAbad);

  fHistClusterMapITSMIskipped = new TH1F("fHistClusterMapITSMIskipped", "N tracks with skip on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMIskipped->Sumw2();
  fHistClusterMapITSMIskipped->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIskipped);
  
  fHistClusterMapITSSAskippedInAcc = new TH1F("fHistClusterMapITSSAskippedInAcc", "N tracks with skip on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAskippedInAcc->Sumw2();
  fHistClusterMapITSSAskippedInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAskippedInAcc);

  fHistClusterMapITSSAskipped = new TH1F("fHistClusterMapITSSAskipped", "N tracks with skip on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAskipped->Sumw2();
  fHistClusterMapITSSAskipped->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAskipped);

  fHistClusterMapITSMIoutinz = new TH1F("fHistClusterMapITSMIoutinz", "N tracks out in z on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMIoutinz->Sumw2();
  fHistClusterMapITSMIoutinz->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIoutinz);
  
  fHistClusterMapITSSAoutinzInAcc = new TH1F("fHistClusterMapITSSAoutinzInAcc", "N tracks with out in z on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAoutinzInAcc->Sumw2();
  fHistClusterMapITSSAoutinzInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAoutinzInAcc);

  fHistClusterMapITSSAoutinz = new TH1F("fHistClusterMapITSSAoutinz", "N tracks with out in z on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAoutinz->Sumw2();
  fHistClusterMapITSSAoutinz->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAoutinz);

  fHistClusterMapITSSAokoutinzbad = new TH1F("fHistClusterMapITSSAokoutinzbad", "N tracks with cluster or bad zone or out in z (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAokoutinzbad->Sumw2();
  fHistClusterMapITSSAokoutinzbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokoutinzbad);

  fHistClusterMapITSMIokoutinzbad = new TH1F("fHistClusterMapITSMIokoutinzbad", "N tracks with cluster or bad zone or out in z (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMIokoutinzbad->Sumw2();
  fHistClusterMapITSMIokoutinzbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIokoutinzbad);

  fHistClusterMapITSSAokoutinzbadInAcc = new TH1F("fHistClusterMapITSSAokoutinzbadInAcc", "N tracks with cluster or bad zone or out in z (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAokoutinzbadInAcc->Sumw2();
  fHistClusterMapITSSAokoutinzbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokoutinzbadInAcc);

  fHistClusterMapITSMInorefit = new TH1F("fHistClusterMapITSMInorefit", "N tracks with norefit on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMInorefit->Sumw2();
  fHistClusterMapITSMInorefit->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMInorefit);
  
  fHistClusterMapITSSAnorefitInAcc = new TH1F("fHistClusterMapITSSAnorefitInAcc", "N tracks with norefit on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAnorefitInAcc->Sumw2();
  fHistClusterMapITSSAnorefitInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnorefitInAcc);

  fHistClusterMapITSSAnorefit = new TH1F("fHistClusterMapITSSAnorefit", "N tracks with norefit on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAnorefit->Sumw2();
  fHistClusterMapITSSAnorefit->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnorefit);

  fHistClusterMapITSMInocls = new TH1F("fHistClusterMapITSMInocls", "N tracks with nocls on Layer (MI); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSMInocls->Sumw2();
  fHistClusterMapITSMInocls->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMInocls);
  
  fHistClusterMapITSSAnoclsInAcc = new TH1F("fHistClusterMapITSSAnoclsInAcc", "N tracks with nocls on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAnoclsInAcc->Sumw2();
  fHistClusterMapITSSAnoclsInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnoclsInAcc);
  
  fHistClusterMapModuleITSSAnoclsInAcc = new TH1F("fHistClusterMapModuleITSSAnoclsInAcc", "N tracks with nocls on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAnoclsInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAnoclsInAcc);

  fHistClusterMapITSSAnocls = new TH1F("fHistClusterMapITSSAnocls", "N tracks with nocls on Layer (SA); Layer; N tracks",6, -0.5, 5.5);
  fHistClusterMapITSSAnocls->Sumw2();
  fHistClusterMapITSSAnocls->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnocls);
  
  fHistPhiTPC = new TH1F("fHistPhiTPC","Azimuthal distribution of TPC tracks; #phi; N tracks",100, -3.1415, 3.1415);
  fHistPhiTPC->Sumw2();
  fHistPhiTPC->SetMinimum(0);
  fOutput->Add(fHistPhiTPC);
  
  fHistPtTPC = new TH1F("fHistPtTPC","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtTPC->Sumw2();
  fHistPtTPC->SetMinimum(0);
  fOutput->Add(fHistPtTPC);
  
  fHistPtITSMI6 = new TH1F("fHistPtITSMI6","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMI6->Sumw2();
  fHistPtITSMI6->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6);
  
  fHistPtITSMI5 = new TH1F("fHistPtITSMI5","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMI5->Sumw2();
  fHistPtITSMI5->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5);
  
  fHistPtITSMI4 = new TH1F("fHistPtITSMI4","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMI4->Sumw2();
  fHistPtITSMI4->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4);
  
  fHistPtITSMI3 = new TH1F("fHistPtITSMI3","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMI3->Sumw2();
  fHistPtITSMI3->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3);
  
  fHistPtITSMI2 = new TH1F("fHistPtITSMI2","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMI2->Sumw2();
  fHistPtITSMI2->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2);
  
  fHistPtITSMISPD = new TH1F("fHistPtITSMISPD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",80, 0, 40);
  fHistPtITSMISPD->Sumw2();
  fHistPtITSMISPD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPD);
  


  // ntuples
  //
  fNtupleESDTracks = new TNtuple("fNtupleESDTracks","tracks","pt:eta:phi:d0:z0:sigmad0:sigmaz0:ptMC:pdgMC:d0MC:d0MCv:z0MCv:sigmad0MCv:sigmaz0MCv:ITSflag");  
  fOutput->Add(fNtupleESDTracks);

  fNtupleITSAlignExtra = new TNtuple("fNtupleITSAlignExtra","ITS alignment checks: extra clusters","layer:x:y:z:dxy:dz:xloc:zloc:npoints");  
  fOutput->Add(fNtupleITSAlignExtra);

  fNtupleITSAlignSPDTracklets = new TNtuple("fNtupleITSAlignSPDTracklets","ITS alignment checks: SPD tracklets wrt SPD vertex","phi:theta:z:dxy:dz");  
  fOutput->Add(fNtupleITSAlignSPDTracklets);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  //if(fESD->GetEventType()!=7) return;

  const AliESDVertex *vertexESD = fESD->GetPrimaryVertexTracks();

  // ***********  MC info ***************
  TArrayF mcVertex(3);
  mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
  Float_t dNchdy=-999.;

  TParticle *part=0;
  AliESDVertex *vertexMC=0;
  AliStack *stack=0;
  if (fReadMC) {
    AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    
    AliHeader* header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    AliGenEventHeader* genHeader = header->GenEventHeader();
    genHeader->PrimaryVertex(mcVertex);


    Int_t ngenpart = (Int_t)stack->GetNtrack();
    //printf("# generated particles = %d\n",ngenpart);
    dNchdy=0;
    for(Int_t ip=0; ip<ngenpart; ip++) {
      part = (TParticle*)stack->Particle(ip);
      // keep only electrons, muons, pions, kaons and protons
      Int_t apdg = TMath::Abs(part->GetPdgCode());
      if(apdg!=11 && apdg!=13 && apdg!=211 && apdg!=321 && apdg!=2212) continue;      
      // reject secondaries
      if(TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]))>0.0010) continue;
      // reject incoming protons
      Double_t energy  = part->Energy();
      if(energy>900.) continue;
      Double_t pz = part->Pz();
      Double_t y = 0.5*TMath::Log((energy+pz+1.e-13)/(energy-pz+1.e-13));
      if(TMath::Abs(y)<1.0) dNchdy += 0.5; // count 1/2 of particles in |y|<1
    }
    //printf("# primary particles = %7.1f\n",dNchdy);
  } 
  // ***********  MC info ***************
  Double_t mcVtxPos[3]={mcVertex[0],mcVertex[1],mcVertex[2]},mcVtxSigma[3]={0,0,0};
  vertexMC = new AliESDVertex(mcVtxPos,mcVtxSigma);

  // ***********  ESD friends ***********
  fESD->SetESDfriend(fESDfriend); //Attach the friend to the ESD
  // ***********  ESD friends ***********
  if(!fESDfriend) printf("no ESD friend\n");

  //
  
  /*  
  // **********  Trigger *****************
  ULong64_t triggerMask;
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);
  
  triggerMask=fESD->GetTriggerMask();
  // MB1: SPDFO || V0L || V0R
  Bool_t eventTriggered = (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right))); 
  //MB2: GFO && V0R
  //triggerMask & spdFO && ((triggerMask&v0left) || (triggerMask&v0right))
  // ************ Trigger ******************
  if(!eventTriggered) return;
  */

  // SPD vertex
  const AliESDVertex *spdv=fESD->GetPrimaryVertexSPD();
 
  
  Int_t ntracks = fESD->GetNumberOfTracks();
  //printf("Tracks # = %d\n",fESD->GetNumberOfTracks());

  fHistNtracks->Fill(ntracks);
  // Post the data already here
  PostData(0, fOutput);

  Int_t idet,status; Float_t xloc,zloc;
  
  // loop on tracks
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *track = fESD->GetTrack(itr);

    // remove kink daughters
    if(track->GetKinkIndex(0)>0) continue;


    Bool_t itsrefit=kFALSE,tpcin=kFALSE,itsfindableAcc=kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit)) itsrefit=kTRUE;
    if ((track->GetStatus() & AliESDtrack::kTPCin)) tpcin=kTRUE;


    //if(TMath::Abs(track->GetD(0,0,0))>1) continue;

    Int_t trkLabel = TMath::Abs(track->GetLabel());
    Int_t nclsITS = track->GetNcls(0);

    Bool_t outInZ=kFALSE;

    for(Int_t layer=0; layer<6; layer++) {
      track->GetITSModuleIndexInfo(layer,idet,status,xloc,zloc);
      if(layer>=2) idet+=240; // add n SPD modules
      if(layer>=4) idet+=260; // add n SPD modules
      if(status==4) outInZ=kTRUE;
      if(tpcin) {
	if(status==1) fHistClusterMapITSMIok->Fill(layer);
	if(status==2) fHistClusterMapITSMIbad->Fill(layer);
	if(status==3) fHistClusterMapITSMIskipped->Fill(layer);
	if(status==4) fHistClusterMapITSMIoutinz->Fill(layer);
	if(status==5) fHistClusterMapITSMInocls->Fill(layer);
	if(status==6) fHistClusterMapITSMInorefit->Fill(layer);
	if(status==1 || status==2 || status==4) fHistClusterMapITSMIokoutinzbad->Fill(layer);
      } else {
	if(status==1) fHistClusterMapITSSAok->Fill(layer);
	if(status==2) fHistClusterMapITSSAbad->Fill(layer);
	if(status==3) fHistClusterMapITSSAskipped->Fill(layer);
	if(status==4) fHistClusterMapITSSAoutinz->Fill(layer);
	if(status==5) fHistClusterMapITSSAnocls->Fill(layer);
	if(status==6) fHistClusterMapITSSAnorefit->Fill(layer);
	if(status==1 || status==2 || status==4) fHistClusterMapITSSAokoutinzbad->Fill(layer);
	if(status==1 && !outInZ) {fHistClusterMapITSSAokInAcc->Fill(layer);fHistClusterMapModuleITSSAokInAcc->Fill(idet);}
	if(status==2 && !outInZ) {fHistClusterMapITSSAbadInAcc->Fill(layer);fHistClusterMapModuleITSSAbadInAcc->Fill(idet);}
	if(status==3 && !outInZ) fHistClusterMapITSSAskippedInAcc->Fill(layer);
	if(status==4 && !outInZ) fHistClusterMapITSSAoutinzInAcc->Fill(layer);
	if(status==5 && !outInZ) {fHistClusterMapITSSAnoclsInAcc->Fill(layer);fHistClusterMapModuleITSSAnoclsInAcc->Fill(idet);}
	if(status==6 && !outInZ) fHistClusterMapITSSAnorefitInAcc->Fill(layer);
	if((status==1 || status==2 || status==4) && !outInZ) fHistClusterMapITSSAokoutinzbadInAcc->Fill(layer);
      }
      if(TESTBIT(track->GetITSClusterMap(),layer)) {
	if(tpcin) {
	  fHistClusterMapITSMI->Fill(layer);
	} else {
	  fHistClusterMapITSSA->Fill(layer);
	  if(!outInZ) fHistClusterMapITSSAInAcc->Fill(layer);
	}
      }
    }  

    // TPC track in ITS acceptance
    if(tpcin && TMath::Abs(track->Eta())<0.9 && track->GetNcls(1)>50) {
      itsfindableAcc=kTRUE;
      fHistPtTPC->Fill(track->Pt());  
      fHistPhiTPC->Fill(track->Phi());  
    }
  
    if(itsfindableAcc) {
      if(nclsITS==6) fHistPtITSMI6->Fill(track->Pt());
      if(nclsITS==5) fHistPtITSMI5->Fill(track->Pt());
      if(nclsITS==4) fHistPtITSMI4->Fill(track->Pt());
      if(nclsITS==3) fHistPtITSMI3->Fill(track->Pt());
      if(nclsITS==2) fHistPtITSMI2->Fill(track->Pt());
      if(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))
	fHistPtITSMISPD->Fill(track->Pt());
    }


    if(tpcin) {
      fHistNclsITSMI->Fill(nclsITS);
    } else {
      fHistNclsITSSA->Fill(nclsITS);
      if(!outInZ) fHistNclsITSSAInAcc->Fill(nclsITS);
    }


    Int_t iITSflag=0; //ITSflag takes the value 0 if the track has no cluster assigned in the SPDs, 1 (2) if one cluster is assigned in SPD1(2), 3 if two clusters are present. Then the same adding 10,20 or 30 for SDD and 100,200 or 300 for SSD

    if(track->HasPointOnITSLayer(0)) iITSflag+=1;
    if(track->HasPointOnITSLayer(1)) iITSflag+=2;
    if(track->HasPointOnITSLayer(2)) iITSflag+=10;
    if(track->HasPointOnITSLayer(3)) iITSflag+=20;
    if(track->HasPointOnITSLayer(4)) iITSflag+=100;
    if(track->HasPointOnITSLayer(5)) iITSflag+=200;

    if(iITSflag==333 && track->GetNcls(0)<6) 
      printf(" ERROR %d   %d\n",track->GetNcls(0),track->GetLabel());
    
    // number of associated ITS clusters
    iITSflag += 1000*track->GetNcls(0);
     
    // number of associated TPC clusters
    iITSflag += 100000*track->GetNcls(1);
     
    // if MC info and is available
    // write the number of ITS clusters produced by this track
    Int_t nITSclsMC=0;
    if(fReadMC && fReadRPLabels) {
      nITSclsMC = NumberOfITSClustersMC(trkLabel);
      if(nITSclsMC>=0) iITSflag += 10000*nITSclsMC;    
    }

    if(!vertexESD) return;
    if(!(vertexESD->GetStatus())) return;

    // impact parameter to VertexTracks
    Double_t d0z0[2],covd0z0[3];
    track->PropagateToDCA(vertexESD,fESD->GetMagneticField(),100.,d0z0,covd0z0);
    if(covd0z0[0]<0. || covd0z0[2]<0.) continue;

    // if MC info is available: get particle properties
    Float_t ptMC=-999.,pdgMC=-999.,d0MC=-999.;
    Double_t d0z0MCv[2]={-999.,-999.},covd0z0MCv[3]={1.,1.,1.};
    if(fReadMC) {
      part = (TParticle*)stack->Particle(trkLabel);
      ptMC=part->Pt();
      pdgMC=part->GetPdgCode();
      d0MC=ParticleImpParMC(part,vertexMC,0.1*fESD->GetMagneticField());
      track->PropagateToDCA(vertexMC,fESD->GetMagneticField(),100.,d0z0MCv,covd0z0MCv);
      if(covd0z0MCv[0]<0. || covd0z0MCv[2]<0.) continue;
      // flag fake tracks
      if(track->GetLabel()<0) iITSflag *= -1;
    }

    Double_t sigmad0MCv=TMath::Sqrt(covd0z0MCv[0]);
    if(!itsrefit) sigmad0MCv *= -1.;

    // fill ntuple with track properties
    if(fFillNtuples && SelectPt(track->Pt())) {
      Float_t fillArray[19]={track->Pt(),track->Eta(),track->Phi(),d0z0[0],d0z0[1],TMath::Sqrt(covd0z0[0]),TMath::Sqrt(covd0z0[2]),ptMC,pdgMC,d0MC,d0z0MCv[0],d0z0MCv[1],sigmad0MCv,TMath::Sqrt(covd0z0MCv[2]),(Float_t)iITSflag};
      fNtupleESDTracks->Fill(fillArray);
    }

    //---------------------------------------------    
    // AliTrackPoints: alignment checks
    // 
    if(!fFillNtuples) continue;
    if(!fESDfriend) continue;


    const AliTrackPointArray *array = track->GetTrackPointArray();
    if(!array) continue;
    AliTrackPoint point;
    Int_t pointOnLayer[6]={0,0,0,0,0,0};
    Int_t indexAssociated[6]={-1,-1,-1,-1,-1,-1},indexExtra=-1;
    Bool_t extra=kFALSE;
    Int_t layerId,layerExtra=-1;
    for(Int_t ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      Float_t r = TMath::Sqrt(point.GetX()*point.GetX()+point.GetY()*point.GetY());

      if(r>3 && r<6) {
	layerId = 0;
      } else if(r>6 && r<8) {
	layerId = 1;
      } else if(r>8 && r<18) {
	layerId = 2;
      } else if(r>18 && r<25) {
	layerId = 3;
      } else if(r>25 && r<38) {
	layerId = 4;
      } else if(r>38 && r<50) {
	layerId = 5;
      } else {
	layerId=100;
      }

      // only ITS points
      if(layerId>5) continue;

      if(!point.IsExtra()) {
	pointOnLayer[layerId]++;
	indexAssociated[layerId]=ipt;
      } else {
	// this is an extra cluster
	extra=kTRUE;
	layerExtra=layerId;
	indexExtra=ipt;
      }
    } // end loop on AliTrackPoints

    TString vtitle = spdv->GetTitle();
    if(!vtitle.Contains("3D")) continue; 

    // SPD tracklet
    if(indexAssociated[0]>=0 && indexAssociated[1]>=0) {
      AliTrackPoint pointSPD1,pointSPD2;
      array->GetPoint(pointSPD1,indexAssociated[0]);
      array->GetPoint(pointSPD2,indexAssociated[1]);
      Float_t phi=TMath::ATan2(pointSPD2.GetY()-pointSPD1.GetY(),pointSPD2.GetX()-pointSPD1.GetX());
      Float_t lambda=TMath::ATan((pointSPD2.GetZ()-pointSPD1.GetZ())/TMath::Sqrt((pointSPD2.GetX()-pointSPD1.GetX())*(pointSPD2.GetX()-pointSPD1.GetX())+(pointSPD2.GetY()-pointSPD1.GetY())*(pointSPD2.GetY()-pointSPD1.GetY())));
      Float_t theta=0.5*TMath::Pi()-lambda;
      TParticle particle(211,0,0,0,0,0,TMath::Cos(phi),TMath::Sin(phi),TMath::Tan(lambda),10.,pointSPD1.GetX(),pointSPD1.GetY(),pointSPD1.GetZ(),0);
      AliESDtrack tracklet(&particle);
      Float_t dz[2];
      // distance to primary SPD (only if 3D and high multiplicity)
      if(spdv->GetNContributors()>10) { 
	//tracklet.GetDZ(spdv->GetXv(),spdv->GetYv(),spdv->GetZv(),0,dz);
	tracklet.GetDZ(-0.07,0.25,spdv->GetZv(),0,dz);
	fNtupleITSAlignSPDTracklets->Fill(phi,theta,0.5*(pointSPD1.GetZ()+pointSPD2.GetZ()),dz[0],dz[1]);
      }
    }

    // distance to extra
    if(extra && spdv->GetNContributors()>4) {
      AliTrackPoint pointExtra,pointAssociated;
      array->GetPoint(pointAssociated,indexAssociated[layerExtra]);
      array->GetPoint(pointExtra,indexExtra);
      Float_t phiExtra = TMath::ATan2(pointExtra.GetY()-spdv->GetYv(),pointExtra.GetX()-spdv->GetXv());
      Float_t phiAssociated = TMath::ATan2(pointAssociated.GetY()-spdv->GetYv(),pointAssociated.GetX()-spdv->GetXv());
      Float_t rExtra = TMath::Sqrt((pointExtra.GetX()-spdv->GetXv())*(pointExtra.GetX()-spdv->GetXv())+(pointExtra.GetY()-spdv->GetYv())*(pointExtra.GetY()-spdv->GetYv()));
      Float_t rAssociated = TMath::Sqrt((pointAssociated.GetX()-spdv->GetXv())*(pointAssociated.GetX()-spdv->GetXv())+(pointAssociated.GetY()-spdv->GetYv())*(pointAssociated.GetY()-spdv->GetYv()));
      Float_t dzExtra[2];
      dzExtra[0] = (phiExtra-phiAssociated)*0.5*(rExtra+rAssociated);
      dzExtra[1] = pointExtra.GetZ()-pointAssociated.GetZ()-(rExtra-rAssociated)*(pointAssociated.GetZ()-spdv->GetZv())/rAssociated;
      Float_t xlocExtra=-100.,zlocExtra=-100.;
      fNtupleITSAlignExtra->Fill(layerExtra,pointExtra.GetX(),pointExtra.GetY(),pointExtra.GetZ(),xlocExtra,zlocExtra,dzExtra[0],dzExtra[1],nclsITS);  
    }
    

  } // end loop on tracks
  
  if(vertexMC) { delete vertexMC; vertexMC=0; }

  return;
}      

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }

  fHistNtracks = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNtracks"));
  fHistNclsITSMI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNclsITSMI"));
  fHistNclsITSSA = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNclsITSSA"));
  fHistClusterMapITSMI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMI"));
  fHistClusterMapITSMIbad = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMIbad"));
  fHistClusterMapITSMIskipped = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMIskipped"));
  fHistClusterMapITSMIoutinz = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMIoutinz"));
  fHistClusterMapITSMInorefit = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMInorefit"));
  fHistClusterMapITSMInocls = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSMInocls"));
  fHistClusterMapITSSA = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSA"));
  fHistClusterMapITSSAbad = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSAbad"));
  fHistClusterMapITSSAskipped = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSAskipped"));
  fHistClusterMapITSSAoutinz = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSAoutinz"));
  fHistClusterMapITSSAnorefit = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSAnorefit"));
  fHistClusterMapITSSAnocls = dynamic_cast<TH1F*>(fOutput->FindObject("fHistClusterMapITSSAnocls"));
  fNtupleESDTracks = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleESDTracks"));
  fNtupleITSAlignExtra = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleITSAlignExtra"));
  fNtupleITSAlignSPDTracklets = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleITSAlignSPDTracklets"));

  return;
}
//---------------------------------------------------------------------------
Int_t AliAnalysisTaskITSTrackingCheck::NumberOfITSClustersMC(Int_t label) const
{
  //
  // Return number of ITS clusters produced by MC particle with given label
  //
  
  AliESDInputHandlerRP *esdHRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdHRP) return -1;
  TTree *cTree = (TTree*)esdHRP->GetTreeR("ITS");
  if(!cTree) return -1;
  TClonesArray *clusters=0;   // new TClonesArray("AliITSRecPoint",10000);
  cTree->SetBranchAddress("ITSRecPoints",&clusters);
  if(!clusters) return -1;

  AliITSRecPoint *c=0;
  Int_t i,n,icl,lay,ilab;
  Int_t ncls[6]={0,0,0,0,0,0};
  Int_t nclstot=0;

  for(i=0; i<2198; i++) {
    cTree->GetEvent(i);
    n=clusters->GetEntriesFast();
    for (icl=0; icl<n; icl++) {
      c=(AliITSRecPoint*)clusters->UncheckedAt(icl);
      lay=c->GetLayer();
      for(ilab=0;ilab<3;ilab++) {
        if(c->GetLabel(ilab)==label) ncls[lay]++;
      }
    }
  }
  for(i=0;i<6;i++) { if(ncls[i]) nclstot++; }

  return nclstot;
    //return label*0;
}
//---------------------------------------------------------------------------
Double_t AliAnalysisTaskITSTrackingCheck::ParticleImpParMC(TParticle *part,
							   AliESDVertex *vert,
							   Double_t bzT) const
{
  //
  // Return the MC value of the impact parameter
  //
 
  Double_t vx=part->Vx()-vert->GetX();
  Double_t vy=part->Vy()-vert->GetY();
      
  Double_t pt=part->Pt();     
  Double_t px=part->Px();     
  Double_t py=part->Py();     
  Double_t charge = (part->GetPdgCode()>0. ? 1. : -1.);
  if(TMath::Abs(part->GetPdgCode())<100) charge*=-1.;

  if(px<0.000001) px=0.000001;     
  Double_t rAnd=((10./2.99792458)*pt/bzT)*100.;
  Double_t center[3],d0;
  center[0]=vx-(1./charge)*rAnd*(py/pt);
  center[1]=vy+(1./charge)*rAnd*(px/pt);
  center[2]=TMath::Sqrt(center[0]*center[0]+center[1]*center[1]);
  d0 = -center[2]+rAnd;

  return d0;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::SelectPt(Double_t pt)
{
  //
  // Keep only tracks in given pt bins
  // 
  Double_t ptlower[10]={0.29,0.49,0.75,0.9,1.9,3.5,6.5, 9.,19.,27.};
  Double_t ptupper[10]={0.31,0.51,0.85,1.1,2.1,4.5,7.5,11.,21.,33.};

  for(Int_t i=0; i<10; i++) {
    if(pt>ptlower[i] && pt<ptupper[i]) {
      fCountsPerPtBin[i]++;
      return kTRUE;
    }
  }
  return kFALSE;
  //return kTRUE;
}




