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


///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskHypertriton3 class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
// This task is optimized for ESDs.root
//
// Author: 
// S. Trogolo, trogolo@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <vector>

#include <TAxis.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHypertriton3.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliVertexerTracks.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#define dcaDP 0.2
#define dcaPIP 0.5
#define dcaDPI 0.5

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHypertriton3)
/* $Id$ */
//________________________________________________________________________
AliAnalysisTaskHypertriton3::AliAnalysisTaskHypertriton3():
AliAnalysisTaskSE("taskHypertriton"),
  fESDevent(0),
  fESDtrackCuts(0x0),
  fESDtrackCutsV0(0x0),
  fPrimaryVertex(0x0),
  fPIDResponse(0),
  fMC(kFALSE),
  fFillTree(kFALSE),
  fCentrality(),
  fCentralityPercentile(),
  fOutput(0),
  fHistCount(0),
  fHistCentralityClass(0),
  fHistCentralityPercentile(0),
  fHistTrigger(0),
  fHistMultiplicity(0),
  fHistZPrimaryVtx(0),
  fHistXPrimaryVtx(0),
  fHistYPrimaryVtx(0),
  fHistChi2perTPCcluster(0),
  fHistTPCpid(0),
  fHistTPCdeusignal(0),
  fHistTPCprosignal(0),
  fHistTPCpionsignal(0),
  fHistTPCantideusignal(0),
  fHistTPCantiprosignal(0),
  fHistTPCpionplussignal(0),
  fHistTOFsignal(0),
  fHistTOFdeusignal(0),
  fHistTOFprosignal(0),
  fHistTOFantideusignal(0),
  fHistTOFantiprosignal(0),
  fHistTOFdeumass(0),
  fHistTOFpromass(0),
  fHistpionTPCcls(0),
  fHistpTpion(0),
  fHistCorrDCAdprimary(0),
  fHistCorrDCApprimary(0),
  fHistCorrDCApiprimary(0),
  fHistDCApiprimary(0),
  fHistDCAdeupro(0),
  fHistDCApiondeu(0),	  
  fHistDCApionpro(0),
  fHistZDecayVtx(0),
  fHistXDecayVtx(0),
  fHistYDecayVtx(0),
  fHistDecayLengthH3L(0),
  fHistDCAXYdeuvtx(0),
  fHistDCAZdeuvtx(0),
  fHistDCAXYprovtx(0),
  fHistDCAZprovtx(0),
  fHistDCAXYpionvtx(0),
  fHistDCAZpionvtx(0),
  fHistCosPointingAngle(0),
  fHistMassHypertriton(0),
  fHistParticle(0x0),
  fHistpionTPCclsMCt(0x0),
  fHistpTpionMCt(0x0),
  fHistpTproMCt(0x0),
  fHistpTdeuMCt(0x0),
  fHistDCAdeuproMCt(0x0),
  fHistDCApiondeuMCt(0x0),
  fHistDCApionproMCt(0x0),
  fHistCorrDCAdprimaryMCt(0x0),
  fHistCorrDCApprimaryMCt(0x0),
  fHistCorrDCApiprimaryMCt(0x0),
  fHistDCApiprimaryMCt(0x0),
  fHistDCApprimaryMCt(0x0),
  fHistDCAdprimaryMCt(0x0),
  fHistDCAXYdeuvtxMCt(0x0),
  fHistDCAZdeuvtxMCt(0x0),
  fHistDCAXYprovtxMCt(0x0),
  fHistDCAZprovtxMCt(0x0),
  fHistDCAXYpionvtxMCt(0x0),
  fHistDCAZpionvtxMCt(0x0),
  fHistZDecayVtxMCt(0x0),
  fHistXDecayVtxMCt(0x0),
  fHistYDecayVtxMCt(0x0),
  fHistMassHypertritonMCt(0x0),
  fTTree(0x0),
  fTchi2NDFdeu(0),
  fTPCclsdeu(0),
  fTPCclsPIDdeu(0),
  fTpTPCdeu(0),
  fTpTdeu(0),
  fTpdeu(0),
  fTTPCnsigmadeu(0),
  fTTOFmassdeu(0),
  fTDCAXYdeuprvtx(0),
  fTDCAZdeuprvtx(0),
  fTchi2NDFpro(0),
  fTPCclspro(0),
  fTPCclsPIDpro(0),
  fTpTPCpro(0),
  fTpTpro(0),
  fTppro(0),
  fTTPCnsigmapro(0),
  fTTOFmasspro(0),
  fTDCAXYproprvtx(0),
  fTDCAZproprvtx(0),
  fTchi2NDFpion(0),
  fTPCclspion(0),
  fTPCclsPIDpion(0),
  fTpTPCpion(0),
  fTpTpion(0),
  fTppion(0),
  fTTPCnsigmapion(0),
  fTDCAXYpioprvtx(0),
  fTDCAZpioprvtx(0),
  fTDCAdp(0),
  fTDCAdpi(0),
  fTDCAppi(0),
  fTDCAXYdvtx(0),
  fTDCAZdvtx(0),
  fTDCAXYpvtx(0),
  fTDCAZpvtx(0),
  fTDCAXYpivtx(0),
  fTDCAZpivtx(0),
  fTDecayLength(0),
  fTCosPA(0),
  fTInvariantMass(0)
{
  //Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  
}


//________________________________________________________________________
AliAnalysisTaskHypertriton3::~AliAnalysisTaskHypertriton3(){
  //Destructor
  if(fOutput){
    delete fOutput;
    fOutput = 0;
  }

  if(fTTree) delete fTTree;

    if(fPIDResponse){
    delete fPIDResponse;
  }

    if(fESDtrackCuts) delete fESDtrackCuts;
    if(fESDtrackCutsV0) delete fESDtrackCutsV0;


} // end of Destructor

//________________________________________________________________________
Double_t AliAnalysisTaskHypertriton3::GetDCAcut(Int_t part, Double_t dca)const{
  Double_t cut;
  if(part == 5){ // p-d vs pi-d
    cut = dcaDPI*(1-(dca/dcaDP));
    return cut;
  }
  else if(part == 4){ // p-d vs p-pi 
    cut = dcaPIP*(1-(dca/dcaDP));
    return cut;
  }
  return -1;
}


//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserCreateOutputObjects(){
  // Create a TList with Histograms
  // Called once
  printf("**************************************************\n");
  printf("AliAnalysisTaskHypertriton3::CreateOutputObjects()\n");
  printf("**************************************************\n");

  //ESD Track cuts
  if(!fESDtrackCuts) fESDtrackCuts = new AliESDtrackCuts();
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);

  //ESD Track cuts V0
  if(!fESDtrackCutsV0) fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsV0->SetMinNClustersTPC(100);
  fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
  fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
  fESDtrackCutsV0->SetPtRange(0.2,1.2);

  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("clistHypertriton");

  fHistCount = new TH1F("fHistCount","Counter histogram",2,-0.5,1.5);
  fHistCount->GetXaxis()->SetBinLabel(1,"Event");
  
  fHistCentralityClass = new TH1F("fHistCentralityClass","Centrality Class; centrality class; entries",11,-0.5,10.5);
  fHistCentralityPercentile = new TH1F("fHistCentralityPercentile","Centrality; centrality percentile; entries",101,-0.5,100.5);
  
  fHistTrigger = new TH1F("fHistTrigger","Trigger statistics",3,-0.5,3.5);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"kMB");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");

  fHistMultiplicity = new TH1F("fHistMultiplicity ", "Multiplicity; multiplicity; entries", 100, 0, 20000);
  
  fHistZPrimaryVtx = new TH1F("fHistZPrimaryVtx","primary vertex - Z coordinate; z_{primary vtx} (cm); entries",40,-10.,10.);

  fHistXPrimaryVtx = new TH1F("fHistXPrimaryVtx","primary vertex - X coordinate; x_{primary vtx} (cm); entries",2000,-1.,1.);
  
  fHistYPrimaryVtx = new TH1F("fHistYPrimaryVtx","primary vertex - Y coordinate; y_{primary vtx} (cm); entries",2000,-1.,1.);

  fHistChi2perTPCcluster = new TH1F("fHistChi2perTPCcluster","#Chi^{2}/TPC cluster distribution; #Chi^{2}/TPCcls; entries",44,-0.5,10.5);
 
  
  //TPC
  fHistTPCpid = new TH2F("fHistTPCpid", "dE/dx for all tracks; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCpid->SetOption("scat");
  fHistTPCpid->SetMarkerStyle(kFullCircle);

  //Hypertriton prongs
  fHistTPCdeusignal = new TH2F("fHistTPCdeusignal", "dE/dx after d PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCdeusignal->SetOption("scat");
  fHistTPCdeusignal->SetMarkerStyle(kFullCircle);

  fHistTPCprosignal = new TH2F("fHistTPCprosignal", "dE/dx after p PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCprosignal->SetOption("scat");
  fHistTPCprosignal->SetMarkerStyle(kFullCircle);

  fHistTPCpionsignal= new TH2F("fHistTPCpionsignal", "dE/dx after  #pi^{-} PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCpionsignal->SetOption("scat");
  fHistTPCpionsignal->SetMarkerStyle(kFullCircle);
  
  //Anti-hypertriton prongs
  fHistTPCantideusignal = new TH2F("fHistTPCantideusignal", "dE/dx after d PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCantideusignal->SetOption("scat");
  fHistTPCantideusignal->SetMarkerStyle(kFullCircle);

  fHistTPCantiprosignal = new TH2F("fHistTPCantiprosignal", "dE/dx after p PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCantiprosignal->SetOption("scat");
  fHistTPCantiprosignal->SetMarkerStyle(kFullCircle);

  fHistTPCpionplussignal= new TH2F("fHistTPCpionplussignal", "dE/dx after  #pi^{+} PID; p (GeV/c); TPC signal", 1000, 0., 10, 300, 0, 2100);
  fHistTPCpionplussignal->SetOption("scat");
  fHistTPCpionplussignal->SetMarkerStyle(kFullCircle);

  //TOF
  
  fHistTOFsignal = new TH2F("fHistTOFsignal","TOF signal; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFdeusignal = new TH2F("fHistTOFdeusignal","#beta vs TPCmom - deuteron; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFprosignal = new TH2F("fHistTOFprosignal","#beta vs TPCmom - proton; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFantideusignal = new TH2F("fHistTOFantideusignal","#beta vs TPCmom - anti-deuteron; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFantiprosignal = new TH2F("fHistTOFantiprosignal","#beta vs TPCmom - anti-proton; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);
  
  fHistTOFdeumass = new TH1F("fHistTOFdeumass","deuteron mass distribution - TOF; mass (GeV/c^{2}); entries",400,0.8,2.8);
  
  fHistTOFpromass = new TH1F("fHistTOFpromass","proton mass distribution - TOF; mass (GeV/c^{2}); entries",200,0.5,1.5);


  fHistpionTPCcls = new TH1F("fHistpionTPCcls","#pi^{-} TPC clusters; TPC clusters; entries",201,-0.5,200.5);
  fHistpTpion = new TH1F("fHistpTpion","pion p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
  fHistCorrDCAdprimary = new TH2F("fHistCorrDCAdprimary","DCA_{PV,xy} vs DCA_{PV,z} - deuteron; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
  fHistCorrDCApprimary = new TH2F("fHistCorrDCApprimary","DCA_{PV,xy} vs DCA_{PV,z} - proton; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
  fHistCorrDCApiprimary = new TH2F("fHistCorrDCApiprimary","DCA_{PV,xy} vs DCA_{PV,z} - pion; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
  fHistDCApiprimary = new TH1F("fHistDCApiprimary","DCA pion-primary vertex; DCA (cm); entries",3200,0.f,80.f);
  
  //DCA prongs
  fHistDCAdeupro = new TH1F("fHistDCAdeupro","DCA d-p tracks;d-p DCA (cm);entries",800,-1.0,5.0);

  fHistDCApiondeu = new TH1F("fHistDCApiondeu","DCA #pi^{-}-d tracks; #pi^{-}-d DCA (cm); entries",800,-1.0,5.0);
  
  fHistDCApionpro = new TH1F("fHistDCApionpro","DCA #pi^{-}-p tracks; #pi^{-}-p DCA (cm); entries",800,-1.0,5.0);

  //Decay histo
  fHistZDecayVtx = new TH1F("fHistZDecayVtx","decay vertex - z coordinate; z_{decay vtx} (cm); entries",800,-20.f,20.f);
  fHistXDecayVtx = new TH1F("fHistXDecayVtx","decay vertex - x coordinate; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistYDecayVtx = new TH1F("fHistYDecayVtx","decay vertex - y coordinate; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistDecayLengthH3L = new TH1F("fHistDecayLengthH3L","decay length ^{3}H_{#Lambda}; decay length (cm); entries",400,0.,100.);

  fHistDCAXYdeuvtx = new TH1F("fHistDCAXYdeuvtx","DCA candidate d-decay vertex - xy coordinate; DCA_{xy} (cm); entries",200,-5.,5.);
  fHistDCAZdeuvtx = new TH1F("fHistDCAZdeuvtx","DCA candidate d-decay vertex - z coordinate; DCA_{z} (cm); entries",200,-10.,10.);
  fHistDCAXYprovtx = new TH1F("fHistDCAXYprovtx","DCA candidate p-decay vertex - xy coordinate; DCA_{xy} (cm); entries",200,-5.,5.);
  fHistDCAZprovtx = new TH1F("fHistDCAZprovtx","DCA candidate p-decay vertex - z coordinate; DCA_{z} (cm); entries",200,-10.,10.);
  fHistDCAXYpionvtx = new TH1F("fHistDCAXYpionvtx","DCA candidate #pi^{-}-decay vertex - xy coordinate; DCA_{xy} (cm); entries",200,-5.,5.);
  fHistDCAZpionvtx = new TH1F("fHistDCAZpionvtx","DCA candidate #pi^{-}-decay vertex - z coordinate; DCA_{z} (cm); entries",200,-10.,10.);
  fHistCosPointingAngle= new TH1F("fHistCosPointingAngle", "Cos pointing angle distribution; cos point angle; entries", 220, -1.1, 1.1);
  fHistMassHypertriton = new TH1F("fHistMassHypertriton", "Invariant mass distribution d+p+#pi^{-};invariant mass d+p+#pi^{-} (GeV/c^{2}); entries ", 400, 2.9, 3.1);
  if(fMC){
    fHistParticle = new TH1F("fHistParticle","Check particle candidate",23,-0.5,22.5);
    fHistParticle->GetXaxis()->SetBinLabel(1,"electron");
    fHistParticle->GetXaxis()->SetBinLabel(2,"positron");
    fHistParticle->GetXaxis()->SetBinLabel(3,"#pi^{+}");
    fHistParticle->GetXaxis()->SetBinLabel(4,"#pi^{-}");
    fHistParticle->GetXaxis()->SetBinLabel(5,"K^{+}");
    fHistParticle->GetXaxis()->SetBinLabel(6,"K^{-}");
    fHistParticle->GetXaxis()->SetBinLabel(7,"proton");
    fHistParticle->GetXaxis()->SetBinLabel(8,"anti-proton");
    fHistParticle->GetXaxis()->SetBinLabel(9,"deuteron");
    fHistParticle->GetXaxis()->SetBinLabel(10,"anti-deuteron");
    fHistParticle->GetXaxis()->SetBinLabel(11,"triton");
    fHistParticle->GetXaxis()->SetBinLabel(12,"anti-triton");
    fHistParticle->GetXaxis()->SetBinLabel(13,"He3");
    fHistParticle->GetXaxis()->SetBinLabel(14,"anti-He3");
    fHistParticle->GetXaxis()->SetBinLabel(15,"#alpha");
    fHistParticle->GetXaxis()->SetBinLabel(16,"anti-#alpha");
    fHistParticle->GetXaxis()->SetBinLabel(17,"^{3}H_{#Lambda}");
    fHistParticle->GetXaxis()->SetBinLabel(18,"anti-^{3}H_{#Lambda}");
    fHistParticle->GetXaxis()->SetBinLabel(19,"^{3}H_{#Lambda} - 2 body");
    fHistParticle->GetXaxis()->SetBinLabel(20,"anti-^{3}H_{#Lambda} - 2 body");
    fHistParticle->GetXaxis()->SetBinLabel(21,"#Lambda - 2 body");
    
    fHistpionTPCclsMCt = new TH1F("fHistpionTPCclsMCt","#pi^{-} TPC clusters - MCtruth; TPC clusters; entries",201,-0.5,200.5);
    fHistpTpionMCt = new TH1F("fHistpTpionMCt","pion p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistpTproMCt = new TH1F("fHistpTproMCt","proton p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistpTdeuMCt = new TH1F("fHistpTdeuMCt","deuteron p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistDCAdeuproMCt = new TH1F("fHistDCAdeuproMCt","DCA d-p tracks - MCtruth;d-p DCA (cm);entries",250,-1.0,1.0);
    fHistDCApiondeuMCt = new TH1F("fHistDCApiondeuMCt","DCA #pi^{-}-d tracks - MCtruth;#pi^{-}-d DCA (cm);entries",250,-1.0,1.0);
    fHistDCApionproMCt = new TH1F("fHistDCApionproMCt","DCA #pi^{-}-p tracks - MCtruth;#pi^{-}-p DCA (cm);entries",250,-1.0,1.0);
    
    fHistCorrDCAdprimaryMCt = new TH2F("fHistCorrDCAdprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - deuteron MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
    fHistCorrDCApprimaryMCt = new TH2F("fHistCorrDCApprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - proton MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
    fHistCorrDCApiprimaryMCt = new TH2F("fHistCorrDCApiprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - pion MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",3200,-40.f,40.f,3200,-40.f,40.f);
    fHistDCApiprimaryMCt = new TH1F("fHistDCApiprimaryMCt","DCA pion-primary vertex; DCA (cm); entries",3200,0.f,80.f);
    fHistDCApprimaryMCt = new TH1F("fHistDCApprimaryMCt","DCA proton-primary vertex; DCA (cm); entries",3200,0.f,80.f);
    fHistDCAdprimaryMCt = new TH1F("fHistDCAdprimaryMCt","DCA deuteron-primary vertex; DCA (cm); entries",3200,0.f,80.f);
    
    fHistDCAXYdeuvtxMCt = new TH1F("fHistDCAXYdeuvtxMCt","DCA candidate d-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZdeuvtxMCt = new TH1F("fHistDCAZdeuvtxMCt","DCA candidate d-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);
    fHistDCAXYprovtxMCt = new TH1F("fHistDCAXYprovtxMCt","DCA candidate p-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZprovtxMCt = new TH1F("fHistDCAZprovtxMCt","DCA candidate p-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);
    fHistDCAXYpionvtxMCt = new TH1F("fHistDCAXYpionvtxMCt","DCA candidate #pi^{-}-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZpionvtxMCt = new TH1F("fHistDCAZpionvtxMCt","DCA candidate #pi^{-}-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);
    
    fHistZDecayVtxMCt = new TH1F("fHistZDecayVtxMCt","decay vertex - z coordinate - MCtruth; z_{decay vtx} (cm); entries",800,-20.f,20.f);
    fHistXDecayVtxMCt = new TH1F("fHistXDecayVtxMCt","decay vertex - x coordinate - MCtruth; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
    fHistYDecayVtxMCt = new TH1F("fHistYDecayVtxMCt","decay vertex - y coordinate - MCtruth; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
    
    fHistMassHypertritonMCt = new TH1F("fHistMassHypertritonMCt","^{3}H_{#Lambda} invariant mass - MCtruth; invariant mass d+p+#pi^{-} (GeV/c^{2}); entries", 400,2.9,3.1);
  }
  
  fOutput->Add(fHistCount);
  fOutput->Add(fHistCentralityClass);
  fOutput->Add(fHistCentralityPercentile);
  fOutput->Add(fHistTrigger);
  fOutput->Add(fHistMultiplicity);
  fOutput->Add(fHistZPrimaryVtx);
  fOutput->Add(fHistXPrimaryVtx);
  fOutput->Add(fHistYPrimaryVtx);
  fOutput->Add(fHistChi2perTPCcluster);
  fOutput->Add(fHistTPCpid);
  fOutput->Add(fHistTPCdeusignal);
  fOutput->Add(fHistTPCprosignal);
  fOutput->Add(fHistTPCpionsignal);
  fOutput->Add(fHistTPCantideusignal);
  fOutput->Add(fHistTPCantiprosignal);
  fOutput->Add(fHistTPCpionplussignal);
  fOutput->Add(fHistTOFsignal);
  fOutput->Add(fHistTOFdeusignal);
  fOutput->Add(fHistTOFprosignal);
  fOutput->Add(fHistTOFantideusignal);
  fOutput->Add(fHistTOFantiprosignal);
  fOutput->Add(fHistTOFdeumass);
  fOutput->Add(fHistTOFpromass);
  fOutput->Add(fHistpionTPCcls);
  fOutput->Add(fHistpTpion);
  fOutput->Add(fHistCorrDCAdprimary);
  fOutput->Add(fHistCorrDCApprimary);
  fOutput->Add(fHistCorrDCApiprimary);
  fOutput->Add(fHistDCApiprimary);
  fOutput->Add(fHistDCAdeupro);
  fOutput->Add(fHistDCApiondeu);	  
  fOutput->Add(fHistDCApionpro);
  fOutput->Add(fHistZDecayVtx);
  fOutput->Add(fHistXDecayVtx);
  fOutput->Add(fHistYDecayVtx);
  fOutput->Add(fHistDecayLengthH3L);
  fOutput->Add(fHistDCAXYdeuvtx);
  fOutput->Add(fHistDCAZdeuvtx);
  fOutput->Add(fHistDCAXYprovtx);
  fOutput->Add(fHistDCAZprovtx);
  fOutput->Add(fHistDCAXYpionvtx);
  fOutput->Add(fHistDCAZpionvtx);
  fOutput->Add(fHistCosPointingAngle);
  fOutput->Add(fHistMassHypertriton);

  if(fMC){
    fOutput->Add(fHistParticle);
    fOutput->Add(fHistpionTPCclsMCt);
    fOutput->Add(fHistpTpionMCt);
    fOutput->Add(fHistpTproMCt);
    fOutput->Add(fHistpTdeuMCt);
    fOutput->Add(fHistDCAdeuproMCt);
    fOutput->Add(fHistDCApiondeuMCt);
    fOutput->Add(fHistDCApionproMCt);
    fOutput->Add(fHistCorrDCAdprimaryMCt);
    fOutput->Add(fHistCorrDCApprimaryMCt);
    fOutput->Add(fHistCorrDCApiprimaryMCt);
    fOutput->Add(fHistDCApiprimaryMCt);
    fOutput->Add(fHistDCApprimaryMCt);
    fOutput->Add(fHistDCAdprimaryMCt);
    fOutput->Add(fHistDCAXYdeuvtxMCt);
    fOutput->Add(fHistDCAZdeuvtxMCt);
    fOutput->Add(fHistDCAXYprovtxMCt);
    fOutput->Add(fHistDCAZprovtxMCt);
    fOutput->Add(fHistDCAXYpionvtxMCt);
    fOutput->Add(fHistDCAZpionvtxMCt);
    fOutput->Add(fHistZDecayVtxMCt);
    fOutput->Add(fHistXDecayVtxMCt);
    fOutput->Add(fHistYDecayVtxMCt);
    fOutput->Add(fHistMassHypertritonMCt);
  }
  
  // Post output data.
  PostData(1,fOutput);
  printf("**************************************************\n");
  printf("end of fOutput\n");
  printf("**************************************************\n");

  if(fFillTree){
  OpenFile(2);
  fTTree = new TTree("hypertriton","hypertriton candidates");
  fTTree->Branch("Chi2NDFdeu",&fTchi2NDFdeu,"Chi2NDFdeu/F");
  fTTree->Branch("TPCclsdeu",&fTPCclsdeu,"TPCclsdeu/s");
  fTTree->Branch("TPCclsPIDdeu",&fTPCclsPIDdeu,"TPCclsPIDdeu/s");
  fTTree->Branch("pTPCdeu",&fTpTPCdeu,"pTPCdeu/F");
  fTTree->Branch("pTdeu",&fTpTdeu,"pTdeu/F");
  fTTree->Branch("pdeu",&fTpdeu,"pdeu/F");
  fTTree->Branch("TPCnsigmadeu",&fTTPCnsigmadeu,"TPCnsigmadeu/F");
  fTTree->Branch("TOFmassdeu",&fTTOFmassdeu,"TOFmassdeu/F");
  fTTree->Branch("DCAxydeuprim",&fTDCAXYdeuprvtx,"DCAxydeuprim/F");
  fTTree->Branch("DCAzdeuprim",&fTDCAZdeuprvtx,"DCAzdeuprim/F");
  fTTree->Branch("Chi2NDFpro",&fTchi2NDFpro,"Chi2NDFpro/F");
  fTTree->Branch("TPCclspro",&fTPCclspro,"TPCclspro/s");
  fTTree->Branch("TPCclsPIDpro",&fTPCclsPIDpro,"TPCclsPIDpro/s");
  fTTree->Branch("pTPCpro",&fTpTPCpro,"pTPCpro/F");
  fTTree->Branch("pTpro",&fTpTpro,"pTpro/F");
  fTTree->Branch("ppro",&fTppro,"ppro/F");
  fTTree->Branch("TPCnsigmapro",&fTTPCnsigmapro,"TPCnsigmapro/F");
  fTTree->Branch("TOFmasspro",&fTTOFmasspro,"TOFmasspro/F");
  fTTree->Branch("DCAxyproprim",&fTDCAXYproprvtx,"DCAxyproprim/F");
  fTTree->Branch("DCAzproprim",&fTDCAZproprvtx,"DCAzproprim/F");
  fTTree->Branch("Chi2NDFpion",&fTchi2NDFpion,"Chi2NDFpion/F");
  fTTree->Branch("TPCclspion",&fTPCclspion,"TPCclspion/s");
  fTTree->Branch("TPCclsPIDpion",&fTPCclsPIDpion,"TPCclsPIDpion/s");
  fTTree->Branch("pTPCpion",&fTpTPCpion,"pTPCpion/F");
  fTTree->Branch("pTpion",&fTpTpion,"pTpion/F");
  fTTree->Branch("ppion",&fTppion,"ppion/F");
  fTTree->Branch("TPCnsigmapion",&fTTPCnsigmapion,"TPCnsigmapion/F");
  fTTree->Branch("DCAxypioprim",&fTDCAXYpioprvtx,"DCAxypioprim/F");
  fTTree->Branch("DCAzpioprim",&fTDCAZpioprvtx,"DCAzpioprim/F");
  fTTree->Branch("DCAdp",&fTDCAdp,"DCAdp/F");
  fTTree->Branch("DCAdpi",&fTDCAdpi,"DCAdpi/F");
  fTTree->Branch("DCAppi",&fTDCAppi,"DCAppi/F");
  fTTree->Branch("DCAXYdeuvtx",&fTDCAXYdvtx,"DCAXYdeuvtx/F");
  fTTree->Branch("DCAZdeuvtx",&fTDCAZdvtx,"DCAZdeuvtx/F");
  fTTree->Branch("DCAXYprovtx",&fTDCAXYpvtx,"DCAXYprovtx/F");
  fTTree->Branch("DCAZprovtx",&fTDCAZpvtx,"DCAZprovtx/F");
  fTTree->Branch("DCAXYpionvtx",&fTDCAXYpivtx,"DCAXYpionvtx/F");
  fTTree->Branch("DCAZpionvtx",&fTDCAZpivtx,"DCAZpionvtx/F");
  fTTree->Branch("DecayLength",&fTDecayLength,"DecayLength/F");
  fTTree->Branch("CosPA",&fTCosPA,"CosPA/F");
  fTTree->Branch("InvariantMass",&fTInvariantMass,"InvariantMass/F");

  fTTree->SetAutoSave(150000000);
  PostData(2,fTTree);
  } else {
    fTTree = new TTree();
  } 

  
  printf("**************************************************\n");
  printf("end of CreateOutputObjects\n");
  printf("**************************************************\n");
  
} // end of UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserExec(Option_t *){
  // Main loop
  // Called for each event
  //----------------------------------------------------------------------------
  // Mass definition
  Double_t pionMass     = TDatabasePDG::Instance()->GetParticle(211)->Mass(); // pion mass = 0.13957 GeV/c2
  Double_t protonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); // proton mass = 0.93827 GeV/c2
  Double_t deuteronMass = AliPID::ParticleMass(AliPID::kDeuteron); // deuteron mass = 1.87561 GeV/c2
  //Double_t helium3Mass = AliPID::ParticleMass(AliPID::kHe3); // helium3 mass = 2.80923 GeV/c2
  
  //PDGCodes
  //particle
  //  Long_t pdgPionPlus = 211;
  Long_t pdgProton = 2212;
  Long_t pdgDeuteron = 1000010020;
  //Long_t pdgHypertriton = 1010010030;
  
  //antiparticle
  Long_t pdgPionMinus = -211;
  //Long_t pdgAntiProton = -2212;
  //Long_t pdgAntiDeuteron = -1000010020;
  //Long_t pdgAntiHypertriton = -1010010030;
  //----------------------------------------------------------------------------


  //==========Define variables==========//
  //===PID loop===//
  Int_t ntracks,label, labelM = 0 ;
  Double_t chi2PerClusterTPC, nClustersTPC=0.;
  Double_t p = 0.;
  Float_t mass, beta, time, time0, gamma = 0.;
  Float_t minPos = 9999.;
  Float_t nsigma[3] = {9999.,9999.,9999.};
  Int_t partID = -1;
  AliESDtrack *track = 0x0;


  //===Combining tracks loop===//
  Int_t ndau = 0;
  Float_t dcapiprim, dcapprim, dcadprim = 0.;
  Float_t dprim[2] = {0.,0.}; //array for the dca from primary vertex coordinates (xy,z) - deuteron
  Float_t dprimc[3] = {0.,0.,0.}; // covariance matrix elements of the TPC only impact parameter
  Float_t pprim[2] = {0.,0.}; //array for the dca from primary vertex coordinates (xy,z) - proton
  Float_t pprimc[3] = {0.,0.,0.}; // covariance matrix elements of the TPC only impact parameter
  Float_t piprim[2] = {0.,0.}; //array for the dca from primary vertex coordinates (xy,z) - pion
  Float_t piprimc[3] = {0.,0.,0.}; // covariance matrix elements of the TPC only impact parameter
  Double_t dlh[3] = {0,0,0}; //array for the coordinates of the decay length
  Double_t dca_dp, dca_dpi, dca_ppi = 0.;
  Double_t dcad[2] = {0.,0.}; // dca between the candidate d,p,pi
  Double_t dcap[2] = {0.,0.}; // and the candidate decay vertex
  Double_t dcapi[2] = {0.,0.}; // dcad[0]= transverse plane coordinate; dcad[1]= z coordinate
  Double_t decayLengthH3L, pTmom, rapidity, pointingAngleH= 0.;
  Double_t lD, lP, lPi = 0;
  Double_t decVt[3] = {0.,0.,0.};
  Bool_t brotherHood = kFALSE;
  TLorentzVector posD, posP, negPi;
  AliESDtrack *trackD = 0x0;
  AliESDtrack *trackP = 0x0;
  AliESDtrack *trackNPi = 0x0;
  
  
  
  //==========ESD event==========
  fESDevent=(AliESDEvent*)InputEvent();
  if(!fESDevent){
    printf("AliAnalysisTaskHypertriton3::Exec(): bad ESD\n");
    return;
  }
  fHistCount->Fill(0);

 //==========MC info==========
  AliStack *stack=0x0;
  AliMCEvent *mcEvent = 0x0;
  if(fMC){//MC info and sample selection
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      printf("ERROR: stack not available\n");
      return;
    }
  } // end of MC info

  
  
  //==========Trigger class==========
  UInt_t maskPhysSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelectedCentral = (maskPhysSel & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (maskPhysSel & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB = (maskPhysSel & AliVEvent::kMB);

  if(isSelectedMB) { // Minimum Bias
    fHistTrigger->Fill(0);
  }

  if(isSelectedCentral) { // Central
    fHistTrigger->Fill(1);
  }
  
  if(isSelectedSemiCentral) { // SemiCentral
    fHistTrigger->Fill(2);
  }

  if(!isSelectedMB && !isSelectedCentral && !isSelectedSemiCentral) return;

  //==========Centrality==========
  if(fESDevent->GetEventSpecie() == 4){ // Event Specie == 4 == PbPb
    AliCentrality *centr=fESDevent->GetCentrality();
    fCentrality = centr->GetCentralityClass10("V0M");
    fCentralityPercentile = centr->GetCentralityPercentile("V0M");
  }
  
  fHistCentralityClass->Fill(fCentrality);
  fHistCentralityPercentile->Fill(fCentralityPercentile);

  if (fCentrality < 0. || fCentrality > 8.) return; //0 bis 80 %
  
  //==========Multiplicity==========
  Int_t refMultTpc = AliESDtrackCuts::GetReferenceMultiplicity(fESDevent, kTRUE);
  fHistMultiplicity->Fill(refMultTpc);
  

  //==========Primary Vertex==========
  fPrimaryVertex = (AliESDVertex*)fESDevent->GetPrimaryVertex();
  if (fPrimaryVertex->GetNContributors()<1) {
      // SPD vertex
    fPrimaryVertex = (AliESDVertex*)fESDevent->GetPrimaryVertexSPD();
      if(fPrimaryVertex->GetNContributors()<1) return;  
  }
  
  if (TMath::Abs(fPrimaryVertex->GetZ()) > 10) return;

  fHistZPrimaryVtx->Fill(fPrimaryVertex->GetZ());
  fHistXPrimaryVtx->Fill(fPrimaryVertex->GetX());
  fHistYPrimaryVtx->Fill(fPrimaryVertex->GetY());

  
  // -------------------------------------------------------
  // Loop for PID on ESD tracks
  // -------------------------------------------------------

  fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  
  vector <AliESDtrack*> cdeuteron;
  vector <AliESDtrack*> cproton;
  vector <AliESDtrack*> cpion;
  
  vector <Float_t> cmassd;
  vector <Float_t> cmassp;
    
  ntracks = fESDevent->GetNumberOfTracks();

  for(Int_t i=0; i < ntracks; i++) {
    track = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(i));
    minPos = 9999.;
    partID = -1;
    nsigma[0] = 9999.;
    nsigma[1] = 9999.;
    nsigma[2] = 9999.;
    
    // Chi2/TPCcls
    nClustersTPC = track->GetTPCclusters(0);
    chi2PerClusterTPC = track->GetTPCchi2()/nClustersTPC;
    fHistChi2perTPCcluster->Fill(chi2PerClusterTPC);

    if(fMC){if(track->GetLabel()<0) continue;}
    if(track->GetID()<0) continue;     

    if(!fESDtrackCuts->AcceptTrack(track)) continue;
    if(!track->GetInnerParam()) continue;
        
    if(track->GetTPCsignalN()<80) continue;
    
    p = track->GetInnerParam()->GetP(); //track->GetTPCmomentum()
    //if(p<0.2) continue;
    
    fHistTPCpid->Fill(p, track->GetTPCsignal());

    // TOF signal
    mass = 0.;
    beta = 0.;  
 
      if(track->GetIntegratedLength() > 350.){
	time0 = fPIDResponse->GetTOFResponse().GetStartTime(track->P());
	time = track->GetTOFsignal() - time0;

	if(time > 0){
	  beta = (track->GetIntegratedLength()) / (2.99792457999999984e-02 * time);
	  gamma = 1/TMath::Sqrt(1-beta*beta);
	  mass = p/TMath::Sqrt(gamma*gamma - 1);

	  fHistTOFsignal->Fill(p,beta);
	} // time > 0	
      } // track->GetIntergratedLength > 350
      
    //Filling PID histo
    label = track->GetLabel();
    
    if(fMC){
      TParticle *tparticleDaughter = stack->Particle(TMath::Abs(label));
      if(track->GetSign()>0){ // track sign +
	if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kDeuteron)) <= 4 && tparticleDaughter->GetPdgCode() == 1000010020){ //d
	  fHistTPCdeusignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.){
	    fHistTOFdeusignal->Fill(p,beta);
	    fHistTOFdeumass->Fill(mass);
	  }
	  cdeuteron.push_back(track);
	  cmassd.push_back(mass); // da spostare nel if length > 350.?
	}
	else if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton)) <= 4 && tparticleDaughter->GetPdgCode() == 2212) { // p
	  fHistTPCprosignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.){
	    fHistTOFprosignal->Fill(p,beta);
	    fHistTOFpromass->Fill(mass);
	  }
	  cproton.push_back(track);
	  cmassp.push_back(mass); // da spostare nel if length > 350.?
	}
	else if (tparticleDaughter->GetPdgCode() == 211) fHistTPCpionplussignal->Fill(p, track->GetTPCsignal());
      } // end of "track sign +"
      else if(track->GetSign()<0){ // track sign -
	if (tparticleDaughter->GetPdgCode() == -1000010020) { //anti-d
	  fHistTPCantideusignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.) fHistTOFantideusignal->Fill(p,beta);
	}
	else if (tparticleDaughter->GetPdgCode() == -2212) { //anti-p
	  fHistTPCantiprosignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.) fHistTOFantiprosignal->Fill(p,beta);
	}
	else if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)) <= 4 && tparticleDaughter->GetPdgCode() == -211) { // pi-
	  fHistTPCpionsignal->Fill(p, track->GetTPCsignal());
	  cpion.push_back(track);
	}
      } // end of "track sign -"
    } // end of MC PID
    
    else {
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)) <= 3) nsigma[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton)) <= 3) nsigma[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kDeuteron)) <= 3) nsigma[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kDeuteron));

      for(Int_t iS = 0; iS < 3; iS++){
	if(nsigma[iS] < minPos){
	  minPos = nsigma[iS];
	  partID = iS;
	}
      }
      
      
      if(track->GetSign()>0){
	if(partID == 2) { //deuteron
	  fHistTPCdeusignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.){
	    fHistTOFdeusignal->Fill(p,beta);
	    fHistTOFdeumass->Fill(mass);
	    cmassd.push_back(mass);
	  }
	  cdeuteron.push_back(track);
	}
	if(partID == 1) { // proton
	  fHistTPCprosignal->Fill(p, track->GetTPCsignal());
	  if(track->GetIntegratedLength() > 350.){
	    fHistTOFprosignal->Fill(p,beta);
	    fHistTOFpromass->Fill(mass);
	    cmassp.push_back(mass);
	  }
	  cproton.push_back(track);
	}
      } // end of "track sign +"
      else if(track->GetSign()<0){
	if(partID == 0) { //pion^-
	  fHistTPCpionsignal->Fill(p, track->GetTPCsignal());
	  cpion.push_back(track);
	}
      }// end of "track sign -"
    } 
  } // end of PID loop
  


  
  // -------------------------------------------------------
  // Loop for Invariant Mass
  // -------------------------------------------------------

  Double_t bz = fESDevent->GetMagneticField();
  Double_t xthiss(0.0);
  Double_t xpp(0.0);

  AliESDVertex *esdV1 = new AliESDVertex(*fPrimaryVertex);
  AliVertexerTracks *vertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
  TObjArray *trkarray = new TObjArray(3);
  AliESDVertex *decayVtx = 0x0;


  //========Set of Starting Cut========//

  Double_t fgDCAPiPVmin = 0.1;
  Double_t fgCosPointingAngle = 0.998;
  Double_t fgDecayLength = 15.;
  Double_t fgPtMother = 10.;
  Double_t fgDCAPiSVxymax = 0.6;
  Double_t fgDCAPiSVzmax = 0.8;
  Double_t fgDCAProSVmax = 0.7;
  Double_t fgDCADeuSVmax = 0.6;
  
  for(UInt_t j=0; j<cdeuteron.size(); j++){ // candidate deuteron loop
    trackD = cdeuteron.at(j);
    
    trackD->GetImpactParameters(dprim,dprimc);
    dcadprim = TMath::Sqrt((dprim[0]*dprim[0])+(dprim[1]*dprim[1]));

    fHistCorrDCAdprimary->Fill(dprim[0],dprim[1]);

    for(UInt_t m=0; m<cproton.size(); m++){ // candidate proton loop
      trackP = cproton.at(m);
	  
      if(fMC) {if(trackD->GetLabel() == trackP->GetLabel()) continue;}

      trackP->GetImpactParameters(pprim,pprimc);
      dcapprim = TMath::Sqrt((pprim[0]*pprim[0])+(pprim[1]*pprim[1]));
      fHistCorrDCApprimary->Fill(pprim[0],pprim[1]);
      
      dca_dp = trackD->GetDCA(trackP,bz,xthiss,xpp);

      if(dca_dp>dcaDP) continue;
      fHistDCAdeupro->Fill(dca_dp);
      
      for(UInt_t s=0; s<cpion.size(); s++ ){ // candidate pion loop
	trackNPi = cpion.at(s);
	brotherHood = kFALSE;
	
	fHistpionTPCcls->Fill(trackNPi->GetTPCclusters(0));
	fHistpTpion->Fill(trackNPi->Pt());

	if(!fESDtrackCutsV0->AcceptTrack(trackNPi)) continue;

	if(fMC){
	if(trackNPi->GetLabel() == trackP->GetLabel()) continue;
	if(trackNPi->GetLabel() == trackD->GetLabel()) continue;
	}

	trackNPi->GetImpactParameters(piprim,piprimc);
	
	dcapiprim = TMath::Sqrt((piprim[0]*piprim[0])+(piprim[1]*piprim[1]));
	       
	fHistDCApiprimary->Fill(dcapiprim);
	fHistCorrDCApiprimary->Fill(piprim[0],piprim[1]);
	
	if(dcapiprim < fgDCAPiPVmin) continue;
	
	dca_dpi = trackNPi->GetDCA(trackD,bz,xthiss,xpp);
	dca_ppi = trackNPi->GetDCA(trackP,bz,xthiss,xpp);

	if(dca_dpi > GetDCAcut(5,dca_dp)) continue;
	if(dca_ppi > GetDCAcut(4,dca_dp)) continue;

	fHistDCApiondeu->Fill(dca_dpi);
	fHistDCApionpro->Fill(dca_ppi);
	
	trkarray->AddAt(trackD,0);
	trkarray->AddAt(trackP,1);
	trkarray->AddAt(trackNPi,2);

	vertexer->SetVtxStart(esdV1);
	decayVtx = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkarray);

	fHistZDecayVtx->Fill(decayVtx->GetZ());
	fHistXDecayVtx->Fill(decayVtx->GetX());
	fHistYDecayVtx->Fill(decayVtx->GetY());
	decVt[0] = decayVtx->GetX();
	decVt[1] = decayVtx->GetY();
	decVt[2] = decayVtx->GetZ();
	
	dlh[0]=fESDevent->GetPrimaryVertex()->GetX() - decayVtx->GetX();
	dlh[1]=fESDevent->GetPrimaryVertex()->GetY() - decayVtx->GetY();
	dlh[2]=fESDevent->GetPrimaryVertex()->GetZ() - decayVtx->GetZ();

	decayLengthH3L = TMath::Sqrt((dlh[0]*dlh[0]) + (dlh[1]*dlh[1]) + (dlh[2]*dlh[2]));
	
	if(decayLengthH3L > fgDecayLength) continue;
	fHistDecayLengthH3L->Fill(decayLengthH3L);
	
	trackD->PropagateToDCA(decayVtx, bz, 10,dcad);
	fHistDCAXYdeuvtx->Fill(dcad[0]);
	fHistDCAZdeuvtx->Fill(dcad[1]);

	if(TMath::Sqrt((dcad[0]*dcad[0])+(dcad[1]*dcad[1])) > fgDCADeuSVmax) continue;
	
	trackP->PropagateToDCA(decayVtx, bz, 10,dcap);
	fHistDCAXYprovtx->Fill(dcap[0]);
	fHistDCAZprovtx->Fill(dcap[1]);

	if(TMath::Sqrt((dcap[0]*dcap[0])+(dcap[1]*dcap[1])) > fgDCAProSVmax) continue;
	
	trackNPi->PropagateToDCA(decayVtx, bz, 10,dcapi);
	if(TMath::Abs(dcapi[0]) > fgDCAPiSVxymax) continue;
	if(TMath::Abs(dcapi[1]) > fgDCAPiSVzmax) continue;
	
	fHistDCAXYpionvtx->Fill(dcapi[0]);
	fHistDCAZpionvtx->Fill(dcapi[1]);
	
	//if (decayVtx) delete decayVtx;

	posD.SetXYZM(trackD->Px(),trackD->Py(),trackD->Pz(),deuteronMass);
	posP.SetXYZM(trackP->Px(),trackP->Py(),trackP->Pz(),protonMass);
	negPi.SetXYZM(trackNPi->Px(),trackNPi->Py(),trackNPi->Pz(),pionMass);

	pTmom = TMath::Sqrt((trackD->Pt()*trackD->Pt())+(trackP->Pt()*trackP->Pt())+(trackNPi->Pt()*trackNPi->Pt()));
	if(pTmom > fgPtMother) continue;
	
	TLorentzVector Hypertriton;
	TVector3 h1;

	Hypertriton=posD+posP+negPi;
	rapidity = Hypertriton.Rapidity();
	h1.SetXYZ(-dlh[0],-dlh[1],-dlh[2]);
	pointingAngleH = Hypertriton.Angle(h1);

	if (TMath::Cos(pointingAngleH) < fgCosPointingAngle) continue;
	fHistCosPointingAngle->Fill(TMath::Cos(pointingAngleH));

	

	fHistMassHypertriton->Fill(Hypertriton.M());

	
	if(fMC){
	  lD = trackD->GetLabel();
	  lP = trackP->GetLabel();
	  lPi = trackNPi->GetLabel();
	  TParticle *tparticleD = stack->Particle(TMath::Abs(lD));
	  labelM = tparticleD->GetFirstMother();
	  TParticle *tparticleMother = stack->Particle(TMath::Abs(labelM));
	  ndau = tparticleMother->GetNDaughters();
	  
	  if(tparticleMother->GetPdgCode() == 1010010030 && ndau == 3){
	    Int_t labelFirstDau = tparticleMother->GetDaughter(0);
	    Int_t labelSecondDau = labelFirstDau + 1;
	    Int_t labelThirdDau = tparticleMother->GetDaughter(1);
	    if(labelFirstDau == lD && labelSecondDau == lP && labelThirdDau == lPi){
	      //signature of dca pion-primary in case of MCtruth
	      dca_dp = -dca_dp;
	      dca_dpi = -dca_dpi;
	      dca_ppi = -dca_ppi;
	      brotherHood = kTRUE;
	    }
	    else brotherHood = kFALSE;
	  }
	  else brotherHood = kFALSE;
	}
	
	if(fFillTree){
	  fTchi2NDFdeu = trackD->GetTPCchi2()/trackD->GetTPCclusters(0);
	  fTPCclsdeu = trackD->GetTPCclusters(0);
	  fTPCclsPIDdeu = trackD->GetTPCsignalN() ;
	  fTpTPCdeu = trackD->GetTPCmomentum();
	  fTpTdeu = trackD->Pt();
	  fTpdeu = trackD->P();
	  fTTPCnsigmadeu = fPIDResponse->NumberOfSigmasTPC(trackD,AliPID::kDeuteron);
	  fTTOFmassdeu = cmassd.at(j) ;
	  fTDCAXYdeuprvtx = dprim[0];
	  fTDCAZdeuprvtx = dprim[1];
	  fTchi2NDFpro = trackP->GetTPCchi2()/trackP->GetTPCclusters(0);
	  fTPCclspro = trackP->GetTPCclusters(0);
	  fTPCclsPIDpro = trackP->GetTPCsignalN();
	  fTpTPCpro = trackP->GetTPCmomentum();
	  fTpTpro = trackP->Pt();
	  fTppro = trackP->P();
	  fTTPCnsigmapro = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kProton);
	  fTTOFmasspro = cmassp.at(m);
	  fTDCAXYproprvtx = pprim[0];
	  fTDCAZproprvtx = pprim[1];
	  fTchi2NDFpion = trackNPi->GetTPCchi2()/trackNPi->GetTPCclusters(0);
	  fTPCclspion = trackNPi->GetTPCclusters(0);
	  fTPCclsPIDpion = trackNPi->GetTPCsignalN();
	  fTpTPCpion = trackNPi->GetTPCmomentum();
	  fTpTpion = trackNPi->Pt();
	  fTppion = trackNPi->P();
	  fTTPCnsigmapion = fPIDResponse->NumberOfSigmasTPC(trackNPi,AliPID::kPion);
	  fTDCAXYpioprvtx = piprim[0];
	  fTDCAZpioprvtx = piprim[1];
	  fTDCAdp = dca_dp;
	  fTDCAdpi = dca_dpi;
	  fTDCAppi = dca_ppi;
	  fTDCAXYdvtx = dcad[0];
	  fTDCAZdvtx = dcad[1];
	  fTDCAXYpvtx = dcap[0];
	  fTDCAZpvtx = dcap[1];
	  fTDCAXYpivtx = dcapi[0];
	  fTDCAZpivtx = dcapi[1];
	  fTDecayLength = decayLengthH3L;
	  fTCosPA = TMath::Cos(pointingAngleH);
	  fTInvariantMass = Hypertriton.M();

	  fTTree->Fill();
	  PostData(2,fTTree);
	} //end of Fill Tree

	//Pure MC part - Hypertriton truth

	if(fMC){
	  if(brotherHood) {// H3L --> d+p+pi 
	    fHistpionTPCclsMCt->Fill(trackNPi->GetTPCclusters(0));
	    fHistpTpionMCt->Fill(trackNPi->Pt());
	    fHistpTproMCt->Fill(trackP->Pt());
	    fHistpTdeuMCt->Fill(trackD->Pt());
	    fHistDCAdeuproMCt->Fill(-dca_dp);
	    fHistDCApiondeuMCt->Fill(-dca_dpi);
	    fHistDCApionproMCt->Fill(-dca_ppi);
	    fHistCorrDCApiprimaryMCt->Fill(piprim[0],piprim[1]);
	    fHistCorrDCApprimaryMCt->Fill(pprim[0],pprim[1]);
	    fHistCorrDCAdprimaryMCt->Fill(dprim[0],dprim[1]);
	    fHistDCApiprimaryMCt->Fill(dcapiprim);
	    fHistDCApprimaryMCt->Fill(dcapprim);
	    fHistDCAdprimaryMCt->Fill(dcadprim);
	    fHistDCAXYdeuvtxMCt->Fill(dcad[0]);
	    fHistDCAZdeuvtxMCt->Fill(dcad[1]);
	    fHistDCAXYprovtxMCt->Fill(dcap[0]);
	    fHistDCAZprovtxMCt->Fill(dcap[1]);
	    fHistDCAXYpionvtxMCt->Fill(dcapi[0]);
	    fHistDCAZpionvtxMCt->Fill(dcapi[1]);
	    fHistZDecayVtxMCt->Fill(decVt[2]);
	    fHistXDecayVtxMCt->Fill(decVt[0]);
	    fHistYDecayVtxMCt->Fill(decVt[1]);
	    fHistMassHypertritonMCt->Fill(Hypertriton.M());
	  }
	} // end of Pure MC part	
      } // end of candidate pion loop
    } // end of candidate proton loop
  }// end of candidate deuteron loop



  if(fMC){
    for(Int_t iMC=0; iMC<stack->GetNtrack(); iMC++){ // check MC stack content
      TParticle *pstack = stack->Particle(iMC);
      if(pstack->GetPdgCode() == 11) fHistParticle->Fill(0); // e-
      if(pstack->GetPdgCode() == -11) fHistParticle->Fill(1); // e+
      if(pstack->GetPdgCode() == 211) fHistParticle->Fill(2); // pi+
      if(pstack->GetPdgCode() == -211) fHistParticle->Fill(3); // pi-
      if(pstack->GetPdgCode() == 321) fHistParticle->Fill(4); // k+
      if(pstack->GetPdgCode() == -321) fHistParticle->Fill(5); // k-
      if(pstack->GetPdgCode() == 2212) fHistParticle->Fill(6); // p
      if(pstack->GetPdgCode() == -2212) fHistParticle->Fill(7); // pbar
      if(pstack->GetPdgCode() == 1000010020) fHistParticle->Fill(8); // d
      if(pstack->GetPdgCode() == -1000010020) fHistParticle->Fill(9); // dbar
      if(pstack->GetPdgCode() == 1000010030) fHistParticle->Fill(10); // t
      if(pstack->GetPdgCode() == -1000010030) fHistParticle->Fill(11); // tbar
      if(pstack->GetPdgCode() == 1000020030) fHistParticle->Fill(12); // He3
      if(pstack->GetPdgCode() == -1000020030) fHistParticle->Fill(13); // He3bar
      if(pstack->GetPdgCode() == 1000020040) fHistParticle->Fill(14); // He4
      if(pstack->GetPdgCode() == -1000020040) fHistParticle->Fill(15); // He4bar
      if(pstack->GetPdgCode() == 1010010030) fHistParticle->Fill(16); // H3L
      if(pstack->GetPdgCode() == -1010010030) fHistParticle->Fill(17); // H3Lbar
      if(pstack->GetPdgCode() == 1010010030 && pstack->GetNDaughters() == 2) fHistParticle->Fill(18); // H3L
      if(pstack->GetPdgCode() == -1010010030 && pstack->GetNDaughters() == 2) fHistParticle->Fill(19); // H3Lbar
      if(pstack->GetPdgCode() == 3122 && pstack->GetNDaughters() == 2) fHistParticle->Fill(20); // Lambda
    } //end of MC content truth check
  }
  
  // Post output data.
  PostData(1,fOutput);
  if(fFillTree) PostData(2,fTTree);
  
  if(vertexer) delete vertexer;
  if(trkarray) delete trkarray;
  if(esdV1) delete esdV1; 
  
} // end of UserExec

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::Terminate(Option_t *){
  // Merge output
  // Called once at the end of the query

  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  printf("end of Terminate");
  return;

} // end of Terminate
