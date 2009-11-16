/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
// Class AliAnalysisTaskSPDdNdEta                                       //
// Analysis task for dN/dEta reconstruction with the SPD                //
//                                                                      //
// Author:  M. Nicassio (INFN Bari)                                     //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it         //
//////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h" 
#include "TH3F.h" 

#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliTrackReference.h"

#include "AliGenEventHeader.h" 
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliAnalysisTaskSPDdNdEta.h"


ClassImp(AliAnalysisTaskSPDdNdEta)

//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::AliAnalysisTaskSPDdNdEta(const char *name) 
  : AliAnalysisTask(name, "SPDdNdEtaTask"), 

  fESD(0), 
  fOutput(0), 
  
  fpythia(kTRUE),
  fCorr(kFALSE), 
  fTrigger(0),

  // Data to be corrected
  fHistSPDRAWMultvsZ(0),
  fHistSPDRAWMultvsZTriggEvts(0),
  fHistSPDRAWEtavsZ(0),

  // Clusters inner layer and tracklets
  fHistSPDmultEtacut(0),
  fHistSPDmult(0),
  fHistSPDeta(0),               
  fHistSPDcl1multEtacutLay1(0),
  fHistSPDcl1mult(0),
  fHistSPDcl1eta(0),
  fHistSPDphi(0),
  fHistSPDcl1phi(0),
  fHistSPDtheta(0),
  fHistSPDcl1theta(0),
  fHistSPDdePhi(0),
  fHistSPDdePhiZ(0),
  fHistSPDdePhi3D(0),
  fHistSPDphivsSPDeta(0),
  fHistSPDcl1phivsSPDcl1eta(0),
  fHistSPDdeTheta(0),

  // SPD vertex 
  fHistSPDvtxAnalysis(0),                 
  fHistSPDvtx3D(0),
  fHistSPDvtxZ(0),
  fHistNcontribSPDvtxvsSPDvtx(0),
  fHistNcontribSPDvtx3D(0),
  fHistNcontribSPDvtxZ(0),
  fHistNcontribSPDvtxall(0),

  // SPD fired chips  
  fHistSPDcl1multvsnFiredChipsLay1(0),
  fHistSPDmultvsnFiredChipsLay1(0),
  fHistSPDmultvsnFiredChipsLay2(0),
  fHistnFiredChipsLay2vsnFiredChipsLay1(0),
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec(0),
  fHistSPDvtxRec(0),

  // Track level correction histograms
  fHistBkgCorrNum(0),   
  fHistBkgCorrDen(0),

  fHistAlgEffNum(0),

  fHistNonDetectableCorrNum(0),  
  fHistNonDetectableCorrDen(0),

  fHistTrackTrigVtxCorrNum(0),   
  fHistTrackTrigCorrDen(0),      

  fHistTrackTrigVtxCorrNumNSD(0), 
  fHistTrackTrigNSD(0),

  // Event level correction histograms
  fHistTrigVtxCorrNum(0),           
  fHistTrigVtxCorrDen(0),           
  
  fHistTrigCorrDen(0),              
  fHistTrigVtxCorrNumNSD(0),
  fHistEvTrigNSD(0),

  // MC distributions
  fHistMCEtavsZTriggMCvtxEvts(0),
  fHistMCEtavsZTriggESDvtxEvts(0),
  fHistMCEtavsZ(0),                

  // MC dN/dEta for each process type
  fHistMCEtaInel(0),
  fHistMCEtaNonDiffractive(0), 
  fHistMCEtaNonSingleDiffractive(0), 
  fHistoProcessType(0),
  fHistoProcessTypeTriggered(0),
  
  // Additional check histos
  fHistContributorsvsMCVtx(0),
  fHistoDetectableNotr(0),
  fHistoDetectabletr(0),
  fHistoNonStoppingTracks(0),
  fHistoDetectedLay1(0),
  fHistoDetectedLay2(0),
  fHistoPt(0),
  fHistoRTRm1(0),
  fHistMCvtx(0),

  // Trigger efficiency vs MC multiplicity 
  fHistMultAllNonDiff(0), 

  fHistMultAllSingleDiff(0),
  fHistMultAllDoubleDiff(0),
  fHistMultTrVtxNonDiff(0),
  fHistMultTrVtxSingleDiff(0),
  fHistMultTrVtxDoubleDiff(0),

  fHistMCEtaNonSingleDiffractiveLargeBin(0) 

{

  // Constructor

  // Define input and output slots here

  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());  

}
//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::~AliAnalysisTaskSPDdNdEta()
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
void AliAnalysisTaskSPDdNdEta::ConnectInputData(Option_t *) 
{

  // Connect ESD 
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));   

  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else fESD = esdH->GetEvent();
  }
  // Disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);  

}

//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::CreateOutputObjects()
{

  // Create histograms
  fOutput = new TList;
  fOutput->SetOwner();

  const Int_t nBinMultCorr = 28;
  Double_t binLimMultCorr[nBinMultCorr+1];
  for (Int_t i = 0; i<nBinMultCorr+1; ++i) {
    if (i<21) binLimMultCorr[i] = (Double_t) i;
    else if (i<27) binLimMultCorr[i] = (Double_t) 20 +5*(i-20);
    else if (i==27) binLimMultCorr[i] = 100;
    else if (i==28) binLimMultCorr[i] = 200; 
  }

  // Event level correction   
  fHistSPDRAWMultvsZ= new TH2F("fHistSPDRAWMultvsZ", "",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  fHistSPDRAWMultvsZ->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZ);

  fHistSPDRAWMultvsZTriggEvts = new TH2F("fHistSPDRAWMultvsZTriggEvts", "",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  fHistSPDRAWMultvsZTriggEvts->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZTriggEvts);

  fHistSPDRAWEtavsZ = new TH2F("fHistSPDRAWEtavsZ", "Tracklet pseudorapidity distribution", 120, -3.,3.,40,-20.,20.);
  fHistSPDRAWEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDRAWEtavsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWEtavsZ->Sumw2();
  fOutput->Add(fHistSPDRAWEtavsZ);
  
  fHistSPDmultEtacut = new TH1F("fHistSPDmultEtacut", "Tracklet multiplicity distribution",201,-0.5,200.5);
  fHistSPDmultEtacut->GetXaxis()->SetTitle("Reconstructed multiplicity (|#eta|<1.5)");
  fHistSPDmultEtacut->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultEtacut);

  fHistSPDmult = new TH1F("fHistSPDmult", "Tracklet multiplicity distribution", 201,-0.5,200.5);
  fHistSPDmult->GetXaxis()->SetTitle("Reconstructed tracklet multiplicity");
  fHistSPDmult->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmult);

  fHistSPDeta = new TH1F("fHistSPDeta", "Tracklet pseudorapidity distribution", 120, -3.,3.);
  fHistSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDeta->GetYaxis()->SetTitle("Entries");
  fHistSPDeta->SetLineColor(kGreen);
  fHistSPDeta->SetLineWidth(3);
  fHistSPDeta->Sumw2();
  fOutput->Add(fHistSPDeta); 

  fHistSPDcl1multEtacutLay1 = new TH1F("fHistSPDcl1multEtacutLay1", "Cluster multiplicity (inner layer)",201,-0.5,200.5);
  fHistSPDcl1multEtacutLay1->GetXaxis()->SetTitle("Cluster multiplicity lay1 (|#eta|<2.)");
  fHistSPDcl1multEtacutLay1->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1multEtacutLay1);

  fHistSPDcl1mult = new TH1F("fHistSPDcl1mult", "Cluster multiplicity (inner layer)",201,-0.5,200.5);
  fHistSPDcl1mult->GetXaxis()->SetTitle("Cluster multiplicity lay1");
  fHistSPDcl1mult->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1mult);

  fHistSPDcl1eta  = new TH1F("fHistSPDcl1eta", "Cluster pseudorapidity (inner layer)", 120, -3.,3.);
  fHistSPDcl1eta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDcl1eta->GetYaxis()->SetTitle("Entries");
  fHistSPDcl1eta->Sumw2();
  fOutput->Add(fHistSPDcl1eta);

  fHistSPDphi = new TH1F("fHistSPDphi", "Tracklet #phi  distribution", 360, 0.,2*TMath::Pi());
  fHistSPDphi->GetXaxis()->SetTitle("#varphi [rad]");
  fHistSPDphi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDphi);

  fHistSPDcl1phi= new TH1F("fHistSPDcl1phi", "Cluster #phi (inner layer) ", 360, 0.,2*TMath::Pi());
  fHistSPDcl1phi->GetXaxis()->SetTitle("#varphi [rad]");
  fHistSPDcl1phi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1phi);

  fHistSPDtheta = new TH1F("fHistSPDtheta", "Tracklet #theta  distribution", 180, 0.,TMath::Pi());
  fHistSPDtheta->GetXaxis()->SetTitle("#theta [rad]");
  fHistSPDtheta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDtheta);

  fHistSPDcl1theta = new TH1F("fHistSPDcl1theta", "Cluster #theta (inner layer)", 180, 0.,TMath::Pi());
  fHistSPDcl1theta->GetXaxis()->SetTitle("#theta [rad]");
  fHistSPDcl1theta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1theta);

  fHistSPDdePhi= new TH1F("fHistSPDdePhi", "Tracklet #Delta#varphi  distribution",400,-0.1,.1);
  fHistSPDdePhi->GetXaxis()->SetTitle("#Delta#varphi [rad]");
  fHistSPDdePhi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDdePhi);

  fHistSPDdePhiZ= new TH1F("fHistSPDdePhiZ", "Tracklet #Delta#varphi  distribution",400,-0.1,.1);
  fOutput->Add(fHistSPDdePhiZ);

  fHistSPDdePhi3D= new TH1F("fHistSPDdePhi3D", "Tracklet #Delta#varphi  distribution",400,-0.1,.1);
  fOutput->Add(fHistSPDdePhi3D);

  fHistSPDphivsSPDeta= new TH2F("fHistSPDphivsSPDeta", "Tracklets - #varphi vs #eta",120,-3.,3,360,0.,2*TMath::Pi());
  fHistSPDphivsSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDphivsSPDeta->GetYaxis()->SetTitle("#varphi [rad]");
  fOutput->Add(fHistSPDphivsSPDeta);

  fHistSPDcl1phivsSPDcl1eta= new TH2F("fHistSPDcl1phivsSPDcl1eta", "Clusters layer1 - #varphi vs #eta",120,-3.,3,360,0.,2*TMath::Pi());
  fHistSPDcl1phivsSPDcl1eta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDcl1phivsSPDcl1eta->GetYaxis()->SetTitle("#varphi [rad]");
  fOutput->Add(fHistSPDcl1phivsSPDcl1eta);

  fHistSPDdeTheta= new TH1F("fHistSPDdeTheta", "Tracklet #Delta#theta distribution",100,-0.05,.05);
  fHistSPDdeTheta->GetXaxis()->SetTitle("#Delta#theta [rad]");
  fHistSPDdeTheta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDdeTheta);


  fHistSPDvtxAnalysis = new TH1F("fHistSPDvtxAnalysis", "SPD vertex  distribution - all events",20,-20.,20.);
  fHistSPDvtxAnalysis->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtxAnalysis->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtxAnalysis);

  fHistSPDvtx3D = new TH3F("fHistSPDvtx3D", "SPD vertex distribution",100,-1.,1.,100,-1.,1.,500,-50.,50.);
  fOutput->Add(fHistSPDvtx3D);

  fHistSPDvtxZ = new TH1F("fHistSPDvtxZ", "SPD vertex distribution",500,-50.,50.);
  fOutput->Add(fHistSPDvtxZ);

  fHistNcontribSPDvtxvsSPDvtx= new TH2F("fHistNcontribSPDvtxvsSPDvtx", " ",100,-50.,50.,10002,-2.,10000.);
  fHistNcontribSPDvtxvsSPDvtx->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistNcontribSPDvtxvsSPDvtx->GetYaxis()->SetTitle("# contributors");
  fOutput->Add(fHistNcontribSPDvtxvsSPDvtx);

  fHistNcontribSPDvtx3D= new TH1F("fHistNcontribSPDvtx3D", "SPD vtx 3D",10002,-2.,10000.);
  fHistNcontribSPDvtx3D->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtx3D->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtx3D);

  fHistNcontribSPDvtxZ= new TH1F("fHistNcontribSPDvtxZ", "SPD vtx Z",10002,-2.,10000.);
  fHistNcontribSPDvtxZ->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtxZ->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtxZ);

  fHistNcontribSPDvtxall= new TH1F("fHistNcontribSPDvtxall", "SPD vtx - all events",10002,-2.,10000.);
  fHistNcontribSPDvtxall->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtxall->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtxall);

  fHistSPDcl1multvsnFiredChipsLay1 = new TH2F("fHistSPDcl1multvsnFiredChipsLay1", "",401,0.,401.,201,-0.5,200.5);
  fHistSPDcl1multvsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistSPDcl1multvsnFiredChipsLay1->GetYaxis()->SetTitle("Cluster lay1 multiplicity");
  fOutput->Add(fHistSPDcl1multvsnFiredChipsLay1);

  fHistSPDmultvsnFiredChipsLay1 = new TH2F("fHistSPDmultvsnFiredChipsLay1","",401,0.,401.,201,-0.5,200.5);
  fHistSPDmultvsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistSPDmultvsnFiredChipsLay1->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvsnFiredChipsLay1);

  fHistSPDmultvsnFiredChipsLay2 = new TH2F("fHistSPDmultvsnFiredChipsLay2","",801,0.,801.,201,-0.5,200.5);
  fHistSPDmultvsnFiredChipsLay2->GetXaxis()->SetTitle("# fired chips lay2");
  fHistSPDmultvsnFiredChipsLay2->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvsnFiredChipsLay2);

  fHistnFiredChipsLay2vsnFiredChipsLay1 = new TH2F("fHistnFiredChipsLay2vsnFiredChipsLay1","",401,0.,401.,801,0.,801.);
  fHistnFiredChipsLay2vsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistnFiredChipsLay2vsnFiredChipsLay1->GetYaxis()->SetTitle("# fired chip lay2");
  fOutput->Add(fHistnFiredChipsLay2vsnFiredChipsLay1);

  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec = new TH2F("fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec","",401,0.,401.,801,0.,801.);
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec->GetXaxis()->SetTitle("# fired chips lay1");
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec->GetYaxis()->SetTitle("# fired chip lay2");
  fOutput->Add(fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec);

  fHistSPDvtxRec = new TH1F("fHistSPDvtxRec", "SPD vertex  distribution",80,-20.,20.);
  fHistSPDvtxRec->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtxRec->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtxRec);


  if (fCorr) {
    fHistBkgCorrNum = new TH2F("fHistBkgCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistBkgCorrNum);

    fHistBkgCorrDen = new TH2F("fHistBkgCorrDen","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistBkgCorrDen);

    fHistAlgEffNum = new TH2F("fHistAlgEffNum","",120,-3.00,3.00,40,-20.,20.);
    fOutput->Add(fHistAlgEffNum);  

    fHistNonDetectableCorrNum = new TH2F("fHistNonDetectableCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistNonDetectableCorrNum);

    fHistNonDetectableCorrDen = new TH2F("fHistNonDetectableCorrDen","",120,-3.00,3.00,40,-20.,20.);
    fOutput->Add(fHistNonDetectableCorrDen);

    fHistTrackTrigVtxCorrNum = new TH2F("fHistTrackTrigVtxCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistTrackTrigVtxCorrNum);

    fHistTrackTrigCorrDen = new TH2F("fHistTrackTrigCorrDen","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistTrackTrigCorrDen);

    fHistTrackTrigVtxCorrNumNSD = new TH2F("fHistTrackTrigVtxCorrNumNSD","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistTrackTrigVtxCorrNumNSD);

    fHistTrackTrigNSD = new TH2F("fHistTrackTrigNSD","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistTrackTrigNSD);


    // Event level correction histograms  
    fHistTrigVtxCorrNum = new TH2F("fHistTrigVtxCorrNum","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistTrigVtxCorrNum);

    fHistTrigVtxCorrDen = new TH2F("fHistTrigVtxCorrDen","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistTrigVtxCorrDen);

    fHistTrigCorrDen = new TH2F("fHistTrigCorrDen","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistTrigCorrDen);

    fHistTrigVtxCorrNumNSD = new TH2F("fHistTrigVtxCorrNumNSD","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistTrigVtxCorrNumNSD);

    fHistEvTrigNSD = new TH2F("fHistEvTrigNSD","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistEvTrigNSD);

    // MC distributions
    fHistMCEtavsZTriggMCvtxEvts = new TH2F("fHistMCEtavsZTriggMCvtxEvts","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZTriggMCvtxEvts->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZTriggMCvtxEvts->GetYaxis()->SetTitle("z_{MCvtx} [cm]");
    fHistMCEtavsZTriggMCvtxEvts->Sumw2();
    fOutput->Add(fHistMCEtavsZTriggMCvtxEvts);

    fHistMCEtavsZTriggESDvtxEvts = new TH2F("fHistMCEtavsZTriggESDvtxEvts","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZTriggESDvtxEvts->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZTriggESDvtxEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
    fHistMCEtavsZTriggESDvtxEvts->Sumw2();
    fOutput->Add(fHistMCEtavsZTriggESDvtxEvts);

    fHistMCEtavsZ = new TH2F("fHistMCEtavsZ","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZ->GetYaxis()->SetTitle("z_{MCvtx} [cm]");
    fHistMCEtavsZ->Sumw2();
    fOutput->Add(fHistMCEtavsZ);

    fHistMCEtaInel = new TH1F("fHistMCEtaInel","",7,-3.5,3.5);
    fHistMCEtaInel->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fOutput->Add(fHistMCEtaInel);

    fHistMCEtaNonDiffractive = new TH1F("fHistMCEtaNonDiffractive","",60,-3.,3.);
    fHistMCEtaNonDiffractive->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fOutput->Add(fHistMCEtaNonDiffractive);

    fHistMCEtaNonSingleDiffractive = new TH1F("fHistMCEtaNonSingleDiffractive","",60,-3.,3.);
    fHistMCEtaNonSingleDiffractive->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fOutput->Add(fHistMCEtaNonSingleDiffractive);


    
    fHistoProcessType = new TH1F("fHistoProcessType","",5,0.,5.);
    fOutput->Add(fHistoProcessType);

    fHistoProcessTypeTriggered = new TH1F("fHistoProcessTypeTriggered","",5,0.,5);
    fOutput->Add(fHistoProcessTypeTriggered);

    fHistContributorsvsMCVtx = new TH2F("fHistContributorsvsMCVtx","",200,-20.,20.,202,-2.,200.);
    fOutput->Add(fHistContributorsvsMCVtx);

    fHistoDetectableNotr = new TH3F("fHistoDetectableNotr","",60,-3.00,3.00,20,-20.,20.,100,0.,10.); 
    fOutput->Add(fHistoDetectableNotr);
 
    fHistoDetectabletr = new TH2F("fHistoDetectabletr","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectabletr);

    fHistoNonStoppingTracks = new TH2F("fHistoNonStoppingTracks","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoNonStoppingTracks);

    fHistoDetectedLay1 = new TH2F("fHistoDetectedLay1","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectedLay1);

    fHistoDetectedLay2 = new TH2F("fHistoDetectedLay2","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectedLay2);

    fHistoPt = new TH1F("fHistoPt","",100,.0,10.);
    fOutput->Add(fHistoPt);

    fHistoRTRm1 = new TH1F("fHistoRTRm1","",10000,0.,5000);
    fOutput->Add(fHistoRTRm1);

    fHistMCvtx = new TH3F("fHistMCvtx", "MC vertex distribution",100,-.5,.5,100,-.5,.5,500,-50.,50.);
    fOutput->Add(fHistMCvtx);

    fHistMultAllNonDiff  = new TH1F("fHistMultAllNonDiff","",20,-0.5,19.5);
    fOutput->Add(fHistMultAllNonDiff);
    fHistMultAllSingleDiff = new TH1F("fHistMultAllSingleDiff","",20,-0.5,19.5);
    fHistMultAllSingleDiff->Sumw2();
    fOutput->Add(fHistMultAllSingleDiff);
    fHistMultAllDoubleDiff = new TH1F("fHistMultAllDoubleDiff","",20,-0.5,19.5);
    fHistMultAllDoubleDiff->Sumw2();
    fOutput->Add(fHistMultAllDoubleDiff);
    fHistMultTrVtxNonDiff = new TH1F("fHistMultTrVtxNonDiff","",20,-0.5,19.5);
    fHistMultTrVtxNonDiff->Sumw2();
    fOutput->Add(fHistMultTrVtxNonDiff);
    fHistMultTrVtxSingleDiff = new TH1F("fHistMultTrVtxSingleDiff","",20,-0.5,19.5);
    fHistMultTrVtxSingleDiff->Sumw2();
    fOutput->Add(fHistMultTrVtxSingleDiff);
    fHistMultTrVtxDoubleDiff = new TH1F("fHistMultTrVtxDoubleDiff","",20,-0.5,19.5);
    fHistMultTrVtxDoubleDiff->Sumw2();
    fOutput->Add(fHistMultTrVtxDoubleDiff);
   
    fHistMCEtaNonSingleDiffractiveLargeBin = new TH1F("fHistMCEtaNonSingleDiffractiveLargeBin","",7,-3.5,3.5);
    fHistMCEtaNonSingleDiffractiveLargeBin->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fOutput->Add(fHistMCEtaNonSingleDiffractiveLargeBin);

  }
}

//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::Exec(Option_t *) 
{
  // Main loop
  // Called for each event


  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // ESD vertex
  Bool_t eventWithVertex = kFALSE;
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  Int_t nContrib = vtxESD->GetNContributors();
  //... check resolution
  Double_t zRes = vtxESD->GetZRes();
  const AliMultiplicity* multESD = fESD->GetMultiplicity();

  // Loading tracklets...
  Int_t multSPD = 0;
  multSPD = multESD->GetNumberOfTracklets();

  AliMultiplicity * mult1 = (AliMultiplicity*)multESD;
  Short_t nFiredChipsLay1 = mult1->GetNumberOfFiredChips(0);
  Short_t nFiredChipsLay2 = mult1->GetNumberOfFiredChips(1);
  Int_t multSPDEtacut = 0;
  Int_t multSPDcl1 = 0;
  Int_t nSingleCl1 = 0;
  Int_t multSPDcl1EtacutLay1 = 0;
  nSingleCl1 = multESD->GetNumberOfSingleClusters();
  multSPDcl1 = nSingleCl1 + multSPD;
  Float_t* recEtaSPDcl1 = new Float_t[multSPD+nSingleCl1];
  if (esdvtx[2]!=0.) fHistSPDvtxRec->Fill(esdvtx[2]);

  if (esdvtx[2]!=0.&&zRes<0.1) eventWithVertex = kTRUE;
  else {
    esdvtx[2] = 0.;
    multSPD = 0;
  }

//  Printf("There are %d tracklets in this event", multSPD);

  // Trigger
  Bool_t eventTriggered = kFALSE;
  ULong64_t triggerMask;
  ULong64_t spdFO   = (1 << 14);
  ULong64_t v0left  = (1 << 11);
  ULong64_t v0right = (1 << 12);

  triggerMask=fESD->GetTriggerMask();

  // No trigger
  if (fTrigger==0) eventTriggered = kTRUE;
  // MB1: GFO || V0OR
  if (fTrigger==1) eventTriggered = triggerMask&spdFO || ((triggerMask&v0left) || (triggerMask&v0right));
  // MB2: GFO && V0OR
  if (fTrigger==2) eventTriggered = triggerMask&spdFO && ((triggerMask&v0left) || (triggerMask&v0right));

  PostData(0, fOutput);


  // Selected events: triggered with vertex
  if (eventTriggered&&eventWithVertex) {

    if (multSPD!=0) fHistSPDvtxAnalysis->Fill(esdvtx[2]);
    if (strcmp(vtxESD->GetTitle(),"vertexer: 3D") == 0) fHistSPDvtx3D->Fill(esdvtx[0],esdvtx[1],esdvtx[2]);
    else fHistSPDvtxZ->Fill(esdvtx[2]);

    for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
      Float_t thetaTr= multESD->GetTheta(itracklet);
      Float_t phiTr= multESD->GetPhi(itracklet);
      Float_t dePhiTr= multESD->GetDeltaPhi(itracklet);
//      Double_t deThetaTr= multESD->GetDeltaTheta(itracklet);
      Float_t recEtaSPD =multESD->GetEta(itracklet);
      recEtaSPDcl1[itracklet] = recEtaSPD;

      fHistSPDeta->Fill(recEtaSPD);
      fHistSPDRAWEtavsZ->Fill(recEtaSPD,esdvtx[2]);
      fHistSPDcl1eta->Fill(recEtaSPD);
      fHistSPDphi->Fill(phiTr);
      fHistSPDcl1phi->Fill(phiTr);
      fHistSPDtheta->Fill(thetaTr);
      fHistSPDcl1theta->Fill(thetaTr);
      fHistSPDdePhi->Fill(dePhiTr);
//      fHistSPDdeTheta->Fill(deThetaTr);
 
      if (strcmp(vtxESD->GetTitle(),"vertexer: Z") == 0) fHistSPDdePhiZ->Fill(dePhiTr);
      if (strcmp(vtxESD->GetTitle(),"vertexer: 3D") == 0) fHistSPDdePhi3D->Fill(dePhiTr);
      fHistSPDphivsSPDeta->Fill(recEtaSPD,phiTr);
      fHistSPDcl1phivsSPDcl1eta->Fill(recEtaSPD,phiTr);

      // Calculate multiplicity in etacut
      if (TMath::Abs(recEtaSPD)<1.5) multSPDEtacut++;
      if (TMath::Abs(recEtaSPD)<2.)  multSPDcl1EtacutLay1++;
    }

    for (Int_t iCl1=0; iCl1<nSingleCl1; ++iCl1) {
      Float_t thetaSingleCl1 = multESD->GetThetaSingle(iCl1);
      // Calculate eta
      Float_t etaSingleCl1 = -TMath::Log(TMath::Tan(thetaSingleCl1/2.));
      Float_t phiSingleCl1 = multESD->GetPhiSingle(iCl1);
      recEtaSPDcl1[iCl1+multSPD] = etaSingleCl1;

      fHistSPDcl1eta->Fill(etaSingleCl1);
      fHistSPDcl1phi->Fill(phiSingleCl1);
      fHistSPDcl1theta->Fill(thetaSingleCl1);
      fHistSPDcl1phivsSPDcl1eta->Fill(etaSingleCl1,phiSingleCl1);
      if (TMath::Abs(etaSingleCl1)<2.) multSPDcl1EtacutLay1++;
    }

    fHistSPDmultEtacut->Fill(multSPDEtacut);
    fHistSPDmult->Fill(multSPD);
    fHistSPDcl1multEtacutLay1->Fill(multSPDcl1EtacutLay1);
    if (strcmp(vtxESD->GetTitle(),"vertexer: 3D") == 0) {
      fHistNcontribSPDvtx3D->Fill(nContrib);
    }
    if (strcmp(vtxESD->GetTitle(),"vertexer: Z") == 0)  {
     fHistNcontribSPDvtxZ->Fill(nContrib);
    }
    fHistNcontribSPDvtxvsSPDvtx->Fill(esdvtx[2],nContrib);
    fHistSPDRAWMultvsZ->Fill(multSPD,esdvtx[2]);
    fHistSPDmultvsnFiredChipsLay1->Fill(nFiredChipsLay1,multSPD);
    fHistSPDmultvsnFiredChipsLay2->Fill(nFiredChipsLay2,multSPD);

  } // End selected events

  if (eventTriggered) {
    fHistSPDRAWMultvsZTriggEvts->Fill(multSPD,esdvtx[2]); 
  }

  fHistSPDcl1mult->Fill(multSPDcl1);
  fHistNcontribSPDvtxall->Fill(nContrib);

  fHistSPDcl1multvsnFiredChipsLay1->Fill(nFiredChipsLay1,multSPDcl1);

  fHistnFiredChipsLay2vsnFiredChipsLay1->Fill(nFiredChipsLay1,nFiredChipsLay2);

  if (esdvtx[2]==0.) {  
    fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec->Fill(nFiredChipsLay1,nFiredChipsLay2);
  }
 
  delete[] recEtaSPDcl1;

  if (fCorr) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
       Printf("ERROR: Could not retrieve MC event");
       return;
    }

    AliStack* stack = mcEvent->Stack();
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
    // MC vertex
    TArrayF vtxMC(3);
    genHeader->PrimaryVertex(vtxMC);

    fHistMCvtx->Fill(vtxMC[0],vtxMC[1],vtxMC[2]);   

    //Adding process type selection
    Int_t processType;
    Bool_t nsdEv = kFALSE;
    Bool_t ndEv = kFALSE; 
    Bool_t sdEv = kFALSE;
    Bool_t ddEv = kFALSE;
    Bool_t inelEv = kFALSE;

    AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);

    AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);

    if (fpythia&&pythiaGenHeader) { 
      processType = pythiaGenHeader->ProcessType();
      if (processType!=91&&processType!=92&&processType!=93&&processType!=94)  ndEv = kTRUE; // non difffractive
      if (processType!=91&&processType!=92&&processType!=93)  nsdEv = kTRUE;                 // non single diffractive
      if (processType==92||processType==93)  sdEv = kTRUE;                                   //single diffractive
      if (processType==94)  ddEv = kTRUE;                                                    //double diffractive
      if (processType!=91)  inelEv = kTRUE;                                                  //inelastic
      
    } else if (!fpythia&&dpmHeader) {
      processType = dpmHeader->ProcessType(); 
      if (processType==1)  ndEv = kTRUE;                                         // non diffractive
      if (processType!=2&&processType!=5&&processType!=6)  nsdEv = kTRUE;        // non single diffractive
      if (processType==5||processType==6)  sdEv = kTRUE;                         // single diffractive
      if (processType==4||processType==7)  ddEv = kTRUE;                         // double diffractive
      if (processType!=2)  inelEv = kTRUE;                                       // inelastic 
    } else if (!pythiaGenHeader&&!dpmHeader) {
      printf("Unknown header type: neither DPMjet nor Pythia. \n");
      return ;
    }

    if (ndEv) {
      fHistoProcessType->Fill(0);
      if (eventTriggered) fHistoProcessTypeTriggered->Fill(0);
    } 
    if (nsdEv) {
      fHistoProcessType->Fill(1);
      if (eventTriggered) fHistoProcessTypeTriggered->Fill(1);
    } else if (sdEv) {
      fHistoProcessType->Fill(3);
      if (eventTriggered) fHistoProcessTypeTriggered->Fill(3);
    }
    if (ddEv) {
      fHistoProcessType->Fill(4);
      if (eventTriggered) fHistoProcessTypeTriggered->Fill(4);
    }
    if (inelEv) {
      fHistoProcessType->Fill(2);                          // inel is always true
      if (eventTriggered) fHistoProcessTypeTriggered->Fill(2);
    } 
 
    // Tracks from MC
    Int_t  multMCCharged = 0;
    Int_t  multMCChargedEtacut = 0;
    Int_t  nMCPart = stack->GetNprimary();
    Float_t* etagen = new Float_t[nMCPart];  
    Float_t* ptgen = new Float_t[nMCPart];
    Int_t* stackIndexOfPrimaryParts = new Int_t[nMCPart];
    Int_t* reconstructedPrimaryPart = new Int_t[nMCPart];
    Int_t* detectedPrimaryPartLay1 = new Int_t[nMCPart];
    Int_t* detectedPrimaryPartLay2 = new Int_t[nMCPart];
    Int_t* detectablePrimaryPart = new Int_t[nMCPart];

    // Loading track references...
    Float_t rminL1 = 3.4;
    Float_t rmaxL1 = 4.4;
    Float_t rminL2 = 6.9;
    Float_t rmaxL2 = 7.9;

    TTree* tRefTree = eventHandler->TreeTR();

    AliTrackReference *tref=0x0;
    mcEvent->ConnectTreeTR(tRefTree);

    // Loop over MC particles
    for (Int_t imc=0; imc<nMCPart; imc++) {
      TParticle* part = stack->Particle(imc);
      Bool_t isPrimary = stack->IsPhysicalPrimary(imc);
      if (!isPrimary)                        continue;
      TParticlePDG* pdgPart = part->GetPDG();
      if (TMath::Abs(pdgPart->Charge())!=3)  continue;
      Float_t theta = part->Theta();
      if (theta==0 || theta==TMath::Pi())    continue;
      Float_t eta = part->Eta();
      Float_t pt = part->Pt();
      etagen[multMCCharged] = eta; 
      ptgen[multMCCharged] = pt;
      stackIndexOfPrimaryParts[multMCCharged] = imc;

      reconstructedPrimaryPart[multMCCharged]=kFALSE;
      detectedPrimaryPartLay1[multMCCharged]=kFALSE;
      detectedPrimaryPartLay2[multMCCharged]=kFALSE;
      detectablePrimaryPart[multMCCharged]=kFALSE;

      fHistoPt->Fill(ptgen[multMCCharged]);
      if (ndEv)                                                    // non difffractive
        fHistMCEtaNonDiffractive->Fill(etagen[multMCCharged]);
      if (nsdEv) {                                                 // non single diffractive
        fHistMCEtaNonSingleDiffractive->Fill(etagen[multMCCharged]);
        fHistMCEtaNonSingleDiffractiveLargeBin->Fill(etagen[multMCCharged]);
      }
      // inel is always true     
      fHistMCEtaInel->Fill(etagen[multMCCharged]); 
      

      AliMCParticle* mcpart = (AliMCParticle*) mcEvent->GetTrack(imc);
      Int_t nref = mcpart->GetNumberOfTrackReferences();

      // Detectable primaries 
      if (nref==0) {
        detectablePrimaryPart[multMCCharged]=kTRUE;
        fHistoDetectableNotr->Fill(etagen[multMCCharged],vtxMC[2],ptgen[multMCCharged]);
      } else if (nref>0) {
        tref = mcpart->GetTrackReference(nref-1);
        if (tref->DetectorId()!=-1) {
          detectablePrimaryPart[multMCCharged]=kTRUE;  
          fHistoNonStoppingTracks->Fill(etagen[multMCCharged],vtxMC[2]); 
          for (Int_t iref=0;iref<nref;iref++) { //since it is detectable, check if it is also detected
            tref = mcpart->GetTrackReference(iref);
            if (tref) {
              if (tref->R()>rminL2&&tref->R()<rmaxL2) {
                if (tref->DetectorId()==0) {
                  detectedPrimaryPartLay2[multMCCharged]=kTRUE;
                  fHistoDetectedLay2->Fill(etagen[multMCCharged],vtxMC[2]);
                  break;
                }
              }
            }
          }
        } else { //last is -1 -> particle disappeared. Where?
          tref = mcpart->GetTrackReference(nref-1);
          fHistoRTRm1->Fill(tref->R());
          if (tref->R()>rmaxL2) {
            detectablePrimaryPart[multMCCharged]=kTRUE;
            fHistoDetectabletr->Fill(etagen[multMCCharged],vtxMC[2]);
            for (Int_t iref=0;iref<nref;iref++) { //since it is detectable, check if it is also detected
              tref = mcpart->GetTrackReference(iref);
              if (tref) {
                if (tref->R()>rminL2&&tref->R()<rmaxL2) {
                  if (tref->DetectorId()==0) {
                    detectedPrimaryPartLay2[multMCCharged]=kTRUE;
                    fHistoDetectedLay2->Fill(etagen[multMCCharged],vtxMC[2]);
                    break;
                  }
                }
              }
            }
          } else if (tref->R()>=rminL2&&tref->R()<=rmaxL2) {
            for (Int_t iref=0;iref<nref;iref++) {
              tref = mcpart->GetTrackReference(iref);
              if (tref) {
                if (tref->R()>rminL2&&tref->R()<rmaxL2) {
                  if (tref->DetectorId()==0) {
                    detectablePrimaryPart[multMCCharged]=kTRUE;
                    detectedPrimaryPartLay2[multMCCharged]=kTRUE;
                    fHistoDetectedLay2->Fill(etagen[multMCCharged],vtxMC[2]);
                    fHistoDetectabletr->Fill(etagen[multMCCharged],vtxMC[2]);
                    break;
                  } 
                }
              }
            }
          }
        }
        // Find out detected prims on each layer
        for (Int_t iref=0; iref<nref; iref++) {
          tref = mcpart->GetTrackReference(iref);
          if (tref->R()>rminL1&&tref->R()<rmaxL1&&tref->DetectorId()==0) {
            detectedPrimaryPartLay1[multMCCharged] = kTRUE;
            fHistoDetectedLay1->Fill(etagen[multMCCharged],vtxMC[2]);
          }
        }
      } 

      multMCCharged++;
      if (TMath::Abs(eta)<1.4) multMCChargedEtacut++;
    } // End of MC particle loop


    if (ndEv) { // non diffractive
      fHistMultAllNonDiff->Fill(multMCChargedEtacut);
      if (eventTriggered&&esdvtx[2]!=0.) fHistMultTrVtxNonDiff->Fill(multMCChargedEtacut);

    }
    if (sdEv) { // single diffractive
      fHistMultAllSingleDiff->Fill(multMCChargedEtacut);
      if (eventTriggered&&esdvtx[2]!=0.) fHistMultTrVtxSingleDiff->Fill(multMCChargedEtacut);
    }
    if (ddEv) { // double diffractive
      fHistMultAllDoubleDiff->Fill(multMCChargedEtacut);
      if (eventTriggered&&esdvtx[2]!=0.) fHistMultTrVtxDoubleDiff->Fill(multMCChargedEtacut);
    }
 
    if (esdvtx[2]==0.) {
      fHistContributorsvsMCVtx->Fill(vtxMC[2],nContrib);
    }
    // Event selection
    if (eventTriggered&&eventWithVertex) {

      for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {

        Int_t labL1 = multESD->GetLabel(itracklet,0);
        Int_t labL2 = multESD->GetLabel(itracklet,1);

        fHistBkgCorrDen->Fill(multESD->GetEta(itracklet),esdvtx[2]);

        if (labL1==labL2) {
          for (Int_t imc=0; imc<multMCCharged; imc++) {
            if (labL1==stackIndexOfPrimaryParts[imc]) {
              if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) {
                reconstructedPrimaryPart[imc]=kTRUE;
                break;
              }
            }
          }
        }
      }

      for (Int_t imc=0; imc<multMCCharged; imc++) {
        if (reconstructedPrimaryPart[imc]) {
          fHistBkgCorrNum->Fill(etagen[imc],vtxMC[2]);
        }
        if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) fHistAlgEffNum->Fill(etagen[imc],vtxMC[2]);

        fHistNonDetectableCorrNum->Fill(etagen[imc],vtxMC[2]);

        if (detectablePrimaryPart[imc]) fHistNonDetectableCorrDen->Fill(etagen[imc],vtxMC[2]);

        fHistMCEtavsZTriggESDvtxEvts->Fill(etagen[imc],esdvtx[2]);
        fHistMCEtavsZTriggMCvtxEvts->Fill(etagen[imc],vtxMC[2]);
      }

      fHistTrigVtxCorrDen->Fill(multSPD,vtxMC[2]);
    } // End of selected events

    if (eventTriggered) {
        fHistTrigCorrDen->Fill(multSPD,vtxMC[2]);
        fHistTrigVtxCorrNum->Fill(multSPD,vtxMC[2]);
      if (nsdEv) {
        fHistEvTrigNSD->Fill(multSPD,vtxMC[2]); //to compute errors
        fHistTrigVtxCorrNumNSD->Fill(multSPD,vtxMC[2]);
      }
      for (Int_t imc=0; imc<multMCCharged; imc++) {
        fHistTrackTrigCorrDen->Fill(etagen[imc],vtxMC[2]);
        if (nsdEv) fHistTrackTrigNSD->Fill(etagen[imc],vtxMC[2]); //to compute errors 
      }
    } else {
       fHistTrigVtxCorrNum->Fill(0.,vtxMC[2]); 
       if (nsdEv) fHistTrigVtxCorrNumNSD->Fill(0.,vtxMC[2]); 

    }

    for (Int_t imc=0; imc<multMCCharged; imc++) {
      fHistTrackTrigVtxCorrNum->Fill(etagen[imc],vtxMC[2]);
      fHistMCEtavsZ->Fill(etagen[imc],vtxMC[2]);
      if (nsdEv) 
        fHistTrackTrigVtxCorrNumNSD->Fill(etagen[imc],vtxMC[2]);
    }

    delete[] etagen;
    delete[] ptgen;
    delete[] stackIndexOfPrimaryParts;
    delete[] reconstructedPrimaryPart;
    delete[] detectedPrimaryPartLay1;
    delete[] detectedPrimaryPartLay2;
    delete[] detectablePrimaryPart;
  }
}      

//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::Terminate(Option_t *) 
{
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {     Printf("ERROR: fOutput not available");
    return;
  }
  
  fHistSPDRAWMultvsZ= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDRAWMultvsZ"));
  fHistSPDRAWMultvsZTriggEvts = dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDRAWMultvsZTriggEvts"));
  fHistSPDRAWEtavsZ = dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDRAWEtavsZ"));

  fHistSPDmultEtacut = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDmultEtacut"));
  fHistSPDmult = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDmult"));
  fHistSPDeta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDeta"));
  fHistSPDcl1multEtacutLay1 = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDcl1multEtacutLay1"));
  fHistSPDcl1mult = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDcl1mult"));
  fHistSPDcl1eta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDcl1eta"));
  fHistSPDphi = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDphi"));
  fHistSPDcl1phi = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDcl1phi"));
  fHistSPDtheta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDtheta"));
  fHistSPDcl1theta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDcl1theta"));
  fHistSPDdePhi = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDdePhi"));
  fHistSPDdePhiZ = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDdePhiZ"));
  fHistSPDdePhi3D = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDdePhi3D"));
  fHistSPDphivsSPDeta= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDphivsSPDeta"));
  fHistSPDcl1phivsSPDcl1eta= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDcl1phivsSPDcl1eta"));
  fHistSPDdeTheta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDdeTheta"));

  fHistSPDvtxAnalysis = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDvtxAnalysis"));
  fHistSPDvtx3D = dynamic_cast<TH3F*> (fOutput->FindObject("fHistSPDvtx3D"));
  fHistSPDvtxZ = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDvtxZ"));
  fHistNcontribSPDvtxvsSPDvtx = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNcontribSPDvtxvsSPDvtx"));
  fHistNcontribSPDvtx3D = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtx3D"));
  fHistNcontribSPDvtxZ = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtxZ"));
  fHistNcontribSPDvtxall = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtxall"));

  fHistSPDcl1multvsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDcl1multvsnFiredChipsLay1"));
  fHistSPDmultvsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDmultvsnFiredChipsLay1"));
  fHistSPDmultvsnFiredChipsLay2= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDmultvsnFiredChipsLay2"));
  fHistnFiredChipsLay2vsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistnFiredChipsLay2vsnFiredChipsLay1"));
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec= dynamic_cast<TH2F*> (fOutput->FindObject("fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec"));

  fHistSPDvtxRec = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDvtxRec"));

  if (fCorr) {

    fHistBkgCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistBkgCorrNum"));
    fHistBkgCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistBkgCorrDen"));

    fHistAlgEffNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistAlgEffNum"));

    fHistNonDetectableCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNonDetectableCorrNum"));
    fHistNonDetectableCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNonDetectableCorrDen"));

    fHistTrackTrigVtxCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigVtxCorrNum"));

    fHistTrackTrigCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigCorrDen")); 
    fHistTrackTrigVtxCorrNumNSD = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigVtxCorrNumNSD"));
    fHistTrackTrigNSD = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigNSD"));


    fHistTrigVtxCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigVtxCorrNum"));
    fHistTrigVtxCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigVtxCorrDen"));

    fHistTrigCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigCorrDen"));

    fHistTrigVtxCorrNumNSD = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigVtxCorrNumNSD"));
    fHistEvTrigNSD = dynamic_cast<TH2F*> (fOutput->FindObject("fHistEvTrigNSD"));

    fHistMCEtavsZTriggMCvtxEvts = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZTriggMCvtxEvts"));
    fHistMCEtavsZTriggESDvtxEvts = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZTriggESDvtxEvts"));
    fHistMCEtavsZ = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZ"));

    fHistMCEtaInel = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMCEtaInel"));
    fHistMCEtaNonDiffractive = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMCEtaNonDiffractive"));
    fHistMCEtaNonSingleDiffractive = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMCEtaNonSingleDiffractive"));

    fHistoProcessType = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoProcessType")); 
    fHistoProcessTypeTriggered = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoProcessTypeTriggered"));

    fHistContributorsvsMCVtx = dynamic_cast<TH2F*> (fOutput->FindObject("fHistContributorsvsMCVtx"));

    fHistoDetectableNotr = dynamic_cast<TH3F*> (fOutput->FindObject("fHistoDetectableNotr"));
    fHistoDetectabletr = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectabletr"));

    fHistoNonStoppingTracks = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoNonStoppingTracks"));

    fHistoDetectedLay1 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectedLay1"));
    fHistoDetectedLay2 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectedLay2"));

    fHistoPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoPt"));

    fHistoRTRm1 = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoRTRm1"));
    fHistMCvtx = dynamic_cast<TH3F*> (fOutput->FindObject("fHistMCvtx"));  
   
    fHistMultAllNonDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultAllNonDiff"));
    fHistMultAllSingleDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultAllSingleDiff"));
    fHistMultAllDoubleDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultAllDoubleDiff"));
    fHistMultTrVtxNonDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultTrVtxNonDiff"));
    fHistMultTrVtxSingleDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultTrVtxSingleDiff"));
    fHistMultTrVtxDoubleDiff = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMultTrVtxDoubleDiff"));

    fHistMCEtaNonSingleDiffractiveLargeBin = dynamic_cast<TH1F*> (fOutput->FindObject("fHistMCEtaNonSingleDiffractiveLargeBin"));
  }
}
