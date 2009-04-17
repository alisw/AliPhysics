#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include <TClonesArray.h>

#include "TH1F.h"
#include "TH2F.h" 
#include "TH3F.h" 
#include "TCanvas.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>

#include <AliMultiplicity.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>

#include "AliAnalysisTaskRL.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include <AliStack.h>
#include <AliLog.h>
#include "AliTrackReference.h"

#include <AliGenEventHeader.h> 
#include "AliAnalysisTaskSPDdNdEta.h"


ClassImp(AliAnalysisTaskSPDdNdEta)

//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::AliAnalysisTaskSPDdNdEta(const char *name) 
  : AliAnalysisTask(name, "SPDdNdEtaTask"), 

  fESD(0), 
  fOutput(0), 

  fCorr(kFALSE), 
  fTrigger(0),

  //Data to be corrected
  fHistSPDRAWMultvsZ(0),
  fHistSPDRAWMultvsZTriggEvts(0),
  fHistSPDRAWEtavsZ(0),

  //Clusters inner layer and tracklets
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

  //SPD vertex 
  fHistSPDvtx(0),                 
  fHistSPDvtx3D(0),
  fHistSPDvtxZ(0),
  fHistNcontribSPDvtxvsSPDvtx(0),
  fHistNcontribSPDvtx3D(0),
  fHistNcontribSPDvtxZ(0),
  fHistNcontribSPDvtxall(0),
  fHistSPDmultvsSPDvtx(0),

  //SPD fired chips  
  fHistSPDcl1multvsnFiredChipsLay1(0),
  fHistSPDmultvsnFiredChipsLay1(0),
  fHistSPDmultvsnFiredChipsLay2(0),
  fHistnFiredChipsLay2vsnFiredChipsLay1(0),
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec(0),

  //Track level correction histograms
  fHistBkgCorrNum(0),   
  fHistBkgCorrDen(0),

  fHistAlgEffNum(0),

  fHistNonDetectableCorrNum(0),  
  fHistNonDetectableCorrDen(0),

  fHistTrackTrigVtxCorrNum(0),   
  fHistTrackTrigCorrDen(0),      

  //Event level correction histograms
  fHistTrigVtxCorrNum(0),           
  fHistTrigVtxCorrDen(0),           
  
  fHistTrigCorrDen(0),              

  // MC distribustions
  fHistMCEtavsZTriggMCvtxEvts(0),
  fHistMCEtavsZTriggESDvtxEvts(0),
  fHistMCEtavsZ(0),                

  // Check histos
  fHistTRradius(0),
  fHistContributorsvsDeVtx(0),
  fHistoDetectableNotr(0),
  fHistoDetectabletr(0),
  fHistoNonStoppingTracks(0),
  fHistoDetectedLay1(0),
  fHistoDetectedLay2(0),
  fHistoPt(0),
  fHistoDetectableTRm1(0),
  fHistoDetectableTR0(0),
  fHistoDetectableTR1(0),
  fHistoDetectableTR2(0),
  fHistoDetectableTR3(0),
  fHistoDetectableTR4(0),
  fHistoDetectableTR5(0),
  fHistoDetectableTR6(0),   
  fHistoRTRm1(0)
  
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
  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);  

}

//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::CreateOutputObjects()
{

  // Create histograms
  fOutput = new TList;
  fOutput->SetOwner();

  fHistSPDRAWMultvsZ= new TH2F("fHistSPDRAWMultvsZ", "",200,0.,200.,20,-20.,20.);
  fHistSPDRAWMultvsZ->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZ);

  fHistSPDRAWMultvsZTriggEvts = new TH2F("fHistSPDRAWMultvsZTriggEvts", "",200,0.,200.,20,-20.,20.);
  fHistSPDRAWMultvsZTriggEvts->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZTriggEvts);

  fHistSPDRAWEtavsZ = new TH2F("fHistSPDRAWEtavsZ", "Tracklet pseudorapidity distribution", 120, -3.,3.,40,-20.,20.);
  fHistSPDRAWEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDRAWEtavsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWEtavsZ->Sumw2();
  fOutput->Add(fHistSPDRAWEtavsZ);

  
  fHistSPDmultEtacut = new TH1F("fHistSPDmultEtacut", "Tracklet multiplicity distribution", 200, 0., 200);
  fHistSPDmultEtacut->GetXaxis()->SetTitle("Reconstructed multiplicity (|#eta|<1.5)");
  fHistSPDmultEtacut->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultEtacut);

  fHistSPDmult = new TH1F("fHistSPDmult", "Tracklet multiplicity distribution", 200, 0., 200);
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

  fHistSPDcl1multEtacutLay1 = new TH1F("fHistSPDcl1multEtacutLay1", "Cluster multiplicity (inner layer)", 200, 0., 200);
  fHistSPDcl1multEtacutLay1->GetXaxis()->SetTitle("Cluster multiplicity lay1 (|#eta|<2.)");
  fHistSPDcl1multEtacutLay1->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1multEtacutLay1);

  fHistSPDcl1mult = new TH1F("fHistSPDcl1mult", "Cluster multiplicity (inner layer)", 200, 0., 200);
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

  fHistSPDtheta = new TH1F("fHistSPDtheta", "Tracklet #theta  distribution", 360, 0.,2*TMath::Pi());
  fHistSPDtheta->GetXaxis()->SetTitle("#theta [rad]");
  fHistSPDtheta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDtheta);

  fHistSPDcl1theta = new TH1F("fHistSPDcl1theta", "Cluster #theta (inner layer)", 360, 0.,2*TMath::Pi());
  fHistSPDcl1theta->GetXaxis()->SetTitle("#theta [rad]");
  fHistSPDcl1theta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDcl1theta);

  fHistSPDdePhi= new TH1F("fHistSPDdePhi", "Tracklet #Delta#varphi  distribution",400,-0.1,.1);
  fHistSPDdePhi->GetXaxis()->SetTitle("#Delta#varphi [rad]");
  fHistSPDdePhi->GetYaxis()->SetTitle("z_{MCvtx}");
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


  fHistSPDvtx = new TH1F("fHistSPDvtx", "SPD vertex  distribution - all events",20,-20.,20.);
  fHistSPDvtx->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtx->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtx);

  fHistSPDvtx3D = new TH3F("fHistSPDvtx3D", "SPD vertex distribution",100,-5.,5.,100,-5.,5.,400,-20.,20.);
  fOutput->Add(fHistSPDvtx3D);

  fHistSPDvtxZ = new TH3F("fHistSPDvtxZ", "SPD vertex distribution",100,-5.,5.,100,-5.,5.,400,-20.,20.);
  fOutput->Add(fHistSPDvtxZ);

  fHistNcontribSPDvtxvsSPDvtx= new TH2F("fHistNcontribSPDvtxvsSPDvtx", " ",100,-50.,50.,202,-2.,200.);
  fHistNcontribSPDvtxvsSPDvtx->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistNcontribSPDvtxvsSPDvtx->GetYaxis()->SetTitle("# contributors");
  fOutput->Add(fHistNcontribSPDvtxvsSPDvtx);

  fHistNcontribSPDvtx3D= new TH1F("fHistNcontribSPDvtx_3D", "SPD vtx 3D",202,-2.,200.);
  fHistNcontribSPDvtx3D->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtx3D->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtx3D);

  fHistNcontribSPDvtxZ= new TH1F("fHistNcontribSPDvtxZ", "SPD vtx Z",202,-2.,200.);
  fHistNcontribSPDvtxZ->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtxZ->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtxZ);

  fHistNcontribSPDvtxall= new TH1F("fHistNcontribSPDvtxall", "SPD vtx - all events",202,-2.,200.);
  fHistNcontribSPDvtxall->GetXaxis()->SetTitle("# contributors");
  fHistNcontribSPDvtxall->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistNcontribSPDvtxall);

  fHistSPDmultvsSPDvtx= new TH2F("fHistSPDmultvsSPDvtx", "",20,-20.,20.,200,0.,200.);
  fHistSPDmultvsSPDvtx->GetXaxis()->SetTitle("z_{recvtx} [cm]");
  fHistSPDmultvsSPDvtx->GetYaxis()->SetTitle("Reconstructed multiplicity");
  fOutput->Add(fHistSPDmultvsSPDvtx);

  fHistSPDcl1multvsnFiredChipsLay1 = new TH2F("fHistSPDcl1multvsnFiredChipsLay1", "",200,0.,200.,200,0.,200.);
  fHistSPDcl1multvsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistSPDcl1multvsnFiredChipsLay1->GetYaxis()->SetTitle("Cluster lay1 multiplicity");
  fOutput->Add(fHistSPDcl1multvsnFiredChipsLay1);

  fHistSPDmultvsnFiredChipsLay1 = new TH2F("fHistSPDmultvsnFiredChipsLay1","",200,0.,200.,200,0.,200.);
  fHistSPDmultvsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistSPDmultvsnFiredChipsLay1->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvsnFiredChipsLay1);

  fHistSPDmultvsnFiredChipsLay2 = new TH2F("fHistSPDmultvsnFiredChipsLay2","",200,0.,200.,200,0.,200.);
  fHistSPDmultvsnFiredChipsLay2->GetXaxis()->SetTitle("# fired chips lay2");
  fHistSPDmultvsnFiredChipsLay2->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvsnFiredChipsLay2);

  fHistnFiredChipsLay2vsnFiredChipsLay1 = new TH2F("fHistnFiredChipsLay2vsnFiredChipsLay1","",200,0.,200.,200,0.,200.);
  fHistnFiredChipsLay2vsnFiredChipsLay1->GetXaxis()->SetTitle("# fired chips lay1");
  fHistnFiredChipsLay2vsnFiredChipsLay1->GetYaxis()->SetTitle("# fired chip lay2");
  fOutput->Add(fHistnFiredChipsLay2vsnFiredChipsLay1);

  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec = new TH2F("fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec","",200,0.,200.,200,0.,200.);
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec->GetXaxis()->SetTitle("# fired chips lay1");
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec->GetYaxis()->SetTitle("# fired chip lay2");
  fOutput->Add(fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec);


  if (fCorr) {
    fHistBkgCorrNum = new TH2F("fHistBkgCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fHistBkgCorrNum->Sumw2();
    fOutput->Add(fHistBkgCorrNum);

    fHistBkgCorrDen = new TH2F("fHistBkgCorrDen","",60,-3.00,3.00,20,-20.,20.);
    fHistBkgCorrDen->Sumw2();
    fOutput->Add(fHistBkgCorrDen);

    fHistAlgEffNum = new TH2F("fHistAlgEffNum","",60,-3.00,3.00,20,-20.,20.);
    fHistAlgEffNum->Sumw2();  
    fOutput->Add(fHistAlgEffNum);  

    fHistNonDetectableCorrNum = new TH2F("fHistNonDetectableCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fHistNonDetectableCorrNum->Sumw2();
    fOutput->Add(fHistNonDetectableCorrNum);

    fHistNonDetectableCorrDen = new TH2F("fHistNonDetectableCorrDen","",60,-3.00,3.00,20,-20.,20.);
    fHistNonDetectableCorrDen->Sumw2();
    fOutput->Add(fHistNonDetectableCorrDen);

    fHistTrackTrigVtxCorrNum = new TH2F("fHistTrackTrigVtxCorrNum","",60,-3.00,3.00,20,-20.,20.);
    fHistTrackTrigVtxCorrNum->Sumw2();
    fOutput->Add(fHistTrackTrigVtxCorrNum);

    fHistTrackTrigCorrDen = new TH2F("fHistTrackTrigCorrDen","",60,-3.00,3.00,20,-20.,20.);
    fHistTrackTrigCorrDen->Sumw2();
    fOutput->Add(fHistTrackTrigCorrDen);

    //Event level correction histograms
    fHistTrigVtxCorrNum = new TH2F("fHistTrigVtxCorrNum","",200,0.,200.,20,-20.,20.);
    fHistTrigVtxCorrNum->Sumw2();
    fOutput->Add(fHistTrigVtxCorrNum);

    fHistTrigVtxCorrDen = new TH2F("fHistTrigVtxCorrDen","",200,0.,200.,20,-20.,20.);
    fHistTrigVtxCorrDen->Sumw2();
    fOutput->Add(fHistTrigVtxCorrDen);

    fHistTrigCorrDen = new TH2F("fHistTrigCorrDen","",200,0.,200.,20,-20.,20.);
    fHistTrigCorrDen->Sumw2();
    fOutput->Add(fHistTrigCorrDen);

    // MC distributions
    fHistMCEtavsZTriggMCvtxEvts = new TH2F("fHistMCEtavsZTriggMCvtxEvts","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZTriggMCvtxEvts->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZTriggMCvtxEvts->GetYaxis()->SetTitle("Entries");
    fHistMCEtavsZTriggMCvtxEvts->Sumw2();
    fOutput->Add(fHistMCEtavsZTriggMCvtxEvts);

    fHistMCEtavsZTriggESDvtxEvts = new TH2F("fHistMCEtavsZTriggESDvtxEvts","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZTriggESDvtxEvts->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZTriggESDvtxEvts->GetYaxis()->SetTitle("Entries");
    fHistMCEtavsZTriggESDvtxEvts->Sumw2();
    fOutput->Add(fHistMCEtavsZTriggESDvtxEvts);

    fHistMCEtavsZ = new TH2F("fHistMCEtavsZ","Generated pseudorapidity distribution",60,-3.00,3.00,20,-20.,20.);
    fHistMCEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
    fHistMCEtavsZ->GetYaxis()->SetTitle("Entries");
    fHistMCEtavsZ->Sumw2();
    fOutput->Add(fHistMCEtavsZ);
   
    fHistTRradius = new TH1F("fHistTRradius","ITS track reference rad",200,0.,10.);
    fOutput->Add(fHistTRradius); 
 
    fHistContributorsvsDeVtx = new TH2F("fHistContributorsvsDeVtx","",200,-20.,20.,22,-2.,20.);
    fOutput->Add(fHistContributorsvsDeVtx);
 
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

    fHistoDetectableTRm1 = new TH2F("fHistoDetectableTRm1","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTRm1);

    fHistoDetectableTR0 = new TH2F("fHistoDetectableTR0","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR0);
 
    fHistoDetectableTR1 = new TH2F("fHistoDetectableTR1","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR1);

    fHistoDetectableTR2 = new TH2F("fHistoDetectableTR2","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR2);
 
    fHistoDetectableTR3 = new TH2F("fHistoDetectableTR3","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR3);

    fHistoDetectableTR4 = new TH2F("fHistoDetectableTR4","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR4);

    fHistoDetectableTR5 = new TH2F("fHistoDetectableTR5","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR5);

    fHistoDetectableTR6 = new TH2F("fHistoDetectableTR6","",60,-3.00,3.00,20,-20.,20.);
    fOutput->Add(fHistoDetectableTR6);
 
    fHistoRTRm1 = new TH1F("fHistoRTRm1","",10000,0.,5000);
    fOutput->Add(fHistoRTRm1);

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
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  Int_t nContrib = vtxESD->GetNContributors();

  const AliMultiplicity* multESD = fESD->GetMultiplicity();

  // Loading tracklets...
  Int_t multSPD = 0;
  multSPD = multESD->GetNumberOfTracklets();

  // Trigger
  ULong64_t triggerMask;
  ULong64_t spdFO   = (1 << 14);
  ULong64_t v0left  = (1 << 11);
  ULong64_t v0right = (1 << 12);

  triggerMask=fESD->GetTriggerMask();

  Bool_t eventTriggered = kFALSE;
  // No trigger
  if (fTrigger==0) eventTriggered = kTRUE;
  // MB1: GFO || V0OR
  if (fTrigger==1) eventTriggered = triggerMask&spdFO || ((triggerMask&v0left) || (triggerMask&v0right));
  // MB2: GFO && V0OR
  if (fTrigger==2) eventTriggered = triggerMask&spdFO && ((triggerMask&v0left) || (triggerMask&v0right));

  PostData(0, fOutput);

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

//  Printf("There are %d tracklets in this event", multSPD);
  if (esdvtx[2]!=0.&&eventTriggered&&multSPD!=0) {

    fHistSPDvtx->Fill(esdvtx[2]);
    if (strcmp(vtxESD->GetTitle(),"vertexer: 3D") == 0) fHistSPDvtx3D->Fill(esdvtx[0],esdvtx[1],esdvtx[2]);
    else fHistSPDvtxZ->Fill(esdvtx[0],esdvtx[1],esdvtx[2]);

    for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
      Float_t thetaTr= multESD->GetTheta(itracklet);
      Float_t phiTr= multESD->GetPhi(itracklet);
      Float_t dePhiTr= multESD->GetDeltaPhi(itracklet);
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
      if (strcmp(vtxESD->GetTitle(),"vertexer: Z") == 0) fHistSPDdePhiZ->Fill(dePhiTr);
      if (strcmp(vtxESD->GetTitle(),"vertexer: 3D") == 0) fHistSPDdePhi3D->Fill(dePhiTr);
      fHistSPDphivsSPDeta->Fill(recEtaSPD,phiTr);
      fHistSPDcl1phivsSPDcl1eta->Fill(recEtaSPD,phiTr);

      // calculate multiplicity in etacut
      if (TMath::Abs(recEtaSPD)<1.5) multSPDEtacut++;
      if (TMath::Abs(recEtaSPD)<2.)  multSPDcl1EtacutLay1++;
    }

    for (Int_t iCl1=0; iCl1<nSingleCl1; ++iCl1) {
      Float_t thetaSingleCl1 = multESD->GetThetaSingle(iCl1);
      //calculate eta
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
    fHistSPDmultvsSPDvtx->Fill(esdvtx[2],multSPD);
    fHistNcontribSPDvtxvsSPDvtx->Fill(esdvtx[2],nContrib);
    fHistSPDRAWMultvsZ->Fill(multSPD,esdvtx[2]);
    fHistSPDmultvsnFiredChipsLay1->Fill(nFiredChipsLay1,multSPD);
    fHistSPDmultvsnFiredChipsLay2->Fill(nFiredChipsLay2,multSPD);

  } //ev selection

  if (eventTriggered) fHistSPDRAWMultvsZTriggEvts->Fill(multSPD,esdvtx[2]);

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
    //MC vertex
    TArrayF vtxMC(3);
    genHeader->PrimaryVertex(vtxMC);


    //Tracks from MC
    Int_t  multMCCharged = 0;
    Int_t  nMCPart = stack->GetNprimary();
    Float_t* eta_gen = new Float_t[nMCPart];  
    Float_t* pt_gen = new Float_t[nMCPart];
    Int_t* stackIndexOfPrimaryParts = new Int_t[nMCPart];
    Int_t* reconstructedPrimaryPart = new Int_t[nMCPart];
    Int_t* detectedPrimaryPartLay1 = new Int_t[nMCPart];
    Int_t* detectedPrimaryPartLay2 = new Int_t[nMCPart];
    Int_t* detectablePrimaryPart = new Int_t[nMCPart];
    Bool_t* stoppingPartAtLay2 = new Bool_t[nMCPart];  

    // Loading track references...
    Float_t rminL1 = 3.4;
    Float_t rmaxL1 = 4.4;
    Float_t rminL2 = 6.9;
    Float_t rmaxL2 = 7.9;

    Int_t* nTrackRefPerPrim = new Int_t[nMCPart];
    Float_t trackRefPerPrim[nMCPart][500][2];

    TTree* tRefTree = eventHandler->TreeTR();

    AliTrackReference *tref=0x0;
    mcEvent->ConnectTreeTR(tRefTree);

    // Loop over MC particles
    for (Int_t imc=0; imc<nMCPart; imc++) {
      TParticle* part = stack->Particle(imc);
      Bool_t IsPrimary = stack->IsPhysicalPrimary(imc);
      if (!IsPrimary)                        continue;
      TParticlePDG* pdgPart = part->GetPDG();
      if (TMath::Abs(pdgPart->Charge())!=3)  continue;
      Float_t theta = part->Theta();
      if (theta==0 || theta==TMath::Pi())    continue;
      Float_t eta = part->Eta();
      Float_t pt = part->Pt();
      eta_gen[multMCCharged] = eta; 
      pt_gen[multMCCharged] = pt;
      stackIndexOfPrimaryParts[multMCCharged] = imc;

      reconstructedPrimaryPart[multMCCharged]=kFALSE;
      detectedPrimaryPartLay1[multMCCharged]=kFALSE;
      detectedPrimaryPartLay2[multMCCharged]=kFALSE;
      detectablePrimaryPart[multMCCharged]=kFALSE;
      stoppingPartAtLay2[multMCCharged]=kFALSE;

      fHistoPt->Fill(pt_gen[multMCCharged]);

      AliMCParticle* mcpart = mcEvent->GetTrack(imc);
      Int_t nref = mcpart->GetNumberOfTrackReferences();
      nTrackRefPerPrim[multMCCharged]=nref;

      Bool_t nonStoppingPart=kFALSE;

      // Loop over all the track refs of the track
      for (Int_t iref=0;iref<nref;iref++) { 
        tref = mcpart->GetTrackReference(iref);
        if (tref) {
          trackRefPerPrim[multMCCharged][iref][0]=tref->R();
          trackRefPerPrim[multMCCharged][iref][1]=(Float_t)tref->DetectorId();
          if (trackRefPerPrim[multMCCharged][iref][1]==0) fHistTRradius->Fill(trackRefPerPrim[multMCCharged][iref][0]); // ITS only 
//        Printf("det ID: %f R: %f", trackRefPerPrim[multMCCharged][iref][1], trackRefPerPrim[multMCCharged][iref][0]);
          if (iref==nref-1&&trackRefPerPrim[multMCCharged][iref][1]!=-1) {
            nonStoppingPart=kTRUE;
          }
        }
      }
      if (nref==0||nonStoppingPart) {
        detectablePrimaryPart[multMCCharged]=kTRUE; //0. if prim doesn't disappear it is detectable
        if (nonStoppingPart) fHistoNonStoppingTracks->Fill(eta_gen[multMCCharged],vtxMC[2]);   
        if (nTrackRefPerPrim[multMCCharged]==0) fHistoDetectableNotr->Fill(eta_gen[multMCCharged],vtxMC[2],pt_gen[multMCCharged]); 
      }

      multMCCharged++;
    } //end of MC particle loop

    // Check histos
    if (multSPD==0 && eventTriggered && esdvtx[2]!=0.) {
      fHistContributorsvsDeVtx->Fill(esdvtx[2]-vtxMC[2],nContrib);
    }

    //Event selection
    if (eventTriggered && esdvtx[2]!=0. && multSPD!=0) {

      for (Int_t imc=0; imc<multMCCharged; imc++) {

        // Find out detected prims on each layer
        for (Int_t iref=0; iref<nTrackRefPerPrim[imc]; iref++) {
          if (trackRefPerPrim[imc][iref][0]>rminL1&&trackRefPerPrim[imc][iref][0]<rmaxL1&&trackRefPerPrim[imc][iref][1]==0) {
            detectedPrimaryPartLay1[imc] = kTRUE;
            fHistoDetectedLay1->Fill(eta_gen[imc],vtxMC[2]);
            break; 
          }
        }
        for (Int_t iref=0; iref<nTrackRefPerPrim[imc]; iref++) {
          if (trackRefPerPrim[imc][iref][0]>rminL2&&trackRefPerPrim[imc][iref][0]<rmaxL2) {
            if (trackRefPerPrim[imc][iref][1]==0) {
              detectedPrimaryPartLay2[imc] = kTRUE;
              fHistoDetectedLay2->Fill(eta_gen[imc],vtxMC[2]);
              detectablePrimaryPart[imc] = kTRUE;  //1. if prim is detected it is also detectable
              break;
            } 
            if (trackRefPerPrim[imc][iref][1]==-1) {
              stoppingPartAtLay2[imc]=kTRUE;
            }
          }
        }
   
        // Find out SPD detectable prims 
        if (!detectablePrimaryPart[imc]) {
          for (Int_t iref=0; iref<nTrackRefPerPrim[imc]; iref++) {
            if (trackRefPerPrim[imc][iref][0]>rmaxL2) { 
              detectablePrimaryPart[imc] = kTRUE; //2. if it has crossed layer 2 it is detectable
              fHistoDetectabletr->Fill(eta_gen[imc],vtxMC[2]);
              if (trackRefPerPrim[imc][iref][1]==-1) {
                fHistoDetectableTRm1->Fill(eta_gen[imc],vtxMC[2]);
                fHistoRTRm1->Fill(trackRefPerPrim[imc][iref][0]);
              }
              else if (trackRefPerPrim[imc][iref][1]==0) fHistoDetectableTR0->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==1) fHistoDetectableTR1->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==2) fHistoDetectableTR2->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==3) fHistoDetectableTR3->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==4) fHistoDetectableTR4->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==5) fHistoDetectableTR5->Fill(eta_gen[imc],vtxMC[2]);
              else if (trackRefPerPrim[imc][iref][1]==6) fHistoDetectableTR6->Fill(eta_gen[imc],vtxMC[2]);
              break;
            }
          }
        }
//        Printf("Primary nr %d", stackIndexOfPrimaryParts[imc]);  
      }
 
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
          fHistBkgCorrNum->Fill(eta_gen[imc],vtxMC[2]);
        }
        if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) fHistAlgEffNum->Fill(eta_gen[imc],vtxMC[2]);

        fHistNonDetectableCorrNum->Fill(eta_gen[imc],vtxMC[2]);

        if (detectablePrimaryPart[imc]) fHistNonDetectableCorrDen->Fill(eta_gen[imc],vtxMC[2]);

        fHistMCEtavsZTriggESDvtxEvts->Fill(eta_gen[imc],esdvtx[2]);
        fHistMCEtavsZTriggMCvtxEvts->Fill(eta_gen[imc],vtxMC[2]);
      }

      fHistTrigVtxCorrDen->Fill(multSPD,vtxMC[2]);

    } //ev selection
 
    if (eventTriggered) {
      fHistTrigCorrDen->Fill(multSPD,vtxMC[2]);
      for (Int_t imc=0; imc<multMCCharged; imc++) {
        fHistTrackTrigCorrDen->Fill(eta_gen[imc],vtxMC[2]); 
      }
    }

    for (Int_t imc=0; imc<multMCCharged; imc++) {
      fHistTrackTrigVtxCorrNum->Fill(eta_gen[imc],vtxMC[2]);
      fHistMCEtavsZ->Fill(eta_gen[imc],vtxMC[2]);
    }

    fHistTrigVtxCorrNum->Fill(multSPD,vtxMC[2]);

    delete[] eta_gen;
    delete[] pt_gen;
    delete[] nTrackRefPerPrim;
    delete[] stackIndexOfPrimaryParts;
    delete[] reconstructedPrimaryPart;
    delete[] detectedPrimaryPartLay1;
    delete[] detectedPrimaryPartLay2;
    delete[] detectablePrimaryPart;
    delete[] stoppingPartAtLay2;
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

  fHistSPDvtx = dynamic_cast<TH1F*> (fOutput->FindObject("fHistSPDvtx"));
  fHistSPDvtx3D = dynamic_cast<TH3F*> (fOutput->FindObject("fHistSPDvtx3D"));
  fHistSPDvtxZ = dynamic_cast<TH3F*> (fOutput->FindObject("fHistSPDvtxZ"));
  fHistNcontribSPDvtxvsSPDvtx = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNcontribSPDvtxvsSPDvtx"));
  fHistNcontribSPDvtx3D = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtx3D"));
  fHistNcontribSPDvtxZ = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtxZ"));
  fHistNcontribSPDvtxall = dynamic_cast<TH1F*> (fOutput->FindObject("fHistNcontribSPDvtxall"));
  fHistSPDmultvsSPDvtx= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDmultvsSPDvtx"));

  fHistSPDcl1multvsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDcl1multvsnFiredChipsLay1"));
  fHistSPDmultvsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDmultvsnFiredChipsLay1"));
  fHistSPDmultvsnFiredChipsLay2= dynamic_cast<TH2F*> (fOutput->FindObject("fHistSPDmultvsnFiredChipsLay2"));
  fHistnFiredChipsLay2vsnFiredChipsLay1= dynamic_cast<TH2F*> (fOutput->FindObject("fHistnFiredChipsLay2vsnFiredChipsLay1"));
  fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec= dynamic_cast<TH2F*> (fOutput->FindObject("fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec"));

  if (fCorr) {
    fHistBkgCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistBkgCorrNum"));
    fHistBkgCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistBkgCorrDen"));

    fHistAlgEffNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistAlgEffNum"));

    fHistNonDetectableCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNonDetectableCorrNum"));
    fHistNonDetectableCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistNonDetectableCorrDen"));

    fHistTrackTrigVtxCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigVtxCorrNum"));

    fHistTrackTrigCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrackTrigCorrDen")); 

    //Event level
    fHistTrigVtxCorrNum = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigVtxCorrNum"));
    fHistTrigVtxCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigVtxCorrDen"));

    fHistTrigCorrDen = dynamic_cast<TH2F*> (fOutput->FindObject("fHistTrigCorrDen"));

    //MC distrib
    fHistMCEtavsZTriggMCvtxEvts = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZTriggMCvtxEvts"));
    fHistMCEtavsZTriggESDvtxEvts = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZTriggESDvtxEvts"));
    fHistMCEtavsZ = dynamic_cast<TH2F*> (fOutput->FindObject("fHistMCEtavsZ"));

    //Check histos  
    fHistTRradius = dynamic_cast<TH1F*> (fOutput->FindObject("fHistTRradius"));

    fHistContributorsvsDeVtx = dynamic_cast<TH2F*> (fOutput->FindObject("fHistContributorsvsDeVtx"));

    fHistoDetectableNotr = dynamic_cast<TH3F*> (fOutput->FindObject("fHistoDetectableNotr"));
    fHistoDetectabletr = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectabletr"));

    fHistoNonStoppingTracks = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoNonStoppingTracks"));

    fHistoDetectedLay1 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectedLay1"));
    fHistoDetectedLay2 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectedLay2"));

    fHistoPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoPt"));

    fHistoDetectableTRm1 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTRm1"));
    fHistoDetectableTR0 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR0"));
    fHistoDetectableTR1 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR1"));
    fHistoDetectableTR2 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR2"));
    fHistoDetectableTR3 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR3"));
    fHistoDetectableTR4 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR4"));
    fHistoDetectableTR5 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR5"));
    fHistoDetectableTR6 = dynamic_cast<TH2F*> (fOutput->FindObject("fHistoDetectableTR6"));
    fHistoRTRm1 = dynamic_cast<TH1F*> (fOutput->FindObject("fHistoRTRm1"));
  }
}
