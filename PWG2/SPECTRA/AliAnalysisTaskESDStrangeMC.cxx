
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

// macro to study V0s, with Monte Carlo information access
// loops over ESD files, and creates AliAODv0
// Author: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliLog.h"

#include "AliAnalysisTaskESDStrangeMC.h"


ClassImp(AliAnalysisTaskESDStrangeMC)

//________________________________________________________________________
AliAnalysisTaskESDStrangeMC::AliAnalysisTaskESDStrangeMC(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fListHist(), 
    fHistPtMC(0),
    fHistMCMultiplicity(0),
    fHistMCPtVsYK0s(0),
    fHistMCPtVsYLambda(0),
    fHistMCPtVsYAntiLambda(0),
    fHistTrackPerEvent(0),
    fHistMCDaughterTrack(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),
    fHistDcaPosToPrimVertex(0),
    fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0),
    fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0),
    fHistDecayLengthV0(0),
    fHistDcaV0Daughters(0),
    fHistChi2(0),
    fHistCosPointAngle(0),
    fHistCosPointAngleZoom(0),
    fHistPtVsYK0s(0),
    fHistPtVsYK0sMI(0),
    fHistPtVsYLambda(0),
    fHistPtVsYLambdaMI(0),
    fHistPtVsYAntiLambda(0),
    fHistPtVsYAntiLambdaMI(0),
    fHistMassK0(0),
    fHistMassK0MI(0),
    fHistMassLambda(0),
    fHistMassLambdaMI(0),
    fHistMassAntiLambda(0),
    fHistMassAntiLambdaMI(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusK0MI(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusLambdaMI(0),
    fHistMassVsRadiusAntiLambda(0),
    fHistMassVsRadiusAntiLambdaMI(0),
    fHistArmenterosPodolanski(0),
    fHistArmenterosPodolanskiMI(0),
    fHistAsMcPtK0(0),
    fHistAsMcPtK0MI(0),
    fHistAsMcPtLambda(0),
    fHistAsMcPtLambdaMI(0),
    fHistAsMcPtAntiLambda(0),
    fHistAsMcPtAntiLambdaMI(0),
    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomK0MI(0),
    fHistAsMcPtZoomLambda(0),
    fHistAsMcPtZoomLambdaMI(0),
    fHistPidMcMassK0(0),
    fHistPidMcMassK0MI(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassLambdaMI(0),
    fHistPidMcMassAntiLambda(0),
    fHistPidMcMassAntiLambdaMI(0),
    fHistAsMcMassK0(0),
    fHistAsMcMassK0MI(0),
    fHistAsMcMassLambda(0),
    fHistAsMcMassLambdaMI(0),
    fHistAsMcMassAntiLambda(0),
    fHistAsMcMassAntiLambdaMI(0),
    fHistAsMcMassVsRadiusK0(0),
    fHistAsMcMassVsRadiusK0MI(0),
    fHistAsMcMassVsRadiusLambda(0),
    fHistAsMcMassVsRadiusLambdaMI(0),
    fHistAsMcMassVsRadiusAntiLambda(0),
    fHistAsMcMassVsRadiusAntiLambdaMI(0),
    fHistAsMcResxK0(0),
    fHistAsMcResyK0(0),
    fHistAsMcReszK0(0),
    fHistAsMcResrVsRadiusK0(0),
    fHistAsMcReszVsRadiusK0(0),
    fHistAsMcResxK0MI(0),
    fHistAsMcResyK0MI(0),
    fHistAsMcReszK0MI(0),
    fHistAsMcResrVsRadiusK0MI(0),
    fHistAsMcReszVsRadiusK0MI(0),
    fHistAsMcResxLambda(0),
    fHistAsMcResyLambda(0),
    fHistAsMcReszLambda(0),
    fHistAsMcResrVsRadiusLambda(0),
    fHistAsMcReszVsRadiusLambda(0),
    fHistAsMcResxLambdaMI(0),
    fHistAsMcResyLambdaMI(0),
    fHistAsMcReszLambdaMI(0),
    fHistAsMcResrVsRadiusLambdaMI(0),
    fHistAsMcReszVsRadiusLambdaMI(0),
    fHistAsMcResxAntiLambda(0),
    fHistAsMcResyAntiLambda(0),
    fHistAsMcReszAntiLambda(0),
    fHistAsMcResrVsRadiusAntiLambda(0),
    fHistAsMcReszVsRadiusAntiLambda(0),
    fHistAsMcResxAntiLambdaMI(0),
    fHistAsMcResyAntiLambdaMI(0),
    fHistAsMcReszAntiLambdaMI(0),
    fHistAsMcResrVsRadiusAntiLambdaMI(0),
    fHistAsMcReszVsRadiusAntiLambdaMI(0)
    
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  //DefineOutput(0, TH1F::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskESDStrangeMC::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    tree->SetBranchStatus("fV0s.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskESDStrangeMC::CreateOutputObjects() 
{
  // Create histograms
  // Called once

  fListHist = new TList();


  //***************
  // MC histograms
  //***************
  fHistPtMC = new TH1F("h1PtMC", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPtMC->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtMC->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtMC->SetMarkerStyle(kFullCircle);
  fListHist->Add(fHistPtMC);

  // Multiplicity
  fHistMCMultiplicity           = new TH1F("h1MCMultiplicity", "MC Multiplicity;Ntracks;Count", 201, -0.5, 200.5);
  fListHist->Add(fHistMCMultiplicity);

  // Pt and rapidity distribution:
  fHistMCPtVsYK0s               = new TH2F("h2MCPtVsYK0s", "K^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYK0s);

  fHistMCPtVsYLambda            = new TH2F("h2MCPtVsYLambda", "#Lambda^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYLambda);

  fHistMCPtVsYAntiLambda        = new TH2F("h2MCPtVsYAntiLambda", "#bar{#Lambda}^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYAntiLambda);
 

  //***********************************
  // Reconstructed particles histograms
  //***********************************

  // multiplicity
  fHistTrackPerEvent           = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",50,0,50);
  fListHist->Add(fHistTrackPerEvent);

  fHistMCDaughterTrack         = new TH1F("h1MCDaughterTrack","Distribution of mc id for daughters;id tags;Counts",15,0,15);
  fListHist->Add(fHistMCDaughterTrack);

  // Primary Vertex:
  fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",40,-1,1);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",40,-1,1);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",60,-5,5);
  fListHist->Add(fHistPrimaryVertexZ);

  // Cut checks:
  fHistDcaPosToPrimVertex      = new TH2F("h2DcaPosToPrimVertex", "Positive V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertex);

  fHistDcaNegToPrimVertex      = new TH2F("h2DcaNegToPrimVertex", "Negative V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertex);

  fHistDcaPosToPrimVertexZoom  = new TH2F("h2DcaPosToPrimVertexZoom", "Positive V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexZoom);

  fHistDcaNegToPrimVertexZoom  = new TH2F("h2DcaNegToPrimVertexZoom", "Negative V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexZoom);

  fHistRadiusV0                = new TH2F("h2RadiusV0", "Radius;Radius(cm);Status",1000,0,100,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0);

  fHistDecayLengthV0           = new TH2F("h2DecayLengthV0", "V0s decay Length;decay length(cm);Status", 200, 0, 100,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0);

  fHistDcaV0Daughters          = new TH2F("h2DcaV0Daughters", "DCA between daughters;dca(cm);Status", 160, 0, 4,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0Daughters);

  fHistChi2                    = new TH2F("h2Chi2", "V0s chi2;chi2;Status", 33, 0, 33,2,-0.5,1.5);
  fListHist->Add(fHistChi2);

  fHistCosPointAngle           = new TH2F("h2CosPointAngle", "Cosine of V0's pointing angle", 100,0,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngle);

  fHistCosPointAngleZoom       = new TH2F("h2CosPointAngleZoom", "Cosine of V0's pointing angle", 100,0.9,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleZoom);

  // Pt and rapidity distribution:
  fHistPtVsYK0s                = new TH2F("h2PtVsYK0s", "K^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0s);
  fHistPtVsYK0sMI              = new TH2F("h2PtVsYK0sMI", "K^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0sMI);

  fHistPtVsYLambda             = new TH2F("h2PtVsYLambda", "#Lambda^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambda);
  fHistPtVsYLambdaMI           = new TH2F("h2PtVsYLambdaMI", "#Lambda^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambdaMI);

  fHistPtVsYAntiLambda         = new TH2F("h2PtVsYAntiLambda", "#bar{#Lambda}^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambda);
  fHistPtVsYAntiLambdaMI       = new TH2F("h2PtVsYAntiLambdaMI", "#bar{#Lambda}^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambdaMI);

  // Mass:
  fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0);
  fHistMassK0MI                 = new TH1F("h1MassK0MI", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0MI);

  fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);
  fHistMassLambdaMI             = new TH1F("h1MassLambdaMI", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambdaMI);

  fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);
  fHistMassAntiLambdaMI         = new TH1F("h1MassAntiLambdaMI", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambdaMI);

  // invariant mass vs radius
  const Double_t radius[10] = {0.0,2.5,2.9,3.9,7.6,15.0,23.9,37.8,42.8,100.0};
  Int_t NbinRadius        = 9;
  Int_t NbinInvMassLambda = 300;

  fHistMassVsRadiusK0           = new TH2F("h2MassVsRadiusK0", "K^{0} reconstructed;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0);

  fHistMassVsRadiusK0MI         = new TH2F("h2MassVsRadiusK0MI", "K^{0} MI reconstructed;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0MI);
  
  fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda reconstructed;radius (cm);M(p#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambda);

  fHistMassVsRadiusLambdaMI     = new TH2F("h2MassVsRadiusLambdaMI", "#Lambda MI reconstructed;radius (cm);M(p#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambdaMI);

  fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} reconstructed;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",NbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambda);

  fHistMassVsRadiusAntiLambdaMI = new TH2F("h2MassVsRadiusAntiLambdaMI", "#bar{#Lambda} reconstructed;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",NbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambdaMI);


  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fHistArmenterosPodolanskiMI   = new TH2F("h2ArmenterosPodolanskiMI","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);


  //********************************
  // Associated particles histograms
  //********************************
  //Pt distribution
  fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtK0);
  fHistAsMcPtK0MI              = new TH1F("h1AsMcPtK0MI", "K^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtK0MI);

  fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtLambda);
  fHistAsMcPtLambdaMI          = new TH1F("h1AsMcPtLambdaMI", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtLambdaMI);

  fHistAsMcPtAntiLambda        = new TH1F("h1AsMcPtAntiLambda", "#bar{#Lambda}^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtAntiLambda);
  fHistAsMcPtAntiLambdaMI      = new TH1F("h1AsMcPtAntiLambdaMI", "#bar{#Lambda}^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtAntiLambdaMI);

  fHistAsMcPtZoomK0            = new TH1F("h1AsMcPtZoomK0", "K^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0);
  fHistAsMcPtZoomK0MI         = new TH1F("h1AsMcPtZoomK0MI", "K^{0} MI candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0MI);

  fHistAsMcPtZoomLambda        = new TH1F("h1AsMcPtZoomLambda", "#Lambda^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambda);
  fHistAsMcPtZoomLambdaMI      = new TH1F("h1AsMcPtZoomLambdaMI", "#Lambda^{0} MI candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambdaMI);


  // Mass
  fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0);
  fHistPidMcMassK0MI           = new TH1F("h1PidMcMassK0MI", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0MI);

  fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambda);
  fHistPidMcMassLambdaMI       = new TH1F("h1PidMcMassLambdaMI", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambdaMI);
  
  fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambda);
  fHistPidMcMassAntiLambdaMI   = new TH1F("h1PidMcMassAntiLambdaMI", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambdaMI);

  fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0);
  fHistAsMcMassK0MI            = new TH1F("h1AsMcMassK0MI", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0MI);
  
  fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambda);
  fHistAsMcMassLambdaMI        = new TH1F("h1AsMcMassLambdaMI", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambdaMI);

  fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambda);
  fHistAsMcMassAntiLambdaMI    = new TH1F("h1AsMcMassAntiLambdaMI", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambdaMI);


  // invariant mass vs radius
  fHistAsMcMassVsRadiusK0           = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0);

  fHistAsMcMassVsRadiusK0MI         = new TH2F("h2AsMcMassVsRadiusK0MI", "K^{0} MI associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",NbinRadius,radius, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0MI);
  
  fHistAsMcMassVsRadiusLambda       = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",NbinRadius,radius, NbinInvMassLambda, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambda);

  fHistAsMcMassVsRadiusLambdaMI     = new TH2F("h2AsMcMassVsRadiusLambdaMI", "#Lambda MI associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",NbinRadius,radius, NbinInvMassLambda, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambdaMI);

  fHistAsMcMassVsRadiusAntiLambda   = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",NbinRadius,radius,NbinInvMassLambda , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambda);
  
  fHistAsMcMassVsRadiusAntiLambdaMI = new TH2F("h2AsMcMassVsRadiusAntiLambdaMI", "#bar{#Lambda} MI associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",NbinRadius,radius,NbinInvMassLambda , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambdaMI);
    

  // Resolution
  fHistAsMcResxK0                   = new TH1F("h1AsMcResxK0", "K^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxK0);
  fHistAsMcResyK0                   = new TH1F("h1AsMcResyK0", "K^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyK0);
  fHistAsMcReszK0                   = new TH1F("h1AsMcReszK0", "K^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszK0);
  fHistAsMcResrVsRadiusK0           = new TH2F("h2AsMcResrVsRadiusK0", "K^{0} associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusK0);
  fHistAsMcReszVsRadiusK0           = new TH2F("h2AsMcReszVsRadiusK0", "K^{0} associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusK0);

  fHistAsMcResxK0MI                 = new TH1F("h1AsMcResxK0MI", "K^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxK0MI);
  fHistAsMcResyK0MI                 = new TH1F("h1AsMcResyK0MI", "K^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyK0MI);
  fHistAsMcReszK0MI                 = new TH1F("h1AsMcReszK0MI", "K^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszK0MI);
  fHistAsMcResrVsRadiusK0MI         = new TH2F("h2AsMcResrVsRadiusK0MI", "K^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusK0MI);
  fHistAsMcReszVsRadiusK0MI         = new TH2F("h2AsMcReszVsRadiusK0MI", "K^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusK0MI);

  fHistAsMcResxLambda             = new TH1F("h1AsMcResxLambda", "#Lambda^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxLambda);
  fHistAsMcResyLambda             = new TH1F("h1AsMcResyLambda", "#Lambda^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyLambda);
  fHistAsMcReszLambda             = new TH1F("h1AsMcReszLambda", "#Lambda^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszLambda);
  fHistAsMcResrVsRadiusLambda     = new TH2F("h2AsMcResrVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusLambda);
  fHistAsMcReszVsRadiusLambda     = new TH2F("h2AsMcReszVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusLambda);

  fHistAsMcResxLambdaMI           = new TH1F("h1AsMcResxLambdaMI", "#Lambda^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxLambdaMI);
  fHistAsMcResyLambdaMI           = new TH1F("h1AsMcResyLambdaMI", "#Lambda^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyLambdaMI);
  fHistAsMcReszLambdaMI           = new TH1F("h1AsMcReszLambdaMI", "#Lambda^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszLambdaMI);
  fHistAsMcResrVsRadiusLambdaMI   = new TH2F("h2AsMcResrVsRadiusLambdaMI", "#Lambda^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusLambdaMI);
  fHistAsMcReszVsRadiusLambdaMI   = new TH2F("h2AsMcReszVsRadiusLambdaMI", "#Lambda^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusLambdaMI);

  fHistAsMcResxAntiLambda         = new TH1F("h1AsMcResxAntiLambda", "#bar{#Lambda}^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxAntiLambda);
  fHistAsMcResyAntiLambda         = new TH1F("h1AsMcResyAntiLambda", "#bar{#Lambda}^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyAntiLambda);
  fHistAsMcReszAntiLambda         = new TH1F("h1AsMcReszAntiLambda", "#bar{#Lambda}^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszAntiLambda);
  fHistAsMcResrVsRadiusAntiLambda = new TH2F("h2AsMcResrVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusAntiLambda);
  fHistAsMcReszVsRadiusAntiLambda = new TH2F("h2AsMcReszVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusAntiLambda);

  fHistAsMcResxAntiLambdaMI         = new TH1F("h1AsMcResxAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxAntiLambdaMI);
  fHistAsMcResyAntiLambdaMI         = new TH1F("h1AsMcResyAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyAntiLambdaMI);
  fHistAsMcReszAntiLambdaMI         = new TH1F("h1AsMcReszAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszAntiLambdaMI);
  fHistAsMcResrVsRadiusAntiLambdaMI = new TH2F("h2AsMcResrVsRadiusAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusAntiLambdaMI);
  fHistAsMcReszVsRadiusAntiLambdaMI = new TH2F("h2AsMcReszVsRadiusAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusAntiLambdaMI);


}

//________________________________________________________________________
void AliAnalysisTaskESDStrangeMC::Exec(Option_t *) 
{
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //**********************************************
  // Check Monte Carlo information access first:
  //**********************************************
  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event

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

  //**********************************************
  // MC loop
  // Called for each event
  //**********************************************

  //Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
  //Printf("MC primary particles: %d", stack->GetNprimary());

  for (Int_t iTracksMC = 0; iTracksMC < mcEvent->GetNumberOfTracks(); iTracksMC++) {
    AliMCParticle* trackMC = mcEvent->GetTrack(iTracksMC);
    if (!trackMC) {
      Printf("ERROR: Could not receive track %d (mc loop)", iTracksMC);
      continue;
    }
    fHistMCMultiplicity->Fill(mcEvent->GetNumberOfTracks());
    fHistPtMC->Fill(trackMC->Pt());
  } //end MC tracks loop 

  Double_t y =999.9;
  for (Int_t iMc = 0; iMc < stack->GetNprimary(); ++iMc) {  
      TParticle *p0 = stack->Particle(iMc);
      if (!p0) {
	Printf("ERROR: particle with label %d not found in stack (mc loop)", iMc);
	continue;
      }
      y=myRap(p0->Energy(),p0->Pz());
      if (p0->GetPdgCode()==310) fHistMCPtVsYK0s->Fill(p0->Pt(),y);
      else if (p0->GetPdgCode()==3122)  fHistMCPtVsYLambda->Fill(p0->Pt(),y);
      else if (p0->GetPdgCode()==-3122) fHistMCPtVsYAntiLambda->Fill(p0->Pt(),y);

    }// end loop MC primary particles


  //*********************************************
  // ESD loop
  // Called for each event
  //*********************************************

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  //Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  if (!(fESD->GetNumberOfTracks())) {
    //Printf("No ESD track in the event");
  return;
  }
  fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());

  // Track loop to fill a pT spectrum
  for (Int_t iesdTracks = 0; iesdTracks < fESD->GetNumberOfTracks(); iesdTracks++) {
    AliESDtrack* esdtrack = fESD->GetTrack(iesdTracks);
    if (!esdtrack) {
      Printf("ERROR: Could not receive track %d", iesdTracks);
      continue;
    }
  }

  // Primary Vertex
  Double_t  PrimaryVtxPosition[3];
  Double_t  PrimaryVtxCov[6];

  const AliESDVertex *primaryVtx = fESD->GetPrimaryVertex();

  primaryVtx->GetXYZ(PrimaryVtxPosition);
  primaryVtx->GetCovMatrix(PrimaryVtxCov); 

  AliAODVertex *primary = new AliAODVertex(PrimaryVtxPosition, PrimaryVtxCov, primaryVtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);

  fHistPrimaryVertexX->Fill(primary->GetX());
  fHistPrimaryVertexY->Fill(primary->GetY());
  fHistPrimaryVertexZ->Fill(primary->GetZ());


  // V0 variables:
  // to get info from ESD files and fill AliAODVertex:
  Float_t   tdcaPosToPrimVertexXYZ[2], tdcaNegToPrimVertexXYZ[2]; // ..[0] = Impact parameter in XY plane and ..[1] = Impact parameter in Z            
  Double_t  tdcaDaughterToPrimVertex[2];                          // ..[0] = Pos and ..[1] = Neg
  Double_t  tdcaV0Daughters     = 0, tdcaV0ToPrimVertex   = 0;
  Double_t  tMomPos[3];
  Double_t  tMomNeg[3];
  Double_t  V0Position[3];
  Double_t  V0Cov[6];

  // to fill AliAODtrack:
  Double_t  TrackP[3];
  Double_t  TrackPosition[3];
  Double_t  TrackcovTr[21];
  Double_t  Trackpid[10];

  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;
  UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg       = 0;
  Int_t    lCheckPIdK0Short     = 0, lCheckMcK0Short      = 0;
  Int_t    lCheckPIdLambda      = 0, lCheckMcLambda       = 0;
  Int_t    lCheckPIdAntiLambda  = 0, lCheckMcAntiLambda   = 0;

  Double_t mcPosX        = 0,  mcPosY  = 0, mcPosZ  = 0;
  Double_t rcPosX        = 0,  rcPosY  = 0, rcPosZ  = 0;
  Double_t rcPosR        = 0,  mcPosR  = 0;
  Double_t cosPointAngle = 0;
  Double_t deltaPos2     = 0;
  Double_t deltaPos[3];

  Int_t    myStatus  = 0;
  Double_t pt        = 0;

  AliAODTrack  *myPosAodTrack  = new AliAODTrack();
  AliAODTrack  *myNegAodTrack  = new AliAODTrack();
  AliAODVertex *myAODVertex    = new AliAODVertex();
  AliAODv0     *myAODv0        = new AliAODv0();

  // V0 loop
  for (Int_t iV0 = 0; iV0 < fESD->GetNumberOfV0s(); iV0++) {
  
    AliESDv0 *v0 = fESD->GetV0(iV0);
    if (!v0) continue;
    myAODv0->ResetV0();

    // AliAODVertex
    v0->GetXYZ(V0Position[0], V0Position[1], V0Position[2]);
    v0->GetPosCov(V0Cov);
    myAODVertex->SetPosition(V0Position[0],V0Position[1],V0Position[2]);
    myAODVertex->SetCovMatrix(V0Cov);
    myAODVertex->SetChi2perNDF(v0->GetChi2V0());
    myAODVertex->SetID((Short_t)iV0);
    myAODVertex->SetParent(primary);
    myAODVertex->SetType(AliAODVertex::kV0);

    // AliAODtrack (V0 Daughters)
    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *TrackPos = fESD->GetTrack(lIndexTrackPos);
    AliESDtrack *TrackNeg = fESD->GetTrack(lIndexTrackNeg);
    if (!TrackPos || !TrackNeg) {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }
    lLabelTrackPos = (UInt_t)TMath::Abs(TrackPos->GetLabel());
    lLabelTrackNeg = (UInt_t)TMath::Abs(TrackNeg->GetLabel());

    myPosAodTrack->SetID((Short_t)(TrackPos->GetID()));  
    myPosAodTrack->SetLabel(lLabelTrackPos);
    TrackPos->GetPxPyPz(TrackP);
    myPosAodTrack->SetP(TrackP);
    TrackPos->GetXYZ(TrackPosition);
    myPosAodTrack->SetPosition(TrackPosition,kFALSE);
    TrackPos->GetCovarianceXYZPxPyPz(TrackcovTr);
    myPosAodTrack->SetCovMatrix(TrackcovTr);
    TrackPos->GetESDpid(Trackpid);
    myPosAodTrack->SetPID(Trackpid);
    myPosAodTrack->SetCharge((Short_t)(TrackPos->Charge()));
    myPosAodTrack->SetITSClusterMap(TrackPos->GetITSClusterMap());
    myPosAodTrack->SetProdVertex(myAODVertex);
    myPosAodTrack->SetUsedForVtxFit(kTRUE);
    myPosAodTrack->SetUsedForPrimVtxFit(kFALSE);
    myPosAodTrack->SetType(AliAODTrack::kSecondary);
    myPosAodTrack->ConvertAliPIDtoAODPID();

    myNegAodTrack->SetID((Short_t)(TrackNeg->GetID()));
    myNegAodTrack->SetLabel(lLabelTrackNeg);
    TrackNeg->GetPxPyPz(TrackP);
    myNegAodTrack->SetP(TrackP);
    TrackNeg->GetXYZ(TrackPosition);
    myNegAodTrack->SetPosition(TrackPosition,kFALSE);
    TrackNeg->GetCovarianceXYZPxPyPz(TrackcovTr);
    myNegAodTrack->SetCovMatrix(TrackcovTr);
    TrackNeg->GetESDpid(Trackpid);
    myNegAodTrack->SetPID(Trackpid);
    myNegAodTrack->SetCharge((Short_t)(TrackNeg->Charge()));
    myNegAodTrack->SetITSClusterMap(TrackPos->GetITSClusterMap());
    myNegAodTrack->SetProdVertex(myAODVertex);
    myNegAodTrack->SetUsedForVtxFit(kTRUE);
    myNegAodTrack->SetUsedForPrimVtxFit(kFALSE);
    myNegAodTrack->SetType(AliAODTrack::kSecondary);
    myNegAodTrack->ConvertAliPIDtoAODPID();
   
    myAODVertex->AddDaughter(myPosAodTrack);
    myAODVertex->AddDaughter(myNegAodTrack);

    // filling myAODv0
    tdcaV0Daughters    = v0->GetDcaV0Daughters();
    tdcaV0ToPrimVertex = v0->GetD(primary->GetX(),primary->GetY(),primary->GetZ());

    if (TrackPos) TrackPos->GetImpactParameters(tdcaPosToPrimVertexXYZ[0],tdcaPosToPrimVertexXYZ[1]);
    if (TrackNeg) TrackNeg->GetImpactParameters(tdcaNegToPrimVertexXYZ[0],tdcaNegToPrimVertexXYZ[1]);
    tdcaDaughterToPrimVertex[0] = TMath::Sqrt(tdcaPosToPrimVertexXYZ[0]*tdcaPosToPrimVertexXYZ[0]+tdcaPosToPrimVertexXYZ[1]*tdcaPosToPrimVertexXYZ[1]);
    tdcaDaughterToPrimVertex[1] = TMath::Sqrt(tdcaNegToPrimVertexXYZ[0]*tdcaNegToPrimVertexXYZ[0]+tdcaNegToPrimVertexXYZ[1]*tdcaNegToPrimVertexXYZ[1]);

    v0->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
    v0->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 

    myAODv0->Fill(myAODVertex, tdcaV0Daughters, tdcaV0ToPrimVertex, tMomPos, tMomNeg, tdcaDaughterToPrimVertex);
    
    // Checking if the v0s is associated with a Monte Carlo particle
    //AliMCParticle  *lPartPos= mcEvent->GetTrack(lLabelTrackPos);
    //AliMCParticle  *lPartNeg= mcEvent->GetTrack(lLabelTrackNeg);
    TParticle  *lPartPos = stack->Particle(lLabelTrackPos);
    TParticle  *lPartNeg = stack->Particle(lLabelTrackNeg);
    Int_t lPartPosMother = lPartPos->GetFirstMother();
    Int_t lPartNegMother = lPartNeg->GetFirstMother();
    
    lCheckPIdK0Short    = 0; lCheckMcK0Short     = 0;
    lCheckPIdLambda     = 0; lCheckMcLambda      = 0;
    lCheckPIdAntiLambda = 0; lCheckMcAntiLambda  = 0;
    
    if( (lPartPosMother==-1) ||	(lPartNegMother==-1) ) {
      fHistMCDaughterTrack->Fill(1);
    }
    else if( (lPartPos->GetPdgCode()==+211) &&
	     (lPartNeg->GetPdgCode()==-211) ) {
      lCheckPIdK0Short    = 1;
      fHistMCDaughterTrack->Fill(3);
      if ( (lPartPosMother==lPartNegMother) &&
	   (stack->Particle(lPartPosMother)->GetPdgCode()==310) ){
	lCheckMcK0Short  = 1;
	mcPosX = lPartPos->Vx();
	mcPosY = lPartPos->Vy();
	mcPosZ = lPartPos->Vz();
	mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      }
    }
    else if( (lPartPos->GetPdgCode()==+2212) &&
	     (lPartNeg->GetPdgCode()==-211) ) {
      lCheckPIdLambda     = 1;
      fHistMCDaughterTrack->Fill(5);
      if ( (lPartPosMother==lPartNegMother) &&
	   (stack->Particle(lPartPosMother)->GetPdgCode()==3122) ){
	lCheckMcLambda  = 1;
	mcPosX = lPartPos->Vx();
	mcPosY = lPartPos->Vy();
	mcPosZ = lPartPos->Vz();
	mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      }
    }
    else if( (lPartPos->GetPdgCode()==211) &&
	     (lPartNeg->GetPdgCode()==-2212) ) {
      lCheckPIdAntiLambda = 1;
      fHistMCDaughterTrack->Fill(7);
      if ( (lPartPosMother==lPartNegMother) &&
	   (stack->Particle(lPartPosMother)->GetPdgCode()==-3122) ){
	lCheckMcAntiLambda  = 1;
	mcPosX = lPartPos->Vx();
	mcPosY = lPartPos->Vy();
	mcPosZ = lPartPos->Vz();
	mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      }
    }

    // Reconstructed V0 position and Cos pointing angle:
    rcPosX   = myAODv0->DecayVertexV0X();
    rcPosY   = myAODv0->DecayVertexV0Y();
    rcPosZ   = myAODv0->DecayVertexV0Z();
    rcPosR   = TMath::Sqrt(rcPosX*rcPosX+rcPosY*rcPosY);

    deltaPos[0]   = rcPosX - PrimaryVtxPosition[0];
    deltaPos[1]   = rcPosY - PrimaryVtxPosition[1];
    deltaPos[2]   = rcPosZ - PrimaryVtxPosition[2];
    deltaPos2     = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
    cosPointAngle = (deltaPos[0]*(myAODv0->MomV0X()) + deltaPos[1]*(myAODv0->MomV0Y()) + deltaPos[2]*(myAODv0->MomV0Z()))/TMath::Sqrt((myAODv0->Ptot2V0())*deltaPos2);
    


    myStatus = v0->GetOnFlyStatus();
    pt     = TMath::Sqrt(myAODv0->Ptot2V0());

    // filling histograms
    fHistDcaPosToPrimVertex->Fill(myAODv0->DcaPosToPrimVertex(),myStatus);
    fHistDcaNegToPrimVertex->Fill(myAODv0->DcaNegToPrimVertex(),myStatus);
    fHistDcaPosToPrimVertexZoom->Fill(myAODv0->DcaPosToPrimVertex(),myStatus);
    fHistDcaNegToPrimVertexZoom->Fill(myAODv0->DcaNegToPrimVertex(),myStatus);
    fHistRadiusV0->Fill(myAODv0->RadiusV0(),myStatus);
    fHistDecayLengthV0->Fill(myAODv0->DecayLengthV0(PrimaryVtxPosition),myStatus);
    fHistDcaV0Daughters->Fill(myAODv0->DcaV0Daughters(),myStatus);
    fHistChi2->Fill(myAODv0->Chi2V0(),myStatus);
    fHistCosPointAngle->Fill(cosPointAngle,myStatus);
    if (cosPointAngle >= 0.9) fHistCosPointAngleZoom->Fill(cosPointAngle,myStatus);
    if (!myStatus) {
      fHistPtVsYK0s->Fill(pt,myAODv0->RapK0Short());
      fHistPtVsYLambda->Fill(pt,myAODv0->RapLambda());
      fHistPtVsYAntiLambda->Fill(pt,myAODv0->RapLambda());
      fHistArmenterosPodolanski->Fill(myAODv0->AlphaV0(),myAODv0->PtArmV0());
    }
    else {
      fHistPtVsYK0sMI->Fill(pt,myAODv0->RapK0Short());
      fHistPtVsYLambdaMI->Fill(pt,myAODv0->RapLambda());
      fHistPtVsYAntiLambdaMI->Fill(pt,myAODv0->RapLambda());
      fHistArmenterosPodolanskiMI->Fill(myAODv0->AlphaV0(),myAODv0->PtArmV0());
    }
    // K0s associated histograms:
    if (TMath::Abs(myAODv0->RapK0Short()) < 1) {
      switch (myStatus){
      case 0 : 
	fHistMassK0->Fill(myAODv0->MassK0Short());
	fHistMassVsRadiusK0->Fill(myAODv0->RadiusV0(),myAODv0->MassK0Short());
	if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(myAODv0->MassK0Short());
	if(lCheckMcK0Short) {
	  fHistAsMcMassK0->Fill(myAODv0->MassK0Short());
	  fHistAsMcPtK0->Fill(pt);
	  if (pt <= 1) fHistAsMcPtZoomK0->Fill(pt);
	  fHistAsMcResxK0->Fill(rcPosX-mcPosX);
	  fHistAsMcResyK0->Fill(rcPosY-mcPosY);
	  fHistAsMcReszK0->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusK0->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusK0->Fill(rcPosR,rcPosZ-mcPosZ);
	}
	break;
	
      case 1 :
	fHistMassK0MI->Fill(myAODv0->MassK0Short());
	fHistMassVsRadiusK0MI->Fill(myAODv0->RadiusV0(),myAODv0->MassK0Short());
	if(lCheckPIdK0Short) fHistPidMcMassK0MI->Fill(myAODv0->MassK0Short());
	if(lCheckMcK0Short) {
	  fHistAsMcMassK0MI->Fill(myAODv0->MassK0Short());
	  fHistAsMcPtK0MI->Fill(pt);
	  if (pt <= 1) fHistAsMcPtZoomK0MI->Fill(pt);
	  fHistAsMcResxK0MI->Fill(rcPosX-mcPosX);
	  fHistAsMcResyK0MI->Fill(rcPosY-mcPosY);
	  fHistAsMcReszK0MI->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusK0MI->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusK0MI->Fill(rcPosR,rcPosZ-mcPosZ);
	}
	break;
	
      }
    }
    // Lambda and AntiLambda associated histograms:
    if (TMath::Abs(myAODv0->RapLambda()) < 1) {
      switch (myStatus){
      case 0 : 
	fHistMassLambda->Fill(myAODv0->MassLambda());
	fHistMassAntiLambda->Fill(myAODv0->MassAntiLambda());
	fHistMassVsRadiusLambda->Fill(myAODv0->RadiusV0(),myAODv0->MassLambda());
	fHistMassVsRadiusAntiLambda->Fill(myAODv0->RadiusV0(),myAODv0->MassAntiLambda());
	if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(myAODv0->MassLambda());
	else if (lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(myAODv0->MassAntiLambda());
	if(lCheckMcLambda) {
	  fHistAsMcMassLambda->Fill(myAODv0->MassLambda());
	  fHistAsMcPtLambda->Fill(pt);
	  if (pt <= 1) fHistAsMcPtZoomLambda->Fill(pt);
	  fHistAsMcMassVsRadiusLambda->Fill(myAODv0->RadiusV0(),myAODv0->MassLambda());
	  fHistAsMcResxLambda->Fill(rcPosX-mcPosX);
	  fHistAsMcResyLambda->Fill(rcPosY-mcPosY);
	  fHistAsMcReszLambda->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusLambda->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusLambda->Fill(rcPosR,rcPosZ-mcPosZ);
	  
	}
	else if(lCheckMcAntiLambda) {
	  fHistAsMcMassAntiLambda->Fill(myAODv0->MassAntiLambda());
	  fHistAsMcPtAntiLambda->Fill(pt);
	  fHistAsMcMassVsRadiusAntiLambda->Fill(myAODv0->RadiusV0(),myAODv0->MassAntiLambda());
	  fHistAsMcResxAntiLambda->Fill(rcPosX-mcPosX);
	  fHistAsMcResyAntiLambda->Fill(rcPosY-mcPosY);
	  fHistAsMcReszAntiLambda->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusAntiLambda->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusAntiLambda->Fill(rcPosR,rcPosZ-mcPosZ);
	}
	break;
	
      case 1 :
	fHistMassLambdaMI->Fill(myAODv0->MassLambda());
	fHistMassAntiLambdaMI->Fill(myAODv0->MassAntiLambda());
	fHistMassVsRadiusLambdaMI->Fill(myAODv0->RadiusV0(),myAODv0->MassLambda());
	fHistMassVsRadiusAntiLambdaMI->Fill(myAODv0->RadiusV0(),myAODv0->MassAntiLambda());
	if(lCheckPIdLambda) fHistPidMcMassLambdaMI->Fill(myAODv0->MassLambda());
	else if (lCheckPIdAntiLambda) fHistPidMcMassAntiLambdaMI->Fill(myAODv0->MassAntiLambda());
	if(lCheckMcLambda) {
	  fHistAsMcMassLambdaMI->Fill(myAODv0->MassLambda());
	  fHistAsMcPtLambdaMI->Fill(pt);
	  fHistAsMcMassVsRadiusLambdaMI->Fill(myAODv0->RadiusV0(),myAODv0->MassLambda());
	  fHistAsMcResxLambdaMI->Fill(rcPosX-mcPosX);
	  fHistAsMcResyLambdaMI->Fill(rcPosY-mcPosY);
	  fHistAsMcReszLambdaMI->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusLambdaMI->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusLambdaMI->Fill(rcPosR,rcPosZ-mcPosZ);
	}
	else if(lCheckMcAntiLambda) {
	  fHistAsMcMassAntiLambdaMI->Fill(myAODv0->MassAntiLambda());
	  fHistAsMcPtAntiLambdaMI->Fill(pt);
	  fHistAsMcMassVsRadiusAntiLambdaMI->Fill(myAODv0->RadiusV0(),myAODv0->MassAntiLambda());
	  fHistAsMcResxAntiLambdaMI->Fill(rcPosX-mcPosX);
	  fHistAsMcResyAntiLambdaMI->Fill(rcPosY-mcPosY);
	  fHistAsMcReszAntiLambdaMI->Fill(rcPosZ-mcPosZ);
	  fHistAsMcResrVsRadiusAntiLambdaMI->Fill(rcPosR,rcPosR-mcPosR);
	  fHistAsMcReszVsRadiusAntiLambdaMI->Fill(rcPosR,rcPosZ-mcPosZ);
	}
	break;
	
      }
    }
 
  } // end V0 loop

  if (primary) delete primary;
  if (myPosAodTrack) delete myPosAodTrack;
  if (myNegAodTrack) delete myNegAodTrack;
  //if (myAODVertex) delete myAODVertex; // don't !!!
  if (myAODv0) delete myAODv0;
  
  
  // Post output data
  PostData(0, fListHist);
}      

//________________________________________________________________________
void AliAnalysisTaskESDStrangeMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  /*
  fHistPtMC = dynamic_cast<TH1F*> (((TList*)GetOutputData(0))->FindObject("h1PtMC"));
  if (!fHistPtMC) {
    Printf("ERROR: fHistPtMC not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskESDStrangeMC","Pt MC",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPtMC->DrawCopy("E");
*/
}

//----------------------------------------------------------------------------

Double_t myRap(Double_t rE, Double_t rPz)
{
  Double_t lRapidity = 1.e30;
  if(rPz && rE && (rPz != rE) && (1.+(2./((rE/rPz)-1.))>0))
    lRapidity = 0.5*log(1.+(2./((rE/rPz)-1.)));
  return lRapidity;
  } 

//----------------------------------------------------------------------------

