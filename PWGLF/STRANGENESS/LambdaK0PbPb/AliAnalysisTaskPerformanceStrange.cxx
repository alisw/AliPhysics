/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//--------------------------------------------------------------------------
//	        AliAnalysisTaskPerformanceStrange class
//    This task is for a performance study of V0 identification.
//                It works with MC info and ESD tree.
//                 Author: Peter Kalinak  pkalinak@cern.ch kalinak@saske.sk
//--------------------------------------------------------------------------

#include <Riostream.h>

#include <stdio.h>
#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"

#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliMultiplicity.h"

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCHeader.h"
#include "AliAODInputHandler.h"

//#include "AliV0vertexer.h"

#include "AliAODMCParticle.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliLog.h"

#include "AliKFVertex.h"
#include "AliVertexerTracks.h"

#include "AliAnalysisTaskPerformanceStrange.h"
#include "AliAnalysisCentralitySelector.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"



ClassImp(AliAnalysisTaskPerformanceStrange)


//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange()
: AliAnalysisTaskSE(), fAnalysisMC(0), fAnalysisType("infoType"),  fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infoCut"),fDown(0),fUp(0), fESD(0), fListHist(0),fCentrSelector(0),fTracksCuts(0),fPIDResponse(0),fQASelector(0), 


    /////// primary vertex
    fHistMCPrimaryVertexX(0),   fHistMCPrimaryVertexY(0),   fHistMCPrimaryVertexZ(0),
    ////// tracks & multiplicity
    fHistPtTracks(0), 
    fHistMCMultiplicityPrimary(0),  fHistMCMultiplicityTracks(0),  fHistTPCTracks(0), 
    ///// Transverse Momentum
    fHistMCPtAllK0s(0),
    fHistMCPtAllLambda(0),   fHistMCPtAllAntiLambda(0),
    fHistMCPtAllXi(0),  fHistMCPtAllAntiXi(0),
    fHistMCPtAllOmega(0), fHistMCPtAllAntiOmega(0),
    /// Rapidity
    fHistMCRapK0s(0),
    fHistMCRapLambda(0), fHistMCRapAntiLambda(0),
    fHistMCRapXi(0), 

    ///// Transverse Momentum primary
    fHistMCPtK0s(0),
    fHistMCPtLambda(0),  fHistMCPtAntiLambda(0),

    ///////////////////////////////////////////
    ////   ESD

    fHistNumberEvents(0),
    fHistTrackPerEvent(0),
    fHistTPCMult(0),
    fHistTrackletPerEvent(0),
    fHistSPDPrimaryVertexZ(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),

    fHistV0Multiplicity(0),
    ///// inv. mass 
    fHistMassK0(0),
    fHistMassLambda(0),
    fHistMassAntiLambda(0),
    fHistMassXi(0),
    fHistMassAntiXi(0),
    fHistMassOmega(0),
    fHistMassAntiOmega(0),

    //inv mass vs PID
    fHistMassXiVsPID(0),
 
    ///////////////////////////////////////
    fHistPtVsMassK0(0),
    fHistPtVsMassLambda(0),
    fHistPtVsMassAntiLambda(0),

    ////////////////////////////////////////

    fHistArmenterosPodolanski(0),
    fHistK0sMassVsLambdaMass(0),
    fHistTPCsigPLambda(0),
    fHistTPCsigPAntiLambda(0),
    fHistNSigmaProton(0),  
  
    /// Associated histos 
    ///rapidity
    fHistAsMcRapK0(0),
    fHistAsMcRapLambda(0),   fHistAsMcRapAntiLambda(0),

    // pt distribution    /////////////////////
    fHistAsMcPtK0(0),
    fHistAsMcPtLambda(0),   fHistAsMcPtAntiLambda(0),

    /////////////////////////////////////
    fHistPidMcMassK0(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassAntiLambda(0),

    ///inv. mass
    fHistAsMcMassK0(0),
    fHistAsMcMassLambda(0), fHistAsMcMassAntiLambda(0),
    ///transv. momentum
    fHistAsMcPtVsMassK0(0),
    fHistAsMcPtVsMassLambda(0), fHistAsMcPtVsMassAntiLambda(0),

    /////background composition
    fHistCompositionXi(0),
    fHistCompositionAntiXi(0),
    fHistCompositionOmega(0),
  fHistCompositionAntiOmega(0),

  fHistMCIndexes(0)
{
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange(const char *name)
  : AliAnalysisTaskSE(name), fAnalysisMC(0), fAnalysisType("infoType"), fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infocut"),fDown(0),fUp(0), fESD(0), fListHist(),fCentrSelector(0), fTracksCuts(0),fPIDResponse(0),fQASelector(0),

    /////// primary vertex
    fHistMCPrimaryVertexX(0),   fHistMCPrimaryVertexY(0),   fHistMCPrimaryVertexZ(0),
    ////// tracks & multiplicity
    fHistPtTracks(0), 
    fHistMCMultiplicityPrimary(0),  fHistMCMultiplicityTracks(0),  fHistTPCTracks(0), 
    ///// Transverse Momentum
    fHistMCPtAllK0s(0),
    fHistMCPtAllLambda(0),   fHistMCPtAllAntiLambda(0),
    fHistMCPtAllXi(0),  fHistMCPtAllAntiXi(0),
    fHistMCPtAllOmega(0), fHistMCPtAllAntiOmega(0),
    /// Rapidity
    fHistMCRapK0s(0),
    fHistMCRapLambda(0), fHistMCRapAntiLambda(0),
    fHistMCRapXi(0), 

    ///// Transverse Momentum primary
    fHistMCPtK0s(0),
    fHistMCPtLambda(0),  fHistMCPtAntiLambda(0),

    ///////////////////////////////////////////
    ////   ESD

    fHistNumberEvents(0),
    fHistTrackPerEvent(0),
    fHistTPCMult(0),
    fHistTrackletPerEvent(0),
    fHistSPDPrimaryVertexZ(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),

    fHistV0Multiplicity(0),
    ///// inv. mass 
    fHistMassK0(0),
    fHistMassLambda(0),
    fHistMassAntiLambda(0),
    fHistMassXi(0),
    fHistMassAntiXi(0),
    fHistMassOmega(0),
    fHistMassAntiOmega(0),

    //inv mass vs PID
    fHistMassXiVsPID(0),
 
    ///////////////////////////////////////
    fHistPtVsMassK0(0),
    fHistPtVsMassLambda(0),
    fHistPtVsMassAntiLambda(0),

    ////////////////////////////////////////

    fHistArmenterosPodolanski(0),
    fHistK0sMassVsLambdaMass(0),
    fHistTPCsigPLambda(0),
    fHistTPCsigPAntiLambda(0),
    fHistNSigmaProton(0),  
  
    /// Associated histos 
    ///rapidity
    fHistAsMcRapK0(0),
    fHistAsMcRapLambda(0),   fHistAsMcRapAntiLambda(0),

    // pt distribution    /////////////////////
    fHistAsMcPtK0(0),
    fHistAsMcPtLambda(0),   fHistAsMcPtAntiLambda(0),

    /////////////////////////////////////
    fHistPidMcMassK0(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassAntiLambda(0),

    ///inv. mass
    fHistAsMcMassK0(0),
    fHistAsMcMassLambda(0), fHistAsMcMassAntiLambda(0),
    ///transv. momentum
    fHistAsMcPtVsMassK0(0),
    fHistAsMcPtVsMassLambda(0), fHistAsMcPtVsMassAntiLambda(0),

    /////background composition
    fHistCompositionXi(0),
    fHistCompositionAntiXi(0),
    fHistCompositionOmega(0),
    fHistCompositionAntiOmega(0),

    fHistMCIndexes(0)
{
  // Constructor

  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliAnalysisCentralitySelector::Class());
  DefineOutput(3, AliESDtrackCuts::Class());
}
AliAnalysisTaskPerformanceStrange::~AliAnalysisTaskPerformanceStrange() {
  //
  // Destructor
  //
  if (fListHist && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fListHist;     fListHist = 0x0;    }
  if (fCentrSelector && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fCentrSelector;    fCentrSelector = 0x0;    }
  if (fTracksCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fTracksCuts;     fTracksCuts = 0x0;    }


}
//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserCreateOutputObjects() 
{

  //******************
  // Create histograms
  //*******************
  fListHist = new TList();
  fListHist->SetOwner();
  //fListHistCuts = new TList();
  //fListHistCuts->SetOwner();

  // Bo: tbd: condition before allocation (i.e. if (!fHistMCMultiplicityPrimary){...} for each histo...

  //***************
  // MC histograms
  //***************
 
  // Primary Vertex:
  fHistMCPrimaryVertexX          = new TH1F("h1MCPrimaryVertexX", "MC Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexX);

  fHistMCPrimaryVertexY          = new TH1F("h1MCPrimaryVertexY", "MC Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexY);

  fHistMCPrimaryVertexZ          = new TH1F("h1MCPrimaryVertexZ", "MC Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistMCPrimaryVertexZ);

  // Pt track distribution

  fHistPtTracks            = new TH1F("h1PtTracks", "Pt tracks;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistPtTracks);

  // Multiplicity
  fHistMCMultiplicityPrimary           = new TH1F("h1MCMultiplicityPrimary", "MC Primary Particles;NPrimary;Count", 201, -0.5, 200.5);
  fListHist->Add(fHistMCMultiplicityPrimary);

  fHistMCMultiplicityTracks            = new TH1F("h1MCMultiplicityTracks", "MC Tracks;Ntracks;Count", 201, -0.5, 200.5);
  fListHist->Add(fHistMCMultiplicityTracks);

  // Pt Distribution of non-primary particles:
  fHistMCPtAllK0s                      = new TH1F("h1MCPtAllK0s", "Non-primary MC K^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllK0s);

  fHistMCPtAllLambda                   = new TH1F("h1MCPtAllLambda", "Non-primary MC #Lambda^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllLambda);

  fHistMCPtAllAntiLambda               = new TH1F("h1MCPtAllAntiLambda", "Non-primary MC #bar{#Lambda}^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllAntiLambda);

  fHistMCPtAllXi               = new TH1F("h1MCPtAllXi", "Non-primary MC #Xi;p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllXi);

  fHistMCPtAllAntiXi               = new TH1F("h1MCPtAllAntiXi", "Non-primary MC #bar{#Xi};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllAntiXi);

  fHistMCPtAllOmega               = new TH1F("h1MCPtAllOmega", "Non-primary MC #Omega;p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllOmega);

  fHistMCPtAllAntiOmega               = new TH1F("h1MCPtAllAntiOmega", "Non-primary MC #bar{#Omega};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllAntiOmega);


  // Rapidity distribution:
  fHistMCRapK0s                 = new TH1F("h1MCRapK0s", "K^{0};y",160,-4,4);
  fListHist->Add(fHistMCRapK0s);

  fHistMCRapLambda              = new TH1F("h1MCRapLambda", "#Lambda;y",160,-4,4);
  fListHist->Add(fHistMCRapLambda);


  fHistMCRapAntiLambda          = new TH1F("h1MCRapAntiLambda", "#bar{#Lambda};y",160,-4,4);
  fListHist->Add(fHistMCRapAntiLambda);


  // Pt distribution:
  fHistMCPtK0s               = new TH1F("h1MCPtK0s", "K^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtK0s);

  fHistMCPtLambda            = new TH1F("h1MCPtLambda", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtLambda);

  fHistMCPtAntiLambda            = new TH1F("h1MCPtAntiLambda", "#AntiLambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtAntiLambda);

  //***********************************
  // Reconstructed particles histograms
  //***********************************

  // Number of events;
  fHistNumberEvents           = new TH1F("h1NumberEvents", "Number of events; index;Number of Events",10,0,10);
  fListHist->Add(fHistNumberEvents);

  // multiplicity
  fHistTrackPerEvent           = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",10000,0,10000);
  fListHist->Add(fHistTrackPerEvent);

  fHistTPCMult           = new TH1F("h1HistTPCMult", "TPC tracks per event;Number of Tracks;Number of Events",10000,0,10000);
  fListHist->Add(fHistTPCMult);


  fHistTrackletPerEvent       = new TH1F("h1TrackletPerEvent", "Number of tracklets;Number of tracklets per events;Number of events",1000,0,1000);
  fListHist->Add(fHistTrackletPerEvent);

  fHistTPCTracks               = new TH1F("h1TPCTracks","Distribution of TPC tracks;Number of TPC tracks:Number of events",1000,0,10000);
  fListHist->Add(fHistTPCTracks);

  // Primary Vertex:
  fHistSPDPrimaryVertexZ          = new TH1F("h1SPDPrimaryVertexZ", "SPD Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistSPDPrimaryVertexZ);

  fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistPrimaryVertexZ);

  //////K0s///////////////// 2D histos: cut vs on fly status////

  // V0 Multiplicity
  if (!fHistV0Multiplicity) {
    if (fCollidingSystems)
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 200, 0, 40000);
    else
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 10, 0, 10); 
    fListHist->Add(fHistV0Multiplicity);
  }


  // Mass:
  fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 200, 0.4, 0.6);
  fListHist->Add(fHistMassK0);

  fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);

  fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);

  fHistMassXi           = new TH1F("h1MassXi", "#Xi candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.25, 1.4);  // 1321.71
  fListHist->Add(fHistMassXi);

  fHistMassAntiXi           = new TH1F("h1MassAntiXi", "#bar{#Xi} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.25, 1.4);
  fListHist->Add(fHistMassAntiXi);

  fHistMassOmega           = new TH1F("h1MassOmega", "#Omega candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.6, 1.75); //1672.45
  fListHist->Add(fHistMassOmega);

  fHistMassAntiOmega           = new TH1F("h1MassAntiOmega", "#bar{#Omega} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.6, 1.75);
  fListHist->Add(fHistMassAntiOmega);

  //PID vs Mass
  fHistMassXiVsPID           = new TH2F("h1MassXiVsPID", "#Xi candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});PID", 150, 1.25, 1.4,5,0,5);  // 1321.71
  fListHist->Add(fHistMassXiVsPID);



  // Pt Vs Mass
  fHistPtVsMassK0               = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",400, 0.4, 0.6,240,0,12);
  fListHist->Add(fHistPtVsMassK0);

  fHistPtVsMassLambda           = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistPtVsMassLambda);
  
  fHistPtVsMassAntiLambda           = new TH2F("h2PtVsMassAntiLambda","#AntiLambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistPtVsMassAntiLambda);

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///Armenteros Podolansky
  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fListHist->Add(fHistArmenterosPodolanski);

  ///Inv. Mass K0s vs Inv. Mass. Lambda
  fHistK0sMassVsLambdaMass      = new TH2F("h2HistK0sMassVsLambdaMass","K^{0} vs #Lambda^{0} candidates; M(#pi^{+}#pi^{-}) (GeV/c^{2}); M(p#pi^{-}) (GeV/c^{2})",200, 0.4, 0.6,140, 1.06, 1.2);
  fListHist->Add(fHistK0sMassVsLambdaMass);

  //dE/dx vs P daughters
  fHistTPCsigPLambda                            = new TH2F("h2TPCsignalVsPLambda","TPC signal Vs p_{t} daughters;  p (GeV/c);TPC signal",1000,0,2,1000,0,1000);
  fListHist->Add(fHistTPCsigPLambda);


  fHistTPCsigPAntiLambda                            = new TH2F("h2TPCsignalVsPAntiLambda","TPC signal Vs p_{t} daughters;  p (GeV/c);TPC signal",1000,0,2,1000,0,1000);
  fListHist->Add(fHistTPCsigPAntiLambda);
 

  fHistNSigmaProton                          =new TH1F("h1NSigmaProton","Number of sigmas for proton;;Count",600,0,6);
  fListHist->Add(fHistNSigmaProton);


  //********************************
  // Associated particles histograms
  //********************************

  // Rap distribution
  fHistAsMcRapK0                = new TH1F("h1AsMcRapK0", "K^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapK0);

  fHistAsMcRapLambda            = new TH1F("h1AsMcRapLambda", "#Lambda^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapLambda);

  fHistAsMcRapAntiLambda        = new TH1F("h1AsMcRapAntiLambda", "#bar{#Lambda}^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapAntiLambda);

  //Pt distribution
  fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtK0);

  fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtLambda);

  fHistAsMcPtAntiLambda            = new TH1F("h1AsMcPtAntiLambda", "#AntiLambda^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtAntiLambda);

  // Mass
  fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0);

  fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambda);
  
  fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambda);
 
  fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0);
  
  fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambda);

  fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambda);

  //Pt versus Mass
  fHistAsMcPtVsMassK0               = new TH2F("h2AsMcPtVsMassK0","K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassK0);

  fHistAsMcPtVsMassLambda           = new TH2F("h2AsMcPtVsMassLambda","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassLambda);

  fHistAsMcPtVsMassAntiLambda       = new TH2F("h2AsMcPtVsMassAntiLambda","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassAntiLambda);

  fHistCompositionXi = new TH1F("h1CompisitionXi","background composition of Xi;Part ID;PID",8000,-4000,4000);
  fListHist->Add(fHistCompositionXi);

  fHistCompositionAntiXi = new TH1F("h1CompisitionAntiXi","background composition of AntiXi;Part ID;Counts",8000,-4000,4000);
  fListHist->Add(fHistCompositionAntiXi);

  fHistCompositionOmega = new TH1F("h1CompisitionOmega","background composition of Omega;Part ID;Counts",8000,-4000,4000);
  fListHist->Add(fHistCompositionOmega);

  fHistCompositionAntiOmega = new TH1F("h1CompisitionAntiOmega","background composition of AntiOmega;Part ID;Counts",8000,-4000,4000);
  fListHist->Add(fHistCompositionAntiOmega);

  fHistMCIndexes = new TH1I("h1MCIndexes","MC Indexes;Index;Counts",210,-10,200);
  fListHist->Add(fHistMCIndexes);

  PostData(1, fListHist);
  PostData(2, fCentrSelector);
  PostData(3, fTracksCuts);
}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliStack* stack = NULL;
  //  TClonesArray *mcArray = NULL;
  TArrayF mcPrimaryVtx;

  fESD=(AliESDEvent *)InputEvent();

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  AliVEvent* lEvent = InputEvent();
  
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }

  ///////
  // PID
  ///////

  if (fUsePID.Contains("withPID")) {
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }

  fHistNumberEvents->Fill(0.5);
  fHistNumberEvents->Fill(1.5);  
  // Centrality selection

  static AliESDtrackCuts * trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1); 
  Bool_t isCentralitySelected = fCentrSelector->IsCentralityBinSelected(fESD,trackCuts);
  if(!isCentralitySelected) return;

  // Remove Events with no tracks
  if (!(fESD->GetNumberOfTracks()))  return;

  //*************************************
  // Cut used:
  //*************************************
      
  // Cut Rapidity:
  Double_t lCutRap  = 0.5;

  // cut Pseudorapidity:

  Double_t lCutPseudorap = 0.8;

  // If PID is used:
  Double_t lLimitPPID    = 0.7;
  Float_t cutNSigmaLowP  = 1E3;
  Float_t cutNSigmaHighP = 1E3;
  if (fUsePID.Contains("withPID")) {
    cutNSigmaLowP  = 4.0;
    cutNSigmaHighP = 4.0;
  }
  // Cut Daughters pt (GeV/c):
  //Double_t cutMinPtDaughter = 0.160;

  // Cut primary vertex:
  Double_t cutPrimVertex = 10.0;

  // cut ctau
  Double_t cutcTauL = 3*7.89;
  Double_t cutcTauK0s = 3*2.68;

  // Min number of TPC clusters:
  // Int_t nbMinTPCclusters = 80;
  // PID flags:
  Int_t LambdaPID = 0;
  Int_t AntiLambdaPID = 0;

  //////******************************************
  ////// Access MC: ******************************
  //////******************************************

  if (fAnalysisMC) {
    if(fAnalysisType == "ESD") {
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
      stack = mcEvent->Stack();
      if (!stack) {
	Printf("ERROR: Could not retrieve stack");
	return;
      }
      
      AliGenEventHeader* mcHeader=mcEvent->GenEventHeader();
      if(!mcHeader) return;
      mcHeader->PrimaryVertex(mcPrimaryVtx);

      if (TMath::Abs(mcPrimaryVtx.At(2)) > cutPrimVertex) return;  /// cut on z of prim. vertex !!!!!!      
    }
    
    /*
    // PID parameters for MC simulations:
    fAlephParameters[0] = 2.15898e+00/50.;
    fAlephParameters[1] = 1.75295e+01;
    fAlephParameters[2] = 3.40030e-09;
    fAlephParameters[3] = 1.96178e+00;
    fAlephParameters[4] = 3.91720e+00; 
    */
  }

  //**********************************************
  // MC loop
  //**********************************************

  //  Double_t lmcPrimVtxR      = 0;

  Int_t lNbMCPrimary        = 0;
  Int_t lNbMCPart           = 0;

  Int_t lPdgcodeCurrentPart = 0;
  Double_t lRapCurrentPart  = 0;
  Double_t lPtCurrentPart   = 0;
  
  Int_t lComeFromSigma      = 0;
  
  // Production Radius
  Double_t mcPosX     = 0.0,  mcPosY      = 0.0,  mcPosZ      = 0.0;
  Double_t mcPosR     = 0.0;

  // Decay Radius
  Double_t mcDecayPosX = 0, mcDecayPosY = 0, mcDecayPosR = 0;

  // current mc particle 's mother
  Int_t iCurrentMother  = 0, lPdgCurrentMother    = 0;
  //  Bool_t lCurrentMotherIsPrimary;

  // variables for multiple reconstruction studies:
  Int_t id0           = 0, id1          = 0;
  Int_t lNtimesReconstructedK0s   = 0, lNtimesReconstructedLambda   = 0, lNtimesReconstructedAntiLambda   = 0;
  // Int_t lNtimesReconstructedK0sMI = 0, lNtimesReconstructedLambdaMI = 0, lNtimesReconstructedAntiLambdaMI = 0;

  ////********************************
  ////Start loop over MC particles****
  ////********************************
  
  if (fAnalysisMC) {

    // Primary vertex
    fHistMCPrimaryVertexX->Fill(mcPrimaryVtx.At(0));
    fHistMCPrimaryVertexY->Fill(mcPrimaryVtx.At(1));
    fHistMCPrimaryVertexZ->Fill(mcPrimaryVtx.At(2));
  
    //    lmcPrimVtxR = TMath::Sqrt(mcPrimaryVtx.At(0)*mcPrimaryVtx.At(0)+mcPrimaryVtx.At(1)*mcPrimaryVtx.At(1));

    if(fAnalysisType == "ESD") {
    
      lNbMCPrimary = stack->GetNprimary(); 
     lNbMCPart    = stack->GetNtrack();
      
      fHistMCMultiplicityPrimary->Fill(lNbMCPrimary);
      fHistMCMultiplicityTracks->Fill(lNbMCPart);
      
      for (Int_t iMc = 0; iMc < (stack->GetNtrack()); iMc++) {  
	TParticle *p0 = stack->Particle(iMc);
	if (!p0) {
	  //Printf("ERROR: particle with label %d not found in stack (mc loop)", iMc);
	  continue;
	}
	lPdgcodeCurrentPart = p0->GetPdgCode();
	
	// Keep only K0s, Lambda and AntiLambda, Xi-, antiXi+ and Omega-, antiOmega+:
	if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) && (lPdgcodeCurrentPart != 3312 ) && (lPdgcodeCurrentPart != -3312) && (lPdgcodeCurrentPart != 3334)  && (lPdgcodeCurrentPart != -3334) ) continue;
	
	lRapCurrentPart   = MyRapidity(p0->Energy(),p0->Pz());
	//lEtaCurrentPart   = p0->Eta();
	lPtCurrentPart    = p0->Pt();
	iCurrentMother    = p0->GetFirstMother();

	//	lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();
	if (iCurrentMother == -1){lPdgCurrentMother=0; } else {lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();}	

	mcPosX = p0->Vx();
	mcPosY = p0->Vy();
	mcPosZ = p0->Vz();
	mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
	
	id0  = p0->GetDaughter(0);
	id1  = p0->GetDaughter(1);

	// Decay Radius and Production Radius
	if ( id0 <= lNbMCPart && id0 > 0 && id1 <= lNbMCPart && id1 > 0) {
	  TParticle *pDaughter0 = stack->Particle(id0);
	  //	  TParticle *pDaughter1 = stack->Particle(id1);
	  //	  lPdgCurrentDaughter0 = pDaughter0->GetPdgCode();
	  //	  lPdgCurrentDaughter1 = pDaughter1->GetPdgCode();
	  
	  mcDecayPosX = pDaughter0->Vx();
	  mcDecayPosY = pDaughter0->Vy();
	  mcDecayPosR = TMath::Sqrt(mcDecayPosX*mcDecayPosX+mcDecayPosY*mcDecayPosY);
	}
	else  {
	  //Printf("ERROR: particle with label %d and/or %d not found in stack (mc loop)", id0,id1);
	  mcDecayPosR = -1.0;
	}
	
	if (lPdgcodeCurrentPart==310)   {
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllK0s->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==3122)  {
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllLambda->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==-3122) {
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==3312) { // Xi
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllXi->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==-3312) { //Anti Xi
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllAntiXi->Fill(lPtCurrentPart);
	}	
	else if (lPdgcodeCurrentPart==3334) {  //Omega
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllXi->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==-3334) { //Anti Omega
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllAntiXi->Fill(lPtCurrentPart);
	}	

	if ( ( ( TMath::Abs(lPdgCurrentMother) == 3212)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3224)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3214)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3114) )
	     && ( iCurrentMother <= lNbMCPrimary )
	     ) lComeFromSigma = 1;
	else lComeFromSigma = 0;
      
	//*********************************************
	// Now keep only primary particles   
	// if ( ( iMc > lNbMCPrimary ) && (!lComeFromSigma) ) continue;
	
	//*************************************//
        // new definition of primary particles //
        //*************************************//
 
	Double_t dx = 0;
	Double_t dy = 0;
	Double_t dz = 0;
	Double_t ProdDistance = 0;

	dx = ( (mcPrimaryVtx.At(0)) - (p0->Vx()) ); 
	dy = ( (mcPrimaryVtx.At(1)) - (p0->Vy()) );
	dz = ( (mcPrimaryVtx.At(2)) - (p0->Vz()) );

	ProdDistance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

	if (ProdDistance > 0.001) continue; // secondary V0
    
	lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;

        // Rap distribution
        if (lPdgcodeCurrentPart==310) {
	  fHistMCRapK0s->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==3122) {
	  fHistMCRapLambda->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==-3122) {
	  fHistMCRapAntiLambda->Fill(lRapCurrentPart);
	}
	/*	
	if (lPdgcodeCurrentPart==3312) {
	  fHistMCRapXi->Fill(lRapCurrentPart);
	}
	*/
	// Rapidity Cut
	if (TMath::Abs(lRapCurrentPart) > lCutRap) continue;
 
	if (lPdgcodeCurrentPart==310) {
	  fHistMCPtK0s->Fill(lPtCurrentPart);
	}
	else 
	  if (lPdgcodeCurrentPart==3122) {
	    fHistMCPtLambda->Fill(lPtCurrentPart);	  
	  }
	  else 
	    if (lPdgcodeCurrentPart==-3122) {
	      fHistMCPtAntiLambda->Fill(lPtCurrentPart);	  
	    }
 
      } // end loop ESD MC
    } // end ESD condition
  } // End Loop over MC condition
  
  //************************************
  // ESD loop 
  //************************************

  //  Double_t  lLambdaMass = 1.115683;  //PDG
  Double_t lPLambda = 0;
  // Double_t lPAntiLambda = 0;
  //Double_t lPK0s = 0;
  Double_t lMagneticField = 999;

  //Multiplcity:
  Int_t    nv0sTot= 0, nv0s = 0;
  Int_t    nCasTot= 0;

  // Qualities
  Double_t v0q = 0, kinCasQual = 0;
  //  Int_t nv0sMI =0;   
  // Variables:
  Double_t  lV0Position[3];
 
  Double_t lDcaPosToPrimVertex = 0;
  Double_t lDcaNegToPrimVertex = 0;
  Double_t lDcaV0Daughters     = 0;
  Double_t lV0cosPointAngle    = 0;
  Double_t lChi2V0             = 0;
  Double_t lV0DecayLength      = 0;
  Double_t lV0tDecayLength     = 0; //transverse decay length in xy plain
  Double_t lV0Radius           = 0;
  //  Double_t lDcaV0ToPrimVertex  = 0;
  Double_t lcTauLambda         = 0;  
  Double_t lcTauAntiLambda     = 0;
  Double_t lcTauK0s            = 0;
  Int_t    lOnFlyStatus        = 0;
  

  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lInvMassXi = 0,lInvMassAntiXi = 0,lInvMassOmega = 0,lInvMassAntiOmega = 0;
  Double_t lPtK0s      = 0, lPtLambda      = 0, lPtAntiLambda      = 0;
  Double_t lRapK0s     = 0, lRapLambda     = 0, lRapAntiLambda     = 0;
  //  Double_t lEtaK0s     = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
  Double_t lAlphaV0      = 0, lPtArmV0       = 0;

  // to study Associated V0s:
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg         = 0;
  UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg         = 0;
  Int_t    lCheckPIdK0Short     = 0, lCheckMcK0Short        = 0;
  Int_t    lCheckPIdLambda      = 0, lCheckMcLambda         = 0;
  Int_t    lCheckPIdAntiLambda  = 0, lCheckMcAntiLambda     = 0;
  Int_t    lCheckSecondaryK0s   = 0, lCheckSecondaryLambda  = 0, lCheckSecondaryAntiLambda  = 0;
  //  Int_t    lCheckGamma          = 0;
  Double_t mcPosMotherX         = 0, mcPosMotherY           = 0, mcPosMotherZ  = 0;
  Double_t mcPosMotherR         = 0;
  //  Double_t mcMotherPt           = 0;

  Int_t lIndexPosMother        = 0;
  Int_t lIndexNegMother        = 0;
  Int_t lIndexMotherOfMother   = 0;
  Int_t lPDGCodePosDaughter    = 0;
  Int_t lPDGCodeNegDaughter    = 0;
  Int_t lPdgcodeMother         = 0;
  Int_t lPdgcodeMotherOfMother = 0;

  // Reconstructed position
  Double_t rcPosXK0s        = 0,  rcPosYK0s        = 0, rcPosZK0s        = 0;
  Double_t rcPosRK0s        = 0;
  Double_t rcPosXLambda     = 0,  rcPosYLambda     = 0, rcPosZLambda     = 0;
  Double_t rcPosRLambda     = 0;
  Double_t rcPosXAntiLambda = 0,  rcPosYAntiLambda = 0, rcPosZAntiLambda = 0;
  Double_t rcPosRAntiLambda = 0;

  // Pt resolution
  Double_t deltaPtK0s  = 0, deltaPtLambda  = 0, deltaPtAntiLambda  = 0;

  // Daughters
  AliESDtrack  *myTrackPos  = NULL;
  AliESDtrack  *myTrackNeg  = NULL;

  // Daughters' momentum:
  Double_t  lMomPos[3] = {999,999,999};
  Double_t  lMomNeg[3] = {999,999,999};

  // Inner Wall parameters:
  Double_t  lMomInnerWallPos =999, lMomInnerWallNeg = 999;

  // PID
  Float_t nSigmaPosPion   = 0;
  Float_t nSigmaNegPion   = 0;

  Float_t nSigmaPosProton = 0;
  Float_t nSigmaNegProton = 0;
  
// Bachelor                                                                                                                                                                  
Bool_t lIsBachelorKaonForTPC = kFALSE;
Bool_t lIsBachelorPionForTPC = kFALSE;

// Negative V0 daughter                                                                                                                                                      
Bool_t lIsNegPionForTPC   = kFALSE;
Bool_t lIsNegProtonForTPC = kFALSE;

// Positive V0 daughter                                                                                                                                                      
Bool_t lIsPosPionForTPC   = kFALSE;
Bool_t lIsPosProtonForTPC = kFALSE;
      
 Double_t lEffMassXi      = 0. ;
 Double_t lChi2Xi         = -1. ;
 Double_t lDcaXiDaughters = -1. ;
 Double_t lXiCosineOfPointingAngle = -1. ;
 Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
 Double_t lXiRadius2D = -1000. ;
 Double_t lXiRadius3D = -1000. ;

 Int_t lIDs = 0;
 Int_t lIndexBachelorMother = 0;
 Int_t lIndexBachelor = 0;
  //  Int_t lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
  Int_t lCheckPIDLambdaPosDaughter     = 0;// lCheckPIDLambdaNegDaughter     = 0;
  //  Int_t lCheckPIDAntiLambdaPosDaughter = 0;
  Int_t lCheckPIDAntiLambdaNegDaughter = 0;
  
  //****************************************************
  // Primary Vertex cuts &
  // Magnetic field and Quality tracks cuts 
  //****************************************************
  Double_t  lPrimaryVtxPosition[3];
  Double_t  lPrimaryVtxCov[6];
  Double_t  lPrimaryVtxChi2 = 999;
  Double_t  lResPrimaryVtxX = 999;
  Double_t  lResPrimaryVtxY = 999;
  Double_t  lResPrimaryVtxZ = 999;
     
  AliAODVertex *myPrimaryVertex = NULL;
  //const AliVVertex *mySPDPrimaryVertex = NULL;
    
  const AliMultiplicity *myMultiplicty = ((AliESDEvent*)fESD)->GetMultiplicity();

  if(fAnalysisType == "ESD") {  
    ////////////////////////////////////////////////////////////////////////////////////
    //////   Best Primary Vertex:  
    if(fCollidingSystems ==0 || fCollidingSystems == 1){ //pp, PbPb Analysis
    const AliESDVertex *myBestPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertex();
    myBestPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertex();
    if (!myBestPrimaryVertex) return;
    if (!myBestPrimaryVertex->GetStatus()) return;

    fHistNumberEvents->Fill(3.5);

    myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
    myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
    if ( ( TMath::Abs(lPrimaryVtxPosition[2]) ) > cutPrimVertex) return ; //// cut on z of prim. vertex!!!!!
    fHistNumberEvents->Fill(4.5);    
    lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
    lResPrimaryVtxX = myBestPrimaryVertex->GetXRes();
    lResPrimaryVtxY = myBestPrimaryVertex->GetYRes();
    lResPrimaryVtxZ = myBestPrimaryVertex->GetZRes();
    // remove TPC-only primary vertex : retain only events with tracking + SPD vertex
    const AliESDVertex *mySPDPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertexSPD();
    if (!mySPDPrimaryVertex) return;
    fHistSPDPrimaryVertexZ->Fill(mySPDPrimaryVertex->GetZ());
    const AliESDVertex *myPrimaryVertexTracking = ((AliESDEvent*)fESD)->GetPrimaryVertexTracks();
    if (!myPrimaryVertexTracking) return;
    if (!mySPDPrimaryVertex->GetStatus() && !myPrimaryVertexTracking->GetStatus() ) return;
    fHistNumberEvents->Fill(5.5);
    fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());   
    myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    if (!myPrimaryVertex) return;
  }
    ///  pPb analysis
    if(fCollidingSystems == 2){  //twiky https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAVertexSelectionStudies 
      const AliESDVertex* trkVtx = ((AliESDEvent*)fESD)->GetPrimaryVertex();
      if (!trkVtx || trkVtx->GetNContributors()<=0) return;
      TString vtxTtl = trkVtx->GetTitle();
      if (!vtxTtl.Contains("VertexerTracks")) return;
      Float_t zvtx = trkVtx->GetZ();
      const AliESDVertex* spdVtx = ((AliESDEvent*)fESD)->GetPrimaryVertexSPD();
      if (spdVtx->GetNContributors()<=0) return;
      TString vtxTyp = spdVtx->GetTitle();
      Double_t cov[6]={0};
      spdVtx->GetCovarianceMatrix(cov);
      Double_t zRes = TMath::Sqrt(cov[5]);
      if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
      if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
      if (TMath::Abs(zvtx) > cutPrimVertex) return;
      trkVtx->GetXYZ(lPrimaryVtxPosition);
    }

    // Number of Tracklets:
    fHistTrackletPerEvent->Fill(myMultiplicty->GetNumberOfTracklets());
    lMagneticField = ((AliESDEvent*)fESD)->GetMagneticField();
    fHistTPCTracks->Fill(AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fESD, kTRUE));
    //////simple chack for multiplicity////////////////////////////////////////////////////////
    Int_t i =0;

    for (Int_t jTracks=0;jTracks<fESD->GetNumberOfTracks();jTracks++){
		
      AliESDtrack* tPCtrack=fESD->GetTrack(jTracks);
      Float_t xy=0;
      Float_t z=0;
      tPCtrack->GetImpactParameters(xy,z);
      if ((fTracksCuts->IsSelected(tPCtrack))&&(xy<1.0)&&(z<1.0)) {i=i+1;}
    }
    fHistTPCMult->Fill(i);
    /////////////////////////////////////////////////////////////////////////////////////////
      }

  fHistPrimaryVertexX->Fill(lPrimaryVtxPosition[0]);
  fHistPrimaryVertexY->Fill(lPrimaryVtxPosition[1]);
  fHistPrimaryVertexZ->Fill(lPrimaryVtxPosition[2]);
  //Double_t lrcPrimVtxR = TMath::Sqrt(lPrimaryVtxPosition[0]*lPrimaryVtxPosition[0]+lPrimaryVtxPosition[0]*lPrimaryVtxPosition[0]);
 
  ////*************************
  //// Cascade loop ***********
  ////*************************
  nCasTot = fESD->GetNumberOfCascades();

  for (Int_t iCas = 0; iCas < nCasTot; iCas++) {

    AliESDcascade *Cas = ((AliESDEvent*)fESD)->GetCascade(iCas);
      if (!Cas) continue;

      lIndexBachelor = Cas->GetBindex();
      fHistMCIndexes->Fill(lIndexBachelor);


      UInt_t lIdxPosXi        = (UInt_t) TMath::Abs( Cas->GetPindex() );
      UInt_t lIdxNegXi        = (UInt_t) TMath::Abs( Cas->GetNindex() );
      UInt_t lBachIdx         = (UInt_t) TMath::Abs( lIndexBachelor );
      // Care track label can be negative in MC production (linked with the track quality)                                                                                 
      // However = normally, not the case for track index ...                                                                                                              

      // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)                                     
      if(lBachIdx == lIdxNegXi) {
	AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
      }
      if(lBachIdx == lIdxPosXi) {
	AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
      }

      lEffMassXi                      = Cas->GetEffMassXi();
      lChi2Xi                         = Cas->GetChi2Xi();
      lDcaXiDaughters                 = Cas->GetDcaXiDaughters();
      lXiCosineOfPointingAngle        = Cas->GetCascadeCosineOfPointingAngle( lPrimaryVtxPosition[0],
									     lPrimaryVtxPosition[1],
									     lPrimaryVtxPosition[2] );
      // Take care : the best available vertex should be used (like in AliCascadeVertexer)                                                                                 

      Cas->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );
      lXiRadius2D    = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
      lXiRadius3D    = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] +  lPosXi[2]*lPosXi[2]);


      AliESDtrack *pTrackXi           = fESD->GetTrack( lIdxPosXi );
      AliESDtrack *nTrackXi           = fESD->GetTrack( lIdxNegXi );
      AliESDtrack *bachTrackXi        = fESD->GetTrack( lBachIdx );
      if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
	AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
	continue;
      }
            
      // Bachelor                                                                                                                                                                  
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;

      // Negative V0 daughter                                                                                                                                                      
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;

      // Positive V0 daughter                                                                                                                                                      
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
      

      if( bachTrackXi->Charge() < 0 ) {
	v0q = 0;
      kinCasQual =  Cas->ChangeMassHypothesis(v0q,3312);
      lInvMassXi = Cas->GetEffMassXi();
	v0q = 0;
      kinCasQual =  Cas->ChangeMassHypothesis(v0q,3334);
      lInvMassOmega = Cas->GetEffMassXi();
      }

      if( bachTrackXi->Charge() > 0 ) {
	v0q = 0;
      kinCasQual =  Cas->ChangeMassHypothesis(v0q,-3312);
      lInvMassAntiXi = Cas->GetEffMassXi();
	v0q = 0;
      kinCasQual =  Cas->ChangeMassHypothesis(v0q,-3334);
      lInvMassAntiOmega = Cas->GetEffMassXi();
      }


      if( lIsBachelorPionForTPC == kTRUE && lIsNegPionForTPC == kTRUE && lIsPosProtonForTPC == kTRUE){
      fHistMassXi->Fill(lInvMassXi);
      if(fAnalysisMC && (lInvMassXi > 1.25) && (lInvMassXi < 1.4)){
	TParticle *Bachelor = stack->Particle(lIndexBachelor);        
		//        lIndexBachelorMother = Bachelor->GetMother(1);
	lIndexBachelorMother = Bachelor->GetFirstMother();
	if(lIndexBachelorMother == -1) continue;
        if(lIndexBachelorMother == 1){  fHistMassXiVsPID->Fill(lInvMassXi,0);fHistCompositionXi->Fill(0); continue;} //combinatorial
	lIDs = stack->Particle(lIndexBachelorMother)->GetPdgCode();
	if(lIDs == 111) { fHistMassXiVsPID->Fill(lInvMassXi,1); } //pi0
	if(lIDs == 213) { fHistMassXiVsPID->Fill(lInvMassXi,2); } //rho+
	if(lIDs == 2212) { fHistMassXiVsPID->Fill(lInvMassXi,3); } //Proton
	if(lIDs == 2112) { fHistMassXiVsPID->Fill(lInvMassXi,4); } //Neutron
        fHistCompositionXi->Fill(lIDs);      
       }

      }
      if( lIsBachelorKaonForTPC == kTRUE && lIsNegPionForTPC == kTRUE && lIsPosProtonForTPC == kTRUE){
      fHistMassOmega->Fill(lInvMassOmega);
      if(fAnalysisMC && lInvMassOmega > 1.6 && lInvMassOmega < 1.75){
	TParticle *Bachelor = stack->Particle(lIndexBachelor);        
		//        lIndexBachelorMother = Bachelor->GetMother(1);
	lIndexBachelorMother = Bachelor->GetFirstMother();
	if(lIndexBachelorMother == -1) continue;
	if(lIndexBachelorMother == 1) { fHistCompositionOmega->Fill(0); continue;}
	lIDs = stack->Particle(lIndexBachelorMother)->GetPdgCode();
      fHistCompositionOmega->Fill(lIDs);      
       }

      }
      if( lIsBachelorPionForTPC == kTRUE && lIsPosPionForTPC == kTRUE && lIsNegProtonForTPC == kTRUE){
      fHistMassAntiXi->Fill(lInvMassAntiXi);
      if(fAnalysisMC && lInvMassAntiXi > 1.25 && lInvMassAntiXi < 1.4){
	TParticle *Bachelor = stack->Particle(lIndexBachelor);        
		//        lIndexBachelorMother = Bachelor->GetMother(1);
	lIndexBachelorMother = Bachelor->GetFirstMother();
	if(lIndexBachelorMother == -1) continue;
	if(lIndexBachelorMother == 1) { fHistCompositionAntiXi->Fill(0); continue;}
	lIDs = stack->Particle(lIndexBachelorMother)->GetPdgCode();
      fHistCompositionAntiXi->Fill(lIDs);      
       }

      }
      if( lIsBachelorKaonForTPC == kTRUE && lIsPosPionForTPC == kTRUE && lIsNegProtonForTPC == kTRUE){
      fHistMassAntiOmega->Fill(lInvMassAntiOmega);
      if(fAnalysisMC && lInvMassAntiOmega > 1.6 && lInvMassAntiOmega < 1.75){
	TParticle *Bachelor = stack->Particle(lIndexBachelor);        
		//        lIndexBachelorMother = Bachelor->GetMother(1);
	lIndexBachelorMother = Bachelor->GetFirstMother();
	if(lIndexBachelorMother == -1) continue;
	if(lIndexBachelorMother == 1) { fHistCompositionAntiOmega->Fill(0); continue;}
	lIDs = stack->Particle(lIndexBachelorMother)->GetPdgCode();
      fHistCompositionAntiOmega->Fill(lIDs);      
       }

      }
           
  
    
  
  } /// end of Cascade loop

  
  ////*************************
  //// V0 loop ****************
  ////*************************
    
  nv0sTot = fESD->GetNumberOfV0s();
  if (!nv0sTot) fHistNumberEvents->Fill(6.5);

  for (Int_t iV0 = 0; iV0 < nv0sTot; iV0++) {
     
    lIndexPosMother     = 0; lIndexNegMother     = 0; lIndexMotherOfMother       = 0;
    lCheckPIdK0Short    = 0; lCheckMcK0Short     = 0; lCheckSecondaryK0s         = 0;
    lCheckPIdLambda     = 0; lCheckMcLambda      = 0; lCheckSecondaryLambda      = 0;
    lCheckPIdAntiLambda = 0; lCheckMcAntiLambda  = 0; lCheckSecondaryAntiLambda  = 0;       
    lComeFromSigma      = 0; //lCheckGamma = 0;
    
    if(fAnalysisType == "ESD") {

      AliESDv0 *v0 = ((AliESDEvent*)fESD)->GetV0(iV0);
      if (!v0) continue;
      //      if ((v0->Pt())<0.6) continue;
     
      // V0's Daughters
      lIndexTrackPos = TMath::Abs(v0->GetPindex());
      lIndexTrackNeg = TMath::Abs(v0->GetNindex());
      AliESDtrack *myTrackPosTest = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);
      AliESDtrack *myTrackNegTest = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);
      if (!myTrackPosTest || !myTrackNegTest) {
	Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
	continue;
      }
      // Remove like-sign
      if ((Int_t)myTrackPosTest->GetSign() == (Int_t)myTrackNegTest->GetSign()){
	continue;
      } 
     
      // VO's main characteristics to check the reconstruction cuts
      lOnFlyStatus       = v0->GetOnFlyStatus();
      lChi2V0            = v0->GetChi2V0();
      lDcaV0Daughters    = v0->GetDcaV0Daughters();
      //      lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
      lV0cosPointAngle   = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);

      v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);

      lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
      lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
		                   TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
		                   TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));
      lV0tDecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
				    TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2));

      if( myTrackPosTest->GetSign() ==1){
	myTrackPos = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);
	myTrackNeg = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);
	// Daughters' momentum;
	v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
	v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
      }
           
      if( myTrackPosTest->GetSign() ==-1){
	myTrackPos = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);
	myTrackNeg = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);
	// Daughters' momentum;
	v0->GetPPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
	v0->GetNPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      }

      lLabelTrackPos = (UInt_t)TMath::Abs(myTrackPos->GetLabel());
      lLabelTrackNeg = (UInt_t)TMath::Abs(myTrackNeg->GetLabel());

      // Inner Wall parameter:
      const AliExternalTrackParam *myInnerWallTrackPos = myTrackPos->GetInnerParam(); 

      if (myInnerWallTrackPos) {lMomInnerWallPos = myInnerWallTrackPos->GetP();} 
      else {continue;}

      const AliExternalTrackParam *myInnerWallTrackNeg = myTrackNeg->GetInnerParam(); 

      if (myInnerWallTrackNeg) {lMomInnerWallNeg = myInnerWallTrackNeg->GetP();}  
      else {continue;}

      // DCA between daughter and Primary Vertex:
      if (myTrackPos) lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      if (myTrackNeg) lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      // Quality tracks cuts:
      if ( !(fTracksCuts->IsSelected(myTrackPos)) || !(fTracksCuts->IsSelected(myTrackNeg)) ) 
	{	  continue;}

      if ( ( TMath::Abs(myTrackPos->Eta()) > lCutPseudorap ) || ( TMath::Abs(myTrackNeg->Eta()) > lCutPseudorap ) ) {continue;}

      // Armenteros variables:
      lAlphaV0      =  v0->AlphaV0();
      lPtArmV0      =  v0->PtArmV0();

      // Pseudorapidity:
      //      lV0Eta = v0->Eta();
      //////////////////////////////////////////////////////////////////////////
      // Invariant mass
      v0->ChangeMassHypothesis(310);
      lInvMassK0s = v0->GetEffMass();
      lPtK0s = v0->Pt();
      // lPzK0s = v0->Pz();

      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      lPtLambda = v0->Pt();
      //lPzLambda = v0->Pz();

      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();
      lPtAntiLambda = v0->Pt();
      //lPzAntiLambda = v0->Pz();

      // Rapidity:
      lRapK0s    = v0->Y(310);
      lRapLambda = v0->Y(3122);
      lRapAntiLambda = v0->Y(-3122);
	
      if (lPtK0s==0) 	{continue;}
      if (lPtLambda==0) 	{continue;}

      if (lPtAntiLambda==0) 	{continue;}

      ///////////////////////////////////////////////////////////////////////      

      // PID  new method July 2011
      if (fUsePID.Contains("withPID")) {
	//	nSigmaPosPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kPion));
	nSigmaPosPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion));
	//	nSigmaNegPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kPion));
	nSigmaNegPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion));
	//	nSigmaPosProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kProton));
	nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton));
	//	nSigmaNegProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kProton));
        nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton));
	
      }
      else {
	nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;
      }
      
      // Monte-Carlo particle associated to reconstructed particles: 
      if (fAnalysisMC) {
	//if (lLabelTrackPos < 0 || lLabelTrackNeg < 0) continue;
	TParticle  *lMCESDPartPos  = stack->Particle(lLabelTrackPos);
	if(!lMCESDPartPos) { 
	  //  Printf("no MC particle for positive and/or negative daughter\n");
	  continue;
	}
	TParticle  *lMCESDPartNeg  = stack->Particle(lLabelTrackNeg);
	if (!lMCESDPartNeg) 	{  continue;}
	lPDGCodePosDaughter = lMCESDPartPos->GetPdgCode();
	lPDGCodeNegDaughter = lMCESDPartNeg->GetPdgCode();
	lIndexPosMother = lMCESDPartPos->GetFirstMother();
	lIndexNegMother = lMCESDPartNeg->GetFirstMother();

	if (lIndexPosMother == -1) {
	  lPdgcodeMother = 0;
	  lIndexMotherOfMother = 0;
	  mcPosX = 0;
	  mcPosY = 0;
	  mcPosZ = 0;
	  mcPosR = 0;
	  mcPosMotherX = 0;
	  mcPosMotherY = 0;
	  mcPosMotherZ = 0;
	  mcPosMotherR = 0;
	  // mcMotherPt = 1;
	}

	else {

	  TParticle  *lMCESDMother    = stack->Particle(lIndexPosMother);
	  if (!lMCESDMother) 	{ continue;}
	  lPdgcodeMother         = lMCESDMother->GetPdgCode();
	  lIndexMotherOfMother   = lMCESDMother->GetFirstMother();
	  if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
	  else {
	    TParticle  *lMCESDMotherOfMother    = stack->Particle(lIndexMotherOfMother);
	    if (!lMCESDMotherOfMother) 	{ continue;}
	    lPdgcodeMotherOfMother = lMCESDMotherOfMother->GetPdgCode();
	  }
	
	  mcPosX = lMCESDPartPos->Vx();
	  mcPosY = lMCESDPartPos->Vy();
	  mcPosZ = lMCESDPartPos->Vz();
	  mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
	  mcPosMotherX = lMCESDMother->Vx();
	  mcPosMotherY = lMCESDMother->Vy();
	  mcPosMotherZ = lMCESDMother->Vz();
	  mcPosMotherR = TMath::Sqrt(mcPosMotherX*mcPosMotherX+mcPosMotherY*mcPosMotherY);
	
	  //  mcMotherPt   = lMCESDMother->Pt();
	}
      }
    } // end ESD condition
    
    // Multiplicity:
    if(!lOnFlyStatus) nv0s++;
    //    else  if(lOnFlyStatus) nv0sMI++;

    // Daughter momentum cut: ! FIX it in case of AOD !
    //if ( (lPtPos  < cutMinPtDaughter ) ||
    //     (lPtNeg  < cutMinPtDaughter )
    //  ) 	{ continue;}
    
    // Look for associated particles:
    if (fAnalysisMC) {
      if( (lIndexPosMother==-1) || (lIndexNegMother==-1) ) {
	//	fHistMCDaughterTrack->Fill(1);
      }
      
      else if( ( (lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211) ) ) 
	{
	  lCheckPIdK0Short    = 1;
	  //	  fHistMCDaughterTrack->Fill(3);
	  if ( (lIndexPosMother==lIndexNegMother) &&
	       (lPdgcodeMother==310) ) {


	    //if (lIndexPosMother <= lNbMCPrimary) lCheckMcK0Short  = 1;
	    //else lCheckSecondaryK0s = 1;

	    Double_t dx = 0;
	    Double_t dy = 0;
	    Double_t dz = 0;
	    Double_t ProdDistance = 0;

	    dx = ( (mcPrimaryVtx.At(0)) - (mcPosMotherX) ); 
	    dy = ( (mcPrimaryVtx.At(1)) - (mcPosMotherY) );
	    dz = ( (mcPrimaryVtx.At(2)) - (mcPosMotherZ) );

	    ProdDistance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

	    if (ProdDistance < 0.001) lCheckMcK0Short  = 1;
	    else lCheckSecondaryK0s = 1;

	  }
	}
      else if( ( (lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211)  )  ) 
	{
	  lCheckPIdLambda     = 1;
	  //	  fHistMCDaughterTrack->Fill(5);
	  if ( (lIndexPosMother==lIndexNegMother) &&
	       (lPdgcodeMother==3122)  ){
	    if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
		 ) lComeFromSigma = 1;
	    else lComeFromSigma = 0; 


	    // if ( (lIndexPosMother <= lNbMCPrimary) || 
	    //     ( ( lIndexPosMother > lNbMCPrimary) && (lComeFromSigma) )
	    //      ) lCheckMcLambda  = 1; 
	    // else lCheckSecondaryLambda    = 1;

	    Double_t dx = 0;
	    Double_t dy = 0;
	    Double_t dz = 0;
	    Double_t ProdDistance = 0;

	    dx = ( (mcPrimaryVtx.At(0)) - (mcPosMotherX) ); 
	    dy = ( (mcPrimaryVtx.At(1)) - (mcPosMotherY) );
	    dz = ( (mcPrimaryVtx.At(2)) - (mcPosMotherZ) );

	    ProdDistance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

	    if (ProdDistance < 0.001) lCheckMcLambda  = 1;
	    else lCheckSecondaryLambda = 1;
	    
	
	  }
	}
      else if( ( (lPDGCodePosDaughter==211)   && (lPDGCodeNegDaughter==-2212) )	) 
	{
	  lCheckPIdAntiLambda = 1;
	  //	  fHistMCDaughterTrack->Fill(7);
	  if ( (lIndexPosMother==lIndexNegMother) &&
	       (lPdgcodeMother==-3122) ) {
	    if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
		 ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
		 ) lComeFromSigma = 1;
	    else lComeFromSigma = 0;  

	    //  if ( (lIndexPosMother <= lNbMCPrimary) || 
	    //       ( ( lIndexPosMother > lNbMCPrimary) && (lComeFromSigma) )
	    //       ) lCheckMcAntiLambda  = 1;
	    //  else lCheckSecondaryAntiLambda = 1;

	    Double_t dx = 0;
	    Double_t dy = 0;
	    Double_t dz = 0;
	    Double_t ProdDistance = 0;

	    dx = ( (mcPrimaryVtx.At(0)) - (mcPosMotherX) ); 
	    dy = ( (mcPrimaryVtx.At(1)) - (mcPosMotherY) );
	    dz = ( (mcPrimaryVtx.At(2)) - (mcPosMotherZ) );

	    ProdDistance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

	    if (ProdDistance < 0.001) lCheckMcAntiLambda = 1;
	    else lCheckSecondaryAntiLambda = 1;

	  }
	}
      
      // Gamma conversion
      //   else if ( (lPDGCodePosDaughter==-11) &&
      //	(lPDGCodeNegDaughter==11) &&
      //	(lPdgcodeMother==22 ) )
      //	lCheckGamma = 1;

    } // end "look for associated particles  
   
    /////////////////////////////////////     
    // PID condition for daughters tracks
    //////////////////////////////////////

    //    lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
    lCheckPIDLambdaPosDaughter     = 0;//, lCheckPIDLambdaNegDaughter     = 0;
    //lCheckPIDAntiLambdaPosDaughter = 0;,
    lCheckPIDAntiLambdaNegDaughter = 0;

    if (lMomInnerWallPos < lLimitPPID) {
      if (nSigmaPosPion < cutNSigmaLowP)   {
	//	lCheckPIDK0sPosDaughter        = 1;
	//lCheckPIDAntiLambdaPosDaughter = 1;
      }
      if (nSigmaPosProton < cutNSigmaLowP) lCheckPIDLambdaPosDaughter    = 1;      
    }

    else if (lMomInnerWallPos > lLimitPPID) {
      if (nSigmaPosPion < cutNSigmaHighP)   {
	//	lCheckPIDK0sPosDaughter        = 1;
	//	lCheckPIDAntiLambdaPosDaughter = 1;
      }
      if (nSigmaPosProton < cutNSigmaHighP) lCheckPIDLambdaPosDaughter    = 1;
    }

    if (lMomInnerWallNeg < lLimitPPID) {
      if (nSigmaNegPion < cutNSigmaLowP)    {
	//	lCheckPIDK0sNegDaughter       = 1;
	//	lCheckPIDLambdaNegDaughter    = 1;
      }
      if (nSigmaNegProton < cutNSigmaLowP)  lCheckPIDAntiLambdaNegDaughter = 1;
      
    }
    else if (lMomInnerWallNeg > lLimitPPID) {
      if (nSigmaNegPion < cutNSigmaHighP)   {
	//	lCheckPIDK0sNegDaughter       = 1;
	//	lCheckPIDLambdaNegDaughter    = 1;
      }
      if (nSigmaNegProton < cutNSigmaHighP) lCheckPIDAntiLambdaNegDaughter = 1;
    }
 

   
    ///////////////values for cuts/////////////////////////////////////////////////////////////////////////////////////////
    if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) || (lDcaV0Daughters > 1.00) || 
	(lV0cosPointAngle < 0.998) || (lV0Radius < 0.9) || (lV0Radius > 100) ) 

      {continue;}
	

    /////////////////////////////////
    //PID for Lambda and AntiLambda
    /////////////////////////////////
    if(fUsePID.Contains("withPID") && (lCheckPIDAntiLambdaNegDaughter==0) && (lCheckPIDLambdaPosDaughter==1)) LambdaPID = 1;
    else LambdaPID =0;
    if(fUsePID.Contains("withPID") && (lCheckPIDLambdaPosDaughter==0) && (lCheckPIDAntiLambdaNegDaughter==1)) AntiLambdaPID = 1;
    else AntiLambdaPID =0;

    //ctau for lambda
    lcTauLambda     = (lV0tDecayLength*lInvMassLambda)/lPtLambda;
    //ctau for antilambda
    lcTauAntiLambda = (lV0tDecayLength*lInvMassAntiLambda)/lPtAntiLambda;
    //ctau for K0s
    lcTauK0s        = (lV0tDecayLength*lInvMassK0s)/lPtK0s;

    //*****************************
    // filling histograms
    //*****************************
    if (lPLambda <1 && lOnFlyStatus==0 ){
      fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
    }

    ////////////////
    //K0s particle
    ////////////////
    if (TMath::Abs(lRapK0s) < lCutRap ) {
      if (lcTauK0s< cutcTauK0s) {

	if (lOnFlyStatus==0){
          fHistMassK0->Fill(lInvMassK0s);
	  fHistPtVsMassK0->Fill(lInvMassK0s,lPtK0s);
            } //end !lOnFlystatus
          } //end ctau cut
    } //end Rap condition

    ///////////////////
    //Lambda particle    
    ///////////////////
    if ((LambdaPID==1 && lMomInnerWallPos <=1 ) || (lMomInnerWallPos > 1) ||  !(fUsePID.Contains("withPID")  )){  
      if (TMath::Abs(lRapLambda) < lCutRap) {
	if (lcTauLambda < cutcTauL){
	  if (lOnFlyStatus==0){
	    fHistMassLambda->Fill(lInvMassLambda);
	    fHistPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
	    if(lPtLambda <=1) fHistNSigmaProton->Fill(nSigmaPosProton);
	  }
	}// end ctau condition
      } //end of Rap condition
    }// end of PID condition

    //////////////////////////////
    // Anti Lambda ///////////////
    //////////////////////////////
    if ((AntiLambdaPID==1 && lMomInnerWallNeg <=1) || (lMomInnerWallNeg>1) ||  !(fUsePID.Contains("withPID"))){  
      if (TMath::Abs(lRapAntiLambda) < lCutRap) {
	if (lcTauAntiLambda < cutcTauL){
	  if (lOnFlyStatus==0){
	    fHistMassAntiLambda->Fill(lInvMassAntiLambda);
	    fHistPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
	  }
	} //end of Rap condition
      } // end of PID condition
    } //end ctau cut




    // Histo versus Rap and armenteros plot
    if (!lOnFlyStatus){
      if (lCheckMcK0Short) fHistAsMcRapK0->Fill(lRapK0s);
      if (lCheckMcLambda) fHistAsMcRapLambda->Fill(lRapLambda);
      if (lCheckMcAntiLambda) fHistAsMcRapLambda->Fill(lRapAntiLambda);
      //      fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
      if ((TMath::Abs(lRapK0s) < lCutRap)&&(TMath::Abs(lRapLambda) < lCutRap)) fHistK0sMassVsLambdaMass->Fill(lInvMassK0s,lInvMassLambda);
    }


    ///////////////////////////////////////////////////    
    // K0s associated histograms in |rap| < lCutRap:
    ///////////////////////////////////////////////////
    if (TMath::Abs(lRapK0s) < lCutRap) {
      switch (lOnFlyStatus){
      case 0 : 
	if (lcTauK0s< cutcTauK0s) {
	  if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(lInvMassK0s);
	  if(lCheckMcK0Short) {
	    fHistAsMcMassK0->Fill(lInvMassK0s);
	    fHistAsMcPtK0->Fill(lPtK0s);
	    fHistAsMcPtVsMassK0->Fill(lInvMassK0s,lPtK0s);
	  }
	}
      } // end rapidity condition
    } // end ctau cut
    
    ///////////////////////////////////////////////////
    // Associated Lambda histograms in |rap| < lCutRap
    ////////////////////////////////////////////////////
     if ((LambdaPID==1 && lMomInnerWallPos <=1) || (lMomInnerWallPos>1) ||  !(fUsePID.Contains("withPID"))){  
    if (TMath::Abs(lRapLambda) < lCutRap) {
      switch (lOnFlyStatus){
      case 0 : 
	if (lcTauLambda < cutcTauL){    
	  if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(lInvMassLambda);
	  if(lCheckMcLambda) {
	    fHistAsMcMassLambda->Fill(lInvMassLambda);
	    fHistAsMcPtLambda->Fill(lPtLambda);
	    //	  fHistCosPointAngleLvsMassVsPtsigL->Fill(lPtLambda,lV0cosPointAngle,lInvMassLambda);
	    fHistAsMcPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
	  }
	}
      } // end rapidity condition
     }// end PID condition
    }// end ctau condition

    ////////////////////////////////////////////////////////
    // Associated AntiLambda histograms in |rap| < lCutRap
    ////////////////////////////////////////////////////////
      if ((AntiLambdaPID==1 && lMomInnerWallNeg <=1) || (lMomInnerWallNeg>1) ||  !(fUsePID.Contains("withPID"))){          
    if (TMath::Abs(lRapAntiLambda) < lCutRap) {
      switch (lOnFlyStatus){
      case 0 : 
	if (lcTauAntiLambda < cutcTauL){
	  if(lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(lInvMassAntiLambda);
	  if(lCheckMcAntiLambda) {
	    fHistAsMcMassAntiLambda->Fill(lInvMassAntiLambda);
	    fHistAsMcPtAntiLambda->Fill(lPtAntiLambda);
	    fHistAsMcPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
	  }
	}
      } // end rapidity condition
     }// end PID condition     
    }// end ctau cut
  
  } // end V0 loop

  fHistV0Multiplicity->Fill(nv0s);
  //  fHistV0MultiplicityMI->Fill(nv0sMI);

  if (fAnalysisType == "ESD") { if(myPrimaryVertex) delete myPrimaryVertex; }
  //  if (fAnalysisType == "ESD") { if(TestTrackCuts) delete TestTrackCuts; }

  
  // Post output data
}      

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::Terminate(Option_t *) 
{/*
 // Draw result to the screen
 // Called once at the end of the query

 TList *cRetrievedList = 0x0;
 cRetrievedList = (TList*)GetOutputData(1);
  
 if(!cRetrievedList){
 AliWarning("ERROR - AliAnalysisTaskPerformanceStrange: output data container list not available\n"); return;
 }
  
  
 fHistV0Multiplicity = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistV0Multiplicity"));
 if (!fHistV0Multiplicity) {
 Printf("ERROR: fHistV0Multiplicity not available");
 return;
 }

 fHistV0MultiplicityMI = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistV0MultiplicityMI"));
 if (!fHistV0MultiplicityMI) {
 Printf("ERROR: fHistV0MultiplicityMI not available");
 return;
 }

 TCanvas *canPerformanceStrange = new TCanvas("AliAnalysisTaskCheckV0","Multiplicity",10,10,510,510);
 canPerformanceStrange->Divide(2,1);
 if (fHistV0Multiplicity->GetMaximum() > 0.) canPerformanceStrange->cd(1)->SetLogy();
 fHistV0Multiplicity->SetMarkerStyle(25);
 fHistV0Multiplicity->DrawCopy("E");
 if (fHistV0MultiplicityMI->GetMaximum() > 0.) canPerformanceStrange->cd(2)->SetLogy();
 fHistV0MultiplicityMI->SetMarkerStyle(24);
 fHistV0MultiplicityMI->DrawCopy("E");
  

 */ 
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskPerformanceStrange::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 
//----------------------------------------------------------------------------

