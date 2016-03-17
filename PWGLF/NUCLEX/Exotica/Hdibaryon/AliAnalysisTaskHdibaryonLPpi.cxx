/**************************************************************************
 * Author : Benjamin Dönigus (benjamin.doenigus@cern.ch)                  *
 *                                                                        *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskHdibaryonLPpi class
//          used to search for the H-Dibaryon in weak 
//          (Lambda Proton Pion) and strong (Lambda Lambda) decays
//-----------------------------------------------------------------

#include "Riostream.h"
#include "TROOT.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TNtuple.h"
#include "TObjString.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliESD.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliGenEventHeader.h"
#include "AliESDInputHandler.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliHeader.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0Cuts.h"
#include "AliAnalysisTaskHdibaryonLPpi.h"

#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TList.h"
#include "TParallelCoord.h"

#include "AliMCParticle.h"
#include "AliGenPythiaEventHeader.h"

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"

#include "THnSparse.h"

#include "AliVertexerTracks.h"

using namespace std;
ClassImp(AliAnalysisTaskHdibaryonLPpi)
//________________________________________________________________________
AliAnalysisTaskHdibaryonLPpi::AliAnalysisTaskHdibaryonLPpi() : AliAnalysisTaskSE()/*AliAnalysisTask(name, ""), fMCEvent(0)*/, fESD(0),   fESDtrackCutsV0(0),
  fESDCutsV0(0),
  fEsdTrackCuts(0),
  fBin(0),
  fEvent(0x0),
  fHistList(0),  
  fHistMassDPi(0), 
  fHistMassLPi(0),
  fHistMassLambdaPi(0),
  fHistMassLambda(0),
  fHistMassLambdaPPi(0),
  fHistMassLambdaP(0),
  fHistMassLambdaK(0),
  fHistMassK0onl(0),
  fHistMassK0offl(0),
  fHistMassK0onlC(0),
  fHistMassK0offlC(0),
  fHistMassPQonl(0), 
  fHistMassPQoffl(0),
  fHistDC(0),
  fHistArmenterosPodolanski(0),
  fHistArmenterosPodolanskiCut(0), 
  fHistHDibaryonInvaMassGen(0),
  fHistHDibaryonInvaMassGenRes(0),
  fHistAntiHDibaryonInvaMassGen(0),
  fHistHDibaryonInvaMassAso(0),
  fHistHDibaryonInvaMassAsoReso(0),
  fHistAntiHDibaryonInvaMassAso(0),
  fHistCheck(0),
  fHistHPointingAngle(0),
  fHistMassH(0),
  fHistMassLambdaFromH(0),
  fHistMassLambdaFromHtLorentz(0),
  fHistMassPpi(0),
  fHistMassPpiReso(0),
  fHistMassLpi(0),
  fHistMassLP(0),
  fHistProtonPIDBb(0),
  fHistPionPIDBb(0),
  fHistProtonPIDLambda(0),
  fHistPionPIDLambda(0),
  fHistMCdcaPvtxDvtx(0),
  fHistMCdcaPvtxLvtx(0),
  fHistMCdcaDvtxLvtx(0),
  fHistMCangleLH(0),
  fHistMCdecayAngle(0),
  fHistMCpointingAngle(0),
  fHistMCap(0),
  fHistMCdcaPvtxDvtxReso(0),
  fHistMCdcaPvtxLvtxReso(0),
  fHistMCdcaDvtxLvtxReso(0),
  fHistMCangleLHReso(0),
  fHistMCdecayAngleReso(0),
  fHistMCpointingAngleReso(0),
  fHistMCapReso(0),
  fHistCentrality(0),
  fHistCentralityAC(0), 
  fHistMultiplicity(0),
  fHistHilf1(0),
  fHistHilf2(0), 
  fHistHilf3(0),
  fHistHilf4(0),
  fHistHilf5(0), 
  fHistHilf6(0),
  fHistPtvsEtaGen(0),
  fHistPtvsEtaAso(0),
  fHistPtvsYGen(0),
  fHistPtvsYAso(0), 
  fHistRap(0),
  fHistCount(0),
  fPIDtpcESD(0),
  fHistTriggerStat(0),
  fHistTriggerStatAfterEventSelection(0), 
  fHistMassHcentMult(0),
  fHistNdim(0)
{
  // DefaultConstructor

}

//________________________________________________________________________
AliAnalysisTaskHdibaryonLPpi::AliAnalysisTaskHdibaryonLPpi(const char *name) : AliAnalysisTaskSE(name)/*AliAnalysisTask(name, ""), fMCEvent(0)*/, fESD(0),   fESDtrackCutsV0(0),
  fESDCutsV0(0),
  fEsdTrackCuts(0),
  fBin(0),
  fEvent(0x0),
  fHistList(0),  
  fHistMassDPi(0), 
  fHistMassLPi(0),
  fHistMassLambdaPi(0),
  fHistMassLambda(0),
  fHistMassLambdaPPi(0),
  fHistMassLambdaP(0),
  fHistMassLambdaK(0),
  fHistMassK0onl(0),
  fHistMassK0offl(0),
  fHistMassK0onlC(0),
  fHistMassK0offlC(0),
  fHistMassPQonl(0), 
  fHistMassPQoffl(0),
  fHistDC(0),
  fHistArmenterosPodolanski(0),
  fHistArmenterosPodolanskiCut(0), 
  fHistHDibaryonInvaMassGen(0),
  fHistHDibaryonInvaMassGenRes(0),
  fHistAntiHDibaryonInvaMassGen(0),
  fHistHDibaryonInvaMassAso(0),
  fHistHDibaryonInvaMassAsoReso(0),
  fHistAntiHDibaryonInvaMassAso(0),
  fHistCheck(0),
  fHistHPointingAngle(0),
  fHistMassH(0),
  fHistMassLambdaFromH(0),
  fHistMassLambdaFromHtLorentz(0),
  fHistMassPpi(0),
  fHistMassPpiReso(0),
  fHistMassLpi(0),
  fHistMassLP(0),
  fHistProtonPIDBb(0),
  fHistPionPIDBb(0),
  fHistProtonPIDLambda(0),
  fHistPionPIDLambda(0),
  fHistMCdcaPvtxDvtx(0),
  fHistMCdcaPvtxLvtx(0),
  fHistMCdcaDvtxLvtx(0),
  fHistMCangleLH(0),
  fHistMCdecayAngle(0),
  fHistMCpointingAngle(0),
  fHistMCap(0),
  fHistMCdcaPvtxDvtxReso(0),
  fHistMCdcaPvtxLvtxReso(0),
  fHistMCdcaDvtxLvtxReso(0),
  fHistMCangleLHReso(0),
  fHistMCdecayAngleReso(0),
  fHistMCpointingAngleReso(0),
  fHistMCapReso(0),
  fHistCentrality(0),
  fHistCentralityAC(0), 
  fHistMultiplicity(0), 
  fHistHilf1(0),
  fHistHilf2(0), 
  fHistHilf3(0),
  fHistHilf4(0),
  fHistHilf5(0), 
  fHistHilf6(0),
  fHistPtvsEtaGen(0),
  fHistPtvsEtaAso(0),
  fHistPtvsYGen(0),
  fHistPtvsYAso(0), 
  fHistRap(0),
  fHistCount(0),
  fPIDtpcESD(0),
  fHistTriggerStat(0),
  fHistTriggerStatAfterEventSelection(0),
  fHistMassHcentMult(0),
  fHistNdim(0)

{
  // Constructor

  // Define input and output slots here
  // Input from a TChain
  DefineInput(0, TChain::Class());
  // Output to TList container
  DefineOutput(1, TList::Class()); //full

  //MC info contol
  if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())SetHasMC();

  //V0 cuts

  fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsV0->SetMinNClustersTPC(80);
  fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
  fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
  fESDtrackCutsV0->SetPtRange(0.2,1.5);
  fESDtrackCutsV0->SetMinDCAToVertexXY(2); //war inzwischen 1 & 3
  fESDtrackCutsV0->SetMinDCAToVertexZ(2); //war inzwischen 1 & 3

  fESDCutsV0 = new AliESDv0Cuts("AliESDCutsV0","AliESDCutsV0");
  fESDCutsV0->SetMaxDcaV0Daughters(1.0);
  fESDCutsV0->SetMinDcaNegToVertex(2); //1.5
    fESDCutsV0->SetMinDcaPosToVertex(2); //1.5

  //ESD Track cuts
  fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");    
  fEsdTrackCuts->SetMinNClustersTPC(80);
  fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fEsdTrackCuts->SetMaxChi2PerClusterTPC(5);
  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
  fEsdTrackCuts->SetEtaRange(-0.9,0.9);
}

//____________________________________________________________
AliAnalysisTaskHdibaryonLPpi::~AliAnalysisTaskHdibaryonLPpi(){
  //
  // Destructor
  //
  if(fHistList){ 
    fHistList->Clear();
    delete fHistList;
  }
  if(fEsdTrackCuts) delete fEsdTrackCuts;
  if(fESDtrackCutsV0) delete fESDtrackCutsV0;
  if(fESDCutsV0) delete fESDCutsV0;

}

//________________________________________________________________________
void AliAnalysisTaskHdibaryonLPpi::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

 fHistList = new TList();
 fHistList->SetOwner();

 fHistMassDPi = new TH1F("fHistMassDPi", "Invariant mass distribution p+#pi^{-} ", 500, 1.0, 1.25);
 fHistMassDPi->GetXaxis()->SetTitle("Invariant mass p+#pi^{-} (GeV/c^{2})");
 fHistMassDPi->GetYaxis()->SetTitle("Entries");
 fHistMassDPi->SetMarkerStyle(kFullCircle); 

 fHistMassLPi = new TH1F("fHistMassLPi", "Offline Invariant mass distribution p+#pi^{-} ", 500, 1.0, 1.25);
 fHistMassLPi->GetXaxis()->SetTitle("Invariant mass p+#pi^{-} (GeV/c^{2})");
 fHistMassLPi->GetYaxis()->SetTitle("Entries");
 fHistMassLPi->SetMarkerStyle(kFullCircle);

 fHistMassLambdaPi = new TH1F("fHistMassLambdaPi", "Invariant mass distribution #Lambda+#pi^{-} ", 500, 1.2, 1.5);
 fHistMassLambdaPi->GetXaxis()->SetTitle("Invariant mass #Lambda+#pi^{-} (GeV/c^{2})");
 fHistMassLambdaPi->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaPi->SetMarkerStyle(kFullCircle);
 
 fHistMassLambda = new TH1F("fHistMassLambda", "Invariant mass distribution of #Lambda for further analyis", 500, 1.0, 1.2);
 fHistMassLambda->GetXaxis()->SetTitle("Invariant mass p+#pi^{+} (GeV/c^{2})");
 fHistMassLambda->GetYaxis()->SetTitle("Entries");
 fHistMassLambda->SetMarkerStyle(kFullCircle);
 
 fHistMassLambdaPPi = new TH1F("fHistMassLambdaPPi", "Invariant mass distribution #Lambdap#pi^{-} ", 300, 2.2, 2.5);
 fHistMassLambdaPPi->GetXaxis()->SetTitle("Invariant mass #Lambdap#pi^{-} (GeV/c^{2})");
 fHistMassLambdaPPi->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaPPi->SetMarkerStyle(kFullCircle);

 fHistMassLambdaP = new TH1F("fHistMassLambdaP", "Invariant mass distribution #Lambdap ", 300, 2.2, 2.5);
 fHistMassLambdaP->GetXaxis()->SetTitle("Invariant mass #Lambdap (GeV/c^{2})");
 fHistMassLambdaP->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaP->SetMarkerStyle(kFullCircle);

 fHistMassLambdaK = new TH1F("fHistMassLambdaK", "Invariant mass distribution #LambdaK^{-} ", 300, 1.6, 1.9);
 fHistMassLambdaK->GetXaxis()->SetTitle("Invariant mass #LambdaK^{-} (GeV/c^{2})");
 fHistMassLambdaK->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaK->SetMarkerStyle(kFullCircle);

 fHistMassK0onl = new TH1F("fHistMassK0onl", "Invariant mass distribution K_{s}^{0} online V0 finder ", 400, 0.2, 1.0);
 fHistMassK0onl->GetXaxis()->SetTitle("Invariant mass #pi^{+}#pi^{-} (GeV/c^{2})");
 fHistMassK0onl->GetYaxis()->SetTitle("Entries");
 fHistMassK0onl->SetMarkerStyle(kFullCircle);

 fHistMassK0offl = new TH1F("fHistMassK0offl", "Invariant mass distribution  K_{s}^{0} offline V0 finder ", 400, 0.2, 1.0);
 fHistMassK0offl->GetXaxis()->SetTitle("Invariant mass  #pi^{+}#pi^{-} (GeV/c^{2})");
 fHistMassK0offl->GetYaxis()->SetTitle("Entries");
 fHistMassK0offl->SetMarkerStyle(kFullCircle);
 
 fHistMassK0onlC = new TH1F("fHistMassK0onlC", "Invariant mass distribution K_{s}^{0} online V0 finder ", 400, 0.2, 1.0);
 fHistMassK0onlC->GetXaxis()->SetTitle("Invariant mass #pi^{+}#pi^{-} (GeV/c^{2})");
 fHistMassK0onlC->GetYaxis()->SetTitle("Entries");
 fHistMassK0onlC->SetMarkerStyle(kFullCircle);

 fHistMassK0offlC = new TH1F("fHistMassK0offlC", "Invariant mass distribution  K_{s}^{0} offline V0 finder ", 400, 0.2, 1.0);
 fHistMassK0offlC->GetXaxis()->SetTitle("Invariant mass  #pi^{+}#pi^{-} (GeV/c^{2})");
 fHistMassK0offlC->GetYaxis()->SetTitle("Entries");
 fHistMassK0offlC->SetMarkerStyle(kFullCircle);

 fHistMassPQonl = new TH1F("fHistMassPQonl", "Invariant mass distribution K_{s}^{0}p using online V0 finder ", 500, 1.3, 2.3);
 fHistMassPQonl->GetXaxis()->SetTitle("Invariant mass K_{s}^{0}p (GeV/c^{2})");
 fHistMassPQonl->GetYaxis()->SetTitle("Entries");
 fHistMassPQonl->SetMarkerStyle(kFullCircle);

 fHistMassPQoffl = new TH1F("fHistMassPQoffl", "Invariant mass distribution  K_{s}^{0}p using offline V0 finder ", 500, 1.3, 2.3);
 fHistMassPQoffl->GetXaxis()->SetTitle("Invariant mass K_{s}^{0}p (GeV/c^{2})");
 fHistMassPQoffl->GetYaxis()->SetTitle("Entries");
 fHistMassPQoffl->SetMarkerStyle(kFullCircle);

 fHistDC = new TH1F("fHistDC", "Proper decay length", 500, 0.0, 25);
 fHistDC->GetXaxis()->SetTitle("c#tau (cm)");
 fHistDC->GetYaxis()->SetTitle("Entries");
 fHistDC->SetMarkerStyle(kFullCircle);
  
 fHistArmenterosPodolanski = new TH2F("fHistArmenterosPodolanski", "Armenteros-Podolanski plot", 200,-1.0,1.0, 500,0,1);
 fHistArmenterosPodolanski->GetXaxis()->SetTitle("#alpha");
 fHistArmenterosPodolanski->GetYaxis()->SetTitle("q_{t}");
 fHistArmenterosPodolanski->SetMarkerStyle(kFullCircle);

 fHistArmenterosPodolanskiCut = new TH2F("fHistArmenterosPodolanskiCut", "Armenteros-Podolanski plot after cut", 200,-1.0,1.0, 500,0,1);
 fHistArmenterosPodolanskiCut->GetXaxis()->SetTitle("#alpha");
 fHistArmenterosPodolanskiCut->GetYaxis()->SetTitle("q_{t}");
 fHistArmenterosPodolanskiCut->SetMarkerStyle(kFullCircle);

 fHistHDibaryonInvaMassGen = new TH1F("fHistHDibaryonInvaMassGen", "Generated  #Lambda p #pi^{-}", 200, 2.1, 2.3);
 fHistHDibaryonInvaMassGen->GetYaxis()->SetTitle("Counts");
 fHistHDibaryonInvaMassGen->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistHDibaryonInvaMassGenRes = new TH1F("fHistHDibaryonInvaMassGenRes", "Generated  #Lambda p #pi^{-} with particles in rapidity!!!", 200, 2.1, 2.3);
 fHistHDibaryonInvaMassGenRes->GetYaxis()->SetTitle("Counts");
 fHistHDibaryonInvaMassGenRes->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistAntiHDibaryonInvaMassGen = new TH1F("fHistAntiHDibaryonInvaMassGen", "Generated  #bar{#Lambda} #bar{p} #pi^{+}", 200, 2.1, 2.3);
 fHistAntiHDibaryonInvaMassGen->GetYaxis()->SetTitle("Counts");
 fHistAntiHDibaryonInvaMassGen->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistHDibaryonInvaMassAso = new TH1F("fHistHDibaryonInvaMassAso", "Associated  #Lambda p #pi^{-}", 200, 2.1, 2.3);
 fHistHDibaryonInvaMassAso->GetYaxis()->SetTitle("Counts");
 fHistHDibaryonInvaMassAso->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistHDibaryonInvaMassAsoReso = new TH1F("fHistHDibaryonInvaMassAsoReso", "Associated  #Lambda p #pi^{-}", 200, 2.1, 2.3);
 fHistHDibaryonInvaMassAsoReso->GetYaxis()->SetTitle("Counts");
 fHistHDibaryonInvaMassAsoReso->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistAntiHDibaryonInvaMassAso = new TH1F("fHistAntiHDibaryonInvaMassAso", "Associated  #bar{#Lambda} #bar{p} #pi^{+}", 200, 2.1, 2.3);
 fHistAntiHDibaryonInvaMassAso->GetYaxis()->SetTitle("Counts");
 fHistAntiHDibaryonInvaMassAso->GetXaxis()->SetTitle("Invariant mass  (GeV/c^{2})");

 fHistCheck = new TH2F("fHistCheck", "Check online/offline", 200, -0.5, 1.5, 200, -0.5, 1.5);
 fHistCheck->GetXaxis()->SetTitle("offline");
 fHistCheck->GetYaxis()->SetTitle("online");
 fHistCheck->SetMarkerStyle(kFullCircle);

 fHistHPointingAngle= new TH1F("fHistHPointingAngle", "Pointing angle distribution for #Lambdap#pi^{-}", 200, 0., 2*TMath::Pi());
 fHistHPointingAngle->GetXaxis()->SetTitle("Pointing angle distribution for #Lambdap#pi^{-}");
 fHistHPointingAngle->GetYaxis()->SetTitle("Entries");
 fHistHPointingAngle->SetMarkerStyle(kFullCircle);

 fHistMassH= new TH1F("fHistMassH", "Invariant mass distribution #Lambdap#pi^{-} (GeV/c^{2}) after pointing angle cut", 3000, 2.2, 2.5);
 fHistMassH->GetXaxis()->SetTitle("Invariant mass #Lambdap#pi^{-} (GeV/c^{2})");
 fHistMassH->GetYaxis()->SetTitle("Entries");
 fHistMassH->SetMarkerStyle(kFullCircle);
 
 fHistMassLambdaFromH= new TH1F("fHistMassLambdaFromH", "Invariant mass distribution #Lambda (GeV/c^{2}) asking for Mother to be H", 300, 1.0, 1.3);
 fHistMassLambdaFromH->GetXaxis()->SetTitle("Invariant mass #Lambda (GeV/c^{2})");
 fHistMassLambdaFromH->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaFromH->SetMarkerStyle(kFullCircle);

 fHistMassLambdaFromHtLorentz= new TH1F("fHistMassLambdaFromHtLorentz", "Invariant mass distribution #Lambda (GeV/c^{2}) asking for Mother to be H", 300, 1.0, 1.3);
 fHistMassLambdaFromHtLorentz->GetXaxis()->SetTitle("Invariant mass #Lambda (GeV/c^{2})");
 fHistMassLambdaFromHtLorentz->GetYaxis()->SetTitle("Entries");
 fHistMassLambdaFromHtLorentz->SetMarkerStyle(kFullCircle);
 
 fHistMassPpi= new TH1F("fHistMassPpi", "Invariant mass distribution of the p#pi^{-} used for combing with Lambda", 300, 1.0, 1.6);
 fHistMassPpi->GetXaxis()->SetTitle("Invariant mass p#pi^{-} (GeV/c^{2})");
 fHistMassPpi->GetYaxis()->SetTitle("Entries");
 fHistMassPpi->SetMarkerStyle(kFullCircle);

 fHistMassPpiReso= new TH1F("fHistMassPpiReso", "Invariant mass distribution of the p#pi^{-} used for combing with Lambda", 300, 1.0, 1.6);
 fHistMassPpiReso->GetXaxis()->SetTitle("Invariant mass p#pi^{-} (GeV/c^{2})");
 fHistMassPpiReso->GetYaxis()->SetTitle("Entries");
 fHistMassPpiReso->SetMarkerStyle(kFullCircle);

 fHistMassLpi= new TH1F("fHistMassLpi", "Invariant mass distribution of the #Lambda#pi^{-} used for combing with p", 300, 1.1, 1.7);
 fHistMassLpi->GetXaxis()->SetTitle("Invariant mass #Lambda#pi^{-} (GeV/c^{2})");
 fHistMassLpi->GetYaxis()->SetTitle("Entries");
 fHistMassLpi->SetMarkerStyle(kFullCircle);

 fHistMassLP= new TH1F("fHistMassLP", "Invariant mass distribution of the #Lambda p used for combing with #pi^{-}", 300, 2.0, 2.3);
 fHistMassLP->GetXaxis()->SetTitle("Invariant mass #Lambda p (GeV/c^{2})");
 fHistMassLP->GetYaxis()->SetTitle("Entries");
 fHistMassLP->SetMarkerStyle(kFullCircle);

 fHistProtonPIDBb = new TH2F("fHistProtonPIDBb", "dE/dx after p PID", 100, 0., 10, 100, 0, 100);
 fHistProtonPIDBb->GetYaxis()->SetTitle("TPC Signal");
 fHistProtonPIDBb->GetXaxis()->SetTitle("P (GeV/c)");
 fHistProtonPIDBb->SetOption("scat");
 fHistProtonPIDBb->SetMarkerStyle(kFullCircle);

 fHistPionPIDBb = new TH2F("fHistPionPIDBb", "dE/dx after K PID", 100, 0., 10, 100, 0, 100);
 fHistPionPIDBb->GetYaxis()->SetTitle("TPC Signal");
 fHistPionPIDBb->GetXaxis()->SetTitle("P (GeV/c)");
 fHistPionPIDBb->SetOption("scat");
 fHistPionPIDBb->SetMarkerStyle(kFullCircle);

 fHistProtonPIDLambda = new TH2F("fHistProtonPIDLambda", "dE/dx after p PID from V0", 100, 0., 10, 100, 0, 100);
 fHistProtonPIDLambda->GetYaxis()->SetTitle("TPC Signal");
 fHistProtonPIDLambda->GetXaxis()->SetTitle("P (GeV/c)");
 fHistProtonPIDLambda->SetOption("scat");
 fHistProtonPIDLambda->SetMarkerStyle(kFullCircle);

 fHistPionPIDLambda = new TH2F("fHistPionPIDLambda", "dE/dx after #pi PID from V0", 100, 0, 10, 100, 0, 100);
 fHistPionPIDLambda->GetYaxis()->SetTitle("TPC Signal");
 fHistPionPIDLambda->GetXaxis()->SetTitle("P (GeV/c)");
 fHistPionPIDLambda->SetOption("scat");
 fHistPionPIDLambda->SetMarkerStyle(kFullCircle);

 fHistMCdcaPvtxDvtx= new TH1F("fHistMCdcaPvtxDvtx", "MC True DCA Primary Vertex - H Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaPvtxDvtx->GetXaxis()->SetTitle("dca prim. vtx- decay vtx (cm)");
 fHistMCdcaPvtxDvtx->GetYaxis()->SetTitle("Entries");
 fHistMCdcaPvtxDvtx->SetMarkerStyle(kFullCircle);

 fHistMCdcaPvtxLvtx= new TH1F("fHistMCdcaPvtxLvtx", "MC True DCA Primary Vertex - Lambda Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaPvtxLvtx->GetXaxis()->SetTitle("dca prim. vtx-#Lambda decay vtx (cm)");
 fHistMCdcaPvtxLvtx->GetYaxis()->SetTitle("Entries");
 fHistMCdcaPvtxLvtx->SetMarkerStyle(kFullCircle);

 fHistMCdcaDvtxLvtx= new TH1F("fHistMCdcaDvtxLvtx", "MC True DCA H Decay Vertex - Lambda Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaDvtxLvtx->GetXaxis()->SetTitle("dca H deacy vtx-#Lambda decay vtx (cm)");
 fHistMCdcaDvtxLvtx->GetYaxis()->SetTitle("Entries");
 fHistMCdcaDvtxLvtx->SetMarkerStyle(kFullCircle);

 fHistMCangleLH= new TH1F("fHistMCangleLH", "MC True Angle between #Lambda and H", 300, 0., 2*TMath::Pi());
 fHistMCangleLH->GetXaxis()->SetTitle("Angle (#Lambda-H)");
 fHistMCangleLH->GetYaxis()->SetTitle("Entries");
 fHistMCangleLH->SetMarkerStyle(kFullCircle);

 fHistMCdecayAngle= new TH1F("fHistMCdecayAngle", "MC True Angle between decay products", 300, 0., 2*TMath::Pi());
 fHistMCdecayAngle->GetXaxis()->SetTitle("Angle (#Lambda-p#pi)");
 fHistMCdecayAngle->GetYaxis()->SetTitle("Entries");
 fHistMCdecayAngle->SetMarkerStyle(kFullCircle);

 fHistMCpointingAngle= new TH1F("fHistMCpointingAngle", "MC True Pointing Angle", 3000, 0., 2*TMath::Pi());
 fHistMCpointingAngle->GetXaxis()->SetTitle("Pointing Angle");
 fHistMCpointingAngle->GetYaxis()->SetTitle("Entries");
 fHistMCpointingAngle->SetMarkerStyle(kFullCircle);

 fHistMCap = new TH2F("fHistMCap", "True MC Armenteros-Podolanski", 200,-1.0,1.0, 500,0,1);
 fHistMCap->GetYaxis()->SetTitle("#alpha");
 fHistMCap->GetXaxis()->SetTitle("q_{t}");
 fHistMCap->SetOption("scat");
 fHistMCap->SetMarkerStyle(kFullCircle);

 fHistMCdcaPvtxDvtxReso= new TH1F("fHistMCdcaPvtxDvtxReso", "MC True DCA Primary Vertex - H Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaPvtxDvtxReso->GetXaxis()->SetTitle("dca prim. vtx- decay vtx (cm)");
 fHistMCdcaPvtxDvtxReso->GetYaxis()->SetTitle("Entries");
 fHistMCdcaPvtxDvtxReso->SetMarkerStyle(kFullCircle);

 fHistMCdcaPvtxLvtxReso= new TH1F("fHistMCdcaPvtxLvtxReso", "MC True DCA Primary Vertex - Lambda Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaPvtxLvtxReso->GetXaxis()->SetTitle("dca prim. vtx-#Lambda decay vtx (cm)");
 fHistMCdcaPvtxLvtxReso->GetYaxis()->SetTitle("Entries");
 fHistMCdcaPvtxLvtxReso->SetMarkerStyle(kFullCircle);

 fHistMCdcaDvtxLvtxReso= new TH1F("fHistMCdcaDvtxLvtxReso", "MC True DCA H Decay Vertex - Lambda Decay Vertex", 300, -0.1, 11.9);
 fHistMCdcaDvtxLvtxReso->GetXaxis()->SetTitle("dca H deacy vtx-#Lambda decay vtx (cm)");
 fHistMCdcaDvtxLvtxReso->GetYaxis()->SetTitle("Entries");
 fHistMCdcaDvtxLvtxReso->SetMarkerStyle(kFullCircle);

 fHistMCangleLHReso= new TH1F("fHistMCangleLHReso", "MC True Angle between #Lambda and H", 300, 0., 2*TMath::Pi());
 fHistMCangleLHReso->GetXaxis()->SetTitle("Angle (#Lambda-H)");
 fHistMCangleLHReso->GetYaxis()->SetTitle("Entries");
 fHistMCangleLHReso->SetMarkerStyle(kFullCircle);

 fHistMCdecayAngleReso= new TH1F("fHistMCdecayAngleReso", "MC True Angle between decay products", 300, 0., 2*TMath::Pi());
 fHistMCdecayAngleReso->GetXaxis()->SetTitle("Angle (#Lambda-p#pi)");
 fHistMCdecayAngleReso->GetYaxis()->SetTitle("Entries");
 fHistMCdecayAngleReso->SetMarkerStyle(kFullCircle);

 fHistMCpointingAngleReso= new TH1F("fHistMCpointingAngleReso", "MC True Pointing Angle", 300, 0., 2*TMath::Pi());
 fHistMCpointingAngleReso->GetXaxis()->SetTitle("Pointing Angle");
 fHistMCpointingAngleReso->GetYaxis()->SetTitle("Entries");
 fHistMCpointingAngleReso->SetMarkerStyle(kFullCircle);

 fHistMCapReso = new TH2F("fHistMCapReso", "True MC Armenteros-Podolanski", 200,-1.0,1.0, 500,0,1);
 fHistMCapReso->GetYaxis()->SetTitle("#alpha");
 fHistMCapReso->GetXaxis()->SetTitle("q_{t}");
 fHistMCapReso->SetOption("scat");
 fHistMCapReso->SetMarkerStyle(kFullCircle);

 fHistCentrality = new TH1F("Centrality ", "Centrality", 11, -0.5, 10.5);
 fHistCentrality ->GetXaxis()->SetTitle("Centrality");
 fHistCentrality ->GetYaxis()->SetTitle("Entries");

 fHistCentralityAC = new TH1F("CentralityAC ", "CentralityAC", 11, -0.5, 10.5);
 fHistCentralityAC ->GetXaxis()->SetTitle("Centrality");
 fHistCentralityAC ->GetYaxis()->SetTitle("Entries");

 fHistMultiplicity = new TH1F("Multiplicity ", "Multiplicity", 100, 0, 20000);
 fHistMultiplicity ->GetXaxis()->SetTitle("Centrality");
 fHistMultiplicity ->GetYaxis()->SetTitle("Entries");

 fHistList->Add(fHistMassDPi);
 fHistList->Add(fHistMassLPi);
 fHistList->Add(fHistMassLambdaPi);
 fHistList->Add(fHistMassLambda);
 fHistList->Add(fHistMassLambdaPPi);
 fHistList->Add(fHistMassLambdaP);
 fHistList->Add(fHistMassLambdaK);
 fHistList->Add(fHistMassK0onl);
 fHistList->Add(fHistMassK0offl);
 fHistList->Add(fHistMassK0onlC);
 fHistList->Add(fHistMassK0offlC);
 fHistList->Add(fHistMassPQonl); 
 fHistList->Add(fHistMassPQoffl);
 fHistList->Add(fHistDC); 
 fHistList->Add(fHistArmenterosPodolanski);
 fHistList->Add(fHistArmenterosPodolanskiCut);
 fHistList->Add(fHistHDibaryonInvaMassGen);
 fHistList->Add(fHistHDibaryonInvaMassGenRes);
 fHistList->Add(fHistAntiHDibaryonInvaMassGen);
 fHistList->Add(fHistHDibaryonInvaMassAso);
 fHistList->Add(fHistHDibaryonInvaMassAsoReso);
 fHistList->Add(fHistAntiHDibaryonInvaMassAso);
 fHistList->Add(fHistCheck);
 fHistList->Add(fHistHPointingAngle);
 fHistList->Add(fHistMassH);
 fHistList->Add(fHistMassPpi);
 fHistList->Add(fHistMassPpiReso);
 fHistList->Add(fHistMassLpi);
 fHistList->Add(fHistMassLP);
 fHistList->Add(fHistMassLambdaFromH);
 fHistList->Add(fHistMassLambdaFromHtLorentz);
 fHistList->Add(fHistProtonPIDBb);
 fHistList->Add(fHistPionPIDBb);
 fHistList->Add(fHistProtonPIDLambda);
 fHistList->Add(fHistPionPIDLambda);
 fHistList->Add(fHistMCdcaPvtxDvtx);
 fHistList->Add(fHistMCdcaPvtxLvtx);
 fHistList->Add(fHistMCdcaDvtxLvtx);
 fHistList->Add(fHistMCangleLH);
 fHistList->Add(fHistMCdecayAngle);
 fHistList->Add(fHistMCpointingAngle);
 fHistList->Add(fHistMCap);
 fHistList->Add(fHistMCdcaPvtxDvtxReso);
 fHistList->Add(fHistMCdcaPvtxLvtxReso);
 fHistList->Add(fHistMCdcaDvtxLvtxReso);
 fHistList->Add(fHistMCangleLHReso);
 fHistList->Add(fHistMCdecayAngleReso);
 fHistList->Add(fHistMCpointingAngleReso);
 fHistList->Add(fHistMCapReso);
 fHistList->Add(fHistCentrality);
 fHistList->Add(fHistCentralityAC);
 fHistList->Add(fHistMultiplicity);
 
 fHistHilf1= new TH1F("fHistHilf1", "HD", 300, 0., 10);
 fHistHilf1->GetXaxis()->SetTitle("");
 fHistHilf1->GetYaxis()->SetTitle("Entries");

 fHistHilf2= new TH1F("fHistHilf2", "HD", 300, 0., 10);
 fHistHilf2->GetXaxis()->SetTitle("");
 fHistHilf2->GetYaxis()->SetTitle("Entries");

 fHistHilf3= new TH1F("fHistHilf3", "HD", 300, 0., 10);
 fHistHilf3->GetXaxis()->SetTitle("");
 fHistHilf3->GetYaxis()->SetTitle("Entries");

 fHistHilf4= new TH1F("fHistHilf4", "HD", 300, 0., 10);
 fHistHilf4->GetXaxis()->SetTitle("");
 fHistHilf4->GetYaxis()->SetTitle("Entries");

 fHistHilf5= new TH1F("fHistHilf5", "HD", 300, 0., 10);
 fHistHilf5->GetXaxis()->SetTitle("");
 fHistHilf5->GetYaxis()->SetTitle("Entries");

 fHistHilf6= new TH1F("fHistHilf6", "HD", 300, 0., 10);
 fHistHilf6->GetXaxis()->SetTitle("");
 fHistHilf6->GetYaxis()->SetTitle("Entries");
 
 fHistPtvsEtaGen = new TH2F("fHistPtvsEtaGen", "p_{t} vs #eta from generated H", 200,0.0,10.0, 200,-1,1);
 fHistPtvsEtaGen->GetXaxis()->SetTitle("p_{t} (GeV/c)");
 fHistPtvsEtaGen->GetYaxis()->SetTitle("#eta");
 fHistPtvsEtaGen->SetOption("scat");
 fHistPtvsEtaGen->SetMarkerStyle(kFullCircle);

 fHistPtvsEtaAso = new TH2F("fHistPtvsEtaAso", "p_{t} vs #eta from associated H", 200,0.0,10.0, 200,-1,1);
 fHistPtvsEtaAso->GetYaxis()->SetTitle("p_{t} (GeV/c)");
 fHistPtvsEtaAso->GetXaxis()->SetTitle("#eta");
 fHistPtvsEtaAso->SetOption("scat");
 fHistPtvsEtaAso->SetMarkerStyle(kFullCircle);

 fHistPtvsYGen = new TH2F("fHistPtvsYGen", "p_{t} vs rapidity from generated H", 200,0.0,10.0, 200,-1,1);
 fHistPtvsYGen->GetXaxis()->SetTitle("p_{t} (GeV/c)");
 fHistPtvsYGen->GetYaxis()->SetTitle("y");
 fHistPtvsYGen->SetOption("scat");
 fHistPtvsYGen->SetMarkerStyle(kFullCircle);

 fHistPtvsYAso = new TH2F("fHistPtvsYAso", "p_{t} vs rapidity from associated H", 200,0.0,10.0, 200,-1,1);
 fHistPtvsYAso->GetXaxis()->SetTitle("p_{t} (GeV/c)");
 fHistPtvsYAso->GetYaxis()->SetTitle("y");
 fHistPtvsYAso->SetOption("scat");
 fHistPtvsYAso->SetMarkerStyle(kFullCircle);

 fHistRap= new TH1F("fHistRap", "Rapidity", 400, -2., 2);
 fHistRap->GetXaxis()->SetTitle("Y");
 fHistRap->GetYaxis()->SetTitle("Entries");
 fHistPtvsEtaAso->SetMarkerStyle(kFullCircle);

 fHistList->Add(fHistHilf1);
 fHistList->Add(fHistHilf2); 
 fHistList->Add(fHistHilf3);
 fHistList->Add(fHistHilf4);
 fHistList->Add(fHistHilf5); 
 fHistList->Add(fHistHilf6);
 fHistList->Add(fHistPtvsEtaGen);
 fHistList->Add(fHistPtvsEtaAso);
 fHistList->Add(fHistPtvsYGen);
 fHistList->Add(fHistPtvsYAso);
 fHistList->Add(fHistRap);

 fHistCount = new TH1F("fHistCount","test",17,0,17);
 fHistCount->GetXaxis()->SetBinLabel(1,"Events");
 fHistCount->GetXaxis()->SetBinLabel(2,"MC All");
 fHistCount->GetXaxis()->SetBinLabel(3,"MC from Primary Vtx");
 fHistCount->GetXaxis()->SetBinLabel(4,"Horst");
 fHistCount->GetXaxis()->SetBinLabel(5,"Lambda Candidates");
 fHistCount->GetXaxis()->SetBinLabel(6,"Sigma Candidates");
 fHistCount->GetXaxis()->SetBinLabel(7,"Horst");
 fHistCount->GetXaxis()->SetBinLabel(8,"Horst");
 fHistCount->GetXaxis()->SetBinLabel(9,"Horst");
 fHistCount->GetXaxis()->SetBinLabel(10,"MC All #bar{Lambda}(1520)s");
 fHistCount->GetXaxis()->SetBinLabel(11,"");
 fHistCount->GetXaxis()->SetBinLabel(12,"H-Dibaryon");
 fHistCount->GetXaxis()->SetBinLabel(13,"Hypertriton 2-Body");
 fHistCount->GetXaxis()->SetBinLabel(14,"Hypertriton 3-Body");
 fHistCount->GetXaxis()->SetBinLabel(15,"");
 fHistCount->GetXaxis()->SetBinLabel(16,"");
 fHistCount->GetXaxis()->SetBinLabel(17,"Lambdas!!!");
 fHistCount->SetStats(0);
 fHistCount->SetFillColor(38);
 fHistList->Add(fHistCount);

 //trigger statistics histogram
  fHistTriggerStat = new TH1F("fHistTriggerStat","Trigger statistics", 4,-0.5, 3.5);
  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral" };
  for ( Int_t ii=0; ii < 3; ii++ )
    fHistTriggerStat->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);

  fHistTriggerStatAfterEventSelection = new TH1F("fHistTriggerStatAfterEventSelection","Trigger statistics after event selection", 4,-0.5, 3.5);
  for ( Int_t ii=0; ii < 3; ii++ )
    fHistTriggerStatAfterEventSelection->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);
  fHistList->Add(fHistTriggerStat);
  fHistList->Add(fHistTriggerStatAfterEventSelection);

  fHistMassHcentMult =  new TH3F("fHistMassHcentMult", "Inv. Mass vs. centrality vs. multiplicity", 100, 2.2, 2.3, 5, 0, 4, 300, 0, 6000);
  fHistMassHcentMult->GetXaxis()->SetTitle("Invariant mass #Lambdap#pi^{-} (GeV/c^{2})"); // inv. mass
  fHistMassHcentMult->GetYaxis()->SetTitle("Centrality"); // triggertype
  fHistMassHcentMult->GetZaxis()->SetTitle("Multiplicity"); // refTPC
  fHistList->Add(fHistMassHcentMult);

  
  const Double_t kz = 2*TMath::Pi();
  Int_t binsD01[16]={  300, 200, 100, 100, 100, 100, 100, 100, 200, 200, 200, 200, 400, 200, 200, 3};
  Double_t xminD01[16]={2.0, 1.0, 0., -1, 0., 0., 0., 0., 0., 0., 0., 0., -1, 0.,  0., 0};
  Double_t xmaxD01[16]={2.3, 1.2, kz, 1, 1, 10, 10, 5, 5, 5, 5, 100, 1, 100,  4000, 1};
  

  fHistNdim = new THnSparseF("fHistNdim","THnS;InvMass, InvMassLambda, pointingAngle, armPoAlpha, armPoQt, pTL, pTH, d0p, d0n, dcaHd, dca, decayL, cosPA, centr, multi, mcinf;InvMassH", 16,binsD01,xminD01,xmaxD01);
  fHistList->Add(fHistNdim);

 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDtpcESD = inputHandler->GetPIDResponse();

// Post output data (if histograms are not used later, PostData is at least called here)
  PostData(1, fHistList);

}

 //________________________________________________________________________
void AliAnalysisTaskHdibaryonLPpi::UserExec(Option_t *) 
{ 
  // Main loop
  // Called for each event

  //define improtant masses
  Double_t pionMass     = 0.13957;
  Double_t protonMass   = 0.93827;

  //define PDGCodes
  Long_t pdgPionPlus          =         211;
  Long_t pdgPionMinus         =        -211;
  Long_t pdgProton            =        2212;
  Long_t pdgAntiProton        =       -2212;
  Long_t pdgLambda            =        3122;
  Long_t pdgAntiLambda        =       -3122;
  Long_t pdgHDibaryon         =  1020000020;
  Long_t pdgAntiHDibaryon     = -1020000020;

  AliStack* stack(NULL);
  if(HasMC()){

    if(!fMCEvent)return;

    AliHeader *head = fMCEvent->Header();
    if(!head)return;
    AliGenPythiaEventHeader *pyheader = (AliGenPythiaEventHeader*)head->GenEventHeader();
    if(!pyheader)return;
    
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
      if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()){
	if(!(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->TreeK()))return;
	if(!(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->TreeTR()))return;
	if(!(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->Init("local")))return;
	stack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
	if(!(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack()->TreeK()))return;
      }
    }
	
    if(!stack)return;
  }

  // -------------------------------------------------------
  // Loop for Inv. Mass via ESD tracks
  // -------------------------------------------------------

  fHistCount->Fill(0);
  
  fESD = dynamic_cast<AliESDEvent *>(fInputEvent);

    if (!fESD) {
      //Printf("ERROR: fESD not available");
      return;
    }


  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if (vertex->GetNContributors()<1) 
    {
      // SPD vertex
      vertex = fESD->GetPrimaryVertexSPD();
       if(vertex->GetNContributors()<1) {
	 return;
       }  
    }
  if (TMath::Abs(vertex->GetZ()) > 10) return;
  
  Int_t centrality = -5;
  Double_t centrPerc = -5;

  if (fESD->GetEventSpecie() == 4) 
    { // PbPb
      AliCentrality *esdCentrality = fESD->GetCentrality();
      centrality = esdCentrality->GetCentralityClass10("V0M"); // centrality percentile determined with V0
      centrPerc = esdCentrality->GetCentralityPercentile("V0M");
//      if (centrality < 0. || centrality > 8.) return; //0 bis 80 %
      if (centrality > 8) return; //0 bis 80 %
      //  cout<<"Centrality: "<< centrality << endl;
    }


  fHistCentrality->Fill(centrality);   

  //*****************//  
  //*   Centrality  *//
  //*****************//
 
  //  Float_t percentile=centrality->GetCentralityPercentile("V0M");

  Bool_t isSelectedCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
 
  Int_t triggertype = 17;

  if(isSelectedCentral){
    fHistTriggerStat->Fill(1);
    triggertype=1;
  }

  if(isSelectedSemiCentral){
    fHistTriggerStat->Fill(2);
    triggertype=2;
  }

  if(isSelectedMB){
    fHistTriggerStat->Fill(0);
    triggertype=3;
  }

  //  if(isSelectedCentral || isSelectedSemiCentral || isSelectedMB){

  //*******************************

  Int_t runNumber = 0;
  //  itrk = 0;
  runNumber = fESD->GetRunNumber();
/*  
  if (!fPIDtpcESD) fPIDtpcESD = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (!fPIDtpcESD) {
    fPIDtpcESD = new AliESDpid(); // HACK FOR MC PBPB --> PLEASE REMOVE AS SOON AS POSSIBLE
    fPIDtpcESD->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
*/
  Double_t pionK=1;
  Double_t pK=1;
 
  TObjArray* listCrossV0   = fESDCutsV0->GetAcceptedV0s(fESD);  
  Int_t nGoodV0s      = listCrossV0->GetEntries();

  const AliESDVertex *esdVer = fESD->GetPrimaryVertex();
  AliESDVertex *esdVer1 = new AliESDVertex(*esdVer);
   
  AliVertexerTracks *vertexer = new AliVertexerTracks(fESD->GetMagneticField());
  TObjArray *trkArray = new TObjArray(2);
  AliVertexerTracks *vertexer1 = new AliVertexerTracks(fESD->GetMagneticField());
  TObjArray *trkArray1 = new TObjArray(2);

  AliKFParticle::SetField(fESD->GetMagneticField());
 	
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));

  Int_t refMultTpc = AliESDtrackCuts::GetReferenceMultiplicity(fESD, kTRUE);
  //cout<<"Multiplicity: "<< refMultTpc << endl;
  fHistMultiplicity->Fill(refMultTpc);
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0}; 
  Double_t mm[3] = {0,0,0};
  Double_t dd[3] = {0,0,0};
  Double_t dd1[3] = {0,0,0};
  const Double_t cProtonMass=TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass=TDatabasePDG::Instance()->GetParticle(11)->Mass();
  const Double_t cLambdaMass=TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t decayLength=0;
  Double_t decayLengthH=0;

  //V0 Loop for Lambda and Anti-Lambda
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s ; iV0MI++) {
    AliESDv0 * fV0MIs = fESD->GetV0(iV0MI);
    Int_t    lOnFlyStatus = 0;

    lOnFlyStatus = fV0MIs->GetOnFlyStatus();
    Double_t lInvMassLambda=0;
    Double_t lInvMassLambdaPi=0;
    Double_t lPtLambda=0;
    Double_t lPzLambda=0;
    Double_t lPLambda=0;
    Int_t onl=0;
    Int_t offl=0;

    TLorentzVector posE;
    TLorentzVector negE;
    TLorentzVector photon;

    if (lOnFlyStatus){
      onl=1; 
      //      return;
    }
    if (!lOnFlyStatus){
      offl=1;
      //return;
    }

    //    fHistMultiplicity->Fill(refMultTpc); 
    fHistCentralityAC->Fill(centrality);
    fHistCheck->Fill(offl,onl);

    AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

    //    if (!
    if (!fEsdTrackCuts->AcceptTrack(trackPosTest)) continue;
    if (!fESDtrackCutsV0->AcceptTrack(trackNegTest)) continue;


      //PID via specific energy loss in the TPC
      //define the arrays for the Bethe-Bloch-Parameters
      Double_t parProton[5]   = {0,0,0,0,0};

      if(runNumber < 166500) //LHC10h
	{
	  parProton[0]   = 1.45802; // ALEPH parameters for protons (pass2)
	  parProton[1]   = 27.4992;
	  parProton[2]   = 4.00313e-15;
	  parProton[3]   = 2.48485;
	  parProton[4]   = 8.31768; 
	}
      
      if(runNumber > 166500) //LHC11h
	{ 
	  parProton[0]   = 1.11243; // ALEPH parameters for protons (pass2)
	  parProton[1]   = 26.1144;
	  parProton[2]   = 4.00313e-15;
	  parProton[3]   = 2.72969;
	  parProton[4]   = 9.15038; 
	}
 
      //Get the total momentum for each track at the inner readout of the TPC
      Double_t ptotN = trackNegTest->GetInnerParam()->GetP();
      Double_t ptotP = trackPosTest->GetInnerParam()->GetP();

      //define expected signals for the various species
      Double_t expSignalPionP = 0;
      Double_t expSignalPionN = 0;
      Double_t expSignalProtonN = 0;
      Double_t expSignalProtonP = 0;

      //for data
      if(!HasMC())
	{
	  expSignalProtonN = AliExternalTrackParam::BetheBlochAleph(ptotN/(protonMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]);
	  expSignalProtonP = AliExternalTrackParam::BetheBlochAleph(ptotP/(protonMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]);
	}
    
      //for MC
      if(HasMC())
	{
	  expSignalPionP = 0.7*AliExternalTrackParam::BetheBlochAleph(ptotP/(pionMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]); 
	  expSignalPionN = 0.7*AliExternalTrackParam::BetheBlochAleph(ptotN/(pionMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]); 

	  expSignalProtonN = 0.65*AliExternalTrackParam::BetheBlochAleph(ptotN/(protonMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]);
	  expSignalProtonP = 0.65*AliExternalTrackParam::BetheBlochAleph(ptotP/(protonMass),parProton[0],parProton[1],parProton[2],parProton[3],parProton[4]);	 
	}

      // PID cut on the nuclei (proton, deuteron, triton, helium3)
      Bool_t corrParticle = kFALSE;

      Bool_t posProton = kFALSE;
      //data
      if(!HasMC())
	{
	  if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackPosTest, AliPID::kProton)) < 3)
	    {
	      posProton = kTRUE;
	      corrParticle = kTRUE;
	    }
	}
      //MC
      if(HasMC())
	{
	  if(//trackPosTest->GetTPCsignal() < 1200 && 
	     TMath::Abs(trackPosTest->GetTPCsignal() - expSignalProtonP)/expSignalProtonP < 0.4)
	    {
	      posProton = kTRUE;
	      corrParticle = kTRUE;
	    }
	}

      Bool_t negProton = kFALSE;
      //data
      if(!HasMC())
	{
	  if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackNegTest, AliPID::kProton)) < 3)
	    {
	      negProton = kTRUE;
	      corrParticle = kTRUE;
	    }
	}
     //MC
      if(HasMC())
	{
	  if(//trackNegTest->GetTPCsignal() < 1200 && 
	     TMath::Abs(trackNegTest->GetTPCsignal() - expSignalProtonN)/expSignalProtonN < 0.4)
	    {
	      negProton = kTRUE;
	      corrParticle = kTRUE;
	    }
	}

     //PID cut for pions
      //data: 3sigma cut on the pions

      Bool_t negPion = kFALSE;     
      Bool_t posPion = kFALSE;

      if (!HasMC())
	{
	  if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackPosTest, AliPID::kPion)) < 4) posPion=kTRUE; //pos daughter has to be a pion
	  if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackNegTest, AliPID::kPion)) < 4) negPion=kTRUE; // negative daughter has to be a pion
	}
      
      //MC: like the nuclei via the specific energyloss in the TPC
      if (HasMC())
	{
	  if (TMath::Abs(trackPosTest->GetTPCsignal() - expSignalPionP)/expSignalPionP < 0.4) posPion=kTRUE;
	  if (TMath::Abs(trackNegTest->GetTPCsignal() - expSignalPionN)/expSignalPionN < 0.4) negPion=kTRUE;
	}

      if (!(posProton==kTRUE)) continue;
      if (!(negPion==kTRUE)) continue;

    //To avoid ghosts

    if( !(trackPosTest->GetStatus() & AliESDtrack::kTPCrefit)){
      continue;
    }

    if( !(trackNegTest->GetStatus() & AliESDtrack::kTPCrefit)){
      continue;
    }

    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0MIs->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0MIs->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter 
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0MIs->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0MIs->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    fV0MIs->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);
    
    Double_t thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    
    Double_t alfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/
      ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) ;
    Double_t qt = vecP.Mag()*sin(thetaP);

    fHistArmenterosPodolanski->Fill(alfa,qt);    //Armenteros-Podolanski calculation

    TLorentzVector k0;
    TLorentzVector k0daugh1;
    TLorentzVector k0daugh2;
    TLorentzVector proton;
    TLorentzVector pq;

    k0daugh1.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    k0daugh2.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    k0=k0daugh1+k0daugh2;

    fV0MIs->ChangeMassHypothesis(3122);
    lInvMassLambda = fV0MIs->GetEffMass();
    lPtLambda = fV0MIs->Pt();
    lPzLambda = fV0MIs->Pz();
    lPLambda = fV0MIs->P();

    trkArray->AddAt(trackPosTest,0);
    trkArray->AddAt(trackNegTest,1);

    vertexer->SetVtxStart(esdVer1);
    AliESDVertex *decayVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);

    dd[0]=fESD->GetPrimaryVertexSPD()->GetX()-decayVertex->GetX();
    dd[1]=fESD->GetPrimaryVertexSPD()->GetY()-decayVertex->GetY();
    dd[2]=fESD->GetPrimaryVertexSPD()->GetZ()-decayVertex->GetZ();
    
    decayLength=sqrt(dd[0]*dd[0]+dd[1]*dd[1]+dd[2]*dd[2]);

      if (decayVertex == NULL) cout << "Lambda decay vtx pointer NULL" << endl;
    if (decayVertex) delete decayVertex;

    TLorentzVector negPio1;
    TLorentzVector posProt1;        
    TLorentzVector posP;
    TLorentzVector posProt;
    TLorentzVector negK;
    TLorentzVector negPio;
    TLorentzVector negPi;
    TLorentzVector omega;
    TLorentzVector threeSum;
    TLorentzVector fourSum;
    TLorentzVector ppK;
    TLorentzVector posPiK;
    TLorentzVector negPiK;
    TLorentzVector kaon;
    TLorentzVector lambda;
    TLorentzVector lambdaH;
    TLorentzVector hDibaryon;
    TVector3 h;
    TVector3 h1;

    Int_t mcStatus=0;

    h.SetXYZ(-dd[0],-dd[1],-dd[2]);


    if (offl==1)fHistMassDPi->Fill(lInvMassLambda);

    if (onl==1){
      fHistMassLPi->Fill(lInvMassLambda);
      
      negE.SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
      posE.SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
      photon=posE+negE;
      
      negPiK.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
      posPiK.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
      kaon=posPiK+negPiK;
      
      negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
      posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
      lambda=negPi+posP;
      lambdaH.SetXYZM(mm[0],mm[1],mm[2],cLambdaMass);

      if (lInvMassLambda>1.1113&&lInvMassLambda<1.1202){
	
	if (!HasMC()){
	  if (qt<-2.21*alfa*alfa+2.945*alfa-0.887) continue;
	  if (qt>-2.21*alfa*alfa+2.945*alfa-0.873) continue;
	  if (photon.M()<0.005) continue;
	  if (kaon.M()>0.495 && kaon.M()<0.500 ) continue;
	}
	

	fHistMassLambda->Fill(lInvMassLambda);
	//
	Bool_t isCorrectlyAssociatedLambda = kFALSE;
	Bool_t isPartialCorrectlyAssociatedLambda = kFALSE;
	Int_t labelAssociatedH=-1; 
	Int_t labelLambda=-1;
	//
	if (HasMC()) {
	  Int_t labelPosTest = trackPosTest->GetLabel();
	  TParticle *tparticleDaughter = stack->Particle(TMath::Abs(labelPosTest));
	  Int_t labelMother = tparticleDaughter->GetFirstMother();
	  TParticle *tparticleMother = stack->Particle(TMath::Abs(labelMother));
	  
	  Int_t labelOma = tparticleMother->GetFirstMother();
	  TParticle *tparticleOma = stack->Particle(TMath::Abs(labelOma));

	  if ((tparticleOma->GetPdgCode() < 0) && TMath::Abs(tparticleMother->GetPdgCode())==pdgLambda){// check mother to be Lambda 
	    Int_t labelFirstDaughter  = tparticleMother->GetDaughter(1);// Proton
	    Int_t labelThirdDaughter  = tparticleMother->GetDaughter(0);// Pion
	    
	    TParticle *tparticleFirstDaughter  = stack->Particle(TMath::Abs(labelFirstDaughter));
	    TParticle *tparticleThirdDaughter  = stack->Particle(TMath::Abs(labelThirdDaughter));
	    
	    if((tparticleFirstDaughter->GetPdgCode() == pdgProton && tparticleThirdDaughter->GetPdgCode()== pdgPionMinus) || 
	       (tparticleFirstDaughter->GetPdgCode() == pdgPionMinus && tparticleThirdDaughter->GetPdgCode()== pdgProton)){ //daughter PDGs
	      isPartialCorrectlyAssociatedLambda = kTRUE;
	      labelLambda=labelMother;
	    }
	  }
	  
	  //H-Dibaryon
	  if(tparticleOma->GetPdgCode() == pdgHDibaryon){ //check grandmother to be H PDG
	    if (TMath::Abs(tparticleMother->GetPdgCode())==pdgLambda){// check mother to be Lambda 
	      Int_t labelFirstDaughter  = tparticleMother->GetDaughter(1);// Proton
	      Int_t labelThirdDaughter  = tparticleMother->GetDaughter(0);// Pion
			
	      TParticle *tparticleFirstDaughter  = stack->Particle(TMath::Abs(labelFirstDaughter));
	      TParticle *tparticleThirdDaughter  = stack->Particle(TMath::Abs(labelThirdDaughter));
	      
	      if((tparticleFirstDaughter->GetPdgCode() == pdgProton && tparticleThirdDaughter->GetPdgCode()== pdgPionMinus) || 
		 (tparticleFirstDaughter->GetPdgCode() == pdgPionMinus && tparticleThirdDaughter->GetPdgCode()== pdgProton)){ //daughter PDGs
		isCorrectlyAssociatedLambda = kTRUE;
		labelAssociatedH=labelOma;
		fHistMassLambdaFromH->Fill(lInvMassLambda);
		fHistMassLambdaFromHtLorentz->Fill(lambda.M());
	      }
	    }
	  }
	}

	fHistProtonPIDLambda->Fill(trackPosTest->GetInnerParam()->GetP(), trackPosTest->GetTPCsignal());
	fHistPionPIDLambda->Fill(trackNegTest->GetInnerParam()->GetP(), trackNegTest->GetTPCsignal());

	//---------------------------------------------------------
	// Proton track loop
	//---------------------------------------------------------
	fHistArmenterosPodolanskiCut->Fill(alfa,qt);

	for (Int_t iTracksP = 0; iTracksP < fESD->GetNumberOfTracks(); iTracksP++) {
	  AliESDtrack* trackP = dynamic_cast<AliESDtrack*> (fESD->GetTrack(iTracksP));
	  if (trackP->GetSign()<0) continue;
	  
	  if (iTracksP==fV0MIs->GetPindex())continue;
	  if (iTracksP==fV0MIs->GetNindex())continue;
		
	  if (!fEsdTrackCuts->AcceptTrack(trackP)) continue;
		
	  AliKFParticle protonKF( *(trackP), 2212);

	  if (!trackP->GetInnerParam()) continue;
	  
	  if (HasMC()) {
	    pK=0.65;
	  }

	  if (!HasMC()){
	    if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackP, AliPID::kProton)) > 3) continue;
	  }		
	  
	  fHistProtonPIDBb->Fill(trackP->GetInnerParam()->GetP(), trackP->GetTPCsignal());
	  
	  posProt.SetXYZM(trackP->Px(),trackP->Py(),trackP->Pz(),cProtonMass);

	  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //Pion Track loop!!!!!!!!!!!!!!!!!!!!!
	  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  for (Int_t iTracksN = iTracksP+1; iTracksN < fESD->GetNumberOfTracks(); iTracksN++) {
	    AliESDtrack* trackN = dynamic_cast<AliESDtrack*> (fESD->GetTrack(iTracksN));

	    if (iTracksN==fV0MIs->GetPindex())continue;
	    if (iTracksN==fV0MIs->GetNindex())continue;
	    if (trackN->GetSign()>0) continue;
	  
	    if (!fEsdTrackCuts->AcceptTrack(trackN)) continue;
	    if (!fESDtrackCutsV0->AcceptTrack(trackN)) continue;

	    Double_t bz = fESD->GetMagneticField();
	   
	    Double_t xthiss(0.0);
	    Double_t xpp(0.0);
	    Double_t dca = trackN->GetDCA(trackP,bz,xthiss,xpp);
	    if (dca>0.5) continue;

	    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
	    posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
	    negPio.SetXYZM(trackN->Px(),trackN->Py(),trackN->Pz(),cPionMass);
	      
	    threeSum=negPi+posP+negPio;
	    lInvMassLambdaPi=threeSum.M();	   
	    
	    fHistDC->Fill(decayLength*lInvMassLambda/lPLambda);
	    
	    AliKFParticle posPionKF( *(trackN) ,-211);
	    
	    if (!trackN->GetInnerParam()) continue;
	    if (HasMC()) {
	      pionK=0.7;
	    }
	    
	  //  if (!HasMC()){
	      if (TMath::Abs(fPIDtpcESD->NumberOfSigmasTPC(trackN, AliPID::kPion)) > 3) continue;
	  //  }
	    fHistPionPIDBb->Fill(trackN->GetInnerParam()->GetP(), trackN->GetTPCsignal());

	    trkArray1->AddAt(trackP,0);
	    trkArray1->AddAt(trackN,1);
	    
	    vertexer1->SetVtxStart(esdVer1);
	    AliESDVertex *decayVertex1 = (AliESDVertex*)vertexer1->VertexForSelectedESDTracks(trkArray1);
	    
	    dd1[0]=fESD->GetPrimaryVertexSPD()->GetX()-decayVertex1->GetX();
	    dd1[1]=fESD->GetPrimaryVertexSPD()->GetY()-decayVertex1->GetY();
	    dd1[2]=fESD->GetPrimaryVertexSPD()->GetZ()-decayVertex1->GetZ();
	    
	    decayLengthH=sqrt(dd1[0]*dd1[0]+dd1[1]*dd1[1]+dd1[2]*dd1[2]);

	      //            Double_t bz = fESD->GetMagneticField();
	    
	    trackP->PropagateToDCA(decayVertex1, bz, 10);
	    trackN->PropagateToDCA(decayVertex1, bz, 10);

	      //	    Double_t xthiss(0.0);
	      //	    Double_t xpp(0.0);
	      //	    Double_t dca = trackN->GetDCA(trackP,bz,xthiss,xpp);

	      if (decayVertex1 == NULL) cout << "Secondary decay vtx pointer NULL" << endl;
	    if (decayVertex1) delete decayVertex1;
	    h1.SetXYZ(-dd1[0],-dd1[1],-dd1[2]);

	    //	    if (dca>1) continue;
	      //	    if (dca>0.1) continue;

	    fourSum=threeSum+posProt;

	    posProt1.SetXYZM(trackP->Px(),trackP->Py(),trackP->Pz(),cProtonMass);
	    negPio1.SetXYZM(trackN->Px(),trackN->Py(),trackN->Pz(),cPionMass);
	    hDibaryon=lambdaH+posProt1+negPio1;
	    Double_t hPointingAngle = hDibaryon.Angle(h);
	    Double_t pointingAngleH = hDibaryon.Angle(h1);
	    Double_t decayAngleH = h.Angle(h1);
	    TVector3 vecDist(dd[0]-dd1[0],dd[1]-dd1[1],dd[2]-dd1[2]);
	    fHistMassLambdaPPi->Fill(hDibaryon.M());
	    fHistHPointingAngle->Fill(pointingAngleH);

            fHistMassHcentMult->Fill(hDibaryon.M(),triggertype,refMultTpc);

	    Double_t rapidity = hDibaryon.Rapidity();
	    if(rapidity > 1.0 || rapidity < -1.0) continue;

	      //Double_t vec[16]={hDibaryon.M(), lInvMassLambda, pointingAngleH, alfa, qt, lPtLambda, hDibaryon.Pt(), posPionKF.GetDistanceFromVertex(primVtx), protonKF.GetDistanceFromVertex(primVtx), dca, protonKF.GetDistanceFromVertex(posPionKF), TMath::Cos(pointingAngleH), centrPerc, refMultTpc, mcStatus};
	      //fHistNdim->Fill(vec);

	    fHistRap->Fill(rapidity);
	    //if (pointingAngleH > 0.1) continue;
	    if (pointingAngleH > 0.05) continue;

	    ///////////////////////////
	    //MC part for Associated H
	    ///////////////////////////

	    if (HasMC() && isCorrectlyAssociatedLambda) {
	      Int_t labelP = trackP->GetLabel();
	      TParticle *tparticleDaughter = stack->Particle(TMath::Abs(labelP));
	      Int_t labelMother = tparticleDaughter->GetFirstMother();
	      TParticle *tparticleMother = stack->Particle(TMath::Abs(labelMother));
	      
	      Int_t labelProton = trackP->GetLabel();
	      Int_t labelPion = trackN->GetLabel();
	      
	      //H-Dibaryon
	      if(tparticleMother->GetPdgCode() == pdgHDibaryon && labelAssociatedH==labelMother){ //check mother PDG
		Int_t labelFirstDaughter  = tparticleMother->GetDaughter(0);
		Int_t labelThirdDaughter  = tparticleMother->GetDaughter(1);
		Int_t labelSecondDaughter = labelFirstDaughter +1;

		TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));
		TParticle *tparticleThirdDaughter  = stack->Particle(TMath::Abs(labelThirdDaughter));
		
		TLorentzVector ppi;
		TLorentzVector lpi;
		TLorentzVector lP;
		
		ppi=posProt+negPio;
		lpi=lambdaH+negPio;
		lP=lambdaH+posProt;

		if((tparticleThirdDaughter->GetPdgCode() == pdgPionMinus && labelPion==labelThirdDaughter)&&(tparticleSecondDaughter->GetPdgCode() == pdgProton||tparticleSecondDaughter->GetPdgCode() == pdgPionMinus) && labelProton==labelSecondDaughter) fHistMassPpi->Fill(ppi.M());
		
		if(tparticleThirdDaughter->GetPdgCode() == pdgPionMinus && labelPion==labelThirdDaughter) fHistMassLpi->Fill(lpi.M());
		
		if(tparticleSecondDaughter->GetPdgCode() == pdgProton && labelProton==labelSecondDaughter) fHistMassLP->Fill(lP.M());
		
		if(tparticleSecondDaughter->GetPdgCode() == pdgProton && labelProton==labelSecondDaughter){//check second daughter PDG
		  if(tparticleThirdDaughter->GetPdgCode() == pdgPionMinus && labelPion==labelThirdDaughter){//check second daughter PDG

		    fHistHDibaryonInvaMassAso->Fill(hDibaryon.M()); 
		    
		    Double_t distance01=vecDist.Mag();
		    fHistMCdcaPvtxDvtx->Fill(decayLengthH);
		    fHistMCdcaPvtxLvtx->Fill(decayLength);
		    fHistMCdcaDvtxLvtx->Fill(distance01);
		    fHistMCangleLH->Fill(hPointingAngle);
		    fHistMCdecayAngle->Fill(decayAngleH);
		    fHistMCpointingAngle->Fill(pointingAngleH);
		    fHistMCap->Fill(alfa,qt);

		    fHistHilf1->Fill(posPionKF.GetDistanceFromVertex(primVtx));
		    fHistHilf2->Fill(protonKF.GetDistanceFromVertex(primVtx));
		    fHistHilf3->Fill(protonKF.GetDistanceFromVertex(posPionKF));
		    fHistHilf6->Fill(dca);
		    fHistPtvsYAso->Fill(hDibaryon.Pt(),hDibaryon.Rapidity());
		    fHistPtvsEtaAso->Fill(hDibaryon.Pt(),hDibaryon.Eta());
		    mcStatus=1;
		  }//end check for third daughter PDG
		}//end check second daughter PDG
	      }//end H-Dibaryon
	    }//end MC
	    
	    //	    cout<<"Trigger: "<<triggertype<<endl;
	    fHistMassH->Fill(hDibaryon.M());
	    fHistMassHcentMult->Fill(hDibaryon.M(),triggertype,refMultTpc);
	    ppK=lambdaH+posProt;
	    fHistMassLambdaP->Fill(ppK.M());

	      //fHistNdim = new THnSparseF("fHistNdim","THnS;InvMass, InvMassLambda, pointingAngle, armPoAlpha, armPoQt, pTL, pTH, d0p, d0n, dcaHd, dca, decayL, cosPA, centr, multi, mcinf;InvMassH", 16,binsD01,xminD01,xmaxD01);

	    Double_t vec[16]={hDibaryon.M(), lInvMassLambda, pointingAngleH, alfa, qt, lPtLambda, hDibaryon.Pt(), posPionKF.GetDistanceFromVertex(primVtx), protonKF.GetDistanceFromVertex(primVtx), dca, protonKF.GetDistanceFromVertex(posPionKF), TMath::Cos(pointingAngleH), centrPerc, static_cast<Double_t>(refMultTpc), static_cast<Double_t>(mcStatus)};
		      fHistNdim->Fill(vec);

	  }
	}
      }
    }    
  }
  
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Pure MC Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Monte Carlo for genenerated particles
  if (HasMC()) //MC loop  
    {

      Int_t stackN = 0;

      Double_t momentumPionGen[3]={0,0,0};
      Double_t momentumNucleonGen[3]={0,0,0};
      Double_t momentumLambdaGen[3]={0,0,0};

      Double_t energyPionGen = 0;
      Double_t energyNucleonGen = 0;
      Double_t energyLambdaGen = 0;

      Double_t transversMomentumMotherGen = 0;
      Double_t longitudinalMomentumMotherGen = 0;
      Double_t totalEnergyMotherGen = 0;
      
      Double_t rapidityGen = 2;

      for(stackN = 0; stackN < stack->GetNtrack(); stackN++) //loop over stack
	{

	  TParticle *tparticleMother = stack->Particle(stackN);

	  if(tparticleMother->GetPdgCode() == pdgLambda) fHistCount->Fill(16); 

	  //H-Dibaryon
	  if(tparticleMother->GetPdgCode() == pdgHDibaryon) //check mother PDG
	    {
	      Int_t labelFirstDaughter  = tparticleMother->GetDaughter(0);
	      Int_t labelThirdDaughter  = tparticleMother->GetDaughter(1);
	      Int_t labelSecondDaughter = labelFirstDaughter +1;

	      TParticle *tparticleFirstDaughter  = stack->Particle(TMath::Abs(labelFirstDaughter));
	      TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));
	      TParticle *tparticleThirdDaughter  = stack->Particle(TMath::Abs(labelThirdDaughter));
	     
	      if(tparticleFirstDaughter->GetPdgCode() == pdgLambda) //check first daughter PDG
		{
		  if(tparticleSecondDaughter->GetPdgCode() == pdgProton)//check second daughter PDG
		    {
		      if(tparticleThirdDaughter->GetPdgCode() == pdgPionMinus)//check second daughter PDG
			{		 		      
			  momentumLambdaGen[0] = tparticleFirstDaughter->Px();
			  momentumLambdaGen[1] = tparticleFirstDaughter->Py();
			  momentumLambdaGen[2] = tparticleFirstDaughter->Pz();
		
			  momentumNucleonGen[0] = tparticleSecondDaughter->Px();
			  momentumNucleonGen[1] = tparticleSecondDaughter->Py();
			  momentumNucleonGen[2] = tparticleSecondDaughter->Pz();

			  momentumPionGen[0] = tparticleThirdDaughter->Px();
			  momentumPionGen[1] = tparticleThirdDaughter->Py();
			  momentumPionGen[2] = tparticleThirdDaughter->Pz();

			  TLorentzVector lorentzVectorLambda;
			  TLorentzVector lorentzVectorProton;
			  TLorentzVector lorentzVectorPion;
			  TLorentzVector lorentzVectorHDibaryon;
			  
			  lorentzVectorLambda.SetXYZM(momentumLambdaGen[0],momentumLambdaGen[1],momentumLambdaGen[2],1.115);
			  lorentzVectorProton.SetXYZM(momentumNucleonGen[0],momentumNucleonGen[1],momentumNucleonGen[2],protonMass);
			  lorentzVectorPion.SetXYZM(momentumPionGen[0],momentumPionGen[1],momentumPionGen[2],pionMass);
			  
			  lorentzVectorHDibaryon = lorentzVectorLambda + lorentzVectorProton + lorentzVectorPion;
			  rapidityGen=lorentzVectorHDibaryon.Rapidity();
			  transversMomentumMotherGen = lorentzVectorHDibaryon.Pt();
			  longitudinalMomentumMotherGen = lorentzVectorHDibaryon.Pz();
			  totalEnergyMotherGen = lorentzVectorHDibaryon.Energy();

			  if(rapidityGen > 1.0 || rapidityGen < -1 ) continue;
			  //lorentzVectorLambda
			  fHistHDibaryonInvaMassGen->Fill(lorentzVectorHDibaryon.M()); 
			  if (lorentzVectorLambda.Rapidity()  > 1.0 || lorentzVectorLambda.Rapidity() < -1) continue;
			  if (lorentzVectorProton.Rapidity()  > 1.0 ||      lorentzVectorProton.Rapidity() < -1) continue;

			  if (lorentzVectorPion.Rapidity()  > 1.0 ||  lorentzVectorPion.Rapidity() < -1) continue;
			  fHistHDibaryonInvaMassGenRes->Fill(lorentzVectorHDibaryon.M());
			  fHistPtvsEtaGen->Fill(lorentzVectorHDibaryon.Pt(),lorentzVectorHDibaryon.Eta());
			  fHistPtvsYGen->Fill(lorentzVectorHDibaryon.Pt(),lorentzVectorHDibaryon.Rapidity());
			  fHistPtvsEtaGen->Fill(lorentzVectorHDibaryon.Pt(),lorentzVectorHDibaryon.Eta());
			  fHistCount->Fill(11);
			}//end of check third daughter PDG
		    }//end of check second daughter PDG
		}//end of check first daughter PDG
	    }//end of H-Dibaryon

	  //Anti-H-Dibaryon
	  if(tparticleMother->GetPdgCode() == pdgAntiHDibaryon) //check mother PDG
	    {
	      Int_t labelFirstDaughter  = tparticleMother->GetDaughter(0);
	      Int_t labelThirdDaughter  = tparticleMother->GetDaughter(1);
	      Int_t labelSecondDaughter = labelFirstDaughter +1;

	      TParticle *tparticleFirstDaughter  = stack->Particle(TMath::Abs(labelFirstDaughter));
	      TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));
	      TParticle *tparticleThirdDaughter  = stack->Particle(TMath::Abs(labelThirdDaughter));

	      if(tparticleFirstDaughter->GetPdgCode() == pdgAntiLambda) //check first daughter PDG
		{
		  if(tparticleSecondDaughter->GetPdgCode() == pdgAntiProton)//check second daughter PDG
		    {
		      if(tparticleThirdDaughter->GetPdgCode() == pdgPionPlus)//check second daughter PDG
			{		 
			  momentumLambdaGen[0] = tparticleFirstDaughter->Px();
			  momentumLambdaGen[1] = tparticleFirstDaughter->Py();
			  momentumLambdaGen[2] = tparticleFirstDaughter->Pz();
		
			  momentumNucleonGen[0] = tparticleSecondDaughter->Px();
			  momentumNucleonGen[1] = tparticleSecondDaughter->Py();
			  momentumNucleonGen[2] = tparticleSecondDaughter->Pz();

			  momentumPionGen[0] = tparticleThirdDaughter->Px();
			  momentumPionGen[1] = tparticleThirdDaughter->Py();
			  momentumPionGen[2] = tparticleThirdDaughter->Pz();
			  
			  energyLambdaGen  = tparticleFirstDaughter->Energy();
			  energyNucleonGen = tparticleSecondDaughter->Energy();
			  energyPionGen    = tparticleThirdDaughter->Energy();
			   
			  TLorentzVector lorentzVectorLambda;
			  TLorentzVector lorentzVectorProton;
			  TLorentzVector lorentzVectorPion;
			  TLorentzVector lorentzVectorHDibaryon;
			  
			  lorentzVectorLambda.SetXYZM(momentumLambdaGen[0],momentumLambdaGen[1],momentumLambdaGen[2],1.115);
			  lorentzVectorProton.SetXYZM(momentumNucleonGen[0],momentumNucleonGen[1],momentumNucleonGen[2],protonMass);
			  lorentzVectorPion.SetXYZM(momentumPionGen[0],momentumPionGen[1],momentumPionGen[2],pionMass);
			  
			  lorentzVectorHDibaryon = lorentzVectorLambda + lorentzVectorProton + lorentzVectorPion;

			  rapidityGen=lorentzVectorHDibaryon.Rapidity();
			  if(rapidityGen > 1.0 || rapidityGen < -1 ) continue;
			  fHistAntiHDibaryonInvaMassGen->Fill(lorentzVectorHDibaryon.M()); 
			}//end of check third daughter PDG
		    }//end of check second daughter PDG
		}//end of check first daughter PDG
	    }//end of Anti-H-Dibaryon
	}      
    }//end MC

  // Post output data.
  PostData(1,fHistList);
  //PostData(0,fHistList);

    if (listCrossV0 == NULL) return;

  if (listCrossV0) delete listCrossV0;
  if (esdVer1) delete esdVer1;
  if (vertexer) delete vertexer;
  if (vertexer1) delete vertexer1;
  if (trkArray) delete trkArray;
  if (trkArray1) delete trkArray1;
}

//________________________________________________________________________
void AliAnalysisTaskHdibaryonLPpi::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

}

