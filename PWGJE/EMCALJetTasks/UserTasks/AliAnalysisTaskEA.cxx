#ifndef ALIANALYSISTASKSE_H

#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TArrayF.h>
#include <TArrayD.h>
#include <TVector2.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliParticleContainer.h"
#include "AliInputEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#endif

#include <string>
#include <time.h>
#include <TRandom3.h>
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliRhoParameter.h"
#include "TVector3.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "AliGenDPMjetEventHeader.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskEA.h"
#include "AliHeader.h" 
#include "AliRunLoader.h"  
#include "AliVVZERO.h"
#include "AliAODZDC.h" 
#include "AliVZDC.h"
#include "AliMultSelection.h"
//#include "AliEmcalDownscaleFactorsOCDB.h"
//#include "AliEmcalAnalysisFactory.h"

using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (7.Oct. 2015)

ClassImp(AliAnalysisTaskEA)
//________________________________________________________________________________________

AliAnalysisTaskEA::AliAnalysisTaskEA(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskEA", kTRUE),  
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyClusterContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fClusterContainerDetLevel(0x0),
fRhoTaskName(""),
fRhoTaskNameMC(""),
fCentralityTree(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsEmcalTrig(0),
fIsHighMultTrig(0),
fCentralityV0A(-1),
fCentralityV0C(-1),
fCentralityCL1(-1),
fCentralityZNA(-1),
fCentralityZNC(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
fVertexer3d(1),
fNTracklets(-1),
fIsV0ATriggered(0),
fIsV0CTriggered(0),
fMultV0A(0.),
fMultV0C(0.),
fMultV0Anorm(0.),
fMultV0Cnorm(0.),
fMultV0AV0Cnorm(0.),
fZEM1Energy(0),
fZEM2Energy(0),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fMC(0),
fHelperClass(0), fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhTrackPhiIncl(0x0), fhTrackEtaIncl(0x0), 
fhJetPhiIncl(0x0), fhJetEtaIncl(0x0),
fhClusterPhiInclMB(0x0), fhClusterEtaInclMB(0x0),
fhClusterPhiInclGA(0x0), fhClusterEtaInclGA(0x0),
fhRhoMB(0x0),
fhRhoHM(0x0),
fhV0AvsV0C(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhFractionOfSecInJet(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),
fhJetPtPartLevelVsJetPtDetLevelZero(0x0),
fhJetPtResolutionVsPtPartLevel(0x0),
fhOneOverPtVsPhiNeg(0x0),
fhOneOverPtVsPhiPos(0x0),
fhSigmaPtOverPtVsPt(0x0),
fhDCAinXVsPt(0x0),
fhDCAinYVsPt(0x0),
fhDCAinXVsPtPhysPrimary(0x0),
fhDCAinYVsPtPhysPrimary(0x0),
fhDCAinXVsPtSecondary(0x0),
fhDCAinYVsPtSecondary(0x0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fnJetChTTBins(0),
fnClusterTTBins(0),
fFillTTree(0),
fSystem(AliAnalysisTaskEA::kpPb),
fFiducialCellCut(0x0),
fMeanV0A(1.),
fMeanV0C(1.),
fPhiCut(TMath::Pi()-0.6),
fRandom(0)
{
   //default constructor

   for(Int_t i=0; i<2; i++) fNClusters[i] = 0;
   for(Int_t i=0; i<8; i++) fRingMultV0[i] = 0;

   for(Int_t i=0; i<5; i++){
      fZNCtower[i] = 0;
      fZPCtower[i] = 0;
      fZNAtower[i] = 0;
      fZPAtower[i] = 0;
      fZNCtowerLG[i] = 0;
      fZPCtowerLG[i] = 0;
      fZNAtowerLG[i] = 0;
      fZPAtowerLG[i] = 0;
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i]   = 0;
      fJetChTT[i]    = 0;
      fClusterTT[i]  = 0;

      fhMultTTHinMB[i] = 0x0;   
      fhMultTTHinHM[i] = 0x0;   
      fhMultTTJinMB[i] = 0x0;  
      fhMultTTJinHM[i] = 0x0;  
      fhMultTTCinMB[i] = 0x0;  
      fhMultTTCinHM[i] = 0x0;  
      fhMultTTCinGA[i] = 0x0; 

      fhTTHinMB[i] = 0x0;
      fhTTHinHM[i] = 0x0;
      fhTTCinMB[i] = 0x0;
      fhTTCinHM[i] = 0x0;
      fhTTCinGA[i] = 0x0;

      fhDrecoilTTHinMB[i] = 0x0;
      fhDrecoilTTHinHM[i] = 0x0;
      fhDrecoilTTCinMB[i] = 0x0;
      fhDrecoilTTCinHM[i] = 0x0;
      fhDrecoilTTCinGA[i] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
      fJetChTTLowPt[i]=-1;
      fJetChTTHighPt[i]=-1;
      fClusterTTLowPt[i]=-1;
      fClusterTTHighPt[i]=-1;
 
      fhV0AvsV0CTTH[i] = 0x0; 
      fhV0AvsV0CTTJ[i] = 0x0;
      fhV0AvsV0CTTCinMB[i] = 0x0;  
      fhV0AvsV0CTTCinGA[i] = 0x0; 
   }
 
   for(Int_t iv=0; iv<fkVtx;iv++){
      fhVertex[iv]=0x0;
      for(Int_t i=0; i<fkTTbins;i++){
         fhVertexTTH[iv][i]=0x0;
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){
      fhCentralityMB[ic] = 0x0;
      fhSignalMB[ic] = 0x0; 
      fhSignalHM[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhCentralityTTH[ic][i] = 0x0;
         fhCentralityTTJ[ic][i] = 0x0;
         fhCentralityTTCinMB[ic][i] = 0x0;
         fhCentralityTTCinGA[ic][i] = 0x0;

         fhSignalTTHinMB[ic][i] = 0x0;
         fhSignalTTHinHM[ic][i] = 0x0;
         fhSignalTTJinMB[ic][i] = 0x0;
         fhSignalTTJinHM[ic][i] = 0x0;
         fhSignalTTCinMB[ic][i] = 0x0;
         fhSignalTTCinHM[ic][i] = 0x0;
         fhSignalTTCinGA[ic][i] = 0x0;
      }
   }

   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMB[i]=0x0;
      fhRhoTTHinHM[i]=0x0;   
      fhRhoTTJinMB[i]=0x0;  
      fhRhoTTJinHM[i]=0x0; 
      fhRhoTTCinMB[i]=0x0; 
      fhRhoTTCinHM[i]=0x0;  
      fhRhoTTCinGA[i]=0x0; 
   }

   fFiducialCellCut = new AliEMCALRecoUtils();
 
   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1; 

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);
   }


   sprintf(fTrigClass,"%s","");
}

//________________________________________________________________________
AliAnalysisTaskEA::AliAnalysisTaskEA(const char *name): 
AliAnalysisTaskEmcalJet(name,kTRUE),  
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyClusterContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fClusterContainerDetLevel(0x0),
fRhoTaskName(""),
fRhoTaskNameMC(""),
fCentralityTree(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsEmcalTrig(0),
fIsHighMultTrig(0),
fCentralityV0A(-1),
fCentralityV0C(-1),
fCentralityCL1(-1),
fCentralityZNA(-1),
fCentralityZNC(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
fVertexer3d(1),
fNTracklets(-1),
fIsV0ATriggered(0),
fIsV0CTriggered(0),
fMultV0A(0.),
fMultV0C(0.),
fMultV0Anorm(0.),
fMultV0Cnorm(0.),
fMultV0AV0Cnorm(0.),
fZEM1Energy(0),
fZEM2Energy(0),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fMC(0),
fHelperClass(0), fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhTrackPhiIncl(0x0), fhTrackEtaIncl(0x0), 
fhJetPhiIncl(0x0), fhJetEtaIncl(0x0), 
fhClusterPhiInclMB(0x0), fhClusterEtaInclMB(0x0),
fhClusterPhiInclGA(0x0), fhClusterEtaInclGA(0x0),
fhRhoMB(0x0),
fhRhoHM(0x0),
fhV0AvsV0C(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhFractionOfSecInJet(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),
fhJetPtPartLevelVsJetPtDetLevelZero(0x0),
fhJetPtResolutionVsPtPartLevel(0x0),
fhOneOverPtVsPhiNeg(0x0),
fhOneOverPtVsPhiPos(0x0),
fhSigmaPtOverPtVsPt(0x0),
fhDCAinXVsPt(0x0),
fhDCAinYVsPt(0x0),
fhDCAinXVsPtPhysPrimary(0x0),
fhDCAinYVsPtPhysPrimary(0x0),
fhDCAinXVsPtSecondary(0x0),
fhDCAinYVsPtSecondary(0x0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fnJetChTTBins(0),
fnClusterTTBins(0),
fFillTTree(0),
fSystem(AliAnalysisTaskEA::kpPb),
fFiducialCellCut(0x0),
fMeanV0A(1.),
fMeanV0C(1.),
fPhiCut(TMath::Pi()-0.6),
fRandom(0)
{
   //Constructor

   for(Int_t i=0; i<2; i++) fNClusters[i] = 0;
   for(Int_t i=0; i<8; i++) fRingMultV0[i] = 0;

   for(Int_t i=0; i<5; i++){
      fZNCtower[i] = 0;
      fZPCtower[i] = 0;
      fZNAtower[i] = 0;
      fZPAtower[i] = 0;
      fZNCtowerLG[i] = 0;
      fZPCtowerLG[i] = 0;
      fZNAtowerLG[i] = 0;
      fZPAtowerLG[i] = 0;
   }
  
   //arrays number of triggers
   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i] = 0;
      fJetChTT[i]  = 0;
      fClusterTT[i]  = 0;

      fhMultTTHinMB[i] = 0x0;   
      fhMultTTHinHM[i] = 0x0;   
      fhMultTTJinMB[i] = 0x0;  
      fhMultTTJinHM[i] = 0x0;  
      fhMultTTCinMB[i] = 0x0;  
      fhMultTTCinHM[i] = 0x0;  
      fhMultTTCinGA[i] = 0x0; 

      fhTTHinMB[i] = 0x0;
      fhTTHinHM[i] = 0x0;
      fhTTCinMB[i] = 0x0;
      fhTTCinHM[i] = 0x0;
      fhTTCinGA[i] = 0x0;

      fhDrecoilTTHinMB[i] = 0x0;
      fhDrecoilTTHinHM[i] = 0x0;
      fhDrecoilTTCinMB[i] = 0x0;
      fhDrecoilTTCinHM[i] = 0x0;
      fhDrecoilTTCinGA[i] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
      fJetChTTLowPt[i]=-1;
      fJetChTTHighPt[i]=-1;
      fClusterTTLowPt[i]=-1;
      fClusterTTHighPt[i]=-1;

      fhV0AvsV0CTTH[i] = 0x0; 
      fhV0AvsV0CTTJ[i] = 0x0;
      fhV0AvsV0CTTCinMB[i] = 0x0;  
      fhV0AvsV0CTTCinGA[i] = 0x0; 
 
   }
 
   for(Int_t iv=0; iv<fkVtx;iv++){
      fhVertex[iv]=0x0;
      for(Int_t i=0; i<fkTTbins;i++){
         fhVertexTTH[iv][i]=0x0;
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){
      fhCentralityMB[ic] = 0x0;
      fhSignalMB[ic] = 0x0; 
      fhSignalHM[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhCentralityTTH[ic][i] = 0x0;
         fhCentralityTTJ[ic][i] = 0x0;
         fhCentralityTTCinMB[ic][i] = 0x0;
         fhCentralityTTCinGA[ic][i] = 0x0;

         fhSignalTTHinMB[ic][i] = 0x0;
         fhSignalTTHinHM[ic][i] = 0x0;
         fhSignalTTJinMB[ic][i] = 0x0;
         fhSignalTTJinHM[ic][i] = 0x0;
         fhSignalTTCinMB[ic][i] = 0x0;
         fhSignalTTCinHM[ic][i] = 0x0;
         fhSignalTTCinGA[ic][i] = 0x0;
      }
   }

   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMB[i]=0x0;
      fhRhoTTHinHM[i]=0x0;   
      fhRhoTTJinMB[i]=0x0;  
      fhRhoTTJinHM[i]=0x0; 
      fhRhoTTCinMB[i]=0x0; 
      fhRhoTTCinHM[i]=0x0;  
      fhRhoTTCinGA[i]=0x0; 
   }

   sprintf(fTrigClass,"%s","");
   //inclusive pT spectrum times the boost function

   fFiducialCellCut = new AliEMCALRecoUtils();

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1; 

      fTTC[i].resize(0); 
      fTTH[i].resize(0); 
      fTTJ[i].resize(0); 
   }

   DefineOutput(1, TList::Class());
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

//_____________________________________________________________________________________
AliAnalysisTaskEA*  AliAnalysisTaskEA::AddTaskEA(
  Int_t       system,  
  const char* jetarrayname, 
  const char* jetarraynameMC, 
  const char* trackarrayname, 
  const char* mcpariclearrayname, 
  const char* clusterarrayname, 
  const char* rhoname, 
  const char* mcrhoname, 
  Double_t    jetRadius, 
  UInt_t      trigger, 
  Int_t       isMC, 
  Double_t    trackEtaWindow,
  Bool_t      useVertexCut,
  Bool_t      usePileUpCut, 
  Double_t    acut,
  Double_t    emcaltofcut, 
  const char* suffix 
){
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================

   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);


   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskEA.cxx", "No analysis manager to connect to.");
      return NULL;
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString myContName("");
   myContName = Form("JetAnalysisR%02d", TMath::Nint(jetRadius*10));
   myContName.Append(suffix);


   AliAnalysisTaskEA *task = new AliAnalysisTaskEA(myContName.Data());

   if(isMC){  //for PYTHIA
      task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer    *trackCont      = 0x0; //detector level track container 
   AliParticleContainer *trackContTrue  = 0x0; //mc particle container
   AliClusterContainer  *clusterCont    = 0x0; //detector level track container 

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks 
   trackCont->SetMinPt(0.15);
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
   trackCont->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
   trackCont->SetAODFilterBits((1 << 8) | (1 << 9));

   if(isMC){
      trackContTrue = task->AddMCParticleContainer(mcpariclearrayname); //particle level MC particles   
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(0.15);
      trackContTrue->SetEtaLimits(-trackEtaWindow,trackEtaWindow);
   }

   clusterCont = task->AddClusterContainer(clusterarrayname);  //detector level tracks 
   clusterCont->SetMinPt(0.3);
   clusterCont->SetExoticCut(1);
   clusterCont->SetClusTimeCut(0, emcaltofcut);

   //   clusterCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
 
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //jet container with detector level tracks
   AliJetContainer *jetContTrue   = 0x0; //jet container with mc particles

   jetContRec   = task->AddJetContainer(jetarrayname,"TPC",jetRadius);

   if(jetContRec) { //DETECTOR LEVEL JET
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);
      jetContRec->SetMinPt(0.150);
      jetContRec->SetMaxTrackPt(1000);
      jetContRec->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
    }

    if(isMC){
      //AKT JETS PARTICLE LEVEL
      jetContTrue = task->AddJetContainer(jetarraynameMC,"TPC",jetRadius);
      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }

   // #### Task configuration 
   task->SetMC(isMC);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetAcceptanceWindows(trackEtaWindow);
   task->SelectCollisionCandidates(trigger);
   task->SetExternalRhoTaskName(rhoname);
   task->SetExternalRhoTaskNameMC(mcrhoname);
   task->SetTrackContainerName(trackarrayname);
   task->SetSystem(system);

 
   task->SetMCParticleContainerName(mcpariclearrayname);
   task->SetClusterContainerName(clusterarrayname);
   task->SetJetContainerName(jetarrayname);
   task->SetMCJetContainerName(jetarraynameMC);
   if(system!=AliAnalysisTaskEA::kpp){
     task->SetUseNewCentralityEstimation(kTRUE);  //CENTRALITY
   }else{
     task->SetUseNewCentralityEstimation(kFALSE);  //CENTRALITY
   }

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));


   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedGATrigger(){
  //EG1 high EMCAL trigger
  TString trigger = fInputEvent->GetFiredTriggerClasses();
  bool passedGammaTrigger = kFALSE;
  //EG1 high EMCAL trigger, EG2 low EMCAL trigger, DG1 high DCAL trigger, DG2 low emcal trigger 
  if(trigger.Contains("EG1")){
     passedGammaTrigger = kTRUE;
  }
  return passedGammaTrigger;

}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedMinBiasTrigger(){
  //minimum bias trigger
  TString trigger = fInputEvent->GetFiredTriggerClasses();
  bool passedTrigger = kFALSE;
  if(trigger.Contains("INT7")){
     passedTrigger = kTRUE;
  }
  return passedTrigger;

}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedHighMultTrigger(){
  //high multiplicity V0M trigger
  TString trigger = fInputEvent->GetFiredTriggerClasses();
  bool passedTrigger = kFALSE;
  if(trigger.Contains("HMV0M")){
     passedTrigger = kTRUE;
  }
  return passedTrigger;

}



//_____________________________________________________________________________________
Double_t AliAnalysisTaskEA::GetExternalRho(Bool_t isMC){

   // Get rho from event using CMS approach
   AliRhoParameter* rho = NULL;
   TString rhoname = (!isMC) ? fRhoTaskName : fRhoTaskNameMC;
   if(!rhoname.IsNull()){
      rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(rhoname.Data()));
      if (!rho) {
        //AliWarningF(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), rhoname.Data());
        return 0.;
      }
   }else{
      //AliWarningF(MSGWARNING("No %s Rho task name provided"), (!isMC ? "DATA" : "MC"));
      return 0.;
   }
   
   return rho->GetVal();
}
//________________________________________________________________________

Bool_t AliAnalysisTaskEA::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION RECONSTRUCTED DATA

   if(!event) return kFALSE;

   //incomplete DAQ events rejection Run2 data 2015 
   // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
   Bool_t bIncompleteDAQ = event->IsIncompleteDAQ();
   if(bIncompleteDAQ){

      fHistEvtSelection->Fill(1.5); // count events with incomplete DAQ
      return kFALSE;
   }
   //___________________________________________________
   //TEST PILE UP
   if(fUsePileUpCut){
      if(!fHelperClass || fHelperClass->IsPileUpEvent(event)){ 
         fHistEvtSelection->Fill(2.5); //count events rejected by pileup
         return kFALSE;
      }

      if(!fHelperClass || fHelperClass->IsSPDClusterVsTrackletBG(event)){
         fHistEvtSelection->Fill(2.5); //count events rejected by pileup
         return kFALSE;
      }
   }
   //BEFORE VERTEX CUT
   fhVertexZall->Fill(event->GetPrimaryVertex()->GetZ()); 
   //___________________________________________________
   //VERTEX CUT

   if(fUseDefaultVertexCut){
      if(!fHelperClass || !fHelperClass->IsVertexSelected2013pA(event)){  //??? USE THIS OR SOMETHING ELSE
         fHistEvtSelection->Fill(3.5); //count events rejected by vertex cut 
         return kFALSE;
      }
   }

   if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > fZVertexCut){
      fHistEvtSelection->Fill(3.5); //count events rejected by vertex cut 
      return kFALSE;
   }
   
   //___________________________________________________
   //AFTER VERTEX CUT
   fhVertexZ->Fill(event->GetPrimaryVertex()->GetZ()); 


  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEA::IsTrackInAcceptance(AliVParticle* track, Bool_t isGen){
   // Check if the track pt and eta range 
   if(!track) return kFALSE;

   if(isGen == kPartLevel){ //particle level MC:   select charged physical primary tracks 
      //Apply only for kine level or MC containers   
      if(!track->Charge()) return kFALSE;
      if(!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()) return kFALSE;
   }
   if(TMath::Abs(track->Eta()) < fTrackEtaWindow){ //APPLY TRACK ETA CUT
      if(track->Pt() > fMinTrackPt){   //APPLY TRACK PT CUT
         return kTRUE;
      }
   }
   return kFALSE;
}
//________________________________________________________________________
void AliAnalysisTaskEA::ExecOnceLocal(){
   // Initialization of jet containers done in  AliAnalysisTaskEmcalJet::ExecOnce()
   //Read arrays of jets and tracks
   fInitializedLocal = kTRUE; 

   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm


   return;
}

//________________________________________________________________________
Int_t  AliAnalysisTaskEA::GetMaxDistanceFromBorder(AliVCluster* cluster){
   //Distance o
   Int_t max = 0;
   
   for(Int_t n=0; n<6; n++){
      fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(n);
      if(fFiducialCellCut->CheckCellFiducialRegion(fGeom, cluster, fCaloCells)){
         max = n;
      }else{
         break;
      }
   }

   return max;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEA::FinalClusterCuts(AliVCluster* cluster){

   //General QA. 
   if( cluster->GetNCells() < 1) return kFALSE;
   
   Int_t disToBad = cluster->GetDistanceToBadChannel();
   if(-1<disToBad && disToBad<2) return kFALSE;
   
   Int_t disToBorder = GetMaxDistanceFromBorder(cluster);
   if(-1<disToBorder && disToBorder<1) return kFALSE;

   Double_t lambda0 =  cluster->GetM02();  //avoid merged clusters   
   if(lambda0 > 0.4) return kFALSE;

   return kTRUE;
}

//________________________________________________________________________
/*
std::string AliAnalysisTaskEA::MatchTrigger(const std::string &triggerstring){
  auto triggerclasses = PWG::EMCAL::Triggerinfo::DecodeTriggerString(triggerstring);
  std::string result;
  for(const auto &t : triggerclasses) {
    // Use CENT cluster for downscaling
    if(t.Triggercluster() != "CENT") continue;
    if(t.Triggerclass().find(fTriggerSelectionString.Data()) == std::string::npos) continue;
    result = t.ExpandClassName();
    break;
  }
  return result;
}*/

//________________________________________________________________________
Bool_t AliAnalysisTaskEA::FillHistograms(){  
   // executed in each event 
   //called in AliAnalysisTaskEmcal::UserExec(Option_t *)
   //   Analyze the event and Fill histograms

   if(!InputEvent()){
      AliError("??? Event pointer == 0 ???");
      return kFALSE;
   }

   //Execute only once:  Get tracks, jets from arrays if not already given 
   if(!fInitializedLocal) ExecOnceLocal(); 

   Double_t jetPtcorr;

   //_________________________________________________________________
   // EVENT SELECTION
   fHistEvtSelection->Fill(0.5); //Count input event

   //Check Reconstructed event vertex and pileup
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec


   fIsMinBiasTrig = kFALSE; //Minimum bias event flag
   if(PassedMinBiasTrigger()){
      fIsMinBiasTrig = kTRUE;
      fHistEvtSelection->Fill(4.5); //Count Accepted input event
   }
  
   fIsEmcalTrig = kFALSE; //EMCAL triggered event flag
   Double_t weight = 1.;
   if(PassedGATrigger()){
     
      fIsEmcalTrig = kTRUE; 
      fHistEvtSelection->Fill(5.5); //Count Accepted input event 

      //read downscaling factor: code from  AliAnalysisTaskEmcalJetEnergySpectrum.cxx
      //weight = 1./PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(MatchTrigger(fInputEvent->GetFiredTriggerClasses().Data()));
      //cout<<"FK:WORK ON THE CODE INCLUDE THE WEIGHT "<<weight<<endl;
   }

   fIsHighMultTrig = kFALSE; //high multiplicity trigger flag
   if(PassedHighMultTrigger()){
      fIsHighMultTrig = kTRUE; //Count Accepted input event
      fHistEvtSelection->Fill(6.5); //Count Accepted input event 
   }


   if(!fIsEmcalTrig && !fIsMinBiasTrig && !fIsHighMultTrig)  return kFALSE; //post data is in UserExec

   
   // END EVENT SELECTION
   //_________________________________________________________________
   //_________________________________________________________________
   //                EVENT PROPERTIES   

   for(int ir=0; ir<8; ir++) fRingMultV0[ir]=0.;

   // ***** Trigger selection
   TString triggerClass = InputEvent()->GetFiredTriggerClasses();
   sprintf(fTrigClass,"%s",triggerClass.Data());


   if(fSystem!=AliAnalysisTaskEA::kpp){ 
     fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
     if(fMultSelection){  
         fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
         fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
         fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
         fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
         fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
      }else{
         fCentralityV0A = -1; 
         fCentralityV0C = -1;
         fCentralityCL1 = -1;
         fCentralityZNA = -1;
         fCentralityZNC = -1;
      }
   }

   const AliVVertex *vertex = InputEvent()->GetPrimaryVertexSPD();
   if(vertex){ 
      fxVertex = vertex->GetX();
      fyVertex = vertex->GetY();
      fzVertex = vertex->GetZ();
      if(vertex->IsFromVertexer3D()) fVertexer3d = kTRUE;
      else fVertexer3d = kFALSE;
   }else{
      fxVertex = 9999.;
      fyVertex = 9999.;
      fzVertex = 9999.;
      fVertexer3d = kFALSE;
   }

   const AliVMultiplicity *mult = InputEvent()->GetMultiplicity();
   if(mult){
      fNTracklets = mult->GetNumberOfTracklets();

      for(Int_t ilay=0; ilay<2; ilay++){
         fNClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
      }
   }else{
      fNTracklets = -9999;
      for(Int_t ilay=0; ilay<2; ilay++){
         fNClusters[ilay] = -9999; 
      }
   }

  

   AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
   if(vzeroAOD){
      fMultV0A = vzeroAOD->GetMTotV0A();
      fMultV0C = vzeroAOD->GetMTotV0C();
      fMultV0Anorm = fMultV0A/fMeanV0A;
      fMultV0Cnorm = fMultV0C/fMeanV0C;
      fMultV0AV0Cnorm = fMultV0Anorm + fMultV0Cnorm;

      fIsV0ATriggered = vzeroAOD->GetV0ADecision();
      fIsV0CTriggered = vzeroAOD->GetV0CDecision();
      
      for(Int_t iRing = 0; iRing < 8; ++iRing){
         for(Int_t i = 0; i < 8; ++i){
            fRingMultV0[iRing] += vzeroAOD->GetMultiplicity(8*iRing+i);
         }
      }
   }else{
      fMultV0A = -1; 
      fMultV0C = -1; 
      fIsV0ATriggered = kFALSE; 
      fIsV0CTriggered = kFALSE; 
      
      for(Int_t iRing = 0; iRing < 8; ++iRing){
         for(Int_t i = 0; i < 8; ++i){
            fRingMultV0[iRing] += 0; 
         }
      }
  }


   AliAODZDC *aodZDC =dynamic_cast<AliAODZDC*> (InputEvent()->GetZDCData());
   if(aodZDC){ 

      fZEM1Energy = (Float_t) (aodZDC->GetZEM1Energy());
      fZEM2Energy = (Float_t) (aodZDC->GetZEM2Energy());
      
      const Double_t* towZNC = aodZDC->GetZNCTowerEnergy();
      const Double_t* towZPC = aodZDC->GetZPCTowerEnergy();
      const Double_t* towZNA = aodZDC->GetZNATowerEnergy();
      const Double_t* towZPA = aodZDC->GetZPATowerEnergy();
      //
      const Double_t* towZNCLG = aodZDC->GetZNCTowerEnergyLR();
      const Double_t* towZPCLG = aodZDC->GetZPCTowerEnergyLR();
      const Double_t* towZNALG = aodZDC->GetZNATowerEnergyLR();
      const Double_t* towZPALG = aodZDC->GetZPATowerEnergyLR();
      //
      for(Int_t it=0; it<5; it++){
         fZNCtower[it] = (Float_t) (towZNC[it]);
         fZPCtower[it] = (Float_t) (towZPC[it]);
         fZNAtower[it] = (Float_t) (towZNA[it]);
         fZPAtower[it] = (Float_t) (towZPA[it]);
         fZNCtowerLG[it] = (Float_t) (towZNCLG[it]);
         fZPCtowerLG[it] = (Float_t) (towZPCLG[it]);
         fZNAtowerLG[it] = (Float_t) (towZNALG[it]);
         fZPAtowerLG[it] = (Float_t) (towZPALG[it]);
      }
   }else{
      fZEM1Energy = -1; 
      fZEM2Energy = -1; 
       for(Int_t it=0; it<5; it++){
         fZNCtower[it] = -1;
         fZPCtower[it] = -1; 
         fZNAtower[it] = -1;
         fZPAtower[it] = -1; 
         fZNCtowerLG[it] = -1;
         fZPCtowerLG[it] = -1; 
         fZNAtowerLG[it] = -1;
         fZPAtowerLG[it] = -1; 
      }
   }




   Double_t rho = GetExternalRho(kDetLevel); //estimated backround pt density
   Double_t rhoMC = 0.;

   //_________________________________________________________________
   //                    JET+TRACK CONTAINERS

   AliEmcalJet  *jet = NULL;  //jet pointer real jet
   AliEmcalJet  *jetMC = NULL;  //jet pointer real jet
   AliVParticle *track = NULL; //jet constituent
   AliVParticle *mcParticle = NULL; //mc particle

   //_________________________________________________________
   //READ  TRACK AND JET CONTAINERS
   //Container operations   http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerIterateTechniques

   fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(fMyTrackContainerName.Data())); //reconstructed particle container 
   fJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyJetContainerName.Data())); //AKT jet



   if(fMC){
      fParticleContainerPartLevel = GetParticleContainer(fMyParticleContainerName.Data()); //reconstructed particle container 
      fJetContainerPartLevel      = GetJetContainer(fMyJetParticleContainerName.Data()); //reconstructed particle container 

      rhoMC = GetExternalRho(kPartLevel); //estimated backround pt density
   }
   //___________________________________________
   //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX
   if(fMC){
      Bool_t bRecPrim = kFALSE;

      if(fParticleContainerPartLevel){

         //pT spectrum of particle level physical primary mc particles
         for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
            mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!mcParticle)  continue; 

            if(IsTrackInAcceptance(mcParticle, kPartLevel)){
               fhPtTrkTruePrimGen->Fill(mcParticle->Pt(), mcParticle->Eta());
            }
         }
      }

      //pT spectrum of detector level physical primary tracks and secondary tracks
      if(fTrkContainerDetLevel && fParticleContainerPartLevel){
         for(auto trkIterator : fTrkContainerDetLevel->accepted_momentum() ){
            track = trkIterator.second;  // Get the pointer to mc particle object
            if(!track)  continue; 

            if(!IsTrackInAcceptance(track, kDetLevel)) continue; //reconstructed level tracks
            bRecPrim = kFALSE; //not yet matched to generator level physical primary

            for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
               mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
               if(!mcParticle)  continue; 

               if(IsTrackInAcceptance(mcParticle, kPartLevel)){
                  if(TMath::Abs(track->GetLabel()) == TMath::Abs(mcParticle->GetLabel())){
                     //has the same label as reconstr track
                  
                     bRecPrim = kTRUE;
                     fhPtTrkTruePrimRec->Fill(mcParticle->Pt(), mcParticle->Eta()); //this is well recontr phys primary
                     break;
                  }//same label with rec particle
               }
            }//loop over gen tracks
            if(!bRecPrim){
               fhPtTrkSecOrFakeRec->Fill(track->Pt(), track->Eta()); //matchnig to phys primary not found, this is fake or second.
            }
         }
      }

      //__________________________________________________________
      //  FILL JET RESPONSE MATRIX

      Double_t jetPtCorrPart = 0.; //particle level jet pt corrected for rho
      Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho

      //Response matrix normalization - spectrum of all generator level jets in acceptance
      if(fJetContainerPartLevel){
         for(auto jetPartIterator : fJetContainerPartLevel->accepted_momentum() ){
            jetMC = jetPartIterator.second;  // Get the pointer to mc particle object
            if(!jetMC)  continue; 

            jetPtCorrPart = jetMC->Pt() - jetMC->Area()*rhoMC;

            fhJetPtPartLevelCorr->Fill(jetPtCorrPart);
            fhJetPtPartLevelZero->Fill(jetMC->Pt());
         }
      }

      //1) Find closest particle level and detector level jets
      //2) Get momentum shift due to fake tracks
      if(fJetContainerDetLevel){
         Double_t  sumAllTrackPtInJet   = 0.;
         Double_t  sumFakeTrackPtInJet  = 0.;


         for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
            jet = jetIterator.second;  // Get the pointer to jet object
            if(!jet)  continue; 
 
            //Get momentum shift due to fake tracks
            sumAllTrackPtInJet  = 0.;
            sumFakeTrackPtInJet  = 0.;

            for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
               track = static_cast<AliVParticle*> (jet->TrackAt(iq, fTrkContainerDetLevel->GetArray()));
               if(!track) continue;
               bRecPrim = kFALSE; //not yet matched to generator level physical primary

               for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
                  mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
                  if(!mcParticle)  continue; 

                  if(IsTrackInAcceptance(mcParticle, kPartLevel)){
                     if(TMath::Abs(track->GetLabel()) == TMath::Abs(mcParticle->GetLabel())){
                        bRecPrim = kTRUE;
                        break;
                     }
                  }
               }
               if(!(bRecPrim && mcParticle &&  mcParticle->IsPhysicalPrimary())){ //this is a fake track

                  sumFakeTrackPtInJet += track->Pt();
               }
               sumAllTrackPtInJet += track->Pt();
            }

            if(sumAllTrackPtInJet>0){
               jetPtCorrDet = jet->Pt() - jet->Area()*rho;
               fhFractionOfSecInJet->Fill(jetPtCorrDet, sumFakeTrackPtInJet/jetPtCorrDet); 
            }

            //Fill Response matrix
            jetMC =  jet->ClosestJet();
            if(!jetMC) continue;
            if(jetMC->Pt()<1e-3) continue;

            jetPtCorrPart =  jetMC->Pt() - jetMC->Area()*rhoMC; 
            jetPtCorrDet  =  jet->Pt() - jet->Area()*rho; 
 
            fhJetPtPartLevelVsJetPtDetLevelCorr->Fill(jetPtCorrDet,jetPtCorrPart); //response matrix
            fhJetPtPartLevelVsJetPtDetLevelZero->Fill(jet->Pt(),jetMC->Pt()); //response matrix

            if(jetPtCorrPart>0){
               fhJetPtResolutionVsPtPartLevel->Fill(jetPtCorrPart,(jetPtCorrDet-jetPtCorrPart)/jetPtCorrPart); //jet pT resolution
            }
         }
      }
   
      //TO BE DONE
      //1) Delta pT for MC particle level events
      //2) Signal for MC particle level events



   }
   //___________________________________________
   //    INCLUSIVE EVENTS (WITHOUT TT REQUIREMENT)

    if(fIsMinBiasTrig){ 
       fhVertex[0]->Fill(fxVertex);
       fhVertex[1]->Fill(fyVertex);
       fhVertex[2]->Fill(fzVertex);

       fhRhoMB->Fill(rho);

       if(fSystem!=AliAnalysisTaskEA::kpp){ 
          fhCentralityMB[fkV0A]->Fill(fCentralityV0A); 
          fhCentralityMB[fkV0C]->Fill(fCentralityV0C); 
          fhCentralityMB[fkSPD]->Fill(fCentralityCL1);
          fhCentralityMB[fkZNA]->Fill(fCentralityZNA);
          fhCentralityMB[fkZNC]->Fill(fCentralityZNC);
       }

       fhSignalMB[fkV0A]->Fill(fMultV0A);
       fhSignalMB[fkV0C]->Fill(fMultV0C);
       fhSignalMB[fkSPD]->Fill(fNTracklets); 
       fhSignalMB[fkZNA]->Fill(fZNAtower[0]); 
       fhSignalMB[fkZNC]->Fill(fZNCtower[0]);
       fhSignalMB[fkV0M]->Fill(fMultV0A+fMultV0C);
       fhSignalMB[fkV0Mnorm]->Fill(fMultV0AV0Cnorm);

       fhV0AvsV0C->Fill(fMultV0C, fMultV0A);
       fhV0AvsSPD->Fill(fNTracklets, fMultV0A);
       fhV0CvsSPD->Fill(fNTracklets,fMultV0C);

    }
 
    if(fIsHighMultTrig){ 
       fhRhoHM->Fill(rho);

       fhSignalHM[fkV0A]->Fill(fMultV0A);
       fhSignalHM[fkV0C]->Fill(fMultV0C);
       fhSignalHM[fkSPD]->Fill(fNTracklets); 
       fhSignalHM[fkZNA]->Fill(fZNAtower[0]); 
       fhSignalHM[fkZNC]->Fill(fZNCtower[0]);
       fhSignalHM[fkV0M]->Fill(fMultV0A+fMultV0C);
       fhSignalHM[fkV0Mnorm]->Fill(fMultV0AV0Cnorm);
   }

   //_________________________________________________________
   //                 TT

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1; 

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);
   }

   TLorentzVector myTT;
   Int_t idx;
   //_________________________________________________________
   //LOOP OVER EMCAL CLUSTERS
   TLorentzVector ph;
   for(Int_t i=0; i<fnClusterTTBins; i++){
      fClusterTT[i] = 0;
   }

   if(fMyClusterContainerName.Data()){
      fClusterContainerDetLevel =  static_cast<AliClusterContainer*> ( GetClusterContainer(fMyClusterContainerName.Data()));
  
 
      for(auto cluster: fClusterContainerDetLevel->accepted()){
         fClusterContainerDetLevel->GetMomentum(ph, cluster);

         if(!FinalClusterCuts(cluster)) continue;

         if(fIsMinBiasTrig){
            //fill some histograms for detector level tracks 
            fhClusterPhiInclMB->Fill(ph.Pt(), ph.Phi());
            fhClusterEtaInclMB->Fill(ph.Pt(), ph.Eta());
         }

         if(fIsEmcalTrig){
            fhClusterPhiInclGA->Fill(ph.Pt(), ph.Phi());
            fhClusterEtaInclGA->Fill(ph.Pt(), ph.Eta());
         }

         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            if(fClusterTTLowPt[igg] < ph.Pt() && ph.Pt() < fClusterTTHighPt[igg]){
               myTT.SetPtEtaPhiM(ph.Pt(),ph.Eta(),ph.Phi(),0.); 
               fTTC[igg].push_back(myTT);
               fClusterTT[igg]++;   // there was a high pt emcal cluster 
            } 
         }
      }

     //chose trigger emcal cluster TT 
     for(Int_t igg=0; igg<fnClusterTTBins; igg++){
        if(fClusterTT[igg]>0){
           fIndexTTC[igg] = fRandom->Integer(fClusterTT[igg]);
        }
     }


 
      if(fIsMinBiasTrig){ 
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         
            fhMultTTCinMB[igg]->Fill(fClusterTT[igg]); 
            
            if(!fClusterTT[igg]) continue;

            fhRhoTTCinMB[igg]->Fill(rho); 
            
            if(fSystem!=AliAnalysisTaskEA::kpp){ 
               fhCentralityTTCinMB[fkV0A][igg]->Fill(fCentralityV0A); 
               fhCentralityTTCinMB[fkV0C][igg]->Fill(fCentralityV0C); 
               fhCentralityTTCinMB[fkSPD][igg]->Fill(fCentralityCL1);
               fhCentralityTTCinMB[fkZNA][igg]->Fill(fCentralityZNA);
               fhCentralityTTCinMB[fkZNC][igg]->Fill(fCentralityZNC);
            }
 
            fhSignalTTCinMB[fkV0A][igg]->Fill(fMultV0A);
            fhSignalTTCinMB[fkV0C][igg]->Fill(fMultV0C);
            fhSignalTTCinMB[fkSPD][igg]->Fill(fNTracklets); 
            fhSignalTTCinMB[fkZNA][igg]->Fill(fZNAtower[0]); 
            fhSignalTTCinMB[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinMB[fkV0M][igg]->Fill(fMultV0A+fMultV0C);
            fhSignalTTCinMB[fkV0Mnorm][igg]->Fill(fMultV0AV0Cnorm);


            fhV0AvsV0CTTCinMB[igg]->Fill(fMultV0C, fMultV0A);


            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){ 
               fhTTCinMB[igg]->Fill(fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhDrecoilTTCinMB[igg]->Fill(jetPtcorr);
                  }
               } 
            }
         }
      }

      if(fIsHighMultTrig){ 
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         
            fhMultTTCinHM[igg]->Fill(fClusterTT[igg]); 
            
            if(!fClusterTT[igg]) continue;

            fhRhoTTCinHM[igg]->Fill(rho); 
 
            fhSignalTTCinHM[fkV0A][igg]->Fill(fMultV0A);
            fhSignalTTCinHM[fkV0C][igg]->Fill(fMultV0C);
            fhSignalTTCinHM[fkSPD][igg]->Fill(fNTracklets); 
            fhSignalTTCinHM[fkZNA][igg]->Fill(fZNAtower[0]); 
            fhSignalTTCinHM[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinHM[fkV0M][igg]->Fill(fMultV0A+fMultV0C);
            fhSignalTTCinHM[fkV0Mnorm][igg]->Fill(fMultV0AV0Cnorm);

            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){ 
               fhTTCinHM[igg]->Fill(fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhDrecoilTTCinHM[igg]->Fill(jetPtcorr);
                  }
               } 
            }
         }
      }

      if(fIsEmcalTrig){ 
     
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         
            fhMultTTCinGA[igg]->Fill(fClusterTT[igg]); 
            
            if(!fClusterTT[igg]) continue;
            
            fhRhoTTCinGA[igg]->Fill(rho); 
            
            if(fSystem!=AliAnalysisTaskEA::kpp){ 
               fhCentralityTTCinGA[fkV0A][igg]->Fill(fCentralityV0A); 
               fhCentralityTTCinGA[fkV0C][igg]->Fill(fCentralityV0C); 
               fhCentralityTTCinGA[fkSPD][igg]->Fill(fCentralityCL1);
               fhCentralityTTCinGA[fkZNA][igg]->Fill(fCentralityZNA);
               fhCentralityTTCinGA[fkZNC][igg]->Fill(fCentralityZNC);
            }

            fhSignalTTCinGA[fkV0A][igg]->Fill(fMultV0A);
            fhSignalTTCinGA[fkV0C][igg]->Fill(fMultV0C);
            fhSignalTTCinGA[fkSPD][igg]->Fill(fNTracklets); 
            fhSignalTTCinGA[fkZNA][igg]->Fill(fZNAtower[0]); 
            fhSignalTTCinGA[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinGA[fkV0M][igg]->Fill(fMultV0A+fMultV0C);
            fhSignalTTCinGA[fkV0Mnorm][igg]->Fill(fMultV0AV0Cnorm);
            
            fhV0AvsV0CTTCinGA[igg]->Fill(fMultV0C, fMultV0A);
 
            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){ 
               fhTTCinGA[igg]->Fill(fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhDrecoilTTCinGA[igg]->Fill(jetPtcorr);
                  }
               } 
            }//trigger exists
         }//loop over TT bins
      }//emcale trigger 
   }//cluster container   

 



   //_________________________________________________________
   //LOOP OVER TRACKS DETECTOR LEVEL + SEARCH FOR HIGH PT HADRON TRIGGER 

   for(Int_t i=0; i<fnHadronTTBins; i++){
      fHadronTT[i] = 0;
   }

   Double_t xyz[50];
   Double_t pxpypz[50];
   Double_t cv[21];
   Int_t label, labelMC;                                     //AID
   Bool_t labelfound=0;
   AliAODMCParticle* particleMC = NULL;             //AID 
   AliAODTrack *trackAOD=NULL ;
 
   for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      track = trackIterator.second;  // Get the full track
      if(!track) continue;

      if(IsTrackInAcceptance(track, kDetLevel)){  


         if(fIsMinBiasTrig){
            //fill some histograms for detector level tracks 
            fhTrackPhiIncl->Fill(track->Pt(), track->Phi());
            fhTrackEtaIncl->Fill(track->Pt(), track->Eta());

            //get sigma pT / pT  
            //Taken from AliEMCalTriggerExtraCuts::CalculateTPCTrackLength
            memset(cv, 0, sizeof(Double_t) * 21); //cleanup arrays
            memset(pxpypz, 0, sizeof(Double_t) * 50);
            memset(xyz, 0, sizeof(Double_t) * 50);

            trackAOD = static_cast <AliAODTrack*>( track);
            if(trackAOD){
               trackAOD->GetXYZ(xyz);
               trackAOD->GetPxPyPz(pxpypz);
               trackAOD->GetCovarianceXYZPxPyPz(cv);
               
               AliExternalTrackParam  par(xyz, pxpypz, cv, trackAOD->Charge());
               fhSigmaPtOverPtVsPt->Fill(trackAOD->Pt(), TMath::Abs(sqrt(par.GetSigma1Pt2())/par.GetSigned1Pt()));
               
               if(trackAOD->Charge()<0){
                  fhOneOverPtVsPhiNeg->Fill(trackAOD->Phi(), 1.0/trackAOD->Pt());
               }else{
                  fhOneOverPtVsPhiPos->Fill(trackAOD->Phi(), 1.0/trackAOD->Pt());
               }
               
               //DCA distributions
               fhDCAinXVsPt->Fill(trackAOD->Pt(), trackAOD->XAtDCA());
               fhDCAinYVsPt->Fill(trackAOD->Pt(), trackAOD->YAtDCA());
               
               if(fMC){
                  label = TMath::Abs(trackAOD->GetLabel());        //AID
                  
                  particleMC = NULL;
                  labelfound=0;
                  for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
                     particleMC  =  static_cast <AliAODMCParticle*>(  mcPartIterator.second);  // Get the pointer to mc particle object
               
                     labelMC = TMath::Abs(particleMC->GetLabel());
                     if(labelMC==label && label > -1){
                        labelfound=1;
                        break;
                     }
                  }
                  if(labelfound && particleMC && particleMC->IsPhysicalPrimary()){
                     fhDCAinXVsPtPhysPrimary->Fill(trackAOD->Pt(), trackAOD->XAtDCA());
                     fhDCAinYVsPtPhysPrimary->Fill(trackAOD->Pt(), trackAOD->YAtDCA());
                  }else{
                     fhDCAinXVsPtSecondary->Fill(trackAOD->Pt(), trackAOD->XAtDCA());
                     fhDCAinYVsPtSecondary->Fill(trackAOD->Pt(), trackAOD->YAtDCA());
                  }//AID
               }
            }
         }

         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(fHadronTTLowPt[itt] < track->Pt() && track->Pt() < fHadronTTHighPt[itt]){
               myTT.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),0.); 
               fTTH[itt].push_back(myTT);
               fHadronTT[itt]++;   // there was a high pt 
            } 
         }
      } 
   }

   //chose trigger hadron TT
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      if(fHadronTT[itt]>0){
         fIndexTTH[itt] = fRandom->Integer(fHadronTT[itt]);
      }
   }

   if(fIsMinBiasTrig){ 
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      
         fhMultTTHinMB[itt]->Fill(fHadronTT[itt]); 
      
         if(!fHadronTT[itt]) continue;
         fhVertexTTH[0][itt]->Fill(fxVertex);
         fhVertexTTH[1][itt]->Fill(fyVertex);
         fhVertexTTH[2][itt]->Fill(fzVertex);
      
         fhRhoTTHinMB[itt]->Fill(rho); 

         if(fSystem!=AliAnalysisTaskEA::kpp){ 
            fhCentralityTTH[fkV0A][itt]->Fill(fCentralityV0A); 
            fhCentralityTTH[fkV0C][itt]->Fill(fCentralityV0C); 
            fhCentralityTTH[fkSPD][itt]->Fill(fCentralityCL1);
            fhCentralityTTH[fkZNA][itt]->Fill(fCentralityZNA);
            fhCentralityTTH[fkZNC][itt]->Fill(fCentralityZNC);
         }
 
         fhSignalTTHinMB[fkV0A][itt]->Fill(fMultV0A);
         fhSignalTTHinMB[fkV0C][itt]->Fill(fMultV0C);
         fhSignalTTHinMB[fkSPD][itt]->Fill(fNTracklets); 
         fhSignalTTHinMB[fkZNA][itt]->Fill(fZNAtower[0]); 
         fhSignalTTHinMB[fkZNC][itt]->Fill(fZNCtower[0]); 
         fhSignalTTHinMB[fkV0M][itt]->Fill(fMultV0A+fMultV0C);
         fhSignalTTHinMB[fkV0Mnorm][itt]->Fill(fMultV0AV0Cnorm);

         fhV0AvsV0CTTH[itt]->Fill(fMultV0C, fMultV0A);

         //hadron trigger
         idx = fIndexTTH[itt];
         if(idx>-1){ 
            fhTTHinMB[itt]->Fill(fTTH[itt][idx].Pt()); //fill trigger track pT
          
            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue; 

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){     
                  //recoil jet
                  jetPtcorr = jet->Pt() - rho*jet->Area();
                  fhDrecoilTTHinMB[itt]->Fill(jetPtcorr);
               }
            } 
         }
      }
   }

   if(fIsHighMultTrig){ 
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      
         fhMultTTHinHM[itt]->Fill(fHadronTT[itt]); 
      
         if(!fHadronTT[itt]) continue;
      
         fhRhoTTHinHM[itt]->Fill(rho); 
 
         fhSignalTTHinHM[fkV0A][itt]->Fill(fMultV0A);
         fhSignalTTHinHM[fkV0C][itt]->Fill(fMultV0C);
         fhSignalTTHinHM[fkSPD][itt]->Fill(fNTracklets); 
         fhSignalTTHinHM[fkZNA][itt]->Fill(fZNAtower[0]); 
         fhSignalTTHinHM[fkZNC][itt]->Fill(fZNCtower[0]); 
         fhSignalTTHinHM[fkV0M][itt]->Fill(fMultV0A+fMultV0C);
         fhSignalTTHinHM[fkV0Mnorm][itt]->Fill(fMultV0AV0Cnorm);

         //hadron trigger
         idx = fIndexTTH[itt];
         if(idx>-1){ 
            fhTTHinHM[itt]->Fill(fTTH[itt][idx].Pt()); //fill trigger track pT
          
            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue; 

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){     
                  //recoil jet
                  jetPtcorr = jet->Pt() - rho*jet->Area();
                  fhDrecoilTTHinHM[itt]->Fill(jetPtcorr);
               }
            } 
         }

      }
   }

   //_________________________________________________________
   //LOOP OVER JETS  DETECTOR LEVEL
 
   for(Int_t i=0; i<fnJetChTTBins; i++){
      fJetChTT[i] = 0;
   }



   for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      jet = jetIterator.second;  // Get the pointer to jet object
      if(!jet)  continue; 
   
      jetPtcorr = jet->Pt() - rho*jet->Area();
      if(fIsMinBiasTrig){
         //fill some histograms for detector level jets 
         fhJetPhiIncl->Fill(jetPtcorr, jet->Phi());
         fhJetEtaIncl->Fill(jetPtcorr, jet->Eta());
      }

      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         if(fJetChTTLowPt[ijj] < jetPtcorr && jetPtcorr < fJetChTTHighPt[ijj]){
            myTT.SetPtEtaPhiM(jetPtcorr, jet->Eta(), jet->Phi(), 0.); 
            fTTJ[ijj].push_back(myTT);
            fJetChTT[ijj]++;   // there was a high pt jet
         } 
      }
   
      //loop over jet constituents
      //for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
      //   track = (AliVParticle*) (jet->TrackAt(iq, fTrkContainerDetLevel->GetArray()));
         //here one can e.g. analyze jet shapes
 
      //}

       
      //you can also find the closest particle level jet to given detector level
      //the mateching betwe particle and detector level jets is done in Tagger task
      //if(fMC){
      //   jetMC = jet->ClosestJet();
      //}
   }

   //chose trigger emcal cluster TT 
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      if(fJetChTT[ijj]>0){
         fIndexTTJ[ijj] = fRandom->Integer(fJetChTT[ijj]);
      }
   }


   if(fIsMinBiasTrig){ 
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      
         fhMultTTJinMB[ijj]->Fill(fJetChTT[ijj]); 
      
         if(!fJetChTT[ijj]) continue; 
       
         fhRhoTTJinMB[ijj]->Fill(rho);

         if(fSystem!=AliAnalysisTaskEA::kpp){ 
            fhCentralityTTJ[fkV0A][ijj]->Fill(fCentralityV0A); 
            fhCentralityTTJ[fkV0C][ijj]->Fill(fCentralityV0C); 
            fhCentralityTTJ[fkSPD][ijj]->Fill(fCentralityCL1);
            fhCentralityTTJ[fkZNA][ijj]->Fill(fCentralityZNA);
            fhCentralityTTJ[fkZNC][ijj]->Fill(fCentralityZNC);
         } 
      
         fhSignalTTJinMB[fkV0A][ijj]->Fill(fMultV0A);
         fhSignalTTJinMB[fkV0C][ijj]->Fill(fMultV0C);
         fhSignalTTJinMB[fkSPD][ijj]->Fill(fNTracklets); 
         fhSignalTTJinMB[fkZNA][ijj]->Fill(fZNAtower[0]); 
         fhSignalTTJinMB[fkZNC][ijj]->Fill(fZNCtower[0]); 
         fhSignalTTJinMB[fkV0M][ijj]->Fill(fMultV0A+fMultV0C);
         fhSignalTTJinMB[fkV0Mnorm][ijj]->Fill(fMultV0AV0Cnorm);
      
         fhV0AvsV0CTTJ[ijj]->Fill(fMultV0C, fMultV0A);

      }
   }

   if(fIsHighMultTrig){ 
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      
         fhMultTTJinHM[ijj]->Fill(fJetChTT[ijj]); 
      
         if(!fJetChTT[ijj]) continue; 
       
         fhRhoTTJinHM[ijj]->Fill(rho);
      
         fhSignalTTJinHM[fkV0A][ijj]->Fill(fMultV0A);
         fhSignalTTJinHM[fkV0C][ijj]->Fill(fMultV0C);
         fhSignalTTJinHM[fkSPD][ijj]->Fill(fNTracklets); 
         fhSignalTTJinHM[fkZNA][ijj]->Fill(fZNAtower[0]); 
         fhSignalTTJinHM[fkZNC][ijj]->Fill(fZNCtower[0]); 
         fhSignalTTJinHM[fkV0M][ijj]->Fill(fMultV0A+fMultV0C);
         fhSignalTTJinHM[fkV0Mnorm][ijj]->Fill(fMultV0AV0Cnorm);
      }
   }

 
   //___________________________________________________________

   if(fFillTTree){ 
      fCentralityTree->Fill();
   }


   return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEA::Terminate(Option_t *){
   //Treminate 
   PostData(1, fOutput);

   // Mandatory
   fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1)); // '1' refers to the output slot
   if(!fOutput) {
      printf("ERROR: Output list not available\n");
      return;
   }
}

//________________________________________________________________________
AliAnalysisTaskEA::~AliAnalysisTaskEA(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fRandom;
   delete fHelperClass;
 
} 
//________________________________________________________________________
void AliAnalysisTaskEA::UserCreateOutputObjects(){
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.
  //fOutput TList defined in the mother class
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name, object;

   fRandom = new TRandom3(0);
   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 7, 0, 7);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"events IN"); //0-1
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"incomplete DAQ (rejected)"); //1-2
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)"); //2-3
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)"); //3-4
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"MB"); //4-5
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"EMCAL"); //5-6
   fHistEvtSelection->GetXaxis()->SetBinLabel(7,"High Mult"); //6-7


   fOutput->Add(fHistEvtSelection);


   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   fhVertexZall =  new TH1F("fhVertexZall","z vertex without cut",40,-20,20);
   fOutput->Add(fhVertexZall); 
 
   fhVertexZ = new TH1F("fhVertexZ","z vertex",40,-20,20);
   fOutput->Add(fhVertexZ);
 
   //-------------------------

   fhTrackEtaIncl = new TH2F("fhTrackEtaIncl","Eta dist inclusive track vs pT", 50,0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhTrackEtaIncl);

   fhTrackPhiIncl = new TH2F("fhTrackPhiIncl","Azim dist tracks vs pT", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2F*) fhTrackPhiIncl);

   fhJetEtaIncl = new TH2F("fhJetEtaIncl","Eta dist inclusive jets vs pTjet", 150, -20, 130, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhJetEtaIncl);

   fhJetPhiIncl = new TH2F("fhJetPhiIncl","Azim dist jets vs pTjet", 60, -20, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2F*) fhJetPhiIncl);

   fhClusterEtaInclMB = new TH2F("fhClusterEtaInclMB","Eta dist inclusive clusters vs pT", 100, 0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhClusterEtaInclMB);

   fhClusterPhiInclMB = new TH2F("fhClusterPhiInclMB","Azim dist clusters vs pT", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2F*) fhClusterPhiInclMB);

   fhClusterEtaInclGA = new TH2F("fhClusterEtaInclGA","Eta dist inclusive clusters vs pT", 100, 0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhClusterEtaInclGA);

   fhClusterPhiInclGA = new TH2F("fhClusterPhiInclGA","Azim dist clusters vs pT", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2F*) fhClusterPhiInclGA);


   //RHO 
   fhRhoMB = new TH1F("hRho","Rho",1000,0,100);
   fOutput->Add((TH1F*) fhRhoMB); 

    name = Form("hRho_HM");
   fhRhoHM = (TH1F*)  fhRhoMB->Clone(name.Data()); 
   fOutput->Add((TH1F*) fhRhoHM); 


   for(Int_t itt=0; itt<fnHadronTTBins;itt++){
      name = Form("%s_MB_TTH%d_%d", fhRhoMB->GetName(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRhoTTHinMB[itt] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTHinMB[itt]); 
   }
   for(Int_t itt=0; itt<fnHadronTTBins;itt++){
      name = Form("%s_HM_TTH%d_%d", fhRhoMB->GetName(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRhoTTHinHM[itt] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events HM with hadron TT
      fOutput->Add((TH1F*) fhRhoTTHinHM[itt]); 
   }
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("%s_MB_TTJ%d_%d", fhRhoMB->GetName(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhRhoTTJinMB[ijj] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTJinMB[ijj]); 
   } 
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("%s_HM_TTJ%d_%d", fhRhoMB->GetName(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhRhoTTJinHM[ijj] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTJinHM[ijj]); 
   } 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("%s_MB_TTC%d_%d", fhRhoMB->GetName(), fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinMB[igg] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTCinMB[igg]); 
   } 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("%s_HM_TTC%d_%d", fhRhoMB->GetName(), fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinHM[igg] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTCinHM[igg]); 
   } 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("%s_GA_TTC%d_%d", fhRhoMB->GetName(), fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinGA[igg] = (TH1F*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1F*) fhRhoTTCinGA[igg]); 
   } 



   //VERTEX
   for(Int_t iv=0; iv<fkVtx;iv++){
      if(iv==0)       fhVertex[iv] = new TH1D("hVertexX","VertexX",100,-1,1);
      else if(iv==1)  fhVertex[iv] = new TH1D("hVertexY","VertexY",100,-1,1);
      else            fhVertex[iv] = new TH1D("hVertexZ","VertexZ",400,-20,20);
      fOutput->Add((TH1D*) fhVertex[iv]); 
   }

   for(Int_t iv=0; iv<fkVtx;iv++){
      for(Int_t itt=0; itt<fnHadronTTBins;itt++){
         name = Form("%s_TTH%d_%d", fhVertex[iv]->GetName(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhVertexTTH[iv][itt] = (TH1D*) fhVertex[iv]->Clone(name.Data()); 
         fOutput->Add((TH1D*) fhVertexTTH[iv][itt]); 
      }
   }


   TString cest[fkCE] = {"V0A", "V0C", "SPD", "ZNA", "ZNC", "V0M", "V0Mnorm"}; //centrality estimators

   if(fSystem!=AliAnalysisTaskEA::kpp){ //not pp
      for(Int_t ic=0; ic<fkCE;ic++){
         name = Form("hCentrality_MB_%s",cest[ic].Data());
         fhCentralityMB[ic] = new TH1D(name.Data(), name.Data(),101,0,101);
         fOutput->Add((TH1D*) fhCentralityMB[ic]); 
      }
      
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hCentrality_MB_%s_TTH%d_%d",cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
            fhCentralityTTH[ic][itt] = new TH1D(name.Data(), name.Data(),101,0,101);
            fOutput->Add((TH1D*) fhCentralityTTH[ic][itt]); 
         }
      }
      
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
            name = Form("hCentrality_MB_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
            fhCentralityTTJ[ic][ijj] = new TH1D(name.Data(), name.Data(),101,0,101);
            fOutput->Add((TH1D*) fhCentralityTTJ[ic][ijj]); 
         }
      }
      
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
            name = Form("hCentrality_MB_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
            fhCentralityTTCinMB[ic][ijj] = new TH1D(name.Data(), name.Data(),101,0,101);
            fOutput->Add((TH1D*) fhCentralityTTCinMB[ic][ijj]); 
         }
      }
      
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
            name = Form("hCentrality_GA_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
            fhCentralityTTCinGA[ic][ijj] = new TH1D(name.Data(), name.Data(),101,0,101);
            fOutput->Add((TH1D*) fhCentralityTTCinGA[ic][ijj]); 
         }
      }
   }//not pp

   TString signal[]={"multV0A", "multV0C", "nTracklets", "znatower0", "znctower0","multV0M","fhNormSumV0AV0C"};
   Float_t signalL[]={0,0,0,0,0,0,0};
   Float_t signalH[]={1000,1000,500,30000,30000,1200,40};
   Int_t signalN[]={1000,1000,500,100,100,1200,400};

   for(Int_t ic=0; ic<fkCE;ic++){ //MB
      name = Form("hSignal_MB_%s", cest[ic].Data());
      fhSignalMB[ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
      fOutput->Add((TH1D*) fhSignalMB[ic]); 
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM
      name = Form("hSignal_HM_%s", cest[ic].Data());
      fhSignalHM[ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
      fOutput->Add((TH1D*) fhSignalHM[ic]); 
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT hadron
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hSignal_MB_%s_TTH%d_%d", cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhSignalTTHinMB[ic][itt] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTHinMB[ic][itt]); 
      }
   }
 
   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT hadron
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hSignal_HM_%s_TTH%d_%d", cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhSignalTTHinHM[ic][itt] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTHinHM[ic][itt]); 
      }
   }
   
   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT jet
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hSignal_MB_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhSignalTTJinMB[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTJinMB[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT jet 
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hSignal_HM_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhSignalTTJinHM[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTJinHM[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT emcal 
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_MB_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinMB[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTCinMB[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT emcal 
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_HM_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinHM[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTCinHM[ic][ijj]); 
      }
   }
 
   for(Int_t ic=0; ic<fkCE;ic++){ //GA && TT emcal
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_GA_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinGA[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTCinGA[ic][ijj]); 
      }
   }


   name = Form("fhV0AvsV0C_MB");
   fhV0AvsV0C = new TH2F(name.Data(),name.Data(),100,0,1000, 100,0,1000);
   fOutput->Add((TH2F*) fhV0AvsV0C); 

   name = Form("fhV0AvsSPD_MB");
   fhV0AvsSPD = new TH2F(name.Data(),name.Data(),100,0,500, 100,0,500);
   fOutput->Add((TH2F*) fhV0AvsSPD);

   name = Form("fhV0CvsSPD_MB");
   fhV0CvsSPD = new TH2F(name.Data(),name.Data(),100,0,500, 100,0,500);
   fOutput->Add((TH2F*) fhV0CvsSPD);

 
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("fhV0AvsV0C_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhV0AvsV0CTTH[itt] = (TH2F*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTH[itt]->SetTitle(name.Data());
      fOutput->Add((TH2F*) fhV0AvsV0CTTH[itt]); 
   } 
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("fhV0AvsV0C_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhV0AvsV0CTTJ[ijj] =  (TH2F*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTJ[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2F*) fhV0AvsV0CTTJ[ijj]); 
   }
   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("fhV0AvsV0C_MB_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhV0AvsV0CTTCinMB[ijj] = (TH2F*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTCinMB[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2F*) fhV0AvsV0CTTCinMB[ijj]); 
   } 
   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("fhV0AvsV0C_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhV0AvsV0CTTCinGA[ijj] = (TH2F*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTCinGA[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2F*) fhV0AvsV0CTTCinGA[ijj]); 
   } 

   //Trigger track candidate multiplicity
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hMultTT_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhMultTTHinMB[itt] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*)  fhMultTTHinMB[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hMultTT_HM_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhMultTTHinHM[itt] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*)  fhMultTTHinHM[itt]); 
   }

   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hMultTT_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhMultTTJinMB[ijj] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*) fhMultTTJinMB[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hMultTT_HM_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhMultTTJinHM[ijj] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*) fhMultTTJinHM[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_MB_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinMB[ijj] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*) fhMultTTCinMB[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_HM_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinHM[ijj] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*) fhMultTTCinHM[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinGA[ijj] = new TH1I(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1I*) fhMultTTCinGA[ijj]); 
   }

   //Trigger track pT spectrum single inclusive
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinMB[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhTTHinMB[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_HM_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinHM[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhTTHinHM[itt]); 
   }

   //Trigger emcal cluster pT spectrum single inclusive
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_MB_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinMB[igg] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhTTCinMB[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_HM_TTC%d_%d",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinHM[igg] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhTTCinHM[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_GA_TTC%d_%d",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinGA[igg] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhTTCinGA[igg]); 
   }

   //DELTA RECOIL SPECTRA
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //! counter of semi-inclusive hadron trigges  in MB
      name = Form("fhDrecoil_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDrecoilTTHinMB[itt] = new TH1D(name.Data(), name.Data(), 200, -20, 180);            
      fOutput->Add((TH1D*) fhDrecoilTTHinMB[itt]); 
   }
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){         //! counter of semi-inclusive hadron trigges  in HM
      name = Form("fhDrecoil_HM_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDrecoilTTHinHM[itt] = (TH1D*) fhDrecoilTTHinMB[itt]->Clone(name.Data()); 
      fOutput->Add((TH1D*) fhDrecoilTTHinHM[itt]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDrecoil_MB_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDrecoilTTCinMB[igg] = (TH1D*) fhDrecoilTTHinMB[0]->Clone(name.Data()); 
      fOutput->Add((TH1D*) fhDrecoilTTCinMB[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDrecoil_HM_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDrecoilTTCinHM[igg] =  (TH1D*) fhDrecoilTTHinMB[0]->Clone(name.Data()); 
      fOutput->Add((TH1D*) fhDrecoilTTCinHM[igg]);   
   }
   
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDrecoil_GA_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDrecoilTTCinGA[igg] =   (TH1D*) fhDrecoilTTHinMB[0]->Clone(name.Data()); 
      fOutput->Add((TH1D*) fhDrecoilTTCinGA[igg]); 
   } 
   if(fMC){
      fhPtTrkTruePrimGen = new TH2D("fhPtTrkTruePrimGen","fhPtTrkTruePrimGen",100,0,100,20,-1,1);
      fOutput->Add((TH2D*) fhPtTrkTruePrimGen); 

      fhPtTrkTruePrimRec = new TH2D("fhPtTrkTruePrimRec","fhPtTrkTruePrimRec",100,0,100,20,-1,1);
      fOutput->Add((TH2D*) fhPtTrkTruePrimRec); 

      fhPtTrkSecOrFakeRec = new TH2D("fhPtTrkSecOrFakeRec","fhPtTrkSecOrFakeRec",100,0,100,20,-1,1);
      fOutput->Add((TH2D*) fhPtTrkSecOrFakeRec); 
      
      fhJetPtPartLevelCorr = new TH1D("fhJetPtPartLevelCorr","fhJetPtPartLevelCorr",270,-20,250);
      fOutput->Add((TH1D*) fhJetPtPartLevelCorr);

      fhJetPtPartLevelZero = new TH1D("fhJetPtPartLevelZero","fhJetPtPartLevelZero",250,0,250);
      fOutput->Add((TH1D*) fhJetPtPartLevelZero);

      fhFractionOfSecInJet = new TH2D("fhFractionOfSecInJet", "Frac of jet pT carried by secondary tracks",50,0,50,210,0,1.05); 
      fOutput->Add((TH2D*) fhFractionOfSecInJet);

      fhJetPtPartLevelVsJetPtDetLevelCorr = new TH2D("fhJetPtPartLevelVsJetPtDetLevelCorr","fhJetPtPartLevelVsJetPtDetLevelCorr",270,-20,250,270,-20,250);
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr); 

      fhJetPtPartLevelVsJetPtDetLevelZero = new TH2D("fhJetPtPartLevelVsJetPtDetLevelCorr","fhJetPtPartLevelVsJetPtDetLevelCorr",250,0,250,250,0,250);
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero);
 
      fhJetPtResolutionVsPtPartLevel = new TH2D("fhJetPtResolutionVsPtPartLevel","fhJetPtResolutionVsPtPartLevel",100,0,100,50,0,2);
      fOutput->Add((TH2D*) fhJetPtResolutionVsPtPartLevel);
   }


    fhOneOverPtVsPhiNeg = new TH2F("fhOneOverPtVsPhiNeg","1/pt versus track phi negative tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
    fOutput->Add((TH2F*) fhOneOverPtVsPhiNeg);
 
    fhOneOverPtVsPhiPos = new TH2F("fhOneOverPtVsPhiPos","1/pt versus track phi positive tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
    fOutput->Add((TH2F*) fhOneOverPtVsPhiPos);
 
    fhSigmaPtOverPtVsPt = new TH2F("fhSigmaPtOverPtVsPt",
                                       "track sigma(1/pt)/ 1/pt vs pt", 100, 0, 100, 250, 0, 1);
    fOutput->Add((TH2F*) fhSigmaPtOverPtVsPt); 

    Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
    Int_t nbins = sizeof(bins)/sizeof(Double_t)-1; //pT binning for DCA distribution

    fhDCAinXVsPt = new TH2F("fhDCAinXVsPt","fhDCAinXVsPt", nbins, bins, 200, -10.,10);
    fOutput->Add((TH2F*) fhDCAinXVsPt); 

    fhDCAinYVsPt = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPt");
    fOutput->Add((TH2F*) fhDCAinYVsPt);

    if(fMC){
       fhDCAinXVsPtPhysPrimary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinXVsPtPhysPrimary");
       fOutput->Add((TH2F*) fhDCAinXVsPtPhysPrimary);

       fhDCAinYVsPtPhysPrimary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPtPhysPrimary"); 
       fOutput->Add((TH2F*) fhDCAinYVsPtPhysPrimary); 
 
       fhDCAinXVsPtSecondary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinXVsPtSecondary");
       fOutput->Add((TH2F*) fhDCAinXVsPtSecondary); 

       fhDCAinYVsPtSecondary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPtSecondary");
       fOutput->Add((TH2F*) fhDCAinYVsPtSecondary); 
    }



   // OUTPUT TREE
   if(fFillTTree){
      fCentralityTree = new TTree("fCentralityTree", "Centrality vs. multiplicity tree");
      //
      fCentralityTree->Branch("trigClass",fTrigClass,"trigClass/C");
      fCentralityTree->Branch("xVertex", &fxVertex,"xVertex/D");
      fCentralityTree->Branch("yVertex", &fyVertex,"yVertex/D");
      fCentralityTree->Branch("zVertex", &fzVertex,"zVertex/D");
      fCentralityTree->Branch("vertexer3d", &fVertexer3d,"vertexer3d/O");
      fCentralityTree->Branch("nTracklets", &fNTracklets,"nTracklets/I");
      fCentralityTree->Branch("nClusters", fNClusters,"nClusters[2]/I");
      //
      fCentralityTree->Branch("isV0ATriggered", &fIsV0ATriggered,"isV0ATriggered/I");
      fCentralityTree->Branch("isV0CTriggered", &fIsV0CTriggered,"isV0CTriggered/I");
      fCentralityTree->Branch("multV0A", &fMultV0A,"multV0A/F");
      fCentralityTree->Branch("multV0C", &fMultV0C,"multV0C/F");
      fCentralityTree->Branch("ringmultV0", fRingMultV0,"ringmultV0[8]/F");

      fCentralityTree->Branch("znctower", fZNCtower, "znctower[5]/F");
      fCentralityTree->Branch("zpctower", fZPCtower, "zpctower[5]/F");
      fCentralityTree->Branch("znatower", fZNAtower, "znatower[5]/F");
      fCentralityTree->Branch("zpatower", fZPAtower, "zpatower[5]/F");
      fCentralityTree->Branch("znctowerLG", fZNCtowerLG, "znctowerLG[5]/F");
      fCentralityTree->Branch("zpctowerLG", fZPCtowerLG, "zpctowerLG[5]/F");
      fCentralityTree->Branch("znatowerLG", fZNAtowerLG, "znatowerLG[5]/F");
      fCentralityTree->Branch("zpatowerLG", fZPAtowerLG, "zpatowerLG[5]/F");
      
      //fCentralityTree->Branch("tdc", fTDCvalues, "tdc[32][4]/I");
      //fCentralityTree->Branch("tdcSum", &fTDCSum, "tdcSum/F");
      //fCentralityTree->Branch("tdcDiff", &fTDCDiff, "tdcDiff/F");
     
      if(fSystem!=AliAnalysisTaskEA::kpp){ 
         fCentralityTree->Branch("centrV0Amult", &fCentralityV0A, "centrV0Amult/F");
         fCentralityTree->Branch("centrV0Cmult", &fCentralityV0C, "centrV0Cmult/F");
         fCentralityTree->Branch("centrSPDclu1", &fCentralityCL1, "centrSPDclu1/F");
         fCentralityTree->Branch("centrZNA", &fCentralityZNA, "centrZNA/F");
         fCentralityTree->Branch("centrZNC", &fCentralityZNC, "centrZNC/F");
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name    = Form("hadronTTbin_%d_%d",fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         object  = name;
         object.Append("/I"); //Number of tracks in given bin
        
         fCentralityTree->Branch(name.Data(), &(fHadronTT[itt]), object.Data());
      }
      
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name    = Form("jetchTTbin_%d_%d",fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         object  = name;
         object.Append("/I"); //Number of jets in given bin
        
         fCentralityTree->Branch(name.Data(), &(fJetChTT[ijj]), object.Data());
      }
      
      fOutput->Add(fCentralityTree);
   } 



   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutput->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
      if(hn){
         hn->Sumw2();
      }
   }
   TH1::AddDirectory(oldStatus);


   PostData(1, fOutput);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEA::RetrieveEventObjects() {
   //
   // retrieve event objects
   //
    if(!AliAnalysisTaskEmcalJet::RetrieveEventObjects())  return kFALSE;
 
   return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEA::Run(){
   // Run analysis code here, if needed. It will be executed before FillHistograms().

   if(!fCaloCells){
      if(fCaloCellsName.IsNull()){
         fCaloCells = InputEvent()->GetEMCALCells();
      }else{
         fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
         if(!fCaloCells) AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data()));
      }
      cout<<"load calo cells"<<endl;
   }


   
   return kTRUE;
}

//________________________________________________________________________
 
