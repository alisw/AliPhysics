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
fCentralityV0M(-1),
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
fMultV0M(0.),
fMultV0Anorm(0.),
fMultV0Cnorm(0.),
fMultV0Mnorm(0.),
fMultV0AV0Cnorm(0.),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Anorm_PartLevel(0.),
fMultV0Cnorm_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
fMultV0AV0Cnorm_PartLevel(0.),
fZEM1Energy(0),
fZEM2Energy(0),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fMC(0),
fHelperClass(0), fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhTrackPhiInclMB(0x0), fhTrackEtaInclMB(0x0),  fhTrackEtaInclHM(0x0), 
fhJetPhiIncl(0x0), fhJetEtaIncl(0x0),
fhClusterPhiInclMB(0x0), fhClusterEtaInclMB(0x0),
fhClusterPhiInclGA(0x0), fhClusterEtaInclGA(0x0),
fhRhoMB(0x0),
fhRhoHM(0x0),
fhRhoMBpart(0x0),
fhV0AvsV0C(0x0),
fhV0MvsV0Mnorm(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhTrackMultMB(0x0),
fhTrackMultHM(0x0),
fhMeanTrackPtMB(0x0),
fhMeanTrackPtHM(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhFractionOfSecInJet(0x0),
fhV0ARunByRunMB(0x0),
fhV0CRunByRunMB(0x0),
fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
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
fSystem(AliAnalysisTaskEA::kpp),
fFiducialCellCut(0x0),
fnRun(1171),
fMeanV0A_PartLevel(12.6936),
fMeanV0C_PartLevel(12.5898),
fMeanV0M_PartLevel(25.2834),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0)
{
   //default constructor

   for(Int_t i=0; i<2; i++) fNClusters[i] = 0;
   for(Int_t i=0; i<8; i++) fRingMultV0[i] = 0;

   for(Int_t i=0; i<2000; i++){
     fMeanV0A[i] = 1.;
     fMeanV0C[i] = 1.;
     fMeanV0M[i] = 1.;
     fRuns[i]    = 0;
   }


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


      fHadronTT_PartLevel[i]   = 0;
      fClusterTT_PartLevel[i]   = 0;

      fhMultTTHinMB[i] = 0x0;   
      fhMultTTHinHM[i] = 0x0;   
      fhMultTTJinMB[i] = 0x0;  
      fhMultTTJinHM[i] = 0x0;  
      fhMultTTCinMB[i] = 0x0;  
      fhMultTTCinHM[i] = 0x0;  
      fhMultTTCinGA[i] = 0x0; 

      fhTTHinMB_V0M[i]      = 0x0;
      fhTTHinMB_CentV0M[i]  = 0x0;
      fhTTHinMB_V0Mnorm1[i] = 0x0;
      fhTTHinMB_V0Mnorm2[i] = 0x0;

      fhTTHinMB_V0M_PartLevel[i]     = 0x0;
      fhTTHinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTHinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhTTHinHM_V0M[i]      = 0x0;
      fhTTHinHM_CentV0M[i]  = 0x0;
      fhTTHinHM_V0Mnorm1[i] = 0x0;
      fhTTHinHM_V0Mnorm2[i] = 0x0;

      fhTTCinMB_V0M[i]      = 0x0;
      fhTTCinMB_CentV0M[i]  = 0x0;
      fhTTCinMB_V0Mnorm1[i] = 0x0;
      fhTTCinMB_V0Mnorm2[i] = 0x0;

      fhTTCinMB_V0M_PartLevel[i]     = 0x0;
      fhTTCinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTCinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhTTCinHM_V0M[i]      = 0x0;
      fhTTCinHM_CentV0M[i]  = 0x0;
      fhTTCinHM_V0Mnorm1[i] = 0x0;
      fhTTCinHM_V0Mnorm2[i] = 0x0;

      fhTTCinGA_V0M[i]      = 0x0;
      fhTTCinGA_CentV0M[i]  = 0x0;
      fhTTCinGA_V0Mnorm1[i] = 0x0;
      fhTTCinGA_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTHinMB_V0M[i]      = 0x0;
      fhRecoilJetPtTTHinMB_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTHinMB_V0M_PartLevel[i]     = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhRecoilJetPtTTHinHM_V0M[i]      = 0x0;
      fhRecoilJetPtTTHinHM_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTHinHM_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTHinHM_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinMB_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinMB_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinMB_V0M_PartLevel[i]     = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhRecoilJetPtTTCinHM_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinHM_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinHM_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinHM_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinGA_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinGA_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinGA_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinGA_V0Mnorm2[i] = 0x0;


 
      fhDeltaPtTTHinMB_RC_CentV0M[i] = 0x0;  
      fhDeltaPtTTHinHM_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinMB_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinHM_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinGA_RC_CentV0M[i] = 0x0;

      fhDeltaPtTTHinMB_RC_V0Mnorm1[i] = 0x0;  
      fhDeltaPtTTHinHM_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinHM_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinGA_RC_V0Mnorm1[i] = 0x0;

      fhDeltaPtTTHinMB_RC_V0Mnorm2[i] = 0x0;  
      fhDeltaPtTTHinHM_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinHM_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinGA_RC_V0Mnorm2[i] = 0x0;


      fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[i] = 0x0;  
      fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[i] = 0x0;  
      fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[i] = 0x0;
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
      fhCentralityHM[ic] = 0x0;
      fhSignalMB[ic] = 0x0; 
      fhSignalHM[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhCentralityTTHinMB[ic][i] = 0x0;
         fhCentralityTTHinHM[ic][i] = 0x0;
         fhCentralityTTJinMB[ic][i] = 0x0;
         fhCentralityTTJinHM[ic][i] = 0x0;
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

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignalMB_PartLevel[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTHinMB_PartLevel[ic][i] = 0x0;
         fhSignalTTCinMB_PartLevel[ic][i] = 0x0;
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

      fhRhoTTHinMBpart[i]=0x0;
      fhRhoTTCinMBpart[i]=0x0; 
   }

   fFiducialCellCut = new AliEMCALRecoUtils();
 
   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1;
      fdeltapT[i]  = 0.; 
      fdeltapT_PartLevel[i]  = 0.; 

      fIndexTTH_PartLevel[i] = -1;
      fIndexTTC_PartLevel[i] = -1;

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      fTTC_PartLevel[i].resize(0);
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
fCentralityV0M(-1),
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
fMultV0M(0.),
fMultV0Anorm(0.),
fMultV0Cnorm(0.),
fMultV0Mnorm(0.),
fMultV0AV0Cnorm(0.),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Anorm_PartLevel(0.),
fMultV0Cnorm_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
fMultV0AV0Cnorm_PartLevel(0.),
fZEM1Energy(0),
fZEM2Energy(0),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fMC(0),
fHelperClass(0), fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhTrackPhiInclMB(0x0), fhTrackEtaInclMB(0x0),  fhTrackEtaInclHM(0x0), 
fhJetPhiIncl(0x0), fhJetEtaIncl(0x0), 
fhClusterPhiInclMB(0x0), fhClusterEtaInclMB(0x0),
fhClusterPhiInclGA(0x0), fhClusterEtaInclGA(0x0),
fhRhoMB(0x0),
fhRhoHM(0x0),
fhRhoMBpart(0x0),
fhV0AvsV0C(0x0),
fhV0MvsV0Mnorm(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhTrackMultMB(0x0),
fhTrackMultHM(0x0),
fhMeanTrackPtMB(0x0),
fhMeanTrackPtHM(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhFractionOfSecInJet(0x0),
fhV0ARunByRunMB(0x0),
fhV0CRunByRunMB(0x0),
fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
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
fSystem(AliAnalysisTaskEA::kpp),
fFiducialCellCut(0x0),
fnRun(1171),
fMeanV0A_PartLevel(12.6936),
fMeanV0C_PartLevel(12.5898),
fMeanV0M_PartLevel(25.2834),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0)
{
   //Constructor

   for(Int_t i=0; i<2; i++) fNClusters[i] = 0;
   for(Int_t i=0; i<8; i++) fRingMultV0[i] = 0;

   for(Int_t i=0; i<2000; i++){ 
      fMeanV0A[i] = 1.;
      fMeanV0C[i] = 1.;
      fMeanV0M[i] = 1.;
      fRuns[i]    = 0;
   }

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

      fHadronTT_PartLevel[i] = 0;
      fClusterTT_PartLevel[i] = 0;

      fhMultTTHinMB[i] = 0x0;   
      fhMultTTHinHM[i] = 0x0;   
      fhMultTTJinMB[i] = 0x0;  
      fhMultTTJinHM[i] = 0x0;  
      fhMultTTCinMB[i] = 0x0;  
      fhMultTTCinHM[i] = 0x0;  
      fhMultTTCinGA[i] = 0x0; 

      fhTTHinMB_V0M[i]      = 0x0;
      fhTTHinMB_CentV0M[i]  = 0x0;
      fhTTHinMB_V0Mnorm1[i] = 0x0;
      fhTTHinMB_V0Mnorm2[i] = 0x0;

      fhTTHinMB_V0M_PartLevel[i]     = 0x0;
      fhTTHinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTHinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhTTHinHM_V0M[i]      = 0x0;
      fhTTHinHM_CentV0M[i]  = 0x0;
      fhTTHinHM_V0Mnorm1[i] = 0x0;
      fhTTHinHM_V0Mnorm2[i] = 0x0;

      fhTTCinMB_V0M[i]      = 0x0;
      fhTTCinMB_CentV0M[i]  = 0x0;
      fhTTCinMB_V0Mnorm1[i] = 0x0;
      fhTTCinMB_V0Mnorm2[i] = 0x0;

      fhTTCinMB_V0M_PartLevel[i]     = 0x0;
      fhTTCinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTCinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhTTCinHM_V0M[i]      = 0x0;
      fhTTCinHM_CentV0M[i]  = 0x0;
      fhTTCinHM_V0Mnorm1[i] = 0x0;
      fhTTCinHM_V0Mnorm2[i] = 0x0;

      fhTTCinGA_V0M[i]      = 0x0;
      fhTTCinGA_CentV0M[i]  = 0x0;
      fhTTCinGA_V0Mnorm1[i] = 0x0;
      fhTTCinGA_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTHinMB_V0M[i]      = 0x0;
      fhRecoilJetPtTTHinMB_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTHinMB_V0M_PartLevel[i]     = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhRecoilJetPtTTHinHM_V0M[i]      = 0x0;
      fhRecoilJetPtTTHinHM_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTHinHM_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTHinHM_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinMB_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinMB_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinMB_V0M_PartLevel[i]     = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[i] = 0x0;

      fhRecoilJetPtTTCinHM_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinHM_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinHM_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinHM_V0Mnorm2[i] = 0x0;

      fhRecoilJetPtTTCinGA_V0M[i]      = 0x0;
      fhRecoilJetPtTTCinGA_CentV0M[i]  = 0x0;
      fhRecoilJetPtTTCinGA_V0Mnorm1[i] = 0x0;
      fhRecoilJetPtTTCinGA_V0Mnorm2[i] = 0x0;

      fhDeltaPtTTHinMB_RC_CentV0M[i] = 0x0;  
      fhDeltaPtTTHinHM_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinMB_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinHM_RC_CentV0M[i] = 0x0;
      fhDeltaPtTTCinGA_RC_CentV0M[i] = 0x0;

      fhDeltaPtTTHinMB_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTHinHM_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinHM_RC_V0Mnorm1[i] = 0x0;
      fhDeltaPtTTCinGA_RC_V0Mnorm1[i] = 0x0;

      fhDeltaPtTTHinMB_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTHinHM_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinHM_RC_V0Mnorm2[i] = 0x0;
      fhDeltaPtTTCinGA_RC_V0Mnorm2[i] = 0x0;


      fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[i] = 0x0;  
      fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[i] = 0x0;  
      fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[i] = 0x0;
      fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[i] = 0x0;
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
      fhCentralityHM[ic] = 0x0;
      fhSignalMB[ic] = 0x0; 
      fhSignalHM[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhCentralityTTHinMB[ic][i] = 0x0;
         fhCentralityTTHinHM[ic][i] = 0x0;
         fhCentralityTTJinMB[ic][i] = 0x0;
         fhCentralityTTJinHM[ic][i] = 0x0;
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

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignalMB_PartLevel[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTHinMB_PartLevel[ic][i] = 0x0;
         fhSignalTTCinMB_PartLevel[ic][i] = 0x0;
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

      fhRhoTTHinMBpart[i]=0x0;
      fhRhoTTCinMBpart[i]=0x0; 
   }

   sprintf(fTrigClass,"%s","");
   //inclusive pT spectrum times the boost function

   fFiducialCellCut = new AliEMCALRecoUtils();

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1;

      fdeltapT[i]  = 0.; 
      fdeltapT_PartLevel[i]  = 0.; 
 
      fIndexTTH_PartLevel[i] = -1;
      fIndexTTC_PartLevel[i] = -1;

      fTTC[i].resize(0); 
      fTTH[i].resize(0); 
      fTTJ[i].resize(0);
 
      fTTH_PartLevel[i].resize(0); 
      fTTC_PartLevel[i].resize(0); 
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
   AliTrackContainer    *trackCont        = 0x0; //detector level track container 
   AliParticleContainer *trackContTrue    = 0x0; //mc particle container for jets
   AliClusterContainer  *clusterCont      = 0x0; //detector level track container 

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks 
   trackCont->SetMinPt(0.15);
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
   trackCont->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
   trackCont->SetAODFilterBits((1 << 8) | (1 << 9));

   if(isMC){
      trackContTrue = task->AddMCParticleContainer(mcpariclearrayname); //particle level MC particles   
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(0.15);
      trackContTrue->SetEtaLimits(-5.1,5.1);
      //trackContTrue->SetEtaLimits(-trackEtaWindow,trackEtaWindow); //V0 eta range
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
   task->SetUseNewCentralityEstimation(kTRUE);  //CENTRALITY

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
  UInt_t triggerMask = fInputHandler->IsEventSelected();
  bool passedGammaTrigger = kFALSE;

  if(triggerMask & AliVEvent::AliVEvent::kEMCEGA){
     //EG1 high EMCAL trigger, EG2 low EMCAL trigger, DG1 high DCAL trigger, DG2 low emcal trigger 
     if(trigger.Contains("EG1")){
        passedGammaTrigger = kTRUE;
     }
  }
  return passedGammaTrigger;

}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedMinBiasTrigger(){
  //minimum bias trigger
  //TString trigger = fInputEvent->GetFiredTriggerClasses();
  bool passedTrigger = kFALSE;

  if(!fMC){  // REAL DATA
     UInt_t triggerMask = fInputHandler->IsEventSelected();
     if(triggerMask & AliVEvent::kINT7){
        passedTrigger = kTRUE;
     }
  }else{   //MONTE CARLO
     passedTrigger = kTRUE;
  }
  return passedTrigger;

}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedHighMultTrigger(){
   //high multiplicity V0M trigger
   UInt_t triggerMask = fInputHandler->IsEventSelected();
   bool passedTrigger = kFALSE;

   if(!fMC){
      if(triggerMask & AliVEvent::kHighMultV0){
         passedTrigger = kTRUE;
      }
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
   TLorentzVector myTT;
   Int_t idx;
   TString name;
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
   // DECIDE WHETHER TO FILL SIGNAL TT OR REFERENCE TT  DEPENDING ON RANDOM  NUMBER  

   fFillSigTT = kTRUE;  
   if( fRandom->Integer(100) < 5) fFillSigTT = kFALSE; 

   //_________________________________________________________________
   //                EVENT PROPERTIES   

   for(int ir=0; ir<8; ir++) fRingMultV0[ir]=0.;

   // ***** Trigger selection
   TString triggerClass = InputEvent()->GetFiredTriggerClasses();
   sprintf(fTrigClass,"%s",triggerClass.Data());


   fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(fMultSelection){
 
      fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      //fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      //fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
      //fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
   }else{
      fCentralityV0A = -1; 
      fCentralityV0C = -1;
      fCentralityV0M = -1;
      fCentralityCL1 = -1;
      fCentralityZNA = -1;
      fCentralityZNC = -1;
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

   Int_t runnumber  = InputEvent()->GetRunNumber();  

   AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
   if(vzeroAOD){
      fMultV0A = vzeroAOD->GetMTotV0A();
      fMultV0C = vzeroAOD->GetMTotV0C();
      fMultV0M = fMultV0A + fMultV0C;

      Long64_t irun = 0;
      if(runnumber < fRuns[0]){
         AliError(Form("%s: RUN NOT FOUND  %d!", GetName(), runnumber)); //index will remain 0
      }else{ 
         irun = TMath::BinarySearch(fnRun, fRuns, runnumber); //index of the given run number
      
         if(fRuns[irun] != runnumber){ 

            AliError(Form("%s: RUN NOT FOUND  %d!", GetName(), runnumber));

            if(runnumber > fRuns[fnRun-1])  irun = fnRun-1;  //index will correspond to the last run in the list
         }
      }

      fMultV0Anorm = fMultV0A/fMeanV0A[irun];
      fMultV0Cnorm = fMultV0C/fMeanV0C[irun];
      fMultV0Mnorm = fMultV0M/fMeanV0M[irun];
      fMultV0AV0Cnorm = 0.5*(fMultV0Anorm + fMultV0Cnorm);

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
      fMultV0M = -1; 
      fMultV0Anorm = -1; 
      fMultV0Cnorm = -1;
      fMultV0Mnorm = -1; 
      fMultV0AV0Cnorm = -1;

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

   fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(fMyTrackContainerName.Data())); //reconstructed track container detector-level 
   fJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyJetContainerName.Data())); //AKT jet  container detector-level 


   if(fMC){
      fParticleContainerPartLevel = GetParticleContainer(fMyParticleContainerName.Data()); // particle container particle-level 
      fJetContainerPartLevel      = GetJetContainer(fMyJetParticleContainerName.Data());   //jet container  particle-level

      rhoMC = GetExternalRho(kPartLevel); //estimated backround pt density
   }
   //________________________________________________________
   //     Find the leading and subleading jets  for estimates  of Delta pt

   Double_t ptLJ=-1, etaLJ=999, phiLJ=0; //leading jet
   Double_t ptSJ=-1, etaSJ=999, phiSJ=0; //subleading jet
   //Exclude 2 leading jets 
   for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                   // trackIterator is a std::map of AliTLorentzVector and AliVTrack
       jet = jetIterator.second;  // Get the pointer to jet object
       if(!jet)  continue; 
    
       if(jet->Pt() > ptLJ){
          ptSJ  = ptLJ;
          etaSJ = etaLJ;
          phiSJ = phiLJ;

          ptLJ  = jet->Pt(); 
          etaLJ = jet->Eta(); 
          phiLJ = jet->Phi(); 
       }else if(jet->Pt() > ptSJ){
          ptSJ  = jet->Pt();
          etaSJ = jet->Eta();
          phiSJ = jet->Phi(); 
       }
   }

   //Exclude 2 leading jets MC 
   Double_t ptLJmc=-1, etaLJmc=999, phiLJmc=0; //leading jet
   Double_t ptSJmc=-1, etaSJmc=999, phiSJmc=0; //subleading jet
   if(fMC){
      for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
                      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
          jet = jetIterator.second;  // Get the pointer to jet object
          if(!jet)  continue; 
       
          if(jet->Pt() > ptLJmc){
             ptSJmc  = ptLJmc;
             etaSJmc = etaLJmc;
             phiSJmc = phiLJmc;
      
             ptLJmc  = jet->Pt(); 
             etaLJmc = jet->Eta(); 
             phiLJmc = jet->Phi(); 
          }else if(jet->Pt() > ptSJmc){
             ptSJmc  = jet->Pt();
             etaSJmc = jet->Eta();
             phiSJmc = jet->Phi(); 
          }
      }
   }
 
   //_________________________________________________________
   //                 TT

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1;
 
      fIndexTTH_PartLevel[i] = -1;
      fIndexTTC_PartLevel[i] = -1;

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      fTTC_PartLevel[i].resize(0);

      fdeltapT[i]  = 0.; 
      fdeltapT_PartLevel[i]  = 0.; 
   }

   TLorentzVector ph;
   for(Int_t i=0; i<fnClusterTTBins; i++){
      fClusterTT[i] = 0;
      fClusterTT_PartLevel[i] = 0;
   }

   for(Int_t i=0; i<fnHadronTTBins; i++){
      fHadronTT[i] = 0;
      fHadronTT_PartLevel[i] = 0;
   }


   //___________________________________________
   //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX

   fMultV0A_PartLevel = 0.;
   fMultV0C_PartLevel = 0.;
   fMultV0M_PartLevel = 0.;

   fMultV0Anorm_PartLevel = 0.;
   fMultV0Cnorm_PartLevel = 0.;
   fMultV0Mnorm_PartLevel = 0.;
   fMultV0AV0Cnorm_PartLevel = 0.;


   if(fMC){
 
      fhRhoMBpart->Fill(rhoMC);

      TClonesArray* arrayMC = 0; // array particles in the MC event
      arrayMC = (TClonesArray*) InputEvent()->FindListObject(AliAODMCParticle::StdBranchName());
      if(!arrayMC){
         AliError("No MC array found!");
         return kFALSE;
      }
      AliAODMCParticle* daughtermc; 

      Bool_t bRecPrim = kFALSE;

      if(fParticleContainerPartLevel){

         //pT spectrum of particle level physical primary mc particles
         for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
            mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!mcParticle)  continue; 

            if(IsTrackInAcceptance(mcParticle, kPartLevel)){
               fhPtTrkTruePrimGen->Fill(mcParticle->Pt(), mcParticle->Eta());

               for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                  if(fHadronTTLowPt[itt] < mcParticle->Pt() && mcParticle->Pt() < fHadronTTHighPt[itt]){
                     myTT.SetPtEtaPhiM(mcParticle->Pt(),mcParticle->Eta(),mcParticle->Phi(),0.); 
                     fTTH_PartLevel[itt].push_back(myTT);
                     fHadronTT_PartLevel[itt]++;   // there was a high pt 
                  }
               }
            } 
          

            if(mcParticle->Charge()){ 
               if((static_cast<AliAODMCParticle*>(mcParticle))->IsPhysicalPrimary()){
                  //get particle level charged particles multiplicities in V0A and V0C
                  if(-3.7 < mcParticle->Eta() && mcParticle->Eta() < -1.7) fMultV0C_PartLevel++;
                  if( 2.8 < mcParticle->Eta() && mcParticle->Eta() < 5.1)  fMultV0A_PartLevel++; 
               } 
            }else{
               //TT Cluster
                if(((static_cast<AliAODMCParticle*>(mcParticle))->IsPhysicalPrimary()) && TMath::Abs(mcParticle->Eta())<0.7){ //EMCAL acceptance
                   if(((static_cast<AliAODMCParticle*>(mcParticle))->GetPdgCode())==22){ //photon
                       //skip photons which have a photon as a daughter particle
                       Int_t d1 = TMath::Abs((static_cast<AliAODMCParticle*> (mcParticle))->GetDaughterLabel(0));
                       Int_t d2 = TMath::Abs((static_cast<AliAODMCParticle*> (mcParticle))->GetDaughterLabel(1)); 
                       Bool_t hasPhotonicDaughter = 0;
                       for(Int_t id=d1;id<=d2; id++){
                          daughtermc = (AliAODMCParticle*) arrayMC->At(id);
                          if(daughtermc->GetPdgCode()==22){ hasPhotonicDaughter=1;  break;}
                       }
                       if(!hasPhotonicDaughter){
                   
                          for(Int_t igg=0; igg<fnClusterTTBins; igg++){
                             if(fClusterTTLowPt[igg] < mcParticle->Pt() && mcParticle->Pt() < fClusterTTHighPt[igg]){
                                         
                                myTT.SetPtEtaPhiM(mcParticle->Pt(),mcParticle->Eta(),mcParticle->Phi(),0.); 
                                fTTC_PartLevel[igg].push_back(myTT);
                                fClusterTT_PartLevel[igg]++;   // there was a high pt emcal cluster 
                             }
                          }
                       }
                    }
                } 
             } 
         }//end of mc particle loop

         //combined V0 multiplicities particle level
         fMultV0M_PartLevel = fMultV0A_PartLevel + fMultV0C_PartLevel;
         fMultV0Anorm_PartLevel = fMultV0A_PartLevel/fMeanV0A_PartLevel;
         fMultV0Cnorm_PartLevel = fMultV0C_PartLevel/fMeanV0C_PartLevel;
         fMultV0Mnorm_PartLevel = fMultV0M_PartLevel/fMeanV0M_PartLevel;
         fMultV0AV0Cnorm_PartLevel = 0.5*(fMultV0Anorm_PartLevel + fMultV0Cnorm_PartLevel);


         //minimum bias particle level

         fhSignalMB_PartLevel[fkV0A]->Fill(fMultV0A_PartLevel);
         fhSignalMB_PartLevel[fkV0C]->Fill(fMultV0C_PartLevel);
         fhSignalMB_PartLevel[fkV0M]->Fill(fMultV0M_PartLevel);
         fhSignalMB_PartLevel[fkV0Mnorm1]->Fill(fMultV0Mnorm_PartLevel);
         fhSignalMB_PartLevel[fkV0Mnorm2]->Fill(fMultV0AV0Cnorm_PartLevel);


         //chose trigger hadron TT   particle level
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(fHadronTT_PartLevel[itt]>0){
               fIndexTTH_PartLevel[itt] = fRandom->Integer(fHadronTT_PartLevel[itt]); 
               idx = fIndexTTH_PartLevel[itt]; 
               
               fdeltapT_PartLevel[itt] = GetDeltaPt(fTTH_PartLevel[itt][idx].Phi(), fTTH_PartLevel[itt][idx].Eta(), phiLJmc, etaLJmc, phiSJmc, etaSJmc, rhoMC, kPartLevel);

               //signal in events with hadron TT   particle level
               fhSignalTTHinMB_PartLevel[fkV0A][itt]->Fill(fMultV0A_PartLevel);
               fhSignalTTHinMB_PartLevel[fkV0C][itt]->Fill(fMultV0C_PartLevel);
               fhSignalTTHinMB_PartLevel[fkV0M][itt]->Fill(fMultV0M_PartLevel);
               fhSignalTTHinMB_PartLevel[fkV0Mnorm1][itt]->Fill(fMultV0Mnorm_PartLevel);
               fhSignalTTHinMB_PartLevel[fkV0Mnorm2][itt]->Fill(fMultV0AV0Cnorm_PartLevel);

               //hadron trigger particle level
              
               if(idx>-1){ 

                  fhRhoTTHinMBpart[itt]->Fill(rhoMC);
                  fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fdeltapT_PartLevel[itt]);
                  fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[itt]->Fill(fMultV0AV0Cnorm_PartLevel, fdeltapT_PartLevel[itt]);


                  if(fFillSigTT && itt==0) continue;  // Do not fill reference 
                  if(!fFillSigTT && itt>0) continue;  // Do not fill signal 

                  fhTTHinMB_V0M_PartLevel[itt]->Fill(fMultV0M_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0M
                  fhTTHinMB_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
                  fhTTHinMB_V0Mnorm2_PartLevel[itt]->Fill(fMultV0AV0Cnorm_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
                
                  //recoil jets  PARTICLE LEVEL
                  for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
                     // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                     jet = jetIterator.second;  // Get the pointer to jet object
                     if(!jet)  continue; 
              
                     if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH_PartLevel[itt][idx].Phi())) > fPhiCut){     
                        //recoil jet hadron trigger
                        jetPtcorr = jet->Pt() - rhoMC*jet->Area();
                        fhRecoilJetPtTTHinMB_V0M_PartLevel[itt]->Fill(fMultV0M_PartLevel, jetPtcorr);
                        fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jetPtcorr);
                        fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[itt]->Fill(fMultV0AV0Cnorm_PartLevel, jetPtcorr);
                     }
                  } 
               }
            }
         }

         //chose trigger emcal cluster TT 
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            if(fClusterTT_PartLevel[igg]>0){ 
               fIndexTTC_PartLevel[igg] = fRandom->Integer(fClusterTT_PartLevel[igg]);
               idx = fIndexTTC_PartLevel[igg];// gamma trigger

               fdeltapT_PartLevel[igg] = GetDeltaPt(fTTC_PartLevel[igg][idx].Phi(), fTTC_PartLevel[igg][idx].Eta(), phiLJmc, etaLJmc, phiSJmc, etaSJmc, rhoMC, kPartLevel);

               //signal in events with hadron TT   particle level
               fhSignalTTCinMB_PartLevel[fkV0A][igg]->Fill(fMultV0A_PartLevel);
               fhSignalTTCinMB_PartLevel[fkV0C][igg]->Fill(fMultV0C_PartLevel);
               fhSignalTTCinMB_PartLevel[fkV0M][igg]->Fill(fMultV0M_PartLevel);
               fhSignalTTCinMB_PartLevel[fkV0Mnorm1][igg]->Fill(fMultV0Mnorm_PartLevel);
               fhSignalTTCinMB_PartLevel[fkV0Mnorm2][igg]->Fill(fMultV0AV0Cnorm_PartLevel);

               if(idx>-1){ 

                  fhRhoTTCinMBpart[igg]->Fill(rhoMC);
                  fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[igg]->Fill(fMultV0Mnorm_PartLevel, fdeltapT_PartLevel[igg]);
                  fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[igg]->Fill(fMultV0AV0Cnorm_PartLevel, fdeltapT_PartLevel[igg]);

                  if(fFillSigTT && igg==0) continue;  // Do not fill reference 
                  if(!fFillSigTT && igg>0) continue;  // Do not fill signal 

                  fhTTCinMB_V0M_PartLevel[igg]->Fill(fMultV0M_PartLevel, fTTC_PartLevel[igg][idx].Pt()); //fill trigger track pT
                  fhTTCinMB_V0Mnorm1_PartLevel[igg]->Fill(fMultV0Mnorm_PartLevel, fTTC_PartLevel[igg][idx].Pt()); //fill trigger track pT
                  fhTTCinMB_V0Mnorm2_PartLevel[igg]->Fill(fMultV0AV0Cnorm_PartLevel, fTTC_PartLevel[igg][idx].Pt()); //fill trigger track pT
 
                  //recoil jets PARTICLE LEVEL
                  for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
                     // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                     jet = jetIterator.second;  // Get the pointer to jet object
                     if(!jet)  continue; 
                 
                     if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC_PartLevel[igg][idx].Phi())) > fPhiCut){     
                        //recoil jet
                        jetPtcorr = jet->Pt() - rhoMC*jet->Area();
                        fhRecoilJetPtTTCinMB_V0M_PartLevel[igg]->Fill(fMultV0M_PartLevel, jetPtcorr);
                        fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[igg]->Fill(fMultV0Mnorm_PartLevel, jetPtcorr);
                        fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[igg]->Fill(fMultV0AV0Cnorm_PartLevel, jetPtcorr);
                     }
                  } 
               }
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
            sumFakeTrackPtInJet = 0.;

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
               fhFractionOfSecInJet->Fill(jetPtCorrDet, sumFakeTrackPtInJet/sumAllTrackPtInJet); 
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
   }
   //___________________________________________
   //    INCLUSIVE EVENTS (WITHOUT TT REQUIREMENT)

   Int_t    trackMult  = 0;
   Double_t sumTrackPt = 0.;

   for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      track = trackIterator.second;  // Get the full track
      if(!track) continue;

      if(IsTrackInAcceptance(track, kDetLevel)){  
         trackMult++;
         sumTrackPt += track->Pt();
      }
   }
   if(trackMult>0){
      sumTrackPt = sumTrackPt/trackMult;
   }

   if(fIsMinBiasTrig){ // run for all MC events   and for real data with min bias trigger
       fhVertex[0]->Fill(fxVertex);
       fhVertex[1]->Fill(fyVertex);
       fhVertex[2]->Fill(fzVertex);

       fhRhoMB->Fill(rho);

       fhCentralityMB[fkV0A]->Fill(fCentralityV0A, fMultV0A); 
       fhCentralityMB[fkV0C]->Fill(fCentralityV0C, fMultV0C); 
       fhCentralityMB[fkV0M]->Fill(fCentralityV0M, fMultV0M); 
       fhCentralityMB[fkV0Mnorm1]->Fill(fCentralityV0M, fMultV0Mnorm); 
       fhCentralityMB[fkV0Mnorm2]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
       //fhCentralityMB[fkSPD]->Fill(fCentralityCL1);
       //fhCentralityMB[fkZNA]->Fill(fCentralityZNA);
       //fhCentralityMB[fkZNC]->Fill(fCentralityZNC);

       fhSignalMB[fkV0A]->Fill(fMultV0A);
       fhSignalMB[fkV0C]->Fill(fMultV0C);
       //fhSignalMB[fkSPD]->Fill(fNTracklets); 
       //fhSignalMB[fkZNA]->Fill(fZNAtower[0]); 
       //fhSignalMB[fkZNC]->Fill(fZNCtower[0]);
       fhSignalMB[fkV0M]->Fill(fMultV0M);
       fhSignalMB[fkV0Mnorm1]->Fill(fMultV0Mnorm);
       fhSignalMB[fkV0Mnorm2]->Fill(fMultV0AV0Cnorm);

       name = Form("%d", runnumber);      
       fhV0ARunByRunMB->Fill(name.Data(), fMultV0A, 1.0);
       fhV0CRunByRunMB->Fill(name.Data(), fMultV0C, 1.0);
       fhV0MRunByRunMB->Fill(name.Data(), fMultV0M, 1.0);
       fhV0MnormRunByRunMB->Fill(name.Data(), fMultV0Mnorm, 1.0);

       fhV0AvsV0C->Fill(fMultV0C, fMultV0A);
       fhV0MvsV0Mnorm->Fill(fMultV0AV0Cnorm, fMultV0M);
       fhV0AvsSPD->Fill(fNTracklets, fMultV0A);
       fhV0CvsSPD->Fill(fNTracklets, fMultV0C);

       fhTrackMultMB->Fill(fCentralityV0M, trackMult); 
       fhMeanTrackPtMB->Fill(fCentralityV0M, sumTrackPt);
    }

    
    if(fIsHighMultTrig){ 
       fhRhoHM->Fill(rho);

       fhCentralityHM[fkV0A]->Fill(fCentralityV0A, fMultV0A); 
       fhCentralityHM[fkV0C]->Fill(fCentralityV0C, fMultV0C); 
       fhCentralityHM[fkV0M]->Fill(fCentralityV0M, fMultV0M); 
       fhCentralityHM[fkV0Mnorm1]->Fill(fCentralityV0M, fMultV0Mnorm); 
       fhCentralityHM[fkV0Mnorm2]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
 
       fhSignalHM[fkV0A]->Fill(fMultV0A);
       fhSignalHM[fkV0C]->Fill(fMultV0C);
       //fhSignalHM[fkSPD]->Fill(fNTracklets); 
       //fhSignalHM[fkZNA]->Fill(fZNAtower[0]); 
       //fhSignalHM[fkZNC]->Fill(fZNCtower[0]);
       fhSignalHM[fkV0M]->Fill(fMultV0M);
       fhSignalHM[fkV0Mnorm1]->Fill(fMultV0Mnorm);
       fhSignalHM[fkV0Mnorm2]->Fill(fMultV0AV0Cnorm);

        
       fhTrackMultHM->Fill(fCentralityV0M, trackMult); 
       fhMeanTrackPtHM->Fill(fCentralityV0M, sumTrackPt);
   }

   //_________________________________________________________
   //LOOP OVER EMCAL CLUSTERS
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
           idx = fIndexTTC[igg];// gamma trigger
           fdeltapT[igg] = GetDeltaPt(fTTC[igg][idx].Phi(), fTTC[igg][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, rho, kDetLevel);
        }
     }


 
      if(fIsMinBiasTrig){ 
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         
            fhMultTTCinMB[igg]->Fill(fClusterTT[igg]); 
            
            if(!fClusterTT[igg]) continue;

            fhRhoTTCinMB[igg]->Fill(rho);

           
            fhCentralityTTCinMB[fkV0A][igg]->Fill(fCentralityV0A, fMultV0A); 
            fhCentralityTTCinMB[fkV0C][igg]->Fill(fCentralityV0C, fMultV0C); 
            fhCentralityTTCinMB[fkV0M][igg]->Fill(fCentralityV0M, fMultV0M); 
            fhCentralityTTCinMB[fkV0Mnorm1][igg]->Fill(fCentralityV0M, fMultV0Mnorm); 
            fhCentralityTTCinMB[fkV0Mnorm2][igg]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
            //fhCentralityTTCinMB[fkSPD][igg]->Fill(fCentralityCL1);
            //fhCentralityTTCinMB[fkZNA][igg]->Fill(fCentralityZNA);
            //fhCentralityTTCinMB[fkZNC][igg]->Fill(fCentralityZNC);
 
            fhSignalTTCinMB[fkV0A][igg]->Fill(fMultV0A);
            fhSignalTTCinMB[fkV0C][igg]->Fill(fMultV0C);
            //fhSignalTTCinMB[fkSPD][igg]->Fill(fNTracklets); 
            //fhSignalTTCinMB[fkZNA][igg]->Fill(fZNAtower[0]); 
            //fhSignalTTCinMB[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinMB[fkV0M][igg]->Fill(fMultV0M);
            fhSignalTTCinMB[fkV0Mnorm1][igg]->Fill(fMultV0Mnorm);
            fhSignalTTCinMB[fkV0Mnorm2][igg]->Fill(fMultV0AV0Cnorm);


            fhV0AvsV0CTTCinMB[igg]->Fill(fMultV0C, fMultV0A);


            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){

               fhDeltaPtTTCinMB_RC_CentV0M[igg]->Fill(fCentralityV0M, fdeltapT[igg]);
               fhDeltaPtTTCinMB_RC_V0Mnorm1[igg]->Fill(fMultV0Mnorm, fdeltapT[igg]);
               fhDeltaPtTTCinMB_RC_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, fdeltapT[igg]);

               if(fFillSigTT && igg==0) continue;  // Do not fill reference 
               if(!fFillSigTT && igg>0) continue;  // Do not fill signal 

               fhTTCinMB_V0M[igg]->Fill(fMultV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinMB_V0Mnorm1[igg]->Fill(fMultV0Mnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinMB_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinMB_CentV0M[igg]->Fill(fCentralityV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhRecoilJetPtTTCinMB_V0M[igg]->Fill(fMultV0M, jetPtcorr);
                     fhRecoilJetPtTTCinMB_V0Mnorm1[igg]->Fill(fMultV0Mnorm, jetPtcorr);
                     fhRecoilJetPtTTCinMB_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, jetPtcorr);
                     fhRecoilJetPtTTCinMB_CentV0M[igg]->Fill(fCentralityV0M, jetPtcorr);
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
            //fhSignalTTCinHM[fkSPD][igg]->Fill(fNTracklets); 
            //fhSignalTTCinHM[fkZNA][igg]->Fill(fZNAtower[0]); 
            //fhSignalTTCinHM[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinHM[fkV0M][igg]->Fill(fMultV0M);
            fhSignalTTCinHM[fkV0Mnorm1][igg]->Fill(fMultV0Mnorm);
            fhSignalTTCinHM[fkV0Mnorm2][igg]->Fill(fMultV0AV0Cnorm);

            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){ 

               fhDeltaPtTTCinHM_RC_CentV0M[igg]->Fill(fCentralityV0M,  fdeltapT[igg]);
               fhDeltaPtTTCinHM_RC_V0Mnorm1[igg]->Fill(fMultV0Mnorm,  fdeltapT[igg]);
               fhDeltaPtTTCinHM_RC_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm,  fdeltapT[igg]);

               if(fFillSigTT && igg==0) continue;  // Do not fill reference 
               if(!fFillSigTT && igg>0) continue;  // Do not fill signal 

               fhTTCinHM_V0M[igg]->Fill(fMultV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinHM_V0Mnorm1[igg]->Fill(fMultV0Mnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinHM_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinHM_CentV0M[igg]->Fill(fCentralityV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhRecoilJetPtTTCinHM_V0M[igg]->Fill(fMultV0M, jetPtcorr);
                     fhRecoilJetPtTTCinHM_V0Mnorm1[igg]->Fill(fMultV0Mnorm, jetPtcorr);
                     fhRecoilJetPtTTCinHM_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, jetPtcorr);
                     fhRecoilJetPtTTCinHM_CentV0M[igg]->Fill(fCentralityV0M, jetPtcorr);
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
            
            fhCentralityTTCinGA[fkV0A][igg]->Fill(fCentralityV0A, fMultV0A); 
            fhCentralityTTCinGA[fkV0C][igg]->Fill(fCentralityV0C, fMultV0C); 
            fhCentralityTTCinGA[fkV0M][igg]->Fill(fCentralityV0M, fMultV0M); 
            fhCentralityTTCinGA[fkV0Mnorm1][igg]->Fill(fCentralityV0M, fMultV0Mnorm); 
            fhCentralityTTCinGA[fkV0Mnorm2][igg]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
            //fhCentralityTTCinGA[fkSPD][igg]->Fill(fCentralityCL1);
            //fhCentralityTTCinGA[fkZNA][igg]->Fill(fCentralityZNA);
            //fhCentralityTTCinGA[fkZNC][igg]->Fill(fCentralityZNC);

            fhSignalTTCinGA[fkV0A][igg]->Fill(fMultV0A);
            fhSignalTTCinGA[fkV0C][igg]->Fill(fMultV0C);
            //fhSignalTTCinGA[fkSPD][igg]->Fill(fNTracklets); 
            //fhSignalTTCinGA[fkZNA][igg]->Fill(fZNAtower[0]); 
            //fhSignalTTCinGA[fkZNC][igg]->Fill(fZNCtower[0]);
            fhSignalTTCinGA[fkV0M][igg]->Fill(fMultV0M);
            fhSignalTTCinGA[fkV0Mnorm1][igg]->Fill(fMultV0Mnorm);
            fhSignalTTCinGA[fkV0Mnorm2][igg]->Fill(fMultV0AV0Cnorm);
            
            fhV0AvsV0CTTCinGA[igg]->Fill(fMultV0C, fMultV0A);
 
            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){

               fhDeltaPtTTCinGA_RC_CentV0M[igg]->Fill(fCentralityV0M, fdeltapT[igg]);
               fhDeltaPtTTCinGA_RC_V0Mnorm1[igg]->Fill(fMultV0Mnorm, fdeltapT[igg]);
               fhDeltaPtTTCinGA_RC_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, fdeltapT[igg]);

               if(fFillSigTT && igg==0) continue;  // Do not fill reference 
               if(!fFillSigTT && igg>0) continue;  // Do not fill signal 

               fhTTCinGA_V0M[igg]->Fill(fMultV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinGA_V0Mnorm1[igg]->Fill(fMultV0Mnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinGA_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, fTTC[igg][idx].Pt()); //fill trigger track pT
               fhTTCinGA_CentV0M[igg]->Fill(fCentralityV0M, fTTC[igg][idx].Pt()); //fill trigger track pT
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
       
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     jetPtcorr = jet->Pt() - rho*jet->Area();
                     fhRecoilJetPtTTCinGA_V0M[igg]->Fill(fMultV0M, jetPtcorr);
                     fhRecoilJetPtTTCinGA_V0Mnorm1[igg]->Fill(fMultV0Mnorm, jetPtcorr);
                     fhRecoilJetPtTTCinGA_V0Mnorm2[igg]->Fill(fMultV0AV0Cnorm, jetPtcorr);
                     fhRecoilJetPtTTCinGA_CentV0M[igg]->Fill(fCentralityV0M, jetPtcorr);
                  }
               } 
            }//trigger exists
         }//loop over TT bins
      }//emcale trigger 
   }//cluster container   

 



   //_________________________________________________________
   //LOOP OVER TRACKS DETECTOR LEVEL + SEARCH FOR HIGH PT HADRON TRIGGER 

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
            fhTrackPhiInclMB->Fill(track->Pt(), track->Phi());
            fhTrackEtaInclMB->Fill(track->Pt(), track->Eta());

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
      fdeltapT[itt] = 0.;
 
      if(fHadronTT[itt]>0){
         fIndexTTH[itt] = fRandom->Integer(fHadronTT[itt]);
         idx = fIndexTTH[itt];
         fdeltapT[itt] = GetDeltaPt(fTTH[itt][idx].Phi(), fTTH[itt][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, rho, kDetLevel);
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

         fhCentralityTTHinMB[fkV0A][itt]->Fill(fCentralityV0A, fMultV0A); 
         fhCentralityTTHinMB[fkV0C][itt]->Fill(fCentralityV0C, fMultV0C);
         fhCentralityTTHinMB[fkV0M][itt]->Fill(fCentralityV0M, fMultV0M); 
         fhCentralityTTHinMB[fkV0Mnorm1][itt]->Fill(fCentralityV0M, fMultV0Mnorm); 
         fhCentralityTTHinMB[fkV0Mnorm2][itt]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
         //fhCentralityTTHinMB[fkSPD][itt]->Fill(fCentralityCL1);
         //fhCentralityTTHinMB[fkZNA][itt]->Fill(fCentralityZNA);
         //fhCentralityTTHinMB[fkZNC][itt]->Fill(fCentralityZNC);
 
         fhSignalTTHinMB[fkV0A][itt]->Fill(fMultV0A);
         fhSignalTTHinMB[fkV0C][itt]->Fill(fMultV0C);
         //fhSignalTTHinMB[fkSPD][itt]->Fill(fNTracklets); 
         //fhSignalTTHinMB[fkZNA][itt]->Fill(fZNAtower[0]); 
         //fhSignalTTHinMB[fkZNC][itt]->Fill(fZNCtower[0]); 
         fhSignalTTHinMB[fkV0M][itt]->Fill(fMultV0M);
         fhSignalTTHinMB[fkV0Mnorm1][itt]->Fill(fMultV0Mnorm);
         fhSignalTTHinMB[fkV0Mnorm2][itt]->Fill(fMultV0AV0Cnorm);

         fhV0AvsV0CTTH[itt]->Fill(fMultV0C, fMultV0A);

         //hadron trigger
         idx = fIndexTTH[itt];
         if(idx>-1){

            fhDeltaPtTTHinMB_RC_CentV0M[itt]->Fill(fCentralityV0M, fdeltapT[itt]);
            fhDeltaPtTTHinMB_RC_V0Mnorm1[itt]->Fill(fMultV0Mnorm, fdeltapT[itt]);         
            fhDeltaPtTTHinMB_RC_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, fdeltapT[itt]);         
 
            if(fFillSigTT && itt==0) continue;  // Do not fill reference 
            if(!fFillSigTT && itt>0) continue;  // Do not fill signal 

            fhTTHinMB_V0M[itt]->Fill(fMultV0M, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0M
            fhTTHinMB_V0Mnorm1[itt]->Fill(fMultV0Mnorm, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            fhTTHinMB_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            fhTTHinMB_CentV0M[itt]->Fill(fCentralityV0M, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0M centrality

            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue; 

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){     
                  //recoil jet
                  jetPtcorr = jet->Pt() - rho*jet->Area();
                  fhRecoilJetPtTTHinMB_V0M[itt]->Fill(fMultV0M, jetPtcorr);
                  fhRecoilJetPtTTHinMB_V0Mnorm1[itt]->Fill(fMultV0Mnorm, jetPtcorr);
                  fhRecoilJetPtTTHinMB_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, jetPtcorr);
                  fhRecoilJetPtTTHinMB_CentV0M[itt]->Fill(fCentralityV0M, jetPtcorr);
               }
            } 
         }
      }
   }

   if(fIsHighMultTrig){
      for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         track = trackIterator.second;  // Get the full track
         if(!track) continue;
   
         if(IsTrackInAcceptance(track, kDetLevel)){  
            //fill some histograms for detector level tracks 
            fhTrackEtaInclHM->Fill(track->Pt(), track->Eta());
         } 
      }
 
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      
         fhMultTTHinHM[itt]->Fill(fHadronTT[itt]); 
      
         if(!fHadronTT[itt]) continue;
      
         fhRhoTTHinHM[itt]->Fill(rho); 
 
         fhCentralityTTHinHM[fkV0A][itt]->Fill(fCentralityV0A, fMultV0A); 
         fhCentralityTTHinHM[fkV0C][itt]->Fill(fCentralityV0C, fMultV0C);
         fhCentralityTTHinHM[fkV0M][itt]->Fill(fCentralityV0M, fMultV0M); 
         fhCentralityTTHinHM[fkV0Mnorm1][itt]->Fill(fCentralityV0M, fMultV0Mnorm); 
         fhCentralityTTHinHM[fkV0Mnorm2][itt]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
 
         fhSignalTTHinHM[fkV0A][itt]->Fill(fMultV0A);
         fhSignalTTHinHM[fkV0C][itt]->Fill(fMultV0C);
         //fhSignalTTHinHM[fkSPD][itt]->Fill(fNTracklets); 
         //fhSignalTTHinHM[fkZNA][itt]->Fill(fZNAtower[0]); 
         //fhSignalTTHinHM[fkZNC][itt]->Fill(fZNCtower[0]); 
         fhSignalTTHinHM[fkV0M][itt]->Fill(fMultV0M);
         fhSignalTTHinHM[fkV0Mnorm1][itt]->Fill(fMultV0Mnorm);
         fhSignalTTHinHM[fkV0Mnorm2][itt]->Fill(fMultV0AV0Cnorm);

         //hadron trigger
         idx = fIndexTTH[itt];
         if(idx>-1){ 

            fhDeltaPtTTHinHM_RC_CentV0M[itt]->Fill(fCentralityV0M, fdeltapT[itt]);
            fhDeltaPtTTHinHM_RC_V0Mnorm1[itt]->Fill(fMultV0Mnorm, fdeltapT[itt]);
            fhDeltaPtTTHinHM_RC_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, fdeltapT[itt]);

            if(fFillSigTT && itt==0) continue;  // Do not fill reference 
            if(!fFillSigTT && itt>0) continue;  // Do not fill signal 

            fhTTHinHM_V0M[itt]->Fill(fMultV0M, fTTH[itt][idx].Pt()); //fill trigger track pT
            fhTTHinHM_V0Mnorm1[itt]->Fill(fMultV0Mnorm, fTTH[itt][idx].Pt()); //fill trigger track pT
            fhTTHinHM_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, fTTH[itt][idx].Pt()); //fill trigger track pT
            fhTTHinHM_CentV0M[itt]->Fill(fCentralityV0M, fTTH[itt][idx].Pt()); //fill trigger track pT
          

            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue; 

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){     
                  //recoil jet
                  jetPtcorr = jet->Pt() - rho*jet->Area();
                  fhRecoilJetPtTTHinHM_V0M[itt]->Fill(fMultV0M, jetPtcorr);
                  fhRecoilJetPtTTHinHM_V0Mnorm1[itt]->Fill(fMultV0Mnorm, jetPtcorr);
                  fhRecoilJetPtTTHinHM_V0Mnorm2[itt]->Fill(fMultV0AV0Cnorm, jetPtcorr);
                  fhRecoilJetPtTTHinHM_CentV0M[itt]->Fill(fCentralityV0M, jetPtcorr);
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

         fhCentralityTTJinMB[fkV0A][ijj]->Fill(fCentralityV0A, fMultV0A); 
         fhCentralityTTJinMB[fkV0C][ijj]->Fill(fCentralityV0C, fMultV0C); 
         fhCentralityTTJinMB[fkV0M][ijj]->Fill(fCentralityV0M, fMultV0M); 
         fhCentralityTTJinMB[fkV0Mnorm1][ijj]->Fill(fCentralityV0M, fMultV0Mnorm); 
         fhCentralityTTJinMB[fkV0Mnorm2][ijj]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 
         //fhCentralityTTJ[fkSPD][ijj]->Fill(fCentralityCL1);
         //fhCentralityTTJ[fkZNA][ijj]->Fill(fCentralityZNA);
         //fhCentralityTTJ[fkZNC][ijj]->Fill(fCentralityZNC);
      
         fhSignalTTJinMB[fkV0A][ijj]->Fill(fMultV0A);
         fhSignalTTJinMB[fkV0C][ijj]->Fill(fMultV0C);
         //fhSignalTTJinMB[fkSPD][ijj]->Fill(fNTracklets); 
         //fhSignalTTJinMB[fkZNA][ijj]->Fill(fZNAtower[0]); 
         //fhSignalTTJinMB[fkZNC][ijj]->Fill(fZNCtower[0]); 
         fhSignalTTJinMB[fkV0M][ijj]->Fill(fMultV0M);
         fhSignalTTJinMB[fkV0Mnorm1][ijj]->Fill(fMultV0Mnorm);
         fhSignalTTJinMB[fkV0Mnorm2][ijj]->Fill(fMultV0AV0Cnorm);
      
         fhV0AvsV0CTTJ[ijj]->Fill(fMultV0C, fMultV0A);

      }
   }

   if(fIsHighMultTrig){ 
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      
         fhMultTTJinHM[ijj]->Fill(fJetChTT[ijj]); 
      
         if(!fJetChTT[ijj]) continue; 
       
         fhRhoTTJinHM[ijj]->Fill(rho);
      
         fhCentralityTTJinHM[fkV0A][ijj]->Fill(fCentralityV0A, fMultV0A); 
         fhCentralityTTJinHM[fkV0C][ijj]->Fill(fCentralityV0C, fMultV0C); 
         fhCentralityTTJinHM[fkV0M][ijj]->Fill(fCentralityV0M, fMultV0M); 
         fhCentralityTTJinHM[fkV0Mnorm1][ijj]->Fill(fCentralityV0M, fMultV0Mnorm); 
         fhCentralityTTJinHM[fkV0Mnorm2][ijj]->Fill(fCentralityV0M, fMultV0AV0Cnorm); 

         fhSignalTTJinHM[fkV0A][ijj]->Fill(fMultV0A);
         fhSignalTTJinHM[fkV0C][ijj]->Fill(fMultV0C);
         //fhSignalTTJinHM[fkSPD][ijj]->Fill(fNTracklets); 
         //fhSignalTTJinHM[fkZNA][ijj]->Fill(fZNAtower[0]); 
         //fhSignalTTJinHM[fkZNC][ijj]->Fill(fZNCtower[0]); 
         fhSignalTTJinHM[fkV0M][ijj]->Fill(fMultV0M);
         fhSignalTTJinHM[fkV0Mnorm1][ijj]->Fill(fMultV0Mnorm);
         fhSignalTTJinHM[fkV0Mnorm2][ijj]->Fill(fMultV0AV0Cnorm);
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
   delete fFiducialCellCut; 
 
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
   fHistEvtSelection = new TH1D("fHistEvtSelection", "event selection", 7, 0, 7);
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
   fhVertexZall =  new TH1D("fhVertexZall","z vertex without cut",40,-20,20);
   fOutput->Add(fhVertexZall); 
 
   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);
   fOutput->Add(fhVertexZ);
 
   //-------------------------

   fhTrackEtaInclMB = new TH2D("fhTrackEtaInclMB","Eta dist inclusive track vs pT MB", 50,0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2D*) fhTrackEtaInclMB);

   fhTrackPhiInclMB = new TH2D("fhTrackPhiInclMB","Azim dist tracks vs pT MB", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2D*) fhTrackPhiInclMB);
   
   fhTrackEtaInclHM = new TH2D("fhTrackEtaInclHM","Eta dist inclusive track vs pT HM", 50,0, 100, 40,-0.9,0.9);
   if(!fMC) fOutput->Add((TH2D*) fhTrackEtaInclHM);

   fhJetEtaIncl = new TH2D("fhJetEtaIncl","Eta dist inclusive jets vs pTjet", 150, -20, 130, 40,-0.9,0.9);
   fOutput->Add((TH2D*) fhJetEtaIncl);

   fhJetPhiIncl = new TH2D("fhJetPhiIncl","Azim dist jets vs pTjet", 60, -20, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2D*) fhJetPhiIncl);

   fhClusterEtaInclMB = new TH2D("fhClusterEtaInclMB","Eta dist inclusive clusters vs pT", 100, 0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2D*) fhClusterEtaInclMB);

   fhClusterPhiInclMB = new TH2D("fhClusterPhiInclMB","Azim dist clusters vs pT", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2D*) fhClusterPhiInclMB);

   fhClusterEtaInclGA = new TH2D("fhClusterEtaInclGA","Eta dist inclusive clusters vs pT", 100, 0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2D*) fhClusterEtaInclGA);

   fhClusterPhiInclGA = new TH2D("fhClusterPhiInclGA","Azim dist clusters vs pT", 50, 0, 100, 50,0,2*TMath::Pi());
   fOutput->Add((TH2D*) fhClusterPhiInclGA);


   //RHO 
   fhRhoMB = new TH1D("hRho_MB","Rho minimum bias det level",1000,0,100);
   fOutput->Add((TH1D*) fhRhoMB); 
 
   name = Form("hRho_HM");
   fhRhoHM = (TH1D*)  fhRhoMB->Clone(name.Data()); 
   fhRhoHM->SetTitle("Rho high multiplicity det level"); 
   if(!fMC) fOutput->Add((TH1D*) fhRhoHM); 


   for(Int_t itt=0; itt<fnHadronTTBins;itt++){
      name = Form("hRho_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRhoTTHinMB[itt] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1D*) fhRhoTTHinMB[itt]); 
   }
   
   for(Int_t itt=0; itt<fnHadronTTBins;itt++){
      name = Form("hRho_HM_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRhoTTHinHM[itt] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events HM with hadron TT
      if(!fMC) fOutput->Add((TH1D*) fhRhoTTHinHM[itt]); 
   }
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hRho_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhRhoTTJinMB[ijj] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1D*) fhRhoTTJinMB[ijj]); 
   } 
   
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hRho_HM_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhRhoTTJinHM[ijj] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      if(!fMC) fOutput->Add((TH1D*) fhRhoTTJinHM[ijj]); 
   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hRho_MB_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinMB[igg] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1D*) fhRhoTTCinMB[igg]); 
   } 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hRho_HM_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinHM[igg] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT      
      if(!fMC) fOutput->Add((TH1D*) fhRhoTTCinHM[igg]); 
   }
    
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hRho_GA_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRhoTTCinGA[igg] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
      fOutput->Add((TH1D*) fhRhoTTCinGA[igg]); 
   } 

   if(fMC){
      name = Form("hRhoMBpart");
      fhRhoMBpart = (TH1D*)  fhRhoMB->Clone(name.Data()); 
      fhRhoMBpart->SetTitle("Rho min bias part level"); 
      fOutput->Add((TH1D*) fhRhoMBpart); 

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){
         name = Form("hRho_MB_TTH%d_%d_part", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRhoTTHinMBpart[itt] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTHinMBpart[itt]); 
      }

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hRho_MB_TTC%d_%d_part", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhRhoTTCinMBpart[igg] = (TH1D*)  fhRhoMB->Clone(name.Data());                      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTCinMBpart[igg]); 
      }

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

   Int_t io = 0; 
   // 	FILTER_p-p_208_LHC16d  
   fRuns[io] = 252235;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252248;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252271;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252310;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252317;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252319;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252322;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252325;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252326;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;
   fRuns[io] = 252330;  fMeanV0A[io] = 50.4618;  fMeanV0C[io] = 77.4992;  fMeanV0M[io] = 127.961;  io++;

   //FILTER_p-p_208_LHC16e
   fRuns[io] = 253437;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253478;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253481;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253482;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253488;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253517;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253529;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253530;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253563;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253589;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;
   fRuns[io] = 253591;  fMeanV0A[io] = 49.9215;  fMeanV0C[io] = 75.7266;  fMeanV0M[io] = 125.648;  io++;

   // 	FILTER_p-p_208_LHC16g 
   fRuns[io] = 254128;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254147;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254149;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254174;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254175;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254178;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254193;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254199;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254204;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254205;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254293;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254302;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254303;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254304;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254330;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254331;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;
   fRuns[io] = 254332;  fMeanV0A[io] = 44.3838;  fMeanV0C[io] = 68.8414;  fMeanV0M[io] = 113.225;  io++;

   //FILTER_p-p_208_LHC16h 
   fRuns[io] = 254604;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254606;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254621;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254629;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254630;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254632;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254640;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254644;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254646;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254648;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254649;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254651;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254652;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254653;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254654;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254983;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 254984;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255079;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255082;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255085;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255086;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255091;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255111;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255154;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255159;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255162;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255167;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255171;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255173;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255174;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255176;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255177;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255180;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255181;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255182;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255240;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255242;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255247;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255248;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255249;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255251;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255252;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255253;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255255;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255256;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255275;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255276;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255280;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255283;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255350;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255351;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255352;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255398;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255402;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255407;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255415;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255418;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255419;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255420;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255421;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255440;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255442;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255447;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255463;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255465;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255466;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;
   fRuns[io] = 255467;  fMeanV0A[io] = 40.8315;  fMeanV0C[io] = 63.1392;  fMeanV0M[io] = 103.971;  io++;

   //FILTER_p-p_208_LHC16i 
   fRuns[io] = 255539;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255540;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255541;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255542;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255543;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255577;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255582;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255583;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255591;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255614;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255615;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255616;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255617;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;
   fRuns[io] = 255618;  fMeanV0A[io] = 36.9083;  fMeanV0C[io] = 58.1123;  fMeanV0M[io] = 95.0206;  io++;

   //FILTER_p-p_208_LHC16j
   fRuns[io] = 256219;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256223;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256227;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256228;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256231;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256281;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256282;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256283;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256284;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256287;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256289;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256290;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256292;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256295;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256297;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256299;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256302;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256307;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256309;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256311;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256356;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256361;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256362;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256363;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256364;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256365;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256366;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256368;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256371;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256372;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256373;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256415;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256417;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;
   fRuns[io] = 256418;  fMeanV0A[io] = 36.0796;  fMeanV0C[io] = 57.1701;  fMeanV0M[io] = 93.2497;  io++;

   // FILTER_p-p_208_LHC16k
   fRuns[io] = 256941;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 256942;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 256944;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257011;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257012;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257021;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257026;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257028;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257077;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257080;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257082;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257084;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257086;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257092;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257095;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257100;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257136;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257137;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257138;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257139;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257141;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257144;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257204;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257206;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257209;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257224;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257260;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257318;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257320;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257322;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257330;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257358;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257364;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257433;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257457;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257468;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257474;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257487;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257488;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257490;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257491;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257492;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257530;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257531;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257537;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257539;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257540;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257541;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257560;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257561;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257562;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257566;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257587;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257588;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257590;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257592;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257594;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257595;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257601;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257604;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257605;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257606;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257630;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257632;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257635;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257636;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257642;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257644;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257682;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257684;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257685;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257687;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257688;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257689;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257691;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257692;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257694;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257697;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257724;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257725;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257727;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257733;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257734;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257735;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257737;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257754;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257757;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257765;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257773;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257797;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257798;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257799;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257800;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257803;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257804;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257850;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257851;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257853;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257855;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257892;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257936;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257937;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257939;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257957;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257960;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257963;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257979;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257986;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257989;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 257992;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258003;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258008;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258012;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258014;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258017;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258019;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258039;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258041;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258042;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258045;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258049;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258053;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258059;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258060;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258062;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258063;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258107;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258108;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258109;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258113;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258114;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258117;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258178;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258197;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258198;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258202;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258203;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258204;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258256;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258257;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258258;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258270;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258271;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258273;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258274;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258278;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258299;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258301;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258302;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258303;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258306;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258307;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258332;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258336;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258359;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258387;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258391;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258393;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258426;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258452;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258454;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258456;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258477;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258499;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;
   fRuns[io] = 258537;  fMeanV0A[io] = 34.0369;  fMeanV0C[io] = 54.9189;  fMeanV0M[io] = 88.9558;  io++;

   // FILTER_p-p_208_LHC16l 
   fRuns[io] = 258919;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 258923;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 258962;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 258964;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259088;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259090;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259091;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259096;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259099;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259117;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259118;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259162;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259204;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259257;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259261;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259263;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259264;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259269;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259270;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259271;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259272;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259273;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259274;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259302;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259303;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259305;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259307;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259334;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259336;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259339;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259340;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259341;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259342;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259378;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259382;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259388;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259389;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259394;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259395;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259396;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259473;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259477;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259747;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259748;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259750;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259751;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259752;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259756;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259781;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259788;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259789;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259822;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259841;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259842;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259860;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259866;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259867;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259868;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;
   fRuns[io] = 259888;  fMeanV0A[io] = 33.7389;  fMeanV0C[io] = 54.8508;  fMeanV0M[io] = 88.5897;  io++;

   // FILTER_p-p_208_LHC16o 
   fRuns[io] = 262424;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262425;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262426;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262428;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262705;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262706;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262708;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262713;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262717;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262719;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262723;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262725;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262727;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262760;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262768;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262776;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262777;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262778;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262841;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262842;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262844;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262847;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262849;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262853;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262855;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 262858;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263331;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263332;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263487;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263490;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263496;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263497;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263529;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263647;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263652;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263654;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263657;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263662;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263663;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263682;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263690;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263691;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263737;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263738;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263739;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263741;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263743;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263744;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263784;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263785;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263786;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263787;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263790;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263792;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263793;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263803;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263810;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263863;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263866;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263905;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263916;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263917;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263920;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263923;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263977;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263978;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263981;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263984;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 263985;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 264033;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;
   fRuns[io] = 264035;  fMeanV0A[io] = 32.4107;  fMeanV0C[io] = 52.7467;  fMeanV0M[io] = 85.1573;  io++;


   // FILTER_p-p_208_LHC16p
   fRuns[io] = 264076;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264078;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264082;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264085;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264086;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264109;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264110;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264129;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264137;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264138;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264139;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264164;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264168;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264188;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264190;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264194;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264197;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264198;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264232;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264233;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264235;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264238;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264259;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264260;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264261;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264262;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264264;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264265;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264266;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264267;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264273;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264277;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264279;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264281;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264305;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264306;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264312;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264336;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264341;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264345;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264346;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;
   fRuns[io] = 264347;  fMeanV0A[io] = 53.9016;  fMeanV0C[io] = 82.9246;  fMeanV0M[io] = 136.826;  io++;

   // FILTER_p-p_208_LHC17c
   fRuns[io] = 270581;  fMeanV0A[io] = 53.2824;  fMeanV0C[io] = 79.1717;  fMeanV0M[io] = 132.454;  io++;
   fRuns[io] = 270661;  fMeanV0A[io] = 53.2824;  fMeanV0C[io] = 79.1717;  fMeanV0M[io] = 132.454;  io++;
   fRuns[io] = 270663;  fMeanV0A[io] = 53.2824;  fMeanV0C[io] = 79.1717;  fMeanV0M[io] = 132.454;  io++;
   fRuns[io] = 270665;  fMeanV0A[io] = 53.2824;  fMeanV0C[io] = 79.1717;  fMeanV0M[io] = 132.454;  io++;
   fRuns[io] = 270667;  fMeanV0A[io] = 53.2824;  fMeanV0C[io] = 79.1717;  fMeanV0M[io] = 132.454;  io++;

   // FILTER_p-p_208_LHC17e
   fRuns[io] = 270822;  fMeanV0A[io] = 53.8262;  fMeanV0C[io] = 79.4514;  fMeanV0M[io] = 133.277;  io++;
   fRuns[io] = 270824;  fMeanV0A[io] = 53.8262;  fMeanV0C[io] = 79.4514;  fMeanV0M[io] = 133.277;  io++;
   fRuns[io] = 270827;  fMeanV0A[io] = 53.8262;  fMeanV0C[io] = 79.4514;  fMeanV0M[io] = 133.277;  io++;
   fRuns[io] = 270828;  fMeanV0A[io] = 53.8262;  fMeanV0C[io] = 79.4514;  fMeanV0M[io] = 133.277;  io++;
   fRuns[io] = 270830;  fMeanV0A[io] = 53.8262;  fMeanV0C[io] = 79.4514;  fMeanV0M[io] = 133.277;  io++;

   // FILTER_p-p_208_LHC17f 
   fRuns[io] = 270854;  fMeanV0A[io] = 53.6570;  fMeanV0C[io] = 78.4221;  fMeanV0M[io] = 132.079;  io++;
   fRuns[io] = 270855;  fMeanV0A[io] = 53.6570;  fMeanV0C[io] = 78.4221;  fMeanV0M[io] = 132.079;  io++;
   fRuns[io] = 270856;  fMeanV0A[io] = 53.6570;  fMeanV0C[io] = 78.4221;  fMeanV0M[io] = 132.079;  io++;
   fRuns[io] = 270861;  fMeanV0A[io] = 53.6570;  fMeanV0C[io] = 78.4221;  fMeanV0M[io] = 132.079;  io++;
   fRuns[io] = 270865;  fMeanV0A[io] = 53.6570;  fMeanV0C[io] = 78.4221;  fMeanV0M[io] = 132.079;  io++;

   // FILTER_p-p_208_LHC17h 
   fRuns[io] = 271870;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 271871;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 271873;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 271874;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 271880;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 271886;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272018;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272020;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272036;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272038;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272039;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272040;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272042;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272076;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272100;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272101;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272123;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272151;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272152;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272153;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272154;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272155;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272156;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272194;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272335;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272340;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272359;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272360;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272388;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272389;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272394;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272395;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272399;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272400;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272411;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272413;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272461;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272462;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272463;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272466;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272468;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272521;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272574;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272575;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272577;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272585;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272607;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272608;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272610;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272620;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272690;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272691;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272712;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272747;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272749;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272760;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272763;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272764;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272782;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272783;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272784;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272828;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272829;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272833;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272834;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272836;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272870;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272871;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272873;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272880;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272903;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272905;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272932;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272933;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272934;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272935;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272939;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272947;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272949;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272976;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272983;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 272985;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273009;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273010;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273077;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273099;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273100;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;
   fRuns[io] = 273103;  fMeanV0A[io] = 51.1986;  fMeanV0C[io] = 74.6724;  fMeanV0M[io] = 125.871;  io++;

   // FILTER_p-p_208_LHC17i 
   fRuns[io] = 273591;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273592;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273593;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273653;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273654;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273824;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273825;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273885;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273886;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273887;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273889;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273918;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273942;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273943;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273946;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273985;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 273986;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274058;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274092;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274094;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274125;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274147;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274148;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274174;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274212;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274232;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274258;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274259;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274263;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274264;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274266;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274268;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274269;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274270;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274271;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274276;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274278;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274280;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274281;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274283;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274329;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274352;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274360;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274363;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274364;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274385;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274386;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274387;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274388;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274389;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274390;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;
   fRuns[io] = 274442;  fMeanV0A[io] = 49.4576;  fMeanV0C[io] = 72.9966;  fMeanV0M[io] = 122.454;  io++;

   // FILTER_p-p_208_LHC17j 
   fRuns[io] = 274593;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274594;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274595;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274596;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274601;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274653;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274657;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274667;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274669;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;
   fRuns[io] = 274671;  fMeanV0A[io] = 49.3295;  fMeanV0C[io] = 72.7080;  fMeanV0M[io] = 122.038;  io++;

   // FILTER_p-p_208_LHC17k
   fRuns[io] = 274690;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274708;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274801;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274802;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274803;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274806;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274815;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274821;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274822;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274877;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274878;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274882;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274886;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274978;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 274979;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275067;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275068;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275073;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275075;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275076;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275149;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275150;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275151;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275173;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275174;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275177;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275180;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275184;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275188;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275239;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275245;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275246;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275247;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275283;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275314;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275322;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275324;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275326;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275328;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275332;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275333;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275360;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275361;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275369;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275372;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275401;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275404;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275406;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275443;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275448;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275452;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275453;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275456;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275457;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275459;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275467;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275471;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275472;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275515;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275558;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275559;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275612;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275617;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275621;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275622;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275623;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275624;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275647;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275648;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275650;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275661;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275664;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 275847;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276097;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276098;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276099;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276102;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276104;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276135;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276140;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276145;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276166;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276169;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276170;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276177;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276178;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276205;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276230;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276257;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276259;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276290;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276292;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276294;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276297;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276302;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276348;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276351;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276435;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276437;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276438;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276439;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276462;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276506;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276507;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;
   fRuns[io] = 276508;  fMeanV0A[io] = 48.2257;  fMeanV0C[io] = 71.2655;  fMeanV0M[io] = 119.491;  io++;

   // FILTER_p-p_208_LHC17l 
   fRuns[io] = 276551;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276552;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276553;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276556;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276557;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276608;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276644;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276670;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276671;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276672;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276674;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276675;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276762;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276916;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276917;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276920;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276967;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276969;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276970;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276971;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 276972;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277015;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277016;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277017;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277037;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277073;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277076;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277079;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277082;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277087;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277091;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277117;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277121;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277155;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277180;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277181;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277182;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277183;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277184;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277188;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277189;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277193;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277194;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277196;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277197;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277256;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277257;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277262;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277293;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277310;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277312;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277314;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277360;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277383;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277384;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277385;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277386;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277389;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277416;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277417;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277418;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277472;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277473;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277476;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277477;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277478;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277479;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277530;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277531;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277534;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277536;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277537;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277574;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277575;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277576;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277577;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277721;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277722;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277723;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277725;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277745;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277746;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277747;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277749;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277794;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277795;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277799;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277800;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277801;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277802;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277805;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277834;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277836;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277841;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277842;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277845;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277847;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277848;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277870;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277876;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277897;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277898;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277899;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277900;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277903;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277904;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277907;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277930;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277952;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277987;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277989;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277991;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 277996;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278121;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278122;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278123;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278126;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278127;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278158;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278164;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278165;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278166;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278167;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278189;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278191;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278215;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;
   fRuns[io] = 278216;  fMeanV0A[io] = 47.3948;  fMeanV0C[io] = 69.9313;  fMeanV0M[io] = 117.326;  io++;

   // FILTER_p-p_208_LHC17m 
   fRuns[io] = 278914;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278915;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278936;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278939;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278941;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278959;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278960;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278963;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278964;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 278999;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279000;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279005;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279007;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279008;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279035;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279036;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279041;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279043;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279044;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279068;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279069;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279073;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279074;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279075;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279106;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279107;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279117;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279118;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279122;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279123;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279130;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279155;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279157;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279199;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279201;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279207;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279208;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279232;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279234;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279235;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279238;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279242;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279264;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279265;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279267;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279268;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279270;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279273;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279274;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279309;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279310;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279312;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279342;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279344;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279348;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279349;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279354;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279355;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279391;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279410;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279435;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279439;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279441;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279483;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279487;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279488;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279491;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279550;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279559;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279630;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279632;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279641;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279642;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279676;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279677;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279679;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279682;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279683;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279684;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279687;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279688;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279689;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279715;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279718;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279719;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279747;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279749;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279773;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279826;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279827;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279830;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279853;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279854;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279855;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 279879;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280051;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280052;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280066;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280107;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280108;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280111;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280114;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280118;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280126;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280131;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280134;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280135;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;
   fRuns[io] = 280140;  fMeanV0A[io] = 46.6462;  fMeanV0C[io] = 68.8063;  fMeanV0M[io] = 115.452;  io++;

   // FILTER_p-p_208_LHC17o 
   fRuns[io] = 280282;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280284;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280285;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280286;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280290;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280310;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280312;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280348;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280349;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280350;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280351;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280374;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280375;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280403;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280405;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280406;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280412;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280415;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280419;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280443;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280445;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280446;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280447;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280448;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280490;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280499;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280518;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280519;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280546;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280547;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280550;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280551;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280574;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280581;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280583;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280613;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280634;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280636;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280637;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280639;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280645;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280647;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280671;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280679;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280681;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280705;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280706;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280729;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280753;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280754;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280755;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280756;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280757;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280761;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280762;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280763;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280764;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280765;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280766;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280767;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280768;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280786;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280787;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280792;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280793;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280842;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280844;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280847;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280848;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280849;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280854;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280856;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280880;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280897;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280936;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280940;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280943;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280947;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280990;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280994;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280996;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280997;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280998;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 280999;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281032;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281033;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281035;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281036;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281060;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281061;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281062;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281080;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281081;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281179;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281180;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281181;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281189;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281190;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281191;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281212;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281213;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281240;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281241;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281242;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281243;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281244;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281271;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281273;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281275;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281277;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281301;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281321;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281415;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281441;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281443;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281444;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281446;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281449;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281450;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281475;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281477;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281509;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281511;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281557;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281562;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281563;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281568;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281569;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281574;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281583;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281592;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281633;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281892;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281893;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281894;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281895;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281915;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281916;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281918;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281920;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281928;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281931;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281932;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281939;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281940;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281953;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281956;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;
   fRuns[io] = 281961;  fMeanV0A[io] = 45.0074;  fMeanV0C[io] = 66.7215;  fMeanV0M[io] = 111.729;  io++;

   // FILTER_p-p_208_LHC17r 
   fRuns[io] = 282528;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282544;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282545;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282546;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282573;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282575;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282579;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282580;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282606;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282607;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282608;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282609;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282618;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282620;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282622;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282629;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282651;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282666;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282667;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282670;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282671;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282673;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282676;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282677;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282700;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282702;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282703;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;
   fRuns[io] = 282704;  fMeanV0A[io] = 44.2918;  fMeanV0C[io] = 65.4901;  fMeanV0M[io] = 109.782;  io++;


   // run by run V0M 
   fnRun = 1171; //the number of runs

   if(fnRun!=io){
      PostData(1, fOutput);
      return;
   }

   if(fMC){   //MC V0 amplitudes do not seem to exhibit period by period variations
      for(Int_t i=0; i<fnRun; i++){
         fMeanV0A[i] = 40.2292;  fMeanV0C[i] = 56.0226;  fMeanV0M[i] = 96.2518;  
      }
   }



 
  fhV0MRunByRunMB = new TH2D("fhV0MRunByRunMB","fhV0MRunByRunMB", fnRun, 0, fnRun, 180,0,1800); 
   for(Int_t ir=0; ir < fnRun; ir++){
      fhV0MRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d",fRuns[ir]));
   } 
   fOutput->Add((TH2D*) fhV0MRunByRunMB);
 
   name = "fhV0ARunByRunMB";
   fhV0ARunByRunMB = (TH2D*)  fhV0MRunByRunMB->Clone(name.Data());
   fhV0ARunByRunMB->SetTitle(name.Data());
   fOutput->Add((TH2D*) fhV0ARunByRunMB);
 
   name = "fhV0CRunByRunMB";
   fhV0CRunByRunMB = (TH2D*)  fhV0MRunByRunMB->Clone(name.Data());
   fhV0CRunByRunMB->SetTitle(name.Data());
   fOutput->Add((TH2D*) fhV0CRunByRunMB);
 
   fhV0MnormRunByRunMB = new TH2D("fhV0MnormRunByRunMB","fhV0MnormRunByRunMB", fnRun, 0, fnRun, 200,0,20); 
   for(Int_t ir=0; ir < fnRun; ir++){
      fhV0MnormRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d",fRuns[ir]));
   } 
   fOutput->Add((TH2D*) fhV0MnormRunByRunMB);
 

   TString cest[] = {"V0A", "V0C", "V0M", "V0Mnorm", "V0ACnorm"}; //centrality estimators

   const Int_t narrV0 = 1700;
   Double_t arrV0[narrV0+1];
   for(Int_t i=0; i<1600; i++){
      arrV0[i]=0.5*i;  //0-800
   }
   for(Int_t i=0; i<=100; i++){
      arrV0[1600+i] = 800 + 10.*i;  //800-1800
   }



   Double_t arrcent[] = {
     0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
     0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
     0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
     0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
     0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
     0.5,0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,  
     1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     2.,2.5,  
     3.,3.5,  
     4.,4.5,  
     5.,5.5,  
     6.,6.5,  
     7.,7.5,  
     8.,8.5,  
     9.,9.5,  
     10,11,12,13,14,15,16,17,18,19,20,
     25,30,35,40,45,50,55,60,65,70,75,80,85,90,100};

    Int_t narrcent = sizeof(arrcent)/sizeof(Double_t)-1;


   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      name = Form("hCentrality_MB_%s",cest[ic].Data());
      fhCentralityMB[ic] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
      fOutput->Add((TH2D*) fhCentralityMB[ic]); 
   }
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      name = Form("hCentrality_MB_%s",cest[ic].Data());
      fhCentralityMB[ic] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 400,0,20);
      fOutput->Add((TH2D*) fhCentralityMB[ic]);
   } 
 

   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      name = Form("hCentrality_HM_%s",cest[ic].Data());
      fhCentralityHM[ic] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
      if(!fMC) fOutput->Add((TH2D*) fhCentralityHM[ic]); 
   }

   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      name = Form("hCentrality_HM_%s",cest[ic].Data());
      fhCentralityHM[ic] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
      if(!fMC) fOutput->Add((TH2D*) fhCentralityHM[ic]); 
   }  
 
   //TTH MB
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hCentrality_MB_%s_TTH%d_%d",cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhCentralityTTHinMB[ic][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         fOutput->Add((TH2D*) fhCentralityTTHinMB[ic][itt]); 
      }
   }
 
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hCentrality_MB_%s_TTH%d_%d",cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhCentralityTTHinMB[ic][itt] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         if(!fMC) fOutput->Add((TH2D*) fhCentralityTTHinMB[ic][itt]); 
      }
   }
   //TTH HM 
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hCentrality_HM_%s_TTH%d_%d",cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhCentralityTTHinHM[ic][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         if(!fMC) fOutput->Add((TH2D*) fhCentralityTTHinHM[ic][itt]); 
      }
   }
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hCentrality_HM_%s_TTH%d_%d",cest[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhCentralityTTHinHM[ic][itt] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         if(!fMC) fOutput->Add((TH2D*) fhCentralityTTHinHM[ic][itt]); 
      } 
   }
   //TTJ MB
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hCentrality_MB_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhCentralityTTJinMB[ic][ijj] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         fOutput->Add((TH2D*) fhCentralityTTJinMB[ic][ijj]); 
      }
   }
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hCentrality_MB_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhCentralityTTJinMB[ic][ijj] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         fOutput->Add((TH2D*) fhCentralityTTJinMB[ic][ijj]); 
      }
   } 
   //TTJ HM
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hCentrality_HM_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhCentralityTTJinHM[ic][ijj] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         if(!fMC) fOutput->Add((TH2D*) fhCentralityTTJinHM[ic][ijj]); 
      }
   }
 
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hCentrality_HM_%s_TTJ%d_%d", cest[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhCentralityTTJinHM[ic][ijj] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         if(!fMC) fOutput->Add((TH2D*) fhCentralityTTJinHM[ic][ijj]); 
      }
   }
   //TTC  MB 
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hCentrality_MB_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhCentralityTTCinMB[ic][ijj] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         fOutput->Add((TH2D*) fhCentralityTTCinMB[ic][ijj]); 
      }
   }
   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hCentrality_MB_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhCentralityTTCinMB[ic][ijj] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         fOutput->Add((TH2D*) fhCentralityTTCinMB[ic][ijj]); 
      }
   }
   //TTC GA
   for(Int_t ic=0; ic<fkV0Mnorm1;ic++){
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hCentrality_GA_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhCentralityTTCinGA[ic][ijj] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         fOutput->Add((TH2D*) fhCentralityTTCinGA[ic][ijj]);
      }
   }

   for(Int_t ic=fkV0Mnorm1; ic<=fkV0Mnorm2; ic++){
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hCentrality_GA_%s_TTC%d_%d", cest[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhCentralityTTCinGA[ic][ijj] = (TH2D*) fhCentralityMB[ic]->Clone(name.Data());
         fOutput->Add((TH2D*) fhCentralityTTCinGA[ic][ijj]);
      }
   }

   //TString signal[]={"multV0A", "multV0C", "nTracklets", "znatower0", "znctower0","multV0M","fhNormSumV0AV0C"};
   //Float_t signalL[]={0,0,0,0,0,0,0};
   //Float_t signalH[]={1000,1000,500,30000,30000,1200,40};
   //Int_t signalN[]={1000,1000,500,100,100,1200,400};
   TString signal[]={"multV0A", "multV0C", "multV0M","multV0Mnorm","multV0ACnorm"};
   Float_t signalL[]={0,0,0,0,0};
   Float_t signalH[]={1000,1000,1800,15,15};
   Int_t   signalN[]={100,100,180,150,150};


   for(Int_t ic=0; ic<fkCE;ic++){ //MB
      name = Form("hSignal_MB_%s", signal[ic].Data());
      fhSignalMB[ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
      fOutput->Add((TH1D*) fhSignalMB[ic]); 
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM
      name = Form("hSignal_HM_%s",  signal[ic].Data());
      fhSignalHM[ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
      if(!fMC) fOutput->Add((TH1D*) fhSignalHM[ic]); 
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT hadron
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hSignal_MB_%s_TTH%d_%d", signal[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhSignalTTHinMB[ic][itt] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTHinMB[ic][itt]); 
      }
   }
 
   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT hadron
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hSignal_HM_%s_TTH%d_%d",  signal[ic].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhSignalTTHinHM[ic][itt] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         if(!fMC) fOutput->Add((TH1D*) fhSignalTTHinHM[ic][itt]); 
      }
   }
   
   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT jet
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hSignal_MB_%s_TTJ%d_%d", signal[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhSignalTTJinMB[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTJinMB[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT jet 
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hSignal_HM_%s_TTJ%d_%d", signal[ic].Data(), fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
         fhSignalTTJinHM[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         if(!fMC) fOutput->Add((TH1D*) fhSignalTTJinHM[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //MB && TT emcal 
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_MB_%s_TTC%d_%d", signal[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinMB[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTCinMB[ic][ijj]); 
      }
   }

   for(Int_t ic=0; ic<fkCE;ic++){ //HM && TT emcal 
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_HM_%s_TTC%d_%d", signal[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinHM[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         if(!fMC) fOutput->Add((TH1D*) fhSignalTTCinHM[ic][ijj]); 
      }
   }
 
   for(Int_t ic=0; ic<fkCE;ic++){ //GA && TT emcal
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("hSignal_GA_%s_TTC%d_%d", signal[ic].Data(), fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhSignalTTCinGA[ic][ijj] = new TH1D(name.Data(),name.Data(),signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignalTTCinGA[ic][ijj]); 
      }
   }

   if(fMC){ //PARTICLE LEVEL SIGNAL DISTRIBUTIONS
      //TString signalmc[]={"multV0A", "multV0C", "multV0M","fhNormSumV0AV0C"};
      //Float_t signalLmc[]={0,0,0,0};
      //Float_t signalHmc[]={500,500,500,40};
      //Int_t signalNmc[]={500,500,500,400};
      TString signalmc[]={"multV0A", "multV0C", "multV0M", "multV0Mnorm","multV0ACnorm"};
      Float_t signalLmc[]={0,0,0,0,0};
      Float_t signalHmc[]={500,500,500,20,20};
      Int_t signalNmc[]={500,500,500,200,200};
      
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_MB_%s_PartLevel", signalmc[ic].Data());
         fhSignalMB_PartLevel[ic] = new TH1D(name.Data(), name.Data(), signalNmc[ic], signalLmc[ic], signalHmc[ic]);
         fOutput->Add((TH1D*) fhSignalMB_PartLevel[ic]); 
      }

      //TT hadron
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_MB_%s_TTH%d_%d_PartLevel", signalmc[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTHinMB_PartLevel[ic][itt] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
            fOutput->Add((TH1D*) fhSignalTTHinMB_PartLevel[ic][itt]); 
         }
      }

      //TT cluster
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            name = Form("hSignal_MB_%s_TTC%d_%d_PartLevel", signalmc[ic].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
            fhSignalTTCinMB_PartLevel[ic][igg] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
            fOutput->Add((TH1D*) fhSignalTTCinMB_PartLevel[ic][igg]); 
         }
      }
   }
 

   name = Form("fhV0AvsV0C_MB");
   fhV0AvsV0C = new TH2D(name.Data(),name.Data(),100,0,1000, 100,0,1000);
   fOutput->Add((TH2D*) fhV0AvsV0C); 

   name = Form("fhV0MvsV0Mnorm_MB");
   fhV0MvsV0Mnorm = new TH2D(name.Data(),name.Data(),100,0,40, 100,0,1200);
   fOutput->Add((TH2D*) fhV0MvsV0Mnorm); 


   name = Form("fhV0AvsSPD_MB");
   fhV0AvsSPD = new TH2D(name.Data(),name.Data(),100,0,500, 100,0,500);
   fOutput->Add((TH2D*) fhV0AvsSPD);

   name = Form("fhV0CvsSPD_MB");
   fhV0CvsSPD = new TH2D(name.Data(),name.Data(),100,0,500, 100,0,500);
   fOutput->Add((TH2D*) fhV0CvsSPD);

 
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("fhV0AvsV0C_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhV0AvsV0CTTH[itt] = (TH2D*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTH[itt]->SetTitle(name.Data());
      fOutput->Add((TH2D*) fhV0AvsV0CTTH[itt]); 
   } 
   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("fhV0AvsV0C_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhV0AvsV0CTTJ[ijj] =  (TH2D*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTJ[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2D*) fhV0AvsV0CTTJ[ijj]); 
   }
   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("fhV0AvsV0C_MB_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhV0AvsV0CTTCinMB[ijj] = (TH2D*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTCinMB[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2D*) fhV0AvsV0CTTCinMB[ijj]); 
   } 
   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("fhV0AvsV0C_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhV0AvsV0CTTCinGA[ijj] = (TH2D*)  fhV0AvsV0C->Clone(name.Data());
      fhV0AvsV0CTTCinGA[ijj]->SetTitle(name.Data()); 
      fOutput->Add((TH2D*) fhV0AvsV0CTTCinGA[ijj]); 
   } 

   //+++++++++++++++++++++++++++++++
   fhTrackMultMB = new TH2D("fhTrackMultMB","fhTrackMultMB", narrcent, arrcent, 1000, 0, 1000); 
   fOutput->Add((TH2D*) fhTrackMultMB); 

   fhTrackMultHM = new TH2D("fhTrackMultHM","fhTrackMultHM", narrcent, arrcent, 1000, 0, 1000); 
   if(!fMC) fOutput->Add((TH1D*) fhTrackMultHM); 

   fhMeanTrackPtMB = new TH2D("fhMeanTrackPtMB","fhMeanTrackPtMB", narrcent, arrcent, 100, 0, 20);
   fOutput->Add((TH1D*) fhMeanTrackPtMB); 

   fhMeanTrackPtHM = new TH2D("fhMeanTrackPtHM","fhMeanTrackPtHM", narrcent, arrcent, 100, 0, 20);
   if(!fMC) fOutput->Add((TH1D*) fhMeanTrackPtHM); 


   //Trigger track candidate multiplicity
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hMultTT_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhMultTTHinMB[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*)  fhMultTTHinMB[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hMultTT_HM_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhMultTTHinHM[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
      if(!fMC) fOutput->Add((TH1D*)  fhMultTTHinHM[itt]); 
   }

   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hMultTT_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhMultTTJinMB[ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhMultTTJinMB[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
      name = Form("hMultTT_HM_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
      fhMultTTJinHM[ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
      if(!fMC) fOutput->Add((TH1D*) fhMultTTJinHM[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_MB_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinMB[ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhMultTTCinMB[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_HM_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinHM[ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
      if(!fMC) fOutput->Add((TH1D*) fhMultTTCinHM[ijj]); 
   }

   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
      name = Form("hMultTT_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
      fhMultTTCinGA[ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
      fOutput->Add((TH1D*) fhMultTTCinGA[ijj]); 
   }

   //Trigger track pT spectrum single inclusive for  MB  versus V0M
   Int_t    nbinsV0M     = 100;
   Double_t maxV0M       = 1000;
   Double_t maxV0Mmc     = 500;
   Int_t    nbinsV0Mnorm = 200;
   Double_t maxV0Mnorm   = 20;
   //Double_t maxCentV0M   = 100;
 
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_MB_TTH%d_%d_V0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinMB_V0M[itt] = new TH2D(name.Data(),name.Data(), nbinsV0M, 0, maxV0M, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTHinMB_V0M[itt]); 
   }

  for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_MB_TTH%d_%d_CentV0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinMB_CentV0M[itt] = new TH2D(name.Data(),name.Data(),  narrcent, arrcent, 1000, 0, 100);
      fOutput->Add((TH2D*) fhTTHinMB_CentV0M[itt]); 
   }


   //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_MB_TTH%d_%d_V0Mnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinMB_V0Mnorm1[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTHinMB_V0Mnorm1[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_MB_TTH%d_%d_V0ACnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinMB_V0Mnorm2[itt] = (TH2D*) fhTTHinMB_V0Mnorm1[itt]->Clone(name.Data());
      fOutput->Add((TH2D*) fhTTHinMB_V0Mnorm2[itt]); 
   }

   if(fMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0M_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTHinMB_V0M_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0M, 0, maxV0Mmc, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTHinMB_V0M_PartLevel[itt]); 
      }
      
      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTHinMB_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTHinMB_V0Mnorm1_PartLevel[itt]); 
      }

      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0ACnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTHinMB_V0Mnorm2_PartLevel[itt] = (TH2D*) fhTTHinMB_V0Mnorm1_PartLevel[itt]->Clone(name.Data());
         fOutput->Add((TH2D*) fhTTHinMB_V0Mnorm2_PartLevel[itt]); 
      }

   }


   //Trigger track pT spectrum single inclusive for  HM  versus V0M
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_HM_TTH%d_%d_V0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinHM_V0M[itt] = new TH2D(name.Data(),name.Data(), nbinsV0M, 0, maxV0M, 100, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTHinHM_V0M[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_HM_TTH%d_%d_CentV0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinHM_CentV0M[itt] = new TH2D(name.Data(),name.Data(),  narrcent, arrcent, 1000, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTHinHM_CentV0M[itt]); 
   }


   //Trigger track pT spectrum single inclusive for  HM  versus V0Mnorm
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_HM_TTH%d_%d_V0Mnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinHM_V0Mnorm1[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTHinHM_V0Mnorm1[itt]); 
   }
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      name = Form("hTT_HM_TTH%d_%d_V0ACnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhTTHinHM_V0Mnorm2[itt] =  (TH2D*) fhTTHinHM_V0Mnorm1[itt]->Clone(name.Data());
      if(!fMC) fOutput->Add((TH2D*) fhTTHinHM_V0Mnorm2[itt]); 
   }



   //TT emcal cluster pT spectrum single inclusive  in MB   with V0M
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_MB_TTC%d_%d_V0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinMB_V0M[igg] = new TH2D(name.Data(),name.Data(), nbinsV0M, 0, maxV0M, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTCinMB_V0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_MB_TTC%d_%d_CentV0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinMB_CentV0M[igg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
      fOutput->Add((TH2D*) fhTTCinMB_CentV0M[igg]); 
   }


   //TT emcal cluster pT spectrum single inclusive  in MB   with V0Mnorm
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_MB_TTC%d_%d_V0Mnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinMB_V0Mnorm1[igg] = new TH2D(name.Data(),name.Data(),  nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTCinMB_V0Mnorm1[igg]); 
   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_MB_TTC%d_%d_V0ACnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinMB_V0Mnorm2[igg] = (TH2D*) fhTTCinMB_V0Mnorm1[igg]->Clone(name.Data());
      fOutput->Add((TH2D*) fhTTCinMB_V0Mnorm2[igg]); 
   }
  
   if(fMC){
      //TT emcal cluster pT spectrum single inclusive  in MB   with V0M
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hTT_MB_TTC%d_%d_V0M_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhTTCinMB_V0M_PartLevel[igg] = new TH2D(name.Data(),name.Data(), nbinsV0M, 0, maxV0Mmc, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTCinMB_V0M_PartLevel[igg]); 
      }
     
      //TT emcal cluster pT spectrum single inclusive  in MB   with V0Mnorm
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hTT_MB_TTC%d_%d_V0Mnorm_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhTTCinMB_V0Mnorm1_PartLevel[igg] = new TH2D(name.Data(),name.Data(),  nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTCinMB_V0Mnorm1_PartLevel[igg]); 
      }

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hTT_MB_TTC%d_%d_V0ACnorm_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhTTCinMB_V0Mnorm2_PartLevel[igg] = (TH2D*) fhTTCinMB_V0Mnorm1_PartLevel[igg]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhTTCinMB_V0Mnorm2_PartLevel[igg]); 
      }
   }


   //TT emcal cluster pT spectrum single inclusive  in HM   with V0M
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_HM_TTC%d_%d_V0M",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinHM_V0M[igg] = new TH2D(name.Data(), name.Data(), nbinsV0M, 0, maxV0M, 100, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTCinHM_V0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_HM_TTC%d_%d_CentV0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinHM_CentV0M[igg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTCinHM_CentV0M[igg]); 
   }


   //TT emcal cluster pT spectrum single inclusive  in HM   with V0Mnorm
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_HM_TTC%d_%d_V0Mnorm",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinHM_V0Mnorm1[igg] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
      if(!fMC) fOutput->Add((TH2D*) fhTTCinHM_V0Mnorm1[igg]); 
   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_HM_TTC%d_%d_V0ACnorm",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinHM_V0Mnorm2[igg] = (TH2D*) fhTTCinHM_V0Mnorm1[igg]->Clone(name.Data());
      if(!fMC) fOutput->Add((TH2D*) fhTTCinHM_V0Mnorm2[igg]); 
   }


   //TT emcal cluster pT spectrum single inclusive  in GA   with V0M
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_GA_TTC%d_%d_V0M",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinGA_V0M[igg] = new TH2D(name.Data(), name.Data(), nbinsV0M, 0, maxV0M, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTCinGA_V0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_GA_TTC%d_%d_CentV0M",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinGA_CentV0M[igg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
      fOutput->Add((TH2D*) fhTTCinGA_CentV0M[igg]); 
   }


   //TT emcal cluster pT spectrum single inclusive  in GA   with V0Mnorm
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_GA_TTC%d_%d_V0Mnorm",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinGA_V0Mnorm1[igg] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
      fOutput->Add((TH2D*) fhTTCinGA_V0Mnorm1[igg]); 
   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("hTT_GA_TTC%d_%d_V0ACnorm",  fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhTTCinGA_V0Mnorm2[igg] = (TH2D*) fhTTCinGA_V0Mnorm1[igg]->Clone(name.Data());
      fOutput->Add((TH2D*) fhTTCinGA_V0Mnorm2[igg]); 
   }



   //RECOIL JET SPECTRA
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinMB_V0M[itt] = new TH2D(name.Data(), name.Data(), nbinsV0M, 0, maxV0M, 200, -20, 180);            
      fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0M[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      name = Form("fhRecoilJetPt_MB_TTH%d_%d_CentV0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinMB_CentV0M[itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);            
      fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_CentV0M[itt]); 
   }


   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
      name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0Mnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinMB_V0Mnorm1[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);            
      fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[itt]); 
   }
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
      name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0ACnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinMB_V0Mnorm2[itt] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[itt]->Clone(name.Data());
      fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm2[itt]); 
   }


   if(fMC){ // particle level 
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
         name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0M_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRecoilJetPtTTHinMB_V0M_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0M, 0, maxV0Mmc, 200, -20, 180);            
         fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0M_PartLevel[itt]); 
      }
      
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);            
         fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[itt]); 
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT   with (V0A/norm+V0C/norm)/2 in MB  
         name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[itt] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[itt]->Clone(name.Data());
         fOutput->Add((TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[itt]); 
      }
   }


   for(Int_t itt=0; itt<fnHadronTTBins; itt++){         //! recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
      name = Form("fhRecoilJetPt_HM_TTH%d_%d_V0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinHM_V0M[itt] = (TH2D*) fhRecoilJetPtTTHinMB_V0M[itt]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTHinHM_V0M[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){         //! recoil jets associated to semi-inclusive hadron TT  in HM  with V0M centrality
      name = Form("fhRecoilJetPt_HM_TTH%d_%d_CentV0M", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinHM_CentV0M[itt] = (TH2D*) fhRecoilJetPtTTHinMB_CentV0M[itt]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTHinHM_CentV0M[itt]); 
   }


   for(Int_t itt=0; itt<fnHadronTTBins; itt++){         //! recoil jets associated to semi-inclusive hadron TT  in HM  with V0Mnorm
      name = Form("fhRecoilJetPt_HM_TTH%d_%d_V0Mnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinHM_V0Mnorm1[itt] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[itt]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTHinHM_V0Mnorm1[itt]); 
   }
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){         //! recoil jets associated to semi-inclusive hadron TT  in HM  with (V0A/norm+V0C/norm)/2
      name = Form("fhRecoilJetPt_HM_TTH%d_%d_V0ACnorm", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhRecoilJetPtTTHinHM_V0Mnorm2[itt] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[itt]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTHinHM_V0Mnorm2[itt]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
      name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinMB_V0M[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M centrality
      name = Form("fhRecoilJetPt_MB_TTC%d_%d_CentV0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinMB_CentV0M[igg] = (TH2D*) fhRecoilJetPtTTHinMB_CentV0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_CentV0M[igg]); 
   }


   for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
      name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0Mnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinMB_V0Mnorm1[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0Mnorm1[igg]); 
   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
      name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0ACnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinMB_V0Mnorm2[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0Mnorm2[igg]); 
   }

   if(fMC){ //particle level
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
         name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0M_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhRecoilJetPtTTCinMB_V0M_PartLevel[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0M_PartLevel[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0M_PartLevel[igg]); 
      }
   
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
         name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0Mnorm_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[igg]); 
      }
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
         name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0ACnorm_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[igg] = (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[igg]); 
      }

   }
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_HM_TTC%d_%d_V0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinHM_V0M[igg] =  (TH2D*) fhRecoilJetPtTTHinMB_V0M[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTCinHM_V0M[igg]);   
   }
 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_HM_TTC%d_%d_CentV0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinHM_CentV0M[igg] =  (TH2D*) fhRecoilJetPtTTHinMB_CentV0M[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTCinHM_CentV0M[igg]);   
   }
 
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_HM_TTC%d_%d_V0Mnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinHM_V0Mnorm1[igg] =  (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTCinHM_V0Mnorm1[igg]);   
   }
    for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_HM_TTC%d_%d_V0ACnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinHM_V0Mnorm2[igg] =  (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhRecoilJetPtTTCinHM_V0Mnorm2[igg]);   
   }
   
   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_GA_TTC%d_%d_V0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinGA_V0M[igg] =   (TH2D*) fhRecoilJetPtTTHinMB_V0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinGA_V0M[igg]); 
   } 

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_GA_TTC%d_%d_CentV0M", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinGA_CentV0M[igg] =   (TH2D*) fhRecoilJetPtTTHinMB_CentV0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinGA_CentV0M[igg]); 
   }


   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_GA_TTC%d_%d_V0Mnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinGA_V0Mnorm1[igg] =   (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinGA_V0Mnorm1[igg]); 
   } 

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhRecoilJetPt_GA_TTC%d_%d_V0ACnorm", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhRecoilJetPtTTCinGA_V0Mnorm2[igg] =   (TH2D*) fhRecoilJetPtTTHinMB_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhRecoilJetPtTTCinGA_V0Mnorm2[igg]); 
   } 


   //delta pT distributions versus V0M CENTRALITY 
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      name = Form("fhDeltaPtTTHinMB_RC_CentV0M_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinMB_RC_CentV0M[itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);            
      fOutput->Add((TH2D*) fhDeltaPtTTHinMB_RC_CentV0M[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
      name = Form("fhDeltaPtTTHinHM_RC_CentV0M_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinHM_RC_CentV0M[itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);            
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTHinHM_RC_CentV0M[itt]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinMB_RC_CentV0M_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinMB_RC_CentV0M[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_CentV0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinMB_RC_CentV0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinHM_RC_CentV0M_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinHM_RC_CentV0M[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_CentV0M[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTCinHM_RC_CentV0M[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinGA_RC_CentV0M_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinGA_RC_CentV0M[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_CentV0M[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinGA_RC_CentV0M[igg]); 
   }

   //delta pT distributions versus V0Mnorm   = V0M/mean V0M
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
      name = Form("fhDeltaPtTTHinMB_RC_V0Mnorm_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinMB_RC_V0Mnorm1[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);            
      fOutput->Add((TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
      name = Form("fhDeltaPtTTHinHM_RC_V0Mnorm_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinHM_RC_V0Mnorm1[itt] = (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTHinHM_RC_V0Mnorm1[itt]); 
   }

    for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinMB_RC_V0Mnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinMB_RC_V0Mnorm1[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinMB_RC_V0Mnorm1[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinHM_RC_V0Mnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinHM_RC_V0Mnorm1[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTCinHM_RC_V0Mnorm1[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinGA_RC_V0Mnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinGA_RC_V0Mnorm1[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinGA_RC_V0Mnorm1[igg]); 
   }
 
   //delta pT distributions versus ( V0A/mean + V0C/mean )/2 
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
      name = Form("fhDeltaPtTTHinMB_RC_V0ACnorm_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinMB_RC_V0Mnorm2[itt] = (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm2[itt]); 
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
      name = Form("fhDeltaPtTTHinHM_RC_V0ACnorm_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
      fhDeltaPtTTHinHM_RC_V0Mnorm2[itt] = (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTHinHM_RC_V0Mnorm2[itt]); 
   }

    for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinMB_RC_V0ACnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinMB_RC_V0Mnorm2[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinMB_RC_V0Mnorm2[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinHM_RC_V0ACnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinHM_RC_V0Mnorm2[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      if(!fMC) fOutput->Add((TH2D*) fhDeltaPtTTCinHM_RC_V0Mnorm2[igg]); 
   }

   for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      name = Form("fhDeltaPtTTCinGA_RC_V0ACnorm_TTC%d_%d", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
      fhDeltaPtTTCinGA_RC_V0Mnorm2[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
      fOutput->Add((TH2D*) fhDeltaPtTTCinGA_RC_V0Mnorm2[igg]); 
   }
   
   if(fMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
         name = Form("fhDeltaPtTTHinMB_RC_V0Mnorm_TTH%d_%d_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[itt] = (TH2D*)  fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[itt]); 
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with (V0A/norm+V0C/norm)/2
         name = Form("fhDeltaPtTTHinMB_RC_V0ACnorm_TTH%d_%d_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[itt] = (TH2D*)  fhDeltaPtTTHinMB_RC_V0Mnorm1[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[itt]); 
      }

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("fhDeltaPtTTCinMB_RC_V0Mnorm_TTC%d_%d_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[igg]); 
      }
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("fhDeltaPtTTCinMB_RC_V0ACnorm_TTC%d_%d_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
         fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[igg] =   (TH2D*) fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[0]->Clone(name.Data()); 
         fOutput->Add((TH2D*) fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[igg]); 
      }

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

      fhJetPtPartLevelVsJetPtDetLevelZero = new TH2D("fhJetPtPartLevelVsJetPtDetLevelZero","fhJetPtPartLevelVsJetPtDetLevelZero",250,0,250,250,0,250);
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero);
 
      fhJetPtResolutionVsPtPartLevel = new TH2D("fhJetPtResolutionVsPtPartLevel","fhJetPtResolutionVsPtPartLevel",100,0,100,50,0,2);
      fOutput->Add((TH2D*) fhJetPtResolutionVsPtPartLevel);
   }


    fhOneOverPtVsPhiNeg = new TH2D("fhOneOverPtVsPhiNeg","1/pt versus track phi negative tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
    fOutput->Add((TH2D*) fhOneOverPtVsPhiNeg);
 
    fhOneOverPtVsPhiPos = new TH2D("fhOneOverPtVsPhiPos","1/pt versus track phi positive tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
    fOutput->Add((TH2D*) fhOneOverPtVsPhiPos);
 
    fhSigmaPtOverPtVsPt = new TH2D("fhSigmaPtOverPtVsPt",
                                       "track sigma(1/pt)/ 1/pt vs pt", 100, 0, 100, 250, 0, 1);
    fOutput->Add((TH2D*) fhSigmaPtOverPtVsPt); 

    Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
    Int_t nbins = sizeof(bins)/sizeof(Double_t)-1; //pT binning for DCA distribution

    fhDCAinXVsPt = new TH2D("fhDCAinXVsPt","fhDCAinXVsPt", nbins, bins, 200, -10.,10);
    fOutput->Add((TH2D*) fhDCAinXVsPt); 

    fhDCAinYVsPt = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPt");
    fOutput->Add((TH2D*) fhDCAinYVsPt);

    if(fMC){
       fhDCAinXVsPtPhysPrimary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinXVsPtPhysPrimary");
       fOutput->Add((TH2D*) fhDCAinXVsPtPhysPrimary);

       fhDCAinYVsPtPhysPrimary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPtPhysPrimary"); 
       fOutput->Add((TH2D*) fhDCAinYVsPtPhysPrimary); 
 
       fhDCAinXVsPtSecondary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinXVsPtSecondary");
       fOutput->Add((TH2D*) fhDCAinXVsPtSecondary); 

       fhDCAinYVsPtSecondary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPtSecondary");
       fOutput->Add((TH2D*) fhDCAinYVsPtSecondary); 
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

Double_t AliAnalysisTaskEA::GetDeltaPt(Double_t phiTT, Double_t etaTT, Double_t phiLJ, Double_t etaLJ, Double_t phiSJ, Double_t etaSJ, Double_t rho, Int_t level){

   Double_t rcEta = fRandom->Uniform( fJetContainerDetLevel->GetJetEtaMin(), fJetContainerDetLevel->GetJetEtaMax());
   Double_t rcPhi = fRandom->Uniform(0, TMath::TwoPi());
   Double_t jetR  = fJetContainerDetLevel->GetJetRadius(); 
   Double_t jetR2  = jetR*jetR; //square of jet R 
   Double_t exclR2 = 4*jetR2;  //rc axis has to be 2R far away from the LJ and SJ jet axis 
   Double_t dphirc=0, detarc=0; //distance of jet from the random cone

   Int_t irc=0;
   Bool_t isCloseLJ = kTRUE;
   Bool_t isCloseSJ = kTRUE;
   Bool_t isCloseTT = kTRUE;
   while(irc<1000){
      isCloseLJ = kTRUE;
      isCloseSJ = kTRUE;
      isCloseTT = kTRUE;
 
      if(etaLJ<10){

          dphirc = TVector2::Phi_mpi_pi(phiLJ - rcPhi);
          detarc = etaLJ - rcEta;

          if( dphirc*dphirc + detarc*detarc >  exclR2 ) isCloseLJ = kFALSE;
      }else{
         isCloseLJ = kFALSE;
      }
    
      if(etaSJ<10){
         dphirc = TVector2::Phi_mpi_pi(phiSJ - rcPhi);
         detarc = etaSJ - rcEta;

         if( dphirc*dphirc + detarc*detarc >  exclR2 ) isCloseSJ = kFALSE;
      }else{
         isCloseSJ = kFALSE; 
      }

      dphirc = TVector2::Phi_mpi_pi(phiTT - rcPhi);
      detarc = etaTT - rcEta;

      if( dphirc*dphirc + detarc*detarc >  exclR2 ) isCloseTT = kFALSE;

      if(!isCloseSJ && !isCloseLJ && !isCloseTT){
         //this random cone is far away from leading and subleading jet
         break;
      }else{ //generate a new random cone position
         rcEta = fRandom->Uniform( fJetContainerDetLevel->GetJetEtaMin(), fJetContainerDetLevel->GetJetEtaMax());
         rcPhi = fRandom->Uniform(0, TMath::TwoPi());
      }
      irc++;
   }

   Double_t sumptrc = 0.;
   AliVParticle *track = NULL; //jet constituent

   if(level == kDetLevel){
      for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         track = trackIterator.second;  // Get the full track
         if(!track) continue;
      
         if(IsTrackInAcceptance(track, kDetLevel)){  
            dphirc = TVector2::Phi_mpi_pi(track->Phi() - rcPhi);
            detarc = track->Eta() - rcEta;
      
            if( dphirc*dphirc + detarc*detarc <  jetR2 ){
                sumptrc +=  track->Pt();
            }
         }
      }
      //Delta pT  sum of momenta in the cone 
   }else{
      for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
         track = mcPartIterator.second;  // Get the pointer to mc particle object
         if(!track)  continue; 

         if(IsTrackInAcceptance(track, kPartLevel)){
            dphirc = TVector2::Phi_mpi_pi(track->Phi() - rcPhi);
            detarc = track->Eta() - rcEta;
      
            if( dphirc*dphirc + detarc*detarc <  jetR2 ){
                sumptrc +=  track->Pt();
            }
         } 
      }
   }

   return ( sumptrc - TMath::Pi()*jetR2*rho);
}

