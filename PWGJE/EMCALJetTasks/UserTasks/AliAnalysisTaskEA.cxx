#ifndef ALIANALYSISTASKSE_H

#include <Riostream.h>
#include <TROOT.h>
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
#include <TH2D.h>
#include <TH3D.h>
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
#include "AliAnalysisEmcalJetHelperEA.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

//#include "AliEmcalDownscaleFactorsOCDB.h"
//#include "AliEmcalAnalysisFactory.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEA)

using namespace PWGJE::EMCALJetTasks;
using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN PP 13 TeV 
// Author Filip Krizek   (8.Aug. 2019)

//________________________________________________________________________________________

AliAnalysisTaskEA::AliAnalysisTaskEA(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskEA", kTRUE),  
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyDetLevelContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyJetDetLevelContainerName(""),
fMyClusterContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName(""),
fMyKTJetDetLevelContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fTrkContainerDetLevelEMB(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fJetContainerDetLevelEMB(0x0),
fClusterContainerDetLevel(0x0),
fKTJetContainerDetLevel(0x0),
fKTJetContainerPartLevel(0x0),
fKTJetContainerDetLevelEMB(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsEmcalTrig(0),
fIsHighMultTrig(0),
fCentralityV0A(-1),
fCentralityV0C(-1),
fCentralityV0M(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
fNTracklets(-1),
fMultV0A(0.),
fMultV0C(0.),
fMultV0M(0.),
fMultV0Mnorm(0.),
fAsymV0M(999),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
fAsymV0M_PartLevel(999),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0), 
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhV0MAssymVsV0Mnorm_PartLevel(0x0),
fhV0AvsV0C(0x0),
//fhV0MvsV0Mnorm(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelZero(0x0),
fhV0ARunByRunMB(0x0),
fhV0CRunByRunMB(0x0),
fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
fhTrackEtaInclEMB(0x0),
fhJetPtPartLevelVsJetPtDetLevelZero(0x0),
fhOneOverPtVsPhiNeg(0x0),
fhOneOverPtVsPhiPos(0x0),
fhSigmaPtOverPtVsPt(0x0),
fhDCAinXVsPt(0x0),
fhDCAinYVsPt(0x0),
fhDCAinXVsPtPhysPrimary(0x0),
fhDCAinYVsPtPhysPrimary(0x0),
fhDCAinXVsPtSecondary(0x0),
fhDCAinYVsPtSecondary(0x0),
fMinFractionShared(0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fnJetChTTBins(0),
fnClusterTTBins(0),
fMode(AliAnalysisTaskEA::kNormal),
fFiducialCellCut(0x0),
fHelperEA(0x0),
fMeanV0M(1.), 
fMeanV0M_PartLevel(1.),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.)                            
{
   //default constructor
   for(Int_t ir=0; ir< kRho; ir++){
      fhRhoMBpart[ir]= 0x0 ;
      fhJetPtAreaV0norm_PartLevel[ir] = 0x0;
      fhJetPtPartLevelVsJetPtDetLevelCorr[ir] = 0x0;
      fhJetPtResolutionVsPtPartLevel[ir] = 0x0;
      fhJetPtPartLevelCorr[ir] = 0x0;
      fhFractionOfSecInJet[ir] = 0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;

      fhTrackMult[itg]=0x0;
      fhMeanTrackPt[itg]=0x0;

      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      fhClusterPhiIncl[itg] = 0x0;
      fhClusterEtaIncl[itg] = 0x0;

      fhTrackPtEtaPhiV0norm[itg] = 0x0;
      fhJetPtEtaPhiV0norm[itg] = 0x0;
      for(Int_t i=0; i<fkTTbins; i++){
         fhJetPtEtaPhiV0normTTH[itg][i] = 0x0;
      } 
 
      for(Int_t ir=0; ir<kRho; ir++){
         fhJetPtAreaV0norm[itg][ir] = 0x0;
         fhRho[itg][ir] = 0x0;
      }
    
      for(Int_t i=0; i<fkTTbins; i++){
         for(Int_t ir=0; ir<kRho; ir++){
            fhRhoTTH[itg][i][ir]=0x0;
            fhRhoTTC[itg][i][ir]=0x0; 
         }
         fhRhoTTJ[itg][i]=0x0;  
      }

      fhSharedJetFraction[itg] = 0x0;
      fhTrialsEMBtot[itg] = 0x0;
      fhXsectionEMBtot[itg] = 0x0;
      fhTrialsEMB[itg] = 0x0;
      fhXsectionEMB[itg] = 0x0;
      fhPtHardEMB[itg] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i]   = 0;
      fJetChTT[i]    = 0;
      fClusterTT[i]  = 0;


      fHadronTT_PartLevel[i]   = 0;
      fClusterTT_PartLevel[i]   = 0;

      //TT
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         fhMultTTH[itg][i] = 0x0;   
         fhMultTTJ[itg][i] = 0x0;  
         fhMultTTC[itg][i] = 0x0;  
         
         //fhTTH_CentV0M[itg][i]  = 0x0;
	 fhTTH_V0Mnorm1[itg][i] = 0x0;
         fhTTH_3D_V0Mnorm1[itg][i] = 0x0;

         //fhTTC_CentV0M[itg][i]  = 0x0;
         fhTTC_V0Mnorm1[itg][i] = 0x0;

         fhV0MAssymVsV0MnormTTH[itg][i] = 0x0;
      }
 
      fhTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTH_3D_V0Mnorm1_PartLevel[i] = 0x0;
               
      fhTTC_V0Mnorm1_PartLevel[i] = 0x0;

      fhV0MAssymVsV0MnormTTH_PartLevel[i] = 0x0;
         
      //RECOIL JET SPECTRA   
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         for(Int_t ir=0; ir<kRho; ir++){ 
            //fhRecoilJetPtTTH_CentV0M[itg][i][ir]  = 0x0;
            fhRecoilJetPtTTH_V0Mnorm1[itg][i][ir] = 0x0;
           
            fhRecoilJetPhiTTH_V0Mnorm1[itg][i][ir] = 0x0;
            fhRecoilJetTTH_V0Mnorm1[itg][i][ir]    = 0x0;
            
            //fhRecoilJetPtTTC_CentV0M[itg][i][ir]  = 0x0;
            fhRecoilJetPtTTC_V0Mnorm1[itg][i][ir] = 0x0;
         }
      }
         
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhRecoilJetPtTTH_V0Mnorm1_PartLevel[i][ir] = 0x0;
         fhRecoilJetPtTTC_V0Mnorm1_PartLevel[i][ir] = 0x0;
            
         fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[i][ir] = 0x0;
         fhRecoilJetTTH_V0Mnorm1_PartLevel[i][ir]    = 0x0;
      }      
   
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         //fhDeltaPtTTH_RC_CentV0M[itg][i] = 0x0;  
         //fhDeltaPtTTC_RC_CentV0M[itg][i] = 0x0;
        
         for(Int_t ir=0; ir<kRho; ir++){ 
            fhDeltaPtTTH_RC_V0Mnorm1[itg][i][ir] = 0x0;  
            fhDeltaPtTTC_RC_V0Mnorm1[itg][i][ir] = 0x0;
         }
      } 
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[i][ir] = 0x0;  
         fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[i][ir] = 0x0;
      }
   
      //remx 
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhJetPtPartLevelCorrTTHdl[i][ir] = 0x0;
         fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[i][ir] = 0x0;
      }   
      //embedding

      for(Int_t itg=kMB; itg<=kGA; itg++){
         for(Int_t ir=0; ir<kRho; ir++){ 
            fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][i][ir] = 0x0;
            fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][i][ir] = 0x0;
            
            fhJetPtPartLevelCorrTTHdl_EMB[itg][i][ir] = 0x0;
            
            fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[itg][i][ir] = 0x0;
            fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[itg][i][ir] = 0x0;
         }
         fhJetPtPartLevelZeroTTHdl_EMB[itg][i] = 0x0;
      }  
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir] = 0x0; 
         fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir] = 0x0; 
      }
   }


   for(Int_t itg=kMB; itg<=kGA; itg++){
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhJetPtPartLevelCorr_EMB[itg][ir] = 0x0;
      }      
      fhJetPtPartLevelZero_EMB[itg] = 0x0;           
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[itg][ir] = 0x0;
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[itg][ir] = 0x0;
      }
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
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      for(Int_t ic=0; ic<fkCE;ic++){
         fhCentrality[itg][ic] = 0x0;
         fhSignal[itg][ic] = 0x0; 
      
         for(Int_t i=0; i<fkTTbins;i++){
            fhCentralityTTH[itg][ic][i] = 0x0;
            fhCentralityTTJ[itg][ic][i] = 0x0;
            fhCentralityTTC[itg][ic][i] = 0x0;
      
            fhSignalTTH[itg][ic][i] = 0x0;
            fhSignalTTJ[itg][ic][i] = 0x0;
            fhSignalTTC[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
         fhSignalTTC_PartLevel[ic][i] = 0x0;
      } 
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhV0MAssymVsV0Mnorm[itg] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins;i++){
      for(Int_t ir=0; ir<kRho;ir++){
         fhRhoTTHinMBpart[i][ir]=0x0;
         fhRhoTTCinMBpart[i][ir]=0x0;
      } 
   }

   fFiducialCellCut = new AliEMCALRecoUtils();
 
   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1;
      for(Int_t ir =0; ir< kRho; ir++){
         fdeltapT[i][ir]  = 0.; 
         fdeltapT_PartLevel[i][ir]  = 0.; 
      }

      fIndexTTH_PartLevel[i] = -1;
      fIndexTTC_PartLevel[i] = -1;

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      fTTC_PartLevel[i].resize(0);
   }

   for(Int_t i=0; i<999; i++){
      frhoveckt[i] = 0.;
      frhovecms[i] = 0.;
   }

   fHelperEA = new PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA();
   fMeanV0M_PartLevel = fHelperEA->GetV0MPartLevel(); 
}

//________________________________________________________________________
AliAnalysisTaskEA::AliAnalysisTaskEA(const char *name): 
AliAnalysisTaskEmcalJet(name,kTRUE),  
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyDetLevelContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyJetDetLevelContainerName(""),
fMyClusterContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName(""),
fMyKTJetDetLevelContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fTrkContainerDetLevelEMB(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fJetContainerDetLevelEMB(0x0),
fClusterContainerDetLevel(0x0),
fKTJetContainerDetLevel(0x0),
fKTJetContainerPartLevel(0x0),
fKTJetContainerDetLevelEMB(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsEmcalTrig(0),
fIsHighMultTrig(0),
fCentralityV0A(-1),
fCentralityV0C(-1),
fCentralityV0M(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
fNTracklets(-1),
fMultV0A(0.),
fMultV0C(0.),
fMultV0M(0.),
fMultV0Mnorm(0.),
fAsymV0M(999),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
fAsymV0M_PartLevel(999),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0), 
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhV0MAssymVsV0Mnorm_PartLevel(0x0),
fhV0AvsV0C(0x0),
//fhV0MvsV0Mnorm(0x0),
fhV0AvsSPD(0x0),
fhV0CvsSPD(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelZero(0x0),
fhV0ARunByRunMB(0x0),
fhV0CRunByRunMB(0x0),
fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
fhTrackEtaInclEMB(0x0),
fhJetPtPartLevelVsJetPtDetLevelZero(0x0),
fhOneOverPtVsPhiNeg(0x0),
fhOneOverPtVsPhiPos(0x0),
fhSigmaPtOverPtVsPt(0x0),
fhDCAinXVsPt(0x0),
fhDCAinYVsPt(0x0),
fhDCAinXVsPtPhysPrimary(0x0),
fhDCAinYVsPtPhysPrimary(0x0),
fhDCAinXVsPtSecondary(0x0),
fhDCAinYVsPtSecondary(0x0),
fMinFractionShared(0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fnJetChTTBins(0),
fnClusterTTBins(0),
fMode(AliAnalysisTaskEA::kNormal),
fFiducialCellCut(0x0),
fHelperEA(0x0),
fMeanV0M(1.),
fMeanV0M_PartLevel(1.),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.)                            
{
   //Constructor
   for(Int_t ir=0; ir<kRho; ir++){
      fhRhoMBpart[ir] = 0x0;
      fhJetPtAreaV0norm_PartLevel[ir] = 0x0;
      fhJetPtPartLevelVsJetPtDetLevelCorr[ir] = 0x0;
      fhJetPtResolutionVsPtPartLevel[ir] = 0x0;
      fhJetPtPartLevelCorr[ir] = 0x0;
      fhFractionOfSecInJet[ir] = 0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;

      fhTrackMult[itg]=0x0;
      fhMeanTrackPt[itg]=0x0;

      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      fhTrackPtEtaPhiV0norm[itg] = 0x0;
      fhJetPtEtaPhiV0norm[itg] = 0x0;
      for(Int_t i=0; i<fkTTbins; i++){
         fhJetPtEtaPhiV0normTTH[itg][i] = 0x0;
      } 

      fhClusterPhiIncl[itg] = 0x0;
      fhClusterEtaIncl[itg] = 0x0;

      for(Int_t ir=0; ir<kRho; ir++){
         fhJetPtAreaV0norm[itg][ir] = 0x0;
         fhRho[itg][ir] = 0x0;
      }
    
      for(Int_t i=0; i<fkTTbins; i++){
         for(Int_t ir=0; ir<kRho; ir++){
            fhRhoTTH[itg][i][ir]=0x0;
            fhRhoTTC[itg][i][ir]=0x0;
         } 
         fhRhoTTJ[itg][i]=0x0;  
      }

      fhSharedJetFraction[itg] = 0x0;
      fhTrialsEMBtot[itg] = 0x0;
      fhXsectionEMBtot[itg] = 0x0;
      fhTrialsEMB[itg] = 0x0;
      fhXsectionEMB[itg] = 0x0;
      fhPtHardEMB[itg] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i]   = 0;
      fJetChTT[i]    = 0;
      fClusterTT[i]  = 0;


      fHadronTT_PartLevel[i]   = 0;
      fClusterTT_PartLevel[i]   = 0;

      //TT
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         fhMultTTH[itg][i] = 0x0;   
         fhMultTTJ[itg][i] = 0x0;  
         fhMultTTC[itg][i] = 0x0;  
         
         //fhTTH_CentV0M[itg][i]  = 0x0;
         fhTTH_V0Mnorm1[itg][i] = 0x0;
         fhTTH_3D_V0Mnorm1[itg][i] = 0x0;

         //fhTTC_CentV0M[itg][i]  = 0x0;
         fhTTC_V0Mnorm1[itg][i] = 0x0;

         fhV0MAssymVsV0MnormTTH[itg][i] = 0x0;
      }
 
      fhTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTH_3D_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTC_V0Mnorm1_PartLevel[i] = 0x0;

      fhV0MAssymVsV0MnormTTH_PartLevel[i] = 0x0;
         
      //RECOIL JET SPECTRA   
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         for(Int_t ir=0; ir<kRho; ir++){ 
            //fhRecoilJetPtTTH_CentV0M[itg][i][ir]  = 0x0;
            fhRecoilJetPtTTH_V0Mnorm1[itg][i][ir] = 0x0;
            
            fhRecoilJetPhiTTH_V0Mnorm1[itg][i][ir] = 0x0;
            fhRecoilJetTTH_V0Mnorm1[itg][i][ir]    = 0x0;
            
            //fhRecoilJetPtTTC_CentV0M[itg][i][ir]  = 0x0;
            fhRecoilJetPtTTC_V0Mnorm1[itg][i][ir] = 0x0;
         }
      }
         
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhRecoilJetPtTTH_V0Mnorm1_PartLevel[i][ir] = 0x0;
         fhRecoilJetPtTTC_V0Mnorm1_PartLevel[i][ir] = 0x0;
         
         fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[i][ir] = 0x0;
         fhRecoilJetTTH_V0Mnorm1_PartLevel[i][ir]    = 0x0;         
      }      
   
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         //fhDeltaPtTTH_RC_CentV0M[itg][i] = 0x0;  
         //fhDeltaPtTTC_RC_CentV0M[itg][i] = 0x0;
        
         for(Int_t ir=0; ir<kRho; ir++){ 
            fhDeltaPtTTH_RC_V0Mnorm1[itg][i][ir] = 0x0;  
            fhDeltaPtTTC_RC_V0Mnorm1[itg][i][ir] = 0x0;
         } 
      } 
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[i][ir] = 0x0;  
         fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[i][ir] = 0x0;
      }
 
      //remx 
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhJetPtPartLevelCorrTTHdl[i][ir] = 0x0;
         fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[i][ir] = 0x0;
      }
   
      //embedding
      for(Int_t itg=kMB; itg<=kGA; itg++){ 
         for(Int_t ir=0; ir<kRho; ir++){ 
            fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][i][ir] = 0x0;
            fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][i][ir] = 0x0;
            
            fhJetPtPartLevelCorrTTHdl_EMB[itg][i][ir] = 0x0;
            
            fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[itg][i][ir] = 0x0;
            fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[itg][i][ir]= 0x0;
         }
         fhJetPtPartLevelZeroTTHdl_EMB[itg][i] = 0x0;
      }  
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      for(Int_t ir=0; ir<kRho; ir++){ 
         fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir] = 0x0; 
         fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir] = 0x0; 
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      for(Int_t ir=0; ir<kRho; ir++){
         fhJetPtPartLevelCorr_EMB[itg][ir] = 0x0;
      }
      fhJetPtPartLevelZero_EMB[itg] = 0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      for(Int_t ir=0; ir<kRho; ir++){
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[itg][ir] = 0x0;
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[itg][ir] = 0x0;
      }
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
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){ 
      for(Int_t ic=0; ic<fkCE;ic++){
         fhCentrality[itg][ic] = 0x0;
         fhSignal[itg][ic] = 0x0; 
      
         for(Int_t i=0; i<fkTTbins;i++){
            fhCentralityTTH[itg][ic][i] = 0x0;
            fhCentralityTTJ[itg][ic][i] = 0x0;
            fhCentralityTTC[itg][ic][i] = 0x0;
      
            fhSignalTTH[itg][ic][i] = 0x0;
            fhSignalTTJ[itg][ic][i] = 0x0;
            fhSignalTTC[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0; 

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
         fhSignalTTC_PartLevel[ic][i] = 0x0;
      } 
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhV0MAssymVsV0Mnorm[itg] = 0x0;
   }


   for(Int_t i=0; i<fkTTbins;i++){
      for(Int_t ir=0; ir<kRho;ir++){

         fhRhoTTHinMBpart[i][ir]=0x0;
         fhRhoTTCinMBpart[i][ir]=0x0; 
      }
   }

   fFiducialCellCut = new AliEMCALRecoUtils();
 
   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      fIndexTTJ[i] = -1;
      for(Int_t ir=0; ir<kRho; ir++){
         fdeltapT[i][ir]  = 0.; 
         fdeltapT_PartLevel[i][ir]  = 0.; 
      }

      fIndexTTH_PartLevel[i] = -1;
      fIndexTTC_PartLevel[i] = -1;

      fTTC[i].resize(0);
      fTTH[i].resize(0);
      fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      fTTC_PartLevel[i].resize(0);
   }

   for(Int_t i=0; i<999; i++){
      frhoveckt[i] = 0.;
      frhovecms[i] = 0.;
   }

   fHelperEA = new PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA();
   fMeanV0M_PartLevel = fHelperEA->GetV0MPartLevel(); 

   DefineOutput(1, TList::Class());
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

//_____________________________________________________________________________________
AliAnalysisTaskEA*  AliAnalysisTaskEA::AddTaskEA(
  Int_t       mode,  
  const char* jetarrayname, 
  const char* jetarraynamePartMC, 
  const char* jetarraynameDetMC, 
  const char* trackarrayname, 
  const char* mcpariclearraynamePartMC, 
  const char* tracknameDetMC, 
  const char* clusterarrayname, 
  const char* ktjetarrayname, 
  const char* ktjetarraynamePartMC,
  const char* ktjetarraynameDetMC, 
  Double_t    jetRadius, 
  UInt_t      trigger, 
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
   Double_t jetRadiuskt   = 0.4;  //for all kt jets use fixed jet radius
   Double_t jetEtaRangekt = TMath::Abs(trackEtaWindow - jetRadiuskt);

   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskEA.cxx", "No analysis manager to connect to.");
      return NULL;
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString myContName("");
   myContName = Form("JetAnalysisR%02d_Acut%02d", TMath::Nint(jetRadius*10), TMath::Nint(acut*10));
   myContName.Append(suffix);
   if(mode == AliAnalysisTaskEA::kEmbedding)  myContName.Append("EMB");

   AliAnalysisTaskEA *task = new AliAnalysisTaskEA(myContName.Data());

   if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding){  //TO BE CHECKED FOR EMBEDDING 
      //for PYTHIA        
      task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer    *trackCont        = 0x0; // detector level track container (or tracks in  combined events when embedding )
   AliParticleContainer *trackContTrue    = 0x0; //mc particle container on  particle level for jets
   AliTrackContainer    *trackContDet     = 0x0; //mc particle container on  detector level for jets (for embedding)
   AliClusterContainer  *clusterCont      = 0x0; //detector level track container 

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks (or combined tracks if embedding) 
   trackCont->SetMinPt(0.15);
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);

   if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding){
      trackContTrue = task->AddMCParticleContainer(mcpariclearraynamePartMC); //particle level MC particles   
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(0.15);
      trackContTrue->SetEtaLimits(-5.1,5.1); //V0 eta range

      if(mode == AliAnalysisTaskEA::kEmbedding) trackContTrue->SetIsEmbedding(kTRUE);
   }

   if(mode == AliAnalysisTaskEA::kEmbedding){
      trackContDet = task->AddTrackContainer(tracknameDetMC);  //detector level pythia tracks when embedding
      trackContDet->SetMinPt(0.15);
      trackContDet->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
      trackContDet->SetIsEmbedding(kTRUE);
   }


   clusterCont = task->AddClusterContainer(clusterarrayname);  //detector level tracks (needs to be checked for embedding) 
   clusterCont->SetMinPt(0.3);
   clusterCont->SetExoticCut(1);
   clusterCont->SetClusTimeCut(0, emcaltofcut);

   //   clusterCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
 
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //AKT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrue   = 0x0; //AKT jet container with mc particle level jets pythia
   AliJetContainer *jetContDet    = 0x0; //AKT jet container used when embedding with mc jets at detector level pythia

   AliJetContainer *jetContRecKT  = 0x0; //KT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrueKT = 0x0; //KT jet container with mc particle level jets pythia
   AliJetContainer *jetContDetKT  = 0x0; //KT jet container used when embedding with mc jets at detector level pythia



   //AKT DETECTOR LEVEL JET    (or combined event jet container when embedding)
   jetContRec   = task->AddJetContainer(jetarrayname,"TPC",jetRadius);

   if(jetContRec) { 
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);
      jetContRec->SetMinPt(0.150);
      jetContRec->SetMaxTrackPt(100.);
      jetContRec->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);

    }

   //KT DETECTOR LEVEL JET    (or combined event jet container when embedding)
   jetContRecKT   = task->AddJetContainer(ktjetarrayname,"TPC",jetRadiuskt);

   if(jetContRecKT) { 
      jetContRecKT->ConnectParticleContainer(trackCont);
      //jetContRecKT->SetPercAreaCut(acut);
      jetContRecKT->SetMinPt(0.);
      jetContRecKT->SetMaxTrackPt(100.);
      jetContRecKT->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);

    }


    if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding){
      //AKT JETS PARTICLE LEVEL
      jetContTrue = task->AddJetContainer(jetarraynamePartMC,"TPC",jetRadius);

      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(100.);
         jetContTrue->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT JETS PARTICLE LEVEL
      jetContTrueKT = task->AddJetContainer(ktjetarraynamePartMC,"TPC",jetRadiuskt);

      if(jetContTrueKT){
         jetContTrueKT->ConnectParticleContainer(trackContTrue);
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetMaxTrackPt(100.);
         jetContTrueKT->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
      }
   }

   if(mode == AliAnalysisTaskEA::kEmbedding){
      //AKT DETECTOR LEVEL JET    (or combined event jet container when embedding)
      jetContDet = task->AddJetContainer(jetarraynameDetMC,"TPC",jetRadius);

      if(jetContDet) { 
         jetContDet->SetPercAreaCut(acut);
         jetContDet->SetMinPt(0.150);
         jetContDet->SetMaxTrackPt(100.);
         jetContDet->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContDet->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT DETECTOR LEVEL JET    (or combined event jet container when embedding)
      jetContDetKT = task->AddJetContainer(ktjetarraynameDetMC,"TPC",jetRadiuskt);

      if(jetContDetKT) { 
         //jetContDetKT->SetPercAreaCut(acut);
         jetContDetKT->SetMinPt(0.0);
         jetContDetKT->SetMaxTrackPt(100.);
         jetContDetKT->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContDetKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
      }
   }

   // #### Task configuration 
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetAcceptanceWindows(trackEtaWindow);
   task->SelectCollisionCandidates(trigger);
   task->SetMode(mode);
 
   task->SetTrackContainerName(trackarrayname);
   task->SetMCParticleContainerName(mcpariclearraynamePartMC); 
   task->SetMCDetLevelContainerName(tracknameDetMC);
   task->SetClusterContainerName(clusterarrayname);

   task->SetJetContainerName(jetarrayname);
   task->SetMCPartJetContainerName(jetarraynamePartMC); 
   task->SetMCDetJetContainerName(jetarraynameDetMC);
   task->SetKTJetContainerName(ktjetarrayname);
   task->SetKTMCPartJetContainerName(ktjetarraynamePartMC);
   task->SetKTMCDetJetContainerName(ktjetarraynameDetMC);

   task->SetJetRadius(jetRadius); 
   task->SetJetAcut(acut);

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

  if(fMode == AliAnalysisTaskEA::kMC) return kFALSE; //MC

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

  bool passedTrigger = kFALSE;

  if(fMode == AliAnalysisTaskEA::kMC){ 
     //mc simulation emulate V0 coincidence trigger
     AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
     if(vzeroAOD){
        if(vzeroAOD->GetMTotV0A() > 0 && vzeroAOD->GetMTotV0C() > 0)  passedTrigger = kTRUE; 
     }
  }else{
     //real data and embedding  take trigger decission from data
     UInt_t triggerMask = fInputHandler->IsEventSelected();
     if(triggerMask & AliVEvent::kINT7){
        passedTrigger = kTRUE;
     }
  }

  return passedTrigger;
}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEA::PassedHighMultTrigger(){
   //high multiplicity V0M trigger

   if(fMode == AliAnalysisTaskEA::kMC) return kFALSE; //MC

   bool passedTrigger = kFALSE;
   UInt_t triggerMask = fInputHandler->IsEventSelected();
   if(triggerMask & AliVEvent::kHighMultV0){
      passedTrigger = kTRUE;
   }

   return passedTrigger;
}
//_____________________________________________________________________________________
void AliAnalysisTaskEA::GetMyRho(AliJetContainer* ktjets, Double_t &rhokt, Double_t &rhocms){

   // Get rho from event using CMS approach
   rhokt  = 0.;
   rhocms = 0.;

   Int_t nJetAcckt = 0;
   Int_t nJetAccms = 0;
   Double_t areaPhysJets = 0.0;
   Double_t areaAllJets  = 0.0;


   Double_t ptLJ=-1;
   Double_t ptSJ=-1;
   AliEmcalJet*  jetLJ = 0x0;
   AliEmcalJet*  jetSJ = 0x0;
   AliEmcalJet*  jet   = 0x0;

   //Exclude 2 leading jets 
   for(auto jetIterator : ktjets->accepted_momentum() ){
                   // trackIterator is a std::map of AliTLorentzVector and AliVTrack
       jet = jetIterator.second;  // Get the pointer to jet object
       if(!jet)  continue; 
    
       if(jet->Pt() > ptLJ){
          ptSJ  = ptLJ;
          jetSJ = jetLJ;

          ptLJ  = jet->Pt(); 
          jetLJ = jet; 
       }else if(jet->Pt() > ptSJ){
          ptSJ  = jet->Pt();
          jetSJ = jet;
       }
   }


   for(auto jetIterator : ktjets->accepted_momentum() ){
                   // trackIterator is a std::map of AliTLorentzVector and AliVTrack
       jet = jetIterator.second;  // Get the pointer to jet object
       if(!jet)  continue; 

       if(jet==jetLJ) continue; //skip two leading kT jets 
       if(jet==jetSJ) continue; 
   
       //standard area based approach
       frhoveckt[nJetAcckt]  = jet->Pt()/jet->Area();
       nJetAcckt++;

       //cms modification
       areaAllJets += jet->Area();

       if(jet->Pt() > 0.1){
          areaPhysJets     += jet->Area();
          frhovecms[nJetAccms]  = jet->Pt()/jet->Area();
          nJetAccms++;
       }
   }

   if(nJetAcckt>0){
      rhokt = TMath::Median(nJetAcckt, frhoveckt);
   }

   if(nJetAccms>0){
      rhocms = TMath::Median(nJetAccms, frhovecms)*(areaPhysJets/areaAllJets);
   }

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

   TLorentzVector myTT;
   Int_t idx;
   TString name;
   Double_t tmparr[4]; 
   Double_t tmparr3[3]; 
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
   //Double_t weight = 1.;
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

   //_________________________________________________________________
   Double_t rho[kRho]    = {0.,0.}; 
   Double_t rhoMC[kRho]  = {0.,0.};  
   Double_t rhoEMB[kRho] = {0.,0.};  
   Double_t dphi     = 999.;
   
   Double_t jetPtCorrPart = 0.; //particle level jet pt corrected for rho
   Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
   Double_t jetPtCorrDetShift1  = 0.;  //detector level jet pt corrected for rho
   Double_t jetPtCorrDetShift2  = 0.;  //detector level jet pt corrected for rho
   Double_t sharedFraction = 0.; //shared fraction between detector level and combined level

   Double_t ptLJ=-1, etaLJ=999, phiLJ=0; //leading jet
   Double_t ptSJ=-1, etaSJ=999, phiSJ=0; //subleading jet
   Double_t ptLJmc=-1, etaLJmc=999, phiLJmc=0; //leading jet
   Double_t ptSJmc=-1, etaSJmc=999, phiSJmc=0; //subleading jet

   //_________________________________________________________________
   //                    JET+TRACK CONTAINERS
  
   AliEmcalJet  *jet = NULL;        //jet pointer real jet 
   AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
   AliEmcalJet  *jetDetMC  = NULL;  //jet pointed detector level MC jet
   AliVParticle *track = NULL; //jet constituent
   AliVParticle *mcParticle = NULL; //mc particle
   //_________________________________________________________________


   Bool_t trigflag[] = {fIsMinBiasTrig, fIsHighMultTrig, fIsEmcalTrig};

   //_________________________________________________________________
   
   
   // END EVENT SELECTION
   //_________________________________________________________________
   // DECIDE WHETHER TO FILL SIGNAL TT OR REFERENCE TT  DEPENDING ON RANDOM  NUMBER  
   
   fFillSigTT = kTRUE;  
   if( fRandom->Integer(100) < 5) fFillSigTT = kFALSE; 
   
   //_________________________________________________________________
   //                EVENT PROPERTIES   
   
   
   fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(fMultSelection){
   
      fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
   }else{
      fCentralityV0A = -1; 
      fCentralityV0C = -1;
      fCentralityV0M = -1;
   }
   
   const AliVVertex *vertex = InputEvent()->GetPrimaryVertexSPD();
   if(vertex){ 
      fxVertex = vertex->GetX();
      fyVertex = vertex->GetY();
      fzVertex = vertex->GetZ();
   }else{
      fxVertex = 9999.;
      fyVertex = 9999.;
      fzVertex = 9999.;
   }
   
   const AliVMultiplicity *mult = InputEvent()->GetMultiplicity();
   if(mult){
      fNTracklets = mult->GetNumberOfTracklets();
   }else{
      fNTracklets = -9999;
   }
   
   Int_t runnumber  = InputEvent()->GetRunNumber();  
   
   AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
   if(vzeroAOD){
      Double_t meanV0A, meanV0C;
   
      fMultV0A = vzeroAOD->GetMTotV0A();
      fMultV0C = vzeroAOD->GetMTotV0C();
      fMultV0M = fMultV0A + fMultV0C;
   
      if(fMode != AliAnalysisTaskEA::kMC){
         fMeanV0M = fHelperEA->GetV0M(runnumber);
         meanV0A  = fHelperEA->GetV0A(runnumber);
         meanV0C  = fHelperEA->GetV0C(runnumber);
      }else{
         fMeanV0M = fHelperEA->GetV0MDetLevel();
         meanV0A  = fHelperEA->GetV0ADetLevel();
         meanV0C  = fHelperEA->GetV0CDetLevel();
      }      
   
      //V0 estimators of event activity normalized per minimum bias activity 
      fMultV0Mnorm = fMultV0M/fMeanV0M;
   
      fAsymV0M = 999;
   
      if(meanV0A>0 && meanV0C>0){
         Double_t normV0A = fMultV0A/meanV0A;
         Double_t normV0C = fMultV0C/meanV0C;
         if((normV0A + normV0C) > 0){
            fAsymV0M = (normV0A - normV0C)/(normV0A + normV0C);
         } 
      }
   }else{
      fMultV0A = -1; 
      fMultV0C = -1; 
      fMultV0M = -1; 
      fMultV0Mnorm = -1; 
      fAsymV0M = 999;
   }
   
  
   
   
   //_________________________________________________________
   //READ  TRACK AND JET CONTAINERS
   //Container operations   http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerIterateTechniques
   
   //fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(fMyTrackContainerName.Data())); //track container detector-level   real data only
   fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(0)); //track container detector-level   real data only
   fJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyJetContainerName.Data()));       //detector-level AKT jets    real data or hybrid event
   fKTJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetContainerName.Data()));     //detector-level KT jets    real data or hybrid event
   
   GetMyRho(fKTJetContainerDetLevel, (Double_t&) rho[krhokt], (Double_t&) rho[krhocms]); //estimated backround pt density
   
   if( fMode != AliAnalysisTaskEA::kNormal){  //particle level particles and jets  for  MC and embedding
      //fParticleContainerPartLevel = GetParticleContainer(fMyParticleContainerName.Data()); //pythia particle level particles 
      fParticleContainerPartLevel = GetParticleContainer(1); //pythia particle level particles 
      fJetContainerPartLevel      = static_cast<AliJetContainer*> (GetJetContainer(fMyJetParticleContainerName.Data()));   //pythia particle level AKT jets
      fKTJetContainerPartLevel    = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetParticleContainerName.Data()));   //pythia particle level KT jets
   
      GetMyRho(fKTJetContainerPartLevel, (Double_t&) rhoMC[krhokt], (Double_t&) rhoMC[krhocms]); //estimated backround pt density
   }
   
   if( fMode == AliAnalysisTaskEA::kEmbedding){ //Detector level pythia  for  embedding
   
      //fTrkContainerDetLevelEMB = static_cast<AliTrackContainer*> (GetTrackContainer(fMyDetLevelContainerName.Data())); //pythia detector level tracks 
      fTrkContainerDetLevelEMB   = static_cast<AliTrackContainer*> (GetTrackContainer(2)); //pythia detector level tracks 
      fJetContainerDetLevelEMB   = static_cast<AliJetContainer*> (GetJetContainer(fMyJetDetLevelContainerName.Data()));  //pythia detector level AKT jets 
      fKTJetContainerDetLevelEMB = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetDetLevelContainerName.Data()));  //pythia detector level KT jets 
   
      GetMyRho(fKTJetContainerDetLevelEMB, (Double_t&) rhoEMB[krhokt], (Double_t&) rhoEMB[krhocms]); //estimated backround pt density
   }
   
   
   //________________________________________________________
   //     Find the leading and subleading jets  for estimates  of Delta pt
   
   
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
   if(fMode != AliAnalysisTaskEA::kNormal){
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
  
      for(Int_t ir=0; ir< kRho; ir++){ 
         fdeltapT[i][ir]  = 0.; 
         fdeltapT_PartLevel[i][ir]  = 0.; 
      }
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
   
   for(Int_t i=0; i<fnJetChTTBins; i++){
      fJetChTT[i] = 0;
   }
   
   Double_t xyz[50];
   Double_t pxpypz[50];
   Double_t cv[21];
   Double_t shift1, shift2;
   
   Int_t label, labelMC;                  
   Bool_t labelfound = 0;
   AliAODMCParticle* particleMC = NULL;    
   AliAODTrack *trackAOD=NULL ;

   Int_t    trackMult  = 0;
   Double_t sumTrackPt = 0.;
   
   if(fIsEmcalTrig || fIsMinBiasTrig || fIsHighMultTrig){  //real data + mc det level + embedding
      //___________________________________________
      //    INCLUSIVE EVENTS (WITHOUT TT REQUIREMENT)
   
  
      for(Int_t itg=kMB; itg<=kGA; itg++){    //@@@
         if(!trigflag[itg]) continue; 
         //events without TT requirement
         for(Int_t ir = 0; ir < kRho; ir++){
            fhRho[itg][ir]->Fill(rho[ir]);
         }
   
          fhCentrality[itg][fkV0A]->Fill(fCentralityV0A, fMultV0A); 
          fhCentrality[itg][fkV0C]->Fill(fCentralityV0C, fMultV0C); 
          fhCentrality[itg][fkV0M]->Fill(fCentralityV0M, fMultV0M); 
          fhCentrality[itg][fkV0Mnorm1]->Fill(fCentralityV0M, fMultV0Mnorm); 
   
          fhSignal[itg][fkV0A]->Fill(fMultV0A);
          fhSignal[itg][fkV0C]->Fill(fMultV0C);
          fhSignal[itg][fkV0M]->Fill(fMultV0M);
          fhSignal[itg][fkV0Mnorm1]->Fill(fMultV0Mnorm);
   
          fhV0MAssymVsV0Mnorm[itg]->Fill(fMultV0Mnorm, fAsymV0M);
      }
   
   
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
   
      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         if(!trigflag[itg]) continue; 
    
         fhTrackMult[itg]->Fill(fMultV0Mnorm, trackMult); 
         fhMeanTrackPt[itg]->Fill(fMultV0Mnorm, sumTrackPt);
      }
   
      if(fIsMinBiasTrig){ // run for all MC events   and for real data with min bias trigger
         fhVertex[0]->Fill(fxVertex);
         fhVertex[1]->Fill(fyVertex);
         fhVertex[2]->Fill(fzVertex);
   
   
         name = Form("%d", runnumber);      
         fhV0ARunByRunMB->Fill(name.Data(), fMultV0A, 1.0);
         fhV0CRunByRunMB->Fill(name.Data(), fMultV0C, 1.0);
         fhV0MRunByRunMB->Fill(name.Data(), fMultV0M, 1.0);
         fhV0MnormRunByRunMB->Fill(name.Data(), fMultV0Mnorm, 1.0);
   
         fhV0AvsV0C->Fill(fMultV0C, fMultV0A);
         fhV0AvsSPD->Fill(fNTracklets, fMultV0A);
         fhV0CvsSPD->Fill(fNTracklets, fMultV0C);
      }
   
   
      //_________________________________________________________
      //LOOP OVER TRACKS DETECTOR LEVEL
    
      for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         track = trackIterator.second;  // Get the full track
         if(!track) continue;
   
         if(IsTrackInAcceptance(track, kDetLevel)){  
   
            for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
               if(!trigflag[itg]) continue; 
               fhTrackPhiIncl[itg]->Fill(track->Pt(), track->Phi());
               fhTrackEtaIncl[itg]->Fill(track->Pt(), track->Eta());
   
               tmparr[0] = track->Pt();
               tmparr[1] = track->Eta();
               tmparr[2] = track->Phi();
               tmparr[3] = fMultV0Mnorm;
               fhTrackPtEtaPhiV0norm[itg]->Fill(tmparr);
            }
         }
      }
   
   
   
      //___________________________________________
      //          EMBEDDED EVENTS     
      if(fMode == kEmbedding){
   
         if(!fIsMinBiasTrig && !fIsHighMultTrig) return kTRUE; //if this is not MB or HM  skip the rest
   
         const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
         double ptHardBin = embeddingHelper->GetPtHardBin();
   
         for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
            if(!trigflag[itg]) continue; 
            fhTrialsEMBtot[itg]->Fill(0.5, embeddingHelper->GetPythiaTrials());
            fhXsectionEMBtot[itg]->Fill(0.5, embeddingHelper->GetPythiaXSection());
        
            fhTrialsEMB[itg]->Fill( ptHardBin, embeddingHelper->GetPythiaTrials());
            fhXsectionEMB[itg]->Fill( ptHardBin, embeddingHelper->GetPythiaXSection());
            fhPtHardEMB[itg]->Fill( embeddingHelper->GetPythiaPtHard());
         }
   
         //Find TT among the PYTHIA Detector level tracks
         // This TT will be used for pythia detector level events as well as for the combined event      
    
         if(fTrkContainerDetLevelEMB){
   
            //detector level pythia mc particles
            for(auto mcDetIterator : fTrkContainerDetLevelEMB->accepted_momentum() ){
               track = mcDetIterator.second;  // Get the pointer to mc particle object
               if(!track)  continue; 
   
               if(IsTrackInAcceptance(track, kDetLevel)){
                  fhTrackEtaInclEMB->Fill(track->Pt(), track->Eta());
   
                  for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                     if(fHadronTTLowPt[itt] < track->Pt() && track->Pt() < fHadronTTHighPt[itt]){
                        myTT.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),0.); 
                        fTTH[itt].push_back(myTT);
                        fHadronTT[itt]++;   // there was a high pt 
                     }
                  }
               }        
            }
   
            //chose trigger hadron TT which will be common for  Detector level pythia and the combined event
            for(Int_t itt=0; itt<fnHadronTTBins; itt++){
               if(fHadronTT[itt]>0){
                  fIndexTTH[itt] = fRandom->Integer(fHadronTT[itt]);
               }
            }
      
            for(Int_t itt=0; itt<fnHadronTTBins; itt++){
               idx = fIndexTTH[itt];//hadron trigger
               if(idx<0) continue;
   
               if(fFillSigTT && itt==0) continue;  // Do not fill reference 
               if(!fFillSigTT && itt>0) continue;  // Do not fill signal 
   
               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!trigflag[itg]) continue; 
                  fhTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
               }
   
               //recoil jets pythia detector level event 
               for(auto jetIterator : fJetContainerDetLevelEMB->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
             
                  dphi = TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi());  

                  for(Int_t ir=0; ir<kRho; ir++){ 
                     jetPtCorrDet = jet->Pt() - rhoEMB[ir]*jet->Area();
                     
                     for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                        if(!trigflag[itg]) continue; 
                        fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi); 
                     }
                  } 
               } 
   
               //recoil jets in the combined event
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
              
                  dphi = TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi()); 
                  for(Int_t ir=0; ir<kRho; ir++){ 
                     jetPtCorrDet = jet->Pt() - rho[ir]*jet->Area();
                     
                     for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                        if(!trigflag[itg]) continue; 
                        fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
                     
                     }
                  }
                  //fill similar disribution but just for jets which have pythia partner
                  jetDetMC =  jet->ClosestJet(); //This is the closes pythia Detector level jet
                  if(jetDetMC){ 

                     dphi = TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[itt][idx].Phi()); 
                     for(Int_t ir=0; ir<kRho; ir++){ 
                        jetPtCorrDet  =  jetDetMC->Pt() - jetDetMC->Area()*rhoEMB[ir]; 
                     
                        for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                           if(!trigflag[itg]) continue; 
                           fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi); 
                        }
                     } 
                  }
               }//jet loop 
            }//TT loop
         }//EMB track container
   
         //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         //                RESPONSE MATRIX FROM EMBEDDED EVENTS
         //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
         //Response matrix normalization.  The matrix will be constructed using 
         //- inclusive generator level jets in acceptance
         //- recoil jets  (that recould from detector level pythia TT)
         if(fJetContainerPartLevel){
            for(auto jetPartIterator : fJetContainerPartLevel->accepted_momentum() ){
               jetPartMC = jetPartIterator.second;  // Get the pointer to mc particle object
               if(!jetPartMC)  continue; 
   
               for(Int_t ir = 0; ir<kRho; ir++){
                  jetPtCorrPart = jetPartMC->Pt() - jetPartMC->Area()*rhoMC[ir];
                  
                  for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                     if(!trigflag[itg]) continue; //inclusive jet spectrum for ReMx normalization 
                     fhJetPtPartLevelCorr_EMB[itg][ir]->Fill(jetPtCorrPart);
                     if(ir==0) fhJetPtPartLevelZero_EMB[itg]->Fill(jetPartMC->Pt());
                  }
                  
                  for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                     idx = fIndexTTH[itt];//hadron trigger
                     if(idx<0) continue;
                  
                     dphi = TVector2::Phi_mpi_pi(jetPartMC->Phi()-fTTH[itt][idx].Phi()); 
                  
                     if(TMath::Abs(dphi) > fPhiCut){    //fill with recoil jets only 
                        for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                           if(!trigflag[itg]) continue; //recoil jet spectrum for ReMx normalization 
                           fhJetPtPartLevelCorrTTHdl_EMB[itg][itt][ir]->Fill(jetPtCorrPart);
                           if(ir==0) fhJetPtPartLevelZeroTTHdl_EMB[itg][itt]->Fill(jetPartMC->Pt());
                        }
                     }
                  }
               }
            }
         }
   
         //FILL 2D RESPONSE MATRIX Find closest particle level and detector level jets  and detector level  combined level
         if(fJetContainerDetLevel){
   
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               jet = jetIterator.second;  // jet on combined level 
               if(!jet)  continue; 
    
               //find closest pythia detector level jet
               jetDetMC =  jet->ClosestJet();
               if(!jetDetMC) continue;
   
               sharedFraction = fJetContainerDetLevel->GetFractionSharedPt(jet); //Check shared momentum fraction
   
               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!trigflag[itg]) continue; //recoil jet spectrum for ReMx normalization 
                  fhSharedJetFraction[itg]->Fill(jetDetMC->Pt(), sharedFraction);
               } 
   
               if(sharedFraction < fMinFractionShared) continue; 
   
               jetPartMC = jetDetMC->ClosestJet(); //This is the closes pythia particle level jet to the pythia detector level jet
               if(!jetPartMC) continue; 
               if(jetPartMC->Pt()<1e-3) continue; //prevents matching with a ghost
   
               for(Int_t ir = 0; ir<kRho; ir++){
                  jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*rhoMC[ir]; 
                  jetPtCorrDet  =  jet->Pt() - jet->Area()*rho[ir]; 
                  
                  for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                     if(!trigflag[itg]) continue; //recoil jet spectrum for ReMx normalization 
                     fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir]->Fill(jetPtCorrDet, jetPtCorrPart); //response matrix
                     fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir]->Fill(jetPtCorrDet, jetPartMC->Pt()); //response matrix
                  }
                  
                  for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                     idx = fIndexTTH[itt];//hadron trigger
                     if(idx<0) continue;
                  
                     dphi = TVector2::Phi_mpi_pi(jetPartMC->Phi()-fTTH[itt][idx].Phi()); 
                  
                     if(TMath::Abs(dphi) > fPhiCut){    //fill with recoil jets only
                        for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                           if(!trigflag[itg]) continue; //recoil jet spectrum for ReMx normalization 
                           fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[itg][itt][ir]->Fill(jetPtCorrDet, jetPtCorrPart);
                           fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[itg][itt][ir]->Fill(jetPtCorrDet, jetPartMC->Pt());
                        }
                     }
                  }
               }
            }
         }
      
         return kTRUE;
      }
   
   
      //____________ INCLUSIVE EVENTS _______________________________
      if(fIsMinBiasTrig){
   
         for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
            // trackIterator is a std::map of AliTLorentzVector and AliVTrack
            track = trackIterator.second;  // Get the full track
            if(!track) continue;
    
            if(IsTrackInAcceptance(track, kDetLevel)){  
    
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
                 
                  //SINGLE TRACK EFFICIENCY AND CONTAMINATION
                  if(fMode == AliAnalysisTaskEA::kMC){
                     label = TMath::Abs(trackAOD->GetLabel()); 
                     
                     particleMC = NULL;
                     labelfound = 0;
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
         }
      }
   
      //______________________________________________________________
      //                       TTH  ANALYSIS 
      //______________________________________________________________
   
      if(fIsMinBiasTrig || fIsHighMultTrig){
         //LOOP SEARCH FOR HIGH PT HADRON TRIGGER IN INCLUSIVE EVENTS 
         for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
            // trackIterator is a std::map of AliTLorentzVector and AliVTrack
            track = trackIterator.second;  // Get the full track
            if(!track) continue;
      
            if(IsTrackInAcceptance(track, kDetLevel)){  
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
            for(Int_t ir=0; ir<kRho; ir++){
               fdeltapT[itt][ir] = 0.;
            }
            
            if(fHadronTT[itt]>0){
               fIndexTTH[itt] = fRandom->Integer(fHadronTT[itt]);
               idx = fIndexTTH[itt];
            
               for(Int_t ir=0; ir<kRho; ir++){
                  fdeltapT[itt][ir] = GetDeltaPt(fTTH[itt][idx].Phi(), fTTH[itt][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, rho[ir], kDetLevel);
               }
            }
         }
      
         for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
            if(!trigflag[itg]) continue; //check which trigger fired
            for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      
               fhMultTTH[itg][itt]->Fill(fHadronTT[itt]); 
      
               if(!fHadronTT[itt]) continue; //check whether there was hadron TT
      
               for(Int_t ir=0; ir<kRho; ir++){
                  fhRhoTTH[itg][itt][ir]->Fill(rho[ir]);
               } 
            }
         }
      
      
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(!fHadronTT[itt]) continue; //analyze events with hadron TT only 
   
            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!trigflag[itg]) continue; //check which trigger fired
       
               fhCentralityTTH[itg][fkV0A][itt]->Fill(fCentralityV0A, fMultV0A); 
               fhCentralityTTH[itg][fkV0C][itt]->Fill(fCentralityV0C, fMultV0C);
               fhCentralityTTH[itg][fkV0M][itt]->Fill(fCentralityV0M, fMultV0M); 
               fhCentralityTTH[itg][fkV0Mnorm1][itt]->Fill(fCentralityV0M, fMultV0Mnorm); 
       
               fhSignalTTH[itg][fkV0A][itt]->Fill(fMultV0A);
               fhSignalTTH[itg][fkV0C][itt]->Fill(fMultV0C);
               fhSignalTTH[itg][fkV0M][itt]->Fill(fMultV0M);
               fhSignalTTH[itg][fkV0Mnorm1][itt]->Fill(fMultV0Mnorm);
   
               fhV0MAssymVsV0MnormTTH[itg][itt]->Fill(fMultV0Mnorm, fAsymV0M);
            }
      
            if(fIsMinBiasTrig){
               fhV0AvsV0CTTH[itt]->Fill(fMultV0C, fMultV0A);
            }
   
            //pick up TTH hadron accoding to the index
            idx = fIndexTTH[itt];
            if(idx>-1){
   
               for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                  if(!trigflag[itg]) continue; //check which trigger fired
      
                  //fhDeltaPtTTH_RC_CentV0M[itg][itt]->Fill(fCentralityV0M, fdeltapT[itt]);
                  for(Int_t ir=0; ir<kRho; ir++){  
                     fhDeltaPtTTH_RC_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, fdeltapT[itt][ir]);        
                  } 
               }
   
               if(fFillSigTT  && itt==0) continue;  // Do not fill reference 
               if(!fFillSigTT && itt>0)  continue;  // Do not fill signal 
    
               for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                  if(!trigflag[itg]) continue; //check which trigger fired
      
                  //fhTTH_CentV0M[itg][itt]->Fill(fCentralityV0M, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0M centrality
                  fhTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm,  fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
                  fhTTH_3D_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, fAsymV0M, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
               }
      
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
      
                  dphi = TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi()); 

                  for(Int_t ir=0; ir<kRho; ir++){  
                     jetPtCorrDet = jet->Pt() - rho[ir]*jet->Area();
                    
                     for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                        if(!trigflag[itg]) continue; //check which trigger fired
                        fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
                    
                        tmparr[0] = fMultV0Mnorm;
                        tmparr[1] = fAsymV0M;
                        tmparr[2] = jetPtCorrDet;
                        tmparr[3] = TMath::Abs(dphi);
                        fhRecoilJetTTH_V0Mnorm1[itg][itt][ir]->Fill(tmparr); 

                        if(ir==0){
                           tmparr[0] = jet->Pt();
                           tmparr[1] = jet->Eta();
                           tmparr[2] = jet->Phi();
                           tmparr[3] = fMultV0Mnorm;
                           fhJetPtEtaPhiV0normTTH[itg][itt]->Fill(tmparr); 
                        }
                     } 
                    
                     if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > fPhiCut){     //select recoil jet
                        for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                           if(!trigflag[itg]) continue; //check which trigger fired
                           //fhRecoilJetPtTTH_CentV0M[itg][itt]->Fill(fCentralityV0M, jetPtCorrDet);
                           fhRecoilJetPtTTH_V0Mnorm1[itg][itt][ir]->Fill(fMultV0Mnorm, jetPtCorrDet);
                    
                           if(itt==0){ //reference with shifted rho
                              if(TMath::Abs( fJetAcut) < 1e-4){  //no area cut                            

                                 if(TMath::Abs(fJetR -0.4) < 1e-3){ //R=0.4 jets 

                                    if(itg==kMB){
                                       if(ir==krhokt){
                                          shift1 =  0.000896396;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.00603773;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = -0.00500762;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.0120839;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }else{
                                       if(ir==krhokt){
                                          shift1 = 0.0179261;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = 0.0322556;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = 0.0117027;  // pT shift of 6,7 w.r.t. 12 20            
                                          shift2 = 0.0207614;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }
                                 }else{  //R=0.2 jets
                                    if(itg==kMB){
                                       if(ir==krhokt){
                                          shift1 = -0.0123206;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.0217714;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = -0.0170195;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.0292258;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }else{
                                       if(ir==krhokt){
                                          shift1 = 0.00595735;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = 0.010721;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = 0.00394271;  // pT shift of 6,7 w.r.t. 12 20            
                                          shift2 = 0.00688605;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }
                                 }
                              }else{ //jet area cut
                                 if(TMath::Abs(fJetR -0.4) < 1e-3){ //R=0.4 jets 

                                    if(itg==kMB){
                                       if(ir==krhokt){
                                          shift1 = -0.0021434;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = 0.000527553;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = -0.00283202;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.00473093;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }else{
                                       if(ir==krhokt){
                                          shift1 = 0.0223449;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = 0.0386102;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = 0.0136156;  // pT shift of 6,7 w.r.t. 12 20            
                                          shift2 = 0.0222889;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }
                                 }else{  //R=0.2 jets
                                    if(itg==kMB){
                                       if(ir==krhokt){
                                          shift1 = -0.00873796;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.0166677;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = -0.0146481;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = -0.0272;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }else{
                                       if(ir==krhokt){
                                          shift1 = 0.00578544;  // pT shift of 6,7 w.r.t. 12 20             
                                          shift2 = 0.0110027;  // pT shift of 6,7 w.r.t. 20 30           
                                       }else{
                                          shift1 = 0.00355317;  // pT shift of 6,7 w.r.t. 12 20            
                                          shift2 = 0.0070802;  // pT shift of 6,7 w.r.t. 20 30           
                                       }
                                    }
                                 }
                              }
                              jetPtCorrDetShift1 = jetPtCorrDet - shift1;
                              jetPtCorrDetShift2 = jetPtCorrDet - shift2;
                    
                              fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[itg][ir]->Fill(fMultV0Mnorm, jetPtCorrDetShift1);
                              fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[itg][ir]->Fill(fMultV0Mnorm, jetPtCorrDetShift2);
                           }
                        }
                     }
                  } 
               }
            }
         }
      }
   
   
      //_________________________________________________________
      //              EMCAL CLUSTERS   TTC
      //_________________________________________________________
      if(fMyClusterContainerName.Data()){
         fClusterContainerDetLevel =  static_cast<AliClusterContainer*> ( GetClusterContainer(fMyClusterContainerName.Data()));
     
         for(auto cluster: fClusterContainerDetLevel->accepted()){
            fClusterContainerDetLevel->GetMomentum(ph, cluster);
   
            if(!FinalClusterCuts(cluster)) continue;
   
            for(Int_t itg = kMB; itg<=kGA; itg++){
               if(!trigflag[itg])  continue;
               fhClusterPhiIncl[itg]->Fill(ph.Pt(), ph.Phi());
               fhClusterEtaIncl[itg]->Fill(ph.Pt(), ph.Eta());
            }
   
   
            for(Int_t igg=0; igg<fnClusterTTBins; igg++){ // seatch for TTC candidates
               if(fClusterTTLowPt[igg] < ph.Pt() && ph.Pt() < fClusterTTHighPt[igg]){
                  myTT.SetPtEtaPhiM(ph.Pt(),ph.Eta(),ph.Phi(),0.); 
                  fTTC[igg].push_back(myTT);
                  fClusterTT[igg]++;   // there was a high pt emcal cluster 
               } 
            }
         }
   
         //chose trigger emcal cluster TTC
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            for(Int_t ir =0; ir< kRho; ir++){
               fdeltapT[igg][ir] = 0; 
            }

            if(fClusterTT[igg]>0){
               fIndexTTC[igg] = fRandom->Integer(fClusterTT[igg]);
               idx = fIndexTTC[igg];// gamma trigger
               for(Int_t ir =0; ir< kRho; ir++){
                  fdeltapT[igg][ir] = GetDeltaPt(fTTC[igg][idx].Phi(), fTTC[igg][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, rho[ir], kDetLevel);
               }
            }
         }
        
         for(Int_t itg=kMB; itg<=kGA; itg++){ //@@@
            if(!trigflag[itg]) continue; //check which trigger fired
            for(Int_t igg=0; igg<fnClusterTTBins; igg++){
               fhMultTTC[itg][igg]->Fill(fClusterTT[igg]); 
        
               if(!fClusterTT[igg]) continue;  //check whether there was TT
        
               for(Int_t ir =0; ir< kRho; ir++){
                  fhRhoTTC[itg][igg][ir]->Fill(rho[ir]);
               }
            }
         }
   
         //  analysis of      TTC   bias events
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
             
            if(!fClusterTT[igg]) continue;
         
            for(Int_t itg=kMB; itg<=kGA; itg++){ //@@@
               if(!trigflag[itg]) continue; //check which trigger fired
               
               fhCentralityTTC[itg][fkV0A][igg]->Fill(fCentralityV0A, fMultV0A); 
               fhCentralityTTC[itg][fkV0C][igg]->Fill(fCentralityV0C, fMultV0C); 
               fhCentralityTTC[itg][fkV0M][igg]->Fill(fCentralityV0M, fMultV0M); 
               fhCentralityTTC[itg][fkV0Mnorm1][igg]->Fill(fCentralityV0M, fMultV0Mnorm); 
         
               fhSignalTTC[itg][fkV0A][igg]->Fill(fMultV0A);
               fhSignalTTC[itg][fkV0C][igg]->Fill(fMultV0C);
               fhSignalTTC[itg][fkV0M][igg]->Fill(fMultV0M);
               fhSignalTTC[itg][fkV0Mnorm1][igg]->Fill(fMultV0Mnorm);
            }
         
            if(fIsMinBiasTrig){
               fhV0AvsV0CTTCinMB[igg]->Fill(fMultV0C, fMultV0A);
            }else if(fIsEmcalTrig){
               fhV0AvsV0CTTCinGA[igg]->Fill(fMultV0C, fMultV0A);
            }
         
            //Recoil jets 
            idx = fIndexTTC[igg];// gamma trigger
            if(idx>-1){
         
               for(Int_t itg=kMB; itg<=kGA; itg++){ //@@@
                  if(!trigflag[itg]) continue; //check which trigger fired
        
                  for(Int_t ir = 0; ir<kRho; ir++){
                     //fhDeltaPtTTC_RC_CentV0M[itg][igg]->Fill(fCentralityV0M, fdeltapT[igg]);
                     fhDeltaPtTTC_RC_V0Mnorm1[itg][igg][ir]->Fill(fMultV0Mnorm, fdeltapT[igg][ir]);
                  }
               }
         
               if(fFillSigTT && igg==0) continue;  // Do not fill reference 
               if(!fFillSigTT && igg>0) continue;  // Do not fill signal 
         
         
               for(Int_t itg=kMB; itg<=kGA; itg++){ //@@@
                  if(!trigflag[itg]) continue; //check which trigger fired
        
                  //fhTTC_CentV0M[itg][igg]->Fill(fCentralityV0M, fTTC[igg][idx].Pt()); //fill TTC trigger track pT
                  fhTTC_V0Mnorm1[itg][igg]->Fill(fMultV0Mnorm, fTTC[igg][idx].Pt()); //fill  TTC trigger track pT
               }
             
               //recoil jets
               for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
                  // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue; 
           
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC[igg][idx].Phi())) > fPhiCut){     
                     //recoil jet
                     for(Int_t ir=0; ir<kRho; ir++){ 
                        jetPtCorrDet = jet->Pt() - rho[ir]*jet->Area();
                        
                        for(Int_t itg=kMB; itg<=kGA; itg++){ //@@@
                           if(!trigflag[itg]) continue; //check which trigger fired
                        
                           //fhRecoilJetPtTTC_CentV0M[itg][igg]->Fill(fCentralityV0M, jetPtCorrDet);
                           fhRecoilJetPtTTC_V0Mnorm1[itg][igg][ir]->Fill(fMultV0Mnorm, jetPtCorrDet);
                        }
                     }
                  } 
               }  
            }
         }
      }//cluster container   
   
   
      //_________________________________________________________
      //      LOOP OVER JETS  DETECTOR LEVEL  TTJ
      //_________________________________________________________
    
      for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         jet = jetIterator.second;  // Get the pointer to jet object
         if(!jet)  continue; 
  
         tmparr[0] = jet->Pt();
         tmparr[1] = jet->Eta();
         tmparr[2] = jet->Phi();
         tmparr[3] = fMultV0Mnorm;

         for(Int_t ir=0; ir<kRho; ir++){ 
            jetPtCorrDet = jet->Pt() - rho[ir]*jet->Area();
            
            tmparr3[0] = jetPtCorrDet;
            tmparr3[1] = jet->Area();
            tmparr3[2] = fMultV0Mnorm;

            
            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!trigflag[itg]) continue; //check which trigger fired
           
               if(ir==0){ 
                  fhJetPhiIncl[itg]->Fill(jetPtCorrDet, jet->Phi());
                  fhJetEtaIncl[itg]->Fill(jetPtCorrDet, jet->Eta());
                  fhJetPtEtaPhiV0norm[itg]->Fill(tmparr);
               }
            
               fhJetPtAreaV0norm[itg][ir]->Fill(tmparr3);
            }
         }

         jetPtCorrDet = jet->Pt() - rho[krhokt]*jet->Area();
   
         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){ //search for TTJ candidates
            if(fJetChTTLowPt[ijj] < jetPtCorrDet && jetPtCorrDet < fJetChTTHighPt[ijj]){
               myTT.SetPtEtaPhiM(jetPtCorrDet, jet->Eta(), jet->Phi(), 0.); 
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
   
   
      for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
         if(!trigflag[itg]) continue; //check which trigger fired
         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
   
            fhMultTTJ[itg][ijj]->Fill(fJetChTT[ijj]); 
   
            if(!fJetChTT[ijj]) continue; //check if there is jet TT 
   
            fhRhoTTJ[itg][ijj]->Fill(rho[krhokt]);
         }
      }
   
      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         
         if(!fJetChTT[ijj]) continue; 
          
         for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
            if(!trigflag[itg]) continue; //check which trigger fired
    
            fhCentralityTTJ[itg][fkV0A][ijj]->Fill(fCentralityV0A, fMultV0A); 
            fhCentralityTTJ[itg][fkV0C][ijj]->Fill(fCentralityV0C, fMultV0C); 
            fhCentralityTTJ[itg][fkV0M][ijj]->Fill(fCentralityV0M, fMultV0M); 
            fhCentralityTTJ[itg][fkV0Mnorm1][ijj]->Fill(fCentralityV0M, fMultV0Mnorm); 
         
            fhSignalTTJ[itg][fkV0A][ijj]->Fill(fMultV0A);
            fhSignalTTJ[itg][fkV0C][ijj]->Fill(fMultV0C);
            fhSignalTTJ[itg][fkV0M][ijj]->Fill(fMultV0M);
            fhSignalTTJ[itg][fkV0Mnorm1][ijj]->Fill(fMultV0Mnorm);
         } 
         
         if(fIsMinBiasTrig){
            fhV0AvsV0CTTJ[ijj]->Fill(fMultV0C, fMultV0A);
         }
      }
   }//trigger selection for real, mc detector level, embedded



   //___________________________________________
   //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX

   fMultV0A_PartLevel = 0.;
   fMultV0C_PartLevel = 0.;
   fMultV0M_PartLevel = 0.;

   fMultV0Mnorm_PartLevel = 0.;
   fAsymV0M_PartLevel = 999.;

   Bool_t isMBpartlevel = 0; // requires coincidence of V0A and V0C on parton level

   if(fMode == AliAnalysisTaskEA::kMC){

      for(Int_t ir=0; ir<kRho; ir++){   
         fhRhoMBpart[ir]->Fill(rhoMC[ir]);
      }

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

            if(mcParticle->Charge()){ 
               if((static_cast<AliAODMCParticle*>(mcParticle))->IsPhysicalPrimary()){
                  //get particle level charged particles multiplicities in V0A and V0C
                  if(-3.7 < mcParticle->Eta() && mcParticle->Eta() < -1.7) fMultV0C_PartLevel++;
                  if( 2.8 < mcParticle->Eta() && mcParticle->Eta() < 5.1)  fMultV0A_PartLevel++; 
               } 
            }
         }

         if(fMultV0A_PartLevel>0 && fMultV0C_PartLevel>0) isMBpartlevel = kTRUE; //Minimum bias trigger on particle level
         //combined V0 multiplicities particle level
         fMultV0M_PartLevel = fMultV0A_PartLevel + fMultV0C_PartLevel;
         //fMultV0Anorm_PartLevel = fMultV0A_PartLevel/fMeanV0A_PartLevel;
         //fMultV0Cnorm_PartLevel = fMultV0C_PartLevel/fMeanV0C_PartLevel;
         fMultV0Mnorm_PartLevel = fMultV0M_PartLevel/fMeanV0M_PartLevel;

         if(isMBpartlevel){
            fHistEvtSelection->Fill(6.5); //Count Accepted input event 

            //pT spectrum of particle level physical primary mc particles
            for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum() ){
               mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
               if(!mcParticle)  continue; 
            
               if(IsTrackInAcceptance(mcParticle, kPartLevel)){
                  fhPtTrkTruePrimGen->Fill(mcParticle->Pt(), mcParticle->Eta(), fMultV0Mnorm_PartLevel);
            
                  for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                     if(fHadronTTLowPt[itt] < mcParticle->Pt() && mcParticle->Pt() < fHadronTTHighPt[itt]){
                        myTT.SetPtEtaPhiM(mcParticle->Pt(),mcParticle->Eta(),mcParticle->Phi(),0.); 
                        fTTH_PartLevel[itt].push_back(myTT);
                        fHadronTT_PartLevel[itt]++;   // there was a high pt 
                     }
                  }
               } 
             
            
               if(!mcParticle->Charge()){ 
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
            


            fhSignal_PartLevel[fkV0A]->Fill(fMultV0A_PartLevel);
            fhSignal_PartLevel[fkV0C]->Fill(fMultV0C_PartLevel);
            fhSignal_PartLevel[fkV0M]->Fill(fMultV0M_PartLevel);
            fhSignal_PartLevel[fkV0Mnorm1]->Fill(fMultV0Mnorm_PartLevel);
         

            Double_t meanV0Apart  =  fHelperEA->GetV0APartLevel(); 
            Double_t meanV0Cpart  =  fHelperEA->GetV0CPartLevel(); 
            if(meanV0Apart > 0 && meanV0Cpart > 0){
               Double_t multV0Anorm_part = fMultV0A_PartLevel/meanV0Apart;
               Double_t multV0Cnorm_part = fMultV0C_PartLevel/meanV0Cpart;
         
               if((multV0Anorm_part + multV0Cnorm_part)>0){
                  fAsymV0M_PartLevel = (multV0Anorm_part - multV0Cnorm_part )/(multV0Anorm_part + multV0Cnorm_part);
               }
            }
         
            fhV0MAssymVsV0Mnorm_PartLevel->Fill(fMultV0Mnorm_PartLevel, fAsymV0M_PartLevel);
         
            //chose trigger hadron TT   particle level
            for(Int_t itt=0; itt<fnHadronTTBins; itt++){
               if(fHadronTT_PartLevel[itt]>0){
                  fIndexTTH_PartLevel[itt] = fRandom->Integer(fHadronTT_PartLevel[itt]); 
                  idx = fIndexTTH_PartLevel[itt]; 
                 
                  for(Int_t ir=0; ir<kRho; ir++){ 
                     fdeltapT_PartLevel[itt][ir] = GetDeltaPt(fTTH_PartLevel[itt][idx].Phi(), fTTH_PartLevel[itt][idx].Eta(), phiLJmc, etaLJmc, phiSJmc, etaSJmc, rhoMC[ir], kPartLevel);
                  }
                  //signal in events with hadron TT   particle level
         
                  fhSignalTTH_PartLevel[fkV0A][itt]->Fill(fMultV0A_PartLevel);
                  fhSignalTTH_PartLevel[fkV0C][itt]->Fill(fMultV0C_PartLevel);
                  fhSignalTTH_PartLevel[fkV0M][itt]->Fill(fMultV0M_PartLevel);
                  fhSignalTTH_PartLevel[fkV0Mnorm1][itt]->Fill(fMultV0Mnorm_PartLevel);
                  
                  fhV0MAssymVsV0MnormTTH_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fAsymV0M_PartLevel);
                  
                  //hadron trigger particle level
                 
                  if(idx>-1){ 
        
                     for(Int_t ir=0; ir<kRho; ir++){ 
                        fhRhoTTHinMBpart[itt][ir]->Fill(rhoMC[ir]);
                        fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt][ir]->Fill(fMultV0Mnorm_PartLevel, fdeltapT_PartLevel[itt][ir]);
                     }
         
                     if(fFillSigTT && itt==0) continue;  // Do not fill reference 
                     if(!fFillSigTT && itt>0) continue;  // Do not fill signal 
         
                     fhTTH_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
                     fhTTH_3D_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fAsymV0M_PartLevel, fTTH_PartLevel[itt][idx].Pt());
                   
                     //recoil jets  PARTICLE LEVEL
                     for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
                        // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                        jet = jetIterator.second;  // Get the pointer to jet object
                        if(!jet)  continue; 
         
                        dphi = TVector2::Phi_mpi_pi(jet->Phi()-fTTH_PartLevel[itt][idx].Phi());                       

                        for(Int_t ir=0; ir<kRho; ir++){ 
                           jetPtCorrDet = jet->Pt() - rhoMC[ir]*jet->Area();
                           
                           fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt][ir]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet, dphi);
                           
                           tmparr[0] = fMultV0Mnorm_PartLevel;
                           tmparr[1] = fAsymV0M_PartLevel;
                           tmparr[2] = jetPtCorrDet;
                           tmparr[3] = TMath::Abs(dphi);
                           fhRecoilJetTTH_V0Mnorm1_PartLevel[itt][ir]->Fill(tmparr); 
                           
                           
                           if(TMath::Abs(dphi) > fPhiCut){     
                              //recoil jet hadron trigger
                           
                              fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt][ir]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet);
                           }
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
        
                  for(Int_t ir=0; ir<kRho; ir++){ 
                     fdeltapT_PartLevel[igg][ir] = GetDeltaPt(fTTC_PartLevel[igg][idx].Phi(), fTTC_PartLevel[igg][idx].Eta(), phiLJmc, etaLJmc, phiSJmc, etaSJmc, rhoMC[ir], kPartLevel);
                  }
         
                  //signal in events with hadron TT   particle level
                  fhSignalTTC_PartLevel[fkV0A][igg]->Fill(fMultV0A_PartLevel);
                  fhSignalTTC_PartLevel[fkV0C][igg]->Fill(fMultV0C_PartLevel);
                  fhSignalTTC_PartLevel[fkV0M][igg]->Fill(fMultV0M_PartLevel);
                  fhSignalTTC_PartLevel[fkV0Mnorm1][igg]->Fill(fMultV0Mnorm_PartLevel);
         
                  if(idx>-1){ 
         
                     for(Int_t ir=0; ir<kRho; ir++){ 
                        fhRhoTTCinMBpart[igg][ir]->Fill(rhoMC[ir]);
                        fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[igg][ir]->Fill(fMultV0Mnorm_PartLevel, fdeltapT_PartLevel[igg][ir]);
                     }

                     if(fFillSigTT && igg==0) continue;  // Do not fill reference 
                     if(!fFillSigTT && igg>0) continue;  // Do not fill signal 
         
                     fhTTC_V0Mnorm1_PartLevel[igg]->Fill(fMultV0Mnorm_PartLevel, fTTC_PartLevel[igg][idx].Pt()); //fill trigger track pT
         
                     //recoil jets PARTICLE LEVEL
                     for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
                        // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                        jet = jetIterator.second;  // Get the pointer to jet object
                        if(!jet)  continue; 
                    
                        if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTC_PartLevel[igg][idx].Phi())) > fPhiCut){     
                           //recoil jet
                           for(Int_t ir=0; ir<kRho; ir++){ 
                              jetPtCorrDet = jet->Pt() - rhoMC[ir]*jet->Area();
                              fhRecoilJetPtTTC_V0Mnorm1_PartLevel[igg][ir]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet);
                           }
                        }
                     } 
                  }
               }
            }
         }
      }


     if(isMBpartlevel && fIsMinBiasTrig){
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
                       fhPtTrkTruePrimRec->Fill(mcParticle->Pt(), mcParticle->Eta(), fMultV0Mnorm_PartLevel); //this is well recontr phys primary
                       break;
                    }//same label with rec particle
                 }
              }//loop over gen tracks
              if(!bRecPrim){
                 fhPtTrkSecOrFakeRec->Fill(track->Pt(), track->Eta(), fMultV0Mnorm_PartLevel); //matchnig to phys primary not found, this is fake or second.
              }
           }
        }
       
        //__________________________________________________________
        //  FILL JET RESPONSE MATRIX
        //__________________________________________________________
       
        //Response matrix normalization - spectrum of all generator level jets in acceptance
        if(fJetContainerPartLevel){
           for(auto jetPartIterator : fJetContainerPartLevel->accepted_momentum() ){
              jetPartMC = jetPartIterator.second;  // Get the pointer to mc particle object
              if(!jetPartMC)  continue; 
      
              for(Int_t ir=0; ir<kRho; ir++){ 
                 jetPtCorrPart = jetPartMC->Pt() - jetPartMC->Area()*rhoMC[ir];
                 
                 fhJetPtPartLevelCorr[ir]->Fill(jetPtCorrPart);
                 if(ir==0) fhJetPtPartLevelZero->Fill(jetPartMC->Pt());
                 
                 
                 tmparr3[0] = jetPtCorrPart;
                 tmparr3[1] = jetPartMC->Area();
                 tmparr3[2] = fMultV0Mnorm_PartLevel;
                 
                 fhJetPtAreaV0norm_PartLevel[ir]->Fill(tmparr3); 
                 
                 for(Int_t itt=0; itt<fnHadronTTBins; itt++){ //event contains a particle level trigger
                    if(fHadronTT[itt]>0){
                       idx = fIndexTTH[itt];
                       dphi = TVector2::Phi_mpi_pi(jetPartMC->Phi()-fTTH[itt][idx].Phi()); 
                 
                       if(TMath::Abs(dphi) > fPhiCut){    //fill with recoil jets only 
                          fhJetPtPartLevelCorrTTHdl[itt][ir]->Fill(jetPtCorrPart);
                       }
                    }
                 }
              }
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
                 for(Int_t ir=0; ir<kRho; ir++){
                    jetPtCorrDet = jet->Pt() - jet->Area()*rho[ir];
                    fhFractionOfSecInJet[ir]->Fill(jetPtCorrDet, sumFakeTrackPtInJet/sumAllTrackPtInJet);
                 }  
              }
       
              //Fill Response matrix
              jetPartMC =  jet->ClosestJet();
              if(!jetPartMC) continue;
              if(jetPartMC->Pt()<1e-3) continue; //prevents matching with a ghost
       
              for(Int_t ir=0; ir<kRho; ir++){
                 jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*rhoMC[ir]; 
                 jetPtCorrDet  =  jet->Pt() - jet->Area()*rho[ir]; 
                 
                 fhJetPtPartLevelVsJetPtDetLevelCorr[ir]->Fill(jetPtCorrDet,jetPtCorrPart); //response matrix
                 if(ir==0) fhJetPtPartLevelVsJetPtDetLevelZero->Fill(jet->Pt(),jetPartMC->Pt()); //response matrix
                 
                 for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                    if(fHadronTT[itt]>0){
                       idx = fIndexTTH[itt];
                       dphi = TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi()); 
                 
                       if(TMath::Abs(dphi) > fPhiCut){    //fill with recoil jets only 
                          fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[itt][ir]->Fill(jetPtCorrDet,jetPtCorrPart); 
                       } 
                    }
                 }
                 
                 
                 if(jetPtCorrPart>0){
                    fhJetPtResolutionVsPtPartLevel[ir]->Fill(jetPtCorrPart,(jetPtCorrDet-jetPtCorrPart)/jetPtCorrPart); //jet pT resolution
                 }
              } 
           }
        }
      }
   }

 
   //___________________________________________________________


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
   delete fHelperEA; 
 
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
   TString rhotype[kRho] = {"kt","cms"};

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


   //Trigger track pT spectrum single inclusive for  MB  versus V0M
   Int_t    nbinsV0M     = 100;
   Double_t maxV0M       = 1000.;
   Double_t maxV0Mmc     = 500.;
   Int_t    nbinsV0Mnorm = 200;
   Double_t maxV0Mnorm   = 20.;
 

   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   fhVertexZall =  new TH1D("fhVertexZall","z vertex without cut",40,-20,20);
   fOutput->Add(fhVertexZall); 
 
   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);
   fOutput->Add(fhVertexZ);
 
   //-------------------------
   TString trig[]={"MB","HM","GA"};


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 

      name   = Form("fhTrackEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive track vs pT %s",trig[itg].Data());
      fhTrackEtaIncl[itg] = new TH2D( name.Data(), object.Data(), 50,0, 100, 40,-0.9,0.9);
      fOutput->Add((TH2D*) fhTrackEtaIncl[itg]);
 
      name   = Form("fhTrackPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist tracks vs pT %s",trig[itg].Data());
      fhTrackPhiIncl[itg] = new TH2D( name.Data(), object.Data(), 50, 0, 100, 50,0,2*TMath::Pi());
      fOutput->Add((TH2D*) fhTrackPhiIncl[itg]);
   }

   if(fMode == AliAnalysisTaskEA::kEmbedding){
     fhTrackEtaInclEMB = (TH2D*) fhTrackEtaIncl[kMB]->Clone("fhTrackEtaInclEMB");
     fOutput->Add((TH2D*) fhTrackEtaInclEMB);
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 

      name   = Form("fhJetEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive jets vs pTjet %s",trig[itg].Data());
      fhJetEtaIncl[itg] = new TH2D(name.Data(),object.Data(), 150, -20, 130, 40,-0.9,0.9);
      fOutput->Add((TH2D*) fhJetEtaIncl[itg]);
 
      name   = Form("fhJetPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist jets vs pTjet %s",trig[itg].Data());
      fhJetPhiIncl[itg] = new TH2D(name.Data(),object.Data(), 60, -20, 100, 50, 0, 2*TMath::Pi());
      fOutput->Add((TH2D*) fhJetPhiIncl[itg]);
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      name   = Form("fhClusterEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive clusters vs pT %s",trig[itg].Data());
      fhClusterEtaIncl[itg] = new TH2D( name.Data(), object.Data(), 100, 0, 100, 40,-0.9,0.9);
      fOutput->Add((TH2D*) fhClusterEtaIncl[itg]);
 
      name   = Form("fhClusterPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist clusters vs pT %s",trig[itg].Data());
      fhClusterPhiIncl[itg] = new TH2D( name.Data(), object.Data(), 50, 0, 100, 50,0,2*TMath::Pi());
      fOutput->Add((TH2D*) fhClusterPhiIncl[itg]);
   }

   const Int_t ktdim = 4;
   Int_t   tbins[ktdim] = { 50,   40,         140, 10};
   Double_t txmin[ktdim] = { 0., -0.9, 0,  0.};  
   Double_t txmax[ktdim] = {50.,  0.9, 2*TMath::Pi(), 10.};  


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      name = Form("fhTrackPtEtaPhiV0norm_%s",trig[itg].Data());

      fhTrackPtEtaPhiV0norm[itg] = new  THnSparseF(name.Data(),"Tracks pt eta phi V0nom", ktdim, tbins, txmin, txmax);
      fOutput->Add((THnSparse*) fhTrackPtEtaPhiV0norm[itg]); 
   } 

   //jets
   const Int_t kjetdim = 4;
   Int_t   jetbins[kjetdim] = {100,   40, 140, 10};
   Double_t jetxmin[kjetdim] = { 0., -0.9,  0,  0.};  
   Double_t jetxmax[kjetdim] = {100.,  0.9, 2*TMath::Pi(), 10.};  


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      name = Form("fhJetPtEtaPhiV0norm_%s",trig[itg].Data());

      fhJetPtEtaPhiV0norm[itg] = new  THnSparseF(name.Data(),"Jet pt eta phi V0nom", kjetdim, jetbins, jetxmin, jetxmax);
      fOutput->Add((THnSparse*) fhJetPtEtaPhiV0norm[itg]); 
   } 

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         name = Form("fhJetPtEtaPhiV0norm_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);

         fhJetPtEtaPhiV0normTTH[itg][itt] = new  THnSparseF(name.Data(),"Jet pt eta phi V0nom", kjetdim, jetbins, jetxmin, jetxmax);
         fOutput->Add((THnSparse*) fhJetPtEtaPhiV0normTTH[itg][itt]); 
      }
   } 



   const Int_t kjdim = 3;
   Int_t   jbins[ktdim]  = {110,  50, 10};
   Double_t jxmin[ktdim] = {-10.,  0,  0.};  
   Double_t jxmax[ktdim] = {100.,  2, 10.};  

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhJetPtAreaV0norm_%s_Rho%s",trig[itg].Data(), rhotype[ir].Data());
         
         fhJetPtAreaV0norm[itg][ir] = new  THnSparseF(name.Data(),"Jet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
         fOutput->Add((THnSparse*) fhJetPtAreaV0norm[itg][ir]); 
      } 
   }
   
   if(fMode == AliAnalysisTaskEA::kMC){
      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhJetPtAreaV0norm_Rho%s_PartLevel", rhotype[ir].Data());
         
         fhJetPtAreaV0norm_PartLevel[ir] = new  THnSparseF(name.Data(),"Part LevelJet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
         fOutput->Add((THnSparse*) fhJetPtAreaV0norm_PartLevel[ir]);
      } 
   } 

   //RHO 
   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      for(Int_t ir=0; ir<kRho; ir++){
         name   = Form("hRho%s_%s",rhotype[ir].Data(),trig[itg].Data());
         object = Form("Rho %s det level %s",rhotype[ir].Data(),trig[itg].Data());
         
         fhRho[itg][ir] = new TH1D( name.Data(), object.Data(),1000,0,100);
         fOutput->Add((TH1D*) fhRho[itg][ir]); 
      }
   } 

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         for(Int_t ir=0; ir<kRho; ir++){
            name = Form("hRho%s_%s_TTH%d_%d", rhotype[ir].Data(), trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhRhoTTH[itg][itt][ir] = (TH1D*)  fhRho[itg][ir]->Clone(name.Data());      //! in events MB with hadron TT
            fOutput->Add((TH1D*) fhRhoTTH[itg][itt][ir]); 
         }
      }

      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){  //JET TT
         name = Form("hRho%s_%s_TTJ%d_%d", rhotype[krhokt].Data(), trig[itg].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
         fhRhoTTJ[itg][ijj] = (TH1D*)  fhRho[itg][krhokt]->Clone(name.Data());                      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTJ[itg][ijj]); 
      } 
   }
 
   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){ //GAMMA TT
         for(Int_t ir=0; ir<kRho; ir++){
            name = Form("hRho%s_%s_TTC%d_%d", rhotype[ir].Data(), trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
            fhRhoTTC[itg][igg][ir] = (TH1D*)  fhRho[itg][ir]->Clone(name.Data());                      //! in events MB with hadron TT
            fOutput->Add((TH1D*) fhRhoTTC[itg][igg][ir]); 
         }
      }
   } 

   if(fMode == AliAnalysisTaskEA::kMC){
      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("hRho%s_MB_part", rhotype[ir].Data());
         fhRhoMBpart[ir] = (TH1D*)  fhRho[kMB][ir]->Clone(name.Data()); 
         fhRhoMBpart[ir]->SetTitle(Form("Rho %s  min bias part level", rhotype[ir].Data())); 
         fOutput->Add((TH1D*) fhRhoMBpart[ir]); 
         
         for(Int_t itt=0; itt<fnHadronTTBins;itt++){
            name = Form("hRho%s_MB_TTH%d_%d_part", rhotype[ir].Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
            fhRhoTTHinMBpart[itt][ir] = (TH1D*)  fhRho[kMB][ir]->Clone(name.Data());                      //! in events MB with hadron TT
            fOutput->Add((TH1D*) fhRhoTTHinMBpart[itt][ir]); 
         }
         
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            name = Form("hRho%s_MB_TTC%d_%d_part", rhotype[ir].Data(), fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
            fhRhoTTCinMBpart[igg][ir] = (TH1D*)  fhRho[kMB][ir]->Clone(name.Data());                      //! in events MB with hadron TT
            fOutput->Add((TH1D*) fhRhoTTCinMBpart[igg][ir]); 
         }
      }
   }

   //VERTEX
   fhVertex[0] = new TH1D("hVertexX","VertexX",100,-1,1);
   fhVertex[1] = new TH1D("hVertexY","VertexY",100,-1,1);
   fhVertex[2] = new TH1D("hVertexZ","VertexZ",400,-20,20);

   for(Int_t iv=0; iv<fkVtx;iv++){
      fOutput->Add((TH1D*) fhVertex[iv]); 
   }


   Int_t nRun = fHelperEA->GetNRuns();

 
  fhV0MRunByRunMB = new TH2D("fhV0MRunByRunMB","fhV0MRunByRunMB", nRun, 0, nRun, 180,0,1800); 
   for(Int_t ir=0; ir < nRun; ir++){
      fhV0MRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d", fHelperEA->GetRun(ir)));
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
 
   fhV0MnormRunByRunMB = new TH2D("fhV0MnormRunByRunMB","fhV0MnormRunByRunMB", nRun, 0, nRun, 200,0,20); 
   for(Int_t ir=0; ir < nRun; ir++){
      fhV0MnormRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d", fHelperEA->GetRun(ir)));
   } 
   fOutput->Add((TH2D*) fhV0MnormRunByRunMB);
 

   //CENTRALITY
   TString cest[] = {"V0A", "V0C", "V0M", "V0Mnorm"}; //centrality estimators

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


   //CENTRALITY
   for(Int_t itg=kMB; itg<=kGA;itg++){
      for(Int_t ic=0; ic<fkCE;ic++){
         if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
         if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

         name = Form("hCentrality_%s_%s", trig[itg].Data(), cest[ic].Data());
         if(ic!=fkV0Mnorm1){
            fhCentrality[itg][ic] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, narrV0, arrV0);
         }else{
            fhCentrality[itg][ic] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 400,0,20);
         }
         fOutput->Add((TH2D*) fhCentrality[itg][ic]); 
      }
   }
 
 
   //CENTRALITY TTH 
   for(Int_t itg=kMB; itg<=kHM;itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
   
      for(Int_t ic=0; ic<fkCE; ic++){
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hCentrality_%s_%s_TTH%d_%d", trig[itg].Data(), cest[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhCentralityTTH[itg][ic][itt] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
            fOutput->Add((TH2D*) fhCentralityTTH[itg][ic][itt]); 
         }
      }
   }
 
   //TTJ MB
   for(Int_t itg=kMB; itg<=kHM;itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
 
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
            name = Form("hCentrality_%s_%s_TTJ%d_%d", trig[itg].Data(), cest[ic].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
            fhCentralityTTJ[itg][ic][ijj] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
            fOutput->Add((TH2D*) fhCentralityTTJ[itg][ic][ijj]); 
         }
      }
   }

   //TTC  MB 
   for(Int_t itg=kMB; itg<=kGA;itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
 
      for(Int_t ic=0; ic<fkCE;ic++){
         for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
            name = Form("hCentrality_%s_%s_TTC%d_%d", trig[itg].Data(), cest[ic].Data(), fClusterTTLowPt[ijj], fClusterTTHighPt[ijj]);
            fhCentralityTTC[itg][ic][ijj] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
            fOutput->Add((TH2D*) fhCentralityTTC[itg][ic][ijj]); 
         }
      }
   }

   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //SIGNAL

   TString signal[]={"multV0A", "multV0C", "multV0M","multV0Mnorm"};
   Float_t signalL[]={0,0,0,0};
   Float_t signalH[]={1000,1000,1800,15};
   Int_t   signalN[]={100,100,180,150};

   for(Int_t itg=kMB; itg<=kGA;itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
 
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_%s_%s", trig[itg].Data(), signal[ic].Data());
         fhSignal[itg][ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignal[itg][ic]); 
      }
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
 
      for(Int_t ic=0; ic<fkCE; ic++){ //TT hadron
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_%s_%s_TTH%d_%d", trig[itg].Data(), signal[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTH[itg][ic][itt] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
            fOutput->Add((TH1D*) fhSignalTTH[itg][ic][itt]); 
         }
      }
 
      for(Int_t ic=0; ic<fkCE; ic++){ //TT jet
         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
            name = Form("hSignal_%s_%s_TTJ%d_%d", trig[itg].Data(), signal[ic].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
            fhSignalTTJ[itg][ic][ijj] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
            fOutput->Add((TH1D*) fhSignalTTJ[itg][ic][ijj]); 
         }
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
 
      for(Int_t ic=0; ic<fkCE; ic++){ //HM && TT jet     
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            name = Form("hSignal_%s_%s_TTC%d_%d",trig[itg].Data(), signal[ic].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
            fhSignalTTC[itg][ic][igg] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
            fOutput->Add((TH1D*) fhSignalTTC[itg][ic][igg]);   
         }
      }
   }



   if(fMode == AliAnalysisTaskEA::kMC){ //PARTICLE LEVEL SIGNAL DISTRIBUTIONS
      TString signalmc[]={"multV0A", "multV0C", "multV0M", "multV0Mnorm"};
      Float_t signalLmc[]={0,0,0,0};
      Float_t signalHmc[]={500,500,500,20};
      Int_t signalNmc[]={500,500,500,200};
      
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_MB_%s_PartLevel", signalmc[ic].Data());
         fhSignal_PartLevel[ic] = new TH1D(name.Data(), name.Data(), signalNmc[ic], signalLmc[ic], signalHmc[ic]);
         fOutput->Add((TH1D*) fhSignal_PartLevel[ic]); 
      }

      //TT hadron
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_MB_%s_TTH%d_%d_PartLevel", signalmc[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTH_PartLevel[ic][itt] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
            fOutput->Add((TH1D*) fhSignalTTH_PartLevel[ic][itt]); 
         }
      }

      //TT cluster
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
            name = Form("hSignal_MB_%s_TTC%d_%d_PartLevel", signalmc[ic].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
            fhSignalTTC_PartLevel[ic][igg] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
            fOutput->Add((TH1D*) fhSignalTTC_PartLevel[ic][igg]); 
         }
      }
   }
 

   // V0 assymetery versus V0norm
   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
   
      name = Form("fhV0MAssymVsV0Mnorm_%s",trig[itg].Data()); 
      fhV0MAssymVsV0Mnorm[itg] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
      fOutput->Add((TH2D*) fhV0MAssymVsV0Mnorm[itg]); 
   }

   if(fMode == AliAnalysisTaskEA::kMC){
      name = Form("fhV0MAssymVsV0Mnorm_MB_PartLevel");
      fhV0MAssymVsV0Mnorm_PartLevel = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
      fOutput->Add((TH2D*) fhV0MAssymVsV0Mnorm_PartLevel); 
   } 
   for(Int_t itg=kMB; itg<=kGA; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
   
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("fhV0MAssymVsV0Mnorm_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhV0MAssymVsV0MnormTTH[itg][itt] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
         fOutput->Add((TH2D*) fhV0MAssymVsV0MnormTTH[itg][itt]); 
      }
   }

   if(fMode == AliAnalysisTaskEA::kMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("fhV0MAssymVsV0Mnorm_MB_TTH%d_%d_PartLevel",  fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhV0MAssymVsV0MnormTTH_PartLevel[itt] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
         fOutput->Add((TH2D*) fhV0MAssymVsV0MnormTTH_PartLevel[itt]);
      } 
   } 



   name = Form("fhV0AvsV0C_MB");
   fhV0AvsV0C = new TH2D(name.Data(),name.Data(),100,0,1000, 100,0,1000);
   fOutput->Add((TH2D*) fhV0AvsV0C); 

   //name = Form("fhV0MvsV0Mnorm_MB");
   //fhV0MvsV0Mnorm = new TH2D(name.Data(),name.Data(),100,0,40, 100,0,1200);
   //fOutput->Add((TH2D*) fhV0MvsV0Mnorm); 


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

    
   if(fMode != AliAnalysisTaskEA::kMC){ 
      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
         name = Form("fhV0AvsV0C_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
         fhV0AvsV0CTTCinGA[ijj] = (TH2D*)  fhV0AvsV0C->Clone(name.Data());
         fhV0AvsV0CTTCinGA[ijj]->SetTitle(name.Data()); 
         fOutput->Add((TH2D*) fhV0AvsV0CTTCinGA[ijj]); 
      } 
   }
   //+++++++++++++++++++++++++++++++
   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      name = Form("fhTrackMult%s", trig[itg].Data());
      fhTrackMult[itg] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 1000, 0, 1000); 
      fOutput->Add((TH2D*) fhTrackMult[itg]); 

      name = Form("fhMeanTrackPt%s", trig[itg].Data());
      fhMeanTrackPt[itg] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 20);
      fOutput->Add((TH1D*) fhMeanTrackPt[itg]); 
   }


   //Trigger track candidate multiplicity
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hMultTT_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhMultTTH[itg][itt] = new TH1D(name.Data(),name.Data(),100,0,100);
         fOutput->Add((TH1D*)  fhMultTTH[itg][itt]); 
      }
   } 

   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTJ
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
         name = Form("hMultTT_%s_TTJ%d_%d", trig[itg].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
         fhMultTTJ[itg][ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
         fOutput->Add((TH1D*) fhMultTTJ[itg][ijj]); 
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTC
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hMultTT_%s_TTC%d_%d", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
         fhMultTTC[itg][igg] = new TH1D(name.Data(),name.Data(),100,0,100);
         fOutput->Add((TH1D*) fhMultTTC[itg][igg]); 
      }
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      //for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      //    name = Form("hTT_%s_TTH%d_%d_CentV0M", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      //    fhTTH_CentV0M[itg][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
      //    fOutput->Add((TH2D*) fhTTH_CentV0M[itg][itt]); 
      //}

      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_%s_TTH%d_%d_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm1[itg][itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm1[itg][itt]); 
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_%s_3D_TTH%d_%d_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhTTH_3D_V0Mnorm1[itg][itt] = new TH3D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 21, -1, 1.1, 100, 0, 100);
         fOutput->Add((TH3D*) fhTTH_3D_V0Mnorm1[itg][itt]); 
      }

   }


   if(fMode == AliAnalysisTaskEA::kMC){
      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm1_PartLevel[itt]); 
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_3D_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTH_3D_V0Mnorm1_PartLevel[itt] = new TH3D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 21, -1, 1.1, 100, 0, 100);
         fOutput->Add((TH3D*) fhTTH_3D_V0Mnorm1_PartLevel[itt]); 
      }
   }


   //TT emcal cluster pT spectrum single inclusive  in MB   with V0M
   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      //   name = Form("hTT_%s_TTC%d_%d_CentV0M", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
      //   fhTTC_CentV0M[itg][igg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
      //   fOutput->Add((TH2D*) fhTTC_CentV0M[itg][igg]); 
      //}
     
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hTT_%s_TTC%d_%d_V0Mnorm", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
         fhTTC_V0Mnorm1[itg][igg] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTC_V0Mnorm1[itg][igg]); 
      }
   }   
  
   if(fMode == AliAnalysisTaskEA::kMC){
      //TT emcal cluster pT spectrum single inclusive  in MB   with V0Mnorm
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         name = Form("hTT_MB_TTC%d_%d_V0Mnorm_PartLevel", fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
         fhTTC_V0Mnorm1_PartLevel[igg] = new TH2D(name.Data(), name.Data(),  nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTC_V0Mnorm1_PartLevel[igg]); 
      }
   }



   //RECOIL JET SPECTRA
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      //for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      //   name = Form("fhRecoilJetPt_%s_TTH%d_%d_CentV0M", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      //   fhRecoilJetPtTTH_CentV0M[itg][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);            
      //   fOutput->Add((TH2D*) fhRecoilJetPtTTH_CentV0M[itg][itt]); 
      //}
      
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         for(Int_t ir=0; ir<kRho; ir++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
            name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhRecoilJetPtTTH_V0Mnorm1[itg][itt][ir] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);            
            fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm1[itg][itt][ir]);
         } 
      }
   }

   //reference with shifted rho
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t ir=0; ir<kRho; ir++){  //TTH
         name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rho%s_rhoShift1", trig[itg].Data(), fHadronTTLowPt[0], fHadronTTHighPt[0], rhotype[ir].Data());
     
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[itg][ir] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[itg][ir]);
      }
      
      for(Int_t ir=0; ir<kRho; ir++){  //TTH
         name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rho%s_rhoShift2", trig[itg].Data(), fHadronTTLowPt[0], fHadronTTHighPt[0], rhotype[ir].Data());
         
         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[itg][ir] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[itg][ir]);
      }
   }
 


   Double_t pTbins3[]   = {-20,-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,140,160,180,200};  
   const Int_t npTbins3 = sizeof(pTbins3)/sizeof(Double_t)-1;

   const Int_t narrPhi=100;
   Double_t arrPhi[narrPhi+1];
   Double_t p = TMath::TwoPi()/narrPhi;
   for(Int_t i=0; i<=narrPhi; i++) arrPhi[i] = -TMath::Pi() + i*p;

   Double_t arrV0Mnorm[nbinsV0Mnorm+1];
   p = maxV0Mnorm/nbinsV0Mnorm;
   for(Int_t i=0; i<=nbinsV0Mnorm; i++) arrV0Mnorm[i] = i*p;


  // dphi of recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);    
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
 
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){ 
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhRecoilJetPhi_%s_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, npTbins3, pTbins3, narrPhi, arrPhi); 
            fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir]);
         } 
      }
   }

   //recoil jet distribution as a function V0norm, V0 assymetery, jet pt, jet |dphi|
   const Int_t rldim = 4;
   Int_t   rlbins[ktdim] = {10,  21, 130, 50};
   Double_t rlmin[ktdim] = { 0., -1, -10,  0.};  
   Double_t rlmax[ktdim] = {10., 1.1, 120, TMath::Pi()};  


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue; 
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue; 

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         for(Int_t ir=0; ir<kRho; ir++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
            name = Form("fhRecoilJet4D_%s_TTH%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            
            fhRecoilJetTTH_V0Mnorm1[itg][itt][ir] = new  THnSparseF(name.Data(),"V0norm, V0 assym, jet pT,  abs(dphi)", rldim, rlbins, rlmin, rlmax);
            fOutput->Add((THnSparse*) fhRecoilJetTTH_V0Mnorm1[itg][itt][ir]);
         }
      } 
   }



   //TTH recoil jet distributions for MC 
   if(fMode == AliAnalysisTaskEA::kMC){ // particle level 
      
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         for(Int_t ir=0; ir<kRho; ir++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
            name = Form("fhRecoilJetPt_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt][ir] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, 200, -20, 180);            
            fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt][ir]);
         } 
      }

      // dphi of recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){ 
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhRecoilJetPhi_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt][ir] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, npTbins3, pTbins3, narrPhi, arrPhi); 
            fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt][ir]);
         } 
      }


      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         for(Int_t ir=0; ir<kRho; ir++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
            name = Form("fhRecoilJet4D_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            
            fhRecoilJetTTH_V0Mnorm1_PartLevel[itt][ir] = new  THnSparseF(name.Data(),"V0norm, V0 assym, jet pT,  abs(dphi)", rldim, rlbins, rlmin, rlmax);
            fOutput->Add((THnSparse*) fhRecoilJetTTH_V0Mnorm1_PartLevel[itt][ir]);
         }
      } 
   }


  
   if(fMode == AliAnalysisTaskEA::kEmbedding){ 
      //! dphi of recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);    
   
      for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){ 
            for(Int_t ir=0; ir<kRho; ir++){ 
        
               //!  filled with any detector level pythia recoil jet 
               name = Form("fhRecoilJetPhi_%s_EMB_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
               fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt][ir] = (TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir]->Clone(name.Data());
               fOutput->Add((TH3D*) fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt][ir]); 
               
               //!  filled  tagged closest detector level pythia recoil jet 
               name = Form("fhRecoilJetPhi_%s_TAG_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
               fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt][ir] = (TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt][ir]->Clone(name.Data());
               fOutput->Add((TH3D*) fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt][ir]);
            } 
         }
      }
   }

   //+++++++++++++++++++++++++++ RECOIL JETS WITH TTC ++++++++++++++++++++++++++++++++++++++++++
   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M centrality
      //   name = Form("fhRecoilJetPt_%s_TTC%d_%d_CentV0M", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
      //   fhRecoilJetPtTTC_CentV0M[itg][igg] = (TH2D*) fhRecoilJetPtTTH_CentV0M[0][0]->Clone(name.Data()); 
      //   fOutput->Add((TH2D*) fhRecoilJetPtTTC_CentV0M[itg][igg]); 
      //}
      
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
         for(Int_t ir=0; ir<kRho; ir++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
            name = Form("fhRecoilJetPt_%s_TTC%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg], rhotype[ir].Data());
            fhRecoilJetPtTTC_V0Mnorm1[itg][igg][ir] = (TH2D*) fhRecoilJetPtTTH_V0Mnorm1[0][0][0]->Clone(name.Data()); 
            fOutput->Add((TH2D*) fhRecoilJetPtTTC_V0Mnorm1[itg][igg][ir]);
         } 
      }
   }

   if(fMode == AliAnalysisTaskEA::kMC){ //particle level
   
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
         for(Int_t ir=0; ir<kRho; ir++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
            name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0Mnorm_Rho%s_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg], rhotype[ir].Data());
            fhRecoilJetPtTTC_V0Mnorm1_PartLevel[igg][ir] = (TH2D*) fhRecoilJetPtTTH_V0Mnorm1_PartLevel[0][0]->Clone(name.Data()); 
            fOutput->Add((TH2D*) fhRecoilJetPtTTC_V0Mnorm1_PartLevel[igg][ir]);
         } 
      }
   }


   //delta pT distributions versus V0M CENTRALITY
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
 
      //for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      //   name = Form("fhDeltaPtTTH_%s_RC_CentV0M_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      //   fhDeltaPtTTH_RC_CentV0M[itg][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);            
      //   fOutput->Add((TH2D*) fhDeltaPtTTH_RC_CentV0M[itg][itt]); 
      //}
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
 
      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){
      //   name = Form("fhDeltaPtTTC_%s_RC_CentV0M_TTC%d_%d", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
      //   fhDeltaPtTTC_RC_CentV0M[itg][igg] = (TH2D*) fhDeltaPtTTH_RC_CentV0M[0][0]->Clone(name.Data()); 
      //   fOutput->Add((TH2D*) fhDeltaPtTTC_RC_CentV0M[itg][igg]); 
      //}
   }

   //delta pT distributions versus V0Mnorm   = V0M/mean V0M
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhDeltaPtTTH_%s_RC_V0Mnorm_TTH%d_%d_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhDeltaPtTTH_RC_V0Mnorm1[itg][itt][ir] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);            
            fOutput->Add((TH2D*) fhDeltaPtTTH_RC_V0Mnorm1[itg][itt][ir]); 
         }
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
 
      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhDeltaPtTTC_%s_RC_V0Mnorm_TTC%d_%d_Rho%s", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg], rhotype[ir].Data());
            fhDeltaPtTTC_RC_V0Mnorm1[itg][igg][ir] =   (TH2D*) fhDeltaPtTTH_RC_V0Mnorm1[0][0][0]->Clone(name.Data()); 
            fOutput->Add((TH2D*) fhDeltaPtTTC_RC_V0Mnorm1[itg][igg][ir]);
         } 
      }
   }

   
   if(fMode == AliAnalysisTaskEA::kMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhDeltaPtTTH_MB_RC_V0Mnorm_TTH%d_%d_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt][ir] = (TH2D*)  fhDeltaPtTTH_RC_V0Mnorm1[0][0][0]->Clone(name.Data()); 
            fOutput->Add((TH2D*) fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt][ir]);
         } 
      }

      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
         for(Int_t ir=0; ir<kRho; ir++){ 
            name = Form("fhDeltaPtTTC_RC_V0Mnorm_TTC%d_%d_Rho%s_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg], rhotype[ir].Data());
            fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[igg][ir] =   (TH2D*) fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[0][0]->Clone(name.Data()); 
            fOutput->Add((TH2D*) fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[igg][ir]);
         } 
      }
   }

   if(fMode == AliAnalysisTaskEA::kMC){
      fhPtTrkTruePrimGen = new TH3D("fhPtTrkTruePrimGen","fhPtTrkTruePrimGen",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkTruePrimGen); 

      fhPtTrkTruePrimRec = new TH3D("fhPtTrkTruePrimRec","fhPtTrkTruePrimRec",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkTruePrimRec); 

      fhPtTrkSecOrFakeRec = new TH3D("fhPtTrkSecOrFakeRec","fhPtTrkSecOrFakeRec",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkSecOrFakeRec); 
     
      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhJetPtPartLevelCorr_Rho%s", rhotype[ir].Data());
         fhJetPtPartLevelCorr[ir] = new TH1D(name.Data(), name.Data(),270,-20,250);
         fOutput->Add((TH1D*) fhJetPtPartLevelCorr[ir]);
      }

      fhJetPtPartLevelZero = new TH1D("fhJetPtPartLevelZero","fhJetPtPartLevelZero",250,0,250);
      fOutput->Add((TH1D*) fhJetPtPartLevelZero);

      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhFractionOfSecInJet_Rho%s", rhotype[ir].Data());
         fhFractionOfSecInJet[ir] = new TH2D(name.Data(), "Frac of jet pT carried by secondary tracks",50,0,50,210,0,1.05); 
         fOutput->Add((TH2D*) fhFractionOfSecInJet[ir]);
      }

      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhJetPtPartLevelVsJetPtDetLevelCorr_Rho%s", rhotype[ir].Data());
         fhJetPtPartLevelVsJetPtDetLevelCorr[ir] = new TH2D(name.Data(), name.Data(),270,-20,250,270,-20,250);
         fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr[ir]);
      } 

      fhJetPtPartLevelVsJetPtDetLevelZero = new TH2D("fhJetPtPartLevelVsJetPtDetLevelZero","fhJetPtPartLevelVsJetPtDetLevelZero",250,0,250,250,0,250);
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero);
 
      for(Int_t ir=0; ir<kRho; ir++){
         name = Form("fhJetPtResolutionVsPtPartLevel_Rho%s", rhotype[ir].Data());
         fhJetPtResolutionVsPtPartLevel[ir] = new TH2D(name.Data(), name.Data(),100,0,100,50,0,2);
         fOutput->Add((TH2D*) fhJetPtResolutionVsPtPartLevel[ir]);
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH   

         for(Int_t ir=0; ir<kRho; ir++){
            name = Form("fhJetPtPartLevelCorr_TTHdl%d_%d_Rho%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhJetPtPartLevelCorrTTHdl[itt][ir] = (TH1D*)  fhJetPtPartLevelCorr[ir]->Clone(name.Data());
            fOutput->Add((TH1D*) fhJetPtPartLevelCorrTTHdl[itt][ir]); //Norm spectrum for detector level TTH
         }

         for(Int_t ir=0; ir<kRho; ir++){
            name = Form("fhJetPtPartLevelVsJetPtDetLevelCorr_TTHdl%d_%d_Rho%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype[ir].Data());
            fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[itt][ir] = (TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr[ir]->Clone(name.Data());
            fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[itt][ir]); //ReMx for detector level TTH
         }
      }
   }

   //+++++++++++++++++++++++++++ EMBEDDING +++++++++++++++++++++++
   if(fMode == AliAnalysisTaskEA::kEmbedding){


      for(Int_t itg=kMB; itg<=kHM; itg++){   //@@@
         //remx normalization spectra
         for(Int_t ir =0; ir<kRho; ir++){
            name = Form("fhJetPtPartLevelCorr_EMB_%s_Rho%s",trig[itg].Data(),rhotype[ir].Data());
            fhJetPtPartLevelCorr_EMB[itg][ir] = new TH1D(name.Data(), name.Data(), 270, -20, 250);
            fOutput->Add((TH1D*) fhJetPtPartLevelCorr_EMB[itg][ir]);
         }
   
         name = Form("fhJetPtPartLevelZero_EMB_%s",trig[itg].Data());
         fhJetPtPartLevelZero_EMB[itg] = new TH1D(name.Data(), name.Data(), 250, 0, 250);
         fOutput->Add((TH1D*) fhJetPtPartLevelZero_EMB[itg]);


         for(Int_t itt=0; itt<fnHadronTTBins; itt++){   //@@@  
            //normalization for response matrix in events with TTH   
            for(Int_t ir=0; ir<kRho; ir++){   //@@@  
               name = Form("%s_TTHdl%d_%d", fhJetPtPartLevelCorr_EMB[itg][ir]->GetName(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
               fhJetPtPartLevelCorrTTHdl_EMB[itg][itt][ir] = (TH1D*) fhJetPtPartLevelCorr_EMB[itg][ir]->Clone(name.Data());    
               fOutput->Add((TH1D*) fhJetPtPartLevelCorrTTHdl_EMB[itg][itt][ir]);
            }

            name = Form("%s_TTHdl%d_%d",fhJetPtPartLevelZero_EMB[itg]->GetName(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhJetPtPartLevelZeroTTHdl_EMB[itg][itt] = (TH1D*) fhJetPtPartLevelZero_EMB[itg]->Clone(name.Data());    
            fOutput->Add((TH1D*) fhJetPtPartLevelZeroTTHdl_EMB[itg][itt]);
         }
      }

      //remx
      for(Int_t itg=kMB; itg<=kHM; itg++){   //@@@
         for(Int_t ir=0; ir<kRho; ir++){  
            name = Form("fhJetPtPartLevelVsJetPtDetLevelCorr_EMB_%s_Rho%s",trig[itg].Data(), rhotype[itg].Data());
            fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir] = new TH2D(name.Data(), name.Data(), 270, -20, 250, 270, -20, 250);
            fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir]);
         }

         for(Int_t ir=0; ir<kRho; ir++){  
            name = Form("fhJetPtPartLevelVsJetPtDetLevelZero_EMB_%s_Rho%s",trig[itg].Data(), rhotype[itg].Data());
            fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir] = new TH2D(name.Data(), name.Data(), 270, -20, 250, 250, 0, 250);        
            fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir]);
         }

         for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH   
            for(Int_t ir=0; ir<kRho; ir++){  
               name = Form("%s_TTHdl%d_%d", fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir]->GetName(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
               fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[itg][itt][ir] = (TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg][ir]->Clone(name.Data());
               fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[itg][itt][ir]);
            }

            for(Int_t ir=0; ir<kRho; ir++){  
               name = Form("%s_TTHdl%d_%d", fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir]->GetName(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
               fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[itg][itt][ir] = (TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg][ir]->Clone(name.Data());
               fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[itg][itt][ir]);
            }
         } 
      }

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhSharedJetFraction%s", trig[itg].Data());
         fhSharedJetFraction[itg] = new TH2D(name.Data(),name.Data(), 40,0,200, 20,0,2);
         fOutput->Add((TH2D*) fhSharedJetFraction[itg]); 
      }

 
      const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
      Int_t nPtHardBins = embeddingHelper->GetNPtHardBins();

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhTrialsEMBtot_%s", trig[itg].Data());
         fhTrialsEMBtot[itg] = new TH1F(name.Data(), name.Data(),  1, 0, 1);
         fOutput->Add((TH1F*) fhTrialsEMBtot[itg]); 
      }
 
      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhXsectionEMBtot_%s", trig[itg].Data());
         fhXsectionEMBtot[itg] = new TProfile(name.Data(), name.Data(),  1, 0, 1);
         fOutput->Add((TProfile*) fhXsectionEMBtot[itg]); 
      }


      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhTrialsEMB_%s", trig[itg].Data());
         fhTrialsEMB[itg] = new TH1F(name.Data(), name.Data(),  nPtHardBins, 0, nPtHardBins);
         fOutput->Add((TH1F*) fhTrialsEMB[itg]); 
      }
 
      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhXsectionEMB_%s", trig[itg].Data());
         fhXsectionEMB[itg] = new TProfile(name.Data(), name.Data(),  nPtHardBins, 0, nPtHardBins);
         fOutput->Add((TProfile*) fhXsectionEMB[itg]); 
      }

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhPtHardEMB_%s", trig[itg].Data());
         //fhPtHardEMB[itg] = new TH1F(name.Data(), name.Data(), fNbins*2, fMinBinPt, fMaxBinPt*4); 
         fhPtHardEMB[itg] = new TH1F(name.Data(), name.Data(), 1000, 0, 1000); 
         fOutput->Add((TProfile*) fhPtHardEMB[itg]); 
      }


   }

   //+++++++++++++++++++++++ MOEMENTU SMEARING HISTOGRAMS  +++++++++++++++++++++++++++++++++++++++++++++++

   fhOneOverPtVsPhiNeg = new TH2D("fhOneOverPtVsPhiNeg","1/pt versus track phi negative tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
   fOutput->Add((TH2D*) fhOneOverPtVsPhiNeg);
 
   fhOneOverPtVsPhiPos = new TH2D("fhOneOverPtVsPhiPos","1/pt versus track phi positive tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);
   fOutput->Add((TH2D*) fhOneOverPtVsPhiPos);
 
   fhSigmaPtOverPtVsPt = new TH2D("fhSigmaPtOverPtVsPt",
                                      "track sigma(1/pt)/ 1/pt vs pt", 100, 0, 100, 250, 0, 1);
   fOutput->Add((TH2D*) fhSigmaPtOverPtVsPt); 

   //+++++++++++++++++++++++ DCA HISTOGRAMS FOR SECONDARY TRACK CONTAMINATION ++++++++++++++++++++++++++++
   Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
   Int_t nbins = sizeof(bins)/sizeof(Double_t)-1; //pT binning for DCA distribution

   fhDCAinXVsPt = new TH2D("fhDCAinXVsPt","fhDCAinXVsPt", nbins, bins, 200, -10.,10);
   fOutput->Add((TH2D*) fhDCAinXVsPt); 

   fhDCAinYVsPt = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPt");
   fOutput->Add((TH2D*) fhDCAinYVsPt);

   if(fMode == AliAnalysisTaskEA::kMC){
      fhDCAinXVsPtPhysPrimary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinXVsPtPhysPrimary");
      fOutput->Add((TH2D*) fhDCAinXVsPtPhysPrimary);

      fhDCAinYVsPtPhysPrimary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPtPhysPrimary"); 
      fOutput->Add((TH2D*) fhDCAinYVsPtPhysPrimary); 
 
      fhDCAinXVsPtSecondary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinXVsPtSecondary");
      fOutput->Add((TH2D*) fhDCAinXVsPtSecondary); 

      fhDCAinYVsPtSecondary = (TH2D*) fhDCAinXVsPt->Clone("fhDCAinYVsPtSecondary");
      fOutput->Add((TH2D*) fhDCAinYVsPtSecondary); 
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

