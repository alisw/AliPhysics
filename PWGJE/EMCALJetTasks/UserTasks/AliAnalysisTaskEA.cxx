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
#include "AliFJWrapper.h"  //EMB_clus

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
fEmbeddPerpendicular(1),  // EMB_clus
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
//fCentralityV0A(-1),
//fCentralityV0C(-1),
fCentralityV0M(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
//fNTracklets(-1),
fMultV0A(0.),
fMultV0C(0.),
fMultV0M(0.),
fMultV0Mnorm(0.),
//fAsymV0M(999),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
//fAsymV0M_PartLevel(999),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0),
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhJetPtAreaV0norm_PartLevel(0x0),
fhRhoMBpart(0x0),
//fhV0MAssymVsV0Mnorm_PartLevel(0x0),
fhSignal_V0M_trueMB_PartLevel(0x0),
//fhV0A_V0C_PartLevel(0x0),
//fhV0A_V0APartLevel(0x0),
//fhV0C_V0CPartLevel(0x0),
//fhV0MvsV0Mnorm(0x0),
//fhV0AvsSPD(0x0),
//fhV0CvsSPD(0x0),
fFastJetWrapper(0x0), // EMB_clus
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),                         //1D unfolding
fhJetPtZeroPartLevelVsJetPtZeroDetLevel(0x0),                     //1D unfolding
fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr(0x0),                   //1D unfolding (added by KA)
fhImpurityInclusive_DetJetPtVsPartJetPtCorr(0x0),                 //Impurity distribution (added by KA)
fhPurity_InclusiveJets(nullptr),                                  //Purity estimation
fhNormalization_Purity(nullptr),                                  //Purity estimation
fhJetEnergy_Angle_Resol(nullptr),                                 //Jet energy and angular resolution
fhInclusiveRespMatrix_HM(nullptr),                                //Check unfolding, HM events
fhMissedEvents_HM(nullptr),                                       //Check unfolding, HM events
fhPhi_JetPtPartLevel_InclusiveJets(0x0),                          //2D unfolding
fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets(0x0),     //2D unfolding
fhPhi_JetPtZeroPartLevel_InclusiveJets(0x0),                      //2D unfolding (added by KA)
fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
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
fhFractionOfSecInJet(0x0),
//fhV0ARunByRunMB(0x0),
//fhV0CRunByRunMB(0x0),
//fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
fhJetPtEvtByEvent(0x0),
//fhRecoilJetPtEvtByEventRandomTT(0x0),
//fhNumberOfHighPtJetsRecoilRandomTTPartLevel(0x0),
fhJetPtEvtByEventPartLevel(0x0),
//fhRecoilJetPtEvtByEventRandomTTPartLevel(0x0),
fhNjetsPt10_15_PartLevel(0x0),
fhNjetsPtAbove15_PartLevel(0x0),
fhNjetsPt10_15_PartLevelV0MnormDetLev(0x0),
fhNjetsPtAbove15_PartLevelV0MnormDetLev(0x0),
fPhiSmearingTrack(0x0),
fPhiSmearingJet(0x0),
fDeltaPhiSmearingTT2030Jet(0x0),
fhTrackEtaInclEMB(0x0),
fMinFractionShared(0),
fZVertexCut(10.0),
fnHadronTTBins(0),
//fnJetChTTBins(0),
//fnClusterTTBins(0),
fMode(AliAnalysisTaskEA::kNormal),
//fFiducialCellCut(0x0),
fHelperEA(0x0),
fMeanV0M(1.),
fMeanV0M_PartLevel(1.),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.),
//fRhoType(0),
kOldV0MC(kFALSE),
fMultFramework(kFALSE),
fRho(0.),
fRhoMC(0.),
fRhoEMB(0.),
fRunnumber(0),
fMaxFacPtHard(0),
fhNjetReMx_V0MnormDetLev_15GeV(0x0),
fhNjetNorm_V0MnormDetLev_15GeV(0x0),
fhNjetReMx_V0MnormDetLev_20GeV(0x0),
fhNjetNorm_V0MnormDetLev_20GeV(0x0)



{
   //default constructor
   //1D respnse matrix from recoil jets  //FF
   for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
      fhRecoilJetPtPartLevelCorr[itt] = NULL;
      fhRecoilJetPtZeroPartLevel[itt] = NULL;

      fhRecoilJetPtPartLevel_CorrespTT[itt] = NULL;                        // Modified by KA
      fhRecoilJetPtZeroPartLevel_CorrespTT[itt] = NULL;                        // Modified by KA

      fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt] = NULL;                   // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt] = NULL;               // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt] = NULL;           // Modified by KA

      fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt] = NULL;         // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt] = NULL;     // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt] = NULL; // Modified by KA

      fhImpurityRecoil_DetJetPtVsPartJetPtCorr[itt] = NULL;
   }
   //2D unfolding
   for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
      fhDeltaPhi_JetPtPartLevel[itt] = NULL;                                             //2D unfolding
      fhDeltaPhi_JetPtPartLevel_CorrespTT[itt] = NULL;                                   //2D unfolding
      fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL;                               //2D unfolding
      fhDeltaPhi_JetPtZero_PartLevel[itt] = NULL;                                        //2D unfolding (added by KA)
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt] = NULL;                   //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt] = NULL;         //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL;     //2D unfolding
      fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL; //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = NULL;               //2D unfolding (added by KA)
      fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = NULL;           //2D unfolding (added by KA)
      // fhNumberOf_ChosenTT_PartLevel[itt] = NULL;
   }

   for(Int_t iq = 0; iq < 4; ++iq){
      fArray_for_filling[iq] = 0.;
   }
   /////////////////


   for(Int_t i=0; i<fkTTbins; i++){
      fhJetPtAreaV0normTTH_PartLevel[i] = 0x0;

      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhJetPtAreaV0normTTH[itg][i] = 0x0;
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      //fhV0A_V0C_V0Mnorm[itg]=0x0;

      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;

      fhTrackMult[itg]=0x0;

      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      //fhClusterPhiIncl[itg] = 0x0;
      //fhClusterEtaIncl[itg] = 0x0;

      fhTrackPtEtaPhiV0norm[itg] = 0x0;
      fhJetPtEtaPhiV0norm[itg] = 0x0;
      for(Int_t i=0; i<fkTTbins; i++){
         fhJetPtEtaPhiV0normTTH[itg][i] = 0x0;
      }

      fhJetPtAreaV0norm[itg] = 0x0;
      fhRho[itg] = 0x0;

      for(Int_t i=0; i<fkTTbins; i++){
         fhRhoTTH[itg][i]=0x0;
         //fhRhoTTC[itg][i]=0x0;
         //fhRhoTTJ[itg][i]=0x0;
      }

      fhSharedJetFraction[itg] = 0x0;
      fhTrialsEMBtot[itg] = 0x0;
      fhXsectionEMBtot[itg] = 0x0;
      //fhTrialsEMB[itg] = 0x0;
      //fhXsectionEMB[itg] = 0x0;
      fhPtHardEMB[itg] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i]   = 0;
      fHadronTT_Labels[i].resize(0); //Modified by KA
      //fJetChTT[i]    = 0;
      //fClusterTT[i]  = 0;


      fHadronTT_PartLevel[i]   = 0;
      //fClusterTT_PartLevel[i]   = 0;

      //TT
      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhMultTTH[itg][i] = 0x0;
         //fhMultTTJ[itg][i] = 0x0;
         //fhMultTTC[itg][i] = 0x0;

         //fhTTH_CentV0M[itg][i]  = 0x0;
         fhTTH_V0Mnorm1[itg][i] = 0x0;
         //fhTTH_3D_V0Mnorm1[itg][i] = 0x0;

         //fhTTC_CentV0M[itg][i]  = 0x0;
         //fhTTC_V0Mnorm1[itg][i] = 0x0;

         //fhV0MAssymVsV0MnormTTH[itg][i] = 0x0;
      }

      fhTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP
      fhTT_Corresp[i] = 0x0; //KA
      //fhTTH_3D_V0Mnorm1_PartLevel[i] = 0x0;

      //fhTTC_V0Mnorm1_PartLevel[i] = 0x0;

      //fhV0MAssymVsV0MnormTTH_PartLevel[i] = 0x0;

      //RECOIL JET SPECTRA
      for(Int_t itg=kMB; itg<=kGA; itg++){
         //fhRecoilJetPtTTH_CentV0M[itg][i]  = 0x0;
         fhRecoilJetPtTTH_V0Mnorm1[itg][i] = 0x0;

         fhRecoilJetPhiTTH_V0Mnorm1[itg][i] = 0x0;
         //fhRecoilJetTTH_V0Mnorm1[itg][i]    = 0x0;

         //fhRecoilJetPtTTC_CentV0M[itg][i]  = 0x0;
         //fhRecoilJetPtTTC_V0Mnorm1[itg][i] = 0x0;
      }

      fhRecoilJetPtTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_TTH_V0MnormDet_PartLevel[i] = NULL; //added by KA
      //fhRecoilJetPtTTC_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhRecoilJetPtTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP

      fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm1_PartLevel[i] = NULL; //added by KA
      //fhRecoilJetTTH_V0Mnorm1_PartLevel[i]    = 0x0;
      fhRecoilJetPhiTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhRecoilJetPhiTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP


      for(Int_t itg=kMB; itg<=kGA; itg++){
         //fhDeltaPtTTH_RC_CentV0M[itg][i] = 0x0;
         //fhDeltaPtTTC_RC_CentV0M[itg][i] = 0x0;

         fhDeltaPtTTH_RC_V0Mnorm1[itg][i] = 0x0;
         //fhDeltaPtTTC_RC_V0Mnorm1[itg][i] = 0x0;

         fhDeltaPtEmbeddPerpendicular[itg][i] = 0x0; // EMB_clus
      }

      fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[i] = 0x0;
      //fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[i] = 0x0;



      //embedding
      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][i] = 0x0;
         //fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][i] = 0x0;
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg] = 0x0;
      fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg] = 0x0;
   }


   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhJetPtPartLevelCorr_EMB[itg] = 0x0;
      fhJetPtPartLevelZero_EMB[itg] = 0x0;
   }

   //for(Int_t itg=kMB; itg<=kGA; itg++){
   //   for(Int_t is=0; is<fkShift; is++){
   //      fhRecoilJetPtTTHref_V0Mnorm1_rhoShift[itg][is] = 0x0;
   //   }
   //}


   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
      //fJetChTTLowPt[i]=-1;
      //fJetChTTHighPt[i]=-1;
      //fClusterTTLowPt[i]=-1;
      //fClusterTTHighPt[i]=-1;

      //fhV0AvsV0CTTH[i] = 0x0;
      //fhV0AvsV0CTTJ[i] = 0x0;
      //fhV0AvsV0CTTCinMB[i] = 0x0;
      //fhV0AvsV0CTTCinGA[i] = 0x0;
   }

   for(Int_t iv=0; iv<fkVtx;iv++){
      fhVertex[iv]=0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhCentrality[itg]= 0x0;
      for(Int_t ic=0; ic<fkCE;ic++){

         fhSignal[itg][ic] = 0x0;

         for(Int_t i=0; i<fkTTbins;i++){
            //fhCentralityTTH[itg][ic][i] = 0x0;
            //fhCentralityTTJ[itg][ic][i] = 0x0;
            //fhCentralityTTC[itg][ic][i] = 0x0;

            fhSignalTTH[itg][ic][i] = 0x0;
            //fhSignalTTJ[itg][ic][i] = 0x0;
            //fhSignalTTC[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0;

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
         //fhSignalTTC_PartLevel[ic][i] = 0x0;
      }
   }

   //for(Int_t itg=kMB; itg<=kGA; itg++){
   //   fhV0MAssymVsV0Mnorm[itg] = 0x0;
   //}

   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMBpart[i]=0x0;
      //fhRhoTTCinMBpart[i]=0x0;
   }

   //fFiducialCellCut = new AliEMCALRecoUtils();

   for(Int_t i=0; i<fkTTbins; i++){
      //fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      //fIndexTTJ[i] = -1;

      fdeltapT[i]  = 0.;
      fdeltapT_PartLevel[i]  = 0.;

      fIndexTTH_PartLevel[i] = -1;
      //fIndexTTC_PartLevel[i] = -1;

      //fTTC[i].resize(0);
      fTTH[i].resize(0);
      //fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      //fTTC_PartLevel[i].resize(0);
   }

   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }



   //JET AND TRACK PT ASYMMETRY
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      for(Int_t itg=kMB; itg<=kHM; itg++){
         //fhJetPtAsymmetryCB[itg][itt]   = NULL;
         //fhTrackPtAsymmetryCB[itg][itt] = NULL;

//         fhNumberOfHighPtJetsCB[itg][itt]     = NULL;
         fhNumberOfHighPtJetsRecoil[itg][itt] = NULL;
//         if(itt==0) fhNumberOfHighPtJetsRecoilRandomTT[itg] = NULL;
      }
      fhRecoilJetPtEvtByEvent[itt] = NULL;

      //fhJetPtAsymmetryCBPartLevel[itt] = NULL;
      //fhTrackPtAsymmetryCBPartLevel[itt] = NULL;
      //fhNumberOfHighPtJetsCBPartLevel[itt] = NULL;
      fhNumberOfHighPtJetsRecoilPartLevel[itt] = NULL;
      fhNumberOfHighPtJetsRecoilPartLevelV0MnormDetLev[itt] = NULL; //FILIP
      fhRecoilJetPtEvtByEventPartLevel[itt] = NULL;

   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhNjetsPt10_15[itg] = NULL;             
      fhNjetsPtAbove15[itg] = NULL;           
   } 

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){  //fk
      fhNotMatchedJetPt[itt] = NULL;
   }

   fHelperEA = new PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA();
   fMeanV0M_PartLevel = fHelperEA->GetV0MPartLevel();
}

//________________________________________________________________________
AliAnalysisTaskEA::AliAnalysisTaskEA(const char *name):
AliAnalysisTaskEmcalJet(name,kTRUE),
fEmbeddPerpendicular(1),  // EMB_clus
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
//fCentralityV0A(-1),
//fCentralityV0C(-1),
fCentralityV0M(-1),
fxVertex(-1),
fyVertex(-1),
fzVertex(-1),
//fNTracklets(-1),
fMultV0A(0.),
fMultV0C(0.),
fMultV0M(0.),
fMultV0Mnorm(0.),
//fAsymV0M(999),
fMultV0A_PartLevel(0.),
fMultV0C_PartLevel(0.),
fMultV0M_PartLevel(0.),
fMultV0Mnorm_PartLevel(0.),
//fAsymV0M_PartLevel(999),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0),
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhJetPtAreaV0norm_PartLevel(0x0),
fhRhoMBpart(0x0),
//fhV0MAssymVsV0Mnorm_PartLevel(0x0),
fhSignal_V0M_trueMB_PartLevel(0x0),
//fhV0A_V0C_PartLevel(0x0),
//fhV0A_V0APartLevel(0x0),
//fhV0C_V0CPartLevel(0x0),
//fhV0MvsV0Mnorm(0x0),
//fhV0AvsSPD(0x0),
//fhV0CvsSPD(0x0),
fFastJetWrapper(0x0),  // EMB_clus
fhPtTrkTruePrimGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkSecOrFakeRec(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),
fhJetPtZeroPartLevelVsJetPtZeroDetLevel(0x0),
fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr(0x0),                   //1D unfolding (added by KA)
fhImpurityInclusive_DetJetPtVsPartJetPtCorr(0x0),                 //Impurity distribution (added by KA)
fhPurity_InclusiveJets(nullptr),                                  //Purity estimation
fhNormalization_Purity(nullptr),                                  //Purity estimation
fhJetEnergy_Angle_Resol(nullptr),                                 //Jet energy and angular resolution
fhInclusiveRespMatrix_HM(nullptr),                                //Check unfolding, HM events
fhMissedEvents_HM(nullptr),                                       //Check unfolding, HM events
fhPhi_JetPtPartLevel_InclusiveJets(0x0),                          //2D unfolding
fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets(0x0),     //2D unfolding
fhPhi_JetPtZeroPartLevel_InclusiveJets(0x0),                      //2D unfolding (added by KA)
fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
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
fhFractionOfSecInJet(0x0),
//fhV0ARunByRunMB(0x0),
//fhV0CRunByRunMB(0x0),
//fhV0MRunByRunMB(0x0),
fhV0MnormRunByRunMB(0x0),
fhJetPtEvtByEvent(0x0),
//fhRecoilJetPtEvtByEventRandomTT(0x0),
//fhNumberOfHighPtJetsRecoilRandomTTPartLevel(0x0),
fhJetPtEvtByEventPartLevel(0x0),
//fhRecoilJetPtEvtByEventRandomTTPartLevel(0x0),
fhNjetsPt10_15_PartLevel(0x0), 
fhNjetsPtAbove15_PartLevel(0x0),
fhNjetsPt10_15_PartLevelV0MnormDetLev(0x0),
fhNjetsPtAbove15_PartLevelV0MnormDetLev(0x0),
fPhiSmearingTrack(0x0),
fPhiSmearingJet(0x0),
fDeltaPhiSmearingTT2030Jet(0x0),
fhTrackEtaInclEMB(0x0),
fMinFractionShared(0),
fZVertexCut(10.0),
fnHadronTTBins(0),
//fnJetChTTBins(0),
//fnClusterTTBins(0),
fMode(AliAnalysisTaskEA::kNormal),
//fFiducialCellCut(0x0),
fHelperEA(0x0),
fMeanV0M(1.),
fMeanV0M_PartLevel(1.),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.),
//fRhoType(0),
kOldV0MC(kFALSE),
fMultFramework(kFALSE),
fRho(0.),
fRhoMC(0.),
fRhoEMB(0.),
fRunnumber(0),
fMaxFacPtHard(0),
fhNjetReMx_V0MnormDetLev_15GeV(0x0),
fhNjetNorm_V0MnormDetLev_15GeV(0x0),
fhNjetReMx_V0MnormDetLev_20GeV(0x0),
fhNjetNorm_V0MnormDetLev_20GeV(0x0)
{
   //Constructor

   //1D respnse matrix from recoil jets  //FF
   for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
      fhRecoilJetPtPartLevelCorr[itt] = NULL;
      fhRecoilJetPtZeroPartLevel[itt] = NULL;

      fhRecoilJetPtPartLevel_CorrespTT[itt] = NULL; // Modified by KA
      fhRecoilJetPtZeroPartLevel_CorrespTT[itt] = NULL; // Modified by KA

      fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt] = NULL;
      fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt] = NULL;
      fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt] = NULL;

      fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt] = NULL; // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt] = NULL; // Modified by KA
      fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt] = NULL; // Modified by KA

      fhImpurityRecoil_DetJetPtVsPartJetPtCorr[itt] = NULL;
   }

   //2D unfolding
   for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
      fhDeltaPhi_JetPtPartLevel[itt] = NULL;                                             //2D unfolding
      fhDeltaPhi_JetPtPartLevel_CorrespTT[itt] = NULL;                                   //2D unfolding
      fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL;                               //2D unfolding
      fhDeltaPhi_JetPtZero_PartLevel[itt] = NULL;                                        //2D unfolding (added by KA)
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt] = NULL;                   //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt] = NULL;         //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL;     //2D unfolding
      fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = NULL; //2D unfolding
      fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = NULL;               //2D unfolding (added by KA)
      fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = NULL;           //2D unfolding (added by KA)
      // fhNumberOf_ChosenTT_PartLevel[itt] = NULL;
   }

   for(Int_t iq = 0; iq < 4; ++iq){
      fArray_for_filling[iq] = 0.;
   }
   /////////////////


   for(Int_t i=0; i<fkTTbins; i++){
      fhJetPtAreaV0normTTH_PartLevel[i] = 0x0;

      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhJetPtAreaV0normTTH[itg][i] = 0x0;
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      //fhV0A_V0C_V0Mnorm[itg]=0x0;

      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;

      fhTrackMult[itg]=0x0;

      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      fhTrackPtEtaPhiV0norm[itg] = 0x0;
      fhJetPtEtaPhiV0norm[itg] = 0x0;
      for(Int_t i=0; i<fkTTbins; i++){
         fhJetPtEtaPhiV0normTTH[itg][i] = 0x0;
      }

      //fhClusterPhiIncl[itg] = 0x0;
      //fhClusterEtaIncl[itg] = 0x0;

      fhJetPtAreaV0norm[itg] = 0x0;
      fhRho[itg] = 0x0;

      for(Int_t i=0; i<fkTTbins; i++){
         fhRhoTTH[itg][i]=0x0;
         //fhRhoTTC[itg][i]=0x0;
         //fhRhoTTJ[itg][i]=0x0;
      }

      fhSharedJetFraction[itg] = 0x0;
      fhTrialsEMBtot[itg] = 0x0;
      fhXsectionEMBtot[itg] = 0x0;
      //fhTrialsEMB[itg] = 0x0;
      //fhXsectionEMB[itg] = 0x0;
      fhPtHardEMB[itg] = 0x0;
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT[i]   = 0;
      fHadronTT_Labels[i].resize(0); // Modified by KA
      //fJetChTT[i]    = 0;
      //fClusterTT[i]  = 0;


      fHadronTT_PartLevel[i]   = 0;
      fHadronTT_Labels_PartLevel[i].resize(0); // Modified by KA
      //fClusterTT_PartLevel[i]   = 0;

      //TT
      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhMultTTH[itg][i] = 0x0;
         //fhMultTTJ[itg][i] = 0x0;
         //fhMultTTC[itg][i] = 0x0;

         //fhTTH_CentV0M[itg][i]  = 0x0;
         fhTTH_V0Mnorm1[itg][i] = 0x0;
         //fhTTH_3D_V0Mnorm1[itg][i] = 0x0;

         //fhTTC_CentV0M[itg][i]  = 0x0;
         //fhTTC_V0Mnorm1[itg][i] = 0x0;

         //fhV0MAssymVsV0MnormTTH[itg][i] = 0x0;
      }

      fhTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP
      fhTT_Corresp[i] = 0x0; //KA
      //fhTTH_3D_V0Mnorm1_PartLevel[i] = 0x0;
      //fhTTC_V0Mnorm1_PartLevel[i] = 0x0;

      //fhV0MAssymVsV0MnormTTH_PartLevel[i] = 0x0;

      //RECOIL JET SPECTRA
      for(Int_t itg=kMB; itg<=kGA; itg++){
         //fhRecoilJetPtTTH_CentV0M[itg][i]  = 0x0;
         fhRecoilJetPtTTH_V0Mnorm1[itg][i] = 0x0;

         fhRecoilJetPhiTTH_V0Mnorm1[itg][i] = 0x0;
         //fhRecoilJetTTH_V0Mnorm1[itg][i]    = 0x0;

         //fhRecoilJetPtTTC_CentV0M[itg][i]  = 0x0;
         //fhRecoilJetPtTTC_V0Mnorm1[itg][i] = 0x0;
      }

      fhRecoilJetPtTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_TTH_V0MnormDet_PartLevel[i] = NULL; //added by KA
      //fhRecoilJetPtTTC_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhRecoilJetPtTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP

      fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm1_PartLevel[i] = NULL; //added by KA
      //fhRecoilJetTTH_V0Mnorm1_PartLevel[i]    = 0x0;
      fhRecoilJetPhiTTH_V0Mnorm_AllMB_PartLevel[i] = 0x0; //FF
      fhRecoilJetPhiTTH_V0MnormDetLev_AllMB_PartLevel[i] = 0x0; //FILIP

      for(Int_t itg=kMB; itg<=kGA; itg++){
         //fhDeltaPtTTH_RC_CentV0M[itg][i] = 0x0;
         //fhDeltaPtTTC_RC_CentV0M[itg][i] = 0x0;

         fhDeltaPtTTH_RC_V0Mnorm1[itg][i] = 0x0;
         //fhDeltaPtTTC_RC_V0Mnorm1[itg][i] = 0x0;

	 fhDeltaPtEmbeddPerpendicular[itg][i] = 0x0; // EMB_clus
      }

      fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[i] = 0x0;
      //fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[i] = 0x0;



      //embedding
      for(Int_t itg=kMB; itg<=kGA; itg++){
         fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][i] = 0x0;
         //fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][i] = 0x0;
      }
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg] = 0x0;
      fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg] = 0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhJetPtPartLevelCorr_EMB[itg] = 0x0;

      fhJetPtPartLevelZero_EMB[itg] = 0x0;
   }

   //for(Int_t itg=kMB; itg<=kGA; itg++){
   //   for(Int_t is=0; is<fkShift; is++){
   //      fhRecoilJetPtTTHref_V0Mnorm1_rhoShift[itg][is] = 0x0;
   //   }
   //}

   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
      //fJetChTTLowPt[i]=-1;
      //fJetChTTHighPt[i]=-1;
      //fClusterTTLowPt[i]=-1;
      //fClusterTTHighPt[i]=-1;

      //fhV0AvsV0CTTH[i] = 0x0;
      //fhV0AvsV0CTTJ[i] = 0x0;
      //fhV0AvsV0CTTCinMB[i] = 0x0;
      //fhV0AvsV0CTTCinGA[i] = 0x0;
   }

   for(Int_t iv=0; iv<fkVtx;iv++){
      fhVertex[iv]=0x0;
   }

   for(Int_t itg=kMB; itg<=kGA; itg++){
      fhCentrality[itg] = 0x0;
      for(Int_t ic=0; ic<fkCE;ic++){

         fhSignal[itg][ic] = 0x0;

         for(Int_t i=0; i<fkTTbins;i++){
            //fhCentralityTTH[itg][ic][i] = 0x0;
            //fhCentralityTTJ[itg][ic][i] = 0x0;
            //fhCentralityTTC[itg][ic][i] = 0x0;

            fhSignalTTH[itg][ic][i] = 0x0;
            //fhSignalTTJ[itg][ic][i] = 0x0;
            //fhSignalTTC[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0;

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
         //fhSignalTTC_PartLevel[ic][i] = 0x0;
      }
   }

   //for(Int_t itg=kMB; itg<=kGA; itg++){
   //   fhV0MAssymVsV0Mnorm[itg] = 0x0;
   //}


   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMBpart[i]=0x0;
      //fhRhoTTCinMBpart[i]=0x0;
   }

   //fFiducialCellCut = new AliEMCALRecoUtils();

   for(Int_t i=0; i<fkTTbins; i++){
      //fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      //fIndexTTJ[i] = -1;

      fdeltapT[i]  = 0.;
      fdeltapT_PartLevel[i]  = 0.;

      fIndexTTH_PartLevel[i] = -1;
      //fIndexTTC_PartLevel[i] = -1;

      //fTTC[i].resize(0);
      fTTH[i].resize(0);
      //fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      //fTTC_PartLevel[i].resize(0);
   }

   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }

   //JET AND TRACK PT ASYMMETRY
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      for(Int_t itg=kMB; itg<=kHM; itg++){
         //fhJetPtAsymmetryCB[itg][itt]   = NULL;
         //fhTrackPtAsymmetryCB[itg][itt] = NULL;

//         fhNumberOfHighPtJetsCB[itg][itt] = NULL;
         fhNumberOfHighPtJetsRecoil[itg][itt] = NULL;
//         if(itt==0) fhNumberOfHighPtJetsRecoilRandomTT[itg] = NULL;
      }

      fhRecoilJetPtEvtByEvent[itt] = NULL;

      //fhJetPtAsymmetryCBPartLevel[itt] = NULL;
      //fhTrackPtAsymmetryCBPartLevel[itt] = NULL;
      //fhNumberOfHighPtJetsCBPartLevel[itt] = NULL;
      fhNumberOfHighPtJetsRecoilPartLevel[itt] = NULL;
      fhNumberOfHighPtJetsRecoilPartLevelV0MnormDetLev[itt] = NULL; //FILIP
      fhRecoilJetPtEvtByEventPartLevel[itt] = NULL;

   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhNjetsPt10_15[itg] = NULL;             
      fhNjetsPtAbove15[itg] = NULL;           
   } 

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){  //fk
      fhNotMatchedJetPt[itt] = NULL;
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

   if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding || mode == AliAnalysisTaskEA::kKine){  //TO BE CHECKED FOR EMBEDDING
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

   if(mode != AliAnalysisTaskEA::kKine){
      trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks (or combined tracks if embedding)
      trackCont->SetMinPt(0.15);
      trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
   }

   if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding){
      trackContTrue = task->AddMCParticleContainer(mcpariclearraynamePartMC); //particle level MC particles
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(1e-3); // KA: old value was 0.15 GeV/c
      trackContTrue->SetEtaLimits(-5.1,5.1); //V0 eta range

      if(mode == AliAnalysisTaskEA::kEmbedding) trackContTrue->SetIsEmbedding(kTRUE);
   }

   if(mode == AliAnalysisTaskEA::kKine){
      trackContTrue = task->AddParticleContainer(mcpariclearraynamePartMC); //particle level MC particles
      trackContTrue->SetMinPt(1e-3); // KA: old value was 0.15 GeV/c
      trackContTrue->SetEtaLimits(-5.1,5.1); //V0 eta range
   }


   if(mode == AliAnalysisTaskEA::kEmbedding){
      trackContDet = task->AddTrackContainer(tracknameDetMC);  //detector level pythia tracks when embedding
      trackContDet->SetMinPt(0.15);
      trackContDet->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
      trackContDet->SetIsEmbedding(kTRUE);
   }


   if(mode != AliAnalysisTaskEA::kKine){
      clusterCont = task->AddClusterContainer(clusterarrayname);  //detector level tracks (needs to be checked for embedding)
      clusterCont->SetMinPt(0.3);
      clusterCont->SetExoticCut(1);
      clusterCont->SetClusTimeCut(0, emcaltofcut);
   }
   //   clusterCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);

   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //AKT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrue   = 0x0; //AKT jet container with mc particle level jets pythia
   AliJetContainer *jetContDet    = 0x0; //AKT jet container used when embedding with mc jets at detector level pythia

   AliJetContainer *jetContRecKT  = 0x0; //KT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrueKT = 0x0; //KT jet container with mc particle level jets pythia
   AliJetContainer *jetContDetKT  = 0x0; //KT jet container used when embedding with mc jets at detector level pythia

   if(mode != AliAnalysisTaskEA::kKine){
      //AKT DETECTOR LEVEL JET    (or combined event jet container when embedding)
      jetContRec   = task->AddJetContainer(jetarrayname,"TPC",jetRadius);

      if(jetContRec){
         jetContRec->ConnectParticleContainer(trackCont);
         jetContRec->SetPercAreaCut(acut);
         jetContRec->SetMinPt(0.150);
         jetContRec->SetMaxTrackPt(100.0);
         jetContRec->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT DETECTOR LEVEL JET    (or combined event jet container when embedding)
      jetContRecKT   = task->AddJetContainer(ktjetarrayname,"TPC",jetRadiuskt);

      if(jetContRecKT){
         jetContRecKT->ConnectParticleContainer(trackCont);
         //jetContRecKT->SetPercAreaCut(acut);
         jetContRecKT->SetMinPt(0.);
         jetContRecKT->SetMaxTrackPt(100.0);
         jetContRecKT->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContRecKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
      }
    }

    if(mode == AliAnalysisTaskEA::kMC || mode == AliAnalysisTaskEA::kEmbedding || mode == AliAnalysisTaskEA::kKine || mode == AliAnalysisTaskEA::kEmbPy){
      //AKT JETS PARTICLE LEVEL
      jetContTrue = task->AddJetContainer(jetarraynamePartMC,"TPC",jetRadius);

      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);
         jetContTrue->SetMinPt(1e-3); // KA: old value is 0.15 GeV/c
         jetContTrue->SetMaxTrackPt(1000.0); // added by K.A.
         jetContTrue->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT JETS PARTICLE LEVEL
      jetContTrueKT = task->AddJetContainer(ktjetarraynamePartMC,"TPC",jetRadiuskt);

      if(jetContTrueKT){
         jetContTrueKT->ConnectParticleContainer(trackContTrue);
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetMaxTrackPt(1000.0); // added by K.A.
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
         jetContDet->SetMaxTrackPt(100.0);
         jetContDet->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContDet->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT DETECTOR LEVEL JET    (or combined event jet container when embedding)
      jetContDetKT = task->AddJetContainer(ktjetarraynameDetMC,"TPC",jetRadiuskt);

      if(jetContDetKT) {
         //jetContDetKT->SetPercAreaCut(acut);
         jetContDetKT->SetMinPt(0.0);
         jetContDetKT->SetMaxTrackPt(100.0);
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

  if(fMode == AliAnalysisTaskEA::kMC)   return kFALSE; //MC
  if(fMode == AliAnalysisTaskEA::kKine) return kFALSE; //MC

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
     /*AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
     if(vzeroAOD){
        if(vzeroAOD->GetMTotV0A() > 0 && vzeroAOD->GetMTotV0C() > 0)  passedTrigger = kTRUE;
     }*/
     passedTrigger = kTRUE;
  }else if(fMode == AliAnalysisTaskEA::kKine){
     passedTrigger = kTRUE;
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

   if(fMode == AliAnalysisTaskEA::kMC)   return kFALSE; //MC
   if(fMode == AliAnalysisTaskEA::kKine) return kFALSE; //MC

   bool passedTrigger = kFALSE;
   UInt_t triggerMask = fInputHandler->IsEventSelected();
   if(triggerMask & AliVEvent::kHighMultV0){
      passedTrigger = kTRUE;
   }

   return passedTrigger;
}
//_____________________________________________________________________________________
Double_t AliAnalysisTaskEA::GetMyRho(AliJetContainer* ktjets){

   // Get rho from event
   Double_t myrho = 0.;


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

   //if(fRhoType == krhokt){ //KT BACKGROUND
      Int_t nJetAcckt = 0;

      for(auto jetIterator : ktjets->accepted_momentum() ){
                      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
          jet = jetIterator.second;  // Get the pointer to jet object
          if(!jet)  continue;

          if(jet==jetLJ) continue; //skip two leading kT jets
          if(jet==jetSJ) continue;

          //standard area based approach
          frhovec[nJetAcckt]  = jet->Pt()/jet->Area();
          nJetAcckt++;
      }

      if(nJetAcckt>0){
         myrho = TMath::Median(nJetAcckt, frhovec);
      }
// }else{ //CMS MODIFICATION
//    Int_t nJetAccms = 0;
//    Double_t areaPhysJets = 0.0;
//    Double_t areaAllJets  = 0.0;
//
//    for(auto jetIterator : ktjets->accepted_momentum() ){
//                    // trackIterator is a std::map of AliTLorentzVector and AliVTrack
//        jet = jetIterator.second;  // Get the pointer to jet object
//        if(!jet)  continue;
//
//        if(jet==jetLJ) continue; //skip two leading kT jets
//        if(jet==jetSJ) continue;
//
//        //cms modification
//        areaAllJets += jet->Area();
//
//        if(jet->Pt() > 0.1){
//           areaPhysJets     += jet->Area();
//           frhovec[nJetAccms]  = jet->Pt()/jet->Area();
//           nJetAccms++;
//        }
//    }
//
//    if(nJetAccms>0){
//       myrho = TMath::Median(nJetAccms, frhovec)*(areaPhysJets/areaAllJets);
//    }
// }

   return  myrho;
}
//________________________________________________________________________

Bool_t AliAnalysisTaskEA::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION RECONSTRUCTED DATA

   if(!event) return kFALSE;
   if(fMode == AliAnalysisTaskEA::kKine)   return kTRUE;

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

      if(event->IsPileupFromSPDInMultBins()){
         fHistEvtSelection->Fill(2.5); //count events rejected by pileup
         return kFALSE;
      }
   }
   //___________________________________________________
   //MULTIPLICITY SELECTIONS

   if(fMode == AliAnalysisTaskEA::kNormal){
      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection ||
         //!fMultSelection->IsEventSelected() ||
         !fMultSelection->GetThisEventIsNotPileup() ||
         !fMultSelection->GetThisEventIsNotPileupInMultBins() ||
	 !fMultSelection->GetThisEventHasNoInconsistentVertices() ||
	 !fMultSelection->GetThisEventPassesTrackletVsCluster()){

         fHistEvtSelection->Fill(7.5); //count events rejected by multiplicity selection
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

   if(fEmbeddPerpendicular){
      fFastJetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper"); //EMB_clus
      fFastJetWrapper->SetAreaType(fastjet::active_area);          //EMB_clus
      fFastJetWrapper->SetGhostArea(0.005);
      fFastJetWrapper->SetR(0.4);                                  //EMB_clus
      fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);    //EMB_clus
      fFastJetWrapper->SetRecombScheme(fastjet::pt_scheme);        //EMB_clus
   }

   return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEA::IsOutlier(){ //FK// whole function
   //Checks that this event is pthard bin outlier
   //inspired by Bool_t AliConvEventCuts::IsJetJetMCEventAccepted

   if(TMath::Abs(fMaxFacPtHard) < 1e-6) return kFALSE; //FK// skip

   TList *genHeaders         = 0x0;
   AliGenEventHeader* gh     = 0;
   Float_t ptHard;
   AliEmcalJet* jetMC = 0x0;

   Bool_t bPythiaHeader = 0; // flag whether pythia header was found

   if(MCEvent()){
      genHeaders = MCEvent()->GetCocktailList(); //get list of MC cocktail headers
   }

   if(genHeaders){
      for(Int_t i = 0; i<genHeaders->GetEntries(); i++){
         gh = (AliGenEventHeader*)genHeaders->At(i);

         AliGenPythiaEventHeader* pyhead= dynamic_cast<AliGenPythiaEventHeader*>(gh); //identify pythia header

         if(pyhead){
            bPythiaHeader = 1;
            ptHard = pyhead->GetPtHard();

            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jetMC = jetIterator.second;  // Get the pointer to jet object
               if(!jetMC)  continue;

               //Compare jet pT and pt Hard
               if(jetMC->Pt() > fMaxFacPtHard * ptHard){
                  fHistEvtSelection->Fill(9.5); // I think idea of that was the following: to count outliers (to show total number of rejected events)
                  return kTRUE;
               }
            }
         }
      }

      if(!bPythiaHeader){ //ptyhia header was not found
          AliWarning("AliAnalysisTaskEA MC header not found");
          fHistEvtSelection->Fill(9.5);
          return kTRUE; //skip the event
      }

      return kFALSE;  //there was not outlier all jets have pT below fMaxFacPtHard * ptHard

   } else {

      fHistEvtSelection->Fill(9.5);
      AliWarning("AliAnalysisTaskEA MC header not found");
      return kTRUE; //MC header not found
   }
}


//________________________________________________________________________
//Int_t  AliAnalysisTaskEA::GetMaxDistanceFromBorder(AliVCluster* cluster){
//   //Distance o
//   Int_t max = 0;
//
//   for(Int_t n=0; n<6; n++){
//      fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(n);
//      if(fFiducialCellCut->CheckCellFiducialRegion(fGeom, cluster, fCaloCells)){
//         max = n;
//      }else{
//         break;
//      }
//   }
//
//   return max;
//}
//
//
////________________________________________________________________________
//Bool_t AliAnalysisTaskEA::FinalClusterCuts(AliVCluster* cluster){
//
//   //General QA.
//   if( cluster->GetNCells() < 1) return kFALSE;
//
//   Int_t disToBad = cluster->GetDistanceToBadChannel();
//   if(-1<disToBad && disToBad<2) return kFALSE;
//
//   Int_t disToBorder = GetMaxDistanceFromBorder(cluster);
//   if(-1<disToBorder && disToBorder<1) return kFALSE;
//
//   Double_t lambda0 =  cluster->GetM02();  //avoid merged clusters
//   if(lambda0 > 0.4) return kFALSE;
//
//   return kTRUE;
//}
//
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

   fRho    = 0.;
   fRhoMC  = 0.;
   fRhoEMB = 0.;

   fRunnumber = 0;

   AliGenEventHeader* mcHeader = NULL;
   AliAODMCHeader* aodMCH = NULL;

   //+++++++++++++++++++++++++++++ check MC z vertex position ++++++++++++++++++++++++++
   if(fMode == AliAnalysisTaskEA::kKine){ //KINE
      AliRunLoader *rl = AliRunLoader::Instance();
      if(rl)  mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(rl->GetHeader()->GenEventHeader());
   }else if(fMode == AliAnalysisTaskEA::kMC){
      if(MCEvent()){

         mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
         if(!mcHeader){
            // Check if AOD
             aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

             if(aodMCH){
                for(UInt_t i = 0; i<aodMCH->GetNCocktailHeaders(); i++){
                  mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
                  if(mcHeader) break;
               }
            }
         }
      }
   }

   if(mcHeader){
      TArrayF pyVtx;
      mcHeader->PrimaryVertex(pyVtx);
      if(TMath::Abs(pyVtx[2]) > fZVertexCut) return kTRUE; //skip events with particle level vertex out of +-10 cm
   }

   //_________________________________________________________________
   // INITIALIZATION
   for(Int_t i=0; i<fkTTbins; i++){
      //fIndexTTC[i] = -1;
      fIndexTTH[i] = -1;
      //fIndexTTJ[i] = -1;

      fIndexTTH_PartLevel[i] = -1;
      //fIndexTTC_PartLevel[i] = -1;

      //fTTC[i].resize(0);
      fTTH[i].resize(0);
      //fTTJ[i].resize(0);

      fTTH_PartLevel[i].resize(0);
      //fTTC_PartLevel[i].resize(0);

      fdeltapT[i]  = 0.;
      fdeltapT_PartLevel[i]  = 0.;
   }


// for(Int_t i=0; i<fnClusterTTBins; i++){
//    fClusterTT[i] = 0;
//    fClusterTT_PartLevel[i] = 0;
// }

   for(Int_t i=0; i<fnHadronTTBins; i++){
      fHadronTT[i] = 0;
      fHadronTT_PartLevel[i] = 0;
   }

// for(Int_t i=0; i<fnJetChTTBins; i++){
//    fJetChTT[i] = 0;
// }

   //_________________________________________________________________

   fhJetPtEvtByEvent->Reset();
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      fhRecoilJetPtEvtByEvent[itt]->Reset();
   }
//   fhRecoilJetPtEvtByEventRandomTT->Reset();

   if(fMode == AliAnalysisTaskEA::kMC || fMode == kKine){
      fhJetPtEvtByEventPartLevel->Reset();
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         fhRecoilJetPtEvtByEventPartLevel[itt]->Reset();
      }
//      fhRecoilJetPtEvtByEventRandomTTPartLevel->Reset();
   }
   //_________________________________________________________________
   // EVENT SELECTION
   fHistEvtSelection->Fill(0.5); //Count input event

   //_________________________________________________________
   //READ  TRACK AND JET CONTAINERS
   //Container operations   http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerIterateTechniques

   if(fMode == AliAnalysisTaskEA::kNormal || fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kEmbedding || fMode == AliAnalysisTaskEA::kEmbPy){
      //fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(fMyTrackContainerName.Data())); //track container detector-level   real data only
      fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(0)); //track container detector-level   real data only
      fJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyJetContainerName.Data()));     //detector-level AKT jets real data or hybrid event
      fKTJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetContainerName.Data())); //detector-level KT jets real data or hybrid event

      fRho = GetMyRho(fKTJetContainerDetLevel); //estimated backround pt density
   }

   if( fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kEmbedding || fMode == kKine ||  fMode == AliAnalysisTaskEA::kEmbPy){  //particle level particles and jets  for  MC and embedding
      //fParticleContainerPartLevel = GetParticleContainer(fMyParticleContainerName.Data()); //pythia particle level particles
      fParticleContainerPartLevel = GetParticleContainer(1); //pythia particle level particles
      fJetContainerPartLevel      = static_cast<AliJetContainer*> (GetJetContainer(fMyJetParticleContainerName.Data()));   //pythia particle level AKT jets
      fKTJetContainerPartLevel    = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetParticleContainerName.Data()));   //pythia particle level KT jets

      fRhoMC = GetMyRho(fKTJetContainerPartLevel); //estimated backround pt density
   }

   if(fMode == AliAnalysisTaskEA::kEmbedding){ //Detector level pythia  for  embedding

      //fTrkContainerDetLevelEMB = static_cast<AliTrackContainer*> (GetTrackContainer(fMyDetLevelContainerName.Data())); //pythia detector level tracks
      fTrkContainerDetLevelEMB   = static_cast<AliTrackContainer*> (GetTrackContainer(2)); //pythia detector level tracks from AOD
      fJetContainerDetLevelEMB   = static_cast<AliJetContainer*> (GetJetContainer(fMyJetDetLevelContainerName.Data()));  //pythia detector level AKT jets
      fKTJetContainerDetLevelEMB = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetDetLevelContainerName.Data()));  //pythia detector level KT jets

      fRhoEMB = GetMyRho(fKTJetContainerDetLevelEMB); //estimated backround pt density
   }

   if(fMode == AliAnalysisTaskEA::kMC  || fMode == AliAnalysisTaskEA::kKine){
      if(IsOutlier()) return kFALSE;
   }
   //_________________________________________________________________
   // DECIDE WHETHER TO FILL SIGNAL TT OR REFERENCE TT  DEPENDING ON RANDOM  NUMBER
   fFillSigTT = kTRUE;
   if(fRandom->Integer(100) < 5) fFillSigTT = kFALSE;
   //________________________________________________________________
   //DATA ANALYSIS PARTICLE LEVEL
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


   fTrigflag[0] = fIsMinBiasTrig;
   fTrigflag[1] = fIsHighMultTrig;
   fTrigflag[2] = fIsEmcalTrig;


   InitEventProperties(); //FILIP
   if(fMode == AliAnalysisTaskEA::kMC  || fMode == AliAnalysisTaskEA::kKine)  FindParticleLevelTT(); // Modified by KA
   if(fMode == AliAnalysisTaskEA::kMC  || fMode == AliAnalysisTaskEA::kKine)  AnalyzeParticleLevel();
   if(fMode == AliAnalysisTaskEA::kKine) return kTRUE;

   //________________________________________________________________
   //DATA ANALYSIS DETECTOR LEVEL

   //Check Reconstructed event vertex and pileup
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec

    //_________________________________________________________________

   
   if(fMode == kEmbPy){                  EmbeddingFromTxtFile(); return kTRUE; }
   if(fMode == kEmbedding){              EmbeddingFromAODFile(); return kTRUE; }

   FindDetectorLevelTT(); //fk Find detector level TT

   if(fMode == AliAnalysisTaskEA::kMC)   FillResponseMatrix();
   if(fMode == AliAnalysisTaskEA::kMC)   FillResponseMatrix2D();

   GeneralTrackProperties();
   AnalyzeRawData();

   return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEA::FindDetectorLevelTT(){
   //fk Find Detector Level TT
   TLorentzVector myTT;
   AliVParticle *track = NULL;

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTH[i] = -1;
      fTTH[i].resize(0);
   }

   for(Int_t i=0; i < fnHadronTTBins; i++){
      fHadronTT[i] = 0;
      fHadronTT_Labels[i].resize(0); // Modified by KA
   }

   if(fIsMinBiasTrig || fIsHighMultTrig){
      if(fTrkContainerDetLevel){

         for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum()){
            // trackIterator is a std::map of AliTLorentzVector and AliVTrack
            track = trackIterator.second;  // Get the full track
            if(!track) continue;

            if(IsTrackInAcceptance(track, kDetLevel)){
               for(Int_t itt=0; itt < fnHadronTTBins; itt++){
                  if(fHadronTTLowPt[itt] < track->Pt() && track->Pt() < fHadronTTHighPt[itt]){
                     fHadronTT_Labels[itt].push_back(TMath::Abs(track->GetLabel())); // Modified by KA
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
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskEA::FindParticleLevelTT(){ // Modified by KA
   //inspired by "FindDetectorLevelTT" function

   // Find Particle Level TT
   TLorentzVector myTT;
   AliVParticle *mcParticle = NULL;

   for(Int_t i = 0; i < fkTTbins; i++){
      fIndexTTH_PartLevel[i] = -1;
      fTTH_PartLevel[i].resize(0);
   }

   for(Int_t i = 0; i < fnHadronTTBins; i++){
      fHadronTT_PartLevel[i]        = 0;
      fHadronTT_Labels_PartLevel[i].resize(0); // Modified by KA
   }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      if(fParticleContainerPartLevel){

         for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum()){
            mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!mcParticle)  continue;

            if(IsTrackInAcceptance(mcParticle, kPartLevel)){
               for(Int_t itt=0; itt<fnHadronTTBins; itt++){
                  if(fHadronTTLowPt[itt] < mcParticle->Pt() && mcParticle->Pt() < fHadronTTHighPt[itt]){
                     fHadronTT_Labels_PartLevel[itt].push_back(TMath::Abs(mcParticle->GetLabel())); // Modified by KA
                     myTT.SetPtEtaPhiM(mcParticle->Pt(), mcParticle->Eta(), mcParticle->Phi(), 0.);
                     fTTH_PartLevel[itt].push_back(myTT);
                     fHadronTT_PartLevel[itt]++; // there was a high pt
                  }
               }
            }
         }

         for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
            if(fHadronTT_PartLevel[itt] > 0){
               fIndexTTH_PartLevel[itt] = fRandom->Integer(fHadronTT_PartLevel[itt]);
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskEA::AnalyzeRawData(){
   //Analyze raw data

   AliEmcalJet  *jet = NULL;        //jet pointer real jet
   AliVParticle *track = NULL; //jet constituent
   Double_t ptLJ=-1, etaLJ=999, phiLJ=0; //leading jet
   Double_t ptSJ=-1, etaSJ=999, phiSJ=0; //subleading jet
   Int_t b1,b2;
   Double_t tmparr[4];
   Double_t tmparr3[3];
   Int_t idx;
   Double_t dphi     = 999.;
   Double_t DeltaPhi_Angle_Abs_0Pi = 0.0;
   Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho

   Int_t njets10_15 = 0;
   Int_t njetsAbove15 = 0;


   //     Find the leading and subleading jets  for estimates  of Delta pt and Exclude 2 leading jets
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

   //______________________________________________________________
   //                       TTH  ANALYSIS
   //______________________________________________________________

   if(fIsMinBiasTrig || fIsHighMultTrig){

      //chose trigger hadron TT
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         fdeltapT[itt] = 0.;

         if(fHadronTT[itt]>0){
            idx = fIndexTTH[itt];

            fdeltapT[itt] = GetDeltaPt(fTTH[itt][idx].Phi(), fTTH[itt][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, fRho, kDetLevel);
         }
      }

      for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
         if(!fTrigflag[itg]) continue; //check which trigger fired
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){

            fhMultTTH[itg][itt]->Fill(fHadronTT[itt]);

            if(!fHadronTT[itt]) continue; //check whether there was hadron TT

            fhRhoTTH[itg][itt]->Fill(fRho);
         }
      }

      for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
         jet = jetIterator.second;
         if(!jet)  continue;
         fhJetPtEvtByEvent->Fill(jet->Pt());  //Fill event by event inclusive jet spectrum
      }


      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         if(!fHadronTT[itt]) continue; //analyze events with hadron TT only

         for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
            if(!fTrigflag[itg]) continue; //check which trigger fired

            fhSignalTTH[itg][fkV0A][itt]->Fill(fMultV0A);
            fhSignalTTH[itg][fkV0C][itt]->Fill(fMultV0C);
            fhSignalTTH[itg][fkV0M][itt]->Fill(fMultV0M);
            fhSignalTTH[itg][fkV0Mnorm1][itt]->Fill(fMultV0Mnorm);

         }


         //pick up TTH hadron accoding to the index
         idx = fIndexTTH[itt];
         if(idx>-1){

            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!fTrigflag[itg]) continue; //check which trigger fired

               fhDeltaPtTTH_RC_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, fdeltapT[itt]);
            }

            if(fFillSigTT  && itt==0) continue;  // Do not fill reference
            if(!fFillSigTT && itt>0)  continue;  // Do not fill signal

            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!fTrigflag[itg]) continue; //check which trigger fired

               fhTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm,  fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            }

            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());
               DeltaPhi_Angle_Abs_0Pi = TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())); // Modified by KA.

               jetPtCorrDet = jet->Pt() - fRho*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                  if(!fTrigflag[itg]) continue; //check which trigger fired
                  fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, DeltaPhi_Angle_Abs_0Pi); // Modified by KA.

                  tmparr[0] = jet->Pt();
                  tmparr[1] = jet->Eta();
                  tmparr[2] = jet->Phi();
                  tmparr[3] = fMultV0Mnorm;
                  fhJetPtEtaPhiV0normTTH[itg][itt]->Fill(tmparr);
               }

               if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > TMath::Pi()/2){ //select recoil hemisphere and count jets
                  fhRecoilJetPtEvtByEvent[itt]->Fill(jet->Pt());
                  if(10.<= jet->Pt() && jet->Pt() < 15.)  njets10_15++;
		  if(15.<=jet->Pt()) njetsAbove15++;
               }
              

               if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > fPhiCut){     //select recoil jet
                  for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                     if(!fTrigflag[itg]) continue; //check which trigger fired
                     fhRecoilJetPtTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet);

                     tmparr3[0] = jetPtCorrDet;
                     tmparr3[1] = jet->Area();
                     tmparr3[2] = fMultV0Mnorm;
                     fhJetPtAreaV0normTTH[itg][itt]->Fill(tmparr3);
                  }
               }
            }

            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!fTrigflag[itg]) continue; //check which trigger fired

               //count number of jets with pT larger than something in recoil region
               tmparr3[2] = fMultV0Mnorm;
               for(Int_t ii = 1; ii<=fhNumberOfHighPtJetsRecoil[itg][itt]->GetAxis(0)->GetNbins(); ii++){
                  tmparr3[0] = fhNumberOfHighPtJetsRecoil[itg][itt]->GetAxis(0)->GetBinLowEdge(ii);
                  b1 = fhRecoilJetPtEvtByEvent[itt]->GetXaxis()->FindBin(tmparr3[0] + 1e-5);
                  b2 = fhRecoilJetPtEvtByEvent[itt]->GetXaxis()->GetNbins()+1;  //include overflow bin
                  tmparr3[1] = fhRecoilJetPtEvtByEvent[itt]->Integral(b1,b2);
                  fhNumberOfHighPtJetsRecoil[itg][itt]->Fill(tmparr3);
               }


               fhNjetsPt10_15[itg]->Fill(fMultV0Mnorm, njets10_15 + 1e-9);
	       fhNjetsPtAbove15[itg]->Fill(fMultV0Mnorm, njetsAbove15 + 1e-9);
            }
         }
      }

      //EMB_clus  DELTA PT FROM TRACK EMBEDDING
      if(fEmbeddPerpendicular){
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(!fHadronTT[itt]) continue; //analyze events with hadron TT only

	    fFastJetWrapper->Clear();
	    Double_t gen_pt = fRandom->Uniform(0.15,20);
	    TLorentzVector lVec;
            lVec.SetPtEtaPhiM(gen_pt, fTTH[itt][fIndexTTH[itt]].Eta(), fTTH[itt][fIndexTTH[itt]].Phi() + TMath::Pi()/2, 0);
            fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E(), -99999);

            for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               track = trackIterator.second;  // Get the full track
               if(!track) continue;

               if(IsTrackInAcceptance(track, kDetLevel)){
                  fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), 1);
	       }
	    }

	    fFastJetWrapper->Run();

	    std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();
            AliVParticle*  hytrk  = NULL;  //track hybrid event jet
            Double_t deltaPtEmb;
            Double_t sumTrkEmbeddedPt=0;

            for(UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet) {
               std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
               sumTrkEmbeddedPt=0;
               for(UInt_t ic = 0; ic < constituents.size(); ++ic){
                  if(constituents[ic].user_index() == -99999){
                     sumTrkEmbeddedPt += constituents[ic].pt();
                     break;
                  }
               }

               if(sumTrkEmbeddedPt>0){
                  deltaPtEmb = jets_incl.at(ijet).pt() - jets_incl.at(ijet).area() * fRho - sumTrkEmbeddedPt;


                  for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                     if(!fTrigflag[itg]) continue; //check which trigger fired

                     fhDeltaPtEmbeddPerpendicular[itg][itt]->Fill(fMultV0Mnorm, sumTrkEmbeddedPt, deltaPtEmb);
		  }
	       }
            }
         }
      }
   }


   //_________________________________________________________
   //      LOOP OVER JETS  DETECTOR LEVEL  TTJ
   //_________________________________________________________

   for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      jet = jetIterator.second;  // Get the pointer to jet object
      if(!jet)  continue;
//
      tmparr[0] = jet->Pt();
      tmparr[1] = jet->Eta();
      tmparr[2] = jet->Phi();
      tmparr[3] = fMultV0Mnorm;

      jetPtCorrDet = jet->Pt() - fRho*jet->Area();

      tmparr3[0] = jetPtCorrDet;
      tmparr3[1] = jet->Area();
      tmparr3[2] = fMultV0Mnorm;


      for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
         if(!fTrigflag[itg]) continue; //check which trigger fired

         fhJetPhiIncl[itg]->Fill(jetPtCorrDet, jet->Phi());
         fhJetEtaIncl[itg]->Fill(jetPtCorrDet, jet->Eta());
         fhJetPtEtaPhiV0norm[itg]->Fill(tmparr);

         fhJetPtAreaV0norm[itg]->Fill(tmparr3);
      }
   }


     return;
}
//_________________________________________________________________
void AliAnalysisTaskEA::InitEventProperties(){
   // EVENT PROPERTIES
   Double_t normV0A = -1;
   Double_t normV0C = -1;
   TString name;
   AliVParticle *track = NULL; //jet constituent
   Int_t    trackMult  = 0;

   if(fMode != AliAnalysisTaskEA::kKine){

      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(fMultSelection){
//
        fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");

        if(fMultFramework){
           normV0A = fMultSelection->GetZ("V0A");
           normV0C = fMultSelection->GetZ("V0C");
           fMultV0Mnorm =  fMultSelection->GetZ("V0M");
        }
      }else{
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

//      const AliVMultiplicity *mult = InputEvent()->GetMultiplicity();
//      if(mult){
//         fNTracklets = mult->GetNumberOfTracklets();
//      }else{
//         fNTracklets = -9999;
//      }

      fRunnumber  = InputEvent()->GetRunNumber();

      AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
      if(vzeroAOD){
         Double_t meanV0A, meanV0C;

         fMultV0A = vzeroAOD->GetMTotV0A();
         fMultV0C = vzeroAOD->GetMTotV0C();
         fMultV0M = fMultV0A + fMultV0C;

         if(fMode != AliAnalysisTaskEA::kMC){
            fMeanV0M = fHelperEA->GetV0M(fRunnumber);
         }else{
            fMeanV0M = fHelperEA->GetV0MDetLevel();
         }

         //V0 estimators of event activity normalized per minimum bias activity
         if(!fMultFramework){
            fMultV0Mnorm = fMultV0M/fMeanV0M;  //either from mult framework of from my analysis
          }

      }else{
         fMultV0A = -1;
         fMultV0C = -1;
         fMultV0M = -1;
         fMultV0Mnorm = -1;
      }
   }

   if((fIsEmcalTrig || fIsMinBiasTrig || fIsHighMultTrig) && (fMode != AliAnalysisTaskEA::kKine)){  //real data + mc det level + embedding

      //___________________________________________
      //    INCLUSIVE EVENTS (WITHOUT TT REQUIREMENT)


      for(Int_t itg=kMB; itg<=kGA; itg++){    //@@@
         if(!fTrigflag[itg]) continue;
         //events without TT requirement
         fhRho[itg]->Fill(fRho);

          fhCentrality[itg]->Fill(fCentralityV0M, fMultV0Mnorm);

          fhSignal[itg][fkV0A]->Fill(fMultV0A);
          fhSignal[itg][fkV0C]->Fill(fMultV0C);
          fhSignal[itg][fkV0M]->Fill(fMultV0M);
          fhSignal[itg][fkV0Mnorm1]->Fill(fMultV0Mnorm);

      }


      for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         track = trackIterator.second;  // Get the full track
         if(!track) continue;

         if(IsTrackInAcceptance(track, kDetLevel)){
            trackMult++;
         }
      }

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         if(!fTrigflag[itg]) continue;

         fhTrackMult[itg]->Fill(fMultV0Mnorm, trackMult);

      }

      if(fIsMinBiasTrig){ // run for all MC events   and for real data with min bias trigger
         fhVertex[0]->Fill(fxVertex);
         fhVertex[1]->Fill(fyVertex);
         fhVertex[2]->Fill(fzVertex);


         name = Form("%d", fRunnumber);
//         fhV0ARunByRunMB->Fill(name.Data(), fMultV0A, 1.0);
//         fhV0CRunByRunMB->Fill(name.Data(), fMultV0C, 1.0);
//         fhV0MRunByRunMB->Fill(name.Data(), fMultV0M, 1.0);
         fhV0MnormRunByRunMB->Fill(name.Data(), fMultV0Mnorm, 1.0);
      }
   }
}
//_________________________________________________________________
void AliAnalysisTaskEA::AnalyzeParticleLevel(){

   TLorentzVector myTT;
   AliEmcalJet  *jet = NULL;        //jet pointer real jet
   AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
   AliVParticle *track = NULL; //jet constituent
   AliVParticle *mcParticle = NULL; //mc particle
   AliAODMCParticle* daughtermc;
   Int_t b1,b2;
   Double_t ptLJmc=-1, etaLJmc=999, phiLJmc=0; //leading jet
   Double_t ptSJmc=-1, etaSJmc=999, phiSJmc=0; //subleading jet
   Int_t idx;
   Double_t tmparr3[3];
   Double_t jetPtCorrDet  = 0.;
   Double_t jetPtCorrPart = 0.;

   Double_t dphi             = 999.;
   Double_t DeltaPhi_Angle_Abs_0Pi = 0.0; // Modified by KA

   fMultV0A_PartLevel = 0.;
   fMultV0C_PartLevel = 0.;
   fMultV0M_PartLevel = 0.;

   fMultV0Mnorm_PartLevel = 0.;

   Bool_t isMBpartlevel = 0; // requires coincidence of V0A and V0C on parton level

   Int_t njets10_15 = 0;
   Int_t njetsAbove15 = 0;


   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
   
      fhRhoMBpart->Fill(fRhoMC);

      //Exclude 2 leading jets MC
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


      if(fParticleContainerPartLevel){
         //pT spectrum of particle level physical primary mc particles
         for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum()){
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
         fMultV0Mnorm_PartLevel = fMultV0M_PartLevel/fMeanV0M_PartLevel;

         fhSignal_V0M_trueMB_PartLevel->Fill(fMultV0M_PartLevel); //Fill V0M multiplicity without the conicidence condition

         if(isMBpartlevel){ //V0 COINCIDENCE
            fHistEvtSelection->Fill(8.5); //Count Accepted input event


            fhSignal_PartLevel[fkV0A]->Fill(fMultV0A_PartLevel);
            fhSignal_PartLevel[fkV0C]->Fill(fMultV0C_PartLevel);
            fhSignal_PartLevel[fkV0M]->Fill(fMultV0M_PartLevel); //V0M multiplicity
            fhSignal_PartLevel[fkV0Mnorm1]->Fill(fMultV0Mnorm_PartLevel);
         }
 
	 for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum()){
            mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!mcParticle)  continue;

            if(IsTrackInAcceptance(mcParticle, kPartLevel)){
               fhPtTrkTruePrimGen->Fill(mcParticle->Pt(), mcParticle->Eta(), fMultV0Mnorm_PartLevel); //Track Efficiency denominator
            }
	 }
     
         //++++++++++++++++++++++++++++++++++++++++++++++++++
         for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
            jet = jetIterator.second;
            if(!jet)  continue;
            fhJetPtEvtByEventPartLevel->Fill(jet->Pt());  //Fill event by event inclusive jet spectrum in CB
         }
         //++++++++++++++++++++++++++++++++++++++++++++++++++
      
         //chose trigger hadron TT   particle level
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(fHadronTT_PartLevel[itt]>0){
               idx = fIndexTTH_PartLevel[itt]; // Initialize at FindParticleLevelTT, KA


               fdeltapT_PartLevel[itt] = GetDeltaPt(fTTH_PartLevel[itt][idx].Phi(), fTTH_PartLevel[itt][idx].Eta(), phiLJmc, etaLJmc, phiSJmc, etaSJmc, fRhoMC, kPartLevel);
               //signal in events with hadron TT   particle level

               if(isMBpartlevel){
                  fhSignalTTH_PartLevel[fkV0A][itt]->Fill(fMultV0A_PartLevel);
                  fhSignalTTH_PartLevel[fkV0C][itt]->Fill(fMultV0C_PartLevel);
                  fhSignalTTH_PartLevel[fkV0M][itt]->Fill(fMultV0M_PartLevel);
                  fhSignalTTH_PartLevel[fkV0Mnorm1][itt]->Fill(fMultV0Mnorm_PartLevel);
               }
	   
               //hadron trigger particle level
               if(idx>-1){

                  if(isMBpartlevel){
                     fhRhoTTHinMBpart[itt]->Fill(fRhoMC);
                     fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fdeltapT_PartLevel[itt]);
                  }
                  if(fFillSigTT && itt==0) continue;  // Do not fill reference
                  if(!fFillSigTT && itt>0) continue;  // Do not fill signal

                  if(isMBpartlevel){ //FF
                     fhTTH_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
                  }
                  fhTTH_V0Mnorm_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, fTTH_PartLevel[itt][idx].Pt()); //FF
                  fhTTH_V0MnormDetLev_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm, fTTH_PartLevel[itt][idx].Pt()); //FILIP
	      
                  //recoil jets  PARTICLE LEVEL
                  for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
                     // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                     jet = jetIterator.second;  // Get the pointer to jet object
                     if(!jet)  continue;

                     dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH_PartLevel[itt][idx].Phi());
                     DeltaPhi_Angle_Abs_0Pi = TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH_PartLevel[itt][idx].Phi())); // Modified by KA.

	             jetPtCorrDet = jet->Pt() - fRhoMC*jet->Area();

                     if(isMBpartlevel){
                        fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet, DeltaPhi_Angle_Abs_0Pi); // Modified by KA
                        fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jet->Pt(), DeltaPhi_Angle_Abs_0Pi); // Modified by KA

                        if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > TMath::Pi()/2){ //select recoil hemisphere and count jets
                           fhRecoilJetPtEvtByEventPartLevel[itt]->Fill(jet->Pt());

			   if(10.<= jet->Pt() && jet->Pt()<15.) njets10_15++;
			   else if(15.<=jet->Pt())              njetsAbove15++;
                        }
                     }

                     fhRecoilJetPhiTTH_V0Mnorm_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet, DeltaPhi_Angle_Abs_0Pi); // Modified by KA
                     fhRecoilJetPhiTTH_V0MnormDetLev_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm, jetPtCorrDet, DeltaPhi_Angle_Abs_0Pi); // FILIP 

                     if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > fPhiCut){
                        //recoil jet hadron trigger
                        if(isMBpartlevel){
                           fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet);

	                           tmparr3[0] = jetPtCorrDet;
	                           tmparr3[1] = jet->Area();
	                           tmparr3[2] = fMultV0Mnorm_PartLevel;
	                           fhJetPtAreaV0normTTH_PartLevel[itt]->Fill(tmparr3);
                        }
                        fhRecoilJetPtTTH_V0Mnorm_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm_PartLevel, jetPtCorrDet);
                        fhRecoilJetPtZero_TTH_V0MnormDet_PartLevel[itt]->Fill(fMultV0Mnorm, jet->Pt()); //FILIP
                        fhRecoilJetPtTTH_V0MnormDetLev_AllMB_PartLevel[itt]->Fill(fMultV0Mnorm, jetPtCorrDet); //FILIP
                     }
                  }
	     
                  //count number of jets with pT larger than something in recoil region
                  for(Int_t ii = 1; ii<=fhNumberOfHighPtJetsRecoilPartLevel[itt]->GetAxis(0)->GetNbins(); ii++){
	             tmparr3[0] = fhNumberOfHighPtJetsRecoilPartLevel[itt]->GetAxis(0)->GetBinLowEdge(ii);
	             b1 = fhRecoilJetPtEvtByEventPartLevel[itt]->GetXaxis()->FindBin(tmparr3[0] + 1e-5);
	             b2 = fhRecoilJetPtEvtByEventPartLevel[itt]->GetXaxis()->GetNbins()+1;  //include overflow bin
	             tmparr3[1] = fhRecoilJetPtEvtByEventPartLevel[itt]->Integral(b1,b2);
                     tmparr3[2] = fMultV0Mnorm_PartLevel;
	             fhNumberOfHighPtJetsRecoilPartLevel[itt]->Fill(tmparr3);

                     tmparr3[2] = fMultV0Mnorm;//FILIP
	             fhNumberOfHighPtJetsRecoilPartLevelV0MnormDetLev[itt]->Fill(tmparr3);//FILIP
                  }

                  fhNjetsPt10_15_PartLevel->Fill(fMultV0Mnorm_PartLevel, njets10_15 + 1e-9);
                  fhNjetsPtAbove15_PartLevel->Fill(fMultV0Mnorm_PartLevel, njetsAbove15 + 1e-9);
                  fhNjetsPt10_15_PartLevelV0MnormDetLev->Fill(fMultV0Mnorm, njets10_15 + 1e-9); 
                  fhNjetsPtAbove15_PartLevelV0MnormDetLev->Fill(fMultV0Mnorm, njetsAbove15 + 1e-9); 
               }
            }
	 }
      }
   }

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fMode != AliAnalysisTaskEA::kKine) {
      //Response matrix normalization - spectrum of all generator level jets in acceptance
      if(fJetContainerPartLevel){
         for(auto jetPartIterator : fJetContainerPartLevel->accepted_momentum() ){
            jetPartMC = jetPartIterator.second;  // Get the pointer to mc particle object
            if(!jetPartMC)  continue;

            jetPtCorrPart = jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;

            fhJetPtPartLevelCorr->Fill(jetPtCorrPart);
            fhJetPtPartLevelZero->Fill(jetPartMC->Pt());


            tmparr3[0] = jetPtCorrPart;
            tmparr3[1] = jetPartMC->Area();
            tmparr3[2] = fMultV0Mnorm_PartLevel;

            fhJetPtAreaV0norm_PartLevel->Fill(tmparr3);

         }
      }
   }
}
//_________________________________________________________________
void AliAnalysisTaskEA::FillResponseMatrix(){

   if(fIsMinBiasTrig && fMode != AliAnalysisTaskEA::kKine){

      //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX
      AliEmcalJet  *jet = NULL;        //jet pointer real jet
      AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
      AliEmcalJet  *jetDetMC  = NULL;  //jet pointed detector level MC jet
      AliVParticle *track = NULL; //jet constituent
      AliVParticle *mcParticle = NULL; //mc particle
      Bool_t bRecPrim = kFALSE;
      Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
      Double_t jetPtCorrPart = 0.;
      Double_t array[5];
      Double_t xphi;

      //pT spectrum of detector level physical primary tracks and secondary tracks
      if(fTrkContainerDetLevel && fParticleContainerPartLevel){
         for(auto trkIterator : fTrkContainerDetLevel->accepted_momentum()){
            track = trkIterator.second;  // Get the pointer to mc particle object
            if(!track)  continue;

            if(!IsTrackInAcceptance(track, kDetLevel)) continue; //reconstructed level tracks
            bRecPrim = kFALSE; //not yet matched to generator level physical primary

            for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum()){
               mcParticle = mcPartIterator.second;  // Get the pointer to mc particle object
               if(!mcParticle)  continue;

               if(IsTrackInAcceptance(mcParticle, kPartLevel)){
                  if(TMath::Abs(track->GetLabel()) == TMath::Abs(mcParticle->GetLabel())){
                     //has the same label as reconstr track

                     bRecPrim = kTRUE;
                     fhPtTrkTruePrimRec->Fill(mcParticle->Pt(), mcParticle->Eta(), fMultV0Mnorm_PartLevel); //this is well recontr phys primary

		     //NEW angular smearing of inclusive tracks
                     fPhiSmearingTrack->Fill( fMultV0Mnorm, track->Pt(), TVector2::Phi_mpi_pi( mcParticle->Phi() - track->Phi())); //NEW  should jet pT be particle level or detector level?  

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


      //1) Find closest particle level and detector level jets
      //2) Get momentum shift due to fake tracks
      if(fJetContainerDetLevel){

         Double_t  sumAllTrackPtInJet   = 0.;
         Double_t  sumFakeTrackPtInJet  = 0.;

         for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
            jet = jetIterator.second;  // Get the pointer to jet object
            if(!jet)  continue;

            fhNormalization_Purity->Fill(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi())), 
                                         jet->Pt() - jet->Area()*fRho);

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
               jetPtCorrDet = jet->Pt() - jet->Area()*fRho;
               fhFractionOfSecInJet->Fill(jetPtCorrDet, sumFakeTrackPtInJet/sumAllTrackPtInJet);
            }

            //Fill Response matrix
            jetPartMC =  jet->ClosestJet();
            
            if(!jetPartMC){               
               fhPurity_InclusiveJets->Fill(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi())), 
                                            jet->Pt() - jet->Area()*fRho);
               continue;
            }

            // KA: Set new relax condition "5e-4" because we consider particle level jets with pT > 1 MeV
            if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

            jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
            jetPtCorrDet  =  jet->Pt() - jet->Area()*fRho;

            fhJetPtPartLevelVsJetPtDetLevelCorr->Fill(jetPtCorrDet,jetPtCorrPart); //response matrix
            fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr->Fill(jetPtCorrDet, jetPartMC->Pt()); //response matrix (added by KA)
            fhJetPtZeroPartLevelVsJetPtZeroDetLevel->Fill(jet->Pt(),jetPartMC->Pt()); //response matrix

            if(TMath::Abs(jetPartMC->Eta()) > 0.5) fhImpurityInclusive_DetJetPtVsPartJetPtCorr->Fill(jetPtCorrDet, jetPtCorrPart); // added by KA

            if(jetPtCorrPart > 0) fhJetPtResolutionVsPtPartLevel->Fill(jetPtCorrPart,(jetPtCorrDet-jetPtCorrPart)/jetPtCorrPart); //jet pT resolution
         }
      }

      //FF fill response matrix from recoil jets
      //Search for TT candidates in particle level physical primary mc particles
      //++++++++++++++++++++++++++++++++++++++++++++++++++
      Int_t idx;
      Int_t idx_PartLevel; // Modified by KA

      //chose trigger hadron TT particle level
      for(Int_t itt=0; itt < fnHadronTTBins; itt++){
         if(fFillSigTT  && itt == 0) continue;  //fk
         if(!fFillSigTT && itt > 0)  continue;  //fk

         if(fHadronTT[itt] > 0){
            idx = fIndexTTH[itt]; //fk initialized in FindDetectorLevelTT

            //recoil jets  PARTICLE LEVEL
            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;
               jetDetMC =  jet->ClosestJet();
               xphi = (jetDetMC) ?  jetDetMC->Phi() : jet->Phi();  //for those cases where there is detector level phi use detectro level phi, else use particle level phi

               if(TMath::Abs(TVector2::Phi_mpi_pi(xphi - fTTH[itt][idx].Phi())) > fPhiCut){  //fk
                  jetPtCorrPart = jet->Pt() - fRhoMC*jet->Area();
                  fhRecoilJetPtPartLevelCorr[itt]->Fill(jetPtCorrPart);
                  fhRecoilJetPtZeroPartLevel[itt]->Fill(jet->Pt());
               }
            }

            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum()){
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               jetPartMC =  jet->ClosestJet();
               if(!jetPartMC){
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){
                     array[0] = jet->Pt();
	                  array[1] = jet->Area();
	                  array[2] = jet->Eta();
	                  array[3] = TVector2::Phi_0_2pi(jet->Phi());
	                  array[4] = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());
                     fhNotMatchedJetPt[itt]->Fill(array);
                  }
                  continue;
               }
               if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

               jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
               jetPtCorrDet  =  jet->Pt() - jet->Area()*fRho;

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){   //fk
                  fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt]->Fill(jetPtCorrDet,jetPtCorrPart);
                  fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt]->Fill(jetPtCorrDet,jetPartMC->Pt());
                  fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt]->Fill(jet->Pt(),jetPartMC->Pt());
                  if(TMath::Abs(jetPartMC->Eta()) > 0.5) fhImpurityRecoil_DetJetPtVsPartJetPtCorr[itt]->Fill(jetPtCorrDet,jetPtCorrPart); // added by KA
               }
            }
         }
      }

      // KA: fill response matrix with jets on particle and detector levels recoiling from corresponding TT
      //Important: TT labels must be the same!!!
      Int_t Label_Det_TT, Label_Part_TT; // Modified by KA

      //FILIP: count number of jets with pT larger than 15 GeV and 20 GeV in particle and detector level events
      Double_t nJetsPt15_PartLevel = 1e-5; //FILIP
      Double_t nJetsPt15_DetLevel  = 1e-5; //FILIP
      Double_t nJetsPt20_PartLevel = 1e-5; //FILIP
      Double_t nJetsPt20_DetLevel  = 1e-5; //FILIP


      for(Int_t iTT = 0; iTT < fnHadronTTBins; iTT++){
         if(fFillSigTT  && iTT == 0) continue;  //fk
         if(!fFillSigTT && iTT > 0)  continue;  //fk

         if(fHadronTT[iTT] > 0 && fHadronTT_PartLevel[iTT] > 0){ // presence of TT on both levels

            Label_Part_TT = 0;
            idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT

            Label_Det_TT = TMath::Abs(fHadronTT_Labels[iTT][idx]); // Label of selected TT on detector level

            for(UInt_t k = 0; k < fHadronTT_Labels_PartLevel[iTT].size(); k++){
               if(Label_Det_TT == fHadronTT_Labels_PartLevel[iTT][k]){  // Look for TT on particle level with same label
                  idx_PartLevel = k;
                  Label_Part_TT = fHadronTT_Labels_PartLevel[iTT][k];
                  break;
               }
            }

            if(!Label_Part_TT) continue;

            fhTT_Corresp[iTT]->Fill(fTTH[iTT][idx].Pt(), fTTH_PartLevel[iTT][idx_PartLevel].Pt()); // Modified by KA
	 
            //recoil jets PARTICLE LEVEL
            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet) continue;

	       //FILIP: select recoil hemisphere and count jets with pT gt 15/20 GeV   on particle level
               if(fFillSigTT){ //just for TT{20,30}
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi() - fTTH_PartLevel[iTT][idx_PartLevel].Phi()) > TMath::Pi()/2)){ 
                     if(jet->Pt() > 15.0)  nJetsPt15_PartLevel++;      
                     if(jet->Pt() > 20.0)  nJetsPt20_PartLevel++;      
	          }
	       }//FILIP

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi() - fTTH_PartLevel[iTT][idx_PartLevel].Phi())) > fPhiCut){  // KA
                  jetPtCorrPart = jet->Pt() - fRhoMC*jet->Area();
                  fhRecoilJetPtPartLevel_CorrespTT[iTT]->Fill(jetPtCorrPart); // Modified by KA
                  fhRecoilJetPtZeroPartLevel_CorrespTT[iTT]->Fill(jet->Pt()); // Modified by KA
               }
            }

            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum()){
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet) continue;

               jetPartMC =  jet->ClosestJet();
               if(!jetPartMC) continue; // IMPORTANT TO ADD
               if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

               jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
               jetPtCorrDet  =  jet->Pt() - jet->Area()*fRho;

	       //FILIP: select recoil hemisphere and count jets with pT gt 15/20 GeV   on detector level 
	       if(fFillSigTT){ //just for TT{20,30}
	          if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[iTT][idx].Phi())) > TMath::Pi()/2){ //select recoil hemisphere and count jets
                     if(jet->Pt() > 15.0)  nJetsPt15_DetLevel++;      
                     if(jet->Pt() > 20.0)  nJetsPt20_DetLevel++;      
                  }//FILIP
	       }

               if(TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi()-fTTH[iTT][idx].Phi())) > fPhiCut && TMath::Abs(TVector2::Phi_mpi_pi(jetPartMC->Phi()-fTTH_PartLevel[iTT][idx_PartLevel].Phi())) > fPhiCut){   // KA: look for events when both matched jets are in the recoil
                  fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[iTT]->Fill(jetPtCorrDet, jetPtCorrPart); // Modified by KA
                  fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[iTT]->Fill(jetPtCorrDet, jetPartMC->Pt()); // Modified by KA
                  fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[iTT]->Fill(jet->Pt(), jetPartMC->Pt()); // Modified by KA
               }
            }

            //FILIP Fill Njet response matrices
            if(fFillSigTT){ //just for TT{20,30}
               fhNjetReMx_V0MnormDetLev_15GeV->Fill(nJetsPt15_DetLevel, nJetsPt15_PartLevel, fMultV0Mnorm); 
               fhNjetNorm_V0MnormDetLev_15GeV->Fill(nJetsPt15_PartLevel, fMultV0Mnorm);  
               
               fhNjetReMx_V0MnormDetLev_20GeV->Fill(nJetsPt20_DetLevel, nJetsPt20_PartLevel, fMultV0Mnorm);
               fhNjetNorm_V0MnormDetLev_20GeV->Fill(nJetsPt20_PartLevel, fMultV0Mnorm);
            }//FILIP

         }
      }


      //NEW -----------------------------------------------------------
      // make difference of delta phi  on detector and particle level
      if(fnHadronTTBins==2){  // there must be 2 trigger track bins, in such a case the index 1 is signal 
         if(fHadronTT[1] > 0){  //there is a TT{20,30} on the detector level
            if(fHadronTT_PartLevel[1] > 0){ // there is also TT{20,30} at particle level

               Int_t idxDetTT = -1;
               Int_t idxPartTT = -1;
               
	       //search for trigger track candidates which have the same MC labels on particle and detenctor level
               for(UInt_t id = 0; id < fTTH[1].size(); id++){
                  for(UInt_t ip = 0; ip < fTTH_PartLevel[1].size(); ip++){
                     if(fHadronTT_Labels_PartLevel[1][ip] == 0) continue;
	             if(fHadronTT_Labels_PartLevel[1][ip] != fHadronTT_Labels[1][id]) continue;
                     idxDetTT  = id;  //found match between detector level and particle level TT
                     idxPartTT = ip;
	             break;
	          }
	       }

	       if(idxDetTT> -1 && idxPartTT > -1){
	          Double_t deltaPhiDet  = -100; //   phi TT - phi jet detector level
	          Double_t deltaPhiPart = -100; //   phi TT - phi jet particle level
                  for(auto jetIterator : fJetContainerDetLevel->accepted_momentum()){
                     jet = jetIterator.second;  // Get the pointer to jet object
                     if(!jet) continue;
                 
                     jetPartMC =  jet->ClosestJet();
                     if(!jetPartMC) continue; // IMPORTANT TO ADD
                     if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost
                 
                     jetPtCorrDet  =  jet->Pt() - jet->Area()*fRho;
                 
                     deltaPhiDet  =  jet->Phi() - fTTH[1][idxDetTT].Phi();
	             deltaPhiPart =  jetPartMC->Phi() - fTTH_PartLevel[1][idxPartTT].Phi();
                     if(TMath::Abs(TVector2::Phi_mpi_pi(deltaPhiDet)) > TMath::Pi()/2){ //recoil jets only
                        fDeltaPhiSmearingTT2030Jet->Fill(fMultV0Mnorm, jetPtCorrDet, TVector2::Phi_mpi_pi(deltaPhiDet - deltaPhiPart));
		     }
	          }
	       }
	    }
         }  	
      }
      //NEW -----------------------------------------------------------
   }
   return;
}
//________________________________________________________________________
void AliAnalysisTaskEA::FillResponseMatrix2D(){
   //Fill 4D response matrix

   if(fIsMinBiasTrig && fMode == AliAnalysisTaskEA::kMC){
      //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX
      AliEmcalJet  *jet             = nullptr; //jet pointer real jet
      AliEmcalJet  *jetDetMC        = nullptr; //jet pointed detector level MC jet
      AliVParticle *track_PartLevel = nullptr; //jet constituent
      AliVParticle *track_DetLevel  = nullptr; //mc particle

      Int_t iterationStep = 0;
      Double_t smearing_Of_PhiAngle = 0;
      Double_t random_PhiAnglePartLevel = 0;
      Double_t random_PhiAngleDetLevel = 0;
      Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
      Double_t jetPtCorrPart = 0.;
      Double_t deltaPhi_angle_ParticleLevel = 0.;
      Double_t deltaPhi_angle_DetLevel = 0.;
      Double_t phi_Angle_Of_TT = 0.; //fk
      Int_t idx = 0; //fk
      Int_t idx_PartLevel; // Modified by KA

      //FILL JET RESPONSE MATRIX
      if(fJetContainerPartLevel){

         //Inclusive jets
         for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
            jet = jetIterator.second;  // Get the pointer to particle level jet object
            if(!jet)  continue;

            // Matching
            jetDetMC =  jet->ClosestJet();

            //Averaging over phi angle by random generation of phi angle PartLvl. Phi angle DetLvl is obtained by adding smearing to PartLvl phi angle
            random_PhiAnglePartLevel = TMath::Pi()*fRandom->Uniform(-1,1);  //(-pi,pi)
            jetPtCorrPart = jet->Pt() - jet->Area()*fRhoMC;

	         if(jet->Pt() > 5e-4){

               if(!jetDetMC){ //No associated MC detector level jet ->  Fill input for the miss function

                  fhPhi_JetPtPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jetPtCorrPart);
                  fhPhi_JetPtZeroPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jet->Pt()); // (added by KA)

                  //HM events
                  if(fMultV0Mnorm  < 4.0) continue; 
                  fhMissedEvents_HM->Fill(fMultV0Mnorm, TMath::Abs(random_PhiAnglePartLevel), jetPtCorrPart);

               } else if(jetDetMC->Pt() < 1e-10){ //associated MC detector level jet is a ghost jet ->  Fill input for the miss function

                  fhPhi_JetPtPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jetPtCorrPart);
                  fhPhi_JetPtZeroPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jet->Pt()); // (added by KA)

                  //HM events
                  if(fMultV0Mnorm  < 4.0) continue; 
                  fhMissedEvents_HM->Fill(fMultV0Mnorm, TMath::Abs(random_PhiAnglePartLevel), jetPtCorrPart);

               } else { //Associated Detector level jet is a physical jet -> fill response matrix
                  smearing_Of_PhiAngle = jetDetMC->Phi() - jet->Phi(); //Smearing of phi angle
                  random_PhiAngleDetLevel = TMath::Abs(TVector2::Phi_mpi_pi(random_PhiAnglePartLevel + smearing_Of_PhiAngle)); //TVector2::Phi_mpi_pi
                  jetPtCorrDet = jetDetMC->Pt() - jetDetMC->Area()*fRho;

                  fArray_for_filling[0] = random_PhiAngleDetLevel;
		            fArray_for_filling[1] = jetPtCorrDet;
		            fArray_for_filling[2] = TMath::Abs(random_PhiAnglePartLevel);
		            fArray_for_filling[3] = jetPtCorrPart;

	               fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->Fill(fArray_for_filling);

                  //Particle level jet pT is not corrected by RhokT
                  fArray_for_filling[3] = jet->Pt();
                  fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->Fill(fArray_for_filling); // added by KA

                  //Jet pT is not corrected by RhokT
                  fArray_for_filling[1] = jetDetMC->Pt();
                  fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->Fill(fArray_for_filling); // added by KA

		            //ANGLUAR SMEARING OF INCLUSIVE JETS
                  fPhiSmearingJet->Fill(fMultV0Mnorm, jetPtCorrDet, TVector2::Phi_mpi_pi(smearing_Of_PhiAngle)); //NEW  should jet pT be particle level or detector level?

                  //------------------------------------------------
                  //Jet energy and angle resolution
                  vector<double> tmp_vector;
                  if(jetPtCorrPart >= 0.0 && jetPtCorrPart <= 100.0){ // Specified range for jet pT. To exclude overflow/underflow bins
                     tmp_vector.push_back(fMultV0Mnorm); // Multiplicity norm, det level
                     tmp_vector.push_back(TVector2::Phi_mpi_pi(smearing_Of_PhiAngle)); // phi anlge resolution
                     tmp_vector.push_back(TVector2::Phi_0_2pi(jet->Phi())); // phi angle, part level jet
                     tmp_vector.push_back((jetPtCorrDet - jetPtCorrPart)/jetPtCorrPart); // jet pT resolution
                     tmp_vector.push_back(jetPtCorrPart); // part level jet pT
                     fhJetEnergy_Angle_Resol->Fill(&tmp_vector[0]);
                  }
                  tmp_vector.clear();

                  //------------------------------------------------
                  //Inclusive response matrix in HM events
                  //To save memory, exclude filling of overflow/underflow bins
                  if(fMultV0Mnorm  < 4.0) continue; 
                     
                  tmp_vector.push_back(fMultV0Mnorm);                         // Multiplicity norm, det level
                  tmp_vector.push_back(random_PhiAngleDetLevel);              // phi anlge det
                  tmp_vector.push_back(jetPtCorrDet);                         // jet pT det
                  tmp_vector.push_back(TMath::Abs(random_PhiAnglePartLevel)); // phi anlge part
                  tmp_vector.push_back(jetPtCorrPart);                        // jet pT part
                  fhInclusiveRespMatrix_HM->Fill(&tmp_vector[0]);
                  tmp_vector.clear();
               }
            }
         }//end inclusive jets

         //TT events
         for(Int_t iTT = 0; iTT < fnHadronTTBins; iTT++){
            if(fHadronTT[iTT] == 0) continue;      //fk skip event, since no TT on detector level
            if(fFillSigTT  && iTT == 0) continue;  //fk
            if(!fFillSigTT && iTT > 0)  continue;  //fk

            idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT
            phi_Angle_Of_TT = fTTH[iTT][idx].Phi(); //fk
            // fhNumberOf_ChosenTT_PartLevel[iTT]->Fill(fTTH_PartLevel[iTT][index_OF_RandomTT].Pt());

            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               jetDetMC = jet->ClosestJet();  // Matched detector level jet

               //For TT events we take Delta phi angle = TT phi angle - Jet phi angle
               deltaPhi_angle_ParticleLevel = TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi() - phi_Angle_Of_TT)); // changed DeltaPhi angle (previous range 0,2*Pi)
               jetPtCorrPart = jet->Pt() - jet->Area()*fRhoMC;

               if(jet->Pt() > 5e-4){

                  if(!jetDetMC){  //no matched detector level jet

                     fhDeltaPhi_JetPtPartLevel[iTT]->Fill(deltaPhi_angle_ParticleLevel, jetPtCorrPart);
                     fhDeltaPhi_JetPtZero_PartLevel[iTT]->Fill(deltaPhi_angle_ParticleLevel, jet->Pt()); // added by KA

                  }else if(jetDetMC->Pt() < 1e-10) {  //matched to a ghost

                     fhDeltaPhi_JetPtPartLevel[iTT]->Fill(deltaPhi_angle_ParticleLevel, jetPtCorrPart);
                     fhDeltaPhi_JetPtZero_PartLevel[iTT]->Fill(deltaPhi_angle_ParticleLevel, jet->Pt()); // added by KA

                  }else{
                     deltaPhi_angle_DetLevel = TMath::Abs(TVector2::Phi_mpi_pi(jetDetMC->Phi() - phi_Angle_Of_TT)); //fk
                     jetPtCorrDet = jetDetMC->Pt() - jetDetMC->Area()*fRho;

                     fArray_for_filling[0] = deltaPhi_angle_DetLevel;
		               fArray_for_filling[1] = jetPtCorrDet;
		               fArray_for_filling[2] = deltaPhi_angle_ParticleLevel;
		               fArray_for_filling[3] = jetPtCorrPart;
                     fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[iTT]->Fill(fArray_for_filling);

                     //Particle level jet pT is NOT corrected for RhokT
		               fArray_for_filling[3] = jet->Pt();
                     fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[iTT]->Fill(fArray_for_filling); // added by KA

                     //Particle level jet pT is NOT corrected for RhokT
		               fArray_for_filling[1] = jetDetMC->Pt();
                     fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[iTT]->Fill(fArray_for_filling); // added by KA
                  }
               }
            }//end recoil jets
         }

         // KA: fill response matrix with jets on particle and detector levels recoiling from corresponding TT
         Int_t Label_Det_TT, Label_Part_TT; // Modified by KA

         for(Int_t iTT = 0; iTT < fnHadronTTBins; iTT++){
            if(fFillSigTT  && iTT == 0) continue;  //fk
            if(!fFillSigTT && iTT > 0)  continue;  //fk

            if(fHadronTT[iTT] > 0 && fHadronTT_PartLevel[iTT] > 0){ // presence of TT on both levels

               Label_Part_TT = 0;
               idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT
               // idx_PartLevel = fIndexTTH_PartLevel[iTT]; //KA initialized in FindParticleLevelTT

               Label_Det_TT = TMath::Abs(fHadronTT_Labels[iTT][idx]); // Label of selected TT on detector level

               for(UInt_t k = 0; k < fHadronTT_Labels_PartLevel[iTT].size(); k++){
                  if(Label_Det_TT == fHadronTT_Labels_PartLevel[iTT][k]){  // Look for TT on particle level with same label
                     idx_PartLevel = k;
                     Label_Part_TT = fHadronTT_Labels_PartLevel[iTT][k];
                     break;
                  }
               }

               if(!Label_Part_TT) continue;

               for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
                  jet = jetIterator.second;  // Get the pointer to jet object
                  if(!jet)  continue;

                  jetDetMC = jet->ClosestJet();  // Matched detector level jet

                  //For TT events we take Delta phi angle = TT phi angle - Jet phi angle
                  deltaPhi_angle_ParticleLevel = TMath::Abs(TVector2::Phi_mpi_pi(jet->Phi() - fTTH_PartLevel[iTT][idx_PartLevel].Phi())); // changed DeltaPhi angle (previous range 0,2*Pi)
                  jetPtCorrPart = jet->Pt() - jet->Area()*fRhoMC;

                  if(jet->Pt() > 5e-4){

                     if(!jetDetMC){  //no matched detector level jet

                        fhDeltaPhi_JetPtPartLevel_CorrespTT[iTT]->Fill(deltaPhi_angle_ParticleLevel, jetPtCorrPart); // Modified by KA
                        fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[iTT]->Fill(deltaPhi_angle_ParticleLevel, jet->Pt()); // Modified by KA

                     }else if(jetDetMC->Pt() < 1e-10) {  //matched to a ghost

                        fhDeltaPhi_JetPtPartLevel_CorrespTT[iTT]->Fill(deltaPhi_angle_ParticleLevel, jetPtCorrPart); // Modified by KA
                        fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[iTT]->Fill(deltaPhi_angle_ParticleLevel, jet->Pt()); // Modified by KA

                     }else{

                        deltaPhi_angle_DetLevel = TMath::Abs(TVector2::Phi_mpi_pi(jetDetMC->Phi() - fTTH[iTT][idx].Phi())); //fk
                        jetPtCorrDet = jetDetMC->Pt() - jetDetMC->Area()*fRho;

                        fArray_for_filling[0] = deltaPhi_angle_DetLevel;
		                  fArray_for_filling[1] = jetPtCorrDet;
		                  fArray_for_filling[2] = deltaPhi_angle_ParticleLevel;
		                  fArray_for_filling[3] = jetPtCorrPart;
                        fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[iTT]->Fill(fArray_for_filling); // Modified by KA

                        //Particle level jet pT is NOT corrected for RhokT
                        fArray_for_filling[3] = jet->Pt();
                        fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[iTT]->Fill(fArray_for_filling); // added by KA

                        //Particle level jet pT is NOT corrected for RhokT
                        fArray_for_filling[1] = jetDetMC->Pt();
                        fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[iTT]->Fill(fArray_for_filling); // added by KA
                     }
                  }
               }//end recoil jets
            }
         }
      }
   }
}
//________________________________________________________________________
void  AliAnalysisTaskEA::EmbeddingFromTxtFile(){
   //   EMBEDDED EVENTS FROM PYTHIA TEXT FILE

   AliEmcalJet  *jet = NULL;        //jet pointer real jet
   AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
   AliEmcalJet  *jetDetMC  = NULL;  //jet pointed detector level MC jet
   AliVParticle *track = NULL; //jet constituent
   AliVParticle *mcParticle = NULL; //mc particle
   Bool_t bRecPrim = kFALSE;
   Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
   Double_t jetPtCorrPart = 0.;
   TLorentzVector myTT;
   Int_t idx = -1;
   Double_t dphi = 999.;

   for(Int_t i=0; i<fkTTbins; i++){
      fTTH[i].resize(0);
   }
   for(Int_t i=0; i<fnHadronTTBins; i++){
      fHadronTT[i] = 0;
   }

   if(fMode == kEmbPy){

      if(!fIsMinBiasTrig && !fIsHighMultTrig) return; //if this is not MB or HM  skip the rest

      if(fParticleContainerPartLevel){

         //detector level pythia mc particles
         for(auto mcDetIterator : fParticleContainerPartLevel->accepted_momentum() ){
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
            fIndexTTH[itt] = -1;
            if(fHadronTT[itt]>0){
               fIndexTTH[itt] = fRandom->Integer(fHadronTT[itt]);
            }
         }

         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            idx = fIndexTTH[itt];//hadron trigger
            if(idx<0) continue;


            for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
               if(!fTrigflag[itg]) continue;
               fhTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            }

            //recoil jets pythia detector level event
            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());

               jetPtCorrDet = jet->Pt() - fRhoMC*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!fTrigflag[itg]) continue;
                  fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
               }
            }

            //recoil jets in the combined event
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());

               jetPtCorrDet = jet->Pt() - fRho*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!fTrigflag[itg]) continue;
                  fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
               }

               //fill similar disribution but just for jets which have pythia partner
//               jetDetMC =  jet->ClosestJet(); //This is the closes pythia Detector level jet
//               if(jetDetMC){
//
//                  dphi = TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[itt][idx].Phi());
//
//                  jetPtCorrDet  =  jetDetMC->Pt() - jetDetMC->Area()*fRhoMC;
//
//                  for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
//                     if(!fTrigflag[itg]) continue;
//                     fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
//                  }
//               }
            }//jet loop
         }//TT loop
      }//EMB track container
   }
   return;
}


//________________________________________________________________________
void  AliAnalysisTaskEA::EmbeddingFromAODFile(){
   //   EMBEDDED EVENTS FROM AOD
   AliEmcalJet  *jet = NULL;        //jet pointer real jet
   AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
   AliEmcalJet  *jetDetMC  = NULL;  //jet pointed detector level MC jet
   AliVParticle *track = NULL; //jet constituent
   AliVParticle *mcParticle = NULL; //mc particle
   TLorentzVector myTT;
   Int_t idx = -1;
   Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
   Double_t jetPtCorrPart = 0.;
   Double_t sharedFraction = 0.; //shared fraction between detector level and combined level
   Double_t dphi = 999.0;

   for(Int_t i=0; i<fkTTbins; i++){
      fTTH[i].resize(0);
   }
   for(Int_t i=0; i<fnHadronTTBins; i++){
      fHadronTT[i] = 0;
   }

   if(fMode == kEmbedding){

      if(!fIsMinBiasTrig && !fIsHighMultTrig) return; //if this is not MB or HM  skip the rest

      const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
      if(!embeddingHelper) return;
      //double ptHardBin = embeddingHelper->GetPtHardBin();

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         if(!fTrigflag[itg]) continue;
         fhTrialsEMBtot[itg]->Fill(0.5, embeddingHelper->GetPythiaTrials());
         fhXsectionEMBtot[itg]->Fill(0.5, embeddingHelper->GetPythiaXSection());

         //fhTrialsEMB[itg]->Fill( ptHardBin, embeddingHelper->GetPythiaTrials());
         //fhXsectionEMB[itg]->Fill( ptHardBin, embeddingHelper->GetPythiaXSection());
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
            fIndexTTH[itt] = -1;
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
               if(!fTrigflag[itg]) continue;
               fhTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            }

            //recoil jets pythia detector level event
            for(auto jetIterator : fJetContainerDetLevelEMB->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());

               jetPtCorrDet = jet->Pt() - fRhoEMB*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!fTrigflag[itg]) continue;
                  fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
               }
            }

            //recoil jets in the combined event
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());

               jetPtCorrDet = jet->Pt() - fRho*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
                  if(!fTrigflag[itg]) continue;
                  fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);

               }

               //fill similar disribution but just for jets which have pythia partner
//               jetDetMC =  jet->ClosestJet(); //This is the closes pythia Detector level jet
//               if(jetDetMC){
//
//                  dphi = TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[itt][idx].Phi());
//
//                  jetPtCorrDet  =  jetDetMC->Pt() - jetDetMC->Area()*fRhoEMB;
//
//                  for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
//                     if(!fTrigflag[itg]) continue;
//                     fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, dphi);
//                  }
//               }
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

            jetPtCorrPart = jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;

            for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
               if(!fTrigflag[itg]) continue; //inclusive jet spectrum for ReMx normalization
               fhJetPtPartLevelCorr_EMB[itg]->Fill(jetPtCorrPart);
               fhJetPtPartLevelZero_EMB[itg]->Fill(jetPartMC->Pt());
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
               if(!fTrigflag[itg]) continue; //recoil jet spectrum for ReMx normalization
               fhSharedJetFraction[itg]->Fill(jetDetMC->Pt(), sharedFraction);
            }

            if(sharedFraction < fMinFractionShared) continue;

            jetPartMC = jetDetMC->ClosestJet(); //This is the closes pythia particle level jet to the pythia detector level jet
            if(!jetPartMC) continue;
            if(jetPartMC->Pt()<1e-3) continue; //prevents matching with a ghost

            jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
            jetPtCorrDet  =  jet->Pt() - jet->Area()*fRho;

            for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
               if(!fTrigflag[itg]) continue; //recoil jet spectrum for ReMx normalization
               fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg]->Fill(jetPtCorrDet, jetPtCorrPart); //response matrix
               fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg]->Fill(jetPtCorrDet, jetPartMC->Pt()); //response matrix
            }
         }
      }
   }

   return;
}
//________________________________________________________________________
void AliAnalysisTaskEA::GeneralTrackProperties(){

   Double_t xyz[50];
   Double_t pxpypz[50];
   Double_t cv[21];
   Double_t tmparr[4];

   Int_t label, labelMC;
   Bool_t labelfound = 0;
   AliAODMCParticle* particleMC = NULL;
   AliVParticle *track = NULL; //jet constituent
   AliAODTrack *trackAOD=NULL ;




   //_________________________________________________________
   //LOOP OVER TRACKS DETECTOR LEVEL

   for(auto trackIterator : fTrkContainerDetLevel->accepted_momentum() ){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      track = trackIterator.second;  // Get the full track
      if(!track) continue;

      if(IsTrackInAcceptance(track, kDetLevel)){

         for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
            if(!fTrigflag[itg]) continue;
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

            trackAOD = dynamic_cast <AliAODTrack*>( track);
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
   //delete fFiducialCellCut;
   delete fHelperEA;
   if(fEmbeddPerpendicular){
      delete fFastJetWrapper; //EMB_clus
   }
}
//________________________________________________________________________
void AliAnalysisTaskEA::UserCreateOutputObjects(){
  // called once to create user defined output objects like histograms, plots etc.
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.
  //fOutput TList defined in the mother class

   if(kOldV0MC){
      fHelperEA->SetV0MeanForMCWithDeltaElectronBug(); //cout set old V0M MC values which suffered from delta electron bug
   }

   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name, object;
   TString rhotype = "kt";
   //if(fRhoType == krhokt)  rhotype ="kt";
   //else rhotype = "cms";

   fRandom = new TRandom3(0);
   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1D("fHistEvtSelection", "event selection", 10, 0, 10);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"events IN"); //0-1
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"incomplete DAQ (rejected)"); //1-2
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)"); //2-3
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)"); //3-4
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"MB"); //4-5
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"EMCAL"); //5-6
   fHistEvtSelection->GetXaxis()->SetBinLabel(7,"High Mult"); //6-7
   fHistEvtSelection->GetXaxis()->SetBinLabel(8,"IsEventSelected"); //7-8
   fHistEvtSelection->GetXaxis()->SetBinLabel(9,"MC MB events"); //8-9
   fHistEvtSelection->GetXaxis()->SetBinLabel(10,"IsOutlier"); //9-10


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
   if(fMode != AliAnalysisTaskEA::kKine)  fOutput->Add(fhVertexZall);

   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);
   if(fMode != AliAnalysisTaskEA::kKine)  fOutput->Add(fhVertexZ);

   //-------------------------
   TString trig[]={"MB","HM","GA"};


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      name   = Form("fhTrackEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive track vs pT %s",trig[itg].Data());
      fhTrackEtaIncl[itg] = new TH2D(name.Data(), object.Data(), 200, 0, 200, 60, -0.9, 0.9); // Modified by KA (old value: 50,0,100)
      fOutput->Add((TH2D*) fhTrackEtaIncl[itg]);

      name   = Form("fhTrackPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist tracks vs pT %s",trig[itg].Data());
      fhTrackPhiIncl[itg] = new TH2D(name.Data(), object.Data(), 200, 0, 200, 50, 0, 2*TMath::Pi()); // Modified by KA (old value: 50,0,100)
      fOutput->Add((TH2D*) fhTrackPhiIncl[itg]);
   }

   if(fMode == AliAnalysisTaskEA::kEmbedding || fMode == AliAnalysisTaskEA::kEmbPy){
     fhTrackEtaInclEMB = (TH2D*) fhTrackEtaIncl[kMB]->Clone("fhTrackEtaInclEMB");
     fOutput->Add((TH2D*) fhTrackEtaInclEMB);
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      name   = Form("fhJetEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive jets vs pTjet %s",trig[itg].Data());
      fhJetEtaIncl[itg] = new TH2D(name.Data(),object.Data(), 150, -20, 130, 60, -0.9, 0.9); // Modified by KA (old value: 40,-0.9,0.9)
      fOutput->Add((TH2D*) fhJetEtaIncl[itg]);

      name   = Form("fhJetPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist jets vs pTjet %s",trig[itg].Data());
      fhJetPhiIncl[itg] = new TH2D(name.Data(),object.Data(), 60, -20, 100, 50, 0, 2*TMath::Pi());
      fOutput->Add((TH2D*) fhJetPhiIncl[itg]);
   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      name   = Form("fhClusterEtaIncl%s",trig[itg].Data());
//      object = Form("Eta dist inclusive clusters vs pT %s",trig[itg].Data());
//      fhClusterEtaIncl[itg] = new TH2D( name.Data(), object.Data(), 100, 0, 100, 40,-0.9,0.9);
//      fOutput->Add((TH2D*) fhClusterEtaIncl[itg]);
//
//      name   = Form("fhClusterPhiIncl%s",trig[itg].Data());
//      object = Form("Azim dist clusters vs pT %s",trig[itg].Data());
//      fhClusterPhiIncl[itg] = new TH2D( name.Data(), object.Data(), 50, 0, 100, 50,0,2*TMath::Pi());
//      fOutput->Add((TH2D*) fhClusterPhiIncl[itg]);
//   }

   const Int_t ktdim = 4;
   Int_t   tbins[ktdim] = { 50,   40, 80, 10};
   Double_t txmin[ktdim] = { 0., -0.9, 0,  0.};
   Double_t txmax[ktdim] = {50.,  0.9, 2*TMath::Pi(), 10.};


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      name = Form("fhTrackPtEtaPhiV0norm_%s",trig[itg].Data());

      fhTrackPtEtaPhiV0norm[itg] = new  THnSparseF(name.Data(),"Tracks pt eta phi V0nom", ktdim, tbins, txmin, txmax);
      fOutput->Add((THnSparse*) fhTrackPtEtaPhiV0norm[itg]);
   }

   //jets
   const Int_t kjetdim = 4;
   Int_t    jetbins[kjetdim] = {100,  60,   80, 10}; // Modified by KA (old value: x,40,x,x)
   Double_t jetxmin[kjetdim] = { 0., -0.9,  0,  0.};
   Double_t jetxmax[kjetdim] = {100.,  0.9, 2*TMath::Pi(), 10.};


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      name = Form("fhJetPtEtaPhiV0norm_%s",trig[itg].Data());

      fhJetPtEtaPhiV0norm[itg] = new  THnSparseF(name.Data(),"Jet pt eta phi V0nom", kjetdim, jetbins, jetxmin, jetxmax);
      fOutput->Add((THnSparse*) fhJetPtEtaPhiV0norm[itg]);
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
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
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      name = Form("fhJetPtAreaV0norm_%s_Rho%s",trig[itg].Data(), rhotype.Data());

      fhJetPtAreaV0norm[itg] = new  THnSparseF(name.Data(),"Jet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
      fOutput->Add((THnSparse*) fhJetPtAreaV0norm[itg]);

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         if(fMode == AliAnalysisTaskEA::kKine) continue;

         name = Form("fhJetPtAreaV0norm_%s_TTH%d_%d_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhJetPtAreaV0normTTH[itg][itt] = new  THnSparseF(name.Data(),"Jet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
         fOutput->Add((THnSparse*) fhJetPtAreaV0normTTH[itg][itt]);
      }
   }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      name = Form("fhJetPtAreaV0norm_Rho%s_PartLevel", rhotype.Data());

      fhJetPtAreaV0norm_PartLevel = new  THnSparseF(name.Data(),"Part LevelJet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
      fOutput->Add((THnSparse*) fhJetPtAreaV0norm_PartLevel);

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         name = Form("fhJetPtAreaV0norm_MB_TTH%d_%d_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhJetPtAreaV0normTTH_PartLevel[itt] = new  THnSparseF(name.Data(),"Part LevelJet V0Mnorm pt eta phi", kjdim, jbins, jxmin, jxmax);
         fOutput->Add((THnSparse*) fhJetPtAreaV0normTTH_PartLevel[itt]);
      }
   }

   //RHO
   for(Int_t itg=kMB; itg<=kGA; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      name   = Form("hRho%s_%s",rhotype.Data(),trig[itg].Data());
      object = Form("Rho %s det level %s",rhotype.Data(),trig[itg].Data());

      fhRho[itg] = new TH1D( name.Data(), object.Data(),1000,0,100);
      fOutput->Add((TH1D*) fhRho[itg]);
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         name = Form("hRho%s_%s_TTH%d_%d", rhotype.Data(), trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRhoTTH[itg][itt] = (TH1D*)  fhRho[itg]->Clone(name.Data());      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTH[itg][itt]);
      }

//      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){  //JET TT
//         name = Form("hRho%s_%s_TTJ%d_%d", rhotype.Data(), trig[itg].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
//         fhRhoTTJ[itg][ijj] = (TH1D*)  fhRho[itg]->Clone(name.Data());                      //! in events MB with hadron TT
//         fOutput->Add((TH1D*) fhRhoTTJ[itg][ijj]);
//      }
   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){ //GAMMA TT
//         name = Form("hRho%s_%s_TTC%d_%d", rhotype.Data(), trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//         fhRhoTTC[itg][igg] = (TH1D*)  fhRho[itg]->Clone(name.Data());                      //! in events MB with hadron TT
//         fOutput->Add((TH1D*) fhRhoTTC[itg][igg]);
//      }
//   }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      name = Form("hRho%s_MB_part", rhotype.Data());
      fhRhoMBpart = new TH1D( name.Data(), name.Data(),1000,0,100);
      fhRhoMBpart->SetTitle(Form("Rho %s  min bias part level", rhotype.Data()));
      fOutput->Add((TH1D*) fhRhoMBpart);

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){
         name = Form("hRho%s_MB_TTH%d_%d_part", rhotype.Data(), fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRhoTTHinMBpart[itt] = (TH1D*) fhRhoMBpart->Clone(name.Data());                      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTHinMBpart[itt]);
      }

//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("hRho%s_MB_TTC%d_%d_part", rhotype.Data(), fClusterTTLowPt[igg],fClusterTTHighPt[igg]);
//         fhRhoTTCinMBpart[igg] = (TH1D*)  fhRhoMBpart->Clone(name.Data());                      //! in events MB with hadron TT
//         fOutput->Add((TH1D*) fhRhoTTCinMBpart[igg]);
//      }
   }

   //VERTEX
   fhVertex[0] = new TH1D("hVertexX","VertexX",100,-1,1);
   fhVertex[1] = new TH1D("hVertexY","VertexY",100,-1,1);
   fhVertex[2] = new TH1D("hVertexZ","VertexZ",400,-20,20);

   for(Int_t iv=0; iv<fkVtx;iv++){
      if(fMode != AliAnalysisTaskEA::kKine) fOutput->Add((TH1D*) fhVertex[iv]);
   }


   Int_t nRun = fHelperEA->GetNRuns();

   if(fMode != AliAnalysisTaskEA::kKine){
//      fhV0MRunByRunMB = new TH2D("fhV0MRunByRunMB","fhV0MRunByRunMB", nRun, 0, nRun, 180,0,1800);
//      for(Int_t ir=0; ir < nRun; ir++){
//         fhV0MRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d", fHelperEA->GetRun(ir)));
//      }
//      fOutput->Add((TH2D*) fhV0MRunByRunMB);
//
//      name = "fhV0ARunByRunMB";
//      fhV0ARunByRunMB = (TH2D*)  fhV0MRunByRunMB->Clone(name.Data());
//      fhV0ARunByRunMB->SetTitle(name.Data());
//      fOutput->Add((TH2D*) fhV0ARunByRunMB);
//
//      name = "fhV0CRunByRunMB";
//      fhV0CRunByRunMB = (TH2D*)  fhV0MRunByRunMB->Clone(name.Data());
//      fhV0CRunByRunMB->SetTitle(name.Data());
//      fOutput->Add((TH2D*) fhV0CRunByRunMB);
//
      fhV0MnormRunByRunMB = new TH2D("fhV0MnormRunByRunMB","fhV0MnormRunByRunMB", nRun, 0, nRun, 200,0,20);
      for(Int_t ir=0; ir < nRun; ir++){
         fhV0MnormRunByRunMB->GetXaxis()->SetBinLabel(ir+1,Form("%d", fHelperEA->GetRun(ir)));
      }
      fOutput->Add((TH2D*) fhV0MnormRunByRunMB);
   }

   //CENTRALITY
//   TString cest[] = {"V0A", "V0C", "V0M", "V0Mnorm"}; //centrality estimators

//   const Int_t narrV0 = 1700;
//   Double_t arrV0[narrV0+1];
//   for(Int_t i=0; i<1600; i++){
//      arrV0[i]=0.5*i;  //0-800
//   }
//   for(Int_t i=0; i<=100; i++){
//      arrV0[1600+i] = 800 + 10.*i;  //800-1800
//   }
//


   Double_t arrcent[] = {
     0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
     0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
     0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
     0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
     0.5,0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
     1.,
     2.,
     3.,
     4.,
     5.,
     6.,
     7.,
     8.,
     9.,
     10,11,12,13,14,15,16,17,18,19,20,
     25,30,35,40,45,50,55,60,65,70,75,80,85,90,100};

   Int_t narrcent = sizeof(arrcent)/sizeof(Double_t)-1;


   //CENTRALITY
   for(Int_t itg=kMB; itg<=kGA;itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      name = Form("hCentrality_%s_V0M", trig[itg].Data());
      fhCentrality[itg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 400,0,20);
      fOutput->Add((TH2D*) fhCentrality[itg]);
   }


//   //CENTRALITY TTH
//   for(Int_t itg=kMB; itg<=kHM;itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t ic=0; ic<fkCE; ic++){
//         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//            name = Form("hCentrality_%s_%s_TTH%d_%d", trig[itg].Data(), cest[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//            fhCentralityTTH[itg][ic][itt] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
//            fOutput->Add((TH2D*) fhCentralityTTH[itg][ic][itt]);
//         }
//      }
//   }
//
//   //TTJ MB
//   for(Int_t itg=kMB; itg<=kHM;itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t ic=0; ic<fkCE;ic++){
//         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
//            name = Form("hCentrality_%s_%s_TTJ%d_%d", trig[itg].Data(), cest[ic].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
//            fhCentralityTTJ[itg][ic][ijj] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
//            fOutput->Add((TH2D*) fhCentralityTTJ[itg][ic][ijj]);
//         }
//      }
//   }
//
//   //TTC  MB
//   for(Int_t itg=kMB; itg<=kGA;itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t ic=0; ic<fkCE;ic++){
//         for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
//            name = Form("hCentrality_%s_%s_TTC%d_%d", trig[itg].Data(), cest[ic].Data(), fClusterTTLowPt[ijj], fClusterTTHighPt[ijj]);
//            fhCentralityTTC[itg][ic][ijj] = (TH2D*) fhCentrality[itg][ic]->Clone(name.Data());
//            fOutput->Add((TH2D*) fhCentralityTTC[itg][ic][ijj]);
//         }
//      }
//   }

   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //SIGNAL

   TString signal[]={"multV0A", "multV0C", "multV0M","multV0Mnorm"};
   Float_t signalL[]={0,0,0,0};
   Float_t signalH[]={1000,1000,1800,15};
   Int_t   signalN[]={100,100,180,150};

   for(Int_t itg=kMB; itg<=kGA;itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;

      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_%s_%s", trig[itg].Data(), signal[ic].Data());
         fhSignal[itg][ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignal[itg][ic]);
      }
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t ic=0; ic<fkCE; ic++){ //TT hadron
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_%s_%s_TTH%d_%d", trig[itg].Data(), signal[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTH[itg][ic][itt] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
            fOutput->Add((TH1D*) fhSignalTTH[itg][ic][itt]);
         }
      }

//      for(Int_t ic=0; ic<fkCE; ic++){ //TT jet
//         for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
//            name = Form("hSignal_%s_%s_TTJ%d_%d", trig[itg].Data(), signal[ic].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
//            fhSignalTTJ[itg][ic][ijj] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
//            fOutput->Add((TH1D*) fhSignalTTJ[itg][ic][ijj]);
//         }
//      }
   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t ic=0; ic<fkCE; ic++){ //HM && TT jet
//         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//            name = Form("hSignal_%s_%s_TTC%d_%d",trig[itg].Data(), signal[ic].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//            fhSignalTTC[itg][ic][igg] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
//            fOutput->Add((TH1D*) fhSignalTTC[itg][ic][igg]);
//         }
//      }
//   }



   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){ //PARTICLE LEVEL SIGNAL DISTRIBUTIONS

      TString signalmc[]={"multV0A", "multV0C", "multV0M", "multV0Mnorm"};
      Float_t signalLmc[]={0,0,0,0};
      Float_t signalHmc[]={500,500,500,20};
      Int_t signalNmc[]={500,500,500,200};

      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_MB_%s_PartLevel", signalmc[ic].Data());
         fhSignal_PartLevel[ic] = new TH1D(name.Data(), name.Data(), signalNmc[ic], signalLmc[ic], signalHmc[ic]);
         fOutput->Add((TH1D*) fhSignal_PartLevel[ic]);
      }

      name = Form("fhSignal_V0M_trueMB_PartLevel");
      fhSignal_V0M_trueMB_PartLevel = new TH1D(name.Data(), name.Data(), 500, 0, 500);
      fOutput->Add((TH1D*) fhSignal_V0M_trueMB_PartLevel);


      //TT hadron
      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_MB_%s_TTH%d_%d_PartLevel", signalmc[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTH_PartLevel[ic][itt] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
            fOutput->Add((TH1D*) fhSignalTTH_PartLevel[ic][itt]);
         }
      }

      //TT cluster
//      for(Int_t ic=0; ic<fkCE;ic++){ //MB
//         for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//            name = Form("hSignal_MB_%s_TTC%d_%d_PartLevel", signalmc[ic].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//            fhSignalTTC_PartLevel[ic][igg] = new TH1D(name.Data(),name.Data(),signalNmc[ic], signalLmc[ic], signalHmc[ic]);
//            fOutput->Add((TH1D*) fhSignalTTC_PartLevel[ic][igg]);
//         }
//      }
   }


   // V0 assymetery versus V0norm
//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      name = Form("fhV0MAssymVsV0Mnorm_%s",trig[itg].Data());
//      fhV0MAssymVsV0Mnorm[itg] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
//      fOutput->Add((TH2D*) fhV0MAssymVsV0Mnorm[itg]);
//   }
//
//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      name = Form("fhV0MAssymVsV0Mnorm_MB_PartLevel");
//      fhV0MAssymVsV0Mnorm_PartLevel = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
//      fOutput->Add((TH2D*) fhV0MAssymVsV0Mnorm_PartLevel);
//   }
//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//         name = Form("fhV0MAssymVsV0Mnorm_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhV0MAssymVsV0MnormTTH[itg][itt] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
//         fOutput->Add((TH2D*) fhV0MAssymVsV0MnormTTH[itg][itt]);
//      }
//   }
//
//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//         name = Form("fhV0MAssymVsV0Mnorm_MB_TTH%d_%d_PartLevel",  fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhV0MAssymVsV0MnormTTH_PartLevel[itt] = new TH2D(name.Data(),name.Data(),10,0,10,21,-1,1.1);
//         fOutput->Add((TH2D*) fhV0MAssymVsV0MnormTTH_PartLevel[itt]);
//      }
//   }
//
//
//   for(Int_t itg=kMB; itg<=kGA; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      name = Form("fhV0A_V0C_V0Mnorm_%s", trig[itg].Data());
//      fhV0A_V0C_V0Mnorm[itg] = new TH3D(name.Data(),name.Data(),100,0,1000, 100,0,1000, 10, 0, 10);
//      fOutput->Add((TH3D*) fhV0A_V0C_V0Mnorm[itg]);
//   }
//
//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      name = Form("fhV0A_V0C_PartLevel_MB");
//      fhV0A_V0C_PartLevel = new TH2D(name.Data(),name.Data(),100,0,200, 100,0,200);
//      fOutput->Add((TH2D*) fhV0A_V0C_PartLevel);
//
//      name = Form("fhV0A_V0APartLevel_MB");
//      fhV0A_V0APartLevel = new TH2D(name.Data(),name.Data(),100,0,1000, 100,0,200);
//      fOutput->Add((TH2D*) fhV0A_V0APartLevel);
//
//      name = Form("fhV0C_V0CPartLevel_MB");
//      fhV0C_V0CPartLevel = (TH2D*)  fhV0A_V0APartLevel->Clone(name.Data());
//      fOutput->Add((TH2D*) fhV0C_V0CPartLevel);
//   }
//
//
   //name = Form("fhV0MvsV0Mnorm_MB");
   //fhV0MvsV0Mnorm = new TH2D(name.Data(),name.Data(),100,0,40, 100,0,1200);
   //fOutput->Add((TH2D*) fhV0MvsV0Mnorm);


//   name = Form("fhV0AvsSPD_MB");
//   fhV0AvsSPD = new TH2D(name.Data(),name.Data(),100,0,500, 100,0,500);
//   if(fMode != AliAnalysisTaskEA::kKine)  fOutput->Add((TH2D*) fhV0AvsSPD);
//
//   name = Form("fhV0CvsSPD_MB");
//   fhV0CvsSPD = new TH2D(name.Data(),name.Data(),100,0,500, 100,0,500);
//   if(fMode != AliAnalysisTaskEA::kKine)  fOutput->Add((TH2D*) fhV0CvsSPD);
//

//   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      name = Form("fhV0AvsV0C_MB_TTH%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
//      fhV0AvsV0CTTH[itt] = new TH2D(name.Data(), name.Data(),100,0,1000,100,0,1000);
//      fhV0AvsV0CTTH[itt]->SetTitle(name.Data());
//      fOutput->Add((TH2D*) fhV0AvsV0CTTH[itt]);
//   }
//   for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      name = Form("fhV0AvsV0C_MB_TTJ%d_%d", fJetChTTLowPt[ijj],fJetChTTHighPt[ijj]);
//      fhV0AvsV0CTTJ[ijj] =  (TH2D*)  fhV0AvsV0CTTH[0]->Clone(name.Data());
//      fhV0AvsV0CTTJ[ijj]->SetTitle(name.Data());
//      fOutput->Add((TH2D*) fhV0AvsV0CTTJ[ijj]);
//   }
//   for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      name = Form("fhV0AvsV0C_MB_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
//      fhV0AvsV0CTTCinMB[ijj] = (TH2D*)  fhV0AvsV0CTTH[0]->Clone(name.Data());
//      fhV0AvsV0CTTCinMB[ijj]->SetTitle(name.Data());
//      fOutput->Add((TH2D*) fhV0AvsV0CTTCinMB[ijj]);
//   }
//
//
//   if(fMode != AliAnalysisTaskEA::kMC && fMode != AliAnalysisTaskEA::kKine){
//      for(Int_t ijj=0; ijj<fnClusterTTBins; ijj++){
//         name = Form("fhV0AvsV0C_GA_TTC%d_%d", fClusterTTLowPt[ijj],fClusterTTHighPt[ijj]);
//         fhV0AvsV0CTTCinGA[ijj] = (TH2D*)  fhV0AvsV0CTTH[0]->Clone(name.Data());
//         fhV0AvsV0CTTCinGA[ijj]->SetTitle(name.Data());
//         fOutput->Add((TH2D*) fhV0AvsV0CTTCinGA[ijj]);
//      }
//   }
   //+++++++++++++++++++++++++++++++
   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      name = Form("fhTrackMult%s", trig[itg].Data());
      fhTrackMult[itg] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 1000, 0, 1000);
      fOutput->Add((TH2D*) fhTrackMult[itg]);
   }


   //Trigger track candidate multiplicity
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hMultTT_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhMultTTH[itg][itt] = new TH1D(name.Data(),name.Data(),100,0,100);
         fOutput->Add((TH1D*)  fhMultTTH[itg][itt]);
      }
   }

//   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTJ
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t ijj=0; ijj<fnJetChTTBins; ijj++){
//         name = Form("hMultTT_%s_TTJ%d_%d", trig[itg].Data(), fJetChTTLowPt[ijj], fJetChTTHighPt[ijj]);
//         fhMultTTJ[itg][ijj] = new TH1D(name.Data(),name.Data(),100,0,100);
//         fOutput->Add((TH1D*) fhMultTTJ[itg][ijj]);
//      }
//   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTC
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("hMultTT_%s_TTC%d_%d", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//         fhMultTTC[itg][igg] = new TH1D(name.Data(),name.Data(),100,0,100);
//         fOutput->Add((TH1D*) fhMultTTC[itg][igg]);
//      }
//   }


   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if(fMode == AliAnalysisTaskEA::kKine) continue;
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

//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//         name = Form("hTT_%s_3D_TTH%d_%d_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhTTH_3D_V0Mnorm1[itg][itt] = new TH3D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 21, -1, 1.1, 100, 0, 100);
//         fOutput->Add((TH3D*) fhTTH_3D_V0Mnorm1[itg][itt]);
//      }

   }


   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm1_PartLevel[itt]);
      }

      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm without V0 coincidence
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0Mnorm_AllMB_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm_AllMB_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm_AllMB_PartLevel[itt]);

	 name = Form("hTT_MB_TTH%d_%d_V0MnormDetLev_AllMB_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]); //FILIP
         fhTTH_V0MnormDetLev_AllMB_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0MnormDetLev_AllMB_PartLevel[itt]);

         //Modified by KA
         name = Form("NumberOf_Corresp_TT%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTT_Corresp[itt] = new TH2D(name.Data(),name.Data(), 100, 0, 100, 100, 0, 100);
         fOutput->Add((TH2D*) fhTT_Corresp[itt]);
      }

//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
//         name = Form("hTT_MB_3D_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
//         fhTTH_3D_V0Mnorm1_PartLevel[itt] = new TH3D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 21, -1, 1.1, 100, 0, 100);
//         fOutput->Add((TH3D*) fhTTH_3D_V0Mnorm1_PartLevel[itt]);
//      }
   }


   //TT emcal cluster pT spectrum single inclusive  in MB   with V0M
//   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//      //   name = Form("hTT_%s_TTC%d_%d_CentV0M", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//      //   fhTTC_CentV0M[itg][igg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 1000, 0, 100);
//      //   fOutput->Add((TH2D*) fhTTC_CentV0M[itg][igg]);
//      //}
//
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("hTT_%s_TTC%d_%d_V0Mnorm", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//         fhTTC_V0Mnorm1[itg][igg] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
//         fOutput->Add((TH2D*) fhTTC_V0Mnorm1[itg][igg]);
//      }
//   }

//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      //TT emcal cluster pT spectrum single inclusive  in MB   with V0Mnorm
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("hTT_MB_TTC%d_%d_V0Mnorm_PartLevel", fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//         fhTTC_V0Mnorm1_PartLevel[igg] = new TH2D(name.Data(), name.Data(),  nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
//         fOutput->Add((TH2D*) fhTTC_V0Mnorm1_PartLevel[igg]);
//      }
//   }

   //reference with shifted rho
//   Int_t shiftMeV;
//   TString sign;
//   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t is=0; is<fkShift; is++){  //TTH
//         shiftMeV = TMath::Nint(1000*( -0.3 + is*0.01));
//         if(shiftMeV<0){
//            sign = "Minus";
//         }else if(shiftMeV>0){
//            sign = "Plus";
//         }
//         name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rho%s_rhoShift%s%dMeV",
//            trig[itg].Data(), fHadronTTLowPt[0], fHadronTTHighPt[0], rhotype.Data(), sign.Data(), TMath::Abs(shiftMeV));
//
//         fhRecoilJetPtTTHref_V0Mnorm1_rhoShift[itg][is] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);
//         fOutput->Add((TH2D*) fhRecoilJetPtTTHref_V0Mnorm1_rhoShift[itg][is]);
//      }
//   }

   // ++++++++++++++++++++++++++++++++++ Binning ++++++++++++++++++++++++++++++++++++
   //Jet pT corrected on RhokT (added by KA)
   const Int_t numberBins_jetpT_RhokT = 270;
   Double_t jetPtRhokT_Bins[numberBins_jetpT_RhokT + 1];
   for(Int_t i = 0; i <= numberBins_jetpT_RhokT; ++i) jetPtRhokT_Bins[i] = i - 20;  //(-20,250)

   //Jet pT not corrected on RhokT (added by KA)
   const Int_t numberBins_jetpT_Zero = 250;
   Double_t jetPtZero_Bins[numberBins_jetpT_Zero + 1];
   for(Int_t i = 0; i <= numberBins_jetpT_Zero; ++i) jetPtZero_Bins[i] = i;  //(0,250)

   //Phi angle
   const Int_t ndeltaPhiBins = 40; // Modified by KA. Previous value was 80
   Double_t deltaPhiBins[ndeltaPhiBins+1];
   Double_t p = TMath::Pi()/ndeltaPhiBins;
   for(Int_t i = 0; i <= ndeltaPhiBins; ++i) deltaPhiBins[i] = i*p;  // Modified by KA. Previous range was (0,2*pi)

   //Phi angle: Inclusive jets
   // Double_t deltaPhiBins_InclusiveJets[2*ndeltaPhiBins+1]; // Modified by KA
   // Double_t step = 0.5*TMath::TwoPi()/ndeltaPhiBins;
   // for(Int_t i = 0; i <= 2*ndeltaPhiBins; i++) deltaPhiBins_InclusiveJets[i] = i*step;

   //V0M normalized
   Double_t arrV0Mnorm[nbinsV0Mnorm+1];
   p = maxV0Mnorm/nbinsV0Mnorm;
   for(Int_t i = 0; i<=nbinsV0Mnorm; i++){
      arrV0Mnorm[i] = i*p; //(0,20.0)
   }

   //Bins for sumTrkEmbeddedPt ("fhDeltaPtEmbeddPerpendicular" histogram)
   const Int_t number_sumTrkEmbeddedPtBins = 20;
   Double_t sumTrkEmbeddedPtBins[number_sumTrkEmbeddedPtBins+1];
   for(Int_t i=0; i<=number_sumTrkEmbeddedPtBins; i++) sumTrkEmbeddedPtBins[i] = i;  //(0,20.0)

   //dphi of recoil jets associated to semi-inclusive hadron TT in MB with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("RecoilJetPhi_%s_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPhiTTH_V0Mnorm1[itg][itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); // added by KA
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]);
      }
   }

   //recoil jet distribution as a function V0norm, V0 assymetery, jet pt, jet |dphi|
//   const Int_t rldim = 4;
//   Int_t   rlbins[ktdim] = {10,  21, 130, 50};
//   Double_t rlmin[ktdim] = { 0., -1, -10,  0.};
//   Double_t rlmax[ktdim] = {10., 1.1, 120, TMath::Pi()};
//
//
//   for(Int_t itg=kMB; itg<=kHM; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
//         name = Form("fhRecoilJet4D_%s_TTH%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
//
//         fhRecoilJetTTH_V0Mnorm1[itg][itt] = new  THnSparseF(name.Data(),"V0norm, V0 assym, jet pT,  abs(dphi)", rldim, rlbins, rlmin, rlmax);
//         fOutput->Add((THnSparse*) fhRecoilJetTTH_V0Mnorm1[itg][itt]);
//      }
//   }

   //RECOIL JET SPECTRA
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      //for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
      //   name = Form("fhRecoilJetPt_%s_TTH%d_%d_CentV0M", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      //   fhRecoilJetPtTTH_CentV0M[itg][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);
      //   fOutput->Add((TH2D*) fhRecoilJetPtTTH_CentV0M[itg][itt]);
      //}

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPtTTH_V0Mnorm1[itg][itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); //changed range: old->(200,   -20, 180) (added by KA)
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm1[itg][itt]);
      }
   }

   //TTH recoil jet distributions for MC
   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){ //! recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         // Jet pT is CORRECTED on RhokT
         name = Form("RecoilJetPt_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // changed range: old-> (200, -20, 180) (added by KA)
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm1_PartLevel[itt]);

         // Jet pT is NOT CORRECTED on RhokT
         name = Form("RecoilJetPtZero_MB_TTH%d_%d_V0MnormDet_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPtZero_TTH_V0MnormDet_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_Zero, jetPtZero_Bins); // range (250, 0, 250) (added by KA)
         fOutput->Add((TH2D*) fhRecoilJetPtZero_TTH_V0MnormDet_PartLevel[itt]);


         name = Form("RecoilJetPt_AllMB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPtTTH_V0Mnorm_AllMB_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // changed range: old-> (200, -20, 180) (added by FF)
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm_AllMB_PartLevel[itt]);

         name = Form("RecoilJetPt_AllMB_TTH%d_%d_V0MnormDetLev_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data()); //FILIP
         fhRecoilJetPtTTH_V0MnormDetLev_AllMB_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // FILIP 
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0MnormDetLev_AllMB_PartLevel[itt]);

      }

      // dphi of recoil jets associated to semi-inclusive hadron TT in MB with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         // Jet pT is CORRECTED on RhokT
         name = Form("RecoilJetPhi_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); // added by KA
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[itt]);

         // Jet pT is NOT CORRECTED on RhokT
         name = Form("RecoilJetPtZero_DeltaPhi_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm1_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_Zero, jetPtZero_Bins, ndeltaPhiBins, deltaPhiBins); // added by KA
         fOutput->Add((TH3D*) fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm1_PartLevel[itt]);

         // Jet pT is CORRECTED on RhokT
         name = Form("RecoilJetPhi_AllMB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhRecoilJetPhiTTH_V0Mnorm_AllMB_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); // added by FF
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm_AllMB_PartLevel[itt]);

         name = Form("RecoilJetPhi_AllMB_TTH%d_%d_V0MnormDetLev_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data()); //FILIP
         fhRecoilJetPhiTTH_V0MnormDetLev_AllMB_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); // FILIP 
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0MnormDetLev_AllMB_PartLevel[itt]);
      }


//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
//         name = Form("fhRecoilJet4D_MB_TTH%d_%d_V0Mnorm_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
//         fhRecoilJetTTH_V0Mnorm1_PartLevel[itt] = new  THnSparseF(name.Data(),"V0norm, V0 assym, jet pT,  abs(dphi)", rldim, rlbins, rlmin, rlmax);
//         fOutput->Add((THnSparse*) fhRecoilJetTTH_V0Mnorm1_PartLevel[itt]);
//      }
   }



   if(fMode == AliAnalysisTaskEA::kEmbedding || fMode == AliAnalysisTaskEA::kEmbPy){
      //! dphi of recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);

      for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){

            //!  filled with any detector level pythia recoil jet
            name = Form("fhRecoilJetPhi_%s_EMB_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
            fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt] = (TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]->Clone(name.Data());
            fOutput->Add((TH3D*) fhRecoilJetPhiTTH_EMB_V0Mnorm1[itg][itt]);

            //!  filled  tagged closest detector level pythia recoil jet
//            name = Form("fhRecoilJetPhi_%s_TAG_TTH%d_%d_Rho%s_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
//            fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt] = (TH3D*) fhRecoilJetPhiTTH_V0Mnorm1[itg][itt]->Clone(name.Data());
//            fOutput->Add((TH3D*) fhRecoilJetPhiTTH_TAG_V0Mnorm1[itg][itt]);
         }
      }
   }

   //+++++++++++++++++++++++++++ RECOIL JETS WITH TTC ++++++++++++++++++++++++++++++++++++++++++
//   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M centrality
//      //   name = Form("fhRecoilJetPt_%s_TTC%d_%d_CentV0M", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//      //   fhRecoilJetPtTTC_CentV0M[itg][igg] = (TH2D*) fhRecoilJetPtTTH_CentV0M[0][0]->Clone(name.Data());
//      //   fOutput->Add((TH2D*) fhRecoilJetPtTTC_CentV0M[itg][igg]);
//      //}
//
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
//         name = Form("fhRecoilJetPt_%s_TTC%d_%d_V0Mnorm_Rho%s", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg], rhotype.Data());
//         fhRecoilJetPtTTC_V0Mnorm1[itg][igg] = (TH2D*) fhRecoilJetPtTTH_V0Mnorm1[0][0]->Clone(name.Data());
//         fOutput->Add((TH2D*) fhRecoilJetPtTTC_V0Mnorm1[itg][igg]);
//      }
//   }

//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){  //! recoil jets associated to semi-inclusive cluster TT  in MB  with V0M
//         name = Form("fhRecoilJetPt_MB_TTC%d_%d_V0Mnorm_Rho%s_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg], rhotype.Data());
//         fhRecoilJetPtTTC_V0Mnorm1_PartLevel[igg] = (TH2D*) fhRecoilJetPtTTH_V0Mnorm1_PartLevel[0]->Clone(name.Data());
//         fOutput->Add((TH2D*) fhRecoilJetPtTTC_V0Mnorm1_PartLevel[igg]);
//      }
//   }


   //delta pT distributions versus V0M CENTRALITY
//   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      //for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0M
//      //   name = Form("fhDeltaPtTTH_%s_RC_CentV0M_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//      //   fhDeltaPtTTH_RC_CentV0M[itg][itt] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 200, -20, 180);
//      //   fOutput->Add((TH2D*) fhDeltaPtTTH_RC_CentV0M[itg][itt]);
//      //}
//   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTH
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      //for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//      //   name = Form("fhDeltaPtTTC_%s_RC_CentV0M_TTC%d_%d", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg]);
//      //   fhDeltaPtTTC_RC_CentV0M[itg][igg] = (TH2D*) fhDeltaPtTTH_RC_CentV0M[0][0]->Clone(name.Data());
//      //   fOutput->Add((TH2D*) fhDeltaPtTTC_RC_CentV0M[itg][igg]);
//      //}
//   }

   //delta pT distributions versus V0Mnorm   = V0M/mean V0M
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhDeltaPtTTH_%s_RC_V0Mnorm_TTH%d_%d_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhDeltaPtTTH_RC_V0Mnorm1[itg][itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // changed range: old-> (200, -20, 180) (added by KA)
         fOutput->Add((TH2D*) fhDeltaPtTTH_RC_V0Mnorm1[itg][itt]);
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhDeltaPtTTH_%s_EMB_V0Mnorm_TTH%d_%d_Rho%s", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhDeltaPtEmbeddPerpendicular[itg][itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, number_sumTrkEmbeddedPtBins, sumTrkEmbeddedPtBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); //Old: (200, -20, 180) added by KA
         fOutput->Add((TH3D*) fhDeltaPtEmbeddPerpendicular[itg][itt]);
      }
   }

//   for(Int_t itg=kMB; itg<=kGA; itg++){  //TTC
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kGA) continue;
//
//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("fhDeltaPtTTC_%s_RC_V0Mnorm_TTC%d_%d_Rho%s", trig[itg].Data(), fClusterTTLowPt[igg], fClusterTTHighPt[igg], rhotype.Data());
//         fhDeltaPtTTC_RC_V0Mnorm1[itg][igg] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);
//         fOutput->Add((TH2D*) fhDeltaPtTTC_RC_V0Mnorm1[itg][igg]);
//      }
//   }


   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in HM  with V0M
         name = Form("fhDeltaPtTTH_MB_RC_V0Mnorm_TTH%d_%d_Rho%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], rhotype.Data());
         fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // old (200, -20.0. 180) added by KA
         fOutput->Add((TH2D*) fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[itt]);
      }

//      for(Int_t igg=0; igg<fnClusterTTBins; igg++){
//         name = Form("fhDeltaPtTTC_RC_V0Mnorm_TTC%d_%d_Rho%s_PartLevel", fClusterTTLowPt[igg],fClusterTTHighPt[igg], rhotype.Data());
//         fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[igg] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 200, -20, 180);
//         fOutput->Add((TH2D*) fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[igg]);
//      }
   }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      fhPtTrkTruePrimGen = new TH3D("fhPtTrkTruePrimGen","fhPtTrkTruePrimGen",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkTruePrimGen);
   }

   if(fMode == AliAnalysisTaskEA::kMC){
      fhPtTrkTruePrimRec = new TH3D("fhPtTrkTruePrimRec","fhPtTrkTruePrimRec",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkTruePrimRec);

      fhPtTrkSecOrFakeRec = new TH3D("fhPtTrkSecOrFakeRec","fhPtTrkSecOrFakeRec",100,0,100,20,-1,1, 10,0,10);
      fOutput->Add((TH3D*) fhPtTrkSecOrFakeRec);

      name = Form("JetPtPartLevel_Rho%s", rhotype.Data());
      fhJetPtPartLevelCorr = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins); //(270. -20, 250) added by KA
      fOutput->Add((TH1D*) fhJetPtPartLevelCorr);

      fhJetPtPartLevelZero = new TH1D("JetPtZeroPartLevel","JetPtZeroPartLevel", numberBins_jetpT_Zero, jetPtZero_Bins); //(250, 0, 250) added by KA
      fOutput->Add((TH1D*) fhJetPtPartLevelZero);

      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //FF Normalization of response matrix filled from recoil jets
         name = Form("RecoilJetPtPartLevel_Rho%s_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtPartLevelCorr[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtPartLevelCorr[itt]);

         //Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevel_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtZeroPartLevel[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtZeroPartLevel[itt]);

         //______________________________________________________
         //KA Normalization of response matrix filled from recoil jets on PARTICLE level
         name = Form("RecoilJetPt_Rho%s_PartLevel_Corresp_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtPartLevel_CorrespTT[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtPartLevel_CorrespTT[itt]);

         //Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevel_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtZeroPartLevel_CorrespTT[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtZeroPartLevel_CorrespTT[itt]);
      }

      name = Form("fhFractionOfSecInJet_Rho%s", rhotype.Data());
      fhFractionOfSecInJet = new TH2D(name.Data(), "Frac of jet pT carried by secondary tracks",50,0,50,210,0,1.05);
      fOutput->Add((TH2D*) fhFractionOfSecInJet);

      //1D unfolding
      name = Form("JetPtPartLevelVsJetPtDetLevel_Rho%s", rhotype.Data());
      fhJetPtPartLevelVsJetPtDetLevelCorr = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // (270, -20.0, 250, 270, -20.0, 270) added by KA
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr);

      fhJetPtZeroPartLevelVsJetPtZeroDetLevel = new TH2D("JetPtZeroPartLevel_Vs_JetPtZeroDetLevel","JetPtZeroPartLevel_Vs_JetPtZeroDetLevel", numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); //(250, 0, 250, 250, 0, 250) added by KA
      fOutput->Add((TH2D*) fhJetPtZeroPartLevelVsJetPtZeroDetLevel);

      name = Form("JetPtZeroPartLevel_Vs_JetPtDetLevel_Rho%s", rhotype.Data());
      fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); // (270. -20, 250, 250, 0, 250) added by KA
      fOutput->Add((TH2D*) fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr);

      name = Form("ImpurityInclusive_DetJetPtVsPartJetPt_Rho%s", rhotype.Data());
      fhImpurityInclusive_DetJetPtVsPartJetPtCorr = new TH2D(name.Data(), "Matched inclusive jets where part. level one outside TPC fid. cut", numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // (270, -20.0, 250, 270, -20.0, 270) added by KA
      fOutput->Add((TH2D*) fhImpurityInclusive_DetJetPtVsPartJetPtCorr);

      name = Form("Purity_InclusiveJets_Rho%s", rhotype.Data());
      fhPurity_InclusiveJets = new TH2D(name.Data(), name.Data(), ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
      fOutput->Add((TH2D*) fhPurity_InclusiveJets);

      fhNormalization_Purity = new TH2D("Normalization_Purity", "Normalization_Purity", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
      fOutput->Add((TH2D*) fhNormalization_Purity);

      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //FF Response matrix filled from recoil jets
	      name = Form("RecoilJetPtPartLevelVsJetPtDetLevel_Rho%s_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); //(270, - 20.0, 250, 270, -20.0, 270) added by FF
         fOutput->Add((TH2D*) fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt]);

         // Part level jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtDetLevel_Rho%s_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); //(270, - 20.0, 250, 250, 0.0, 250) added by FF
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt]);

         //  Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); //(250, 0.0, 250, 250, 0.0, 250) added by FF
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt]);

         //____________________________________________________________
         //KA Reponse matrix filled from recoil jets wrt corresponding TT
	      name = Form("RecoilJetPtPartLevelVsJetPtDetLevel_Rho%s_Corresp_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); //(270, - 20.0, 250, 270, -20.0, 270) added by KA
         fOutput->Add((TH2D*) fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt]);

         // Part level jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtDetLevel_Rho%s_Corresp_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); //(270, - 20.0, 250, 250, 0.0, 250) added by KA
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt]);

         //  Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); //(250, 0.0, 250, 250, 0.0, 250) added by KA
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt]);

         //____________________________________________________________
         name = Form("ImpurityRecoil_DetJetPtVsPartJetPt_Rho%s_TTH%d_%d", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]); //FF
         fhImpurityRecoil_DetJetPtVsPartJetPtCorr[itt] = new TH2D(name.Data(), "Matched recoil jets where part. level ones outside TPC fid. cut", numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); //(270, - 20.0, 250, 270, -20.0, 270) added by KA
         fOutput->Add((TH2D*) fhImpurityRecoil_DetJetPtVsPartJetPtCorr[itt]);
      }

      name = Form("fhJetPtResolutionVsPtPartLevel_Rho%s", rhotype.Data());
      fhJetPtResolutionVsPtPartLevel = new TH2D(name.Data(), name.Data(),100,0,100,50,0,2);
      fOutput->Add((TH2D*) fhJetPtResolutionVsPtPartLevel);
   }

   if(fMode == AliAnalysisTaskEA::kMC){

         const Int_t knmjdim = 5; //pt, Area, eta, phi, delta phi wrt TT
         Int_t   nmjbins[knmjdim]  = {  50, 20,   12, 20,               16   };
         Double_t nmjxmin[knmjdim] = {  0.,  0,  -0.6, 0,               TMath::Pi()-0.8};
         Double_t nmjxmax[knmjdim] = {100.,  1,   0.6,TMath::TwoPi(),   TMath::Pi()+0.8};

         for(Int_t itt=0; itt<fnHadronTTBins; itt++){  //fk
            name = Form("fhNotMatchedJetPt_%s_TTH%d_%d", trig[kMB].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhNotMatchedJetPt[itt] = new THnSparseF(name.Data(),"Not matched Recoil jet  pt area eta phi deltaphi", knmjdim, nmjbins, nmjxmin, nmjxmax);

            fOutput->Add((THnSparse*) fhNotMatchedJetPt[itt]);
         }
      }


   //2D unfolding -------------------------------
   if(fMode == AliAnalysisTaskEA::kMC){

      //Auxiliary variables for Sparse object

      // Particle Level jet pT is CORRECTED on Rhokt
      const Int_t fNumberOfDimensions = 4;
      const Int_t fSparseBinsNumber_JetPtRhokT[fNumberOfDimensions] = {ndeltaPhiBins, numberBins_jetpT_RhokT, ndeltaPhiBins, numberBins_jetpT_RhokT};

      //Inclusive jets
      // Missed events
      name = Form("MissedEvents_Phi_JetPt_Rho%s_PartLevel_InclusiveJets_%s", rhotype.Data(), trig[kMB].Data());
      fhPhi_JetPtPartLevel_InclusiveJets = new TH2D (name.Data(), "Missed events phi vs inclusive jet pT RhokT part level",  ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
      fOutput->Add((TH2D*) fhPhi_JetPtPartLevel_InclusiveJets);

      //Sparse object as basis for RooUnfoldResponse object
      name = Form("Phi_JetPtRhokT_DetLevel_Vs_Phi_JetPt_Rho%s_PartLevel_InclusiveJets_%s", rhotype.Data(), trig[kMB].Data());
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets = new THnSparseD (name.Data(), "Phi vs Jet pT RhokT for filling response matrix with inclusive jets", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
      fOutput->Add((THnSparse*) fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets);

      //Jets from events with TT
      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //Missed events for RM initialization, K.A.
         name = Form("MissedEvents_DeltaPhi_JetPt_Rho%s_PartLevel_TTH%d_%d_%s", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtPartLevel[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT RhokT part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtPartLevel[itt]);

         // Modified by KA
         name = Form("MissedEvents_DeltaPhi_JetPt_Rho%s_PartLevel_Corresp_TTH%d_%d_%s", rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtPartLevel_CorrespTT[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT RhokT part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtPartLevel_CorrespTT[itt]);

         //Sparse object as basis for RooUnfoldResponse object
         name = Form("DeltaPhi_JetPt_Rho%s_DetLevel_Vs_DeltaPhi_JetPt_Rho%s_PartLevel_TTH%d_%d_%s", rhotype.Data(), rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT RhokT for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]);

         // Modified by KA
         name = Form("DeltaPhi_JetPt_Rho%s_DetLevel_Vs_DeltaPhi_JetPt_Rho%s_PartLevel_Corresp_TTH%d_%d_%s", rhotype.Data(), rhotype.Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT RhokT for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]);

         //Count TT on Particle level
         // name = Form("fhNumberOf_ChosenTTH%d_%d_%s_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         // fhNumberOf_ChosenTT_PartLevel[itt] = new TH1D (name.Data(), "Number of Chosen TT Particle Level", 100, fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         // fOutput->Add((TH1D*) fhNumberOf_ChosenTT_PartLevel[itt]);
      }


      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      // Particle Level jet pT is NOT CORRECTED on Rhokt (added by KA)
      const Int_t fSparseBinsNumber_PartLevel_JetPtZero[fNumberOfDimensions] = {ndeltaPhiBins, numberBins_jetpT_RhokT, ndeltaPhiBins, numberBins_jetpT_Zero};
      const Int_t fSparseBinsNumber_JetPtZero[fNumberOfDimensions]           = {ndeltaPhiBins, numberBins_jetpT_Zero, ndeltaPhiBins, numberBins_jetpT_Zero};

      //Inclusive jets
      //Missed events
      name = Form("MissedEvents_Phi_JetPtZeroPartLevel_InclusiveJets_%s", trig[kMB].Data());
      fhPhi_JetPtZeroPartLevel_InclusiveJets = new TH2D (name.Data(), "Missed events phi vs inclusive jet pT zero part level",  ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_Zero, jetPtZero_Bins); // added by KA
      fOutput->Add((TH2D*) fhPhi_JetPtZeroPartLevel_InclusiveJets);

      // Sparse object as basis for RooUnfoldResponse object. Part level jet pT is not corrected on Rhokt
      name = Form("Phi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets_%s", trig[kMB].Data());
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets = new THnSparseD (name.Data(), "Phi vs Jet pT zero for filling response matrix with inclusive jets", fNumberOfDimensions, fSparseBinsNumber_PartLevel_JetPtZero, NULL, NULL); // added by KA
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(3, jetPtZero_Bins);  //Jet pT particle level
      fOutput->Add((THnSparse*) fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets);

      // Sparse object as basis for RooUnfoldResponse object. Jet pT is not corrected on Rhokt
      name = Form("Phi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets_%s", trig[kMB].Data());
      fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets = new THnSparseD (name.Data(), "Phi vs Jet pT zero for filling response matrix with inclusive jets", fNumberOfDimensions, fSparseBinsNumber_JetPtZero, NULL, NULL); // added by KA
      fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(0, deltaPhiBins);   //Delta phi detector level
      fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(1, jetPtZero_Bins); //Jet pT detector level
      fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(2, deltaPhiBins);   //Delta phi particle level
      fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets->SetBinEdges(3, jetPtZero_Bins); //Jet pT particle level
      fOutput->Add((THnSparse*) fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets);

      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //Missed events for RM initialization, K.A. Part level jet pT is not corrected on RhokT
         name = Form("MissedEvents_DeltaPhi_JetPtZero_PartLevel_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtZero_PartLevel[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT zero part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_Zero, jetPtZero_Bins); // added by KA
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtZero_PartLevel[itt]);

         //Sparse object as basis for RooUnfoldResponse object. Part level jet pT is NOT corrected on RhokT
         name = Form("DeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT zero for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_PartLevel_JetPtZero, NULL, NULL); // added by KA
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(3, jetPtZero_Bins);  //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]);

         //Sparse object as basis for RooUnfoldResponse object. Jet pT is NOT corrected on RhokT
         name = Form("DeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT zero for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtZero, NULL, NULL); // added by KA
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(0, deltaPhiBins);   //Delta phi detector level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(1, jetPtZero_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(2, deltaPhiBins);   //Delta phi particle level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]->SetBinEdges(3, jetPtZero_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[itt]);

         //_________________________________________________________
         // Response matrix constructed w.r.t. corresponding TT

         //Missed events for RM initialization, K.A. Part level jet pT is not corrected on RhokT
         name = Form("MissedEvents_DeltaPhi_JetPtZero_PartLevel_Corresp_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT zero part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_Zero, jetPtZero_Bins); // added by KA
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]);

         //Sparse object as basis for RooUnfoldResponse object. Part level jet pT is NOT corrected on RhokT
         name = Form("DeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_Corresp_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT zero for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_PartLevel_JetPtZero, NULL, NULL); // added by KA
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(3, jetPtZero_Bins);  //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]);

         //Sparse object as basis for RooUnfoldResponse object. Jet pT is NOT corrected on RhokT
         name = Form("DeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_Corresp_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT zero for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtZero, NULL, NULL); // added by KA
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(0, deltaPhiBins);   //Delta phi detector level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(1, jetPtZero_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(2, deltaPhiBins);   //Delta phi particle level
         fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]->SetBinEdges(3, jetPtZero_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[itt]);
      }


      //--------------------------------------------------------------------------------------
      // Jet energy and angular resolution
      constexpr int numOfDim = 5; //Mult, phi resol, phi, pT resol, pT
      vector<double> mult_axis {0.0, 4.0, 4.4, 5.0, 9.0, 15.0};
      const int nBins = mult_axis.size() - 1;
      vector<int> numOfBins {nBins, 40, 40, 50, 100};

      fhJetEnergy_Angle_Resol = new THnSparseD("JetEnergy_Angle_Resol","JetEnergy_Angle_Resol", numOfDim, &numOfBins[0], nullptr, nullptr);
      fhJetEnergy_Angle_Resol->GetAxis(0)->Set(numOfBins[0], &mult_axis[0]);
      fhJetEnergy_Angle_Resol->GetAxis(1)->Set(numOfBins[1], -1.0, 1.0);
      fhJetEnergy_Angle_Resol->GetAxis(2)->Set(numOfBins[2], 0.0, 2*TMath::Pi());
      fhJetEnergy_Angle_Resol->GetAxis(3)->Set(numOfBins[3], -1.0, 1.5);
      fhJetEnergy_Angle_Resol->GetAxis(4)->Set(numOfBins[4], 0.0, 100.0);
      fOutput->Add((THnSparse*) fhJetEnergy_Angle_Resol);

      //--------------------------------------------------------------------------------------
      // Inclusive response matrix in HM events
      mult_axis = {4.0, 4.4, 5.0, 9.0};
      numOfBins = {(int) mult_axis.size() - 1, ndeltaPhiBins, numberBins_jetpT_RhokT, ndeltaPhiBins, numberBins_jetpT_RhokT};

      fhInclusiveRespMatrix_HM = new THnSparseD("InclusiveRespMatrix_HM","InclusiveRespMatrix_HM", numOfDim, &numOfBins[0], nullptr, nullptr);
      fhInclusiveRespMatrix_HM->SetBinEdges(0, &mult_axis[0]);
      fhInclusiveRespMatrix_HM->SetBinEdges(1, deltaPhiBins);
      fhInclusiveRespMatrix_HM->SetBinEdges(2, jetPtRhokT_Bins);
      fhInclusiveRespMatrix_HM->SetBinEdges(3, deltaPhiBins);
      fhInclusiveRespMatrix_HM->SetBinEdges(4, jetPtRhokT_Bins);
      fOutput->Add((THnSparse*) fhInclusiveRespMatrix_HM);

      // Missed events
      fhMissedEvents_HM = new TH3D ("MissedEvents_InclusiveJets_HM", "MissedEvents_InclusiveJets_HM", numOfBins[0], &mult_axis[0], ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
      fOutput->Add((TH3D*) fhMissedEvents_HM);

   }
   //2D unfolding -------------------------------


   //JET PT ASYMMETRY
//   for(Int_t itg=kMB; itg<=kHM; itg++){
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//         name = Form("fhJetPtAsymmetryCB_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhJetPtAsymmetryCB[itg][itt] = new TH2D(name.Data(),name.Data(),10,0,10, 21,-1,1.1);
//         fOutput->Add((TH2D*) fhJetPtAsymmetryCB[itg][itt]); //ReMx for detector level TTH
//      }
//   }

//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//         name = Form("fhJetPtAsymmetryCB_MB_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhJetPtAsymmetryCBPartLevel[itt] = new TH2D(name.Data(),name.Data(),10,0,10, 21,-1,1.1);
//         fOutput->Add((TH2D*) fhJetPtAsymmetryCBPartLevel[itt]); //ReMx for detector level TTH
//      }
//   }

   //TRACK PT ASYMMETRY
//   for(Int_t itg=kMB; itg<=kHM; itg++){
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//         name = Form("fhTrackPtAsymmetryCB_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhTrackPtAsymmetryCB[itg][itt] = new TH2D(name.Data(),name.Data(),10,0,10, 21,-1,1.1);
//         fOutput->Add((TH2D*) fhTrackPtAsymmetryCB[itg][itt]); //ReMx for detector level TTH
//      }
//   }


//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//         name = Form("fhTrackPtAsymmetryCB_MB_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//         fhTrackPtAsymmetryCBPartLevel[itt] = new TH2D(name.Data(),name.Data(),10,0,10, 21,-1,1.1);
//         fOutput->Add((TH2D*) fhTrackPtAsymmetryCBPartLevel[itt]); //ReMx for detector level TTH
//      }
//   }


   //Auxiliary jet pT spectra filled event by event which will not go to output
   fhJetPtEvtByEvent       = new TH1D("fhJetPtEvtByEvent","fhJetPtEvtByEvent", 100,0,100);
//   fhRecoilJetPtEvtByEventRandomTT = new TH1D("fhRecoilJetPtEvtByEventRandomTT","fhRecoilJetPtEvtByEventRandomTT", 100,0,100);

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      fhJetPtEvtByEventPartLevel = new TH1D("fhJetPtEvtByEventPartLevel","fhJetPtEvtByEventPartLevel", 100,0,100);
//      fhRecoilJetPtEvtByEventRandomTTPartLevel = new TH1D("fhRecoilJetPtEvtByEventRandomTTPartLevel","fhRecoilJetPtEvtByEventRandomTTPartLevel", 100,0,100);
   }

   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      name = Form("fhRecoilJetPtEvtByEvent_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      fhRecoilJetPtEvtByEvent[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
   }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
         name = Form("fhRecoilJetPtEvtByEvent_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtEvtByEventPartLevel[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
      }
   }

   //Count multiplicity of high-pT jets per event
   const Int_t khighptjetdim = 3;
   Int_t    highptjetbins[khighptjetdim] = {40, 50, 10};
   Double_t highptjetxmin[khighptjetdim] = { 0.,  0,  0.};
   Double_t highptjetxmax[khighptjetdim] = {40., 50, 10.};

//   for(Int_t itg=kMB; itg<=kHM; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;
//
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//
//         name = Form("fhNumberOfHighPtJetsCB_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//
//         fhNumberOfHighPtJetsCB[itg][itt] = new  THnSparseF(name.Data(),"Number of jets with pt above X", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
//         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsCB[itg][itt]);
//      }
//   }

//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
//
//         name = Form("fhNumberOfHighPtJetsCB_MB_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
//
//         fhNumberOfHighPtJetsCBPartLevel[itt] = new  THnSparseF(name.Data(),"Number of jets with pt above X particle level", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
//         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsCBPartLevel[itt]);
//      }
//   }


   //Count multiplicity of recoil high-pT jets per event
   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
         name = Form("fhNumberOfHighPtJetsRecoil_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);

         fhNumberOfHighPtJetsRecoil[itg][itt] = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoil[itg][itt]);
      }
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if(fMode == AliAnalysisTaskEA::kKine) continue;
      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;


      // For TTsig number of recoi jets with pt 10-15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPt10_15_%s", trig[itg].Data());
      fhNjetsPt10_15[itg] = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);  
      fOutput->Add((TH2D*) fhNjetsPt10_15[itg]); 

      // For TTsig number of recoi jets with pt above 15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPtAbove15_%s",trig[itg].Data());    
      fhNjetsPtAbove15[itg] = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);    
      fOutput->Add((TH2D*) fhNjetsPtAbove15[itg]);    
   }


//   for(Int_t itg=kMB; itg<=kHM; itg++){
//      if(fMode == AliAnalysisTaskEA::kKine) continue;
//      if((fMode == AliAnalysisTaskEA::kMC) && itg == kHM) continue;

 //     name = Form("fhNumberOfHighPtJetsRecoil_%s_randomTTH", trig[itg].Data());
//      fhNumberOfHighPtJetsRecoilRandomTT[itg] = new THnSparseF(name.Data(),"Number of recoil jets with pt above X (random TT)", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
 //     fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoilRandomTT[itg]);
 //  }

   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){  //  V0Mnorm from particle level
         name = Form("fhNumberOfHighPtJetsRecoil_MB_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);

         fhNumberOfHighPtJetsRecoilPartLevel[itt] = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X particle level", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoilPartLevel[itt]);
      }

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){  //FILIP    V0Mnorm from detector level
         name = Form("fhNumberOfHighPtJetsRecoil_MB_TTH%d_%d_PartLevel_V0MnormDetLev", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhNumberOfHighPtJetsRecoilPartLevelV0MnormDetLev[itt] = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X particle level", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoilPartLevelV0MnormDetLev[itt]);
      }//FILIP

      // For TTsig number of recoi jets with pt 10-15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPt10_15_PartLevel");
      fhNjetsPt10_15_PartLevel = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);  
      fOutput->Add((TH2D*) fhNjetsPt10_15_PartLevel); 

      // For TTsig number of recoi jets with pt above 15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPtAbove15_PartLevel");    
      fhNjetsPtAbove15_PartLevel = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);    
      fOutput->Add((TH2D*) fhNjetsPtAbove15_PartLevel);   

      // For TTsig number of recoi jets with pt 10-15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPt10_15_PartLevelV0MnormDetLev");
      fhNjetsPt10_15_PartLevelV0MnormDetLev = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);  
      fOutput->Add((TH2D*) fhNjetsPt10_15_PartLevelV0MnormDetLev); 

      // For TTsig number of recoi jets with pt above 15 GeV/c vs V0Mnorm
      name = Form("fhNjetsPtAbove15_PartLevelV0MnormDetLev");    
      fhNjetsPtAbove15_PartLevelV0MnormDetLev = new TH2D(name.Data(),name.Data(),20,0,20,10,0,10);    
      fOutput->Add((TH2D*) fhNjetsPtAbove15_PartLevelV0MnormDetLev);    

   }

//   if(fMode == AliAnalysisTaskEA::kMC || fMode == AliAnalysisTaskEA::kKine){
//      name = Form("fhNumberOfHighPtJetsRecoil_MB_randomTTH_PartLevel");
//      fhNumberOfHighPtJetsRecoilRandomTTPartLevel = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X particle level (random TT)", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
//      fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoilRandomTTPartLevel);
//  }
   //+++++++++++++++++++++++++++ EMBEDDING +++++++++++++++++++++++
   if(fMode == AliAnalysisTaskEA::kEmbedding){


      for(Int_t itg=kMB; itg<=kHM; itg++){   //@@@
         //remx normalization spectra
         name = Form("fhJetPtPartLevelCorr_EMB_%s_Rho%s",trig[itg].Data(),rhotype.Data());
         fhJetPtPartLevelCorr_EMB[itg] = new TH1D(name.Data(), name.Data(), 270, -20, 250);
         fOutput->Add((TH1D*) fhJetPtPartLevelCorr_EMB[itg]);

         name = Form("fhJetPtPartLevelZero_EMB_%s",trig[itg].Data());
         fhJetPtPartLevelZero_EMB[itg] = new TH1D(name.Data(), name.Data(), 250, 0, 250);
         fOutput->Add((TH1D*) fhJetPtPartLevelZero_EMB[itg]);
      }

      //remx
      for(Int_t itg=kMB; itg<=kHM; itg++){   //@@@
         name = Form("fhJetPtPartLevelVsJetPtDetLevelCorr_EMB_%s_Rho%s",trig[itg].Data(), rhotype.Data());
         fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg] = new TH2D(name.Data(), name.Data(), 270, -20, 250, 270, -20, 250);
         fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[itg]);

         name = Form("fhJetPtPartLevelVsJetPtDetLevelZero_EMB_%s_Rho%s",trig[itg].Data(), rhotype.Data());
         fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg] = new TH2D(name.Data(), name.Data(), 270, -20, 250, 250, 0, 250);
         fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelZero_EMB[itg]);

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


 //     for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
 //        name = Form("fhTrialsEMB_%s", trig[itg].Data());
 //        fhTrialsEMB[itg] = new TH1F(name.Data(), name.Data(),  nPtHardBins, 0, nPtHardBins);
 //        fOutput->Add((TH1F*) fhTrialsEMB[itg]);
 //     }
 //
 //     for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
 //        name = Form("fhXsectionEMB_%s", trig[itg].Data());
 //        fhXsectionEMB[itg] = new TProfile(name.Data(), name.Data(),  nPtHardBins, 0, nPtHardBins);
 //        fOutput->Add((TProfile*) fhXsectionEMB[itg]);
 //     }

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         name = Form("fhPtHardEMB_%s", trig[itg].Data());
         //fhPtHardEMB[itg] = new TH1F(name.Data(), name.Data(), fNbins*2, fMinBinPt, fMaxBinPt*4);
         fhPtHardEMB[itg] = new TH1F(name.Data(), name.Data(), 1000, 0, 1000);
         fOutput->Add((TProfile*) fhPtHardEMB[itg]);
      }


   }

   //+++++++++++++++++++++++ MOEMENTU SMEARING HISTOGRAMS  +++++++++++++++++++++++++++++++++++++++++++++++
   if(fMode != AliAnalysisTaskEA::kKine){
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
   }

   if(fMode == AliAnalysisTaskEA::kMC){ //FILIP
      Double_t NjetBins [] = {0,1,2,3,4,20};  
      Int_t nNjetBins = sizeof(NjetBins)/sizeof(Double_t) - 1;  
      Double_t NV0MnormBins [] = {0,4,5,9,20};  
      Int_t nNV0MnormBins = sizeof(NV0MnormBins)/sizeof(Double_t) - 1;  
		   
      fhNjetReMx_V0MnormDetLev_15GeV = new TH3D("fhNjetReMx_V0MnormDetLev_15GeV","fhNjetReMx_V0MnormDetLev_15GeV",
        nNjetBins,NjetBins,nNjetBins,NjetBins,nNV0MnormBins,NV0MnormBins);

      fOutput->Add((TH3D*) fhNjetReMx_V0MnormDetLev_15GeV);

      fhNjetNorm_V0MnormDetLev_15GeV = new TH2D("fhNjetNorm_V0MnormDetLev_15GeV","fhNjetNorm_V0MnormDetLev_15GeV", nNjetBins,NjetBins,nNV0MnormBins,NV0MnormBins);
      fOutput->Add((TH2D*) fhNjetNorm_V0MnormDetLev_15GeV);

      fhNjetReMx_V0MnormDetLev_20GeV = new TH3D("fhNjetReMx_V0MnormDetLev_20GeV","fhNjetReMx_V0MnormDetLev_20GeV",
        nNjetBins,NjetBins,nNjetBins,NjetBins,nNV0MnormBins,NV0MnormBins);
      fOutput->Add((TH3D*) fhNjetReMx_V0MnormDetLev_20GeV);

      fhNjetNorm_V0MnormDetLev_20GeV = new TH2D("fhNjetNorm_V0MnormDetLev_20GeV","fhNjetNorm_V0MnormDetLev_20GeV", nNjetBins,NjetBins,nNV0MnormBins,NV0MnormBins);
      fOutput->Add((TH2D*) fhNjetNorm_V0MnormDetLev_20GeV); 
   }//FILIP


   if(fMode == AliAnalysisTaskEA::kMC){ 
   //NEW
      Double_t myNV0MnormBins [] = {0,4,4.1,4.2,4.3,4.4,4.5,5,9,20};  //  V0M/<V0M> bins 
      Int_t nmyNV0MnormBins = sizeof(myNV0MnormBins)/sizeof(Double_t) - 1;  
	
      Double_t myTTPtBins [] = {0,6,7,10,15,20,30,100}; //trigger track pT bins 
      Int_t nmyTTPtBins = sizeof(myTTPtBins)/sizeof(Double_t) - 1;  

      Double_t myJetPtBins [] = {0,5,10,15,20,25,30,35,40,45,50,60,80,100}; //trigger track pT bins 
      Int_t nmyJetPtBins = sizeof(myJetPtBins)/sizeof(Double_t) - 1;  

      const Int_t nmyPhiBins = 800; 
      Double_t myPhiBins[nmyPhiBins+1];
      Double_t p = TMath::Pi()/nmyPhiBins;
      for(Int_t i = 0; i <= nmyPhiBins; i++) myPhiBins[i] = -TMath::Pi()/2 + i*p;  // range was (-pi/2,pi/2)


      fPhiSmearingTrack = new TH3D("fPhiSmearingTrack","fPhiSmearingTrack", nmyNV0MnormBins, myNV0MnormBins, nmyTTPtBins, myTTPtBins, nmyPhiBins, myPhiBins); 
      fOutput->Add((TH3D*) fPhiSmearingTrack); 

      fPhiSmearingJet = new TH3D("fPhiSmearingJet","fPhiSmearingJet", nmyNV0MnormBins, myNV0MnormBins,  nmyJetPtBins,  myJetPtBins, nmyPhiBins, myPhiBins); 
      fOutput->Add((TH3D*) fPhiSmearingJet); 

      fDeltaPhiSmearingTT2030Jet = new TH3D("fDeltaPhiSmearingTT2030Jet","fDeltaPhiSmearingTT2030Jet", nmyNV0MnormBins, myNV0MnormBins,  nmyJetPtBins,  myJetPtBins, nmyPhiBins, myPhiBins); 
      fOutput->Add((TH3D*) fDeltaPhiSmearingTT2030Jet); 
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

//   if(!fCaloCells){
//      if(fCaloCellsName.IsNull()){
//         fCaloCells = InputEvent()->GetEMCALCells();
//      }else{
//         fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
//         if(!fCaloCells) AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data()));
//      }
//      cout<<"load calo cells"<<endl;
//   }



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
