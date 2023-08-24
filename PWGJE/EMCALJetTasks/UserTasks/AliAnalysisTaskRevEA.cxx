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
#include "AliAnalysisTaskRevEA.h"
#include "AliHeader.h"
#include "AliRunLoader.h"
#include "AliVVZERO.h"
#include "AliAODZDC.h"
#include "AliVZDC.h"
#include "AliMultSelection.h"

//#include "AliEmcalDownscaleFactorsOCDB.h"
//#include "AliEmcalAnalysisFactory.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskRevEA)

using namespace PWGJE::EMCALJetTasks;
using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN PP 13 TeV
// Author Filip Krizek   (8.Aug. 2019)

//________________________________________________________________________________________

AliAnalysisTaskRevEA::AliAnalysisTaskRevEA():
AliAnalysisTaskEmcalJet("AliAnalysisTaskRevEA", kTRUE),
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyDetLevelContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fKTJetContainerDetLevel(0x0),
fKTJetContainerPartLevel(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsHighMultTrig(0),
fCentralityV0M(-1),
fMultV0Mnorm(0.),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0),
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhRhoMBpart(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),                         //1D unfolding
fhJetPtZeroPartLevelVsJetPtZeroDetLevel(0x0),                         //1D unfolding
fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr(0x0),                   //1D unfolding (added by KA)
fhPhi_JetPtPartLevel_InclusiveJets(0x0),                          //2D unfolding
fhPhi_JetPtZeroPartLevel_InclusiveJets(0x0),                      //2D unfolding (added by KA)
fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets(0x0),     //2D unfolding
fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhJetPtResolutionVsPtPartLevel(0x0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fMode(AliAnalysisTaskRevEA::kNormal),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.),
fRho(0.),
fRhoMC(0.),
fMaxFacPtHard(0),
fhNjetReMx_V0MnormDetLev_15GeV(0),
fhNjetNorm_V0MnormDetLev_15GeV(0),
fhNjetReMx_V0MnormDetLev_20GeV(0),
fhNjetNorm_V0MnormDetLev_20GeV(0)
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
   }

   for(Int_t iq = 0; iq < 4; ++iq){
      fArray_for_filling[iq] = 0.;
   }
   /////////////////



   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;

      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      fhRho[itg] = 0x0;

      for(Int_t i=0; i<fkTTbins; i++){
         fhRhoTTH[itg][i]=0x0;
      }

   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT_Labels[i].resize(0); //Modified by KA

      //TT
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhMultTTH[itg][i] = 0x0;

         fhTTH_V0Mnorm[itg][i] = 0x0;
      }

      fhTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhTT_Corresp[i] = 0x0; //KA

      //RECOIL JET SPECTRA
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhRecoilJetPtTTH_V0Mnorm[itg][i] = 0x0;
         fhRecoilJetPhiTTH_V0Mnorm[itg][i] = 0x0;
      }

      fhRecoilJetPtTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[i] = NULL; //added by KA

      fhRecoilJetPhiTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[i] = NULL; //added by KA


      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhDeltaPtTTH_RC_V0Mnorm[itg][i] = 0x0;
      }
   }


   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhCentrality[itg]= 0x0;
      for(Int_t ic=0; ic<fkCE;ic++){
         fhSignal[itg][ic] = 0x0;

         for(Int_t i=0; i<fkTTbins;i++){
            fhSignalTTH[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0;

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
      }
   }


   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMBpart[i]=0x0;
   }


   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTH[i] = -1;
      fIndexTTH_PartLevel[i] = -1;

      fTTH[i].resize(0);
      fTTH_PartLevel[i].resize(0);

      fdeltapT[i]  = 0.;
   }

   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }



   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhNumberOfHighPtJetsRecoil[itg][itt] = NULL;
      }
      fhRecoilJetPtEvtByEvent[itt] = NULL;

      fhNumberOfHighPtJetsRecoilPartLevel[itt] = NULL;
      fhRecoilJetPtEvtByEventPartLevel[itt] = NULL;
   }


}

//________________________________________________________________________
AliAnalysisTaskRevEA::AliAnalysisTaskRevEA(const char *name):
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyDetLevelContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName(""),
fTrkContainerDetLevel(0x0),
fParticleContainerPartLevel(0x0),
fJetContainerDetLevel(0x0),
fJetContainerPartLevel(0x0),
fKTJetContainerDetLevel(0x0),
fKTJetContainerPartLevel(0x0),
fMultSelection(0x0),
fIsMinBiasTrig(0),
fIsHighMultTrig(0),
fCentralityV0M(-1),
fMultV0Mnorm(0.),
fTrackEtaWindow(0.9),
fMinTrackPt(0.150),
fHelperClass(0),
fInitializedLocal(0),
fHistEvtSelection(0x0),
fhVertexZall(0x0),
fhVertexZ(0x0),
fhRhoMBpart(0x0),
fhJetPtPartLevelCorr(0x0),
fhJetPtPartLevelZero(0x0),
fhJetPtPartLevelVsJetPtDetLevelCorr(0x0),                         //1D unfolding
fhJetPtZeroPartLevelVsJetPtZeroDetLevel(0x0),                         //1D unfolding
fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr(0x0),                   //1D unfolding (added by KA)
fhPhi_JetPtPartLevel_InclusiveJets(0x0),                          //2D unfolding
fhPhi_JetPtZeroPartLevel_InclusiveJets(0x0),                      //2D unfolding (added by KA)
fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets(0x0),     //2D unfolding
fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets(0x0), //2D unfolding (added by KA)
fhJetPtResolutionVsPtPartLevel(0x0),
fZVertexCut(10.0),
fnHadronTTBins(0),
fMode(AliAnalysisTaskRevEA::kNormal),
fFillSigTT(1),
fPhiCut(TMath::Pi()-0.6),
fRandom(0),
fJetR(0.4),
fJetAcut(0.),
fRho(0.),
fRhoMC(0.),
fMaxFacPtHard(0),
fhNjetReMx_V0MnormDetLev_15GeV(0),
fhNjetNorm_V0MnormDetLev_15GeV(0),
fhNjetReMx_V0MnormDetLev_20GeV(0),
fhNjetNorm_V0MnormDetLev_20GeV(0)
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
   }

   for(Int_t iq = 0; iq < 4; ++iq){
      fArray_for_filling[iq] = 0.;
   }
   /////////////////



   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhTrackPhiIncl[itg]=0x0;
      fhTrackEtaIncl[itg]=0x0;
      fhJetPhiIncl[itg]=0x0;
      fhJetEtaIncl[itg]=0x0;

      fhRho[itg] = 0x0;

      for(Int_t i=0; i<fkTTbins; i++){
         fhRhoTTH[itg][i]=0x0;
      }
   }

   for(Int_t i=0; i<fkTTbins; i++){
      fHadronTT_Labels[i].resize(0); // Modified by KA
      fHadronTT_Labels_PartLevel[i].resize(0); // Modified by KA

      //TT
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhMultTTH[itg][i] = 0x0;
         fhTTH_V0Mnorm[itg][i] = 0x0;
      }

      fhTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhTT_Corresp[i] = 0x0; //KA

      //RECOIL JET SPECTRA
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhRecoilJetPtTTH_V0Mnorm[itg][i] = 0x0;
         fhRecoilJetPhiTTH_V0Mnorm[itg][i] = 0x0;
      }

      fhRecoilJetPtTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[i] = NULL; //added by KA

      fhRecoilJetPhiTTH_V0Mnorm_PartLevel[i] = 0x0;
      fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[i] = NULL; //added by KA

      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhDeltaPtTTH_RC_V0Mnorm[itg][i] = 0x0;
      }
   }



   for(Int_t i=0; i<fkTTbins;i++){
      fHadronTTLowPt[i]=-1;
      fHadronTTHighPt[i]=-1;
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      fhCentrality[itg] = 0x0;
      for(Int_t ic=0; ic<fkCE;ic++){
         fhSignal[itg][ic] = 0x0;

         for(Int_t i=0; i<fkTTbins;i++){
            fhSignalTTH[itg][ic][i] = 0x0;
         }
      }
   }

   //particle level
   for(Int_t ic=0; ic<fkCE;ic++){
      fhSignal_PartLevel[ic] = 0x0;

      for(Int_t i=0; i<fkTTbins;i++){
         fhSignalTTH_PartLevel[ic][i] = 0x0;
      }
   }


   for(Int_t i=0; i<fkTTbins;i++){
      fhRhoTTHinMBpart[i]=0x0;
   }


   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTH[i] = -1;
      fIndexTTH_PartLevel[i] = -1;

      fTTH[i].resize(0);
      fTTH_PartLevel[i].resize(0);

      fdeltapT[i]  = 0.;
   }

   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }

   //JET AND TRACK PT ASYMMETRY
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      for(Int_t itg=kMB; itg<=kHM; itg++){
         fhNumberOfHighPtJetsRecoil[itg][itt] = NULL;
      }
      fhRecoilJetPtEvtByEvent[itt] = NULL;
      fhNumberOfHighPtJetsRecoilPartLevel[itt] = NULL;
      fhRecoilJetPtEvtByEventPartLevel[itt] = NULL;
   }


   DefineOutput(1, TList::Class());
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

//_____________________________________________________________________________________
AliAnalysisTaskRevEA*  AliAnalysisTaskRevEA::AddTaskRevEA(
  Int_t       mode,
  const char* jetarrayname,
  const char* jetarraynamePartMC,
  const char* trackarrayname,
  const char* mcpariclearraynamePartMC,
  const char* ktjetarrayname,
  const char* ktjetarraynamePartMC,
  Double_t    jetRadius,
  UInt_t      trigger,
  Double_t    trackEtaWindow,
  Bool_t      useVertexCut,
  Bool_t      usePileUpCut,
  Double_t    acut,
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
      ::Error("AliAnalysisTaskRevEA.cxx", "No analysis manager to connect to.");
      return NULL;
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString myContName("");
   myContName = Form("JetAnalysisR%02d_Acut%02d", TMath::Nint(jetRadius*10), TMath::Nint(acut*10));
   myContName.Append(suffix);

   AliAnalysisTaskRevEA *task = new AliAnalysisTaskRevEA(myContName.Data());

   if(mode == AliAnalysisTaskRevEA::kMC){  
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

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks (or combined tracks if embedding)
   trackCont->SetMinPt(0.15);
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);

   if(mode == AliAnalysisTaskRevEA::kMC){
      trackContTrue = task->AddMCParticleContainer(mcpariclearraynamePartMC); //particle level MC particles
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(1e-3); // KA: old value was 0.15 GeV/c
      trackContTrue->SetEtaLimits(-5.1,5.1); //V0 eta range
   }

   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //AKT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrue   = 0x0; //AKT jet container with mc particle level jets pythia

   AliJetContainer *jetContRecKT  = 0x0; //KT jet container with detector level tracks   or combined event jets after embedding
   AliJetContainer *jetContTrueKT = 0x0; //KT jet container with mc particle level jets pythia

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
      jetContRecKT->SetMinPt(0.);
      jetContRecKT->SetMaxTrackPt(100.0);
      jetContRecKT->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
   }

   if(mode == AliAnalysisTaskRevEA::kMC){
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
         jetContTrueKT->SetMaxTrackPt(1000.0); 
         jetContTrueKT->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
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

   task->SetJetContainerName(jetarrayname);
   task->SetMCPartJetContainerName(jetarraynamePartMC);

   task->SetKTJetContainerName(ktjetarrayname);
   task->SetKTMCPartJetContainerName(ktjetarraynamePartMC);

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
Bool_t AliAnalysisTaskRevEA::PassedMinBiasTrigger(){
  //minimum bias trigger

  bool passedTrigger = kFALSE;

  if(fMode == AliAnalysisTaskRevEA::kMC){
     //mc simulation emulate V0 coincidence trigger
    
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
Bool_t AliAnalysisTaskRevEA::PassedHighMultTrigger(){
   //high multiplicity V0M trigger

   if(fMode == AliAnalysisTaskRevEA::kMC)   return kFALSE; //MC

   bool passedTrigger = kFALSE;
   UInt_t triggerMask = fInputHandler->IsEventSelected();
   if(triggerMask & AliVEvent::kHighMultV0){
      passedTrigger = kTRUE;
   }

   return passedTrigger;
}
//_____________________________________________________________________________________
Double_t AliAnalysisTaskRevEA::GetMyRho(AliJetContainer* ktjets){

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

      Int_t nJetAcckt = 0;

      for(auto jetIterator : ktjets->accepted_momentum() ){
                      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
          jet = jetIterator.second;  // Get the pointer to jet object
          if(!jet)  continue;

          if(jet==jetLJ) continue; //skip two leading kT jets
          if(jet==jetSJ) continue;

          //standard area based approach
	  if(jet->Area() > 0.2){
             frhovec[nJetAcckt]  = jet->Pt()/jet->Area();
             nJetAcckt++;
	  }
      }

      if(nJetAcckt>0){
         myrho = TMath::Median(nJetAcckt, frhovec);
      }


   return  myrho;
}
//________________________________________________________________________

Bool_t AliAnalysisTaskRevEA::IsEventInAcceptance(AliVEvent* event){
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

      if(event->IsPileupFromSPDInMultBins()){
         fHistEvtSelection->Fill(2.5); //count events rejected by pileup
         return kFALSE;
      }
   }
   //___________________________________________________
   //MULTIPLICITY SELECTIONS

   if(fMode == AliAnalysisTaskRevEA::kNormal){
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
Bool_t AliAnalysisTaskRevEA::IsTrackInAcceptance(AliVParticle* track, Bool_t isGen){
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
void AliAnalysisTaskRevEA::ExecOnceLocal(){
   // Initialization of jet containers done in  AliAnalysisTaskEmcalJet::ExecOnce()
   //Read arrays of jets and tracks
   fInitializedLocal = kTRUE;

   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm

   return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRevEA::IsOutlier(){ //FK// whole function
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
          AliWarning("AliAnalysisTaskRevEA MC header not found");
          fHistEvtSelection->Fill(9.5);
          return kTRUE; //skip the event
      }

      return kFALSE;  //there was not outlier all jets have pT below fMaxFacPtHard * ptHard

   } else {

      fHistEvtSelection->Fill(9.5);
      AliWarning("AliAnalysisTaskRevEA MC header not found");
      return kTRUE; //MC header not found
   }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskRevEA::FillHistograms(){
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


   AliGenEventHeader* mcHeader = NULL;
   AliAODMCHeader* aodMCH = NULL;

   //+++++++++++++++++++++++++++++ check MC z vertex position ++++++++++++++++++++++++++
   if(fMode == AliAnalysisTaskRevEA::kMC){
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

   //cut on vertex z in particle level
   if(mcHeader){
      TArrayF pyVtx;
      mcHeader->PrimaryVertex(pyVtx);
      if(TMath::Abs(pyVtx[2]) > fZVertexCut) return kTRUE; //skip events with particle level vertex out of +-10 cm
   }

   //cut on vertex z in detector level
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec

   //_________________________________________________________________
   // GET V0M/<V0M> from centrality framework use the same number for detector and particle level
   fCentralityV0M = -1;
   fMultV0Mnorm   = -1;
   
   fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(fMultSelection){
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      fMultV0Mnorm   =  fMultSelection->GetZ("V0M");
   }else{
      return kFALSE;
   }

   //_________________________________________________________________
   // INITIALIZATION
   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTH[i] = -1;
      fIndexTTH_PartLevel[i] = -1;

      fTTH[i].resize(0);
      fTTH_PartLevel[i].resize(0);

      fdeltapT[i]  = 0.;
   }


   //_________________________________________________________________
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){
      fhRecoilJetPtEvtByEvent[itt]->Reset();
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         fhRecoilJetPtEvtByEventPartLevel[itt]->Reset();
      }
   }
   //_________________________________________________________________
   // EVENT SELECTION
   fHistEvtSelection->Fill(0.5); //Count input event

   //_________________________________________________________
   //READ  TRACK AND JET CONTAINERS
   //Container operations   http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerIterateTechniques

   if(fMode == AliAnalysisTaskRevEA::kNormal || fMode == AliAnalysisTaskRevEA::kMC){
      //fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(fMyTrackContainerName.Data())); //track container detector-level   real data only
      fTrkContainerDetLevel = static_cast<AliTrackContainer*> (GetTrackContainer(0)); //track container detector-level   real data only
      fJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyJetContainerName.Data()));     //detector-level AKT jets real data or hybrid event
      fKTJetContainerDetLevel = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetContainerName.Data())); //detector-level KT jets real data or hybrid event

      fRho = GetMyRho(fKTJetContainerDetLevel); //estimated backround pt density
   }

   if( fMode == AliAnalysisTaskRevEA::kMC){  //particle level particles and jets  for  MC and embedding
      //fParticleContainerPartLevel = GetParticleContainer(fMyParticleContainerName.Data()); //pythia particle level particles
      fParticleContainerPartLevel = GetParticleContainer(1); //pythia particle level particles
      fJetContainerPartLevel      = static_cast<AliJetContainer*> (GetJetContainer(fMyJetParticleContainerName.Data()));   //pythia particle level AKT jets
      fKTJetContainerPartLevel    = static_cast<AliJetContainer*> (GetJetContainer(fMyKTJetParticleContainerName.Data()));   //pythia particle level KT jets

      fRhoMC = GetMyRho(fKTJetContainerPartLevel); //estimated backround pt density
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){
      if(IsOutlier()) return kFALSE;
   }

   //_________________________________________________________________
   // DECIDE WHETHER TO FILL SIGNAL TT OR REFERENCE TT  DEPENDING ON RANDOM  NUMBER
   fFillSigTT = kTRUE;
   if(fRandom->Integer(100) < 5) fFillSigTT = kFALSE;
   //________________________________________________________________
   //DATA ANALYSIS PARTICLE LEVEL

   if(fMode == AliAnalysisTaskRevEA::kMC){  
      FindParticleLevelTT(); 
      AnalyzeParticleLevel();
   }
   //________________________________________________________________
   //DATA ANALYSIS DETECTOR LEVEL

   //Check Reconstructed event vertex and pileup

   fIsMinBiasTrig = kFALSE; //Minimum bias event flag
   if(PassedMinBiasTrigger()){
      fIsMinBiasTrig = kTRUE;
      fHistEvtSelection->Fill(4.5); //Count Accepted input event
   }

   fIsHighMultTrig = kFALSE; //high multiplicity trigger flag
   if(PassedHighMultTrigger()){
      fIsHighMultTrig = kTRUE; //Count Accepted input event
      fHistEvtSelection->Fill(6.5); //Count Accepted input event
   }

   //_________________________________________________________________


   fTrigflag[0] = fIsMinBiasTrig;
   fTrigflag[1] = fIsHighMultTrig;

   //_________________________________________________________________

   InitEventProperties();

   FindDetectorLevelTT(); //fk Find detector level TT

   if(fMode == AliAnalysisTaskRevEA::kMC){   
      FillResponseMatrix();
      FillResponseMatrix2D();
   }

   GeneralTrackProperties();
   AnalyzeRawData();

   return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskRevEA::FindDetectorLevelTT(){
   //fk Find Detector Level TT
   TLorentzVector myTT;
   AliVParticle *track = NULL;

   for(Int_t i=0; i<fkTTbins; i++){
      fIndexTTH[i] = -1;
      fTTH[i].resize(0);
   }

   for(Int_t i=0; i < fnHadronTTBins; i++){
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
                  }
               }
            }
         }
         //chose trigger hadron TT
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(fTTH[itt].size()>0){
               fIndexTTH[itt] = fRandom->Integer((Int_t) (fTTH[itt].size()));
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskRevEA::FindParticleLevelTT(){ // Modified by KA
   //inspired by "FindDetectorLevelTT" function

   // Find Particle Level TT
   TLorentzVector myTT;
   AliVParticle *mcParticle = NULL;

   for(Int_t i = 0; i < fkTTbins; i++){
      fIndexTTH_PartLevel[i] = -1;
      fTTH_PartLevel[i].resize(0);
   }

   for(Int_t i = 0; i < fnHadronTTBins; i++){
      fHadronTT_Labels_PartLevel[i].resize(0); // Modified by KA
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){
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
                  }
               }
            }
         }

         for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
            if(fTTH_PartLevel[itt].size() > 0){
	       fIndexTTH_PartLevel[itt] = fRandom->Integer((Int_t) (fTTH_PartLevel[itt].size()));
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskRevEA::AnalyzeRawData(){
   //Analyze raw data

   AliEmcalJet  *jet = NULL;        //jet pointer real jet
   Double_t ptLJ=-1, etaLJ=999, phiLJ=0; //leading jet
   Double_t ptSJ=-1, etaSJ=999, phiSJ=0; //subleading jet
   Int_t b1,b2;
   Double_t tmparr3[3];
   Int_t idx;
   Double_t dphi     = 999.;
   Double_t deltaPhi_Angle_Abs_0Pi = 0.0;
   Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho



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

      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //chose trigger hadron TT
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         fdeltapT[itt] = 0.;

         if(fTTH[itt].size()>0){
            idx = fIndexTTH[itt];

            fdeltapT[itt] = GetDeltaPt(fTTH[itt][idx].Phi(), fTTH[itt][idx].Eta(), phiLJ, etaLJ, phiSJ, etaSJ, fRho, kDetLevel);
         }
      }

      for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
         if(!fTrigflag[itg]) continue; //check which trigger fired
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){

            fhMultTTH[itg][itt]->Fill(fTTH[itt].size());

            if(fTTH[itt].size()==0) continue; //check whether there was hadron TT

            fhRhoTTH[itg][itt]->Fill(fRho);
         }
      }



      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         if(fTTH[itt].size()==0) continue; //analyze events with hadron TT only

         for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
            if(!fTrigflag[itg]) continue; //check which trigger fired

            fhSignalTTH[itg][fkV0Mnorm][itt]->Fill(fMultV0Mnorm);
         }


         //pick up TTH hadron accoding to the index
         idx = fIndexTTH[itt];
         if(idx>-1){
//
            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!fTrigflag[itg]) continue; //check which trigger fired

               fhDeltaPtTTH_RC_V0Mnorm[itg][itt]->Fill(fMultV0Mnorm, fdeltapT[itt]);
            }

            if(fFillSigTT  && itt==0) continue;  // Do not fill reference
            if(!fFillSigTT && itt>0)  continue;  // Do not fill signal

            for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
               if(!fTrigflag[itg]) continue; //check which trigger fired

               fhTTH_V0Mnorm[itg][itt]->Fill(fMultV0Mnorm,  fTTH[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm
            }

            //recoil jets
            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jet = jetIterator.second;  // Get the pointer to jet object
               if(!jet)  continue;

               dphi = TVector2::Phi_0_2pi(jet->Phi()-fTTH[itt][idx].Phi());
               deltaPhi_Angle_Abs_0Pi = TMath::Abs(TVector2::Phi_mpi_pi(dphi)); // Modified by KA.

               jetPtCorrDet = jet->Pt() - fRho*jet->Area();

               for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                  if(!fTrigflag[itg]) continue; //check which trigger fired
                  fhRecoilJetPhiTTH_V0Mnorm[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet, deltaPhi_Angle_Abs_0Pi); // Modified by KA.
               }

               if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > TMath::Pi()/2){ //select recoil hemisphere and count jets
                  fhRecoilJetPtEvtByEvent[itt]->Fill(jet->Pt());
               }
               
               if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > fPhiCut){     //select recoil jet
                  for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
                     if(!fTrigflag[itg]) continue; //check which trigger fired
                     fhRecoilJetPtTTH_V0Mnorm[itg][itt]->Fill(fMultV0Mnorm, jetPtCorrDet);
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
      for(Int_t itg=kMB; itg<=kHM; itg++){ //@@@
         if(!fTrigflag[itg]) continue; //check which trigger fired

         fhJetPhiIncl[itg]->Fill(jet->Pt(), jet->Phi());
         fhJetEtaIncl[itg]->Fill(jet->Pt(), jet->Eta());
      }
   }

     return;
}
//_________________________________________________________________
void AliAnalysisTaskRevEA::InitEventProperties(){
   // EVENT PROPERTIES

   if((fIsMinBiasTrig || fIsHighMultTrig)){  //real data + mc det level 

      //___________________________________________
      //    INCLUSIVE EVENTS (WITHOUT TT REQUIREMENT)

      for(Int_t itg=kMB; itg<=kHM; itg++){    //@@@
         if(!fTrigflag[itg]) continue;
         //events without TT requirement
         fhRho[itg]->Fill(fRho);

         fhCentrality[itg]->Fill(fCentralityV0M, fMultV0Mnorm);
         fhSignal[itg][fkV0Mnorm]->Fill(fMultV0Mnorm);
      }
   }
}
//_________________________________________________________________
void AliAnalysisTaskRevEA::AnalyzeParticleLevel(){

   TLorentzVector myTT;
   AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
   AliVParticle *mcParticle = NULL; //mc particle
   Int_t b1,b2;
   Int_t idx;
   Double_t tmparr3[3];
   Double_t jetPtCorrPart = 0.;

   Double_t dphi             = 999.;
   Double_t deltaPhi_Angle_Abs_0Pi = 0.0; // Modified by KA


   if(fMode == AliAnalysisTaskRevEA::kMC){

      fHistEvtSelection->Fill(8.5); //Count Accepted input event

      fhRhoMBpart->Fill(fRhoMC);
      fhSignal_PartLevel[fkV0Mnorm]->Fill(fMultV0Mnorm);


      if(fParticleContainerPartLevel){
         //chose trigger hadron TT   particle level
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            if(fTTH_PartLevel[itt].size()>0){
               idx = fIndexTTH_PartLevel[itt]; // Initialize at FindParticleLevelTT, KA

               //signal in events with hadron TT   particle level

               fhSignalTTH_PartLevel[fkV0Mnorm][itt]->Fill(fMultV0Mnorm);

               //hadron trigger particle level
               if(idx>-1){

                  fhRhoTTHinMBpart[itt]->Fill(fRhoMC);

                  if(fFillSigTT && itt==0) continue;  // Do not fill reference
                  if(!fFillSigTT && itt>0) continue;  // Do not fill signal

                  fhTTH_V0Mnorm_PartLevel[itt]->Fill(fMultV0Mnorm, fTTH_PartLevel[itt][idx].Pt()); //fill trigger track pT for given V0Mnorm

                  //recoil jets  PARTICLE LEVEL
                  for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
                     // trackIterator is a std::map of AliTLorentzVector and AliVTrack
                     jetPartMC = jetIterator.second;  // Get the pointer to jet object
                     if(!jetPartMC)  continue;

                     dphi = TVector2::Phi_0_2pi(jetPartMC->Phi()-fTTH_PartLevel[itt][idx].Phi());
                     deltaPhi_Angle_Abs_0Pi = TMath::Abs(TVector2::Phi_mpi_pi(dphi)); // Modified by KA.

	             jetPtCorrPart = jetPartMC->Pt() - fRhoMC*jetPartMC->Area();

                     fhRecoilJetPhiTTH_V0Mnorm_PartLevel[itt]->Fill(fMultV0Mnorm,  jetPtCorrPart, deltaPhi_Angle_Abs_0Pi); // Modified by KA
                     fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[itt]->Fill(fMultV0Mnorm, jetPartMC->Pt(), deltaPhi_Angle_Abs_0Pi); // Modified by KA

                     if(deltaPhi_Angle_Abs_0Pi > TMath::Pi()/2){ //select recoil hemisphere and count jets
                        fhRecoilJetPtEvtByEventPartLevel[itt]->Fill(jetPartMC->Pt());
                     }

                     if(TMath::Abs(TVector2::Phi_mpi_pi(dphi)) > fPhiCut){
                        //recoil jet hadron trigger
                        fhRecoilJetPtTTH_V0Mnorm_PartLevel[itt]->Fill(fMultV0Mnorm,  jetPtCorrPart);
                        fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[itt]->Fill(fMultV0Mnorm, jetPartMC->Pt());
                     }
                  }
               }

               //count number of jets with pT larger than something in recoil region
               tmparr3[2] = fMultV0Mnorm;
               for(Int_t ii = 1; ii<=fhNumberOfHighPtJetsRecoilPartLevel[itt]->GetAxis(0)->GetNbins(); ii++){
	          tmparr3[0] = fhNumberOfHighPtJetsRecoilPartLevel[itt]->GetAxis(0)->GetBinLowEdge(ii);
	          b1 = fhRecoilJetPtEvtByEventPartLevel[itt]->GetXaxis()->FindBin(tmparr3[0] + 1e-5);
	          b2 = fhRecoilJetPtEvtByEventPartLevel[itt]->GetXaxis()->GetNbins()+1;  //include overflow bin
	          tmparr3[1] = fhRecoilJetPtEvtByEventPartLevel[itt]->Integral(b1,b2);
	          fhNumberOfHighPtJetsRecoilPartLevel[itt]->Fill(tmparr3);
               }
            }
         }
      }
    
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //Response matrix normalization - spectrum of all generator level jets in acceptance
      if(fJetContainerPartLevel){
         for(auto jetPartIterator : fJetContainerPartLevel->accepted_momentum() ){
            jetPartMC = jetPartIterator.second;  // Get the pointer to mc particle object
            if(!jetPartMC)  continue;
    
            jetPtCorrPart = jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
    
            fhJetPtPartLevelCorr->Fill(jetPtCorrPart);
            fhJetPtPartLevelZero->Fill(jetPartMC->Pt());
         }
      }
   }
}
//_________________________________________________________________
void AliAnalysisTaskRevEA::FillResponseMatrix(){

   if(fIsMinBiasTrig){

      //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX
      AliEmcalJet  *jetPartMC = NULL;  //jet pointer particle level MC jet
      AliEmcalJet  *jetDetMC  = NULL;  //jet pointed detector level MC jet
      AliVParticle *track = NULL; //jet constituent
      AliVParticle *mcParticle = NULL; //mc particle
      Bool_t bRecPrim = kFALSE;
      Double_t jetPtCorrDet  = 0.;  //detector level jet pt corrected for rho
      Double_t jetPtCorrPart = 0.;
      Double_t array[5];
      Double_t xphi;

      //__________________________________________________________
      //  FILL JET RESPONSE MATRIX
      //__________________________________________________________


      //1) Find closest particle level and detector level jets
      if(fJetContainerDetLevel){

         for(auto jetIterator : fJetContainerDetLevel->accepted_momentum() ){
            jetDetMC = jetIterator.second;  // Get the pointer to jet object
            if(!jetDetMC)  continue;

            //Fill Response matrix
            jetPartMC =  jetDetMC->ClosestJet();
            if(!jetPartMC) continue;

            // KA: Set new relax condition "5e-4" because we consider particle level jets with pT > 1 MeV
            if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

            jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
            jetPtCorrDet  =  jetDetMC->Pt()  - jetDetMC->Area()*fRho;

            fhJetPtPartLevelVsJetPtDetLevelCorr->Fill(jetPtCorrDet,jetPtCorrPart); //response matrix
            fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr->Fill(jetPtCorrDet, jetPartMC->Pt()); //response matrix (added by KA)
            fhJetPtZeroPartLevelVsJetPtZeroDetLevel->Fill(jetDetMC->Pt(),jetPartMC->Pt()); //response matrix

            fhJetPtResolutionVsPtPartLevel->Fill(jetPartMC->Pt(),(jetDetMC->Pt()-jetPartMC->Pt())/jetPartMC->Pt()); //jet pT resolution
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

         if(fTTH[itt].size() > 0){
            idx = fIndexTTH[itt]; //fk initialized in FindDetectorLevelTT

            //recoil jets  PARTICLE LEVEL
            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
               // trackIterator is a std::map of AliTLorentzVector and AliVTrack
               jetPartMC = jetIterator.second;  // Get the pointer to jet object
               if(!jetPartMC)  continue;
               jetDetMC =  jetPartMC->ClosestJet();
               xphi = (jetDetMC) ?  jetDetMC->Phi() : jetPartMC->Phi();  //for those cases where there is detector level phi use detectro level phi, else use particle level phi

               if(TMath::Abs(TVector2::Phi_mpi_pi(xphi - fTTH[itt][idx].Phi())) > fPhiCut){  //fk
                  jetPtCorrPart = jetPartMC->Pt() - fRhoMC*jetPartMC->Area();
                  fhRecoilJetPtPartLevelCorr[itt]->Fill(jetPtCorrPart);
                  fhRecoilJetPtZeroPartLevel[itt]->Fill(jetPartMC->Pt());
               }
            }

            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum()){
               jetDetMC = jetIterator.second;  // Get the pointer to jet object
               if(!jetDetMC)  continue;

               jetPartMC =  jetDetMC->ClosestJet();
               if(!jetPartMC){
                  continue;
               }
               if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

               jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
               jetPtCorrDet  =  jetDetMC->Pt() -  jetDetMC->Area()*fRho;

               if(TMath::Abs(TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[itt][idx].Phi())) > fPhiCut){   //fk
                  fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt]->Fill(jetPtCorrDet,jetPtCorrPart);
                  fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt]->Fill(jetPtCorrDet,jetPartMC->Pt());
                  fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt]->Fill(jetDetMC->Pt(),jetPartMC->Pt());
               }
            }
         }
      }

      // KA: fill response matrix with jets on particle and detector levels recoiling from corresponding TT
      //Important: TT labels must be the same!!!
      Int_t Label_Det_TT, Label_Part_TT; // Modified by KA

      //count number of jets with pT larger than 15 GeV and 20 GeV in particle and detector level events
      Double_t nJetsPt15_PartLevel = 1e-5; 
      Double_t nJetsPt15_DetLevel  = 1e-5; 
      Double_t nJetsPt20_PartLevel = 1e-5; 
      Double_t nJetsPt20_DetLevel  = 1e-5; 


      for(Int_t iTT = 0; iTT < fnHadronTTBins; iTT++){
         if(fFillSigTT  && iTT == 0) continue;  //fk
         if(!fFillSigTT && iTT > 0)  continue;  //fk

         if(fTTH[iTT].size() > 0 && fTTH_PartLevel[iTT].size() > 0){ // presence of TT on both levels

            Label_Part_TT = -1;
            idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT

            Label_Det_TT = TMath::Abs(fHadronTT_Labels[iTT][idx]); // Label of selected TT on detector level

            for(Int_t k = 0; k < (Int_t)(fHadronTT_Labels_PartLevel[iTT].size()); k++){
               if(Label_Det_TT == fHadronTT_Labels_PartLevel[iTT][k]){  // Look for TT on particle level with same label
                  idx_PartLevel = k;
                  Label_Part_TT = fHadronTT_Labels_PartLevel[iTT][k];
                  break;
               }
            }

            if(Label_Part_TT<0) continue;

            fhTT_Corresp[iTT]->Fill(fTTH[iTT][idx].Pt(), fTTH_PartLevel[iTT][idx_PartLevel].Pt()); // Modified by KA

            //recoil jets PARTICLE LEVEL
            for(auto jetIterator : fJetContainerPartLevel->accepted_momentum()){
               jetPartMC = jetIterator.second;  // Get the pointer to jet object
               if(!jetPartMC) continue;
               
                if(fFillSigTT){ //just for TT{20,30}
                   if(TMath::Abs(TVector2::Phi_mpi_pi(jetPartMC->Phi() - fTTH_PartLevel[iTT][idx_PartLevel].Phi()) > TMath::Pi()/2)){
                      if(jetPartMC->Pt() > 15.0)  nJetsPt15_PartLevel++;
                      if(jetPartMC->Pt() > 20.0)  nJetsPt20_PartLevel++;
                   }
                }

               if(TMath::Abs(TVector2::Phi_mpi_pi(jetPartMC->Phi() - fTTH_PartLevel[iTT][idx_PartLevel].Phi())) > fPhiCut){  // KA
                  jetPtCorrPart = jetPartMC->Pt() - fRhoMC*jetPartMC->Area();
                  fhRecoilJetPtPartLevel_CorrespTT[iTT]->Fill(jetPtCorrPart); // Modified by KA
                  fhRecoilJetPtZeroPartLevel_CorrespTT[iTT]->Fill(jetPartMC->Pt()); // Modified by KA
               }
            }

            for(auto jetIterator : fJetContainerDetLevel->accepted_momentum()){
               jetDetMC = jetIterator.second;  // Get the pointer to jet object
               if(!jetDetMC) continue;

               jetPartMC = jetDetMC->ClosestJet();
               if(!jetPartMC) continue; // IMPORTANT TO ADD
               if(jetPartMC->Pt() < 5e-4) continue; //prevents matching with a ghost

               jetPtCorrPart =  jetPartMC->Pt() - jetPartMC->Area()*fRhoMC;
               jetPtCorrDet  =  jetDetMC->Pt()  - jetDetMC->Area()*fRho;

               if(fFillSigTT){ //just for TT{20,30}
                  if(TMath::Abs(TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[iTT][idx].Phi())) > TMath::Pi()/2){ //select recoil hemisphere and count jets
                     if(jetDetMC->Pt() > 15.0)  nJetsPt15_DetLevel++;
                     if(jetDetMC->Pt() > 20.0)  nJetsPt20_DetLevel++;
                  }
               }

               if(TMath::Abs(TVector2::Phi_mpi_pi(jetDetMC->Phi()-fTTH[iTT][idx].Phi())) > fPhiCut && TMath::Abs(TVector2::Phi_mpi_pi(jetPartMC->Phi()-fTTH_PartLevel[iTT][idx_PartLevel].Phi())) > fPhiCut){   // KA: look for events when both matched jets are in the recoil
                  fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[iTT]->Fill(jetPtCorrDet, jetPtCorrPart); // Modified by KA
                  fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[iTT]->Fill(jetPtCorrDet, jetPartMC->Pt()); // Modified by KA
                  fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[iTT]->Fill(jetDetMC->Pt(), jetPartMC->Pt()); // Modified by KA
               }
            }

            //Fill Njet response matrices
            if(fFillSigTT){ //just for TT{20,30}
               fhNjetReMx_V0MnormDetLev_15GeV->Fill(nJetsPt15_DetLevel, nJetsPt15_PartLevel, fMultV0Mnorm);
               fhNjetNorm_V0MnormDetLev_15GeV->Fill(nJetsPt15_PartLevel, fMultV0Mnorm);

               fhNjetReMx_V0MnormDetLev_20GeV->Fill(nJetsPt20_DetLevel, nJetsPt20_PartLevel, fMultV0Mnorm);
               fhNjetNorm_V0MnormDetLev_20GeV->Fill(nJetsPt20_PartLevel, fMultV0Mnorm);
            }

         }
      }
   }
   return;
}
//________________________________________________________________________
void AliAnalysisTaskRevEA::FillResponseMatrix2D(){
   //Fill 4D response matrix

   if(fIsMinBiasTrig && fMode == AliAnalysisTaskRevEA::kMC){
      //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX
      AliEmcalJet  *jet             = NULL; //jet pointer real jet
      AliEmcalJet  *jetDetMC        = NULL; //jet pointed detector level MC jet
      AliVParticle *track_PartLevel = NULL; //jet constituent
      AliVParticle *track_DetLevel  = NULL; //mc particle

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

               } else if(jetDetMC->Pt() < 1e-10){ //associated MC detector level jet is a ghost jet ->  Fill input for the miss function

                  fhPhi_JetPtPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jetPtCorrPart);
                  fhPhi_JetPtZeroPartLevel_InclusiveJets->Fill(TMath::Abs(random_PhiAnglePartLevel), jet->Pt()); // (added by KA)

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
               }
            }
         }//end inclusive jets

         //TT events
         for(Int_t iTT = 0; iTT < fnHadronTTBins; iTT++){
            if(fTTH[iTT].size() == 0) continue;      //fk skip event, since no TT on detector level
            if(fFillSigTT  && iTT == 0) continue;  //fk
            if(!fFillSigTT && iTT > 0)  continue;  //fk

            idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT
            phi_Angle_Of_TT = fTTH[iTT][idx].Phi(); //fk

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

            if(fTTH[iTT].size() > 0 && fTTH_PartLevel[iTT].size() > 0){ // presence of TT on both levels

               Label_Part_TT = 0;
               idx = fIndexTTH[iTT]; //fk initialized in FindDetectorLevelTT
               // idx_PartLevel = fIndexTTH_PartLevel[iTT]; //KA initialized in FindParticleLevelTT

               Label_Det_TT = TMath::Abs(fHadronTT_Labels[iTT][idx]); // Label of selected TT on detector level

               for(Int_t k = 0; k < (Int_t) (fHadronTT_Labels_PartLevel[iTT].size()); k++){
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
void AliAnalysisTaskRevEA::GeneralTrackProperties(){

   //basic properties of tracks 
   AliVParticle *track = NULL; //jet constituent

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
         }
      }
   }

}
//________________________________________________________________________
void AliAnalysisTaskRevEA::Terminate(Option_t *){
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
AliAnalysisTaskRevEA::~AliAnalysisTaskRevEA(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fRandom;
   delete fHelperClass;
}
//________________________________________________________________________
void AliAnalysisTaskRevEA::UserCreateOutputObjects(){
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
   fHistEvtSelection = new TH1D("fHistEvtSelection", "event selection", 10, 0, 10);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"events IN"); //0-1
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"incomplete DAQ (rejected)"); //1-2
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)"); //2-3
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)"); //3-4
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"MB"); //4-5
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,""); //5-6
   fHistEvtSelection->GetXaxis()->SetBinLabel(7,"High Mult"); //6-7
   fHistEvtSelection->GetXaxis()->SetBinLabel(8,"IsEventSelected"); //7-8
   fHistEvtSelection->GetXaxis()->SetBinLabel(9,"MC MB events"); //8-9
   fHistEvtSelection->GetXaxis()->SetBinLabel(10,"IsOutlier"); //9-10


   fOutput->Add(fHistEvtSelection);


   //Trigger track pT spectrum single inclusive for  MB  versus V0M
   Int_t    nbinsV0M     = 100;
   Double_t maxV0M       = 1000.;
   Int_t    nbinsV0Mnorm = 200;
   Double_t maxV0Mnorm   = 20.;

   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   fhVertexZall =  new TH1D("fhVertexZall","z vertex without cut",40,-20,20);
   fOutput->Add(fhVertexZall);

   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);
   fOutput->Add(fhVertexZ);

   //-------------------------
   TString trig[]={"MB","HM"};

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      name   = Form("fhTrackEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive track vs pT %s",trig[itg].Data());
      fhTrackEtaIncl[itg] = new TH2D(name.Data(), object.Data(), 200, 0, 200, 60, -0.9, 0.9); // Modified by KA (old value: 50,0,100)
      fOutput->Add((TH2D*) fhTrackEtaIncl[itg]);

      name   = Form("fhTrackPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist tracks vs pT %s",trig[itg].Data());
      fhTrackPhiIncl[itg] = new TH2D(name.Data(), object.Data(), 200, 0, 200, 50, 0, 2*TMath::Pi()); // Modified by KA (old value: 50,0,100)
      fOutput->Add((TH2D*) fhTrackPhiIncl[itg]);
   }


   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      name   = Form("fhJetEtaIncl%s",trig[itg].Data());
      object = Form("Eta dist inclusive jets vs pTjet %s",trig[itg].Data());
      fhJetEtaIncl[itg] = new TH2D(name.Data(),object.Data(), 100, 0, 100, 60, -0.9, 0.9); // Modified by KA (old value: 40,-0.9,0.9)
      fOutput->Add((TH2D*) fhJetEtaIncl[itg]);

      name   = Form("fhJetPhiIncl%s",trig[itg].Data());
      object = Form("Azim dist jets vs pTjet %s",trig[itg].Data());
      fhJetPhiIncl[itg] = new TH2D(name.Data(),object.Data(), 100, 0, 100, 50, 0, 2*TMath::Pi());
      fOutput->Add((TH2D*) fhJetPhiIncl[itg]);
   }

   //RHO
   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      name   = Form("hRhokt_%s",trig[itg].Data());
      object = Form("Rho kt det level %s",trig[itg].Data());

      fhRho[itg] = new TH1D( name.Data(), object.Data(),1000,0,100);
      fOutput->Add((TH1D*) fhRho[itg]);
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){ //HADRON TT
         name = Form("hRhokt_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRhoTTH[itg][itt] = (TH1D*)  fhRho[itg]->Clone(name.Data());      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTH[itg][itt]);
      }
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){
      name = Form("hRhokt_MB_part");
      fhRhoMBpart = new TH1D( name.Data(), name.Data(),1000,0,100);
      fhRhoMBpart->SetTitle(Form("Rho kt  min bias part level"));
      fOutput->Add((TH1D*) fhRhoMBpart);

      for(Int_t itt=0; itt<fnHadronTTBins;itt++){
         name = Form("hRhokt_MB_TTH%d_%d_part", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhRhoTTHinMBpart[itt] = (TH1D*) fhRhoMBpart->Clone(name.Data());                      //! in events MB with hadron TT
         fOutput->Add((TH1D*) fhRhoTTHinMBpart[itt]);
      }
   }


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
   for(Int_t itg=kMB; itg<=kHM;itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      name = Form("hCentrality_%s_V0Mnorm", trig[itg].Data());
      fhCentrality[itg] = new TH2D(name.Data(), name.Data(), narrcent, arrcent, 400,0,20);
      fOutput->Add((TH2D*) fhCentrality[itg]);
   }


   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //SIGNAL

   TString signal[]={"V0Mnorm"};
   Float_t signalL[]={0};
   Float_t signalH[]={15};
   Int_t   signalN[]={150};

   for(Int_t itg=kMB; itg<=kHM;itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t ic=0; ic<fkCE;ic++){ //MB
         name = Form("hSignal_%s_%s", trig[itg].Data(), signal[ic].Data());
         fhSignal[itg][ic] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
         fOutput->Add((TH1D*) fhSignal[itg][ic]);
      }
   }

   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t ic=0; ic<fkCE; ic++){ //TT hadron
         for(Int_t itt=0; itt<fnHadronTTBins; itt++){
            name = Form("hSignal_%s_%s_TTH%d_%d", trig[itg].Data(), signal[ic].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
            fhSignalTTH[itg][ic][itt] = new TH1D(name.Data(), name.Data(), signalN[ic], signalL[ic], signalH[ic]);
            fOutput->Add((TH1D*) fhSignalTTH[itg][ic][itt]);
         }
      }
   }


   if(fMode == AliAnalysisTaskRevEA::kMC){ //PARTICLE LEVEL SIGNAL DISTRIBUTIONS

      TString signalmc[]={"V0Mnorm"};
      Float_t signalLmc[]={0};
      Float_t signalHmc[]={20};
      Int_t signalNmc[]={200};

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
   }


   //+++++++++++++++++++++++++++++++


   //Trigger track candidate multiplicity
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hMultTT_%s_TTH%d_%d", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhMultTTH[itg][itt] = new TH1D(name.Data(),name.Data(),100,0,100);
         fOutput->Add((TH1D*)  fhMultTTH[itg][itt]);
      }
   }



   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;


      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_%s_TTH%d_%d_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm[itg][itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm[itg][itt]);
      }
   }


   if(fMode == AliAnalysisTaskRevEA::kMC){
      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("hTT_MB_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTTH_V0Mnorm_PartLevel[itt] = new TH2D(name.Data(),name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, 100, 0, 100);
         fOutput->Add((TH2D*) fhTTH_V0Mnorm_PartLevel[itt]);
      }

      //Trigger track pT spectrum single inclusive for MB  versus  V0Mnorm without V0 coincidence
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){

         name = Form("NumberOf_Corresp_TT%d_%d", fHadronTTLowPt[itt],fHadronTTHighPt[itt]);
         fhTT_Corresp[itt] = new TH2D(name.Data(),name.Data(), 100, 0, 100, 100, 0, 100);
         fOutput->Add((TH2D*) fhTT_Corresp[itt]);
      }
   }



   // ++++++++++++++++++++++++++++++++++ Binning ++++++++++++++++++++++++++++++++++++
   //Jet pT corrected on RhokT 
   const Int_t numberBins_jetpT_RhokT = 270;
   Double_t jetPtRhokT_Bins[numberBins_jetpT_RhokT + 1];
   for(Int_t i = 0; i <= numberBins_jetpT_RhokT; i++) jetPtRhokT_Bins[i] = i - 20;  //(-20,250)

   //Jet pT not corrected on RhokT 
   const Int_t numberBins_jetpT_Zero = 250;
   Double_t jetPtZero_Bins[numberBins_jetpT_Zero + 1];
   for(Int_t i = 0; i <= numberBins_jetpT_Zero; i++) jetPtZero_Bins[i] = i;  //(0,250)

   //Phi angle
   const Int_t ndeltaPhiBins = 40; 
   Double_t deltaPhiBins[ndeltaPhiBins+1];
   Double_t p = TMath::Pi()/ndeltaPhiBins;
   for(Int_t i = 0; i <= ndeltaPhiBins; i++) deltaPhiBins[i] = i*p;  

   //Phi angle: Inclusive jets
   Double_t deltaPhiBins_InclusiveJets[2*ndeltaPhiBins+1]; 
   Double_t step = 0.5*TMath::TwoPi()/ndeltaPhiBins;
   for(Int_t i = 0; i <= 2*ndeltaPhiBins; i++) deltaPhiBins_InclusiveJets[i] = i*step;

   //V0M normalized
   Double_t arrV0Mnorm[nbinsV0Mnorm+1];
   p = maxV0Mnorm/nbinsV0Mnorm;
   for(Int_t i = 0; i<=nbinsV0Mnorm; i++){
      arrV0Mnorm[i] = i*p; //(0,20.0)
   }


   //dphi of recoil jets associated to semi-inclusive hadron TT in MB with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         name = Form("RecoilJetPhi_%s_TTH%d_%d_Rhokt_V0Mnorm", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPhiTTH_V0Mnorm[itg][itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); 
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm[itg][itt]);
      }
   }


   //RECOIL JET SPECTRA
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhRecoilJetPt_%s_TTH%d_%d_V0Mnorm_Rhokt", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtTTH_V0Mnorm[itg][itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm[itg][itt]);
      }
   }

   //TTH recoil jet distributions for MC
   if(fMode == AliAnalysisTaskRevEA::kMC){

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){ //! recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         // Jet pT is CORRECTED on RhokT
         name = Form("RecoilJetPt_MB_TTH%d_%d_V0Mnorm_Rhokt_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtTTH_V0Mnorm_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtTTH_V0Mnorm_PartLevel[itt]);

         // Jet pT is NOT CORRECTED on RhokT
         name = Form("RecoilJetPtZero_MB_TTH%d_%d_V0Mnorm_Rhokt_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0., maxV0Mnorm, numberBins_jetpT_Zero, jetPtZero_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[itt]);


      }

      // dphi of recoil jets associated to semi-inclusive hadron TT in MB with V0Mnorm (fMultV0Mnorm, jetPtCorrDet, dphi);
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){
         // Jet pT is CORRECTED on RhokT
         name = Form("RecoilJetPhi_MB_TTH%d_%d_V0Mnorm_Rhokt_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPhiTTH_V0Mnorm_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins, ndeltaPhiBins, deltaPhiBins); // added by KA
         fOutput->Add((TH3D*) fhRecoilJetPhiTTH_V0Mnorm_PartLevel[itt]);

         // Jet pT is NOT CORRECTED on RhokT
         name = Form("RecoilJetPtZero_DeltaPhi_TTH%d_%d_V0Mnorm_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[itt] = new TH3D(name.Data(), name.Data(), nbinsV0Mnorm, arrV0Mnorm, numberBins_jetpT_Zero, jetPtZero_Bins, ndeltaPhiBins, deltaPhiBins); // added by KA
         fOutput->Add((TH3D*) fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[itt]);
      }
   }


   //+++++++++++++++++++++++++++ RECOIL JETS WITH TTC ++++++++++++++++++++++++++++++++++++++++++

   //delta pT distributions versus V0Mnorm   = V0M/mean V0M
   for(Int_t itg=kMB; itg<=kHM; itg++){  //TTH
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){        //!  recoil jets associated to semi-inclusive hadron TT  in MB  with V0Mnorm
         name = Form("fhDeltaPtTTH_%s_RC_V0Mnorm_TTH%d_%d_Rhokt", trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhDeltaPtTTH_RC_V0Mnorm[itg][itt] = new TH2D(name.Data(), name.Data(), nbinsV0Mnorm, 0, maxV0Mnorm, numberBins_jetpT_RhokT, jetPtRhokT_Bins); // changed range: old-> (200, -20, 180) (added by KA)
         fOutput->Add((TH2D*) fhDeltaPtTTH_RC_V0Mnorm[itg][itt]);
      }
   }



   if(fMode == AliAnalysisTaskRevEA::kMC){
      name = Form("JetPtPartLevel_Rhokt");
      fhJetPtPartLevelCorr = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
      fOutput->Add((TH1D*) fhJetPtPartLevelCorr);

      fhJetPtPartLevelZero = new TH1D("JetPtZeroPartLevel","JetPtZeroPartLevel", numberBins_jetpT_Zero, jetPtZero_Bins); 
      fOutput->Add((TH1D*) fhJetPtPartLevelZero);

      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //Normalization of response matrix filled from recoil jets
         name = Form("RecoilJetPtPartLevel_Rhokt_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtPartLevelCorr[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtPartLevelCorr[itt]);

         //Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevel_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
         fhRecoilJetPtZeroPartLevel[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtZeroPartLevel[itt]);

         //______________________________________________________
         //Normalization of response matrix filled from recoil jets on PARTICLE level
         name = Form("RecoilJetPt_Rhokt_PartLevel_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtPartLevel_CorrespTT[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtPartLevel_CorrespTT[itt]);

         //Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevel_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtZeroPartLevel_CorrespTT[itt] = new TH1D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins);
         fOutput->Add((TH1D*) fhRecoilJetPtZeroPartLevel_CorrespTT[itt]);
      }


      //1D unfolding
      name = Form("JetPtPartLevelVsJetPtDetLevel_Rhokt");
      fhJetPtPartLevelVsJetPtDetLevelCorr = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
      fOutput->Add((TH2D*) fhJetPtPartLevelVsJetPtDetLevelCorr);

      fhJetPtZeroPartLevelVsJetPtZeroDetLevel = new TH2D("JetPtZeroPartLevel_Vs_JetPtZeroDetLevel","JetPtZeroPartLevel_Vs_JetPtZeroDetLevel", numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
      fOutput->Add((TH2D*) fhJetPtZeroPartLevelVsJetPtZeroDetLevel);

      name = Form("JetPtZeroPartLevel_Vs_JetPtDetLevel_Rhokt");
      fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
      fOutput->Add((TH2D*) fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr);


      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //FF Response matrix filled from recoil jets
	      name = Form("RecoilJetPtPartLevelVsJetPtDetLevel_Rhokt_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[itt]);

         // Part level jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtDetLevel_Rhokt_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[itt]);

         //  Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[itt]);

         //____________________________________________________________
         //KA Reponse matrix filled from recoil jets wrt corresponding TT
	      name = Form("RecoilJetPtPartLevelVsJetPtDetLevel_Rhokt_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_RhokT, jetPtRhokT_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[itt]);

         // Part level jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtDetLevel_Rhokt_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_RhokT, jetPtRhokT_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[itt]);

         //  Jet pT not corrected on RhokT
         name = Form("RecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_Corresp_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]); 
         fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt] = new TH2D(name.Data(), name.Data(), numberBins_jetpT_Zero, jetPtZero_Bins, numberBins_jetpT_Zero, jetPtZero_Bins); 
         fOutput->Add((TH2D*) fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[itt]);
      }

      name = Form("fhJetPtResolutionVsPtPartLevel_Rhokt");
      fhJetPtResolutionVsPtPartLevel = new TH2D(name.Data(), name.Data(),100,0,100,50,0,2);
      fOutput->Add((TH2D*) fhJetPtResolutionVsPtPartLevel);
   }


   //2D unfolding -------------------------------
   if(fMode == AliAnalysisTaskRevEA::kMC){

      //Auxiliary variables for Sparse object

      // Particle Level jet pT is CORRECTED on Rhokt
      const Int_t fNumberOfDimensions = 4;
      const Int_t fSparseBinsNumber_JetPtRhokT[fNumberOfDimensions] = {ndeltaPhiBins, numberBins_jetpT_RhokT, ndeltaPhiBins, numberBins_jetpT_RhokT};

      //Inclusive jets
      // Missed events
      name = Form("MissedEvents_Phi_JetPt_Rhokt_PartLevel_InclusiveJets_%s", trig[kMB].Data());
      fhPhi_JetPtPartLevel_InclusiveJets = new TH2D (name.Data(), "Missed events phi vs inclusive jet pT RhokT part level",  ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
      fOutput->Add((TH2D*) fhPhi_JetPtPartLevel_InclusiveJets);

      //Sparse object as basis for RooUnfoldResponse object
      name = Form("Phi_JetPtRhokT_DetLevel_Vs_Phi_JetPt_Rhokt_PartLevel_InclusiveJets_%s", trig[kMB].Data());
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets = new THnSparseD (name.Data(), "Phi vs Jet pT RhokT for filling response matrix with inclusive jets", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
      fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
      fOutput->Add((THnSparse*) fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets);

      //Jets from events with TT
      for(Int_t itt = 0; itt < fnHadronTTBins; itt++){
         //Missed events for RM initialization, K.A.
         name = Form("MissedEvents_DeltaPhi_JetPt_Rhokt_PartLevel_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtPartLevel[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT RhokT part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtPartLevel[itt]);

         // Modified by KA
         name = Form("MissedEvents_DeltaPhi_JetPt_Rhokt_PartLevel_Corresp_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data()); // first-> event trigger; second-> TT bin
         fhDeltaPhi_JetPtPartLevel_CorrespTT[itt] = new TH2D (name.Data(), "Missed events delta phi vs jet pT RhokT part level", ndeltaPhiBins, deltaPhiBins, numberBins_jetpT_RhokT, jetPtRhokT_Bins);
         fOutput->Add((TH2D*) fhDeltaPhi_JetPtPartLevel_CorrespTT[itt]);

         //Sparse object as basis for RooUnfoldResponse object
         name = Form("DeltaPhi_JetPt_Rhokt_DetLevel_Vs_DeltaPhi_JetPt_Rhokt_PartLevel_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT RhokT for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[itt]);

         // Modified by KA
         name = Form("DeltaPhi_JetPt_Rhokt_DetLevel_Vs_DeltaPhi_JetPt_Rhokt_PartLevel_Corresp_TTH%d_%d_%s", fHadronTTLowPt[itt], fHadronTTHighPt[itt], trig[kMB].Data());
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt] = new THnSparseD (name.Data(), "Delta phi vs jet pT RhokT for filling response matrix", fNumberOfDimensions, fSparseBinsNumber_JetPtRhokT, NULL, NULL);
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(0, deltaPhiBins);    //Delta phi detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(1, jetPtRhokT_Bins); //Jet pT detector level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(2, deltaPhiBins);    //Delta phi particle level
         fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]->SetBinEdges(3, jetPtRhokT_Bins); //Jet pT particle level
         fOutput->Add((THnSparse*) fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[itt]);
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
   }
   //end of 2D unfolding -------------------------------

   //Auxiliary jet pT spectra filled event by event which will not go to output
   for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
      name = Form("fhRecoilJetPtEvtByEvent_TTH%d_%d", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);
      fhRecoilJetPtEvtByEvent[itt] = new TH1D(name.Data(),name.Data(),100,0,100);
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){
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


   //Count multiplicity of recoil high-pT jets per event
   for(Int_t itg=kMB; itg<=kHM; itg++){
      if((fMode == AliAnalysisTaskRevEA::kMC) && itg == kHM) continue;

      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
         name = Form("fhNumberOfHighPtJetsRecoil_%s_TTH%d_%d",trig[itg].Data(), fHadronTTLowPt[itt], fHadronTTHighPt[itt]);

         fhNumberOfHighPtJetsRecoil[itg][itt] = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoil[itg][itt]);
      }
   }


   if(fMode == AliAnalysisTaskRevEA::kMC){
      for(Int_t itt=0; itt<fnHadronTTBins; itt++){    //response matrix in events with TTH
         name = Form("fhNumberOfHighPtJetsRecoil_MB_TTH%d_%d_PartLevel", fHadronTTLowPt[itt], fHadronTTHighPt[itt]);

         fhNumberOfHighPtJetsRecoilPartLevel[itt] = new  THnSparseF(name.Data(),"Number of recoil jets with pt above X particle level", khighptjetdim, highptjetbins, highptjetxmin, highptjetxmax);
         fOutput->Add((THnSparse*) fhNumberOfHighPtJetsRecoilPartLevel[itt]);
      }
   }

   if(fMode == AliAnalysisTaskRevEA::kMC){ 
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
Bool_t AliAnalysisTaskRevEA::RetrieveEventObjects() {
   //
   // retrieve event objects
   //
    if(!AliAnalysisTaskEmcalJet::RetrieveEventObjects())  return kFALSE;

   return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskRevEA::Run(){
   // Run analysis code here, if needed. It will be executed before FillHistograms().


   return kTRUE;
}

//________________________________________________________________________

Double_t AliAnalysisTaskRevEA::GetDeltaPt(Double_t phiTT, Double_t etaTT, Double_t phiLJ, Double_t etaLJ, Double_t phiSJ, Double_t etaSJ, Double_t rho, Int_t level){

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
