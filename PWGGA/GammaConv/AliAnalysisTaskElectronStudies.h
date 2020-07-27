#include "AliAODConversionPhoton.h"
#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliCaloPhotonCuts.h"
#include "AliCalorimeterUtils.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"
#include "AliAODCaloCluster.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliAODMCParticle.h"
#include "TVector3.h"
#include "AliTrackerBase.h"
#include "TLorentzVector.h"
#include "AliRhoParameter.h"
#include "AliCaloTrackMatcher.h"

#ifndef AliAnalysisTaskElectronStudies_cxx
#define AliAnalysisTaskElectronStudies_cxx

class AliAnalysisTaskElectronStudies : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskElectronStudies();
    AliAnalysisTaskElectronStudies(const char *name);
    virtual ~AliAnalysisTaskElectronStudies();
    AliCalorimeterUtils * GetCaloUtils()     { if (!fCaloUtils) fCaloUtils = new AliCalorimeterUtils() ; 
                                             return fCaloUtils     ; }
    virtual void   UserCreateOutputObjects  ();
    virtual Bool_t Notify                   ();
    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    virtual void   UserExec                 ( Option_t *option );
    virtual void   Terminate                ( Option_t * );
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
    void SetIsMC(Int_t setting) {fIsMC = setting;}
    void SetEventCuts                       ( AliConvEventCuts* conversionCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fEventCuts=conversionCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetClusterCutsEMC                  ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fClusterCutsEMC=clusterCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }                                                                                                                                             
    void SetConvCuts                  ( AliConversionPhotonCuts* convCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fConvCuts=convCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                              }          
    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}
    void SetEtaMatching(Double_t p0,Double_t p1=0.,Double_t p2 =0){
        fMatchingParamsEta[0] = p0;
        fMatchingParamsEta[1] = p1;
        fMatchingParamsEta[2] = p2;
    }
    void SetPhiMatching(Double_t p0,Double_t p1=0.,Double_t p2 =0){
        fMatchingParamsPhi[0] = p0;
        fMatchingParamsPhi[1] = p1;
        fMatchingParamsPhi[2] = p2;
    }
  protected:
    AliVEvent*                  fInputEvent;                //!<!
    AliMCEvent*                 fMCEvent;                   //!<!
    Double_t                    fWeightJetJetMC;            //
    TList*                      fOutputList;                //!<!
    TTree*                      fAnalysisTree;              //!<!
    Int_t                       fIsMC;                      //
    Bool_t                      fIsHeavyIon;                //
    AliV0ReaderV1*              fV0Reader;        //!<! V0Reader for basic conversion photon selection
    TString                     fV0ReaderName;    
    AliPIDResponse*             fPIDResponse;       //!<! V0Reader for basic conversion photon selection       
    TClonesArray*               fReaderGammas;     //!<! array with photon from fV0Reader                      //
    TClonesArray*               fAODMCTrackArray;    // storage of track array
   
    AliEMCALGeometry*           fGeomEMCAL;    // pointer to EMCAL geometry
    
    // cuts and setting
    TString                     fCorrTaskSetting;           //
    AliConvEventCuts*           fEventCuts;                 // event cuts
    AliCaloPhotonCuts*          fClusterCutsEMC;            // emc cluster cuts used for signal clusters (clusters that are stored to tree)
    AliConversionPhotonCuts*    fConvCuts;                  // Cuts used by the V0Reader
    AliCalorimeterUtils*        fCaloUtils;
    
    // Track cuts
    Int_t                       fMinClsTPC;  // 
    Double_t                    fChi2PerClsTPC;  // 
    Int_t                       fMinClsITS;  // 
    Double_t                    fEtaCut;  // 
    Double_t                    fPtCut;  // 
    Double_t                    fYMCCut;  // 

    Double_t                    fMatchingParamsPhi[3];// [0] + (pt + [1])^[2]
    Double_t                    fMatchingParamsEta[3];//

    // histos
    TH1F*                       fHistoNEvents;   // 
    TH1F*                       fHistoNEventsWOWeight;   // 

    Long64_t                    fTreeBuffSize;           ///< allowed uncompressed buffer size per tree
    Long64_t                    fMemCountAOD;            //!<! accumulated tree size before AutoSave

    Int_t                       fTrackMatcherRunningMode; // CaloTrackMatcher running mode

    // tree
    Float_t fBuffer_ClusterE; 
    Float_t fBuffer_ClusterEta; 
    Float_t fBuffer_ClusterPhi; 
    Float_t fBuffer_ClusterM02; 
    Float_t fBuffer_Track_Pt; 
    Float_t fBuffer_Track_P; 
    Float_t fBuffer_Track_Eta; 
    Float_t fBuffer_Track_Phi; 
    Float_t fBuffer_Track_NSigmaElec; 
    Int_t  fBuffer_Track_IsFromV0; 
    Float_t fBuffer_MC_True_Cluster_E; 
    Float_t fBuffer_MC_True_Track_E; 
    Float_t fBuffer_MC_True_Track_Pt; 
    Float_t fBuffer_MC_True_Track_P; 
    Int_t fBuffer_MC_Track_Is_Electron; 
  

  private:
    ULong64_t GetUniqueEventID      ( AliVHeader *header);
    void CountTracks                ();
    void ResetBuffer();
    void ProcessMCCaloPhoton(AliAODCaloCluster* clus,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag);
    void ProcessCaloPhotons();
    Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
    void ProcessMatchedTrack(AliAODTrack* track, AliAODCaloCluster* clus, Bool_t isV0);
    void ProcessMCParticles();
    void ProcessTrackMatching(AliAODCaloCluster* clus);
    Bool_t IsMatchedWithConv(AliAODCaloCluster* clus, AliCaloPhotonCuts* cuts);
    Bool_t IsSameTrack(Int_t id1, Int_t id2); // check if GetID() of both tracks points to same base track
    Bool_t IsInEMCalAcceptance(AliAODConversionPhoton *photon); // check if conv photon is in EMC acc
    Bool_t IsInEMCalAcceptance(AliAODMCParticle *part); // check if mcpart is in emc acceptance
    Int_t CheckClustersForMCContribution(Int_t mclabel, TClonesArray *vclus);
    Int_t CheckConvForMCContribution(Int_t mclabel, TClonesArray *vconv);
    void RelabelAODPhotonCandidates(Bool_t mode);

    AliAnalysisTaskElectronStudies(const AliAnalysisTaskElectronStudies&); // Prevent copy-construction
    AliAnalysisTaskElectronStudies& operator=(const AliAnalysisTaskElectronStudies&); // Prevent assignment  
    ClassDef(AliAnalysisTaskElectronStudies, 1);
};

#endif

