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
#include <limits>

#ifndef AliAnalysisTaskElectronStudies_cxx
#define AliAnalysisTaskElectronStudies_cxx

typedef struct {
  UShort_t ClusterE, ClusterM02, ClusterM20,ClusterNCells,Track_E, Track_PonEMCal;
  UShort_t MC_True_Cluster_E, MC_True_Track_E;
  Short_t Track_NSigmaElec,Track_Charge,Track_dEta,Track_dPhi, Track_Px, Track_Py, Track_Pz,MC_True_Track_Px, MC_True_Track_Py, MC_True_Track_Pz;
  Bool_t Track_IsFromV0,MC_Track_Is_Electron,MC_Cluster_Is_Electron;
  UShort_t MC_ClusterTrack_Same_Electron;
  Int_t MC_True_Track_MotherPDG;
  UShort_t matchType, minR, isoE;
} treeWriteContainer;
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
   void SetTMCuts                  ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fTMCuts=clusterCuts           ;
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
    void SetMinClsTPC(Int_t mincls){
        fMinClsTPC = mincls;
    }
    void SetMinFracClsTPC(Double_t mincls){
        fMinFracClsTPC = mincls;
    }
    void SetMinClsITS(Int_t mincls){
        fMinClsITS = mincls;
    }
    void SetChi2PerClsTPC(Double_t p){
        fChi2PerClsTPC = p;
    }
    void SetEtaCut(Double_t p){
        fEtaCut = p;
    }
    void SetMinPtCut(Double_t p){
        fPtCut = p;
    }
    void SetMinNSigmaElec(Double_t nsigma){
        fMinNsigmaElec = nsigma;
    }
    void SetMaxNSigmaElec(Double_t nsigma){
        fMaxNsigmaElec = nsigma;
    }
    void SetMaxDCA(Double_t DCAxy, Double_t DCAz){
        fMaxDCAxy = DCAxy;
        fMaxDCAz = DCAz;
    }
    void SetUseRTrackMatching(Bool_t b){
        fUseRTrackMatching = b;
    }
    void SetMaxIsoRadius(Float_t r){
        fIsoMaxRadius = r;
    }
    void SetRTrackMatching(Double_t r){
        SetUseRTrackMatching(kTRUE);
        fRTrackMatching = r;
    }
    void SetTrackMatcherName(TString s){
        fTrackMatcherName = s;
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
    AliEMCALRecoUtils*          fEMCalRecoUtils;
    // cuts and setting
    TString                     fCorrTaskSetting;           //
    AliConvEventCuts*           fEventCuts;                 // event cuts
    AliCaloPhotonCuts*          fClusterCutsEMC;            // emc cluster cuts used for signal clusters (clusters that are stored to tree)
    AliCaloPhotonCuts*          fTMCuts;                 // only used for track matching
    AliConversionPhotonCuts*    fConvCuts;                  // Cuts used by the V0Reader
    AliCalorimeterUtils*        fCaloUtils;
    
    // Track cuts
    Int_t                       fMinClsTPC;  // 
    Double_t                    fMinFracClsTPC;  // 
    Double_t                    fChi2PerClsTPC;  // 
    Int_t                       fMinClsITS;  // 
    Double_t                    fEtaCut;  // 
    Double_t                    fPtCut;  // 
    Double_t                    fYMCCut;  // 
    Double_t                    fMinNsigmaElec; // only needed for signal histos
    Double_t                    fMaxNsigmaElec; // only needed for signal histos
    
    Double_t                    fMaxDCAxy; // 
    Double_t                    fMaxDCAz; // 


    Double_t                    fMatchingParamsPhi[3];// [0] + (pt + [1])^[2]
    Double_t                    fMatchingParamsEta[3];//

    Bool_t                      fUseRTrackMatching;//
    Double_t                    fRTrackMatching;//

    // histos
    TH1F*                       fHistoNEvents;   // 
    TH1F*                       fHistoNEventsWOWeight;   // 
  
    TH1F*                       fPtElectronTrack;
    TH1F*                       fPtElectronTrackInEmcalAcc;
    TH1F*                       fTruePtCluster;
    TH1F*                       fTruePtElectronCluster;
    TH1F*                       fTruePtElectronClusterMatchedWithTrack;

    TH1F*                       fTruePtElectronTrack;
    TH1F*                       fTruePtElectronTrackInEmcalAcc;

    TH1F*                       fGenPtElectrons;
    TH1F*                       fGenPtElectronsInEmcalAcc;
    
    
    TH2F*                       fTrackPvsPOnSurface;
    TH2F*                       fTrackPvsPOnSurfaceOwn;
    TH2F*                       fTrackRefPvsR;
    TH1F*                       fTrackPOnSurface;
    TH1F*                       fTrackPOnSurfaceOwn;
    TH1F*                       fTrackPOnSurfaceTrue;

    Long64_t                    fTreeBuffSize;           ///< allowed uncompressed buffer size per tree
    Long64_t                    fMemCountAOD;            //!<! accumulated tree size before AutoSave

    Int_t                       fTrackMatcherRunningMode; // CaloTrackMatcher running mode

    Float_t                     fIsoMaxRadius; //
    Float_t                     fConversionTrackMatchR; // maximum dR at which conv photon is considered matched
    // tree
    UShort_t              fBuffer_NPrimaryTracks;
    UShort_t              fBuffer_NClus;
    Bool_t fBuffer_IsProblem; // if true, some conversion ran into limits
    std::vector<UShort_t> fBuffer_ClusterE;     //!<! array buffer
    std::vector<UShort_t> fBuffer_ClusterM02; 
    std::vector<UShort_t> fBuffer_ClusterM20; 
    std::vector<UShort_t> fBuffer_ClusterNCells; 
    std::vector<UShort_t> fBuffer_Track_E; // default is always closest
    std::vector<Short_t> fBuffer_Track_Px; // default is always closest
    std::vector<Short_t> fBuffer_Track_Py; 
    std::vector<Short_t> fBuffer_Track_Pz; 
    std::vector<UShort_t> fBuffer_Track_PonEMCal; 
    std::vector<Short_t> fBuffer_Track_Charge; 
    std::vector<Short_t> fBuffer_Track_dEta; 
    std::vector<Short_t> fBuffer_Track_dPhi; 
    std::vector<Short_t> fBuffer_Track_NSigmaElec; 
    std::vector<Bool_t>  fBuffer_Track_IsFromV0; 
    std::vector<UShort_t>  fBuffer_Track_ClosestR; 
    std::vector<UShort_t>  fBuffer_Track_ChargedIso; 
    std::vector<UShort_t>  fBuffer_MatchType; 

    std::vector<UShort_t> fBuffer_MC_True_Cluster_E; 
    std::vector<UShort_t> fBuffer_MC_True_Track_E; 
    std::vector<Short_t> fBuffer_MC_True_Track_Px; 
    std::vector<Short_t> fBuffer_MC_True_Track_Py; 
    std::vector<Short_t> fBuffer_MC_True_Track_Pz; 
    std::vector<Int_t> fBuffer_MC_True_Track_MotherPDG; 
    std::vector<Bool_t> fBuffer_MC_Track_Is_Electron; 
    std::vector<Bool_t> fBuffer_MC_Cluster_Is_Electron; 
    std::vector<UShort_t> fBuffer_MC_ClusterTrack_Same_Electron; 
    Float_t fBuffer_MC_JetJetWeight; 
    AliCaloTrackMatcher* fTrackMatcher;
    TString  fTrackMatcherName; // track matcher name used for cut histos etc
  private:
    ULong64_t GetUniqueEventID      ( AliVHeader *header);
    void CountTracks                ();
    void ResetBuffer();
    void ProcessMCCaloPhoton(AliAODCaloCluster* clus,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag);
    void ProcessCaloPhotons();
    void ProcessTracks(); // only needed for track effi
    void ProcessTracksESD(); // only needed for track effi
    Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
    void ProcessMatchedTrack(AliAODTrack* track, AliAODCaloCluster* clus, Bool_t isV0, Float_t dEtaV0 = 99, Float_t dPhiV0 = 99);
    void ProcessMCParticles();
    void ProcessTrackMatching(AliAODCaloCluster* clus);
    Bool_t IsMatchedWithConv(AliAODCaloCluster* clus, AliCaloPhotonCuts* cuts);
    Bool_t IsSameTrack(Int_t id1, Int_t id2); // check if GetID() of both tracks points to same base track
    Bool_t IsInEMCalAcceptance(AliAODConversionPhoton *photon); // check if conv photon is in EMC acc
    Bool_t IsInEMCalAcceptance(AliAODMCParticle *part); // check if mcpart is in emc acceptance
    Bool_t IsInEMCalAcceptance(AliAODTrack *part); // check if mcpart is in emc acceptance
    Int_t CheckClustersForMCContribution(Int_t mclabel, TClonesArray *vclus);
    Int_t CheckConvForMCContribution(Int_t mclabel, TClonesArray *vconv);
    void RelabelAODPhotonCandidates(Bool_t mode);
    Short_t ConvertToShort(Float_t input, Int_t scale);
    Short_t ConvertToShort(Double_t input, Int_t scale);
    UShort_t ConvertToUShort(Float_t input, Int_t scale);
    UShort_t ConvertToUShort(Double_t input, Int_t scale);
    void PushToVectors(treeWriteContainer input);
    std::pair<Double_t,Double_t> ProcessChargedIsolation(AliAODTrack* track);
    AliAnalysisTaskElectronStudies(const AliAnalysisTaskElectronStudies&); // Prevent copy-construction
    AliAnalysisTaskElectronStudies& operator=(const AliAnalysisTaskElectronStudies&); // Prevent assignment  
    ClassDef(AliAnalysisTaskElectronStudies, 14);

};

const Int_t kShortScaleLow = 100; // used for ocnversion
const Int_t kShortScaleMiddle = 1000;
const Int_t kShortScaleHigh = 10000;
#endif

