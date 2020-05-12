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

#ifndef AliAnalysisTaskGammaIsoTree_cxx
#define AliAnalysisTaskGammaIsoTree_cxx

typedef struct {
  Double32_t pVtxX,pVtxY,pVtxZ; // primary vertex 
  Int_t runnumber,numberESDtracks;
  Double_t rho;
} dEvtHeader;

typedef struct {
  Double32_t pVtxX,pVtxY,pVtxZ; // primary vertex 
  Int_t runnumber,numberESDtracks;
  Float_t weightJJ;
  Double_t rho;
} mcEvtHeader;

// currently not used, but can be if later simplification is needed
typedef struct {
  TLorentzVector p; 
  Short_t clustertype;
  Int_t nCells, nLM;
  Double_t ToF;
  Double32_t m02,m20,disp;
  Int_t matchedTrackInd[20];
  Int_t nmbMatchedTracks;
  Long_t fCaloPhotonMCLabels[50];
} lightCluster;

typedef struct {
  Double32_t isoRawCharged[2],isoRawNeutral[2], isoCell[2]; // storage for two isolation radii each
  Int_t isTagged; //0 : no 1:withOtherConv 2: withOtherCluster 3: both
} isoInfo;


class AliAnalysisTaskGammaIsoTree : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskGammaIsoTree();
    AliAnalysisTaskGammaIsoTree(const char *name);
    virtual ~AliAnalysisTaskGammaIsoTree();
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
    void SetClusterCutsBackgroundEMC             ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                           fClusterCutsBackgroundEMC=clusterCuts           ;
                                                                                           fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetClusterCutsPHOS                  ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fClusterCutsPHOS=clusterCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          } 
    void SetConvCuts                  ( AliConversionPhotonCuts* convCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fConvCuts=convCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }                                                                                                                                                                        
    void SetYCutMC(Double_t y) {fYMCCut = y;}
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

    void SetEOverP(Double_t p0){ fMatchingEOverP = p0;}
    void SetDoBackgroundTrackMatching(Bool_t p){ fDoBackgroundTrackMatching = p;}

    void SetDoTrackIso(Bool_t p0){ fDoTrackIsolation = p0;}
    void SetTrackIsoR(Float_t r1, Float_t r2){ fTrackIsolationR[0] = r1; fTrackIsolationR[1] = r2;}
    void SetDoNeutralIso(Bool_t p0){ fDoNeutralIsolation = p0;}
    void SetNeutralIsoR(Float_t r1, Float_t r2){ fNeutralIsolationR[0] = r1; fNeutralIsolationR[1] = r2;}

    void SetRhoOutName(TString s){fRhoOutName = s;}
    void SetBuffSize(Long64_t buff){fTreeBuffSize = buff;}

    void SetPi0TaggingWindow(Double_t min,Double_t max=0.){
        fPi0TaggingWindow[0] = min;
        fPi0TaggingWindow[1] = max;
    }
    void SetEtaTaggingWindow(Double_t min,Double_t max=0.){
        fEtaTaggingWindow[0] = min;
        fEtaTaggingWindow[1] = max;
    }

    void SetSaveConversions(Bool_t b){
        fSaveConversions = b;
    }
    void SetSaveEMCClusters(Bool_t b){
        fSaveEMCClusters = b;
    }
    void SetSavePHOSClusters(Bool_t b){
        fSavePHOSClusters = b;
    }
    void SetSaveTracks(Bool_t b){
        fSavePHOSClusters = b;
    }

  protected:
    AliVEvent*                  fInputEvent;                //
    AliMCEvent*                 fMCEvent;                   //
    Double_t                    fWeightJetJetMC;            //
    TList*                      fOutputList;                //
    TTree*                      fAnalysisTree;              //
    Int_t                       fIsMC;
    Bool_t                      fIsHeavyIon;                //
    AliV0ReaderV1*              fV0Reader;        //!<! V0Reader for basic conversion photon selection
    TString                     fV0ReaderName;       ///< Name of the V0 reader
    TClonesArray*               fReaderGammas;     //!<! array with photon from fV0Reader                      //
    std::vector<AliAODConversionPhoton*> fConversionCandidates;  // stores conv candidates of event that fulfill cuts
    std::vector<AliAODCaloCluster*> fClusterEMCalCandidates;   // stores emcal clusters that fulfill cuts
    std::vector<AliAODCaloCluster*> fClusterEMCalCandidatesBackground;  // vector containing clusters used for tagging and isolation, for internal use only
    std::vector<AliAODCaloCluster*> fClusterPHOSCandidates;  // stores phos clusters that fulfill cuts
    //std::vector<lightCluster> fClusterEMCalCandidates;  //
    //std::vector<lightCluster> fClusterPHOSCandidates;  //
    std::vector<AliAODTrack*>   fTracks;  //
    std::vector<AliAODMCParticle*>   fMCParticles;  // stores mc particles
    std::vector<Short_t>        fIDMatchedEMCTrack; // ID of up to 5 tracks per cluster, where index of vector corresponds to emc candidates index
    std::vector<Short_t>        fIDMatchedEMCTrackBackground; // ID of up to 5 tracks per cluster, where index of vector corresponds to emc candidates index
    dEvtHeader                  fDataEvtHeader; // storage for general event properties
    mcEvtHeader                 fMCEvtHeader;   // storage for MC event properties
    std::vector<isoInfo>        fConvIsoInfo;   // storage for isolation info of conv photons, following same ordering as fConversionCandidates
    std::vector<isoInfo>        fCaloIsoInfo;   // storage for isolation of EMC clusters, following same ordering as fConversionCandidates
   
    AliEMCALGeometry*           fGeomEMCAL;                  // pointer to EMCAL geometry
    
    // cuts and setting
    TString                     fCorrTaskSetting;
    AliConvEventCuts*           fEventCuts;                 // event cuts
    AliCaloPhotonCuts*          fClusterCutsEMC;            // emc cluster cuts used for signal clusters (clusters that are stored to tree)
    AliCaloPhotonCuts*          fClusterCutsBackgroundEMC;  // emc cluster cuts used for background clusters (used for tagging and isolation, not stored)
    AliCaloPhotonCuts*          fClusterCutsPHOS;           // phos cluster cuts
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
    Double_t                    fMatchingEOverP; //
    Bool_t                      fDoBackgroundTrackMatching; // should track matching be applied for background clusters (tagging and iso)

    Bool_t                      fDoTrackIsolation; //
    Float_t                     fTrackIsolationR[2];  //

    Bool_t                      fDoNeutralIsolation; //
    Float_t                     fNeutralIsolationR[2]; //
    
    Double_t                    fPi0TaggingWindow[2];    // inv mass window used for pi0 tagging
    Double_t                    fEtaTaggingWindow[2];    // inv mass window used for eta tagging
    
    
    Bool_t                      fSaveConversions; //
    Bool_t                      fSaveEMCClusters; //
    Bool_t                      fSavePHOSClusters; //
    Bool_t                      fSaveTracks; //
    // histos
    TH1F*                       fHistoNEvents;   // 
    TH1F*                       fHistoNEventsWOWeight;   // 
    TH1F*                       fHistoChargedIso;   // 

    TString                     fRhoOutName; // 

    Long64_t                    fTreeBuffSize;           ///< allowed uncompressed buffer size per tree
    Long64_t                    fMemCountAOD;            //!<! accumulated tree size before AutoSave

  private:
    ULong64_t GetUniqueEventID      ( AliVHeader *header);
    void CountTracks                ();
    void ResetBuffer();
    void ProcessConversionPhotons();
    void ProcessCaloPhotons();
    Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
    void ProcessTracks();
    void ProcessMCParticles();
    Int_t ProcessTrackMatching(AliAODCaloCluster* clus, std::vector<AliAODTrack*> tracks);
    void ProcessChargedIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]);
    void ProcessChargedIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]);
    void ProcessNeutralIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]);
    void ProcessCellIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]);
    void ProcessCellIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]);
    void ProcessNeutralIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]);
    Int_t ProcessTagging(AliAODConversionPhoton* photon);
    Int_t ProcessTagging(AliAODCaloCluster* cluster);
    void ReduceTrackInfo();
    void RelabelAODPhotonCandidates(Bool_t mode);

    AliAnalysisTaskGammaIsoTree(const AliAnalysisTaskGammaIsoTree&); // Prevent copy-construction
    AliAnalysisTaskGammaIsoTree& operator=(const AliAnalysisTaskGammaIsoTree&); // Prevent assignment  
    ClassDef(AliAnalysisTaskGammaIsoTree, 5);
};

#endif

