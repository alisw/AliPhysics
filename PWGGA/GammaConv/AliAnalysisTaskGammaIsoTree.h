#include "AliAODConversionPhoton.h"
#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
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
  UInt_t evtType;
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

// small class for extra information on clusters
class AliExtraClusterInfoHelper : public TObject{
    public:
      Short_t nLM,matchedTrackIndex;
      Float_t exoticEFrac;
      AliExtraClusterInfoHelper(Short_t pNLM, Short_t pmatchedTrackIndex,Float_t pexoticEFrac)
       : nLM(pNLM),
         matchedTrackIndex(pmatchedTrackIndex),
         exoticEFrac(pexoticEFrac)
      {}

      AliExtraClusterInfoHelper() : AliExtraClusterInfoHelper(0,0,0) {}
      Bool_t isMatched(){
        if(matchedTrackIndex!=-1){
          return kTRUE;
        } else{
          return kFALSE;
        }
      }

      // copy constructor
      AliExtraClusterInfoHelper(const AliExtraClusterInfoHelper & original) : nLM(original.nLM),matchedTrackIndex(original.matchedTrackIndex),exoticEFrac(original.exoticEFrac){}
      ClassDef(AliExtraClusterInfoHelper,1);
};

// small class for details on isolation and tagging
class AliIsoInfoHelper : public TObject {
  public:
    Double32_t isoRawCharged[2],isoRawNeutral[2], isoCell[2]; // storage for two isolation radii each
    Int_t isTagged; //0 : no 1:withOtherConv 2: withOtherCluster 3: both
    
    AliIsoInfoHelper() 
    : isTagged(0)
    {
        isoRawCharged[0] = -1;
        isoRawCharged[1] = -1;
        isoRawNeutral[0] = -1;
        isoRawNeutral[1] = -1;
        isoCell[0] = -1;
        isoCell[1] = -1;
    }
    
    AliIsoInfoHelper(Double32_t isoRawCh[2], Double32_t isoRawNeut[2],Double32_t isoC[2], Int_t tagged) 
     : isTagged(tagged)
    {
        isoRawCharged[0] = isoRawCh[0];
        isoRawCharged[1] = isoRawCh[1];
        isoRawNeutral[0] = isoRawNeut[0];
        isoRawNeutral[1] = isoRawNeut[1];
        isoCell[0] = isoC[0];
        isoCell[1] = isoC[1];
    }
    // copy construction
    AliIsoInfoHelper(const AliIsoInfoHelper & original) : isTagged(original.isTagged)
      {
        isoRawCharged[0] = original.isoRawCharged[0];
        isoRawCharged[1] = original.isoRawCharged[1];
        isoRawNeutral[0] = original.isoRawNeutral[0];
        isoRawNeutral[1] = original.isoRawNeutral[1];
        isoCell[0] = original.isoCell[0];
        isoCell[1] = original.isoCell[1];
      }
    ClassDef(AliIsoInfoHelper,1);
};


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
    void SetAllowOverlapHeaders( Bool_t allowOverlapHeader ) {fAllowOverlapHeaders = allowOverlapHeader;}
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
    void SetClusterCutsIsolationEMC             ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                           fClusterCutsIsolationEMC=clusterCuts           ;
                                                                                           fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetClusterCutsTaggingEMC             ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                           fClusterCutsTaggingEMC=clusterCuts           ;
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
    void SetMinClsTPC(Int_t mincls){
        fMinClsTPC = mincls;
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
    void SetDoOwnTrackMatching(Bool_t p){ fDoOwnTrackMatching = p;}

    void SetDoTrackIso(Bool_t p0){ fDoTrackIsolation = p0;}
    void SetTrackIsoR(vector<Float_t> rvec){ fTrackIsolationR = rvec;}
    void SetTrackIsoE(vector<Double_t> evec){ fTrackIsolationE = evec;}
    void SetDoNeutralIso(Bool_t p0){ fDoNeutralIsolation = p0;}
    void SetNeutralIsoR(vector<Float_t> rvec){ fNeutralIsolationR = rvec;}
    void SetNeutralIsoE(vector<Double_t> evec){ fNeutralIsolationE = evec;}
    
    void SetDoCellIso(Bool_t p0){ fDoCellIsolation = p0;}
    void SetRhoOutName(TString s){fRhoOutName = s;}
    void SetBuffSize(Long64_t buff){fTreeBuffSize = buff;}
    void SetDoTagging(Bool_t p0){ fDoTagging = p0;}
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
        fSaveTracks = b;
    }
    void SetUseHistograms(Bool_t b){
        fUseHistograms = b;
    }
    void SetUseTree(Int_t b){
        fUseTree = b;
    }
    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

    void SetAntiIsolationE(Double_t emin, Double_t emax){
        fAntiIsolationE[0] = emin;
        fAntiIsolationE[1] = emax;
    }
    void SetSignalMinM02(Double_t m02){
      fMinM02 = m02;
    }
    void SetSignalMaxM02(Double_t m02){
      fMaxM02 = m02;
    }
  protected:
    AliVEvent*                  fInputEvent;                //!<!
    AliMCEvent*                 fMCEvent;                   //!<!
    Double_t                    fWeightJetJetMC;            //
    TList*                      fOutputList;                //!<!
    TList*                      fConvFolderRec;                //!<!
    TList*                      fConvFolderTrue;                //!<!
    TList*                      fCaloFolderRec;                //!<!
    TList*                      fCaloFolderTrue;                //!<!
    TList*                      fGeneralFolder;                //!<!
    TList*                      fQAFolder;                //!<!
    TList*                      fGeneratorFolder;                //!<!
    TTree*                      fAnalysisTree;              //!<!
    Int_t                       fIsMC;                      //
    Bool_t                      fIsHeavyIon;                //
    AliV0ReaderV1*              fV0Reader;        //!<! V0Reader for basic conversion photon selection
    TString                     fV0ReaderName;    
    TClonesArray*               fReaderGammas;     //!<! array with photon from fV0Reader                      //
    TClonesArray* fConversionCandidates;   //!<! stores conv candidates of event that fulfill cuts
    TClonesArray* fClusterEMCalCandidates;    //!<! stores emcal clusters that fulfill cuts
    TClonesArray* fClusterEMCalCandidatesIsolation;   //!<! vector containing clusters used for isolation, for internal use only
    TClonesArray* fClusterEMCalCandidatesTagging;   //!<! vector containing clusters used for tagging, internal use only
    TClonesArray* fClusterPHOSCandidates;   //!<! stores phos clusters that fulfill cuts
    TClonesArray* fTracks;   //!<!
    TClonesArray* fMCParticles;   //!<! stores mc particles
    TClonesArray* fAODMCTrackArray;    // storage of track array
    TClonesArray* fExtraClusterInfo;  //!<! ID of up to 5 tracks per cluster, where index of vector corresponds to emc candidates index
    TClonesArray* fExtraClusterInfoBackground;  //!<! ID of up to 5 tracks per cluster, where index of vector corresponds to emc candidates index
    dEvtHeader                  fDataEvtHeader;  //!<! storage for general event properties
    mcEvtHeader                 fMCEvtHeader;    //!<! storage for MC event properties
    TClonesArray*        fConvIsoInfo;    //!<! storage for isolation info of conv photons, following same ordering as fConversionCandidates
    TClonesArray*        fCaloIsoInfo;    //!<! storage for isolation of EMC clusters, following same ordering as fConversionCandidates
   
    AliEMCALGeometry*           fGeomEMCAL;    // pointer to EMCAL geometry
    
    // cuts and setting
    TString                     fCorrTaskSetting;           //
    AliConvEventCuts*           fEventCuts;                 // event cuts
    AliCaloPhotonCuts*          fClusterCutsEMC;            // emc cluster cuts used for signal clusters (clusters that are stored to tree)
    AliCaloPhotonCuts*          fClusterCutsIsolationEMC;  // emc cluster cuts used for background clusters (used for tagging and isolation, not stored)
    AliCaloPhotonCuts*          fClusterCutsTaggingEMC;  // emc cluster cuts used for background clusters (used for tagging and isolation, not stored)
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
    Bool_t                      fDoOwnTrackMatching; // flag to enable own track matching instead of track matching provided by AliCaloPhotonCuts

    Bool_t                      fDoTrackIsolation; //
    std::vector<Float_t>        fTrackIsolationR;  //
    std::vector<Double_t>       fTrackIsolationE;  //

    Bool_t                      fDoNeutralIsolation; //
    std::vector<Float_t>        fNeutralIsolationR; //
    std::vector<Double_t>       fNeutralIsolationE; //
    Bool_t                      fDoCellIsolation; //
    
    Bool_t                      fDoTagging; //
    Double_t                    fPi0TaggingWindow[2];    // inv mass window used for pi0 tagging
    Double_t                    fEtaTaggingWindow[2];    // inv mass window used for eta tagging
        
    Bool_t                      fSaveConversions; //
    Bool_t                      fSaveEMCClusters; //
    Bool_t                      fSavePHOSClusters; //
    Bool_t                      fSaveTracks; //

    Bool_t                      fUseHistograms; // if activated, histograms will be used instead of a tree
    Int_t                       fUseTree; // 0 no tree, 1, light tree, 2 full tree
    // histos
    TH1F*                       fHistoNEvents;   // 
    TH1F*                       fHistoNEventsWOWeight;   // 
    TH1F*                       fHistoChargedIso;   // 
    TH2F*                       fHistoTaggingPCMPCM;   // 
    TH2F*                       fHistoTaggingPCMEMC;   // 
    TH2F*                       fHistoTaggingEMCPCM;   // 
    TH2F*                       fHistoTaggingEMCEMC;   // 

    //
    // ─── CONVERSION HISTOS ───────────────────────────────────────────
    //

    TH1F*                       fConvPt;//
    TH1F*                       fConvPtBeforeAcc;//

    TH1F*                       fConvPtTaggedCalo; //
    TH1F*                       fConvPtTaggedAsDecayCalo; //
    TH2F*                       fConvIsoRawCharged[5]; //
    TH2F*                       fConvIsoRawNeutral[5]; //
    TH2F*                       fConvIsoRawFull[5]; //
    TH2F*                       fConvIsoCell[5]; //
    TH2F*                       fConvIsoCorr[5]; //
    TH1F*                       fConvRho; //
    TH1F*                       fConvRhoTimesArea; //
    // True conv histos
    TH1F*                       fConvTruePt; //
    TH1F*                       fConvTruePtPrimary; //
    TH1F*                       fConvTruePtDecay; //
    TH1F*                       fConvTruePtDecayFoundOtherInCluster; //
    TH1F*                       fConvTruePtDecayOtherInAcc; //
    TH1F*                       fConvTruePtDecayOtherInAccAboveMinEnergy; //
    TH1F*                       fConvTruePtTaggedCalo; //
    TH1F*                       fConvTruePtTaggedAsDecayCalo; //
    TH1F*                       fConvTrueRecPt; //
    TH1F*                       fConvTrueRecPtPrimary; //
    TH1F*                       fConvTrueRecPtDecay; //
    TH1F*                       fConvTrueRecPtDecayFoundOtherInCluster; //
    TH1F*                       fConvTrueRecPtDecayOtherInAcc; //
    TH1F*                       fConvTrueRecPtDecayOtherInAccAboveMinEnergy; //
    TH1F*                       fConvTrueRecPtTaggedCalo; //
    TH1F*                       fConvTrueRecPtTaggedAsDecayCalo; //
    TH2F*                       fConvTrueIsoRawCharged[5]; //
    TH2F*                       fConvTrueIsoRawNeutral[5]; //
    TH2F*                       fConvTrueIsoRawFull[5]; //
    TH2F*                       fConvTrueIsoCorr[5]; //
    TH2F*                       fConvTrueIsoCell[5]; //
    TH2F*                       fConvTrueIsoRawCharged_FromDecay[5]; //
    TH2F*                       fConvTrueIsoRawNeutral_FromDecay[5]; //
    TH2F*                       fConvTrueIsoRawFull_FromDecay[5]; //
    TH2F*                       fConvTrueIsoCell_FromDecay[5]; //
    TH2F*                       fConvTrueIsoRawCharged_FromDirect[5]; //
    TH2F*                       fConvTrueIsoRawNeutral_FromDirect[5]; //
    TH2F*                       fConvTrueIsoRawFull_FromDirect[5]; //
    TH2F*                       fConvTrueIsoCell_FromDirect[5]; //

    TH1F*                       fConvPtIsoCharged[5][5]; // R , Emin
    TH1F*                       fConvPtIsoNeutral[5][5]; //
    TH1F*                       fConvPtIsoFull[5][5]; //
    TH1F*                       fConvPtIsoCell[5][5]; //
    TH1F*                       fConvPtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fConvPtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fConvPtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fConvPtTaggedCaloIsoCell[5][5]; //

    TH1F*                       fConvTruePtIsoCharged[5][5]; // R , Emin
    TH1F*                       fConvTruePtIsoNeutral[5][5]; //
    TH1F*                       fConvTruePtIsoFull[5][5]; //
    TH1F*                       fConvTruePtIsoCell[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoCell[5][5]; //
    TH1F*                       fConvTruePtIsoChargedFromDirect[5][5]; // R , Emin
    TH1F*                       fConvTruePtIsoNeutralFromDirect[5][5]; //
    TH1F*                       fConvTruePtIsoFullFromDirect[5][5]; //
    TH1F*                       fConvTruePtIsoCellFromDirect[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoChargedFromDirect[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoNeutralFromDirect[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoFullFromDirect[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloIsoCellFromDirect[5][5]; //

    TH1F*                       fConvTrueRecPtIsoCharged[5][5]; // R , Emin
    TH1F*                       fConvTrueRecPtIsoNeutral[5][5]; //
    TH1F*                       fConvTrueRecPtIsoFull[5][5]; //
    TH1F*                       fConvTrueRecPtIsoCell[5][5]; //
    TH1F*                       fConvTrueRecPtIsoChargedFromDirect[5][5]; // R , Emin
    TH1F*                       fConvTrueRecPtIsoNeutralFromDirect[5][5]; //
    TH1F*                       fConvTrueRecPtIsoFullFromDirect[5][5]; //
    TH1F*                       fConvTrueRecPtIsoCellFromDirect[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloIsoCell[5][5]; //

    TH1F*                       fConvTruePtMCIsoCharged[5][5]; // R , Emin
    TH1F*                       fConvTruePtMCIsoNeutral[5][5]; //
    TH1F*                       fConvTruePtMCIsoFull[5][5]; //
    TH1F*                       fConvTruePtMCIsoCell[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloMCIsoCharged[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloMCIsoNeutral[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloMCIsoFull[5][5]; //
    TH1F*                       fConvTruePtTaggedCaloMCIsoCell[5][5]; //

    TH1F*                       fConvTrueRecPtMCIsoCharged[5][5]; // R , Emin
    TH1F*                       fConvTrueRecPtMCIsoNeutral[5][5]; //
    TH1F*                       fConvTrueRecPtMCIsoFull[5][5]; //
    TH1F*                       fConvTrueRecPtMCIsoCell[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloMCIsoCharged[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloMCIsoNeutral[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloMCIsoFull[5][5]; //
    TH1F*                       fConvTrueRecPtTaggedCaloMCIsoCell[5][5]; //

    TH2F*                       fConvInvMass; // inv mass PCM-EMC
    TH2F*                       fConvInvMassChargedIsolated[5][5]; 
    TH2F*                       fConvInvMassAntiChargedIsolated[5];
    TH2F*                       fConvInvMassNeutralIsolated[5][5]; 
    TH2F*                       fConvInvMassAntiNeutralIsolated[5];
    TH2F*                       fConvInvMassCellIsolated[5][5]; 
    TH2F*                       fConvInvMassAntiCellIsolated[5];
    TH2F*                       fConvInvMassFullIsolated[5][5]; 
    TH2F*                       fConvInvMassAntiFullIsolated[5];
    TH2F*                       fConvTrueInvMass; // inv mass of true photon PCM-EMC
    TH2F*                       fConvTrueInvMass_FromDecay;
    TH2F*                       fConvTrueInvMass_FromDirect;

    //
    // ─── CALO HISTOS ─────────────────────────────────────────────────
    //
    TH1F*                       fCaloPt;//
    TH1F*                       fCaloPtBeforeAcc;//

    TH1F*                       fCaloE;//

    TH1F*                       fCaloPtTaggedCalo; //
    TH1F*                       fCaloPtTaggedAsDecayCalo; //
    TH2F*                       fCaloIsoRawCharged[5]; //
    TH2F*                       fCaloIsoRawNeutral[5]; //
    TH2F*                       fCaloIsoRawFull[5]; //
    TH2F*                       fCaloIsoCell[5]; //
    TH2F*                       fCaloIsoCorr[5]; //
    TH1F*                       fCaloRho; //
    TH1F*                       fCaloRhoTimesArea; //
    // True conv histos
    TH1F*                       fCaloTruePt; //
    TH1F*                       fCaloTruePtPrimary; //
    TH1F*                       fCaloTruePtDecay; //
    TH1F*                       fCaloTruePtDecayFoundOtherInCluster; //
    TH1F*                       fCaloTruePtDecayOtherInAcc; //
    TH1F*                       fCaloTruePtDecayOtherInAccAboveMinEnergy; //
    TH1F*                       fCaloTruePtTaggedCalo; //
    TH1F*                       fCaloTruePtTaggedAsDecayCalo; //
    TH1F*                       fCaloTrueRecPt; //
    TH1F*                       fCaloTrueRecPtPrimary; //
    TH1F*                       fCaloTrueRecPtDecay; //
    TH1F*                       fCaloTrueRecPtDecayFoundOtherInCluster; //
    TH1F*                       fCaloTrueRecPtDecayOtherInAcc; //
    TH1F*                       fCaloTrueRecPtDecayOtherInAccAboveMinEnergy; //
    TH1F*                       fCaloTrueRecPtTaggedCalo; //
    TH1F*                       fCaloTrueRecPtTaggedAsDecayCalo; //
    TH2F*                       fCaloTrueIsoRawCharged[5]; //
    TH2F*                       fCaloTrueIsoRawNeutral[5]; //
    TH2F*                       fCaloTrueIsoRawFull[5]; //
    TH2F*                       fCaloTrueIsoCorr[5]; //
    TH2F*                       fCaloTrueIsoCell[5]; //
    TH2F*                       fCaloTrueIsoRawCharged_FromDecay[5]; //
    TH2F*                       fCaloTrueIsoRawNeutral_FromDecay[5]; //
    TH2F*                       fCaloTrueIsoRawFull_FromDecay[5]; //
    TH2F*                       fCaloTrueIsoCell_FromDecay[5]; //
    TH2F*                       fCaloTrueIsoRawCharged_FromDirect[5]; //
    TH2F*                       fCaloTrueIsoRawNeutral_FromDirect[5]; //
    TH2F*                       fCaloTrueIsoRawFull_FromDirect[5]; //
    TH2F*                       fCaloTrueIsoCell_FromDirect[5]; //

    TH1F*                       fCaloPtIsoCharged[5][5]; // R , Emin
    TH1F*                       fCaloPtIsoNeutral[5][5]; //
    TH1F*                       fCaloPtIsoFull[5][5]; //
    TH1F*                       fCaloPtIsoCell[5][5]; //
    TH1F*                       fCaloPtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fCaloPtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fCaloPtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fCaloPtTaggedCaloIsoCell[5][5]; //

    TH1F*                       fCaloTruePtIsoCharged[5][5]; // R , Emin
    TH1F*                       fCaloTruePtIsoNeutral[5][5]; //
    TH1F*                       fCaloTruePtIsoFull[5][5]; //
    TH1F*                       fCaloTruePtIsoCell[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoCell[5][5]; //
    TH1F*                       fCaloTruePtIsoChargedFromDirect[5][5]; // R , Emin
    TH1F*                       fCaloTruePtIsoNeutralFromDirect[5][5]; //
    TH1F*                       fCaloTruePtIsoFullFromDirect[5][5]; //
    TH1F*                       fCaloTruePtIsoCellFromDirect[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoChargedFromDirect[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoNeutralFromDirect[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoFullFromDirect[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloIsoCellFromDirect[5][5]; //

    TH1F*                       fCaloTrueRecPtIsoCharged[5][5]; // R , Emin
    TH1F*                       fCaloTrueRecPtIsoNeutral[5][5]; //
    TH1F*                       fCaloTrueRecPtIsoFull[5][5]; //
    TH1F*                       fCaloTrueRecPtIsoCell[5][5]; //
    TH1F*                       fCaloTrueRecPtIsoChargedFromDirect[5][5]; // R , Emin
    TH1F*                       fCaloTrueRecPtIsoNeutralFromDirect[5][5]; //
    TH1F*                       fCaloTrueRecPtIsoFullFromDirect[5][5]; //
    TH1F*                       fCaloTrueRecPtIsoCellFromDirect[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloIsoCharged[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloIsoNeutral[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloIsoFull[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloIsoCell[5][5]; //


    TH1F*                       fCaloTruePtMCIsoCharged[5][5]; // R , Emin
    TH1F*                       fCaloTruePtMCIsoNeutral[5][5]; //
    TH1F*                       fCaloTruePtMCIsoFull[5][5]; //
    TH1F*                       fCaloTruePtMCIsoCell[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloMCIsoCharged[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloMCIsoNeutral[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloMCIsoFull[5][5]; //
    TH1F*                       fCaloTruePtTaggedCaloMCIsoCell[5][5]; //

    TH1F*                       fCaloTrueRecPtMCIsoCharged[5][5]; // R , Emin
    TH1F*                       fCaloTrueRecPtMCIsoNeutral[5][5]; //
    TH1F*                       fCaloTrueRecPtMCIsoFull[5][5]; //
    TH1F*                       fCaloTrueRecPtMCIsoCell[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloMCIsoCharged[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloMCIsoNeutral[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloMCIsoFull[5][5]; //
    TH1F*                       fCaloTrueRecPtTaggedCaloMCIsoCell[5][5]; //

    TH2F*                       fCaloInvMass; // inv mass EMC-EMC
    TH2F*                       fCaloInvMassChargedIsolated[5][5]; //
    TH2F*                       fCaloInvMassAntiChargedIsolated[5];//
    TH2F*                       fCaloInvMassNeutralIsolated[5][5]; //
    TH2F*                       fCaloInvMassAntiNeutralIsolated[5];//
    TH2F*                       fCaloInvMassCellIsolated[5][5]; //
    TH2F*                       fCaloInvMassAntiCellIsolated[5];//
    TH2F*                       fCaloInvMassFullIsolated[5][5]; //
    TH2F*                       fCaloInvMassAntiFullIsolated[5];//
    TH2F*                       fCaloTrueInvMass; // inv mass of true photon EMC-EMC
    TH2F*                       fCaloTrueInvMass_FromDecay;//
    TH2F*                       fCaloTrueInvMass_FromDirect;//

    THnSparseF*                       fCaloM02; // inv mass EMC-EMC
    THnSparseF*                       fCaloM02ChargedIsolated[5][5]; //
    THnSparseF*                       fCaloM02AntiChargedIsolated[5];//
    THnSparseF*                       fCaloM02NeutralIsolated[5][5]; //
    THnSparseF*                       fCaloM02AntiNeutralIsolated[5];//
    THnSparseF*                       fCaloM02CellIsolated[5][5]; //
    THnSparseF*                       fCaloM02AntiCellIsolated[5];//
    THnSparseF*                       fCaloM02FullIsolated[5][5]; //
    THnSparseF*                       fCaloM02AntiFullIsolated[5];//
    THnSparseF*                       fCaloTrueM02; // inv mass of true photon EMC-EMC
    THnSparseF*                       fCaloTrueM02_FromDecay;//
    THnSparseF*                       fCaloTrueM02_FromDirect;//

    //
    // ─── GENERATOR LEVEL HISTOS ──────────────────────────────────────
    //

    TH1I*                       fHistoMCHeaders;                                      //! array of histos for header names
    TH1F*                       fGenPhotonPt;//
    TH1F*                       fGenPhotonPt_FromDecay;//
    TH1F*                       fGenPhotonPt_FromDirect;//
    TH1F*                       fGenPhotonPtInEMCalAcc;//
    TH1F*                       fGenPhotonPtInEMCalAcc_FromDecay;//
    TH1F*                       fGenPhotonPtInEMCalAcc_FromDirect;//
    TH1F*                       fGenPhotonPtFoundNormCluster;//
    TH1F*                       fGenPhotonPtFoundTaggingCluster;//
    TH1F*                       fGenPhotonPtFoundIsoCluster;//
    TH2F*                       fGenPhotonEFoundNoClusterVsCellE;//
    TH1F*                       fGenPi0Pt;//
    TH1F*                       fGenPi0PtInEMCalAcc;//
    TH1F*                       fGenPi0PtInEMCalAcc_BothGammaInEMCal;//
    TH1F*                       fGenPi0PtInEMCalAcc_BothGammaInClusters;//

    TString                     fRhoOutName; // 

    Long64_t                    fTreeBuffSize;           ///< allowed uncompressed buffer size per tree
    Long64_t                    fMemCountAOD;            //!<! accumulated tree size before AutoSave

    Int_t                       fTrackMatcherRunningMode; // CaloTrackMatcher running mode

    Double_t                    fAntiIsolationE[2];
    Double_t                    fMinM02; // min m02 for signal clusters (separate from normal cuts to allow purity estimation)
    Double_t                    fMaxM02; // max m02 for signal clusters (separate from normal cuts to allow purity estimation)
    
    // MC cluster & headers 
    Bool_t                fIsFromDesiredHeader;                                 // flag for MC headers
    Bool_t                fIsOverlappingWithOtherHeader;                        // flag for particles in MC overlapping between headers
    Bool_t                fAllowOverlapHeaders;                                 // enable overlapping headers for cluster selection
    // //
    // // ─── FOR LIGHT TREE ──────────────────────────────────────────────
    // //

    // Double32_t                  fBuffCaloPt; 
    // Double32_t                  fBuffCaloChargedIso; 
    // Double32_t                  fBuffCaloNeutralIso; 
    // Double32_t                  fBuffCaloFullIso; 
    // Double32_t                  fBuffCaloCellIso; 
    // Bool_t                      fBuffCaloIsDecay; 
    // Bool_t                      fBuffCaloIsDirect; 
    // Bool_t                      fBuffCaloIsPhoton; 
    // Double32_t                  fBuffCaloM02; 
    // Double32_t                  fBuffCaloSubClusInvMass; 


  private:
    ULong64_t GetUniqueEventID      ( AliVHeader *header);
    void CountTracks                ();
    void ResetBuffer();
    void ProcessConversionPhotons();
    void ProcessMCConversionPhoton(AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag);
    void ProcessMCCaloPhoton(AliAODCaloCluster* clus,AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight);
    void ProcessCaloPhotons();
    Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
    void ProcessTracks();
    void ProcessMCParticles();
    Int_t ProcessTrackMatching(AliAODCaloCluster* clus, TClonesArray* tracks);
    vector<Double32_t> ProcessChargedIsolation(AliAODConversionPhoton* photon);
    vector<Double32_t> ProcessChargedIsolation(AliAODCaloCluster* cluster);
    vector<Double32_t> ProcessNeutralIsolation(AliAODConversionPhoton* photon);
    vector<Double32_t> ProcessCellIsolation(AliAODConversionPhoton* photon);
    vector<Double32_t> ProcessCellIsolation(AliAODCaloCluster* cluster);
    vector<Double32_t> ProcessNeutralIsolation(AliAODCaloCluster* cluster);
    vector<Double32_t> ProcessMCIsolation(Int_t mclabel);
    Int_t ProcessTagging(AliAODConversionPhoton* photon);
    Int_t ProcessTagging(AliAODCaloCluster* cluster);
    void ReduceTrackInfo();
    void RelabelAODPhotonCandidates(Bool_t mode);
    void FillConversionHistos(AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag);
    void FillCaloHistosPurity(AliAODCaloCluster* clus, AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight);
    void FillCaloHistos(AliAODCaloCluster* clus, AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight);
    Float_t GetExoticEnergyFraction(AliVCluster *cluster, AliVEvent *event);
    Bool_t IsMatchedWithConv(AliAODCaloCluster* clus, AliCaloPhotonCuts* cuts);
    Bool_t IsSameTrack(Int_t id1, Int_t id2); // check if GetID() of both tracks points to same base track
    Bool_t IsInEMCalAcceptance(AliAODConversionPhoton *photon); // check if conv photon is in EMC acc
    Bool_t IsInEMCalAcceptance(AliAODMCParticle *part); // check if mcpart is in emc acceptance
    Bool_t IsTrueConversionPhoton(AliAODConversionPhoton *photon);
    Int_t GetConvPhotonMCLabel(AliAODConversionPhoton *photon);
    Bool_t IsDecayPhoton(AliAODMCParticle *mcphoton);
    Bool_t IsDecayPhoton(AliAODConversionPhoton *photon);
    Int_t CheckClustersForMCContribution(Int_t mclabel, TClonesArray *vclus);
    Int_t CheckConvForMCContribution(Int_t mclabel, TClonesArray *vconv);
    AliAnalysisTaskGammaIsoTree(const AliAnalysisTaskGammaIsoTree&); // Prevent copy-construction
    AliAnalysisTaskGammaIsoTree& operator=(const AliAnalysisTaskGammaIsoTree&); // Prevent assignment  
    ClassDef(AliAnalysisTaskGammaIsoTree, 22);
};

#endif

