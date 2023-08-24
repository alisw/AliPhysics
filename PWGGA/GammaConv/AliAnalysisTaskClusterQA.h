#ifndef AliAnalysisClusterQA_cxx
#define AliAnalysisClusterQA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskClusterQA : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskClusterQA();
    AliAnalysisTaskClusterQA(const char *name);
    virtual ~AliAnalysisTaskClusterQA();

    virtual void   UserCreateOutputObjects  ();
    virtual Bool_t Notify                   ();
    virtual void   UserExec                 ( Option_t *option );
    virtual void   Terminate                ( Option_t * );

    void SetV0Reader                        ( AliV0ReaderV1 *v0Reader )                   { fV0Reader=v0Reader                  ; }
    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
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
    // void SetMesonCuts                     ( AliConversionMesonCuts* mesonCuts,
    //                                           Bool_t IsHeavyIon )                         {
    //                                                                                         fMesonCuts=mesonCuts           ;
    //                                                                                         fIsHeavyIon = IsHeavyIon            ;
    //                                                                                       }
    void FillType                           ( Double_t fillTree,
                                              Bool_t fillHistorams)                       {
                                                                                            ffillTree = fillTree                ;
                                                                                            ffillHistograms = fillHistorams     ;
                                                                                          }
    void SetIsMC                            ( Int_t isMC )                                { fIsMC                 = isMC        ; }
    void SetDoAdditionalHistos              ( Bool_t val )                                { fSaveAdditionalHistos = val        ; }
    void SetSaveEventProperties             ( Bool_t val  )                               { fSaveEventProperties  = val         ; }
    void SetSaveClusterCells                ( Bool_t val  )                               { fSaveCells            = val         ; }
    void SetSaveSurroundingCells            ( Bool_t val  )                               { fSaveSurroundingCells = val         ; }
    void SetSaveSurroundingTracks           ( Bool_t val  )                               { fSaveTracks           = val         ; }
    void SetMaxConeRadius                   ( Float_t val  )                              { fConeRadius           = val         ; }
    void SetMinTrackPt                      ( Float_t val  )                              { fMinTrackPt           = val         ; }
    void SetMinClusterEnergy                ( Float_t val  )                              { fMinClusterEnergy     = val         ; }
    void SetMinMaxNLMCut                    ( Int_t valmin, Int_t valmax  )               { fMinNLMCut            = valmin      ;
                                                                                            fMaxNLMCut            = valmax      ; }
    void SetSaveMCInformation               ( Bool_t val  )                               { fSaveMCInformation    = val         ; }
    void SetEventwiseClusterOutput          ( Bool_t val )                                { fSaveEventsInVector   = val         ; }
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
    Int_t       FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event);
    void GetRowAndColumnFromAbsCellID(Int_t cellIndex, Int_t& row, Int_t& column);
    Int_t MakePhotonCandidates(AliVCluster* clus, AliVCaloCells* cells, Long_t indexCluster);
    void ProcessTracksAndMatching(AliVCluster* clus, Long_t indexCluster);
    Int_t  GetMCClusterFlag(AliVCluster* clus, AliVCaloCells* cells);
    Float_t GetCentrality(AliVEvent *event);
    //  Int_t       GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event);
    // Int_t       GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event,  Int_t *absCellIdList, Float_t* maxEList);
  private:

    AliAnalysisTaskClusterQA     ( const AliAnalysisTaskClusterQA& ); // Prevent copy-construction
    AliAnalysisTaskClusterQA &operator=( const AliAnalysisTaskClusterQA& ); // Prevent assignment

    ULong64_t GetUniqueEventID      ( AliVHeader *header);
    void ProcessQATreeCluster       ( AliVEvent *event, AliVCluster* cluster, Long_t indexCluster);
    void ProcessQA                  ( AliAODConversionPhoton *gamma );
    void RelabelAODPhotonCandidates ( Bool_t mode );
    void ProcessTrueQAESD           ( AliAODConversionPhoton *TruePhotonCandidate,
                                      AliESDtrack *elec,
                                      AliESDtrack *posi );
    void ProcessTrueQAAOD           ( AliAODConversionPhoton *TruePhotonCandidate,
                                      AliAODTrack *elec,
                                      AliAODTrack *posi );
    Int_t ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster,
                                        AliAODConversionPhoton *TrueSubClusterCandidate1,
                                        AliAODConversionPhoton *TrueSubClusterCandidate2);
    Int_t ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster,
                                        AliAODConversionPhoton *TrueSubClusterCandidate1,
                                        AliAODConversionPhoton *TrueSubClusterCandidate2);
    UInt_t IsTruePhotonESD          ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t IsTruePhotonAOD          ( AliAODConversionPhoton *TruePhotonCandidate );
    void CountTracks                ();
    void SetLogBinningXTH2          ( TH2* histoRebin );
    void ResetBuffer();
    void ResetBufferVectors();

  protected:
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    TClonesArray*               fReaderGammas;          // Array with conversion photons selected by V0Reader Cut
    AliPIDResponse*             fPIDResponse;                         ///< PID response
    TString                     fCorrTaskSetting;
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliCaloPhotonCuts*          fClusterCutsEMC;               // Cuts used by the V0Reader
    // AliConversionMesonCuts*     fMesonCuts;                  // MesonCutObject
    Int_t                       fMinNLMCut;                   ///< save MC information
    Int_t                       fMaxNLMCut;                   ///< save MC information
    AliVEvent*                  fInputEvent;                //
    AliMCEvent*                 fMCEvent;                   //
    Double_t                    fWeightJetJetMC;                                      // weight for Jet-Jet MC
    AliEMCALGeometry*           fGeomEMCAL;                     // pointer to EMCAL geometry
    TTree*                      fClusterTree;                    //
    Bool_t                      fIsHeavyIon;                //
    Double_t                    ffillTree;                  //
    Bool_t                      ffillHistograms;            //
    TList*                      fOutputList;                //
    Int_t                       fIsMC;                      //
    Bool_t                      fCorrectForNonlinearity;                      //

    // Save flags
    Bool_t          fSaveEventProperties;                 ///< save general event properties (centrality etc.)
    Bool_t          fSaveCells;                           ///< save arrays of cluster cells
    Bool_t          fSaveSurroundingCells;                        ///< save arrays of all cells in event
    Bool_t          fSaveTracks;                        ///< save arrays of all cells in event
    Float_t          fConeRadius;                        ///< save arrays of all cells in event
    Float_t          fMinTrackPt;                        ///< save arrays of all cells in event
    Float_t          fMinClusterEnergy;                        ///< save arrays of all cells in event
    Bool_t          fSaveMCInformation;                   ///< save MC information
    Bool_t          fSaveAdditionalHistos;                   ///< save MC information
    Bool_t          fSaveEventsInVector;                    ///< save cluster information in event vectors information

    // Option flags
    std::vector<Float_t> fExtractionPercentages;          ///< Percentages which will be extracted for a given pT bin
    std::vector<Float_t> fExtractionPercentagePtBins;     ///< pT-bins associated with fExtractionPercentages

    // Buffers that will be added to the tree
    Float_t         fBuffer_EventWeight;                     //!<! array buffer
    Float_t         fBuffer_ClusterE;                     //!<! array buffer
    Float_t         fBuffer_ClusterPhi;                   //!<! array buffer
    Float_t         fBuffer_ClusterEta;                   //!<! array buffer
    Bool_t          fBuffer_ClusterIsEMCAL;                   //!<! array buffer
    Int_t           fBuffer_ClusterSupMod;                   //!<! array buffer
    Int_t           fBuffer_MC_Cluster_Flag;                   //!<! array buffer
    Int_t           fBuffer_ClusterNumCells;              //!<! array buffer
    // Int_t           fBuffer_ClusterNLM;              //!<! array buffer
    Int_t           fBuffer_LeadingCell_ID;              //!<! array buffer
    Float_t         fBuffer_LeadingCell_E;              //!<! array buffer
    Float_t         fBuffer_LeadingCell_Eta;              //!<! array buffer
    Float_t         fBuffer_LeadingCell_Phi;              //!<! array buffer
    Float_t         fBuffer_ClusterM02;              //!<! array buffer
    Float_t         fBuffer_ClusterM20;              //!<! array buffer

    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Multiplicity;             //!<! array buffer
    Int_t           fBuffer_Event_NumActiveCells;          //!<! array buffer

    // Int_t*          fBuffer_ClusterNLM_ID;                      //!<! array buffer
    // Float_t*          fBuffer_ClusterNLM_E;                      //!<! array buffer
    Int_t*          fBuffer_Cells_ID;                      //!<! array buffer
    Float_t*        fBuffer_Cells_E;                      //!<! array buffer
    Float_t*        fBuffer_Cells_RelativeEta;                      //!<! array buffer
    Float_t*        fBuffer_Cells_RelativePhi;                      //!<! array buffer

    Int_t           fBuffer_Surrounding_NCells;                //!<! array buffer
    Int_t*          fBuffer_Surrounding_Cells_ID;                //!<! array buffer
    Float_t*          fBuffer_Surrounding_Cells_Time;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_R;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_E;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_RelativeEta;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_RelativePhi;              //!<! array buffer
    Int_t           fBuffer_Surrounding_NTracks;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_R;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_Pt;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_P;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_nSigdEdxE;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_RelativeEta;              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_RelativePhi;              //!<! array buffer
    Bool_t*         fBuffer_Surrounding_Tracks_V0Flag;              //!<! array buffer
    Float_t*         fBuffer_Surrounding_Tracks_TOF;              //!<! array buffer

    Int_t           fBuffer_Cluster_MC_Label;              //!<! array buffer
    Int_t           fBuffer_Mother_MC_Label;              //!<! array buffer
    Float_t         fBuffer_Cluster_MC_EFracFirstLabel;              //!<! array buffer
    Float_t         fBuffer_Cluster_MC_TrueEFirstLabel;              //!<! array buffer
    Int_t         fBuffer_Cluster_MC_FirstLabel;              //!<! array buffer
    Float_t         fBuffer_Cluster_MC_EFracLeadingPi0;              //!<! array buffer
    Float_t         fBuffer_Cluster_MC_LeadingPi0_Pt;              //!<! array buffer
    Float_t         fBuffer_Cluster_MC_LeadingPi0_E;              //!<! array buffer



    // vector buffers for storing eventwise information

    std::vector<Float_t>    fVBuffer_Cluster_E;                         //!<! vector buffer
    std::vector<Float_t>    fVBuffer_Cluster_Eta;                         //!<! vector buffer
    std::vector<Float_t>    fVBuffer_Cluster_Phi;                         //!<! vector buffer
    std::vector<Bool_t>     fVBuffer_Cluster_isEMCal;                        //!<! vector buffer
    std::vector<Int_t>      fVTrueNeutralPionDaughterIndex;                    //!<! vector buffer   store the MC stack ID of mother pi0 for true information


    ClassDef(AliAnalysisTaskClusterQA, 16);
};

const Int_t kMaxActiveCells = 18000;
const Int_t kMaxNTracks = 4000;

#endif
