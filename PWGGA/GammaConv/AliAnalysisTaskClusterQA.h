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
    void SetClusterCutsEMC                     ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fClusterCutsEMC=clusterCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetClusterCutsDMC                     ( AliCaloPhotonCuts* clusterCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fClusterCutsDMC=clusterCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetMesonCuts                     ( AliConversionMesonCuts* mesonCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fMesonCuts=mesonCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }

    void FillType                           ( Double_t fillTree,
                                              Bool_t fillHistorams)                       {
                                                                                            ffillTree = fillTree                ;
                                                                                            ffillHistograms = fillHistorams     ;
                                                                                          }
    void SetIsMC                            ( Bool_t isMC )                               { fIsMC                 = isMC        ; }
    void SetSaveEventProperties             ( Bool_t val  )                               { fSaveEventProperties  = val         ; }
    void SetSaveClusterCells                ( Bool_t val  )                               { fSaveCells            = val         ; }
    void SetSaveSurroundingCells            ( Bool_t val  )                               { fSaveSurroundingCells         = val         ; }
    void SetNSurroundingCells               ( Int_t val  )                                { fNSurroundingCells    = val         ; }
    void SetMinMaxNLMCut                    ( Int_t valmin, Int_t valmax  )               { fMinNLMCut    = valmin              ;
                                                                                            fMaxNLMCut    = valmax              ; }
    void SetSaveMCInformation               ( Bool_t val  )                               { fSaveMCInformation    = val         ; }
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
    Int_t       FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event);
    void GetRowAndColumnFromAbsCellID(Int_t cellIndex, Int_t& row, Int_t& column);
    Int_t MakePhotonCandidates(AliVCluster* clus, AliVCaloCells* cells, Long_t indexCluster);
    Int_t  GetMCClusterFlag(AliVCluster* clus, AliVCaloCells* cells);
  private:
        
    AliAnalysisTaskClusterQA     ( const AliAnalysisTaskClusterQA& ); // Prevent copy-construction
    AliAnalysisTaskClusterQA &operator=( const AliAnalysisTaskClusterQA& ); // Prevent assignment

    void ProcessQATreeCluster       ( AliVEvent *event, AliVCluster* cluster, Long_t indexCluster);
    void ProcessQA                  ( AliAODConversionPhoton *gamma );
    void RelabelAODPhotonCandidates ( Bool_t mode );
    void ProcessTrueQAESD           ( AliAODConversionPhoton *TruePhotonCandidate,
                                      AliESDtrack *elec,
                                      AliESDtrack *posi );
    void ProcessTrueQAAOD           ( AliAODConversionPhoton *TruePhotonCandidate,
                                      AliAODTrack *elec,
                                      AliAODTrack *posi );
    Int_t ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, Float_t m02,
                                        AliAODConversionPhoton *TrueSubClusterCandidate1,
                                        AliAODConversionPhoton *TrueSubClusterCandidate2);
    UInt_t IsTruePhotonESD          ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t IsTruePhotonAOD          ( AliAODConversionPhoton *TruePhotonCandidate );
    void CountTracks                ();
    void SetLogBinningXTH2          ( TH2* histoRebin );
        
  protected:
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    TString                     fCorrTaskSetting;
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliCaloPhotonCuts*          fClusterCutsEMC;               // Cuts used by the V0Reader
    AliCaloPhotonCuts*          fClusterCutsDMC;               // Cuts used by the V0Reader
    AliConversionMesonCuts*          fMesonCuts;               // Cuts used by the V0Reader
    Int_t           fMinNLMCut;                   ///< save MC information
    Int_t           fMaxNLMCut;                   ///< save MC information
    AliVEvent*                  fInputEvent;                //
    AliMCEvent*                 fMCEvent;                   //
    Double_t                    fWeightJetJetMC;                                      // weight for Jet-Jet MC
    AliEMCALGeometry*           fGeomEMCAL;                     // pointer to EMCAL geometry
    TTree*                      fClusterTree;                    //
    Bool_t                      fIsHeavyIon;                //
    Double_t                    ffillTree;                  //
    Bool_t                      ffillHistograms;            //
    TList*                      fOutputList;                //
    Bool_t                      fIsMC;                      //
    Bool_t                      fCorrectForNonlinearity;                      //
    
    // Save flags
    Bool_t          fSaveEventProperties;                 ///< save general event properties (centrality etc.)
    Bool_t          fSaveCells;                           ///< save arrays of cluster cells
    Bool_t          fSaveSurroundingCells;                        ///< save arrays of all cells in event
    Bool_t          fSaveMCInformation;                   ///< save MC information
    Int_t           fNSurroundingCells;                   ///< save MC information

    // Option flags
    std::vector<Float_t> fExtractionPercentages;          ///< Percentages which will be extracted for a given pT bin
    std::vector<Float_t> fExtractionPercentagePtBins;     ///< pT-bins associated with fExtractionPercentages
    
    // Buffers that will be added to the tree
    Float_t         fBuffer_ClusterE;                     //!<! array buffer
    Float_t         fBuffer_ClusterPhi;                   //!<! array buffer
    Float_t         fBuffer_ClusterEta;                   //!<! array buffer
    Bool_t          fBuffer_ClusterIsEMCAL;                   //!<! array buffer
    Short_t         fBuffer_MC_Cluster_Flag;                   //!<! array buffer
    Short_t         fBuffer_ClusterNumCells;              //!<! array buffer
    Short_t         fBuffer_LeadingCell_ID;              //!<! array buffer
    Short_t         fBuffer_LeadingCell_Row;              //!<! array buffer
    Short_t         fBuffer_LeadingCell_Column;              //!<! array buffer
    Float_t         fBuffer_ClusterM02;              //!<! array buffer
    Float_t         fBuffer_ClusterM20;              //!<! array buffer

    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Multiplicity;             //!<! array buffer
    Short_t         fBuffer_Event_NumActiveCells;          //!<! array buffer

    Short_t*        fBuffer_Cells_ID;                      //!<! array buffer
    Float_t*        fBuffer_Cells_E;                      //!<! array buffer
    Short_t*        fBuffer_Cells_RelativeRow;                      //!<! array buffer
    Short_t*        fBuffer_Cells_RelativeColumn;                      //!<! array buffer

    Short_t*        fBuffer_Surrounding_Cells_ID;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_E;                //!<! array buffer
    Short_t*        fBuffer_Surrounding_Cells_RelativeRow;              //!<! array buffer
    Short_t*        fBuffer_Surrounding_Cells_RelativeColumn;              //!<! array buffer
    
    
    Short_t        fBuffer_Cluster_MC_Label;              //!<! array buffer

    
    ClassDef(AliAnalysisTaskClusterQA, 1);
};

static std::vector<Float_t> DEFAULT_VECTOR_FLOAT;
static std::vector<Short_t> DEFAULT_VECTOR_SHORT;
const Int_t kMaxActiveCells = 500;

#endif

