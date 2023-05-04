#ifndef AliAnalysisCorrelationTree_cxx
#define AliAnalysisCorrelationTree_cxx

class AliPIDResponse;


#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliCaloPhotonCuts.h"
#include "AliConversionPhotonCuts.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskCorrelationTree : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskCorrelationTree();
    AliAnalysisTaskCorrelationTree(const char *name);
    virtual ~AliAnalysisTaskCorrelationTree();

    virtual void   UserCreateOutputObjects  ();
    virtual Bool_t Notify                   ();
    virtual void   UserExec                 ( Option_t *option );
    virtual void   Terminate                ( Option_t * );

    void SetV0Reader                        ( AliV0ReaderV1 *v0Reader )                   { fV0Reader=v0Reader                  ; }
    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetConversionCuts                  ( AliConversionPhotonCuts* conversionCuts, 
                                              Bool_t IsHeavyIon )                         {
                                                                                            fConversionCuts=conversionCuts      ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetEventCuts                       ( AliConvEventCuts* conversionCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fEventCuts=conversionCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetCaloCuts                       ( AliCaloPhotonCuts* caloCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fClusterCuts=caloCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }

    void SetCorrectionTaskSetting(TString setting)        { fCorrTaskSetting = setting                                                                    ; }

    void SetIsMC                            ( Bool_t isMC )                               { fIsMC = isMC                        ; }

  private:
        
    AliAnalysisTaskCorrelationTree     ( const AliAnalysisTaskCorrelationTree& ); // Prevent copy-construction
    AliAnalysisTaskCorrelationTree &operator=( const AliAnalysisTaskCorrelationTree& ); // Prevent assignment

    void ProcessQATree              ( AliAODConversionPhoton *gamma );
    void ProcessQA                  ( AliAODConversionPhoton *gamma );
    void RelabelAODPhotonCandidates ( Bool_t mode );
    void ProcessMuons (  );
    void ProcessElectrons (  );
    void ProcessClusters (  );
    void ProcessTrueQAAOD           ( AliAODConversionPhoton *TruePhotonCandidate, 
                                      AliAODTrack *elec,
                                      AliAODTrack *posi );
    UInt_t IsTruePhotonAOD          ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t GetTrueMotherInfoAOD     ( AliAODConversionPhoton *TruePhotonCandidate );
    void ResetBuffer();
        
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliVEvent*                  fInputEvent;                //
    AliMCEvent*                 fMCEvent;                   //
      TClonesArray*   fAODMCTrackArray;                   // pointer to track array
      TClonesArray*           farrClustersProcess;                                // Cluster array
    TString                 fCorrTaskSetting;
    AliCaloPhotonCuts*           fClusterCuts;                 // Cuts used by the V0Reader

    TTree*                      fAnalysisTree;                    //
    Bool_t                      fIsHeavyIon;                //
    TList*                      fOutputList;                //
    AliPIDResponse 			*fPidResponse;
    TH1F*                  fHistoNEvents;                                      //! array of histos with event information
    Int_t                       fBuffer_NContributors; //
    Float_t                     fBuffer_RunNumber; //
    Int_t                       fBuffer_VertexZ; //
    Int_t                       fBuffer_NEventTriggers; //
    Int_t*                      fBuffer_EventTrigger; //
    Int_t                       fBuffer_NElectronCandidates; //
    Float_t*                    fBuffer_ElectronCandidate_E; //
    Float_t*                    fBuffer_ElectronCandidate_Px; //
    Float_t*                    fBuffer_ElectronCandidate_Py; //
    Float_t*                    fBuffer_ElectronCandidate_Pz; //
    Float_t*                    fBuffer_ElectronCandidate_PropEta; //
    Float_t*                    fBuffer_ElectronCandidate_PropPhi; //
    Float_t*                    fBuffer_ElectronCandidate_Charge; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaElecTPC; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaElecTOF; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaPionTPC; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaPionTOF; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaKaonTPC; //
    Float_t*                    fBuffer_ElectronCandidate_NSigmaProtonTPC; //
    Float_t*                     fBuffer_ElectronCandidate_MC_E;                      //
    Float_t*                     fBuffer_ElectronCandidate_MC_Px;                      //
    Float_t*                     fBuffer_ElectronCandidate_MC_Py;                      //
    Float_t*                     fBuffer_ElectronCandidate_MC_Pz;                      //
    Int_t*                     fBuffer_ElectronCandidate_MC_PDG;                      //
    Float_t*                      fBuffer_ElectronCandidate_MC_Mother_E;
    Float_t*                      fBuffer_ElectronCandidate_MC_Mother_Px;
    Float_t*                      fBuffer_ElectronCandidate_MC_Mother_Py;
    Float_t*                      fBuffer_ElectronCandidate_MC_Mother_Pz;
    Int_t*                      fBuffer_ElectronCandidate_MC_Mother_PDG;
    Int_t*                      fBuffer_ElectronCandidate_MC_GrandMother_PDG;
    Float_t*                      fBuffer_ElectronCandidate_MC_GrandMother_Pt;
    Int_t*                      fBuffer_ElectronCandidate_MC_GGrandMother_PDG;
    Float_t*                      fBuffer_ElectronCandidate_MC_GGrandMother_Pt;
    Int_t                       fBuffer_NMuonCandidates; //
    Float_t*                    fBuffer_MuonCandidate_E; //
    Float_t*                    fBuffer_MuonCandidate_Px; //
    Float_t*                    fBuffer_MuonCandidate_Py; //
    Float_t*                    fBuffer_MuonCandidate_Pz; //
    Float_t*                    fBuffer_MuonCandidate_RAbsEnd; //
    Float_t*                    fBuffer_MuonCandidate_Chi2NDF; //
    Float_t*                    fBuffer_MuonCandidate_DCA; //
    Int_t*                      fBuffer_MuonCandidate_Charge; //
    Int_t*                      fBuffer_MuonCandidate_MatchTrigger; //
    Float_t*                     fBuffer_MuonCandidate_MC_E;                      //
    Float_t*                     fBuffer_MuonCandidate_MC_Px;                      //
    Float_t*                     fBuffer_MuonCandidate_MC_Py;                      //
    Float_t*                     fBuffer_MuonCandidate_MC_Pz;                      //
    Int_t*                     fBuffer_MuonCandidate_MC_PDG;                      //
    Float_t*                      fBuffer_MuonCandidate_MC_Mother_E;
    Float_t*                      fBuffer_MuonCandidate_MC_Mother_Px;
    Float_t*                      fBuffer_MuonCandidate_MC_Mother_Py;
    Float_t*                      fBuffer_MuonCandidate_MC_Mother_Pz;
    Int_t*                      fBuffer_MuonCandidate_MC_Mother_PDG;
    Int_t*                      fBuffer_MuonCandidate_MC_GrandMother_PDG;
    Float_t*                      fBuffer_MuonCandidate_MC_GrandMother_Pt;
    Int_t*                      fBuffer_MuonCandidate_MC_GGrandMother_PDG;
    Float_t*                      fBuffer_MuonCandidate_MC_GGrandMother_Pt;
    Int_t                       fBuffer_NClusterCandidates; //
    Float_t*                    fBuffer_ClusterCandidate_E; //
    Float_t*                    fBuffer_ClusterCandidate_Px; //
    Float_t*                    fBuffer_ClusterCandidate_Py; //
    Float_t*                    fBuffer_ClusterCandidate_Pz; //
    Float_t*                    fBuffer_ClusterCandidate_Eta; //
    Float_t*                    fBuffer_ClusterCandidate_Phi; //
    Float_t*                      fBuffer_ClusterCandidate_M02; //
    Float_t*                     fBuffer_ClusterCandidate_MC_E;                      //
    Float_t*                     fBuffer_ClusterCandidate_MC_Px;                      //
    Float_t*                     fBuffer_ClusterCandidate_MC_Py;                      //
    Float_t*                     fBuffer_ClusterCandidate_MC_Pz;                      //
    Int_t*                     fBuffer_ClusterCandidate_MC_PDG;                      //
    Float_t*                      fBuffer_ClusterCandidate_MC_Mother_E;
    Float_t*                      fBuffer_ClusterCandidate_MC_Mother_Px;
    Float_t*                      fBuffer_ClusterCandidate_MC_Mother_Py;
    Float_t*                      fBuffer_ClusterCandidate_MC_Mother_Pz;
    Int_t*                      fBuffer_ClusterCandidate_MC_Mother_PDG;
    Int_t*                      fBuffer_ClusterCandidate_MC_GrandMother_PDG;
    Bool_t                      fIsMC;                      //

    const Int_t kMaxTracks = 200;
    const Int_t kMaxTriggers = 20;
    
    ClassDef(AliAnalysisTaskCorrelationTree, 1);
};

#endif
