#ifndef AliAnalysisConvCaloTree_cxx
#define AliAnalysisConvCaloTree_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
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


class AliAnalysisTaskConvCaloTree : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskConvCaloTree();
    AliAnalysisTaskConvCaloTree(const char *name);
    virtual ~AliAnalysisTaskConvCaloTree();

    virtual void   UserCreateOutputObjects  ();
    virtual Bool_t Notify                   ();
    virtual void   UserExec                 ( Option_t *option );
    virtual void   Terminate                ( Option_t * );

    void SetV0Reader                        ( AliV0ReaderV1 *v0Reader )                   { fV0Reader           = v0Reader     ; }
    void SetV0ReaderName(TString name)                                                    { fV0ReaderName       = name         ; }
    void SetEventCuts                       ( AliConvEventCuts* conversionCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fEventCuts          = conversionCuts;
                                                                                            fIsHeavyIon         = IsHeavyIon   ;
                                                                                          }
    void SetClusterCutsEMC                  ( AliCaloPhotonCuts* clusterCuts)             { fClusterCutsEMC     = clusterCuts  ; }
    void SetClusterCutsPHOS                 ( AliCaloPhotonCuts* clusterCuts)             { fClusterCutsPHOS    = clusterCuts  ; }
    void SetConversionCuts                  ( AliConversionPhotonCuts* convCuts)          { fConversionCuts     = convCuts     ; }
    void SetMesonCuts                       ( AliConversionMesonCuts* mesonCuts)          { fMesonCuts          = mesonCuts    ; }
    void SetIsMC                            ( Int_t isMC )                                { fIsMC               = isMC         ; }
    void SetSaveMCInformation               ( Bool_t val  )                               { fSaveMCInformation  = val          ; }
    void SetSaveClusters                    ( Bool_t val  )                               { fSaveClusters       = val          ; }
    void SetSaveConversions                 ( Bool_t val  )                               { fSaveConversions    = val          ; }
    void SetSaveTracks                      ( Int_t  val  )                               { fSaveTracks         = val          ; }

    void SetMinTrackPt                      ( Double_t tmp )                              { fMinTrackPt         = tmp          ; }

    void SetCorrectionTaskSetting           (TString setting)                             { fCorrTaskSetting    = setting      ; }

  private:

    AliAnalysisTaskConvCaloTree     ( const AliAnalysisTaskConvCaloTree& ); // Prevent copy-construction
    AliAnalysisTaskConvCaloTree &operator=( const AliAnalysisTaskConvCaloTree& ); // Prevent assignment

    void ProcessClustersAOD();
    void ProcessConversionsAOD();
    void ProcessTracksAOD();
    void ProcessAODMCParticles();
    void FillConversionsTree();
    void RelabelAODPhotonCandidates(Bool_t mode);
    void ResetBuffer();
    void ResetBufferVectors();


  protected:
    AliV0ReaderV1*                  fV0Reader;                      ///< V0 Reader
    TString                         fV0ReaderName;                  ///< V0 Reader Name
    TClonesArray*                   fReaderGammas;                  ///< Array with conversion photons selected by V0Reader Cut
    AliPIDResponse*                 fPIDResponse;                   ///< PID response
    TString                         fCorrTaskSetting;               ///<
    AliVEvent*                      fInputEvent;                    ///< pointer to event
    AliMCEvent*                     fMCEvent;                       ///< pointer to MC event
    Double_t                        fWeightJetJetMC;                ///< weight for Jet-Jet MC
    AliConvEventCuts*               fEventCuts;                     ///< Event cuts
    AliCaloPhotonCuts*              fClusterCutsEMC;                ///< cluster cuts for EDC clusters
    AliCaloPhotonCuts*              fClusterCutsPHOS;               ///< cluster cuts for PHOS clusters
    AliCaloPhotonCuts*              fClusterCuts;                   ///< cluster cuts for either PHOS or EDC clusters
    AliConversionPhotonCuts*        fConversionCuts;                ///< conversion cuts
    AliConversionMesonCuts*         fMesonCuts;                     ///< meson cuts
    TTree*                          fPhotonTree;                    ///<
    Bool_t                          fIsHeavyIon;                    ///<
    TClonesArray*                   fAODMCTrackArray;               ///< AOD track array
    TClonesArray*                   farrClustersProcess;            ///< cluster array
    TList*                          fOutputList;                    ///<
    Int_t                           fIsMC;                          ///<
    std::vector<Int_t>              fMCEventPos;                    ///< AOD relabeling
    std::vector<Int_t>              fMCEventNeg;                    ///< AOD relabeling
    std::vector<Int_t>              fESDArrayPos;                   ///< AOD relabeling
    std::vector<Int_t>              fESDArrayNeg;                   ///< AOD relabeling

    Bool_t                          fSaveMCInformation;             ///< flag to decide to save MC information
    Bool_t                          fSaveClusters;                  ///< flag to decide to save clusters
    Bool_t                          fSaveConversions;               ///< flag to decide to save conversions
    Int_t                           fSaveTracks;                    ///< flag to decide to save tracks

    // track cuts
    Double_t                        fMinTrackPt;                                   //!<! minimum required track pT

    TList*                          fInfoList;                                     //!<! list to store MC and event info
    // histograms for InfoList
    TH1F*                           fHistoNEvents;                                 //!<! number of events
    TH1F*                           fHistoMCPi0Pt;                                 //!<! secondaries in acceptance
    TH1F*                           fHistoMCEtaPt;                                 //!<! secondaries in acceptance

    // Buffers that will be added to the tree
    Float_t                         fBuffer_EventWeight;                           //!<! array buffer to store event weights
    Float_t                         fBuffer_Event_Vertex_Z;                        //!<! array buffer store Vertex Z information


    std::vector<Float_t>            fVBuffer_Cluster_E;                            //!<! vector buffer  Calo (PHOS or EDC) cluster energy
    std::vector<Float_t>            fVBuffer_Cluster_Eta;                          //!<! vector buffer  Calo (PHOS or EDC) cluster eta
    std::vector<Float_t>            fVBuffer_Cluster_Phi;                          //!<! vector buffer  Calo (PHOS or EDC) cluster phi
    std::vector<unsigned int>       fVTrueClusterPi0DaughterIndex;                 //!<! vector buffer   store the MC stack ID of mother pi0 for true information
    std::vector<unsigned int>       fVTrueClusterEtaDaughterIndex;                 //!<! vector buffer   store the MC stack ID of mother eta for true information

    std::vector<Float_t>            fVBuffer_Conv_px;                              //!<! vector buffer  Conversion photon px
    std::vector<Float_t>            fVBuffer_Conv_py;                              //!<! vector buffer  Conversion photon py
    std::vector<Float_t>            fVBuffer_Conv_pz;                              //!<! vector buffer  Conversion photon pz
    std::vector<Float_t>            fVBuffer_Elec1etaCalo;                         //!<! vector buffer  eta of electron 1 from photon conversion on Calo surface
    std::vector<Float_t>            fVBuffer_Elec2etaCalo;                         //!<! vector buffer  eta of electron 2 from photon conversion on Calo surface
    std::vector<Float_t>            fVBuffer_Elec1phiCalo;                         //!<! vector buffer  phi of electron 1 from photon conversion on Calo surface
    std::vector<Float_t>            fVBuffer_Elec2phiCalo;                         //!<! vector buffer  phi of electron 2 from photon conversion on Calo surface
    std::vector<unsigned int>       fVTrueConvPi0DaughterIndex;                    //!<! vector buffer   store the MC stack ID of mother pi0 for true information
    std::vector<unsigned int>       fVTrueConvEtaDaughterIndex;                    //!<! vector buffer   store the MC stack ID of mother pi0 for true information

    std::vector<Short_t>                 fVBuffer_Track_px;                         //!<! vector buffer: track px (*100)
    std::vector<Short_t>                 fVBuffer_Track_py;                         //!<! vector buffer: track py (*100)
    std::vector<Short_t>                 fVBuffer_Track_pz;                         //!<! vector buffer: track pz (*100)
    std::vector<Short_t>                 fVBuffer_Track_P;                          //!<! vector buffer: track P  (*100)
    std::vector<Short_t>                 fVBuffer_Track_Calo_eta;                   //!<! vector buffer: track eta on Calo surface (*10000)
    std::vector<UShort_t>                fVBuffer_Track_Calo_phi;                   //!<! vector buffer: track eta on Calo surface (*10000)

    ClassDef(AliAnalysisTaskConvCaloTree, 1);
};


#endif
