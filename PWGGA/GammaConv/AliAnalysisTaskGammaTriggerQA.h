#ifndef ALIANLYSISTASKGAMMATRIGGERQA_cxx
#define ALIANLYSISTASKGAMMATRIGGERQA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include <vector>
#include <map>

class AliAnalysisTaskGammaTriggerQA : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaTriggerQA();
    AliAnalysisTaskGammaTriggerQA(const char *name);
    virtual ~AliAnalysisTaskGammaTriggerQA();

    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
    void InitBack();

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;
    }

    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}

    // switches for additional analysis streams or outputs
    void SetLightOutput(Bool_t flag){fDoLightOutput = flag;}
    void SetDetailedQAFlag (Int_t flag) {fQADetailed  = flag;}

    // Setting the cut lists for the conversion photons
    void SetEventCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fEventCutArray = CutArray;
    }

    // Setting the cut lists for the calo photons
    void SetCaloCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fClusterCutArray = CutArray;
    }

    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

  protected:
    AliV0ReaderV1*        fV0Reader;                                            // basic photon Selection Task
    TString               fV0ReaderName;
    TString               fCorrTaskSetting;
    AliVEvent*            fInputEvent;                                          // current event
    AliMCEvent*           fMCEvent;                                             // corresponding MC event
    TList**               fCutFolder;                                           // Array of lists for containers belonging to cut
    TList**               fESDList;                                             // Array of lists with histograms with reconstructed properties
    TList*                fOutputContainer;                                     // Output container
    TList*                fClusterCandidates;                                   //! current list of cluster candidates
    TList*                fEventCutArray;                                       // List with Event Cuts
    AliConvEventCuts*     fEventCuts;                                           // EventCutObject
    TList*                fClusterCutArray;                                     // List with Cluster Cuts
    AliCaloPhotonCuts*    fCaloPhotonCuts;                                      // CaloPhotonCutObject

    // histograms for rec photon clusters
    TH1F**                fHistoClusGammaPt;                                    //! array of histos with cluster, pt
    TH1F**                fHistoClusGammaE;                                     //! array of histos with cluster, E

    // event histograms
    TH1F**                fHistoNEvents;                                        //! array of histos with event information
    TH1F**                fHistoNEventsWOWeight;                                //! array of histos with event information without event weights
    TH1F**                fHistoNGoodESDTracks;                                 //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                fHistoCent;                                           //! array of histos with centrality slices
    TH1F**                fHistoVertexZ;                                        //! array of histos with vertex z distribution for selected events
    TH1F**                fHistoNGammaCandidates;                               //! array of histos with number of gamma candidates per event
    TH1F**                fHistoNGammaCandidatesBasic;                          //! array of histos with number of gamma candidates per event for basic cluster cut
    TH2F**                fHistoNGoodESDTracksVsNGammaCandidates;               //! array of histos with number of good tracks vs gamma candidates
    TH2F**                fHistoSPDClusterTrackletBackground;                   //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                fHistoNV0Tracks;                                      //! array of histos with V0 counts
    TH1F**                fHistoNV0Trigger;                                     //! array of histos with V0 trigger
    TH2F**                fHistoNV0TriggerTracks;                               //! array of histos with V0 trigger vs tracks
    TProfile**            fProfileEtaShift;                                     //! array of profiles with eta shift
    TProfile**            fProfileJetJetXSection;                               //! array of profiles with xsection for jetjet
    TH1F**                fHistoJetJetNTrials;                                  //! array of histos with ntrials for jetjet

    // tree
    TList**               fTreeList;                                            // array for lists with Tree
    TTree**               fTreeTriggInfo;                                       //! Array of lists with tree
    Float_t               fCent;                                                // cent
    Short_t               fT0Trigg;                                             // T0 trigg info
    UInt_t                fV0Mult;                                              // offline V0 multiplicity
    UInt_t                fV0Trigg;                                             // online V0 multiplicity
    UInt_t                fTPCMult;                                             // TPC multiplicity
    UInt_t                fSPDHit;                                              // SPD hit multiplicity
    UInt_t                fSPDTracklet;                                         // SPD tracklet multiplicity
    Float_t               fZVertex;                                             // Z vertex position

    // additional variables
    Double_t              fEventPlaneAngle;                                     // EventPlaneAngle
    TRandom3              fRandom;                                              // random
    Int_t                 fnCuts;                                               // number of cuts to be analysed in parallel
    Int_t                 fiCut;                                                // current cut
    Int_t                 fIsHeavyIon;                                          // switch for pp = 0, PbPb = 1, pPb = 2
    Bool_t                fDoLightOutput;                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Int_t                 fQADetailed;                                          // switch for detailed QA
    Int_t                 fIsMC;                                                // flag for MC information
    Double_t              fWeightJetJetMC;                                      // weight for Jet-Jet MC
    Int_t                 fNCurrentClusterBasic;                                // current number of cluster without minE

  private:
    AliAnalysisTaskGammaTriggerQA(const AliAnalysisTaskGammaTriggerQA&);                  // Prevent copy-construction
    AliAnalysisTaskGammaTriggerQA &operator=(const AliAnalysisTaskGammaTriggerQA&);       // Prevent assignment

    ClassDef(AliAnalysisTaskGammaTriggerQA, 3);
};

#endif
