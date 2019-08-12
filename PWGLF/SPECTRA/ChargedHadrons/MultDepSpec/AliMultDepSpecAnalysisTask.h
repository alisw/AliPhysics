/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define MAX_ALLOWED_MULT_BINS 500
#define MAX_HISTO_DIM 4
#define PRECISION 1e-6


#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TRandom3.h"

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"


#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliESDZDC.h"

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliVHeader.h"



#include "AlidNdPtTools.h"
#include "AliMCSpectraWeights.h"

#include <iostream>
using std::string;
using std::vector;
using std::array;

class AliMultDepSpecAnalysisTask : public AliAnalysisTaskSE {
  public:
    AliMultDepSpecAnalysisTask();
    AliMultDepSpecAnalysisTask(const char *name);
    virtual ~AliMultDepSpecAnalysisTask();

    static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 119, Bool_t useDataDrivenCorrections = kFALSE, string pccTrainOutputPath = "",  Int_t pccSysFlag = 0,  Int_t secSysFlag = 0);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t*);

    // Setters
    void SetTriggerMask(UInt_t triggermask)  {fTriggerMask = triggermask;}
    void SetCutMode(Int_t cutMode){fCutMode = cutMode;}
    void SetIsMC(Bool_t isMC = kTRUE){fIsMC = isMC;}
    void SetUseESD(){fIsESD = kTRUE;}
    void SetUseAOD(){fIsESD = kFALSE;}
    void SetMCSpectraWeights(AliMCSpectraWeights* mcSpectraWeights){fMCSpectraWeights = mcSpectraWeights;}
    void SetUseDataDrivenCorrections(Bool_t useDataDrivenCorrections = kTRUE){fMCUseDataDrivenCorrections = useDataDrivenCorrections;}
    void SetUseZDCCut(Bool_t useZDC){fUseZDCCut = useZDC;}
    void SetOverridePbPbEventCuts(Bool_t overridePbPbEventCuts){fOverridePbPbEventCuts = overridePbPbEventCuts;}

    void SetSecScalingSysFlag(Int_t sysFlag = 0){fMCSecScalingSysFlag = sysFlag;}

    // Binning
    void SetBinsPt(Int_t nBins, Double_t* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    void SetBinsEta(Int_t nBins, Double_t* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    void SetBinsMult(Int_t nBins, Double_t* binEdges){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(nBins+1,binEdges);}
    void SetBinsCent(Int_t nBins, Double_t* binEdges){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(nBins+1,binEdges);}
    void SetBinsZv(Int_t nBins, Double_t* binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}
    void SetBinsPtReso(Int_t nBins, Double_t* binEdges){if(fBinsPtReso) delete fBinsPtReso; fBinsPtReso = new TArrayD(nBins+1,binEdges);}
    void SetBinsMult(vector<Int_t> multSteps, vector<Int_t> multBinWidth);
    void SetBinsMult(Int_t maxMult);

    // Acceptance cuts
    void SetMinEta(Double_t minEta){fMinEta = minEta;}
    void SetMaxEta(Double_t maxEta){fMaxEta = maxEta;}
    void SetMinPt(Double_t minPt){fMinPt = minPt;}
    void SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
    void SetMaxZv(Double_t maxZv)  {fMaxZv = maxZv;}
    void SetMaxCent(Double_t maxCent)  {fUseCent = kTRUE; fMaxCent = maxCent;}
    void SetMinCent(Double_t minCent)  {fUseCent = kTRUE; fMinCent = minCent;}

  private:

    //
    TList*              fOutputList;		  //!<! Output list
    AliEventCuts        fEventCuts;       //!<! Event cuts
    AliESDtrackCuts*    fESDtrackCuts;    //!<! Track cuts
    TRandom3*           fRand;            //!<! Random generator
    AliMCSpectraWeights* fMCSpectraWeights;            //-> MC spectra weights object

    Int_t               fCutMode;         ///< ID of track cut variation (100=default)
    Bool_t              fIsESD;			      ///< Flag for ESD usage
    Bool_t              fIsMC;            ///< Flag for MC usage
    Bool_t              fUseCent;         ///< Flag for Centrality usage
    Bool_t              fUseZDCCut;         ///< Flag for zdc cut usage
    Bool_t              fOverridePbPbEventCuts;         ///< override centrality cut in PbPb
    Bool_t              fMCUseDataDrivenCorrections; ///< Flag for data driven corrections usage
    Int_t               fMCSecScalingSysFlag; ///< Flag for secondary scaling systematics 0: nominal, -1,1 variations
    // Cuts
    UInt_t                fTriggerMask;   ///< Trigger mask
    Double_t              fMinEta;        ///< Minimum eta cut
    Double_t              fMaxEta;        ///< Maximum eta cut
    Double_t              fMinPt;			    ///< Minimum pT cut
    Double_t              fMaxPt;			    ///< Maximum pT cut
    Double_t              fMaxZv;			    ///< Maximum absolute z vertex cut
    Double_t              fMinCent;       ///< Minimum centrality
    Double_t              fMaxCent;       ///< Maximum centrality

    // Binning
    TArrayD*             fBinsEventCuts;       ///< Array of bins for event cuts
    TArrayD*             fBinsMult;       ///< Array of bins in multiplicity
    TArrayD*             fBinsCent;       ///< Array of bins in centrality
    TArrayD*             fBinsPt;			    ///< Array of bins in pt
    TArrayD*             fBinsEta;		    ///< Array of bins in eta
    TArrayD*             fBinsZv;			    ///< Array of bins in Zv (Z-position of primary vtx)
    TArrayD*             fBinsPtReso;     ///< Array of bins for relative pt resoulution

    // Output Histograms
    THnSparseF* fHistEventSelection;      //!<! Histogram of event selection
    THnSparseF* fHistEvents;              //!<! Histogram of measured event distribution
    THnSparseF* fHistTracks;              //!<! Histogram of measured tracks
    THnSparseF* fHistRelPtReso;           //!<! Histogram of relatvie pT resolution from covariance matrix

    THnSparseF* fHistMCEventEfficiency;   //!<! Histogram of selelcted events vs Nch
    THnSparseF* fHistMCEventEfficiencyScaled; //!<! Histogram of selelcted events vs Nch

    THnSparseF* fHistMCRelPtReso;         //!<! Histogram of relative pt resolution from mc
    THnSparseF* fHistMCMultCorrelMatrix;  //!<! Histogram of multilicity correlation
    THnSparseF* fHistMCPtCorrelMatrix;    //!<! Histogram of pT correlation
    THnSparseF* fHistMCEtaCorrelMatrix;   //!<! Histogram of eta correlation
    THnSparseF* fHistMCPrimTrue;          //!<! Histogram of generated primaries
    THnSparseF* fHistMCPrimMeas;          //!<! Histogram of measured primaries
    THnSparseF* fHistMCSecMeas;           //!<! Histogram of measured secondaries
    THnSparseF* fHistMCEdgeContam;        //!<! Histogram of tracks from particles out of acceptance
    TH1D* fHistMCDoubleCountig;           //!<! Histogram to track double counting



    THnSparseF* fHistMCEventsScaled;  //!<! Histogram
    THnSparseF* fHistMCTracksScaled;  //!<! Histogram
    THnSparseF* fHistMCMultCorrelMatrixScaled; //!<! Histogram of scaled multilicity correlation
    THnSparseF* fHistMCPrimTrueScaled;  //!<! Histogram
    THnSparseF* fHistMCPrimMeasScaled;  //!<! Histogram
    THnSparseF* fHistMCSecMeasScaled; //!<! Histogram
    THnSparseF* fHistMCEdgeContamScaled;  //!<! Histogram

    THnSparseF* fHistMCMultMeasScaleEffect;   //!<! Histogram
    THnSparseF* fHistMCMultTrueScaleEffect;   //!<! Histogram

    // event related properties
    AliVEvent*          fEvent;			      //!<! Event object
    AliMCEvent*         fMCEvent;         //!<! MC event
    Double_t            fCent;            //!<! measured V0M centrality
    Double_t            fMultMeas;        //!<! measured multiplicity
    Double_t            fMultTrue;        //!<! true multiplicity
    Double_t            fMultMeasScaled;        //!<! measured multiplicity adjusted to data
    Double_t            fMultTrueScaled;        //!<! true multiplicity adjusted to data

    Int_t                           fRunNumber;                 //!<! run n
    Int_t                           fEventNumberInFile;         //!<! event number in file
    UInt_t                          fTimeStamp;                 //!<! event time stamp

    // track related properties
    Double_t                        fPt;                        //!<! track pT
    Double_t                        fEta;                       //!<! track Eta
    Double_t                        fSigmaPt;                  //!<! sigma(pT)/pT

    Double_t                        fMCPt;                      //!<! mc pt
    Double_t                        fMCEta;                     //!<! mc eta
    Int_t                           fMCLabel;                   //!<! mc label
    Bool_t                          fIsParticleInAcceptance;     //!<! particle in acceptance
    Bool_t                          fMCIsPhysicalPrimary;       //!<! is physical primary?
    Bool_t                          fMCIsCharged;               //!<! is charged?
    Bool_t                          fMCIsChargedPrimary;        //!<! is charged primary?
    Bool_t                          fMCIsChargedSecondary;        //!<! is charged secondary?

    Double_t                        fMCParticleWeight;          //!<! scaling factor of particle to match data
    Double_t                        fMCSecScaleWeight;          //!<! scaling factor of secondary to match data
    Int_t                           fNRepetitions;               //!<! how often to repeat this particle to match data
    Bool_t                          fUseRandomSeed;              ///<  use a random seed or a deterministic one (default)

    //UE: enumerator for region kTowards, kAway, kTransverse
    // external setter to select region

    // Tracking functions
    void InitESDTrackCuts();
    Bool_t AcceptTrackQuality(AliVTrack* track);
    Double_t GetCentrality(AliVEvent* event);

    // Data driven correction related functions
    Double_t GetSecScalingFactor(AliMCParticle* particle);
    Double_t GetParticleWeight(AliMCParticle* particle);

    Bool_t InitEvent();
    Bool_t InitTrack(AliVTrack* track);
    Bool_t InitParticle(AliMCParticle* particle);

    void LoopMeas(Bool_t count = kFALSE);
    void LoopTrue(Bool_t count = kFALSE);

    void FillEventHistos();
    void FillMeasTrackHistos();
    void FillMeasParticleHistos();
    void FillTrueParticleHistos();

    void FillMeasScaledTrackHistos();
    void FillTrueScaledParticleHistos();
    void FillMeasScaledParticleHistos();

    Int_t GetNRepetitons(Double_t scalingFactor);
    UInt_t GetSeed();

    // Histogramming functions
    THnSparseF* CreateHistogram(const string& name, const vector<string>& axes);
    TH1D* CreateLogHistogram(const string& name);
    TArrayD* GetBinEdges(const string& axisName);
    inline void FillHisto(THnSparseF* histo, const array<Double_t, MAX_HISTO_DIM>& values);
    inline void FillLogHisto(TH1D* logHist, const string& entry);
    string GetAxisTitle(const string& axisName);
    void SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins);

    AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&); // not implemented
    AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliMultDepSpecAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
