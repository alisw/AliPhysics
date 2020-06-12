/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define MAX_ALLOWED_MULT_BINS 500
#define MAX_HISTO_DIM 4
#define PRECISION 1e-6

#include "TSystem.h"
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
#include "AliAODTrack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliVHeader.h"

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

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t*);
    virtual void   Terminate(Option_t*);

    // Setters
    void SetTriggerMask(unsigned int triggermask)  {fTriggerMask = triggermask;}
    void SetIsMC(bool isMC = true){fIsMC = isMC;}
    void SetIsAOD(bool isAOD = true){fIsESD = !isAOD;}
    void SetUseDataDrivenCorrections(bool useDDC = true){fMCUseDDC = useDDC;}
    void SetUseZDCCut(bool useZDC){fUseZDCCut = useZDC;}
    void SetOverridePbPbEventCuts(bool overridePbPbEventCuts){fOverridePbPbEventCuts = overridePbPbEventCuts;}

    // Binning
    void SetBinsPt(int nBins, double* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    void SetBinsEta(int nBins, double* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    void SetBinsMult(int nBins, double* binEdges){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(nBins+1,binEdges);}
    void SetBinsZv(int nBins, double* binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}
    void SetBinsPtReso(int nBins, double* binEdges){if(fBinsPtReso) delete fBinsPtReso; fBinsPtReso = new TArrayD(nBins+1,binEdges);}
    void SetBinsMult(vector<int> multSteps, vector<int> multBinWidth);
    void SetBinsMult(int maxMult);

    // Acceptance cuts
    void SetMinEta(double minEta)   {fMinEta = minEta;}
    void SetMaxEta(double maxEta)   {fMaxEta = maxEta;}
    void SetMinPt(double minPt)     {fMinPt = minPt;}
    void SetMaxPt(double maxPt)     {fMaxPt = maxPt;}
    void SetMaxZv(double maxZv)     {fMaxZv = maxZv;}
    void SetMaxCent(double maxCent) {fMaxCent = maxCent;}
    void SetMinCent(double minCent) {fMinCent = minCent;}

    //void SetMinV0Mult(double minV0Mult)  {fMinV0Mult = minV0Mult;}

    // Configure this object for a train run
    static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(string dataSet, TString options, int cutModeLow = 100, int cutModeHigh = 119, bool isMC = false);
    void SaveTrainMetadata();
    bool InitTask(bool isMC, bool isAOD, string dataSet, TString options, int cutMode = 100);
    bool SetupTask(string dataSet, TString options);


  private:
  
    TList*              fOutputList;		        //!<! Output list
    AliEventCuts        fEventCuts;             //!<! Event cuts
    AliESDtrackCuts*    fTrackCuts;             //-> Track cuts
    TRandom3*           fRand;                  //!<! Random generator

    string            fTrainMetadata;           ///<  metadata of the train run used to generate the output

    bool              fIsESD;			              ///< Flag for ESD usage
    bool              fIsMC;                    ///< Flag for MC usage
    bool              fUseZDCCut;               ///< Flag for zdc cut usage
    bool              fOverridePbPbEventCuts;   ///< override centrality cut in PbPb
    bool              fMCUseDDC;                ///< Flag for data driven corrections usage
    // Cuts
    unsigned int        fTriggerMask;   ///< Trigger mask
    double              fMinEta;        ///< Minimum eta cut
    double              fMaxEta;        ///< Maximum eta cut
    double              fMinPt;			    ///< Minimum pT cut
    double              fMaxPt;			    ///< Maximum pT cut
    double              fMaxZv;			    ///< Maximum absolute z vertex cut
    double              fMinCent;       ///< Minimum centrality
    double              fMaxCent;       ///< Maximum centrality

    // Binning
    TArrayD*             fBinsEventCuts;  ///< Array of bins for event cuts
    TArrayD*             fBinsMult;       ///< Array of bins in central barrel track multiplicity
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

    THnSparseF* fHistMCRelPtReso;         //!<! Histogram of relative pt resolution from mc
    THnSparseF* fHistMCMultCorrelMatrix;  //!<! Histogram of multilicity correlation
    THnSparseF* fHistMCPtCorrelMatrix;    //!<! Histogram of pT correlation
    THnSparseF* fHistMCEtaCorrelMatrix;   //!<! Histogram of eta correlation
    THnSparseF* fHistMCPrimTrue;          //!<! Histogram of generated primaries
    THnSparseF* fHistMCPrimMeas;          //!<! Histogram of measured primaries
    THnSparseF* fHistMCSecMeas;           //!<! Histogram of measured secondaries

    // event related properties
    AliVEvent*          fEvent;			      //!<! Event object
    AliMCEvent*         fMCEvent;         //!<! MC event
    double              fMultMeas;        //!<! measured central barrel track multiplicity
    double              fMultTrue;        //!<! true multiplicity

    int                           fRunNumber;                 //!<! run number
    unsigned long                 fEventNumber;               //!<! event number
    unsigned int                  fTimeStamp;                 //!<! event time stamp

    // track related properties
    double                        fPt;                         //!<! track pT
    double                        fEta;                        //!<! track Eta
    double                        fSigmaPt;                    //!<! sigma(pT)/pT

    double                        fMCPt;                       //!<! mc pt
    double                        fMCEta;                      //!<! mc eta
    int                           fMCLabel;                    //!<! mc label
    bool                          fIsParticleInAcceptance;     //!<! particle in acceptance

    bool                          fMCIsChargedPrimary;          //!<! is charged primary?
    bool                          fMCIsChargedSecDecay;         //!<! is charged secondary from decay?
    bool                          fMCIsChargedSecMat;           //!<! is charged secondary from material?
    bool                          fMCIsChargedSecondary;        //!<! is charged secondary?

    double                        fMCParticleWeight;            //!<! scaling factor of particle to match data
    double                        fMCSecScaleWeight;            //!<! scaling factor of secondary to match data
    int                           fNRepetitions;                //!<! how often to repeat this particle to match data
    bool                          fUseRandomSeed;               ///<  use a random seed or a deterministic one (default)

    // Tracking functions
    void InitTrackCuts();
    bool AcceptTrackQuality(AliVTrack* track);
    double GetCentrality(AliVEvent* event);

    // Data driven correction related functions
    double GetSecScalingFactor(AliVParticle* particle);
    double GetParticleWeight(AliVParticle* particle);

    bool InitEvent();
    bool InitTrack(AliVTrack* track);

    template<typename Particle_t>
    bool InitParticle(Particle_t* particle);

    void LoopMeas(bool count = false);
    void LoopTrue(bool count = false);

    void FillEventHistos();
    void FillMeasTrackHistos();
    void FillMeasParticleHistos();
    void FillTrueParticleHistos();

    int GetNRepetitons(double scalingFactor);
    unsigned long GetSeed();

    // Histogramming functions
    THnSparseF* CreateHistogram(const string& name, const vector<string>& axes);
    TH1D* CreateLogHistogram(const string& name);
    TArrayD* GetBinEdges(const string& axisName);
    inline void FillHisto(THnSparseF* histo, const array<double, MAX_HISTO_DIM>& values);
    inline void FillLogHisto(TH1D* logHist, const string& entry);
    string GetAxisTitle(const string& axisName);
    void SetFixedBinEdges(double* array, double lowerEdge, double upperEdge, int nBins);

    AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&); // not implemented
    AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliMultDepSpecAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
