/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define MAX_HISTO_DIM 4
#define PRECISION 1e-6


#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"

#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"

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
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t*);

    // Setters
    void SetTriggerMask(UInt_t triggermask)  {fTriggerMask = triggermask;}
    void SetCutMode(Int_t cutMode){fCutMode = cutMode;}
    void SetUseMC(Bool_t useMC = kTRUE){fIsMC = useMC;}
    void SetUseESD(){fIsESD = kTRUE;}
    void SetUseAOD(){fIsESD = kFALSE;}

    // Binning
    void SetBinsPt(Int_t nBins, Double_t* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    void SetBinsEta(Int_t nBins, Double_t* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    void SetBinsMult(Int_t nBins, Double_t* binEdges){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(nBins+1,binEdges);}
    void SetBinsCent(Int_t nBins, Double_t* binEdges){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(nBins+1,binEdges);}
    void SetBinsZv(Int_t nBins, Double_t* binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}
    void SetBinsPtReso(Int_t nBins, Double_t* binEdges){if(fBinsPtReso) delete fBinsPtReso; fBinsPtReso = new TArrayD(nBins+1,binEdges);}

    // Acceptance cuts
    void SetMinEta(Double_t minEta){fMinEta = minEta;}
    void SetMaxEta(Double_t maxEta){fMaxEta = maxEta;}
    void SetMinPt(Double_t minPt){fMinPt = minPt;}
    void SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}

  private:

    TList*              fOutputList;		  //!<! Output list
    AliVEvent*          fEvent;			      //!<! Event object
    AliMCEvent*         fMCEvent;         //!<! MC event
    AliEventCuts        fEventCuts;       //!<! Event cuts
    AliESDtrackCuts*    fESDtrackCuts;    //!<! Track cuts

    Int_t               fCutMode;         ///< ID of track cut variation (100=default)
    Bool_t              fIsESD;			      ///< Flag for ESD usage
    Bool_t              fIsMC;            ///< Flag for MC usage

    // Acceptance cuts for tracks
    UInt_t                fTriggerMask;   ///< Trigger mask
    Double_t              fMinEta;        ///< Minimum eta cut
    Double_t              fMaxEta;        ///< Maximum eta cut
    Double_t              fMinPt;			    ///< Minimum pT cut
    Double_t              fMaxPt;			    ///< Maximum pT cut

    // Binning
    TArrayD*             fBinsMult;       ///< Array of bins in multiplicity
    TArrayD*             fBinsCent;       ///< Array of bins in centrality
    TArrayD*             fBinsPt;			    ///< Array of bins in pt
    TArrayD*             fBinsEta;		    ///< Array of bins in eta
    TArrayD*             fBinsZv;			    ///< Array of bins in Zv (Z-position of primary vtx)
    TArrayD*             fBinsPtReso;     ///< Array of bins for relative pt resoulution


    // Output Histograms
    TH1F* fHistEventSelection;            //!<! Histogram for event counting
    THnSparseF* fHistEvents;              //!<! Histogram of measured event distribution
    THnSparseF* fHistTracks;              //!<! Histogram of measured tracks
    THnSparseF* fHistRelPtReso;           //!<! Histogram of relatvie pT resolution from covariance matrix

    THnSparseF* fHistMCRelPtReso;         //!<! Histogram of relative pt resolution from mc
    THnSparseF* fHistMCMultCorrelMatrix;  //!<! Histogram of multilicity correlation
    THnSparseF* fHistMCPtCorrelMatrix;    //!<! Histogram of pT correlation
    THnSparseF* fHistMCEtaCorrelMatrix;   //!<! Histogram of eta correlation
    THnSparseF* fHistMCPrimTrue;          //!<! Histogram of generated primaries
    THnSparseF* fHistMCPrimMeas;          //!<! Histogram of measured primaries
    THnSparseF* fHistMCSecMeas;           //!<! Histogram of measured secondaries


    void InitESDTrackCuts();
    Bool_t AcceptKinematics(AliVParticle* particle);
    Bool_t AcceptTrackQuality(AliVTrack* track);
    Double_t GetCentrality(AliVEvent* event);
    Bool_t IsChargedPrimary(Int_t mcLabel);
    void SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins);

    THnSparseF* CreateHistogram(string name, vector<string> axes);
    TArrayD* GetBinEdges(string& axisName);
    inline void FillHisto(THnSparseF* histo, array<Double_t, MAX_HISTO_DIM> values);
    string GetAxisTitle(string& axisName);

    AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&); // not implemented
    AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliMultDepSpecAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
