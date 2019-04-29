/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

class TParticle;
class AliESDEvent;
class AliVEvent;
class AliESDtrackCuts;

#include "THn.h"
#include "THnSparse.h"
#include "TF1.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

#include <iostream>
using std::string;
using std::vector;
using std::array;


class AliMultDepSpecAnalysisTask : public AliAnalysisTaskSE {
  public:
    static constexpr Int_t MAX_HISTO_DIM = 6;
    static constexpr Double_t PRECISION = 1e-6;
    AliMultDepSpecAnalysisTask();
    AliMultDepSpecAnalysisTask(const char *name);
    virtual ~AliMultDepSpecAnalysisTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t*);


    void SetTriggerMask(UInt_t triggermask)  {fTriggerMask = triggermask;}
    UInt_t GetTriggerMask()  {return fTriggerMask;}

    void SetCutMode(Int_t cutMode){fCutMode = cutMode;}
    // Setters
    void SetUseMC(Bool_t useMC = kTRUE){fIsMC = useMC;}
    void SetUseESD(){fIsESD = kTRUE;}
    void SetUseAOD(){fIsESD = kFALSE;}

    // Binning
    void SetBinsPt(TArrayD* bins){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(*bins);}
    void SetBinsPt(Int_t nBins, Double_t* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}

    void SetBinsEta(TArrayD* bins){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(*bins);}
    void SetBinsEta(Int_t nBins, Double_t* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}

    void SetBinsMult(TArrayD* bins){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(*bins);}
    void SetBinsMult(Int_t nBins, Double_t* binEdges){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(nBins+1,binEdges);}

    void SetBinsCent(TArrayD* bins){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(*bins);}
    void SetBinsCent(Int_t nBins, Double_t* binEdges){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(nBins+1,binEdges);}

    void SetBinsZv(TArrayD* bins){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(*bins);}
    void SetBinsZv(Int_t nBins, Double_t* binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}

    void SetBinsPtReso(Int_t nBins, Double_t* binEdges){if(fBinsPtReso) delete fBinsPtReso; fBinsPtReso = new TArrayD(nBins+1,binEdges);}

    // Acceptance cuts
    void SetMinEta(Double_t minEta){fMinEta = minEta;}
    void SetMaxEta(Double_t maxEta){fMaxEta = maxEta;}
    void SetMinPt(Double_t minPt){fMinPt = minPt;}
    void SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}


    Bool_t AcceptKinematics(AliVParticle* particle);
    Bool_t AcceptTrackQuality(AliVTrack* track);


    Double_t GetCentrality(AliVEvent* event);

    Bool_t IsChargedPrimary(Int_t mcLabel);

    void InitESDTrackCuts();
    void InitdNdPtEventCuts();
    void SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins);

    AliEventCuts fEventCuts; /// Event cuts

  private:

    TList*              fOutputList;		//!<! Output list
    AliVEvent*          fEvent;			    //!<! Event object
    AliMCEvent*         fMCEvent;       //!<! MC event

    AliESDtrackCuts*    fESDtrackCuts;
    Int_t               fCutMode;
    Bool_t              fIsESD;			    ///< Flag for ESD usage
    Bool_t              fIsMC;			    ///< Flag for MC usage

    // Acceptance cuts for tracks
    UInt_t                fTriggerMask;   // trigger mask
    Double_t              fMinEta;			///< Minimum eta cut
    Double_t              fMaxEta;			///< Maximum eta cut
    Double_t              fMinPt;			  ///< Minimum pT cut
    Double_t              fMaxPt;			  ///< Maximum pT cut

    // Binning
    TArrayD*             fBinsMult;		///< Array of bins in multiplicity
    TArrayD*             fBinsCent;		///< Array of bins in centrality
    TArrayD*             fBinsPt;			///< Array of bins in pt
    TArrayD*             fBinsEta;		///< Array of bins in eta
    TArrayD*             fBinsZv;			///< Array of bins in Zv (Z-position of primary vtx)
    TArrayD*             fBinsPtReso;			   ///< Array of bins for relative pt resoulution


    // Output Histograms

    TH1F* fHistEventSelection;    //!<! Histogram for triggered events and events with vertex
    THnSparseF* fHistEvents;      //!<! Histogram of measured events
    THnSparseF* fHistTracks;
    THnSparseF* fHistRelPtReso;

    THnSparseF* fHistMCRelPtReso;
    THnSparseF* fHistMCMultCorrelMatrix;
    THnSparseF* fHistMCPtCorrelMatrix;
    THnSparseF* fHistMCEtaCorrelMatrix;
    THnSparseF* fHistMCPrimTrue;
    THnSparseF* fHistMCPrimMeas;
    THnSparseF* fHistMCSecMeas;


    THnSparseF* CreateHistogram(string name, vector<string> axes);
    TArrayD* GetBinEdges(string& axisName);
    inline void FillHisto(THnSparseF* histo, array<Double_t, MAX_HISTO_DIM> values);
    string GetAxisTitle(string& axisName);
    TH1D* CreateLogHistogram(const char* logHistName);
    void Log(TH1D* logHist, const char* stage);

    AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&); // not implemented
    AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliMultDepSpecAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
