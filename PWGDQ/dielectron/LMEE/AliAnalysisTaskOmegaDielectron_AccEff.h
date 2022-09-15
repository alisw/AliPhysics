/// \classAliAnalysisTaskOmegaDielectron_AccEfff
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliAnalysisTaskOmegaDielectron_AccEff_cxx // header guard in case of multiple includes
#define AliAnalysisTaskOmegaDielectron_AccEff_cxx

#define MAX_HISTO_DIM 4
#define PRECISION 1e-6


#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "AliPID.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"


#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"

#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"

#include <iostream>
#include <vector>
using namespace std;
using std::string;
using std::vector;
using std::array;

class AliAnalysisTaskOmegaDielectron_AccEff : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskOmegaDielectron_AccEff();
    AliAnalysisTaskOmegaDielectron_AccEff(const char *name);
    virtual ~AliAnalysisTaskOmegaDielectron_AccEff();

    static AliAnalysisTaskOmegaDielectron_AccEff* AddTaskMultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 121);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t*);

    // Setters
    void SetTriggerMask(UInt_t triggermask)  {fTriggerMask = triggermask;}
    void SetCutMode(Int_t cutMode){fCutMode = cutMode;}
    void SetIsMC(Bool_t isMC = kTRUE){fIsMC = isMC;}
    void SetUseESD(){fIsESD = kTRUE;}
    void SetUseAOD(){fIsESD = kFALSE;}

    // Binning
    void SetBinsPt(Int_t nBins, Double_t* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    void SetBinsEta(Int_t nBins, Double_t* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    void SetBinsy(Int_t nBins, Double_t* binEdges){if(fBinsy) delete fBinsy; fBinsy = new TArrayD(nBins+1,binEdges);}
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

    TList*              fOutputList;		  //!<! Output list
    AliVEvent*          fEvent;			      //!<! Event object
    AliMCEvent*         fMCEvent;         //!<! MC event
    AliEventCuts        fEventCuts;       //!<! Event cuts
    AliESDtrackCuts*    fESDtrackCuts;    //!<! Track cuts

    Int_t               fCutMode;         ///< ID of track cut variation (100=default)
    Bool_t              fIsESD;			      ///< Flag for ESD usage
    Bool_t              fIsMC;            ///< Flag for MC usage
    Bool_t              fUseCent;         ///< Flag for Centrality usage

    // Acceptance cuts for tracks
    UInt_t                fTriggerMask;   ///< Trigger mask
    Double_t              fMinEta;        ///< Minimum eta cut
    Double_t              fMaxEta;        ///< Maximum eta cut
    Double_t              fMinPt;			    ///< Minimum pT cut
    Double_t              fMaxPt;			    ///< Maximum pT cut
    Double_t              fMaxZv;			    ///< Maximum absolute z vertex cut
    Double_t              fMinCent;       ///< Minimum centrality
    Double_t              fMaxCent;       ///< Maximum centrality

    // Binning
    TArrayD*             fBinsMult;       ///< Array of bins in multiplicity
    TArrayD*             fBinsCent;       ///< Array of bins in centrality
    TArrayD*             fBinsPt;			    ///< Array of bins in pt
    TArrayD*             fBinsEta;		    ///< Array of bins in eta
    TArrayD*             fBinsy;		    ///< Array of bins in y
    TArrayD*             fBinsZv;			    ///< Array of bins in Zv (Z-position of primary vtx)
    TArrayD*             fBinsPtReso;     ///< Array of bins for relative pt resoulution

    // pdg codes:
    Int_t                felectron_pdg;   //!
    Int_t                fpositron_pdg;   //!
    Int_t                fmother_pdg;     //!

    //storage vectors:
    vector<AliAODTrack *> v_elec_true_omega;    //! array of strings containing the electron track from a true omega dielectron decay
    vector<AliAODTrack *> v_posi_true_omega;    //! array of strings containing the positron track from a true omega dielectron decay
    vector<Int_t> v_elec_motherID_true_omega;    //! array of strings containing the electron track from a true omega dielectron decay
    vector<Int_t> v_posi_motherID_true_omega;    //! array of strings containing the positron track from a true omega dielectron decay



    // Output Histograms
    TH1F* fHistMC_Omegas_Rapidity;        //!<! Histogram for event counting
    TH1F* fHistEventSelection;            //!<! Histogram for event counting

    TH1F* fHistMC_ele1_posi2_OmegaDielDeacay;         //!<! Histogram of relative pt resolution from mc
    TH2F* fHistMC_Omegas_gen;          //!<! Histogram of generated primaries
    TH2F* fHistMC_Omegas_gen_DaughtersinAcc;          //!<! Histogram of generated primaries

    TH1F* fHist_rec_true_Ele_Omegas_Mothers;          //!<! Histogram of generated primaries
    TH1F* fHist_rec_true_Pos_Omegas_Mothers;          //!<! Histogram of generated primaries
    TH2D* fHist_rec_true_Dielec;          //!<! Histogram of generated primaries


    void InitESDTrackCuts();
    Bool_t AcceptKinematics(AliVParticle* particle);
    Bool_t EtaCut(AliVParticle* particle);
    Bool_t AcceptTrackQuality(AliVTrack* track);
    Double_t GetCentrality(AliVEvent* event);
    Bool_t IsChargedPrimary(Int_t mcLabel);
    Bool_t CheckDielectronDecay(AliMCParticle *particle);
    Bool_t CheckDielectronDecay_DaughterinAcc(AliMCParticle *particle);
    void SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins);

    THnSparseF* CreateHistogram(string name, vector<string> axes);
    TArrayD* GetBinEdges(string& axisName);
    inline void FillHisto(THnSparseF* histo, array<Double_t, MAX_HISTO_DIM> values);
    string GetAxisTitle(string& axisName);

    AliAnalysisTaskOmegaDielectron_AccEff(const AliAnalysisTaskOmegaDielectron_AccEff&); // not implemented
    AliAnalysisTaskOmegaDielectron_AccEff& operator=(const AliAnalysisTaskOmegaDielectron_AccEff&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskOmegaDielectron_AccEff, 1); // example of analysis
    /// \endcond
};

#endif
