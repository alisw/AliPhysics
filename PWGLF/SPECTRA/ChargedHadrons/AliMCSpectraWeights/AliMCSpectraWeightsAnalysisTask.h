/*!
   \file AliMCSpectraWeightsAnalysisTask.cxx
   \brief A minimal analysis task for AliMCSpectraWeights
   \author Patrick Huhn
   \date 25/10/2019
*/
#ifndef ALIMCSPECTRAWEIGHTSANALYSISTASK_H
#define ALIMCSPECTRAWEIGHTSANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#include "THn.h"
#include "TH3F.h"
class TArrayD;
class TString;
class AliMCSpectraWeights;
class AliStack;
class AliMCEvent;
class AliVEvent;
class TList;


/**
 * \class AliMCSpectraWeightsAnalysisTask
 * \brief Minimal analysis task for AliMCSpectraWeights
 */
class AliMCSpectraWeightsAnalysisTask : public AliAnalysisTaskSE {
  public:
    AliMCSpectraWeightsAnalysisTask();
    AliMCSpectraWeightsAnalysisTask(const char* name);
    virtual ~AliMCSpectraWeightsAnalysisTask();

    #ifdef __CLING__
      //C++ 11 feature
    AliMCSpectraWeightsAnalysisTask(const AliMCSpectraWeightsAnalysisTask&) = delete;
    AliMCSpectraWeightsAnalysisTask& operator=(const AliMCSpectraWeightsAnalysisTask&) = delete;
    #endif

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option) {}

    //------ Setter ---------
    void SetDebugLevel(Int_t level) { fDebugLevel = level; }
    // void SetCollisionSystem(const char* sys) {fstCollisionSystem=sys;}
    // void SetMCTrainOutputPath(const char* path) {fstMCTrainOutput=path;}
    void SetMCSpectraWeightObject(AliMCSpectraWeights *obj) {fMCSpectraWeights=obj;}
    void SetTriggerMask(UInt_t triggermask)  { fTriggerMask = triggermask; }
    /// Set the flag for the use of MC.
    void SetUseMC(Bool_t useMC=kTRUE){fIsMC=useMC;}
    /// Set the flag for the use of ESD \c (fIsESD=kTRUE)
    void SetUseESD(){fIsESD=kTRUE;}
    /// Set the flag for the use of AOD \c (fIsESD=kFALSE)
    void SetUseAOD(){fIsESD=kFALSE;}
    /// Set the flag for the particle composition
    void SetUseMultiplicity(Bool_t useMult){fUseMultiplicity = useMult;}
    void SetUseCentrality(){fUseMultiplicity = kFALSE;}
    // Binning
    /// Set bins in Pt using a TArrayD
    void SetBinsPt(TArrayD *bins){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(*bins);}
    /// Set bins in Pt using number of bins and array of bin edges
    void SetBinsPt(Int_t nBins, Double_t *binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Eta using a TArrayD
    void SetBinsEta(TArrayD *bins){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(*bins);}
    /// Set bins in Eta using number of bins and array of bin edges
    void SetBinsEta(Int_t nBins, Double_t *binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Multiplicity/Centrality using a TArrayD
    void SetBinsMultCent(TArrayD *bins){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(*bins);}
    /// Set bins in Multiplicity/Centrality number of bins and array of bin edges
    void SetBinsMultCent(Int_t nBins, Double_t *binEdges){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Zv using a TArrayD
    void SetBinsZv(TArrayD *bins){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(*bins);}
    /// Set bins in Zv using number of bins and array of bin edges
    void SetBinsZv(Int_t nBins, Double_t *binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}

    //------ Gettter --------
    UInt_t GetTriggerMask()  { return fTriggerMask; }
    TH3F* GetHistMCPartCorr() const {return fHistMCPartCorr;}
    THnF* GetHistMCGenPrimTrack() const {return fHistMCGenPrimTrack;}
    // Binning
    TArrayD* GetBinsPt() {return fBinsPt;}
    TArrayD* GetBinsEta(){return fBinsEta;}
    TArrayD* GetBinsMultCent() {return fBinsMultCent;}
    TArrayD* GetBinsZv() {return fBinsZv;}
    Double_t GetEventMultCent(AliVEvent *event);

  private:
    Int_t         fDebugLevel;           ///!< Debug level
    TList         *fOutputList;		//!<! Output list
    AliVEvent     *fEvent;			//!<! Event object (AliVEvent)
    AliMCEvent    *fMCEvent;		//!<! MC event
    AliStack      *fMCStack;		//!<! MC stack
    Bool_t        fIsESD;//
    Bool_t        fIsMC;//
    Bool_t        fUseMultiplicity;//
    UInt_t        fTriggerMask;    // trigger mask

    //Particle composition
    AliMCSpectraWeights *fMCSpectraWeights;//->
    TH3F                 *fHistMCPartCorr;//!
    THnF                 *fHistMCGenPrimTrack;//!
    TH3F                 *fHistMCFractions;//!
    TH3F                 *fHistDataFractions;//!
    TH3F                 *fHistMCWeights;//!

    //binning
    TArrayD     	*fBinsMultCent;		///< Array of bins in multiplicity or centrality
    TArrayD     	*fBinsPt;			///< Array of bins in pt
    TArrayD     	*fBinsEta;		///< Array of bins in eta
    TArrayD     	*fBinsZv;			///< Array of bins in Zv (Z-position of primary vtx)

    #if defined (__CINT__)
      AliMCSpectraWeightsAnalysisTask(const AliMCSpectraWeightsAnalysisTask&);
      AliMCSpectraWeightsAnalysisTask& operator=(const AliMCSpectraWeightsAnalysisTask&);
    #endif

    ClassDef(AliMCSpectraWeightsAnalysisTask, 1);
};

#endif /* ALIMCSPECTRAWEIGHTSANALYSISTASK_H */
