#ifndef ALIANALYSISTASKPIDMCEffCont_cxx
#define ALIANALYSISTASKPIDMCEffCont_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance
// Function for single particles and pairs
//
// ---------------------------------------------------------------------

class TList;
class TH1F;
class TH3D;
class TH2F;
class TString;
class AliAODEvent;
class AliAODInputHandler;
class TH2D;

#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPIDMCEffCont : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskPIDMCEffCont() : AliAnalysisTaskSE(),
    fAOD(0),
    fArrayMC(0),
    fQAList(0),
    fOutputList(0),
    fHistEventStats(0),
    fHistCentrality(0),
    fHistNMult(0),
    fHistVz(0),
    fHistPt(0),
    fHistPhi(0),
    fHistEta(0),
    fHistNSigmaTPCvsPtbeforePID(0),
    fHistNSigmaTPCvsPtafterPID(0),
    fHistContaminationSecondariesPlus(0),
    fHistContaminationSecondariesMinus(0),
    fHistContaminationPrimariesPlus(0),
    fHistContaminationPrimariesMinus(0),
    fHistGeneratedEtaPtPhiPlus(0),
    fHistSurvivedEtaPtPhiPlus(0),
    fHistGeneratedEtaPtPhiMinus(0),
    fHistSurvivedEtaPtPhiMinus(0),
    fUseCentrality(kFALSE),
    fCentralityEstimator("V0M"),
    fCentralityPercentileMin(0.0),
    fCentralityPercentileMax(5.0),
    fInjectedSignals(kFALSE),
    fElectronRejection(kFALSE),
    fElectronOnlyRejection(kFALSE),
    fElectronRejectionNSigma(-1.),
    fElectronRejectionMinPt(0.),
    fElectronRejectionMaxPt(1000.),
    fVxMax(3.0),
    fVyMax(3.0),
    fVzMax(10.),
    fAODTrackCutBit(128),
    fMinNumberOfTPCClusters(80),
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0),
    fMaxDCAz(3.0),
    fPtMin(0.0),
    fPtMax(20.0),
    fEtaMin(-0.8),
    fEtaMax(0.8),
    fUsePIDNsigma(kFALSE),
    fUsePIDPDG(kFALSE),
    fDCAextended(kFALSE),
    fPIDResponse(0),
    fPIDNSigma(3),
    fParticleOfInterest(kPion),
    fHistDCAXYptprimminus(0),
    fHistDCAXYptprimplus(0),
    fHistDCAXYptsecminusweak(0),
    fHistDCAXYptsecplusweak(0),
    fHistDCAXYptsecminusmat(0),
    fHistDCAXYptsecplusmat(0),
    fHistDCAXYptchargedminus(0),
    fHistDCAXYptchargedplus(0),
    fHistDCAXYptprimminus_ext(0),
    fHistDCAXYptprimplus_ext(0),
    fHistDCAXYptsecminusweak_ext(0),
    fHistDCAXYptsecplusweak_ext(0),
    fHistDCAXYptsecminusmat_ext(0),
    fHistDCAXYptsecplusmat_ext(0),
    fHistDCAXYptchargedminus_ext(0),
    fHistDCAXYptchargedplus_ext(0){}
    AliAnalysisTaskPIDMCEffCont(const char *name);
    virtual ~AliAnalysisTaskPIDMCEffCont() {}
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    Bool_t   IsLabelUsed(TArrayI array, Int_t label);
    
    void SetAODtrackCutBit(Int_t bit){
        fAODTrackCutBit = bit;
    }
    void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
        fVxMax = vx;
        fVyMax = vy;
        fVzMax = vz;
    }
    
    //Centrality
    void UseCentrality() { fUseCentrality = kTRUE;}
    void SetCentralityEstimator(const char* centralityEstimator) {
        fCentralityEstimator = centralityEstimator;}
    void SetCentralityPercentileRange(Float_t min, Float_t max) {
        fCentralityPercentileMin=min;
        fCentralityPercentileMax=max;
    }
    
    //Injected signals
    void SetRejectInjectedSignals() {fInjectedSignals = kTRUE;}
    
    // electron rejection
    void SetElectronRejection(Double_t gMaxNSigma){
        fElectronRejection = kTRUE;
        fElectronRejectionNSigma = gMaxNSigma;
    }
    
    void SetElectronOnlyRejection(Double_t gMaxNSigma){
        fElectronRejection       = kTRUE;
        fElectronOnlyRejection   = kTRUE;
        fElectronRejectionNSigma = gMaxNSigma;
    }
    
    void SetElectronRejectionPt(Double_t minPt,Double_t maxPt){
        fElectronRejectionMinPt  = minPt;
        fElectronRejectionMaxPt  = maxPt;
    }
    
    //Track cuts
    void SetMinNumberOfTPCClusters(Double_t min) {
        fMinNumberOfTPCClusters = min;}
    void SetMaxChi2PerTPCCluster(Double_t max) {
        fMaxChi2PerTPCCluster = max;}
    void SetMaxDCAxy(Double_t max) {
        fMaxDCAxy = max;}
    void SetMaxDCAz(Double_t max) {
        fMaxDCAz = max;}

    void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax) {
        fEtaMin  = etamin; fEtaMax  = etamax;
        fPtMin  = ptmin; fPtMax  = ptmax;
        
    }
    
    void UsePIDNsigma(){fUsePIDNsigma = kTRUE;}
    
    void UsePIDPDG(){fUsePIDPDG = kTRUE;}
    
    void SetNSigmaPID(Int_t nsigma) {
        
        fPIDNSigma = nsigma;
    }
  
    void USEextendedDCA() {fDCAextended = kTRUE;}

    enum kParticleOfInterest { kMuon, kElectron, kPion, kKaon, kProton };
    
    void setParticleType(kParticleOfInterest ptype){
        fParticleOfInterest = ptype;
    }
    
private:
    AliAODEvent* fAOD; //! AOD object
    TClonesArray *fArrayMC; //! array of MC particles
    TList       *fQAList; //! QA list
    TList       *fOutputList; //! Output list
    
    // QA histograms
    TH1F        *fHistEventStats; //!event stats
    TH1F        *fHistCentrality; //!centrality
    TH1F        *fHistNMult; //! nmult
    TH2F        *fHistVz;//!
    TH2F        *fHistNSigmaTPCvsPtbeforePID;//TPC nsigma vs pT before PID cuts (QA histogram)
    TH2F        *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)
    
    TH1F 	*fHistPt;
    TH1F 	*fHistPhi;
    TH1F	*fHistEta;

    // output histograms
    TH3D        *fHistContaminationSecondariesPlus;//!
    TH3D        *fHistContaminationSecondariesMinus;//!
    TH3D        *fHistContaminationPrimariesPlus;//!
    TH3D        *fHistContaminationPrimariesMinus;//!
    
    // output histograms (single particles)
    TH3D        *fHistGeneratedEtaPtPhiPlus;//!correction map for positives (generated)
    TH3D        *fHistSurvivedEtaPtPhiPlus;//!correction map positives (survived)
    
    TH3D        *fHistGeneratedEtaPtPhiMinus;//!correction map for negatives (generated)
    TH3D        *fHistSurvivedEtaPtPhiMinus;//!correction map negatives (survived)
    
    TH3F        *fHistDCAXYptprimminus;
    TH3F        *fHistDCAXYptprimplus;
    TH3F        *fHistDCAXYptsecminusweak;
    TH3F        *fHistDCAXYptsecplusweak;
    TH3F        *fHistDCAXYptsecminusmat;
    TH3F        *fHistDCAXYptsecplusmat;
    TH3F        *fHistDCAXYptchargedminus;
    TH3F        *fHistDCAXYptchargedplus;
    TH3F        *fHistDCAXYptprimminus_ext;
    TH3F        *fHistDCAXYptprimplus_ext;
    TH3F        *fHistDCAXYptsecminusweak_ext;
    TH3F        *fHistDCAXYptsecplusweak_ext;
    TH3F        *fHistDCAXYptsecminusmat_ext;
    TH3F        *fHistDCAXYptsecplusmat_ext;
    TH3F        *fHistDCAXYptchargedminus_ext;
    TH3F        *fHistDCAXYptchargedplus_ext;

    Bool_t  fUseCentrality;// Bool_t use centrality or not
    TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
    Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile
    
    Bool_t fInjectedSignals;//Flag for using the rejection of injected signals
    Bool_t fUsePIDNsigma;
    Bool_t fUsePIDPDG;
    Bool_t fDCAextended;

    AliPIDResponse *fPIDResponse;     //! PID response object
    Bool_t   fElectronRejection;//flag to use electron rejection
    Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
    Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
    Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
    Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)
    
    Double_t fVxMax;//vxmax
    Double_t fVyMax;//vymax
    Double_t fVzMax;//vzmax
    
    Int_t fAODTrackCutBit;//track cut bit from track selection (only used for AODs)
    
    Double_t fMinNumberOfTPCClusters;//!
    Double_t fMaxChi2PerTPCCluster;//!
    Double_t fMaxDCAxy, fMaxDCAz;//!
    Double_t fPtMin, fPtMax;//!
    Double_t fEtaMin, fEtaMax;//!
    
    Int_t fPIDNSigma;
    kParticleOfInterest fParticleOfInterest;
    
    AliAnalysisTaskPIDMCEffCont(const AliAnalysisTaskPIDMCEffCont&); // not implemented
    AliAnalysisTaskPIDMCEffCont& operator=(const AliAnalysisTaskPIDMCEffCont&); // not implemented
    
    ClassDef(AliAnalysisTaskPIDMCEffCont, 1); // example of analysis
};

#endif
