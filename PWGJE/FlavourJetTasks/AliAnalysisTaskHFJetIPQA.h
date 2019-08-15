#ifndef ALIANALYSISTASKJETIPQA_H
#define ALIANALYSISTASKJETIPQA_H
#include "AliAnalysisTaskEmcalJet.h"
#include "TGraph.h"
#include "THistManager.h"
#include <memory>
class AliEmcalJet;
class AliAODVertex;
class AliAODTrack;
class TList;
class TH1D;
class THistManager;
class AliPIDResponse;
class AliHFEpidBayes;
class TTree;
class TH2D;
class TCanvas;
class TParticle;
class TClonesArray;
class AliAODMCParticle;
class AliMCEvent;
class AliESDEvent;
class AliESDtrack;
class TGraph;
class AliAnalysisUtils;
class TRandom3;
class AliTriggerAnalysis;
class THnSparse;
class AliOADBContainer;
class AliEmcalList;
class AliVertexerTracks;
class AliPID;
class TGraph;
class AliVParticle;
class AliPIDCombined;
class AliEmcalJetFinder;
#include "AliFJWrapper.h"


#include "TMatrixD.h"
#include "TF1.h"
#include "AliESDtrackCuts.h"
#include <vector>
#include <utility>
#include <map>

//Define helper classes to cleanup analysis

class AliAnalysisTaskHFJetIPQA: public AliAnalysisTaskEmcalJet
{
public:
    //STATIC ENUM DEFINITIONS
    enum EPileup {kNoPileupSelection,kRejectPileupEvent,kRejectTracksFromPileupVertex};
    enum ERejBits {kNotSelTrigger,kNoVertex,kTooFewVtxContrib,kVertexChi2NDF,kZVtxOutFid,kPileupSPD,kOutsideCentrality,kVertexZContrib,kPhysicsSelection,kNoContributors,kDeltaVertexZ,kNoVertexTracks,kVertexZResolution,kMVPileup,kSPDClusterCut,kZVtxSPDOutFid,
                   kBadDiamondXDistance,kBadDiamondYDistance,kBadDiamondZDistance};
    enum TTypeImpPar {kXY,kXYSig,kXYZ,kXYZSig,kXYZSigmaOnly,kZSig};
    enum EParticleType  {bPi0=111,bEta=221,bEtaPrime=331,bPhi=333,bRho=113,bOmega=223,bSigma0=3212,bK0s=310,bLambda=3122,bPi=211,bProton=2212,bKaon=321,bOmegaBaryon=3334,
                         bAntiOmegaBaryon=-3334,bXiBaryon=3312,bAntiXiBaryon=-3312,bD0=421,bDPlus=411,bDStarPlus=413,bDSPlus=431,bK0l=130,bSigmaMinus = 3112,bSigmaPlus = 3222,bRhoPlus=213,
                         bBPlus = 521,bB0 = 511,bLambdaB =5122,bLambdaC=4122,bBStarPlus=523,bK0S892 = 313,bK0S892plus = 323};
    enum EParticleArrayIdx
    {bIdxPi0=0,bIdxEta=1,bIdxEtaPrime=2,bIdxPhi=3,bIdxRho=4,bIdxOmega=5,bIdxK0s=6,bIdxLambda=7,bIdxPi=8,bIdxProton=9,bIdxKaon=10,bIdxD0=11,bIdxDPlus=12,
        bIdxDStarPlus=13,bIdxDSPlus=14,bIdxLambdaC=15,bIdxBPlus = 16,bIdxB0 = 17,bIdxLambdaB = 18,bIdxBStarPlus=19,bIdxK0S892=20,bIdxK0S892plus=21,bIdxOmegaBaryon=22,bIdxXiBaryon=23,bIdxK0l=24,bIdxRest=25,bIdxSigmaPlus=26,bIdxSigmaMinus=27};

    enum bCuts{
        bAnalysisCut_SigmaDiamond=0,
        bAnalysisCut_Sigma_Z=1,
        bAnalysisCut_Sigma_Y=2,
        bAnalysisCut_RelError_Z=3,
        bAnalysisCut_RelError_Y=4,
        bAnalysisCut_NContibutors=5,
        bAnalysisCut_MaxVtxZ=6,
        bAnalysisCut_DCAJetTrack=7,
        bAnalysisCut_MaxDCA_Z=8,
        bAnalysisCut_MaxDCA_XY=9,
        bAnalysisCut_MaxDecayLength=10,
        bAnalysisCut_Z_Chi2perNDF=11,
        bAnalysisCut_MinTrackPt=12,
        bAnalysisCut_MinTrackPtMC=13,
        bAnalysisCut_MinTPCClus=14,
        bAnalysisCut_MinJetPt=15,
        bAnalysisCut_MaxJetPt=16,
        bAnalysisCut_MinJetEta=17,
        bAnalysisCut_MaxJetEta=18,
    };

    //UTILITY STRUCT DEFINITIONS
    struct SJetIpPati {
        SJetIpPati(Double_t v1, Double_t v2, Bool_t b, Bool_t c,Int_t tl,Double_t pt): first(v1),second(v2),is_electron(b),is_fromB(c),trackLabel(tl),trackpt(pt){}
        Double_t first; // to be compatible with std::pair
        Double_t second;// to be compatible with std::pair
        Bool_t   is_electron; // added for electron contribution check
        Bool_t   is_fromB; // added for electron contribution check
        Int_t trackLabel=-1 ;
        Double_t trackpt=-99;
    };
    struct SQuarks {
        SQuarks(Double_t pt, int label, int pdg, double deta): first(pt),second(label),HasPDG(pdg),HasDETA(deta){}
        Double_t first; // to be compatible with std::pair
        int second;// to be compatible with std::pair
        int HasPDG; // added for electron contribution check
        double HasDETA;
    };
    //FUNCTION DEFINITIONS
    AliAnalysisTaskHFJetIPQA();
    AliAnalysisTaskHFJetIPQA(const char *name);
    AliAnalysisTaskHFJetIPQA(const AliAnalysisTaskHFJetIPQA&); // not implemented
    AliAnalysisTaskHFJetIPQA& operator=(const AliAnalysisTaskHFJetIPQA&); // not implemented
    virtual ~AliAnalysisTaskHFJetIPQA(){;}
    virtual void   UserCreateOutputObjects();
    virtual void   UserExecOnce();
    virtual void   Terminate(Option_t *option="");
    virtual Bool_t Run();
    virtual Bool_t IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits);



    void SetESDCuts (AliESDtrackCuts  *cuts =NULL){fESDTrackCut =  new AliESDtrackCuts(*cuts);}
    void SetUseMonteCarloWeighingLinus(TH1F *Pi0 ,TH1F *Eta,TH1F *EtaP,TH1F *Rho,TH1F *Phi,TH1F *Omega,TH1F *K0s,TH1F *Lambda,TH1F *ChargedPi,
                                       TH1F *ChargedKaon,TH1F *Proton,TH1F *D0,TH1F *DPlus,TH1F *DStarPlus,
                                       TH1F *DSPlus,TH1F *LambdaC,TH1F *BPlus,TH1F *B0,TH1F *LambdaB,TH1F *BStarPlus);
    void SetFlukaFactor(TGraph* GraphOmega, TGraph* GraphXi, TGraph* K0Star, TGraph* Phi);
    void localtoglobal(double alpha, double *local, double *global);
    Bool_t FillTrackHistograms(AliVTrack * track, double * dca , double *cov,double weight);
   // void EventwiseCleanup();
    AliVParticle * GetVParticleMother(AliVParticle *part);
    Bool_t IsPhysicalPrimary(AliVParticle *part);
    void SetDefaultAnalysisCuts();
    void ChangeDefaultCutTo(AliAnalysisTaskHFJetIPQA::bCuts cutname, Double_t newcutvalue);
    Double_t GetLocalAlphaAOD(AliAODTrack *track);
    Double_t GetTrackCurvature(AliAODTrack *track);
    AliAODVertex *RemoveDaughtersFromPrimaryVtx(const AliVTrack * const track);
    void GetMaxImpactParameterCutR(const AliVTrack * const track, Double_t &maximpactRcut);
    Double_t GetLocalThetaAOD(AliAODTrack *track);
    Bool_t IsVertexSelected(const AliVVertex *vertex);
    Bool_t GetImpactParameter(const AliAODTrack *track, const AliAODEvent *event, Double_t *dca, Double_t *cov, Double_t *XYZatDCA);
    AliExternalTrackParam GetExternalParamFromJet(const AliEmcalJet *jet, const AliAODEvent *event);
    Bool_t GetImpactParameterWrtToJet(const AliAODTrack *track, const AliAODEvent *event, const AliEmcalJet *jet, Double_t *dca, Double_t *cov, Double_t *XYZatDCA, Double_t &jetsign);
    Bool_t getJetVtxMass( AliEmcalJet *jet, double &value);
    void SetJetRadius(Double_t fJetRadRead){fJetRadius=fJetRadRead;}
    int GetMCTruth(AliAODTrack *track, int &motherpdg);
    bool GetPIDCombined(AliAODTrack * track, double *prob, int &nDetectors, UInt_t &usedDet , AliPID::EParticleType &MostProbablePID, bool setTrackPID );
    void setFProductionNumberPtHard(Int_t value=-1)
    {
        fProductionNumberPtHard = value;
    }
    void RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont);   //Based on AliAnalysisTaskEmcalQGTagging::RecursiveParents
    void StoreDaughters(AliAODMCParticle* part, int kFirstMotherLabel);
    int DoJetPartonMatching(const AliEmcalJet *jet, int partonlabel);
    void DoFlavourVectorFill(int kJetOrigin, double partonpt, double jetpt, int kPartonsInJet, double deta);
    int DoFlavourDecision();
    Bool_t IsParton(int pdg);
    Double_t CalculateJetProb(AliEmcalJet *jet);
    Double_t CalculatePSTrack(Double_t sign, Double_t significance, Double_t trackPt, Int_t trclass);
    Double_t CalculatePSTrackPID(Double_t sign, Double_t significance, Double_t trackPt, Int_t trclass, Int_t species);

    //Setter for turning Corrections off and on
    void SmearTrack(AliAODTrack *track);
    void setFRunSmearing(Bool_t value){fRunSmearing = value;}
    void setFDoMCCorrection(Bool_t value){fDoMCCorrection=value;}
    void setFDoUnderlyingEventSub(Bool_t value){fDoUnderlyingEventSub=value;}
    void setfDoFlavourMatching(Bool_t value){fDoFlavourMatching=value;}
    void setfDaughterRadius(Double_t value){fDaughtersRadius=value;}
    void setfNoJetConstituents(Int_t value){fNoJetConstituents=value;}
    void setIsPythia(int i=0){if(i==6)fPythia6=kTRUE;if(i==8)fPythia8=kTRUE;}



    Bool_t IsTrackAcceptedJP(AliVTrack *track, Int_t n);
    bool IsFromElectron(AliAODTrack *track);
    bool IsFromProton(AliAODTrack *track);
    bool IsFromKaon(AliAODTrack *track);
    bool IsFromPion(AliAODTrack *track);
    Bool_t SetResFunctionPID(const char * filename);
    Double_t getFMCglobalDCAxyShift() const;
    void setFMCglobalDCAxyShift(const Double_t &value);
    Double_t getFVertexRecalcMinPt() const;
    void setFVertexRecalcMinPt(const Double_t &value);
    void setFMCglobalDCASmear(const Double_t value);
    void setFParam_Smear_Sigma(Double_t value)
    {
        fParam_Smear_Sigma = value;
    }
    void setFParam_Smear_Mean(Double_t value)
    {
        fParam_Smear_Mean = value;
    }

    void setGlobalVertex(Bool_t value)
    {
        fGlobalVertex = value;
    }

   void setDoNotCheckIsPhysicalPrimary(Bool_t value)
    {
        fDoNotCheckIsPhysicalPrimary = value;
    }
     void setDoJetProb(Bool_t value)
    {
        fDoJetProb = value;
    }
    void useTreeForCorrelations(Bool_t value){
        fUseTreeForCorrelations = value;
    }
    //virtual Bool_t IsEventSelected();
    void FillCorrelations(bool bn[3], double v[3], double jetpt);
    void setFFillCorrelations(const Bool_t &value);
    virtual void SetPtHardBin(Int_t b){ fSelectPtHardBin = b;}
    void SetHardCutoff(Double_t t)                            {fHardCutOff = t;}



protected:
    TH1D *fh1dTracksAccepeted; //!
    TH1D *fh1dCuts; //!
    TH2D *fh2dManifoldParton; //!
    TH2D *fh2dLightNotContrib; //!
    TH2D *fh2dCharmNotContrib; //!
    TH2D *fh2dBottomNotContrib; //!
    TH2D *fh2dLightNMatch; //!
    TH2D *fh2dCharmNMatch; //!
    TH2D *fh2dBottomNMatch; //!
    TH2D *fh2dLightDeta; //!
    TH2D *fh2dCharmDeta; //!
    TH2D *fh2dBottomDeta; //!

    THnSparse *fHLundIterative;//!       iterative declustering


private:
    THistManager         fHistManager    ;///< Histogram manager
    const AliAODVertex * fEventVertex;//!
    AliPIDResponse *fPidResponse ;//!
    AliEmcalJet *  GetPerpendicularPseudoJet (AliEmcalJet*jet_in  , bool rev );
    void GetOutOfJetParticleComposition(AliEmcalJet * jet, int flavour);
    void FillParticleCompositionSpectra(AliEmcalJet * jet,const char * histname );
    void FillParticleCompositionEvent();
    void DoJetLoop(); //jet matching function 2/4
    void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Int_t matching=0);
    void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
    void SmearTrackHybrid(AliVTrack * track);
    void FillHist(const char * name,Double_t x ,Double_t w);
    void FillHist(const char * name,Double_t x, Double_t y,Double_t w);
    void IncHist(const char * name,Int_t bin);
    void SubtractMean (Double_t val[2],AliVTrack *track);
    Bool_t IsTrackAccepted(AliVTrack* track,Int_t n=6);
    Bool_t MatchJetsGeometricDefault(); //jet matching function 1/4
    Double_t GetMonteCarloCorrectionFactor(AliVTrack *track, Int_t &pCorr_indx, double &ppt);
    Double_t GetWeightFactor( AliVTrack * mcpart,Int_t &pCorr_indx, double &ppt);
    Bool_t ParticleIsPossibleSource(Int_t pdg);
    Bool_t IsSelectionParticle( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
    Bool_t IsSelectionParticleALICE( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
    Bool_t IsSelectionParticleStrange( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
    Bool_t IsSelectionParticleMeson( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
    Bool_t IsSelectionParticleOmegaXiSigmaP( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx );
    Bool_t IsSecondaryFromWeakDecay( AliVParticle * particle ) ;
    Bool_t IsTruePrimary	(AliVParticle * mcpart);
    Bool_t GetBMesonWeight( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
    Bool_t IsPromptDMeson(AliVParticle * part );
    Bool_t IsPromptBMeson(AliVParticle * part );
    Double_t GetValImpactParameter(TTypeImpPar type, Double_t *impar, Double_t *cov);
    static Bool_t mysort(const SJetIpPati& i, const SJetIpPati& j);
    static Bool_t myquarksort(const SQuarks& i, const SQuarks& j);
    Int_t IsMCJetPartonFast(const AliEmcalJet *jet,  Double_t radius,Bool_t &is_udg);
    Int_t GetRunNr(AliVEvent * event){return event->GetRunNumber();}
    Double_t GetPtCorrected(const AliEmcalJet* jet);
    Double_t GetPtCorrectedMC(const AliEmcalJet *jet);
    void PrintSettings();


    //Functions to allow jet probability/TC System 8 efficiency estimation
    Bool_t IsJetTaggedTC(int n =0 ,double thres = 0.1);
    Bool_t IsJetTaggedJetProb(double thresProb = 0.90);
    TH1 *  AddHistogramm(const char * name,const char * title,Int_t x,Double_t xlow,Double_t xhigh, Int_t y=0,Double_t ylow=0,Double_t yhigh=0);
    TH1D * GetHist1D(const char * name){return (TH1D*)fOutput->FindObject(name);}
    TH2D * GetHist2D(const char * name){return (TH2D*)fOutput->FindObject(name);}


private:
    //Booleans for different corrections
    Bool_t   fRunSmearing;//
    Bool_t   fUsePIDJetProb;//
    Bool_t   fDoMCCorrection;//  Bool to turn on/off MC correction. Take care: some histograms may still be influenced by weighting.
    Bool_t   fDoUnderlyingEventSub;//
    Bool_t   fDoFlavourMatching;//
    Bool_t   fPythia6;//
    Bool_t   fPythia8;//

    Bool_t   fFillCorrelations;//
    Double_t fParam_Smear_Sigma;//
    Double_t fParam_Smear_Mean;//
    Bool_t   fGlobalVertex;//
    Bool_t fDoNotCheckIsPhysicalPrimary;//
    Bool_t fDoJetProb;
    TGraph * fGraphMean;//!
    TGraph * fGraphSigmaData;//!
    TGraph * fGraphSigmaMC;//!
    TGraph * fGraphXi; //!
    TGraph * fGraphOmega;
    TGraph * fK0Star; //!
    TGraph * fPhi; //!
    TGraph * fGeant3FlukaProton;//!
    TGraph * fGeant3FlukaAntiProton;//!
    TGraph * fGeant3FlukaLambda;//!
    TGraph * fGeant3FlukaAntiLambda;//!
    TGraph * fGeant3FlukaKMinus;//!
    //! \brief cCuts
    TCanvas *cCuts; //
    //! \brief fMCArray
    TClonesArray     *fMCArray;//!
    AliMCEvent       *fMCEvent;//!
    AliESDtrackCuts  *fESDTrackCut;//
    AliVertexerTracks *fVertexer;//!
    Bool_t fMcEvtSampled;//
    Double_t fBackgroundFactorLinus[21][498]; //[21][498]FineBinned correction factors up 0.1-25 GeV/c first value below last above 0.05 binwidth
    std::vector <Double_t > fPUdsgJet;//!
    std::vector <Double_t > fPSJet;//!
    std::vector <Double_t > fPCJet;//!
    std::vector <Double_t > fPBJet;//!
    std::vector <Double_t > fJetCont;//!
    std::vector <SQuarks > fHardProcess;//!
    std::vector <SQuarks > fQuarkVec;//!
    std::map<int, int> daughtermother;//!

    TGraph fResolutionFunction[200];//[200]<-
    Double_t fAnalysisCuts[19]; // /Additional (to ESD track cut or AOD filter bits) analysis cuts.
    AliPIDCombined *fCombined ;//!
    Float_t fXsectionWeightingFactor;//
    Int_t   fProductionNumberPtHard;//
    Double_t fJetRadius;//
    Double_t fDaughtersRadius;//
    Int_t fNoJetConstituents;//
    Double_t fMCglobalDCAxyShift;//
    Double_t fMCglobalDCASmear;//
    Double_t fVertexRecalcMinPt;//
    Double_t fHardCutOff;//
//Event mixing for correlation study
    Double_t fn1_mix ;
    Double_t fn2_mix ;
    Double_t fn3_mix ;
    Bool_t   fIsMixSignalReady_n1;
    Bool_t   fIsMixSignalReady_n2;
    Bool_t   fIsMixSignalReady_n3;
    Bool_t   fIsSameEvent_n1;
    Bool_t   fIsSameEvent_n2;
    Bool_t   fIsSameEvent_n3;
    Bool_t   fUseTreeForCorrelations;
    TTree *  fCorrelationCrossCheck;//!
    Float_t  fTREE_n1;
    Float_t  fTREE_n2;
    Float_t  fTREE_n3;
    Float_t  fTREE_pt;


    void SetMixDCA(int n , Double_t v){
    if(n==1){
            if(fIsMixSignalReady_n1) return;
            fn1_mix = v;
            fIsMixSignalReady_n1 = kTRUE;
            fIsSameEvent_n1 = kTRUE;

        }
    else if(n==2){
            if(fIsMixSignalReady_n2) return;
            fn2_mix = v;
            fIsMixSignalReady_n2 = kTRUE;
            fIsSameEvent_n2 = kTRUE;

        }
    else if(n==3){
            if(fIsMixSignalReady_n3) return;
            fn3_mix = v;
            fIsMixSignalReady_n3 = kTRUE;
            fIsSameEvent_n3 = kTRUE;
        }
    }

    Bool_t GetMixDCA(int n , double &v){
        if(n==1){
                if (!fIsMixSignalReady_n1 || fIsSameEvent_n1) return kFALSE;
                v= fn1_mix;
                fIsMixSignalReady_n1 = kFALSE;
            }
        else if(n==2){
                if (!fIsMixSignalReady_n2|| fIsSameEvent_n2) return kFALSE;
                v = fn2_mix;
                fIsMixSignalReady_n2 = kFALSE;
            }
        else if(n==3){
                if (!fIsMixSignalReady_n3|| fIsSameEvent_n3) return kFALSE;
                v = fn3_mix;
                fIsMixSignalReady_n3 = kFALSE;
            }
    return kTRUE;
    }







    ClassDef(AliAnalysisTaskHFJetIPQA, 36)
};

#endif


