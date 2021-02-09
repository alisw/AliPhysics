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
#include <TDatabasePDG.h>
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
        bAnalysisCut_MinITSLayersHit=15,
        bAnalysisCut_MinTrackChi2=16,
        bAnalysisCut_MinJetPt=17,
        bAnalysisCut_MaxJetPt=18,
        bAnalysisCut_MinJetEta=19,
        bAnalysisCut_MaxJetEta=20,
        bAnalysisCut_MinJetArea=21,
        bAnalysisCut_HasSPD=22,
        bAnalysisCut_HasSDD=23,
        bAnalysisCut_HasSSD=24,
        bAnalysisCut_KinkCand=25,
        bAnalysisCut_HasTPCrefit=26,
        bAnalysisCut_HasITSrefit=27,
        bAnalysisCut_PtHardAndJetPtFactor=28,
        bAnalysisCut_MinNewVertexContrib=29,
        bAnalysisCut_SDz=30,
        bAnalysisCut_SDbeta=31,
        bAnalysisCut_MaxIPLNJP=32
    };

    enum V0Cuts{
        DaughMaxEta,
        DaughMinPt,
        MinDCADaughWrtPV,
        MaxDCADaughvsDaugh,
        IsTPCRefitOn,
        DoPosNoTPCClusters,
        MinNoCrossedTPCRows,
        NoCrossedOverNoTPCClustersMin,
        NoCrossedOverNoTPCClustersMax,
        IsKinkCand,

        MaxV0Eta,
        MaxV0Rap,
        MaxSigmadEdxTPC,
        MinDecayRadius,
        MaxDecayRadius,
        MaxCosPALambda,
        MinCosPAK0,
        MaxLifeTimeK0,
        MaxLifeTimeLambda,
        DoArmenteros,
        DoMassWindow,
        InvarMassWindowK0,
        InvarMassWindowLambda,

        fAV0Cut,
        fBV0Cut,
        fCV0Cut
    };

    enum TCTagType{
        TCNo,
        TCIPSigPtDep,
        TCIPFixedPt,
        TCIPSigPtDepVarNTemps
    };

    enum ProbTagType{
        ProbNo,
        ProbJP,
        ProblnJP
    };

    enum V0TagType{
        V0No,
        V0Rec,
        V0MC,
        V0TrueRec
    };

    enum V0Type{
        V0Untagged,
        V0K0s,
        V0Lambda,
        V0AntiLambda
    };

    enum TemplateFlavour{
        Unid,
        UDSG,
        C,
        B,
        UDSGV0,
        CV0
    };


    //UTILITY STRUCT DEFINITIONS
    struct SJetIpPati {
        SJetIpPati(Double_t v1, Double_t v2, Int_t isv0, Bool_t c,Int_t tl,Double_t pt, Int_t v0mcid): first(v1),second(v2),is_V0(isv0),is_fromB(c),trackLabel(tl),trackpt(pt),iv0MCID(v0mcid){}
        Double_t first; // to be compatible with std::pair
        Double_t second;// to be compatible with std::pair
        Int_t   is_V0; // added for electron contribution check
        Bool_t   is_fromB; // added for electron contribution check
        Int_t trackLabel;
        Double_t trackpt;
        Int_t iv0MCID;
    };

    struct SV0Daugh {
        SV0Daugh(): fPt(0), fEta(0), iCharge(0), iCrossedTPC(0), iNoTPCCluster(0), fDCAtoPV(0), bTPCRefitOn(kFALSE), bIsKink(kFALSE){}
          double fPt;
          double fEta;
          int iCharge;
          int iCrossedTPC;
          int iNoTPCCluster;
          double fDCAtoPV;
          bool bTPCRefitOn;
          bool bIsKink;

          void Reset() {memset(this,0, sizeof(*this));}
          void Print() const;
    };

    struct SV0Cand {
        SV0Cand():
        bOnFly(0),
        fDCAV0DaughvsDaugh(0),
        fPA(0),
        fDecayRadius(0),
        fLifetimeK0(0),
        fLifetimeLambda(0),
        fEta(0),
        fPt(0),
        fRapK0(0),
        fRapLambda(0),
        fDecayLength3D(0),
        fDecayLength2D(0),
        fArmenterosAlpha(0),
        fArmenterosPt(0),
        fMassK0(0),
        fMassLambda(0),
        fMassAntilambda(0),
        fSigmaPosPion(0),
        fSigmaPosProton(0),
        fSigmaNegPion(0),
        fSigmaNegProton(0),
        bDaughsMissing(0),

        bIsCandidateK0s (kTRUE), // candidate for K0s
        bIsCandidateLambda (kTRUE), // candidate for Lambda
        bIsCandidateALambda (kTRUE), // candidate for anti-Lambda
        bIsInPeakK0s (kFALSE), // candidate within the K0s mass peak
        bIsInPeakLambda (kFALSE), // candidate within the Lambda mass peak
        bIsInPeakALambda (kFALSE), // candidate within the anti-Lambda mass peak
        bIsInConeJet (kFALSE), // candidate within the jet cones
        bIsInConePerp (kFALSE), // candidate within a perpendicular cone
        bIsInConeRnd (kFALSE), // candidate within the random cone
        bIsInConeMed (kFALSE), // candidate within the median-cluster cone
         bIsOutsideCones (kFALSE) // candidate outside excluded cones
        {}

        bool bOnFly;
        double fDCAV0DaughvsDaugh;
        double fPA;
        double fDecayRadius;
        double fLifetimeK0;
        double fLifetimeLambda;
        double fEta;
        double fPt;
        double fRapK0;
        double fRapLambda;
        double fDecayLength3D;
        double fDecayLength2D;
        double fArmenterosAlpha;
        double fArmenterosPt;
        double fMassK0;
        double fMassLambda;
        double fMassAntilambda;
        double fSigmaPosPion;
        double fSigmaPosProton;
        double fSigmaNegPion;
        double fSigmaNegProton;
        bool bDaughsMissing;

        Bool_t bIsCandidateK0s ; // candidate for K0s
        Bool_t bIsCandidateLambda ; // candidate for Lambda
        Bool_t bIsCandidateALambda ; // candidate for anti-Lambda
        Bool_t bIsInPeakK0s ; // candidate within the K0s mass peak
        Bool_t bIsInPeakLambda ; // candidate within the Lambda mass peak
        Bool_t bIsInPeakALambda ; // candidate within the anti-Lambda mass peak
        Bool_t bIsInConeJet ; // candidate within the jet cones
        Bool_t bIsInConePerp ; // candidate within a perpendicular cone
        Bool_t bIsInConeRnd ; // candidate within the random cone
        Bool_t bIsInConeMed ; // candidate within the median-cluster cone
        Bool_t bIsOutsideCones ; // candidate outside excluded cones

        void Reset() {memset(this,0, sizeof(*this)); bIsCandidateK0s=bIsCandidateLambda=bIsCandidateALambda=kTRUE;}
        void Print() const;
    };

    //_________________________
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

    //__________________________
    //basic stuff
    void localtoglobal(double alpha, double *local, double *global);
   // void EventwiseCleanup();
    AliVParticle * GetVParticleMother(AliVParticle *part);
    //Double_t GetLocalAlphaAOD(AliAODTrack *track);
    //Double_t GetTrackCurvature(AliAODTrack *track);
    //Double_t GetLocalThetaAOD(AliAODTrack *track);
    Bool_t getJetVtxMass( AliEmcalJet *jet, double &value);
    void SetJetRadius(Double_t fJetRadRead){fJetRadius=fJetRadRead;}

    //int GetMCTruth(AliAODTrack *track, int &motherpdg);
    bool GetPIDCombined(AliAODTrack * track, double *prob, int &nDetectors, UInt_t &usedDet , AliPID::EParticleType &MostProbablePID, bool setTrackPID );
    void setFProductionNumberPtHard(Int_t value=-1)
    {
        fProductionNumberPtHard = value;
    }
    Int_t IsInVector(const vector<Int_t>& vec, Int_t iLabel, TString sFunc);
    Bool_t IsParton(int pdg);
    Bool_t IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax);
    Int_t NDaughterInCone(std::vector<Int_t>& vecDaughLabels, const AliEmcalJet* jet, const AliAODEvent* event, Double_t dRMax, Double_t& ipsig);

    void GetLowUpperBinNo(int &iLowerBin, int &iUpperBin, double min, double max, TString type, Int_t iN);

    void PrintAllTreeVars();

    //____________________________
    //Cuts
    void SetESDCuts (AliESDtrackCuts  *cuts =NULL){fESDTrackCut =  new AliESDtrackCuts(*cuts);}
    void SetDefaultAnalysisCuts();
    void SetDefaultV0Cuts();
    void DefaultInitTreeVars();
    Bool_t IsPhysicalPrimary(AliVParticle *part);
    void ChangeDefaultCutTo(AliAnalysisTaskHFJetIPQA::bCuts cutname, Double_t newcutvalue);
    void GetMaxImpactParameterCutR(const AliVTrack * const track, Double_t &maximpactRcut);

    Bool_t IsTrackAccepted(const AliVTrack* track, int jetflavour);
    Bool_t IsDCAAccepted(double decaylength, double ipwrtjet, Double_t * dca, int jetflavour);
    Bool_t IsEventAccepted(AliAODEvent *ev);

    void GetV0Properties(SV0Cand*&  sV0, AliAODv0* &v0);
    void GetV0DaughProperties(SV0Daugh* & sTrack,AliAODv0* &v0, bool isPos);
    void FillV0Candidates(Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut);
    Int_t IsV0Daughter(const AliAODEvent* fAODIn,const AliAODTrack* track, Int_t iTrack);
    void SelectV0Candidates(const AliAODEvent *fAODIn);
    void IdentifyRecV0PDG(Double_t fMassK0, Double_t fMassLambda, Double_t fMassAntiLambda, Bool_t& isK0, Bool_t& IsLambda, Bool_t& IsAntiLambda, TString sIsCalledBy="");
    Bool_t PerformV0AcceptanceCuts(Double_t V0pt, Double_t V0y, Double_t V0PosDaughpt, Double_t V0PosDaughEta,Double_t V0NegDaughpt, Double_t V0NegDaughEta);
    Bool_t PerformV0MCAcceptanceCuts(const AliAODMCParticle* pAODMother, AliAODMCParticle* pAODPosDaugh,AliAODMCParticle* pAODNegDaugh,Bool_t& bV0MCIsK0s,Bool_t& bV0MCIsLambda,Bool_t& bV0MCIsALambda);
    //void GetGeneratedV0();
    Int_t GetGenV0Jets(const AliEmcalJet* jetgen, const AliAODEvent* event, const std::vector<Int_t>& iTrackLabels, const std::vector<Double_t>& fTrackRecIPs, const std::vector<Double_t>& fTrackRecPts, Int_t fGenJetFlavour, Bool_t **kTagDec, Double_t fLNJP);
    Int_t FindAllV0Daughters(AliAODMCParticle* pAOD, const AliAODEvent* event, const AliEmcalJet* jetgen, const vector<Int_t>& iTrackLabels, const vector<Double_t>& fTrackRecIPs,Int_t iCount, Int_t iLevel);
    void GetGenV0DaughterIP(AliAODMCParticle *pAOD, const AliEmcalJet* jetgen, const AliAODEvent* event, const vector<Int_t>& iTrackLabels, const vector<Double_t>& fTrackRecIPs, Int_t& iInVectorInxMaxIP);
    //AliAODMCParticle* GetMCTrack( const AliAODTrack* track);
    AliAODMCParticle* GetMCTrack(int iLabel);
    int GetV0MCVeto(const AliAODEvent* fAODIn, AliAODv0* v0, Int_t tracklabel);
    void FillV0EfficiencyHists(int isV0, int & jetflavour, double jetpt, bool &isV0Jet);

    void FillCandidateJet(Int_t CutIndex, Int_t JetFlavor);
    bool IsFromElectron(AliAODTrack *track);
    bool IsFromProton(AliAODTrack *track);
    bool IsFromKaon(AliAODTrack *track);
    bool IsFromPion(AliAODTrack *track);


    //_____________________________
    //Impact Parameter Generation
    Bool_t GetImpactParameter(const AliAODTrack *track, const AliAODEvent *event, Double_t *dca, Double_t *cov, Double_t *XYZatDCA);
    Bool_t GetMCIP(const AliAODMCParticle* track,const AliAODEvent *event, const AliEmcalJet* jetgen, Double_t& ipsig);
    Double_t GetIPSign(Double_t *XYZatDCA, Double_t* jetp,Double_t* VxVyVz);
    AliExternalTrackParam GetExternalParamFromJet(const AliEmcalJet *jet, const AliAODEvent *event);
    Bool_t GetImpactParameterWrtToJet(const AliAODTrack *track, const AliAODEvent *event, const AliEmcalJet *jet, Double_t *dca, Double_t *cov, Double_t *XYZatDCA, Double_t &jetsign, int jetflavour);
    int DetermineUnsuitableVtxTracks(int *skipped, AliAODEvent * const aod, AliVTrack * const track);
    void DetermineIPVars(std::vector<AliAnalysisTaskHFJetIPQA::SJetIpPati>& sImpParXY, std::vector<AliAnalysisTaskHFJetIPQA::SJetIpPati> sImpParXYSig, std::vector<Float_t> &ipvalsig, std::vector<Float_t> &ipval, Int_t& HasGoodIPTracks);
    //______________________________
    //Corrections
    double DoUESubtraction(AliJetContainer* &jetcongen, AliJetContainer* &jetconrec, AliEmcalJet* &jetrec, double jetpt);
    AliAODVertex *RemoveDaughtersFromPrimaryVtx(const AliVTrack * const track);

    //_______________________________
    //Filling Histograms
    Bool_t FillTrackHistograms(AliVTrack * track, double * dca , double *cov,double weight);
    void FillRecHistograms(Int_t jetflavour, Double_t recjetpt, Double_t fJetGenPt,Double_t fJetRecEta, Double_t fJetGenEta, Double_t fJetRecPhi, Int_t fUnfoldFracCalc);
    void FillGenHistograms(Int_t jetflavour,Double_t jetgenpt, Int_t fUnfoldFracCalc);
    Bool_t PerformGenLevAcceptanceCuts(Double_t fJetGenEta);
    void FillTaggedJetPtDistribution(bool** kTagDec, double jetpt);

    //________________________________
    //Setters
    void SmearTrack(AliAODTrack *track);
    void setFRunSmearing(Bool_t value){fRunSmearing = value;}
    void setFDoMCCorrection(Bool_t value){fDoMCCorrection=value;}
    void setFDoUnderlyingEventSub(Bool_t value){fDoUnderlyingEventSub=value;}
    void setfDoFlavourMatching(Bool_t value){fDoFlavourMatching=value;}
    void setV0Cut(int iCut,double value){fV0Cuts[iCut]=value;}
    void setAnalysisCuts(bCuts cut, bool cutvalue){fAnalysisCuts[cut]=cutvalue;}
    void setAnalysisCuts(bCuts cut, int cutvalue){fAnalysisCuts[cut]=cutvalue;}
    void setAnalysisCuts(bCuts cut, double cutvalue){fAnalysisCuts[cut]=cutvalue;}

    Double_t getFMCglobalDCAxyShift() const;
    void setFMCglobalDCAxyShift(const Double_t &value);
    Double_t getFVertexRecalcMinPt() const;
    void setFVertexRecalcMinPt(const Double_t &value);
    void setFMCglobalDCASmear(const Double_t value);
    void setFParam_Smear_Sigma(Double_t value){fParam_Smear_Sigma = value;}
    void setFParam_Smear_Mean(Double_t value){fParam_Smear_Mean = value;}
    void setGlobalVertex(Bool_t value){fGlobalVertex = value;}
    void setDoNotCheckIsPhysicalPrimary(Bool_t value){fDoNotCheckIsPhysicalPrimary = value;}
    void setDoJetProb(Bool_t value){fDoJetProb = value;}
    void setDoTCTagging(Int_t value) {fDoTCTagging=value;}
    void setDoProbTagging(Int_t value) {fDoProbTagging=value;}
    void setDoMCEffs(Bool_t value){fDoMCEffs=value;}

    void setTrackIPvsPtValues(double fav0cut, double fbv0cut, double fcv0cut){fV0Cuts[fAV0Cut]=fav0cut;fV0Cuts[fBV0Cut]=fbv0cut;fV0Cuts[fCV0Cut]=fcv0cut;}
    void setfDaughterRadius(Double_t value){fDaughtersRadius=value;}
    void setfNoJetConstituents(Int_t value){fNoJetConstituents=value;}
    void setfNThresholds(Int_t value){fNThresholds=value; printf("Setting threshold value=%i\n",fNThresholds);}
    void setfUserSignificance(Bool_t value){fUseSignificance=value;}
    void SetTagSettings(int iTagSetting, int ntracktypes=-1);
    void SetfUnfoldPseudoDataFrac(int frac){fUnfoldPseudeDataFrac=frac;}
    void setfResponseMode(bool value){fResponseMode=value;}
    void setTaskName(const char* name){sTaskName=Form("%s",name);}

    //_____________________________
    //Lund Plane
    void RecursiveParents(const AliEmcalJet *fJet,const AliJetContainer *fJetCont, std::vector<Int_t> &fGroomedJetConst);   //Based on AliAnalysisTaskEmcalQGTagging::RecursiveParents

    //_____________________________
    //Track Counting
    enum TaggingType{
          Full,
          Single1st,
          Single2nd,
          Single3rd,
          Double,
          Triple,
    };

    void DoTCTagging(Float_t jetpt, Int_t nGoodIPTracks, const std::vector<Float_t>& ipval, Bool_t **kTagDec);
    void DoProbTagging(double probval, double jetpt, bool** kTagDec);
    void SetTCThresholds(TObjArray** &threshs);
    void SetProbThresholds(TObjArray** &threshs);
    void ReadProbvsIPLookup(TObjArray *&oLookup);
    void ReadThresholdHists(TString PathToThresholds, TString taskname, int nTCThresh, int iTagSetting, int ntracktypesprob);
    void setTagLevel(int taglevel){kTagLevel=taglevel;}
    void setTCThresholdPtFixed(double value){fTCThresholdPtFixed=value;};

    //________________________________
    //Probability Tagging
    Float_t GetTrackProbability(Float_t jetpt, Int_t nGoodIPTracks, const std::vector<Float_t>& ipval);
    void GetDeltaRij(const AliAODTrack* track, const AliEmcalJet *jet);
    void setDoLundPlane(Bool_t dolundplane){fDoLundPlane=dolundplane;}
    Float_t IntegrateIP(Float_t jetpt, Float_t IP, Int_t iN);

    void useTreeForCorrelations(Bool_t value){fUseTreeForCorrelations = value;}
    //virtual Bool_t IsEventSelected();
    void FillCorrelations(bool bn[3], double v[3], double jetpt);
    void setFFillCorrelations(const Bool_t &value);
    virtual void SetPtHardBin(Int_t b){ fSelectPtHardBin = b;}
    void SetHardCutoff(Double_t t)                            {fHardCutOff = t;}


public:
    AliEventCuts fEventCuts;

private:
    THistManager         fHistManager    ;///< Histogram manager
    AliAODVertex * fEventVertex;//!
    AliPIDResponse *fPidResponse ;//!
    AliEmcalJet *  GetPerpendicularPseudoJet (AliEmcalJet*jet_in  , bool rev );

    TTree* tJetTree;      //! tree containing jet properties
    Float_t fJetRecPt;
    Float_t fJetArea;
    Float_t fMatchedJetPt;
    Float_t fJetProb;
    Float_t fMeanLNKt;
    Float_t fMeanTheta;
    Float_t fMeanLNKtSD;
    Float_t fMeanThetaSD;
    Float_t fJetMass;
    Int_t fJetFlavour;
    Int_t nTracks;
    Int_t fNEvent;
    bool bMatched;
    Int_t bIsTrueGenV0Jet;
    Float_t fTrackIPs[40];
    Float_t fTrackIPSigs[40];
    Float_t fTrackProb[40];
    Float_t fTrackChi2OverNDF[40];
    Float_t fTrackPt[40];
    Float_t fDeltaRij[40];
    Float_t fV0MotherPt[40];
    Float_t fV0MotherPtMC[40];
    Float_t fV0MotherEta[40];
    Float_t fV0MotherEtaMC[40];
    Int_t iTrackITSHits[40];
    Int_t iV0MCID[40];
    Int_t iV0RecID[40];
    Int_t bTrackIsV0[40];
    Bool_t bPassedSD[40];
    //Bool_t bFull[30];
    Bool_t bSingle1st[30];
    //Bool_t bSingle2nd[30];
    //Bool_t bSingle3rd[30];
    Bool_t bDouble[30];
    //Bool_t bTriple[30];

    void FillParticleCompositionSpectra(AliEmcalJet * jet,const char * histname );
    void DoJetLoop(); //jet matching function 2/4
    void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Int_t matching=0);
    void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
    void SmearTrackHybrid(AliVTrack * track);
    void FillHist(const char * name,Double_t x ,Double_t w);
    void FillHist(const char * name,Double_t x, Double_t y,Double_t w);
    void IncHist(const char * name,Int_t bin);
    void SubtractMean (Double_t val[2],AliVTrack *track);
    Bool_t MatchJetsGeometricDefault(); //jet matching function 1/4
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
    Int_t IsMCJetPartonFast(const AliEmcalJet *jet,  Double_t radius,Bool_t &is_udg);
    Int_t GetRunNr(AliVEvent * event){return event->GetRunNumber();}
    Double_t GetPtCorrected(const AliEmcalJet* jet);
    //Double_t GetPtCorrectedMC(const AliEmcalJet *jet);
    void PrintSettings();
    void PrintV0Settings();


    //Functions to allow jet probability/TC System 8 efficiency estimation
    Bool_t IsJetTaggedTC(int n =0 ,double thres = 0.1);
    Bool_t IsJetTaggedJetProb(double thresProb = 0.90);
    TH1 *  AddHistogramm(const char * name,const char * title,Int_t x,Double_t xlow,Double_t xhigh, Int_t y=0,Double_t ylow=0,Double_t yhigh=0);
    TH1D * GetHist1D(const char * name){return (TH1D*)fOutput->FindObject(name);}
    TH2D * GetHist2D(const char * name){return (TH2D*)fOutput->FindObject(name);}


private:
    AliJetContainer*  jetconrec;
    AliJetContainer*  jetcongen;

    //___________________
    //Booleans for settings
    Bool_t   fRunSmearing;//
    Bool_t   fUsePIDJetProb;//
    Bool_t   fDoMCCorrection;//  Bool to turn on/off MC correction. Take care: some histograms may still be influenced by weighting.
    Bool_t   fDoUnderlyingEventSub;//

    Bool_t   fDoFlavourMatching;//
    Double_t fParam_Smear_Sigma;//
    Double_t fParam_Smear_Mean;//
    Bool_t   fGlobalVertex;//
    Bool_t fDoNotCheckIsPhysicalPrimary;//
    Bool_t fDoJetProb; //
    Bool_t   fFillCorrelations;//
    Bool_t fDoLundPlane;//
    Int_t fDoTCTagging;//  //0: no TC tagging, 1: IP Significance tagging, 2: IP tagging, fixed threshold
    Int_t fDoProbTagging;//  //0: no probability tagging, 1: use JP for tagging, 2: use lnJP for tagging
    Bool_t fDoMCEffs;//
    Bool_t fUseSignificance;//
    Bool_t fResponseMode;//

    //_____________________
    //variables
    int kTagLevel; //1: accept single splittings, 2: accept only 2+3, 3: accept only 3 for track counting algorithm
    std::vector<double > fFracs;
    Float_t fXsectionWeightingFactor;//
    Int_t   fProductionNumberPtHard;//
    Int_t fNThresholds;//
    Int_t fNTrackTypes;//
    std::vector<TString> sTemplateFlavour;
    Int_t fUnfoldPseudeDataFrac;//
    TString sTaskName;//

    //______________________
    //Cuts
    Double_t fJetRadius;//
    Double_t fDaughtersRadius;//
    Int_t fNoJetConstituents;//
    Double_t fTCThresholdPtFixed; //

    //*********************************
    //Histograms

    //__________________________
    //Histograms for track counting
    std::vector<TH1D*> h1DThresholdsFirst; //0-> single probability, 1-> double probability, 2-> tripple probability
    std::vector<TH1D*> h1DThresholdsSecond; //
    std::vector<TH1D*> h1DThresholdsThird;//

    //_____________________________
    //Histograms for probability tagging
    std::vector<TH2D*> h2DProbLookup;//
    std::vector<TH1D*> h1DProbThresholds;//

    //______________________________
    //Cut Histograms
    TCanvas *cCuts; //

    TH1D *fh1DCutInclusive;
    TH1D *fh1dCutudg;
    TH1D *fh1dCutc;
    TH1D *fh1dCutb;
    TH1D *fh1dCuts;

    TH1D *fh1dTracksAccepeted; //!
    TH1D* fh1dCutsPrinted;//!

    THnSparse *fHLundIterative;//!       iterative declustering

    THnSparse* fhnV0InJetK0s; //! V0 in jet cones, in a centrality bin, m_V0; pt_V0; eta_V0; pt_jet
    THnSparse* fhnV0InJetLambda; //!
    THnSparse* fhnV0InJetALambda; //!
    TH1D* fh1V0CounterCentK0s; //! number of K0s candidates after various cuts
    TH1D* fh1V0CounterCentLambda; //! number of Lambda candidates after various cuts
    TH1D* fh1V0CounterCentALambda; //! number of ALambda candidates after various cuts
    TH2D* fh2dKshortMassVsPt; //!
    TH2D* fh2dLamdaMassVsPt; //!
    TH2D* fh2dAnLamdaMassVsPt; //!

    TH1D* h1DV0FalseRec; //!
    TH1D* h1DV0TrueRec; //!
    TH1D* h1DV0TrueDataDef; //!
    TH1D* h1DV0TrueMCDef; //!

    TH1D* fh1dKshortPtMC;//!
    TH1D* fh1dLamdaPtMC;//!
    TH1D* fh1dAnLamdaPtMC;//!
    TH1D* fh1dKshortEtaMC;//!
    TH1D* fh1dLamdaEtaMC;//!
    TH1D* fh1dAnLamdaEtaMC;//!
    TH2D *fh2dKshortPtVsJetPtMC;//!
    TH2D *fh2dLamdaPtVsJetPtMC;//!
    TH2D *fh2dAnLamdaPtVsJetPtMC;//!

    //________________________________
    //vectors
    TClonesArray     *fMCArray;//!
    AliMCEvent       *fMCEvent;//!
    AliESDtrackCuts  *fESDTrackCut;//
    AliVertexerTracks *fVertexer;//!
    TClonesArray* fV0CandidateArray;//!
    Bool_t fMcEvtSampled;//
    Double_t fBackgroundFactorLinus[21][498]; //[21][498]FineBinned correction factors up 0.1-25 GeV/c first value below last above 0.05 binwidth
    std::vector <Double_t > fPUdsgJet;//!
    std::vector <Double_t > fPSJet;//!
    std::vector <Double_t > fPCJet;//!
    std::vector <Double_t > fPBJet;//!
    std::vector <Double_t > fJetCont;//!
    std::map<int, int> daughtermother;//!

    TGraph fResolutionFunction[200];//[200]<-
    Double_t fAnalysisCuts[33]; // /Additional (to ESD track cut or AOD filter bits) analysis cuts.
    Double_t fV0Cuts[26];

    AliPIDCombined *fCombined ;//!

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

   /* Bool_t GetMixDCA(int n , double &v){
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
    }*/

   ClassDef(AliAnalysisTaskHFJetIPQA, 79)
};

#endif


