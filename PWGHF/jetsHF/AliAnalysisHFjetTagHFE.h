#ifndef ALIANALYSISTASKHFJETTAGHFE_H
#define ALIANALYSISTASKHFJETTAGHFE_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class TClonesArray;
class TArrayF;
class AliAODEvent;  // sample
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliPIDResponse;
class AliAODMCHeader;
class AliAODMCParticle;  // sample
class AliMultSelection;
class TRandom;

#include "TClonesArray.h"
#include "TObjArray.h"
#include "TObject.h"
//#include "AliAODMCParticle"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFQnVectorHandler.h"
#include "AliQnCorrectionsManager.h"
#include "TProfile.h"

class AliAnalysisHFjetTagHFE : public AliAnalysisTaskEmcalJet {
  public:
    // for Event Plane.
    enum EventPlaneMeth { kTPC, kTPCVZERO, kVZERO, kVZEROA, kVZEROC, kPosTPCVZERO, kNegTPCVZERO };  // Event Plane to be calculated in the task
    enum FlowMethod { kEP, kSP, kEvShapeEP, kEvShapeSP, kEPVsMass, kEvShapeEPVsMass };              // Event Plane, Scalar Product or Event Shape Engeneering methods
    enum EventPlaneDet { kNone = -1, kFullTPC, kPosTPC, kNegTPC, kFullV0, kV0A, kV0C };
    enum q2Method { kq2TPC, kq2PosTPC, kq2NegTPC, kq2VZERO, kq2VZEROA, kq2VZEROC };

    AliAnalysisHFjetTagHFE();
    AliAnalysisHFjetTagHFE(const char *name);
    virtual ~AliAnalysisHFjetTagHFE();

    void UserCreateOutputObjects();
    void Terminate(Option_t *option);

    void SetCentralityMimHFEjet(Int_t centMim) { fcentMim = centMim; };
    void SetCentralityMaxHFEjet(Int_t centMax) { fcentMax = centMax; };
    void SetDebugHFEjet(Bool_t dbHFEj) { idbHFEj = dbHFEj; };
    void SetHybridTrack(Bool_t Hybrid) { iHybrid = Hybrid; };
    void SetOccCorr(Bool_t OccCorr) { iOccCorr = OccCorr; };
    void SetPPcoll(Bool_t ppcoll) { ippcoll = ppcoll; };
    void SetMinSig(Double_t mimSig) { fmimSig = mimSig; };
    void SetMinEop(Double_t mimEop) { fmimEop = mimEop; };
    void SetMinM20(Double_t mimM20) { fmimM20 = mimM20; };
    void SetMaxM20(Double_t maxM20) { fmaxM20 = maxM20; };
    void SetJetEtaCut(Double_t JetEtaCut) { fJetEtaCut = JetEtaCut; };
    void SetEleEtaCut(Double_t EleEtaCut) { fEleEtaCut = EleEtaCut; };
    void SetMCdata(Bool_t mcData) { fmcData = mcData; };
    void SetInvMassCut0(Double_t InvmassCut) { fInvmassCut = InvmassCut; };
    void SetInvMassCut1(Double_t ptAssocut) { fptAssocut = ptAssocut; };
    void SetMCcorr(Bool_t MCcorr) { iMCcorr = MCcorr; };
    void SetDCApTweight(Bool_t DCApTweight) { iDCApTweight = DCApTweight; };
    // void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    // void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1 = flagTr1; };
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2 = flagTr2; };
    void SetMCeta(Bool_t MCEtaFull) { iMCEtaFull = MCEtaFull; };
    void SetSS(Bool_t SSlong) { iSSlong = SSlong; };
    void SetPtHardMax(Double_t PtHardMax) { fPtHardMax = PtHardMax; };
    void SetUEflow(Double_t UEflow){fUEflow = UEflow;};

    /////////////////////////
    ////// Event Plane //////
    /////////////////////////
    void SetTenderTaskName(TString name) { fTenderTaskName = name; }
    void SetHarmonic(int harmonic) { fHarmonic = harmonic; }
    void SetCalibraionType(int calibtype) { fCalibType = calibtype; }
    void SetNormMethod(int norm) { fNormMethod = norm; }
    void SetOADBFileName(TString filename) {
        fOADBFileName = filename;
        cout << "++++++++++++ " << filename << endl;
    }
    void SetFlowMethod(int meth) { fFlowMethod = meth; }
    void SetEventPlaneDetector(int det) { fEvPlaneDet = det; }
    void SetSubEventDetectors(int detsubA, int detsubB) {
        fSubEvDetA = detsubA;
        fSubEvDetB = detsubB;
    }
    void SetqnMethod(int qnmethod) { fqnMeth = qnmethod; }
    void SetqnPercentileSelection(TString splinesfilepath) {
        fPercentileqn = true;
        fqnSplineFileName = splinesfilepath;
    }
    void SetTPCHalvesEtaGap(double etagap = 0.2) { fEtaGapInTPCHalves = etagap; }

    // Correct Ntracklet
    void SetWeightNtrkl(TH1D *hWeight) {
        if (fweightNtrkl) delete fweightNtrkl;
        fweightNtrkl = new TH1D(*hWeight);
    }

    void SetMultiProfileLHC16(TProfile *hprof) { fMultiEstimatorAvg = new TProfile(*hprof); }
    TProfile *GetEstimatorHistogram(const AliAODEvent *fAOD);
    void SetNref(Double_t nref) { fNref = nref; };
    void CountNch();

  protected:
    void ExecOnce();
    Bool_t FillHistograms();
    Bool_t Run();
    void CheckClusTrackMatching();

    // for My Analysis.
    double GetDeltaPsiSubInRange(double psi1, double psi2);

    void GetMainQnVectorInfo(double &mainPsin, double &mainMultQn, double mainQn[2], double &SubAPsin, double &SubAMultQn, double SubAQn[2], double &SubBPsin, double &SubBMultQn, double SubBQn[2], AliHFQnVectorHandler *HFQnVectorHandler);

    bool LoadSplinesForqnPercentile();

    AliVEvent *fVevent;  //! event object
    AliMultSelection *fMultSelection;
    TClonesArray *ftrack;
    TClonesArray *fCaloClusters;
    AliAODMCHeader *fMCheader;
    AliPIDResponse *fpidResponse;  //! pid response

    Float_t fcentMim;  // mim. centrality
    Float_t fcentMax;  // max. centrality
    Bool_t idbHFEj;
    Bool_t iHybrid;
    Bool_t iOccCorr;
    Bool_t ippcoll;
    Bool_t fEMCEG1;       // EMcal Threshold EG1
    Bool_t fEMCEG2;       // EMcal Threshold EG2
    Double_t fmimSig;     // max. centrality
    Double_t fmimEop;     // max. centrality
    Double_t fmimM20;     // max. centrality
    Double_t fmaxM20;     // max. centrality
    Double_t fJetEtaCut;  // max. centrality
    Double_t fEleEtaCut;  // max. centrality
    Double_t fInvmassCut;
    Double_t fptAssocut;
    Double_t fNref;
    Int_t Nch;
    Double_t correctednAcc;
    Bool_t fmcData;
    Bool_t iMCcorr;
    Bool_t iDCApTweight;
    Bool_t iMCEtaFull;
    Bool_t iSSlong;
    Double_t fPtHardMax;
    Double_t fUEflow;
    Int_t NembMCpi0;
    Int_t NembMCeta;
    Int_t NembEPOS;
    Int_t NpureMCproc;
    Double_t GetCorrectedNtrackletsD(TProfile *estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult);  //

    // for Event Plane.
    TString fTenderTaskName;  // name of tender task needed to get the calibrated Qn vectors.
    Int_t fHarmonic;
    Int_t fCalibType;
    Int_t fNormMethod;
    TString fOADBFileName;
    Int_t fFlowMethod;
    Int_t fEvPlaneDet;
    Int_t fSubEvDetA;
    Int_t fSubEvDetB;
    Int_t fqnMeth;
    bool fPercentileqn;
    TString fqnSplineFileName;
    bool fLoadedSplines;
    double fEtaGapInTPCHalves;
    double fScalProdLimit;
    TList *fqnSplinesList[6];

    // General histograms
    TH1 **fHistTracksPt;            //! Track pt spectrum
    TH1 **fHistClustersPt;          //! Cluster pt spectrum
    TH1 **fHistLeadingJetPt;        //! Leading jet pt spectrum
    TH2 **fHistJetsPhiEta;          //! Phi-Eta distribution of jets
    TH2 **fHistJetsPtArea;          //! Jet pt vs. area
    TH2 **fHistJetsPtLeadHad;       //! Jet pt vs. leading hadron
    TH2 **fHistJetsCorrPtArea;      //! Jet pt - bkg vs. area
    TH3 *fHistPtDEtaDPhiTrackClus;  //! track pt, delta eta, delta phi to matched cluster
    TH3 *fHistPtDEtaDPhiClusTrack;  //! cluster pt, delta eta, delta phi to matched track

    TH1F *fHistClustDx;  //!
    TH1F *fHistClustDz;  //!

    // for Event Plane.
    TH2F *fEPcorV0AV0C;    // EP correration of V0A and V0C.
    TH2F *fEPcorV0ATPC;    // EP correration of V0A and TPC.
    TH2F *fEPcorV0CTPC;    // EP correration of V0C and TPC.
    TH2F *fsubV0AV0Ccos2;  // histogram of cos(2*Delta Psi between V0A and V0C).
    TH2F *fsubV0ATPCcos2;  // histogram of cos(2*Delta Psi between V0A and TPC).
    TH2F *fsubV0CTPCcos2;  // histogram of cos(2*Delta Psi between V0C and TPC).
    
    // Background flow subtraction
    TProfile *fUEv2;
    TProfile *fUEv2_sys0;
    TProfile *fUEv2_sys1;
    TProfile *fUEv2_sys2;
    TF1 *fUE;

    // for yield
    TH2F *fJetPhicos2_ele;    // c
    TH2F *fJetPhisin2_ele;    //
    TH1F *fInPlane_eJet;      // In-plane yield of Jet tagged as electron
    TH1F *fOutPlane_eJet;     // Out-of-plane yield of Jet tagged as electron
    TH2F *fInPlane_eJet_2D;   // In-plane yield of Jet tagged as electron
    TH2F *fOutPlane_eJet_2D;  // Out-of-plane yield of Jet tagged as electron

    TH2F *fJetPhicos2_HF;
    TH2F *fJetPhisin2_HF;
    TH1F *fInPlane_HFjet;      // In-plane yield of Jet tagged as HF
    TH1F *fOutPlane_HFjet;     // Out-of-plane yield of Jet tagged as HF
    TH2F *fInPlane_HFjet_2D;   // In-plane yield of Jet tagged as HF
    TH2F *fOutPlane_HFjet_2D;  // Out-of-plane yield of Jet tagged as HF

    TH2F *fJetPhicos2_phoLS;   // Jet matched with photonic
    TH2F *fJetPhisin2_phoLS;   //
    TH1F *fInPlane_phoLS;      //
    TH1F *fOutPlane_phoLS;     //
    TH2F *fInPlane_phoLS_2D;   //
    TH2F *fOutPlane_phoLS_2D;  //

    TH2F *fJetPhicos2_phoULS;   //
    TH2F *fJetPhisin2_phoULS;   //
    TH1F *fInPlane_phoULS;      //
    TH1F *fOutPlane_phoULS;     //
    TH2F *fInPlane_phoULS_2D;   //
    TH2F *fOutPlane_phoULS_2D;  //

    TH2F *fJetPhicos2_Had;      //
    TH2F *fJetPhisin2_Had;      //
    TH1F *fInPlane_Hadjet;      //
    TH1F *fOutPlane_Hadjet;     //
    TH2F *fInPlane_Hadjet_2D;   //
    TH2F *fOutPlane_Hadjet_2D;  //
    TH2F *fHistJetBGflow;


    TH1F *fHistMultCent;  //!
    TH2F *fHistZcorr;     //!
    TH1F *fHistCent;      //!
    TH2F *fHistTPCnSigma;
    TH2F *fHistTPCnSigma_ele;
    TH2F *fHistTPCnSigma_had;
    TH2F *fHistTPCnSigma_eMC;
    TH2F *fHistdEdx;
    TH2F *fHistTPCnSigma_p;
    TH2F *fHistTPCnSigma_pele;
    TH2F *fHistTPCnSigma_phad;
    TH2F *fM20;
    TH2F *fM20_ele;
    TH2F *fM20_had;
    TH2F *fM20_2;
    TH2F *fHistEopNsig;
    TH2F *fHistEop;
    TH2F *fHistEopHFE;
    TH2F *fHistEopHad;
    TH2F *fHistEopIso;
    TH2F *fHistEopHFjet;
    TH2F *fHistNsigHFjet;
    TH2F *fHistEtaHFjet;
    TH1F *fHistJetOrg;
    TH2F *fHistJetOrgArea;
    TH1F *fHistJetBG;
    TH1F *fHistJetSub;
    TH1F *fHisteJetOrg;
    TH1F *fHisteJetBG;
    TH1F *fHisteJetSub;
    TH1F *fHistIncEle;
    TH1F *fHistIncEle2;
    TH1F *fHistIncEleInJet0;
    TH1F *fHistIncEleInJet1;
    TH1F *fHistHfEleMC;
    TH1F *fHistHfEleMCreco;
    TH2F *fHistHfEleMCiso;
    TH2F *fHistEleiso;
    TH1F *fHistEle_wISO;
    TH1F *fHistEle_woISO;
    TH1F *fHistPhoEleMC;
    TH1F *fHistPhoEleMCpi0;
    TH1F *fHistPhoEleMCeta;
    TH1F *fHistPhoEleMCreco;
    TH1F *fHistPhoEleMCrecopi0;
    TH1F *fHistPhoEleMCrecoeta;
    TH1F *fHistMCorgPi0;
    TH1F *fHistMCorgEta;
    TH2F *fHistIncjet;
    TH2F *fHistIncjetFrac;
    TH2F *fHistIncjetOrg;
    TH2F *fHistIncjetBG;
    TH2F *fHistHFjet;
    TH1F *fHistHFdijet;
    TH2F *fHistHFdijetCorr;
    TH2F *fHistULSjet;
    TH2F *fHistLSjet;
    TH2F *fHistHadjet;
    TH2F *fHistHFjet_DCA;
    TH2F *fHistULSjet_DCA;
    TH2F *fHistLSjet_DCA;
    TH2F *fHistHadjet_DCA;
    THnSparse *fHistHFjet_ridge;
    TH2F *fHistHFjetOrder;
    TH2F *fHistDiJetPhi;
    TH2F *fHistDiJetMomBalance;
    TH2F *fHistDiJetMomBalance_All;
    TH2F *fHistDiJetPhi_MC;
    TH2F *fHistDiJetMomBalance_MC;
    TH2F *fHistQjet;
    // TH2F                        *fHistQjet_mult;
    // TH2F                        *fHistGjet_mult;
    TH1F *fHistQjet_mult;
    TH1F *fHistGjet_mult;
    TH2F *fHistJetWidthIncjet;
    TH2F *fHistJetWidthQjet;
    TH2F *fHistIncjetCont;
    TH2F *fHistQjetCont;
    TH2F *fHistIncjetR;
    TH2F *fHistQjetR;
    TH2F *fHistIncjetPhi;
    TH2F *fHistQjetPhi;
    TH2F *fHistIncjetPtD;
    TH2F *fHistQjetPtD;
    TH2F *fInvmassULS;
    TH2F *fInvmassLS;
    TH2F *fInvmassHFuls;
    TH2F *fInvmassHFls;
    TH1F *fLxy_ls;
    TH1F *fLxy_uls;
    TH2D *feJetCorr;
    THnSparse *HFjetCorr0;
    THnSparse *HFjetCorr1;
    THnSparse *HFjetCorr2;
    THnSparse *HFjetCorr3;
    THnSparse *HFjetCorrUE;
    THnSparse *HFjetCorrUE_bgfrac;
    THnSparse *HFjetParticle;
    TH2D *HFjetDCA_c;
    TH2D *HFjetDCA_b;
    TH2D *HFjetDCA_b_FONLL;
    TH2D *HFjetDCA_Dp;
    TH2D *HFjetDCA_Dz;
    TH2D *HFjetDCA_Ds;
    TH2D *HFjetDCA_Lc;
    TH2D *HFjetDCA_Dp_FONLL;
    TH2D *HFjetDCA_Dz_FONLL;
    TH2D *HFjetDCA_Ds_FONLL;
    TH2D *HFjetDCA_Lc_FONLL;
    TH1F *fQAHistJetPhi;
    TH1F *fQAHistTrPhiJet;
    TH1F *fQAHistTrPhi;
    TH1F *fQAHistNits;
    TH2F *fQAHistEleDCAxy;
    TH2F *fQAHistEleDCAz;
    TH2F *fQAHistTrPtJet;
    TH1F *fHistClustE;
    TH1F *fHistClustEtime;
    TH2F *fEMCClsEtaPhi;
    TH1F *fHistRho;
    TH1F *fHistBGfrac;
    TH1F *fHistBGfracHFEev;
    TH1F *fHistBGrandHFEev;
    TH2D *fHistNtrBGfrac;
    TH1F *fHistUE_org;
    TH1F *fHistUE_true;
    TH1F *fHistUE_reco;
    TH2D *fHistJetEnergyReso;
    TH2D *fHistNmatchJet;
    THnSparse *fHistJetEtaCorr0;
    THnSparse *fHistJetEtaCorr1;
    THnSparse *fHistJetEtaCorr2;
    TH1D *fHistDp_POWHEG;
    TH1D *fHistDz_POWHEG;
    TH1D *fHistDs_POWHEG;
    TH1D *fHistLc_POWHEG;
    TH1D *fHistB_POWHEG;
    TH1D *fHistHFEinJet;
    TH1D *fHistHadroninJet;
    TF1 *fPi0Weight;
    TF1 *fEtaWeight;
    TF1 *fpythia_b;
    TF1 *fpowheg_b;
    TF1 *fpythia_c;
    TF1 *fpowheg_c;
    TF1 *fFONLL_D;
    TF1 *fFONLL_Lc;
    TF1 *fFONLL_B;
    TH2F *fzvtx_Ntrkl;
    TH2F *fzvtx_Ntrkl_Corr;
    TH2F *fNchNtr;
    TF1 *fCorrZvtx;
    TF1 *fCorrNtrkl;
    TH1F *fzvtx_Corr;
    TH1F *fNtrkl_Corr;
    TH2F *fzvtx_Nch;
    TH2F *fNchNtr_Corr;
    TH1D *fweightNtrkl;
    THnSparse *fNtrklcorr_HFjet;
    THnSparse *fNtrklcorr_ULSjet;
    THnSparse *fNtrklcorr_LSjet;
    THnSparse *fNtrklcorr_Hadjet;
    THnSparse *fNtrklEopHFE;
    THnSparse *fNtrklEopHad;
    TH1D *fHistphoPi0MC;
    TH1D *fHistphoEtaMC;
    TH2D *fNtrklRhoarea;
    THnSparse *fHistPtfracB;
    THnSparse *fHistPtfracD;
    TString fbgfracFile;
    TH1D *fDelta_pT;

    // for Event Plane.
    TH1F *fHistEvPlane[6];     // histograms of Event Plane angle.
    TH2F *fEvPlaneVsCentr[6];  // histograms of EP angle vs. centrality.
    TH2F *fEPResolVsCentr[3];  // histograms of Resolution(cosine correlation) vs. centrality

    TRandom *generator;

    AliJetContainer *fJetsCont;              //! Jets
    AliJetContainer *fJetsContPart;          //! Jets particle
    AliParticleContainer *fTracksCont;       //! Tracks
    AliClusterContainer *fCaloClustersCont;  //! Clusters
    Bool_t tagHFjet(AliEmcalJet *jet, double *epT, int MCpid, double &maxpT_e, double &ptUS);
    void MakePriorPbPb(AliEmcalJet *jet, double *epT);
    Double_t ReduceJetEnergyScale(AliEmcalJet *jetC, double *epT, double effval);
    // void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec);
    Double_t CalRandomCone(Double_t HFjetPhi[], Double_t HFjetEta[], Double_t HFjetArea);
    Bool_t isHeavyFlavour(int Mompdg, Bool_t &ich, Bool_t &ibe);
    Bool_t isPhotonic(int Mompdg);
    // void MakeParticleLevelJet(THnSparse *pJet);
    void MakeParticleLevelJet(Double_t &pthard);
    // void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
    // void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};
    void GetFakeHadronJet(Double_t pthad, Double_t *hpTarray, Double_t &rho, Double_t PsinV0A);

    Double_t CalOccCorrection();
    void FindMother(AliAODMCParticle *part, int &label, int &pid, double &ptmom);
    Double_t IsolationCut(Int_t itrack, AliVTrack *track, Double_t TrackPt, Double_t MatchPhi, Double_t MatchEta, Double_t MatchclE);
    virtual void IsolationTrackBase(Int_t itrack, AliVTrack *track, Double_t MatchclE, Double_t &IsoEnergyTrack, Int_t &NtrackCone);
    Double_t CalJetWidth(AliEmcalJet *jetC, TH2F *htmp0, TH2F *htmp1, TH2F *htmp2, TH2F *htmp3);
    void dJetHad(AliAODTrack *asstrack, Double_t jetpT, Double_t phijet, Double_t etajet, Bool_t uls, Bool_t ls);

  private:
    // TClonesArray  *ftrack;

    AliAODEvent *fAOD;
    TClonesArray *fMCarray;
    AliAODMCParticle *fMCparticle;
    AliAODMCParticle *fMCparticleMother;

    // Bool_t fmcData;

    AliAnalysisHFjetTagHFE(const AliAnalysisHFjetTagHFE &);             // not implemented
    AliAnalysisHFjetTagHFE &operator=(const AliAnalysisHFjetTagHFE &);  // not implemented

    //=====Multiplicity===========
    TProfile *fMultiEstimatorAvg;

    // ClassDef(AliAnalysisHFjetTagHFE, 7) // jet sample analysis task
    ClassDef(AliAnalysisHFjetTagHFE, 12)  // jet sample analysis task
};

#endif
