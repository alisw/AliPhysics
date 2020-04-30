#ifndef AliAnalysisTaskeeCor_cxx
#define AliAnalysisTaskeeCor_cxx

class TH1F;
class TH1I;
class TH3;
class AliMCEvent;
class TList;
class AliAODEvent;
class AliPIDResponse;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskeeCor : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskeeCor();
    AliAnalysisTaskeeCor(const char *name);
    virtual ~AliAnalysisTaskeeCor();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetTrackCuts(AliESDtrackCuts* const cuts){ fTrackCuts = cuts;}
    void SetMinPtCut(Double_t ptMin){ fPtMinCut = ptMin;}
    void SetEleTPCcuts(Double_t TPCmin, Double_t TPCmax){ fTPCmin = TPCmin; fTPCmax = TPCmax;}
    void SetEleTOFcuts(Double_t TOFmin, Double_t TOFmax){ fTOFmin = TOFmin; fTOFmax = TOFmax;}
    void SetPionTPCrej(Double_t TPCmin, Double_t TPCmax){ fTPCminPio = TPCmin; fTPCmaxPio = TPCmax;}
    void SetProtTPCrej(Double_t TPCmin, Double_t TPCmax){ fTPCminPro = TPCmin; fTPCmaxPro = TPCmax;}
    void SetKaonTPCrej(Double_t TPCmin, Double_t TPCmax){ fTPCminKao = TPCmin; fTPCmaxKao = TPCmax;}
    void SetRecabPID(Bool_t use){ fRecabPID = use;}
    void SetTPCmeanPIDmap(TH2F *h1){ fTPCmean  = (TH2F*)h1->Clone(); if (!fTPCmean)  printf("ERROR -> COULD NOT ADD!\n");}
    void SetTPCwidthPIDmap(TH2F *h2){fTPCwidth = (TH2F*)h2->Clone(); if (!fTPCwidth) printf("ERROR -> COULD NOT ADD!\n");}
    void SetTOFmeanPIDmap(TH2F *h3){ fTOFmean  = (TH2F*)h3->Clone(); if (!fTOFmean)  printf("ERROR -> COULD NOT ADD!\n");}
    void SetTOFwidthPIDmap(TH2F *h4){fTOFwidth = (TH2F*)h4->Clone(); if (!fTOFwidth) printf("ERROR -> COULD NOT ADD!\n");}
    void SetSmearing(Bool_t use){ fUseSmearing = use;}
    void SetDCASmearingByMath(Bool_t use){fDCASmearingByMath = use;}
    void SetDCASmearingByMaps(Bool_t use){fDCASmearingByMaps = use;}
    void SetDCASmearingByPars(Bool_t use, Bool_t smrMethod){fDCASmearingByPars = use; fDCAparSmr = smrMethod;}
    void SetPtSmrMap(TH2F *s1){ fPtSmr  = (TH2F*)s1->Clone(); if (!fPtSmr)  printf("ERROR -> COULD NOT ADD!\n");}
    void SetEtaSmrMap(TH2F *s2){fEtaSmr = (TH2F*)s2->Clone(); if (!fEtaSmr) printf("ERROR -> COULD NOT ADD!\n");}
    void SetPhiEleSmrMap(TH2F *s3){fPhiEleSmr = (TH2F*)s3->Clone(); if (!fPhiEleSmr) printf("ERROR -> COULD NOT ADD!\n");}
    void SetPhiPosSmrMap(TH2F *s4){fPhiPosSmr = (TH2F*)s4->Clone(); if (!fPhiPosSmr) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCASmrMap0(TH2F *s5){fDCASmrMap0 = (TH2F*)s5->Clone(); if (!fDCASmrMap0) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCASmrMap1(TH2F *s6){fDCASmrMap1 = (TH2F*)s6->Clone(); if (!fDCASmrMap1) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCASmrParCen(TH1D *s7){fDCASmrCen = (TH1D*)s7->Clone(); if (!fDCASmrCen) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCASmrParSig(TH1D *s8){fDCASmrSig = (TH1D*)s8->Clone(); if (!fDCASmrSig) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCASmrParMax(TH1D *s9){fDCASmrMax = (TH1D*)s9->Clone(); if (!fDCASmrMax) printf("ERROR -> COULD NOT ADD!\n");}
    void SetDCAmapsFromMC(Bool_t isMC){fDCASmrFromMC = isMC;}
    void FillPIDrecMaps(Bool_t use){ fFillPIDrecMaps = use;}
    void FillSmrMaps(Bool_t use){ fFillSmrMaps = use;}

    // simple struct for electrons
    struct LMEEparticle {
        Double_t genP;
        Double_t genPt;
        Double_t genTheta;
        Double_t genEta;
        Double_t genPhi;
        Double_t genDCAxy;
        Double_t genDCAxySig;
        Double_t recP;
        Double_t recPt;
        Double_t recTheta;
        Double_t recEta;
        Double_t recPhi;
        Double_t recDCAxy;
        Double_t recDCAxySig;
        Double_t charge;
        Int_t label;
        Int_t mlabel;
        Int_t mPDG;
        Int_t grmlabel;
        Int_t grmPDG;
        Bool_t passedCuts;
        Bool_t isCharm;
        Bool_t isBeauty;
        TLorentzVector genLv;
        TLorentzVector recLv;
        // functions
        LMEEparticle() : genP(-99.),genPt(-99.),genTheta(-99.),genEta(-99.),genPhi(-99.),genDCAxy(-99.),genDCAxySig(-99.),
                         recP(-99.),recPt(-99.),recTheta(-99.),recEta(-99.),recPhi(-99.),recDCAxy(-99.),recDCAxySig(-99.),
                         mlabel(-1),mPDG(-1),grmlabel(-1),grmPDG(-1),passedCuts(kFALSE),genLv(),recLv()
        {
        }
        void MakeGenLV(){ genLv.SetPtEtaPhiM(genPt,genEta,genPhi,0.0005109989); }
        void MakeRecLV(){ recLv.SetPtEtaPhiM(recPt,recEta,recPhi,0.0005109989); }
    };

private:

	Bool_t 		IsCharmedEle(Int_t label);
	Bool_t 		IsBeautyEle(Int_t label);
	void	 	  GetRecalibrationPID(Double_t mom, Double_t eta, Double_t *meanTPC, Double_t *widthTPC, Double_t *meanTOF, Double_t *widthTOF);
	Double_t	GetPtSmr(Double_t pt);
	Double_t	GetEtaSmr(Double_t pt);
	Double_t	GetPhiSmr(Double_t pt, Double_t q);
	Double_t	GetDCASmr(Double_t pt);
	Double_t	GetDCASmrByPar(Double_t pt);
	Bool_t 		GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0);

    // transient members are not streamed
    AliMCEvent  *mcEvent;         //! MC event
    AliAODEvent *aodEvent;        //! AOD event
    TList       *fOutputList;     //! Output list
    TList       *fOutputListGen;  //! Output list for generated
    TList       *fOutputListRec;  //! Output list for reconstructed
    TList       *fOutputListGenCC;  //! Output list for generated
    TList       *fOutputListRecCC;  //! Output list for reconstructed
    TList       *fOutputListGenBB;  //! Output list for generated
    TList       *fOutputListRecBB;  //! Output list for reconstructed

    // Event histograms
    TH1I        *fEventStat;        //! Event statistics

    // generated
    TH1F        *fHistEleGenPt;     //!
    TH1F        *fHistPosGenPt;     //!
    TH2F        *fHistEleGenEtaPhi; //!
    TH2F        *fHistPosGenEtaPhi; //!

    TH1F        *fHistEleGenPhi;     //!
    TH1F        *fHistPosGenPhi;     //!
    TH1F        *fHistEleGenDCAxy; //!
    TH1F        *fHistPosGenDCAxy; //!

    TH1I        *fMotherPDGCC;        //!
    TH1F        *fHistEleGenPtCC;     //!
    TH1F        *fHistPosGenPtCC;     //!
    TH2F        *fHistEleGenEtaPhiCC; //!
    TH2F        *fHistPosGenEtaPhiCC; //!
    TH1I        *fMotherPDGBB;        //!
    TH1F        *fHistEleGenPtBB;     //!
    TH1F        *fHistPosGenPtBB;     //!
    TH2F        *fHistEleGenEtaPhiBB; //!
    TH2F        *fHistPosGenEtaPhiBB; //!

    TH1F        *fHistPairsGenPt;   //!
    TH1F        *fHistPairsGenMass;   //!
    TH1F        *fHistPairsGenDphi;   //!
    TH2F        *fHistPairsGenMassPt;   //!
    TH2F        *fHistPairsGenPhiPt;   //!
    TH2F        *fHistPairsGenMassPt2;   //!
    TH2F        *fHistPairsGenPhiPt2;   //!
    TH3F        *fHistPairsGenPtMasPhi;   //!

    TH1F        *fHistPairsGenPtCC;   //!
    TH1F        *fHistPairsGenMassCC;   //!
    TH1F        *fHistPairsGenDphiCC;   //!
    TH2F        *fHistPairsGenMassPtCC;   //!
    TH2F        *fHistPairsGenPhiPtCC;   //!
    TH2F        *fHistPairsGenMassPt2CC;   //!
    TH2F        *fHistPairsGenPhiPt2CC;   //!
    TH3F        *fHistPairsGenPtMasPhiCC;   //!

    TH1F        *fHistPairsGenPtBB;   //!
    TH1F        *fHistPairsGenMassBB;   //!
    TH1F        *fHistPairsGenDphiBB;   //!
    TH2F        *fHistPairsGenMassPtBB;   //!
    TH2F        *fHistPairsGenPhiPtBB;   //!
    TH2F        *fHistPairsGenMassPt2BB;   //!
    TH2F        *fHistPairsGenPhiPt2BB;   //!
    TH3F        *fHistPairsGenPtMasPhiBB;   //!

    TH1F        *fHistPairsGenPtLSBB;   //!
    TH1F        *fHistPairsGenMassLSBB;   //!
    TH1F        *fHistPairsGenDphiLSBB;   //!
    TH2F        *fHistPairsGenMassPtLSBB;   //!
    TH2F        *fHistPairsGenPhiPtLSBB;   //!
    TH2F        *fHistPairsGenMassPt2LSBB;   //!
    TH2F        *fHistPairsGenPhiPt2LSBB;   //!
    TH3F        *fHistPairsGenPtMasPhiLSBB;   //!



    TH2F        *fHistPairsGen;   //!

    TH2F        *fHistPairsCCgen;   //!
    TH2F        *fHistPairsBBgen;   //!
    TH2F        *fHistPairsBBgenLS;   //!

    TH2F        *fHistPairsCCgen2;   //! ready to be used
    TH2F        *fHistPairsBBgen2;   //! ready to be used
    TH2F        *fHistPairsBBgenLS2;   //! ready to be used

    // reconstructed
    TH1F        *fHistEleRecPt;     //!
    TH1F        *fHistPosRecPt;     //!
    TH2F        *fHistEleRecEtaPhi; //!
    TH2F        *fHistPosRecEtaPhi; //!

    TH1F        *fHistEleRecPhi;     //!
    TH1F        *fHistPosRecPhi;     //!
    TH1F        *fHistEleRecDCAxy; //!
    TH1F        *fHistPosRecDCAxy; //!
    TH1F        *fHistEleRecDCAxySig; //!
    TH1F        *fHistPosRecDCAxySig; //!
    TH2F        *fHistEleRecPtDCAxySig; //!
    TH2F        *fHistPosRecPtDCAxySig; //!

    TH1F        *fHistEleRecDCAxyCC; //!
    TH1F        *fHistPosRecDCAxyCC; //!
    TH1F        *fHistEleRecDCAxyBB; //!
    TH1F        *fHistPosRecDCAxyBB; //!

    TH1F        *fHistEleRecPtCC;     //!
    TH1F        *fHistPosRecPtCC;     //!
    TH2F        *fHistEleRecEtaPhiCC; //!
    TH2F        *fHistPosRecEtaPhiCC; //!
    TH1F        *fHistEleRecPtBB;     //!
    TH1F        *fHistPosRecPtBB;     //!
    TH2F        *fHistEleRecEtaPhiBB; //!
    TH2F        *fHistPosRecEtaPhiBB; //!

    TH3F        *TPCnSigmaEle_Eta_P_lin; //!
    TH3F        *TOFnSigmaEle_Eta_P_lin; //!

    TH2F        *TPCnSigmaEle_Eta_P_lin2; //!
    TH2F        *TOFnSigmaEle_Eta_P_lin2; //!

    // tracking info
    TH1I        *fHistNclsTPC;      //!
    TH1I        *fHistNclsSTPC;     //!
    TH1I        *fHistNcrTPC;       //!
    TH1I        *fHistNclsITS;      //!
    TH1F        *fHistDCAxy;        //!
    TH1F        *fHistDCAz;         //!
    TH1F        *fHistChi2perNDF;   //!

    // PID info
    TH2F        *fHistTPCnSigmaEle; //!
    TH2F        *fHistTPCnSigmaPio; //!
    TH2F        *fHistTOFnSigmaEle; //!

    // pairs
    TH1F        *fHistPairsRecPt;   //!
    TH1F        *fHistPairsRecMass;   //!
    TH1F        *fHistPairsRecDphi;   //!
    TH2F        *fHistPairsRecMassPt;   //!
    TH2F        *fHistPairsRecPhiPt;   //!
    TH2F        *fHistPairsRecMassPt2;   //!
    TH2F        *fHistPairsRecPhiPt2;   //!
    TH3F        *fHistPairsRecPtMasPhi;   //!

    TH1F        *fHistPairsRecPtCC;   //!
    TH1F        *fHistPairsRecMassCC;   //!
    TH1F        *fHistPairsRecDphiCC;   //!
    TH2F        *fHistPairsRecMassPtCC;   //!
    TH2F        *fHistPairsRecPhiPtCC;   //!
    TH2F        *fHistPairsRecMassPt2CC;   //!
    TH2F        *fHistPairsRecPhiPt2CC;   //!
    TH3F        *fHistPairsRecPtMasPhiCC;   //!

    TH1F        *fHistPairsRecPtBB;   //!
    TH1F        *fHistPairsRecMassBB;   //!
    TH1F        *fHistPairsRecDphiBB;   //!
    TH2F        *fHistPairsRecMassPtBB;   //!
    TH2F        *fHistPairsRecPhiPtBB;   //!
    TH2F        *fHistPairsRecMassPt2BB;   //!
    TH2F        *fHistPairsRecPhiPt2BB;   //!
    TH3F        *fHistPairsRecPtMasPhiBB;   //!

    TH1F        *fHistPairsRecPtLSBB;   //!
    TH1F        *fHistPairsRecMassLSBB;   //!
    TH1F        *fHistPairsRecDphiLSBB;   //!
    TH2F        *fHistPairsRecMassPtLSBB;   //!
    TH2F        *fHistPairsRecPhiPtLSBB;   //!
    TH2F        *fHistPairsRecMassPt2LSBB;   //!
    TH2F        *fHistPairsRecPhiPt2LSBB;   //!
    TH3F        *fHistPairsRecPtMasPhiLSBB;   //!

    TH2F        *fHistPairsRec;   //! rec dielectrons
    TH1F        *fHistPairsDCAxy;   //! rec dielectrons
    TH1F        *fHistPairsDCAxySig;   //! rec dielectrons
    TH2F        *fHistPairsPtDCAxySig;   //! rec dielectrons
    TH1F        *fHistPairsDCAxyCC;   //! rec dielectrons
    TH1F        *fHistPairsDCAxySigCC;   //! rec dielectrons
    TH2F        *fHistPairsPtDCAxySigCC;   //! rec dielectrons
    TH1F        *fHistPairsDCAxyBB;   //! rec dielectrons
    TH1F        *fHistPairsDCAxySigBB;   //! rec dielectrons
    TH2F        *fHistPairsPtDCAxySigBB;   //! rec dielectrons
    TH1F        *fHistPairsDCAxyLSBB;   //! rec dielectrons
    TH1F        *fHistPairsDCAxySigLSBB;   //! rec dielectrons
    TH2F        *fHistPairsPtDCAxySigLSBB;   //! rec dielectrons
    
    TH1F        *fHistPairsGenDCAxy;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxySig;   //! rec dielectrons
    TH2F        *fHistPairsGenPtDCAxySig;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxyCC;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxySigCC;   //! rec dielectrons
    TH2F        *fHistPairsGenPtDCAxySigCC;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxyBB;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxySigBB;   //! rec dielectrons
    TH2F        *fHistPairsGenPtDCAxySigBB;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxyLSBB;   //! rec dielectrons
    TH1F        *fHistPairsGenDCAxySigLSBB;   //! rec dielectrons
    TH2F        *fHistPairsGenPtDCAxySigLSBB;   //! rec dielectrons

    TH2F        *fHistPairsCCrec;   //! rec dielectrons from cc
    TH2F        *fHistPairsBBrec;   //! rec dielectrons from bb
    TH2F        *fHistPairsBBrecLS;   //! rec dielectrons from bb

    TH2F        *fHistPairsCCrec2;   //! rec dielectrons from cc
    TH2F        *fHistPairsBBrec2;   //! rec dielectrons from bb
    TH2F        *fHistPairsBBrecLS2;   //! rec dielectrons from bb

    // Smearing maps
    TH2F        *fPRec_Gen;  //!
    TH2F        *fPtRec_Gen;  //!
    TH2F        *fPtRecOverGen;  //!
    TH2F        *fEtaRec_Gen; //!
    TH2F        *fElePhiRec_Gen; //!
    TH2F        *fPosPhiRec_Gen; //!
    TH2F        *fPairDCARec_Gen; //!
    TH2F        *fEleDCARes; //!
    TH2F        *fPosDCARes; //!
    TH2F        *fEleDCARes2; //!
    TH2F        *fPosDCARes2; //!
    TH2F        *fEleDCARes3; //!
    TH2F        *fPosDCARes3; //!

    AliPIDResponse *fPIDResponse;   //! PID response object
    // persistent members are streamed (copied/stored)
    AliESDtrackCuts *fTrackCuts; // Track cuts
    Double_t		fPtMinCut; //Gen Track cut
    Bool_t			fRecabPID;//
    TH2F			*fTPCmean; //
    TH2F			*fTPCwidth; //
    TH2F			*fTOFmean; //
    TH2F			*fTOFwidth; //
    Bool_t			fUseSmearing;//
    Bool_t			fDCASmearingByMath;//
    Bool_t			fDCASmearingByMaps;//
    Bool_t			fDCASmearingByPars;//
    Bool_t			fDCAparSmr;//
    Bool_t			fFillPIDrecMaps;//
    Bool_t			fFillSmrMaps;//
    Bool_t			fDCASmrFromMC;//
    TH2F			*fPtSmr; //
    TH2F			*fEtaSmr; //
    TH2F			*fPhiEleSmr; //
    TH2F			*fPhiPosSmr; //
    TH2F			*fDCASmrMap0; //
    TH2F			*fDCASmrMap1; //
    TH1D			*fDCASmrCen; //
    TH1D			*fDCASmrSig; //
    TH1D			*fDCASmrMax; //
    Double_t		fTPCmin; //
    Double_t		fTPCmax; //
    Double_t		fTOFmin; //
    Double_t		fTOFmax; //
    Double_t		fTPCminPio; //
    Double_t		fTPCmaxPio; //
    Double_t		fTPCminPro; //
    Double_t		fTPCmaxPro; //
    Double_t		fTPCminKao; //
    Double_t		fTPCmaxKao; //

    enum {kPythiaCC=0, kPythiaBB, kPythiaB, kOther, kNbinsEvent};

    AliAnalysisTaskeeCor(const AliAnalysisTaskeeCor&); // not implemented
    AliAnalysisTaskeeCor& operator=(const AliAnalysisTaskeeCor&); // not implemented

    ClassDef(AliAnalysisTaskeeCor, 1);
};

#endif
