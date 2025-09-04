#ifndef AliAnalysisTaskEmcalHfeTagging_H
#define AliAnalysisTaskEmcalHfeTagging_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class TList;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAODEvent;
class AliJetContainer;
class AliParticleContainer;
class AliPIDResponse;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliEventCuts;
class AliMultiEventInputHandler;
class AliAODMCHeader;
class AliAnalysisUtils;

#include "AliAnalysisTaskEmcalJet.h"
#include <THnSparse.h>



class AliAnalysisTaskEmcalHfeTagging : public AliAnalysisTaskEmcalJet {
public:
    
    enum JetShapeType {
        kMCTrue = 0,   // generated jets only
        kTrueDet =1,  // detector and generated jets
        kData   = 2,  // raw data
        kDetEmb = 3,  //detector embedded jets
        kDetEmbPart = 4,
        kPythiaDef = 5,
        kDetEmbPartPythia=6,
        kGenOnTheFly = 7
    };
    enum JetShapeSub {
        kNoSub = 0,
        kConstSub = 1,
        kDerivSub = 2
    };
    enum JetSelectionType {
        kInclusive = 0,
        kRecoil = 1
    };
    
    enum DerivSubtrOrder {
        kSecondOrder = 0,
        kFirstOrder = 1
    };
    
    AliAnalysisTaskEmcalHfeTagging();
    AliAnalysisTaskEmcalHfeTagging(const char *name);
    virtual ~AliAnalysisTaskEmcalHfeTagging();
    
    void                                UserCreateOutputObjects();
    void                                Terminate(Option_t *option);
    
    //Setters
    void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
    void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ; }
    void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
    void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
    void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ; }
    void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
    void SetRMatching(Float_t f)                              { fRMatching = f ;}
    void SetSelectShapes(Int_t c)                             {fSelectedShapes = c;}
    void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
    void SetAngularWindowRecoilJet (Float_t t)                {fangWindowRecoil = t; }
    Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
    Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
    void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
    void SetOneConstSelectionOn(Bool_t t)                     { fOneConstSelectOn =t;}
    void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
    void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
    void SetSemigoodCorrect(Int_t yesno)                      {fSemigoodCorrect=yesno;}
    void SetHolePos(Float_t poshole)                          { fHolePos = poshole; }
    void SetHoleWidth(Float_t holewidth)                      { fHoleWidth = holewidth; }
    void SetDerivativeSubtractionOrder(Int_t c)               { fDerivSubtrOrder = c;}
    void SetMCweight(Int_t c)                                 { fMCweight = c;}
    void SetAssPtCut(Double_t d)                              { fAssPtCut = d;}
    void SetITSncut(Int_t c)                                  { fITSncut = c;}
    void SetAssTPCnCut(Int_t c)                               { fAssTPCnCut = c;}
    void SetTPCnCut(Int_t c)                                  { fTPCnCut = c;}
    void SetSigmaTOFcut(Double_t d)                           { fSigmaTOFcut = d;}
    void SetSigmaTPCcutLowPt(Double_t d)                      { fSigmaTPCcutLowPt = d;}
    void SetSigmaTPCcutHighPt(Double_t d)                     { fSigmaTPCcutHighPt = d;}
    void SetSigmTPCcutExcElec(Double_t d)                     { fSigmTPCcutExcElec = d;}
    void SetDcaXYcut(Double_t d)                              { fDcaXYcut = d;}
    void SetDcaZcut(Double_t d)                               { fDcaZcut = d;}
    void SetIMcut(Double_t d)                                 { fIMcut = d;}
    void SetEtaCut(Double_t d)                                { fEtaCut = d;}
    void SetMinEoPcut(Double_t d)                             { fMinEoPcut = d;}
    void SetMaxEoPcut(Double_t d)                             { fMaxEoPcut = d;}
    void SetM20cut(Double_t d)                                { fM20cut = d;}
    void SetMinM02cut(Double_t d)                             { fMinM02cut = d;}
    void SetMaxM02cut(Double_t d)                             { fMaxM02cut = d;}
    void SetMinPtTPC(Double_t d)                              { fMinPtTPC = d;}
    void SetMaxPtTPC(Double_t d)                              { fMaxPtTPC = d;}
    void SetMinPtEMCal(Double_t d)                            { fMinPtEMCal = d;}
    void SetMaxPtEMCal(Double_t d)                            { fMaxPtEMCal = d;}
    void SetMinPtSemiInclusive(Double_t d)                    { fMinPtSemiInclusive = d;}
    
protected:
    Bool_t                              RetrieveEventObjects();
    Bool_t                              Run();
    Bool_t                              FillHistograms();
    
    void                                GetNumberOfElectrons(AliEmcalJet *jet,Int_t jetContNb, Int_t nMother, Double_t listMother[],  Int_t &nIncElec,  Int_t &nPhotElec, Double_t &pElec, Double_t &ptElec, Bool_t &hasElec);
    void                                GetNumberOfTrueElectrons(AliEmcalJet *jet,Int_t jetContNb, Int_t nMother, Double_t listMother[], Int_t &nTrueElec, Int_t &nTrueHFElec, Double_t &ptTrueHFElec);
    void                                GetWeightAndDecay(AliAODMCParticle *particle, Int_t &decay, Double_t &weight);
    Int_t                               GetNumberOfPairs(AliEmcalJet *jet,AliAODTrack *track,const AliVVertex *pVtx, Int_t nMother, Double_t listMother[],Int_t decay, Double_t weight);
    Bool_t                              IsFromHFdecay(AliAODMCParticle *particle);
    Bool_t                              IsFromLMdecay(AliAODMCParticle *particle);
    Bool_t                              IsPrimary(AliAODMCParticle *particle);
    Bool_t                              HasMother(AliAODMCParticle *particle);
    Double_t                            GetPi0weight(Double_t mcPi0pT) const;
    Double_t                            GetEtaweight(Double_t mcEtapT) const;
    Bool_t                              InclElecTrackCuts(const AliVVertex *pVtx,AliAODTrack *ietrack, Int_t nMother, Double_t listMother[]);
    Bool_t                              PhotElecTrackCuts(const AliVVertex *pVtx,AliAODTrack *aetrack, Int_t nMother, Double_t listMother[]);
    Float_t                             GetJetMass(AliEmcalJet *jet,Int_t jetContNb);
    Float_t                             Angularity(AliEmcalJet *jet, Int_t jetContNb);
    Float_t                             GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb);
    Float_t                             PTD(AliEmcalJet *jet, Int_t jetContNb);
    Float_t                             GetJetpTD(AliEmcalJet *jet, Int_t jetContNb);
    Float_t                             GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb);
	Int_t 								MaxPtBinForSemiInclusiveJet(AliEmcalJet *jet, Int_t jetContNb);
	Int_t								SelectTrigger(Float_t minpT, Float_t maxpT);
    Double_t                            RelativePhi(Double_t mphi, Double_t vphi);
	Double_t 							AngularDifference(AliVParticle* jet1, AliVParticle* jet2);
	AliEmcalJet*						GetClosestOnOtherJetContainer(AliEmcalJet* jet1, AliJetContainer* othercontainer);
	Double_t							GetFractionSharedPtBetweenJets(AliEmcalJet* jet1, AliEmcalJet* jetmatched);
    
    AliAODEvent                         *fAOD;                  //! AOD object
    AliVEvent                           *fVevent;               //! VEvent
    AliPIDResponse                      *fpidResponse;          //! PID response
    TClonesArray                        *fTracksTender;         //Tender for tracks
    const AliVVertex                    *pVtx;                  // z vertex
    const AliVVertex                    *spdVtx;                //spd z vertez
    
    AliMCEvent                          *fMC;                   //! MC object
    AliStack                            *fStack;                //! stack
    AliAODMCParticle                    *fMCparticle;           //! MC particle
    TClonesArray                        *fMCarray;              //! MC array
    AliAODMCHeader                      *fMCheader;
    
	// Binning and tree size constants
	static constexpr Int_t 				nbins_ept = 33;			 // number of bins for electron pt
	static constexpr Int_t				nbins_eptTPC = 18;       // number of bins for electron pt using TPC
	static constexpr Int_t 				nbins_jetpt = 6;   		 // number of bins for jet pt for substructure
	static constexpr Int_t 				nbins_g = 14;			 // number of bins for angularity
	static constexpr Int_t 				nbins_ptd = 14;			 // number of bins for momentum dispersion
	static constexpr Int_t 				nbins_gsys = 6;			 // number of bins for angularity systematic
	static constexpr Int_t 				nbins_ptdsys = 7;		 // number of bins for momentum dispersion systematic
	static constexpr Int_t 				nbins_ptauxiliary = 59;  // number of bins for auxiliary pt
    static constexpr Int_t 				nbranches = 21;		     // number of branches in output tree
    static constexpr Int_t				ncutsefflowpt = 5;		 // cuts for HFE efficiency calculation for TPC + TOF
    static constexpr Int_t				ncutseffhighpt = 7;		 // cuts for HFE efficiency calculation for EMCal + TPC
    static constexpr Int_t				ncutseffsubs = 2;	     // cuts for HFE efficiency calculation for substructure
    static constexpr Int_t				ncutssemiincl = 6;	     // cuts for constituent semi inclusive jets (leading track pt)
    
    Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted.
    Float_t                             fMinFractionShared;      // only fill histos for jets if shared fraction larger than X
    JetShapeType                        fJetShapeType;           // jet type to be used
    JetShapeSub                         fJetShapeSub;            // jet subtraction to be used
    JetSelectionType                    fJetSelection;           // Jet selection: inclusive/recoil jet
    Float_t                             fShapesVar[nbranches];   // jet shapes used for the tagging
    Float_t                             fPtThreshold;
    Float_t                             fRMatching;
    Int_t                               fSelectedShapes;         //chose set of shapes
    Float_t                             fminpTTrig;              //min - max pT for trigger particle in case of recoil jet
    Float_t                             fmaxpTTrig;
    Float_t                             fangWindowRecoil;        //angular window for btb recoil analysis
    Int_t                               fSemigoodCorrect;        //if==1 we run over semigood runs
    Float_t                             fHolePos;                //position in radians of the bad TPC sector
    Float_t                             fHoleWidth;              //width of the hole in radians
    Bool_t                              fCentSelectOn;           // switch on/off centrality selection
    Float_t                             fCentMin;                // min centrality value
    Float_t                             fCentMax;                // max centrality value
    Bool_t                              fOneConstSelectOn;       // switch on/off one constituent selection
    Int_t                               fDerivSubtrOrder;
    Int_t                               fMCweight;               //0: no weight, 1: enhanced MC, 2: MB MC
    Int_t                               fRunNumber;              // run number
    
    Double_t                            fAssPtCut;               // pt cut for associated electron
    Int_t                               fITSncut;                // ITC number of clusters for tagged electrons
    Int_t                               fAssTPCnCut;             // TPC number of clusters for associated electron
    Int_t                               fTPCnCut;                // TPC number of clusters for tagged electron
    Bool_t                              fAssITSrefitCut;         // ITS refit for associated electron
    Double_t                            fSigmaTOFcut;            // sigma TOF cut |sigma_TOF|< cut
    Double_t                            fSigmaTPCcutLowPt;       // sigma TPC cut for low pt electron identification
    Double_t                            fSigmaTPCcutHighPt;      // sigma TPC cut for high pt electron identification
    Double_t                            fSigmTPCcutExcElec;      // sigma TPC for exclusion of electron
    Double_t                            fDcaXYcut;               // DCA_xy cut
    Double_t                            fDcaZcut;                // DCA_xy cut
    Double_t                            fIMcut;                  // invariant mass cut
    Double_t                            fEtaCut;                 // eta cut of inclusive electrons
    Double_t                            fMinEoPcut;              // minimum value of the E/p cut
    Double_t                            fMaxEoPcut;              // maximum value of the E/p cut
    Double_t                            fM20cut;                 // maximum value of the M20 cut
    Double_t                            fMinM02cut;                 // minimum value of the M02 cut
    Double_t                            fMaxM02cut;                 // maximum value of the M02 cut
    Double_t                            fMinPtTPC;               // minimum pt for the TPC analysis
    Double_t                            fMaxPtTPC;               // maximum pt for the TPC analysis
    Double_t                            fMinPtEMCal;             // minimum pt for the EMCal analysis
    Double_t                            fMaxPtEMCal;             // maximum pt for the EMCal analysis
    Double_t                            fMinPtSemiInclusive;     // minimum pt for constituent for semi-inclusive distributions
    
    TH1F                                *fNeventV0;
    TH1F                                *fNeventT0;
    TH2F                                *fh2ResponseUW;     
    TH2F                                *fh2ResponseW;
    TH2F                                *fPhiJetCorr6;
    TH2F                                *fPhiJetCorr7;
    TH2F                                *fEtaJetCorr6;
    TH2F                                *fEtaJetCorr7;
    TH2F                                *fPtJetCorr;
    TH1F                                *fPtJet;
    TH1F                                *fPtGenJet;
    TH2F                                *fPhiJet;
    TH2F                                *fEtaJet;
    TH2F                                *fEtaPhiJet;
    TH2F                                *fAreaJet;
    TH2F                                *fJetProbDensityDetPart;
    TH2F                                *fJetProbDensityPartDet;
    TH2F                                *fNbOfConstvspT;
    TH2F                                *fnTPCnTOFnocut;
    TH2F                                *fnTPCnocutP;
    TH2F                                *fnTOFnocutP;
    TH2F                                *fnTPCcutP;
    TH2F                                *fnTPCcutPt;
    TH1F                                *fnTPCSigma[nbins_jetpt][nbins_eptTPC];
    TH2F                                *fnULSmLSpairsPerElectron;
    TH2F                                *fInvmassLS[nbins_jetpt];
    TH2F                                *fInvmassULS[nbins_jetpt];
    TH1F                                *fnPartPerJet;
    TH1F                                *fnElecOverPartPerJet;
    TH1F                                *fnInclElecPerJet;
    TH1F                                *fnPhotElecPerJet;
    TH1F                                *fnIncSubPhotElecPerJet;
    TH1F                                *fnTrueElecPerJet;
    TH1F                                *fnTrueHFElecPerJet;
    TH1F                                *fnTruePElecPerJet;
    TH2F                                *fnTrueElecPerJetPt;
    TH2F                                *fnTrueHFElecPerJetPt;
    TH2F                                *fnTruePElecPerJetPt;
    TH1F                                *fPi0PtGen;
    TH1F                                *fPi0PtEnh;
    TH1F                                *fEtaPtGen;
    TH1F                                *fEtaPtEnh;
    TH1F                                *fGenHfePt;
    TH1F                                *fGenPePt;
    TH1F                                *fRecPEAng[nbins_jetpt][nbins_gsys];
    TH1F                                *fTotPEAng[nbins_jetpt][nbins_gsys];
    TH1F                                *fRecPEDisp[nbins_jetpt][nbins_ptdsys];
    TH1F                                *fTotPEDisp[nbins_jetpt][nbins_ptdsys];
    TH1F                                *fULSptAng[nbins_jetpt][nbins_gsys];
    TH1F                                *fLSptAng[nbins_jetpt][nbins_gsys];
    TH1F                                *fULSptDisp[nbins_jetpt][nbins_ptdsys];
    TH1F                                *fLSptDisp[nbins_jetpt][nbins_ptdsys];
    TH2F                                *fPtP;
    TH1F                                *fptJetIE;               // pT of jets containing IE
    TH1F                                *fptJetPE;               // pT of jets containing PE
    TH1F                                *fptJetHFE;              // pT of jets containing HFE
    TH1F                                *fptJetHadron;           // pT of jets not containing electrons
    TH1F                                *fptRecPE;
    TH1F                                *fptTruePE;
    TH1F                                *fptTrueHFEeffTPCTOF[ncutsefflowpt];
    TH3F                                *fptTrueHFEeffTPCTOFang[ncutseffsubs];
    TH3F                                *fptTrueHFEeffTPCTOFdisp[ncutseffsubs];
    TH1F                                *fptTrueHFEeffEMCal[ncutseffhighpt];
    TH3F                                *fptTrueHFEeffEMCalang[ncutseffsubs];
    TH3F                                *fptTrueHFEeffEMCaldisp[ncutseffsubs];
    TH1F                                *fPtTrack;
    TH2F                                *fPhiTrack;
    TH2F                                *fEtaTrack;
    TH2F                                *fEtaPhiTrack;
    TH2F                                *fPhiRecElecTPC;
    TH2F                                *fEtaRecElecTPC;
    TH2F                                *fEtaPhiRecElecTPC;
    TH2F                                *fPhiRecElecEMCal;
    TH2F                                *fEtaRecElecEMCal;
    TH2F                                *fEtaPhiRecElecEMCal;
    TH2F                                *fPhiTrueElec;
    TH2F                                *fEtaTrueElec;
    TH2F                                *fEtaPhiTrueElec;
    TH2F                                *fnEovPelecNoTPCcut;
    TH2F                                *fnEovPelecTPCcut;
    TH2F                                *fnEovPelecEMCalcut;
    TH2F                                *fnEovPelecTPCEMCalcut;
    TH2F                                *fnEovPelecTPCsscut[nbins_jetpt];
    TH2F                                *fnM20backg;
    TH2F                                *fnM02backg;
    TH2F                                *fnEovPbackg;
    TH2F                                *fnEovPbackgEMCalcut;
    TH2F                                *fnClsE;
    TH2F                                *fnM20;
    TH2F                                *fnM02;
    TH2F                                *fnClsTime;
    TH2F                                *fnM20TrueElec;
    TH2F                                *fnM02TrueElec;
    TH3F                                *fnShowerShapeTrueElec;
    TH2F                                *fnEovPTrueElecnocut;
    TH2F                                *fnEovPTrueElecTPCcut;
    TH2F                                *fnEovPTrueElecEMCalcut;
    TH2F                                *fnEovPTrueElecTPCEMCalcut;
    TH2F                                *fnM20TrueBkg;
    TH2F                                *fnM02TrueBkg;
    TH3F                                *fnShowerShapeTrueBkg;
    TH2F                                *fnEovPTrueBkgnocut;
    TH2F                                *fnEovPTrueBkgTPCcut;
    TH2F                                *fnEovPTrueBkgEMCalcut;
    TH2F                                *fnEovPTrueBkgTPCEMCalcut;
    TH2F                                *fAngULS;
    TH2F                                *fAngLS;
    TH2F                                *fAngChargPart;
    TH2F                                *fAngHadron;
    TH2F                                *fAngIncElec;
    TH2F                                *fAngPhotElec;
    TH2F                                *fAngElecFromD;
    TH2F                                *fAngElecFromB;
    TH2F                                *fAngElecFromDFromB;
    TH2F                                *fAngD;
    TH2F                                *fAngB;
    TH2F                                *fAngCharm;
    TH2F                                *fAngBeauty;
    TH2F                                *fAngQuark;
    TH2F                                *fAngGluon;
    TH2F                                *fDispULS;
    TH2F                                *fDispLS;
    TH2F                                *fDispChargPart;
    TH2F                                *fDispHadron;
    TH2F                                *fDispIncElec;
    TH2F                                *fDispPhotElec;
    TH2F                                *fDispElecFromD;
    TH2F                                *fDispElecFromB;
    TH2F                                *fDispElecFromDFromB;
    TH2F                                *fDispD;
    TH2F                                *fDispB;
    TH2F                                *fDispCharm;
    TH2F                                *fDispBeauty;
    TH2F                                *fDispQuark;
    TH2F                                *fDispGluon;

	// (Semi-)inclusive observables and response matrices
	TH1F								*fPtSemiInclJet[ncutssemiincl];
	TH2F								*fAngSemiInclJet[ncutssemiincl];
	TH2F								*fDispSemiInclJet[ncutssemiincl];
	THnSparseF							*fRMPtSemiInclJet[ncutssemiincl];
	THnSparseF						    *fRMAngSemiInclJet[ncutssemiincl];
	THnSparseF						    *fRMDispSemiInclJet[ncutssemiincl];
	// Unweighted response matrices
	THnSparseF							*fRMUWPtSemiInclJet[ncutssemiincl];
	THnSparseF						    *fRMUWAngSemiInclJet[ncutssemiincl];
	THnSparseF						    *fRMUWDispSemiInclJet[ncutssemiincl];

	// No electron methodology
	TH2F								*fptJetNoElectrons;
	TH2F								*fptExtendedJetNoElectrons;
	TH3F								*fAngJetNoElectrons;
	TH3F								*fDispJetNoElectrons;
    
	TTree                               *fTreeObservableTagging;            // Tree with tagging variables subtracted MC or true MC or raw
    

private:
    AliAnalysisTaskEmcalHfeTagging(const AliAnalysisTaskEmcalHfeTagging&);            // not implemented
    AliAnalysisTaskEmcalHfeTagging &operator=(const AliAnalysisTaskEmcalHfeTagging&); // not implemented
    
    ClassDef(AliAnalysisTaskEmcalHfeTagging, 6)
};
#endif

