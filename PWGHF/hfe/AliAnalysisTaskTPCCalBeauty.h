//
//  AliAnalysisTaskTPCCalBeauty.h
//
//
//  Created by Erin Gauger
//
//

#ifndef AliAnalysisTaskTPCCalBeauty_cxx
#define AliAnalysisTaskTPCCalBeauty_cxx

#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
//#include "AliSelectNonHFE.h"
#include "AliAODMCParticle.h"
#include "TF1.h"

class THnSparse;
class AliMultSelection;
class AliAODMCParticle;
class AliAODMCHeader;

class AliAnalysisTaskTPCCalBeauty : public AliAnalysisTaskSE
{
public:
    //two class constructors
    AliAnalysisTaskTPCCalBeauty();
    AliAnalysisTaskTPCCalBeauty(const char *name);
    //class destructor
    virtual         ~AliAnalysisTaskTPCCalBeauty();
    //called once at beginning of runtime
    virtual void    UserCreateOutputObjects();
    //called for each event
    virtual void    UserExec(Option_t *option);
    //called at end of analysis
    virtual void    Terminate(Option_t *option);
    
    Bool_t          GetTenderSwitch() {return fUseTender;};
    void            SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
    
    Bool_t          GetEMCalTriggerEG1() { return fEMCEG1; };
    void            SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; };
    Bool_t          GetEMCalTriggerDG1() { return fDCalDG1; };
    void            SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; };
    
    void            SetCentSelection(Bool_t applyCent) {fApplyCentrality = applyCent;};
    void            SetFillSprs(Bool_t fillElecSprs) {fFlagFillSprs = fillElecSprs;};
    void            SetMC(Bool_t setMCfill) {fFlagFillMCHistos = setMCfill;};
    void            SetSSCut(Double_t maxM20Cut) {fMaxM20Cut = maxM20Cut;};
    void            UseLongAxis(Bool_t applyM02) {fApplyM02Cut = applyM02;}
    void            SetPileUpCut1(Bool_t EnablePileupCut1){fEnablePileupCut1 = EnablePileupCut1;};
    void            SetPileUpCut2(Bool_t EnablePileupRejVZEROTPCout){fEnablePileupRejVZEROTPCout = EnablePileupRejVZEROTPCout;};
    void            SetEoPShift(Bool_t EnableShiftEoP){fShiftEoP = EnableShiftEoP;};
    void            SetTimeCut(Bool_t EnableTimeCut){fApplyTimeCut = EnableTimeCut;};
    void            SetITSLayer(Int_t EnableLayer){fApplyITSLayer = EnableLayer;};
    void            SetTPCnCrossRows(Int_t nTpcCrossRowCut){fNCrossRows = nTpcCrossRowCut;};
    void            SetITSChi2(Int_t itsChi2Cut){fItsChi2 = itsChi2Cut;};
    
    void            SetEoP(Double_t eopCut) {fMinEoPCut = eopCut;};
    void            SetEoPMax(Double_t eopCutMax) {fMaxEoPCut = eopCutMax;};
    void            SetNSig(Double_t nSigCut) {fMinNSigCut = nSigCut;};
    
    void            SetNSigAsso(Double_t nSigAssoCut) {fMinNSigAssoCut = nSigAssoCut;};
    void            SetPtAsso(Double_t ptAssoCut) {fMinPtAssoCut = ptAssoCut;};
    
    void            SetDCABinSize(Double_t dcaBinning) {fDCABinSize = dcaBinning;};
    void            SetStackLoop(Bool_t runStackLoop) {fFlagRunStackLoop = runStackLoop;};
    void            SetTPCClus(Int_t nTPCclusters) {fNclusTPC = nTPCclusters;};
    void            SetDCAzCut(Double_t dcaZCut) {fDCAzCut = dcaZCut;};
    void            SetMassCut(Double_t minM, Double_t maxM) {
                        fMinMass = minM;
                        fMaxMass = maxM;
                    };
    void            SetEtaCut(Double_t etaMin, Double_t etaMax) {
                        fMinEta = etaMin;
                        fMaxEta = etaMax;
                    }
    void            SetAssoDCACut(Double_t assoXYdca, Double_t assoZdca) {
        fAssoDCAxy = assoXYdca;
        fAssoDCAz = assoZdca;
    };
    void            SetAssoTPCClus(Int_t nAssoTPCclusters) {fAssoTPCnCls = nAssoTPCclusters;};
    void            SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
    void            SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
    void            SetTrkMatch(Double_t maxTrkMatch) {fTrkMatch = maxTrkMatch;};
    
    void            SetCentralitySelection(Double_t centMin, Double_t centMax){fCentralityMin = centMin; fCentralityMax = centMax;};
    Double_t        CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    
    Bool_t          GetNMCPartProduced();
    void            GetPi0EtaWeight(THnSparse *SparseWeight);
    
    void            GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
    void            FindMother(AliAODMCParticle* part, Int_t &fpidSort, Bool_t &kEmbEta, Bool_t &kEmbPi0, Bool_t &kHijing, Double_t &momPt, Double_t &momGamma, Double_t &momTime);
    void            InvMassCheckData(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Double_t fWeight);
    void            InvMassCheckMC(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Bool_t kHijing, Bool_t kEmbEta, Bool_t kEmbPi0, Bool_t &kFlagReco, Double_t fWeight, Int_t fpidSort, Double_t prodRadius);
    void            FillULSSparse(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Double_t fWeight,Double_t prodRadius,Int_t pidSort);
    
    void            SetHadronEoPCut(Bool_t hadronEopCut) {fApplyHadEoPCut = hadronEopCut;};
    void            SetVtxZCut(Double_t zVertexCut) {fVtxZCut = zVertexCut;};
    void            SetDCAxyCut(Double_t dcaXYCut) {fDCAxyCut = dcaXYCut;};
    //void            SetBmesonTauWeight(TF2 *BPlus, TF2 *B0, TF2 *Bs);
    void            SetTauWeight(Bool_t useDLWeight) {fUseTauWeight = useDLWeight;};
    
private:
    AliAODEvent         *fAOD;           //! input event
    AliAODMCHeader      *fMCHeader;      //!
    TClonesArray        *fMCarray;       //! MC array
    AliAODMCParticle    *fMCparticle;    //! MC particle
    AliPIDResponse      *fpidResponse;   //! pid response
    TList               *fOutputList;    //! output list
    
    AliMultSelection    *fMultSelection; //! need to get centrality
    Double_t            fCentrality;     //!
    Double_t            fCentralityMin;  // set event centrality min
    Double_t            fCentralityMax;  // set event centrality max
    
    Bool_t              fEMCEG1;         // EMCal Threshold EG1
    Bool_t              fDCalDG1;        // DCal Threshold DG1
    //Bool_t              fFlagApplySSCut; //switch to turn on SS cut
    Double_t            fMaxM20Cut;      // set eID M20 cut
    Bool_t              fApplyM02Cut;    // apply M02 instead of M20 cut
    Bool_t              fEnablePileupCut1; //add additional pile-up cuts
    Bool_t              fEnablePileupRejVZEROTPCout; //add additional pile-up cuts
    Bool_t              fShiftEoP; //apply E/p shift in MC
    Bool_t              fApplyTimeCut; //apply cluster timing cuts
    Int_t              fApplyITSLayer; //0=kAny, 1=kFirst, 2=kBoth
    Int_t              fNCrossRows; //set N of TPC crossed rows
    Int_t              fItsChi2; //set max ITS chi-2. if negative, does not apply any cut
    Double_t            fMinEoPCut;      // set min eID E/p cut
    Double_t            fMaxEoPCut;      // set max eID E/p cut
    Double_t            fMinNSigCut;     // set eID nSig cut
    
    Double_t            fMinNSigAssoCut; // set nSig cut for asso track
    Double_t            fMinPtAssoCut;   // set min pt for asso track
    
    Double_t            fDCABinSize;     // set DCA bin size for systematics
    Bool_t              fApplyCentrality; //apply centrality cut
    Bool_t              fFlagFillSprs; //switch to fill electron eid sparse
    Bool_t              fFlagFillMCHistos; // switch to fill histos that require MC pid
    Bool_t              fFlagRunStackLoop; //switch to run stack loop to get D & B meson species info
    Int_t               fNclusTPC;       // set number of TPC clusters
    Double_t            fDCAzCut;        // set DCA z cut
    Double_t            fMinMass;        // set min inv mass
    Double_t            fMaxMass;        // set max inv mass
    Double_t            fMinEta;        // set min inv mass
    Double_t            fMaxEta;        // set max inv mass
    
    Double_t            fAssoDCAxy;     //asso track dcaXY cut
    Double_t            fAssoDCAz;      //asso track dcaZ cut
    Int_t               fAssoTPCnCls;   //TPC nCls cut for asso track
    Bool_t              fFlagClsTypeEMC; // switch to select EMC clusters
    Bool_t              fFlagClsTypeDCAL;// switch to select DCAL clusters
    Double_t            fTrkMatch; //set distance to cluster
    
    Bool_t              fUseTender;      // switch to add tender
    Bool_t              fApplyHadEoPCut; // switch to apply E/p cut to hadrons
    Double_t            fVtxZCut;        // z vertex cut (in cm)
    Double_t            fDCAxyCut;       // max DCAxy cut
    Bool_t              fFlagULS;        // flag ULS
    Bool_t              fFlagLS;         // flag LS
    
    TH1F                *fNevents;       //! no of events
    TH1F                *fVtX;           //! vertex x
    TH1F                *fVtY;           //! vertex y
    TH1F                *fVtZ;           //! vertex z
    
    TH1F                *fTrkPtB4TC;     //! track pT before track cuts
    
    TH2F                *fDCAxyz;       //! DCAz vs. DCAxy
    
    TH1F                *fTrkPt;         //! track pT
    TH1F                *fTrkP;          //! track p
    TH1F                *fTrkClsPhi;     //! track and cluster delta phi
    TH1F                *fTrkClsEta;     //! track and cluster delta eta
    TH1F                *fClsPhi;         //! cluster phi
    TH1F                *fClsEta;         //! cluster eta
    //TH1F                *fClsEamDCal;    //! cluster energy after matching to DCal
    //TH1F                *fClsEamEMCal;   //! cluster energy after matching to EMCal
    //TH1F                *fClsEAll;   //! cluster energy of all track-matched clusters
    TH1F                *fClsE;   //! cluster energy of EMCal/DCal depending on subfolder
    TH1F                *fClsEnoTimeCut; //! cluster energy after matching w/o cluster time cut
    //TH1F                *fClsEamElecEMC;   //! cluster energy of e- after matching to EMCal
    //TH1F                *fClsEamElecDC;   //! cluster energy of e- after matching to EMCal
    TH1F                *fTrkPhi;        //! track phi after track matching
    TH1F                *fTrkEta;        //! track eta after track matching
    TH1F                *fdEdx;          //! track dEdx
    TH2F                *fnSigma;          //! track dEdx
    TH2F                *fnSigmaAftTrkMatch;  //! track dEdx after matching cal
    TH1F                *fCentCheck;     //! event centrality
    TH1F                *fTrigCheck;     //! checking trigger used
    TH1F                *fITSLayerCheck;  //! checking ITS layer requirement
    TH2F                *fEMCTrkMatch;   //! plots distance of cluster from closest track
    
    THnSparse           *fBasicSprs;    //! Sparse with pT, DCA, phi, eta, charge, electron bool, trk match bool
    Double_t            *fvalueBasic;   //!Electron info
    
    TH1F                *fInvmassLS;     //! Plots LS mass dist
    TH1F                *fInvmassULS;    //! Plots ULS mass dist
    
    //TH1F                *fInvmassLSWeightEnhEta;   //! Plots LS mass dist
    //TH1F                *fInvmassULSWeightEnhEta;  //!
    //TH1F                *fInvmassLSWeightEnhPi0;   //!
    //TH1F                *fInvmassULSWeightEnhPi0;  //!
    //TH1F                *fInvmassLSHijingEta;      //!
    //TH1F                *fInvmassULSHijingEta;     //!
    //TH1F                *fInvmassLSHijingPi0;      //!
    //TH1F                *fInvmassULSHijingPi0;     //!
    //TH1F                *fInvmassLSHijingPhoton;   //!
    //TH1F                *fInvmassULSHijingPhoton;  //!
    //TH1F                *fInvmassLSEnhPhoton;         //!
    //TH1F                *fInvmassULSEnhPhoton;        //!
    
    TH2F                *fULSdcaBelow;   //! ULS electron DCA vs. pT, m<0.1
    TH2F                *fLSdcaBelow;    //! LS electron DCA vs. pT, m<0.1
    TH2F                *fULSdcaBelowWeight;   //! ULS electron DCA vs. pT, m<0.1
    TH2F                *fLSdcaBelowWeight;    //! LS electron DCA vs. pT, m<0.1
    
    TH1F                *fLSWeightEnhEta;     //! LS for Weighted enhanced eta
    TH1F                *fULSWeightEnhEta; //! ULS for Weighted enhanced eta
    TH1F                *fLSWeightEnhPi0;  //! LS for Weighted enhanced pi0
    TH1F                *fULSWeightEnhPi0; //! ULS for Weighted enhanced pi0
    TH1F                *fLSHijingEta;     //! LS for hijing eta
    TH1F                *fULSHijingEta;    //! ULS for hijing eta
    TH1F                *fLSHijingPi0;     //! LS for hijing pi0
    TH1F                *fULSHijingPi0;    //! ULS for hijing pi0
    TH1F                *fLSHijingPhoton;  //! LS for hijing photon
    TH1F                *fULSHijingPhoton; //! ULS for hijing photon
    TH1F                *fLSEnhPhoton;     //! LS for all photon e
    TH1F                *fULSEnhPhoton;    //! ULS for all photon e
    
    //TH2F                *fPhotonicDCA;   //! Photonic DCA using MC PID
    TH2F                *fInclElecDCA;   //! Inclusive electron DCA vs. pT
    TH2F                *fnSigaftEoPCut;   //! DCA after eID E/p cut
    //TH2F                *fnSigaftSysEoPCut;   //! DCA after eID systematic E/p cut
    TH2F                *fnSigaftM20EoPCut;   //! DCA after eID M20+E/p cut
    //TH2F                *fnSigaftSysM20EoPCut;   //! DCA after Sys eID M20 + E/p cut
    //TH2F                *fInclElecDCAnoSign;   //! Inclusive electron DCA vs. pT, no sign
    //TH2F                *fElecEoPnoSig;  //! Elec EoP w/o sigma cut
    
    //TH2F                *fInclElecEoPnoShift;   //! Inclusive electron EoP vs. pT w/o E/p shift
    //TH2F                *fHadronEoPnoShift;     //! Hadron EoP vs. pT w/o E/p shift
    TH2F                *fInclElecEoP;   //! Inclusive electron EoP vs. pT
    TH2F                *fInclElecEoPNoM20;   //! Inclusive electron EoP vs. pT
    //TH2F                *fTPCElecEoP;   //! EoP vs. pT, -0.1<nsig<3 cut
    TH2F                *fHadronEoP;     //! Hadron EoP vs. pT
    TH2F                *fHadronEoPNoM20;     //! Hadron EoP vs. pT
    TH2F                *fHadronDCA;     //! Hadron DCA vs. pT
    //TH2F                *fHadronCamDCAHij;     //! Hadron DCA vs. pT, no sign, no E/p cut
    //TH2F                *fHadronCamDCA;     //! Hadron DCA vs. pT, no sign, no E/p cut
    
    TF1                 *fPi0Weight;    //! Function to weight enhanced pi0
    TF1                 *fEtaWeight;    //! Function to weight enhanced eta
    //TF1                 *fPi0EtaWeight; //! Function to weight enhanced eta+pi0
    //TH1F                *fDWeight; //!
    TH1F                *fDWeightNew; //!
    TH1F                *fDWeightVar1; //!
    TH1F                *fDWeightVar2; //!
    TH1F                *fDPlusWeightVar1; //!
    TH1F                *fDsWeightVar1; //!
    TH1F                *fLcWeightVar1; //!
    TH1F                *fLcWeightVar2; //!
    //TH1F                *fBWeight; //!
    TH1F                *fBWeightNew; //!
    TH1F                *fBWeightVar1; //!
    TH1F                *fBWeightVar2; //!
    TF1                 *fBPlusTauWeight; //!
    TF1                 *fB0TauWeight; //!
    TF1                 *fBsTauWeight; //!
    TF1                 *fDPlusTauWeight; //!
    TF1                 *fD0TauWeight; //!
    TF1                 *fDsTauWeight; //!
    Bool_t              fUseTauWeight; //Switch to apply delay length weight
    
    
    Double_t            fWeight;        //!
    
    //TH2F                *fPi0DCA;           //! Pi0 DCA vs. pT
    //TH2F                *fEtaDCA;           //! Eta DCA vs. pT
    
    TH2F                *fEnhEtaDCA;           //!
    TH1F                *fEnhEtaWeightedPt;    //!
    TH2F                *fEnhPi0DCA;           //!
    TH1F                *fEnhPi0WeightedPt;    //!
    TH2F                *fEtaHijingDCA;        //!
    TH1F                *fEtaHijingPt;         //!
    TH2F                *fPi0HijingDCA;        //!
    TH1F                *fPi0HijingPt;         //!
    
    TH2F                *fPhotonHijingDCA;     //!
    TH1F                *fPhotonHijingPt;      //!
    TH2F                *fEnhPhotonDCA;        //!
    TH1F                *fEnhPhotonWeightedPt; //!
    
    TH2F                *fPhotonHijingTagDCA;     //!
    TH2F                *fEnhPhotonTagDCA;        //!
    
    TH3F                *fComboNumWeight;       //!
    TH3F                *fComboNumNoWeight;     //!
    TH3F                *fComboDenomWeight;     //!
    TH3F                *fComboDenomNoWeight;   //!
    
    //TH1F                *fDMesonPDG; //! plots abs(pdg) of D mesons in the stack
    TH1F                *fD0MesonPt;  //!
    TH1F                *fD0MesonFromDStarPt; //!
    TH1F                *fDPlusMesonPt; //!
    TH1F                *fDsMesonPt;//!
    TH1F                *fDStarMesonPt;//!
    TH1F                *fAllDMesonPt; //!
    
    TH1F                *fLambdaCPt; //!
    TH1F                *fD0MesonPtWeight;  //!
    TH1F                *fLambdaCPtWeight; //!
    TH1F                *fDPlusMesonPtWeight; //!
    TH1F                *fDsMesonPtWeight; //!
    TH1F                *fEtaCPt; //!
    TH1F                *fCBaryonPt; //!
    
    TH1F                *fBMesonPt; //!
    //TH1F                *fBMesonPtATLAS; //!
    //TH1F                *fBPlusPtATLAS; //!
    //TH1F                *fBMesonPtCMS; //!
    //TH1F                *fBPlusPtCMS; //!
    //TH1F                *fBMesonPtLHCb; //!
    //TH1F                *fBPlusPtLHCb; //!
    TH1F                *fBBaryonPt; //!
    
    TH1F                *fBMesonElecPt; //!
    TH1F                *fBBaryonElecPt; //!
    
    //TH2F    *fPromptD0DCAWeight; //!
    //TH2F    *fD0FromDStarDCAWeight; //!
    //TH2F    *fPromptD0DCANoWeight; //!
    //TH2F    *fD0FromDStarDCANoWeight; //!
    
    Int_t               fNtotMCpart;     //! N of total MC particles produced by generator
    Int_t               fNpureMC;        //! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0;      //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta;      //! N > fNembMCeta = particles from eta generator
    
    THnSparse           *fSprsPi0EtaWeightCal;  //! Sparse for pi0,eta weight calc
    THnSparse           *fSprsTemplatesNoWeight;  //! Sparse for templates
    THnSparse           *fSprsTemplatesWeight;  //! Sparse for templates
    THnSparse           *fSprsTemplatesWeightVar1;  //! Sparse for templates
    THnSparse           *fSprsTemplatesWeightVar2;  //! Sparse for templates
    THnSparse           *fSprsClosureTest;  //! Sparse for templates
    THnSparse           *fSprsClosureTestWeight;  //! Sparse for templates
    
    THnSparse           *fSprsULSdca; //! Sparse with ULS info
    THnSparse           *fSprsULSdcaWeight; //! Sparse with ULS info
    THnSparse           *fSprsLSdca; //! Sparse with LS info
    THnSparse           *fSprsLSdcaWeight; //! Sparse with LS info
    
    //TH2F                *fDTemplateWeight; //!
    //TH2F                *fDTemplateNoWeight; //!
    //TH2F                *fDTemplateWeightNew; //!
    //TH2F                *fDTemplateWeightVar1; //!
    //TH2F                *fDTemplateWeightVar2; //!
    
    //TH2F                *fBTemplateWeight; //!
    //TH2F                *fBTemplateNoWeight; //!
    //TH2F                *fBTemplateWeightNew; //!
    //TH2F                *fBTemplateWeightVar1; //!
    //TH2F                *fBTemplateWeightVar2; //!
    
    TH1F                *fAllElecStack; //!
    TH1F                *fHFElecStack; //!
    TH1F                *fBElecStack; //!
    
    TH1F                *fAllElecStackDiffPID; //!
    TH1F                *fDElecStackDiffPID; //!
    TH1F                *fBElecStackDiffPID; //!
    
    TH1F                *fElecTPCTrk; //!
    TH1F                *fHFElecTPCTrk; //!
    TH1F                *fBElecTPCTrk; //!
    
    TH1F                *fElecAftTrkCuts; //!
    TH1F                *fHFElecAftTrkCuts; //!
    TH1F                *fBElecAftTrkCuts; //!
    
    TH1F                *fElecAftLooseTrkCuts; //!
    TH1F                *fHFElecAftLooseTrkCuts; //!
    TH1F                *fBElecAftLooseTrkCuts; //!
    
    /*TH1F                *fElecAftLooseTrkCutsDiffPID; //!
    TH1F                *fDElecAftLooseTrkCutsDiffPID; //!
    TH1F                *fBElecAftLooseTrkCutsDiffPID; //!*/
    
    TH1F                *fElecAftTrkMatch; //!
    TH1F                *fHFElecAftTrkMatch; //!
    TH1F                *fBElecAftTrkMatch; //!
    
    TH1F                *fElecAftTPCeID; //!
    TH1F                *fHFElecAftTPCeID; //!
    TH1F                *fBElecAftTPCeID; //!
    
    TH1F                *fElecAftEMCeID; //!
    TH1F                *fHFElecAftEMCeID; //!
    TH1F                *fBElecAftEMCeID; //!
    
    TH1F                *fElecAftEoP; //!
    TH1F                *fHFElecAftEoP; //!
    TH1F                *fBElecAftEoP; //!
    
    THnSparse           *fElectronSprs;  //! Sparse with electron cut parameters
    //Double_t            *fvalueElectron;//!Electron info
    
    //Double_t            *fvalueElectron; //! Electron info
    
    AliAnalysisTaskTPCCalBeauty(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    AliAnalysisTaskTPCCalBeauty& operator=(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    ClassDef(AliAnalysisTaskTPCCalBeauty, 1);
};

#endif
