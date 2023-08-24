
#ifndef AliAnalysisTaskEHCorrel_h
#define AliAnalysisTaskEHCorrel_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////
//                                                                    //
//Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation in Run 2//
//                                                                    //
//  Author: Deepa Thomas (University of Texas at Austin)              //
//          Ravindra Singh (IIT Indore)                               //
//          Amanda Flores (University of Texas at Austin)             //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliEventPoolManager;
class AliMultSelection;
class AliAnalysisUtils;
class AliGenEventHeader;

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"

class AliAnalysisTaskEHCorrel : public AliAnalysisTaskSE {
  public:

    enum EnhanceSigOrNot {kMB,kEnhance};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    AliAnalysisTaskEHCorrel();
    AliAnalysisTaskEHCorrel(const char *name);
    virtual ~AliAnalysisTaskEHCorrel();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    Bool_t  PassEventSelect(AliVEvent *fVevent);
    Bool_t  PassAddtionalPileUpCuts();
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    void    CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    Bool_t  PassTrackCuts(AliAODTrack *atrack);
    void    GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
    Bool_t  PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack);
    Bool_t  PassHadronCuts(AliAODTrack *HadTrack);
    void    HadronInfo();
    //void    HadronInfo(Int_t itrack);
    //void    ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl, TH2F *HisDphi);
    void    ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl);
    void    ElectronHadCorrelNoPartner(Int_t itrack, Int_t jtrack, AliVTrack *track, THnSparse *SparseEHCorrlNoPartner);
    void    EMCalClusterInfo();
    void    SelectNonHFElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS);
    TObjArray* CloneAndReduceTrackList();
    void    MixedEvent(AliVTrack *track, THnSparse *SparseMixEHCorrl);

    void    SetCentralitySelection(Double_t centMin, Double_t centMax) {fCentralityMin = centMin; fCentralityMax = centMax;};
    void    SetEleEtaCuts(Double_t Min, Double_t Max) {fEtaCutEleMin = Min; fEtaCutEleMax = Max;};
    void    SetMinTPCNCrossRElec(Int_t MinNCrossRE) {fTPCNCrossRElec = MinNCrossRE;};
    void    SetMinRatioTPCNCrossRElec(Double_t MinRatioNCrossRE) {fRatioTPCNCrossRElec = MinRatioNCrossRE;};
    void    SetMinITSNClsElec(Int_t MinNClsE) {fITSNClsElec = MinNClsE;};
    void    SetTPCnsigCut(Double_t nsigMin, Double_t nsigMax) {fTPCnSigmaMin = nsigMin; fTPCnSigmaMax= nsigMax;};
    void    SetM02Cut(Double_t m02Min, Double_t m02Max) {fM02Min = m02Min; fM02Max = m02Max;};
    void    SetM20Cut(Double_t m20Min, Double_t m20Max) {fM20Min = m20Min; fM20Max = m20Max;};
    void    SetEovPCut(Double_t eovpMin, Double_t eovpMax) {fEovPMin = eovpMin; fEovPMax = eovpMax;};
    void    SetHadronCutCase(Int_t hadCutCase) {fHadCutCase = hadCutCase;};
    void    SetTriggerElePtCut(Bool_t trigElePtcut) {fTrigElePtCut = trigElePtcut;};
    void    GetVtxZCentralityBin();
    Bool_t  GetTenderSwitch() {return fUseTender;};
    Double_t GetElecEffi(AliVTrack *track);

    void    SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};

    void    SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
    void    SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};

    void    SetPartnerEleMinTPCNCls(Int_t MinNClsPE) {fTPCNClsPartnerE = MinNClsPE;};
    void    SetPartnerEleMinPt(Double_t PtPE) {fPartElePt = PtPE;};
    void    SetInvmassCut(Double_t invmasscut) {fInvmassCut = invmasscut;};

    void    SetHadMinTPCNCrossR(Int_t MinRatioNCrossRHad) {fTPCNCrossRHad = MinRatioNCrossRHad;};
    void    SetHadMinRatioTPCNCrossR(Double_t MinNCrossRHad) {fRatioTPCNCrossRHad = MinNCrossRHad;};
    void    SetHadSPDkAny(Bool_t HadSPDkAny) {fFlagHadSPDkAny = HadSPDkAny;};
    void    SetHadLargeITSNCls(Bool_t HadLargITSNCls) {fFlagHadITSNCls = HadLargITSNCls;};

    void    SetHadFiducialCut(Bool_t HadFiducialCut) {fFlagHadFiducialCut = HadFiducialCut;};
    void    SetHadPosEtaOnly(Bool_t HadPosEtaOnly) {fFlagHadPosEtaOnly = HadPosEtaOnly;};
    void    SetHadNegEtaOnly(Bool_t HadNegEtaOnly) {fFlagHadNegEtaOnly = HadNegEtaOnly;};
    void    SetHadEtaCuts(Double_t Min, Double_t Max) {fEtaCutHadMin = Min; fEtaCutHadMax = Max;};

    void    SwitchMECorrec(Bool_t FillMECorr){fFlagFillMECorr = FillMECorr;};
    void    SetMEBinChange(Bool_t MEBinChange) {fFlagMEBinChange = MEBinChange;};
    void    SetElecSPDkFirst(Bool_t EleSPDkFirst) {fFlagEleSPDkFirst = EleSPDkFirst;};

    void    SetEMCClsTimeCut(Bool_t EMCtimeCut) {fEMCClsTimeCut = EMCtimeCut;};

    void    SetAdditionalPileUpCuts(Bool_t addpilupcuts) {fApplyAddPileUpCuts = addpilupcuts;};

    void    SetElecEffi(Bool_t applyElectronEffi){fApplyElectronEffi = applyElectronEffi;};
    void    SetElectronEffiMap(TH1 *HistElecEffi) {fHistElecEffi = (TH1D*)HistElecEffi->Clone();}

    void    IsPbPb(Bool_t isPbPb) {fIsPbPb = isPbPb;};
    void    Ispp(Bool_t ispp) {fIspp = ispp;};
    void    IspPb(Bool_t ispPb) {fIspPb = ispPb;};
    void    IsPASS2weight(Bool_t pPbpass2weight) {fpPbPASS2weight = pPbpass2weight;};
    void    IsMC(Bool_t isMC) {fIsMC = isMC;};

    void    IsClustersOn(Bool_t fSwitch) {fEMCalClusOn = fSwitch;};
    void    IsHadInfoOn(Bool_t fSwitch) {fHadronInfoOn = fSwitch;};
    void    IsTracksOn(Bool_t fSwitch) {fTracksOn = fSwitch;};

    void    SwitchFillEHCorrel(Bool_t fSwitch){fFillEHCorrel = fSwitch;};
    void    SetNDeltaPhiBins(Int_t nbins){fNDelPhiBins = nbins;};

    void    GetHadronTrackingEfficiency();
    void    SwitchHadTrackEffi(Bool_t fSwitch) {fCalcHadronTrackEffi = fSwitch;};

    void    SetNonHFEEffi(Bool_t fSwitch) {fCalculateNonHFEEffi = fSwitch;};
    void    SetWeightCal(Bool_t fSwitch) {fCalPi0EtaWeight = fSwitch;};
    void    GetPi0EtaWeight(THnSparse *SparseWeight);
    void    GetPi0EtaWeightPbPb2018(THnSparse *SparseWeight);
    Bool_t  GetNonHFEEffiDenom(AliVTrack *track);
    Bool_t  GetNonHFEEffiDenomPbPb2018(AliVTrack *track);
    Bool_t  GetNonHFEEffiRecoTag(AliVTrack *track);
    Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
    Int_t   GetPi0EtaType(AliAODMCParticle *part);
    
    void    RemovePileUpInMCGen(Bool_t fSwitch) {fRemovePileUpinMCGen = fSwitch;};
    
    void    SetTrackMatchPar(Double_t deltaEta, Double_t deltaPhi){fDeltaEta = deltaEta; fDeltaPhi = deltaPhi;};

  private:
    Bool_t GetNMCPartProduced();
    Bool_t GetNMCPartProducedPbPb2018();

    AliVEvent           *fVevent;//!V event object
    AliAODEvent         *fAOD;//!AOD object
    const AliVVertex    *fpVtx; //!
    AliPIDResponse      *fpidResponse; //!pid response
    AliMultSelection    *fMultSelection;//!

    TH1D                *fHistElecEffi;//

    TClonesArray        *fTracks_tender;//Tender tracks
    TClonesArray        *fCaloClusters_tender;//Tender cluster
    Bool_t              fApplyAddPileUpCuts;//
    Bool_t              fUseTender;// switch to add tender
    Double_t            fCentrality;//!
    Double_t            fCentralityMin;//
    Double_t            fCentralityMax;//
    Double_t            fMultiplicity;//!
    Bool_t              fEMCEG1;//
    Bool_t              fEMCEG2;//
    Bool_t              fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t              fFlagClsTypeDCAL;//switch to select DCAL clusters
    Int_t               fTPCNCrossRElec;// track TPC NClusters
    Double_t            fRatioTPCNCrossRElec;//
    Bool_t              fFlagEleSPDkFirst;//
    Double_t            fEtaCutEleMin;// Electron track Eta cut min
    Double_t            fEtaCutEleMax;// Electron track Eta cut max
    Double_t            fDeltaEta;//
    Double_t            fDeltaPhi;//
    Double_t            fTPCnSigma;//!
    Double_t            fTPCnSigmaMin;//
    Double_t            fTPCnSigmaMax;//
    Double_t            fM02Min;//
    Double_t            fM02Max;//
    Double_t            fM20Min;//
    Double_t            fM20Max;//
    Double_t            fEovPMin;//
    Double_t            fEovPMax;//
    Int_t               fTPCNCrossRHad;// Had track TPC NClusters
    Double_t            fRatioTPCNCrossRHad;//
    Double_t            fEtaCutHadMin;// Had track Eta cut min
    Double_t            fEtaCutHadMax;// Had track Eta cut max
    Int_t               fITSNClsElec;//
    Int_t               fTPCNClsPartnerE;//
    Double_t            fPartElePt;//
    Double_t            fInvmassCut;//
    Bool_t              fFlagHadSPDkAny;//
    Bool_t              fFlagHadITSNCls;//
    Bool_t              fFlagHadFiducialCut;//
    Bool_t              fFlagHadPosEtaOnly;//
    Bool_t              fFlagHadNegEtaOnly;//
    Double_t            fTPCnSigmaHadMin;//
    Double_t            fTPCnSigmaHadMax;//
    Int_t               fHadCutCase;//
    AliEventPoolManager *fPoolMgr;//!
    Bool_t              fTrigElePtCut;//
    Int_t               fNEle;//!
    Double_t            fVtxZBin;//!
    Double_t            fCentBin;//!
    Bool_t              fFlagFillMECorr;// Switch for doing the ME correction
    Bool_t              fFlagMEBinChange;//
    Bool_t              fIsPbPb;//
    Bool_t              fIspp;//
    Bool_t              fIspPb;//
    Bool_t              fpPbPASS2weight;//
    Bool_t              fEMCClsTimeCut;//
    TClonesArray        *fMCarray;//!
    AliAODMCHeader      *fMCHeader;//!
    Bool_t              fIsMC;// Is MC
    Bool_t              fEMCalClusOn;//
    Bool_t              fHadronInfoOn;//
    Bool_t              fTracksOn;//
    Bool_t              fApplyElectronEffi;//
    Double_t            fEffi;//!
    Double_t            fWeight;//!
    Bool_t              fCalcHadronTrackEffi;//
    Bool_t              fFillEHCorrel;//
    Bool_t              fRemovePileUpinMCGen;//

    //Non-HFE
    Bool_t              fCalculateNonHFEEffi;//
    Bool_t              fCalPi0EtaWeight;//
    Bool_t              fIsFrmEmbPi0;//!
    Bool_t              fIsFrmEmbEta;//!
    Int_t               ftype;//!
    Double_t            fWeightPi0;//!
    Double_t            fWeightEta;//!
    Int_t                fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    Int_t               fNembMCpileup; //! N > fNembMCpileup = PileupNoGenTrig Cocktail Header
    
    TF1                 *fFuncPtDepEta;//!
    TF1                 *fFuncPtDepPhi;//!

    Int_t               fNDelPhiBins; //number of bins to be used for deltaphi distribution
    TList               *fOutputList;//!output list
    TH1F                *fNevents;//!no of events
    TH1F                *fVtxZ;//!
    TH1F                *fVtxX;//!
    TH1F                *fVtxY;//!
    TH1F                *fCentralityNoPass;//!
    TH1F                *fCentralityPass;//!
    TH1F                *fMultiplicityNoPass;//!
    TH1F                *fMultiplicityPass;//!
    TH2F                *fCentMultiplicityNoPass;//!
    TH2F                *fCentMultiplicityPass;//!
    TH1F                *fHistClustE;//!
    TH2F                *fEMCClsEtaPhi;//!
    TH2F                *fHistoNCells;//!
    TH2F                *fHistoTimeEMC;//!

    TH1F                *fNegTrkIDPt;//!
    TH1F                *fTrkPt;//!
    TH1F                *fTrketa;//!
    TH1F                *fTrkphi;//!
    TH2F                *fdEdx;//!
    TH2F                *fTPCnsig;//!
    TH1F                *fTrkNClsF;//!
    TH1F                *fTrkTPCNCrossRows;//!
    TH1F                *fTrkRatCrossRowNclus;//!
    TH1F                *fHistPtMatch;//!
    TH2F                *fEMCTrkMatch;//!
    TH1F                *fEMCTrkPt;//!
    TH1F                *fEMCTrketa;//!
    TH1F                *fEMCTrkphi;//!
    TH2F                *fEMCTPCnsig;//!
    TH1F                *fClsEAftMatch;//!
    TH2F                *fClsEtaPhiAftMatch;//!
    TH2F                *fHistNsigEop;//!
    TH2F                *fM20EovP;//!
    TH2F                *fM02EovP;//!
    TH2F                *fHistEop;//!
    TH2F                *fM20;//!
    TH2F                *fM02;//!
    TH2F                *fHistEop_AftEID;//!
    TH1F                *fInclsElecPt;//!
    TH1F                *fNElecInEvt;//!
    TH2F                *fHadEop;//!
    TH1F                *fHadPt_AftEID;//!
    TH1F                *fULSElecPt;//!
    TH1F                *fLSElecPt;//!
    TH1F                *fTagULSElecPt;//!
    TH1F                *fTagLSElecPt;//!

    TH2F                *fHadronPhiPt;//!
    TH1F                *fHadronPhi;//!
    TH1F                *fHadronPhiTPChalf;//!
    TH1F                *fHadronPt;//!

    TH1F                *fInvmassLS;//!
    TH1F                *fInvmassULS;//!
    TH2F                *fInvmassLSPt;//!
    TH2F                *fInvmassULSPt;//!

    TH1F                *fNoMixedEvents;//!
    TH2F                *fMixStatCent;//!
    TH2F                *fMixStatVtxZ;//!
    TH2F                *fMixStatCentVtxz;//!

    //TH2F                *fHisHadDphi;//!
    //TH2F                *fHisIncEDphi;//!
    //TH2F                *fHisLSDphi;//!
    //TH2F                *fHisULSDphi;//!

    TH1F                *fTrackPhyPrimAll;//!
    TH1F                *fTrackEffiDenomPt;//!
    TH1F                *fTrackPhyPrimAllPDG;//
    TH1F                *fTrackEffiNumHadTrkPt;//!
    TH1F                *fTrackPhyPrimPDGCut;//!

    TH1F                *fRealInclsElecPt;//!
    TH1F                *fNonHFeTrkPt;//!
    TH1F                *fEtaeEmbWeightTrkPt;//!
    TH1F                *fPi0eEmbWeightTrkPt;//!
    TH1F                *fNonHFeEmbWeightTrkPt;//!
    TH1F                *fNonHFeEmbTrkPt;//!
    TF1                 *fPi0Weight;//!
    TF1                 *fEtaWeight;//!
    TH1F                *fRecoNonHFeTrkPt;//!
    TH1F                *fRecoNonHFeEmbTrkPt;//!
    TH1F                *fRecoNonHFeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0eEmbWeightTrkPt;//!
    TH1F                *fRecoEtaeEmbWeightTrkPt;//!

    THnSparse           *fSprsAllHadHCorrl;//!
    THnSparse           *fSprsMixAllHadHCorrl;//!
    THnSparse           *fSprsHadHCorrl;//!
    THnSparse           *fSprsInclusiveEHCorrl;//!
    THnSparse           *fSprsLSEHCorrl;//!
    THnSparse           *fSprsULSEHCorrl;//!
    THnSparse           *fSprsLSNoPartnerEHCorrl;//!
    THnSparse           *fSprsULSNoPartnerEHCorrl;//!
    THnSparse           *fSprsTagULSEHCorrl;//!
    THnSparse           *fSprsTagLSEHCorrl;//!
    THnSparse           *fSprsMixHadHCorrl;//!
    THnSparse           *fSprsMixInclusiveEHCorrl;//!
    THnSparse           *fSprsMixLSEHCorrl;//!
    THnSparse           *fSprsMixULSEHCorrl;//!
    THnSparse           *fSprsMixTagULSEHCorrl;//!
    THnSparse           *fSprsMixTagLSEHCorrl;//!
    THnSparse           *fSprsPi0EtaWeightCal;//! weight cal

    AliAnalysisTaskEHCorrel(const AliAnalysisTaskEHCorrel&); // not implemented
    AliAnalysisTaskEHCorrel& operator=(const AliAnalysisTaskEHCorrel&); // not implemented

    ClassDef(AliAnalysisTaskEHCorrel, 2); //!example of analysis
};

class AliehDPhiBasicParticle : public AliVParticle
{
  public:

    AliehDPhiBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliehDPhiBasicParticle() {}

    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

    virtual Short_t Charge()      const { return fCharge; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge

    ClassDef( AliehDPhiBasicParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
#endif
