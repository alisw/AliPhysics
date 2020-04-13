#ifndef AliAnalysisTaskHFETPCTOFMultiplicity_H
#define AliAnalysisTaskHFETPCTOFMultiplicity_H
class TH1F;
class TH2F;
class TLorentzVector;
class THnSparse;
class TRandom3;

class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliAODMCParticle; // sample
class AliEMCALTriggerPatchInfo;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliMultSelection;
class AliSelectNonHFE;


#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "TProfile.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliVEvent.h"
#include "AliSelectNonHFE.h"




class AliAnalysisTaskHFETPCTOFMultiplicity : public AliAnalysisTaskSE
{
public:
    
    enum HijingOrNot {kHijing,kElse};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    
    AliAnalysisTaskHFETPCTOFMultiplicity();
    AliAnalysisTaskHFETPCTOFMultiplicity(const char *name);
    virtual                 ~AliAnalysisTaskHFETPCTOFMultiplicity();
    virtual void         Init();
    virtual void           LocalInit() {Init();}
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    
    
    
    //-----------------selections cuts----------------------------------------------------------------
    Bool_t          PassEventSelect(AliAODEvent *fAOD); //event selection cuts
    Bool_t        Passtrackcuts(AliAODTrack *atrack); // track selection cuts
    //-----------------getters------------------------------------------------------------------------
    Bool_t          IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
    Bool_t          GetNMCPartProduced();
    Int_t           GetNcharged(); // n charge generated
    Int_t           GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray);
    Int_t           GetPi0EtaType(AliAODMCParticle *part);
    Double_t        GetTrackletsMeanCorrection(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult); /*const*/
    void            GetPi0EtaWeight(THnSparse *SparseWeight);
    Double_t        Beta(AliAODTrack *track);
    Bool_t          GetNonHFEEffiDenomGenPurMC(AliVTrack *track, Int_t correctednAcc1);
    Bool_t          GetNonHFEEffiRecoTagGenPurMC(AliVTrack *track, Int_t correctednAcc1);
    Bool_t          GetNonHFEEffiULSLSGenPurMC(AliAODTrack *track, AliAODTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass, Int_t correctednAcc1);
    //-----------------setters---------------------------------------------------------------- -------
    void            SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Int_t correctednAcc1);
    void             SetReferenceMultiplicity(Double_t multi){fRefMult=multi;};
    void             SetMultiProfileLHC16s(TProfile * hprof){
        if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
        fMultEstimatorAvg[0]=new TProfile(*hprof);
    }
    void             SetMultiProfileLHC16r(TProfile * hprof){
        if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
        fMultEstimatorAvg[1]=new TProfile(*hprof);
    }
    
    void             SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;};
    void            SwitchPi0EtaWeightCalc(Bool_t fSwitch) {fCalculateWeight = fSwitch;};
    
    void             SetNcontVCut(Int_t NcontV){fCutNcontV=NcontV;}
    void             SetEtaRange(Double_t Etarange){fCutTrackEta=Etarange;}
    void             SetMaxTPCCluster(Int_t MaxTPCclus){fCutTPCMaxCls=MaxTPCclus;}
    void             SetNTPCCluster(Int_t TPCNclus){fCutTPCNCls=TPCNclus;}
    void             SetNITSCluster(Int_t ITSNclus){fCutITSNCls=ITSNclus;}
    void             SetTOFNSigmaCut(Double_t CutNsigmaTOF){fCutNsigmaTOF = CutNsigmaTOF;}
    
    void             SetDCACut(Double_t DCAxyCut,Double_t DCAzCut){
        fCutDCAxy=DCAxyCut;
        fCutDCAz=DCAzCut;
    }
    
    void SetHitsOnSPDLayers(Bool_t SPDBoth,Bool_t SPDAny, Bool_t SPDFirst)    {
        fSPDBoth=SPDBoth;
        if(!SPDBoth)fSPDAny=SPDAny;
        if(!SPDBoth && !SPDAny)fSPDFirst=SPDFirst;
    }
    //------------Setters for PID electron------------------
    void             SetTPCnsigma(Double_t TPCNSigMin,Double_t TPCNSigMax){
        fCutNsigmaEMin=TPCNSigMin;
        fCutNsigmaEMax=TPCNSigMax;
    }
    
    
    
    //------------Setters for photonic electron pair------------------
    void SetInvMassCut(Double_t InvmassCut){fCutInvmass=InvmassCut;}
    void SetAssoTPCclus(Int_t AssoTPCCluster){fAssoTPCCluster=AssoTPCCluster;}
    void SetAssoITSclus(Int_t AssoITSCluster){fAssoITSCluster=AssoITSCluster;}
    void SetAssoITSrefit(Bool_t AssoITSRefit){fAssoITSRefit= AssoITSRefit;}
    void SetAssopTMin(Double_t AssoEPt){fCutAssoEPt = AssoEPt;}
    void SetAssoEtarange(Double_t AssoEEtarange){fCutAssoEEta=AssoEEtarange;}
    void SetAssoTPCnsig(Double_t AssoENsigma){fCutAssoENsigma=AssoENsigma;}
    
    //---------------------------------------------------------------------------------
private:
    
    //--------------------Event Cut------------------
    Int_t         fCutNcontV;
    //--------------------Track Cut------------------
    Int_t            fCutTPCMaxCls;
    Int_t            fCutTPCchi2perNDF;
    Int_t         fCutTPCNCls;
    Int_t         fCutITSNCls;
    Double_t         fCutDCAxy;
    Double_t         fCutDCAz;
    Double_t         fCutTrackEta;
    Double_t         fCutNsigmaTOF;//!
    
    Bool_t fSPDBoth;
    Bool_t fSPDAny;
    Bool_t fSPDFirst; 
    //--------------------PID Cut------------------
    
    Double_t         fCutNsigmaEMin;
    Double_t         fCutNsigmaEMax;
    
    //--------------------Loose cuts for photonic electron pair------------------
    Int_t         fAssoTPCCluster;//!
    Int_t         fAssoITSCluster;//!
    Double_t         fCutAssoEPt;
    Double_t         fCutAssoEEta;
    Double_t         fCutAssoENsigma;
    Bool_t         fAssoITSRefit;//!
    //--------------------Mass Cut for photonic electron pair------------------
    Double_t         fCutInvmass;
    
    //---------------------------------------------------------------------------------
    AliAODEvent*          fAOD;//! input event
    AliPIDResponse*       fpidResponse;//!pid response
    
    AliSelectNonHFE*      fNonHFE; //!
    AliAODMCParticle*        fMCparticle;
    AliAODMCHeader*        fMCHeader;//!
    TClonesArray*            fMCArray;//! MC array
    TProfile*             GetEstimatorHistogram(const AliAODEvent *fAOD);
    
    Double_t              fTPCnSigma;//!
    Bool_t                 fRejectPUFromSPD;
    Bool_t                fCalculateElectronEffi;//
    
    Bool_t                  fCalculateWeight;//
    Int_t                   fNTotMCpart; //! N of total MC particles produced by generator
    Int_t                  fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t                   fNembMCpi0; //! N > fNMCpi0 = particles from pi0 generator
    Int_t                   fNembMCeta; //! N > fNMCeta = particles from eta generator
    
    TList*                fOutputList;//! output list
    TList*        fListProfiles; // list of profile histos for z-vtx correction
    
    TH1F*            fNevents;//! no of events
    TH1F*            fVtxZ;//! Vertex z
    TH1F*            fVtxX;//! Vertex x
    TH1F*            fVtxY;//! Vertex y
    TH2F*            fTPCdEdx;//! dedx vs pt
    TH2F*            fTPCnsigma;//! TPC Nsigma
    TH2F*            fTOFnsigma;//! TPC Nsigma
    TH1F*            fTrkPt;//! track pt
    TH1F*            fTrketa;//! track eta
    TH1F*            fTrkphi;//! track phi
    TF1*          fLandau;//!
    TF1*          fErr;//!
    TF1*          fNewHadFunc;//!
    TH2F*            fHadCont_Landau; //!
    TH2F*            fHadCont_Err; //!
    TH2F*            fPt_incl_e_Landau; //!
    TH2F*            fPt_incl_e_Err; //!
    
    TH2F*            fTPCdEdxVsPTOFCut;//!
    TH2F*            fTPCnSigmaVsPTOFCut;//!
    TH2F*         fTPCnSigmaProtonsVsPTOFCut;//!
    TH2F*         fTPCnSigmaKaonsVsPTOFCut;//!
    
    
    TH1F*                 fInvmassLS;//!
    TH1F*                 fInvmassULS;//!
    TH2F*                 fInvmassLSPt;//!
    TH2F*                 fInvmassULSPt;//!
    TH1F*                 fULSElecPt;//!
    TH1F*                 fLSElecPt;//!
    //MC histograms-------------------------------------------
    
    
    Bool_t                fReadMC;    //flag for access to MC
    Bool_t                fIsFrmPi0;//!
    Bool_t                fIsFrmEta;//!
    Int_t                 ftype;//!
    Double_t              fWeight;//!
    Double_t            fWeightPi0;//!
    Double_t            fWeightEta;//!
    
    TH1F*                 fInclsElecPt;//!
    TH1F*                 fInclsElecPtAll;//!
    TH1F*                 fInclsElecPtReco;//!
    
    TH2F*                 fHFElecPtAll;//!
    TH2F*                 fHFElecPtReco_wtrkcuts;//!
    TH2F*                 fHFElecPtReco_wTPCPID;//!
    TH2F*                 fHFElecPtReco_wtrkCalocuts;//!
    TH2F*                 fHFElecPtReco_wTPCTOFPID;//!
    TF1*                  fPi0Weight;//!
    TF1*                  fEtaWeight;//!
    TF1*                  fPi0Weight2;//!
    TF1*                  fEtaWeight2;//!
    
    TH2F*                    fNonHFeTrkPt;//!
    TH2F*                    fNonHFeWeightTrkPt;//!
    TH2F*                    fPi0eWeightTrkPt;//!
    TH2F*                    fEtaeWeightTrkPt;//!
    TH2F*                    fPi0EtaEleTrkPt;//!
    TH2F*                    fRecoPi0EtaEleTrkPt;//!
    
    TH1F*                    fNonHFePairInvmassLS;//!
    TH1F*                    fNonHFePairInvmassULS;//!
    
    TH2F*                    fRecoLSeTrkPt;//!
    TH2F*                    fRecoLSeWeightTrkPt;//!
    TH2F*                    fRecoPi0LSeWeightTrkPt;//!
    TH2F*                    fRecoEtaLSeWeightTrkPt;//!
    
    TH2F*                    fRecoULSeTrkPt;//!
    TH2F*                    fRecoULSeWeightTrkPt;//!
    TH2F*                    fRecoPi0ULSeWeightTrkPt;//!
    TH2F*                    fRecoEtaULSeWeightTrkPt;//!
    
    TH2F*                    fRecoNonHFeTrkPt;//!
    TH2F*                    fRecoNonHFeWeightTrkPt;//!
    TH2F*                    fRecoPi0eWeightTrkPt;//!
    TH2F*                    fRecoEtaeWeightTrkPt;//!
    
    
    THnSparse*        fSprsPi0EtaWeightCal;//!
    
    //--------------------------------------------------------
    TProfile*        fMultEstimatorAvg[2];
    Double_t         fRefMult;
    TRandom3*        gRandom;//!< random number generator
    
    THnSparse*             fSparseElectron; //! Electron information
    THnSparse*             fSparseLSElectron; //! Electron information
    THnSparse*             fSparseULSElectron; //! Electron information
    Double_t*             fvaluePHElectron;//!
    Double_t*             fvalueElectron; //!
    
    
    THnSparse*             fSparseMulti; //! Multiplicity information
    Double_t*             fvalueMulti; //!
    
    
    
    
    
    
    AliAnalysisTaskHFETPCTOFMultiplicity(const AliAnalysisTaskHFETPCTOFMultiplicity& ); // not implemented
    AliAnalysisTaskHFETPCTOFMultiplicity& operator=(const AliAnalysisTaskHFETPCTOFMultiplicity& ); // not implemented
    
    ClassDef(AliAnalysisTaskHFETPCTOFMultiplicity,1);
};

#endif
