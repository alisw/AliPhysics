#ifndef ALICONVERSIONPHOTONCUTS_H
#define ALICONVERSIONPHOTONCUTS_H

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"


class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliPIDResponse;
class AliKFVertex;
class TH1F;
class TH2F;
class TF1;
class TProfile;
class TProfile2D;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;

/**
 * @class AliConversionPhotonCuts
 * @brief Class handling all kinds of selection cuts for Gamma Conversion analysis
 * @author Friederike Bock
 * @ingroup GammaConv
 *
 * The cut configuration is set as a string with an 19 digit number.
 * Each digit in the string corresponds to a certain cut type, while
 * its values represent the cut values. The cut configuration is listed here:
 *
 * | Position in the cut string (from the end) | Cut type               |
 * |-------------------------------------------|------------------------|
 * |                  0                        | V0FinderType           |
 * |                  1                        | EtaCut                 |
 * |                  2                        | MinRCut                |
 * |                  3                        | EtaForPhiCut           |
 * |                  4                        | MinPhiCut              |
 * |                  5                        | MaxPhiCut              |
 * |                  6                        | SinglePtCut            |
 * |                  7                        | ClsTPCCut              |
 * |                  8                        | ededxSigmaCut          |
 * |                  9                        | pidedxSigmaCut         |
 * |                  10                       | piMomdedxSigmaCut      |
 * |                  11                       | piMaxMomdedxSigmaCut   |
 * |                  12                       | LowPRejectionSigmaCut  |
 * |                  13                       | TOFelectronPID         |
 * |                  14                       | ITSelectronPID         |
 * |                  15                       | TRDelectronPID         |
 * |                  16                       | QtMaxCut               |
 * |                  17                       | Chi2GammaCut           |
 * |                  18                       | PsiPair                |
 * |                  19                       | DoPhotonAsymmetryCut   |
 * |                  20                       | CosinePointingAngle    |
 * |                  21                       | SharedElectronCuts     |
 * |                  22                       | RejectToCloseV0s       |
 * |                  23                       | DcaRPrimVtx            |
 * |                  24                       | DcaZPrimVtx            |
 * |                  25                       | EvetPlane              |
*/


class AliConversionPhotonCuts : public AliAnalysisCuts {

  public:
    enum cutIds {
      kv0FinderType,
      ketaCut,
      kRCut,
      kEtaForPhiSector,
      kMinPhiSector,
      kMaxPhiSector,
      ksinglePtCut,
      kclsTPCCut,
      kededxSigmaCut,
      kpidedxSigmaCut,
      kpiMomdedxSigmaCut,
      kpiMaxMomdedxSigmaCut,
      kLowPRejectionSigmaCut,
      kTOFelectronPID,
      kITSelectronPID,
      kTRDelectronPID,
      kQtMaxCut,
      kchi2GammaCut,
      kPsiPair,
      kdoPhotonAsymmetryCut,
      kCosPAngle,
      kElecShare,
      kToCloseV0s,
      kDcaRPrimVtx,
      kDcaZPrimVtx,
      kInPlaneOutOfPlane,
      kNCuts
    };

    enum photonCuts {
        kPhotonIn=0,
        kOnFly,
        kNoTracks,
        kNoV0,
        kTrackCuts,
        kdEdxCuts,
        kConvPointFail,
        kPhotonCuts,
        kEventPlane,
        kPhotonOut
    };


    Bool_t SetCutIds(TString cutString);
    Int_t fCuts[kNCuts];
    Bool_t SetCut(cutIds cutID, Int_t cut);
    Bool_t UpdateCutString();

    static const char * fgkCutNames[kNCuts];

    Double_t GetCosineOfPointingAngle(const AliConversionPhotonBase * photon, AliVEvent * event) const;
    Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
    void FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);
    void SetPreSelectionCutFlag(Bool_t preSelFlag){fPreSelCut = preSelFlag;}

    AliConversionPhotonCuts(const char *name="V0Cuts", const char * title="V0 Cuts");
    AliConversionPhotonCuts(const AliConversionPhotonCuts&);
    AliConversionPhotonCuts& operator=(const AliConversionPhotonCuts&);

    virtual ~AliConversionPhotonCuts();                            //virtual destructor

    static AliConversionPhotonCuts * GetStandardCuts2010PbPb();
    static AliConversionPhotonCuts * GetStandardCuts2010pp();

    Bool_t InitPIDResponse();
    void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
    AliPIDResponse * GetPIDResponse() { return fPIDResponse;}


    virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

    TString GetCutNumber();

    Float_t GetKappaTPC(AliConversionPhotonBase *gamma, AliVEvent *event);
    Bool_t GetBDTVariableValues(AliConversionPhotonBase *gamma, AliVEvent *event, Float_t* values);

    // Cut Selection
    Bool_t TrackIsSelected(AliConversionPhotonBase * photon, AliVEvent  * event);
    Bool_t PhotonIsSelected(AliConversionPhotonBase * photon, AliVEvent  * event);
    Bool_t PhotonIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent,Bool_t checkForConvertedGamma=kTRUE);
    Bool_t PhotonIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray,Bool_t checkForConvertedGamma=kTRUE);
    Bool_t PhotonIsSelectedMCAODESD(AliDalitzAODESDMC *particle,AliDalitzEventMC *mcEvent,Bool_t checkForConvertedGamma) const;
    //Bool_t ElectronIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent);
    Bool_t TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack);
    //Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE);
    //Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliMCEvent *mcEvent, Bool_t bMCDaughtersInAcceptance=kFALSE);

    void PrintCuts();
    void PrintCutsWithValues();

    void SetLightOutput( Int_t flag ){fDoLightOutput = flag; return;}
    void SetPlotTrackPID( Bool_t flag ){fDoPlotTrackPID = flag; return;}
    void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
    void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
    TList *GetCutHistograms(){return fHistograms;}
    void FillPhotonCutIndex(Int_t photoncut){if(fHistoCutIndex)fHistoCutIndex->Fill(photoncut);}
    void FillV0EtaBeforedEdxCuts(Float_t v0Eta){if(fHistoEtaDistV0s)fHistoEtaDistV0s->Fill(v0Eta);}
    void FillV0EtaAfterdEdxCuts(Float_t v0Eta){if(fHistoEtaDistV0sAfterdEdxCuts)fHistoEtaDistV0sAfterdEdxCuts->Fill(v0Eta);}

    void SetV0ReaderName(TString name){fV0ReaderName = name; return;}
    void SetProcessAODCheck(Bool_t flag){fProcessAODCheck = flag; return;}

    AliVTrack * GetTrack(AliVEvent * event, Int_t label);
    AliESDtrack *GetESDTrack(AliESDEvent * event, Int_t label);

    ///Cut functions
    Bool_t SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex);
    Bool_t SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex);
    Bool_t AcceptanceCuts(AliConversionPhotonBase *photon);
    Bool_t AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg);
    Bool_t PhiSectorCut(AliConversionPhotonBase * photon);
    //   Bool_t dEdxCuts(AliVTrack * track);
    Bool_t dEdxCuts(AliVTrack * track, AliConversionPhotonBase * photon);
    Bool_t KappaCuts(AliConversionPhotonBase * photon,AliVEvent *event);
    Bool_t ArmenterosQtCut(AliConversionPhotonBase *photon);
    Bool_t AsymmetryCut(AliConversionPhotonBase *photon,AliVEvent *event);
    Bool_t PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event);
    Bool_t SelectV0Finder(Bool_t onfly){
      if(onfly == fUseOnFlyV0Finder) return kTRUE;
      else return kFALSE;
    }
    Bool_t PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event);
    Bool_t CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event);
    Bool_t PsiPairCut(const AliConversionPhotonBase * photon) const;
    Bool_t CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const;
    Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
    Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);

    UChar_t DeterminePhotonQualityAOD(AliAODConversionPhoton*, AliVEvent*);
    UChar_t DeterminePhotonQualityTRD(AliAODConversionPhoton*, AliVEvent*);
    UChar_t DeterminePhotonQualityTOF(AliAODConversionPhoton*, AliVEvent*);
    Bool_t InPlaneOutOfPlaneCut(Double_t photonPhi, Double_t eventPlaneAngle = -100, Bool_t fill = kTRUE);
    Int_t GetInPlaneOutOfPlaneCut(){return fInPlaneOutOfPlane;}

    // Set Individual Cuts
    Bool_t SetRCut(Int_t RCut);
    Bool_t SetV0Finder(Int_t v0FinderType);
    Bool_t SetChi2GammaCut(Int_t chi2GammaCut);
    Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
    Bool_t SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut);
    Bool_t SetSinglePtCut(Int_t singlePtCut);
    Bool_t SetTPCClusterCut(Int_t clsTPCCut);
    Bool_t SetEtaCut(Int_t etaCut);
    Bool_t SetEtaForPhiCut(Int_t etaPhiCut);
    Bool_t SetMinPhiSectorCut(Int_t minPhiCut);
    Bool_t SetMaxPhiSectorCut(Int_t maxPhiCut);
    Bool_t SetMinMomPiondEdxCut(Int_t piMinMomdedxSigmaCut);
    Bool_t SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut);
    Bool_t SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut);
    Bool_t SetQtMaxCut(Int_t QtMaxCut);
    Bool_t SetTOFElectronPIDCut(Int_t TOFelectronPID);
    Bool_t SetTRDElectronCut(Int_t TRDElectronCut);
    Bool_t SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut);
    Bool_t SetCosPAngleCut(Int_t cosCut);
    Bool_t SetPsiPairCut(Int_t psiCut);
    Bool_t SetSharedElectronCut(Int_t sharedElec);
    Bool_t SetToCloseV0sCut(Int_t toClose);
    Bool_t SetDCARPhotonPrimVtxCut(Int_t DCARPhotonPrimVtx);
    Bool_t SetDCAZPhotonPrimVtxCut(Int_t DCAZPhotonPrimVtx);
    Bool_t SetInPlaneOutOfPlane(Int_t inOutPlane);
    Bool_t SetKappaTPCCut(Int_t kappaCut);
    void SetIsHeavyIon(Int_t isHeavyIon){fIsHeavyIon=isHeavyIon;}
    Int_t GetFirstTPCRow(Double_t radius);
    void GetPhiRegions();
    Bool_t SetITSElectronPIDCut(Int_t ITSelectronPID);
    Bool_t SetTRDElectronPIDCut(Int_t TRDelectronPID);


    // Request Flags
    Bool_t UseElecSharingCut(){return fDoSharedElecCut;}
    Bool_t UseToCloseV0sCut(){return fDoToCloseV0sCut;}
    Double_t GetEtaCut(){return fEtaCut;}
    Double_t GetSingleElectronPtCut(){return fSinglePtCut;}
    void SetDodEdxSigmaCut(Bool_t k=kTRUE)  {fDodEdxSigmaCut=k;}
    void SetSwitchToKappaInsteadOfNSigdEdxTPC(Bool_t k=kTRUE) {fSwitchToKappa=k;}
    void SetDoElecDeDxPostCalibration(Bool_t k=kTRUE)  {fDoElecDeDxPostCalibration=k;}
    Bool_t GetMaterialBudgetWeightsInitialized() {return fMaterialBudgetWeightsInitialized;}
    Bool_t InitializeMaterialBudgetWeights(Int_t flag, TString filename);
    Float_t GetMaterialBudgetCorrectingWeightForTrueGamma(AliAODConversionPhoton* gamma, Double_t magField);

    Int_t GetV0FinderSameSign(){return fUseOnFlyV0FinderSameSign;}
    Bool_t GetUseBDTPhotonCuts(){return fUseBDTPhotonCuts;}
    Int_t GetDoElecDeDxPostCalibration() {return fDoElecDeDxPostCalibration;}
    void  SetElecDeDxPostCalibrationCustomFile(TString filename){fFileNameElecDeDxPostCalibration = filename; return;};
    Bool_t InitializeElecDeDxPostCalibration(TString filename);
    Bool_t LoadElecDeDxPostCalibration(Int_t runNumber);
    Double_t GetCorrectedElectronTPCResponse(Short_t charge,Double_t nsig,Double_t P,Double_t Eta,Double_t TPCCl, Double_t R);
    void ForceTPCRecalibrationAsFunctionOfConvR(){fIsRecalibDepTPCCl = kFALSE;}

  protected:
    TList*            fHistograms;                          ///< List of QA histograms
    AliPIDResponse*   fPIDResponse;                         ///< PID response

    Int_t            fDoLightOutput;                       ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t            fDoPlotTrackPID;                       ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    TString           fV0ReaderName;                        ///< Name of the V0 reader

    //cuts
    Double_t          fMaxR;                                ///< r cut
    Double_t          fMinR;                                ///< r cut
    Double_t          fEtaCut;                              ///< eta cut
    Double_t          fEtaCutMin;                           ///< eta cut
    Float_t           fEtaForPhiCutMin;                     ///< eta cut for phi sector selection
    Float_t           fEtaForPhiCutMax;                     ///< eta cut for phi sector selection
    Float_t           fMinPhiCut;                           ///< phi sector cut
    Float_t           fMaxPhiCut;                           ///< phi sector cut
    Int_t             fDoShrinkTPCAcceptance;               ///< Flag for shrinking the TPC acceptance due to different reasons
    Double_t          fPtCut;                               ///< pt cut
    Double_t          fSinglePtCut;                         ///< pt cut for electron/positron
    Double_t          fSinglePtCut2;                        ///< second pt cut for electron/positron if asymmetric cut is chosen
    Bool_t            fDoAsymPtCut;                         ///< Flag for setting asymmetric pT cut on electron/positron
    Double_t          fMaxZ;                                ///< z cut
    Double_t          fMinClsTPC;                           ///< minimum clusters in the TPC
    Double_t          fMinClsTPCToF;                        ///< minimum clusters to findable clusters
    Double_t          fMaxTPCChi2NDF;                       ///< maximum TPC track Chi2 per NDF
    Double_t          fLineCutZRSlope;                      ///< linecut
    Double_t          fLineCutZValue;                       ///< linecut
    Double_t          fLineCutZRSlopeMin;                   ///< linecut
    Double_t          fLineCutZValueMin;                    ///< linecut
    Double_t          fChi2CutConversion;                   ///< chi2cut
    Double_t          fChi2CutConversionExpFunc;            ///< chi2cutexpfunction
    Double_t          fPIDProbabilityCutNegativeParticle;   ///<
    Double_t          fPIDProbabilityCutPositiveParticle;   ///<
    Bool_t            fDodEdxSigmaCut;                      ///< flag to use the dEdxCut based on sigmas
    Bool_t            fDoTOFsigmaCut;                       ///< flag to use TOF pid cut RRnewTOF
    Double_t          fPIDTRDEfficiency;                    ///< required electron efficiency for TRD PID
    Bool_t            fDoTRDPID;                            ///< flag to use TRD pid
    Double_t          fPIDnSigmaAboveElectronLine;          ///< sigma cut
    Double_t          fPIDnSigmaBelowElectronLine;          ///< sigma cut
    Double_t          fTofPIDnSigmaAboveElectronLine;       ///< sigma cut RRnewTOF
    Double_t          fTofPIDnSigmaBelowElectronLine;       ///< sigma cut RRnewTOF
    Double_t          fPIDnSigmaAbovePionLine;              ///< sigma cut
    Double_t          fPIDnSigmaAbovePionLineHighPt;        ///< sigma cut
    Double_t          fPIDMinPnSigmaAbovePionLine;          ///< sigma cut
    Double_t          fPIDMaxPnSigmaAbovePionLine;          ///< sigma cut
    Double_t          fDoKaonRejectionLowP;                 ///< Kaon rejection at low p
    Double_t          fDoProtonRejectionLowP;               ///< Proton rejection at low p
    Double_t          fDoPionRejectionLowP;                 ///< Pion rejection at low p
    Double_t          fPIDnSigmaAtLowPAroundKaonLine;       ///< sigma cut
    Double_t          fPIDnSigmaAtLowPAroundProtonLine;     ///< sigma cut
    Double_t          fPIDnSigmaAtLowPAroundPionLine;       ///< sigma cut
    Double_t          fPIDMinPKaonRejectionLowP;            ///< Momentum limit to apply kaon rejection
    Double_t          fPIDMinPProtonRejectionLowP;          ///< Momentum limit to apply proton rejection
    Double_t          fPIDMinPPionRejectionLowP;            ///< Momentum limit to apply proton rejection
    Int_t             fDoQtGammaSelection;                  ///< Select gammas using qtMax
    Bool_t            fDo2DQt;                              ///< Select gammas using ellipse cut
    Double_t          fQtMax;                               ///< Maximum Qt from Armenteros to select Gammas
    Double_t          fQtPtMax;                             ///< Maximum Qt vs Pt param from Armenteros to select Gammas
    Double_t          fNSigmaMass;                          ///< nsigma cut
    Bool_t            fUseEtaMinCut;                        ///< flag
    Bool_t            fUseOnFlyV0Finder;                    ///< flag
    Int_t             fUseOnFlyV0FinderSameSign;            ///< int to set same sign pairing
    Bool_t            fUseBDTPhotonCuts;                    ///< flag
    Bool_t            fDoPhotonAsymmetryCut;                ///< flag to use the PhotonAsymetryCut
    Bool_t            fDoPhotonPDependentAsymCut;           ///< flag to use the PhotonAsymetryCut with P dependent cut
    TF1 *             fFAsymmetryCut;                       ///<
    Double_t          fMinPPhotonAsymmetryCut;              ///< Min Momentum for Asymmetry Cut
    Double_t          fMinPhotonAsymmetry;                  ///< Asymmetry Cut
    Double_t          fMaxPhotonAsymmetry;                  ///< Asymmetry Cut
    Bool_t            fUseCorrectedTPCClsInfo;              ///< flag to use corrected tpc cl info
    Bool_t            fUseTOFpid;                           ///< flag to use tof pid
    Bool_t            fUseTOFtiming;                        ///< flag to use tof timing information
    Double_t          fTOFtimeMin;                          ///< minimum TOF time cut on conversion leg
    Double_t          fTOFtimeMax;                          ///< maximum TOF time cut on conversion leg
    Bool_t            fTOFtimingBothLegs;                   ///< flag to use tof timing on both or either one photon leg
    Bool_t            fUseTOFpidMomRange;                   ///< flag to use TOF nSigma cut only above a certain track momentum
    Double_t          fTofPIDMinMom;                        ///< track momentum threshold for TOF nSigma cut
    Double_t          fTofPIDMaxMom;                        ///< track momentum threshold for TOF nSigma cut
    Float_t           fOpeningAngle;                        ///< min opening angle for meson
    Float_t           fPsiPairCut;                          ///<
    Int_t             fDo2DPsiPairChi2;                     ///<
    Bool_t            fIncludeRejectedPsiPair;              ///<
    Float_t           fCosPAngleCut;                        ///<
    Bool_t            fDoToCloseV0sCut;                     ///<
    Double_t          fminV0Dist;                           ///<
    Bool_t            fDoSharedElecCut;                     ///<
    Bool_t            fDoPhotonQualitySelectionCut;         ///<
    Bool_t            fDoPhotonQualityRejectionCut;         ///<
    Int_t             fPhotonQualityCut;                    ///<
    Int_t             fPhotonQualityCutTRD;                 ///<
    Int_t             fPhotonQualityCutTOF;                 ///<
    TRandom3          fRandom;                              ///<
    Int_t             fElectronArraySize;                   ///< Size of electron array
    Int_t*            fElectronLabelArray;                  //[fElectronArraySize]
    Double_t          fDCAZPrimVtxCut;                      ///< cut value for the maximum distance in Z between the photon & the primary vertex [cm]
    Double_t          fDCARPrimVtxCut;                      ///< cut value for the maximum distance in R between the photon & the primary vertex [cm]
    Int_t             fInPlaneOutOfPlane;                   ///< In-Plane Out-Of Plane Analysis
    Float_t           fConversionPointXArray;               ///< Array with conversion Point x
    Float_t           fConversionPointYArray;               ///< Array with conversion Point y
    Float_t           fConversionPointZArray;               ///< Array with conversion Point z
    TObjString*       fCutString;                           ///< cut number used for analysis
    TString           fCutStringRead;                       ///<
    Int_t             fIsHeavyIon;                          ///< flag for pp (0), PbPb (1), pPb (2)
    Bool_t            fUseITSpid;                           ///< flag to use tof pid
    Double_t          fITSPIDnSigmaAboveElectronLine;       ///< sigma cut RRnewTOF
    Double_t          fITSPIDnSigmaBelowElectronLine;       ///< sigma cut RRnewTOF
    Double_t          fMaxPtPIDITS;                         ///< max pt for ITS PID
    Double_t          fTRDPIDAboveCut;                      ///< TRD cut range
    Double_t          fTRDPIDBelowCut;                      ///< TRD cut range
    Bool_t            fDoDoubleCountingCut;                 ///< Flag to reject double counting
    Double_t          fMinRDC;                              ///< Min R for Double Counting Cut
    Double_t          fDeltaR;                              ///< Delta R for Double Counting Cut
    Double_t          fOpenAngle;                           ///< Opening Angle for Double Counting Cut
    Bool_t            fSwitchToKappa;                       ///< switches from standard dEdx nSigma TPC cuts to Kappa TPC
    Float_t           fKappaMinCut;                         ///< maximum Kappa cut
    Float_t           fKappaMaxCut;                         ///< maximum Kappa cut
    Bool_t            fDoElecDeDxPostCalibration;           ///<
    Bool_t            fIsRecalibDepTPCCl;                   ///<
    // Histograms
    TH1F*             fHistoEtaDistV0s;                     ///< eta-distribution of all V0s after Finder selection
    TH1F*             fHistoEtaDistV0sAfterdEdxCuts;        ///< eta-distribution of all V0s after Finder selection after dEdx cuts
    TH2F*             fHistodEdxCuts;                       ///< bookkeeping for dEdx cuts
    TH2F*             fHistoTPCdEdxbefore;                  ///< TPC dEdx before cuts
    TH2F*             fHistoTPCdEdxafter;                   ///< TPC dEdx after cuts
    TH2F*             fHistoTPCdEdxSigbefore;               ///< TPC Sigma dEdx before cuts
    TH2F*             fHistoTPCdEdxSigafter;                ///< TPC Sigm dEdx after cuts
    TH1F*             fHistoTPCChi2NDFBefore;               ///< TPC track Chi2 before cuts
    TH1F*             fHistoTPCChi2NDFAfter;                ///< TPC track Chi2 after cuts
    TH2F*             fHistoTPCChi2NDF2D;                   ///< TPC neg vs pos track Chi2 before cuts
    TH2F*             fHistoKappaafter;                     ///< Kappa vs photon pt after cuts
    TH2F*             fHistoTOFbefore;                      ///< TOF before cuts
    TH2F*             fHistoTOFSigbefore;                   ///< TOF Sigma before cuts
    TH2F*             fHistoTOFSigafter;                    ///< TOF Sigma after cuts
    TH2F*             fHistoITSSigbefore;                   ///< ITS Sigma before cuts
    TH2F*             fHistoITSSigafter;                    ///< ITS Sigma after cuts
    TH2F*             fHistoPsiPairDeltaPhiafter;           ///< TOF Sigma after cuts
    TH1F*             fHistoTrackCuts;                      ///< bookkeeping for track cuts
    TH2F*             fHistoTrackPID;                      ///< bookkeeping for track pid
    TH2F*             fHistoPhotonCuts;                     ///< bookkeeping for photon specific cuts
    TH1F*             fHistoInvMassbefore;                  ///< e+e- inv mass distribution before cuts
    TH2F*             fHistoArmenterosbefore;               ///< armenteros podolanski plot before cuts
    TH1F*             fHistoInvMassafter;                   ///< e+e- inv mass distribution after cuts
    TH2F*             fHistoArmenterosafter;                ///< armenteros podolanski plot after cuts
    TH2F*             fHistoAsymmetrybefore;                 ///< asymmetry plot after cuts
    TH2F*             fHistoAsymmetryafter;                 ///< asymmetry plot after cuts
    TH2F*             fHistoAcceptanceCuts;                 ///< bookkeeping for acceptance cuts
    TH1F*             fHistoCutIndex;                       ///< bookkeeping for cuts
    TH2F*             fHistoTOFtimeVSMomentum;              ///< TOF timing (ns) versus e+- momentum
    TH1F*             fHistoEventPlanePhi;                  ///< EventPlaneAngle Minus Photon Angle
    Bool_t            fPreSelCut;                           ///< Flag for preselection cut used in V0Reader
    Bool_t            fProcessAODCheck;                     ///< Flag for processing check for AOD to be contained in AliAODs.root and AliAODGammaConversion.root
    Bool_t            fMaterialBudgetWeightsInitialized;    ///< weights for conversions photons due due deviating material budget in MC compared to data
    TProfile2D*       fProfileContainingMaterialBudgetWeights;
    TString           fFileNameElecDeDxPostCalibration;     ///< name of recalibration file (if no special non-OADB is required)
    Bool_t            fElecDeDxPostCalibrationInitialized;  ///< flag to check that initialization worked
    Int_t             fRecalibCurrentRun;                   ///< runnumber for correct loading of recalib from OADB
    Int_t             fnRBins;                              //
    TH2S**            fHistoEleMapRecalib;  //[fnRBins]
    TH2S**            fHistoPosMapRecalib;  //[fnRBins]
    Double_t          fGoodRegionCMin;                      ///< regions WITHOUT strong space charge distortions on C side
    Double_t          fGoodRegionAMin;                      ///<
    Double_t          fBadRegionCMin;                       ///< regions WITH strong space charge distortions on C side
    Double_t          fBadRegionAMin;                       ///<
    Double_t          fGoodRegionCMax;                      ///<
    Double_t          fGoodRegionAMax;                      ///<
    Double_t          fBadRegionCMax;                       ///<
    Double_t          fBadRegionAMax;                       ///<
    Double_t          fExcludeMinR;                         ///< r cut exclude region
    Double_t          fExcludeMaxR;                         ///< r cut exclude region

  private:
    /// \cond CLASSIMP
    ClassDef(AliConversionPhotonCuts,36)
    /// \endcond
};

#endif
