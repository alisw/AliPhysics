#ifndef ALICONVK0LAMBDACUTS_H
#define ALICONVK0LAMBDACUTS_H
#include <TObjString.h>
#include "AliAODpidUtil.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "AliAODv0.h"


class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliPIDResponse;
class AliGAKFVertex;
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
 * @class AliConvK0LambdaCuts (adopted from AliConversionPhotonCuts)
 * @brief Class handling all kinds of selection cuts for K0 and Lambda V0s
 * @author Joshua Koenig
 * @ingroup GammaConvBase
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


class AliConvK0LambdaCuts : public AliAnalysisCuts {

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

    AliConvK0LambdaCuts(const char *name="K0sLambdaCuts", const char * title="K0sLambdaCuts");
    AliConvK0LambdaCuts(const AliConvK0LambdaCuts&);
    AliConvK0LambdaCuts& operator=(const AliConvK0LambdaCuts&);
    virtual ~AliConvK0LambdaCuts();

    virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

    // bool SetCutIds(TString cutString);
    int fCuts[kNCuts];
    bool SetCut(cutIds cutID, int cut);
    bool UpdateCutString();

    static const char * fgkCutNames[kNCuts];

    bool InitializeCutsFromCutString(const TString analysisCutSelection);

    bool InitPIDResponse();
    void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
    AliPIDResponse * GetPIDResponse() { return fPIDResponse;}

    TString GetCutNumber();

    // Cut Selection
    void SetLightOutput( int flag ){fDoLightOutput = flag; return;}
    void SetFillCutHistograms(TString name="",bool preCut = kTRUE) { InitCutHistograms(name, preCut); }
    void InitCutHistograms(TString name="",bool preCut = kTRUE);
    TList *GetCutHistograms(){return fHistograms;}
    void PrintCuts();
    void CallSumw2ForLists(TList* l);

    void SetIsMC(int mc) {fIsMC = mc;}
    void SetDoQA(bool qa) {fDoQA = qa;}

    void InitArmPodRefHistos(const char * filename);

    bool IsK0sLambdaAccepted(AliAODv0 *v0, int mesonPDGCode, double weight);

    // Set Individual Cuts
    bool SetRCut(int RCut);
    bool SetV0Finder(int v0FinderType);
    bool SetChi2V0Cut(int chi2GammaCut);
    bool SetTPCdEdxCut(int ededxSigmaCut);
    bool SetSinglePtCut(int singlePtCut);
    bool SetTPCClusterCut(int clsTPCCut);
    bool SetEtaCut(int etaCut);
    bool SetEtaForPhiCut(int etaPhiCut);
    bool SetMinPhiSectorCut(int minPhiCut);
    bool SetMaxPhiSectorCut(int maxPhiCut);
    bool SetMinMomPiondEdxCut(int piMinMomdedxSigmaCut);
    bool SetArmenterosQTCut(int QtMaxCut);
    bool SetTOFElectronPIDCut(int TOFelectronPID);
    bool SetAsymmetryCut(int doPhotonAsymmetryCut);
    bool SetPhotonRDepPtCut(int doPhotonRDepPtCut);
    bool SetCosPAngleCut(int cosCut);
    bool SetPsiPairCut(int psiCut);
    bool SetToCloseV0sCut(int toClose);
    bool SetDCAPrimVtxCut(int DCARPrimVtx);
    bool SetDCAZPrimVtxCut(int DCAZPrimVtx);
    void SetIsHeavyIon(int isHeavyIon){fIsHeavyIon=isHeavyIon;}
    bool SetITSElectronPIDCut(int ITSelectronPID);
    bool SetTRDElectronPIDCut(int TRDelectronPID);


    // Cut functions
    bool DoV0ReaderTypeCut(bool status) const;
    bool DoEtaCut(double etaV0) const;
    bool DoRCut(AliAODv0 *v0) const;
    bool DoSinglePtCut(const AliAODTrack *trNeg, const AliAODTrack *trPos) const;
    bool DoTPCClusCut(const AliAODTrack *trNeg, const AliAODTrack *trPos) const;
    bool DodEdXCut(int mode, const AliAODTrack *trNeg, const AliAODTrack *trPos);
    bool DoChi2Cut(double chi2) const;
    bool DoCosPACut(double cosPA) const;
    bool DoDCAToPrimVtxCut(double dca) const;
    bool DoArmenterosQtCut(AliAODv0* v0, int mode) const;

    void FillTrueHistograms(AliAODv0 *v0, int mode, double weight);

  private:
    TList*            fHistograms;                          ///< List of QA histograms
    AliPIDResponse*   fPIDResponse;                         ///< PID response

    int             fIsMC;
    int             fDoLightOutput;                       ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    bool            fDoQA;                                ///< switch for running QA

    // //cuts
    double          fMaxR;                                ///< r cut
    double          fMinR;                                ///< r cut
    double          fEtaCut;                              ///< eta cut
    double          fEtaCutMin;                           ///< eta cut
    double          fSinglePtCutPos;                         ///< pt cut for positive track
    double          fSinglePtCutNeg;                         ///< pt cut for negative track
    double          fMinClsTPC;                           ///< minimum clusters in the TPC
    double          fLineCutZRSlope;                      ///< linecut
    double          fLineCutZValue;                       ///< linecut
    double          fLineCutZRSlopeMin;                   ///< linecut
    double          fLineCutZValueMin;                    ///< linecut
    double          fChi2CutV0;                   ///< chi2cut
    double          fPIDnSigmaAbove;          ///< sigma cut
    double          fPIDnSigmaBelow;          ///< sigma cut
    double          fTofPIDnSigmaAbove;       ///< sigma cut RRnewTOF
    double          fTofPIDnSigmaBelow;       ///< sigma cut RRnewTOF
    bool            fDoArmenteros1DCuts;
    bool            fDoArmenteros2DCuts;                              ///< Select gammas using ellipse cut
    double          maxDevNegArmPod2D;
    double          maxDevPosArmPod2D;
    double          maxRangeHist2DArmPod;
    double          fQtMax;                               ///< Maximum Qt from Armenteros to select Gammas
    double          fQtMin;                             ///< Maximum Qt vs Pt param from Armenteros to select Gammas
    double          fAlphaMin;
    double          fAlphaMax;
    bool            fUseOnFlyV0Finder;                    ///< flag
    float           fCosPAngleCut;                        ///<
    double          fDCAZPrimVtxCut;                      ///< cut value for the maximum distance in Z between the photon & the primary vertex [cm]
    double          fDCAPrimVtxCut;                      ///< cut value for the maximum distance in R between the photon & the primary vertex [cm]
    TObjString*       fCutString;                           ///< cut number used for analysis
    TString           fCutStringRead;                       ///<
    int             fIsHeavyIon;                          ///< flag for pp (0), PbPb (1), pPb (2)
    double          fExcludeMinR;                         ///< r cut exclude region
    double          fExcludeMaxR;                         ///< r cut exclude region

    // cut histograms
    TH2F*           fHistArmPodRefK0s;                      //
    TH2F*           fHistArmPodRefLambda;                   //
    TH2F*           fHistArmPodRefAntiLambda;               //


    // Histograms
    TH2F*             fHistoCutIndex;                       ///< histogram with individual cuts vs.v0 pt   
    TH2F*             fHistoArmenterosbefore;               ///< armenteros podolanski plot before cuts
    TH2F*             fHistoArmenterosafter;                ///< armenteros podolanski plot after cuts
    TH2F*             fHistoChi2before;                     ///< Chi2 distribution before cuts
    TH2F*             fHistoChi2after;                      ///< Chi2 distribution after cuts
    TH2F*             fHistoArmenterosTrue;                 ///< armenteros podolanski plot for true particles
    TH2F*             fHistoNSigmaPosTrackTrue;             ///< NSigma TPC of positive tracks for true particles
    TH2F*             fHistoNSigmaNegTrackTrue;             ///< NSigma TPC of negative tracks for true particles

    /// \cond CLASSIMP
    ClassDef(AliConvK0LambdaCuts,2)
    /// \endcond
};

#endif
