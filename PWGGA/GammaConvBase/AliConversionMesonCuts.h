#ifndef ALICONVERSIONMESONCUTS_H
#define ALICONVERSIONMESONCUTS_H

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
#include "AliAODMCParticle.h"
#include "AliCaloPhotonCuts.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


/**
 * @class AliConversionMesonCuts
 * @brief Class handling all kinds of selection cuts for Gamma Conversion analysis
 * @author Svein Lindal
 * @author Daniel Lohner
 * @ingroup GammaConv
 *
 * The cut configuration is set as a string with an 19 digit number.
 * Each digit in the string corresponds to a certain cut type, while
 * its values represent the cut values. The cut configuration is listed here:
 *
 * | Position in the cut string (from the end) | Cut type                 |
 * |-------------------------------------------|--------------------------|
 * |                  0                        | MesonKind                |
 * |                  1                        | BackgroundScheme         |
 * |                  2                        | NumberOfBGEvents         |
 * |                  3                        | DegreesForRotationMethod |
 * |                  4                        | RapidityMesonCut         |
 * |                  5                        | PtCut                    |
 * |                  6                        | AlphaMesonCut            |
 * |                  7                        | SelectionWindow          |
 * |                  8                        | SharedElectronCuts       |
 * |                  9                        | RejectToCloseV0s         |
 * |                  10                       | UseMCPSmearing           |
 * |                  11                       | DcaGammaGamma            |
 * |                  12                       | DcaRPrimVtx              |
 * |                  13                       | DcaZPrimVtx              |
 * |                  14                       | MinOpanMesonCut          |
 * |                  15                       | MaxOpanMesonCut          |
*/


class AliConversionMesonCuts : public AliAnalysisCuts {

  public:


    enum cutIds {
      kMesonKind,
      kBackgroundScheme,
      kNumberOfBGEvents,
      kDegreesForRotationMethod,
      kRapidityMesonCut,
      kPtCut,
      kalphaMesonCut,
      kSelectionCut,
      kElecShare,
      kToCloseV0s,
      kuseMCPSmearing,
      kDcaGammaGamma,
      kDcaRPrimVtx,
      kDcaZPrimVtx,
      kMinOpanMesonCut,
      kMaxOpanMesonCut,
      kNCuts
    };

    Bool_t  SetCutIds(TString cutString);
    Int_t   fCuts[kNCuts];
    Bool_t  SetCut(cutIds cutID, Int_t cut);
    Bool_t  UpdateCutString();

    static const char * fgkCutNames[kNCuts];

    Bool_t  InitializeCutsFromCutString(const TString analysisCutSelection);
    void    FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);

    AliConversionMesonCuts(const char *name="MesonCuts", const char * title="Meson Cuts");
    AliConversionMesonCuts(const AliConversionMesonCuts&);
    AliConversionMesonCuts& operator=(const AliConversionMesonCuts&);

    virtual ~AliConversionMesonCuts();                            //virtual destructor

    virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
    Bool_t MesonIsSelectedByMassCut (AliAODConversionMother *meson, Int_t nominalRange);
    Bool_t ArmenterosLikeQtCut(Double_t alpha, Double_t qT);

    TString GetCutNumber();

    // Cut Selection
    Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE, Double_t fRapidityShift=0., Int_t leadingCellID1 = 0, Int_t leadingCellID2 = 0, Char_t recoMeth1 = 0, Char_t  recoMeth2 = 0);
    Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliMCEvent *mcEvent, Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedAODMC(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedMCAODESD(AliDalitzAODESDMC *fMCMother,AliDalitzEventMC *mcEvent, Double_t fRapidityShift=0.) const;
    Bool_t MesonIsSelectedMCDalitz(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma,Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedAODMCDalitz(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma,Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedMCDalitzAODESD(AliDalitzAODESDMC* fMCMother,AliDalitzEventMC *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma,Double_t fRapidityShift=0.) const;
    Bool_t MesonIsSelectedMCEtaPiPlPiMiGamma(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedAODMCEtaPiPlPiMiGamma(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCPiPlPiMiEta(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedAODMCPiPlPiMiEta(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCPiPlPiMiPiZero(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedAODMCPiPlPiMiPiZero(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCPiZeroGamma(TParticle *fMCMother, AliMCEvent *mcEvent, Int_t &labelNeutPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedAODMCPiZeroGamma(AliAODMCParticle *fMCMother, TClonesArray *AODMCArray, Int_t &labelNeutPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCChiC(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &, Int_t &, Int_t &, Double_t fRapidityShift=0. );
    Bool_t MesonIsSelectedAODMCChiC(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &, Int_t &, Int_t &, Double_t fRapidityShift=0. );
    Bool_t MesonIsSelectedMCChiCAODESD(AliDalitzAODESDMC* fMCMother,AliDalitzEventMC *mcEvent, Int_t &, Int_t &, Int_t &, Double_t fRapidityShift=0. ) const;
    Bool_t MesonIsSelectedPiZeroGammaAngle(AliAODConversionMother *omega, AliAODConversionMother *pi0, AliAODConversionPhoton *gamma,
                                           Bool_t DoPiZeroAngleCut, TF1 *maxfit, Double_t lowerFactor, Double_t upperFactor);
    Bool_t MesonIsSelectedPiZeroGammaOAC(AliAODConversionMother *omega, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1, AliAODConversionPhoton *gamma2);
    void   PrintCuts();
    void   PrintCutsWithValues();

    void    SetLightOutput( Int_t flag ){fDoLightOutput = flag; return;}
    void    SetRunningMode(Int_t mode){fMode = mode; return;}
    void    InitCutHistograms(TString name="",Bool_t additionalHists=kFALSE);
    void    SetFillCutHistograms(TString name=""){if(!fHistograms){InitCutHistograms(name);};}
    TList   *GetCutHistograms(){return fHistograms;}
    void    SmearParticle(AliAODConversionPhoton * photon);
    void    SmearVirtualPhoton(AliAODConversionPhoton* photon);
    TLorentzVector SmearElectron(TLorentzVector particle);

    void    SetDefaultSmearing(Double_t p0, Double_t p1, Double_t p2){fUseMCPSmearing=1;fPBremSmearing=p0;fPSigSmearing=p1;fPSigSmearingCte=p2;return;}

    //Cut functions
    Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
    Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);

    void SetCaloMesonCutsObject(AliCaloPhotonCuts* cuts){fCaloPhotonCuts = cuts;}

    // Set Individual Cuts
    Bool_t SetMinPtCut(Int_t PtCut);
    Bool_t SetMesonKind(Int_t mesonKind);
    Bool_t SetSelectionWindowCut(Int_t selectionCut);
    Bool_t SetSelectionWindowMergedCut(Int_t selectionCut);
    Bool_t SetSelectionWindowCutPtDep(Int_t selectionCut);
    Bool_t SetAlphaMesonCut(Int_t alphaMesonCut);
    Bool_t SetAlphaMesonMergedCut(Int_t alphaMesonCut);
    Bool_t SetRapidityMesonCut(Int_t RapidityMesonCut);
    Bool_t SetBackgroundScheme(Int_t BackgroundScheme);
    Bool_t SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod);
    Bool_t SetNumberOfBGEvents(Int_t NumberOfBGEvents);
    Bool_t SetMCPSmearing(Int_t useMCPSmearing);
    Bool_t SetSharedElectronCut(Int_t sharedElec);
    Bool_t SetToCloseV0sCut(Int_t toClose);
    Bool_t SetDCAGammaGammaCut(Int_t DCAGammaGamma);
    Bool_t SetDCAZMesonPrimVtxCut(Int_t DCAZMesonPrimVtx);
    Bool_t SetDCARMesonPrimVtxCut(Int_t DCARMesonPrimVtx);
    void   SetOpeningAngleCut(Float_t OpeningAngle){fOpeningAngle = OpeningAngle;}
    Bool_t SetMinOpanMesonCut(Int_t minOpanMesonCut);
    Bool_t SetMaxOpanMesonCut(Int_t maxOpanMesonCut);
    void   SetEnableOpeningAngleCut (Bool_t isOn) {fEnableMinOpeningAngleCut = isOn;}
    void   SetIsMergedClusterCut(Int_t merged)                { fIsMergedClusterCut = merged; return;}
    void   SetUsePtDepSelectionWindow(Bool_t ptdep)                { fUsePtDepSelectionWindow = ptdep; return;}
    Bool_t GetUsePtDepSelectionWindow(Bool_t ptdep)                { return fUsePtDepSelectionWindow;}
    Int_t  GetIsMergedClusterCut()                            { return fIsMergedClusterCut;}
    Double_t GetRapidityCutValueMin()                            { return fRapidityCutMesonMin; }
    Double_t GetRapidityCutValueMax()                            { return fRapidityCutMesonMax; }
    void   SetEnableOmegaAPlikeCut(Bool_t DoOmegaAPlikeCut) {fEnableOmegaAPlikeCut = DoOmegaAPlikeCut;}

    Float_t FunctionMinMassCut(Float_t e);
    Float_t FunctionMaxMassCut(Float_t e);

    // Request Flags
    Bool_t   UseRotationMethod(){return fUseRotationMethodInBG;}
    Bool_t   UsePtmaxMethod(){return fUsePtmaxMethodForBG;}
    Bool_t   UseTrackMultiplicity(){return fUseTrackMultiplicityForBG;}
    Int_t    GetNumberOfBGEvents(){return fNumberOfBGEvents;}
    Int_t    NDegreesRotation(){return fNDegreeRotationPMForBG;}
    Bool_t   DoBGCalculation(){return fDoBG;}
    Bool_t   DoBGProbability(){return fDoBGProbability;}
    Bool_t   DoConvCaloMixing(){return fDoConvCaloMixing;}
    Bool_t   DoSectorMixing(){return fDoSectorMixing;}
    Bool_t   DoSectorJetMixing(){return fDoSectorJetMixing;}
    Bool_t   DoJetMixing(){return fDoJetMixing;}
    Bool_t   DoJetRotateMixing() {return fDoJetRotateMixing;}
    Bool_t   DoJetPtMixing() {return fDoJetPtMixing;}
    Bool_t   DoSphericityMixing(){return fDoSphericityMixing;}
    Int_t    DoGammaSwappForBg(){return fDoGammaSwappForBg;}
    Bool_t   DoWeightingInSwappBg(){return fDoWeightingInSwappBg;}
    Int_t    GammaSwappMethodBg(){return fGammaSwappMethodBg;}
    Int_t    GetNumberOfSwappsForBg(){return fNumberOfSwappsForBg;}
    Bool_t   DoJetAnalysis(){return fDoJetAnalysis;}
    Bool_t   DoJetQA(){return fDoJetQA;}
    Int_t    DoOutOfJet(){return fDoOutOfJet;}
    Bool_t   DoIsolatedAnalysis(){return fDoIsolatedAnalysis;}
    Bool_t   DoHighPtHadronAnalysis(){return fDoHighPtHadronAnalysis;}
    Bool_t   UseElecSharingCut(){return fDoSharedElecCut;}
    Bool_t   UseToCloseV0sCut(){return fDoToCloseV0sCut;}
    Bool_t   UseMCPSmearing(){if (fUseMCPSmearing > 0) return kTRUE; else return kFALSE;}
    Int_t    BackgroundHandlerType(){return fBackgroundHandler;}
    Double_t GetSelectionLow() const {return fSelectionLow;}
    Double_t GetSelectionHigh() const {return fSelectionHigh;}
    Double_t GetAcceptMassFlag() const {return fAcceptMesonMass;}
    Double_t GetMinPt() const {return fMinPt;}
    Bool_t   UseLikeSignMixing() {return fBackgroundUseLikeSign;}
    Bool_t   UseSidebandMixing() {return fBackgroundUseSideband;}
    Bool_t   UseSidebandMixingBothSides() {return fBackgroundUseSidebandBothSides;}
    Double_t GetSidebandMixingLow() const {return fSidebandMixingLow;}
    Double_t GetSidebandMixingHigh() const {return fSidebandMixingHigh;}
    Double_t GetSidebandMixingLeftLow() const {return fSidebandMixingLeftLow;}
    Double_t GetSidebandMixingLeftHigh() const {return fSidebandMixingLeftHigh;}
    Double_t GetSidebandMixingRightLow() const {return fSidebandMixingRightLow;}
    Double_t GetSidebandMixingRightHigh() const {return fSidebandMixingRightHigh;}
    Int_t    GetBackgroundMode() const {return fBackgroundMode;}
    Bool_t   DoGammaMinEnergyCut() const {return fDoGammaMinEnergyCut;}
    Int_t    GetNDaughterEnergyCut() const {return fNDaughterEnergyCut;}
    Int_t    GetSingleDaughterMinE() const {return fSingleDaughterMinE;}
    Bool_t   UseGammaSelection() const{return fUseGammaSelection;}

  protected:
    TRandom3    fRandom;                        ///<
    AliCaloPhotonCuts* fCaloPhotonCuts;         ///< CaloPhotonCutObject belonging to same main task
    TList*      fHistograms;                    ///< List of QA histograms
    TObjString* fCutString;                     ///< cut number used for analysis
    TString     fCutStringRead;
    // Histograms
    TH2F*       fHistoMesonCuts;                ///< bookkeeping for meson cuts
    TH2F*       fHistoMesonBGCuts;              ///< bookkeeping for meson bg cuts
    TH1F*       fHistoDCAGGMesonBefore;         ///<
    TH1F*       fHistoDCAZMesonPrimVtxBefore;   ///<
    TH1F*       fHistoDCARMesonPrimVtxBefore;   ///<
    TH1F*       fHistoDCAGGMesonAfter;          ///<
    TH2F*       fHistoDCAZMesonPrimVtxAfter;    ///<
    TH1F*       fHistoDCARMesonPrimVtxAfter;    ///<
    TH1F*       fHistoInvMassBefore;            ///<
    TH1F*       fHistoInvMassAfter;             ///<

    TF1*        fBrem;                          ///<
    TF1*        fFAlphaCut;                     ///<
    TF1*        fFMinOpanCut;                   ///<
    TF1*        fFMaxOpanCut;                   ///<

    Double_t    fMaxR;                          ///< max r cut
    Double_t    fMinPt;                         ///< min pT cut
    Double_t    fSelectionLow;                  ///< lower meson inv mass window for further selection
    Double_t    fSelectionHigh;                 ///< higher meson inv mass window for further selection
    Double_t    fSelectionNSigmaLow;            ///< N of sigma for ptdep selection window cut min
    Double_t    fSelectionNSigmaHigh;           ///< N of sigma for ptdep selection window cut max
    Double_t    fAlphaMinCutMeson;              ///< min value for meson alpha cut
    Double_t    fAlphaCutMeson;                 ///< max value for meson alpha cut
    Double_t    fRapidityCutMesonMin;           ///< min value for meson rapidity
    Double_t    fRapidityCutMesonMax;           ///< max value for meson rapidity
    Double_t    fMinV0Dist;                     ///
    UShort_t    fAPlikeSigma;                   ///< sigma range for the lower bound of the AP like cut
    Double_t    fMesonQualityMin;               ///
    Double_t    fPBremSmearing;                 ///
    Double_t    fPSigSmearing;                  ///
    Double_t    fPSigSmearingCte;               ///
    Double_t    fDCAGammaGammaCut;              ///< cut value for the maximum distance between the two photons [cm]
    Double_t    fDCAZMesonPrimVtxCut;           ///< cut value for the maximum distance in Z between the production point of the Meson & the primary vertex [cm]
    Double_t    fDCARMesonPrimVtxCut;           ///< cut value for the maximum distance in R between the production point of the Meson & the primary vertex [cm]
    Double_t    fMinOpanCutMeson;               ///<
    Double_t    fMaxOpanCutMeson;               ///<
    Double_t    fSidebandMixingLow;             ///<
    Double_t    fSidebandMixingHigh;            ///<
    Double_t    fSidebandMixingLeftLow;         ///<
    Double_t    fSidebandMixingLeftHigh;        ///<
    Double_t    fSidebandMixingRightLow;        ///<
    Double_t    fSidebandMixingRightHigh;       ///<

    Float_t     fOpeningAngle;                  ///< min opening angle for meson

    Int_t       fMode;                          ///< running mode of ConversionMesonCuts to select different sets of cut parameters for different running modes
    Int_t       fMesonKind;                     ///<
    Int_t       fIsMergedClusterCut;            ///< flag for merged cluster and di cluster analysis
    Int_t       fUsePtDepSelectionWindow;       ///< flag for usage of pT dependent selection window cut
    Int_t       fUseGammaSelection;             ///< flag for usage of gamma candidate selection e.g. gamma pairs which are in Pi0 mass region are not used as direct gammas for omega reco
    Int_t       fSelectionWindowCut;            ///< selection window for merged ana in mass
    Int_t       fNDegreeRotationPMForBG;        ///<
    Int_t       fNumberOfBGEvents;              ///<
    Int_t       fElectronLabelArraySize;        ///<
    Int_t*      fElectronLabelArray;            //[fElectronLabelArraySize] Array with elec/pos v0 label
    Int_t       fBackgroundHandler;             ///<
    Int_t       fMassParamFunction;             ///< flag to set the functions that should be used to paramterize NDM mass and width

    Int_t      fDoLightOutput;                 ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t      fDoMinPtCut;                    ///< do min pT cut
    Bool_t      fEnableMassCut;                 ///< flag to enable mass cut
    Bool_t      fAcceptMesonMass;               ///< flag to distinguish rejecting and accepting meson mass window for further analysis
    Bool_t      fUseRotationMethodInBG;         ///< flag to apply rotation method for meson bg estimation
    Bool_t      fUsePtmaxMethodForBG;           ///< flag to apply Ptmax method
    Bool_t      fDoBG;                          ///< flag to intialize BG
    Bool_t      fDoBGProbability;               ///< flag to use probability method for meson bg estimation
    Bool_t      fDoConvCaloMixing;              ///< flag to use enable convcalo mixing in addition to caloconv mixing
    Bool_t      fDoSectorMixing;                ///< flag to enable Sectormixing for meson bg estimation
    Bool_t      fDoSectorJetMixing;             ///< flag to enable Sectormixing with jets for meson bg estimation
    Bool_t      fDoJetMixing;                   ///< flag to enable mixing by cluster distance to jet axis
    Bool_t      fDoJetRotateMixing;             ///< flag to enable mixing by rotating calorimeter
    Bool_t      fDoJetPtMixing;                 ///< flag to enbale mixing by jet pt bins
    Bool_t      fDoSphericityMixing;            ///< flag to enable Sphericitymixing for meson bg estimation
    Bool_t      fUseTrackMultiplicityForBG;     ///< flag to use track multiplicity for meson bg estimation (else V0 mult)
    Int_t       fDoGammaSwappForBg;             ///< flag to use cluster swapping for background estimation
    Bool_t      fDoWeightingInSwappBg;          ///< flag to use multiplicity weighting for cluster swapping for background estimation
    Int_t       fGammaSwappMethodBg;            ///< flag to switch between different methods for cluster swapping: 0= 90 degree; 1=random angle
    Int_t       fNumberOfSwappsForBg;           ///< flag to enable multiple rotations for 1 photon pair for cluster swapping Bg
    Bool_t      fEnableMinOpeningAngleCut;      ///< flag to enable min opening angle cut
    Bool_t      fEnableOneCellDistCut;          ///< flag to enable 1 cell dist cut
    Bool_t      fAllowCombOnlyInSameRecMethod;  ///< flag to disable inv mass pairing among different calo's
    Bool_t      fDoToCloseV0sCut;               ///<
    Bool_t      fDoSharedElecCut;               ///<
    Bool_t      fDoMesonQualitySelection;       ///< flag to enable the meson selection based on the quality.
    Int_t       fUseMCPSmearing;                ///< flag
    Bool_t      fAlphaPtDepCut;                 ///<
    Bool_t      fDCAGammaGammaCutOn;            ///< cut flag for the maximum distance between the two photons
    Bool_t      fDCAZMesonPrimVtxCutOn;         ///< cut flag for the maximum distance in Z between the production point of the Meson & the primary vertex
    Bool_t      fDCARMesonPrimVtxCutOn;         ///< cut flag for the maximum distance in R between the production point of the Meson & the primary vertex
    Bool_t      fMinOpanPtDepCut;               ///<
    Bool_t      fMaxOpanPtDepCut;               ///<
    Bool_t      fBackgroundUseSideband;         ///< enable sideband mixing on one side of NDM
    Bool_t      fBackgroundUseSidebandBothSides;///< enable sideband mixing on both sides NDM
    Bool_t      fBackgroundUseLikeSign;         ///< enable likesign mixing
    Int_t       fBackgroundMode;                ///< default is 4: all pions from different event
    Bool_t      fDoJetAnalysis;                 ///< switch to run a jet analysis
    Bool_t      fDoJetQA;                       ///< switch to run a jet QA analysis
    Int_t       fDoOutOfJet;                    ///< switch to Analyse mesons out of jet (0 = switched off, 1 = all mesons out of jet, 2 = mesons on away side of jet, 3 = mesons in "donut shape" around jet)
    Bool_t      fDoIsolatedAnalysis;            ///< switch to run a isolated pi0 analysis
    Bool_t      fDoHighPtHadronAnalysis;        ///< switch to run a pi0 analysis with a high pt hadron in the event
    Bool_t      fEnableOmegaAPlikeCut;          ///< falg to enable the overloaded to close to V0 cut as cut inside an AP like plot

    Bool_t      fDoGammaMinEnergyCut;           ///< if enabled, at least fNDaughterEnergyCut daughter contributing to neutral meson need to fulfill fMinSingleDaughterE
    Int_t       fNDaughterEnergyCut;            ///< if above is enabled, at least fNDaughterEnergyCut daughter contributing to neutral meson needs to fulfill fMinSingleDaughterE
    Float_t     fSingleDaughterMinE;            ///< if above is enabled, at least fNDaughterEnergyCut daughter contributing to neutral meson needs to fulfill fMinSingleDaughterE

  private:

    /// \cond CLASSIMP
    ClassDef(AliConversionMesonCuts,46)
    /// \endcond
};


#endif
