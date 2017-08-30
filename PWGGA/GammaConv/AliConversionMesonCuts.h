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
 * |                  5                        | RCut                     |
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
      kRCut,
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
    virtual Bool_t CheckWhetherInMassRange(Double_t mass){ return mass>fSelectionLow && mass < fSelectionHigh;}

    TString GetCutNumber();

    // Cut Selection
    Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE, Double_t fRapidityShift=0., Int_t leadingCellID1 = 0, Int_t leadingCellID2 = 0);
    Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliMCEvent *mcEvent, Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedAODMC(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedMCDalitz(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma,Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedAODMCDalitz(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma,Double_t fRapidityShift=0.);
    Bool_t MesonIsSelectedMCEtaPiPlPiMiGamma(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCPiPlPiMiPiZero(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCPiZeroGamma(TParticle *fMCMother, AliMCEvent *mcEvent, Int_t &labelNeutPion, Int_t &labelGamma, Double_t fRapidityShift=0);
    Bool_t MesonIsSelectedMCChiC(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &, Int_t &, Int_t &, Double_t fRapidityShift=0. );
    Bool_t MesonIsSelectedPiZeroGammaAngle(AliAODConversionMother *omega, AliAODConversionMother *pi0, AliAODConversionPhoton *gamma, Bool_t DoPiZeroAngleCut, TF1 *maxfit, Double_t lowerFactor, Double_t upperFactor);
    void   PrintCuts();
    void   PrintCutsWithValues();
    
    void    SetLightOutput( Bool_t flag ){fDoLightOutput = flag; return;}
    void    SetRunningMode(Int_t mode){fMode = mode; return;}
    void    InitCutHistograms(TString name="",Bool_t additionalHists=kFALSE);
    void    SetFillCutHistograms(TString name=""){if(!fHistograms){InitCutHistograms(name);};}
    TList   *GetCutHistograms(){return fHistograms;}
    void    SmearParticle(AliAODConversionPhoton * photon);
    void    SmearVirtualPhoton(AliAODConversionPhoton* photon);
    TLorentzVector SmearElectron(TLorentzVector particle);

    void    SetDefaultSmearing(Double_t p0, Double_t p1, Double_t p2){fUseMCPSmearing=1.;fPBremSmearing=p0;fPSigSmearing=p1;fPSigSmearingCte=p2;return;}

    //Cut functions
    Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
    Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);

    void SetCaloMesonCutsObject(AliCaloPhotonCuts* cuts){fCaloPhotonCuts = cuts;}

    // Set Individual Cuts
    Bool_t SetRCut(Int_t RCut);
    Bool_t SetMesonKind(Int_t mesonKind);
    Bool_t SetSelectionWindowCut(Int_t selectionCut);
    Bool_t SetSelectionWindowMergedCut(Int_t selectionCut);
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
    Int_t  GetIsMergedClusterCut()                            { return fIsMergedClusterCut;}
    Double_t GetRapidityCutValue()                            { return fRapidityCutMeson; }

    Float_t FunctionMinMassCut(Float_t e);
    Float_t FunctionMaxMassCut(Float_t e);
    
    // Request Flags
    Bool_t   UseRotationMethod(){return fUseRotationMethodInBG;}
    Bool_t   UsePtmaxMethod(){return fUsePtmaxMethodForBG;}
    Bool_t   UseTrackMultiplicity(){return fUseTrackMultiplicityForBG;}
    Int_t    GetNumberOfBGEvents(){return fNumberOfBGEvents;}
    Int_t    NDegreesRotation(){return fnDegreeRotationPMForBG;}
    Bool_t   DoBGCalculation(){return fDoBG;}
    Bool_t   DoBGProbability(){return fdoBGProbability;}
    Bool_t   UseElecSharingCut(){return fDoSharedElecCut;}
    Bool_t   UseToCloseV0sCut(){return fDoToCloseV0sCut;}
    Bool_t   UseMCPSmearing(){return fUseMCPSmearing;}
    Int_t    BackgroundHandlerType(){return fBackgroundHandler;}
    Double_t GetSelectionLow() const {return fSelectionLow;}
    Double_t GetSelectionHigh() const {return fSelectionHigh;}
    
  protected:
    TList*    fHistograms;				 ///< List of QA histograms
    Bool_t    fDoLightOutput;             ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Int_t     fMode;                      ///< running mode of ConversionMesonCuts to select different sets of cut parameters for different running modes

    AliCaloPhotonCuts* fCaloPhotonCuts;   ///< CaloPhotonCutObject belonging to same main task

    //cuts
    Int_t     fMesonKind;				 ///<
    Int_t     fIsMergedClusterCut;        ///< flag for merged cluster and di cluster analysis
    Double_t  fMaxR;                      ///< r cut
    Bool_t    fEnableMassCut;             ///< flag to enable mass cut
    Double_t  fSelectionLow;              ///< lower meson inv mass window for further selection
    Double_t  fSelectionHigh;             ///< higher meson inv mass window for further selection
    Int_t     fSelectionWindowCut;        ///< selection window for merged ana in mass
    Double_t  fAlphaMinCutMeson;          ///< min value for meson alpha cut
    Double_t  fAlphaCutMeson;             ///< max value for meson alpha cut
    Double_t  fRapidityCutMeson;          ///< max value for meson rapidity
    Bool_t    fUseRotationMethodInBG;     ///< flag to apply rotation method for meson bg estimation
    Bool_t    fUsePtmaxMethodForBG;       ///< flag to apply Ptmax method
    Bool_t    fDoBG;                      ///< flag to intialize BG
    Bool_t    fdoBGProbability;           ///< flag to use probability method for meson bg estimation
    Bool_t    fUseTrackMultiplicityForBG; ///< flag to use track multiplicity for meson bg estimation (else V0 mult)
    Int_t     fnDegreeRotationPMForBG;    ///<
    Int_t     fNumberOfBGEvents;          ///<
    Float_t   fOpeningAngle;              ///< min opening angle for meson
    Bool_t    fEnableMinOpeningAngleCut;  ///< flag to enable min opening angle cut
    Bool_t    fEnableOneCellDistCut;      ///< flag to enable 1 cell dist cut
    Bool_t    fDoToCloseV0sCut;           //
    Double_t  fminV0Dist;                 //
    Bool_t    fDoSharedElecCut;           //
    Bool_t    fUseMCPSmearing;            // flag
    Double_t  fPBremSmearing;             //
    Double_t  fPSigSmearing;              //
    Double_t  fPSigSmearingCte;           //
    TF1*      fBrem;                      //
    TRandom3  fRandom;                    //
    TF1*      fFAlphaCut;                 //
    Bool_t    fAlphaPtDepCut;             //
    Int_t     fElectronLabelArraySize;    //
    Int_t*    fElectronLabelArray;        //[fElectronLabelArraySize] Array with elec/pos v0 label
    Double_t  fDCAGammaGammaCut;          ///< cut value for the maximum distance between the two photons [cm]
    Double_t  fDCAZMesonPrimVtxCut;       ///< cut value for the maximum distance in Z between the production point of the Meson & the primary vertex [cm]
    Double_t  fDCARMesonPrimVtxCut;       ///< cut value for the maximum distance in R between the production point of the Meson & the primary vertex [cm]
    Bool_t    fDCAGammaGammaCutOn;        ///< cut flag for the maximum distance between the two photons
    Bool_t    fDCAZMesonPrimVtxCutOn;     ///< cut flag for the maximum distance in Z between the production point of the Meson & the primary vertex
    Bool_t    fDCARMesonPrimVtxCutOn;     ///< cut flag for the maximum distance in R between the production point of the Meson & the primary vertex
    Double_t  fMinOpanCutMeson;           ///<
    TF1*      fFMinOpanCut;               ///<
    Bool_t    fMinOpanPtDepCut;           ///<
    Double_t  fMaxOpanCutMeson;           ///<
    TF1*      fFMaxOpanCut;               ///<
    Bool_t    fMaxOpanPtDepCut;           ///<
    Int_t     fBackgroundHandler;         ///<
    
    // Histograms
    TObjString* fCutString;                     ///< cut number used for analysis
    TString     fCutStringRead;
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

  private:

    /// \cond CLASSIMP
    ClassDef(AliConversionMesonCuts,19)
    /// \endcond
};


#endif
