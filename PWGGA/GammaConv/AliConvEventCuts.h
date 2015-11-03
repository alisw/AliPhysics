#ifndef ALICONVEVENTCUTS_H
#define ALICONVEVENTCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner                                    *

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "TRandom3.h"
#include "AliVCaloTrigger.h"

class AliESDEvent;
class AliAODEvent;
class TH1F;
class TH2F;
class TF1;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;
class AliEmcalTriggerPatchInfo;

using namespace std;

class AliConvEventCuts : public AliAnalysisCuts {
    
    public: 
      enum cutIds {
        kisHeavyIon,                  
        kCentralityMin,               
        kCentralityMax,               
        kSelectSpecialTriggerAlias,                 
        kSelectSubTriggerClass,             
        kremovePileUp,                
        kExtraSignals, 
        kVertex,
        kNCuts
      };

      enum TriggerTypeEMCAL {
        kND       = -1,  //not defined
        kJ1       = 1,
        kJ2       = 2,
        kG1       = 3,
        kG2       = 4,
        kL0       = 5,
      };

      AliConvEventCuts(const char *name="EventCuts", const char * title="Event Cuts");
      AliConvEventCuts(const AliConvEventCuts&);
      AliConvEventCuts& operator=(const AliConvEventCuts&);

      virtual ~AliConvEventCuts();                            //virtual destructor

  //    static AliConvEventCuts * GetStandardCuts2010PbPb();
  //    static AliConvEventCuts * GetStandardCuts2010pp();

      Int_t     fCuts[kNCuts];
      Bool_t    UpdateCutString();
      static const char * fgkCutNames[kNCuts];

      // Seters
      Bool_t    SetCutIds (TString cutString); 
      Bool_t    SetCut (cutIds cutID, Int_t cut);
      Bool_t    SetIsHeavyIon (Int_t isHeavyIon);
      Bool_t    SetCentralityMax (Int_t centralityBin);
      Bool_t    SetCentralityMin (Int_t centralityBin);
      Bool_t    SetRemovePileUp (Int_t removePileUp);  
      Bool_t    SetMultiplicityMethod (Int_t multiplicityMethod);
      Bool_t    SetSelectSpecialTrigger (Int_t selectSpecialTrigger);
      Bool_t    SetSelectSubTriggerClass (Int_t selectSpecialSubTriggerClass);
      Bool_t    SetRejectExtraSignalsCut (Int_t extraSignal);
      Bool_t    SetVertexCut(Int_t vertexCut);
      void    SetTriggerMimicking(Bool_t value)                                     { fMimicTrigger = value                                     ; 
                                                                                      if(value)AliInfo("enabled trigger mimicking")             ; }
      void    SetTriggerOverlapRejecion (Bool_t value)                              { fRejectTriggerOverlap = value                             ; 
                                                                                      if(value)AliInfo("enabled trigger overlap rejection")     ; }

      void    SetV0ReaderName (TString name)                                        { fV0ReaderName = name                                      ; }
      void    SetAddedSignalPDGCode (Int_t addedSignalPDGcode)                      { fAddedSignalPDGCode = addedSignalPDGcode                  ; }
      void    SetPreSelectionCutFlag (Bool_t preSelFlag)                            { fPreSelCut = preSelFlag                                   ; }   
      void    SetCaloTriggerPatchInfoName(const char *n)                            { fCaloTriggerPatchInfoName = n                             ; }
      void    SetCaloTriggersName(const char *n)                                    { fCaloTriggersName  = n                                    ; }
      void    SetAcceptedHeader(TList *HeaderList)                                  { fHeaderList = HeaderList                                  ; }   
      void    SetFillCutHistograms( TString name="",
                                    Bool_t preCut = kTRUE)                          { if(!fHistograms){ InitCutHistograms(name,preCut);}        ; }
      void    SetEtaShift(Double_t etaShift)                                        { fEtaShift = etaShift                                      ; } // Eta shift Setting
      void    SetEtaShift(TString pPbOrPbp)                                         { Double_t etaShift = 0.0                                   ;
                                                                                      if(!pPbOrPbp.CompareTo("pPb"))      etaShift = -0.465     ;
                                                                                      else if(!pPbOrPbp.CompareTo("Pbp")) etaShift =  0.465     ;
                                                                                      fEtaShift = etaShift                                      ; }
  //    void GetHistoCentralityFlattening
      void    SetUseWeightFlatCentralityFromFile( Int_t doFlattening = 1,
                              TString pathC="$ALICE_PHYSICS/PWGGA/GammaConv/InterpValuesAndFlattening.root",
                              TString histoCentNotFlat="")
                                                                                    { 
                                                                                      AliInfo(Form("enabled centrality flattening with weights from file: %s",pathC.Data()));
                                                                                      fDoCentralityFlat = doFlattening;
                                                                                      fPathWeightsFlatCent=pathC;
                                                                                      fNameHistoNotFlatCentrality = histoCentNotFlat;
                                                                                    }
      void    SetUseReweightingWithHistogramFromFile( Bool_t pi0reweight=kTRUE, 
                                Bool_t etareweight=kFALSE, 
                                Bool_t k0sreweight=kFALSE, 
                                                              TString path="$ALICE_PHYSICS/PWGGA/GammaConv/MCSpectraInput.root",
                                TString histoNamePi0 = "", 
                                TString histoNameEta = "", 
                                TString histoNameK0s = "",
                                TString fitNamePi0 = "", 
                                TString fitNameEta = "", 
                                TString fitNameK0s ="" ) 
                                                                                    {
                                                                                      AliInfo(Form("enabled reweighting for: pi0 : %i, eta: %i, K0s: %i",pi0reweight, etareweight, k0sreweight));
                                                                                      fDoReweightHistoMCPi0 = pi0reweight                       ; 
                                                                                      fDoReweightHistoMCEta = etareweight                       ; 
                                                                                      fDoReweightHistoMCK0s = k0sreweight                       ; 
                                                                                      fPathTrFReweighting=path                                  ;
                                                                                      fNameHistoReweightingPi0 =histoNamePi0                    ;
                                                                                      fNameHistoReweightingEta =histoNameEta                    ;
                                                                                      fNameHistoReweightingK0s =histoNameK0s                    ; 
                                                                                      fNameFitDataPi0 =fitNamePi0                               ;
                                                                                      fNameFitDataEta =fitNameEta                               ;
                                                                                      fNameFitDataK0s =fitNameK0s                               ; 
                                                                                    }
      void    SetMaxFacPtHard(Float_t value)                                        { fMaxFacPtHard = value                                     ; 
                                                                                      AliInfo(Form("maximum factor between pt hard and jet put to: %2.2f",fMaxFacPtHard));
                                                                                    }  
      
      // Geters
      TString   GetCutNumber();
      TString*  GetFoundHeader()                                                    { return fGeneratorNames                                    ; }
      Int_t     GetEventQuality()                                                   { return fEventQuality                                      ; }
      Bool_t    GetIsFromPileup()                                                   { return fRemovePileUp                                      ; }
      void      GetCentralityRange(Double_t range[2])                               { range[0]=10*fCentralityMin                                ;
                                                                                      range[1]=10*fCentralityMax                                ; }
      TList*    GetCutHistograms()                                                  { return fHistograms                                        ; }
      Int_t     GetMultiplicityMethod()                                             { return fMultiplicityMethod                                ; }
      Int_t     GetSignalRejection()                                                { return fRejectExtraSignals                                ; }
      Int_t     GetNAcceptedHeaders()                                               { return fnHeaders                                          ; }
      TString * GetAcceptedHeaderNames()                                            { return fGeneratorNames                                    ; }
      Int_t *   GetAcceptedHeaderStart()                                            { return fNotRejectedStart                                  ; }
      Int_t *   GetAcceptedHeaderEnd()                                              { return fNotRejectedEnd                                    ; }
      Int_t     GetAcceptedHeaderStart(Int_t headernumber)                          { if (headernumber < fnHeaders) 
                                                                                        return fNotRejectedStart[headernumber]                  ;
                                                                                      else 
                                                                                        return -1                                               ;
                                                                                    }
      Int_t     GetAcceptedHeaderEnd(Int_t headernumber)                            { if (headernumber < fnHeaders) 
                                                                                        return fNotRejectedEnd[headernumber]                    ; 
                                                                                      else 
                                                                                        return -1                                               ;
                                                                                    }
      TList*    GetAcceptedHeader()                                                 { return fHeaderList                                        ; }
      Int_t     GetNumberOfContributorsVtx(AliVEvent *event);
      Double_t  GetEtaShift()                                                       { return fEtaShift                                          ; }
      Bool_t    GetDoEtaShift()                                                     { return fDoEtaShift                                        ; }
      TString   GetSpecialTriggerName()                                             { return fSpecialTriggerName                                ; }
      AliEmcalTriggerPatchInfo   *GetMainTriggerPatch();
      ULong_t   GetTriggerList();
      Float_t   GetWeightForCentralityFlattening(AliVEvent *InputEvent = 0x0);
      Float_t   GetWeightForMeson(TString period, Int_t index, AliStack *MCStack, AliVEvent *InputEvent = 0x0);
      Float_t   GetCentrality(AliVEvent *event);
      void      GetCorrectEtaShiftFromPeriod(TString periodName);
      void      GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *MCEvent); 
      TClonesArray*     GetArrayFromEvent(AliVEvent* fInputEvent, const char *name, const char *clname=0);
      
      Bool_t    InitializeCutsFromCutString(const TString analysisCutSelection);
      void      SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kAny) {
                                                                                      fOfflineTriggerMask = offlineTriggerMask                  ;
                                                                                      fTriggerSelectedManually = kTRUE                          ;
                                                                                    }
      void    SelectSpecialTrigger( UInt_t offlineTriggerMask = AliVEvent::kAny, 
                                    TString TriggerClassName = "AliVEvent::kAny" ) {
                                                                                      fOfflineTriggerMask = offlineTriggerMask                  ;
                                                                                      fSpecialTriggerName = TriggerClassName                    ;
                                                                                      AliInfo(fSpecialTriggerName)                              ;
        
                                                                                    }   

      virtual   Bool_t IsSelected(TObject* /*obj*/)                                 { return kTRUE                                              ; }
      virtual   Bool_t IsSelected(TList* /*list*/)                                  { return kTRUE                                              ; }

      
      // Cut Selection
      Bool_t    EventIsSelected(  AliVEvent *fInputEvent, 
                                  AliVEvent *fMCEvent);
      Int_t     IsEventAcceptedByCut( AliConvEventCuts *ReaderCuts, 
                                      AliVEvent *InputEvent, 
                                      AliMCEvent *MCEvent, 
                                      Int_t isHeavyIon, 
                                      Bool_t isEMCALAnalysis);
        
      void    PrintCuts();
      void    PrintCutsWithValues();
      void    InitCutHistograms(  TString name="",
                                  Bool_t preCut = kTRUE);
      
      ///Cut functions
      Int_t   IsParticleFromBGEvent(  Int_t index, 
                                      AliStack *MCStack, 
                                      AliVEvent *InputEvent = 0x0);
      
      void    LoadWeightingFlatCentralityFromFile ();
      
      void    LoadReweightingHistosMCFromFile ();

      // Event Cuts
      Bool_t    IsCentralitySelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent = NULL);
      Bool_t    VertexZCut(AliVEvent *fInputEvent);
      Bool_t    IsJetJetMCEventAccepted(AliVEvent *MCEvent, Double_t& weight);
      Float_t   GetPtHard(AliVEvent *MCEvent);
      void      GetXSectionAndNTrials(AliVEvent *MCEvent, Float_t &XSection, Float_t &NTrials);
      Float_t   GetMaxPtJet()                                                       { return fMaxPtJetMC                                        ; }
      Bool_t    MimicTrigger( AliVEvent *fInputEvent, 
                              Bool_t isMC );
      Bool_t    IsTriggerSelected(  AliVEvent *fInputEvent, 
                                    Bool_t isMC);
      Bool_t    HasV0AND()                                                          { return fHasV0AND                                          ; }
      Bool_t    IsSDDFired()                                                        { return fIsSDDFired                                        ; }
      Int_t     IsSpecialTrigger()                                                  { return fSpecialTrigger                                    ; }
      Int_t     IsSpecialSubTrigger()                                               { return fSpecialSubTrigger                                 ; }
      void      InitializeEMCALTrigger( AliVEvent *fInputEvent);
      Bool_t    HasTriggerType(TriggerTypeEMCAL t);
      
      // Request Flags
      Int_t     IsHeavyIon()                                                        { return fIsHeavyIon                                        ; }
      void      DoEtaShift(Bool_t doEtaShift)                                       { fDoEtaShift = doEtaShift                                  ; }
      
      //MC particle flags - determine whether particle is primary or secondary
      Bool_t    IsConversionPrimaryESD( AliStack *MCStack,
                                        UInt_t stackpos, 
                                        Double_t prodVtxX, 
                                        Double_t prodVtxY,
                                        Double_t prodVtxZ);
      Bool_t    IsConversionPrimaryAOD( AliVEvent *fInputEvent, 
                                        AliAODMCParticle* AODMCParticle,
                                        Double_t prodVtxX, 
                                        Double_t prodVtxY,
                                        Double_t prodVtxZ);
      
      Int_t     SecondaryClassificationPhoton(  TParticle *particle,
                                                AliStack* fMCStack, 
                                                Bool_t isConversion );
      Int_t     SecondaryClassificationPhotonAOD( AliAODMCParticle *particle,
                                                  TClonesArray *aodmcArray, 
                                                  Bool_t isConversion );

      
    protected:
      TList*                      fHistograms;
      TList*                      fHeaderList;

      Int_t                       fEventQuality;                          // EventQuality
      //cuts
      Int_t                       fIsHeavyIon;                            // flag for heavy ion
      Int_t                       fDetectorCentrality;                    // centrality detecotor V0M or CL1
      Int_t                       fModCentralityClass;                    // allows to select smaller centrality classes
      Bool_t                      fEnableVertexCut;                       // enable vertex cut
      Double_t                    fMaxVertexZ;                            // max z offset of vertex
      Int_t                       fCentralityMin;                         // centrality selection lower bin value
      Int_t                       fCentralityMax;                         // centrality selection upper bin value
      Int_t                       fMultiplicityMethod;                    // selected multiplicity method
      Int_t                       fSpecialTrigger;                        // flag
      Int_t                       fSpecialSubTrigger;                     // flag
      Bool_t                      fRemovePileUp;                          // flag
      Int_t                       fRejectExtraSignals;                    //
      UInt_t                      fOfflineTriggerMask;                    // Task processes collision candidates only
      Bool_t                      fHasV0AND;                              // V0AND Offline Trigger
      Bool_t                      fIsSDDFired;                            // SDD FIRED to select with SDD events
      TRandom3                    fRandom;                                //
      Int_t                       fnHeaders;                              // Number of Headers
      Int_t*                      fNotRejectedStart;                      //[fnHeaders]
      Int_t*                      fNotRejectedEnd;                        //[fnHeaders]
      TString*                    fGeneratorNames;                        //[fnHeaders]
      TObjString*                 fCutString;                             // cut number used for analysis
      AliAnalysisUtils*           fUtils;
      Double_t                    fEtaShift;
      Bool_t                      fDoEtaShift;                            // Flag for Etashift
      Int_t                       fDoCentralityFlat;
      TString                     fNameHistoNotFlatCentrality;
      TString                     fPathWeightsFlatCent;
      Bool_t                      fDoReweightHistoMCPi0;                  // Flag for reweighting Pi0 input with histogram
      Bool_t                      fDoReweightHistoMCEta;                  // Flag for reweighting Eta input with histogram
      Bool_t                      fDoReweightHistoMCK0s;                  // Flag for reweighting K0s input with histogram
      TString                     fPathTrFReweighting;                    // Path for file used in reweighting
      TString                     fNameHistoReweightingPi0;               // Histogram name for reweighting Pi0
      TString                     fNameHistoReweightingEta;               // Histogram name for reweighting Eta
      TString                     fNameHistoReweightingK0s;               // Histogram name for reweighting K0s
      TString                     fNameFitDataPi0;                        // Fit name for fit to spectrum of pi0s in Data
      TString                     fNameFitDataEta;                        // Fit name for fit to spectrum of etas in Data
      TString                     fNameFitDataK0s;                        // Fit name for fit to spectrum of k0s in Data
      // Histograms
      TH1F*                       fHistoEventCuts;                        // bookkeeping for event selection cuts
      TH1F*                       hCentrality;                            // centrality distribution for selected events
      TH1D*                       hCentralityNotFlat;                     // centrality distribution loaded for cent. flattening
      //TH2F*                      hCentralityVsNumberOfPrimaryTracks;    // centrality distribution for selected events
      TH1F*                       hVertexZ;                               // vertex z distribution for selected events
      TH1F*                       hTriggerClass;                          // fired offline trigger class
      TH1F*                       hTriggerClassSelected;                  // selected fired offline trigger class
      TH1F*                       hTriggerClassesCorrelated;              // selected trigger class correlation with others
      TH1D*                       hReweightMCHistPi0;                     // histogram input for reweighting Pi0
      TH1D*                       hReweightMCHistEta;                     // histogram input for reweighting Eta
      TH1D*                       hReweightMCHistK0s;                     // histogram input for reweighting K0s
      TF1*                        fFitDataPi0;                            // fit to pi0 spectrum in Data
      TF1*                        fFitDataEta;                            // fit to eta spectrum in Data
      TF1*                        fFitDataK0s;                            // fit to K0s spectrum in Data
      Int_t                       fAddedSignalPDGCode;
      Bool_t                      fPreSelCut;                             // Flag for preselection cut used in V0Reader
      Bool_t                      fTriggerSelectedManually;               // Flag for manual trigger selection
      TString                     fSpecialTriggerName;                    // Name of the Special Triggers
      TString                     fSpecialSubTriggerName;                 // Name of the Special Triggers
      Int_t                       fNSpecialSubTriggerOptions;
      TH2F*                       hSPDClusterTrackletBackground;          // SPD tracklets vs SPD clusters for background-correction
      // trigger information
      TString                fV0ReaderName;                               // Name of V0Reader
      AliVCaloTrigger*            fCaloTriggers;                          //! calo triggers
      TClonesArray*               fTriggerPatchInfo;                      //! trigger patch info array
      AliEmcalTriggerPatchInfo *  fMainTriggerPatchEMCAL;                 // main trigger patch, will be cached after first call
      TString                     fCaloTriggersName;                      // name of calo triggers collection
      TString                     fCaloTriggerPatchInfoName;              // trigger patch info array name
      ULong_t                     fTriggersEMCAL;                         // list of fired EMCAL triggers
      ULong_t                     fTriggersEMCALSelected;                 // list of accepted triggers
      Bool_t                      fEMCALTrigInitialized;                  // EMCAL triggers initialized
      // Primary secondary distinction
      Double_t                    fSecProdBoundary;                       // 3D radius of production (cm) for primary-secodary distinction
      Float_t                     fMaxPtJetMC;                            // maximum jet pt in event
      Float_t                     fMaxFacPtHard;                          // maximum factor between maximum jet pt and pt hard generated
      Float_t                     fMaxFacPtHardSingleParticle;            // maximum factor between maximum single particle pt (pi0/eta) and pt hard generated
      Bool_t                      fMimicTrigger;                          // enable trigger mimiking
      Bool_t                      fRejectTriggerOverlap;                  // enable trigger overlap rejections
    private:

      ClassDef(AliConvEventCuts,15)
};


#endif
