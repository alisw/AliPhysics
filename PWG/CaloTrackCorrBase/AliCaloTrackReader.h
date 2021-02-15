#ifndef ALICALOTRACKREADER_H
#define ALICALOTRACKREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloTrackReader
/// \ingroup CaloTrackCorrelationsBase
/// \brief Base class for event, clusters and tracks filtering and preparation for the analysis.
///
/// Base class for accessing/filtering data, MonteCarlo, ESD or AOD, of PHOS, EMCAL and 
/// Central Barrel Tracking detectors. It filters de events and detector input
/// depending on different selection criteria, like kinematical restrictions, 
/// goodness of the event or cluster/track, etc.
/// 
/// Mother class of 
///  * AliCaloTrackESDReader: Fills ESD data in 3 TObjArrays (PHOS, EMCAL, CTS)
///  * AliCaloTrackMCReader : Fills Kinematics data in 3 TObjArrays (PHOS, EMCAL, CTS)
///  * AliCaloTrackAODReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS)   
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include <TObject.h> 
#include <TString.h>
class TObjArray ; 
class TTree ;
class TArrayI ;
class TObjString;
#include <TRandom3.h>

//--- ANALYSIS system ---
#include "AliVEvent.h"
class AliVCaloCells;
class AliHeader; 
class AliGenEventHeader; 
class AliGenPythiaEventHeader; 
class AliAODEvent;
class AliMCEvent;
class AliMixedEvent;
class AliAODMCHeader;
class AliCentrality;
class AliMultSelection;
class AliESDtrackCuts;
//class AliTriggerAnalysis;
class AliEventplane;
class AliVCluster;
#include "AliLog.h"
#include "AliEventCuts.h"
//#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

// --- CaloTrackCorr / EMCAL ---
#include "AliFiducialCut.h"
class AliCalorimeterUtils;
#include "AliAnaWeights.h"
#include "AliMCAnalysisUtils.h"

class AliCaloTrackReader : public TObject {

public: 
  
                  AliCaloTrackReader() ; // ctor
  virtual        ~AliCaloTrackReader() ; // virtual dtor
  void            DeletePointers();
  
  //--------------------------------
  // General methods
  //--------------------------------

  virtual void    Init();
    
  virtual void    InitParameters();
  
  virtual void    Print(const Option_t * opt) const;
  
  virtual void    ResetLists();

  virtual Int_t   GetDebug()                         const { return fDebug                 ; }
  virtual void    SetDebug(Int_t d)                        { fDebug = d                    ; }
  
  enum inputDataType {kESD, kAOD, kMC};
  virtual Int_t   GetDataType()                      const { return fDataType              ; }
  virtual void    SetDataType(Int_t data )                 { fDataType = data              ; }

  virtual Int_t   GetEventNumber()                   const { return fEventNumber           ; }
	
  virtual TObjString *  GetListOfParameters() ;
  
  TString         GetTaskName()                      const { return fTaskName              ; }
  void            SetTaskName(TString name)                { fTaskName = name              ; }
    
  //---------------------------------------
  // Input/output event setters and getters
  //---------------------------------------
  
  virtual void    SetInputEvent(AliVEvent* input) ;
  virtual void    SetOutputEvent(AliAODEvent*  aod)        { fOutputEvent = aod            ; }
  virtual void    SetMC(AliMCEvent* mc) ;           
  virtual void    SetInputOutputMCEvent(AliVEvent* /*esd*/, AliAODEvent* /*aod*/, AliMCEvent* /*mc*/) { ; }
  
  //
  // Embedded events
  //
  /// Reject clusters without MC label (reject background)
  Bool_t           IsEmbeddedClusterSelectionOn()    const { return fSelectEmbeddedClusters   ; }
  void             SwitchOnEmbeddedClustersSelection()     { fSelectEmbeddedClusters = kTRUE  ; }
  void             SwitchOffEmbeddedClustersSelection()    { fSelectEmbeddedClusters = kFALSE ; }
  
  /// Use as input for the analysis the MCEvent() or the InputEvent() from embedded MC signal.
  /// \param useMCEvt: recover not the standard MCEvent() but an external MC embedded event
  /// \param useInputEvt: recover not the standard InputEvent() but an external embedded input event
  void         UseEmbeddedEvent(Bool_t useMCEvt, Bool_t useInputEvt) { 
    fEmbeddedEvent[0] = useMCEvt ; fEmbeddedEvent[1] = useInputEvt ; }

  Bool_t        IsEmbeddedMCEventUsed   () { return fEmbeddedEvent[0] ; }
  Bool_t        IsEmbeddedInputEventUsed() { return fEmbeddedEvent[1] ; }
  
//  AliVCaloCells * GetEMCALCellsExternalEvent() { 
//    return AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent()->GetEMCALCells() ; }
//  AliVCluster   * GetCaloClusterExternalEvent(Int_t icluster) { 
//    return AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent()->GetCaloCluster(icluster) ; }
//  Int_t           GetNumberOfCaloClustersExternalEvent() { 
//     return AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent()->GetNumberOfCaloClusters() ; }
//  
  // Delta AODs
  
  virtual TList * GetAODBranchList()                 const { return fAODBranchList         ; }
  void            SetDeltaAODFileName(TString name)        { fDeltaAODFileName = name      ; }
  TString         GetDeltaAODFileName()              const { return fDeltaAODFileName      ; }
  void            SwitchOnWriteDeltaAOD()                  { fWriteOutputDeltaAOD = kTRUE  ; }
  void            SwitchOffWriteDeltaAOD()                 { fWriteOutputDeltaAOD = kFALSE ; }
  Bool_t          WriteDeltaAODToFile()              const { return fWriteOutputDeltaAOD   ; } 
  
  virtual TList * GetCreateControlHistograms() ;
  void            SetControlHistogramEnergyBinning(Int_t nBins, Float_t emin, Float_t emax)
  { fEnergyHistogramNbins = nBins ; fEnergyHistogramLimit[0] = emin; fEnergyHistogramLimit[1] = emax ; }
  
  Bool_t          IsHistoCentDependentOn()           const {   return fHistoCentDependent  ; }
  void            SwitchOffHistoCentDependent()            {   fHistoCentDependent = kFALSE; }
  void            SwitchOnHistoCentDependent()             {   fHistoCentDependent = kTRUE ; }
  
  void            SwitchOffHistoPtDependent()              {   fHistoPtDependent = kFALSE  ; }
  void            SwitchOnHistoPtDependent()               {   fHistoPtDependent = kTRUE   ; }
  
  //------------------------------------------------------------
  // Clusters/Tracks arrays filtering/filling methods and switchs 
  //------------------------------------------------------------
  
  // detector identificator enum, used here and in AliAnaCaloTrackBaseClass and derived classes
  enum detector { kEMCAL = AliFiducialCut::kEMCAL, kPHOS = AliFiducialCut::kPHOS,
                  kCTS   = AliFiducialCut::kCTS  , kDCAL = AliFiducialCut::kDCAL,
                  kDCALPHOS = AliFiducialCut::kDCALPHOS } ;
  
  /// Smearing function enum.
  enum smearingFunction     { kNoSmearing, kSmearingLandau, kSmearingLandauShift           } ;

  // Minimum pt setters and getters
  
  Float_t          GetEMCALPtMin()                   const { return fEMCALPtMin            ; }
  Float_t          GetPHOSPtMin()                    const { return fPHOSPtMin             ; }
  Float_t          GetCTSPtMin()                     const { return fCTSPtMin              ; }
  Float_t          GetEMCALPtMax()                   const { return fEMCALPtMax            ; }
  Float_t          GetPHOSPtMax()                    const { return fPHOSPtMax             ; }
  Float_t          GetCTSPtMax()                     const { return fCTSPtMax              ; }
  
  void             SetEMCALPtMin(Float_t  pt)              { fEMCALPtMin = pt              ; }
  void             SetPHOSPtMin (Float_t  pt)              { fPHOSPtMin  = pt              ; }
  void             SetCTSPtMin  (Float_t  pt)              { fCTSPtMin   = pt              ; }  
  void             SetEMCALPtMax(Float_t  pt)              { fEMCALPtMax = pt              ; }
  void             SetPHOSPtMax (Float_t  pt)              { fPHOSPtMax  = pt              ; }
  void             SetCTSPtMax  (Float_t  pt)              { fCTSPtMax   = pt              ; }  
   
  Float_t          GetEMCALEMin()                    const { return GetEMCALPtMin()        ; }
  Float_t          GetPHOSEMin()                     const { return GetPHOSPtMin()         ; }
  Float_t          GetEMCALEMax()                    const { return GetEMCALPtMax()        ; }
  Float_t          GetPHOSEMax()                     const { return GetPHOSPtMax()         ; }
  
  void             SetEMCALEMin (Float_t  en)              { SetEMCALPtMin(en)             ; }
  void             SetPHOSEMin  (Float_t  en)              { SetPHOSPtMin (en)             ; }
  void             SetEMCALEMax (Float_t  en)              { SetEMCALPtMax(en)             ; }
  void             SetPHOSEMax  (Float_t  en)              { SetPHOSPtMax (en)             ; }
  
  virtual Int_t    GetTrackID(AliVTrack* track) ;
  
  // Distance to bad channels cut
  
  Float_t          GetEMCALBadChannelMinDist()       const { return fEMCALBadChMinDist     ; }
  Float_t          GetPHOSBadChannelMinDist()        const { return fPHOSBadChMinDist      ; }

  void             SetEMCALBadChannelMinDist(Float_t di)   { fEMCALBadChMinDist = di       ; }
  void             SetPHOSBadChannelMinDist (Float_t di)   { fPHOSBadChMinDist  = di       ; }

  // Number of cells in cluster cut
  
  Int_t            GetEMCALNCellsCut()               const { return fEMCALNCellsCut       ; }
  Int_t            GetPHOSNCellsCut()                const { return fPHOSNCellsCut        ; }
    
  void             SetEMCALNCellsCut(Int_t nc)             { fEMCALNCellsCut = nc         ; }
  void             SetPHOSNCellsCut (Int_t nc)             { fPHOSNCellsCut  = nc         ; }
  
  Float_t          GetEMCALEnDepNCellsCutEnMin()  const { return fEMCALNCellsCutEnDepEnMin; }
  Float_t          GetEMCALEnDepNCellsCutConst()  const { return fEMCALNCellsCutEnDepConstant; }
  Float_t          GetEMCALEnDepNCellsCutSlope()  const { return fEMCALNCellsCutEnDepSlope; }
  void             SetEMCALEnDepNCellsCut(Float_t emin, Float_t constant, Float_t slope)             
  { fEMCALNCellsCutEnDepEnMin = emin ;  fEMCALNCellsCutEnDepConstant = constant; fEMCALNCellsCutEnDepSlope = slope ;}
  
  // Track DCA cut
  
  Bool_t           AcceptDCA(Float_t pt, Float_t dca);
  Double_t         GetTrackDCACut(Int_t i)           const { if(i >= 0 && i < 3 ) return fTrackDCACut[i] ;
                                                             else return -999              ; }
  
  void             SetTrackDCACut(Int_t i, Float_t cut)    { if(i >= 0 && i < 3 )
                                                             fTrackDCACut[i] = cut         ; }
  
  void             SwitchOnUseTrackDCACut()                { fUseTrackDCACut = kTRUE       ; }
  void             SwitchOffUseTrackDCACut()               { fUseTrackDCACut = kFALSE      ; }
  Bool_t           IsDCACutOn()                      const { return fUseTrackDCACut        ; }
  
  // Time cut
  
  Double_t         GetTrackTimeCutMin()              const { return fTrackTimeCutMin       ; }
  Double_t         GetTrackTimeCutMax()              const { return fTrackTimeCutMax       ; }
  
  void             SetTrackTimeCut(Double_t a, Double_t b) { fTrackTimeCutMin = a ;
                                                             fTrackTimeCutMax = b          ; } // ns
  
  void             SwitchOnUseTrackTimeCut()               { fUseTrackTimeCut = kTRUE      ;  fAccessTrackTOF  = kTRUE ; }
  void             SwitchOffUseTrackTimeCut()              { fUseTrackTimeCut = kFALSE     ; }

  void             SwitchOnAccessTrackTimeCut()            { fAccessTrackTOF  = kTRUE      ; }
  void             SwitchOffAccessTrackTimeCut()           { fAccessTrackTOF  = kFALSE     ; }
  Bool_t           IsAccessToTrackTimeOn()           const { return fAccessTrackTOF        ; }

  
  Double_t         GetEMCALTimeCutMin()              const { return fEMCALTimeCutMin       ; }
  Double_t         GetEMCALTimeCutMax()              const { return fEMCALTimeCutMax       ; }	

  Bool_t           IsInTimeWindow(Double_t tof, Float_t energy)  const ;
  
  void             SetEMCALTimeCut(Double_t a, Double_t b) { fEMCALTimeCutMin = a ; 
                                                             fEMCALTimeCutMax = b          ; } // ns
  
  void             SetEMCALParametrizedMinTimeCut(Int_t i, Float_t par) { fEMCALParamTimeCutMin[i] = par ; } 
  void             SetEMCALParametrizedMaxTimeCut(Int_t i, Float_t par) { fEMCALParamTimeCutMax[i] = par ; }
  
  void             SwitchOnUseEMCALTimeCut()               { fUseEMCALTimeCut = kTRUE      ; }
  void             SwitchOffUseEMCALTimeCut()              { fUseEMCALTimeCut = kFALSE     ; }
  
  void             SwitchOnUseParametrizedTimeCut()        { fUseParamTimeCut = kTRUE      ; }
  void             SwitchOffUseParametrizedTimeCut()       { fUseParamTimeCut = kFALSE     ; }

  Float_t          GetEMCALHighEnergyNdiffCut()     const  { return fEMCALHighEnergyNdiffCut;}
  void             SetEMCALHighEnergyNdiffCut(Float_t en)  { fEMCALHighEnergyNdiffCut = en ; }
 
  Float_t          GetEMCALMinCellEnNdiffCut()      const  { return fEMCALMinCellEnNdiffCut; }
  void             SetEMCALMinCellEnNdiffCut(Float_t en)   { fEMCALMinCellEnNdiffCut  = en ; }
  
  // Fidutial cuts
  
  virtual AliFiducialCut * GetFiducialCut()                { 
                    if(!fFiducialCut) fFiducialCut = new AliFiducialCut(); 
                    return  fFiducialCut                                                   ; }
  virtual void     SetFiducialCut(AliFiducialCut * fc)     { fFiducialCut = fc           ; }
  virtual Bool_t   IsFiducialCutOn()                 const { return fCheckFidCut           ; }
  virtual void     SwitchOnFiducialCut()                   { fCheckFidCut = kTRUE          ; 
                                                             fFiducialCut = new AliFiducialCut() ; }
  virtual void     SwitchOffFiducialCut()                  { fCheckFidCut = kFALSE         ; }
    
  // Cluster/track/cells switchs
  
  Bool_t           IsCTSSwitchedOn()                 const { return fFillCTS               ; }
  void             SwitchOnCTS()                           { fFillCTS = kTRUE              ; }
  void             SwitchOffCTS()                          { fFillCTS = kFALSE             ; }

  Bool_t           IsEMCALSwitchedOn()               const { return fFillEMCAL             ; }
  void             SwitchOnEMCAL()                         { fFillEMCAL = kTRUE            ; }
  void             SwitchOffEMCAL()                        { fFillEMCAL = kFALSE           ; }

  Bool_t           IsDCALSwitchedOn()                const { return fFillDCAL              ; }
  void             SwitchOnDCAL()                          { fFillDCAL = kTRUE             ; }
  void             SwitchOffDCAL()                         { fFillDCAL = kFALSE            ; }
  
  Bool_t           IsPHOSSwitchedOn()                const { return fFillPHOS              ; }
  void             SwitchOnPHOS()                          { fFillPHOS = kTRUE             ; }
  void             SwitchOffPHOS()                         { fFillPHOS = kFALSE            ; }

  Bool_t           IsEMCALCellsSwitchedOn()          const { return fFillEMCALCells        ; }
  void             SwitchOnEMCALCells()                    { fFillEMCALCells = kTRUE       ; }
  void             SwitchOffEMCALCells()                   { fFillEMCALCells = kFALSE      ; }

  Bool_t           IsPHOSCellsSwitchedOn()           const { return fFillPHOSCells         ; }
  void             SwitchOnPHOSCells()                     { fFillPHOSCells = kTRUE        ; }
  void             SwitchOffPHOSCells()                    { fFillPHOSCells = kFALSE       ; }

  Bool_t           AreClustersRecalculated()         const { return fRecalculateClusters   ; }
  void             SwitchOnClusterRecalculation()          { fRecalculateClusters = kTRUE  ; }
  void             SwitchOffClusterRecalculation()         { fRecalculateClusters = kFALSE ; }  

  void             SwitchOnClusterELinearityCorrection()   { fCorrectELinearity = kTRUE    ; }
  void             SwitchOffClusterELinearityCorrection()  { fCorrectELinearity = kFALSE   ; }

  void             SwitchOnClusterEScalePerSMCorrection()  { fScaleEPerSM = kTRUE          ; }
  void             SwitchOffClusterEScalePerSMCorrection() { fScaleEPerSM = kFALSE         ; }
  void             SetScaleFactorPerSM(Int_t ism, Float_t factor)          
                                                           { if ( ism < 22 && ism >= 0 ) fScaleFactorPerSM[ism] = factor ; }

  // Shower shape smearing function
  
  void             SetSmearingFunction(Int_t smfu)         { fSmearingFunction = smfu      ; }
  Int_t            GetSmearingFunction()             const { return fSmearingFunction      ; }
  
  Bool_t           IsShowerShapeSmeared()            const { return fSmearShowerShape      ; }
  void             SwitchOnShowerShapeSmearing()           { fSmearShowerShape = kTRUE     ; }
  void             SwitchOffShowerShapeSmearing()          { fSmearShowerShape = kFALSE    ; }
  
  void             SetShowerShapeSmearWidth(Float_t w )    { fSmearShowerShapeWidth = w    ; }
  
  void             SetSmearingNLMRange(Int_t mi, Int_t ma) { fSmearNLMMin = mi ; fSmearNLMMax = ma ; }
  
  // Filling/ filtering / detector information access methods
  
  virtual Bool_t   FillInputEvent(Int_t iEntry, const char *currentFileName)  ;
  virtual void     FillInputCTS() ;
  virtual void     FillInputCTSSelectTrack(AliVTrack * track, Int_t itrack, Bool_t & bc0) ;
  virtual void     FillInputEMCAL() ;
  virtual void     FillInputEMCALSelectCluster(AliVCluster * clus, Int_t iclus) ;
  virtual void     FillInputPHOS() ;
  virtual void     FillInputEMCALCells() ;
  virtual void     FillInputPHOSCells() ;
  virtual void     FillInputVZERO() ;  
  
  Int_t            GetV0Signal(Int_t i)              const { return fV0ADC[i]               ; }
  Int_t            GetV0Multiplicity(Int_t i)        const { return fV0Mul[i]               ; }
  
  void             SetEMCALClusterListName(TString name)   { fEMCALClustersListName = name  ; }
  TString          GetEMCALClusterListName()         const { return fEMCALClustersListName  ; }

  void             SetEMCALCellsListName(TString name)     { fEMCALCellsListName = name     ; }
  TString          GetEMCALCellsListName()           const { return fEMCALCellsListName     ; }
  
  // Arrays with clusters/track/cells access method
  
  virtual TObjArray*     GetCTSTracks()              const { return fCTSTracks              ; }
  virtual TObjArray*     GetEMCALClusters()          const { return fEMCALClusters          ; }
  virtual TObjArray*     GetDCALClusters()           const { return fDCALClusters           ; }
  virtual TObjArray*     GetPHOSClusters()           const { return fPHOSClusters           ; }
  virtual AliVCaloCells* GetEMCALCells()             const { return fEMCALCells             ; }
  virtual AliVCaloCells* GetPHOSCells()              const { return fPHOSCells              ; }
  
  //-------------------------------------
  // Event/track selection methods
  //-------------------------------------
  
  void             AcceptFastClusterEvents()               { fAcceptFastCluster     = kTRUE  ; } 
  void             RejectFastClusterEvents()               { fAcceptFastCluster     = kFALSE ; }  
  Bool_t           IsFastClusterAccepted()           const { return fAcceptFastCluster       ; }   
  
  Bool_t           AcceptEventWithTriggerBit(UInt_t trigFired);
  Bool_t           RejectEventWithTriggerBit(UInt_t trigFired);
  void             SetAcceptEventsWithBit(UInt_t bit)      { Int_t n = fAcceptEventsWithBit.GetSize();
                                                             fAcceptEventsWithBit.Set(n+1);
                                                             fAcceptEventsWithBit.AddAt(bit,n) ; }
  
  void             SetRejectEventsWithBit(UInt_t bit)      { Int_t n = fRejectEventsWithBit.GetSize();
                                                             fRejectEventsWithBit.Set(n+1);
                                                             fRejectEventsWithBit.AddAt(bit,n) ; }
  /// Activate removal of LED events depending on number of cells in SM
  /// \param opt: 1- default, check only SM3, 2- or larger check all SMs
  Bool_t           IsLEDEventRemoved()               const { return fRemoveLEDEvents         ; }   
  Bool_t           RejectLEDEvents();
  
  void             SwitchOnLEDEventsRemoval(Int_t opt = 1) { fRemoveLEDEvents       = opt    ; }
  void             SwitchOffLEDEventsRemoval()             { fRemoveLEDEvents       = 0      ; }
  
  void             SetLEDHighEnergyCutSM(Float_t e)        { fLEDHighEnergyCutSM         = e ; }   
  void             SetLEDHighNCellsCutSM(Int_t   n)        { fLEDHighNCellsCutSM         = n ; }   
  void             SetLEDLowEnergyCutSM3(Float_t e)        { fLEDLowEnergyCutSM3         = e ; }   
  void             SetLEDLowNCellsCutSM3(Int_t   n)        { fLEDLowNCellsCutSM3         = n ; }   
  void             SetLEDMinCellEnergy  (Float_t e)        { fLEDMinCellEnergy           = e ; }   
  void             SetLEDMaxCellEnergy  (Float_t e)        { fLEDMaxCellEnergy           = e ; }  
  
  void             SwitchOnLEDStripEventsRemoval()         { fRemoveLEDStripEvents  = kTRUE  ; }
  void             SwitchOffLEDStripEventsRemoval()        { fRemoveLEDStripEvents  = kFALSE ; }
  
  void             SetLEDStripHighEnergyCutSM(Float_t eFull, Float_t eThird)   
                                                           { fLEDHighEnergyCutStrip[0] = eFull ; fLEDHighEnergyCutStrip[1] = eThird ; }   
  void             SetLEDStripHighNCellsCutSM(Int_t   nFull, Int_t   nThird) 
                                                           { fLEDHighNCellsCutStrip[0] = nFull ; fLEDHighNCellsCutStrip[0] = nThird ; }   
  void             SetLEDStripLowEnergyCutSM3(Float_t e)   { fLEDLowEnergyCutSM3Strip    = e ; }   
  void             SetLEDStripLowNCellsCutSM3(Int_t   n)   { fLEDLowNCellsCutSM3Strip    = n ; }   
  void             SetLEDEventMaxNumberOfStrips(Int_t n)   { fLEDEventMaxNumberOfStrips  = n ; }   

  void             SetFiredTriggerClassName(TString name)  { fFiredTriggerClassName = name   ; }
  TString          GetFiredTriggerClassName()        const { return fFiredTriggerClassName   ; }
  TString          GetFiredTriggerClasses()          const { return GetInputEvent()->GetFiredTriggerClasses() ; }
  
  
  // Event selection when mixed event is used
  
  UInt_t           GetEventTriggerMaskInput()         const { return fEventTriggerMaskInput    ; }
  void             SetEventTriggerMaskInput(UInt_t evtTrig = AliVEvent::kAny) 
                                                           { fEventTriggerMaskInput = evtTrig ; }
  
  UInt_t           GetEventTriggerMask()             const { return fEventTriggerMask        ; }
  void             SetEventTriggerMask(UInt_t evtTrig = AliVEvent::kAny) 
                                                           { fEventTriggerMask = evtTrig     ; }
  UInt_t           GetMixEventTriggerMask()          const { return fMixEventTriggerMask     ; }
  void             SetMixEventTriggerMask(UInt_t evtTrig = AliVEvent::kAnyINT)
                                                           { fMixEventTriggerMask = evtTrig  ; }
  Bool_t           IsEventTriggerAtSEOn()            const { return fEventTriggerAtSE        ; }
  void             SwitchOnEventTriggerAtSE()              { fEventTriggerAtSE      = kTRUE  ; }
  void             SwitchOffEventTriggerAtSE()             { fEventTriggerAtSE      = kFALSE ; }
		
  // EMCal Triggered events selection, studies
  
  TArrayI          GetTriggerPatches(Int_t tmin, Int_t tmax);
  void             MatchTriggerCluster(TArrayI patches);
  void             SetEMCALTriggerThresholds();
  
  Bool_t           CheckEventTriggers();
  
  Bool_t           IsExoticEvent()                   const { return fIsExoticEvent           ; }
  Bool_t           IsBadCellTriggerEvent()           const { return fIsBadCellEvent          ; }
  Bool_t           IsBadMaxCellTriggerEvent()        const { return fIsBadMaxCellEvent       ; }
  Bool_t           IsTriggerMatched()                const { return fIsTriggerMatch          ; }
  Bool_t           IsTriggerMatchedOpenCuts(Int_t i) const { return fIsTriggerMatchOpenCut[i]; }
  
  Int_t            GetTriggerClusterBC()             const { return fTriggerClusterBC        ; }
  Int_t            GetTriggerClusterIndex()          const { return fTriggerClusterIndex     ; }
  Int_t            GetTriggerClusterId()             const { return fTriggerClusterId        ; }
  
  Float_t          GetEventTriggerL0Threshold()      const { return fTriggerL0EventThreshold ; }
  void             SetEventTriggerL0Threshold(Float_t tr)  { fTriggerL0EventThreshold   = tr ; }
  Float_t          GetEventTriggerL1Threshold()      const { return fTriggerL1EventThreshold ; }
  void             SetEventTriggerL1Threshold(Float_t tr)  { fTriggerL1EventThreshold   = tr ; fTriggerL1EventThresholdFix = kTRUE; }

  void             SetEventTriggerL1Bit(Int_t ega, Int_t eje) { fBitEGA   = ega ; fBitEJE = eje; }
  
  void             SetTriggerPatchTimeWindow(Int_t min, Int_t max) { fTriggerPatchTimeWindow[0] = min ;
                                                                     fTriggerPatchTimeWindow[1] = max ; }
  
  Bool_t           AreBadTriggerEventsFromTriggerMakerRemoved()      const { return fRemoveBadTriggerEventsFromEMCalTriggerMaker     ; }
  void             SwitchOffBadTriggerEventsFromTriggerMakerRemoval()      { fRemoveBadTriggerEventsFromEMCalTriggerMaker   = kFALSE ; }
  void             SwitchOnBadTriggerEventsFromTriggerMakerRemoval()       { fRemoveBadTriggerEventsFromEMCalTriggerMaker   = kTRUE  ; }
  
  void             SetEMCalTriggerMakerDecisionContainerName(TString name) { fEMCalTriggerMakerDecissionContainerName = name ; }
  TString          GetEMCalTriggerMakerDecisionContainerName()       const { return fEMCalTriggerMakerDecissionContainerName ; }
   
  Bool_t           AreBadTriggerEventsRemoved()      const { return fRemoveBadTriggerEvents     ; }
  void             SwitchOffBadTriggerEventsRemoval()      { fRemoveBadTriggerEvents   = kFALSE ; }
  void             SwitchOnBadTriggerEventsRemoval()       { fRemoveBadTriggerEvents   = kTRUE  ; }

  Bool_t           AreUnMatchedTriggerEventsRemoved()const { return fRemoveUnMatchedTriggers    ; }
  void             SwitchOffUnMatchedTriggerEventsRemoval(){ fRemoveUnMatchedTriggers  = kFALSE ; }
  void             SwitchOnUnMatchedTriggerEventsRemoval() { fRemoveUnMatchedTriggers  = kTRUE  ; }
  
  Bool_t           IsTriggerPatchMatchedToCluster()  const { return fTriggerPatchClusterMatch   ; }
  void             SwitchOffTriggerPatchMatching()         { fTriggerPatchClusterMatch = kFALSE ; }
  void             SwitchOnTriggerPatchMatching()          { fTriggerPatchClusterMatch = kTRUE  ; }
  
  Bool_t           IsTriggerClusterTimeRecal()       const { return fTriggerClusterTimeRecal    ; }
  void             SwitchOnTriggerClusterTimeRecal ()      { fTriggerClusterTimeRecal  = kTRUE  ; }
  void             SwitchOffTriggerClusterTimeRecal()      { fTriggerClusterTimeRecal  = kFALSE ; }
  
  void             SetEventTriggerBit(UInt_t mask);
  Bool_t           IsEventMinimumBias()              const { return fEventTrigMinBias        ; }
  Bool_t           IsEventCentral()                  const { return fEventTrigCentral        ; }
  Bool_t           IsEventSemiCentral()              const { return fEventTrigSemiCentral    ; }
	
  Bool_t           IsEventEMCALL0()                  const { return fEventTrigEMCALL0        ; }
  Bool_t           IsEventEMCALL1Gamma1()            const { return fEventTrigEMCALL1Gamma1  ; }
  Bool_t           IsEventEMCALL1Gamma2()            const { return fEventTrigEMCALL1Gamma2  ; }
  Bool_t           IsEventEMCALL1Jet1()              const { return fEventTrigEMCALL1Jet1    ; }
  Bool_t           IsEventEMCALL1Jet2()              const { return fEventTrigEMCALL1Jet2    ; }
  Bool_t           IsEventEMCALL1Gamma()             const { return (fEventTrigEMCALL1Gamma1 || fEventTrigEMCALL1Gamma2) ; }
  Bool_t           IsEventEMCALL1Jet()               const { return (fEventTrigEMCALL1Jet1   || fEventTrigEMCALL1Jet2  ) ; }
  Bool_t           IsEventEMCALL1()                  const { return (IsEventEMCALL1Gamma()   || IsEventEMCALL1Jet()    ) ; }
  
  Bool_t           IsEventDCALL0()                   const { return fEventTrigDCALL0        ; }
  Bool_t           IsEventDCALL1Gamma1()             const { return fEventTrigDCALL1Gamma1  ; }
  Bool_t           IsEventDCALL1Gamma2()             const { return fEventTrigDCALL1Gamma2  ; }
  Bool_t           IsEventDCALL1Jet1()               const { return fEventTrigDCALL1Jet1    ; }
  Bool_t           IsEventDCALL1Jet2()               const { return fEventTrigDCALL1Jet2    ; }
  Bool_t           IsEventDCALL1Gamma()              const { return (fEventTrigDCALL1Gamma1 || fEventTrigDCALL1Gamma2) ; }
  Bool_t           IsEventDCALL1Jet()                const { return (fEventTrigDCALL1Jet1   || fEventTrigDCALL1Jet2  ) ; }
  Bool_t           IsEventDCALL1()                   const { return (IsEventDCALL1Gamma()   || IsEventDCALL1Jet()    ) ; }
  
  Bool_t           IsEventMinimumBiasCaloOnly()      const { return fEventTrigMinBiasCaloOnly        ; }
  Bool_t           IsEventEMCALL0CaloOnly()          const { return fEventTrigEMCALL0CaloOnly        ; }
  Bool_t           IsEventEMCALL1Gamma1CaloOnly()    const { return fEventTrigEMCALL1Gamma1CaloOnly  ; }
  Bool_t           IsEventEMCALL1Gamma2CaloOnly()    const { return fEventTrigEMCALL1Gamma2CaloOnly  ; }
  Bool_t           IsEventEMCALL1Jet1CaloOnly()      const { return fEventTrigEMCALL1Jet1CaloOnly    ; }
  Bool_t           IsEventEMCALL1Jet2CaloOnly()      const { return fEventTrigEMCALL1Jet2CaloOnly    ; }
  Bool_t           IsEventEMCALL1GammaCaloOnly()     const { return (fEventTrigEMCALL1Gamma1CaloOnly || fEventTrigEMCALL1Gamma2CaloOnly) ; }
  Bool_t           IsEventEMCALL1JetCaloOnly()       const { return (fEventTrigEMCALL1Jet1CaloOnly   || fEventTrigEMCALL1Jet2CaloOnly  ) ; }
  Bool_t           IsEventEMCALL1CaloOnly()          const { return (IsEventEMCALL1GammaCaloOnly()   || IsEventEMCALL1JetCaloOnly()    ) ; }
  
  Bool_t           IsEventDCALL0CaloOnly()           const { return fEventTrigDCALL0CaloOnly        ; }
  Bool_t           IsEventDCALL1Gamma1CaloOnly()     const { return fEventTrigDCALL1Gamma1CaloOnly  ; }
  Bool_t           IsEventDCALL1Gamma2CaloOnly()     const { return fEventTrigDCALL1Gamma2CaloOnly  ; }
  Bool_t           IsEventDCALL1Jet1CaloOnly()       const { return fEventTrigDCALL1Jet1CaloOnly    ; }
  Bool_t           IsEventDCALL1Jet2CaloOnly()       const { return fEventTrigDCALL1Jet2CaloOnly    ; }
  Bool_t           IsEventDCALL1GammaCaloOnly()      const { return (fEventTrigDCALL1Gamma1CaloOnly || fEventTrigDCALL1Gamma2CaloOnly) ; }
  Bool_t           IsEventDCALL1JetCaloOnly()        const { return (fEventTrigDCALL1Jet1CaloOnly   || fEventTrigDCALL1Jet2CaloOnly  ) ; }
  Bool_t           IsEventDCALL1CaloOnly()           const { return (IsEventDCALL1GammaCaloOnly()   || IsEventDCALL1JetCaloOnly()    ) ; }
  
  void             SwitchOnEMCALEventRejectionL1HighWithL1Low()  { fRejectEMCalTriggerEventsL1HighWithL1Low = kTRUE  ; }
  void             SwitchOffEMCALEventRejectionL1HighWithL1Low() { fRejectEMCalTriggerEventsL1HighWithL1Low = kFALSE ; }
	
  void             SwitchOnRemoveCentralityTriggerOutliers()     { fRemoveCentralityTriggerOutliers = kTRUE  ; }
  void             SwitchOffRemoveCentralityTriggerOutliers()    { fRemoveCentralityTriggerOutliers = kFALSE ; }
  
  // Other event rejections criteria
  
  void             SwitchOnPileUpEventRejection()          { fDoPileUpEventRejection= kTRUE  ; }
  void             SwitchOffPileUpEventRejection()         { fDoPileUpEventRejection= kFALSE ; }
  Bool_t           IsPileUpEventRejectionDone()      const { return fDoPileUpEventRejection  ; }
  
  void             SwitchOnV0ANDSelection()                { fDoV0ANDEventSelection = kTRUE  ; }
  void             SwitchOffV0ANDSelection()               { fDoV0ANDEventSelection = kFALSE ; }
  Bool_t           IsV0ANDEventSelectionDone()       const { return fDoV0ANDEventSelection   ; } 

  void             SwitchOnVertexBCEventSelection()        { fDoVertexBCEventSelection = kTRUE  ; }
  void             SwitchOffVertexBCEventSelection()       { fDoVertexBCEventSelection = kFALSE ; }
  Bool_t           IsVertexBCEventSelectionDone()    const { return fDoVertexBCEventSelection   ; }
  
  void             SwitchOnPrimaryVertexSelection()        { fUseEventsWithPrimaryVertex = kTRUE  ; }
  void             SwitchOffPrimaryVertexSelection()       { fUseEventsWithPrimaryVertex = kFALSE ; }
  Bool_t           IsPrimaryVertexSelectionDone()    const { return fUseEventsWithPrimaryVertex   ; } 
  
  void             SwitchOnRejectNoTrackEvents()           { fDoRejectNoTrackEvents = kTRUE  ; }
  void             SwitchOffRejectNoTrackEvents()          { fDoRejectNoTrackEvents = kFALSE ; }
  Bool_t           IsEventWithNoTrackRejectionDone() const { return fDoRejectNoTrackEvents   ; }

  void             UseEventCutsClass(Bool_t use)           { fUseEventCutsClass = use  ; }
  AliEventCuts    &GetEventCutsClass()                     { return fEventCuts         ; }
  
  // Time Stamp
  
  Double_t         GetRunTimeStampMin()              const { return fTimeStampRunMin         ; }
  Double_t         GetRunTimeStampMax()              const { return fTimeStampRunMax         ; }
  
  void             SetRunTimeStamp(Double_t a, Double_t b) { fTimeStampRunMin = a            ;
                                                             fTimeStampRunMax = b            ; } // seconds
  
  Float_t          GetEventTimeStampFractionMin()    const { return fTimeStampEventFracMin   ; }
  Float_t          GetEventTimeStampFractionMax()    const { return fTimeStampEventFracMax   ; }
  
  void             SetEventTimeStampFraction(Float_t a, Float_t b) { fTimeStampEventFracMin = a ;
                                                                     fTimeStampEventFracMax = b ; }
  
  void             SwitchOnSelectEventTimeStamp()          { fTimeStampEventSelect = kTRUE   ; }
  void             SwitchOffSelectEventTimeStamp()         { fTimeStampEventSelect = kFALSE  ; }
  
  Bool_t           IsSelectEventTimeStampOn()              { return  fTimeStampEventSelect   ; }

  // Time Stamp CTP corrected
    
  Double_t         GetTimeStampEventCTPBCCorrMin()   const { return fTimeStampEventCTPBCCorrMin ; }
  Double_t         GetTimeStampEventCTPBCCorrMax()   const { return fTimeStampEventCTPBCCorrMax ; }
  
  void             SetTimeStampEventCTPBCCorrRange(Double_t a, Double_t b) { 
                                                             fTimeStampEventCTPBCCorrMin = a    ;
                                                             fTimeStampEventCTPBCCorrMax = b    ; } // seconds
  
  void             SwitchOnExcludeEventTimeCTPBCCorrStamp() { fTimeStampEventCTPBCCorrExclude = kTRUE   ; }
  void             SwitchOffExcludeEventTimeCTPBCCorrStamp(){ fTimeStampEventCTPBCCorrExclude = kFALSE  ; }
  
  Bool_t           IsExcludeEventTimeStampCTPBCCorrOn()     { return  fTimeStampEventCTPBCCorrExclude ; }

  
  // Event tagging as pile-up
  
  Bool_t           IsPileUpFromSPD()               const ;
  Bool_t           IsPileUpFromEMCal()             const ;
  Bool_t           IsPileUpFromSPDAndEMCal()       const ;
  Bool_t           IsPileUpFromSPDOrEMCal()        const ;
  Bool_t           IsPileUpFromSPDAndNotEMCal()    const ;
  Bool_t           IsPileUpFromEMCalAndNotSPD()    const ;
  Bool_t           IsPileUpFromNotSPDAndNotEMCal() const ;

  void             SetPileUpParamForSPD  (Int_t i, Double_t param)
                                                           { fPileUpParamSPD[i]  = param  ; }
  void             SetPileUpParamForEMCal(Int_t param)     { fNPileUpClustersCut = param  ; }
  
  Int_t            GetNPileUpClusters()                    { return  fNPileUpClusters     ; }
  Int_t            GetNNonPileUpClusters()                 { return  fNNonPileUpClusters  ; }
  
  Int_t            GetEMCalEventBC(Int_t bc)     const     { if(bc >=0 && bc < 19) return  fEMCalBCEvent   [bc] ; else return 0 ; }
  Int_t            GetTrackEventBC(Int_t bc)     const     { if(bc >=0 && bc < 19) return  fTrackBCEvent   [bc] ; else return 0 ; }
  Int_t            GetEMCalEventBCcut(Int_t bc)  const     { if(bc >=0 && bc < 19) return  fEMCalBCEventCut[bc] ; else return 0 ; }
  Int_t            GetTrackEventBCcut(Int_t bc)  const     { if(bc >=0 && bc < 19) return  fTrackBCEventCut[bc] ; else return 0 ; }

  void             SetEMCalEventBC(Int_t bc)               { if(bc >=0 && bc < 19) fEMCalBCEvent   [bc] = 1 ; }
  void             SetTrackEventBC(Int_t bc)               { if(bc >=0 && bc < 19) fTrackBCEvent   [bc] = 1 ; }
  void             SetEMCalEventBCcut(Int_t bc)            { if(bc >=0 && bc < 19) fEMCalBCEventCut[bc] = 1 ; }
  void             SetTrackEventBCcut(Int_t bc)            { if(bc >=0 && bc < 19) fTrackBCEventCut[bc] = 1 ; }

  Int_t            GetVertexBC(const AliVVertex * vtx);
  Int_t            GetVertexBC()                  const    { return fVertexBC              ; }
  void             SwitchOnRecalculateVertexBC()           { fRecalculateVertexBC = kTRUE  ; fAccessTrackTOF  = kTRUE ; }
  void             SwitchOffRecalculateVertexBC()          { fRecalculateVertexBC = kFALSE ; }
  
  // Track selection
  
  ULong_t          GetTrackStatus()                  const { return fTrackStatus          ; }
  void             SetTrackStatus(ULong_t bit)             { fTrackStatus = bit           ; }

  virtual Bool_t   SelectTrack(AliVTrack* , Double_t*)     { return kFALSE                ; } // See AOD/ESD reader
  
  void             SwitchOnTrackHitSPDSelection()          { fSelectSPDHitTracks = kTRUE  ; }
  void             SwitchOffTrackHitSPDSelection()         { fSelectSPDHitTracks = kFALSE ; }
  
  Int_t            GetMinimumITSclusters()           const { return fSelectMinITSclusters ; }
  void             SetMinimumITSclusters(Int_t min)        { fSelectMinITSclusters = min  ; }
  
  Float_t          GetMaximumChi2PerITScluster()      const { return fSelectMaxChi2PerITScluster ; }
  void             SetMaximumChi2PerITScluster(Float_t max) { fSelectMaxChi2PerITScluster = max  ; }

  Int_t            GetMinimumTPCclusters()           const { return fSelectMinTPCclusters ; }
  void             SetMinimumTPCclusters(Int_t min)        { fSelectMinTPCclusters = min  ; }

  Float_t          GetMaximumChi2PerTPCcluster()      const { return fSelectMaxChi2PerTPCcluster ; }
  void             SetMaximumChi2PerTPCcluster(Float_t max) { fSelectMaxChi2PerTPCcluster = max  ; }

  Int_t            GetTrackMultiplicity(Int_t cut=0) const 
  {  if(cut < 10)  return fTrackMult [cut] ; else return 0 ; }
  Float_t          GetTrackSumPt(Int_t cut=0) const 
  {  if(cut < 10)  return fTrackSumPt[cut] ; else return 0 ; }

  void             SetTrackMultiplicityNPtCut(Float_t ncut){ fTrackMultNPtCut = ncut      ; }
  Int_t            GetTrackMultiplicityNPtCut() const      { return fTrackMultNPtCut      ; }

  void             SetTrackMultiplicityPtCut(Int_t cut, Float_t  pt) {  if(cut < 10)  fTrackMultPtCut[cut] = pt; }
  Float_t          GetTrackMultiplicityPtCut(Int_t cut=0) const 
  {  if(cut < 10)  return fTrackMultPtCut[cut] ;  else return 0 ; }
 
  Float_t          GetTrackMultiplicityEtaCut()      const { return fTrackMultEtaCut      ; }
  void             SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut = eta       ; }		
  
  // Virtual for AliCaloTrackAODReader
  
  virtual ULong_t  GetTrackFilterMask()               const { return 0 ; }
  virtual void     SetTrackFilterMask(ULong_t)              { ; }
  
  virtual ULong_t  GetTrackFilterMaskComplementary()  const { return 0 ; }
  virtual void     SetTrackFilterMaskComplementary(ULong_t) { ; }
  
  virtual void     SwitchOnAODHybridTrackSelection()        { ; }
  virtual void     SwitchOffAODHybridTrackSelection()       { ; }
  
  virtual void     SwitchOnAODPrimaryTrackSelection()       { ; }
  virtual void     SwitchOffAODPrimaryTrackSelection()      { ; }
  
  virtual void     SwitchOnAODTrackSharedClusterSelection() { ; }
  virtual void     SwitchOffAODTrackSharedClusterSelection(){ ; }
  
  virtual void     SetTPCSharedClusterFraction(Float_t)     { ; }
  virtual Float_t  GetTPCSharedClusterFraction() const      { return 0 ; }

  // Virtual for AliCaloTrackESDReader
  
  virtual AliESDtrackCuts* GetTrackCuts()              const { return 0 ; }
  virtual AliESDtrackCuts* GetTrackComplementaryCuts() const { return 0 ; }
  
  virtual void     SetTrackCuts(AliESDtrackCuts *)               { ; }
  virtual void     SetTrackComplementaryCuts(AliESDtrackCuts *)  { ; }
  
  virtual void     SwitchOnConstrainTrackToVertex()        { ; }
  virtual void     SwitchOffConstrainTrackToVertex()       { ; }

  // Events species selection
  
  void             AnalyzeOnlyLEDEvents()                  { fEventType  = 8      ; }
  void             AnalyzeOnlyPhysicsEvents()              { fEventType  = 7      ; }
  void             AnalyzeOnlyEventsOfType(Int_t specie)   { fEventType  = specie ; }
  
  //-------------------------------
  // Vertex methods
  //-------------------------------

  virtual void      GetVertex(Double_t v[3])         const ;
  virtual Double_t* GetVertex(Int_t evtIndex)        const { return fVertex[evtIndex]            ; }
  virtual void      GetVertex(Double_t vertex[3],    const Int_t evtIndex) const ;
  virtual void      FillVertexArray();
  virtual Bool_t    CheckForPrimaryVertex()          const { return kTRUE                        ; } // algorithm in ESD/AOD Readers
  virtual Float_t   GetZvertexCut()                  const { return fZvtxCut                     ; } // cut on vertex position
  virtual void      SetZvertexCut(Float_t zcut=10.)        { fZvtxCut=zcut                       ; } // cut on vertex position

  //--------------------------
  // Centrality / Event Plane
  //--------------------------

  virtual AliCentrality*    GetCentrality()          const ;
  
  virtual AliMultSelection* GetMultSelCen()          const ;

  void             SetMultiplicityWithPhysSel( Bool_t ps ) { fMultWithEventSel  = ps             ; }
  
  virtual void     SwitchOnAliCentrality ()                { fUseAliCentrality  = kTRUE          ; }
  virtual void     SwitchOffAliCentrality()                { fUseAliCentrality  = kFALSE         ; }
  
  virtual void     SetCentralityClass(TString name)        { fCentralityClass   = name           ; }
  virtual void     SetCentralityOpt(Int_t opt)             { fCentralityOpt     = opt            ; }
  virtual TString  GetCentralityClass()              const { return fCentralityClass             ; }
  virtual Int_t    GetCentralityOpt()                const { return fCentralityOpt               ; }
  virtual Float_t  GetEventCentralityF()             const ;
  virtual Int_t    GetEventCentrality()              const { return TMath::Floor(GetEventCentralityF()); }
  virtual void     SetCentralityBin(Int_t min, Int_t max) //Set the centrality bin to select the event. If used, then need to get percentile
                                                           { fCentralityBin[0]=min; fCentralityBin[1]=max;  
                                                             if(min>=0 && max > 0) fCentralityOpt = 100 ; }
  virtual Float_t  GetCentralityBin(Int_t i)         const { if(i < 0 || i > 1) return 0 ; 
                                                             else return fCentralityBin[i]              ; }
  
  virtual AliEventplane* GetEventPlane()             const { if(fDataType!=kMC) return fInputEvent->GetEventplane() ;
                                                             else               return 0x0       ; }  
  virtual Double_t       GetEventPlaneAngle()        const ;          
  virtual void           SetEventPlaneMethod(TString m)    { fEventPlaneMethod = m               ; }
  virtual TString        GetEventPlaneMethod()       const { return fEventPlaneMethod            ; }

  //--------------------
  // Mixing
  //--------------------

  Int_t   GetLastCaloMixedEvent()                    const { return fLastMixedCaloEvent          ; }
  Int_t   GetLastTracksMixedEvent ()                 const { return fLastMixedTracksEvent        ; }
  
  TList * GetListWithMixedEventsForCalo  (Int_t bi)  const { if(fListMixedCaloEvents)   return fListMixedCaloEvents  [bi] ; else return 0 ; }
  TList * GetListWithMixedEventsForTracks(Int_t bi)  const { if(fListMixedTracksEvents) return fListMixedTracksEvents[bi] ; else return 0 ; }
   
  Bool_t  ListWithMixedEventsForCaloExists()         const { if(fListMixedCaloEvents) return kTRUE  ;
                                                             else                     return kFALSE ; }

  Bool_t  ListWithMixedEventsForTracksExists()       const { if(fListMixedTracksEvents) return kTRUE  ;
                                                             else                       return kFALSE ; }
  
  void    SetLastCaloMixedEvent  (Int_t e)                 { fLastMixedCaloEvent    = e          ; }
  void    SetLastTracksMixedEvent(Int_t e)                 { fLastMixedTracksEvent  = e          ; }
  
  void    SetListWithMixedEventsForCalo (TList ** l)       { 
            if(fListMixedCaloEvents) printf("AliCaloTrackReader::SetListWithMixedEventsForCalo() - Track Mixing event list already set, nothing done\n");
            else                        fListMixedCaloEvents    = l ; }
  
  void    SetListWithMixedEventsForTracks(TList ** l)      { 
            if(fListMixedTracksEvents)  printf("AliCaloTrackReader::SetListWithMixedEventsForTracks() - Calorimeter Mixing event list already set, nothing done\n");
            else                        fListMixedTracksEvents  = l ; }
  
  //-------------------------------------
  // Other methods
  //-------------------------------------
  
  AliCalorimeterUtils * GetCaloUtils()               const { return fCaloUtils                       ; }
  void                  SetCaloUtils(AliCalorimeterUtils * caloutils)
                                                           { fCaloUtils = caloutils                  ; }
    
  Double_t              GetEventWeight()             const { return fEventWeight                     ; }
  AliAnaWeights       * GetWeightUtils()                   { if ( !fWeightUtils ) fWeightUtils = new AliAnaWeights() ;
                                                             return               fWeightUtils       ; }

  virtual Double_t      GetBField()                  const { return fInputEvent->GetMagneticField()  ; }
  
  /// Shift phi angle in case of negative value 360 degrees. Example TLorenzVector::Phi defined in -pi to pi
  Float_t               GetPhi  (Float_t phi)        const { if ( phi < 0 ) phi += TMath::TwoPi() ; return phi ; }
  
  Float_t               DegToRad(Float_t deg)        const { deg *= TMath::DegToRad(); return deg ; }
  
  Float_t               RadToDeg(Float_t rad)        const { rad *= TMath::RadToDeg(); return rad ; }

  virtual AliMCAnalysisUtils * GetMCAnalysisUtils()        { return           fMCUtils ; } 
  virtual void                 SetMCAnalysisUtils(AliMCAnalysisUtils * mcutils) { 
                                                             if (  fMCUtils ) delete fMCUtils; 
                                                             fMCUtils = mcutils ; }

  //------------------------------------------------
  // MC analysis specific methods
  //-------------------------------------------------
  
  // Kinematics and galice.root available
  
  virtual AliHeader*         GetHeader()            const ;
  virtual AliGenEventHeader* GetGenEventHeader()    const { return fGenEventHeader       ; }
  virtual AliGenPythiaEventHeader* GetGenPythiaEventHeader() 
                                                    const { return fGenPythiaEventHeader ; }
  
  virtual void SwitchOnPythiaEventHeaderUse()             { fCheckPythiaEventHeader = kTRUE ; }
  virtual void SwitchOffPythiaEventHeaderUse()            { fCheckPythiaEventHeader = kFALSE ; }
  virtual Bool_t IsPythiaEventHeaderUsed()          const { return fCheckPythiaEventHeader ; }
  
  // See implementation in AOD and ESD readers
  
  // Filtered kinematics in AOD
  
  virtual TClonesArray*     GetAODMCParticles()      const ;
  virtual AliAODMCHeader*   GetAODMCHeader   ()      const ;
  
  virtual AliVEvent*        GetInputEvent()          const { return fInputEvent            ; }
  virtual AliVEvent*        GetOriginalInputEvent()  const { return 0x0                    ; }
  virtual AliAODEvent*      GetOutputEvent()         const { return fOutputEvent           ; }
  virtual AliMCEvent*       GetMC()                  const { return fMC                    ; }
  virtual AliMixedEvent*    GetMixedEvent()          const { return fMixedEvent            ; }
  virtual Int_t             GetNMixedEvent()         const { return fNMixedEvent           ; } 
  
  void             SwitchOnStack()                         { AliError("Obsolete, remove this setting in AddTask") ; }
  void             SwitchOffStack()                        { AliError("Obsolete, remove this setting in AddTask") ; }
  void             SwitchOnAODMCParticles()                { AliError("Obsolete, remove this setting in AddTask") ; }
  void             SwitchOffAODMCParticles()               { AliError("Obsolete, remove this setting in AddTask") ; }

  void             RemapMCLabelForAODs(Int_t &label);
  
  // Select generated events, depending on comparison of pT hard and jets
    
  virtual Bool_t   ComparePtHardAndJetPt(Int_t process, TString processName) ;
  virtual Bool_t   IsPtHardAndJetPtComparisonSet()       const { return  fComparePtHardAndJetPt   ; }
  virtual void     SetPtHardAndJetPtComparison(Bool_t compare) { fComparePtHardAndJetPt = compare ; }	
  virtual Float_t  GetPtHardAndJetFactor()               const { return  fPtHardAndJetPtFactor    ; }
  virtual void     SetPtHardAndJetPtFactor(Float_t factor)     { fPtHardAndJetPtFactor = factor   ; }		
  
  virtual Bool_t   ComparePtHardAndClusterPt(Int_t process, TString processName) ;
  virtual Bool_t   IsPtHardAndClusterPtComparisonSet()       const { return  fComparePtHardAndClusterPt   ; }
  virtual void     SetPtHardAndClusterPtComparison(Bool_t compare) { fComparePtHardAndClusterPt = compare ; }	
  virtual Float_t  GetPtHardAndClusterFactor()               const { return  fPtHardAndClusterPtFactor    ; }
  virtual void     SetPtHardAndClusterPtFactor(Float_t factor)     { fPtHardAndClusterPtFactor = factor   ; }		

  virtual Bool_t   ComparePtHardAndPromptPhotonPt(Int_t process, TString processName) ;
  virtual Bool_t   IsPtHardAndPromptPhotonPtComparisonSet()       const { return  fComparePtHardAndPromptPhotonPt   ; }
  virtual void     SetPtHardAndPromptPhotonPtComparison(Bool_t compare) { fComparePtHardAndPromptPhotonPt = compare ; }
  virtual Float_t  GetPtHardAndPromptPhotonFactor()               const { return  fPtHardAndPromptPhotonPtFactor    ; }
  virtual void     SetPtHardAndPromptPhotonPtFactor(Float_t factor)     { fPtHardAndPromptPhotonPtFactor = factor   ; }

  // Select particles or clusters depending on generator
  virtual void     SetNumberOfMCGeneratorsToAccept(Int_t nGen) 
  { fNMCGenerToAccept = nGen ; 
    if      ( nGen > 5 ) fNMCGenerToAccept = 5 ; 
    else if ( nGen < 0 ) fNMCGenerToAccept = 0 ; }
  virtual Int_t    GetNumberOfMCGeneratorsToAccept()         const { return fNMCGenerToAccept ; } 
  
  virtual void     SetNameOfMCGeneratorsToAccept(Int_t ig, TString name) 
  { if ( ig < 5 && ig >= 0 ) fMCGenerToAccept[ig] = name ; }  
  virtual void     SetIndexOfMCGeneratorsToAccept(Int_t ig, Int_t index) 
  { if ( ig < 5 && ig >= 0 ) fMCGenerIndexToAccept[ig] = index ; }  
  virtual TString GetNameOfMCGeneratorsToAccept(Int_t ig)   const { return fMCGenerToAccept[ig] ; }
  virtual Int_t   GetIndexOfMCGeneratorsToAccept(Int_t ig)  const { return fMCGenerIndexToAccept[ig] ; }
  
  Bool_t           AcceptParticleMCLabel(Int_t mcLabel)     const ;
  Int_t            GetCocktailGeneratorAndIndex(Int_t index, TString & nameGen) const ;
  TString          GetGeneratorNameAndIndex(Int_t index, Int_t & genIndex) const ;

  virtual void     SetNameOfMCEventHederGeneratorToAccept(TString name) { fMCGenerEventHeaderToAccept = name ; }
  virtual TString  GetNameOfMCEventHederGeneratorToAccept()       const { return fMCGenerEventHeaderToAccept ; }
  
  void   SwitchOnMCPromptPhotonsSelection()         { fAcceptMCPromptPhotonOnly    = kTRUE ; }
  void   SwitchOffMCPromptPhotonsSelection()        { fAcceptMCPromptPhotonOnly    = kFALSE; }
  Bool_t AreMCPromptPhotonsSelected()         const { return fAcceptMCPromptPhotonOnly     ; }
  void   SwitchOnMCFragmentationPhotonsRejection()  { fRejectMCFragmentationPhoton = kTRUE ; }
  void   SwitchOffMCFragmentationPhotonsRejection() { fRejectMCFragmentationPhoton = kFALSE; }
  Bool_t AreMCFragmentationPhotonsRejected()  const { return fRejectMCFragmentationPhoton  ; }

  // MC reader methods, declared there to allow compilation, they are only used in the MC reader
  
  virtual void AddNeutralParticlesArray(TArrayI & /*array*/) { ; }  
  virtual void AddChargedParticlesArray(TArrayI & /*array*/) { ; } 
  virtual void AddStatusArray(TArrayI & /*array*/)           { ; }
  
  virtual void SwitchOnPi0Decay()                            { ; } 
  virtual void SwitchOffPi0Decay()                           { ; } 
  virtual void SwitchOnStatusSelection()                     { ; }
  virtual void SwitchOffStatusSelection()                    { ; }
  virtual void SwitchOnOverlapCheck()                        { ; }
  virtual void SwitchOffOverlapCheck()                       { ; }
  virtual void SwitchOnOnlyGeneratorParticles()              { ; }
  virtual void SwitchOffOnlyGeneratorParticles()             { ; }

  virtual void SetEMCALOverlapAngle(Float_t /*angle*/)       { ; }
  virtual void SetPHOSOverlapAngle(Float_t /*angle*/)        { ; }

  //-------------
  // Jets
  //-------------
  
  Bool_t       IsNonStandardJetsSwitchedOn()           const { return fFillInputNonStandardJetBranch   ; }
  void         SwitchOnNonStandardJets()                     { fFillInputNonStandardJetBranch = kTRUE  ; }
  void         SwitchOffNonStandardJets()                    { fFillInputNonStandardJetBranch = kFALSE ; }
  
  Bool_t       IsBackgroundJetsSwitchedOn()           const { return fFillInputBackgroundJetBranch   ; }
  void         SwitchOnBackgroundJets()                     { fFillInputBackgroundJetBranch = kTRUE  ; }
  void         SwitchOffBackgroundJets()                    { fFillInputBackgroundJetBranch = kFALSE ; }

  virtual void FillInputNonStandardJets() ;
  virtual TClonesArray* GetNonStandardJets()            const { return fNonStandardJets                 ; }
  virtual void SetInputNonStandardJetBranchName(TString name) { fInputNonStandardJetBranchName   = name ; }
  virtual TString GetInputNonStandardJetBranchName()          { return fInputNonStandardJetBranchName   ; }
  
  virtual void FillInputBackgroundJets() ;
  //  virtual AliAODJetEventBackground* GetBackgroundJets() const { return fBackgroundJets                 ; }
  virtual TClonesArray* GetBackgroundJets() const { return fBackgroundJets                 ; }
  virtual void SetInputBackgroundJetBranchName(TString name) { fInputBackgroundJetBranchName   = name ; }
  virtual TString GetInputBackgroundJetBranchName()          { return fInputBackgroundJetBranchName   ; }

  //------------------
  // PAR runs
  //------------------
  Bool_t IsParRun()                                     const { return fParRun   ; }     
  void SwitchOnParRun()                                       { fParRun = kTRUE  ; }
  void SwitchOffParRun()                                      { fParRun = kFALSE ; }
  
 protected:
  
  Int_t	           fEventNumber;                   ///<  Event number.
  Int_t            fDataType ;                     ///<  Select MC: Kinematics, Data: ESD/AOD, MCData: Both.
  Int_t            fDebug;                         ///<  Debugging level.
  AliFiducialCut * fFiducialCut;                   ///<  Acceptance cuts.
  Bool_t           fCheckFidCut ;                  ///<  Do analysis for clusters in defined region.         

  Bool_t           fComparePtHardAndJetPt;         ///<  In MonteCarlo, jet events, reject fake events with wrong jet energy.
  Float_t          fPtHardAndJetPtFactor;          ///<  Factor between ptHard and jet pT to reject/accept event.

  Bool_t           fComparePtHardAndClusterPt;     ///<  In MonteCarlo, jet events, reject events with too large cluster energy.
  Float_t          fPtHardAndClusterPtFactor;      ///<  Factor between ptHard and cluster pT to reject/accept event.
  
  Bool_t           fComparePtHardAndPromptPhotonPt;///<  In MonteCarlo, prompt photon events, reject fake events with wrong jet energy.
  Float_t          fPtHardAndPromptPhotonPtFactor; ///<  Factor between ptHard and prompt photon pT to reject/accept event.

  Float_t          fCTSPtMin;                      ///<  pT Threshold on charged particles. 
  Float_t          fEMCALPtMin;                    ///<  pT Threshold on emcal clusters.
  Float_t          fPHOSPtMin;                     ///<  pT Threshold on phos clusters.
  Float_t          fCTSPtMax;                      ///<  pT Threshold on charged particles.
  Float_t          fEMCALPtMax;                    ///<  pT Threshold on emcal clusters.
  Float_t          fPHOSPtMax;                     ///<  pT Threshold on phos clusters.
  
  Float_t          fEMCALBadChMinDist ;            ///<  Minimal distance to bad channel to accept cluster in EMCal, cell units
  Float_t          fPHOSBadChMinDist ;             ///<  Minimal distance to bad channel to accept cluster in PHOS, cm

  Int_t            fEMCALNCellsCut ;               ///<  Accept for the analysis EMCAL clusters with more than fNCellsCut cells
  Float_t          fEMCALNCellsCutEnDepEnMin ;     ///<  Minimum cluster energy to apply energy dependent N cell cut in EMCal
  Float_t          fEMCALNCellsCutEnDepConstant ;  ///<  Constant value of energy dependent N cell cut in EMCal
  Float_t          fEMCALNCellsCutEnDepSlope ;     ///<  Slope value of energy depedent N cell cut in EMCal
  Int_t            fPHOSNCellsCut ;                ///<  Accept for the analysis PHOS clusters with more than fNCellsCut cells
  
  Float_t          fEMCALHighEnergyNdiffCut;       ///<  Minimum energy for which the cut on n diff T-Card = 0 is applied
  Float_t          fEMCALMinCellEnNdiffCut;        ///<  Minimum energy of cells used counting n diff in T-Card 
  
  Bool_t           fUseEMCALTimeCut;               ///<  Do time cut selection.
  Bool_t           fUseParamTimeCut;               ///<  Use simple or parametrized time cut.
  Bool_t           fUseTrackTimeCut;               ///<  Do time cut selection.
  Bool_t           fAccessTrackTOF;                ///<  Access the track TOF, in case of problems when accessing GetTOFBunchCrossing.
  Double_t         fEMCALTimeCutMin;               ///<  Remove clusters/cells with time smaller than this value, in ns.
  Double_t         fEMCALTimeCutMax;               ///<  Remove clusters/cells with time larger than this value, in ns.
  Float_t          fEMCALParamTimeCutMin[4];       ///<  Remove clusters/cells with time smaller than parametrized value, in ns.
  Double_t         fEMCALParamTimeCutMax[4];       ///<  Remove clusters/cells with time larger than parametrized value, in ns.
  Double_t         fTrackTimeCutMin;               ///<  Remove tracks with time smaller than this value, in ns.
  Double_t         fTrackTimeCutMax;               ///<  Remove tracks with time larger than this value, in ns.
  
  Bool_t           fUseTrackDCACut;                ///<  Do DCA selection.
  Double_t         fTrackDCACut[3];                ///<  Remove tracks with DCA larger than cut, parameters of function stored here.

  /// List with AOD branches created and needed in analysis.
  TList          * fAODBranchList ;                //-> 
  
  /// Temporal array with tracks.
  TObjArray      * fCTSTracks ;                    //-> 
  
  /// Temporal array with EMCAL CaloClusters.
  TObjArray      * fEMCALClusters ;                //-> 

  /// Temporal array with DCAL CaloClusters, not needed in the normal case, use just EMCal array with DCal limits.
  TObjArray      * fDCALClusters ;                 //-> 
  
  /// Temporal array with PHOS  CaloClusters.
  TObjArray      * fPHOSClusters ;                 //-> 
  
  AliVCaloCells  * fEMCALCells ;                   //!<! Temporal array with EMCAL AliVCaloCells.
  AliVCaloCells  * fPHOSCells ;                    //!<! Temporal array with PHOS  AliVCaloCells.

  AliVEvent      * fInputEvent;                    //!<! pointer to esd or aod input.
  AliAODEvent    * fOutputEvent;                   //!<! pointer to aod output.
  AliMCEvent     * fMC;                            //!<! Monte Carlo Event Handler.  

  Bool_t           fEmbeddedEvent[2];              ///< Data and MC events embedded with AliAnalysisTaskEmcalEmbeddingHelper
  Bool_t           fSelectEmbeddedClusters;        ///<  Use only simulated clusters that come from embedding.

  Bool_t           fFillCTS;                       ///<  Use data from CTS.
  Bool_t           fFillEMCAL;                     ///<  Use data from EMCAL.
  Bool_t           fFillDCAL;                      ///<  Use data from DCAL, not needed in the normal case, use just EMCal array with DCal limits.
  Bool_t           fFillPHOS;                      ///<  Use data from PHOS.
  Bool_t           fFillEMCALCells;                ///<  Use data from EMCAL.
  Bool_t           fFillPHOSCells;                 ///<  Use data from PHOS.
  Bool_t           fRecalculateClusters;           ///<  Correct clusters, recalculate them if recalibration parameters is given.
  Bool_t           fCorrectELinearity;             ///<  Correct cluster linearity, always on.
  
  Bool_t           fScaleEPerSM ;                  ///<  Scale cluster energy by a constant factor, depending on SM 
  Float_t          fScaleFactorPerSM[22];          ///<  Scale factor depending on SM number to be applied to cluster energy
  
  Bool_t           fSmearShowerShape;              ///<  Smear shower shape (use in MC).
  Float_t          fSmearShowerShapeWidth;         ///<  Smear shower shape landau function "width" (use in MC).
  TRandom3         fRandom ;                       //!<! Random generator.
  Int_t            fSmearingFunction;              ///<  Choice of smearing function. 0 no smearing. 1 smearing from Gustavo (Landau center at 0). 2 smearing from Astrid (Landau center at 0.05). See enum smearingFunction 
  Int_t            fSmearNLMMin ;                  ///< Do smearing for clusters with at least this value 
  Int_t            fSmearNLMMax ;                  ///< Do smearing for clusters with at maximum this value
  
  // Track selection and counting
  ULong_t          fTrackStatus        ;           ///<  Track selection bit, select tracks refitted in TPC, ITS ...
  Bool_t           fSelectSPDHitTracks ;           ///<  Ensure that track hits SPD layers, AOD
  Int_t            fSelectMinITSclusters;          ///< Ensure track as at least this number of ITS clusters, AOD
  Float_t          fSelectMaxChi2PerITScluster;    ///< Ensure track as less than this Chi2/nCls of ITS, AOD 
  Int_t            fSelectMinTPCclusters;          ///< Ensure track as at least this number of TPC clusters, AOD
  Float_t          fSelectMaxChi2PerTPCcluster;    ///< Ensure track as less than this Chi2/nCls of TPC, AOD
  
  Int_t            fTrackMult[10]      ;           ///<  Track multiplicity, count for different pT cuts
  Float_t          fTrackSumPt[10]     ;           ///<  Track sum pT, count for different pT cuts
  Int_t            fTrackMultNPtCut    ;           ///<  Track multiplicty, number of pt cuts
  Float_t          fTrackMultPtCut[10] ;           ///<  Track multiplicity and sum pt cuts list
  Float_t          fTrackMultEtaCut    ;           ///<  Track multiplicity eta cut.
  	
  TString          fDeltaAODFileName   ;           ///<  Delta AOD file name.
  TString          fFiredTriggerClassName;         ///<  Name of trigger event type used to do the analysis.

  // Trigger bit
  UInt_t           fEventTriggerMaskInput ;        ///<  Event input trigger mask
  UInt_t           fEventTriggerMask ;             ///<  Select this triggerered event.
  UInt_t           fMixEventTriggerMask ;          ///<  Select this triggerered event for mixing, tipically kMB or kINT7.
  Bool_t           fEventTriggerAtSE;              ///<  Select triggered event at SE base task or here.
  
  Bool_t           fEventTrigMinBias ;             ///<  Event is min bias on its name, it should correspond to AliVEvent::kMB, AliVEvent::kINT7.
  Bool_t           fEventTrigCentral ;             ///<  Event is AliVEvent::kCentral on its name,  it should correspond to PbPb.
  Bool_t           fEventTrigSemiCentral ;         ///<  Event is AliVEvent::kSemiCentral on its name,  it should correspond to PbPb.
  
  Bool_t           fEventTrigEMCALL0 ;             ///<  Event is EMCal L0 on its name, it should correspond to AliVEvent::kEMC7, AliVEvent::kEMC1.
  Bool_t           fEventTrigEMCALL1Gamma1 ;       ///<  Event is EMCal L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigEMCALL1Gamma2 ;       ///<  Event is EMCal L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigEMCALL1Jet1 ;         ///<  Event is EMCal L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigEMCALL1Jet2 ;         ///<  Event is EMCal L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA.
	Bool_t           fEventTrigDCALL0 ;              ///<  Event is DCal L0 on its name, it should correspond to AliVEvent::kEMC7, AliVEvent::kEMC1.
  Bool_t           fEventTrigDCALL1Gamma1 ;        ///<  Event is DCal L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigDCALL1Gamma2 ;        ///<  Event is DCal L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigDCALL1Jet1 ;          ///<  Event is DCal L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA.
  Bool_t           fEventTrigDCALL1Jet2 ;          ///<  Event is DCal L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA.
  
  Bool_t           fEventTrigMinBiasCaloOnly ;       ///<  Event is min bias on its name, it should correspond to AliVEvent::kCaloOnly.
  Bool_t           fEventTrigEMCALL0CaloOnly ;       ///<  Event is EMCal L0 on its name, it should correspond to AliVEvent::kCaloOnly.
  Bool_t           fEventTrigEMCALL1Gamma1CaloOnly ; ///<  Event is EMCal L1-Gamma, threshold 1 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigEMCALL1Gamma2CaloOnly ; ///<  Event is EMCal L1-Gamma, threshold 2 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigEMCALL1Jet1CaloOnly ;   ///<  Event is EMCal L1-Gamma, threshold 1 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigEMCALL1Jet2CaloOnly ;   ///<  Event is EMCal L1-Gamma, threshold 2 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigDCALL0CaloOnly ;        ///<  Event is DCal L0 on its name, it should correspond to AliVEvent::kCaloOnly.
  Bool_t           fEventTrigDCALL1Gamma1CaloOnly ;  ///<  Event is DCal L1-Gamma, threshold 1 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigDCALL1Gamma2CaloOnly ;  ///<  Event is DCal L1-Gamma, threshold 2 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigDCALL1Jet1CaloOnly ;    ///<  Event is DCal L1-Gamma, threshold 1 on its name,  it should correspond kCaloOnly.
  Bool_t           fEventTrigDCALL1Jet2CaloOnly ;    ///<  Event is DCal L1-Gamma, threshold 2 on its name,  it should correspond kCaloOnly.
  
  Int_t            fBitEGA;                        ///<  Trigger bit on VCaloTrigger for EGA.
  Int_t            fBitEJE;                        ///<  Trigger bit on VCaloTrigger for EJE.
	
  Int_t            fEventType;                     ///<  Set the event species: 7 physics, 0 MC, 8 LED (not useful now)

  TString          fTaskName;                      ///<  Name of task that executes the analysis.
	
  AliCalorimeterUtils * fCaloUtils ;               ///<  Pointer to AliCalorimeterUtils.
  AliMCAnalysisUtils  * fMCUtils;                  ///<  MonteCarlo Analysis utils. Initialized in SetMC()

  AliAnaWeights  * fWeightUtils ;                  ///<  Pointer to AliAnaWeights.
  Double_t         fEventWeight ;                  ///<  Weight assigned to the event when filling histograms.
    
  AliMixedEvent  * fMixedEvent  ;                  //!<! Mixed event object. This class is not the owner.
  Int_t            fNMixedEvent ;                  ///<  Number of events in mixed event buffer.
  Double_t      ** fVertex      ;                  //!<! Vertex array 3 dim for each mixed event buffer.
  
  AliEventCuts     fEventCuts;                     ///< Event selection utility
  Bool_t           fUseEventCutsClass;             ///< Use AliEventCuts class 
  
  TList **         fListMixedTracksEvents;         //!<! Container for tracks stored for different events, used in case of own mixing, set in analysis class.
  TList **         fListMixedCaloEvents  ;         //!<! Container for photon stored for different events, used in case of own mixing, set in analysis class.
  Int_t            fLastMixedTracksEvent ;         ///<  Temporary container with the last event added to the mixing list for tracks.
  Int_t            fLastMixedCaloEvent   ;         ///<  Temporary container with the last event added to the mixing list for photons.
   
  Bool_t           fWriteOutputDeltaAOD;           ///<  Write the created delta AOD objects into file.  
  
  Int_t            fV0ADC[2]    ;                  ///<  Integrated V0 signal.
  Int_t            fV0Mul[2]    ;                  ///<  Integrated V0 Multiplicity.

  TString          fEMCALClustersListName;         ///<  Alternative list of clusters produced elsewhere and not from InputEvent.
  TString          fEMCALCellsListName;            ///<  Alternative list of cells produced elsewhere and not from InputEvent.
  
  //  Event selection
  
  Float_t          fZvtxCut ;	                     ///<  Cut on vertex position.
  Bool_t           fAcceptFastCluster;             ///<  Accept events from fast cluster, exclude these events for LHC11a.
 
  // LED events
  Int_t            fRemoveLEDEvents;               ///<  Remove events where LED was wrongly firing - only EMCAL LHC11a for this equal to 1, generalized to any SM for larger
  Float_t          fLEDHighEnergyCutSM;            ///<  SM is too active if energy above this value, likely LED event 
  Int_t            fLEDHighNCellsCutSM;            ///<  SM is too active if n cells above this value, likely LED event 
  Float_t          fLEDLowEnergyCutSM3;            ///<  SM3 low activity if energy below this value, check activity on other SM for LED event (Run2)
  Int_t            fLEDLowNCellsCutSM3;            ///<  SM3 low activity if n cells below this value, check activity on other SM LED event (Run2)
  Float_t          fLEDMinCellEnergy;              ///<  Count or sum cells energy above this value to determine if event had LEDs
  Float_t          fLEDMaxCellEnergy;              ///<  Count or sum cells energy below this value to determine if event had LEDs

  Int_t            fRemoveLEDStripEvents;          ///<  Remove events where an LED strip or more was wrongly firing - only EMCAL 
  Int_t            fLEDEventMaxNumberOfStrips;     ///<  Cut on events with a number of too active strips
  Float_t          fLEDHighEnergyCutStrip[2];      ///<  SM strip is too active if energy above this value, likely LED event. [0] Full SM, [1] 1/3 SM 
  Int_t            fLEDHighNCellsCutStrip[2];      ///<  SM strip is too active if n cells above this value, likely LED event. [0] Full SM, [1] 1/3 SM  
  Float_t          fLEDLowEnergyCutSM3Strip;       ///<  SM3 strip low activity if energy below this value, check activity on other SM for LED event (Run2)
  Int_t            fLEDLowNCellsCutSM3Strip;       ///<  SM3 strip low activity if n cells below this value, check activity on other SM LED event (Run2)
 
  // Triggered event selection
  
  Bool_t           fRemoveBadTriggerEventsFromEMCalTriggerMaker;  ///<  Remove triggered events because recalculated trigger was bad
  TString          fEMCalTriggerMakerDecissionContainerName;      ///<  Name of container with trigger decission
  
  Bool_t           fRemoveBadTriggerEvents;        ///<  Remove triggered events because trigger was exotic, bad, or out of BC.
  Bool_t           fTriggerPatchClusterMatch;      ///<  Search for the trigger patch and check if associated cluster was the trigger.
  Int_t            fTriggerPatchTimeWindow[2];     ///<  Trigger patch selection window.
  
  Float_t          fTriggerL0EventThreshold;       ///<  L0 Threshold to look for triggered events, set outside.
  Float_t          fTriggerL1EventThreshold;       ///<  L1 Threshold to look for triggered events, set in data.
  Bool_t           fTriggerL1EventThresholdFix;    ///<  L1 Threshold is fix and set outside.
  
  Int_t            fTriggerClusterBC;              ///<  Event triggered by a cluster in BC -5 0 to 5.
  Int_t            fTriggerClusterIndex;           ///<  Index in clusters array of trigger cluster.
  Int_t            fTriggerClusterId;              ///<  Id of trigger cluster (cluster->GetID()).
  Bool_t           fIsExoticEvent;                 ///<  Exotic trigger event flag.
  Bool_t           fIsBadCellEvent;                ///<  Bad cell triggered event flag, any cell in cluster is bad.
  Bool_t           fIsBadMaxCellEvent;             ///<  Bad cell triggered event flag, only max energy cell is bad.
  Bool_t           fIsTriggerMatch;                ///<  Could match the event to a trigger patch?
  Bool_t           fIsTriggerMatchOpenCut[3];      ///<  Could not match the event to a trigger patch?, retry opening cuts.
  Bool_t           fTriggerClusterTimeRecal;       ///<  In case cluster already calibrated, do not try to recalibrate even if recalib on in AliEMCALRecoUtils.
  Bool_t           fRemoveUnMatchedTriggers;       ///<  Analyze events where trigger patch and cluster where found or not.
  
  
  Bool_t           fDoPileUpEventRejection;        ///<  Select pile-up events by SPD.
  Bool_t           fDoV0ANDEventSelection;         ///<  Select events depending on V0AND.
  Bool_t           fDoVertexBCEventSelection;      ///<  Select events with vertex on BC=0 or -100.
  Bool_t           fDoRejectNoTrackEvents;         ///<  Reject events with no selected tracks in event.
  Bool_t           fUseEventsWithPrimaryVertex ;   ///<  Select events with primary vertex.
//AliTriggerAnalysis* fTriggerAnalysis;            ///<  Access to trigger selection algorithm for V0AND calculation.
  
  Bool_t           fTimeStampEventSelect;          ///<  Select events within a fraction of data taking time.
  Float_t          fTimeStampEventFracMin;         ///<  Minimum value of time stamp fraction event.
  Float_t          fTimeStampEventFracMax;         ///<  Maximum value of time stamp fraction event.
  Double_t         fTimeStampRunMin;               ///<  Minimum value of time stamp in run.
  Double_t         fTimeStampRunMax;               ///<  Maximum value of time stamp in run.
  
  Bool_t           fTimeStampEventCTPBCCorrExclude; ///<  Activate event selection within a range of data taking time CTP corrected. ESD only.
  Double_t         fTimeStampEventCTPBCCorrMin;    ///<  Minimum value of time stamp corrected by CTP in run.
  Double_t         fTimeStampEventCTPBCCorrMax;    ///<  Maximum value of time stamp corrected by CTP in run.
  
  ///< Parameters to pass to method IsPileupFromSPD:
  ///< Int_t minContributors, Double_t minZdist, Double_t nSigmaZdist,Double_t nSigmaDiamXY,Double_t nSigmaDiamZ
  Double_t         fPileUpParamSPD[5];
    
  // Pile-up in EMCal
  
  Int_t            fNPileUpClusters;               ///<  Number of clusters with time avobe 20 ns.
  Int_t            fNNonPileUpClusters;            ///<  Number of clusters with time below 20 ns.
  Int_t            fNPileUpClustersCut;            ///<  Cut to select event as pile-up.
  Int_t            fEMCalBCEvent[19];              ///<  Fill one entry per event if there is a cluster in a given BC.
  Int_t            fEMCalBCEventCut[19];           ///<  Fill one entry per event if there is a cluster in a given BC, depend on cluster E, acceptance cut.
  Int_t            fTrackBCEvent[19];              ///<  Fill one entry per event if there is a track in a given BC.
  Int_t            fTrackBCEventCut[19];           ///<  Fill one entry per event if there is a track in a given BC, depend on track pT, acceptance cut
  Int_t            fVertexBC;                      ///<  Vertex BC.
  Bool_t           fRecalculateVertexBC;           ///<  Recalculate vertex BC from tracks pointing to vertex.
  
  // Centrality/Event plane
  Bool_t           fUseAliCentrality;              ///<  Select as centrality estimator AliCentrality (Run1) or AliMultSelection (Run1 and Run2)
  Bool_t           fMultWithEventSel;              ///<  Embedded event selection in multiplicity task activated
  TString          fCentralityClass;               ///<  Name of selected centrality class.     
  Int_t            fCentralityOpt;                 ///<  Option for the returned value of the centrality, possible options 5, 10, 100.
  Int_t            fCentralityBin[2];              ///<  Minimum and maximum value of the centrality for the analysis.
  TString          fEventPlaneMethod;              ///<  Name of event plane method, by default "Q".

  // Jets
  Bool_t           fFillInputNonStandardJetBranch; ///<  Flag to use data from non standard jets.
  TClonesArray *   fNonStandardJets;               //!<! Temporal array with jets.
  TString          fInputNonStandardJetBranchName; ///<  Name of non standard jet branch.
  Bool_t           fFillInputBackgroundJetBranch;  ///<  Flag to use data from background jets.
  //  AliAODJetEventBackground * fBackgroundJets;      //!<! Background jets.
  TClonesArray * fBackgroundJets;      //!<! Background jets.
  TString          fInputBackgroundJetBranchName;  ///<  Name of background jet branch.

  TArrayI          fAcceptEventsWithBit;           ///<  Accept events if trigger bit is on.
  TArrayI          fRejectEventsWithBit;           ///<  Reject events if trigger bit is on.

  Bool_t           fRejectEMCalTriggerEventsL1HighWithL1Low; ///< Reject events L1 high threshold events also triggered by low threshold.
  Bool_t           fRemoveCentralityTriggerOutliers; ///< Reject events from centrality triggers out of expected ranges (PbPb 2011,2018)
  
  TLorentzVector   fMomentum;                      //!<! Temporal TLorentzVector container, avoid declaration of TLorentzVectors per event.

  // Handle runs affected by PAR
  Bool_t           fParRun;                        ///<  Flag set true when run affected by PAR
  Short_t          fCurrentParIndex;               //!<! temporal PAR number based on event global to get L1 phase correction in PAR runs
  
  // cut control histograms
  
  TList *          fOutputContainer;               //!<! Output container with cut control histograms.
  TH2F  *          fhEMCALClusterEtaPhi;           //!<! Control histogram on EMCAL clusters acceptance, before fiducial cuts
  TH2F  *          fhEMCALClusterEtaPhiFidCut;     //!<! Control histogram on EMCAL clusters acceptance, after fiducial cuts
  TH2F  *          fhEMCALClusterDisToBadE;        //!<! Control histogram on EMCAL clusters distance to bad channels
  TH2F  *          fhEMCALClusterTimeE;            //!<! Control histogram on EMCAL timing
  TH1F  *          fhEMCALClusterCutsE         [9];//!<! Control histogram on the different EMCal cluster selection cuts, E
  TH1F  *          fhEMCALClusterCutsESignal   [9];//!<! Control histogram on the different EMCal cluster selection cuts, E. Embedded signal clusters.
  TH2F  *          fhEMCALClusterCutsECen      [9];//!<! Control histogram on the different EMCal cluster selection cuts, E vs centrality
  TH2F  *          fhEMCALClusterCutsECenSignal[9];//!<! Control histogram on the different EMCal cluster selection cuts, E vs centrality. Embedded signal clusters.
  TH1F  *          fhPHOSClusterCutsE [7];         //!<! Control histogram on the different PHOS cluster selection cuts, E
  TH1F  *          fhCTSTrackCutsPt   [6];         //!<! Control histogram on the different CTS tracks selection cuts, pT
  TH1F  *          fhCTSTrackCutsPtSignal[6];      //!<! Control histogram on the different CTS tracks selection cuts, pT. Embedded signal
  TH2F  *          fhCTSTrackCutsPtCen[6];         //!<! Control histogram on the different CTS tracks selection cuts, pT vs centrality
  TH2F  *          fhCTSTrackCutsPtCenSignal[6];   //!<! Control histogram on the different CTS tracks selection cuts, pT vs centrality. Embedded signal
  TH1F  *          fhEMCALClusterBadTrigger;       //!<! Control histogram on clusters E on bad triggered events
  TH1F  *          fhCentralityBadTrigger;         //!<! Control histogram on event centrality for bad triggered events
  TH2F  *          fhEMCALClusterCentralityBadTrigger; //!<! Control histogram on clusters E vs centrality on bad triggered events

  TH2F  *          fhEMCALNSumEnCellsPerSM;        //!<! Control histogram of LED events rejection
  TH2F  *          fhEMCALNSumEnCellsPerSMAfter;   //!<! Control histogram of LED events rejection, after cut
  TH2F  *          fhEMCALNSumEnCellsPerSMAfterStripCut; //!<! Control histogram of LED events rejection, after LED strip rejection
  TH2F  *          fhEMCALNSumEnCellsPerStrip;     //!<! Control histogram of LED events on strips rejection, after LED SM rejection
  TH2F  *          fhEMCALNSumEnCellsPerStripAfter;//!<! Control histogram of LED events on strips rejection, after strip LED and SM rejection

  TH1F  *          fhPtHardPtJetPtRatio;           //!<! Control histogram fraction of generated jet pT over pT hard
  TH1F  *          fhPtHardPromptPhotonPtRatio;    //!<! Control histogram fraction of generated prompt photon Pt over pT hard
  TH1F  *          fhPtHardEnClusterRatio;         //!<! Control histogram fraction of generated calorimeter cluster energy over pT hard
  TH2F  *          fhPtHardEnClusterCenRatio;      //!<! Control histogram fraction of generated calorimeter cluster energy over pT hard  vs centrality

  Float_t          fEnergyHistogramLimit[2];       ///<  Binning of the control histograms, number of bins
  Int_t            fEnergyHistogramNbins ;         ///<  Binning of the control histograms, min and max window
  Bool_t           fHistoCentDependent;            ///< Fill centrality dependent of some histograms  
  Bool_t           fHistoPtDependent;              ///< Fill control histograms with Pt not E 
  
  TH1I  *          fhNEventsAfterCut;              //!<! Each bin represents number of events resulting after a given selection cut: vertex, trigger, ...  

  // MC labels to accept
  Int_t            fNMCGenerToAccept;              ///<  Number of MC generators that should not be included in analysis
  TString          fMCGenerToAccept[5];            ///<  List with name of generators that should not be included
  Int_t            fMCGenerIndexToAccept[5];       ///<  List with index of generators that should not be included

  TString          fMCGenerEventHeaderToAccept;    ///<  Accept events that contain at least this event header name

  AliGenEventHeader       * fGenEventHeader;       //!<! Event header
  AliGenPythiaEventHeader * fGenPythiaEventHeader; //!<! Event header casted to pythia
  Bool_t                    fCheckPythiaEventHeader; ///< Switch on/off recovery of the Pythia event header
  
  Bool_t           fAcceptMCPromptPhotonOnly ;     ///< Accept in the analysis task (AliiAnaPhoton) only cluster from prompt photons, to be used in gamma-jet simulations
  Bool_t           fRejectMCFragmentationPhoton ;  ///< Reject in the analysis task (AliAnaPhoton) clusters from fragmentation photons, to be used in jet-jet simulations

  /// Copy constructor not implemented.
  AliCaloTrackReader(              const AliCaloTrackReader & r) ; 
  
  /// Assignment operator not implemented.
  AliCaloTrackReader & operator = (const AliCaloTrackReader & r) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliCaloTrackReader,95) ;
  /// \endcond

} ;

#endif //ALICALOTRACKREADER_H



