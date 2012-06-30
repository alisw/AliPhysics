#ifndef ALICALOTRACKREADER_H
#define ALICALOTRACKREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors.
// Not all MC particles/tracks/clusters are kept, some kinematical restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader : Fills Kinematics data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS)   
// -- Author: Gustavo Conesa (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TObject.h> 
#include <TString.h>
class TObjArray ; 
class TTree ;
class TArrayI ;

//--- ANALYSIS system ---
#include "AliVEvent.h"
class AliVCaloCells;
class AliStack; 
class AliHeader; 
class AliGenEventHeader; 
class AliAODEvent;
class AliMCEvent;
class AliMixedEvent;
class AliAODMCHeader;
class AliESDtrackCuts;
class AliCentrality;
class AliTriggerAnalysis;
class AliEventplane;
class AliVCluster;

// --- CaloTrackCorr / EMCAL ---
#include "AliFiducialCut.h"
class AliEMCALRecoUtils;
class AliCalorimeterUtils;

class AliCaloTrackReader : public TObject {

public: 
  
  AliCaloTrackReader() ;          // ctor
  virtual ~AliCaloTrackReader() ; // virtual dtor
  
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
	
  TString         GetTaskName()                      const { return fTaskName              ; }
  void            SetTaskName(TString name)                { fTaskName = name              ; }
    
  //---------------------------------------
  //Input/output event setters and getters
  //---------------------------------------
  virtual void    SetInputEvent(AliVEvent* const input) ;
  virtual void    SetOutputEvent(AliAODEvent* const aod)   { fOutputEvent = aod            ; }
  virtual void    SetMC(AliMCEvent* const mc)              { fMC          = mc             ; }
  virtual void    SetInputOutputMCEvent(AliVEvent* /*esd*/, AliAODEvent* /*aod*/, AliMCEvent* /*mc*/) { ; }
  
  // Delta AODs
  virtual TList * GetAODBranchList()                 const { return fAODBranchList         ; }
  void            SetDeltaAODFileName(TString name )       { fDeltaAODFileName = name      ; }
  TString         GetDeltaAODFileName()              const { return fDeltaAODFileName      ; }
  void            SwitchOnWriteDeltaAOD()                  { fWriteOutputDeltaAOD = kTRUE  ; }
  void            SwitchOffWriteDeltaAOD()                 { fWriteOutputDeltaAOD = kFALSE ; }
  Bool_t          WriteDeltaAODToFile()              const { return fWriteOutputDeltaAOD   ; } 
  
  //------------------------------------------------------------
  //Clusters/Tracks arrays filtering/filling methods and switchs 
  //------------------------------------------------------------
  
  //Minimum pt setters and getters 
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
  
  void             SetEMCALEMin (Float_t  e)               { SetEMCALPtMin(e)              ; }
  void             SetPHOSEMin  (Float_t  e)               { SetPHOSPtMin (e)              ; }
  void             SetEMCALEMax (Float_t  e)               { SetEMCALPtMax(e)              ; }
  void             SetPHOSEMax  (Float_t  e)               { SetPHOSPtMax (e)              ; }
  
  Double_t         GetEMCALTimeCutMin()              const { return fEMCALTimeCutMin        ; }
  Double_t         GetEMCALTimeCutMax()              const { return fEMCALTimeCutMax        ; }	

  void             SetEMCALTimeCut(Double_t a, Double_t b) { fEMCALTimeCutMin = a; 
                                                             fEMCALTimeCutMax = b          ; } // ns
  // Fidutial cuts  
  virtual AliFiducialCut * GetFiducialCut()                { 
                    if(!fFiducialCut) fFiducialCut = new AliFiducialCut(); 
                    return  fFiducialCut                                                   ; }
  virtual void     SetFiducialCut(AliFiducialCut * const fc) { fFiducialCut = fc           ; }
  virtual Bool_t   IsFiducialCutOn()                 const { return fCheckFidCut           ; }
  virtual void     SwitchOnFiducialCut()                   { fCheckFidCut = kTRUE          ; 
                                                             fFiducialCut = new AliFiducialCut() ; }
  virtual void     SwitchOffFiducialCut()                  { fCheckFidCut = kFALSE         ; }
  
  // Cluster origin
  Bool_t           IsEMCALCluster(AliVCluster *clus) const;
  Bool_t           IsPHOSCluster (AliVCluster *clus) const;
  //Patch for cluster origin for Old AODs implementation
  void             SwitchOnOldAODs()                       { fOldAOD = kTRUE               ; }
  void             SwitchOffOldAODs()                      { fOldAOD = kFALSE              ; }
  
  // Cluster/track/cells switchs
  Bool_t           IsCTSSwitchedOn()                 const { return fFillCTS               ; }
  void             SwitchOnCTS()                           { fFillCTS = kTRUE              ; }
  void             SwitchOffCTS()                          { fFillCTS = kFALSE             ; }

  Bool_t           IsEMCALSwitchedOn()               const { return fFillEMCAL             ; }
  void             SwitchOnEMCAL()                         { fFillEMCAL = kTRUE            ; }
  void             SwitchOffEMCAL()                        { fFillEMCAL = kFALSE           ; }

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
  
  Bool_t           IsEmbeddedClusterSelectionOn()    const { return fSelectEmbeddedClusters   ; }
  void             SwitchOnEmbeddedClustersSelection()     { fSelectEmbeddedClusters = kTRUE  ; }
  void             SwitchOffEmbeddedClustersSelection()    { fSelectEmbeddedClusters = kFALSE ; }
  
  // Filling/ filtering / detector information access methods
  virtual Bool_t   FillInputEvent(const Int_t iEntry, const char *currentFileName)  ;
  virtual void     FillInputCTS() ;
  virtual void     FillInputEMCAL() ;
  virtual void     FillInputEMCALAlgorithm(AliVCluster * clus, const Int_t iclus) ;
  virtual void     FillInputPHOS() ;
  virtual void     FillInputEMCALCells() ;
  virtual void     FillInputPHOSCells() ;
  virtual void     FillInputVZERO() ;  
  
  Int_t            GetV0Signal(Int_t i)              const { return fV0ADC[i]               ; }
  Int_t            GetV0Multiplicity(Int_t i)        const { return fV0Mul[i]               ; }
  
  void             SetEMCALClusterListName(TString &name)  { fEMCALClustersListName = name  ; }
  TString          GetEMCALClusterListName()         const { return fEMCALClustersListName  ; }

  // Arrayes with clusters/track/cells access method
  virtual TObjArray*     GetCTSTracks()              const { return fCTSTracks              ; }
  virtual TObjArray*     GetEMCALClusters()          const { return fEMCALClusters          ; }
  virtual TObjArray*     GetPHOSClusters()           const { return fPHOSClusters           ; }
  virtual AliVCaloCells* GetEMCALCells()             const { return fEMCALCells             ; }
  virtual AliVCaloCells* GetPHOSCells()              const { return fPHOSCells              ; }
   
  //-------------------------------------
  // Event/track selection methods
  //-------------------------------------
  
  void             AcceptFastClusterEvents()               { fAcceptFastCluster     = kTRUE  ; } 
  void             RejectFastClusterEvents()               { fAcceptFastCluster     = kFALSE ; }  
  Bool_t           IsFastClusterAccepted()           const { return fAcceptFastCluster       ; }   
  
  void             SwitchOnLEDEventsRemoval()              { fRemoveLEDEvents       = kTRUE  ; }
  void             SwitchOffLEDEventsRemoval()             { fRemoveLEDEvents       = kFALSE ; } 
  Bool_t           IsLEDEventRemoved()               const { return fRemoveLEDEvents         ; }   

  void             SetFiredTriggerClassName(TString name ) { fFiredTriggerClassName = name   ; }
  TString          GetFiredTriggerClassName()        const { return fFiredTriggerClassName   ; }
  TString          GetFiredTriggerClasses() ;               
  
  UInt_t           GetEventTriggerMask()             const { return fEventTriggerMask        ; }
  void             SetEventTriggerMaks(UInt_t evtTrig = AliVEvent::kAny) 
                                                           { fEventTriggerMask = evtTrig     ; }

  UInt_t           GetMixEventTriggerMask()             const { return fMixEventTriggerMask  ; }
  void             SetMixEventTriggerMaks(UInt_t evtTrig = AliVEvent::kAnyINT) 
                                                           { fMixEventTriggerMask = evtTrig  ; }
  
  
  Bool_t           IsEventTriggerAtSEOn()            const { return fEventTriggerAtSE        ; }
  void             SwitchOnEventTriggerAtSE()              { fEventTriggerAtSE      = kTRUE  ; }
  void             SwitchOffEventTriggerAtSE()             { fEventTriggerAtSE      = kFALSE ; }
  
  void             SwitchOnEventSelection()                { fDoEventSelection      = kTRUE  ; }
  void             SwitchOffEventSelection()               { fDoEventSelection      = kFALSE ; }
  Bool_t           IsEventSelectionDone()            const { return fDoEventSelection        ; } 
  
  void             SwitchOnV0ANDSelection()                { fDoV0ANDEventSelection = kTRUE  ; }
  void             SwitchOffV0ANDSelection()               { fDoV0ANDEventSelection = kFALSE ; }
  Bool_t           IsV0ANDEventSelectionDone()       const { return fDoV0ANDEventSelection   ; } 

  void             SwitchOnPrimaryVertexSelection()        { fUseEventsWithPrimaryVertex = kTRUE  ; }
  void             SwitchOffPrimaryVertexSelection()       { fUseEventsWithPrimaryVertex = kFALSE ; }
  Bool_t           IsPrimaryVertexSelectionDone()    const { return fUseEventsWithPrimaryVertex   ; } 
  
  // Track selection
  ULong_t          GetTrackStatus()                  const { return fTrackStatus       ; }
  void             SetTrackStatus(ULong_t bit)             { fTrackStatus = bit        ; }		

  ULong_t          GetTrackFilterMask()              const {return fTrackFilterMask    ; }
  void             SetTrackFilterMask(ULong_t bit)         { fTrackFilterMask = bit    ; }		
  
  AliESDtrackCuts* GetTrackCuts()                    const { return fESDtrackCuts      ; }
  void             SetTrackCuts(AliESDtrackCuts * cuts)    ;

  void             SwitchOnConstrainTrackToVertex()        { fConstrainTrack = kTRUE   ; } 
  void             SwitchOffConstrainTrackToVertex()       { fConstrainTrack = kFALSE  ; }      
  
  void             SwitchOnAODHybridTrackSelection()       { fSelectHybridTracks = kTRUE  ; } 
  void             SwitchOffAODHybridTrackSelection()      { fSelectHybridTracks = kFALSE ; }      
  
  Int_t            GetTrackMultiplicity()            const { return fTrackMult         ; }
  Float_t          GetTrackMultiplicityEtaCut()      const { return fTrackMultEtaCut   ; }
  void             SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut = eta    ; }		
  
  // Calorimeter specific and patches
  void             AnalyzeOnlyLED()                        { fAnaLED = kTRUE           ; }
  void             AnalyzeOnlyPhysics()                    { fAnaLED = kFALSE          ; }
  
  void             SwitchOnCaloFilterPatch()               { fCaloFilterPatch = kTRUE  ; 
                                                             fFillCTS = kFALSE         ; }
  void             SwitchOffCaloFilterPatch()              { fCaloFilterPatch = kFALSE ; }
  Bool_t           IsCaloFilterPatchOn()             const { 
                    if(fDataType == kAOD) { return fCaloFilterPatch ; } 
                    else                  { return kFALSE           ; }                  }
  	
  //-------------------------------
  //Vertex methods
  //-------------------------------
  virtual void      GetVertex(Double_t v[3])         const ;
  virtual Double_t* GetVertex(const Int_t evtIndex)  const { return fVertex[evtIndex]            ; }
  virtual void      GetVertex(Double_t vertex[3],    const Int_t evtIndex) const ;
  virtual void      FillVertexArray();
  virtual Bool_t    CheckForPrimaryVertex();
  virtual Float_t   GetZvertexCut()                  const { return fZvtxCut                     ; } //cut on vertex position  
  virtual void      SetZvertexCut(Float_t zcut=10.)        { fZvtxCut=zcut                       ; } //cut on vertex position

  //--------------------------
  // Centrality / Event Plane
  //--------------------------
  virtual AliCentrality* GetCentrality()             const { return fInputEvent->GetCentrality() ; } //Look in AOD reader, different there
  virtual void     SetCentralityClass(TString name)        { fCentralityClass   = name           ; }
  virtual void     SetCentralityOpt(Int_t opt)             { fCentralityOpt     = opt            ; }
  virtual TString  GetCentralityClass()              const { return fCentralityClass             ; }
  virtual Int_t    GetCentralityOpt()                const { return fCentralityOpt               ; }
  virtual Int_t    GetEventCentrality()              const ;
  virtual void     SetCentralityBin(Int_t min, Int_t max) //Set the centrality bin to select the event. If used, then need to get percentile
                                                           { fCentralityBin[0]=min; fCentralityBin[1]=max;  
                                                             if(min>=0 && max > 0) fCentralityOpt = 100 ; }
  virtual Float_t  GetCentralityBin(Int_t i)         const { if(i < 0 || i > 1) return 0 ; 
                                                             else return fCentralityBin[i]              ; }
  
  virtual AliEventplane* GetEventPlane()             const { return fInputEvent->GetEventplane() ; }  
  virtual Double_t       GetEventPlaneAngle()        const ;          
  virtual void           SetEventPlaneMethod(TString m)    { fEventPlaneMethod = m               ; }
  virtual TString        GetEventPlaneMethod()       const { return fEventPlaneMethod            ; }

  //--------------------
  // Mixing
  //--------------------

  Int_t   GetLastCaloMixedEvent()                    const { return fLastMixedCaloEvent          ; }
  Int_t   GetLastTracksMixedEvent ()                 const { return fLastMixedTracksEvent        ; }
  
  TList * GetListWithMixedEventsForCalo  (Int_t bi)  const { if(fListMixedCaloEvents)   return fListMixedCaloEvents[bi]     ; else return 0 ; }
  TList * GetListWithMixedEventsForTracks(Int_t bi)  const { if(fListMixedTracksEvents) return fListMixedTracksEvents [bi]  ; else return 0 ; }  
   
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
  AliCalorimeterUtils * GetCaloUtils()               const { return fCaloUtils                   ; }
  void             SetCaloUtils(AliCalorimeterUtils * caloutils)  { fCaloUtils = caloutils       ; }  
  
  virtual Double_t GetBField()                       const { return fInputEvent->GetMagneticField()  ; } 
  
  void    SetImportGeometryFromFile(Bool_t import, 
                                    TString path = "")     { 
                                                             fImportGeometryFromFile = import    ; 
                                                             fImportGeometryFilePath = path      ; }      
  
  //------------------------------------------------
  // MC analysis specific methods
  //-------------------------------------------------
  
  //Kinematics and galice.root available 
  virtual AliStack*          GetStack()              const ;
  virtual AliHeader*         GetHeader()             const ;
  virtual AliGenEventHeader* GetGenEventHeader()     const ;
  
  //Filtered kinematics in AOD	
  virtual TClonesArray*     GetAODMCParticles(Int_t input = 0) const ;
  virtual AliAODMCHeader*   GetAODMCHeader(Int_t input = 0)    const ;
	
  virtual AliVEvent*        GetInputEvent()          const { return fInputEvent            ; }
  virtual AliVEvent*        GetOriginalInputEvent()  const { return 0x0                    ; }
  virtual AliAODEvent*      GetOutputEvent()         const { return fOutputEvent           ; }
  virtual AliMCEvent*       GetMC()                  const { return fMC                    ; }
  virtual AliMixedEvent*    GetMixedEvent()          const { return fMixedEvent            ; }
  virtual Int_t             GetNMixedEvent()         const { return fNMixedEvent           ; } 
  
  void             SwitchOnStack()                         { fReadStack          = kTRUE   ; }
  void             SwitchOffStack()                        { fReadStack          = kFALSE  ; }
  void             SwitchOnAODMCParticles()                { fReadAODMCParticles = kTRUE   ; }
  void             SwitchOffAODMCParticles()               { fReadAODMCParticles = kFALSE  ; }
  Bool_t           ReadStack()                       const { return fReadStack             ; }
  Bool_t           ReadAODMCParticles()              const { return fReadAODMCParticles    ; }
	
  //Select generated events, depending on comparison of pT hard and jets.
  virtual Bool_t   ComparePtHardAndJetPt() ;
  virtual Bool_t   IsPtHardAndJetPtComparisonSet()       const { return  fComparePtHardAndJetPt   ; }
  virtual void     SetPtHardAndJetPtComparison(Bool_t compare) { fComparePtHardAndJetPt = compare ; }	
  virtual Float_t  GetPtHardAndJetFactor()               const { return  fPtHardAndJetPtFactor    ; }
  virtual void     SetPtHardAndJetPtFactor(Float_t factor)     { fPtHardAndJetPtFactor = factor   ; }		
  
  //MC reader methods, declared there to allow compilation, they are only used in the MC reader:
  
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

  
 protected:
  Int_t	           fEventNumber;            // Event number
  Int_t            fDataType ;              // Select MC:Kinematics, Data:ESD/AOD, MCData:Both
  Int_t            fDebug;                  // Debugging level
  AliFiducialCut * fFiducialCut;            //! Acceptance cuts
  Bool_t           fCheckFidCut ;           // Do analysis for clusters in defined region         

  Bool_t           fComparePtHardAndJetPt;  // In MonteCarlo, jet events, reject fake events with wrong jet energy.
  Float_t          fPtHardAndJetPtFactor;   // Factor between ptHard and jet pT to reject/accept event.

  Float_t          fCTSPtMin;               // pT Threshold on charged particles 
  Float_t          fEMCALPtMin;             // pT Threshold on emcal clusters
  Float_t          fPHOSPtMin;              // pT Threshold on phos clusters
  Float_t          fCTSPtMax;               // pT Threshold on charged particles 
  Float_t          fEMCALPtMax;             // pT Threshold on emcal clusters
  Float_t          fPHOSPtMax;              // pT Threshold on phos clusters
  Double_t         fEMCALTimeCutMin;        // Remove clusters/cells with time smaller than this value, in ns
  Double_t         fEMCALTimeCutMax;        // Remove clusters/cells with time larger than this value, in ns
   
  TList          * fAODBranchList ;         //-> List with AOD branches created and needed in analysis  
  TObjArray      * fCTSTracks ;             //-> temporal array with tracks
  TObjArray      * fEMCALClusters ;         //-> temporal array with EMCAL CaloClusters
  TObjArray      * fPHOSClusters ;          //-> temporal array with PHOS  CaloClusters
  AliVCaloCells  * fEMCALCells ;            //! temporal array with EMCAL CaloCells
  AliVCaloCells  * fPHOSCells ;             //! temporal array with PHOS  CaloCells

  AliVEvent      * fInputEvent;             //! pointer to esd or aod input
  AliAODEvent    * fOutputEvent;            //! pointer to aod output
  AliMCEvent     * fMC;                     //! Monte Carlo Event Handler  

  Bool_t           fFillCTS;                // use data from CTS
  Bool_t           fFillEMCAL;              // use data from EMCAL
  Bool_t           fFillPHOS;               // use data from PHOS
  Bool_t           fFillEMCALCells;         // use data from EMCAL
  Bool_t           fFillPHOSCells;          // use data from PHOS
  Bool_t           fRecalculateClusters;    // Correct clusters, recalculate them if recalibration parameters is given
  Bool_t           fSelectEmbeddedClusters; // Use only simulated clusters that come from embedding.
  
  ULong_t          fTrackStatus        ;    // Track selection bit, select tracks refitted in TPC, ITS ...
  ULong_t          fTrackFilterMask    ;    // Track selection bit, for AODs (any difference with track status?)
  AliESDtrackCuts *fESDtrackCuts       ;    // Track cut 
  Bool_t           fConstrainTrack     ;    // Constrain Track to vertex
  Bool_t           fSelectHybridTracks ;    // Select CTS tracks of type hybrid (only for AODs)
  Int_t            fTrackMult          ;    // Track multiplicity
  Float_t          fTrackMultEtaCut    ;    // Track multiplicity eta cut
  Bool_t           fReadStack          ;    // Access kine information from stack
  Bool_t	         fReadAODMCParticles ;    // Access kine information from filtered AOD MC particles
	
  TString          fDeltaAODFileName   ;    // Delta AOD file name
  TString          fFiredTriggerClassName;  // Name of trigger event type used to do the analysis

  UInt_t           fEventTriggerMask ;      // select this triggerered event
  UInt_t           fMixEventTriggerMask ;   // select this triggerered event for mixing, tipically kMB or kAnyINT
  Bool_t           fEventTriggerAtSE;       // select triggered event at SE base task or here
  
  Bool_t           fAnaLED;                 // Analyze LED data only.

  TString          fTaskName;               // Name of task that executes the analysis
	
  AliCalorimeterUtils * fCaloUtils ;        //  Pointer to CalorimeterUtils

  AliMixedEvent  * fMixedEvent  ;           //! mixed event object. This class is not the owner
  Int_t            fNMixedEvent ;           // number of events in mixed event buffer
  Double_t      ** fVertex      ;           //! vertex array 3 dim for each mixed event buffer
  
  TList **         fListMixedTracksEvents ; //! Container for tracks stored for different events, used in case of own mixing, set in analysis class
  TList **         fListMixedCaloEvents;    //! Container for photon stored for different events, used in case of own mixing, set in analysis class
  Int_t            fLastMixedTracksEvent  ; //  Temporary container with the last event added to the mixing list for tracks
  Int_t            fLastMixedCaloEvent ;    //  Temporary container with the last event added to the mixing list for photons
  
  Bool_t           fWriteOutputDeltaAOD;    // Write the created delta AOD objects into file  
	Bool_t           fOldAOD;                 // Old AODs, before revision 4.20
  
  Int_t            fV0ADC[2]    ;           // Integrated V0 signal
  Int_t            fV0Mul[2]    ;           // Integrated V0 Multiplicity

  Bool_t           fCaloFilterPatch;        // CaloFilter patch
  TString          fEMCALClustersListName;  // Alternative list of clusters produced elsewhere and not from InputEvent
  
  // Event selection
  Float_t          fZvtxCut ;	                   // Cut on vertex position
  Bool_t           fAcceptFastCluster;           // Accept events from fast cluster, exclude these events for LHC11a
  Bool_t           fRemoveLEDEvents;             // Remove events where LED was wrongly firing - EMCAL LHC11a
  Bool_t           fDoEventSelection;            // Select events depending on V0, pileup, vertex well reconstructed, at least 1 track ...
  Bool_t           fDoV0ANDEventSelection;       // Select events depending on V0, fDoEventSelection should be on
  Bool_t           fUseEventsWithPrimaryVertex ; // Select events with primary vertex
  AliTriggerAnalysis* fTriggerAnalysis;          // Access to trigger selection algorithm for V0AND calculation
  
  //Centrality/Event plane
  TString          fCentralityClass;        // Name of selected centrality class     
  Int_t            fCentralityOpt;          // Option for the returned value of the centrality, possible options 5, 10, 100
  Int_t            fCentralityBin[2];       // Minimum and maximum value of the centrality for the analysis
  TString          fEventPlaneMethod;       // Name of event plane method, by default "Q"
  
  Bool_t           fImportGeometryFromFile; // Import geometry settings in geometry.root file
  TString          fImportGeometryFilePath; // path fo geometry.root file

  AliCaloTrackReader(              const AliCaloTrackReader & r) ; // cpy ctor
  AliCaloTrackReader & operator = (const AliCaloTrackReader & r) ; // cpy assignment
  
  ClassDef(AliCaloTrackReader,42)
  
} ;


#endif //ALICALOTRACKREADER_H



