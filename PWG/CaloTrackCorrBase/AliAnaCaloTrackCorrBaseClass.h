#ifndef ALIANACALOTRACKCORRBASECLASS_H
#define ALIANACALOTRACKCORRBASECLASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Base class for CaloTrackCorr analysis algorithms
//-- Author: Gustavo Conesa (INFN-LNF, LPSC-Grenoble)
//-- Add the possibility for event selection analysis based on 
//   vertex and multiplicity bins (Yaxian Mao, 10/10/2010)
//
//_________________________________________________________________________

#include <cstdlib>

//ROOT
#include <TList.h> 
#include <TObject.h>
class TClonesArray ;
class TObjArray ;
class TObjString;

//Analysis
#include "AliCaloTrackReader.h"   
#include "AliCaloPID.h"
#include "AliFiducialCut.h"
#include "AliIsolationCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliCalorimeterUtils.h" 
#include "AliHistogramRanges.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMixedEvent.h" 
class AliVCaloCells;
class AliStack ; 
class AliHeader ; 
class AliGenEventHeader ; 
class AliEMCALGeometry;
class AliPHOSGeoUtils;
class AliCentrality;
class AliEventplane;

class AliAnaCaloTrackCorrBaseClass : public TObject {
	
public:   
  AliAnaCaloTrackCorrBaseClass() ;          // default ctor
  virtual ~AliAnaCaloTrackCorrBaseClass() ; // virtual dtor
  
  //General methods, to be declared in deriving classes if needed
  
  virtual TList *        GetCreateOutputObjects()               { return (new TList)          ; }
  
  virtual void           Init() {;}
  virtual void           InitParameters() ;
   
  virtual void           MakeAnalysisFillAOD()                  { ; }
  
  virtual void           MakeAnalysisFillHistograms()           { ; }  
    
  virtual void           Print(const Option_t * ) const ;
  
  virtual void           Terminate(TList * /*outputList*/)      { ; }

  //Histograms, cuts 
	
  virtual void           AddToHistogramsName(TString add)       { fAddToHistogramsName = add  ; }  
  virtual TString        GetAddedHistogramsStringToName() const { return fAddToHistogramsName ; }
  
  virtual TObjString *   GetAnalysisCuts()                      { return 0x0                  ; }
  virtual TString	       GetBaseParametersList();
  
  //Getters, setters
  virtual Int_t          GetDebug()                       const { return fDebug               ; }
  virtual void           SetDebug(Int_t d)                      { fDebug = d                  ; }
  
  virtual Int_t          GetEventNumber() const ;
   	
  //Centrality
  virtual AliCentrality* GetCentrality()                   const { return fReader->GetCentrality()       ; }
  virtual Int_t          GetEventCentrality()              const { return fReader->GetEventCentrality()  ; }
	
  //Event plane
  virtual AliEventplane* GetEventPlane()                   const { return fReader->GetEventPlane()       ; }           
  virtual TString        GetEventPlaneMethod()             const { return fReader->GetEventPlaneMethod() ; }
  
  //AOD branch
  virtual void           AddAODParticle(AliAODPWG4Particle part) ;
  
  virtual void           ConnectInputOutputAODBranches();
  
  virtual TClonesArray * GetCreateOutputAODBranch() ;
  
  virtual TString        GetInputAODName()                 const { return fInputAODName  ; }
  virtual void           SetInputAODName(TString name)           { fInputAODName = name  ; }	
  
  virtual TString        GetOutputAODName()                const { return fOutputAODName ; }
  virtual void           SetOutputAODName(TString name)          { fNewAOD = kTRUE ; fOutputAODName = name; }
  
  virtual Bool_t         NewOutputAOD()                    const { return fNewAOD        ; }
  
  virtual TString        GetOutputAODClassName()           const { return fOutputAODClassName ; }
  virtual void           SetOutputAODClassName(TString name)     { fOutputAODClassName = name ; }
  
  virtual TString        GetAODObjArrayName()              const { return fAODObjArrayName ; }
  virtual void           SetAODObjArrayName(TString name)        { fAODObjArrayName = name ; }
  
  virtual TClonesArray * GetInputAODBranch()               const { return fInputAODBranch  ; }
  virtual TClonesArray * GetOutputAODBranch()              const { if(fNewAOD) return fOutputAODBranch; else return fInputAODBranch ; }
  virtual TClonesArray * GetAODBranch(TString aodBranchName) const ;
	
  //Track cluster arrays access methods
  virtual TClonesArray*  GetAODCaloClusters()              const ; // Output AOD clusters, not used?
  virtual TClonesArray*  GetAODTracks()                    const ; // Output AOD tracks,   not used?
  virtual AliVCaloCells* GetPHOSCells()                    const { return fReader->GetPHOSCells()  ; }
  virtual AliVCaloCells* GetEMCALCells()                   const { return fReader->GetEMCALCells() ; }
  virtual TObjArray*     GetCTSTracks()                    const ;
  virtual TObjArray*     GetEMCALClusters()                const ;
  virtual TObjArray*     GetPHOSClusters()                 const ;
  
	
  // Common analysis switchs 
  
  virtual Bool_t         IsDataMC()                        const { return fDataMC                ; }
  virtual void           SwitchOnDataMC()                        { fDataMC = kTRUE ; if(!fMCUtils)fMCUtils = new AliMCAnalysisUtils();}
  virtual void           SwitchOffDataMC()                       { fDataMC = kFALSE              ; }
  
  virtual Bool_t         IsFiducialCutOn()                 const { return fCheckFidCut           ; }
  virtual void           SwitchOnFiducialCut()                   { fCheckFidCut = kTRUE;  if(!fFidCut)fFidCut = new AliFiducialCut();}
  virtual void           SwitchOffFiducialCut()                  { fCheckFidCut = kFALSE         ; }
  
  virtual Bool_t         IsCaloPIDOn()                     const { return fCheckCaloPID          ; }
  virtual void           SwitchOnCaloPID()                       { fCheckCaloPID = kTRUE; if(!fCaloPID)fCaloPID = new AliCaloPID();}
  virtual void           SwitchOffCaloPID()                      { fCheckCaloPID = kFALSE        ; }
  
  virtual Bool_t         MakePlotsOn()                     const { return fMakePlots        ; }
  virtual void           SwitchOnPlotsMaking()                   { fMakePlots = kTRUE       ; }
  virtual void           SwitchOffPlotsMaking()                  { fMakePlots = kFALSE      ; }
  
  virtual Bool_t         DoEventSelect()                   const { return fUseSelectEvent   ; }   // Do correlation analysis with different event buffers

  virtual void           SwitchOnEventSelection()                { fUseSelectEvent = kTRUE  ; }
  virtual void           SwitchOffEventSelection()               { fUseSelectEvent = kFALSE ; } 
  
  // Cluster energy/momentum cut
  
  virtual Float_t        GetMaxPt()                        const { return fMaxPt ; }
  virtual Float_t        GetMinPt()                        const { return fMinPt ; }
  virtual void           SetMaxPt(Float_t pt)                    { fMaxPt = pt   ; }
  virtual void           SetMinPt(Float_t pt)                    { fMinPt = pt   ; }
  virtual void           SetPtCutRange(Double_t mi, Double_t ma) { fMaxPt = ma;  fMinPt=mi; }
  
  virtual Float_t        GetMaxEnergy()                    const { return fMaxPt ; }
  virtual Float_t        GetMinEnergy()                    const { return fMinPt ; }
  virtual void           SetMaxEnergy(Float_t e)                 { fMaxPt = e    ; }
  virtual void           SetMinEnergy(Float_t e)                 { fMinPt = e    ; }
  virtual void           SetEnergyCutRange(Double_t mi, Double_t ma) { fMaxPt = ma;   fMinPt = mi; }
  
  //Cluster Pairs Time cut  
  virtual void           SetPairTimeCut(Float_t t)               { fPairTimeCut  = t   ; } //ns
  virtual Float_t        GetPairTimeCut()                  const { return fPairTimeCut ; } //ns
  
  //Getters / Setters for parameters of event buffers
  
  virtual Int_t          GetMultiBin()                     const { return fMultiBin  ; } //number of bins in Multiplicity 
  virtual Int_t          GetNZvertBin()                    const { return fNZvertBin ; } //number of bins in vertex   
  virtual Int_t          GetNRPBin()                       const { return fNrpBin    ; } //number of bins in reaction plain 
  virtual Int_t          GetNCentrBin()                    const { return fNCentrBin ; } //number of bins in centrality
  virtual Int_t          GetNMaxEvMix()                    const { return fNmaxMixEv ; } //maximal number of events for mixin
  virtual Float_t        GetZvertexCut()                   const { return GetReader()->GetZvertexCut();} //cut on vertex position  
  virtual Int_t          GetMaxMulti()                     const { return fMaxMulti  ; }  
  virtual Int_t          GetMinMulti()                     const { return fMinMulti  ; }  
  
  virtual void           SetMultiBin(Int_t n=1)                  { fMultiBin  = n ;} //number of bins in Multiplicity  
  virtual void           SetNZvertBin(Int_t n=1)                 { fNZvertBin = n ;} //number of bins for vertex position
  virtual void           SetNRPBin(Int_t n=1)                    { fNrpBin    = n ;} //number of bins in reaction plain  
  virtual void           SetNCentrBin(Int_t n=1)                 { fNCentrBin = n ;} //number of bins in centrality 
  virtual void           SetNMaxEvMix(Int_t n=20)                { fNmaxMixEv = n ;} //maximal number of events for mixing
  virtual void           SetMultiplicity(Int_t multimin, Int_t multimax) {fMinMulti = multimin ; fMaxMulti = multimax ; }
  
  //Mixed event  
  virtual Int_t           CheckMixedEventVertex(const Int_t caloLabel, const Int_t trackLabel) ;
  virtual AliMixedEvent * GetMixedEvent()                        { return GetReader()->GetMixedEvent()  ; } 
  virtual Int_t           GetNMixedEvent()                 const { return GetReader()->GetNMixedEvent() ; } 
  
  //Vertex methods
  virtual void           GetVertex(Double_t vertex[3])     const { GetReader()->GetVertex(vertex)       ; } 
  virtual Double_t*      GetVertex(const Int_t evtIndex)   const { return GetReader()->GetVertex(evtIndex) ; } 
  virtual void           GetVertex(Double_t vertex[3], const Int_t evtIndex) const { GetReader()->GetVertex(vertex,evtIndex) ; } 

  
  //MULTIPLICITY
  
  virtual Int_t GetTrackMultiplicity()                     const { return fReader->GetTrackMultiplicity() ; }
  
  //VZERO
  
  virtual Int_t GetV0Signal(Int_t i )                      const { return fReader->GetV0Signal(i)         ; }
  
  virtual Int_t GetV0Multiplicity(Int_t i )                const { return fReader->GetV0Multiplicity(i)   ; }
  
  
  //MC event acces methods
  virtual AliStack *                 GetMCStack()          const ;
  
  virtual AliHeader*                 GetMCHeader()         const ;
  
  virtual AliGenEventHeader        * GetMCGenEventHeader() const ;
  
  //Analysis helpers classes pointers setters and getters
  
  virtual AliCaloPID               * GetCaloPID()                { if(!fCaloPID) fCaloPID = new AliCaloPID();               return  fCaloPID ; }
  
  virtual AliCalorimeterUtils      * GetCaloUtils()        const { return fCaloUtils            ; }
  
  virtual AliFiducialCut           * GetFiducialCut()            { if(!fFidCut)  fFidCut  = new AliFiducialCut();           return  fFidCut  ; }
  
  virtual AliHistogramRanges       * GetHistogramRanges()        { if(!fHisto)   fHisto   = new AliHistogramRanges();       return  fHisto   ; }
  
  virtual AliIsolationCut          * GetIsolationCut()           { if(!fIC)      fIC      = new AliIsolationCut();          return  fIC      ; }
  
  virtual AliMCAnalysisUtils       * GetMCAnalysisUtils()        { if(!fMCUtils) fMCUtils = new AliMCAnalysisUtils();       return  fMCUtils ; }
  
  virtual AliNeutralMesonSelection * GetNeutralMesonSelection()  { if(!fNMS)     fNMS     = new AliNeutralMesonSelection(); return  fNMS     ; }
  
  virtual AliCaloTrackReader       * GetReader()           const { return fReader                        ; }
  
  virtual AliEMCALGeometry         * GetEMCALGeometry()    const { return fCaloUtils->GetEMCALGeometry() ; }
  
  virtual AliPHOSGeoUtils          * GetPHOSGeometry()     const { return fCaloUtils->GetPHOSGeometry()  ; }

  virtual void                       SetCaloPID(AliCaloPID * const pid)                             { delete fCaloPID; fCaloPID = pid     ; }
  
  virtual void                       SetCaloUtils(AliCalorimeterUtils * caloutils)                  { fCaloUtils = caloutils              ; }	
  
  virtual void                       SetFiducialCut(AliFiducialCut * const fc)                      { delete fFidCut;  fFidCut  = fc      ; }
  
  virtual void                       SetHistogramRanges(AliHistogramRanges * const hr)              { delete fHisto;   fHisto   = hr      ; }
  
  virtual void                       SetIsolationCut(AliIsolationCut * const ic)                    { delete fIC;      fIC      = ic      ; }
  
  virtual void                       SetMCAnalysisUtils(AliMCAnalysisUtils * const mcutils)         { delete fMCUtils; fMCUtils = mcutils ; }	
  
  virtual void                       SetNeutralMesonSelection(AliNeutralMesonSelection * const nms) { delete fNMS;     fNMS     = nms     ; }
  
  virtual void                       SetReader(AliCaloTrackReader * const reader)                   { fReader = reader                    ; }
  
  //Calorimeter specific access methods and calculations
    
  virtual Bool_t         IsTrackMatched(AliVCluster * cluster, AliVEvent* event) {
   return GetCaloPID()->IsTrackMatched(cluster, fCaloUtils, event) ; } 
  
  virtual Int_t          GetModuleNumberCellIndexes(const Int_t absId, const TString calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const {
	  return fCaloUtils->GetModuleNumberCellIndexes(absId, calo, icol, irow,iRCU) ; }
  
  virtual Int_t          GetModuleNumber(AliAODPWG4Particle * part) const {
	  return fCaloUtils->GetModuleNumber(part, fReader->GetInputEvent()) ; }
  
  virtual Int_t          GetModuleNumber(AliVCluster * cluster)     const {
	  return fCaloUtils->GetModuleNumber(cluster) ; }
  
  virtual AliVCluster*   FindCluster(TObjArray* clusters, const Int_t id, Int_t & iclus, const Int_t first=0) ;

private:    
  
  Bool_t                     fDataMC ;             // Flag to access MC data when using ESD or AOD     
  Int_t                      fDebug ;              // Debug level
  Bool_t                     fCheckFidCut ;        // Do analysis for clusters in defined region         
  Bool_t                     fCheckCaloPID ;       // Do analysis for calorimeters
  Bool_t                     fRecalculateCaloPID ; // Recalculate PID or use PID weights in calorimeters
  Float_t                    fMinPt ;              // Maximum pt of (trigger) particles in the analysis
  Float_t                    fMaxPt ;              // Minimum pt of (trigger) particles in the analysis
  Float_t                    fPairTimeCut;         // Maximum difference between time of cluster pairs (ns)
  Int_t                      fMultiBin ;	         // Number of bins in event container for multiplicity
  Int_t                      fNZvertBin ;	         // Number of bins in event container for vertex position
  Int_t                      fNrpBin ;	           // Number of bins in event container for reaction plain
  Int_t                      fNCentrBin ;	         // Number of bins in event container for centrality
  Int_t                      fNmaxMixEv ;	         // Maximal number of events stored in buffer for mixing
  Int_t                      fMaxMulti ;           // Maximum multiplicity of particles in the analysis
  Int_t                      fMinMulti ;           // Maximum multiplicity of particles in the analysis
  Bool_t                     fUseSelectEvent ;     // Select events based on multiplicity and vertex cuts
  Bool_t                     fMakePlots   ;        // Print plots
    
  TClonesArray*              fInputAODBranch ;     //! Selected input particles branch
  TString                    fInputAODName ;       //  Name of input AOD branch;
  TClonesArray*              fOutputAODBranch ;    //! Selected output particles branch
  Bool_t                     fNewAOD ;             //  Flag, new aod branch added to the analysis or not.
  TString                    fOutputAODName ;      //  Name of output AOD branch;
  TString                    fOutputAODClassName;  //  Type of aod objects to be stored in the TClonesArray (AliAODPWG4Particle, AliAODPWG4ParticleCorrelation ...)	
  TString                    fAODObjArrayName ;    //  Name of ref array kept in a TList in AliAODParticleCorrelation with clusters or track references.
  TString                    fAddToHistogramsName; //  Add this string to histograms name
  
  //Analysis helper classes access pointers
  AliCaloPID               * fCaloPID;             // PID calculation
  AliCalorimeterUtils      * fCaloUtils ;          // Pointer to CalorimeterUtils
  AliFiducialCut           * fFidCut;              // Acceptance cuts
  AliHistogramRanges       * fHisto ;              // Histogram ranges container
  AliIsolationCut          * fIC;                  // Isolation cut 
  AliMCAnalysisUtils       * fMCUtils;             // MonteCarlo Analysis utils 
  AliNeutralMesonSelection * fNMS;                 // Neutral Meson Selection
  AliCaloTrackReader       * fReader;              // Acces to ESD/AOD/MC data

  AliAnaCaloTrackCorrBaseClass(              const AliAnaCaloTrackCorrBaseClass & bc) ; // cpy ctor
  AliAnaCaloTrackCorrBaseClass & operator = (const AliAnaCaloTrackCorrBaseClass & bc) ; // cpy assignment
  
  ClassDef(AliAnaCaloTrackCorrBaseClass,20)
} ;


#endif //ALIANACALOTRACKCORRBASECLASS_H



