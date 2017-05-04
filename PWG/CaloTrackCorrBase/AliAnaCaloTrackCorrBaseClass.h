#ifndef ALIANACALOTRACKCORRBASECLASS_H
#define ALIANACALOTRACKCORRBASECLASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaCaloTrackCorrBaseClass
/// \ingroup CaloTrackCorrelationsBase
/// \brief Base class for CaloTrackCorr analysis algorithms
/// 
/// All analysis classes in the package CaloTrackCorrelations
/// must derive from this class.
/// Common data members and methods are defined here
/// The main methods of the class that must be declared in the 
/// daughter analysis classes are
/// * *GetCreateOutputObjects()*: All histograms are declared here 
/// * *MakeAnalysisFillAOD()*: If a particle object is created/filtered, it has to be done here
/// * *MakeAnalysisFillHistograms()*: The loop on clusters/tracks/filtered particles is done here, histograms filled
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
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
class AliMCEvent ; 
class AliHeader ; 
class AliGenEventHeader ; 
class AliEMCALGeometry;
class AliPHOSGeoUtils;
class AliCentrality;
class AliMultSelection;
class AliEventplane;
#include "AliAnalysisManager.h"
#include "AliLog.h"

// Jets
class AliAODJetEventBackground;

class AliAnaCaloTrackCorrBaseClass : public TObject {
	
public:   
  
  AliAnaCaloTrackCorrBaseClass() ;          // default ctor
  virtual ~AliAnaCaloTrackCorrBaseClass() ; // virtual dtor
  
  // General methods, to be declared in deriving classes if needed
  
  virtual TList *        GetCreateOutputObjects()               { return (new TList)          ; }
  
  virtual void           Init()                                 { ; }
  virtual void           InitDebug()      ;
  virtual void           InitParameters() ;
  virtual void           InitCaloParameters() ;
  
  virtual void           FillEventMixPool()                     { ; }

  virtual void           MakeAnalysisFillAOD()                  { ; }
  
  virtual void           MakeAnalysisFillHistograms()           { ; }  
    
  virtual void           Print(const Option_t * ) const ;
  
  virtual void           Terminate(TList * /*outputList*/)      { ; }

  // Histograms, cuts 
	
  virtual void           AddToHistogramsName(TString add)       { fAddToHistogramsName = add  ; }
  virtual TString        GetAddedHistogramsStringToName() const { return fAddToHistogramsName ; }
  
  virtual TObjString *   GetAnalysisCuts()                      { return 0x0                  ; }
  virtual TString	       GetBaseParametersList();
  
  // Getters, setters
  
  virtual Int_t          GetDebug()                       const { return fDebug               ; }
  virtual void           SetDebug(Int_t d)                      { fDebug = d                  ; }
  
  virtual Int_t          GetEventNumber() const ;
  
  // Centrality, multiplicity selection
  
  virtual Int_t             GetTrackMultiplicity()        const { return fReader->GetTrackMultiplicity() ; }
  virtual AliCentrality*    GetCentrality()               const { return fReader->GetCentrality() ; }
  virtual AliMultSelection* GetMultSelCen()               const { return fReader->GetMultSelCen() ; }
  virtual Int_t             GetEventCentrality()          const { if(fUseTrackMultBins)
                                                                        return GetTrackMultiplicity();
                                                                   else return fReader->GetEventCentrality(); }
  
  // Event plane
  
  virtual AliEventplane* GetEventPlane()                   const { return fReader->GetEventPlane()       ; }           
  virtual Double_t       GetEventPlaneAngle()              const { return fReader->GetEventPlaneAngle()  ; }           
  virtual TString        GetEventPlaneMethod()             const { return fReader->GetEventPlaneMethod() ; }
  
  // AOD branch
  
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
  virtual TClonesArray * GetAODBranch(const TString & aodBranchName) const ;
	
  // Track cluster arrays access methods
  
  virtual TClonesArray*  GetAODCaloClusters()              const ; // Output AOD clusters, not used?
  virtual TClonesArray*  GetAODTracks()                    const ; // Output AOD tracks,   not used?
  virtual AliVCaloCells* GetPHOSCells()                    const { return fReader->GetPHOSCells()  ; }
  virtual AliVCaloCells* GetEMCALCells()                   const { return fReader->GetEMCALCells() ; }
  virtual TObjArray*     GetCTSTracks()                    const ;
  virtual TObjArray*     GetEMCALClusters()                const ;
  virtual TObjArray*     GetPHOSClusters()                 const ;
  
  // Jets
  
  virtual TClonesArray*  GetNonStandardJets()              const { return fReader->GetNonStandardJets() ;}
  virtual AliAODJetEventBackground*  GetBackgroundJets()   const { return fReader->GetBackgroundJets() ;}

  // Common analysis switchs 
  
  /// Set the tag identifing the main detector used in the analysis 
  enum detector 
  { 
    kEMCAL    = AliFiducialCut::kEMCAL,   /// EMCAL (+DCAL) 
    kPHOS     = AliFiducialCut::kPHOS,    /// PHOS 
    kCTS      = AliFiducialCut::kCTS,     /// Tracking
    kDCAL     = AliFiducialCut::kDCAL,    /// DCal, not used so far, just in case
    kDCALPHOS = AliFiducialCut::kDCALPHOS /// DCal+PHOS, not used so far, just in case
  } ;

  virtual Int_t          GetCalorimeter()                 const  { return fCalorimeter          ; }
  virtual TString        GetCalorimeterString()           const  { return fCalorimeterString    ; }
  virtual void           SetCalorimeter(TString & calo);
  virtual void           SetCalorimeter(Int_t calo) ;

  virtual Bool_t         IsDataMC()                        const { return fDataMC                ; }
  virtual void           SwitchOnDataMC()                        { fDataMC = kTRUE ;
                                                                   if(!fMCUtils) fMCUtils = new AliMCAnalysisUtils() ; }
  virtual void           SwitchOffDataMC()                       { fDataMC = kFALSE              ; }
  
  virtual Bool_t         IsFiducialCutOn()                 const { return fCheckFidCut           ; }
  virtual void           SwitchOnFiducialCut()                   { fCheckFidCut = kTRUE ;
                                                                   if(!fFidCut)  fFidCut   = new AliFiducialCut()     ; }
  virtual void           SwitchOffFiducialCut()                  { fCheckFidCut = kFALSE         ; }

  virtual Bool_t         IsRealCaloAcceptanceOn()          const { return fCheckRealCaloAcc      ; }
  virtual void           SwitchOnRealCaloAcceptance()            { fCheckRealCaloAcc = kTRUE;  }
  virtual void           SwitchOffRealCaloAcceptance()           { fCheckRealCaloAcc = kFALSE    ; }
  
  virtual Bool_t         IsCaloPIDOn()                     const { return fCheckCaloPID          ; }
  virtual void           SwitchOnCaloPID()                       { fCheckCaloPID = kTRUE ;
                                                                   if(!fCaloPID)  fCaloPID = new AliCaloPID()         ; }
  virtual void           SwitchOffCaloPID()                      { fCheckCaloPID = kFALSE        ; }
  
  virtual Bool_t         MakePlotsOn()                     const { return fMakePlots             ; }
  virtual void           SwitchOnPlotsMaking()                   { fMakePlots = kTRUE            ; }
  virtual void           SwitchOffPlotsMaking()                  { fMakePlots = kFALSE           ; }
  
  virtual Bool_t         IsPileUpAnalysisOn()              const { return fFillPileUpHistograms  ; }
  virtual void           SwitchOnFillPileUpHistograms()          { fFillPileUpHistograms = kTRUE ; }
  virtual void           SwitchOffFillPileUpHistograms()         { fFillPileUpHistograms = kFALSE; }
  
  virtual Bool_t         IsHighMultiplicityAnalysisOn()     const { return fFillHighMultHistograms   ; }
  virtual void           SwitchOnFillHighMultiplicityHistograms() { fFillHighMultHistograms = kTRUE  ; }
  virtual void           SwitchOffFillHighMultiplicityHistograms(){ fFillHighMultHistograms = kFALSE ; }

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
  
  // Cluster Pairs Time cut  
  
  virtual void           SetPairTimeCut(Float_t t)               { fPairTimeCut  = t   ; } /// Time cut in ns
  virtual Float_t        GetPairTimeCut()                  const { return fPairTimeCut ; } /// Time cut in ns
  
  // Number of TRD modules in front of EMCAL (year <=2012)
  
  Int_t                  GetFirstSMCoveredByTRD()          const { return fTRDSMCovered ; }
  void                   SetFirstSMCoveredByTRD(Int_t n)         { fTRDSMCovered    = n ; }
  
  // Getters / Setters for parameters of event buffers
  
  virtual Int_t          GetNZvertBin()                    const { return fNZvertBin ; }    /// Number of bins in vertex
  virtual Int_t          GetNRPBin()                       const { return fNrpBin    ; }    /// Number of bins in reaction plain
  virtual Int_t          GetNCentrBin()                    const { return fNCentrBin ; }    /// Number of bins in centrality
  virtual Int_t          GetNTrackMultBin()                const { return GetNCentrBin(); } /// Number of bins in track multiplicity
  virtual Int_t          GetNMaxEvMix()                    const { return fNmaxMixEv ; }    /// Maximal number of events for mixin
  virtual Float_t        GetZvertexCut()                   const { return GetReader()->GetZvertexCut();} /// Cut on vertex position
  virtual Int_t          GetTrackMultiplicityBin()         const ;
  virtual Int_t          GetEventCentralityBin()           const ;
  virtual Int_t          GetEventRPBin()                   const ;
  virtual Int_t          GetEventVzBin()                   const ;
  virtual Int_t          GetEventMixBin()                  const ;
  virtual Int_t          GetEventMixBin(Int_t iCen, Int_t iVz, Int_t iRP) const;
  
  virtual Double_t       GetEventWeight()                  const { return GetReader()->GetEventWeight()      ; }
  virtual Double_t       GetParticlePtWeight(Float_t pt, Int_t pdg, TString genName, Int_t igen) 
    const { return (GetReader()->GetWeightUtils())->GetParticlePtWeight(pt, pdg, genName, igen)  ; }
  
  virtual void           SetNZvertBin(Int_t n = 1 )              { fNZvertBin = n ; if(n < 1) fNZvertBin = 1 ; } /// Number of bins for vertex position
  virtual void           SetNRPBin   (Int_t n = 1 )              { fNrpBin    = n ; if(n < 1) fNrpBin    = 1 ; } /// Number of bins in reaction plain
  virtual void           SetNCentrBin(Int_t n = 1 )              { fNCentrBin = n ; if(n < 1) fNCentrBin = 1 ; } /// Number of bins in centrality
  virtual void           SetNTrackMultBin(Int_t n = 1 )          { SetNCentrBin(n); }                            /// Number of bins in track multiplicity
  virtual void           SetNMaxEvMix(Int_t n = 20)              { fNmaxMixEv = n ; if(n < 1) fNmaxMixEv = 1 ; } /// Maximal number of events for mixing
  virtual void           SetTrackMultiplicityBin(Int_t bin, Int_t mult) { if(bin < 20) fTrackMultBins[bin] = mult ; }
  
  virtual void           SwitchOnTrackMultBins()                 { fUseTrackMultBins = kTRUE  ; }
  virtual void           SwitchOffTrackMultBins()                { fUseTrackMultBins = kFALSE ; }
  
  virtual void           SwitchOnOwnMix()                        { fDoOwnMix         = kTRUE  ; }
  virtual void           SwitchOffOwnMix()                       { fDoOwnMix         = kFALSE ; }
  
  virtual Bool_t         DoOwnMix()                        const { return fDoOwnMix           ; }
  virtual Bool_t         UseTrackMultBins()                const { return fUseTrackMultBins   ; }

  // Mixed event 
  
  virtual Int_t           CheckMixedEventVertex(Int_t caloLabel, Int_t trackLabel) ;
  
  virtual AliMixedEvent * GetMixedEvent()                  const { return GetReader()->GetMixedEvent()     ; }
  
  virtual Int_t           GetNMixedEvent()                 const { return GetReader()->GetNMixedEvent()    ; } 
  
  // Vertex methods
  
  virtual void           GetVertex(Double_t vertex[3])     const { GetReader()->GetVertex(vertex)          ; } 
  
  virtual Double_t*      GetVertex(Int_t evtIndex)         const { return GetReader()->GetVertex(evtIndex) ; }
  
  virtual void           GetVertex(Double_t vertex[3],
                                   Int_t evtIndex)         const { GetReader()->GetVertex(vertex,evtIndex) ; }
  // VZERO
  
  virtual Int_t          GetV0Signal(Int_t i )             const { return fReader->GetV0Signal(i)          ; }
  
  virtual Int_t          GetV0Multiplicity(Int_t i )       const { return fReader->GetV0Multiplicity(i)    ; }
  
  
  // Some utilities dealing with cluster angles
  
  /// Shift phi angle in case of negative value 360 degrees. Example TLorenzVector::Phi defined in -pi to pi
  Float_t                GetPhi  (Float_t phi)             const { if ( phi < 0 ) phi += TMath::TwoPi() ; return phi ; }
  
  Float_t                DegToRad(Float_t deg)             const { deg *= TMath::DegToRad(); return deg ; }
  
  Float_t                RadToDeg(Float_t rad)             const { rad *= TMath::RadToDeg(); return rad ; }
  
  // Calorimeter specific access methods and calculations
  
  virtual Bool_t         IsTrackMatched(AliVCluster * cluster, AliVEvent* event) 
  { return GetCaloPID()->IsTrackMatched(cluster, fCaloUtils, event)             ; } 
  
  virtual Int_t          GetModuleNumberCellIndexes(Int_t absId, Int_t calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const 
  { return fCaloUtils->GetModuleNumberCellIndexes(absId, calo, icol, irow,iRCU) ; }

  virtual Int_t          GetModuleNumberCellIndexesAbsCaloMap(Int_t absId, Int_t calo, 
                                                              Int_t & icol, Int_t & irow, Int_t &iRCU,
                                                              Int_t & icolAbs, Int_t & irowAbs) const 
  { return fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(absId, calo, icol, irow,iRCU,icolAbs,irowAbs) ; }
  
  virtual Int_t          GetModuleNumber(AliAODPWG4Particle * part) const 
  { return fCaloUtils->GetModuleNumber(part, fReader->GetInputEvent())          ; }
  
  virtual Int_t          GetModuleNumber(AliVCluster * cluster)     const 
  { return fCaloUtils->GetModuleNumber(cluster)                                 ; }
  
  virtual AliVCluster*   FindCluster(TObjArray* clusters, Int_t clId, Int_t & iclus, Int_t first = 0 ) ;

  
  // MC event acces methods
  
  virtual AliMCEvent *               GetMC()               const ;
  
  virtual AliHeader*                 GetMCHeader()         const ;
  
  virtual AliGenEventHeader        * GetMCGenEventHeader() const ;
  
  // Analysis helpers classes pointers setters and getters
  
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

  virtual void                       SetCaloPID(AliCaloPID * pid)                                   { delete fCaloPID; fCaloPID = pid     ; }
  
  virtual void                       SetCaloUtils(AliCalorimeterUtils * caloutils)                  { fCaloUtils = caloutils              ; }	
  
  virtual void                       SetFiducialCut(AliFiducialCut * fc)                            { delete fFidCut;  fFidCut  = fc      ; }
  
  virtual void                       SetHistogramRanges(AliHistogramRanges * hr)                    { delete fHisto;   fHisto   = hr      ; }
  
  virtual void                       SetIsolationCut(AliIsolationCut * ic)                          { delete fIC;      fIC      = ic      ; }
  
  virtual void                       SetMCAnalysisUtils(AliMCAnalysisUtils * mcutils)               { delete fMCUtils; fMCUtils = mcutils ; }
  
  virtual void                       SetNeutralMesonSelection(AliNeutralMesonSelection * const nms) { delete fNMS;     fNMS     = nms     ; }
  
  virtual void                       SetReader(AliCaloTrackReader * reader)                         { fReader = reader                    ; }
  
  // Cocktail generator studies
  
  void                               SwitchOnStudyClusterOverlapsPerGenerator()   { fStudyClusterOverlapsPerGenerator = kTRUE  ; }
  void                               SwitchOffStudyClusterOverlapsPerGenerator()  { fStudyClusterOverlapsPerGenerator = kFALSE ; }  
  Bool_t                             IsStudyClusterOverlapsPerGeneratorOn() const {  return fStudyClusterOverlapsPerGenerator  ; }  
  
  void                               SetNCocktailGenNamesToCheck(Int_t nb)        { fNCocktailGenNames = nb   ; }
  Int_t                              GetNCocktailGenNamesToCheck()          const { return fNCocktailGenNames ; }
  
  void                               SetCocktailGenNameToCheck(Int_t i, TString v){ if(i < 10) fCocktailGenNames[i] = v   ; }
  TString                            GetCocktailGenNameToCheck(Int_t i)     const { if(i < 10) return fCocktailGenNames[i]; 
    else       return ""                  ; }
  void                               SetCocktailGenIndexToCheck(Int_t i, Int_t v) { if(i < 10) fCocktailGenIndeces[i] = v ; }
  Int_t                              GetCocktailGenIndexToCheck(Int_t i)    const { if(i < 10) return fCocktailGenIndeces[i]; 
    else       return -1                  ; }
  
  Int_t                              GetCocktailGeneratorBackgroundTag(AliVCluster * clus, Int_t mctag,
                                                                       TString & genName   , Int_t & index,
                                                                       TString & genNameBkg, Int_t & indexBkg);
  
protected:

  Int_t                      fNModules    ;        ///<  Number of EMCAL/PHOS modules to use in analysis, set in CaloUtils
  Int_t                      fFirstModule ;        ///<  First EMCAL/PHOS module, set in CaloUtils or depending fidutial cuts
  Int_t                      fLastModule  ;        ///<  Last EMCAL/PHOS module, set in CaloUtils or depending fidutial cuts
  Int_t                      fNRCU        ;        ///<  Number of EMCAL/PHOS RCU
  Int_t                      fNMaxCols    ;        ///<  Number of EMCAL/PHOS columns per SM
  Int_t                      fNMaxRows    ;        ///<  Number of EMCAL/PHOS rows per SM
  Int_t                      fNMaxColsFull;        ///<  Number of EMCAL/PHOS columns full detector
  Int_t                      fNMaxRowsFull;        ///<  Number of EMCAL/PHOS rows full detector
  Int_t                      fNMaxRowsFullMax;     ///<  First of EMCAL/PHOS rows full detector
  Int_t                      fNMaxRowsFullMin;     ///<  Last of EMCAL/PHOS rows full detector

private:    
  
  Bool_t                     fDataMC ;             ///< Flag to access MC data when using ESD or AOD.    
  Int_t                      fDebug ;              ///< Debug level.
  Int_t                      fCalorimeter ;        ///< Calorimeter selection.
  TString                    fCalorimeterString ;  ///< Calorimeter selection.
  Bool_t                     fCheckFidCut ;        ///< Do analysis for clusters in defined region.
  Bool_t                     fCheckRealCaloAcc ;   ///< When analysis of MC particle kinematics, check their hit in Calorimeter in Real Geometry or use AliFiducialCut.
  Bool_t                     fCheckCaloPID ;       ///< Do analysis for calorimeters.
  Bool_t                     fRecalculateCaloPID ; ///< Recalculate PID or use PID weights in calorimeters.
  Float_t                    fMinPt ;              ///< Maximum pt of (trigger) particles in the analysis.
  Float_t                    fMaxPt ;              ///< Minimum pt of (trigger) particles in the analysis.
  Float_t                    fPairTimeCut;         ///< Maximum difference between time of cluster pairs (ns).
  Int_t                      fTRDSMCovered;        ///< From which SM EMCal is covered by TRD.
  
  Int_t                      fNZvertBin ;	         ///< Number of bins in event container for vertex position.
  Int_t                      fNrpBin ;	           ///< Number of bins in event container for reaction plain.
  Int_t                      fNCentrBin ;	         ///< Number of bins in event container for centrality.
  Int_t                      fNmaxMixEv ;	         ///< Maximal number of events stored in buffer for mixing.
  Bool_t                     fDoOwnMix;            ///< Do combinatorial background not the one provided by the frame.
  Bool_t                     fUseTrackMultBins;    ///< Use track multiplicity and not centrality bins in mixing.
  Int_t                      fTrackMultBins[20];   ///< Multiplicity bins limits. Number of bins set with SetNTrackMult() that calls SetNCentrBin().
  Bool_t                     fFillPileUpHistograms;   ///< Fill pile-up related histograms.
  Bool_t                     fFillHighMultHistograms; ///< Histograms with centrality and event plane for triggers pT.
  Bool_t                     fMakePlots   ;        ///< Print plots.
    
  TClonesArray*              fInputAODBranch ;     //!<! Selected input particles branch.
  TString                    fInputAODName ;       ///<  Name of input AOD branch.
  TClonesArray*              fOutputAODBranch ;    //!<! Selected output particles branch.
  Bool_t                     fNewAOD ;             ///<  Flag, new aod branch added to the analysis or not.
  TString                    fOutputAODName ;      ///<  Name of output AOD branch.
  TString                    fOutputAODClassName;  ///<  Type of aod objects to be stored in the TClonesArray (AliAODPWG4Particle, AliAODPWG4ParticleCorrelation ...).	
  TString                    fAODObjArrayName ;    ///<  Name of ref array kept in a TList in AliAODParticleCorrelation with clusters or track. references.
  TString                    fAddToHistogramsName; ///<  Add this string to histograms name.
  
  // Analysis helper classes access pointers
  AliCaloPID               * fCaloPID;             ///< PID calculation utils.
  AliCalorimeterUtils      * fCaloUtils ;          ///< Pointer to Calorimeter Utils.
  AliFiducialCut           * fFidCut;              ///< Acceptance cuts detector dependent.
  AliHistogramRanges       * fHisto ;              ///< Histogram ranges container.
  AliIsolationCut          * fIC;                  ///< Isolation cut utils. 
  AliMCAnalysisUtils       * fMCUtils;             ///< MonteCarlo Analysis utils. 
  AliNeutralMesonSelection * fNMS;                 ///< Neutral Meson Selection utities.
  AliCaloTrackReader       * fReader;              ///< Access to ESD/AOD/MC data and other utilities.

  // Cocktail generator studies
  Bool_t                     fStudyClusterOverlapsPerGenerator; ///<  In case of coctail generators, check the content of the cluster
  Int_t                      fNCocktailGenNames;                ///<  Number of generators to study
  TString                    fCocktailGenNames[10];             ///<  Array with name of generators to study, first must be always empty
  Int_t                      fCocktailGenIndeces[10];           ///<  Array with indeces of generators to study
  
  /// Copy constructor not implemented.
  AliAnaCaloTrackCorrBaseClass(              const AliAnaCaloTrackCorrBaseClass & bc) ; 
  
  /// Assignment operator not implemented.
  AliAnaCaloTrackCorrBaseClass & operator = (const AliAnaCaloTrackCorrBaseClass & bc) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliAnaCaloTrackCorrBaseClass,29) ;
  /// \endcond

} ;

#endif //ALIANACALOTRACKCORRBASECLASS_H



