#ifndef ALIANAPARTCORRBASECLASS_H
#define ALIANAPARTCORRBASECLASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

  //_________________________________________________________________________
  // Base class for analysis algorithms
  //-- Author: Gustavo Conesa (INFN-LNF)
//-Add the possibality for event selection analysis based on vertex and multiplicity bins (Yaxian Mao, 10/10/2010)
#include <cstdlib>

  //ROOT
class TClonesArray ;
class TObjArray ;
#include <TList.h> 
#include <TObject.h>
class TObjString;

  //Analysis
class AliVCaloCells;
#include "AliCaloTrackReader.h"   
#include "AliCaloPID.h"
#include "AliFiducialCut.h"
#include "AliIsolationCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliCalorimeterUtils.h" 
class AliStack ; 
class AliHeader ; 
class AliGenEventHeader ; 
#include "AliAODPWG4ParticleCorrelation.h"
class AliEMCALGeoUtils;
class AliPHOSGeoUtils;
#include "AliMixedEvent.h" 

class AliAnaPartCorrBaseClass : public TObject {
	
public:   
  AliAnaPartCorrBaseClass() ; // default ctor
  virtual ~AliAnaPartCorrBaseClass() ; //virtual dtor
  
private:
  AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & g) ; // cpy ctor
  AliAnaPartCorrBaseClass & operator = (const AliAnaPartCorrBaseClass & g) ;//cpy assignment
  
public:

  virtual void AddAODParticle(AliAODPWG4Particle part) ;
  
  virtual void ConnectInputOutputAODBranches();
  
  virtual TList * GetCreateOutputObjects()      { return (new TList) ;}
	
  virtual void AddToHistogramsName(TString add) { fAddToHistogramsName = add; }  
  virtual TString GetAddedHistogramsStringToName() {return fAddToHistogramsName ;}
  
  virtual void Init() {;}
  virtual void InitParameters() ;
  
  virtual void Print(const Option_t * ) const ;
  
  virtual void MakeAnalysisFillAOD()  {;}
  
  virtual void MakeAnalysisFillHistograms() {;}
  	
  virtual TObjString * GetAnalysisCuts() {return 0x0;}
	
  virtual Int_t GetDebug() const  { return fDebug ; }
  virtual void SetDebug(Int_t d)   { fDebug = d ; }
  
  virtual Int_t GetEventNumber() const ;
  
  virtual AliCaloTrackReader * GetReader() const {return fReader ; }
  virtual void SetReader(AliCaloTrackReader * const reader) { fReader = reader ; }
  
  Int_t GetTrackMultiplicity() const {return fReader->GetTrackMultiplicity();}
  
  //Calorimeter helper class access methods
  AliEMCALGeoUtils *  GetEMCALGeometry() const { return fCaloUtils->GetEMCALGeometry(); }
  AliPHOSGeoUtils  *  GetPHOSGeometry()  const { return fCaloUtils->GetPHOSGeometry() ; }
  
  Int_t GetModuleNumberCellIndexes(const Int_t absId, const TString calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const {
	  return fCaloUtils->GetModuleNumberCellIndexes(absId, calo, icol, irow,iRCU);}
  Int_t GetModuleNumber(AliAODPWG4Particle * part) const {
	  return fCaloUtils->GetModuleNumber(part, fReader->GetInputEvent());}
  Int_t GetModuleNumber(AliVCluster * cluster) const {
	  return fCaloUtils->GetModuleNumber(cluster);}
 	
  virtual void Terminate(TList * /*outputList*/) {;}
	
  //analysis AOD branch
  virtual TClonesArray * GetCreateOutputAODBranch() ;
  virtual TString GetInputAODName() const {return fInputAODName ; }
  virtual void SetInputAODName(TString name)   { fInputAODName = name; }	
  virtual TString GetOutputAODName()  const {return fOutputAODName ; }
  virtual void SetOutputAODName(TString name)   { fNewAOD = kTRUE ; fOutputAODName = name; }
  virtual Bool_t NewOutputAOD() const {return fNewAOD;}
  virtual TString GetOutputAODClassName() const {return fOutputAODClassName;}
  virtual void SetOutputAODClassName(TString name) {fOutputAODClassName = name; }
  virtual AliCalorimeterUtils * GetCaloUtils() const {return fCaloUtils ; }
  void SetCaloUtils(AliCalorimeterUtils * caloutils) { fCaloUtils = caloutils ; }	
	
  virtual TString GetAODObjArrayName() const {return fAODObjArrayName;}
  virtual void SetAODObjArrayName(TString name) {fAODObjArrayName = name; }
  
  virtual TClonesArray* GetInputAODBranch() const {return fInputAODBranch ;}
  virtual TClonesArray* GetOutputAODBranch() const {if(fNewAOD) return fOutputAODBranch; else return fInputAODBranch ;}
  virtual TClonesArray* GetAODBranch(TString aodBranchName) const ;
	
  virtual TClonesArray* GetAODCaloClusters() const ;
  virtual TClonesArray* GetAODTracks() const ;	
  virtual AliVCaloCells* GetPHOSCells()  const {return fReader->GetPHOSCells()  ;}
  virtual AliVCaloCells* GetEMCALCells() const {return fReader->GetEMCALCells() ;}

  virtual TObjArray* GetAODCTS() const ;
  virtual TObjArray* GetAODEMCAL() const ;
  virtual TObjArray* GetAODPHOS() const ;
  
  virtual TString	GetBaseParametersList();
    
  virtual AliStack * GetMCStack() const ;
  virtual AliHeader* GetMCHeader() const ;
  virtual AliGenEventHeader* GetMCGenEventHeader() const ;
  
  //Analysis helpers classes pointers setters and getters
  virtual AliCaloPID * GetCaloPID() {if(!fCaloPID) fCaloPID = new AliCaloPID(); return  fCaloPID ;}
  virtual void SetCaloPID(AliCaloPID * const pid) { fCaloPID = pid ;}
  
  virtual AliFiducialCut * GetFiducialCut() {if(!fFidCut) fFidCut = new AliFiducialCut(); return  fFidCut ;}
  virtual void SetFiducialCut(AliFiducialCut * const fc) { fFidCut = fc ;}
  
  virtual AliIsolationCut * GetIsolationCut() {if(!fIC) fIC = new AliIsolationCut();  return  fIC ;}
  virtual void SetIsolationCut(AliIsolationCut * const ic) { fIC = ic ;}
  
  virtual AliMCAnalysisUtils * GetMCAnalysisUtils()  {if(!fMCUtils) fMCUtils = new AliMCAnalysisUtils(); return  fMCUtils ;}
  virtual void SetMCAnalysisUtils(AliMCAnalysisUtils * const mcutils) { fMCUtils = mcutils ;}	
  
  virtual AliNeutralMesonSelection * GetNeutralMesonSelection()  {if(!fNMS) fNMS = new AliNeutralMesonSelection(); return  fNMS ;}
  virtual void SetNeutralMesonSelection(AliNeutralMesonSelection * const nms) { fNMS = nms ;}
	
  virtual Bool_t     IsDataMC()       {return fDataMC ; }
  virtual void SwitchOnDataMC()       {fDataMC = kTRUE ; if(!fMCUtils)fMCUtils = new AliMCAnalysisUtils();}
  virtual void SwitchOffDataMC()      {fDataMC = kFALSE ; }
  
  virtual Bool_t IsFiducialCutOn()       { return fCheckFidCut ; }
  virtual void SwitchOnFiducialCut()     { fCheckFidCut = kTRUE;  if(!fFidCut)fFidCut = new AliFiducialCut();}
  virtual void SwitchOffFiducialCut()    { fCheckFidCut = kFALSE;}
  
  virtual Bool_t IsCaloPIDOn()       { return fCheckCaloPID ; }
  virtual void SwitchOnCaloPID()     { fCheckCaloPID = kTRUE; if(!fCaloPID)fCaloPID = new AliCaloPID();}
  virtual void SwitchOffCaloPID()    { fCheckCaloPID = kFALSE;}
  
  virtual Bool_t IsCaloPIDRecalculationOn()       { return fRecalculateCaloPID ; }
  virtual void SwitchOnCaloPIDRecalculation()     { fRecalculateCaloPID  = kTRUE;}
  virtual void SwitchOffCaloPIDRecalculation()    { fRecalculateCaloPID  = kFALSE;}
  
  virtual Float_t    GetMaxPt()         const {return fMaxPt ; }
  virtual Float_t    GetMinPt()         const {return fMinPt ; }
  virtual void SetMaxPt(Float_t pt)           {fMaxPt = pt ; }
  virtual void SetMinPt(Float_t pt)           {fMinPt = pt ; }
  virtual void SetPtCutRange(Double_t ptmin, Double_t ptmax)
  {  fMaxPt=ptmax;   fMinPt=ptmin;}
  //Setters for parameters of event buffers
  virtual void SetMultiBin(Int_t n=1) {fMultiBin=n ;} //number of bins in Multiplicity  
  virtual void SetNZvertBin(Int_t n=1) {fNZvertBin=n ;} //number of bins for vertex position
  virtual void SetNRPBin(Int_t n=1)    {fNrpBin=n ;}    //number of bins in reaction plain  
  virtual void SetZvertexCut(Float_t zcut=40.){fZvtxCut=zcut ;} //cut on vertex position
  virtual void SetMultiplicity(Int_t multimin, Int_t multimax) {fMinMulti = multimin ; fMaxMulti = multimax ; }
  virtual void SwitchOnEventSelection()    {fUseSelectEvent = kTRUE ; }
  virtual void SwitchOffEventSelection()   {fUseSelectEvent = kFALSE ; } 
  //Getters for event selection
  virtual Int_t GetMultiBin()  const       {return fMultiBin ;} //number of bins in Multiplicity 
  virtual Int_t GetNZvertBin() const       {return fNZvertBin ;} //number of bins in vertex   
  virtual Int_t GetNRPBin()    const       {return fNrpBin ;}    //number of bins in reaction plain 
  //Getters for event selection
  virtual Float_t GetZvertexCut() const {return fZvtxCut ;} //cut on vertex position  
  virtual Int_t GetMaxMulti()     const {return fMaxMulti  ; }  
  virtual Int_t GetMinMulti()     const {return fMinMulti  ; }  
  
  // Do correlation analysis with different event buffers
  virtual Bool_t DoEventSelect() const {return fUseSelectEvent ; }
  
  //Histogrammes setters and getters
  //Pt, Energy 
  virtual void SetHistoPtRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPtBins = n ;
    fHistoPtMax = max ;
    fHistoPtMin = min ;
  }
  
  virtual Int_t   GetHistoPtBins()  const { return fHistoPtBins; }
  virtual Float_t GetHistoPtMin()   const { return fHistoPtMin ; }
  virtual Float_t GetHistoPtMax()   const { return fHistoPtMax ; }
  
    //Azimuthal angle
  virtual void SetHistoPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPhiBins  = n ;
    fHistoPhiMax   = max ;
    fHistoPhiMin   = min ;
  }
  
  virtual Int_t   GetHistoPhiBins()  const { return fHistoPhiBins; }
  virtual Float_t GetHistoPhiMin()   const { return fHistoPhiMin ; }
  virtual Float_t GetHistoPhiMax()   const { return fHistoPhiMax ; }
  
  //Pseudorapidity-rapidity
  virtual void SetHistoEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoEtaBins = n ;
    fHistoEtaMax  = max ;
    fHistoEtaMin  = min ;
  }
  
  virtual Int_t   GetHistoEtaBins()  const { return fHistoEtaBins; }
  virtual Float_t GetHistoEtaMin()   const { return fHistoEtaMin ; }
  virtual Float_t GetHistoEtaMax()   const { return fHistoEtaMax ; }
  
  //Mass
  virtual void SetHistoMassRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoMassBins = n ;
    fHistoMassMax  = max ;
    fHistoMassMin  = min ;
  }
	
  virtual Int_t   GetHistoMassBins()  const { return fHistoMassBins ; }
  virtual Float_t GetHistoMassMin()   const { return fHistoMassMin ; }
  virtual Float_t GetHistoMassMax()   const { return fHistoMassMax ; }
	
  //Asymetry
  virtual void SetHistoAsymmetryRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoAsymBins = n ;
    fHistoAsymMax  = max ;
    fHistoAsymMin  = min ;
  }
	
  virtual Int_t   GetHistoAsymmetryBins()  const { return fHistoAsymBins ; }
  virtual Float_t GetHistoAsymmetryMin()   const { return fHistoAsymMin ; }
  virtual Float_t GetHistoAsymmetryMax()   const { return fHistoAsymMax ; }	
  
  virtual AliMixedEvent * GetMixedEvent()          { return GetReader()->GetMixedEvent() ; } 
  virtual Int_t           GetNMixedEvent()   const { return GetReader()->GetNMixedEvent() ; } 
  
  virtual void      GetVertex(Double_t vertex[3])   const { GetReader()->GetVertex(vertex) ; } 
  virtual void      GetVertex(Double_t vertex[3],const Int_t evtIndex) const { GetReader()->GetVertex(vertex,evtIndex) ; } 
  virtual Double_t* GetVertex(const Int_t evtIndex) const { return GetReader()->GetVertex(evtIndex) ; } 

	virtual Bool_t IsTrackMatched(const AliVCluster * cluster) const { return fCaloPID->IsTrackMatched(cluster); } 
  
private:    
  
  Bool_t  fDataMC ;             // Flag to access MC data when using ESD or AOD     
  Int_t   fDebug ;              // Debug level
  Bool_t  fCheckFidCut ;        // Do analysis for clusters in defined region         
  Bool_t  fCheckCaloPID ;       // Do analysis for calorimeters
  Bool_t  fRecalculateCaloPID ; // Recalculate PID or use PID weights in calorimeters
  Float_t fMinPt ;              // Maximum pt of (trigger) particles in the analysis
  Float_t fMaxPt ;              // Minimum pt of (trigger) particles in the analysis
  Int_t    fMultiBin ;	  // Number of bins in event container for multiplicity
  Int_t    fNZvertBin ;	  // Number of bins in event container for vertex position
  Int_t    fNrpBin ;	    // Number of bins in event container for reaction plain
  Float_t  fZvtxCut ;	    // Cut on vertex position  
  Int_t    fMaxMulti ;              // Maximum multiplicity of particles in the analysis
  Int_t    fMinMulti ;              // Maximum multiplicity of particles in the analysis
  Bool_t   fUseSelectEvent ; // Select events based on multiplicity and vertex cuts
  
	
  AliCaloTrackReader * fReader; // Acces to ESD/AOD/MC data
  
  TClonesArray* fInputAODBranch ;    //! Selected input particles branch
  TString       fInputAODName ;      //  Name of input AOD branch;
  TClonesArray* fOutputAODBranch ;   //! Selected output particles branch
  Bool_t        fNewAOD ;            //  Flag, new aod branch added to the analysis or not.
  TString       fOutputAODName ;     //  Name of output AOD branch;
  TString       fOutputAODClassName; //  Type of aod objects to be stored in the TClonesArray (AliAODPWG4Particle, AliAODPWG4ParticleCorrelation ...)	
  TString       fAODObjArrayName ;   // Name of ref array kept in a TList in AliAODParticleCorrelation with clusters or track references.
  TString       fAddToHistogramsName;// Add this string to histograms name
  
  //Analysis helper classes access pointers
  AliCaloPID               * fCaloPID; //! PID calculation
  AliFiducialCut           * fFidCut;  //! Acceptance cuts
  AliIsolationCut          * fIC;      //! Isolation cut 
  AliMCAnalysisUtils       * fMCUtils; //! MonteCarlo Analysis utils 
  AliNeutralMesonSelection * fNMS;     //! Neutral Meson Selection
  AliCalorimeterUtils      * fCaloUtils ; //  Pointer to CalorimeterUtils

  //Histograms binning and range    
  Int_t   fHistoPtBins   ;  // Number of bins in pt axis
  Float_t fHistoPtMax    ;  // Maximum value of pt histogram range
  Float_t fHistoPtMin    ;  // Minimum value of pt histogram range
  Int_t   fHistoPhiBins  ;  // Number of bins in phi axis
  Float_t fHistoPhiMax   ;  // Maximum value of phi histogram range
  Float_t fHistoPhiMin   ;  // Minimum value of phi histogram range
  Int_t   fHistoEtaBins  ;  // Number of bins in eta axis
  Float_t fHistoEtaMax   ;  // Maximum value of eta histogram range
  Float_t fHistoEtaMin   ;  // Minimum value of eta histogram range
  Int_t   fHistoMassBins ;  // Number of bins in mass axis
  Float_t fHistoMassMax  ;  // Maximum value of mass histogram range
  Float_t fHistoMassMin  ;  // Minimum value of mass histogram range
  Int_t   fHistoAsymBins ;  // Number of bins in asymmetry axis
  Float_t fHistoAsymMax  ;  // Maximum value of asymmetry histogram range
  Float_t fHistoAsymMin  ;  // Minimum value of asymmetry histogram range
  
  ClassDef(AliAnaPartCorrBaseClass,10)
} ;


#endif //ALIANAPARTCORRBASECLASS_H



