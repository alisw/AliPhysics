#ifndef ALIANAPARTCORRBASECLASS_H
#define ALIANAPARTCORRBASECLASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Base class for analysis algorithms
//-- Author: Gustavo Conesa (INFN-LNF)

#include <cstdlib>

//ROOT
class TClonesArray ;
class TRefArray ;
#include <TList.h>
#include <TObject.h>

//Analysis
class AliAODCaloCluster;
class AliAODCaloCells;
#include "AliAODPWG4Particle.h"
class AliCaloTrackReader ;   
class AliCaloPID ;
class AliFidutialCut ;
class AliIsolationCut ;
class AliMCAnalysisUtils ;
class AliNeutralMesonSelection ;
class AliStack ; 
class AliHeader ; 
class AliGenEventHeader ; 


class AliAnaPartCorrBaseClass : public TObject {
	
public: 
	
  AliAnaPartCorrBaseClass() ; // default ctor
  AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & g) ; // cpy ctor
  AliAnaPartCorrBaseClass & operator = (const AliAnaPartCorrBaseClass & g) ;//cpy assignment
  virtual ~AliAnaPartCorrBaseClass() ; //virtual dtor
  
//	virtual void AddAODCaloCluster(AliAODCaloCluster calo) ;
  virtual void AddAODParticle(AliAODPWG4Particle pc) ;
//
//	virtual void ConnectAODCaloClusters();
  virtual void ConnectAODPHOSCells();
  virtual void ConnectAODEMCALCells();
  virtual void ConnectInputOutputAODBranches();
  
  virtual TList * GetCreateOutputObjects() { return (new TList) ;}
  
  virtual void Init() {;}
  virtual void InitParameters() ;
  
  virtual void Print(const Option_t * ) const ;
  
  virtual void MakeAnalysisFillAOD()  {;}
  
  virtual void MakeAnalysisFillHistograms() {;}
  
  virtual Int_t GetDebug() const  { return fDebug ; }
  virtual void SetDebug(Int_t d)   { fDebug = d ; }
  
  virtual Int_t GetEventNumber() const ;
  
  virtual AliCaloTrackReader * GetReader() const {return fReader ; }
  virtual void SetReader(AliCaloTrackReader * reader) { fReader = reader ; }
  
  virtual void Terminate() {;}
  
  //analysis AOD branch
  virtual TClonesArray * GetCreateOutputAODBranch() ;
  virtual TString GetInputAODName() const {return fInputAODName ; }
  virtual void SetInputAODName(TString name)   { fInputAODName = name; }	
  virtual TString GetOutputAODName()  const {return fOutputAODName ; }
  virtual void SetOutputAODName(TString name)   { fNewAOD = kTRUE ; fOutputAODName = name; }
  virtual Bool_t NewOutputAOD() const {return fNewAOD;}
  virtual TString GetOutputAODClassName() const {return fOutputAODClassName;}
  virtual void SetOutputAODClassName(TString name) {fOutputAODClassName = name; }
  
  virtual TClonesArray* GetInputAODBranch() const {return fInputAODBranch ;}
  virtual TClonesArray* GetOutputAODBranch() const {return fOutputAODBranch ;}

//	virtual TClonesArray* GetAODCaloClusters() const {return fAODCaloClusters ;}
  virtual TClonesArray* GetAODCaloClusters() const ;
  virtual TClonesArray* GetAODTracks() const ;	
  virtual AliAODCaloCells* GetAODCaloCells() const {return fAODCaloCells ;}
  
  virtual TRefArray* GetAODCTS() const ;
  virtual TRefArray* GetAODEMCAL() const ;
  virtual TRefArray* GetAODPHOS() const ;
  
  virtual TString	GetBaseParametersList();
  
  virtual TNamed * GetEMCALCells() const ;
  virtual TNamed * GetPHOSCells() const ;
  
  virtual AliStack * GetMCStack() const ;
  virtual AliHeader* GetMCHeader() const ;
  virtual AliGenEventHeader* GetMCGenEventHeader() const ;
  
  //Analysis helpers classes pointers setters and getters
  virtual AliCaloPID * GetCaloPID() const {return  fCaloPID ;}
  virtual void SetCaloPID(AliCaloPID * pid) { fCaloPID = pid ;}
  
  virtual AliFidutialCut * GetFidutialCut() const {return  fFidCut ;}
  virtual void SetFidutialCut(AliFidutialCut * fc) { fFidCut = fc ;}
  
  virtual AliIsolationCut * GetIsolationCut() const {return  fIC ;}
  virtual void SetIsolationCut(AliIsolationCut * fc) { fIC = fc ;}
  
  virtual AliMCAnalysisUtils * GetMCAnalysisUtils() const {return  fMCUtils ;}
  virtual void SetMCAnalysisUtils(AliMCAnalysisUtils * mcutils) { fMCUtils = mcutils ;}	
  
  virtual AliNeutralMesonSelection * GetNeutralMesonSelection() const {return  fNMS ;}
  virtual void SetNeutralMesonSelection(AliNeutralMesonSelection * nms) { fNMS = nms ;}
  
  virtual Bool_t     IsDataMC() const {return fDataMC ; }
  virtual void SwitchOnDataMC()    {fDataMC = kTRUE ; }
  virtual void SwitchOffDataMC()    {fDataMC = kFALSE ; }
  
  virtual Bool_t IsFidutialCutOn() const {return fCheckFidCut ; }
  virtual void SwitchOnFidutialCut() { fCheckFidCut = kTRUE;}
  virtual void SwitchOffFidutialCut() { fCheckFidCut = kFALSE;}
  
  virtual Bool_t IsCaloPIDOn() const {return fCheckCaloPID ; }
  virtual void SwitchOnCaloPID() { fCheckCaloPID = kTRUE;}
  virtual void SwitchOffCaloPID() { fCheckCaloPID = kFALSE;}
  
  virtual Bool_t IsCaloPIDRecalculationOn() const {return fRecalculateCaloPID ; }
  virtual void SwitchOnCaloPIDRecalculation() { fRecalculateCaloPID  = kTRUE;}
  virtual void SwitchOffCaloPIDRecalculation() { fRecalculateCaloPID  = kFALSE;}
  
  virtual Float_t    GetMaxPt()         const {return fMaxPt ; }
  virtual Float_t    GetMinPt()         const {return fMinPt ; }
  virtual void SetMaxPt(Float_t pt)              {fMaxPt = pt ; }
  virtual void SetMinPt(Float_t pt)              {fMinPt = pt ; }
  void SetPtCutRange(Double_t ptmin, Double_t ptmax)
  {  fMaxPt=ptmax;   fMinPt=ptmin;}
  
  //Histogrammes setters and getters
  virtual void SetHistoPtRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtBins = n ;
    fHistoPtMax = max ;
    fHistoPtMin = min ;
  }
  
  Int_t   GetHistoNPtBins() const { return fHistoNPtBins ; }
  Float_t GetHistoPtMin()   const { return fHistoPtMin ; }
  Float_t GetHistoPtMax()   const { return fHistoPtMax ; }
  
  virtual void SetHistoPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPhiBins = n ;
    fHistoPhiMax = max ;
    fHistoPhiMin = min ;
  }
  
  Int_t   GetHistoNPhiBins() const { return fHistoNPhiBins ; }
  Float_t GetHistoPhiMin()   const { return fHistoPhiMin ; }
  Float_t GetHistoPhiMax()   const { return fHistoPhiMax ; }
  
  virtual void SetHistoEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNEtaBins = n ;
    fHistoEtaMax = max ;
    fHistoEtaMin = min ;
  }
  
  Int_t   GetHistoNEtaBins() const { return fHistoNEtaBins ; }
  Float_t GetHistoEtaMin()   const { return fHistoEtaMin ; }
  Float_t GetHistoEtaMax()   const { return fHistoEtaMax ; }
  
  
 private:    
  
  Bool_t  fDataMC ;             // Flag to access MC data when using ESD or AOD     
  Int_t   fDebug ;              // Debug level
  Bool_t  fCheckFidCut ;        // Do analysis for clusters in defined region         
  Bool_t  fCheckCaloPID ;       // Do analysis for calorimeters
  Bool_t  fRecalculateCaloPID ; // Recalculate PID or use PID weights in calorimeters
  Float_t fMinPt ;              // Maximum pt of (trigger) particles in the analysis
  Float_t fMaxPt ;              // Minimum pt of (trigger) particles in the analysis
  
  AliCaloTrackReader * fReader; // Acces to ESD/AOD/MC data
  
  TClonesArray* fInputAODBranch ;    //! Selected input particles branch
  TString       fInputAODName ;      //  Name of input AOD branch;
  TClonesArray* fOutputAODBranch ;   //! Selected output particles branch
  Bool_t        fNewAOD ;            //  Flag, new aod branch added to the analysis or not.
  TString       fOutputAODName ;     //  Name of output AOD branch;
  TString       fOutputAODClassName; //  Type of aod objects to be stored in the TClonesArray (AliAODPWG4Particle, AliAODPWG4ParticleCorrelation ...)	
  
  //TClonesArray* fAODCaloClusters ;  //! selected PHOS/EMCAL CaloClusters
  AliAODCaloCells * fAODCaloCells ; //! selected PHOS/EMCAL CaloCells
  
  //Analysis helper classes access pointers
  AliCaloPID               * fCaloPID; // PID calculation
  AliFidutialCut           * fFidCut;  // Acceptance cuts
  AliIsolationCut          * fIC;      // Isolation cut 
  AliMCAnalysisUtils       * fMCUtils; // MonteCarlo Analysis utils 
  AliNeutralMesonSelection * fNMS;     // Neutral Meson Selection
  
  //Histograms binning and range    
  Int_t   fHistoNPtBins ;  // Number of bins in pt axis
  Float_t fHistoPtMax ;    // Maximum value of pt histogram range
  Float_t fHistoPtMin ;    // Minimum value of pt histogram range
  Int_t   fHistoNPhiBins ; // Number of bins in phi axis
  Float_t fHistoPhiMax ;   // Maximum value of phi histogram range
  Float_t fHistoPhiMin ;   // Minimum value of phi histogram range
  Int_t   fHistoNEtaBins ; // Number of bins in eta axis
  Float_t fHistoEtaMax ;   // Maximum value of eta histogram range
  Float_t fHistoEtaMin ;   // Minimum value of eta histogram range
  
  ClassDef(AliAnaPartCorrBaseClass,3)
    } ;


#endif //ALIANAPARTCORRBASECLASS_H



