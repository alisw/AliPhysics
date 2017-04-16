#ifndef ALIANALYSISTASKCALOFILTER_H
#define ALIANALYSISTASKCALOFILTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliAnalysisTaskCaloFilter
/// \ingroup EMCALPerformance 
/// \brief Filter Calorimeter ESDs into AODs
///
/// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
/// PHOS or both, creating the corresponing AODCaloClusters
/// and AODCaloCells.
///
/// Also AODs are filtered. The main idea is to create lightweight AOD output with calorimeter
/// information mainly but also storing tracks so that full analysis like
/// correlations can be done.
///
/// Events can be filtered requiring hits in the EMCal or PHOS. Tracks are filtered
/// depending on track bits. MC particle info is also filtered.
/// Fill also the AODHeader information and the vertex.
/// Fill tracks if requested, option of only hybrid or all.
///
/// Option to store pure Minimum Bias events without event selection.
///
/// Copy of AliAnalysisTaskESDfilter.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

#include "AliAnalysisTaskSE.h"
class AliEMCALRecoUtils;
class AliEMCALGeometry;
class AliESDEvent;
class AliAODEvent;

class AliAnalysisTaskCaloFilter : public AliAnalysisTaskSE
{
 public:
  
  AliAnalysisTaskCaloFilter();
  AliAnalysisTaskCaloFilter(const char* name);
  virtual ~AliAnalysisTaskCaloFilter() ;
    
  // * General analysis frame methods *
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() { Init() ; }
  virtual void   UserExec(Option_t *option);
        
  // * Task own methods *
  
  Bool_t  AcceptEvent() ;
  
  Bool_t  AcceptEventEMCAL();

  Bool_t  AcceptEventPHOS();
  
  Bool_t  AcceptEventTrack();
  
  Bool_t  AcceptEventVertex();
  
  Bool_t  CheckForPrimaryVertex();
  Bool_t  CheckForPrimaryVertexInESDs();
  Bool_t  CheckForPrimaryVertexInAODs();
  
  void    CorrectionsInEMCAL();
  
  void    FillAODHeader();

  void    FillAODCaloCells();
  
  void    FillAODCaloClusters();
  
  void    FillAODCaloTrigger();

  void    FillAODMCParticles();
  
  void    FillAODTracks();

  void    FillAODv0s();
  
  void    FillAODVertices();
  
  void    FillAODVZERO();
  
  void    PrintInfo();
  
  // * Task settings *
  
  // Geometry methods
  
  void    SetEMCALGeometryName(TString name)                  { fEMCALGeoName = name        ; }
  TString    EMCALGeometryName()                       const  { return fEMCALGeoName        ; }
  
  void    SwitchOnLoadOwnEMCALGeometryMatrices()              { fLoadEMCALMatrices = kTRUE  ; }
  void    SwitchOffLoadOwnEMCALGeometryMatrices()             { fLoadEMCALMatrices = kFALSE ; }
  void    SetEMCALGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i]    = m      ; }
  
  //void    SwitchOnLoadOwnPHOSGeometryMatrices()               { fLoadPHOSMatrices = kTRUE  ; }
  //void    SwitchOffLoadOwnPHOSGeometryMatrices()              { fLoadPHOSMatrices = kFALSE ; }
  //void    SetPHOSGeometryMatrixInSM(TGeoHMatrix* m, Int_t i)  { fPHOSMatrix[i]    = m      ; }
    
  void    SwitchOnCheckEventPrimaryVertex()       { fCheckEventVertex = kTRUE   ; }
  void    SwitchOffCheckEventPrimaryVertex()      { fCheckEventVertex = kFALSE  ; }

  void    SwitchOnFillAODFile()                   { fFillAODFile = kTRUE        ; }
  void    SwitchOffFillAODFile()                  { fFillAODFile = kFALSE       ; }

  void    SwitchOnFillMCParticles()               { fFillMCParticles = kTRUE    ; }
  void    SwitchOffFillMCParticles()              { fFillMCParticles = kFALSE   ; }
  
  void    SwitchOnFillTracks()                    { fFillTracks  = kTRUE        ; }
  void    SwitchOffFillTracks()                   { fFillTracks  = kFALSE       ; }
  
  void    SwitchOnFillHybridTracks()              { fFillTracks  = kTRUE        ;
                                                    fFillHybridTracks  = kTRUE  ; }
  void    SwitchOffFillHybridTracks()             { fFillHybridTracks  = kFALSE ; }
  
  void    SwitchOnFillv0s()                       { fFillv0s     = kTRUE        ; }
  void    SwitchOffFillv0s()                      { fFillv0s     = kFALSE       ; }
  
  void    SwitchOnFillVZERO()                     { fFillVZERO   = kTRUE        ; }
  void    SwitchOffFillVZERO()                    { fFillVZERO   = kFALSE       ; }
  
  void    SwitchOnFillAllVertices()               { fFillAllVertices = kTRUE    ; }
  void    SwitchOffFillAllVertices()              { fFillAllVertices = kFALSE   ; }
  
  /// Enumerate the options to filter:
  /// * kBoth: Both EMCal and PHOS,
  /// * kEMCAL
  /// * kPHOS
  enum    caloFilter { kBoth = 0, kEMCAL = 1, kPHOS = 2 } ;
  
  void    SetCaloFilter(Int_t calo)               { fCaloFilter = calo          ; }
  TString GetCaloFilter()                  const  { return fCaloFilter          ; }  

  /// Select which detector involved in event selection: PHOS, EMCAL or tracking
  /// you can filter EMCal or PHOS clusters, but the event selection is independent.
  void    SetEventSelection(Bool_t emcal, Bool_t phos, Bool_t track) {
    fEventSelection[0] = emcal ; fEventSelection[1] = phos ; fEventSelection[2] = track ; }
  
  void    SwitchOnAcceptAllMBEvent()              { fAcceptAllMBEvent = kTRUE   ; }
  void    SwitchOffAcceptAllMBEvent()             { fAcceptAllMBEvent = kFALSE  ; }

  void    SetMBTriggerMask(UInt_t mask)           { fMBTriggerMask    = mask    ; }
  
  void    SetEMCALRecoUtils(AliEMCALRecoUtils* ru){ fEMCALRecoUtils = ru        ; }
  AliEMCALRecoUtils* GetEMCALRecoUtils()   const  { return fEMCALRecoUtils      ; }

  void    SwitchOnClusterCorrection()             { fCorrect = kTRUE            ; }
  void    SwitchOffClusterCorrection()            { fCorrect = kFALSE           ; }
  
  void    SetConfigFileName(TString name)         { fConfigName = name          ; }
  
  void    SetEMCALEnergyCut(Float_t cut)          { fEMCALEnergyCut = cut       ; }
  Float_t GetEMCALEnergyCut()               const { return fEMCALEnergyCut      ; }
  void    SetEMCALNcellsCut(Int_t cut)            { fEMCALNcellsCut = cut       ; }
  Int_t   GetEMCALNcellsCut()               const { return fEMCALNcellsCut      ; }
  
  void    SetPHOSEnergyCut(Float_t cut)           { fPHOSEnergyCut = cut        ; }
  Float_t GetPHOSEnergyCut()                const { return fPHOSEnergyCut       ; }
  void    SetPHOSNcellsCut(Int_t cut)             { fPHOSNcellsCut = cut        ; }
  Int_t   GetPHOSNcellsCut()                const { return fPHOSNcellsCut       ; }
  
  void    SetTrackPtCut(Float_t cut)              { fTrackPtCut = cut           ; }
  Float_t GetTrackPtCut()                   const { return fTrackPtCut          ; }
  
  void    SetVzCut(Float_t cut)                   { fVzCut = cut                ; }
  Float_t GetVzCut()                        const { return fVzCut               ; }
  
  
private:
    
  Int_t               fCaloFilter;        ///<  Calorimeter to filter: kBoth, kEMCAL, kPHOS.
  Bool_t              fEventSelection[3]; ///<  Define which detector is used to select the event: {EMCAL,PHOS,Tracks}.
  Bool_t              fAcceptAllMBEvent;  ///<  Do not select the MB events with same Event selection cuts as other triggers.
  UInt_t              fMBTriggerMask;     ///<  Define the mask for MB events, it should be kMB, but not always defined, use kAnyINT instead.
  Int_t               fCorrect;           ///<  Recalibrate or recalculate different cluster parameters, only for EMCal.
  
  //EMCAL specific
  AliEMCALGeometry  * fEMCALGeo;          //!<! EMCAL geometry.
  TString             fEMCALGeoName;      ///<  Name of geometry to use.
  AliEMCALRecoUtils * fEMCALRecoUtils;    ///<  Pointer to EMCAL utilities for clusterization.
  
  //Geometry
  Bool_t              fLoadEMCALMatrices; ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs.
  TGeoHMatrix       * fEMCALMatrix[22];   ///<  Geometry matrices with alignments.
  //Bool_t            fLoadPHOSMatrices;  ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs.
  //TGeoHMatrix *     fPHOSMatrix[5];     ///<  Geometry matrices with alignments.
  Bool_t              fGeoMatrixSet;      ///<  Set geometry matrices only once, for the first event.
  
  TString             fConfigName;        ///<  Name of analysis configuration file.
  
  Bool_t              fFillAODFile;       ///<  Fill the output AOD file with clusters.
  Bool_t              fFillMCParticles;   ///<  Fill the output AOD file with MC particles.
  Bool_t              fFillTracks;        ///<  Fill the output AOD file with tracks.
  Bool_t              fFillHybridTracks;  ///<  Fill the output AOD file with hybrid tracks, only when fFillTracks = kTRUE.

  Bool_t              fFillAllVertices;   ///<  Fill the output AOD file with all vertices.
  Bool_t              fFillv0s;           ///<  Fill the output AOD file with v0s.
  Bool_t              fFillVZERO;         ///<  Fill the output AOD file with VZERO input.

  Float_t             fEMCALEnergyCut;    ///<  At least an EMCAL cluster with this energy in the event.
  Int_t               fEMCALNcellsCut;    ///<  At least an EMCAL cluster with fNCellsCut cells over fEnergyCut.

  Float_t             fPHOSEnergyCut;     ///<  At least a PHOS cluster with this energy in the event.
  Int_t               fPHOSNcellsCut;     ///<  At least a PHOS cluster with fNCellsCut cells over fEnergyCut.
  
  Float_t             fTrackPtCut;        ///<  Select events with at least a track with this pT.
  
  Float_t             fVzCut;             ///<  At least events with vertex within cut.
  
  Bool_t              fCheckEventVertex;  ///<  Check the primary vertex of the event or not.

  AliVEvent*          fEvent;             //!<! Event pointer.
  AliESDEvent*        fESDEvent;          //!<! ESD event pointer (cast of fEvent).
  AliAODEvent*        fAODEvent;          //!<! AOD event pointer (cast of fEvent).
  
  /// Copy constructor not implemented.
  AliAnalysisTaskCaloFilter(           const AliAnalysisTaskCaloFilter&) ;
    
  /// Assignment operator not implemented.
  AliAnalysisTaskCaloFilter& operator=(const AliAnalysisTaskCaloFilter&) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCaloFilter, 12) ;
  /// \endcond

};

#endif
