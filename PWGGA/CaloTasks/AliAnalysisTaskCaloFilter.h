#ifndef ALIANALYSISTASKCALOFILTER_H
#define ALIANALYSISTASKCALOFILTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Fill also the AODHeader information and the vertex.
// Fill tracks if requested
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

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
    
  //General analysis frame methods
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() { Init() ; }
  virtual void   UserExec(Option_t *option);
        
  // Task own methods
  
  Bool_t  AcceptEvent() ;
  
  Bool_t  AcceptEventEMCAL();

  Bool_t  AcceptEventPHOS();
  
  Bool_t  AcceptEventTrack();
  
  Bool_t  AcceptEventVertex();
  
  Bool_t  CheckForPrimaryVertex();
  
  void    CorrectionsInEMCAL();
  
  void    FillAODHeader();

  void    FillAODCaloCells();
  
  void    FillAODCaloClusters();
  
  void    FillAODCaloTrigger();

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
    
  void    SwitchOnFillAODFile()                   { fFillAODFile = kTRUE        ; }
  void    SwitchOffFillAODFile()                  { fFillAODFile = kFALSE       ; }

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
  
  enum    caloFilter {kBoth = 0, kEMCAL = 1, kPHOS=2};
  
  void    SetCaloFilter(Int_t calo)               { fCaloFilter = calo          ; }
  TString GetCaloFilter()                  const  { return fCaloFilter          ; }  
  
  void    SetEventSelection(Bool_t emcal, Bool_t phos, Bool_t track) { 
    // Detector involved in event selection
    fEventSelection[0] = emcal ; fEventSelection[1] = phos ; fEventSelection[2] = track ; }
  
  void    SwitchOnAcceptAllMBEvent()              { fAcceptAllMBEvent = kTRUE   ; }
  void    SwitchOffAcceptAllMBEvent()             { fAcceptAllMBEvent = kFALSE  ; }

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
    
  Int_t               fCaloFilter;        // Calorimeter to filter
  Bool_t              fEventSelection[3]; // Define which detector is used to select the event
  Bool_t              fAcceptAllMBEvent;    // Do not select the MB events with same cuts as other triggers
  Int_t               fCorrect;           // Recalibrate or recalculate different cluster parameters
  
  //EMCAL specific
  AliEMCALGeometry  * fEMCALGeo;          //! EMCAL geometry
  TString             fEMCALGeoName;      // Name of geometry to use.
  AliEMCALRecoUtils * fEMCALRecoUtils;    // Pointer to EMCAL utilities for clusterization
  
  //Geometry
  Bool_t              fLoadEMCALMatrices; // Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  TGeoHMatrix       * fEMCALMatrix[12];   // Geometry matrices with alignments
  //Bool_t            fLoadPHOSMatrices;  // Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  //TGeoHMatrix *     fPHOSMatrix[5];     // Geometry matrices with alignments
  Bool_t              fGeoMatrixSet;      // Set geometry matrices only once, for the first event.   
  
  TString             fConfigName;        // Name of analysis configuration file
  
  Bool_t              fFillAODFile;       // Fill the output AOD file with clusters 
  Bool_t              fFillTracks;        // Fill tracks
  Bool_t              fFillHybridTracks;  // Fill hybrid tracks

  Bool_t              fFillAllVertices;   // Fill all vertices
  Bool_t              fFillv0s;           // Fill v0s
  Bool_t              fFillVZERO;         // Fill VZERO

  Float_t             fEMCALEnergyCut;    //  At least an EMCAL cluster with this energy in the event
  Int_t               fEMCALNcellsCut;    //  At least an EMCAL cluster with fNCellsCut cells over fEnergyCut

  Float_t             fPHOSEnergyCut;     //  At least a PHOS cluster with this energy in the event
  Int_t               fPHOSNcellsCut;     //  At least a PHOS cluster with fNCellsCut cells over fEnergyCut
  
  Float_t             fTrackPtCut;        //  At least a track with this pT in the event
  
  Float_t             fVzCut;             //  At least events with vertex within cut
  
  AliVEvent*          fEvent;             //! event pointer
  AliESDEvent*        fESDEvent;          //! ESD event pointer
  AliAODEvent*        fAODEvent;          //! AOD event pointer

  
  AliAnalysisTaskCaloFilter(           const AliAnalysisTaskCaloFilter&);
  AliAnalysisTaskCaloFilter& operator=(const AliAnalysisTaskCaloFilter&);
  
  ClassDef(AliAnalysisTaskCaloFilter, 9); // Analysis task for standard ESD filtering
  
};

#endif
