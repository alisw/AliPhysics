#ifndef ALIANALYSISTASKCALOFILTER_H
#define ALIANALYSISTASKCALOFILTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskCaloFilter.h  $ */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Needed for calorimeter calibration.
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

class TList;

#include "AliAnalysisTaskSE.h"
class AliEMCALRecoUtils;
class AliEMCALGeometry;
class AliESDtrackCuts;
class AliTriggerAnalysis;
class TNtuple;
class AliAnalysisTaskCaloFilter : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCaloFilter();
  AliAnalysisTaskCaloFilter(const char* name);
  virtual ~AliAnalysisTaskCaloFilter() ;
  
private:
  AliAnalysisTaskCaloFilter(const AliAnalysisTaskCaloFilter&);
  AliAnalysisTaskCaloFilter& operator=(const AliAnalysisTaskCaloFilter&);
  
public:
  //General analysis frame methods
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() { Init() ; }
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
      
  //Geometry methods
  void SetEMCALGeometryName(TString name)                  { fEMCALGeoName = name        ; }
  TString EMCALGeometryName()                       const  { return fEMCALGeoName        ; }
  void SwitchOnLoadOwnEMCALGeometryMatrices()              { fLoadEMCALMatrices = kTRUE  ; }
  void SwitchOffLoadOwnEMCALGeometryMatrices()             { fLoadEMCALMatrices = kFALSE ; }
  void SetEMCALGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i]    = m      ; }
  //void SwitchOnLoadOwnPHOSGeometryMatrices()               { fLoadPHOSMatrices = kTRUE  ; }
  //void SwitchOffLoadOwnPHOSGeometryMatrices()              { fLoadPHOSMatrices = kFALSE ; }
  //void SetPHOSGeometryMatrixInSM(TGeoHMatrix* m, Int_t i)  { fPHOSMatrix[i]    = m      ; }
  
  //Task settings
  void    FillAODFile(Bool_t yesno)               { fFillAODFile = yesno    ; }
  enum    caloFilter {kBoth = 0, kEMCAL = 1, kPHOS=2};
  void    SetCaloFilter(Int_t calo)               { fCaloFilter = calo      ; }
  TString GetCaloFilter()                  const  { return fCaloFilter      ; }  
  
  void SetEMCALRecoUtils(AliEMCALRecoUtils * ru)  { fEMCALRecoUtils = ru    ; }
  AliEMCALRecoUtils* GetEMCALRecoUtils()   const  { return fEMCALRecoUtils  ; }

  void    SwitchOnClusterCorrection()             { fCorrect = kTRUE        ; }
  void    SwitchOffClusterCorrection()            { fCorrect = kFALSE       ; }
  
  //Event selection
  AliESDtrackCuts* GetTrackCuts()          const  { return fESDtrackCuts    ; }
  void    SetTrackCuts(AliESDtrackCuts * cuts)    { fESDtrackCuts = cuts    ; }		  
  Float_t GetTrackMultiplicityEtaCut()     const  { return fTrackMultEtaCut ; }
  void    SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut = eta  ; }		
  virtual Bool_t CheckForPrimaryVertex();

  void    PrintInfo();
  
  void    SetConfigFileName(TString name)         { fConfigName = name      ; }
  
private:
  
  //TList* fCuts ;      //! List with analysis cuts
  Int_t               fCaloFilter;        // Calorimeter to filter
  Int_t               fCorrect;           // Recalibrate or recalculate different cluster parameters
  //EMCAL specific
  AliEMCALGeometry  * fEMCALGeo;          //! EMCAL geometry
  TString             fEMCALGeoName;      // Name of geometry to use.
  AliEMCALRecoUtils * fEMCALRecoUtils;    // Pointer to EMCAL utilities for clusterization

  AliESDtrackCuts   * fESDtrackCuts;      // Track cut  
  AliTriggerAnalysis* fTriggerAnalysis;   // Access to trigger selection algorithm for V0AND calculation
  Float_t             fTrackMultEtaCut;   // Track multiplicity eta cut
  
  //Geometry
  Bool_t              fLoadEMCALMatrices; // Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  TGeoHMatrix       * fEMCALMatrix[10];   // Geometry matrices with alignments
  //Bool_t            fLoadPHOSMatrices;  // Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  //TGeoHMatrix *     fPHOSMatrix[5];     // Geometry matrices with alignments
  Bool_t              fGeoMatrixSet;      // Set geometry matrices only once, for the first event.   

  TNtuple           * fEventNtuple;       // NTuple with event parameters 
  
  TString             fConfigName;        // Name of analysis configuration file
  Bool_t              fFillAODFile;       // Fill the output AOD file with clusters 
  
  ClassDef(AliAnalysisTaskCaloFilter, 6); // Analysis task for standard ESD filtering
};

#endif
