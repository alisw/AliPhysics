#ifndef ALIANALYSISTASKCALOFILTER_H
#define ALIANALYSISTASKCALOFILTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Keep tracks if requested
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
class AliEMCALRecoUtils;
class AliEMCALGeometry;

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
        
  Bool_t  AcceptEventEMCAL();
  
  Bool_t  AcceptEventVertex();
  
  Bool_t  CheckForPrimaryVertex();
  
  void    CorrectionsInEMCAL();
  
  void    FillAODHeader();

  void    FillAODCaloCells();
  
  void    FillAODCaloClusters();
  
  void    FillAODTracks();
  
  void    FillAODVertices();
  
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
    
  void    SwitchOnFillAODFile()                   { fFillAODFile = kTRUE    ; }
  void    SwitchOffFillAODFile()                  { fFillAODFile = kFALSE   ; }

  void    SwitchOnFillTracks()                    { fFillTracks  = kTRUE    ; }
  void    SwitchOffFillTracks()                   { fFillTracks  = kFALSE   ; }
  
  enum    caloFilter {kBoth = 0, kEMCAL = 1, kPHOS=2};
  void    SetCaloFilter(Int_t calo)               { fCaloFilter = calo      ; }
  TString GetCaloFilter()                  const  { return fCaloFilter      ; }  
  
  void    SetEMCALRecoUtils(AliEMCALRecoUtils* ru){ fEMCALRecoUtils = ru    ; }
  AliEMCALRecoUtils* GetEMCALRecoUtils()   const  { return fEMCALRecoUtils  ; }

  void    SwitchOnClusterCorrection()             { fCorrect = kTRUE        ; }
  void    SwitchOffClusterCorrection()            { fCorrect = kFALSE       ; }
  
  void    SetConfigFileName(TString name)         { fConfigName = name      ; }
  
  void    SetEnergyCut(Float_t cut)               { fEnergyCut = cut        ; }
  Float_t GetEnergyCut()                    const { return fEnergyCut       ; }
  void    SetNcellsCut(Int_t cut)                 { fNcellsCut = cut        ; }
  Int_t   GetNcellsCut()                    const { return fNcellsCut       ; }
  void    SetVzCut(Float_t cut)                   { fVzCut = cut            ; }
  Int_t   GetVzCut()                        const { return fVzCut           ; }
  
  
private:
    
  Int_t               fCaloFilter;        // Calorimeter to filter
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
  
  Float_t             fEnergyCut;         //  At least a cluster with this energy in the event
  Int_t               fNcellsCut;         //  At least a cluster with fNCellsCut cells over fEnergyCut
  Float_t             fVzCut;             //  At least events with vertex within cut
  
  AliAnalysisTaskCaloFilter(           const AliAnalysisTaskCaloFilter&);
  AliAnalysisTaskCaloFilter& operator=(const AliAnalysisTaskCaloFilter&);
  
  ClassDef(AliAnalysisTaskCaloFilter, 8); // Analysis task for standard ESD filtering
  
};

#endif
