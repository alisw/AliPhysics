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
  virtual void   UserCreateOutputObjects();
  //virtual void   Init();
  //virtual void   LocalInit() ;
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  
  enum caloFilter {kBoth = 0, kEMCAL = 1, kPHOS=2};
  void SetCaloFilter(Int_t calo) {fCaloFilter = calo;}
  TString GetCaloFilter() const  {return fCaloFilter;}

  void SetGeometryName(TString name) { fEMCALGeoName = name ; }
  TString GeometryName() const       { return fEMCALGeoName ; }
  
  void SetEMCALRecoUtils(AliEMCALRecoUtils * ru) {fEMCALRecoUtils = ru;}
  AliEMCALRecoUtils* GetEMCALRecoUtils() const   {return fEMCALRecoUtils;}

  void SwitchOnClusterCorrection()  {fCorrect = kTRUE ;}
  void SwitchOffClusterCorrection() {fCorrect = kFALSE;}

  AliESDtrackCuts* GetTrackCuts()          const  { return fESDtrackCuts    ; }
  void    SetTrackCuts(AliESDtrackCuts * cuts)    { fESDtrackCuts = cuts    ; }		  
  Float_t GetTrackMultiplicityEtaCut()     const  { return fTrackMultEtaCut ; }
  void    SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut = eta  ; }		
  
  
  void PrintInfo();
  
private:
  
  //TList* fCuts ;      //! List with analysis cuts
  Int_t  fCaloFilter; // Calorimeter to filter
  Int_t  fCorrect;    // Recalibrate or recalculate different cluster parameters
  //EMCAL specific
  AliEMCALGeometry  * fEMCALGeo;       //! EMCAL geometry
  TString             fEMCALGeoName;   // Name of geometry to use.
  AliEMCALRecoUtils * fEMCALRecoUtils; // Pointer to EMCAL utilities for clusterization

  AliESDtrackCuts *fESDtrackCuts       ; // Track cut  
  Float_t          fTrackMultEtaCut    ; // Track multiplicity eta cut
  
  
  ClassDef(AliAnalysisTaskCaloFilter, 3); // Analysis task for standard ESD filtering
};

#endif
