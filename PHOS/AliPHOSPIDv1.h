#ifndef ALIPHOSPIDV1_H
#define ALIPHOSPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
// Identification is based on information from CPV and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH), Gustavo Conesa.

// --- ROOT system ---
//class TFormula ;
class TVector3 ;
class TMatrixD ;
class TPrincipal ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;

#include "AliPHOSPID.h"

class  AliPHOSPIDv1 : public AliPHOSPID {
  
 public:
  
  AliPHOSPIDv1() ;          // ctor            
  AliPHOSPIDv1(const char* headerFile, const char * tsBranch = "Default", const Bool_t toSplit=kFALSE) ;
  
  virtual ~AliPHOSPIDv1() ; // dtor
  
  virtual void Exec(Option_t * option) ;
  //  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  //  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual const Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
  
  virtual void Print(Option_t * option) const {}
  void Print() ; 
  
  //Get files that contain the PCA
  const TString GetPrincipalFile( )const {return fFileName ;}
  const TString GetPrincipalFilePar( )const {return fFileNamePar ;}
  
  // Set and Get all parameters necessary in the PID depending on the 
  // custer energy and Purity-Efficiency point (possible options "HIGH 
  // EFFICIENCY" "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options 
  // changing EFFICIENCY by PURITY)
  void SetCpvtoEmcDistanceCut(Float_t Cluster_En, TString Eff_Pur, Float_t cut)  ; 
  void SetTimeGate(Float_t Cluster_En, TString Eff_Pur, Float_t gate)  ; 
 
  const Double_t GetCpvtoEmcDistanceCut(const Float_t Cluster_En, const TString Eff_Pur) const ;
  const Double_t GetTimeGate(const Float_t Cluster_En, const TString Eff_Pur)  const;
 
  void SetEllipseParameter(TString Param, Int_t i, Double_t par)  ;    
  const Double_t GetParameterToCalculateEllipse(const TString Param, const Int_t i) const  ;     
  const Double_t GetEllipseParameter(const TString Param, Float_t E) const;
  //Get and Set energy calibration parameters
  
  void  SetCalibrationParameter(Int_t Param,Double_t param);
  const Double_t GetCalibrationParameter(const Int_t i) const;
  const Double_t GetCalibratedEnergy(const Float_t e)   const; //Calibrates energy.
  
  //  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  //  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v1" ; }  
  
 private:
  
  const TString BranchName() const ; 
  virtual void Init() ;
  virtual void InitParameters() ;
  void     MakeRecParticles(void ) ;
  // Relative Distance CPV-EMC
  const Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; 
  const Int_t    GetPrincipalSign(const Double_t* P, const Int_t eff_pur, const Float_t E)const ; //Principal cut
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;
  // Gives in wich cluster energy range is the event
  const Int_t     GetClusterOption(const Float_t Cluster_En) const;
  // Gives the Efficiency-Purity point.
  const Int_t    GetEffPurOption(const TString Eff_Pur)const ;
  virtual  void WriteRecParticles(Int_t event) ; 
  void     SetParameters() ; //Fills the matrix of parameters
  
 

 private:

  Bool_t                 fDefaultInit;        //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString    fFileName ;  // File that contains the Principal file for analysis
  TString    fFileNamePar ;// File that contains the parameters for analysis
  
  //  TString    fFrom ;              // name of Recpoints and TrackSegments 
  //  TString    fHeaderFileName ;    // file name with event header
  //  TString    fTrackSegmentsTitle; // branch name with track segments
  //  TString    fRecPointsTitle ;    // branch name with rec points
  //  TString    fRecParticlesTitle ; // branch name with rec particles
 
  Int_t                      fNEvent ;            //! current event number
  //  AliPHOSClusterizer *       fClusterizer ;       //! clusterizer
  //  AliPHOSTrackSegmentMaker * fTSMaker ;           //! track segment maker

  TPrincipal *               fPrincipal ;         //! TPrincipal from pca file 
  Double_t *                 fX ;                 //! Principal data 
  Double_t *                 fP ;                 //! Principal eigenvalues
  Int_t                      fRecParticlesInRun ; //! Total number of recparticles in one run
  TMatrixD *                 fParameters;         //! Matrix of identification Parameters

  ClassDef( AliPHOSPIDv1,6)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
