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
  AliPHOSPIDv1(const char* headerFile, const char * tsBranch = "Default", const char * from = 0) ;
   
  virtual ~AliPHOSPIDv1() ; // dtor

  virtual void Exec(Option_t * option) ;
  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual const Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
 
  virtual void Print(Option_t * option) const {}
  void Print() ; 
  // Get CpvtoEmcDistanceCut and TimeGate parameters depending on the custer energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" "MEDIUM EFFICIENCY" "LOW  
  // EFFICIENCY" and 3 more options changing EFFICIENCY by PURITY)
  Double_t GetCpvtoEmcDistanceCut(const Float_t Cluster_En, const TString Eff_Pur)  ;
  Double_t GetTimeGate(const Float_t Cluster_En, const TString Eff_Pur)  ;

  //Get files that contain the PCA
  const TString GetPrincipalFile5( )const {return fFileName5 ;}
  const TString GetPrincipalFilePar5( )const {return fFileNamePar5 ;}
  const TString GetPrincipalFile100( )const {return fFileName100 ;}
  const TString GetPrincipalFilePar100( )const {return fFileNamePar100 ;}

  // Set all parameters necessary in the PID depending on the custer energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" "MEDIUM EFFICIENCY" "LOW  
  // EFFICIENCY" and 3 more options changing EFFICIENCY by PURITY)
  void SetCpvtoEmcDistanceCut(Float_t Cluster_En, TString Eff_Pur, Float_t cut)  ; 
  void SetTimeGate(Float_t Cluster_En, TString Eff_Pur, Float_t gate)  ; 
  void SetEllipseXCenter(Float_t Cluster_En, TString Eff_Pur, Float_t x)  ;    
  void SetEllipseYCenter(Float_t Cluster_En, TString Eff_Pur, Float_t y)  ;   
  void SetEllipseAParameter(Float_t Cluster_En, TString Eff_Pur, Float_t a)  ;  
  void SetEllipseBParameter(Float_t Cluster_En, TString Eff_Pur, Float_t b)  ; 
  void SetEllipseAngle(Float_t Cluster_En, TString Eff_Pur, Float_t angle)  ;    
  void SetEllipseParameters(Float_t Cluster_En, TString Eff_Pur, Float_t x, Float_t y,Float_t a, Float_t b,Float_t angle) ;  
 
  //Get and Set energy calibration parameters
  Float_t GetACalParameter() {return fACalParameter ;}
  Float_t GetBCalParameter() {return fBCalParameter ;}
  Float_t GetCCalParameter() {return fCCalParameter ;}
  void SetACalParameter(Float_t a) { fACalParameter = a ;}
  void SetBCalParameter(Float_t b) { fBCalParameter = b ;}
  void SetCCalParameter(Float_t c) { fCCalParameter = c ;}


  Float_t GetEnergyAnalysisCut() {return  fEnergyAnalysisCut ;}
  void SetEnergyAnalysisCut(Float_t e) {  fEnergyAnalysisCut = e ;}

  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v1" ; }  
  
 private:

  const TString BranchName() const ; 
  virtual void Init() ;
  virtual void InitParameters() ;
  void     MakeRecParticles(void ) ;
  // Relative Distance CPV-EMC
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; 
  Int_t    GetPrincipalSign(Double_t* P, Int_t ell, Int_t eff_pur)const ; //Principal cut
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;
  // Gives in wich cluster energy range is the event
  void     GetClusterOption(const Float_t Cluster_En,const Bool_t range) ;
  // Gives the Efficiency-Purity point.
  Int_t    GetEffPurOption(const TString Eff_Pur)const ;
  virtual  void WriteRecParticles(Int_t event) ; 
  void     SetParameters() ; //Fills the matrix of parameters
  //Selects principal and parameters file in function of energy range.
  void     GetAnalysisParameters(Float_t ClusterEn) ;
  Double_t CalibratedEnergy(Float_t e)  ; //Calibrates energy.
 

 private:

  Bool_t                 fDefaultInit;        //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString    fFileName5 ;     // File that contains the Principal file for analysis from 0.5 to 5 GeV
  TString    fFileName100 ;   // File that contains the Principal file for analysis from 0.5 to 100 GeV
  TString    fFileNamePar5 ;  // File that contains the parameters for analysis from 0.5 to 5 GeV
  TString    fFileNamePar100 ;// File that contains the parameters for analysis from 0.5 to 100 GeV
  
 
  TString    fFrom ;              // name of Recpoints and TrackSegments 
  TString    fHeaderFileName ;    // file name with event header
  TString    fTrackSegmentsTitle; // branch name with track segments
  TString    fRecPointsTitle ;    // branch name with rec points
  TString    fRecParticlesTitle ; // branch name with rec particles
 
  Int_t                      fNEvent ;            //! current event number
  AliPHOSClusterizer *       fClusterizer ;       //! clusterizer
  AliPHOSTrackSegmentMaker * fTSMaker ;           //! track segment maker

  TPrincipal *               fPrincipal5 ;        //! TPrincipal from fFileName5  
  TPrincipal *               fPrincipal100 ;      //! TPrincipal from fFileName100 
  TPrincipal *               fPrincipal ;         //! TPrincipal copy 
  Double_t *                 fX ;                 //! Principal data 
  Double_t *                 fP ;                 //! Principal eigenvalues
 
  Int_t                      fRecParticlesInRun ; //! Total number of recparticles in one run
 
  TMatrixD *                 fParameters5 ;       //! Matrix of identification Parameters 0.5 to 5 GeV
  TMatrixD *                 fParameters100 ;     //! Matrix of identification Parameters 5-100 GeV
  TMatrixD *                 fParameters;         //! Matrix copy of identification Parameters
  Float_t                    fEnergyAnalysisCut;   // Energy to change from one PCA to the other.
  Int_t                      fCluster;            // Cluster energy range to choose parameters
  Int_t                      fClusterrcpv;        // Cluster energy range to choos rcpv parameters
  Int_t                      fMatrixExtraRow;     // Different size of the parameters file. Depends on range

  Float_t   fACalParameter ;// A parameter energy calibration Encal=A+B*En+C*En^2
  Float_t   fBCalParameter ;// B parameter energy calibration Encal=A+B*En+C*En^2
  Float_t   fCCalParameter ;// B parameter energy calibration Encal=A+B*En+C*En^2

  ClassDef( AliPHOSPIDv1,5)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
