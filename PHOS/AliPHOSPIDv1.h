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
  Double_t GetCpvtoEmcDistanceCut(const Float_t Cluster_En, const TString Eff_Pur)const  ;
  Double_t GetTimeGate(const Float_t Cluster_En, const TString Eff_Pur)const  ;
  const TString GetPrincipalFile( )const {return fFileName ;}
  const TString GetPrincipalFilePar( )const {return fFileNamePar ;}

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
  void SetParameters(TString OptFileName) ; //Fills the matrix of parameters
 
  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v1" ; }  
  
 private:

  const TString AliPHOSPIDv1::BranchName() const ; 
  virtual void Init() ;
  virtual void InitParameters() ;
  void     MakeRecParticles(void ) ;
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; // Relative Distance CPV-EMC
  Int_t    GetPrincipalSign(Double_t* P, Int_t ell, Int_t eff_pur)const ; //Principal cut
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;
  //void     SetParameters() ; //Fills the matrix of parameters
  Int_t    GetClusterOption(const Float_t Cluster_En,const Bool_t rcpv)const ;// Gives in wich cluster energy range is the event
  Int_t    GetEffPurOption(const TString Eff_Pur)const ;// Gives the Efficiency-Purity point.
  virtual void WriteRecParticles(Int_t event) ; 
  void  SetPrincipalFileOptions(TString OptFileName) ;

 private:

  TString                fFileName ;          // Name of the file which contains the Principal file
  TString                fFileNamePar ;       // Name of the file which contains the parameters
  Int_t                  fCluster ;           // Parameter that depens on the energy range to select the parameters of the file       
  TString                fOptFileName ;       // choose de parameters and principal file
  TString                fFrom ;              // name of Recpoints and TrackSegments 
  TString                fHeaderFileName ;    // file name with event header
  TString                fTrackSegmentsTitle; // branch name with track segments
  TString                fRecPointsTitle ;    // branch name with rec points
  TString                fRecParticlesTitle ; // branch name with rec particles
 
  Int_t                      fNEvent ;            //! current event number
  AliPHOSClusterizer *       fClusterizer ;       //! clusterizer
  AliPHOSTrackSegmentMaker * fTSMaker ;           //! track segment maker
  TPrincipal *               fPrincipal ;         //! TPrincipal from fFileName 
  Int_t                      fRecParticlesInRun ; //! Total number of recparticles in one run
  Double_t *                 fX ;                 //! Principal data 
  Double_t *                 fP ;                 //! Principal eigenvalues
  TMatrixD *                 fParameters ;        //! Matrix of all identification Parameters

  ClassDef( AliPHOSPIDv1,3)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
