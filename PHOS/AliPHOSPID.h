#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class                             
//  of identified particles                
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TTask.h" 
class TFormula ;
class TClonesArray ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSGeometry ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;

class AliPHOSPID : public TTask {

 public:

  AliPHOSPID() ;          // ctor            
  AliPHOSPID(const char* headerFile,const char * name, const Bool_t toSplit) ;
  virtual ~AliPHOSPID() ; // dtor

  virtual void Exec(Option_t * option) { cout << "AliPHOSPID::Exec not define " << endl ; }
  //  virtual char * GetRecParticlesBranch()const  { cout << "AliPHOSPID::GetRecParticlesBranch not defined " << endl ; return 0 ; }
  //  virtual char * GetTrackSegmentsBranch()const { cout << "AliPHOSPID::GetTrackSegmentsBranch not defined " << endl ; return 0 ; } 
  virtual const Int_t GetRecParticlesInRun()  const { cout << "AliPHOSPID:::GetRecParticlesInRun not defined " << endl ; return 0 ;} 
  virtual void Print(Option_t * option) const { cout << "AliPHOSPID::Print not defined " << endl ;}
  //virtual void PlotDispersionCuts()const = 0;
  //virtual void SetIdentificationMethod(char * option) = 0 ;
  //virtual void SetShowerProfileCut(char *  formula) = 0  ; 
  //virtual void SetDispersionCut(Float_t cut) = 0  ;   
  virtual void SetCpvtoEmcDistanceCut(Float_t Cluster_En, TString Eff_Pur,Float_t cut ) { cout << "AliPHOSPID::SetCpvtoEmcDistanceCut not defined " << endl ;}
  virtual void SetTimeGate(Float_t Cluster_En, TString Eff_Pur, Float_t gate) { cout << "AliPHOSPID::SetTimeGate not defined " << endl ; }
  //  virtual void SetTrackSegmentsBranch(const char* title) { cout << "AliPHOSPID::Exec not define " << endl ; }
  //  virtual void SetRecParticlesBranch (const char* title) { cout << "AliPHOSPID::SetTecParticlesBranch not defined " << endl ; }
  //  virtual void SetSplitFile(const TString splitFileName = "PHOS.RecData.root") const ; 
  virtual const char * Version() const { cout << "AliPHOSPID::Version not defined " << endl ; return 0 ; }  
  virtual void WriteRecParticles(Int_t event) { cout << "AliPHOSPID::WriteRecParticles not defined " << endl ; }

private: 
  virtual void Init() { cout << "AliPHOSPID::Init not define " << endl ; } 

protected:

  TFile * fSplitFile ;             //! file in which RecParticles will eventually be stored
  Bool_t  fToSplit   ;             //! do we in the split mode  
  ClassDef(AliPHOSPID,1)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
