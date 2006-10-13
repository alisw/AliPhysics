#ifndef ALIEMCALV0_H
#define ALIEMCALV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                      
         */
/* $Id$ */

//_________________________________________________________________________
// Implementation version v0 of EMCAL Manager class 
//*--                  
//*-- Author: Yves Schutz (SUBATECH)
//*-- and   : Sahal Yacoob (LBL / UCT) 
//          : Aleksei Pavlinov (WSU)     SHASHLYK

// --- ROOT system ---

class TFile;
class TList;
class TNode;
class AliEMCALShishKebabTrd1Module;

// --- AliRoot header files ---
#include "AliEMCAL.h"
#include "TGeoManager.h"
#include <TArrayF.h>
//class AliEMCALGeometry ; 

class AliEMCALv0 : public AliEMCAL {

 public:

  AliEMCALv0();
  AliEMCALv0(const char *name, const char *title="") ;
  virtual ~AliEMCALv0(){} 

  AliEMCALv0(const AliEMCALv0 & emcal);

  AliEMCALv0 & operator = (const AliEMCALv0  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this; 
  }
 
  using AliEMCAL::AddHit;

  virtual void  AddAlignableVolumes() const;

  virtual void BuildGeometry();// creates the geometry for the ROOT display
  TNode *BuildGeometryOfWSUC();  // WSUC - test environment
  virtual void CreateGeometry() ;// creates the geometry for GEANT
  virtual void   Init(void) ;                                       // does nothing
  virtual Int_t  IsVersion(void) const { 
    // Gives the version number 
    return 0 ; 
  }
  virtual const TString Version(void) const{ 
    // As above
    return TString("v0") ; 
  }
  
  // ShishKebab 
  void CreateShishKebabGeometry();
  void CreateSmod(const char* mother="XEN1");
  void CreateEmod(const char* mother="SMOD", const char* child="EMOD");
  // TRD1
  void Trd1Tower3X3(const double parSCM0[5]);
  void Trd1Tower4X4();
  void PbInTrap(const double parTRAP[11], TString n);
  // TRD2 - 1th design
  void Scm0InTrd2(const AliEMCALGeometry * g, const Double_t emodPar[5], Double_t parSCM0[5]);
  void Division2X2InScm0(const AliEMCALGeometry * g, const Double_t parSCM0[5]);
  void PbInTrapForTrd2(const double *parTRAP, TString name);
  // TRD2 - 2th design
  void PbmoInTrd2(const AliEMCALGeometry * g, const Double_t emodPar[5], Double_t parPBMO[5]);
  void Division2X2InPbmo(const AliEMCALGeometry * g, const Double_t parPBMO[5]);

  TList  *GetShishKebabModules() {return fShishKebabModules;}
  AliEMCALShishKebabTrd1Module *GetShishKebabModule(Int_t neta=0);

 protected:
  TList   *fShishKebabModules; //! list of modules
 private:
  TArrayF  fEnvelop1;          //! parameters of EMCAL envelop for TRD1(2) case 
  Int_t    fIdRotm;            //! number of rotation matrix (working variable)
  Int_t   *fIdTmedArr;         //! fIdtmed->GetArray() - 1599;
  Double_t fSampleWidth;       //! sample width = double(g->GetECPbRadThick()+g->GetECScintThick());
  Double_t fSmodPar0;          //! x size of super module  
  Double_t fSmodPar1;          //! y size of super module  
  Double_t fSmodPar2;          //! z size of super module  
  Double_t fParEMOD[5];        //! parameters of EMCAL module (TRD1,2)

  ClassDef(AliEMCALv0,3) // Implementation of EMCAL manager class for midrapidity barrel layout between 80 and 180(190) degrees 
    
    };
    
#endif // AliEMCALV0_H
