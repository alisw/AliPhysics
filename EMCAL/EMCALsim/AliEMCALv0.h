#ifndef ALIEMCALV0_H
#define ALIEMCALV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALv0
/// \ingroup EMCALsim
/// \brief EMCal simulation manager class v0
///
/// Implementation version v0 of EMCAL Manager class
/// An object of this class does not produce digits.
/// It is the one to use if you do not want to produce outputs in TREEH or TREED
///
/// This class places a Geometry of the EMCAL in the ALICE Detector as defined in AliEMCALGeometry.cxx
///
/// This class contains old geometry generation of EMCal prototypes (WSUC, 3x3 modules, 1x1 modules).
/// Should the removal be considered?
///
/// WARNING: Do not use for full detector simulations, use v2.
///
/// \author Yves Schutz (CCIN2P3)
/// \author Sahal Yacoob (LBL /UCT)
/// \author Alexei Pavlinov (WSU)
///
//_________________________________________________________________________

// --- ROOT system ---
class TFile;
class TList;
#include "TGeoManager.h"
#include <TArrayF.h>

// --- AliRoot header files ---
class AliEMCALShishKebabTrd1Module;
class AliEMCALSpaceFrame;
#include "AliEMCAL.h"

class AliEMCALv0 : public AliEMCAL 
{
public:
  
  AliEMCALv0();
  AliEMCALv0(const char *name, const char *title="",const Bool_t checkGeoAndRun = kTRUE) ;
  virtual ~AliEMCALv0(){} 
  
  using AliEMCAL::AddHit;
  
  virtual void   AddAlignableVolumes()        const;
  virtual void   AddAlignableVolumesInALICE() const;
  virtual void   AddAlignableVolumesInWSUC()  const;
  
  virtual void   CreateGeometry() ;// creates the geometry for GEANT
  virtual void   Init(void) ;      // does nothing
  
  // Gives the version number 
  virtual       Int_t  IsVersion(void) const { return 0 ; }
  virtual const TString  Version(void) const { return TString("v0") ; }
  
  // ShishKebab 
  void CreateShishKebabGeometry();
  void CreateSmod(const char* mother="XEN1");
  void CreateEmod(const char* mother="SMOD", const char* child="EMOD");
  void CreateAlFrontPlate(const char* mother="EMOD", const char* child="ALFP");
  
  // TRD1
  void Trd1Tower3X3(const Double_t *parSCM0);
  void PbInTrap(const Double_t parTRAP[11], TString n);
  
  // 1X1 case - Nov 22, 2006
  void Trd1Tower1X1(Double_t *parSCM0);
  void PbInTrd1(const Double_t *parTrd1, TString n);
  
  TList                        *GetShishKebabModules() const {return fShishKebabModules;}
  AliEMCALShishKebabTrd1Module *GetShishKebabModule(Int_t neta=0);
  
protected:
  
  TList   *fShishKebabModules; //!<! list of modules
  
private:
  
  TArrayF  fEnvelop1;          //!<! parameters of EMCAL envelop for TRD1(2) case 
  Int_t    fIdRotm;            //!<! number of rotation matrix (working variable)
  Int_t   *fIdTmedArr;         //!<! fIdtmed->GetArray() - 1599;
  
  Double_t fSampleWidth;       //!<! sample width = double(g->GetECPbRadThick()+g->GetECScintThick());
  Double_t fSmodPar0;          //!<! x size of super module  
  Double_t fSmodPar1;          //!<! y size of super module  
  Double_t fSmodPar2;          //!<! z size of super module  
  Double_t fInnerEdge;         //!<! Inner edge of DCAL super module 
  Double_t fParEMOD[5];        //!<! parameters of EMCAL module (TRD1,2)
  
  AliEMCALSpaceFrame* fCalFrame; ///< EMCAL Space frame object
  
  AliEMCALv0              (const AliEMCALv0 & emcal);
  AliEMCALv0 & operator = (const AliEMCALv0 & /*rvalue*/);
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALv0,4) ;
  /// \endcond
  
};
    
#endif // AliEMCALV0_H
