#ifndef ALIPHOSV0_H
#define ALIPHOSV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.44  2006/09/27 19:55:57  kharlov
 * Alignment object with symbolic volume names are introduced
 *
 * Revision 1.43  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// Layout EMC + CPV  has name IHEP
//*--                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

class TFile;
class TFolder;

// --- AliRoot header files ---
#include "AliPHOS.h"

class AliPHOSv0 : public AliPHOS {

 public:

  AliPHOSv0() {}
  AliPHOSv0(const char *name, const char *title="") ;
  virtual ~AliPHOSv0(void){
    // dtor
  } 

//    virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits ) {
  //this function is not a final-overrider for AliPHOS::AddHit, to
  //supress warning, I use using-declaration :)
  using AliPHOS::AddHit;
  virtual void   AddHit( Int_t, Int_t, Int_t, Int_t, Float_t*) {
    // useless since there are no hits
    Fatal("AddHit", "not to be used with v0") ;
  }
  virtual void   BuildGeometry(void) ;             // creates the geometry for the ROOT display
  void           BuildGeometryforEMC(void) ;      // creates the PHOS geometry for the ROOT display
  //  void           BuildGeometryforPPSD(void) ;      // creates the PPSD geometry for the ROOT display
  void           BuildGeometryforCPV(void) ;       // creates the CPV  geometry for the ROOT display
  virtual void   CreateGeometry(void) ;            // creates the geometry for GEANT
  void           CreateGeometryforEMC(void) ;     // creates the PHOS geometry for GEANT
  //  void           CreateGeometryforPPSD(void) ;     // creates the PPSD geometry for GEANT
  void           CreateGeometryforCPV(void) ;      // creates the CPV  geometry for GEANT
  void           CreateGeometryforSupport(void) ;  // creates the Support geometry for GEANT
  virtual void   AddAlignableVolumes() const;      // define sym.names for alignable volumes

  virtual Float_t ZMin() const;                    // overall dimension of the module (min)
  virtual Float_t ZMax() const;                    // overall dimension of the module (max)

  virtual void   Init(void) ;                      // does nothing
  virtual Int_t  IsVersion(void) const { 
    // Gives the version number 
    return 0 ; 
  }
  virtual const TString Version(void)const { 
    // As above
    return TString("v0") ; 
  }
  
  
 private:
  AliPHOSv0(AliPHOSv0 & phos);
  AliPHOSv0 & operator = (const AliPHOSv0 & /*rvalue*/);

  ClassDef(AliPHOSv0,1)  // Implementation of PHOS manager class for layout EMC+PPSD
    
    };
    
#endif // AliPHOSV0_H
