#ifndef ALIEMCALGETTER_H
#define ALIEMCALGETTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A singleton that returns various objects 
//  Should be used on the analysis stage to avoid confusing between different
//  branches of reconstruction tree: e.g. reading RecPoints and TS made from 
//  another set of RecPoints.
//  At this stage the Getter class handles only Hits, Digits, and SDigits.
//  The objects are retrived from folders.  
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSGetter    


// --- ROOT system ---
#include "TClonesArray.h"
#include "TFolder.h"  
#include "TTree.h"
class TString ;
class TParticle ;
class TTask ;

// --- Standard library ---
#include <stdlib.h>
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliEMCALv1.h" 
class AliEMCALGeometry ;
class AliEMCALHit ;
class AliEMCALDigit ;
class AliEMCALDigitizer ;
class AliEMCALSDigitizer ;

class AliEMCALGetter : public TObject {
  
 public:
  
  AliEMCALGetter(){ 
    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    cerr << "ERROR: AliPHOGetter is a singleton default ctor not callable" << endl ;
    abort() ; 
  } 
  AliEMCALGetter(const AliEMCALGetter & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    abort() ; 
  } 
  
  virtual ~AliEMCALGetter() ; 
  
  Bool_t PostHits(void ) const ;  
  Bool_t PostSDigits(      const char * name,  const char * file = 0) const ;  
  Bool_t PostDigits(       const char * name ) const ;  

  Bool_t PostSDigitizer (AliEMCALSDigitizer * sdigitizer) const ;  
  Bool_t PostSDigitizer ( const char * name, const char * file ) const ;  
  Bool_t PostDigitizer (AliEMCALDigitizer * digitizer) const ;  
  Bool_t PostDigitizer  ( const char * name) const ;  
  

  void   Event(const Int_t event, const char * opt = "HSD") ;    
  void   Track(Int_t itrack) ;

  //Method to be used when digitizing under AliRunDigitizer, who opens all files etc.
  void   ReadTreeS(TTree * treeS,Int_t input) ;
  
  Int_t  EventNumber()       { return (Int_t) gAlice->GetEvNumber() ; }
  Int_t  MaxEvent()          { return (Int_t) gAlice->TreeE()->GetEntries() ; }
  static AliEMCALGetter * GetInstance(const char* headerFile,
				     const char* branchTitle = "Default" ) ; 
  static AliEMCALGetter *   GetInstance() ; 

  const AliEMCALv0 *         EMCAL()  ;  
  const  AliEMCALGeometry * EMCALGeometry() ; 

  // Hits
        TClonesArray *  Hits(void) const { return (TClonesArray*)(ReturnO("Hits")) ; }

  // SDigits
        TClonesArray *  SDigits(const char * name = 0, const char * file=0) const 
	                              { return (TClonesArray*)(ReturnO("SDigits", name, file)) ; }

   AliEMCALSDigitizer *  SDigitizer(const char * name =0) const 
                                      { return ((AliEMCALSDigitizer*)(ReturnT("SDigitizer", name))) ; }

  // Digits
        TClonesArray *  Digits(const char * name = 0)   const 
                             { return (TClonesArray*)(ReturnO("Digits", name)) ; }
    AliEMCALDigitizer *  Digitizer(const char * name =0) const 
                             { return (AliEMCALDigitizer*)(ReturnT("Digitizer", name)) ; }

  // Primaries
  const TParticle *           Primary(Int_t index) const ;
  const Int_t                 NPrimaries()const { return fNPrimaries; }


  AliEMCALGetter & operator = (const AliEMCALGetter & ) {
    // assignement operator requested by coding convention, but not needed
    abort() ;
    return *this ; 
  }
  
  TFolder * SDigitsFolder() { return dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ; }

 private:

  AliEMCALGetter(const char* headerFile, const char* branchTitle ="Default") ; 
  void CreateWhiteBoard() const ; 
  const TObject * ReturnO(TString what, TString name=0, TString file=0) const ; 
  const TTask * ReturnT(TString what,TString name=0) const ; 
  void DefineBranchTitles(char* branch, char* branchTitle) ;
  void ReadTreeD() ;
  void ReadTreeH() ;
  void ReadTreeS(Int_t event) ;
  void ReadPrimaries() ;

  void * HitsRef(void) const ;
  void * SDigitsRef(const char * name, const char * file = 0 ) const;
  void * DigitsRef (const char * name)   const ;

  void * SDigitizerRef (const char * name) const ; 
  void * DigitizerRef  (const char * name) const ; 

 private:

  TString        fHeaderFile ;        //! File in which gAlice lives
  TString        fBranchTitle ;       //!
  TString        fDigitsTitle ;       //!
  TString        fSDigitsTitle ;      //!

  Int_t          fDebug ;             // Debug level

  Int_t          fNPrimaries ;        //! # of primaries
  
  TObjArray *    fPrimaries ;         //! list of lists of primaries-for the case of mixing

  TFolder *      fHitsFolder ;        //!Folder that contains the Hits 
  TFolder *      fSDigitsFolder ;     //!Folder that contains the SDigits 
  TFolder *      fDigitsFolder ;      //!Folder that contains the Digits 
  TFolder *      fTasksFolder ;       //!Folder that contains the Tasks (sdigitizer, digitizer, reconstructioner)
  TFolder *      fModuleFolder ;     //!
  static AliEMCALGetter * fgObjGetter; // pointer to the unique instance of the singleton 

  ClassDef(AliEMCALGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliEMCALGETTER_H
