#ifndef ALIEMCALGETTER_H
#define ALIEMCALGETTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
//________________________________________________________________________
//  A singleton that returns various objects 
//  Should be used on the analysis stage to avoid confusing between different
//  branches of reconstruction tree: e.g. reading RecPoints and TS made from 
//  another set of RecPoints.
// 
//  The objects are retrived from folders.  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    

// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
 
// --- ROOT system ---
#include "TClonesArray.h"
#include "TFolder.h"  
#include "TTree.h"
#include "TFile.h"
class TString ;
class TParticle ;
class TTask ;

// --- Standard library ---
#include <stdlib.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliEMCAL.h" 
#include "AliEMCALHit.h" 
#include "AliEMCALDigit.h" 
#include "AliEMCALTowerRecPoint.h"
class AliEMCALGeometry ;
#include "AliEMCALDigitizer.h" 
#include "AliEMCALSDigitizer.h"
class AliEMCALClusterizer ;

class AliEMCALGetter : public TObject {
  
 public:
  
  AliEMCALGetter(){  // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    Fatal("ctor", "singleton: default ctor not callable") ;
  } 
  AliEMCALGetter(const AliEMCALGetter & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    Fatal("cpy ctor", "not implemented") ;  
  } 
  
  AliEMCALGetter & operator = (const AliEMCALGetter & ) {
    // assignement operator requested by coding convention, but not needed
    Fatal("operator =", "not implemented") ;  
    return *this ; 
  }
  virtual ~AliEMCALGetter() ; 
  
  //=========== Instantiators ================  
  static AliEMCALGetter * GetInstance(const char* headerFile,
				      const char* branchTitle = "Default", 
				      const Bool_t toSplit = kFALSE ) ; 
  static AliEMCALGetter * GetInstance() ; 
 
  //=========== General information about run ============== 
  const Int_t  MaxEvent() const    { return static_cast<Int_t>(gAlice->TreeE()->GetEntries()) ; }
  const Int_t  EventNumber() const { return static_cast<Int_t>(gAlice->GetEvNumber()) ; }
  const Bool_t BranchExists(const TString recName) const ;

  //========== EMCALGeometry and EMCAL ============= 
  const AliEMCAL *   EMCAL() ;  
  AliEMCALGeometry * EMCALGeometry() ; 

  //========== Methods to read somethig from file ==========
  void   Event(const Int_t event, const char * opt = "HSDRP") ;    
  void   Track(const Int_t itrack) ;
  void   ReadTreeS(TTree * treeS,Int_t input) ; //Method to be used when 
                                                //digitizing is under the conytrol of AliRunDigitizer, 
                                                //which opens all files etc.
  //-----------------now getter's data--------------------------------------

  //=========== Primaries ============
  TTree *           TreeK(TString filename="") ; 
  TClonesArray *    Primaries(void) const { return dynamic_cast<TClonesArray*>(ReturnO("Primaries")) ; }
  const TParticle * Primary(Int_t index) const;
  const Int_t       NPrimaries()const { return fNPrimaries; }
  const TParticle * Secondary(TParticle * p, Int_t index=1) const ;  

  //=========== Hits =================
  TTree *               TreeH(TString filename="") ; 
  const TClonesArray *  Hits(void) { return dynamic_cast<const TClonesArray*>(ReturnO("Hits")) ; }
  const AliEMCALHit *   Hit(Int_t index) { return dynamic_cast<const AliEMCALHit*>(Hits()->At(index) );}

  //=========== SDigits ==============
  TTree *         TreeS(TString filename="") ; 
  TClonesArray *  SDigits(const char * name = 0, const char * file=0) { 
    return dynamic_cast<TClonesArray*>(ReturnO("SDigits", name, file)) ; }
  const AliEMCALDigit *  SDigit(Int_t index) { return static_cast<const AliEMCALDigit*>(SDigits()->At(index)) ;}
  const AliEMCALSDigitizer *  SDigitizer(const char * name =0) const { 
    return (const AliEMCALSDigitizer*)(ReturnT("SDigitizer", name)) ; // here static or dynamic cast does not work ! why ?
  }

  //========== Digits ================
  TTree *         TreeD(TString filename="") ; 
  TClonesArray *  Digits(const char * name = 0)const  { 
    return dynamic_cast<TClonesArray*>(ReturnO("Digits", name)) ; }
  const AliEMCALDigit *  Digit(Int_t index) { return static_cast<const AliEMCALDigit *>(Digits()->At(index)) ;}
  const AliEMCALDigitizer *  Digitizer(const char * name =0) const { 
    return (const AliEMCALDigitizer*)(ReturnT("Digitizer", name)) ; }
      
  //========== RecPoints =============
  TObjArray * TowerRecPoints(const char * name = 0) const { 
    return (dynamic_cast<TObjArray*>(ReturnO("TowerRecPoints", name))) ; }
  const AliEMCALTowerRecPoint *  TowerRecPoint(Int_t index) { return static_cast<const AliEMCALTowerRecPoint *>(TowerRecPoints()->At(index)) ;}
  TObjArray * PreShowerRecPoints(const char * name = 0) const { 
    return (dynamic_cast<TObjArray*>(ReturnO("PreShowerRecPoints", name))) ; }
  const AliEMCALClusterizer * Clusterizer (const char * name =0) const { 
    return (const AliEMCALClusterizer*)(ReturnT("Clusterizer", name)) ;// here static or dynamic cast does not work ! why ?
  }

  //-----------------Auxiliary methods: cleaners-----------------
  void RemoveTask(TString opt, TString name) const ;
  void RemoveObjects(TString opt, TString name) const ; 
  void RemoveSDigits() const ; 

 //-----------------Auxiliary methods: miscellana-----------------
  void CloseFile() ;  
  const TFolder * Folder(const TString what) const ;
  const Bool_t HasFailed(void) const {return fFailed ;} 
  void ListBranches(Int_t event=0) const ;
  void NewBranch(TString name, Int_t event = 0) ; 
  Bool_t NewFile(TString name) ;
  TFolder * SDigitsFolder() { return dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ; }
  void SetDebug(Int_t level) {fDebug = level;} // Set debug level
  void SetRecParticlesTitle(const TString title) { fRecParticlesTitle = title ; }

  //------------Auxiliary methods: Posters--------------------
  const Bool_t PostPrimaries(void ) const ;  
  const Bool_t PostHits(void ) const ;  
  const Bool_t PostSDigits(      const char * name,  const char * file = 0) const ;  
  const Bool_t PostDigits(       const char * name ) const ;  
  const Bool_t PostRecPoints(    const char * name ) const ;  
  const Bool_t PostClusterizer( const char * name) const ;  
  const Bool_t PostClusterizer(AliEMCALClusterizer * clu) const ;  
  const Bool_t PostSDigitizer (AliEMCALSDigitizer * sdigitizer) const ;  
  const Bool_t PostSDigitizer ( const char * name, const char * file ) const ;  
  const Bool_t PostDigitizer (AliEMCALDigitizer * digitizer) const ;  
  const Bool_t PostDigitizer  ( const char * name) const ;  

private:

  AliEMCALGetter(const char* headerFile, const char* branchTitle ="Default", const Bool_t toSplit = kFALSE) ; 
  TObject * ReturnO(TString what, TString name=0, TString file=0) const ; 
  const TTask * ReturnT(TString what,TString name=0) const ; 
  void DefineBranchTitles(char* branch, char* branchTitle) ;
  Int_t ReadTreeD(const Int_t event) ;
  Int_t ReadTreeH(void) ;
  Int_t ReadTreeR(const Int_t event) ;
  Int_t ReadTreeS(const Int_t event) ;
  void ReadTreeQA(void) ;
  void ReadPrimaries(void) ;
  void CleanWhiteBoard(void) ;
  void CloseSplitFiles(void) ;
  void SetTitle(const char * title) ;

  TObject ** PrimariesRef(void) const ;
  TObject ** HitsRef(void) const ;
  TObject ** SDigitsRef(const char * name, const char * file = 0 ) const;
  TObject ** DigitsRef (const char * name)   const ;
  TObject ** TowerRecPointsRef (const char * name) const ;
  TObject ** PreShowerRecPointsRef (const char * name) const ;
  TObject ** SDigitizerRef (const char * name) const ; 
  TObject ** DigitizerRef  (const char * name) const ; 
  TObject ** ClusterizerRef(const char * name) const ; 

 private:

  static TFile * fFile ;              //!
  Bool_t         fToSplit ;           //! Do we work in the split mode
  TString        fHeaderFile ;        //! File in which gAlice lives
  TString        fBranchTitle ;       //!
  TString        fTrackSegmentsTitle ;//! 
  TString        fTrackSegmentsFileName ;//! 
  TString        fRecPointsTitle ;    //!
  TString        fRecPointsFileName ;    //!
  TString        fRecParticlesTitle ; //!
  TString        fRecParticlesFileName ; //!
  TString        fDigitsTitle ;       //!TDirectory tempo(gDirectory) 
  TString        fDigitsFileName ;    //! TDirectory tempo(gDirectory)  
  TString        fSDigitsTitle ;      //!
  TString        fSDigitsFileName ;      //!
  Bool_t         fFailed ;            //! true if file is not opened and/or galice not found
  Int_t          fDebug ;             // Debug level
  Int_t          fNPrimaries ;        //! # of primaries 
  TObjArray *    fPrimaries ;         //! list of lists of primaries-for the case of mixing
  TFolder *      fModuleFolder ;      //!Folder that contains the modules 
  TFolder *      fPrimariesFolder ;   //!Folder that contains the Primary Particles 
  TFolder *      fHitsFolder ;        //!Folder that contains the Hits 
  TFolder *      fSDigitsFolder ;     //!Folder that contains the SDigits 
  TFolder *      fDigitsFolder ;      //!Folder that contains the Digits 
  TFolder *      fRecoFolder ;        //!Folder that contains the reconstructed objects (RecPoints, TrackSegments, RecParticles) 
  TFolder *      fQAFolder ;          //!Folder that contains the QA objects  
  TFolder *      fTasksFolder ;       //!Folder that contains the Tasks (sdigitizer, digitizer, reconstructioner)
 
  static AliEMCALGetter * fgObjGetter; // pointer to the unique instance of the singleton 

  ClassDef(AliEMCALGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliEMCALGETTER_H

