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
// 
//  The objects are retrived from folders.  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    


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
//class AliEMCALEmcRecPoint ;
//class AliEMCALCpvRecPoint ;
//class AliEMCALClusterizer ;
//class AliEMCALTrackSegment ;
//class AliEMCALTrackSegmentMaker ;
//class AliEMCALRecParticle ;
//class AliEMCALPID ;

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
  //Bool_t PostRecPoints(    const char * name ) const ;  
  //Bool_t PostTrackSegments(const char * name) const ;  
  //Bool_t PostRecParticles( const char * name) const ;  

  //Bool_t PostClusterizer( const char * name) const ;  
  //Bool_t PostClusterizer(AliEMCALClusterizer * clu) const ;  
 Bool_t PostSDigitizer (AliEMCALSDigitizer * sdigitizer) const ;  
  Bool_t PostSDigitizer ( const char * name, const char * file ) const ;  
  Bool_t PostDigitizer (AliEMCALDigitizer * digitizer) const ;  
  Bool_t PostDigitizer  ( const char * name) const ;  
  //Bool_t PostTrackSegmentMaker(AliEMCALTrackSegmentMaker * tsm) const ;  
  //Bool_t PostTrackSegmentMaker(const char * name ) const ;  
  //Bool_t PostPID  (AliEMCALPID * pid) const ;  
  //Bool_t PostPID  (const char * name ) const ;  
  //Bool_t PostQA   (void) const ;
  

  void   Event(const Int_t event, const char * opt = "HSD") ;    
  void   Track(Int_t itrack) ;

  //Method to be used when digitizing under AliRunDigitizer, who opens all files etc.
  void   ReadTreeS(TTree * treeS,Int_t input) ;
  
  Int_t  EventNumber()       { return (Int_t) gAlice->GetEvNumber() ; }
  Int_t  MaxEvent()          { return (Int_t) gAlice->TreeE()->GetEntries() ; }
  static AliEMCALGetter * GetInstance(const char* headerFile,
				     const char* branchTitle = "Default" ) ; 
  static AliEMCALGetter *   GetInstance() ; 

  const AliEMCALv1 *         EMCAL()  ;  
  const  AliEMCALGeometry * EMCALGeometry() ; 
   // Alarms
  //TFolder * Alarms() const { return (TFolder*)(ReturnO("Alarms", 0)) ; }
  //TObjArray *  Alarms(const char * name ) const { return (TObjArray*)(ReturnO("Alarms", name)) ; }

  // QA Tasks
  //TTask * QATasks(const char * name = 0) const { return (TTask*)(ReturnT("QATasks", name)) ; }

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

  // RecPoints
  //TObjArray * EmcRecPoints(const char * name = 0) const { 
    //          return (TObjArray*)(ReturnO("EmcRecPoints", name)) ; }
  //TObjArray * CpvRecPoints(const char * name = 0) const { 
    //          return (TObjArray*)(ReturnO("CpvRecPoints", name)) ; }

  //AliEMCALClusterizer * Clusterizer (const char * name =0) const 
    //          { return (AliEMCALClusterizer*)(ReturnT("Clusterizer", name)) ; }

  // TrackSegments
  //TClonesArray * TrackSegments(const char * name = 0) const 
   //                { return (TClonesArray*)(ReturnO("TrackSegments", name)) ; }
  //AliEMCALTrackSegmentMaker * TrackSegmentMaker (const char * name =0) const 
    //               { return (AliEMCALTrackSegmentMaker*)(ReturnT("TrackSegmentMaker", name)) ; }

  // RecParticles
  //TClonesArray * RecParticles(const char * name = 0) const  
    //               { return (TClonesArray*)(ReturnO("RecParticles", name)) ; }
    //AliEMCALPID * PID(const char * name =0) const 
      //             { return (AliEMCALPID*)(ReturnT("PID", name)) ; }

  // Primaries
  const TParticle *           Primary(Int_t index) const ;
  const Int_t                 NPrimaries()const { return fNPrimaries; }

  void  SetDebug(Int_t level) {fDebug = level;} // Set debug level

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
  //void ReadTreeR() ;
  void ReadTreeS(Int_t event) ;
  //void ReadTreeQA() ;
  void ReadPrimaries() ;

  void * HitsRef(void) const ;
  void * SDigitsRef(const char * name, const char * file = 0 ) const;
  void * DigitsRef (const char * name)   const ;
  //void * EmcRecPointsRef (const char * name) const ;
  //void * CpvRecPointsRef (const char * name) const ;
  //void * TrackSegmentsRef(const char * name)   const ;
  //void * RecParticlesRef (const char * name)   const ;
  //void * AlarmsRef (void)   const ;

  void * SDigitizerRef (const char * name) const ; 
  void * DigitizerRef  (const char * name) const ; 
  //void * ClusterizerRef(const char * name) const ; 
  //void * TSMakerRef    (const char * name) const ; 
  //void * PIDRef        (const char * name) const ; 

 private:

  TString        fHeaderFile ;        //! File in which gAlice lives
  TString        fBranchTitle ;       //!
  //TString        fTrackSegmentsTitle ;//! 
  //TString        fRecPointsTitle ;    //!
  //TString        fRecParticlesTitle ; //!
  TString        fDigitsTitle ;       //!
  TString        fSDigitsTitle ;      //!

  Int_t          fDebug ;             // Debug level

  Int_t          fNPrimaries ;        //! # of primaries
  
  TObjArray *    fPrimaries ;         //! list of lists of primaries-for the case of mixing

  TFolder *      fModuleFolder ;      //!Folder that contains the modules 
  TFolder *      fHitsFolder ;        //!Folder that contains the Hits 
  TFolder *      fSDigitsFolder ;     //!Folder that contains the SDigits 
  TFolder *      fDigitsFolder ;      //!Folder that contains the Digits 
  //TFolder *      fRecoFolder ;        //!Folder that contains the reconstructed objects (RecPoints, TrackSegments, RecParticles) 
  //TFolder *      fQAFolder ;          //!Folder that contains the QA objects  
  TFolder *      fTasksFolder ;       //!Folder that contains the Tasks (sdigitizer, digitizer, reconstructioner)
 
  static AliEMCALGetter * fgObjGetter; // pointer to the unique instance of the singleton 

  ClassDef(AliEMCALGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliEMCALGETTER_H
