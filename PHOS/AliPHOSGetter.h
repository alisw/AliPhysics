#ifndef ALIPHOSGETTER_H
#define ALIPHOSGETTER_H
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
#include "TFile.h"
class TString ;
class TParticle ;
class TTask ;

// --- Standard library ---
#include <stdlib.h>
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOS.h" 
#include "AliPHOSHit.h" 
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"
class AliPHOSGeometry ;
class AliPHOSDigitizer ;
class AliPHOSSDigitizer ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSPID ;

class AliPHOSGetter : public TObject {
  
 public:
  
  AliPHOSGetter(){    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    cerr << "ERROR: AliPHOGetter is a singleton default ctor not callable" << endl ;
    abort() ; 
  } 
  AliPHOSGetter(const AliPHOSGetter & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    abort() ; 
  } 
  
  virtual ~AliPHOSGetter() ; 
  
  void CloseFile() ;  
  const TFolder * Folder(const TString what) const ;
  void ListBranches(Int_t event=0) const ;
  void NewBranch(TString name, Int_t event = 0) ; 
  Bool_t NewFile(TString name) ;
  const Bool_t HasFailed() const { return fFailed ; }
  Bool_t PostPrimaries(void ) const ;  
  Bool_t PostHits(void ) const ;  
  Bool_t PostSDigits(      const char * name,  const char * file = 0) const ;  
  Bool_t PostDigits(       const char * name ) const ;  
  Bool_t PostRecPoints(    const char * name ) const ;  
  Bool_t PostTrackSegments(const char * name) const ;  
  Bool_t PostRecParticles( const char * name) const ;  

  Bool_t PostClusterizer( const char * name) const ;  
  Bool_t PostClusterizer(AliPHOSClusterizer * clu) const ;  
  Bool_t PostSDigitizer (AliPHOSSDigitizer * sdigitizer) const ;  
  Bool_t PostSDigitizer ( const char * name, const char * file ) const ;  
  Bool_t PostDigitizer (AliPHOSDigitizer * digitizer) const ;  
  Bool_t PostDigitizer  ( const char * name) const ;  
  Bool_t PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tsm) const ;  
  Bool_t PostTrackSegmentMaker(const char * name ) const ;  
  Bool_t PostPID  (AliPHOSPID * pid) const ;  
  Bool_t PostPID  (const char * name ) const ;  
  Bool_t PostQA   (void) const ;
  

  void   Event(const Int_t event, const char * opt = "HSDRP") ;    
  void   Track(Int_t itrack) ;

  //Method to be used when digitizing under AliRunDigitizer, who opens all files etc.
  void   ReadTreeS(TTree * treeS,Int_t input) ;
  
  Int_t  EventNumber()       { return (Int_t) gAlice->GetEvNumber() ; }
  Int_t  MaxEvent()          { return (Int_t) gAlice->TreeE()->GetEntries() ; }
  static AliPHOSGetter * GetInstance(const char* headerFile,
				     const char* branchTitle = "Default" ) ; 
  static AliPHOSGetter *   GetInstance() ; 

  const AliPHOS *         PHOS()  ;  
  const  AliPHOSGeometry * PHOSGeometry() ; 
   // Alarms
  TFolder * Alarms() const { return (TFolder*)(ReturnO("Alarms", 0)) ; }
  TObjArray *  Alarms(const char * name ) const { return (TObjArray*)(ReturnO("Alarms", name)) ; }

  // QA Tasks
  TTask * QATasks(const char * name = 0) const { return (TTask*)(ReturnT("QATasks", name)) ; }

  // Primaries
  TClonesArray *  Primaries(void) const { return (TClonesArray*)(ReturnO("Primaries")) ; }
  // Hits
  const TClonesArray *  Hits(void) { return static_cast<const TClonesArray*>(ReturnO("Hits")) ; }
  const AliPHOSHit * Hit(Int_t index)  { return static_cast<const AliPHOSHit*>(Hits()->At(index) );}
  
  // SDigits
  TClonesArray *  SDigits(const char * name = 0, const char * file=0) { 
    return static_cast<TClonesArray*>(ReturnO("SDigits", name, file)) ; 
  }
  const AliPHOSDigit *  SDigit(Int_t index) { return static_cast<const AliPHOSDigit*>(SDigits()->At(index)) ;}
  
  AliPHOSSDigitizer *  SDigitizer(const char * name =0) const { 
    return ((AliPHOSSDigitizer*)(ReturnT("SDigitizer", name))) ; 
  }
  
  // Digits
  TClonesArray *  Digits(const char * name = 0)const  { 
    return static_cast<TClonesArray*>(ReturnO("Digits", name)) ; 
  }
  const AliPHOSDigit *  Digit(Int_t index) { return static_cast<const AliPHOSDigit *>(Digits()->At(index)) ;}
  AliPHOSDigitizer *  Digitizer(const char * name =0) const { 
    return (AliPHOSDigitizer*)(ReturnT("Digitizer", name)) ; 
  }
  
  // RecPoints
  TObjArray * EmcRecPoints(const char * name = 0) {
    return static_cast<TObjArray*>(ReturnO("EmcRecPoints", name)) ; 
  }
  TObjArray * CpvRecPoints(const char * name = 0) { 
    return static_cast<TObjArray*>(ReturnO("CpvRecPoints", name)) ; 
  }
  const AliPHOSEmcRecPoint * EmcRecPoint(Int_t index) { 
    return static_cast<const AliPHOSEmcRecPoint *>(EmcRecPoints()->At(index)) ;
  }
  const AliPHOSCpvRecPoint * CpvRecPoint(Int_t index) { 
    return static_cast<const AliPHOSCpvRecPoint *>(CpvRecPoints()->At(index)) ;  
  }
    
  AliPHOSClusterizer * Clusterizer (const char * name =0) const { 
    return (AliPHOSClusterizer*)(ReturnT("Clusterizer", name)) ; 
  }
  
  // TrackSegments
  TClonesArray * TrackSegments(const char * name = 0) { 
    return static_cast<TClonesArray*>(ReturnO("TrackSegments", name)) ; 
  }
  const AliPHOSTrackSegment * TrackSegment(Int_t index) { 
    return static_cast<const AliPHOSTrackSegment*>(TrackSegments()->At(index)) ; 
  }
  AliPHOSTrackSegmentMaker * TrackSegmentMaker (const char * name =0) const { 
    return (AliPHOSTrackSegmentMaker*)(ReturnT("TrackSegmentMaker", name)) ; 
  }
  
  // RecParticles
  TClonesArray * RecParticles(const char * name = 0) { 
    return static_cast<TClonesArray*>(ReturnO("RecParticles", name)) ; 
  }
  const AliPHOSRecParticle * RecParticle(Int_t index) { 
    return static_cast<const AliPHOSRecParticle*>(RecParticles()->At(index)) ; 
  }
  AliPHOSPID * PID(const char * name =0) const { 
    return (AliPHOSPID*)(ReturnT("PID", name)) ; 
  }
  
  // Primaries
  const TParticle *           Primary(Int_t index) const ;
  const Int_t                 NPrimaries()const { return fNPrimaries; }
  const TParticle *           Secondary(TParticle * p, Int_t index=1) const ;
  
  void  RemoveTask(TString opt, TString name) const ;
  void  RemoveObjects(TString opt, TString name) const ;
  void  RemoveSDigits() const ;
  void  SetDebug(Int_t level) {fDebug = level;} // Set debug level
  
  AliPHOSGetter & operator = (const AliPHOSGetter & ) {
    // assignement operator requested by coding convention, but not needed
    abort() ;
    return *this ; 
  }
  
  TFolder * SDigitsFolder() { return dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ; }

  void SetRecParticlesTitle(const TString title) { fRecParticlesTitle = title ; }
  
private:
  
  AliPHOSGetter(const char* headerFile, const char* branchTitle ="Default") ; 
  TObject * ReturnO(TString what, TString name=0, TString file=0) const ; 
  const TTask * ReturnT(TString what,TString name=0) const ; 
  void DefineBranchTitles(char* branch, char* branchTitle) ;
  Int_t ReadTreeD() ;
  Int_t ReadTreeH() ;
  Int_t ReadTreeR(Bool_t any=kFALSE) ;
  Int_t ReadTreeS(Int_t event) ;
  void ReadTreeQA() ;
  void ReadPrimaries() ;

  TObject** PrimariesRef(void) const ;
  TObject** HitsRef(void) const ;
  TObject** SDigitsRef(const char * name, const char * file = 0 ) const;
  TObject** DigitsRef (const char * name)   const ;
  TObject** EmcRecPointsRef (const char * name) const ;
  TObject** CpvRecPointsRef (const char * name) const ;
  TObject** TrackSegmentsRef(const char * name)   const ;
  TObject** RecParticlesRef (const char * name)   const ;
  TObject** AlarmsRef (void)   const ;

  TObject** SDigitizerRef (const char * name) const ; 
  TObject** DigitizerRef  (const char * name) const ; 
  TObject** ClusterizerRef(const char * name) const ; 
  TObject** TSMakerRef    (const char * name) const ; 
  TObject** PIDRef        (const char * name) const ; 

 private:

  static TFile *        fFile;               //! 
  TString        fHeaderFile ;        //! File in which gAlice lives
  TString        fBranchTitle ;       //!
  TString        fTrackSegmentsTitle ;//! 
  TString        fRecPointsTitle ;    //!
  TString        fRecParticlesTitle ; //!
  TString        fDigitsTitle ;       //! TDirectory tempo(gDirectory) ; 

  TString        fSDigitsTitle ;      //!

  Bool_t         fFailed ;            //! set if file not opend or galice not found
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
 
  static AliPHOSGetter * fgObjGetter; // pointer to the unique instance of the singleton 

  ClassDef(AliPHOSGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliPHOSGETTER_H
