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
    cerr << "ERROR: AliPHOSGetter is a singleton default ctor not callable" << endl ;
    abort() ; 
  } 
  AliPHOSGetter(const AliPHOSGetter & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    abort() ; 
  } 
  
  AliPHOSGetter & operator = (const AliPHOSGetter & ) {
    // assignement operator requested by coding convention, but not needed
    abort() ;
    return *this ; 
  }
  virtual ~AliPHOSGetter() ; 
  
  //=========== Instantiators ================
  static AliPHOSGetter * GetInstance(const char* headerFile,
				     const char* branchTitle = "Default",
                                     const Bool_t toSplit = kFALSE ) ; 
  static AliPHOSGetter * GetInstance() ; 
  
  //=========== General information about run ==============
  const Int_t  MaxEvent() const    { return static_cast<Int_t>(gAlice->TreeE()->GetEntries()) ; }
  const Int_t  EventNumber() const { return static_cast<Int_t>(gAlice->GetEvNumber()) ; }
  const Bool_t BranchExists(const TString tree) const ; 
  
  //========== PHOSGeometry and PHOS ============= 
  const AliPHOS *         PHOS() ;  
  const AliPHOSGeometry * PHOSGeometry() ; 
  
  //========== Methods to read something from file ==========
  void   Event(const Int_t event, const char * opt = "HSDRP") ;    
  void   Track(const Int_t itrack) ;
  void   ReadTreeS(TTree * treeS,Int_t input) ; //Method to be used when 
                                                //digitizing is under the control ofAliRunDigitizer, 
                                                //which opens all files etc.
  //========== Alarms ======================
  TFolder * Alarms() const { return dynamic_cast<TFolder*>(ReturnO("Alarms", 0)) ; }
  const TObjArray *  Alarms(const char * name ) const { return dynamic_cast<const TObjArray*>(ReturnO("Alarms", name)) ; }
  const TTask * QATasks(const char * name = 0) const { return ReturnT("QATasks", name) ; }
  
  //-----------------now getter's data--------------------------------------
  
  //=========== Primaries ============
  TTree *           TreeK(TString filename="") ; 
  TClonesArray *    Primaries(void) const { return dynamic_cast<TClonesArray*>(ReturnO("Primaries")) ; }
  const TParticle * Primary(Int_t index) const ;
  const Int_t       NPrimaries()const { return fNPrimaries; }
  const TParticle * Secondary(TParticle * p, Int_t index=1) const ;  
  
  //=========== Hits =================
  TTree *               TreeH(TString filename="") ; 
  const TClonesArray *  Hits(void) { return dynamic_cast<const TClonesArray*>(ReturnO("Hits")) ; }
  const AliPHOSHit *    Hit(Int_t index)  { return dynamic_cast<const AliPHOSHit*>(Hits()->At(index) );}
  
  //=========== SDigits ==============
  TTree *                    TreeS(TString filename="") ; 
  TClonesArray *             SDigits(const char * name = 0, const char * file=0) { 
    return dynamic_cast<TClonesArray*>(ReturnO("SDigits", name, file)) ;   }
  //const AliPHOSDigit *  SDigit(Int_t index) { return static_cast<const AliPHOSDigit *>(SDigits()->At(index)) ;} !!! why no such method ?
  const AliPHOSSDigitizer *  SDigitizer(const char * name =0) const { 
    return (const AliPHOSSDigitizer *) ReturnT("SDigitizer", name) ;   // here static or dynamic cast does not work ! why ?
  }
  
  //========== Digits ================
  TTree *                   TreeD(TString filename="") ; 
  TClonesArray *            Digits(const char * name = 0)const  { 
    return dynamic_cast<TClonesArray*>(ReturnO("Digits", name)) ;   }
  //const AliPHOSDigit *  Digit(Int_t index) { return static_cast<const AliPHOSDigit *>(Digits()->At(index)) ;} !!! why no such method ?
  const AliPHOSDigitizer *  Digitizer(const char * name = 0) const { 
    return (const AliPHOSDigitizer*)(ReturnT("Digitizer", name)) ;   }
  
  //========== RecPoints =============
  TObjArray *                EmcRecPoints(const char * name = 0) {
    return dynamic_cast<TObjArray*>(ReturnO("EmcRecPoints", name)) ;   }
  //const AliPHOSEmcRecPoint *  EmcRecPoint(Int_t index) { return static_cast<const AliPHOSEmcRecPoint *>(EmcRecPoints()->At(index)) ;} !!! why no such method ?
  TObjArray *                CpvRecPoints(const char * name = 0) { 
    return dynamic_cast<TObjArray*>(ReturnO("CpvRecPoints", name)) ;   }    
  const AliPHOSClusterizer * Clusterizer (const char * name =0) const { 
    return (const AliPHOSClusterizer*)(ReturnT("Clusterizer", name)) ;   // here static or dynamic cast does not work ! why ?
  }
  
  //========== TrackSegments ==========
  TClonesArray * TrackSegments(const char * name = 0) { 
    return static_cast<TClonesArray*>(ReturnO("TrackSegments", name)) ;   }
  const AliPHOSTrackSegmentMaker * TrackSegmentMaker (const char * name =0) const { 
    return (const AliPHOSTrackSegmentMaker*)(ReturnT("TrackSegmentMaker", name)) ;   }
  
  //========== RecParticles ===========
  TClonesArray * RecParticles(const char * name = 0) { 
    return static_cast<TClonesArray*>(ReturnO("RecParticles", name)) ;   }
  const AliPHOSPID * PID(const char * name =0) const { 
    return (const AliPHOSPID*)(ReturnT("PID", name)) ; } // here static or dynamic cast does not work ! why ? 
    
  //-----------------Auxiliary methods: cleaners-----------------
  void  RemoveTask(TString opt, TString name) const ;
  void  RemoveObjects(TString opt, TString name) const ;
  void  RemoveSDigits() const ;  

  //----------------Auxiliary methods: miscellana----------------
  void CloseFile() ;  
  const TFolder * Folder(const TString what) const ;
  const Bool_t HasFailed() const { return fFailed ; }
  void ListBranches(Int_t event=0) const ;
  void NewBranch(TString name, Int_t event = 0) ; 
  Bool_t NewFile(TString name) ;
  TFolder * SDigitsFolder() { return dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ; }  
  void SetDebug(Int_t level) {fDebug = level;} // Set debug level
  void SetRecParticlesTitle(const TString title) { fRecParticlesTitle = title ; }
  
  //------------Auxiliary methods: Posters--------------------
  const Bool_t PostPrimaries(void ) const ;  
  const Bool_t PostHits(void ) const ;  
  const Bool_t PostSDigits(      const char * name,  const char * file = 0) const ;  
  const Bool_t PostDigits(       const char * name ) const ;  
  const Bool_t PostRecPoints(    const char * name ) const ;  
  const Bool_t PostTrackSegments(const char * name) const ;  
  const Bool_t PostRecParticles( const char * name) const ;  
  const Bool_t PostClusterizer( const char * name) const ;  
  const Bool_t PostClusterizer(AliPHOSClusterizer * clu) const ;  
  const Bool_t PostSDigitizer (AliPHOSSDigitizer * sdigitizer) const ;  
  const Bool_t PostSDigitizer ( const char * name, const char * file ) const ;  
  const Bool_t PostDigitizer (AliPHOSDigitizer * digitizer) const ;  
  const Bool_t PostDigitizer  ( const char * name) const ;  
  const Bool_t PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tsm) const ;  
  const Bool_t PostTrackSegmentMaker(const char * name ) const ;  
  const Bool_t PostPID  (AliPHOSPID * pid) const ;  
  const Bool_t PostPID  (const char * name ) const ;  
  const Bool_t PostQA   (void) const ;

private:
  
  AliPHOSGetter(const char* headerFile, const char* branchTitle ="Default", const Bool_t toSplit = kFALSE) ; 
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

  static TFile * fFile;           //! 
  Bool_t         fToSplit ;              //! Do we work in the split mode
  TString        fHeaderFile ;           //! File in which gAlice lives
  TString        fBranchTitle ;          //!
  TString        fTrackSegmentsTitle ;   //! 
  TString        fTrackSegmentsFileName ;//! 
  TString        fRecPointsTitle ;       //!
  TString        fRecPointsFileName ;    //!
  TString        fRecParticlesTitle ;    //!
  TString        fRecParticlesFileName ; //!
  TString        fDigitsTitle ;          //! TDirectory tempo(gDirectory) 
  TString        fDigitsFileName ;       //! TDirectory tempo(gDirectory) 
  TString        fSDigitsTitle ;         //!
  TString        fSDigitsFileName ;      //!
  Bool_t         fFailed ;            //! set if file not opend or galice not found
  Int_t          fDebug ;             //! Debug level
  AliRun *       fAlice ;             //! needed to read TreeK if in an other file than fHeaderFile
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
