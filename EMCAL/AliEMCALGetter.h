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
#include "TObject.h"  
#include "TClonesArray.h" 
// #include "TFolder.h"  
// #include "TTree.h"
// #include "TFile.h"
// class TString ;
 class TParticle ;
// class TTask ;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h" 

// #include "AliRun.h"
class AliEMCAL ; 
#include "AliEMCALHit.h"   
// #include "AliEMCALDigit.h"
// #include "AliEMCALRecPoint.h"
// #include "AliEMCALRecPoint.h"
// #include "AliEMCALTrackSegment.h"
// #include "AliEMCALRecParticle.h"
class AliEMCALGeometry ;
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h" 
// class AliEMCALClusterizer ;
// class AliEMCALTrackSegmentMaker ;
// class AliEMCALPID ;
// //class AliEMCALCalibrationDB ;
// class AliEMCALConTableDB ;
class AliEMCALBeamTestEvent ;

#include "AliEMCALLoader.h" 

class AliEMCALGetter : public TObject {
  
 public:  
  AliEMCALGetter(){    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    Fatal("ctor", "AliEMCALGetter is a singleton default ctor not callable") ;
  } 
  AliEMCALGetter(const AliEMCALGetter & obj):TObject(obj) {
    // cpy ctor requested by Coding Convention 
    Fatal("cpy ctor", "not implemented") ;
  } 
  
  AliEMCALGetter & operator = (const AliEMCALGetter & ) {
    // assignement operator requested by coding convention, but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  virtual ~AliEMCALGetter() ; 
  
  //=========== Instantiators ================
  static AliEMCALGetter * Instance(const char* headerFile,
				  const char* version = AliConfig::fgkDefaultEventFolderName,
				  Option_t * openingOption = "READ" ) ; 
  static AliEMCALGetter * Instance() ; 

  static void Print() ; 
  
//   //=========== General information about run ==============
  Bool_t IsLoaded(const TString tree) const { return fLoadingStatus.Contains(tree) ; } 
  void   SetLoaded(const TString tree) { fLoadingStatus += tree ; } 

  Int_t  MaxEvent() const ; 
  Int_t  EventNumber() const ; 
  Bool_t VersionExists(TString & opt) const ; 
  UShort_t EventPattern(void) const ; 
  Float_t  BeamEnergy(void) const ;
  
//   //========== EMCALGeometry and EMCAL ============= 
  AliEMCAL *         EMCAL() const  ;  
  AliEMCALGeometry * EMCALGeometry() const ; 
  
//   //========== Methods to read something from file ==========
  void   Event(const Int_t event, const char * opt = "HSDRP") ;    
  void   Track(const Int_t itrack) ;
  
//   //-----------------now getter's data--------------------------------------
// //  AliEMCALCalibrationDB * CalibrationDB(){return  fcdb; }
// //  void ReadCalibrationDB(const char * name, const char * filename) ;

  //=========== Primaries ============
//   TTree *           TreeK(TString filename="") ; 
  TClonesArray *    Primaries(void)  ;
  TParticle * Primary(Int_t index) const ;
  Int_t       NPrimaries()const { return fNPrimaries; }
  TParticle * Secondary(const TParticle * p, const Int_t index=1) const ;  
  
//   //=========== Hits =================
//   TTree *               TreeH(TString filename="") ; 
  TClonesArray *  Hits(void)  ; 
  AliEMCALHit *    Hit(const Int_t index) { return dynamic_cast<AliEMCALHit*>(Hits()->At(index) );}
  TTree *         TreeH() const ; 
  
  //=========== SDigits ==============
  TClonesArray *      SDigits() ;  
  AliEMCALDigit *      SDigit(const Int_t index) { return static_cast<AliEMCALDigit *>(SDigits()->At(index)) ;} 
  TTree *             TreeS() const ; 
  AliEMCALSDigitizer * SDigitizer() ;  

  TString             GetSDigitsFileName() { return EmcalLoader()->GetSDigitsFileName() ; }  
  Int_t               LoadSDigits(Option_t* opt="") { return EmcalLoader()->LoadSDigits(opt) ; }
  Int_t               LoadSDigitizer(Option_t* opt=""){ return  EmcalLoader()->LoadSDigitizer(opt) ; }
  Int_t               WriteSDigits(Option_t* opt="") { return EmcalLoader()->WriteSDigits(opt) ; }
  Int_t               WriteSDigitizer(Option_t* opt=""){
    return  EmcalLoader()->WriteSDigitizer(opt) ; }
  
  //========== Digits ================
  TClonesArray * Digits() ;
  AliEMCALDigit * Digit(const Int_t index) { return static_cast<AliEMCALDigit *>(Digits()->At(index)) ;} 
  TTree *        TreeD() const ; 
  AliEMCALDigitizer * Digitizer() ;
  TString             GetDigitsFileName() { return EmcalLoader()->GetDigitsFileName() ; }  
  Int_t               LoadDigits(Option_t* opt="") { return EmcalLoader()->LoadDigits(opt) ; }
  Int_t               LoadDigitizer(Option_t* opt=""){
    return  EmcalLoader()->LoadDigitizer(opt) ; }
  Int_t               WriteDigits(Option_t* opt="") { return EmcalLoader()->WriteDigits(opt) ; }
  Int_t               WriteDigitizer(Option_t* opt=""){
    return  EmcalLoader()->WriteDigitizer(opt) ; }
  
  //========== RecPoints =============
  TObjArray *             PRERecPoints() ;
  AliEMCALRecPoint *      PRERecPoint(const Int_t index) { return static_cast<AliEMCALRecPoint *>(PRERecPoints()->At(index)) ;} 
  TObjArray *             ECARecPoints() ;
  AliEMCALTowerRecPoint * ECARecPoint(const Int_t index) { return static_cast<AliEMCALTowerRecPoint *>(ECARecPoints()->At(index)) ;} 
  TObjArray *             HCARecPoints() ;
  AliEMCALTowerRecPoint * HCARecPoint(const Int_t index) { return static_cast<AliEMCALTowerRecPoint *>(HCARecPoints()->At(index)) ;} 
  TTree *                 TreeR() const ;
  AliEMCALClusterizer *   Clusterizer()  ;
  TString                 GetRecPointsFileName() { return EmcalLoader()->GetRecPointsFileName() ; } 
  Int_t                   LoadRecPoints(Option_t* opt="") { return EmcalLoader()->LoadRecPoints(opt) ; }
  Int_t                   LoadClusterizer(Option_t* opt=""){
    return  EmcalLoader()->LoadClusterizer(opt) ; }
  Int_t                   WriteRecPoints(Option_t* opt="") { return EmcalLoader()->WriteRecPoints(opt) ; }
  Int_t                   WriteClusterizer(Option_t* opt=""){
    return  EmcalLoader()->WriteClusterizer(opt) ; }

  //========== TrackSegments   TClonesArray * TrackSegments(const char * name = 0) { 
  TClonesArray *           TrackSegments() ;
  AliEMCALTrackSegment *  TrackSegments(const Int_t index) { return static_cast<AliEMCALTrackSegment *>(TrackSegments()->At(index)) ;} 
  TTree *               TreeT() const ;
  AliEMCALTrackSegmentMaker * TrackSegmentMaker() ;
  TString               GetTracksFileName() { return EmcalLoader()->GetTracksFileName() ; } 
  Int_t                 LoadTracks(Option_t* opt="") { return EmcalLoader()->LoadTracks(opt) ; }
  Int_t                 LoadTrackSegementMaker(Option_t* opt=""){
    return  EmcalLoader()->LoadTrackSegmentMaker(opt) ; }
  Int_t                 WriteTracks(Option_t* opt="") { return EmcalLoader()->WriteTracks(opt) ; }
  Int_t                 WriteTrackSegmentMaker(Option_t* opt=""){
    return  EmcalLoader()->WriteTracker(opt) ; }
  //========== RecParticles ===========

  TClonesArray *         RecParticles() ;
  AliEMCALRecParticle *   RecPaticles(const Int_t index) { return static_cast<AliEMCALRecParticle *>(RecParticles()->At(index)) ;} 
  TTree *               TreeP() const ;
  AliEMCALPID * PID() ;
  TString               GetRecParticlesFileName() { return EmcalLoader()->GetRecParticlesFileName() ; } 
  Int_t                 LoadRecParticles(Option_t* opt="") { return EmcalLoader()->LoadRecParticles(opt) ; }
  Int_t                 LoadPID(Option_t* opt=""){
    return  EmcalLoader()->LoadPID(opt) ; }
  Int_t                 WriteRecParticles(Option_t* opt="") { return EmcalLoader()->WriteRecParticles(opt) ; }
  Int_t                 WritePID(Option_t* opt=""){
    return  EmcalLoader()->WritePID(opt) ; }


  void SetDebug(Int_t level) {fgDebug = level;} // Set debug level 
  void PostClusterizer(AliEMCALClusterizer * clu) 
    const{EmcalLoader()->PostClusterizer(clu) ; }
  void PostPID(AliEMCALPID * pid) 
    const{EmcalLoader()->PostPID(pid) ; }
  void PostTrackSegmentMaker(AliEMCALTrackSegmentMaker * tr) 
    const{EmcalLoader()->PostTrackSegmentMaker(tr) ; }
  void PostSDigitizer (AliEMCALSDigitizer * sdigitizer) 
    const {EmcalLoader()->PostSDigitizer(sdigitizer);}    
  void PostDigitizer (AliEMCALDigitizer * digitizer)    
    const {EmcalLoader()->PostDigitizer(dynamic_cast<AliDigitizer *>(digitizer));}

  TString Version() const  { return EmcalLoader()->GetTitle() ; } 
  AliEMCALLoader * EmcalLoader() const { return  fgEmcalLoader ; }
  
private:
  
  AliEMCALGetter(const char* headerFile,
		const char* version = AliConfig::fgkDefaultEventFolderName,
		Option_t * openingOption = "READ") ;

  Int_t ReadTreeD(void) ;
  Int_t ReadTreeH(void) ;
  Int_t ReadTreeR(void) ;
  Int_t ReadTreeT(void) ;
  Int_t ReadTreeS(void) ;
  Int_t ReadTreeP(void) ;


  void ReadPrimaries(void) ;

private:

//   static TFile * fgFile;           //! 

//  AliEMCALBeamTestEvent * fBTE ;           //! Header if BeamTest Event

  static Int_t          fgDebug ;             //! Debug level

  TString           fLoadingStatus ;     //! tells which trees are loaded
  Int_t             fNPrimaries ;        //! # of primaries  
  TClonesArray *    fPrimaries ;         //! list of lists of primaries

//  AliEMCALCalibrationDB * fcdb ;       //!
   
  static AliEMCALLoader * fgEmcalLoader ;
  static AliEMCALGetter * fgObjGetter; // pointer to the unique instance of the singleton 
  
  enum EDataTypes{kHits,kSDigits,kDigits,kRecPoints,kTracks,kNDataTypes};


  ClassDef(AliEMCALGetter,2)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliEMCALGETTER_H
