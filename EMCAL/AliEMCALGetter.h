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
class TParticle ;
class TTree ; 
class TGraph ;
class TF1 ;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h" 

// #include "AliRun.h"
#include "AliEMCALHit.h"   
#include "AliEMCALRecParticle.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h" 
#include "AliEMCALLoader.h" 
class AliEMCAL ; 
class AliEMCALGeometry ; 
class AliEMCALClusterizer ; 
class AliEMCALRecPoint ; 
class AliEMCALTrackSegmentMaker ;
class AliEMCALTrackSegment ; 
class AliEMCALPID ; 
class AliEMCALBeamTestEvent ;


class AliEMCALGetter : public TObject {
  
 public:  
  AliEMCALGetter() {    
    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
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
				  const char* version = AliConfig::GetDefaultEventFolderName(),
				  Option_t * openingOption = "READ" ) ; 
  static AliEMCALGetter * Instance() ; 
  void   OpenFile(const char* headerFile,
		  const char* version = AliConfig::GetDefaultEventFolderName(),
		  Option_t * openingOption = "READ" );

  static void Print() ; 
  
//   //=========== General information about run ==============
  Bool_t IsLoaded(TString tree) const { return fLoadingStatus.Contains(tree) ; } 
  void   SetLoaded(TString tree) { fLoadingStatus += tree ; } 

  Int_t  MaxEvent() const ; 
  Int_t  EventNumber() const ; 
  Bool_t VersionExists(TString & opt) const ; 
  UShort_t EventPattern(void) const ; 
  Float_t  BeamEnergy(void) const ;
  
//   //========== EMCALGeometry and EMCAL ============= 
  AliEMCAL *         EMCAL() const  ;  
  AliEMCALGeometry * EMCALGeometry() const ; 
  
//   //========== Methods to read something from file ==========
  void   Event(Int_t event, const char * opt = "HSDRP") ;    
  void   Track(Int_t itrack) ;
  
//   //-----------------now getter's data--------------------------------------
// //  AliEMCALCalibrationDB * CalibrationDB(){return  fcdb; }
// //  void ReadCalibrationDB(const char * name, const char * filename) ;

  //=========== Primaries ============
//   TTree *           TreeK(TString filename="") ; 
  TParticle * Primary(Int_t index) const ;
  Int_t       NPrimaries() const;
  TParticle * Secondary(const TParticle * p, Int_t index=1) const ;  
  
//   //=========== Hits =================
//   TTree *               TreeH(TString filename="") ; 
  TClonesArray *  Hits(void) const  ; 
  AliEMCALHit *    Hit(Int_t index) const { return dynamic_cast<AliEMCALHit*>(Hits()->At(index) );}
  TTree *         TreeH() const ; 
  
  //=========== SDigits ==============
  TClonesArray *      SDigits() const ;  
  AliEMCALDigit *      SDigit(Int_t index) const { return static_cast<AliEMCALDigit *>(SDigits()->At(index)) ;} 
  TTree *             TreeS() const ; 
  AliEMCALSDigitizer * SDigitizer() ;  

  TString             GetSDigitsFileName() const { return EmcalLoader()->GetSDigitsFileName() ; }  
  Int_t               LoadSDigits(Option_t* opt="") const { return EmcalLoader()->LoadSDigits(opt) ; }
  Int_t               LoadSDigitizer(Option_t* opt="") const { return  EmcalLoader()->LoadSDigitizer(opt) ; }
  Int_t               WriteSDigits(Option_t* opt="") const { return EmcalLoader()->WriteSDigits(opt) ; }
  Int_t               WriteSDigitizer(Option_t* opt="") const {
    return  EmcalLoader()->WriteSDigitizer(opt) ; }
  
  //========== Digits ================
  TClonesArray * Digits() const ;
  AliEMCALDigit * Digit(Int_t index) const { return static_cast<AliEMCALDigit *>(Digits()->At(index)) ;} 
  TTree *        TreeD() const ; 
  AliEMCALDigitizer * Digitizer() ;
  TString             GetDigitsFileName() const { return EmcalLoader()->GetDigitsFileName() ; }  
  Int_t               LoadDigits(Option_t* opt="") const { return EmcalLoader()->LoadDigits(opt) ; }
  Int_t               LoadDigitizer(Option_t* opt="") const {
    return  EmcalLoader()->LoadDigitizer(opt) ; }
  Int_t               WriteDigits(Option_t* opt="") const { return EmcalLoader()->WriteDigits(opt) ; }
  Int_t               WriteDigitizer(Option_t* opt="") const {
    return  EmcalLoader()->WriteDigitizer(opt) ; }
  
  //========== RecPoints =============
  TObjArray *             ECARecPoints() const;
  AliEMCALRecPoint * ECARecPoint(Int_t index) const{ return static_cast<AliEMCALRecPoint *>(ECARecPoints()->At(index)) ;}    
  TTree *                 TreeR() const ;
  AliEMCALClusterizer *   Clusterizer()  ;
  TString                 GetRecPointsFileName() const { return EmcalLoader()->GetRecPointsFileName() ; } 
  Int_t                   LoadRecPoints(Option_t* opt="") const { return EmcalLoader()->LoadRecPoints(opt) ; }
  Int_t                   LoadClusterizer(Option_t* opt="") const {
    return  EmcalLoader()->LoadClusterizer(opt) ; }
  Int_t                   WriteRecPoints(Option_t* opt="") const { return EmcalLoader()->WriteRecPoints(opt) ; }
  Int_t                   WriteClusterizer(Option_t* opt="") const {
    return  EmcalLoader()->WriteClusterizer(opt) ; }

  //========== TrackSegments   TClonesArray * TrackSegments(const char * name = 0) { 
  TClonesArray *           TrackSegments() const ;
  AliEMCALTrackSegment *  TrackSegments(Int_t index) const { return static_cast<AliEMCALTrackSegment *>(TrackSegments()->At(index)) ;} 
  TTree *               TreeT() const ;
  AliEMCALTrackSegmentMaker * TrackSegmentMaker() ;
  TString               GetTracksFileName() const { return EmcalLoader()->GetTracksFileName() ; } 
  Int_t                 LoadTracks(Option_t* opt="") const { return EmcalLoader()->LoadTracks(opt) ; }
  Int_t                 LoadTrackSegementMaker(Option_t* opt="") const {
    return  EmcalLoader()->LoadTrackSegmentMaker(opt) ; }
  Int_t                 WriteTracks(Option_t* opt="") const { return EmcalLoader()->WriteTracks(opt) ; }
  Int_t                 WriteTrackSegmentMaker(Option_t* opt="") const {
    return  EmcalLoader()->WriteTracker(opt) ; }
  //========== RecParticles ===========

  TClonesArray *         RecParticles() const ;
  AliEMCALRecParticle *  RecParticle(Int_t index) const { return static_cast<AliEMCALRecParticle *>(RecParticles()->At(index)) ;} 
  TTree *               TreeP() const ;
  AliEMCALPID * PID() ;
  TString               GetRecParticlesFileName() const { return EmcalLoader()->GetRecParticlesFileName() ; } 
  Int_t                 LoadRecParticles(Option_t* opt="") const { return EmcalLoader()->LoadRecParticles(opt) ; }
  Int_t                 LoadPID(Option_t* opt="") const {
    return  EmcalLoader()->LoadPID(opt) ; }
  Int_t                 WriteRecParticles(Option_t* opt="") const { return EmcalLoader()->WriteRecParticles(opt) ; }
  Int_t                 WritePID(Option_t* opt="") const {
    return  EmcalLoader()->WritePID(opt) ; }

  //========== Raw ===========
  Int_t ReadRaw(Int_t event) ; 

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
  AliEMCALLoader * EmcalLoader() const { return fgEmcalLoader; }
  void Reset() ;  

private:
  
  AliEMCALGetter(const char* headerFile,
		 const char* version = AliConfig::GetDefaultEventFolderName(),
		 Option_t * openingOption = "READ") ;


  Int_t ReadTreeD(void) ;
  Int_t ReadTreeH(void) ;
  Int_t ReadTreeR(void) ;
  Int_t ReadTreeT(void) ;
  Int_t ReadTreeS(void) ;
  Int_t ReadTreeP(void) ;


  void ReadPrimaries(void) ;
 
  void FitRaw(Bool_t lowGainFlag, TGraph * gLowGain, TGraph * gHighGain, TF1* signalF, Int_t & amp, Double_t & time) ; 

private:

//   static TFile * fgFile;           //! 

//  AliEMCALBeamTestEvent * fBTE ;           //! Header if BeamTest Event

  static Int_t          fgDebug ;             //! Debug level

  TString           fLoadingStatus ;     //! tells which trees are loaded
  static TString           fVersion;            //! stores the current folder name

//  AliEMCALCalibrationDB * fcdb ;       //!
   
  static AliEMCALLoader * fgEmcalLoader ; // pointer to EMCAL Loader
  static AliEMCALGetter * fgObjGetter; // pointer to the unique instance of the singleton 
  
  enum EDataTypes{kHits,kSDigits,kDigits,kRecPoints,kTracks,kNDataTypes};


  ClassDef(AliEMCALGetter,5)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliEMCALGETTER_H
