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
#include "TObject.h"  
class TParticle ;
class TTree ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h" 
#include "AliPHOSLoader.h" 
#include "AliPHOSHit.h" 
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h" 
#include "AliPHOSCpvRecPoint.h" 
#include "AliPHOSTrackSegment.h" 
#include "AliPHOSRecParticle.h" 
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
class AliPHOS ;  
class AliPHOSGeometry ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;  
class AliPHOSPID ; 
class AliPHOSBeamTestEvent ;
class AliESD ; 

class AliPHOSGetter : public TObject {
  
public:  
  AliPHOSGetter(){    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    Fatal("ctor", "AliPHOSGetter is a singleton default ctor not callable") ;
  } 
  AliPHOSGetter(const AliPHOSGetter & obj) : TObject(obj) {
    // cpy ctor requested by Coding Convention 
    Fatal("cpy ctor", "not implemented") ;
  } 
  
  AliPHOSGetter & operator = (const AliPHOSGetter & ) {
    // assignement operator requested by coding convention, but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  virtual ~AliPHOSGetter() ; 
  
  //=========== Instantiators ================
  static AliPHOSGetter * Instance(const char* headerFile,
				  const char* version = AliConfig::GetDefaultEventFolderName(),
				  Option_t * openingOption = "READ" ) ; 
  static AliPHOSGetter * Instance() ; 
  
  static void Print() ; 
  
  //=========== General information about run ==============
  Bool_t IsLoaded(TString tree) const { return fLoadingStatus.Contains(tree) ; } 
  void   SetLoaded(TString tree) { fLoadingStatus += tree ; } 
  
  Int_t  MaxEvent() const ; 
  Int_t  EventNumber() const ; 
  Bool_t VersionExists(TString & opt) const ; 
  UShort_t EventPattern(void) const ; 
  Float_t  BeamEnergy(void) const ;
  
  //========== PHOSGeometry and PHOS ============= 
  AliPHOS *         PHOS() const  ;  
  AliPHOSGeometry * PHOSGeometry() const ; 
  
  //========== Methods to read something from file ==========
  void   Event(Int_t event, const char * opt = "HSDRTP") ;    
  void   Track(Int_t itrack) ;
  
  //-----------------now getter's data--------------------------------------
  //  AliPHOSCalibrationDB * CalibrationDB(){return  fcdb; }
  //  void ReadCalibrationDB(const char * name, const char * filename) ;
  
  //=========== Primaries ============
  TClonesArray *    Primaries(void) ;
  TParticle * Primary(Int_t index) const ;
  Int_t       NPrimaries()const { return fNPrimaries; }
  TParticle * Secondary(const TParticle * p, Int_t index=1) const ;  
  
  //=========== Hits =================
  TClonesArray *  Hits(void)  ; 
  AliPHOSHit *    Hit(Int_t index) { return dynamic_cast<AliPHOSHit*>(Hits()->At(index) );}
  TTree *         TreeH() const ; 
  
  //=========== SDigits ==============
  TClonesArray *      SDigits() ;  
  AliPHOSDigit *      SDigit(Int_t index) { return static_cast<AliPHOSDigit *>(SDigits()->At(index)) ;} 
  TTree *             TreeS() const ; 
  AliPHOSSDigitizer * SDigitizer() ;  
  
  TString             GetSDigitsFileName() const { return PhosLoader()->GetSDigitsFileName() ; }  
  Int_t               LoadSDigits(Option_t* opt="") const { return PhosLoader()->LoadSDigits(opt) ; }
  Int_t               LoadSDigitizer(Option_t* opt="") const { return  PhosLoader()->LoadSDigitizer(opt) ; }
  Int_t               WriteSDigits(Option_t* opt="") const  { return PhosLoader()->WriteSDigits(opt) ; }
  Int_t               WriteSDigitizer(Option_t* opt="") const {
    return  PhosLoader()->WriteSDigitizer(opt) ; }
  
  //========== Digits ================
  TClonesArray * Digits() ;
  AliPHOSDigit * Digit(Int_t index) { return static_cast<AliPHOSDigit *>(Digits()->At(index)) ;} 
  TTree *        TreeD() const ; 
  AliPHOSDigitizer * Digitizer() ;
  TString             GetDigitsFileName() const { return PhosLoader()->GetDigitsFileName() ; }  
  Int_t               LoadDigits(Option_t* opt="") const { return PhosLoader()->LoadDigits(opt) ; }
  Int_t               LoadDigitizer(Option_t* opt="") const {
    return  PhosLoader()->LoadDigitizer(opt) ; }
  Int_t               WriteDigits(Option_t* opt="") const { return PhosLoader()->WriteDigits(opt) ; }
  Int_t               WriteDigitizer(Option_t* opt="") const {
    return  PhosLoader()->WriteDigitizer(opt) ; }
  
  //========== RecPoints =============
  TObjArray *           EmcRecPoints() ;
  AliPHOSEmcRecPoint *  EmcRecPoint(Int_t index) { return static_cast<AliPHOSEmcRecPoint *>(EmcRecPoints()->At(index)) ;} 
  TObjArray *           CpvRecPoints() ; 
  AliPHOSCpvRecPoint *  CpvRecPoint(Int_t index) { return static_cast<AliPHOSCpvRecPoint *>(CpvRecPoints()->At(index)) ;} 
  TTree *               TreeR() const ;
  AliPHOSClusterizer * Clusterizer() ;
  TString               GetRecPointsFileName() const { return PhosLoader()->GetRecPointsFileName() ; } 
  Int_t                 LoadRecPoints(Option_t* opt="") const { return PhosLoader()->LoadRecPoints(opt) ; }
  Int_t                 LoadClusterizer(Option_t* opt="") const {
    return  PhosLoader()->LoadClusterizer(opt) ; }
  Int_t                 WriteRecPoints(Option_t* opt="") const { return PhosLoader()->WriteRecPoints(opt) ; }
  Int_t                 WriteClusterizer(Option_t* opt="") const {
    return  PhosLoader()->WriteClusterizer(opt) ; }
  
  //========== TrackSegments   TClonesArray * TrackSegments(const char * name = 0) { 
  TClonesArray *           TrackSegments() ;
  AliPHOSTrackSegment *  TrackSegment(Int_t index) { return static_cast<AliPHOSTrackSegment *>(TrackSegments()->At(index)) ;} 
  TTree *               TreeT() const ;
  AliPHOSTrackSegmentMaker * TrackSegmentMaker() ;
  TString               GetTracksFileName() const { return PhosLoader()->GetTracksFileName() ; } 
  Int_t                 LoadTracks(Option_t* opt="") const { return PhosLoader()->LoadTracks(opt) ; }
  Int_t                 LoadTrackSegementMaker(Option_t* opt="") const {
    return  PhosLoader()->LoadTrackSegmentMaker(opt) ; }
  Int_t                 WriteTracks(Option_t* opt="") const { return PhosLoader()->WriteTracks(opt) ; }
  Int_t                 WriteTrackSegmentMaker(Option_t* opt="") const {
    return  PhosLoader()->WriteTracker(opt) ; }
  
  //========== RecParticles ===========
  TClonesArray *         RecParticles() ;
  AliPHOSRecParticle *   RecParticle(Int_t index) { return static_cast<AliPHOSRecParticle *>(RecParticles()->At(index)) ;} 
  TTree *               TreeP() const ;
  AliPHOSPID * PID() ;
  TString               GetRecParticlesFileName() const { return PhosLoader()->GetRecParticlesFileName() ; } 
  Int_t                 LoadRecParticles(Option_t* opt="") const { return PhosLoader()->LoadRecParticles(opt) ; }
  Int_t                 LoadPID(Option_t* opt="") const {
    return  PhosLoader()->LoadPID(opt) ; }
  Int_t                 WriteRecParticles(Option_t* opt="") const { return PhosLoader()->WriteRecParticles(opt) ; }
  Int_t                 WritePID(Option_t* opt="") const {
    return  PhosLoader()->WritePID(opt) ; }
  
  void SetDebug(Int_t level) {fgDebug = level;} // Set debug level 
  void PostClusterizer(AliPHOSClusterizer * clu) 
    const{PhosLoader()->PostClusterizer(clu) ; }
  void PostPID(AliPHOSPID * pid) 
    const{PhosLoader()->PostPID(pid) ; }
  void PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tr) 
    const{PhosLoader()->PostTrackSegmentMaker(tr) ; }
  void PostSDigitizer (AliPHOSSDigitizer * sdigitizer) 
    const {PhosLoader()->PostSDigitizer(sdigitizer);}    
  void PostDigitizer (AliPHOSDigitizer * digitizer)    
    const {PhosLoader()->PostDigitizer(dynamic_cast<AliDigitizer *>(digitizer));}
  
  TString Version() const  { return PhosLoader()->GetTitle() ; } 
  AliPHOSLoader * PhosLoader() const { return  fgPhosLoader ; }
  void Reset() ;
  
  AliESD * ESD() const { return fESD ; }
  
private:
  
  AliPHOSGetter(const char* headerFile,
		const char* version = AliConfig::GetDefaultEventFolderName(),
		Option_t * openingOption = "READ") ;
  
  Int_t ReadTreeD(void) ;
  Int_t ReadTreeH(void) ;
  Int_t ReadTreeR(void) ;
  Int_t ReadTreeT(void) ;
  Int_t ReadTreeS(void) ;
  Int_t ReadTreeP(void) ;

  Int_t ReadTreeE(Int_t event) ;    
  Bool_t OpenESDFile() ;
  void ReadPrimaries(void) ;
  
private:
  
  AliPHOSBeamTestEvent * fBTE ;           //! Header if BeamTest Event
  static Int_t          fgDebug ;             //! Debug level
  
  TString           fLoadingStatus ;     //! tells which trees are loaded
  Int_t             fNPrimaries ;        //! # of primaries  
  TClonesArray *    fPrimaries ;         //! list of lists of primaries
  TFile *           fESDFile ;           //! ESD file
  TString           fESDFileName ;       //! ESD File Name
  AliESD *          fESD ;               //! ESD object
  TTree *           fESDTree ;           //! ESD Tree

  //  AliPHOSCalibrationDB * fcdb ;       //!
  
  static AliPHOSLoader * fgPhosLoader ; // the loader for the NewIO
  static AliPHOSGetter * fgObjGetter; // pointer to the unique instance of the singleton 
  
  enum EDataTypes{kHits,kSDigits,kDigits,kRecPoints,kTracks,kNDataTypes};
  
  
  ClassDef(AliPHOSGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 
    
    };

#endif // AliPHOSGETTER_H
