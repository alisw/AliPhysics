#ifndef ALIPHOSGETTERLIGHT_H
#define ALIPHOSGETTERLIGHT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class to mask AliPHOSGetter in "on flight" reconstruction (i.e. without 
//  writing to file and creation full ALICE data structure)      
//                  
//*-- Author: D.Peressounko (RRC KI)


// --- ROOT system ---
class TClonesArray ;
class TObjArray ;

// --- Standard library ---
 
// --- AliRoot header files ---
#include "AliPHOSGetter.h" 
#include "AliPHOSGeometry.h"
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSPID ;

class AliPHOSGetterLight : public AliPHOSGetter {

public:
  AliPHOSGetterLight() ;          // ctor

  virtual ~AliPHOSGetterLight() ; // dtor

  static AliPHOSGetterLight * Instance(const char* /*headerFile*/,
				  const char* version = AliConfig::GetDefaultEventFolderName(),
				  Option_t * openingOption = "READ" ) ; 
  static AliPHOSGetterLight * Instance() ; 

  //=========== General information about run ==============
  virtual Bool_t IsLoaded(TString /*tree*/) const { Error("IsLoaded","NotDefined"); return kFALSE ; } 
  virtual void   SetLoaded(TString /*tree*/) { Error("SetLoaded","NotDefined"); } 
  
  virtual Int_t  MaxEvent() const {return 1 ;} //always "read" event 1 
  virtual Int_t  EventNumber() const {return 0; }  //always the same event 
  virtual Bool_t VersionExists(TString & /*opt*/) const {return kFALSE;} 
  virtual UShort_t EventPattern(void) const {return 0;} 
  virtual Float_t  BeamEnergy(void) const {return 10.;}
  
  //========== PHOSGeometry and PHOS ============= 
  virtual AliPHOS *         PHOS() const { Error("PHOS()","NotDefined"); return 0 ;}  
  virtual AliPHOSGeometry * PHOSGeometry() const {AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance("IHEP","IHEP") ; return geom ;} 
  
  //========== Methods to read something from file ==========
  virtual void   Event(Int_t /*event*/, const char * /*opt = "HSDRTP"*/){} //Use data already in memory    
  virtual void   Track(Int_t /*itrack*/) { Error("Track()","NotDefined");}
 
  
  //-----------------now getter's data--------------------------------------
  virtual AliPHOSCalibrationDB * CalibrationDB(){return  fcdb; }
  virtual void SetCalibrationDB(AliPHOSCalibrationDB * cdb) {fcdb = cdb ;}
  
  //=========== Primaries ============
  virtual TClonesArray *    Primaries(void) {Error("Primaries()","NotDefined"); return 0;}
  virtual TParticle * Primary(Int_t /*index*/) const {Error("Primary()","NotDefined"); return 0;}
  virtual Int_t       NPrimaries()const {Error("NPrimaries()","NotDefined"); return 0; }
  virtual TParticle * Secondary(const TParticle * /*p*/, Int_t /*index*/) const {Error("Secondary()","NotDefined"); return 0; }  
  
  //=========== Hits =================
  virtual TClonesArray *  Hits(void) {Error("Hits()","Method not defined") ; return 0;}   
  virtual AliPHOSHit *    Hit(Int_t /*index*/) {Error("Hit()","Method not defined"); return 0 ;}
  virtual TTree *         TreeH() const {Error("TreeH","Method not defined") ; return 0;}  
  
  //=========== SDigits ==============
  virtual TClonesArray *      SDigits(){Error("Sdigits","Method not defined") ; return 0;}   
  virtual AliPHOSDigit *      SDigit(Int_t /*index*/) { Error("SDigit","Method not defined") ;return 0 ;} 
  virtual TTree *             TreeS() const {Error("TresS","Method not defined") ; return 0 ;} 
  virtual AliPHOSSDigitizer * SDigitizer() {Error("SDigitizer","Method not defined") ; return 0 ;}
  
  virtual TString             GetSDigitsFileName() const { Error("OnFlight","Method not defined") ; return "" ; }  
  virtual Int_t               LoadSDigits(Option_t* /*opt=""*/) const { Error("LoadSDigits","Method not defined") ;return 0 ; }
  virtual Int_t               LoadSDigitizer(Option_t* /*opt=""*/) const { Error("LoadSdigitizer","Method not defined") ; return  0 ; }
  virtual Int_t               WriteSDigits(Option_t* /*opt=""*/) const  { Error("WriteSDigits","Method not defined") ; return 0 ; }
  virtual Int_t               WriteSDigitizer(Option_t* /*opt=""*/) const {Error("WriteSDigitizer","Method not defined") ; return  0 ; }
  
  //========== Digits ================
  virtual TClonesArray * Digits(){return fDigits ; }
  virtual AliPHOSDigit * Digit(Int_t index) { return static_cast<AliPHOSDigit *>(fDigits->At(index)) ;} 
  virtual TTree *        TreeD() const {Error("TreeD","Method not defined") ; return 0;}  
  virtual AliPHOSDigitizer * Digitizer(){Error("Digitizer","Method not defined") ; return 0;}
  virtual TString             GetDigitsFileName() const { Error("GetDigitsFileName","Method not defined") ; return "" ; }  
  virtual Int_t               LoadDigits(Option_t* /*opt=""*/) const {Error("LoadDigits","Method not defined") ; return 0 ; }
  virtual Int_t               LoadDigitizer(Option_t* /*opt=""*/) const {Error("LoadDigitizer","Method not defined") ; return 0;}
  virtual Int_t               WriteDigits(Option_t* /*opt=""*/) const { Error("WriteDigits","Method not defined") ; return 0; }
  virtual Int_t               WriteDigitizer(Option_t* /*opt=""*/) const {Error("WriteDigitizer","Method not defined"); return 0 ;}

  //Methods to distinguish raw and simulated digits
  virtual Bool_t              IsRawDigits(void) const {return fRawDigits;}
  virtual void                SetRawDigits(Bool_t isRaw = kTRUE){fRawDigits = isRaw;}
  
  //========== RecPoints =============
  virtual TObjArray *           EmcRecPoints(){return fEmcRecPoints ;}
  virtual AliPHOSEmcRecPoint *  EmcRecPoint(Int_t index) { return static_cast<AliPHOSEmcRecPoint *>(fEmcRecPoints->At(index)) ;} 
  virtual TObjArray *           CpvRecPoints(){return fCpvRecPoints ;} 
  virtual AliPHOSCpvRecPoint *  CpvRecPoint(Int_t index) { return static_cast<AliPHOSCpvRecPoint *>(fCpvRecPoints->At(index)) ;} 
  virtual TTree *               TreeR() const {Error("TreeR","Method not defined") ; return 0;}
  virtual AliPHOSClusterizer * Clusterizer() { return fClusterizer;}
  virtual TString               GetRecPointsFileName() const { Error("RecPointsFileName","Method not defined") ;return "" ; } 
  virtual Int_t                 LoadRecPoints(Option_t* /*opt=""*/) const {Error("LoadRecPoints","Method not defined") ; return 0; }
  virtual Int_t                 LoadClusterizer(Option_t* /*opt=""*/) const {Error("LoadClusterizer","Method not defined") ; return 0 ;} 
  virtual Int_t                 WriteRecPoints(Option_t* /*opt=""*/) const {Error("WriteRecPoints","Method not defined") ; return 0; }
  virtual Int_t                 WriteClusterizer(Option_t* /*opt=""*/) const {Error("WriteClusterizer","Method not defined"); return 0 ;}
  
  //========== TrackSegments   TClonesArray * TrackSegments(const char * name = 0) { 
  virtual TClonesArray *        TrackSegments(){return fTS ;} ;
  virtual AliPHOSTrackSegment * TrackSegment(Int_t index) { return static_cast<AliPHOSTrackSegment *>(fTS->At(index)) ;} 
  virtual TTree *               TreeT() const {Error("TreeT","Method not defined") ; return 0 ;} 
  virtual AliPHOSTrackSegmentMaker * TrackSegmentMaker(){ return fTSM ;}
  virtual TString               GetTracksFileName() const { Error("GetTSFileName","Method not defiled") ; return "" ; } 
  virtual Int_t                 LoadTracks(Option_t* /*opt=""*/) const { Error("LoadTracks","Method not defined") ;return 0; }
  virtual Int_t                 LoadTrackSegementMaker(Option_t* /*opt=""*/) const {Error("LoadTSMaker","Method noe defined") ; return 0 ;} 
  virtual Int_t                 WriteTracks(Option_t* /*opt=""*/) const { Error("WriteTracks","Method not defined") ; return 0; }
  virtual Int_t                 WriteTrackSegmentMaker(Option_t* /*opt=""*/) const {Error("WriteTSM","Method not defined") ; return 0 ;}
  
  //========== RecParticles ===========
  virtual TClonesArray *        RecParticles(){ return fRP;} 
  virtual AliPHOSRecParticle *  RecParticle(Int_t index) { return static_cast<AliPHOSRecParticle *>(fRP->At(index)) ;} 
  virtual TTree *               TreeP() const {Error("TreeP","Method net defined"); return 0 ;} 
  virtual AliPHOSPID *          PID(){return fPID ;} 
  virtual TString               GetRecParticlesFileName() const { Error("GetRPFileName","Method not defined") ; return "" ; } 
  virtual Int_t                 LoadRecParticles(Option_t* /*opt=""*/) const { Error("LoadRP","Method not defined") ;return 0 ; }
  virtual Int_t                 LoadPID(Option_t* /*opt=""*/) const {Error("LoadPID","Method not defined") ;return 0 ; }
  virtual Int_t                 WriteRecParticles(Option_t* /*opt=""*/) const { Error("WriteRP","Method not defined"); return 0 ; }
  virtual Int_t                 WritePID(Option_t* /*opt=""*/) const {Error("WritePID","Method not defined"); return 0 ; } 

  //========== Raw ===========
  //  virtual Int_t ReadRaw(Int_t event) ; 

  //  virtual void SetDebug(Int_t level) {fgDebug = level;} // Set debug level 
  virtual void PostClusterizer(AliPHOSClusterizer * clu)const{((AliPHOSGetterLight*)fgObjGetter)->fClusterizer = clu;} 
  virtual void PostPID(AliPHOSPID * pid)const{((AliPHOSGetterLight*)fgObjGetter)->fPID = pid ;} 
  virtual void PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tr)const{((AliPHOSGetterLight*)fgObjGetter)->fTSM = tr ;} 
  //virtual void PostSDigitizer (AliPHOSSDigitizer * sdigitizer) 
  // const {PhosLoader()->PostSDigitizer(sdigitizer);}    
  //virtual void PostDigitizer (AliPHOSDigitizer * digitizer)    
  // const {PhosLoader()->PostDigitizer(dynamic_cast<AliDigitizer *>(digitizer));}
  
  virtual TString Version() const  { return "OnFlight" ; } 
  virtual AliPHOSLoader * PhosLoader() const { Error("PhosLoader","Method not defined") ;return 0 ; }
  virtual void Reset(){} 
  
  virtual AliESD * ESD() const { return 0 ; }

private:
  
  AliPHOSGetterLight(const char* headerFile,
		const char* version = AliConfig::GetDefaultEventFolderName(),
		Option_t * openingOption = "READ") ;

private :

  TClonesArray * fDigits ;
  TObjArray    * fEmcRecPoints ;
  TObjArray    * fCpvRecPoints ;
  TClonesArray * fTS ;
  TClonesArray * fRP ;

  AliPHOSCalibrationDB * fcdb ;

  AliPHOSClusterizer       * fClusterizer ; 
  AliPHOSTrackSegmentMaker * fTSM ;
  AliPHOSPID               * fPID ;

  Bool_t         fRawDigits ;

  ClassDef(AliPHOSGetterLight,1)  // Getter for \"on flyght\" reconstruction 

};

#endif // AliPHOSGETTERLIGHT_H
