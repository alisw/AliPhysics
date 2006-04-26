#ifndef ALIPHOSGETTERLIGHT_H
#define ALIPHOSGETTERLIGHT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class to substitute AliPHOSGetter in "on flight" reconstruction (i.e. without 
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
  AliPHOSGetterLight(const AliPHOSGetterLight & obj) : AliPHOSGetter(obj) {
    // cpy ctor requested by Coding Convention 
    Fatal("cpy ctor", "not implemented") ;
  } 
  
  AliPHOSGetterLight & operator = (const AliPHOSGetterLight & ) {
    // assignement operator requested by coding convention, but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }

  virtual ~AliPHOSGetterLight() ; // dtor

  //method normally used for creation of this class
  static AliPHOSGetterLight * Instance(const char* /*headerFile*/,
				  const char* version = AliConfig::GetDefaultEventFolderName(),
				  Option_t * openingOption = "READ" ) ; 
  static AliPHOSGetterLight * Instance() ; 

  //=========== General information about run ==============
  
  virtual Int_t  MaxEvent() const {return 1 ;} //always "read" event 1 
  virtual Int_t  EventNumber() const {return 0; }  //always the same event 
  virtual Bool_t VersionExists(TString & /*opt*/) const {return kFALSE;} 
  virtual UShort_t EventPattern(void) const {return 0;}  //not needed in on-flight reconstruction
  virtual Float_t  BeamEnergy(void) const {return 10.;} //not needed in on-flight reconstruction
  
  //========== PHOSGeometry and PHOS ============= 
  //Dummy function not necessary for on-flight reconstruction, but has to be overloaded
  virtual AliPHOSGeometry * PHOSGeometry() const {return AliPHOSGeometry::GetInstance("IHEP","IHEP") ; } //Create if necessary geom
  
  //========== Methods to read something from file ==========
  virtual void   Event(Int_t /*event*/, const char * /*opt = "HSDRTP"*/){} //Use data already in memory    
 
  
  //-----------------now getter's data--------------------------------------
  
  //========== Digits ================
  virtual TClonesArray * Digits() const {return fDigits ; }
  virtual AliPHOSDigit * Digit(Int_t index) const { return static_cast<AliPHOSDigit *>(fDigits->At(index)) ;} 
  //  virtual AliPHOSDigitizer * Digitizer(){Error("Digitizer","Method not defined") ; return 0;}
  
  //========== RecPoints =============
  virtual TObjArray *           EmcRecPoints() const {return fEmcRecPoints ;}
  virtual AliPHOSEmcRecPoint *  EmcRecPoint(Int_t index) const { return static_cast<AliPHOSEmcRecPoint *>(fEmcRecPoints->At(index)) ;} 
  virtual TObjArray *           CpvRecPoints() const {return fCpvRecPoints ;} 
  virtual AliPHOSCpvRecPoint *  CpvRecPoint(Int_t index) const { return static_cast<AliPHOSCpvRecPoint *>(fCpvRecPoints->At(index)) ;} 
  virtual AliPHOSClusterizer * Clusterizer() { return fClusterizer;}
  
  //========== TrackSegments   TClonesArray * TrackSegments(const char * name = 0) { 
  virtual TClonesArray *        TrackSegments() const {return fTS ;} ;
  virtual AliPHOSTrackSegment * TrackSegment(Int_t index) const { return static_cast<AliPHOSTrackSegment *>(fTS->At(index)) ;} 
  virtual AliPHOSTrackSegmentMaker * TrackSegmentMaker(){ return fTSM ;}
  
  //========== RecParticles ===========
  virtual TClonesArray *        RecParticles() const { return fRP;} 
  virtual AliPHOSRecParticle *  RecParticle(Int_t index) const { return static_cast<AliPHOSRecParticle *>(fRP->At(index)) ;} 
  virtual AliPHOSPID *          PID(){return fPID ;} 

  //========== Raw ===========
  //  virtual Int_t ReadRaw(Int_t event) ; 

  virtual void PostClusterizer(AliPHOSClusterizer * clu)const{((AliPHOSGetterLight*)fgObjGetter)->fClusterizer = clu;} 
  virtual void PostPID(AliPHOSPID * pid)const{((AliPHOSGetterLight*)fgObjGetter)->fPID = pid ;} 
  virtual void PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tr)const{((AliPHOSGetterLight*)fgObjGetter)->fTSM = tr ;} 
  virtual TString Version() const  { return "OnFlight" ; } 
  virtual void Reset(){} 
  
  virtual AliESD * ESD() const { return 0 ; }

private:
  
  AliPHOSGetterLight(const char* headerFile,
		const char* version = AliConfig::GetDefaultEventFolderName(),
		Option_t * openingOption = "READ") ;

private :

  TClonesArray * fDigits ;        //Digits container for current event
  TObjArray    * fEmcRecPoints ;  //EmcRecPoints container for current event
  TObjArray    * fCpvRecPoints ;  //CPV RecPoints container for current event
  TClonesArray * fTS ;            //TrackSegments container for currect event
  TClonesArray * fRP ;            //Rec Particles conatiner for currect event

  //  AliPHOSCalibrationDB * fcdb ;   //Pointer to calibration database

  AliPHOSClusterizer       * fClusterizer ; //Pointer to clusterizer 
  AliPHOSTrackSegmentMaker * fTSM ;         //Pointer to TrackSegmentMaker
  AliPHOSPID               * fPID ;         //Pointer to PIDMaker

  //  Bool_t         fRawDigits ;    //Do we reconstruct raw digits

  ClassDef(AliPHOSGetterLight,1)  // Getter for \"on flyght\" reconstruction 

};

#endif // AliPHOSGETTERLIGHT_H
