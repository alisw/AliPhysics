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
#include <stdlib.h>
#include <iostream.h>

class TString ;
class TParticle ;
class TTask ;
class TFolder ; 

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOS ;
class AliPHOSGeometry ;
class AliPHOSHit ;
class AliPHOSDigit ;
class AliPHOSDigitizer ;
class AliPHOSSDigitizer ;
class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegment ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSRecParticle ;
class AliPHOSPID ;

class AliPHOSGetter : public TObject {
  
 public:
  
  AliPHOSGetter(){ 
    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    cerr << "ERROR: AliPHOGetter is a singleton default ctor not callable" << endl ;
    abort() ; 
  } 
  AliPHOSGetter(const AliPHOSGetter & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    abort() ; 
  } 
  
  virtual ~AliPHOSGetter(){
    // dtor
  }
  
  void Post(const char * file, const char * opt, const char * name = 0, const Int_t event=-1) const ;  
  void  Event(Int_t event) ; // reads event from file 
  //     Int_t EventNumber(){ return fEvent; }
  //     Int_t MaxEvent()   { return fMaxEvent;}
  static AliPHOSGetter * GetInstance(const char* headerFile,
				     const char* branchTitle = "No Name" ) ; 
  static AliPHOSGetter * GetInstance() ; 
  const AliPHOS * PHOS() const ;  
  const AliPHOSGeometry * PHOSGeometry() const  ; 
 
  // Alarms
  TFolder * Alarms(const char * name = 0) const { return (TFolder*)(ReturnO("Alarms", name)) ; }

  // Hits
  TClonesArray * Hits(const char * name = 0) const { return (TClonesArray*)(ReturnO("Hits", name)) ; }
  const AliPHOSHit * Hit(Int_t index)
    {if( Hits() ) return (AliPHOSHit*)(Hits()->At(index)); 
      return 0 ;} 
  const Int_t NHits() 
    {if( Hits() ) return Hits()->GetEntriesFast(); 
     return 0 ;} 
  // SDigits
  TClonesArray * SDigits(const char * name = 0, const char * file=0) const { 
    return (TClonesArray*)(ReturnO("SDigits", name, file)) ; }
  const AliPHOSDigit *        SDigit(Int_t index)
    {if( SDigits() ) return (AliPHOSDigit*)(SDigits()->At(index)); 
      return 0 ;} 
  const Int_t                 NSDigits() 
    {if( SDigits() ) return SDigits()->GetEntriesFast(); 
     return 0 ;} 
  AliPHOSSDigitizer *   SDigitizer(const char * name =0) const { return ((AliPHOSSDigitizer*)(ReturnT("SDigitizer", name))) ; }
  // Digits
  TClonesArray *        Digits(const char * name = 0) const { return (TClonesArray*)(ReturnO("Digits", name)) ; }
  const AliPHOSDigit *        Digit(Int_t index) 
                          {if( Digits() ) return (AliPHOSDigit*)Digits()->At(index); 
			  return 0 ;} 
  const Int_t                 NDigits() 
                          {if( Digits() ) return Digits()->GetEntriesFast(); 
			  return 0 ;} 
  AliPHOSDigitizer *    Digitizer(const char * name =0) const { return (AliPHOSDigitizer*)(ReturnT("Digitizer", name)) ; }
  // RecPoints
  TObjArray * EmcRecPoints(const char * name = 0, const char * file=0) const { 
    return (TObjArray*)(ReturnO("EmcRecPoints", name, file)) ; }
  const AliPHOSEmcRecPoint *  EmcRecPoint(Int_t index) 
    { if(EmcRecPoints()) return (AliPHOSEmcRecPoint*)EmcRecPoints()->At(index); return 0 ;} 
  const Int_t NEmcRecPoints() 
    { if(EmcRecPoints()) return EmcRecPoints()->GetEntriesFast(); return 0 ;} 
  TObjArray * CpvRecPoints(const char * name = 0, const char * file=0) const { 
    return (TObjArray*)(ReturnO("CpvRecPoints", name, file)) ; }
  const AliPHOSRecPoint * CpvRecPoint(Int_t index) 
    { if(CpvRecPoints()) return (AliPHOSRecPoint*)CpvRecPoints()->At(index); return 0 ;}  
  const Int_t NCpvRecPoints() 
    { if(CpvRecPoints()) return CpvRecPoints()->GetEntriesFast(); return 0 ;}  
  AliPHOSClusterizer * Clusterizer (const char * name =0) const 
    { return (AliPHOSClusterizer*)(ReturnT("Clusterizer", name)) ; }
  // TrackSegments
  TClonesArray * TrackSegments(const char * name = 0, const char * file=0) const { 
    return (TClonesArray*)(ReturnO("TrackSegments", name, file)) ; }
  const AliPHOSTrackSegment * TrackSegment(Int_t index) 
    { if(TrackSegments()) return (AliPHOSTrackSegment*)TrackSegments()->At(index); return 0 ;} 
  const Int_t NTrackSegments() 
    { if(TrackSegments()) return TrackSegments()->GetEntriesFast(); return 0 ;} 
  AliPHOSTrackSegmentMaker * TrackSegmentMaker (const char * name =0) const 
    { return (AliPHOSTrackSegmentMaker*)(ReturnT("TrackSegmentMaker", name)) ; }
  // RecParticles
  TClonesArray * RecParticles(const char * name = 0, const char * file=0) const { 
    return (TClonesArray*)(ReturnO("RecParticles", name, file)) ; }
  const AliPHOSRecParticle * RecParticle(Int_t index) 
    { if(RecParticles()) return (AliPHOSRecParticle*)RecParticles()->At(index); return 0 ;} 
  const Int_t NRecParticles() 
    { if(RecParticles()) return RecParticles()->GetEntriesFast(); return 0 ;} 
  AliPHOSPID * PID(const char * name =0) const 
    { return (AliPHOSPID*)(ReturnT("PID", name)) ; }
  // Primaries
  const TParticle *           Primary(Int_t index) const ;
  const Int_t                 NPrimaries()const { return fNPrimaries; }


  AliPHOSGetter & operator = (const AliPHOSGetter & ) {
    // assignement operator requested by coding convention
    // but not needed
    abort() ;
    return *this ; 
  }
  
 private:

  AliPHOSGetter(const char* headerFile, const char* branchTitle =0) ; 
  void CreateWhiteBoard() const ; 
  const TObject * ReturnO(TString what, TString name=0, TString file=0) const ; 
  const TTask * ReturnT(TString what,TString name=0) const ; 
  void DefineBranchTitles(char* branch, char* branchTitle) ;
  void ReadTreeD() ;
  void ReadTreeH() ;
  void ReadTreeQA() ;
  void ReadTreeR() ;
  void ReadTreeS() ;
  void ReadPrimaries() ;

 private:

  TString        fHeaderFile ;    //!
  TString        fTrackSegmentsTitle ;//!
  TString        fRecPointsTitle ;//!
  TString        fRecParticlesTitle ;//!
  TString        fDigitsTitle ;   //!
  TString        fSDigitsTitle ;  //!

  Int_t          fNPrimaries ;    //! # of primaries
  
  TObjArray *    fPrimaries ;     //! list of lists of primaries-for the case of mixing

  static AliPHOSGetter * fgObjGetter ; // pointer to the unique instance of the singleton 

  ClassDef(AliPHOSGetter,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliPHOSGETTER_H
