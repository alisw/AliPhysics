#ifndef ALIPHOSINDEXTOOBJECT_H
#define ALIPHOSINDEXTOOBJECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A singleton that returns objects from their indexes. 
//  Should be used on the analysis stage to avoid confusing between different
//  branches of reconstruction tree: e.g. reading RecPoints and TS made from 
//  another set of RecPoints.
// 
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    


// --- ROOT system ---
#include "TClonesArray.h"
#include <stdlib.h>

class TString ;
class TParticle ;

// --- Standard library ---

// --- AliRoot header files ---

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

class AliPHOSIndexToObject : public TObject {
  
 public:
  
  AliPHOSIndexToObject(){ 
    // ctor: this is a singleton, the ctor should never be called but cint needs it as public
    abort() ; 
  } 
  AliPHOSIndexToObject(const AliPHOSIndexToObject & obj) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    abort() ; 
  } 
  
  virtual ~AliPHOSIndexToObject(){
    // dtor
  }
  
  void  GetEvent(Int_t event) ; // reads event from file 
  Int_t GetEventNumber(){ return fEvent; }
  Int_t GetMaxEvent()   { return fMaxEvent;}
  
  static AliPHOSIndexToObject * GetInstance(const char* headerFile,const char* branch = "PHOSRP",
					    const char* branchTitle =0 ) ; 
  static AliPHOSIndexToObject * GetInstance() ; 
  
  AliPHOSDigit *        GimeSDigit(Int_t index)
    {if(fSDigits) return (AliPHOSDigit*)fSDigits->At(index); 
     return 0 ;} 
  Int_t                 GimeNSDigits() 
    {if(fSDigits) return fSDigits->GetEntriesFast(); 
     return 0 ;} 
  AliPHOSDigit *        GimeDigit(Int_t index) 
                          {if(fDigits) return (AliPHOSDigit*)fDigits->At(index); 
			  return 0 ;} 
  Int_t                 GimeNDigits() 
                          {if(fDigits) return fDigits->GetEntriesFast(); 
			  return 0 ;} 
  TParticle *           GimePrimary(Int_t index) ;
  Int_t                 GimeNPrimaries()const { return fNPrimaries; }
  AliPHOSRecParticle *  GimeRecParticle(Int_t index) 
                          { if(fRecParticles) return (AliPHOSRecParticle*)fRecParticles->At(index); 
			  return 0 ;} 
  Int_t                 GimeNRecParticles() 
                          { if(fRecParticles) return fRecParticles->GetEntriesFast(); 
			  return 0 ;} 
  AliPHOSEmcRecPoint *  GimeEmcRecPoint(Int_t index) 
                          { if(fEmcRecPoints) return (AliPHOSEmcRecPoint*)fEmcRecPoints->At(index);
			  return 0 ;} 
  Int_t                 GimeNEmcRecPoints() 
                          { if(fEmcRecPoints) return fEmcRecPoints->GetEntriesFast();
			    return 0 ;} 
  AliPHOSRecPoint *     GimeCpvRecPoint(Int_t index) 
                          { if(fCpvRecPoints) return (AliPHOSRecPoint*)fCpvRecPoints->At(index); 
			    return 0 ;}  
  Int_t                 GimeNCpvRecPoints() 
                          { if(fCpvRecPoints) return fCpvRecPoints->GetEntriesFast(); 
			    return 0 ;}  
  AliPHOSTrackSegment * GimeTrackSegment(Int_t index) 
                          { if(fTS) return (AliPHOSTrackSegment*)fTS->At(index); 
			    return 0 ;} 
  Int_t                 GimeNTrackSegments() 
                          { if(fTS) return fTS->GetEntriesFast(); 
			    return 0 ;} 
  AliPHOSSDigitizer*    GimeSDigitizer() 
                          { return fSDigitizer; } 
  AliPHOSDigitizer*     GimeDigitizer() 
                          { return fDigitizer; } 
  AliPHOSClusterizer*   GimeClusterizer() 
                          { return fClusterizer; } 
  AliPHOSTrackSegmentMaker* GimeTSMaker() 
                          { return fTSMaker; } 
  AliPHOSPID *          GimePID() 
                          { return fPID; } 
  
  AliPHOSIndexToObject & operator = (const AliPHOSIndexToObject & ) {
    // assignement operator requested by coding convention
    // but not needed
    abort() ;
    return *this ; 
  }
  
 private:
  
  AliPHOSIndexToObject(const char* headerFile,const char* branch = "RecParticles",
					      const char* branchTitle =0) ; 
  void DefineBranchTitles(const char* branch,const char* branchTitle) ;
  void ReadTreeR() ;
  void ReadTreeD() ;
  void ReadTreeS() ;
  void ReadPrimaries() ;

 private:

  TString        fHeaderFile ;    //!
  TString        fRPTitle ;       //!
  TString        fTSTitle ;       //!
  TString        fRecPointsTitle ;//!
  TString        fDigitsTitle ;   //!
  TString        fSDigitsTitle ;  //!

  Int_t          fEvent ;        
  Int_t          fMaxEvent ;      //! Number of events in currect headers file
  Int_t          fNPrimaries ;    //! # of primaries
  
  TClonesArray * fSDigits ;       //! list of Sdigits
  TClonesArray * fDigits ;        //! list of digits
  TObjArray *    fEmcRecPoints ;  //! emc rec points
  TObjArray *    fCpvRecPoints ;  //! CPV/PPSD rec points
  TClonesArray * fTS ;            //! TrackSegments
  TClonesArray * fRecParticles ;  //! RecParticles
  TObjArray *    fPrimaries ;     //! list of lists of primaries-for the case of mixing

  AliPHOSSDigitizer * fSDigitizer ;     //!
  AliPHOSDigitizer  * fDigitizer ;      //!
  AliPHOSClusterizer* fClusterizer ;    //!
  AliPHOSTrackSegmentMaker * fTSMaker ; //!
  AliPHOSPID *        fPID ;            //!

  static AliPHOSIndexToObject * fgObjGetter ; // pointer to the unique instance of the singleton 

  ClassDef(AliPHOSIndexToObject,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliPHOSINDEXTOOBJECT_H
