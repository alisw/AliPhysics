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
class TString ;
class TParticle ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSDigit ;
class AliPHOSDigitizer ;
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

  static AliPHOSIndexToObject * GetInstance(char* headerFile,char* branch = "PHOSRP",
					    char* branchTitle =0 ) ; 
  static AliPHOSIndexToObject * GetInstance() ; 
  
  AliPHOSDigit *        GimeDigit(Int_t index) 
                          {if(fDigits) return (AliPHOSDigit*)fDigits->At(index); 
			   else        return 0 ;} 
  TParticle *           GimePrimary(Int_t index) ;
  AliPHOSRecParticle *  GimeRecParticle(Int_t index) 
                          { if(fRecParticles) return (AliPHOSRecParticle*)fRecParticles->At(index); 
			    else return 0 ;} 
  AliPHOSEmcRecPoint *  GimeEmcRecPoint(Int_t index) 
                          { if(fEmcRecPoints) return (AliPHOSEmcRecPoint*)fEmcRecPoints->At(index);
			    else return 0 ;} 
  AliPHOSRecPoint *     GimeCpvRecPoint(Int_t index) 
                          { if(fCpvRecPoints) return (AliPHOSRecPoint*)fCpvRecPoints->At(index); 
			    else return 0 ;}  
  AliPHOSTrackSegment * GimeTrackSegment(Int_t index) 
                          { if(fTS) return (AliPHOSTrackSegment*)fTS->At(index); 
			    else return 0 ;} 
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
  
  AliPHOSIndexToObject(char* headerFile,char* branch = "RecParticles",
					    char* branchTitle =0) ; 
  Bool_t ReadRecParticles(char * branchTitle) ;
  Bool_t ReadTS(char * branchTitle) ;
  Bool_t ReadRecPoints(char * branchTitle) ;
  Bool_t ReadDigits(char * branchTitle) ;
  Bool_t ReadPrimaries() ;


 private:

  TString        fHeaderFile ; 

  Int_t          fEvent ;         //! Current event
  Int_t          fMaxEvent ;      //! Number of events in currect headers file

  TClonesArray * fDigits ;        //! list of digits
  TObjArray *    fEmcRecPoints ;  //! emc rec points
  TObjArray *    fCpvRecPoints ;  //! CPV/PPSD rec points
  TClonesArray * fTS ;            //! TrackSegments
  TClonesArray * fRecParticles ;  //! RecParticles
  TObjArray *    fPrimaries ;     //! list of lists of primaries-for the case of mixing

  AliPHOSDigitizer  * fDigitizer ;
  AliPHOSClusterizer* fClusterizer ;
  AliPHOSTrackSegmentMaker * fTSMaker ;
  AliPHOSPID *        fPID ;

  static AliPHOSIndexToObject * fgObjGetter ; // pointer to the unique instance of the singleton 

  ClassDef(AliPHOSIndexToObject,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliPHOSINDEXTOOBJECT_H
