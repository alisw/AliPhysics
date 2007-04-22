#ifndef ALITRDTRIGGERHLT_H
#define ALITRDTRIGGERHLT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "AliLog.h"
#include "AliTRDgeometry.h"
#include "AliTRDtrigger.h"

/* class TTree; */
/* class TClonesArray; */
/* class TObjArray; */

class AliRunLoader;
class AliRawReader;
class AliRawReaderMemory;

//class AliTRDtrigger : public TNamed {
class AliTRDtriggerHLT : public AliTRDtrigger 
{
 public:  
  AliTRDtriggerHLT();
  AliTRDtriggerHLT(const Text_t* name, const Text_t* title);
  AliTRDtriggerHLT(const AliTRDtriggerHLT &p);   
  virtual         ~AliTRDtriggerHLT();
  AliTRDtriggerHLT   &operator=(const AliTRDtriggerHLT &p); 

  virtual Bool_t   TestTracklet(Int_t det, Int_t row, Int_t seed, Int_t n);
  virtual Bool_t   TreeTracklets(Int_t idet);
  virtual Bool_t   IsTreeOwner() const {return fTreeCreatedHere;}
  TTree *          GetTrackletTree() {return fTrackletTree;}
  virtual Bool_t   ResetTree();
  virtual void     Copy(TObject &p) const;
  virtual void     Init();
  virtual Bool_t   ReadDigits(AliRawReaderMemory* rawReader);
  virtual Bool_t   MakeTracklets(Bool_t makeTracks = kFALSE);
/*   virtual void     MakeTracks(Int_t det); */

 protected:

  Bool_t            fTreeCreatedHere; //flag indicating that AliTRDtriggerHLT has created the cluster tree

 private:

  //functions should not be used on HLT:
  void     SetRunLoader(AliRunLoader *rl)                
  { 
    fRunLoader = rl;    
    AliInfo("Not to be used!");    
  }

  //Bool_t   WriteTracklets(Int_t det);

  Bool_t   ReadDigits()
  {
    AliInfo("Not to be used!");
    return kFALSE;
  }

  Bool_t   ReadTracklets(AliRunLoader *rl)
  {
    AliInfo("Not to be used!");
    return kFALSE;
  }

  Bool_t   Open(const Char_t *name, Int_t nEvent = 0)
  {
    AliInfo("Not to be used!");
    return kFALSE;
  }

  virtual Bool_t   ReadDigits(AliRawReader* rawReader)
  {
    AliInfo("Not to be used!");
    return kFALSE;
  }
      
  ClassDef(AliTRDtriggerHLT,1)                                    //  TRD trigger class

};

#endif
