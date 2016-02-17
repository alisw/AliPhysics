#ifndef ALIVCALOTRIGGER_H
#define ALIVCALOTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//
//   Virtual class to access calorimeter 
//   (EMCAL, PHOS, PMD, FMD) trigger data
//   Author: Salvatore Aiola
//
//-------------------------------------------------------------------------

#include <TNamed.h>

class AliVCaloTrigger : public TNamed 
{
public:

  AliVCaloTrigger(): TNamed() {;}
  AliVCaloTrigger(const char* name, const char* title) : TNamed(name, title) {;}
  AliVCaloTrigger(const AliVCaloTrigger& ctrig);
  virtual ~AliVCaloTrigger() {;}
  AliVCaloTrigger& operator=(const AliVCaloTrigger& ctrig);
	
  virtual Bool_t       IsEmpty()                                             = 0;
  virtual void         Reset()                                               = 0;
  virtual void         Allocate(Int_t /*size*/)                              = 0;
  virtual void         DeAllocate()                                          = 0;
  
  virtual Bool_t       Add(Int_t /*col*/, Int_t /*row*/, 
                           Float_t /*amp*/, Float_t /*time*/, 
                           Int_t* /*trgtimes*/, Int_t /*ntrgtimes*/, 
                           Int_t /*trgts*/, Int_t /*trgbits*/)               = 0;
  
  virtual Bool_t       Add(Int_t /*col*/, Int_t /*row*/, 
                           Float_t /*amp*/, Float_t /*time*/, 
                           Int_t* /*trgtimes*/, Int_t /*ntrgtimes*/, 
                           Int_t /*trgts*/, Int_t /*subr*/, Int_t /*trgbit*/)= 0;

  
  virtual void         SetL1Threshold(Int_t /*i*/, Int_t /*thr*/)            = 0;
  virtual void         SetL1Threshold(Int_t /*i*/, Int_t /*j*/, Int_t /*th*/)= 0;
  
  virtual void         SetL1V0(const Int_t* /*v*/)                           = 0;
  virtual void         SetL1V0(Int_t /*i*/, const Int_t* /*v*/)              = 0;
  
  virtual void         SetL1FrameMask(Int_t /*m*/)                           = 0;
  virtual void         SetL1FrameMask(Int_t /*i*/, Int_t /*m*/)              = 0;

  virtual void         GetPosition(Int_t& /*col*/, Int_t& /*row*/)    const  = 0;
  virtual void         GetAmplitude(Float_t& /*amp*/)                 const  = 0;
  virtual void         GetTime(Float_t& /*time*/)                     const  = 0;
  
  virtual void         GetTriggerBits(Int_t& /*bits*/)                const  = 0;
  virtual void         GetNL0Times(Int_t& /*ntimes*/)                 const  = 0;
  virtual void         GetL0Times(Int_t* /*times*/)                   const  = 0;
  virtual Int_t        GetEntries()                                   const  = 0;

  virtual void         GetL1TimeSum(Int_t& /*timesum*/)               const  = 0;
  virtual Int_t        GetL1TimeSum()                                 const  = 0;
  
  virtual void         GetL1SubRegion(  Int_t& /*subreg*/)            const  = 0;
  virtual Int_t        GetL1SubRegion()                               const  = 0;
  
  virtual Int_t        GetL1Threshold(Int_t /*i*/)                    const  = 0;
  virtual Int_t        GetL1Threshold(Int_t /*i*/, Int_t /*j*/)       const  = 0;
  
  virtual Int_t        GetL1V0(Int_t /*i*/)                           const  = 0;
  virtual Int_t        GetL1V0(Int_t /*i*/, Int_t /*j*/)              const  = 0;
  
  virtual Int_t        GetL1FrameMask()                               const  = 0;
  virtual Int_t        GetL1FrameMask(Int_t /*i*/)                    const  = 0;
 
  virtual Int_t        GetMedian(Int_t /*i*/)                         const  = 0;
  
  virtual Int_t        GetTriggerBitWord()                            const  = 0;
  virtual void         GetTriggerBitWord(Int_t& /*bw*/ )              const  = 0;

  virtual Bool_t       Next()                                                = 0;
  virtual void         Copy(TObject& obj)                             const     ;

  virtual void         Print(const Option_t* /*opt*/)                 const  = 0;
  
private:

  ClassDef(AliVCaloTrigger, 0)
};
#endif //ALIVCALOTRIGGER_H

