#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H
 
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 

#include <TObject.h>

/* $Id$ */

/////////////////////////////////////////////////////////////
//  Base settings for the ITS response classes.            //  
//  The data member of this class are static and set once  //
//  for all the modules.                                   //    
///////////////////////////////////////////////////////////// 

class AliITSresponseSDD : public TObject {
  public:

    AliITSresponseSDD();
    virtual ~AliITSresponseSDD(){};


    virtual void SetTimeOffset(Float_t to){fTimeOffset = to;}
    virtual Float_t GetTimeOffset()const {return fTimeOffset;}
    static Float_t DefaultTimeOffset() {return fgkTimeOffsetDefault;}

    virtual void SetADC2keV(Float_t conv){fADC2keV=conv;}
    virtual Float_t GetADC2keV()const {return fADC2keV;}
    static Float_t DefaulttADC2keV() {return fgkADC2keVDefault;}

    static Float_t GetCarlosRXClockPeriod() {return fgkCarlosRXClockPeriod;}
 

 protected:

    static const Float_t fgkTimeOffsetDefault; // default for fTimeOffset
    static const Float_t fgkADC2keVDefault; // default for fADC2keV
    static const Float_t fgkCarlosRXClockPeriod;  // clock period for CarlosRX

    Float_t  fTimeOffset;     // Time offset due to electronic delays 
    Float_t  fADC2keV;        // Conversion factor from ADC to keV

 private:

   AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
   AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.

    ClassDef(AliITSresponseSDD,15) 
    
    };
#endif
