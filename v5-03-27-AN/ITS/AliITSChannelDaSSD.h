#ifndef ALIITSCHANNELDASSD_H
#define ALIITSCHANNELDASSD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*                                                                        */
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides storage container ITS SSD channel callibration data
/// used by DA. 
///
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"


class AliITSChannelDaSSD : public TObject {
  public :
    AliITSChannelDaSSD();
    explicit AliITSChannelDaSSD(const UShort_t stripId);
    AliITSChannelDaSSD(const UShort_t stripId, const Long_t eventsnumber);
    AliITSChannelDaSSD(const AliITSChannelDaSSD& strip);
    AliITSChannelDaSSD& operator = (const AliITSChannelDaSSD& strip);
    virtual ~AliITSChannelDaSSD();
    
    UShort_t  GetStripId() const { return fStripId; }
    Long_t    GetEventsNumber() const { return fEventsNumber; }
    Short_t*  GetSignal()  const { return fSignal; }
    Short_t   GetSignal(const Long_t eventnumber)  const 
                           { return (eventnumber < fEventsNumber && fSignal) ? *(fSignal+eventnumber) : fgkDefaultSignal; }
    
    Float_t  GetPedestal() const { return fPedestal; }
    Float_t  GetNoise()    const { return fNoise;    }
    Float_t  GetNoiseCM()  const { return fNoiseCM;  }
    Long_t   GetOverflowNumber() const { return fNOverflowEv; }
    Bool_t   CheckNoise()   const { return (fNoise < fgkUndefinedValue) ? kTRUE : kFALSE; }
    Bool_t   CheckNoiseCM() const { return (fNoiseCM < fgkUndefinedValue) ? kTRUE : kFALSE; }

    Bool_t    SetEvenetsNumber(const Long_t eventsnumber);
    Bool_t    SetSignal(const Long_t eventnumber, const Short_t signal);
    void      SetStripId(const UShort_t stripId) { fStripId = stripId; }

    void      SetPedestal(Float_t pedestal) { fPedestal = pedestal; }
    void      SetNoise(Float_t noise)   { fNoise = noise; }
    void      SetNoiseCM(Float_t noise) { fNoiseCM = noise; }
    void      SetOverflowNumber(Long_t ovn) {fNOverflowEv = ovn; }
    void      ClearNoise() { if (fNoise < fgkUndefinedValue) fNoise = fgkUndefinedValue; }
    void      ClearNoiseCM() { if (fNoiseCM < fgkUndefinedValue) fNoiseCM = fgkUndefinedValue; }
    void      ClearPedestal() { if (fPedestal < fgkUndefinedValue) fPedestal = fgkUndefinedValue; }
    void      ClearSignal() { if (fSignal) for (Int_t i = 0; i < fEventsNumber; i++) 
                                             fSignal[i] = 0x100 * fgkDefaultSignal + fgkDefaultSignal; }
    void      DeleteSignal() { if (fSignal) { delete [] fSignal; fSignal = NULL; fEventsNumber = 0;} }

    static  Short_t  GetOverflowConst()  { return fgkSignalOverflow;  }
    static  Short_t  GetUnderflowConst() { return fgkSignalUnderflow; }
    static  Float_t  GetUndefinedValue() { return fgkUndefinedValue;  }
    static  Short_t  GetMaxStripIdConst(){ return fgkMaxStripId;      }

  protected :
    static const Short_t fgkMinStripId;         // minimum strip id
    static const Short_t fgkMaxStripId;         // maximum strip id

    static const Short_t  fgkSignalOverflow;     // ADC overflow value
    static const Short_t  fgkSignalUnderflow;    // ADC underflow value
    static const UShort_t fgkDefaultSignal;      // initialization value for fNoise, fPedestal, fSignal[i]
    static const Float_t  fgkUndefinedValue;     // constant value which indicates that fNoise, fPedestal is undefined
  
    UShort_t          fStripId;             //  channel = strip number within SSD module 0-1535
    Long_t            fEventsNumber;        //  number of events for fSignal memory allocation
    Short_t          *fSignal;              //! array of signal data
    
    Float_t           fPedestal;            //  pedestal
    Float_t           fNoise;               //  noise
    Float_t           fNoiseCM;             //  noise with CM correction
    Long_t            fNOverflowEv;         //  Number of events which exceed the pedestal calculation threshold
		      
    ClassDef(AliITSChannelDaSSD, 3)
};

#endif
