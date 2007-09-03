#ifndef ALIITSCHANNELDASSD_H
#define ALIITSCHANNELDASSD_H


#include <iostream>
#include "TMath.h"
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
    Short_t*  GetSignal(const Long_t eventnumber)  const 
                                { return (eventnumber < fEventsNumber && fSignal) ? (fSignal+eventnumber) : NULL; }
    
    Float_t  GetPedestal() const { return fPedestal; }
    Float_t  GetNoise()    const { return fNoise; }
    Bool_t   CheckNoise()  const { return (fNoise < fgkUndefinedValue) ? kTRUE : kFALSE; }

    UChar_t  GetFeeOffset() const { return (fPedestal < fgkUndefinedValue) ? (static_cast<UChar_t>(TMath::Nint(fPedestal))) : 0; }
    UChar_t  GetFeeZs()     const { return (fNoise < fgkUndefinedValue) ? 
                                                    (static_cast<UChar_t>(fZsThresholdFactor * TMath::Nint(fNoise))) : 0;      }
    Bool_t    SetEvenetsNumber(const Long_t eventsnumber);
    Bool_t    SetSignal(const Long_t eventnumber, const Short_t signal);
    void      SetStripId(const UShort_t stripId) { fStripId = stripId; }

    void      SetPedestal(Float_t pedestal) { fPedestal = pedestal; }
    void      SetNoise(Float_t noise) { fNoise = noise; }
    void      ClearNoise() { if (fNoise < fgkUndefinedValue) fNoise = fgkUndefinedValue; }
    void      ClearPedestal() { if (fPedestal < fgkUndefinedValue) fPedestal = fgkUndefinedValue; }
    void      ClearSignal() { if (fSignal) for (Int_t i = 0; i < fEventsNumber; i++) 
                                             fSignal[i] = 0x100 * fgkDefaultSignal + fgkDefaultSignal; }
    void      DeleteSignal() { if (fSignal) { delete [] fSignal; fSignal = NULL; fEventsNumber = 0;} }

    static  Short_t  GetOverflowConst()  { return fgkSignalOverflow;  }
    static  Short_t  GetUnderflowConst() { return fgkSignalUnderflow; }
    static  Float_t  GetUndefinedValue() { return fgkUndefinedValue;  }

  protected :
  
    static const UShort_t fgkMinStripId = 0;
    static const UShort_t fgkMaxStripId = 1535;

    static const Short_t  fgkSignalOverflow  = 2047;
    static const Short_t  fgkSignalUnderflow = 2048;
    static const UShort_t fgkDefaultSignal   = 0x7F;
    static const Float_t  fgkUndefinedValue;
 
    UShort_t          fStripId;             //  channel = strip number within SSD module 0-1535
    Long_t            fEventsNumber;        //  number of events for fSignal memory allocation
    Short_t          *fSignal;              //! 
    
    Float_t           fPedestal;
    Float_t           fNoise;
 
    Float_t           fZsThresholdFactor;   //  factor for zero suppression threshold

    ClassDef(AliITSChannelDaSSD, 1)
};

#endif
