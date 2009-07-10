#ifndef ALIPHOSRAWFITTERV0_H
#define ALIPHOSRAWFITTERV0_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts the signal parameters (energy, time, quality)
// from ALTRO samples. Energy is in ADC counts, time is in time bin units.
// A coarse algorithm is applied.

class TArrayI;
class AliPHOSCalibData;

class AliPHOSRawFitterv0 : public TObject 
{

public:

  AliPHOSRawFitterv0();
  AliPHOSRawFitterv0(const AliPHOSRawFitterv0& rawFitterv0);
  AliPHOSRawFitterv0& operator = (const AliPHOSRawFitterv0& rawFitterv0);
  virtual ~AliPHOSRawFitterv0();

  void SubtractPedestals(Bool_t subtract) {fPedSubtract  = subtract;}
  void SetAmpOffset     (Int_t extPed=5)  {fAmpOffset    = extPed ;}
  void SetAmpThreshold  (Int_t thr=5)     {fAmpThreshold = thr ;}
  void SetNBunches(const Int_t nBunches) { fNBunches = nBunches; }
  void SetChannelGeo(const Int_t module, const Int_t cellX,
		     const Int_t cellZ,  const Int_t caloFlag);

  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);
  Double_t GetEnergy()        const { return fEnergy;      }
  Double_t GetTime()          const { return fTime;        }
  Double_t GetSignalQuality() const { return fQuality;     }
  Double_t GetPedestalRMS()   const { return fPedestalRMS; }
  Int_t    GetModule()        const { return fModule;      }
  Int_t    GetCellX()         const { return fCellX;       }
  Int_t    GetCellZ()         const { return fCellZ;       }
  Int_t    GetCaloFlag()      const { return fCaloFlag;    }
  Bool_t   IsOverflow()       const { return fOverflow;    }

  void SetCalibData(AliPHOSCalibData * cdata){ fCalibData=cdata ;}

protected:   
  
  Int_t    fModule;         // PHOS module number
  Int_t    fCellX;          // cell number along X-axis
  Int_t    fCellZ;          // cell number along Z-axis
  Int_t    fCaloFlag;       // 0=LG, 1=HG, 2=TRU
  Int_t    fNBunches;       // number of bunches in a signal
  Bool_t   fPedSubtract;    // pedestals subtraction (kTRUE="yes")
  Double_t fEnergy;         // "digit" energy
  Double_t fTime;           // "digit" time
  Double_t fQuality ;       // sample quality
  Double_t fPedestalRMS;    // calciulated RMS of pedestal (non-ZS runs)
  Int_t    fAmpOffset ;     // pedestal offset from ALTRO chips
  Int_t    fAmpThreshold ;  // zero suppression threshold from ALTRO chips
  Bool_t   fOverflow ;      // kTRUE is the signal overflows
  AliPHOSCalibData * fCalibData ;   //! Calibration database if avalable

  ClassDef(AliPHOSRawFitterv0,2)
};

#endif
