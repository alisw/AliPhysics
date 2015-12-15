#ifndef ALIESDCALOTRIGGER_H
#define ALIESDCALOTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 


Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliVCaloTrigger.h"

class TArrayI;

class AliESDCaloTrigger : public AliVCaloTrigger 
{
public:
	         AliESDCaloTrigger();
	         AliESDCaloTrigger(const AliESDCaloTrigger& ctrig);
	virtual ~AliESDCaloTrigger();
	
	AliESDCaloTrigger& operator=(const AliESDCaloTrigger& ctrig);
	
	Bool_t  IsEmpty() {return (fNEntries == 0);}

	virtual void Reset() {fCurrent = -1;}

	void    Allocate(Int_t size);
	void    DeAllocate(        ); 
	
        Bool_t  Add(Int_t col, Int_t row, Float_t amp, Float_t time, Int_t trgtimes[], Int_t ntrgtimes, Int_t trgts, Int_t trgbits);
        Bool_t  Add(Int_t col, Int_t row, Float_t amp, Float_t time, Int_t trgtimes[], Int_t ntrgtimes, Int_t trgts, Int_t subra, Int_t trgbits);
	
	void    SetL1Threshold(Int_t i, Int_t thr) {fL1Threshold[i] = thr;}
        void    SetL1Threshold(Int_t i, Int_t j, Int_t thr) {if (i) fL1DCALThreshold[j] = thr; else fL1Threshold[j] = thr;}
  
	void    SetL1V0(const Int_t* v) {for (int i = 0; i < 2; i++) fL1V0[i] = v[i];}
        void    SetL1V0(Int_t i, const Int_t* v) {
          if (i) {for (int j = 0; j < 2; j++) fL1DCALV0[j] = v[j];} else {for (int j = 0; j < 2; j++) fL1V0[j] = v[j];}
        }
 
        void    SetL1FrameMask(Int_t m) {fL1FrameMask = m;}
        void    SetL1FrameMask(Int_t i, Int_t m) {if (i) fL1DCALFrameMask = m; else fL1FrameMask = m;}
  
	void    SetTriggerBitWord(Int_t w) {fTriggerBitWord = w;}
	void    SetMedian(Int_t i, Int_t m) {fMedian[i] = m;}

	void    GetPosition(       Int_t& col, Int_t& row           ) const;
	void    GetAmplitude(      Float_t& amp                     ) const;
	void    GetTime(           Float_t& time                    ) const;
	
	void    GetTriggerBits(    Int_t& bits                      ) const;
	void    GetNL0Times(       Int_t& ntimes                    ) const;
	void    GetL0Times(        Int_t  times[]                   ) const;
	Int_t   GetEntries(                                         ) const {return fNEntries;}

        void    GetL1TimeSum(      Int_t& timesum                   ) const;
        Int_t   GetL1TimeSum(                                       ) const;
  
        void    GetL1SubRegion(    Int_t& subreg                    ) const;
        Int_t   GetL1SubRegion(                                     ) const;
  
	Int_t   GetL1Threshold(    Int_t  i                         ) const {return fL1Threshold[i];}
	Int_t   GetL1Threshold(    Int_t  i, Int_t j                ) const {return ((i)?fL1DCALThreshold[j]:fL1Threshold[j]);}
  
	Int_t   GetL1V0(           Int_t  i                         ) const {return fL1V0[i];}
        Int_t   GetL1V0(           Int_t  i, Int_t  j               ) const {return ((i)?fL1DCALV0[j]:fL1V0[j]);}
  
        Int_t   GetL1FrameMask(                                     ) const {return fL1FrameMask;}
        Int_t   GetL1FrameMask(    Int_t  i                         ) const {return ((i)?fL1DCALFrameMask:fL1FrameMask);}
  
        Int_t   GetMedian(         Int_t  i                         ) const {return fMedian[i];}
  
        Int_t   GetTriggerBitWord(                                  ) const {return fTriggerBitWord;}
        void    GetTriggerBitWord( Int_t& bw                        ) const {bw = fTriggerBitWord;}
        
	virtual Bool_t Next();

	virtual void Copy(TObject& obj) const;
	
	virtual void Print(const Option_t* opt) const;
	
private:

        Int_t    fNEntries;
        Int_t    fCurrent;

	Int_t*   fColumn;             // [fNEntries]
	Int_t*   fRow;                // [fNEntries]
	Float_t* fAmplitude;          // [fNEntries]
	Float_t* fTime;               // [fNEntries]
	Int_t*   fNL0Times;           // [fNEntries]
	TArrayI* fL0Times;            //
	Int_t*   fL1TimeSum;          // [fNEntries]
	Int_t*   fTriggerBits;        // [fNEntries]
	
	Int_t    fL1Threshold[4];     // L1 thresholds from raw data
	Int_t    fL1V0[2];            // L1 threshold components
	Int_t    fL1FrameMask;        // Validation flag for L1 data
  
        Int_t    fL1DCALThreshold[4]; // L1 thresholds from raw data
        Int_t*   fL1SubRegion;        // [fNEntries]
        Int_t    fL1DCALFrameMask;    // Validation flag for L1 data
        Int_t    fMedian[2];          // Background median
        Int_t    fTriggerBitWord;     // Trigger bit word
        Int_t    fL1DCALV0[2];        // L1 threshold components
	
	ClassDef(AliESDCaloTrigger, 8)
};
#endif

