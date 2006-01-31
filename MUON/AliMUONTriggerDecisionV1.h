/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#ifndef AliMUONTriggerDecisionV1V1_H
#define AliMUONTriggerDecisionV1V1_H

/// \ingroup sim
/// \class AliMUONTriggerDecisionV1
/// \brief MUON trigger decision class 
///
/// interim solution (but needed so far)

#ifndef ROOT_TTask
#include "TTask.h"
#endif

#ifndef ROOT_TArrayI
#include "TArrayI.h"
#endif

class AliMUONData;
class AliMUON;
class TObjArray;

class AliMUONTriggerDecisionV1 : public TTask 
{
 public:
  AliMUONTriggerDecisionV1(AliMUONData* data); // constructor
  AliMUONTriggerDecisionV1(); // constructor
  virtual ~AliMUONTriggerDecisionV1();                  // destructor
  
  void Trigger();
  void ResetBit();
  void SetBit();
  void SetBitUpDownY();

  void TrigX(Int_t ch1q[16], Int_t ch2q[16], Int_t ch3q[32], Int_t ch4q[32], 
	     Int_t coinc44, Int_t minDevStrip[5], Int_t minDev[5]);
  void Sort2x5(Int_t dev1[6], Int_t dev2[6],
	       Int_t minDev[6], Int_t &dev1GTdev2);
  void TrigY(Int_t y1[16], Int_t y2[16], Int_t y3[16], Int_t y4[16],
	     Int_t y3u[16], Int_t y3d[16], Int_t y4u[16], Int_t y4d[16],
	     Int_t x2m, Int_t x2ud, Int_t orMud[2], Int_t resetMid, 
	     Int_t coinc44, Int_t coordY[5]);
  void LocalTrigger(Int_t icirc, Int_t minDevStrip[5], 
		    Int_t minDev[5], Int_t coordY[5], 
		    Int_t &iTrigger);    
  void GlobalTrigger();

  void Exec(Option_t* opt="");

  // print-debug
  void PrintBitPatXInput(Int_t icirc);
  void PrintBitPatYInput(Int_t icirc);
  void PrintLocalOutput(Int_t minDevStrip[5], Int_t minDev[5], 
			Int_t coordY[5]);

  // return member data information
  Int_t GetITrigger(Int_t icirc) const;
  Int_t GetStripX11(Int_t icirc) const;
  Int_t GetDev(Int_t icirc) const;
  Int_t GetStripY11(Int_t icirc) const;
  void GetLutOutput(Int_t icirc, Int_t lpt[2], Int_t hpt[2], Int_t apt[2]) const;
  void GetGlobalTrigger(Int_t singlePlus[3], Int_t singleMinus[3], 
			Int_t singleUndef[3], Int_t pairUnlike[3], 
			Int_t pairLike[3]) const;  


  ClassDef(AliMUONTriggerDecisionV1,1) // Trigger Decision class


protected:

  AliMUONTriggerDecisionV1(const AliMUONTriggerDecisionV1& rhs);
  AliMUONTriggerDecisionV1& operator=(const AliMUONTriggerDecisionV1& rhs);

  void ClearDigitNumbers();

  void DigitFiredCircuit(
                Int_t circuit, Int_t cathode,
                Int_t chamber, Int_t digit
        );

  // Global Trigger information [0] : Low pt, [1] : High pt, [2] : All pt 
  Int_t fGlobalSinglePlus[3];  // tot num of single plus
  Int_t fGlobalSingleMinus[3]; // tot num of single minus
  Int_t fGlobalSingleUndef[3]; // tot num of single undefined
  Int_t fGlobalPairUnlike[3];  // tot num of unlike-sign pairs
  Int_t fGlobalPairLike[3];    // tot num of like-sign pairs

  // Local Trigger information
  Int_t fTrigger[234];  // fTrigger = 0 : no trigger, 1 : trigger
  Int_t fStripX11[234];  // X strip in MC11 which triggers
  Int_t fDev[234];       // deviation which triggers
  Int_t fStripY11[234];  // Y strip in MC11 which triggers
  Int_t fLutLpt[234][2]; // Local Trigger info Low pt
  Int_t fLutHpt[234][2]; // Local Trigger info High pt
  Int_t fLutApt[234][2]; // Local Trigger info All pt

  // bit pattern
  Int_t fXbit11[234][16]; // bit pattern XM11
  Int_t fXbit12[234][16]; // bit pattern XM12
  Int_t fXbit21[234][32]; // bit pattern XM21
  Int_t fXbit22[234][32]; // bit pattern XM22
  Int_t fYbit11[234][16]; // bit pattern YM11
  Int_t fYbit12[234][16]; // bit pattern YM12
  Int_t fYbit21[234][16]; // bit pattern YM21
  Int_t fYbit22[234][16]; // bit pattern YM22

  Int_t fYbit21U[234][16]; // bit pattern YM21 Up
  Int_t fYbit22U[234][16]; // bit pattern YM22 Up
  Int_t fYbit21D[234][16]; // bit pattern YM21 Down
  Int_t fYbit22D[234][16]; // bit pattern YM22 Down

  TArrayI fDigitNumbers[234];  //! The digit number that fired a circuit.
  
  TObjArray*     fTriggerCircuit;     //! List of Trigger Circuit
  AliMUONData*   fMUONData;           //! Data container for MUON subsystem 
  AliMUON*       fMUON;               //! pointer to MUON  
};
#endif







