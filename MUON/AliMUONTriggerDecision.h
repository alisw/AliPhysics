/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
#ifndef ALIMUONTRIGGERDECISION_H
#define ALIMUONTRIGGERDECISION_H
////////////////////////////////////////////////
//  MUON Trigger Decision Class               //
////////////////////////////////////////////////
#include "TObject.h"

class AliMUONHitMapA1;
class TF1;
class TClonesArray;
class AliMUONSegmentation;
class AliMUONResponse;

class AliMUONTriggerDecision :
public TObject {
 public:
  AliMUONTriggerDecision(Int_t iprint);         // constructor
  ~AliMUONTriggerDecision();                  // destructor
  
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

  // print-debug
  void PrintBitPatXInput(Int_t icirc);
  void PrintBitPatYInput(Int_t icirc);
  void PrintLocalOutput(Int_t minDevStrip[5], Int_t minDev[5], 
			Int_t coordY[5]);

  // return member data information
  Int_t GetITrigger(Int_t icirc);
  Int_t GetStripX11(Int_t icirc);
  Int_t GetDev(Int_t icirc);
  Int_t GetStripY11(Int_t icirc);
  void GetLutOutput(Int_t icirc, Int_t lpt[2], Int_t hpt[2], Int_t apt[2]);
  void GetGlobalTrigger(Int_t singlePlus[3], Int_t singleMinus[3], 
			Int_t singleUndef[3], Int_t pairUnlike[3], 
			Int_t pairLike[3]);  
  
// Add a new Local Trigger
  // virtual void AddLocalTrigger(const AliMUONLocalTrigger);
//  Return pointer to Local Triggers
  //  TClonesArray* LocalTriggers(){return fLocalTriggers;}

  ClassDef(AliMUONTriggerDecision,1) // Trigger Decision class

    protected:     
  Int_t fDebug;               // print option     

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

  // ???
  //  TClonesArray* fLocalTriggers;   // Local Triggers
  // Int_t fNLocalTriggers;          // Number of Local Triggers

};
#endif







