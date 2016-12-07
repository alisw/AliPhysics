#ifndef ALITOFTRIGGER_H
#define ALITOFTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//                                             //
//      TOF Trigger Detector Class             //
//                                             //
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"
#include "AliLog.h"
#include "TTree.h"

class AliTOFrawData;
class AliTOFTriggerMask;

class AliTOFTrigger : public AliTriggerDetector
{
 public:
  AliTOFTrigger();  // constructor
  AliTOFTrigger(Int_t HighMultTh, Int_t ppMBTh, Int_t MultiMuonTh, Int_t UPTh, Float_t deltaminpsi, Float_t deltamaxpsi, Float_t deltaminro, Float_t deltamaxro, Int_t stripWindow,Float_t startTimeWindow=0.0,Float_t widthTimeWindow=25.); //constructor with parameters
  virtual ~AliTOFTrigger();  // destructor
  virtual void    CreateInputs();
  virtual void    Trigger();
  Int_t   GetHighMultTh() const {return fHighMultTh;}
  Int_t   GetppMBTh() const {return fppMBTh;}
  Int_t   GetMultiMuonTh() const {return fMultiMuonTh;}
  Int_t   GetUPTh() const {return fUPTh;}
  Float_t Getdeltaminpsi() const {return fdeltaminpsi;}
  Float_t Getdeltamaxpsi() const {return fdeltamaxpsi;}
  Float_t Getdeltaminro() const {return fdeltaminro;}
  Float_t Getdeltamaxro() const {return fdeltamaxro;}
  Int_t  GetstripWindow() const {return fstripWindow;}

  static void LoadActiveMask(); // Load active channel trigger mask
  void GetMapMatrix(Bool_t map[][24]) const;
  void GetMap(Bool_t **map) const;
  //void PrintMap(); // to be checked because of warning problems
  void GetTRDmapMatrix(Bool_t map[][8]) const;
  void GetTRDmap(Bool_t **map) const;
  Bool_t GetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,Int_t iTDC, Int_t iCH);
  Bool_t GetBit(Int_t *detind);
  void SetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,Int_t iTDC, Int_t iCH);
  void SetBit(Int_t *detind);
  void ResetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,Int_t iTDC, Int_t iCH);
  void ResetBit(Int_t *detind);


  void   SetHighMultTh(Int_t HighMultTh){fHighMultTh = HighMultTh;}
  void   SetppMBTh(Int_t ppMBTh){fppMBTh = ppMBTh;}
  void   SetMultiMuonTh(Int_t MultiMuonTh){fMultiMuonTh = MultiMuonTh;}
  void   SetUPTh(Int_t UPTh){fUPTh = UPTh;}
  void   Setdeltaminpsi(Float_t deltaminpsi){fdeltaminpsi = deltaminpsi;}
  void   Setdeltamaxpsi(Float_t deltamaxpsi){fdeltamaxpsi = deltamaxpsi;}
  void   Setdeltaminro(Float_t deltaminro){fdeltaminro = deltaminro;}
  void   Setdeltamaxro(Float_t deltamaxro){fdeltamaxro = deltamaxro;}
  void   SetstripWindow(Int_t stripWindow){fstripWindow = stripWindow;}

  Bool_t Return(Int_t i){if(i==0) return fSel1;
			 else if(i==1) return fSel2;
                         else if(i==2) return fSel3;
                         else if(i==3) return fSel4;
			 else { AliWarning(Form(" Index out of range: %d not in [0,3]",i)); return kFALSE; }
			};

  Float_t GetStartTimeWindow() const {return fStartTimeHit;}; // in ns
  Float_t GetTimeWidthWindow() const {return fTimeWidthTrigger;}; // in ns
  void    SetStartTimeWindow(Float_t val) {fStartTimeHit = val;}; // in ns
  void    SetTimeWidthWindow(Float_t val) {fTimeWidthTrigger = val;}; // in ns
  
  Int_t GetNumberOfCrateOn(){return fNCrateOn;}; 
  Int_t GetNumberOfMaxipadOn(){return fNMaxipadOn;}; 
  Int_t GetNumberOfMaxipadOnAll(){return fNMaxipadOnAll;}; 
  Bool_t *GetLTMarray(){return fLTMarray;};
  void   CreateCTTMMatrix();
  void   CreateLTMMatrix();
  void   CreateLTMMatrixFromDigits();
  void   CreateLTMMatrixFromRaw(AliRawReader *fRawReader);

  static AliTOFTriggerMask *GetTOFTriggerMap() {return fTOFTrigMap;}
  static void PrepareTOFMapFromRaw(AliRawReader *fRawReader,Int_t deltaBC=13600);
  static void PrepareTOFMapFromDigit(TTree *treeD, Float_t startTimeHit=0, Float_t timeWidthTrigger=25);
 private:

  enum{
    kNLTM = 72,          //Number of LTM
    kNLTMchannels = 48,  //Number of channels in a LTM
    kNCTTM = 36,         //Number of CTTM per TOF side
    kNCTTMchannels = 24,  //Number of channels in a CTTM
    kNLTMtoTRDchannels = 8  //Number of channels in a CTTM
  };

  static AliTOFTriggerMask *fTOFTrigMap; // class with the TOF trigger map
  static AliTOFTriggerMask *fTOFTrigMask; // class with the TOF trigger mask

  static Int_t fgFromTriggertoDCS[72]; // dcs to trigger mapping

  AliTOFTrigger& operator=(const AliTOFTrigger &/*source*/); // ass. op.
  AliTOFTrigger(const AliTOFTrigger & tr);

  void    GetCTTMIndex(Int_t *detind, Int_t *indexCTTM);
  static void    GetLTMIndex(const Int_t * const detind, Int_t *LTMIndex);
  Bool_t  fLTMmatrix[kNLTM][kNLTMchannels];         //LTM matrix 
  Bool_t  fLTMarray[kNCTTM];        //LTM array for UPpurposes
  Bool_t  fCTTMmatrixFront[kNCTTM][kNCTTMchannels];//CTTM matrix for TOP FPGA 
  Bool_t  fCTTMmatrixBack[kNCTTM][kNCTTMchannels]; //CTTM matrix for BOTTOM FPGA
  Int_t   fHighMultTh;             //threshold for High Multiplicity trigger
  Int_t   fppMBTh;                 //threshold for pp Minimum Bias trigger
  Int_t   fMultiMuonTh;            //threshold for Multi Muon trigger 
  Int_t   fUPTh;                   //threshold for Ultra-Per coll trigger 
  Float_t fdeltaminpsi;            //min delta phi for J/psi decay (UP trigger)
  Float_t fdeltamaxpsi;            //max delta phi for J/psi decay (UP trigger)
  Float_t fdeltaminro;             //min delta phi for ro decay (UP trigger)
  Float_t fdeltamaxro;             //max delta phi for ro decay (UP trigger)
  Int_t   fstripWindow;            //strip window for triggering

  Bool_t fSel1,fSel2,fSel3,fSel4; // ppMB, PbPbMB2, PbPbMB3, PbPbUP

  UInt_t fPowerMask[kNCTTMchannels+1];  // mask for 24 TDC channels

  Int_t fNCrateOn; // number of crate fired
  Int_t fNMaxipadOn; // number of Maxipad fired
  Int_t fNMaxipadOnAll; // number of Maxipad fired w/o TDC dead mask

  // aggiungere larghezza finestra temporale e tempo0 in ns
  Float_t fStartTimeHit;      // time window start after channel equalization (subtraction of the minimal time per channel default 0 ns)
  Float_t fTimeWidthTrigger;  // time window width (default 25 ns)
   
  ClassDef(AliTOFTrigger,3)  // TOF Trigger Detector class
};
#endif

