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

class AliTOFrawData;

class AliTOFTrigger : public AliTriggerDetector
{
 public:
  AliTOFTrigger();  // constructor
  AliTOFTrigger(Int_t HighMultTh, Int_t ppMBTh, Int_t MultiMuonTh, Int_t UPTh, Float_t deltaminpsi, Float_t deltamaxpsi, Float_t deltaminro, Float_t deltamaxro, Int_t stripWindow); //constructor with parameters
  AliTOFTrigger(const AliTOFTrigger & tr);
  virtual ~AliTOFTrigger(){}  // destructor
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

  void GetMap(Bool_t **map);
  //void PrintMap(); // to be checked because of warning problems
  void GetTRDmap(Bool_t **map);
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

  void   CreateCTTMMatrix();
  void   CreateLTMMatrix();
  void   CreateLTMMatrixFromDigits();
  void   CreateLTMMatrixFromRaw(AliRawReader *fRawReader);
 private:

  enum{
    kNLTM = 72,          //Number of LTM
    kNLTMchannels = 48,  //Number of channels in a LTM
    kNCTTM = 36,         //Number of CTTM per TOF side
    kNCTTMchannels = 24,  //Number of channels in a CTTM
    kNLTMtoTRDchannels = 8  //Number of channels in a CTTM
  };

  void    GetCTTMIndex(Int_t *detind, Int_t *indexCTTM);
  void    GetLTMIndex(Int_t *detind, Int_t *LTMIndex);
  Bool_t  fLTMmatrix[kNLTM][kNLTMchannels];         //LTM matrix  
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

  ClassDef(AliTOFTrigger,1)  // TOF Trigger Detector class
};
#endif
