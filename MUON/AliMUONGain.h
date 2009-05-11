#ifndef ALIMUONGain_H
#define ALIMUONGain_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONDA
/// \brief Implementation of the pedestal and gain computing
/// 
//  Author: Alberto Baldisseri, JL Charvet (05/05/2009)

#include "AliMUONPedestal.h"

/* // global variables */
/* const Int_t kNChannels = AliMpConstants::ManuNofChannels(); */
/* const Int_t kADCMax    = 4095; */

class AliMUONGain : public AliMUONPedestal
{
  public:
    AliMUONGain();
    virtual ~AliMUONGain();

    TString WriteGainData(Int_t bp, Int_t manu, Int_t ch, Double_t p1, Double_t p2, Int_t threshold, Int_t q);
    TString WriteGainHeader();
    void MakePedStoreForGain(TString flatfile);
    void MakeGainStore(TString flatfile); 
    ///
    void SetAliRootDataFileName() {sprintf(fRootDataFileName,"MUONTRKGAINda_data.root");}
    ///
    Char_t* GetRootDataFileName() {return fRootDataFileName;}
    TString WriteDummyHeader();
    ///
    void SetAliInjCharge(Int_t charge) {fInjCharge = charge;}
    ///
    void SetAliPrintLevel(Int_t pri) {fPrintLevel = pri;}
    ///
    void SetAliInit(Int_t ini) {fnInit = ini;}
    ///
    void SetAliEntries(Int_t ent) {fnEntries = ent;}
    ///
    void SetAliNbpf1(Int_t nf1) {fnbpf1 = nf1;}
    ///
    void SetAliPlotLevel(Int_t plo) {fPlotLevel = plo;}

  private:
    Int_t fInjCharge; ///<
    Char_t fRootDataFileName[256]; ///<
    Int_t fnInit; ///<
    Int_t fnEntries; ///<
    Int_t fnbpf1; ///<
    Int_t fPrintLevel; ///< 
    Int_t fPlotLevel; ///< 
	
  ClassDef(AliMUONGain,1) // 
};

#endif
