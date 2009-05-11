#ifndef ALIMUONPEDESTAL_H
#define ALIMUONPEDESTAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONPedestal
/// \brief Implementation of the pedestal computing
/// 
//  Author Alberto Baldisseri, JL Charvet 

#include <TObject.h>

#include "AliMpConstants.h"

class AliMUONVStore;

class TTimeStamp;
class THashTable;

// global variables
const Int_t kNChannels = AliMpConstants::ManuNofChannels();
const Int_t kADCMax    = 4095;

class AliMUONPedestal : public TObject
{
  public:
    AliMUONPedestal();
    virtual ~AliMUONPedestal();
    
    ///
    void SetAliNEvents(Int_t events) {fNEvents = events;}
    ///
    void SetAliRunNumber(Int_t run) {fRunNumber = run;}
    ///
    void SetAliNChannel(Int_t nch) {fNChannel = nch;}
    ///
    void SetAlifilcout(ofstream* stream) {fFilcout = stream;}
    ///
    TTimeStamp* GetDate() {return fDate;}
    ///
    THashTable* GetErrorBuspatchTable() {return fErrorBuspatchTable;}
    ///
    Char_t* GetHistoFileName() {return fHistoFileName;}
    void MakePed(Int_t bp,Int_t manu,Int_t ch,Int_t charge);
    void MakePedStore(TString flatfile);
    TString WritePedData(Int_t bp, Int_t manu, Int_t ch, Double_t mean, Double_t sigma);
    TString WritePedHeader();

    ///
    void SetprefixDA(char* folder) {sprintf(fprefixDA,"%s",folder);}
    ///
    void SetAliIndex(Int_t ind) {fIndex = ind;}

  protected:
    Int_t fNEvents; ///<
    Int_t fRunNumber; ///<
    Int_t fNChannel; ///<
    Int_t fNManu; ///<
    THashTable* fErrorBuspatchTable; ///< Table for buspatches with parity errors 
    TTimeStamp* fDate; ///<
    ofstream* fFilcout; ///<
    Char_t fHistoFileName[256]; ///<
    AliMUONVStore* fPedestalStore; ///<
    Int_t fIndex; ///<
    Char_t fprefixDA[256]; ///<

  private:
    /// Not implemented
    AliMUONPedestal(const AliMUONPedestal& rhs);
    /// Not implemented
    AliMUONPedestal& operator = (const AliMUONPedestal& rhs);

  ClassDef(AliMUONPedestal,1) // 
};

#endif
