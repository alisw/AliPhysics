#ifndef AliPHOSCPVDA1_H
#define AliPHOSCPVDA1_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSCpvDA1 accumulates histograms with amplitudes per CPV channel.
// It is intended to run at DAQ or HLT computers.
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1.h"
#include "TFile.h"

class AliPHOSCpvDA1 : public TNamed {
  
 public:
  
  AliPHOSCpvDA1(Int_t module);
  AliPHOSCpvDA1(const AliPHOSCpvDA1& );
  AliPHOSCpvDA1& operator= (const AliPHOSCpvDA1& );
  ~AliPHOSCpvDA1();
  
  void  FillHistograms(Float_t e[128][56]);
  Int_t GetModule() { return fMod; }
  void  UpdateHistoFile();
  
 private:

  TFile* fHistoFile;            // root file to store histograms in
  TH1F*  fCharge[128][56];      // charge deposited on CPV pads
  Int_t  fMod;                  // PHOS module number (0..4)
  
  ClassDef(AliPHOSCpvDA1,1)

};

#endif
