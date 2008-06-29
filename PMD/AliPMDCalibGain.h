#ifndef ALIPMDCALIBGAIN_H
#define ALIPMDCALIBGAIN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

class TH1F;  
class AliRawReader;

class AliPMDCalibGain : public TObject
{
 public:
  AliPMDCalibGain() ;              // ctor
  AliPMDCalibGain(const AliPMDCalibGain &pmdcalibgain);  // copy constructor
  AliPMDCalibGain &operator=(const AliPMDCalibGain &pmdcalibgain); // assignment op

  virtual ~AliPMDCalibGain() ;           // dtor

  Int_t ExtractPedestal();               // Assign the pedestal values
  void  ReadIntermediateFile();          // reads the intermediate file
  void  WriteIntermediateFile();         // writes the intermediate file

  Bool_t ProcessEvent(AliRawReader *rawReader, TObjArray *pmdddlcont);  //Looks for iso cells

  void Analyse(TTree *gaintree);
  
 private:

  enum
      {
	  kDet    = 2,   // Number of Planes
	  kMaxSMN = 24,  // Number of Modules
	  kMaxRow = 48,  // Number of Rows
	  kMaxCol = 96   // Number of Columns
      };

  Float_t fDetCount[kDet];                             //counter detector wise
  Float_t fDetIso[kDet];
  Float_t fSMIso[kDet][kMaxSMN];
  Float_t fSMCount[kDet][kMaxSMN];                     // counter
  Float_t fCellIso[kDet][kMaxSMN][kMaxRow][kMaxCol];   // adc of iso cells
  Float_t fCellCount[kDet][kMaxSMN][kMaxRow][kMaxCol]; // counter

  Float_t fPedMeanRMS[kDet][kMaxSMN][kMaxRow][kMaxCol];// Pedestal Mean
  FILE    *fpw;         // write the intermediate file

ClassDef(AliPMDCalibGain,4)        // description 
};
#endif // ALIPMDCALIBGAIN_H
