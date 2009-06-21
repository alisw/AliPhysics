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

  Int_t ExtractPedestal(const Char_t *rootFile);    // pedestal 
  Int_t ExtractHotChannel(const Char_t *rootFile);  // Hotchannel root file 
  void  ReadTempFile(const Char_t *tempFile);       // read inter file
  void  WriteTempFile(const Char_t *tempFile);      // write inter file

  Bool_t ProcessEvent(AliRawReader *rawReader, TObjArray *pmdddlcont);  //Looks for iso cells

  void Analyse(TTree *gaintree, TTree *meantree);
  void FindHotCell(TTree *hottree, Float_t xvar); // finds hot cell
  
 private:

  enum
      {
	  kDet    = 2,   // Number of Planes
	  kMaxSMN = 24,  // Number of Modules per plane
	  kMaxRow = 48,  // Number of Rows
	  kMaxCol = 96   // Number of Columns
      };

  Float_t fSMIso[kDet][kMaxSMN];
  Float_t fSMCount[kDet][kMaxSMN];                     // counter
  Float_t fCellIso[kDet][kMaxSMN][kMaxRow][kMaxCol];   // adc of iso cells
  Float_t fCellCount[kDet][kMaxSMN][kMaxRow][kMaxCol]; // counter of iso cell
  Float_t fNhitCell[kDet][kMaxSMN][kMaxRow][kMaxCol];  // counter
  Float_t fPedMeanRMS[kDet][kMaxSMN][kMaxRow][kMaxCol];// Pedestal Mean
  Float_t fHotFlag[kDet][kMaxSMN][kMaxRow][kMaxCol];   // HotChannel Flag 
   
  Float_t fCountSm[kDet][kMaxSMN];     // event counter for each module
  Float_t fTempnhit[kDet][kMaxSMN];    // hit frequency of each module
  Float_t fTempnhitSq[kDet][kMaxSMN];  // square of hit freq. of each mod.
  FILE    *fpw;                        // write the temp file

ClassDef(AliPMDCalibGain,7)            // description 
};
#endif // ALIPMDCALIBGAIN_H
