#ifndef ALIPMDDDLRAWDATA_H
#define ALIPMDDDLRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDDigitization.h, Version 00        //
//                                                     //
//  Date   : September 20 2002                         //
//                                                     //
//-----------------------------------------------------//

#include <TObject.h>

class TClonesArray;
class TTree;

class AliPMDdigit;

class AliPMDDDLRawData:public TObject
{
 public:

  AliPMDDDLRawData();
  virtual ~AliPMDDDLRawData();

  void WritePMDRawData(TTree *treeD);
  void GetUMDigitsData(TTree *treeD, Int_t imodule, Int_t ium, Int_t ddlno,
		       Int_t & totword, UInt_t *buffer);
  void TransformS2H(Int_t smn, Int_t &irow, Int_t &icol);
  void GetMCMCh(Int_t ddlno, Int_t um, Int_t row, Int_t col,
		UInt_t &mcmno, UInt_t &chno);


 protected:

  TClonesArray *fDigits;    //! List of digits
  AliPMDdigit  *fPMDdigit;  //! Pointer to digits

  ClassDef(AliPMDDDLRawData,3)    // To make RAW Data
};
#endif

