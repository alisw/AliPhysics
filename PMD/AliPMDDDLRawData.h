#ifndef ALIPMDDDLRAWDATA_H
#define ALIPMDDDLRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : AliPMDDDLRawData.h, Version 01       //
//                                                     //
//  Date   : June 20 2006                              //
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
  AliPMDDDLRawData (const AliPMDDDLRawData &ddlraw);  // copy constructor
  AliPMDDDLRawData &operator=(const AliPMDDDLRawData &ddlraw); // assignment op

  virtual ~AliPMDDDLRawData();

  void WritePMDRawData(TTree *treeD);
  void GetUMDigitsData(TTree *treeD, Int_t imodule, Int_t ddlno,
		       Int_t *contentsBus, UInt_t busPatch[][1536]);
  void TransformS2H(Int_t smn, Int_t &irow, Int_t &icol);
  void GetMCMCh(Int_t ddlno, Int_t row, Int_t col,
		Int_t beginPatchBus, Int_t endPatchBus,
		Int_t *mcmperBus,
		Int_t *startRowBus, Int_t *startColBus,
		Int_t *endRowBus, Int_t *endColBus,
		Int_t & busno, UInt_t &mcmno, UInt_t &chno);

 protected:

  TClonesArray *fDigits;    //! List of digits
  //  AliPMDdigit  *fPMDdigit;  //! Pointer to digits

  ClassDef(AliPMDDDLRawData,6)    // To make RAW Data
};
#endif

