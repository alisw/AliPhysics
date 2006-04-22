/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpBusPatch.h,v 1.2 2006/03/17 11:35:58 ivana Exp $

/// \ingroup management
/// \class AliMpBusPatch
/// \brief Class that manages the maps buspatch<>DDL<>DE 
///
/// Calculates also the maximum DSP and buspatch numbers for a given DE
///
/// Author: Ch. Finck; Subatech Nantes

#ifndef ALI_MP_BUSPATCH_H
#define ALI_MP_BUSPATCH_H

#include <TObject.h>

#include <TExMap.h>

class TArrayI;

class AliMpBusPatch : public TObject
{

 public:

  AliMpBusPatch();
  AliMpBusPatch(const AliMpBusPatch& src);
  virtual ~AliMpBusPatch();

  // operators  
  AliMpBusPatch& operator = (const AliMpBusPatch& src) ;
  
  // methods
  void ReadBusPatchFile();
  void GetDspInfo(Int_t iCh, Int_t& iDspMax, Int_t* iBusPerDSP) const;

  Int_t    GetDEfromBus(Int_t busPatchId);
  TArrayI* GetBusfromDE(Int_t idDE);
  Int_t    GetDDLfromBus(Int_t busPatchId);

 private:

  TExMap fDetElemIdToBusPatch;       //! Map from idDE to BusPatch   
  TExMap fBusPatchToDetElem;         //! Map from BusPatch to idDE
  TExMap fBusPatchToDDL;             //! Map from BusPatch to iDDL

  Int_t fMaxBusPerCh[10];            //! max buspatch number per chamber

  ClassDef(AliMpBusPatch,1) //utility class for the motif type
};


#endif //ALI_MP_BUSPATCH_H
