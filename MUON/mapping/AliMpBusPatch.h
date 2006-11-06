/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpBusPatch.h,v 1.5 2006/05/24 13:58:16 ivana Exp $

/// \ingroup management
/// \class AliMpBusPatch
/// \brief Class that manages the maps buspatch<>DDL<>DE 
///
/// Calculates also the maximum DSP and buspatch numbers for a given DE
///
/// \author Ch. Finck; Subatech Nantes

#ifndef ALI_MP_BUSPATCH_H
#define ALI_MP_BUSPATCH_H

#include <TObject.h>

#include <TExMap.h>
#include <TArrayI.h>

class TExMapIter;

class AliMpBusPatch : public TObject
{

 public:

  AliMpBusPatch();
  virtual ~AliMpBusPatch();

  // methods
  void ReadBusPatchFile();
  void GetDspInfo(Int_t iDDL, Int_t& iDspMax, Int_t* iBusPerDSP) const;

  Int_t    GetDEfromBus(Int_t busPatchId);
  TArrayI* GetBusfromDE(Int_t idDE);
  Int_t    GetDDLfromBus(Int_t busPatchId);
  void     AddBus(Int_t nDDL, Int_t busPatch);
  void     AddDetElem(Int_t nDDL, Int_t detElem);
  void     Sort();
  TArrayI  GetBusInDDL(Int_t nDDL) {return fBusInDDL[nDDL];}
  TArrayI  GetDeInDDL(Int_t nDDL) {return fDeInDDL[nDDL];}
  Int_t    NextBusInDDL(Int_t iDDL);
  void     ResetBusItr(Int_t iDDL);
  TExMapIter   GetBusItr() {return TExMapIter(&fBusPatchToDDL);}

 private:
  AliMpBusPatch(const AliMpBusPatch& src);
  AliMpBusPatch& operator = (const AliMpBusPatch& src) ;

  TExMap fDetElemIdToBusPatch;       //!< Map from idDE to BusPatch   
  TExMap fBusPatchToDetElem;         //!< Map from BusPatch to idDE
  TExMap fBusPatchToDDL;             //!< Map from BusPatch to iDDL

  TArrayI fBusInDDL[20];             //!< buspatch array per DDL
  Int_t   fBusItr[20];               //!< buspatch in DDL iterator
  TArrayI fDeInDDL[20];              //!< detElem array per DDL


  void Sort(TArrayI& arr, Int_t start, Int_t end);

  ClassDef(AliMpBusPatch,1) //utility class for the motif type
};


#endif //ALI_MP_BUSPATCH_H
