// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TPCData_H
#define ALIEVE_TPCData_H

#include <TEveUtil.h>

#include <TObject.h>

#include <vector>

class TTree;
class AliTPCRawStream;
class AliTPCRawStreamOld;


class AliEveTPCSectorData;

class AliEveTPCData : public TObject, public TEveRefCnt
{
protected:
  std::vector<AliEveTPCSectorData*>  fSectors;
  Int_t                        fSectorBlockSize;
  Short_t                      fLoadThreshold;
  Short_t                      fLoadPedestal;
  Bool_t                       fAutoPedestal;

public:
  AliEveTPCData();
  virtual ~AliEveTPCData();

  void CreateSector(Int_t sector);
  void CreateAllSectors();
  void DropAllSectors();
  void DeleteAllSectors();

  AliEveTPCSectorData* GetSectorData(Int_t sector, Bool_t spawnSectors=kFALSE);

  Int_t GetSectorBlockSize()   const { return fSectorBlockSize; }
  void  SetSectorBlockSize(Int_t bs) { fSectorBlockSize = bs; }

  Short_t GetLoadThreshold()     const { return fLoadThreshold; }
  void    SetLoadThreshold(Short_t lt) { fLoadThreshold = lt; }

  Short_t GetLoadPedestal()     const { return fLoadPedestal; }
  void    SetLoadPedestal(Short_t lp) { fLoadPedestal = lp; }

  Bool_t GetAutoPedestal()     const { return fAutoPedestal; }
  void   SetAutoPedestal(Bool_t ap)  { fAutoPedestal = ap; }

  void LoadDigits(TTree* tree, Bool_t spawnSectors=kTRUE);
  void LoadRaw(AliTPCRawStream& input, Bool_t spawnSectors=kTRUE, Bool_t warn=kFALSE);

  ClassDef(AliEveTPCData, 1); // Manages TPC data for an event.
}; // endclass AliEveTPCData

#endif
