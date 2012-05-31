// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCData_H
#define AliEveTPCData_H

#include <TEveUtil.h>

#include <TObject.h>

#include <vector>

class TTree;
class AliTPCRawStreamV3;

class AliEveTPCSectorData;

//------------------------------------------------------------------------------
// AliEveTPCData
//
// Container for TPC data for all sectors.

class AliEveTPCData : public TObject, public TEveRefCnt
{
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
  void LoadRaw(AliTPCRawStreamV3& input, Bool_t spawnSectors=kTRUE, Bool_t warn=kFALSE);

protected:
  std::vector<AliEveTPCSectorData*>  fSectors; // Vector of sector-data.
  Int_t                      fSectorBlockSize; // Block-size of sector-data.
  Short_t                    fLoadThreshold;   // Threshold at load-time.
  Short_t                    fLoadPedestal;    // Pedestal at load-time.
  Bool_t                     fAutoPedestal;    // If true determine pedestals automatically for each pad.

private:
  AliEveTPCData(const AliEveTPCData&);            // Not implemented
  AliEveTPCData& operator=(const AliEveTPCData&); // Not implemented

  ClassDef(AliEveTPCData, 0); // Manages TPC data for an event.
};

#endif
