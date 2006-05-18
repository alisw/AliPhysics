// $Header$

#ifndef ALIEVE_TPCData_H
#define ALIEVE_TPCData_H

#include <Reve/Reve.h>

#include <TObject.h>

#include <vector>

class TTree;
class AliTPCRawStream;
class AliTPCRawStreamOld;

namespace Alieve {

class TPCSectorData;

class TPCData : public TObject, public Reve::ReferenceCount
{
protected:
  std::vector<TPCSectorData*>  fSectors;
  Int_t                        fSectorBlockSize;
  Short_t                      fLoadThreshold;
  Short_t                      fLoadPedestal;

public:
  TPCData();
  virtual ~TPCData();

  void CreateSector(Int_t sector);
  void CreateAllSectors();

  TPCSectorData* GetSectorData(Int_t sector, Bool_t spawnSectors=kFALSE);

  Int_t GetSectorBlockSize()   const { return fSectorBlockSize; }
  void  SetSectorBlockSize(Int_t bs) { fSectorBlockSize = bs; }

  Short_t GetLoadThreshold()     const { return fLoadThreshold; }
  void    SetLoadThreshold(Short_t lt) { fLoadThreshold = lt; }

  Short_t GetLoadPedestal()     const { return fLoadPedestal; }
  void    SetLoadPedestal(Short_t lp) { fLoadPedestal = lp; }

  void LoadDigits(TTree* tree, Bool_t spawnSectors=kTRUE);
  void LoadRaw(AliTPCRawStream&    input, Bool_t spawnSectors=kTRUE);
  void LoadRaw(AliTPCRawStreamOld& input, Bool_t spawnSectors=kTRUE, Bool_t warn=kFALSE);

  ClassDef(TPCData, 1); // Manages TPC data for an event.
}; // endclass TPCData

}

#endif
