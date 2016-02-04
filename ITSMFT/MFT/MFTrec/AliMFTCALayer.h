#ifndef AliMFTCALayer_H
#define AliMFTCALayer_H

#include "TObject.h"
#include "TClonesArray.h"
#include "AliMFTCALadder.h"

class AliMFTCACell;
class AliMFTCAHit;

class AliMFTCALayer : public TObject {
  
public:
  
  AliMFTCALayer();
  ~AliMFTCALayer() {};
  
  AliMFTCALayer (const AliMFTCALayer &layer);
  AliMFTCALayer &operator=(const AliMFTCALayer&);
  
  AliMFTCALadder *GetLadder(Int_t nl) { return (AliMFTCALadder*)fLadders->At(nl); }
  AliMFTCALadder *GetLadderID(Int_t id);
  
  virtual void Clear(const Option_t *);
  void ClearCells();
  
  void SetID(Int_t id) { fID = id; }
  const Int_t GetID() { return fID; }
  const Int_t GetNhits() { return fNhits; }
  const Int_t GetNcells() { return fNcells; }
  AliMFTCACell *AddCell();
  AliMFTCACell *GetCell(Int_t nc) { return (AliMFTCACell*)fCells->At(nc); }
  AliMFTCAHit  *AddHit();
  AliMFTCAHit  *GetHit(Int_t nh)  { return (AliMFTCAHit*)fHits->At(nh);   }
  
private:
  
  Int_t fID;            // Identifier
  Int_t fNhits;          // Number of hits
  Int_t fNcells;         // Number of cells: track segments with downstream
                         // end in this layer
  Int_t fNladders;       // number of ladders
  
  TClonesArray *fCells;     //! Array of cells
  TClonesArray *fHits;      //! Array of hits
  
  TClonesArray *fLadders;   //! array of ladders
  
  ClassDef(AliMFTCALayer,1);
  
};
#endif
