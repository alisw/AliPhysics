#ifndef AliMFTCALadder_H
#define AliMFTCALadder_H

#include "TObject.h"

#include "AliMFTCACell.h"
#include "AliMFTCAHit.h"

class TClonesArray;

//_________________________________________________________________________________
class AliMFTCALadder : public TObject {
  
public:
  
  AliMFTCALadder();
  ~AliMFTCALadder() {};
  
  AliMFTCALadder (const AliMFTCALadder &ladder);
  AliMFTCALadder &operator=(const AliMFTCALadder&);

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
  
  Int_t     fID;             // Identifier
  Int_t     fNhits;          // Number of hits
  Int_t     fNcells;         // Number of cells: track segments with downstream end in
                             // this layer
  
  TClonesArray *fCells;     //! Array of cells
  TClonesArray *fHits;      //! Array of hits
  
  ClassDef(AliMFTCALadder,1);
  
};
#endif
