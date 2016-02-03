#ifndef AliMFTCARoad_H
#define AliMFTCARoad_H

#include "TObject.h"
#include "TClonesArray.h"
#include "AliMFTConstants.h"
#include "AliMFTCAHit.h"
#include "AliMFTCACell.h"

//_________________________________________________________________________________
class AliMFTCARoad : public TObject {
  
public:
  
  AliMFTCARoad();
  ~AliMFTCARoad() {};
  
  AliMFTCARoad (const AliMFTCARoad &road);
  AliMFTCARoad &operator=(const AliMFTCARoad&);
  
  void SetID(Int_t id) { fID = id; }
  const Int_t GetID() { return fID; }
  const Int_t GetNhits() { return fNhits; }
  void AddHit(AliMFTCAHit *hit);
  const Int_t GetNhitsInLayer(Int_t nl) { return fNhitsInLayer[nl]; }
  AliMFTCAHit *GetHitInLayer(Int_t nl, Int_t nh) { return (AliMFTCAHit*)fHitsInLayer[nl]->At(nh); }
  void SetNHitSta(Int_t n) { fNHitSta = n; }
  const Int_t GetNHitSta() { return fNHitSta; }
  void SetLength(Int_t l1, Int_t l2) {
    fLength = l2/2-l1/2;
    fLayer1 = l1;
    fLayer2 = l2;
  }
  const Int_t GetLength() { return fLength; }
  void SetGood() { fIsGood = kTRUE; }
  const Bool_t IsGood() { return fIsGood; }
  const Int_t GetLayer1() { return fLayer1; }
  const Int_t GetLayer2() { return fLayer2; }
  void AddCell(AliMFTCACell *cell);
  const Int_t GetNcellsInLayer(Int_t nl) { return fNcellsInLayer[nl]; }
  AliMFTCACell *GetCellInLayer(Int_t nl, Int_t nc) { return (AliMFTCACell*)fCellsInLayer[nl]->At(nc); }
  AliMFTCACell *GetCellByGID(Int_t gid);
  
  virtual void Clear(const Option_t *);
  
private:
  static const Int_t fNDetMax = AliMFTConstants::fNMaxPlanes;

  Int_t     fID;             // Identifier
  Int_t     fNhits;          // Number of hits
  Int_t     fNHitSta;        // Number of hit stations (or disks)
  Int_t     fLength;         // distance between the first and the last planes
  Bool_t    fIsGood;         // has the minimum requested number of stations
  Int_t fNhitsInLayer[fNDetMax]; // Number of hits per layer
  Int_t fLayer1;                 // first layer
  Int_t fLayer2;                 // last layer
  Int_t fNcellsInLayer[fNDetMax]; // Number of cells per layer
  
  TClonesArray *fHitsInLayer[fNDetMax];   //! Array of hits per layer
  TClonesArray *fCellsInLayer[fNDetMax];  //! Array of cells per layer
  
  ClassDef(AliMFTCARoad,1);
  
};

#endif

