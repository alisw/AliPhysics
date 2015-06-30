
#include "TClonesArray.h"

#include "AliMFTCARoad.h"

ClassImp(AliMFTCARoad)

//___________________________________________________________________________
AliMFTCARoad::AliMFTCARoad() :
TObject(),
fID(-1),
fNhits(0),
fNHitSta(0),
fLength(0),
fIsGood(kFALSE),
fNhitsInLayer(),
fLayer1(-1),
fLayer2(-1),
fNcellsInLayer()
{
  
  //fHits   = new TClonesArray("AliMFTCAHit",   1000);
  
  for (Int_t i = 0; i < fNDetMax; i++) {
    fHitsInLayer[i] = new TClonesArray("AliMFTCAHit",100);
    fNhitsInLayer[i] = 0;
    fCellsInLayer[i] = new TClonesArray("AliMFTCACell",10);
    fNcellsInLayer[i] = 0;
  }
  
}

//___________________________________________________________________________
void AliMFTCARoad::Clear(const Option_t *)
{
  
  //if (fHits) fHits->Clear("C");
  
  for (Int_t i = 0; i < fNDetMax; i++) {
    if (fHitsInLayer[i]) fHitsInLayer[i]->Clear("C");
    if (fCellsInLayer[i]) fCellsInLayer[i]->Clear("C");
    fNhitsInLayer[i] = 0;
    fNcellsInLayer[i] = 0;
  }
  
  fNhits   = 0;
  fNHitSta = 0;
  fLength  = 0;
  fIsGood  = kFALSE;
  fLayer1 = -1;
  fLayer2 = -1;
  
}

//___________________________________________________________________________
void AliMFTCARoad::AddHit(AliMFTCAHit *hit)
{
  
  //new ((*fHits)[fNhits++]) AliMFTCAHit(*hit);
  fNhits++;
  Int_t layer = hit->GetLayer();
  new ((*fHitsInLayer[layer])[fNhitsInLayer[layer]++]) AliMFTCAHit(*hit);
  
}

//___________________________________________________________________________
void AliMFTCARoad::AddCell(AliMFTCACell *cell) {
  
  Int_t layer = cell->GetLayers()[0];
  new ((*fCellsInLayer[layer])[fNcellsInLayer[layer]++]) AliMFTCACell(*cell);
  
}
//___________________________________________________________________________
AliMFTCACell *AliMFTCARoad::GetCellByGID(Int_t gid) {
  
  AliMFTCACell *cell;
  
  for (Int_t iL = 0; iL < (fNDetMax-1); iL++) {
    for (Int_t iC = 0; iC < GetNcellsInLayer(iL); iC++) {
      cell = GetCellInLayer(iL,iC);
      if (gid == cell->GetGID()) return cell;
    }
  }
  
  return 0;
  
}
