#include "AliMFTCALayer.h"

#include "AliMFTCALadder.h"
#include "AliMFTCACell.h"
#include "AliMFTCAHit.h"

ClassImp(AliMFTCALayer)

//___________________________________________________________________________
AliMFTCALayer::AliMFTCALayer() :
TObject(),
fID(-1),
fNhits(0),
fNcells(0),
fNladders(0)
{
  
  fCells  = new TClonesArray("AliMFTCACell",  1000);
  fHits   = new TClonesArray("AliMFTCAHit",   1000);
  fLadders  = new TClonesArray("AliMFTCALadder",  100);
  
}

//___________________________________________________________________________
void AliMFTCALayer::Clear(const Option_t *) {
  
  if (fCells) fCells->Clear("C");
  if (fHits)  fHits->Clear("C");
  
  fNhits  =  0;
  fNcells =  0;
  
}

//___________________________________________________________________________
void AliMFTCALayer::ClearCells()
{
  
  if (fCells) fCells->Clear("C");
  
  fNcells =  0;
  
}

//___________________________________________________________________________
AliMFTCACell *AliMFTCALayer::AddCell()
{
  
  new ((*fCells)[fNcells++]) AliMFTCACell();
  AliMFTCACell *cell = (AliMFTCACell*)fCells->At(fCells->GetLast());
  
  return cell;
  
}

//___________________________________________________________________________
AliMFTCAHit *AliMFTCALayer::AddHit()
{
  
  new ((*fHits)[fNhits++]) AliMFTCAHit();
  AliMFTCAHit *hit = (AliMFTCAHit*)fHits->At(fHits->GetLast());
  
  return hit;
  
}

//___________________________________________________________________________
AliMFTCALadder *AliMFTCALayer::GetLadderID(Int_t id) {
  
  AliMFTCALadder *ladder;
  Int_t i = 0;
  for (; i < fNladders; i++) {
    ladder = (AliMFTCALadder*)fLadders->At(i);
    if (ladder->GetID() == id) break;
  }
  if (i == fNladders) {
    new ((*fLadders)[fNladders++]) AliMFTCALadder();
    ladder = (AliMFTCALadder*)fLadders->At(fLadders->GetLast());
    ladder->SetID(id);
  }
  
  return ladder;
  
}
