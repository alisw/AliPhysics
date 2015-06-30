
#include "TClonesArray.h"
#include "AliMFTCALadder.h"


ClassImp(AliMFTCALadder)


//___________________________________________________________________________
AliMFTCALadder::AliMFTCALadder() :
TObject(),
fID(-1),
fNhits(0),
fNcells(0)
{
  
  fCells  = new TClonesArray("AliMFTCACell",  1000);
  fHits   = new TClonesArray("AliMFTCAHit",   1000);
  
}

//___________________________________________________________________________
void AliMFTCALadder::Clear(const Option_t *)
{
  
  if (fCells) fCells->Clear("C");
  if (fHits)  fHits->Clear("C");
  
  fNhits  =  0;
  fNcells =  0;
  
}

//___________________________________________________________________________
void AliMFTCALadder::ClearCells()
{
  
  if (fCells) fCells->Clear("C");
  
  fNcells =  0;
  
}

//___________________________________________________________________________
AliMFTCACell *AliMFTCALadder::AddCell()
{
  
  new ((*fCells)[fNcells++]) AliMFTCACell();
  AliMFTCACell *cell = (AliMFTCACell*)fCells->At(fCells->GetLast());
  
  return cell;
  
}

//___________________________________________________________________________
AliMFTCAHit *AliMFTCALadder::AddHit()
{
  
  new ((*fHits)[fNhits++]) AliMFTCAHit();
  AliMFTCAHit *hit = (AliMFTCAHit*)fHits->At(fHits->GetLast());
  
  return hit;
  
}
