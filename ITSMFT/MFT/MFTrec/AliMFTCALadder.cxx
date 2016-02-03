
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
AliMFTCALadder::AliMFTCALadder(const AliMFTCALadder &ladder) :
TObject(ladder),
fID(ladder.fID),
fNhits(ladder.fNhits),
fNcells(ladder.fNcells)
{

  // copy constructor
  
  fCells  = new TClonesArray("AliMFTCACell",  ladder.fNcells);
  fHits   = new TClonesArray("AliMFTCAHit",   ladder.fNhits);

  AliMFTCACell *caCell;
  for (Int_t i = 0; i < ladder.fNcells; i++) {
    caCell = (AliMFTCACell*)(ladder.fCells->At(i));
    new ((*fCells)[i]) AliMFTCACell(*caCell);
  }
  AliMFTCAHit *caHit;
  for (Int_t i = 0; i < ladder.fNhits; i++) {
    caHit = (AliMFTCAHit*)(ladder.fHits->At(i));
    new ((*fHits)[i]) AliMFTCAHit(*caHit);
  }

}

//___________________________________________________________________________
AliMFTCALadder& AliMFTCALadder::operator=(const AliMFTCALadder& ladder) 
{

  // assignment operator

  // check assignement to self
  if (this == &ladder) return *this;

  TObject::operator=(ladder);

  fID = ladder.fID;
  fNhits = ladder.fNhits;
  fNcells = ladder.fNcells;

  fCells  = new TClonesArray("AliMFTCACell",  ladder.fNcells);
  fHits   = new TClonesArray("AliMFTCAHit",   ladder.fNhits);

  AliMFTCACell *caCell;
  for (Int_t i = 0; i < ladder.fNcells; i++) {
    caCell = (AliMFTCACell*)(ladder.fCells->At(i));
    new ((*fCells)[i]) AliMFTCACell(*caCell);
  }
  AliMFTCAHit *caHit;
  for (Int_t i = 0; i < ladder.fNhits; i++) {
    caHit = (AliMFTCAHit*)(ladder.fHits->At(i));
    new ((*fHits)[i]) AliMFTCAHit(*caHit);
  }

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
