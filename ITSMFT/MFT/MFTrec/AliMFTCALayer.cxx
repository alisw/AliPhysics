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
AliMFTCALayer::AliMFTCALayer(const AliMFTCALayer &layer) :
TObject(layer),
fID(layer.fID),
fNhits(layer.fNhits),
fNcells(layer.fNcells),
fNladders(layer.fNladders)
{
  
  // copy constructor
  
  fCells  = new TClonesArray("AliMFTCACell",  layer.fNcells);
  fHits   = new TClonesArray("AliMFTCAHit",   layer.fNhits);
  fLadders  = new TClonesArray("AliMFTCALadder",  layer.fNladders);
  
  AliMFTCALadder *caLadder;
  for (Int_t i = 0; i < layer.fNladders; i++) {
    caLadder = (AliMFTCALadder*)(layer.fLadders->At(i));
    new ((*fLadders)[i]) AliMFTCALadder(*caLadder);
  }
  AliMFTCACell *caCell;
  for (Int_t i = 0; i < layer.fNcells; i++) {
    caCell = (AliMFTCACell*)(layer.fCells->At(i));
    new ((*fCells)[i]) AliMFTCACell(*caCell);
  }
  AliMFTCAHit *caHit;
  for (Int_t i = 0; i < layer.fNhits; i++) {
    caHit = (AliMFTCAHit*)(layer.fHits->At(i));
    new ((*fHits)[i]) AliMFTCAHit(*caHit);
  }

}

//___________________________________________________________________________
AliMFTCALayer& AliMFTCALayer::operator=(const AliMFTCALayer& layer) 
{

  // assignment operator

  // check assignement to self
  if (this == &layer) return *this;

  TObject::operator=(layer);

  fID = layer.fID;
  fNhits = layer.fNhits;
  fNcells = layer.fNcells;
  fNladders = layer.fNladders;

  fCells  = new TClonesArray("AliMFTCACell",  layer.fNcells);
  fHits   = new TClonesArray("AliMFTCAHit",   layer.fNhits);
  fLadders  = new TClonesArray("AliMFTCALadder",  layer.fNladders);
  
  AliMFTCALadder *caLadder;
  for (Int_t i = 0; i < layer.fNladders; i++) {
    caLadder = (AliMFTCALadder*)(layer.fLadders->At(i));
    new ((*fLadders)[i]) AliMFTCALadder(*caLadder);
  }
  AliMFTCACell *caCell;
  for (Int_t i = 0; i < layer.fNcells; i++) {
    caCell = (AliMFTCACell*)(layer.fCells->At(i));
    new ((*fCells)[i]) AliMFTCACell(*caCell);
  }
  AliMFTCAHit *caHit;
  for (Int_t i = 0; i < layer.fNhits; i++) {
    caHit = (AliMFTCAHit*)(layer.fHits->At(i));
    new ((*fHits)[i]) AliMFTCAHit(*caHit);
  }

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
