// $Header$

#include "Plex.h"

using namespace Reve;

//______________________________________________________________________
// Plex
//
// A group of containers with chunked memory allocation.
//


//______________________________________________________________________
// VoidCPlex
//
// Non-structured (Void) Complete Plex.
// Allocation chunk can accommodate fN atoms of byte-size fS each.
// The chunks themselves are TArrayCs and are stored in a std::vector<TArrayC*>.
// Holes in the structure are not supported, neither is removal of atoms.
// The structure can be Refit() to occupy a single contiguous array.
//

void VoidCPlex::ReleaseChunks()
{
  for (Int_t i=0; i<fVecSize; ++i)
    delete fChunks[i];
  fChunks.clear();
}

VoidCPlex::VoidCPlex() : 
  fS(0), fN(0),
  fSize(0), fVecSize(0), fCapacity(0)
{}

VoidCPlex::VoidCPlex(Int_t atom_size, Int_t chunk_size) :
  fS(atom_size), fN(chunk_size),
  fSize(0), fVecSize(0), fCapacity(0)
{}

VoidCPlex::~VoidCPlex()
{
  ReleaseChunks();
}

/**************************************************************************/

void VoidCPlex::Reset(Int_t atom_size, Int_t chunk_size)
{
  ReleaseChunks();
  fS = atom_size;
  fN = chunk_size;
  fSize = fVecSize = fCapacity = 0;
}

void VoidCPlex::Refit()
{
  if (fSize == 0 || (fVecSize == 1 && fSize == fCapacity))
    return;

  TArrayC* one = new TArrayC(fS*fSize);
  Char_t*  pos = one->fArray;
  for (Int_t i=0; i<fVecSize; ++i)
  {
    Int_t size = fS * NAtoms(i);
    memcpy(pos, fChunks[i]->fArray, size);
    pos += size;
  }
  ReleaseChunks();
  fN = fCapacity = fSize;
  fVecSize = 1;
  fChunks.push_back(one);
}

/**************************************************************************/

Char_t* VoidCPlex::NewChunk()
{
  fChunks.push_back(new TArrayC(fS*fN));
  ++fVecSize;
  fCapacity += fN;
  return fChunks.back()->fArray;
}

/**************************************************************************/
/**************************************************************************/
