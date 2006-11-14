// $Header$

#ifndef REVE_PLEX_H
#define REVE_PLEX_H

#include <Reve/Reve.h>

#include <TObject.h>
#include <TArrayC.h>

#include <vector>

namespace Reve {

/**************************************************************************/
// VoidCPlex
/**************************************************************************/

class VoidCPlex
{
private:
  VoidCPlex(const VoidCPlex&);            // Not implemented
  VoidCPlex& operator=(const VoidCPlex&); // Not implemented

protected:
  Int_t fS; // Size of atom
  Int_t fN; // Number of atoms in a chunk

  Int_t fSize;     // Size of container, number of atoms
  Int_t fVecSize;  // Number of allocated chunks
  Int_t fCapacity; // Available capacity within the chunks

  std::vector<TArrayC*> fChunks; // Memory blocks

  void ReleaseChunks();

public:
  VoidCPlex();
  VoidCPlex(Int_t atom_size, Int_t chunk_size);
  virtual ~VoidCPlex();

  void Reset(Int_t atom_size, Int_t chunk_size);
  void Refit();

  Int_t    S() const { return fS; }
  Int_t    N() const { return fN; }
  
  Int_t    Size()     const { return fSize; }
  Int_t    VecSize()  const { return fVecSize; }
  Int_t    Capacity() const { return fCapacity; }

  Char_t* Atom(Int_t idx)   { return fChunks[idx/fN]->fArray + idx%fN*fS; }
  Char_t* Chunk(Int_t chk)  { return fChunks[chk]->fArray; }
  Int_t   NAtoms(Int_t chk) { return (chk < fVecSize-1) ? fN : (fSize-1)%fN + 1; }

  Char_t* NewAtom();
  Char_t* NewChunk();

  ClassDef(VoidCPlex, 1)
};

inline Char_t* VoidCPlex::NewAtom()
{
  Char_t *a = (fSize >= fCapacity) ? NewChunk() : Atom(fSize);
  ++fSize;
  return a;
}


/**************************************************************************/
// Some-class CPlex
/**************************************************************************/

template<class T>
class CPlex : public VoidCPlex
{
private:
  CPlex(const CPlex&);            // Not implemented
  CPlex& operator=(const CPlex&); // Not implemented

public:
  CPlex()                 : VoidCPlex() {}
  CPlex(Int_t chunk_size) : VoidCPlex(sizeof(T), chunk_size) {}
  virtual ~CPlex() {}

  void Reset(Int_t chunk_size) { Reset(sizeof(T), chunk_size); }

  T* At(Int_t idx)  { return reinterpret_cast<T*>(Atom(idx)); }
  T& Ref(Int_t idx) { return *At(idx); }

  ClassDef(CPlex, 1)
}; // endclass CPlex

}

#endif
