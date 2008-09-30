//
// Class AliRsnPairMgr
//
// A collection of pairs for an analysis.
// The function of this collection is just for purposes of well-sorting
// the analyzed pairs into upper-level groups, in the case of a wide
// analysis containing many resonances at once, or different settings for the same one.
//
// Each PairMgr will result in a separate list of histograms, which
// can be seen as a folder in the output file, whose name is given by this object.
//
// author: M. Vala (email: martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliRsnPairMgr.h"

ClassImp(AliRsnPairMgr)

//_____________________________________________________________________________
AliRsnPairMgr::AliRsnPairMgr(const char*name) :
    TNamed(name,name),
    fPairs(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPairMgr::~AliRsnPairMgr()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnPairMgr::AddPair(AliRsnPair *pair)
{
//
// Adds pair
//

  fPairs.Add((AliRsnPair*)pair);
}

//_____________________________________________________________________________
void AliRsnPairMgr::PrintPairs()
{
//
// Prints all pairs
//

  AliRsnPair *pair = 0;
  TObjArrayIter next(&fPairs);
  while ((pair = (AliRsnPair*)next())) pair->Print();
}
