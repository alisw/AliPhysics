//
// Class AliRsnPairManager
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
#include "AliRsnPair.h"

#include "AliRsnPairManager.h"

ClassImp(AliRsnPairManager)

//_____________________________________________________________________________
AliRsnPairManager::AliRsnPairManager(const char*name) :
    AliRsnVManager(name)
{
//
// Default constructor
//

  AliDebug(AliLog::kDebug +2, "<-");
  AliDebug(AliLog::kDebug +2, "->");
}

//_____________________________________________________________________________
void AliRsnPairManager::Add(TObject* objPair)
{
//
// Adds a new AliRsnPair to the list owned by this object.
//

  AliDebug(AliLog::kDebug+2, "<-");
  AliRsnPair *pair = dynamic_cast<AliRsnPair*>(objPair);

  if (!pair) {
    AliWarning(Form("Pair is %p. Skipping ...", pair));
    return;
  }

  AliDebug(AliLog::kDebug+1, Form("Adding %s [%d entries] ...", pair->GetPairName().Data(), fArray.GetEntries()));
  fArray.Add((AliRsnPair*)pair);

  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
void AliRsnPairManager::AddPair(AliRsnPair* pair)
{
//
// Adds a new AliRsnPair to the list owned by this object.
//

  Add(pair);
}


//_____________________________________________________________________________
void AliRsnPairManager::Print(Option_t* /*dummy*/) const
{
//
// Overload of TObject::Print() method.
// With respect to the other print method, adds a title string.
//

  AliDebug(AliLog::kDebug+2,"<-");

  AliInfo(Form("\t\t======== Pair Manager %s ========", GetName()));
  PrintArray();

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPairManager::PrintArray() const
{
//
// Prints all pairs
//

  AliDebug(AliLog::kDebug+2,"<-");

  AliRsnPair *pair = 0;
  TObjArrayIter next(&fArray);
  while ((pair = (AliRsnPair*)next())) pair->Print();

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPairManager::InitAllPairs(TList* list)
{
//
// Initialize all pairs, and builds a TList of histograms
// which are created by each of them, in order to link it
// to the output handler in the AnalysisTasks.
//

  AliDebug(AliLog::kDebug+2, "<-");

//   TList *list = new TList();
//   list->SetName(GetName());
//   list->SetOwner();

  AliRsnPair *pair = 0;
  TObjArrayIter next(&fArray);

  Int_t i = 0;
  while ((pair = (AliRsnPair*)next())) {
    if (!pair) continue;
    AliDebug(AliLog::kDebug+1, Form("InitAllPairs of the PairManager(%s) [%d] ...", pair->GetPairName().Data(), i++));
    pair->GenerateHistograms(GetName(),list);
  }

  AliDebug(AliLog::kDebug+2, "->");
//   return list;
}

//_____________________________________________________________________________
void AliRsnPairManager::ProcessAllPairs
(AliRsnPIDIndex *pidIndexes1, AliRsnEvent *ev1, AliRsnPIDIndex *pidIndexes2, AliRsnEvent *ev2)
{
//
// Processes one (single-event analysis) or two (event-mixing) events
// to fill histograms in all stored pairs.
//

  AliDebug(AliLog::kDebug+2, "<-");

  AliRsnPair *pair = 0;
  TObjArrayIter next(&fArray);

  Int_t i=0;
  while ((pair = (AliRsnPair*)next())) {
    if (!pair) continue;
    AliDebug(AliLog::kDebug+1, Form("ProcessAllPairs of the PairManager(%s) [%d] ...", pair->GetPairName().Data(), i++));
    pair->LoopPair(pidIndexes1, ev1, pidIndexes2, ev2);
  }

  AliDebug(AliLog::kDebug+2, "->");
}
