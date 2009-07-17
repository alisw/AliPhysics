//
// Class AliRsnVManager
//
// Base "manager" class.
// It is built in order to manage a list of objects which share
// the same level in the work flow of the analysis.
// This base class contains simply the list of "child" objects
// and the methods to add objects to the list or retrieve the list.
//
// author     : M. Vala       [martin.vala@cern.ch]
// revised by : A. Pulvirenti [alberto.pulvirenti@ct.infn.it]
//

#include "AliLog.h"

#include "AliRsnVManager.h"

ClassImp(AliRsnVManager)

//_____________________________________________________________________________
AliRsnVManager::AliRsnVManager(const char*name) :
    TNamed(name, name),
    fArray(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnVManager::~AliRsnVManager()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnVManager::Add(TObject*const obj)
{
//
// Add a new object in the list.
//

  fArray.Add((TObject*)obj);
}

//_____________________________________________________________________________
void AliRsnVManager::Print(Option_t* /*dummy*/) const
{
//
// Overload of the standard TObject::Print() method.
//

  PrintArray();
}

//_____________________________________________________________________________
void AliRsnVManager::PrintArray() const
{
//
// Calls the "Print()" method of all objects
// stored in the list, to print their informations.
//

  TObject *obj = 0;
  TObjArrayIter next(&fArray);
  while ((obj = (TObject*)next())) obj->Print();
}
