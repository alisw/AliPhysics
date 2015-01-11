#include "AliEmcalJetUtility.h"

ClassImp(AliEmcalJetUtility)

//______________________________________________________________________________
AliEmcalJetUtility::AliEmcalJetUtility() :
TNamed(),
  fJetTask(0),
  fInit(kFALSE)
{
  // Dummy constructor.

}

//______________________________________________________________________________
AliEmcalJetUtility::AliEmcalJetUtility(const char* name) :
  TNamed(name, name),
  fJetTask(0),
  fInit(kFALSE)
{
  // Default constructor.
}

//______________________________________________________________________________
AliEmcalJetUtility::AliEmcalJetUtility(const AliEmcalJetUtility &other) :
  TNamed(other),
  fJetTask(other.fJetTask),
  fInit(other.fInit)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliEmcalJetUtility& AliEmcalJetUtility::operator=(const AliEmcalJetUtility &other)
{
  // Assignment
  if (&other == this) return *this;
  TNamed::operator=(other);
  fJetTask = other.fJetTask;
  fInit = other.fInit;
  return *this;
}
