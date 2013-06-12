//
// Emcal Container Base class
//
// Author: M. Verweij

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include "AliEmcalJet.h"
#include "AliVEvent.h"
#include "AliLog.h"

#include "AliEmcalContainer.h"

ClassImp(AliEmcalContainer)

//________________________________________________________________________
AliEmcalContainer::AliEmcalContainer():
  TNamed("AliEmcalContainer","AliEmcalContainer"),
  fClArray(0),
  fClArrayName()
{
  // Default constructor.

}

//________________________________________________________________________
AliEmcalContainer::AliEmcalContainer(const char *name):
  TNamed(name,name),
  fClArray(0),
  fClArrayName()
{
  // Standard constructor.

}

//________________________________________________________________________
AliEmcalContainer::~AliEmcalContainer()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalContainer::SetArray(AliVEvent *event, const char *clname) {

  // Get array from event.

  if (!fClArrayName.IsNull() && !fClArray) {
    fClArray = dynamic_cast<TClonesArray*>(event->FindListObject(fClArrayName));
    if (!fClArray) {
      AliWarning(Form("%s: Could not retrieve array with name %s!", GetName(), fClArrayName.Data())); 
      return;
    }
  } else {
    return;
  }

  if (!clname)
    return;

  TString objname(fClArray->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom(clname)) {
    AliWarning(Form("%s: Objects of type %s in %s are not inherited from %s!", 
                    GetName(), cls.GetName(), fClArrayName.Data(), clname)); 
    return;
  }
  return;

 
}

