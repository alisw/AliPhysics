//
// Emcal Container Base class
//
// Author: M. Verweij


#include <TClonesArray.h>

#include "AliEmcalJet.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliNamedArrayI.h"

#include "AliEmcalContainer.h"

ClassImp(AliEmcalContainer)

//________________________________________________________________________
AliEmcalContainer::AliEmcalContainer():
  TNamed("AliEmcalContainer","AliEmcalContainer"),
  fClArrayName(),
  fClassName(),
  fIsParticleLevel(kFALSE),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
AliEmcalContainer::AliEmcalContainer(const char *name):
  TNamed(name,name),
  fClArrayName(),
  fClassName(),
  fIsParticleLevel(kFALSE),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
void AliEmcalContainer::SetArray(AliVEvent *event) 
{
  // Get array from event.

  const AliVVertex *vertex = event->GetPrimaryVertex();
  if (vertex) vertex->GetXYZ(fVertex);

  if (!fClArrayName.IsNull() && !fClArray) {
    fClArray = dynamic_cast<TClonesArray*>(event->FindListObject(fClArrayName));
    if (!fClArray) {
      AliError(Form("%s: Could not retrieve array with name %s!", GetName(), fClArrayName.Data())); 
      return;
    }
  } else {
    return;
  }

  if (!fClassName.IsNull()) {
    TString objname(fClArray->GetClass()->GetName());
    TClass cls(objname);
    if (!cls.InheritsFrom(fClassName)) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from %s!", 
		    GetName(), cls.GetName(), fClArrayName.Data(), fClassName.Data())); 
      fClArray = 0;
    }
  }

  fLabelMap = dynamic_cast<AliNamedArrayI*>(event->FindListObject(fClArrayName + "_Map"));
}

//________________________________________________________________________
Int_t AliEmcalContainer::GetIndexFromLabel(Int_t lab) const
{ 
  if (fLabelMap) {
    if (lab < fLabelMap->GetSize()) {
      return fLabelMap->At(lab); 
    }
    else {
      AliDebug(3,Form("%s_AliEmcalContainer::GetIndexFromLabel - Label not found in the map, returning -1...",fClArrayName.Data()));
      return -1;
    }
  }
  else {
    AliDebug(3,Form("%s_AliEmcalContainer::GetIndexFromLabel - No index-label map found, returning label...",fClArrayName.Data()));
    return lab; 
  }
}
