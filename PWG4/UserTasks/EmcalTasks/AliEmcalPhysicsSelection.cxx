/* $Id$ */

#include "AliEmcalPhysicsSelection.h"
#include "AliESDEvent.h"
#include <AliLog.h>

ClassImp(AliEmcalPhysicsSelection)

//__________________________________________________________________________________________________
UInt_t AliEmcalPhysicsSelection::GetSelectionMask(const TObject* obj) 
{ 
  const AliESDEvent *ev = static_cast<const AliESDEvent*>(obj);
  UInt_t res = IsCollisionCandidate(ev); 

  if (fExcludeFastOnly) {
    if (res & AliVEvent::kFastOnly) {
      Int_t ncells = ev->GetEMCALCells()->GetNumberOfCells();
      if (ncells>0) {
        AliFatal(Form("Number of cells %d, even though EMCAL should not be in fast only partition.",ncells));
      }
      return 0;
    }
  }

  return res;
}
