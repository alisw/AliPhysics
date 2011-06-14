//
// This cut selects the AliRsnDaughter objects pointing
// to tracks with a well defined true particle species,
// defined through its PDG code or species according the
// enumeration defined in AliRsnDaughter class.
// ---
// Using this cut on data results in no tracks passing it.
//

#include "AliRsnCutTrue.h"

ClassImp(AliRsnCutTrue)

//__________________________________________________________________________________________________
AliRsnCutTrue::AliRsnCutTrue(const char *name, Int_t pdg) :
   AliRsnCut(name, AliRsnTarget::kDaughter, pdg)
{
//
// Constructor version #1:
// pass directly the PDG code
//
}

//__________________________________________________________________________________________________
AliRsnCutTrue::AliRsnCutTrue(const char *name, AliRsnDaughter::ESpecies species) :
   AliRsnCut(name, AliRsnTarget::kDaughter, AliRsnDaughter::SpeciesPDG(species))
{
//
// Constructor version #2:
// pass the species from AliRsnDaughter enum, which is converted into PDG code
//
}

//__________________________________________________________________________________________________
AliRsnCutTrue::AliRsnCutTrue(const AliRsnCutTrue &copy) :
   AliRsnCut(copy)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnCutTrue& AliRsnCutTrue::operator=(const AliRsnCutTrue &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   return (*this);
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutTrue::IsSelected(TObject *obj)
{
//
// Check:
// if the MC reference is present, recover PDG
// and check if it matches the required one, in absolute value.
//

   // convert target
   if (!TargetOK(obj)) return kFALSE;
   
   // check if MC is present
   if (!fDaughter->GetRefMC()) {
      AliError("Cannot check cut 'AliRsnCutTrue' without MC information");
      return kFALSE;
   }
   
   // compare PDG
   fCutValueI = fDaughter->GetPDGAbs();
   return OkValueI();
}

   
   
