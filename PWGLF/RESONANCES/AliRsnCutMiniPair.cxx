//
// This cut definition works on mini-pairs for the 'mini' subpackage.
// Since cuts on mini-pairs can usually be just kinematic selections
// or kinematic comparisons between 4-momenta of daughters, they are all
// implemented in this class, by the use of an enumeration which allows
// the user to choose between all possibilities.
//

#include <Riostream.h>

#include "AliRsnMiniPair.h"
#include "AliRsnCutMiniPair.h"

ClassImp(AliRsnCutMiniPair)

//__________________________________________________________________________________________________
AliRsnCutMiniPair::AliRsnCutMiniPair(const char *name, EType type) :
   AliRsnCut(name, AliRsnTarget::kMother),
   fType(type)
{
//
// Constructor
//
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutMiniPair::IsSelected(TObject *obj)
{
//
// Global check
//

   AliRsnMiniPair *pair = dynamic_cast<AliRsnMiniPair *>(obj);
   if (!pair) {
      AliError("This cut applies only to mini-pairs");
      return kFALSE;
   }

   switch (fType) {
      case kRapidityRange:
         fCutValueD = pair->Y(0);
         return OkRangeD();
      case kRapidityRangeMC:
         fCutValueD = pair->Y(1);
         return OkRangeD();
      case kRapidityFiducialRegion:
         fCutValueD = pair->Y(0);
         fRefPtValueD = pair->Pt(0);
         return OkRangeD();
      case kRapidityFiducialRegionMC:
         fCutValueD = pair->Y(1);
         fRefPtValueD = pair->Pt(1);
         return OkRangeD();  
      case kMomentumComparison:
         AliWarning("TODO: implement this");
         return kTRUE;
      case kDCAproduct:
         fCutValueD = pair->DCAProduct();
         return OkRangeD();
      case kDeltaCosRange:
         fCutValueD = pair->DeltaCos(kFALSE);
         return OkRangeD();
      case kPhiVRange:
         fCutValueD = pair->PhiV(kFALSE);
         return OkRangeD();
      case kContainsV0Daughter:
         return pair->ContainsV0Daughter();
      case kMassRange:
         fCutValueD = pair->InvMass(0);
         return OkRangeD();
      case kAsymRange:
         fCutValueD = pair->PairAsymmetry(0);
         return OkRangeD();
      case kAsymRangeMC:
         fCutValueD = pair->PairAsymmetry(1);
         return OkRangeD();
      case kPseudorapidityRange:
         fCutValueD = pair->Eta(0);
         return OkRangeD();
      case kPseudorapidityRangeMC:
         fCutValueD = pair->Eta(1);
         return OkRangeD();
      case kPassesOOBPileupCut:
         return pair->PassesOOBPileupCut();

   default:
         AliWarning("Undefined enum value");
         return kTRUE;
   }
}
