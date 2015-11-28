#include "TString.h"
#include "TList.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliObservableBase.h"
#include "AliEventClassifierBase.h"
#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliObservableBase)

// Constructors
AliObservableBase::AliObservableBase()
: TNamed() {

};
AliObservableBase::AliObservableBase(const char* name, const char* title)
  : TNamed(name, title) {
  // asign valid pdgs; detour for cxx98
  Int_t pdgs[] = {111, 211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
  fSafePdgCodes.assign(&pdgs[0], &pdgs[0] + (sizeof(pdgs) / sizeof(Int_t)));
}
