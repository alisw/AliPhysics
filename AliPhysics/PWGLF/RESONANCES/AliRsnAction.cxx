//
// *** Class AliRsnAction ***
//
//  Base class for Action
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Jan Musinsky (jan.musinsky@cern.ch)
//

#include <TObjArray.h>

#include "AliRsnAction.h"

ClassImp(AliRsnAction)

//__________________________________________________________________________________________________
AliRsnAction::AliRsnAction(const char *name,const char *title) : TNamed(name,title)
{
//
// Default constructor
//

}

//__________________________________________________________________________________________________
AliRsnAction::AliRsnAction(const AliRsnAction &copy) : TNamed(copy)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnAction &AliRsnAction::operator=(const AliRsnAction &copy)
{
//
// Assignment constructor
//
   TNamed::operator=(copy);
   if (this == &copy)
      return *this;

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnAction::~AliRsnAction()
{
//
// Destructor
//
}

//__________________________________________________________________________________________________
Bool_t AliRsnAction::InitAction(TList */*outList*/,TObjArray */*objects*/)
{
//
// Init of Action
//

   return kTRUE;

}

//__________________________________________________________________________________________________
Bool_t AliRsnAction::ExecAction(TObjArray */*objects*/)
{
//
// Execute of action
//

   return kTRUE;

}
