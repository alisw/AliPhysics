//
// *** Class AliRsnAction ***
//
//  Base class for Action
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Jan Musinsky (jan.musinsky@cern.ch)
//

#ifndef ALIRSNACTION_H
#define ALIRSNACTION_H

#include <TNamed.h>

class AliRsnAction : public TNamed {
public:
   AliRsnAction(const char* name="noName",const char *title="No Title");
   AliRsnAction(const AliRsnAction &copy);
   AliRsnAction &operator=(const AliRsnAction &copy);
   virtual ~AliRsnAction();
   
   virtual Bool_t InitAction(TList *outList=0,TObjArray *objects=0);
   virtual Bool_t ExecAction(TObjArray *objects=0);

protected:
   
   ClassDef(AliRsnAction, 1)
};

#endif
