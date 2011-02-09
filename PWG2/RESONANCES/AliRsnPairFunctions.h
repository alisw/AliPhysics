//
// *** Class AliRsnPairFunctions ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnPairFunctions_H
#define AliRsnPairFunctions_H

#include "AliRsnPair.h"

class TH1;
class TH2;
class TList;
class TArrayI;

class AliRsnEvent;
class AliRsnCutSet;
class AliRsnFunction;
class AliRsnValue;

class AliRsnPairFunctions : public AliRsnPair {
public:

   AliRsnPairFunctions(const char *name = "default", AliRsnPairDef *def = 0);
   AliRsnPairFunctions(const AliRsnPairFunctions &copy);
   AliRsnPairFunctions& operator=(const AliRsnPairFunctions&);
   ~AliRsnPairFunctions();

   void         AddFunction(AliRsnFunction*const fcn);
   TList*       GenerateHistograms(const char *prefix = "", TList *list = 0);
   virtual void Compute();
   virtual void Init(const char *prefix, TList *list);

protected:

   TClonesArray   fFunctions;    // a list of functions which generate histograms

   ClassDef(AliRsnPairFunctions, 2)
};

#endif

