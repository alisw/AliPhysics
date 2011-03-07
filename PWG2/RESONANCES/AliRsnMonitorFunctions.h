//
// *** Class AliRsnMonitorFunctions ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnMonitorFunctions_H
#define AliRsnMonitorFunctions_H

#include "AliRsnMonitor.h"

class AliRsnFunction;

class AliRsnMonitorFunctions : public AliRsnMonitor {
public:

   AliRsnMonitorFunctions(const char *name = "default", AliRsnDaughterDef *def = 0);
   AliRsnMonitorFunctions(const AliRsnMonitorFunctions &copy);
   AliRsnMonitorFunctions& operator=(const AliRsnMonitorFunctions&);
   ~AliRsnMonitorFunctions();

   void           AddFunction(AliRsnFunction* const fcn);
   TList*         GenerateHistograms(const char *prefix = "", TList *list = 0);
   virtual void   Compute();
   virtual void   Init(const char *prefix, TList *list);

protected:

   TClonesArray   fFunctions;    // a list of functions which generate histograms

   ClassDef(AliRsnMonitorFunctions, 1)
};

#endif

