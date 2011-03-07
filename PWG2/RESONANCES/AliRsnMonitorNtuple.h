//
// *** Class AliRsnMonitorNtuple ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnMonitorNtuple_H
#define AliRsnMonitorNtuple_H

#include "AliRsnMonitor.h"

class TList;
class TNtuple;

class AliRsnMonitorNtuple : public AliRsnMonitor {
public:

   AliRsnMonitorNtuple(const char *name = "default", AliRsnDaughterDef *def = 0);
   AliRsnMonitorNtuple(const AliRsnMonitorNtuple &copy);
   AliRsnMonitorNtuple& operator=(const AliRsnMonitorNtuple&);
   ~AliRsnMonitorNtuple();

   Bool_t         AddValue(AliRsnValue*const val);
   void           GenerateNtuple(const char *prefix = "", TList *list = 0);
   virtual void   Compute();
   virtual void   Init(const char *prefix, TList *list);

private:

   TClonesArray  fValues;  // single values computed from analyzed objects
   TNtuple      *fNtuple;  // ntuple computed with values

   ClassDef(AliRsnMonitorNtuple, 2)
};

#endif

