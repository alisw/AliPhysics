//
// *** Class AliRsnTrainManager ***
//
//  Base class for Action
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Jan Musinsky (jan.musinsky@cern.ch)
//

#ifndef ALIRSNTRAINMANAGER_H
#define ALIRSNTRAINMANAGER_H

#include <TNamed.h>

class TMap;

class AliRsnTrainManager : public TNamed {
public:
   AliRsnTrainManager(const char *name="RsnTrainManager",const char *title="Resonances Train Manager");
   AliRsnTrainManager(const AliRsnTrainManager &copy);
   AliRsnTrainManager &operator=(const AliRsnTrainManager &copy);
   virtual ~AliRsnTrainManager();

   static void         SetGlobalStr(const char *key, const char *value,Bool_t verbose=kTRUE);
   static void         SetGlobalInt(const char *key, Int_t value,Bool_t verbose=kTRUE);
   static void         SetGlobalDbl(const char *key, Double_t value,Bool_t verbose=kTRUE);
   static void         SetGlobalObj(const char *key, TObject *value,Bool_t verbose=kTRUE);

   static const char  *GetGlobalStr(const char *key, Bool_t &valid);
   static Int_t        GetGlobalInt(const char *key, Bool_t &valid);
   static Double_t     GetGlobalDbl(const char *key, Bool_t &valid);
   static TObject     *GetGlobalObj(const char *key, Bool_t &valid);

   virtual void     Print(Option_t *option="") const;

   TMap               *GetGlobals() { return fGlobals; }
   static AliRsnTrainManager *GetRsnTrainManager() { return fgRsnTrainManager; }

private:
   TMap                      *fGlobals; // Map with global variables
   static AliRsnTrainManager *fgRsnTrainManager;

   ClassDef(AliRsnTrainManager, 1)
};

#endif
