//
// Class AliMixInputHandlerInfo
//
// AliMixInputHandlerInfo is interface with mixed
// input handlers
//
// author:
//        Martin Vala (martin.vala@cern.ch)
//
#ifndef ALIMIXINPUTHANDLERINFO_H
#define ALIMIXINPUTHANDLERINFO_H
#include <TArrayI.h>
#include <TNamed.h>

class TTree;
class TChain;
class TChainElement;
class AliInputEventHandler;
class AliMixInputHandlerInfo : public TNamed {

public:
   AliMixInputHandlerInfo(const char *name = "defautlTree", const char *title = "Defautl tree");
   virtual ~AliMixInputHandlerInfo();
   TChain *GetChain();

   void AddChain(TChain *chain);
//     void AddTreeToChain(TTree *tree);
   void AddTreeToChain(const char *path);

   void PrepareEntry(TChainElement *te, Long64_t entry, AliInputEventHandler *eh, Option_t *opt);

   void SetZeroEntryNumber(Long64_t num) { fZeroEntryNumber = num; }
   TChainElement *GetEntryInTree(Long64_t &entry);
   Long64_t      GetEntries();

private:
   TChain    *fChain;              // current chain
   TArrayI   fChainEntriesArray;   // array of entries of every chaing
   Long64_t  fZeroEntryNumber;     // zero entry number (will be used when we will delete not needed chains)
   Bool_t    fNeedNotify;          // flag if Notify is needed for current input handler

   AliMixInputHandlerInfo(const AliMixInputHandlerInfo &handler);
   AliMixInputHandlerInfo &operator=(const AliMixInputHandlerInfo &handler);

   ClassDef(AliMixInputHandlerInfo, 1); // Mix Input Handler info
};

#endif // ALIMIXINPUTHANDLERINFO_H
