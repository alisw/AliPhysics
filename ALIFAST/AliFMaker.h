#ifndef AliFMaker_H
#define AliFMaker_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast virtual base class for Makers                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif
#ifndef ROOT_TClonesArray
#include <TClonesArray.h>
#endif

class TList;
class TBrowser;
class TChain;

class AliFMaker : public TNamed {

protected:

   Bool_t         fIsClonable;  //!True if Maker objects are clonable
   Int_t          fSave;        // = 1 if m-Maker to be saved in the Tree
   TObject       *fFruits;      //Pointer to maker fruits (result)
   TObject       *fClones;      //Pointer to clones of fruits
   TString        fBranchName;  //Name of branch (if any)
   TList         *fHistograms;  //Pointer to list supporting Maker histograms

public:
                  AliFMaker();
                  AliFMaker(const char *name, const char *title);
   virtual       ~AliFMaker();
   virtual void   Browse(TBrowser *b);
   virtual void   Clear(Option_t *option="");
   virtual void   Draw(Option_t *option="");
   virtual void   Finish();
   TList         *Histograms() {return fHistograms;}
   virtual void   Init();
   Bool_t         IsFolder() {return kTRUE;}
   TObject       *Fruit() {return fFruits;}
   TClonesArray  *Fruits() {return (TClonesArray*)fFruits;}
   TObject       *Clones() {return fClones;}
   virtual void   FillClone();
   virtual void   Make() = 0;
   virtual void   PrintInfo();
   virtual void   MakeBranch();
   virtual void   Save(Int_t save=1) {fSave = save;}
   virtual void   SetChainAddress(TChain *chain);

   ClassDef(AliFMaker, 1)   //AliFast virtual base class for Makers
};

#endif













