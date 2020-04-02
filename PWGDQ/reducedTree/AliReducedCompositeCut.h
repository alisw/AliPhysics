// Class for cutting on ALICE AliReducedVarManager information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   06/02/2018

#ifndef ALIREDUCEDCOMPOSITECUT_H
#define ALIREDUCEDCOMPOSITECUT_H

#include "AliReducedInfoCut.h"
#include "TList.h"

//_________________________________________________________________________
class AliReducedCompositeCut : public AliReducedInfoCut {

public:
   AliReducedCompositeCut(Bool_t useAND=kTRUE);
   AliReducedCompositeCut(const Char_t* name, const Char_t* title, Bool_t useAND=kTRUE);
   virtual ~AliReducedCompositeCut();

   void AddCut(AliReducedInfoCut* cut) {fCuts.Add(cut);};
   
   Bool_t GetUseAND() const {return fOptionUseAND;}
   Int_t    GetNCuts() const {return fCuts.GetEntries();}
   
   virtual Bool_t IsSelected(TObject* obj);
   virtual Bool_t IsSelected(Float_t* values);
   virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
protected:
   Bool_t fOptionUseAND;            // true (default): apply AND on all cuts; false: use OR
   TList   fCuts;                            // list of cuts
    
   AliReducedCompositeCut(const AliReducedCompositeCut &c);
   AliReducedCompositeCut& operator= (const AliReducedCompositeCut &c);
  
   ClassDef(AliReducedCompositeCut,2);
};

#endif
