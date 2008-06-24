#ifndef ALIRSNMVCUTMGR_H
#define ALIRSNMVCUTMGR_H

#include "TNamed.h"
#include "TObjArray.h"

#include "AliRsnCut.h"
// #include "AliRsnCutSet.h"
// class AliRsnCut;
class AliRsnCutSet;

/**
  @author Martin Vala <Martin.Vala@cern.ch>
*/
class AliRsnCutMgr : public TNamed
{
  public:

//     enum ECutSetType {
//       kParticle= 0,
//       kPair,
//       kMixEventFinderCut,
//       kLastCutSetIndex
//     };


    AliRsnCutMgr();
    AliRsnCutMgr(const char *name, const char* title);

    ~AliRsnCutMgr();

    void SetCutSet(AliRsnCut::ECutSetType type,AliRsnCutSet* cutset);
    AliRsnCutSet* GetCutSet(AliRsnCut::ECutSetType type) { return fCutSets[type];}

    Bool_t IsSelected(AliRsnCut::ECutSetType type,TObject *obj);

  private:

    AliRsnCutMgr(const AliRsnCutMgr &cut):TNamed(cut) {}
    AliRsnCutMgr& operator=(const AliRsnCutMgr& /*cut*/) {return *this;}

    AliRsnCutSet *fCutSets[AliRsnCut::kLastCutSetIndex];

    ClassDef ( AliRsnCutMgr,1 );
};

#endif
