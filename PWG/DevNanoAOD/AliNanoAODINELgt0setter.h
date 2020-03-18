#ifndef AliNanoAODINELgt0setter_h
#define AliNanoAODINELgt0setter_h

#include "AliNanoAODCustomSetter.h"
#include "AliNanoAODHeader.h"
#include "AliAODEvent.h"
#include "AliMultSelectionTask.h"


class AliNanoAODINELgt0setter : public AliNanoAODCustomSetter
{
public:
  AliNanoAODINELgt0setter(const char * name = "AliNanoAODINELgt0setter") : AliNanoAODCustomSetter(name), fIndex(-1) {;}
  virtual ~AliNanoAODINELgt0setter() {;}
  inline virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head , TString varListHeader  ) {
    if (fIndex == -1)
      fIndex = varListHeader.Contains("cstINELgt0") ? head->GetVarIndex("cstINELgt0") : -2;
    if (fIndex >= 0)
      head->SetVar(fIndex, AliMultSelectionTask::IsINELgtZERO(event) ? 1 : 0);
  }
  virtual void SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack) {
      // This function will not be used in this class.
      return;
  };

  protected:
  bool fGoodToGo;
  int fIndex;

  ClassDef(AliNanoAODINELgt0setter, 1)
};



#endif /* AliNanoAODINELgt0setter_h */
