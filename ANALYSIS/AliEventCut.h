#ifndef ALIEVENTCUT_H
#define ALIEVENTCUT_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
//
//
//
///////////////////////////////////////////////////////////

#include "TObject.h"

class AliAOD;
class TObjArray;

class AliEventCut: public TObject
{
  public: 
    AliEventCut();
    virtual ~AliEventCut();
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
    
  protected:
    TObjArray* fBaseCuts;
  private:
    ClassDef(AliEventCut,1)
};

#endif
