#ifndef ALIBASEEVENTCUT_H
#define ALIBASEEVENTCUT_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliBaseEventCut
//
// Base class for cauts that checks only one event property
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include "TObject.h"

class AliAOD;

class AliBaseEventCut: public TObject
{
  public: 
    AliBaseEventCut();
    virtual ~AliBaseEventCut(){}
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
  protected:
    virtual Double_t GetValue(AliAOD* aod) const = 0;
    
    Double_t fMin;//Minimum value
    Double_t fMax;//Maximum value
  private:
    ClassDef(AliBaseEventCut,1)
};

#endif
