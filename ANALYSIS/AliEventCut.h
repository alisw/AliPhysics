#ifndef ALIEVENTCUT_H
#define ALIEVENTCUT_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliEventCut
//
// Event cut. It has list of base event cuts. 
// Each of base event cut checks only one property.
// Logical base cuts also exists that point to other base cuts.
// Using them one can build complicated cut with binary tree structure
//
///////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>
#include "AliEventBaseCut.h"

class AliAOD;
enum EEventCutProperty;

class AliEventCut: public TObject
{
  public: 
    AliEventCut();
    AliEventCut(const AliEventCut& in);
    virtual ~AliEventCut();
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
    void           AddBasePartCut(AliEventBaseCut* ebcut);

    void           SetNChargedRange(Int_t min,Int_t max, Double_t etamin = -10.0,Double_t etamax = 10.0);
    
  protected:
    AliEventBaseCut* FindCut(EEventCutProperty prop);
    
    TObjArray fBaseCuts;
  private:
    ClassDef(AliEventCut,1)
};

class AliEventEmptyCut: public TObject
{
  public: 
    AliEventEmptyCut(){}
    virtual ~AliEventEmptyCut(){}
    
    Bool_t Pass(AliAOD* aod) const {return kFALSE;}//always accept
    
  protected:
  private:
    ClassDef(AliEventEmptyCut,1)
};

#endif
