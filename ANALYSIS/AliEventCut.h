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
#include "AliBaseEventCut.h"

class AliAOD;

class AliEventCut: public TObject
{
  public: 
    AliEventCut();
    AliEventCut(const AliEventCut& in);
    virtual ~AliEventCut();
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
    void AddBasePartCut(AliBaseEventCut* ebcut){fBaseCuts.Add(ebcut);}
    
  protected:
    TObjArray fBaseCuts;
  private:
    ClassDef(AliEventCut,1)
};

class AliEmptyEventCut: public TObject
{
  public: 
    AliEmptyEventCut(){}
    virtual ~AliEmptyEventCut(){}
    
    Bool_t Pass(AliAOD* aod) const {return kFALSE;}//always accept
    
  protected:
  private:
    ClassDef(AliEmptyEventCut,1)
};

#endif
