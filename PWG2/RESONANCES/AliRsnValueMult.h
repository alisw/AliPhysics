//
// Class AliRsnValueMult
//
// Inherits from AliRsnValue, and computer multiplicity of the event
// in several ways
//
// Author: A. Pulvirenti
// Email : alberto.pulvirenti@ct.infn.it
//

#ifndef ALIRSNVALUEMULT
#define ALIRSNVALUEMULT

#include "AliRsnValue.h"
#include "AliESDtrackCuts.h"

class AliRsnValueMult : public AliRsnValue
{
  public:
  
    enum EMode
    {
      kESDcuts,
      kNTracks,
      kNTracklets
    };
    
    AliRsnValueMult();
    AliRsnValueMult(const char *name, EValueType type, Int_t n = 0, Double_t min = 0.0, Double_t max = 0.0);
    AliRsnValueMult(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
    AliRsnValueMult(const char *name, EValueType type, Int_t n, Double_t *array);
    AliRsnValueMult(const AliRsnValueMult& copy) : AliRsnValue(copy),fMode(copy.fMode),fESDcuts(copy.fESDcuts) {}
    AliRsnValueMult& operator=(const AliRsnValueMult& copy) {AliRsnValue::operator=(copy); fMode=copy.fMode; fESDcuts=copy.fESDcuts;return (*this);}
    virtual ~AliRsnValueMult() { }
    
    AliESDtrackCuts* GetCuts() {return &fESDcuts;}
    EMode            GetMode() {return fMode;}
    void             SetMode(EMode mode) {fMode = mode;}
    
    virtual Bool_t   Eval(AliRsnMother * const mother, AliRsnPairDef * const pairDef, AliRsnEvent * const event);
    virtual Bool_t   Eval(AliRsnDaughter * const daughter, AliRsnEvent * const event);
    virtual void     Print(Option_t *option = "") const;
  
  protected:
  
    EMode           fMode;     // chosen method to compute multiplicity 
    AliESDtrackCuts fESDcuts;  // ESD track cuts necessary for one of the methods
    
  private:
  
    ClassDef(AliRsnValueMult,1)
};

#endif
