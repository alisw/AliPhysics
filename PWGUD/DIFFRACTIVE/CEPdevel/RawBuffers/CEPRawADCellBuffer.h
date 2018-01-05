#ifndef CEPRawADCellBuffer_H
#define CEPRawADCellBuffer_H

#include "TObject.h"
#include "AliESDAD.h"

class CEPRawADCellBuffer : public TObject 
{

  private:
    /// Number of AD cells
    UInt_t          fNCells;        // this is a fixed size (therefore hardcoded in constructor)
    /// Multiplicity in cell
    Float_t         fMult[16];
    /// ADC charge in cell 
    Float_t         fAdc[16];
    /// Leading time measured by TDC
    Float_t         fTime[16];
    /// Final decision based on average time of channel 
    /// from AliVVZERO:
    ///     enum Decision { kV0Invalid = -1, kV0Empty = 0, kV0BB, kV0BG, kV0Fake }
    Int_t           fDecisionADA;
    Int_t           fDecisionADC;


  public:
                    CEPRawADCellBuffer();
                    ~CEPRawADCellBuffer();
    /// Modifiers
    void            Reset();
    
 
    /// Setter functions
    void            SetADMultiplicity (Float_t mult,   UInt_t i)  { fMult[i]     = mult;      }
    void            SetADCharge       (Float_t adc,    UInt_t i)  { fAdc[i]      = adc;       }
    void            SetADTime         (Float_t ADtime, UInt_t i)  { fTime[i]     = ADtime;    }
    void            SetADDecision_A   (Int_t decisionA)           { fDecisionADA = decisionA; }
    void            SetADDecision_C   (Int_t decisionC)           { fDecisionADC = decisionC; }

    /// Global Setter
    void            SetADVariables(AliESDAD* ADobj);

    /// Accessors
    Float_t         GetADMultiplicity(Int_t i)  const;
    Float_t         GetADCharge(Int_t i)        const;
    Float_t         GetADTime(Int_t i)          const;
    Int_t           GetADADecision_A()          const { return fDecisionADA; }
    Int_t           GetADADecision_C()          const { return fDecisionADC; }

    /// Return number of cells
    UInt_t         GetNCells()                 const { return fNCells; }


    ClassDef(CEPRawADCellBuffer, 1)     // CEP raw track buffer
};

#endif
