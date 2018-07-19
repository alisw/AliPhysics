#ifndef CEPRawV0Buffer_H
#define CEPRawV0Buffer_H

#include "AliESDVZERO.h"

class CEPRawV0Buffer : public TObject 
{

  private:
    /// Number of V0 cells
    UInt_t              fNCells;            // this is fixed (therefore hardcoded in constructor)
    /// Multiplicity in cell
    Float_t             fMult[64];
    /// ADC charge in cell of V0
    Float_t             fAdc[64];
    /// Leading time measured by TDC
    Float_t             fTime[64];
    /// signal width in given cell
    Float_t             fSigWidth[64];
    /// Number of fired PM in V0A & V0C
    Short_t             fNPMV0A;
    Short_t             fNPMV0C;

  public:
                        CEPRawV0Buffer();
    virtual             ~CEPRawV0Buffer() {};
    /// Modifiers
    void                Reset();
 
    /// Setter functions
    void                SetV0Multiplicity (Float_t mult,   UInt_t i)  { fMult[i]     = mult;   }
    void                SetV0Charge       (Float_t adc,    UInt_t i)  { fAdc[i]      = adc;    }
    void                SetV0Time         (Float_t ADtime, UInt_t i)  { fTime[i]     = ADtime; }
    void                SetV0Width        (Float_t sWidth, UInt_t i)  { fSigWidth[i] = sWidth; }
    void                SetV0NPMV0A       (Short_t NpmV0A)            { fNPMV0A      = NpmV0A; }
    void                SetV0NPMV0C       (Short_t NpmV0C)            { fNPMV0C      = NpmV0C; }

    /// Global Setter
    void                SetV0Variables(AliESDVZERO* V0Dobj);

    /// Accessors
    Float_t             GetV0Multiplicity(Int_t i)  const;
    Float_t             GetV0Charge(Int_t i)        const;
    Float_t             GetV0Time(Int_t i)          const;
    Float_t             GetV0Width(Int_t i)         const;
    Short_t             GetV0NumberPM_A()           const { return fNPMV0C; }
    Short_t             GetV0NumberPM_C()           const { return fNPMV0C; }
    /// Total accessors
    Float_t             GetV0TotalMultiplicity()    const;
    Float_t             GetV0TotalCharge()          const;
    Float_t             GetV0TotalTime()            const;
    Float_t             GetV0TotalSignalWidth()     const;

    /// Return number of cells
    UInt_t              GetNCells()                 const { return fNCells; }


    ClassDef(CEPRawV0Buffer,1)     // CEP raw track buffer
};

#endif
