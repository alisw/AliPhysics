#ifndef CEPRawFMDBuffer_H
#define CEPRawFMDBuffer_H

#include "TObject.h"
#include "AliESDFMD.h"

class CEPRawFMDBuffer : public TObject 
{

  private:
    /// Number of sectors
    UShort_t            fNCells;           
    /// Mulitplicity in sector
    Float_t             fMult[140];

  public:
                        CEPRawFMDBuffer();
                        CEPRawFMDBuffer(const CEPRawFMDBuffer& fb);
                        ~CEPRawFMDBuffer();
    // Modifiers
    /// Reset the member variables
    void                Reset();
    void                Copy(TObject &obj) const;
    CEPRawFMDBuffer &   operator=(const CEPRawFMDBuffer& source);
    
    // Setter functions 
    /// Set the cell multiplicity 
    void                SetFMDCellMultiplicity(Float_t mult, UInt_t i); 

    /// Global variable setter
    void                SetFMDVariables(AliESDFMD* FMDObj);
    
    // Accessors
    /// Returns the number of cells/sectors
    UInt_t              GetFMDnCells()                   const { return fNCells;  }
    /// Returns the global multiplicity in all fmd cells
    Float_t             GetFMDTotalMultiplicity()        const;
    /// Returns the amplitde of a given sector 
    Float_t             GetFMDCellMultiplicity(UInt_t i) const;


    ClassDef(CEPRawFMDBuffer,1)
};

#endif
