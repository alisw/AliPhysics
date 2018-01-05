#ifndef CEPRawFMDBuffer_H
#define CEPRawFMDBuffer_H

#include "TObject.h"
#include "AliESDFMD.h"

class CEPRawFMDBuffer : public TObject 
{

  private:
    /// Number of final cells is equal to the number of sectors
    UShort_t        fNCells;           
    /// Mulitplicity in certain cell
    Float_t *       fMult;           

  public:
                    CEPRawFMDBuffer();
                    ~CEPRawFMDBuffer();
    // Modifiers
    /// Reset the member variables
    void            Reset();
    
    // Setter functions 
    /// Set the number of FMD cells to fill
    void            SetFMDnCells           (UInt_t nCells) { fNCells = nCells; }
    /// Set the cell multiplicity 
    void            SetFMDCellMultiplicity (Float_t mult, UInt_t i); 

    /// Global variable setter
    void            SetFMDVariables    (AliESDFMD* FMDObj,   UInt_t nCells=140);
    
    // Accessors
    /// Returns the number of cells/sectors
    UInt_t          GetFMDnCells()                   const { return fNCells;  }
    /// Returns the amplitde of a given sector 
    Float_t         GetFMDCellMultiplicity(UInt_t i) const;


    ClassDef(CEPRawFMDBuffer,1)
};

#endif
