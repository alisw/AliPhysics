#ifndef CEPRawCaloBuffer_H
#define CEPRawCaloBuffer_H

#include "AliESDCaloCells.h"

class CEPRawCaloBuffer : public TObject {

  private:
    /// Number of AD cells
    UInt_t          fNCells;            // Get this with GetNumberOfCells from
    /// Multiplicity in cell
    Double_t *      fAmplitude;
    /// Cell time (s)
    Double_t *      fTime;
    /// MC Label for particle which deposited most energy in cell
    Int_t *         fCellMCLabel;

  public:
                    CEPRawCaloBuffer();
                    ~CEPRawCaloBuffer();
    /// Modifiers
    void            Reset();

 
    /// Setter functions
    void            SetCaloCellAmplitude (Double_t amp,  UInt_t i);
    void            SetCaloCellTime      (Double_t tme,  UInt_t i); 
    /// Set MC information
    void            SetCaloCellMCLabel   (Int_t MCl,     UInt_t i); 

    /// Global Setter
    void            SetCaloVariables(AliESDCaloCells* CellsObj);

    /// Accessors
    Float_t         GetCaloCellAmplitude(UInt_t i)     const;
    Float_t         GetCaloCellTime(UInt_t i)          const;
    Int_t           GetCaloCellMCLabel(UInt_t i)       const;

    /// Return number of cells
    UInt_t         GetNCells()                       const { return fNCells; }


    ClassDef(CEPRawCaloBuffer,1)     // CEP raw track buffer
};

#endif
