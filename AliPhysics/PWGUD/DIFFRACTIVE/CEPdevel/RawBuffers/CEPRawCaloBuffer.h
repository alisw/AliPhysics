#ifndef CEPRawCaloBuffer_H
#define CEPRawCaloBuffer_H

#include "AliESDCaloCells.h"

class CEPRawCaloBuffer : public TObject {

  private:
    /// Number of AD cells
    // we need the [fNCells]. This tells the streamer
    // how many space it has to allocate
    UInt_t              fNCells;            // number of cells, this is a variable number
    /// Multiplicity in cell
    Double_t*           fAmplitude;         //[fNCells]
    /// Cell time (s)
    Double_t*           fTime;              //[fNCells]
    /// MC Label for particle which deposited most energy in cell
    Int_t*              fCellMCLabel;       //[fNCells]

  public:
                        CEPRawCaloBuffer();
                        // necessary as we have pointer variables
                        CEPRawCaloBuffer(const CEPRawCaloBuffer& cb);   
    virtual             ~CEPRawCaloBuffer();
    /// Modifiers
    void                Reset();

    CEPRawCaloBuffer &  operator=(const CEPRawCaloBuffer& source);
      
    /// Setter functions
    void                SetCaloCellAmplitude (Double_t amp,  UInt_t i);
    void                SetCaloCellTime      (Double_t tme,  UInt_t i); 
    /// Set MC information
    void                SetCaloCellMCLabel   (Int_t MCl,     UInt_t i); 

    /// Global Setter
    void                SetCaloVariables(AliESDCaloCells* CellsObj);

    /// Accessors
    Float_t             GetCaloCellAmplitude(UInt_t i)  const;
    Float_t             GetCaloCellTime(UInt_t i)       const;
    Int_t               GetCaloCellMCLabel(UInt_t i)    const;

    /// Sum Accessors
    Float_t             GetCaloTotalAmplitude()         const;  
    Float_t             GetCaloTotalTime()              const;  

    /// Return number of cells
    UInt_t              GetNCells()                     const { return fNCells; }


    ClassDef(CEPRawCaloBuffer,1)     // CEP raw track buffer
};

#endif
