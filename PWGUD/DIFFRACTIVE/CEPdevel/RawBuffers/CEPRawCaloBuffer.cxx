/// ////////////////////////////////////////////////////////////////////////////
///
/// CEP Calorimeter buffer
///
/// structure to hold information on the calorimeter
///
// ____________________________________________________________________________
#include "CEPRawCaloBuffer.h"
#include "CEPTrackBuffer.h"
#include <stdexcept>
#include <stdlib.h> 

ClassImp(CEPRawCaloBuffer)

// ____________________________________________________________________________
CEPRawCaloBuffer::CEPRawCaloBuffer()
  : TObject()
  , fNCells      (0)
  , fAmplitude   (0x0)
  , fCellMCLabel (0x0)
  , fTime        (0x0)
{

}

// ____________________________________________________________________________
CEPRawCaloBuffer::~CEPRawCaloBuffer()
{
    if (fAmplitude)     { delete [] fAmplitude;   fAmplitude = 0x0;   }
    if (fCellMCLabel)   { delete [] fCellMCLabel; fCellMCLabel = 0x0; }
    if (fTime)          { delete [] fTime;        fTime = 0x0;        }
}

// ____________________________________________________________________________
void CEPRawCaloBuffer::Reset()
{
    if (fAmplitude)     { delete [] fAmplitude;   fAmplitude = 0x0;   }
    if (fCellMCLabel)   { delete [] fCellMCLabel; fCellMCLabel = 0x0; }
    if (fTime)          { delete [] fTime;        fTime = 0x0;        }
    fNCells      = 0;
}

// ____________________________________________________________________________
void CEPRawCaloBuffer::SetCaloCellAmplitude(Double_t amp,  UInt_t i)
{
    try {
        if (i>=fNCells) throw std::out_of_range("Calo cell array: index out of range"); 
        fAmplitude[i] = amp;
    } catch (const std::out_of_range& e) {
        printf("%s\n",e.what());
        exit(EXIT_FAILURE);
    }
}

// ____________________________________________________________________________
void CEPRawCaloBuffer::SetCaloCellTime(Double_t tme,  UInt_t i)
{
    try {
        if (i>=fNCells) throw std::out_of_range("Calo cell array: index out of range"); 
        fTime[i] = tme;
    } catch (const std::out_of_range& e) {
        printf("%s\n",e.what());
        exit(EXIT_FAILURE);
    }    
}

void CEPRawCaloBuffer::SetCaloCellMCLabel(Int_t MCl,  UInt_t i)
{
    try {
        if (i>=fNCells) throw std::out_of_range("Calo cell array: index out of range"); 
        fCellMCLabel[i] = MCl;
    } catch (const std::out_of_range& e) {
        printf("%s\n",e.what());
        exit(EXIT_FAILURE);
    }    
}

// ____________________________________________________________________________
Float_t CEPRawCaloBuffer::GetCaloCellAmplitude(UInt_t i) const
{
    return (i<fNCells) ? fAmplitude[i] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Int_t CEPRawCaloBuffer::GetCaloCellMCLabel(UInt_t i) const
{
    return (i<fNCells) ? fCellMCLabel[i] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawCaloBuffer::GetCaloCellTime(UInt_t i) const
{
    return (i<fNCells) ? fTime[i] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
void CEPRawCaloBuffer::SetCaloVariables(AliESDCaloCells* CellsObj)
{
    // Global member variable setter
    fNCells = CellsObj->GetNumberOfCells();
    fAmplitude   = new Double_t[fNCells]; 
    fCellMCLabel = new Int_t[fNCells]; 
    fTime        = new Double_t[fNCells]; 
    for (UInt_t i(0); i<fNCells; i++) 
    {  
        this->SetCaloCellAmplitude(CellsObj->GetAmplitude(i), i);
        this->SetCaloCellTime(CellsObj->GetTime(i), i);
        this->SetCaloCellMCLabel(CellsObj->GetCellMCLabel(i), i);
    }
}

// ____________________________________________________________________________
