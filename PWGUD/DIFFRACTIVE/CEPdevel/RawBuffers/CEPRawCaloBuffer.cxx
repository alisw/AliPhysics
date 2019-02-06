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
  , fNCells(0)
  , fAmplitude(0x0)
  , fCellMCLabel(0x0)
  , fTime(0x0)
{

}

// as we have pointer member variables we have to add a few things to our class
//  - a copy constructor
//  - overwrite the assignment operator
//  - overwrite the Copy function of TObject
///
/// Copy constructor
///
// ____________________________________________________________________________
CEPRawCaloBuffer::CEPRawCaloBuffer(const CEPRawCaloBuffer& cb)
  : TObject(cb)
  , fNCells(cb.fNCells)
  , fAmplitude(0x0)
  , fCellMCLabel(0x0)
  , fTime(0x0)
{
    fAmplitude   = new Double_t[fNCells];
    fCellMCLabel = new Int_t[fNCells];
    fTime        = new Double_t[fNCells];

    for (UInt_t i(0); i<fNCells; i++){
        fAmplitude[i]   = cb.fAmplitude[i]; 
        fCellMCLabel[i] = cb.fCellMCLabel[i]; 
        fTime[i]        = cb.fTime[i]; 
    }

}

///
/// Assignment operator.
///
//__________________________________________________________________________
CEPRawCaloBuffer & CEPRawCaloBuffer::operator =(const CEPRawCaloBuffer& source)  
{
    if (&source == this) return *this;
    TObject::operator=(source);

    this->Reset();
  
    fNCells = source.fNCells;
    if (fNCells>0){
        fAmplitude   = new Double_t[fNCells];
        fTime        = new Double_t[fNCells];
        fCellMCLabel = new Int_t[fNCells];
        for (UInt_t i(0); i<fNCells; i++){
            fAmplitude[i]   = source.fAmplitude[i];
            fTime[i]        = source.fTime[i];
            fCellMCLabel[i] = source.fCellMCLabel[i];
        }
    }

    return *this;
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
    fNCells = 0;
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
Float_t CEPRawCaloBuffer::GetCaloTotalAmplitude() const
{
    Float_t totalAmpl(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalAmpl += fAmplitude[i];
    }

    return totalAmpl;
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
Float_t CEPRawCaloBuffer::GetCaloTotalTime() const
{
    Float_t totalTime(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalTime += fTime[i];
    }

    return totalTime;
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
