// ////////////////////////////////////////////////////////////////////////////
//
// CEP AD buffer containing information on the AD cells
//
// structure to hold track information
//
// ____________________________________________________________________________
#include "CEPRawADCellBuffer.h"
#include "CEPTrackBuffer.h"

ClassImp(CEPRawADCellBuffer)

// ____________________________________________________________________________
CEPRawADCellBuffer::CEPRawADCellBuffer()
  : TObject()
  , fNCells(16)  
  , fDecisionADA (CEPTrackBuffer::kdumval)
  , fDecisionADC (CEPTrackBuffer::kdumval)
{

  this->Reset();

}

// ____________________________________________________________________________
CEPRawADCellBuffer::~CEPRawADCellBuffer()
{
  
}

// ____________________________________________________________________________
void CEPRawADCellBuffer::Reset()
{
    fNCells      = 16;
    fDecisionADA = CEPTrackBuffer::kdumval;
    fDecisionADC = CEPTrackBuffer::kdumval;
    for (UInt_t i(0); i<fNCells; i++){
        fMult[i] = 0.0; 
        fAdc[i]  = 0.0; 
        fTime[i] = 0.0; 
    }
}

// ____________________________________________________________________________
Float_t CEPRawADCellBuffer::GetADMultiplicity(Int_t i) const
{
    return std::abs(i)<fNCells ? fMult[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawADCellBuffer::GetADCharge(Int_t i) const
{
    return std::abs(i)<fNCells ? fAdc[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawADCellBuffer::GetADTime(Int_t i) const
{
    return std::abs(i)<fNCells ? fTime[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
void CEPRawADCellBuffer::SetADVariables(AliESDAD* ADobj)
{
    // Global member variable setter
    this->SetADDecision_A(ADobj->GetADADecision());
    this->SetADDecision_C(ADobj->GetADCDecision());
    for (UInt_t i(0); i<fNCells; i++) 
    {  
        this->SetADMultiplicity(ADobj->GetMultiplicity(i), i);
        this->SetADCharge(ADobj->GetAdc(i), i);
        this->SetADTime(ADobj->GetTime(i), i);
    } 
}

// ____________________________________________________________________________
