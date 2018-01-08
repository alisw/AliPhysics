// ////////////////////////////////////////////////////////////////////////////
//
// CEP AD buffer containing information on the AD cells
//
// structure to hold track information
//
// ____________________________________________________________________________
#include "CEPRawADBuffer.h"
#include "CEPTrackBuffer.h"

ClassImp(CEPRawADBuffer)

// ____________________________________________________________________________
CEPRawADBuffer::CEPRawADBuffer()
  : TObject()
  , fNCells(16)  
  , fDecisionADA (CEPTrackBuffer::kdumval)
  , fDecisionADC (CEPTrackBuffer::kdumval)
{
    this->Reset();
}

// ____________________________________________________________________________
CEPRawADBuffer::CEPRawADBuffer(const CEPRawADBuffer& adb)
  : TObject(adb)
  , fNCells(adb.fNCells)  
  , fDecisionADA (adb.fDecisionADA)
  , fDecisionADC (adb.fDecisionADC)
{
    for (UInt_t i(0); i<fNCells; i++) {
        fMult[i] = adb.fMult[i]; 
        fAdc[i]  = adb.fAdc[i]; 
        fTime[i] = adb.fTime[i]; 
    }
}

// ____________________________________________________________________________
CEPRawADBuffer& CEPRawADBuffer::operator=(const CEPRawADBuffer& source)
{
    if(this==&source) return *this;
    TObject::operator=(source);
    // Assignment operator
    fNCells      = source.fNCells;
    fDecisionADA = source.fDecisionADA;
    fDecisionADC = source.fDecisionADC;

    for(UInt_t i(0); i<fNCells; i++) {
       fMult[i] = source.fMult[i];
       fAdc[i]  = source.fAdc[i];
       fTime[i] = source.fTime[i];
    }

    return *this;
}

// ____________________________________________________________________________
void CEPRawADBuffer::Copy(TObject &obj) const 
{
    if (this==&obj) return;
    CEPRawADBuffer *robj = dynamic_cast<CEPRawADBuffer*>(&obj);
    if (!robj) return; 
    *robj = *this;
}

// ____________________________________________________________________________
CEPRawADBuffer::~CEPRawADBuffer()
{
  
}

// ____________________________________________________________________________
void CEPRawADBuffer::Reset()
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
Float_t CEPRawADBuffer::GetADMultiplicity(Int_t i) const
{
    return std::abs(i)<fNCells ? fMult[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawADBuffer::GetADTotalMultiplicity() const
{
    Float_t totalMult(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalMult += fMult[i]; 
    }

    return totalMult;
}

// ____________________________________________________________________________
Float_t CEPRawADBuffer::GetADCharge(Int_t i) const
{
    return std::abs(i)<fNCells ? fAdc[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawADBuffer::GetADTotalCharge() const
{
    Float_t totalCharge(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalCharge += fAdc[i]; 
    }

    return totalCharge;
}

// ____________________________________________________________________________
Float_t CEPRawADBuffer::GetADTime(Int_t i) const
{
    return std::abs(i)<fNCells ? fTime[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawADBuffer::GetADTotalTime() const
{
    Float_t totalTime(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalTime += fTime[i]; 
    }

    return totalTime;
}


// ____________________________________________________________________________
void CEPRawADBuffer::SetADVariables(AliESDAD* ADobj)
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
