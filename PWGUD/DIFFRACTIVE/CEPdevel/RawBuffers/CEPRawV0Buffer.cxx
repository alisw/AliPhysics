/// ////////////////////////////////////////////////////////////////////////////
//
// CEP V0 buffer
//
// structure to hold information on cell level for V0
//
// ____________________________________________________________________________
#include "CEPTrackBuffer.h"
#include "CEPRawV0Buffer.h"

ClassImp(CEPRawV0Buffer)

// ____________________________________________________________________________
CEPRawV0Buffer::CEPRawV0Buffer()
  : TObject()
  , fNCells(64)
  , fNPMV0A(CEPTrackBuffer::kdumval)
  , fNPMV0C(CEPTrackBuffer::kdumval)
{
    this->Reset();
}

// ____________________________________________________________________________
CEPRawV0Buffer::CEPRawV0Buffer(const CEPRawV0Buffer& vb)
  : TObject(vb)
  , fNCells(vb.fNCells)
  , fNPMV0A(vb.fNPMV0A)
  , fNPMV0C(vb.fNPMV0C)
{
    for (UInt_t i(0); i<fNCells; i++){
        fMult[i]     = vb.fMult[i]; 
        fAdc[i]      = vb.fAdc[i]; 
        fTime[i]     = vb.fTime[i]; 
        fSigWidth[i] = vb.fSigWidth[i]; 
    }
}

// ____________________________________________________________________________
CEPRawV0Buffer::~CEPRawV0Buffer()
{
  
}

// ____________________________________________________________________________
CEPRawV0Buffer& CEPRawV0Buffer::operator=(const CEPRawV0Buffer& source)
{
    if(this==&source) return *this;
    TObject::operator=(source);
    // Assignment operator
    fNCells = source.fNCells;
    fNPMV0A = source.fNPMV0A;
    fNPMV0C = source.fNPMV0C;

    for(UInt_t i(0); i<fNCells; i++) {
       fMult[i]     = source.fMult[i];
       fAdc[i]      = source.fAdc[i];
       fTime[i]     = source.fTime[i];
       fSigWidth[i] = source.fSigWidth[i];
    }

    return *this;
}

// ____________________________________________________________________________
void CEPRawV0Buffer::Copy(TObject &obj) const 
{
    if (this==&obj) return;
    CEPRawV0Buffer *robj = dynamic_cast<CEPRawV0Buffer*>(&obj);
    if (!robj) return; 
    *robj = *this;
}

// ____________________________________________________________________________
void CEPRawV0Buffer::Reset()
{
    fNCells   = 64;
    fNPMV0A   = CEPTrackBuffer::kdumval;
    fNPMV0C   = CEPTrackBuffer::kdumval;
    for (UInt_t i(0); i<fNCells; i++){
        fMult[i]     = 0.0; 
        fAdc[i]      = 0.0; 
        fTime[i]     = 0.0; 
        fSigWidth[i] = 0.0; 
    }
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0Multiplicity(Int_t i) const
{
    return std::abs(i)<fNCells ? fMult[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0TotalMultiplicity() const
{
    Float_t totalMult(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalMult += fMult[i]; 
    }

    return totalMult;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0Charge(Int_t i) const
{
    return std::abs(i)<fNCells ? fAdc[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0TotalCharge() const
{
    Float_t totalCharge(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalCharge += fAdc[i]; 
    }

    return totalCharge;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0Time(Int_t i) const
{
    return std::abs(i)<fNCells ? fTime[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0TotalTime() const
{
    Float_t totalTime(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalTime += fTime[i]; 
    }

    return totalTime;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0Width(Int_t i) const
{
    return std::abs(i)<fNCells ? fSigWidth[std::abs(i)] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawV0Buffer::GetV0TotalSignalWidth() const
{
    Float_t totalWidth(0.0);
    for (UInt_t i(0); i<fNCells; i++){
        totalWidth += fSigWidth[i]; 
    }

    return totalWidth;
}

// ____________________________________________________________________________
void CEPRawV0Buffer::SetV0Variables(AliESDVZERO* V0obj)
{
    // Global member variable setter
    this->SetV0NPMV0A(V0obj->GetNbPMV0A());
    this->SetV0NPMV0C(V0obj->GetNbPMV0C());
    for (UInt_t i(0); i<fNCells; i++) 
    {  
        this->SetV0Multiplicity(V0obj->GetMultiplicity(i), i);
        this->SetV0Charge(V0obj->GetAdc(i), i);
        this->SetV0Time(V0obj->GetTime(i), i);
        this->SetV0Width(V0obj->GetWidth(i), i);
    } 
}

// ____________________________________________________________________________
