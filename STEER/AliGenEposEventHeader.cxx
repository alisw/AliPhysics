/*
 * AliGenEposEventHeader.h
 * 
 * Header for EPOS generated event.
 *
 *      Author: Piotr Ostrowski
 */

#include "AliGenEposEventHeader.h"
#include "TGenerator.h"

ClassImp(AliGenEposEventHeader)

AliGenEposEventHeader::AliGenEposEventHeader(const char* name):
    AliGenEventHeader(name),
    fBimevt(0),
    fPhievt(0),
    fKolevt(0),
    fKoievt(0),
    fPmxevt(0),
    fEgyevt(0),
    fNpjevt(0),
    fNtgevt(0),
    fNpnevt(0),
    fNppevt(0),
    fNtnevt(0),
    fNtpevt(0),
    fJpnevt(0),
    fJppevt(0),
    fJtnevt(0),
    fJtpevt(0),
    fXbjevt(0),
    fQsqevt(0),
    fNglevt(0),
    fZppevt(0),
    fZptevt(0)
{

}

AliGenEposEventHeader::AliGenEposEventHeader() :     fBimevt(0),
    fPhievt(0),
    fKolevt(0),
    fKoievt(0),
    fPmxevt(0),
    fEgyevt(0),
    fNpjevt(0),
    fNtgevt(0),
    fNpnevt(0),
    fNppevt(0),
    fNtnevt(0),
    fNtpevt(0),
    fJpnevt(0),
    fJppevt(0),
    fJtnevt(0),
    fJtpevt(0),
    fXbjevt(0),
    fQsqevt(0),
    fNglevt(0),
    fZppevt(0),
    fZptevt(0)
{

}


Float_t AliGenEposEventHeader::GetBimevt() { return fBimevt; }
Float_t AliGenEposEventHeader::GetPhievt() { return fPhievt; }
Int_t AliGenEposEventHeader::GetKolevt()   { return fKolevt; }
Int_t AliGenEposEventHeader::GetKoievt()   { return fKoievt; }
Float_t AliGenEposEventHeader::GetPmxevt() { return fPmxevt; }
Float_t AliGenEposEventHeader::GetEgyevt() { return fEgyevt; }
Int_t AliGenEposEventHeader::GetNpjevt()   { return fNpjevt; }
Int_t AliGenEposEventHeader::GetNtgevt()   { return fNtgevt; }
Int_t AliGenEposEventHeader::GetNpnevt()   { return fNpnevt; }
Int_t AliGenEposEventHeader::GetNppevt()   { return fNppevt; }
Int_t AliGenEposEventHeader::GetNtnevt()   { return fNtnevt; }
Int_t AliGenEposEventHeader::GetNtpevt()   { return fNtpevt; }
Int_t AliGenEposEventHeader::GetJpnevt()   { return fJpnevt; }
Int_t AliGenEposEventHeader::GetJppevt()   { return fJppevt; }
Int_t AliGenEposEventHeader::GetJtnevt()   { return fJtnevt; }
Int_t AliGenEposEventHeader::GetJtpevt()   { return fJtpevt; }
Float_t AliGenEposEventHeader::GetXbjevt() { return fXbjevt; }
Float_t AliGenEposEventHeader::GetQsqevt() { return fQsqevt; }
Int_t AliGenEposEventHeader::GetNglevt()   { return fNglevt; }
Float_t AliGenEposEventHeader::GetZppevt() { return fZppevt; }
Float_t AliGenEposEventHeader::GetZptevt() { return fZptevt; }

void AliGenEposEventHeader::SetBimevt(Float_t value) { fBimevt = value; }
void AliGenEposEventHeader::SetPhievt(Float_t value) { fPhievt = value; }
void AliGenEposEventHeader::SetKolevt(Int_t value)   { fKolevt = value; }
void AliGenEposEventHeader::SetKoievt(Int_t value)   { fKoievt = value; }
void AliGenEposEventHeader::SetPmxevt(Float_t value) { fPmxevt = value; }
void AliGenEposEventHeader::SetEgyevt(Float_t value) { fEgyevt = value; }
void AliGenEposEventHeader::SetNpjevt(Int_t value)   { fNpjevt = value; }
void AliGenEposEventHeader::SetNtgevt(Int_t value)   { fNtgevt = value; }
void AliGenEposEventHeader::SetNpnevt(Int_t value)   { fNpnevt = value; }
void AliGenEposEventHeader::SetNppevt(Int_t value)   { fNppevt = value; }
void AliGenEposEventHeader::SetNtnevt(Int_t value)   { fNtnevt = value; }
void AliGenEposEventHeader::SetNtpevt(Int_t value)   { fNtpevt = value; }
void AliGenEposEventHeader::SetJpnevt(Int_t value)   { fJpnevt = value; }
void AliGenEposEventHeader::SetJppevt(Int_t value)   { fJppevt = value; }
void AliGenEposEventHeader::SetJtnevt(Int_t value)   { fJtnevt = value; }
void AliGenEposEventHeader::SetJtpevt(Int_t value)   { fJtpevt = value; }
void AliGenEposEventHeader::SetXbjevt(Float_t value) { fXbjevt = value; }
void AliGenEposEventHeader::SetQsqevt(Float_t value) { fQsqevt = value; }
void AliGenEposEventHeader::SetNglevt(Int_t value)   { fNglevt = value; }
void AliGenEposEventHeader::SetZppevt(Float_t value) { fZppevt = value; }
void AliGenEposEventHeader::SetZptevt(Float_t value) { fZptevt = value; }

