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


void AliGenEposEventHeader::FillInternalFields(TGenerator */*epos*/) {
/*
    fBimevt = epos->GetBimevt();
    fPhievt = epos->GetPhievt();
    fKolevt = epos->GetKolevt();
    fKoievt = epos->GetKoievt();
    fPmxevt = epos->GetPmxevt();
    fEgyevt = epos->GetEgyevt();
    fNpjevt = epos->GetNpjevt();
    fNtgevt = epos->GetNtgevt();
    fNpnevt = epos->GetNpnevt();
    fNppevt = epos->GetNppevt();
    fNtnevt = epos->GetNtnevt();
    fNtpevt = epos->GetNtpevt();
    fJpnevt = epos->GetJpnevt();
    fJppevt = epos->GetJppevt();
    fJtnevt = epos->GetJtnevt();
    fJtpevt = epos->GetJtpevt();
    fXbjevt = epos->GetXbjevt();
    fQsqevt = epos->GetQsqevt();
    fNglevt = epos->GetNglevt();
    fZppevt = epos->GetZppevt();
    fZptevt = epos->GetZptevt();
*/
}

Float_t AliGenEposEventHeader::GetBimevt() { return fBimevt; }
Float_t AliGenEposEventHeader::GetPhievt() { return fPhievt; }
Int_t AliGenEposEventHeader::GetKolevt() { return fKolevt; }
Int_t AliGenEposEventHeader::GetKoievt() { return fKoievt; }
Float_t AliGenEposEventHeader::GetPmxevt() { return fPmxevt; }
Float_t AliGenEposEventHeader::GetEgyevt() { return fEgyevt; }
Int_t AliGenEposEventHeader::GetNpjevt() { return fNpjevt; }
Int_t AliGenEposEventHeader::GetNtgevt() { return fNtgevt; }
Int_t AliGenEposEventHeader::GetNpnevt() { return fNpnevt; }
Int_t AliGenEposEventHeader::GetNppevt() { return fNppevt; }
Int_t AliGenEposEventHeader::GetNtnevt() { return fNtnevt; }
Int_t AliGenEposEventHeader::GetNtpevt() { return fNtpevt; }
Int_t AliGenEposEventHeader::GetJpnevt() { return fJpnevt; }
Int_t AliGenEposEventHeader::GetJppevt() { return fJppevt; }
Int_t AliGenEposEventHeader::GetJtnevt() { return fJtnevt; }
Int_t AliGenEposEventHeader::GetJtpevt() { return fJtpevt; }
Float_t AliGenEposEventHeader::GetXbjevt() { return fXbjevt; }
Float_t AliGenEposEventHeader::GetQsqevt() { return fQsqevt; }
Int_t AliGenEposEventHeader::GetNglevt() { return fNglevt; }
Float_t AliGenEposEventHeader::GetZppevt() { return fZppevt; }
Float_t AliGenEposEventHeader::GetZptevt() { return fZptevt; }

