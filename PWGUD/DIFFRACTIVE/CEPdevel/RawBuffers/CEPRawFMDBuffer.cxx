// ////////////////////////////////////////////////////////////////////////////
//
// CEP FMD buffer
//
// structure to hold raw FMD information
//
// ____________________________________________________________________________
#include "CEPTrackBuffer.h"
#include "CEPRawFMDBuffer.h"
#include <stdexcept>
#include <stdlib.h> 
#include <iostream>
#include "AliTriggerAnalysis.h"

ClassImp(CEPRawFMDBuffer)

// ____________________________________________________________________________
CEPRawFMDBuffer::CEPRawFMDBuffer()
  : TObject()
  , fNCells(140)
{
    this->Reset();
}

// ____________________________________________________________________________
void CEPRawFMDBuffer::Reset()
{
    fNCells     = 140;
    for (UInt_t i(0); i<fNCells; i++){
        fMult[i] = 0.0; 
    }
}

// ____________________________________________________________________________
void CEPRawFMDBuffer::SetFMDCellMultiplicity(Float_t mult, UInt_t i)
{
    try {
        if (i>=fNCells) throw std::out_of_range("FMD cell array: index out of range"); 
        fMult[i] = mult;
    } catch (const std::out_of_range& e) {
        printf("%s\n",e.what());
        exit(EXIT_FAILURE);
    }
}

// ____________________________________________________________________________
Float_t CEPRawFMDBuffer::GetFMDCellMultiplicity(UInt_t i) const
{
    return (i<fNCells) ? fMult[i] : CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
Float_t CEPRawFMDBuffer::GetFMDTotalMultiplicity() const
{
    Float_t totalMult(0.0);
    for (UInt_t i(0); i<fNCells; i++) {
        totalMult += fMult[i];
    }

    return totalMult;
}

// ____________________________________________________________________________
void CEPRawFMDBuffer::SetFMDVariables(AliESDFMD* FMDobj)
{
    UInt_t cellIdx(0);
    Float_t sumMult(0.0);
   
    AliTriggerAnalysis::AliceSide side;
    Int_t detFrom, detTo;

    // loop over A and C side
    for (Int_t ii(0); ii<2; ii++) {

	    if (ii == 0) side = AliTriggerAnalysis::kASide;
    	if (ii == 1) side = AliTriggerAnalysis::kCSide;

    	detFrom = (side == AliTriggerAnalysis::kASide) ? 1 : 3;
    	detTo   = (side == AliTriggerAnalysis::kASide) ? 2 : 3;

    	Float_t totalMult = 0;
    	for (UShort_t det=detFrom;det<=detTo;det++) {
      	    Int_t nRings = (det == 1 ? 1 : 2);
      	    for (UShort_t ir = 0; ir < nRings; ir++) {
                Char_t   ring = (ir == 0 ? 'I' : 'O');
                UShort_t nsec = (ir == 0 ? 20  : 40);
                UShort_t nstr = (ir == 0 ? 512 : 256);
                for (UShort_t sec =0; sec < nsec;  sec++) {
                    for (UShort_t strip = 0; strip < nstr; strip++) {
                        Float_t mult = FMDobj->Multiplicity(det,ring,sec,strip);
                        if (mult == AliESDFMD::kInvalidMult) continue;
                        sumMult += mult;
                    }
                    // fill the cell after sector
                    this->SetFMDCellMultiplicity(sumMult, cellIdx); 
                    cellIdx++;
                    sumMult=0.0;
                }
            }
    	}
    }

}

// ____________________________________________________________________________
