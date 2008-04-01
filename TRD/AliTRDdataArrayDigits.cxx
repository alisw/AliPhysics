#include <TROOT.h>

#include "AliLog.h"
#include "AliTRDarrayS.h"
#include "Cal/AliTRDCalPadStatus.h"

#include "AliTRDdataArrayDigits.h"

ClassImp(AliTRDdataArrayDigits)

//_____________________________________________________________________________
AliTRDdataArrayDigits::AliTRDdataArrayDigits(Int_t nrow, Int_t ncol, Int_t ntime): AliTRDdataArrayS(nrow, ncol, ntime)
{
  //
  // Constructor
  //
}

//_____________________________________________________________________________
Short_t AliTRDdataArrayDigits::GetDataUnchecked(Int_t row, Int_t col, Int_t time) const
{
  // 
  // Get the pad status without boundary checking
  // (taking pad maskng into account)
  //
  Short_t value = GetDataFast(GetIdx1Unchecked(row,col),time);
  // Be aware of Manipulations introduced by pad masking in the RawReader
  // Only output the manipulated Value
  CLRBIT(value, 10);
  CLRBIT(value, 11);
  CLRBIT(value, 12);
  return value;
}

//_____________________________________________________________________________
Short_t AliTRDdataArrayDigits::GetData(Int_t row, Int_t col, Int_t time) const
{
  //
  // Get the pad signal
  // (taking pad masking into account)
  if ((row >= 0) && (col >= 0) && (time >= 0)) {
    Int_t idx1 = GetIdx1(row,col);
    if ((idx1 >= 0) && (time < fNdim2)) {
      Short_t value = 0;
      if (fBufType == 0)
	{
	  value = GetDataFast(idx1,time);
	}
      if (fBufType == 1) 
	{
	  value = GetData1(idx1,time);
	}
      // Awareness of the bit masking
      CLRBIT(value, 10);
      CLRBIT(value, 11);
      CLRBIT(value, 12);
      return value;
    }
    else {
      if (idx1 >= 0) {
        AliError(Form("time %d out of bounds (size: %d, this: 0x%08x)"
                     ,time,fNdim2,this));
      }
    }
  }

  return -1;
}


//_____________________________________________________________________________
Int_t AliTRDdataArrayDigits::GetOverThreshold(Short_t threshold)
{
  // 
  // Returns the number of entries over threshold
  // (taking pad masking into account)
  //
 
  if ((fElements == 0) || (fElements->GetSize() <= 0))
    return 0;
 
  Int_t over = 0;
  Short_t value = 0;

  for (Bool_t cont = First(); cont == kTRUE; cont = Next()) {
    if ((fCurrentIdx1 < 0) || (fCurrentIdx1 >= fNdim1)) continue;
    if ((fCurrentIdx2 < 0) || (fCurrentIdx2 >= fNdim2)) continue;
    value = fElements->At(fCurrentIndex);
    CLRBIT(value, 10);
    CLRBIT(value, 11);
    CLRBIT(value, 12);
    if (value > threshold) over++;
  }

  return over;

}

//_____________________________________________________________________________
UChar_t AliTRDdataArrayDigits::GetPadStatus(Int_t row, Int_t col, Int_t time) const
{
  // 
  // Returns the pad status stored in the pad signal
  //
  // Output is a UChar_t value
  // Status Codes:
  //               Noisy Masking:           2
  //               Bridged Left Masking     8
  //               Bridged Right Masking    8
  //               Not Connected Masking    32
  //
  UChar_t padstatus = 0;
  Short_t signal = GetDataFast(GetIdx1Unchecked(row,col),time);
  if(signal > 0 && TESTBIT(signal, 10)){
    if(TESTBIT(signal, 11))
      if(TESTBIT(signal, 12))
	padstatus = AliTRDCalPadStatus::kPadBridgedRight;
      else
	padstatus = AliTRDCalPadStatus::kNotConnected;
    else
      if(TESTBIT(signal, 12))
	padstatus = AliTRDCalPadStatus::kPadBridgedLeft;
      else
	padstatus = AliTRDCalPadStatus::kMasked;
  }
  return padstatus;
}

//_____________________________________________________________________________
void AliTRDdataArrayDigits::SetPadStatus(Int_t row, Int_t col, Int_t time, UChar_t status)
{
  //
  // Setting the pad status into the signal using the Bits 10 to 14 
  // (currently used: 10 to 12)
  //
  // Input codes (Unsigned char):
  //               Noisy Masking:           2
  //               Bridged Left Masking     8
  //               Bridged Right Masking    8
  //               Not Connected Masking    32
  //
  // Status codes: Any masking:             Bit 10(1)
  //               Noisy masking:           Bit 11(0), Bit 12(0)
  //               No Connection masking:   Bit 11(1), Bit 12(0)
  //               Bridged Left masking:    Bit 11(0), Bit 12(1)
  //               Bridged Right masking:   Bit 11(1), Bit 12(1)
  // 
  Short_t signal = GetDataFast(GetIdx1Unchecked(row, col), time);
  // only set the Pad Status if the signal is > 0
  if(signal > 0){
    switch(status)
      {
      case AliTRDCalPadStatus::kMasked:
	SETBIT(signal, 10);
	CLRBIT(signal, 11);
	CLRBIT(signal, 12);
	break;
      case AliTRDCalPadStatus::kNotConnected:
	SETBIT(signal, 10);
	SETBIT(signal, 11);
	CLRBIT(signal, 12);
      break;
      case AliTRDCalPadStatus::kPadBridgedLeft:
	SETBIT(signal, 10);
	CLRBIT(signal, 11);
	SETBIT(signal, 12);
	break;
      case AliTRDCalPadStatus::kPadBridgedRight:
	SETBIT(signal, 10);
	SETBIT(signal, 11);
	SETBIT(signal, 12);
      default:
	CLRBIT(signal, 10);
	CLRBIT(signal, 11);
	CLRBIT(signal, 12);
      }
    SetDataUnchecked(row, col, time, signal);
  }
}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayDigits::IsPadCorrupted(Int_t row, Int_t col, Int_t time)
{
  // 
  // Checks if the pad has any masking as corrupted (Bit 10 in signal set)
  //
  Short_t signal = GetDataFast(GetIdx1Unchecked(row, col), time);
  return (signal > 0 && TESTBIT(signal, 10)) ? kTRUE : kFALSE;
}

