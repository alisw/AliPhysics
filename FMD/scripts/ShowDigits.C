//
// Script to digit multiplicity information to std::cout. 
//
#include <TH1F.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>

class ShowDigits : public AliFMDInputDigits
{
  TH1F* fHist;
  Int_t det
  ShowDigits(Int_t det, const char* file="galice.root") 
    : AliFMDInputDigits(file), fDet(det)
  {
    fHist = new TH1F("digitData", "Digit Data", 128, 0, 1024);
  }
  Bool_t ProcessDigit(AliFMDDigit* digit) 
  {
    if (digit->Counts() > 12) digit->Print();
    if (digit->Detector() == det) fHist->Fill(digit->Counts());
    return kTRUE;
  }
  Bool_t Finish() 
  {
    fHist->Draw();
    return kTRUE;
  }
}
//____________________________________________________________________
//
// EOF
//
