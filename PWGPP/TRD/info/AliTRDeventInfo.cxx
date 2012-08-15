////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Event info for TRD performance train                                  //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TMath.h"

#include "AliESDHeader.h"
#include "AliESDRun.h"

#include "AliTRDeventInfo.h"

ClassImp(AliTRDeventInfo)
Int_t const AliTRDeventInfo::fgkMultBin[AliTRDeventInfo::kCentralityClasses-1] = {
  700, 1400, 2100, 2800/*, 3500*/
};
Float_t const AliTRDeventInfo::fgkCentBin[AliTRDeventInfo::kCentralityClasses-1] = {
  0.1, 0.2, 0.5, 0.8/*, 1.*/
};
//____________________________________________________________________
AliTRDeventInfo::AliTRDeventInfo():
  TObject()
  ,fHeader(NULL)
  ,fRun(NULL)
  ,fCentrality(-1)
  ,fMult(-1)
{
  //
  // Default Constructor
  // 
  SetBit(kOwner, 0);
}

//____________________________________________________________________
AliTRDeventInfo::AliTRDeventInfo(AliESDHeader *header, AliESDRun *run):
  TObject()
  ,fHeader(header)
  ,fRun(run)
  ,fCentrality(-1)
  ,fMult(-1)
{
  //
  // Constructor with Arguments
  //
  SetBit(kOwner, 0);
//  fHeader->Print();
/*  for(Int_t ilevel(0); ilevel<3; ilevel++){
    printf("L%d :: ", ilevel);
    Int_t itrig(0); TString tn;
    do{
      tn = fHeader->GetTriggerInputName(itrig++, ilevel);
      printf("%s ", tn.Data());
    } while(tn.CompareTo(""));
    printf("\n");
  }*/
}

//____________________________________________________________________
AliTRDeventInfo::AliTRDeventInfo(const AliTRDeventInfo &info):
  TObject()
  ,fHeader(info.fHeader)
  ,fRun(info.fRun)
  ,fCentrality(info.fCentrality)
  ,fMult(info.fMult)
{
  //
  // Copy Constructor
  // Flat Copy
  // 
  SetBit(kOwner, 0);
}

//____________________________________________________________________
AliTRDeventInfo& AliTRDeventInfo::operator=(const AliTRDeventInfo& info){
  //
  // Operator=
  // Flat Copy
  //
  if(this == &info) return *this;
  fHeader     = info.fHeader;
  fRun        = info.fRun;
  fCentrality = info.fCentrality;
  fMult       = info.fMult;
  SetBit(kOwner, 0);
  return *this;
}

//____________________________________________________________________
AliTRDeventInfo::~AliTRDeventInfo(){
  //
  // Destructor
  // Delete the entries if it is the Owner
  //
  Delete("");
}

//____________________________________________________________________
void AliTRDeventInfo::Delete(const Option_t *){
  //
  // Delete the Object
  // Delete the entries if it is the Owner
  // 
  if(IsOwner()){
    if(fHeader) delete fHeader;
    if(fRun) delete fRun;
  };
  fHeader = NULL;
  fRun = NULL;
}

//____________________________________________________________________
void AliTRDeventInfo::SetOwner()
{
  // Do deep copy
  
  SetBit(kOwner, 1);
  fHeader = new AliESDHeader(*fHeader);
  fRun = new AliESDRun(*fRun);
}

//____________________________________________________________________
UShort_t  AliTRDeventInfo::GetBunchFill() const
{
  // wrapper
  return fHeader->GetBunchCrossNumber();
}

//____________________________________________________________________
Int_t  AliTRDeventInfo::GetCentralityBin(Float_t cenPer)
{
// calculate centrality bin
  for(Int_t icen(0); icen<kCentralityClasses-1; icen++){
    if(cenPer<fgkCentBin[icen]) return icen;
  }
  return -1;
}

//____________________________________________________________________
TString  AliTRDeventInfo::GetFiredTriggerClasses()
{
// copy of similar function from AliESDEvent
  return (fRun&&fHeader)?fRun->GetFiredTriggerClasses(fHeader->GetTriggerMask()):"";
}

//____________________________________________________________________
Int_t  AliTRDeventInfo::GetMultiplicityBin(Int_t n)
{
// calculate centrality bin
  for(Int_t im(0); im<kCentralityClasses-1; im++){
    if(n<fgkMultBin[im]) return kCentralityClasses-im-1;
  }
  return 0;
}

//____________________________________________________________________
void AliTRDeventInfo::GetListOfIsolatedBunches(TH1D* hbc, Int_t bunchSpacing) {
  //
  // Find the isolated bunch crossings
  //
  // from I.Arsene first implementation in AliTRDcheckESD

  Int_t nBunches(0);
  Bool_t isIsolated[kLHCbunches]; memset(isIsolated, 0, kLHCbunches*sizeof(Bool_t));
  for(Int_t bcBin(1); bcBin<=hbc->GetNbinsX(); bcBin++) {
    Int_t bc(TMath::Nint(hbc->GetBinCenter(bcBin)));
    if(bc<0 || bc>=kLHCbunches) continue; // outside LHC range
    Double_t entries(hbc->GetBinContent(bcBin));
    if(entries<1.) continue;             // not filled
    
    // check isolation
    Bool_t kFOUND(kTRUE);
    for(Int_t ibc = TMath::Max(1,bcBin-bunchSpacing); ibc <= TMath::Min(Int_t(kLHCbunches), bcBin+bunchSpacing); ibc++) {
      if(ibc==bcBin) continue;
      if(hbc->GetBinContent(ibc)>0.) {
        kFOUND = kFALSE;
        break;
      }
    }
    isIsolated[bc] = kFOUND;
    if(kFOUND) nBunches++;
  }   // end loop over BC bins

  printf("Isolated bunches [%d] @ min bunch spacing [%d]:\n{", nBunches, bunchSpacing);
  for(Int_t ibc(0); ibc<kLHCbunches; ++ibc) if(isIsolated[ibc]) printf("%4d, ", ibc);
  printf("};\n");
}
