////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with only one step.                                    //
////////////////////////////////////////////////////////////

#include <TFile.h>
#include "AliITSOnlineSPDscanSingle.h"
#include "AliITSOnlineSPDscanInfo.h"

ClassImp(AliITSOnlineSPDscanSingle)

AliITSOnlineSPDscanSingle::AliITSOnlineSPDscanSingle(Char_t *fileName) {
  // constructor
  sprintf(fFileName,"%s",fileName);
  // look for a previously saved info object 
  // (if file not found create a new one and return, else read)
  FILE* fp0 = fopen(fFileName, "r");
  if (fp0 == NULL) {
    fScanInfo = new AliITSOnlineSPDscanInfo();
    fFile = new TFile(fFileName, "RECREATE");
    fWrite=kTRUE;
  }
  else {
    fclose(fp0);
    fFile = new TFile(fFileName, "READ");
    fWrite=kFALSE;
    fFile->GetObject("SPDscanInfo", fScanInfo);
  }
  Init();
  if (GetNSteps()==0) {
    AddScanStep(); // this is supposedly the only step for this object
  }
}
AliITSOnlineSPDscanSingle::~AliITSOnlineSPDscanSingle() {}

// call the corresponding methods in SPDscan with nsi=0 ******************
void AliITSOnlineSPDscanSingle::SetTriggers(UInt_t val)
  {AliITSOnlineSPDscan::SetTriggers(0,val);}
void AliITSOnlineSPDscanSingle::SetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val)
  {AliITSOnlineSPDscan::SetHits(0,hs,chipi,coli,rowi,val);}
void AliITSOnlineSPDscanSingle::IncrementTriggers()
  {AliITSOnlineSPDscan::IncrementTriggers(0);}
void AliITSOnlineSPDscanSingle::IncrementHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi)
  {AliITSOnlineSPDscan::IncrementHits(0,hs,chipi,coli,rowi);}
void AliITSOnlineSPDscanSingle::SetHitEvents(UInt_t hs, UInt_t chipi, UInt_t val)
  {AliITSOnlineSPDscan::SetHitEvents(0,hs,chipi,val);}
void AliITSOnlineSPDscanSingle::SetHitEventsTot(UInt_t hs, UInt_t val)
  {AliITSOnlineSPDscan::SetHitEventsTot(0,hs,val);}
void AliITSOnlineSPDscanSingle::IncrementHitEvents(UInt_t hs, UInt_t chipi)
  {AliITSOnlineSPDscan::IncrementHitEvents(0,hs,chipi);}
void AliITSOnlineSPDscanSingle::IncrementHitEventsTot(UInt_t hs)
  {AliITSOnlineSPDscan::IncrementHitEventsTot(0,hs);}
UInt_t AliITSOnlineSPDscanSingle::GetTriggers()
  {return AliITSOnlineSPDscan::GetTriggers(0);}
UInt_t AliITSOnlineSPDscanSingle::GetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi)
  {return AliITSOnlineSPDscan::GetHits(0,hs,chipi,coli,rowi);}
Float_t AliITSOnlineSPDscanSingle::GetHitsEfficiency(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi)
  {return AliITSOnlineSPDscan::GetHitsEfficiency(0,hs,chipi,coli,rowi);}
Float_t AliITSOnlineSPDscanSingle::GetHitsEfficiencyError(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi)
  {return AliITSOnlineSPDscan::GetHitsEfficiencyError(0,hs,chipi,coli,rowi);}
UInt_t AliITSOnlineSPDscanSingle::GetHitEvents(UInt_t hs, UInt_t chipi)
  {return AliITSOnlineSPDscan::GetHitEvents(0,hs,chipi);}
UInt_t AliITSOnlineSPDscanSingle::GetHitEventsTot(UInt_t hs)
  {return AliITSOnlineSPDscan::GetHitEventsTot(0,hs);}
Float_t AliITSOnlineSPDscanSingle::GetHitEventsEfficiency(UInt_t hs, UInt_t chipi)
  {return AliITSOnlineSPDscan::GetHitEventsEfficiency(0,hs,chipi);}
Float_t AliITSOnlineSPDscanSingle::GetHitEventsTotEfficiency(UInt_t hs)
  {return AliITSOnlineSPDscan::GetHitEventsTotEfficiency(0,hs);}
Float_t AliITSOnlineSPDscanSingle::GetHitEventsEfficiencyError(UInt_t hs, UInt_t chipi)
  {return AliITSOnlineSPDscan::GetHitEventsEfficiencyError(0,hs,chipi);}
Float_t AliITSOnlineSPDscanSingle::GetHitEventsTotEfficiencyError(UInt_t hs)
  {return AliITSOnlineSPDscan::GetHitEventsTotEfficiencyError(0,hs);}
Float_t AliITSOnlineSPDscanSingle::GetAverageMultiplicity(UInt_t hs, UInt_t chipi)
  {return AliITSOnlineSPDscan::GetAverageMultiplicity(0,hs,chipi);}
Float_t AliITSOnlineSPDscanSingle::GetAverageMultiplicityTot(UInt_t hs)
  {return AliITSOnlineSPDscan::GetAverageMultiplicityTot(0,hs);}
