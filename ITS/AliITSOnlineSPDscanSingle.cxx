/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with only one step.                                    //
////////////////////////////////////////////////////////////


#include <TFile.h>
#include "AliITSOnlineSPDscanSingle.h"
#include "AliITSOnlineSPDscanInfo.h"

AliITSOnlineSPDscanSingle::AliITSOnlineSPDscanSingle(const Char_t *fileName, Bool_t readFromGridFile) {
  // constructor
  fFileName=fileName;
  fModified=kFALSE;
  fInfoModified=kFALSE;
  // look for a previously saved info object 
  // (if file not found create a new one and return, else read)

  Bool_t bRead = readFromGridFile;

  if (!bRead) {
    FILE* fp0 = fopen(fFileName.Data(), "r");
    if (fp0 != NULL) {
      bRead=kTRUE;
      fclose(fp0);
    }
  }

  if (bRead) { // open file for reading
    fFile = TFile::Open(fFileName.Data(), "READ");
    if (fFile==NULL) { // grid file not found, create new local default file
      printf("ERROR: AliITSOnlineSPDscan: File %s not found! Creating 'test999.root' file instead\n",fFileName.Data());
      // create default empty file:
      fFileName = "test999.root";
      fScanInfo = new AliITSOnlineSPDscanInfo();
      fInfoModified=kTRUE;
      fFile = new TFile(fFileName.Data(), "RECREATE");
      fWrite=kTRUE;
    }
    else { // read from file (grid or local)
      fWrite=kFALSE;
      fFile->GetObject("AliITSOnlineSPDscanInfo", fScanInfo);
    }
  }
  else { // create new local file
    fScanInfo = new AliITSOnlineSPDscanInfo();
    fInfoModified=kTRUE;
    fFile = new TFile(fFileName.Data(), "RECREATE");
    fWrite=kTRUE;
  }

  Init();
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
