/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Log$ */

//_________________________________________________________________________
//    Calibration coefficients
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALCalibCoefs.h"

#include "AliRun.h"
#include "AliEMCALCalibData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliEMCALGeometry.h"

#include "TH1F.h"

TableClassImpl(AliEMCALCalibCoefs,calibCoef)

// Get initial Calib Data from DB
AliEMCALCalibCoefs* AliEMCALCalibCoefs::GetCalibTableFromDb(const char *tn)
{
  // Initial cc with decalibration
  char* dbString  = "local:///data/r22b/ALICE/PROD/CALIBRATION_May_2007/PI0/PDSF/10GEV/DECALIB/DeCalibDB"; 
  AliEMCALCalibData* caldata = (AliEMCALCalibData*)
    (AliCDBManager::Instance()->GetStorage(dbString)
     ->Get("EMCAL/Calib/Data",gAlice->GetRunNumber())->GetObject());
  if(caldata == 0) return 0;

  AliEMCALCalibCoefs *tab = new  AliEMCALCalibCoefs(tn);
  tab->SetCalibMethod(AliEMCALCalibCoefs::kMC);

  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();
  for(Int_t id=0; id<g->GetNCells(); id++){
    Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, iphiCell=0, ietaCell=0;
    g->GetCellIndex(id, nSupMod, nModule, nIphi, nIeta);
    g->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphiCell, ietaCell);

    calibCoef r;
    r.absId = id;
    r.cc    = caldata->GetADCchannel(nSupMod, ietaCell, iphiCell);
    r.eCc   = 0.0;

    tab->AddAt(&r);
  }
  tab->Purge();
  tab->SetCalibMethod(AliEMCALCalibCoefs::kPI0);

  return tab;
}

TH1F* AliEMCALCalibCoefs::GetHistOfCalibTableFromDb(const char *tn)
{
  // First SM only
  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance("");

  AliEMCALCalibCoefs* tab = GetCalibTableFromDb(tn);
  if(tab==0) return 0;

  TH1F *h = new TH1F("hCCfirst", " cc first  (in MeV)", 70, 12., 19.);
  calibCoef *r;
  for(Int_t i=0; i<tab->GetNRows(); i++){
    r = tab->GetTable(i);
    if(i>=1152) break;
    h->Fill(r->cc*1000.); 
  }
  delete tab;
  return h;
}

calibCoef *AliEMCALCalibCoefs::GetRow(const int absId)
{
  calibCoef *r=0;
  for(int id=0; id<GetNRows(); id++){
    r = GetTable(id);
    if(r->absId == absId) return r;
  }
  return 0;
}

void AliEMCALCalibCoefs::PrintTable()
{
  printf(" Table : %s : nrows %i \n", GetName(), int(GetNRows()));
  for(int i=0; i<GetNRows(); i++) PrintTable(i);
}

void AliEMCALCalibCoefs::PrintTable(const Int_t i)
{
  if(i>=GetNRows()) return;
  printf("row %i \n", i);
  PrintRec(GetTable(i));
}

void AliEMCALCalibCoefs::PrintRec(calibCoef* r)
{
  if(r==0) return;
  printf(" abs Id %5.5i cc %7.6f eCc %7.6f \n", r->absId, r->cc, r->eCc);
}
