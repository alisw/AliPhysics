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
/* History of cvs commits:
 *
* $Log$
* Revision 1.3  2007/10/16 14:36:39  pavlinov
* fixed code violation (almost)
*
* Revision 1.2  2007/09/11 19:38:15  pavlinov
* added pi0 calibration, linearity, shower profile
* co: warning: `$Log' is obsolescent; use ` * $Log'.
* Revision 1.1  2007/08/08 15:58:01  hristov
*/


//_________________________________________________________________________
//   Calibration coefficients
//   The place for holding all information after calibration
//   - not only calibration coefficients.
//   For flash Adc we should keep at least ratio of high/low gains 
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
#include <TString.h>

ClassImp(AliEMCALCalibCoef) 
// ----------------------------------------------------------------------
AliEMCALCalibCoef::AliEMCALCalibCoef() : fAbsId(-1), fCc(-1), fECc(-1)
{ 
  // default constructor
}

// ----------------------------------------------------------------------
AliEMCALCalibCoef::AliEMCALCalibCoef(const Int_t id, const Double_t c, const Double_t ec) :
fAbsId(id), fCc(c), fECc(ec)
{
  // Oct 16, 2007
}


ClassImp(AliEMCALCalibCoefs)

// ----------------------------------------------------------------------
AliEMCALCalibCoefs::AliEMCALCalibCoefs() : TNamed("",""), fTable(0), fCurrentInd(0), fCalibMethod(0)
{
  // default constructor
}

// ----------------------------------------------------------------------
AliEMCALCalibCoefs::AliEMCALCalibCoefs(const AliEMCALCalibCoefs& coefs) 
  : TNamed(coefs), fTable(coefs.fTable), 
    fCurrentInd(coefs.fCurrentInd), fCalibMethod(coefs.fCalibMethod)
{
  // copy constructor
}

// ----------------------------------------------------------------------
AliEMCALCalibCoefs::AliEMCALCalibCoefs(const char* name, const Int_t nrow) : TNamed(name,"table of cell information") , fTable(0), fCurrentInd(0), fCalibMethod(0)
{
  // Oct 16, 2007
  fTable = new TObjArray(nrow);
}

// ----------------------------------------------------------------------
void AliEMCALCalibCoefs::AddAt(AliEMCALCalibCoef* r)
{
  // Oct 16, 2007
  (*fTable)[fCurrentInd] = new AliEMCALCalibCoef(*r);
  fCurrentInd++;
}

// ----------------------------------------------------------------------
AliEMCALCalibCoefs::~AliEMCALCalibCoefs()
{
  // Destructor
  if(fTable) {
    fTable->Delete();
    delete fTable;
  }
}

// ----------------------------------------------------------------------
AliEMCALCalibCoef* AliEMCALCalibCoefs::GetTable(Int_t i) const
{
  // Oct 16, 2007
  return (AliEMCALCalibCoef*)fTable->At(i);
}

// ----------------------------------------------------------------------
// Get initial Calib Data from DB
AliEMCALCalibCoefs* AliEMCALCalibCoefs::GetCalibTableFromDb(const char *tn, AliEMCALCalibData **calData)
{ 
  //
  // See ~/macros/ALICE/sim.C  for choice of CDB
  // Get calib. table which was used than calculated rec.points
  // 
  static const char *calibType = "EMCAL/Calib/*";
  static const char *calibTypeData = "EMCAL/Calib/Data";
  // Initial cc
  calData[0] = 0;

  AliCDBManager* man = AliCDBManager::Instance(); 
  AliCDBStorage* specificStorage = man->GetSpecificStorage(calibType);

  AliEMCALCalibData* caldata = (AliEMCALCalibData*)
    specificStorage->Get(calibTypeData, AliRunLoader::Instance()->GetRunNumber())->GetObject();
  if(caldata == 0) return 0;

  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();

  AliEMCALCalibCoefs *tab = new  AliEMCALCalibCoefs(tn,g->GetNCells());
  tab->SetCalibMethod(AliEMCALCalibCoefs::kMC);

  for(Int_t id=0; id<g->GetNCells(); id++){
    Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, iphiCell=0, ietaCell=0;
    g->GetCellIndex(id, nSupMod, nModule, nIphi, nIeta);
    g->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphiCell, ietaCell);

    AliEMCALCalibCoef r;
    r.fAbsId = id;
    r.fCc    = caldata->GetADCchannel(nSupMod, ietaCell, iphiCell); // column(ietaCell) : row(iphiCell)
    r.fECc   = 0.0;

    tab->AddAt(&r);
  }
  // tab->Purge();
  tab->SetCalibMethod(AliEMCALCalibCoefs::kPI0);
  calData[0] = caldata; 

  printf("\n <I> AliEMCALCalibCoefs::GetCalibTableFromDb \n name |%s| title |%s| | storage |%s|\n", 
	 caldata->GetName(), caldata->GetTitle(),  specificStorage->GetURI().Data());

  return tab;
}

// ----------------------------------------------------------------------
TH1F* AliEMCALCalibCoefs::GetHistOfCalibTableFromDb(const char *tn)
{
  // First SM only
  AliEMCALCalibData *calData[1]; // unused here

  AliEMCALCalibCoefs* tab = GetCalibTableFromDb(tn, calData);
  if(tab==0) return 0;

  TH1F *h = new TH1F("hCCfirst", " fCc first  (in MeV)", 70, 12., 19.);
  AliEMCALCalibCoef *r;
  for(Int_t i=0; i<tab->GetSize(); i++){
    r = tab->GetTable(i);
    if(i>=1152) break;
    h->Fill(r->fCc*1000.); // GEV ->MeV
  }
  delete tab;
  return h;
}

// ----------------------------------------------------------------------
AliEMCALCalibData* AliEMCALCalibCoefs::GetCalibTableForDb(const AliEMCALCalibCoefs *tab, const char* dbLocation,
const char* coment)
{
  // See EMCAL/macros/CalibrationDB/AliEMCALSetCDB.C
  // Define and save to CDB 
  printf("<I> AliEMCALCalibCoefs::GetCalibTableForDb started \n");
  AliEMCALCalibData* caldata=0;
  if(tab==0) return caldata;

  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  const char* objFormat  = "";
  caldata = new AliEMCALCalibData("EMCAL");
  caldata->SetTitle(coment);

  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();

  for(int id=0; id<tab->GetSize(); id++){
    Float_t ped=0.;

    Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, iphiCell=0, ietaCell=0;
    g->GetCellIndex(id, nSupMod, nModule, nIphi, nIeta);
    g->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphiCell, ietaCell);

    AliEMCALCalibCoef *r = tab->GetTable(id);
    // ietaCell - column; iphiCell - row
    caldata->SetADCchannel (nSupMod, ietaCell, iphiCell, r->fCc);
    caldata->SetADCpedestal(nSupMod, ietaCell, iphiCell, ped);
  }
  printf("<I> Fill AliEMCALCalibData table \n");
  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Aleksei Pavlinov");

  AliCDBManager* man = AliCDBManager::Instance();
  if(man == 0) {
    printf("<E> AliEMCALCalibCoefs::GetCalibTableForDb : define AliCDBManager, NO saving  !! \n"); 
  } else {
    //printf("<I> AliCDBManager %p \n", man); 
    AliCDBId id("EMCAL/Calib/Data",firstRun,lastRun); // create in EMCAL/Calib/Data DBFolder 
    TString dBFolder(dbLocation);
    AliCDBStorage* loc = man->GetStorage(dBFolder.Data());
    loc->Put(caldata, id, &md);
  }

  return caldata;
}

// ----------------------------------------------------------------------
AliEMCALCalibCoef *AliEMCALCalibCoefs::GetRow(const int absId)
{
  // Oct 16, 2007
  AliEMCALCalibCoef *r=0;
  for(int id=0; id<fTable->GetSize(); id++){
    r = GetTable(id);
    if(r->fAbsId == absId) return r;
  }
  return 0;
}

// ----------------------------------------------------------------------
void AliEMCALCalibCoefs::PrintTable()
{
  // Oct 16, 2007
  printf(" Table : %s : nrows %i \n", GetName(), int(fTable->GetSize()));
  for(int i=0; i<fTable->GetSize(); i++) PrintTable(i);
}

// ----------------------------------------------------------------------
void AliEMCALCalibCoefs::PrintTable(const Int_t i)
{
  // Oct 16, 2007
  if(i>=fTable->GetSize()) return;
  printf("row %i \n", i);
  PrintRec(GetTable(i));
}

// ----------------------------------------------------------------------
void AliEMCALCalibCoefs::PrintRec(AliEMCALCalibCoef* r)
{
  // Oct 16, 2007
  if(r==0) return;
  printf(" abs Id %5.5i fCc %7.6f fECc %7.6f \n", r->fAbsId, r->fCc, r->fECc);
}
