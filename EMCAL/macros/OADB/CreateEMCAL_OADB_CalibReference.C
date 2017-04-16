///
/// \file CreateEMCAL_OADB_CalibReference.C
/// \ingroup EMCAL_OADB
/// \brief Create OADB file with reference calibration parameters from OCDB files
///
/// Create update the OADB for EMCAL reference calibration parameters with already existing
/// parameters in the OCDB.
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS ???
/// \author Marcel Figueredo, <marcel.figueredo@cern.ch>, Sao Paulo
///

#if !defined(__CINT__)
#include <TH2F.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TSystem.h>

#include "AliEMCALGeometry.h"
#include "AliEMCALCalibData.h"

#include "AliCDBEntry.h"
#include "AliOADBContainer.h"
#endif

///
/// Recover an existing histogram with recalibration factors  
/// in order to get consistent binning
///
/// \param run: reference run number
///
TH2F* GetHistogramExample(int run=153323)
{ 
  /*
   For 2012-13:  /alice/data/2012/OCDB/EMCAL/Calib/Data/Run177115_999999999_v2_s0.root
   For 2011:  /alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v3_s0.root
   For 2010: /alice/data/2010/OCDB/EMCAL/Calib/Data/Run113461_999999999_v3_s0.root
   */
  
  TString pathOADB = "$ALICE_ROOT/OADB/EMCAL";
  AliOADBContainer *contRF=new AliOADBContainer("");
  contRF->InitFromFile(Form("%s/EMCALRecalib.root",pathOADB.Data()),"AliEMCALRecalib");
  
  TH2F *h;
  TString pass = "pass1";
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(run);
  if(!recal)
  {
    printf("Energy recalibration OADB not found  3\n");
    return 0x0;
  }
  
  TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
  if(!recalpass)
  {
    printf("Energy recalibration OADB not found 2\n");
    return 0x0;
  }
  
  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if(!recalib)
  {
    printf("Energy recalibration OADB not found 1\n");
    return 0x0;
  }
  
  h = (TH2F*)recalib->FindObject("EMCALRecalFactors_SM1");
  h->SetName("hbuffer");
  h->Reset();
  
  return h;
}

///
/// Recover an existing histogram with recalibration factors  
/// in order to get consistent binning
///
/// \param run: reference run number
/// \param name: OCDB file name and path
/// \param arrname: Name inside OADB for array of SM histograms
/// \param divFactor:
///
TObjArray *GetArrayEnergy(Int_t run=153323,
                          TString name="Run113461_999999999_v3_s0.root",
                          TString arrname="array",Float_t divFactor=1)
{
  TString pathOADB = "$ALICE_ROOT/OADB/EMCAL";
  gSystem->Load("libOADB");  
  
  //gROOT->LoadMacro("LoadLibraries.C");
  //LoadLibraries();  
  
  // **** Loading the root files with Recalibration Factors:
  TFile* f	= new TFile(name);
  
  // Instantiate EMCAL geometry for the first time
  AliEMCALGeometry *geom;
  if     (run < 140000) geom = new AliEMCALGeometry("EMCAL_FIRSTYEARV1","EMCAL"); // 2010
  else if(run < 198000) geom = new AliEMCALGeometry("EMCAL_COMPLETEV1","EMCAL");  // 2011-2012-2013
  else                  geom = new AliEMCALGeometry("EMCAL_COMPLETE12SMV1_DCAL_8SM","EMCAL"); // Run2
  
  const Int_t nSM = static_cast<const Int_t> (geom->GetNumberOfSuperModules());
  TObjArray *array = new TObjArray(nSM);
  array->SetName(arrname);
  TH2F *hbuffer = GetHistogramExample(153323); //Just for the binning
  
  TH2F *h[12];
  hbuffer->SetName("hb");
  for(int i=0;i<nSM;i++)
  {
    h[i] = (TH2F*)hbuffer->Clone(Form("EMCALRecalFactors_SM%d",i));
    h[i]->SetTitle(Form("EMCALRecalFactors_SM%d",i));
  }
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  AliEMCALCalibData* cparam = (AliEMCALCalibData*) cdb->GetObject();	
  
  Int_t nDiff = 0;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t paramOCDB = -1;
    Float_t paramOADB = -1;
    
    if(cparam) 
    {
      paramOCDB = cparam->GetADCchannel(iSM,iCol,iRow);
      paramOADB = divFactor!=1?paramOCDB/divFactor:paramOCDB;
      //printf("\n OCDB=%f and OADB=%f",paramOCDB,paramOADB);
    }
    
    if(h[iSM]) h[iSM]->SetBinContent(iCol,iRow,paramOADB);
  }
  
  for(int i=0;i<nSM;i++)
  {
    array->Add(h[i]);
  }
  
  delete geom;
  delete cparam;
  delete cdb;
  
  f->Close();
  
  return array;  
}

/// 
/// Create OADB Container for EMCal Calibration reference objects
/// Main method
///
void CreateEMCAL_OADB_CalibReference()
{  
  TObjArray *array12 = GetArrayEnergy(178000,"Run177115_999999999_v2_s0.root","Recalib",1); //last one is division factor
  TObjArray *array10 = GetArrayEnergy(113461,"Run113461_999999999_v3_s0.root","Recalib",1); 
  TObjArray *array11 = GetArrayEnergy(153323,"Run144484_999999999_v3_s0.root","Recalib",1);
  
  // Creating Container 
  AliOADBContainer* con = new AliOADBContainer("AliEMCALRecalib");
  
  con->AddDefaultObject(*&array10);
  con->AddDefaultObject(*&array11);
  con->AddDefaultObject(*&array12);
  
  // Appending objects to their specific Run number
  con->AppendObject(*&array10,113461,139517);
  con->AppendObject(*&array11,144484,170593);
  con->AppendObject(*&array12,177115,197692);
  
  con->WriteToFile("EMCALCalibOnlineRef.root");
}
