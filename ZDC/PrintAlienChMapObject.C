#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TFile.h>
#include "AliCDBEntry.h"
#include "AliCDBGrid.h"
#include "AliCDBId.h"
#include "AliCDBLocal.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBPath.h"
#include "AliCDBRunRange.h"
#include "AliCDBStorage.h"
#include "AliDCSValue.h"
#include "AliZDC.h"
#include "AliZDCv3.h"
#include "AliZDCPedestals.h"
#include "AliZDCEnCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliZDCChMap.h"
#include "AliZDCLaserCalib.h"
#include "AliZDCMBCalib.h"
#include "AliZDCDataDCS.h"

#endif


void PrintAlienChMapObject(TString histoFName="Run109690_999999999_v5_s0.root ",
	Int_t alienOCDB=1)
{
  TGrid::Connect("alien:",0,0,"t");

  if (alienOCDB==1)histoFName = histoFName.Prepend("alien:///alice/data/2010/OCDB/ZDC/Calib/ChMap/");
  else if(alienOCDB==2) histoFName = histoFName.Prepend("alien:///alice/simulation/2008/v4-15-Release/Ideal/ZDC/Calib/ChMap/");
  else if(alienOCDB==3) histoFName = histoFName.Prepend("alien:///alice/simulation/2008/v4-15-Release/Full/ZDC/Calib/ChMap/");
  else if(alienOCDB==4) histoFName = histoFName.Prepend("alien:///alice/simulation/2008/v4-15-Release/Residual/ZDC/Calib/ChMap/");
  printf(" \n\n Printing object %s from alien OCDB \n\n",histoFName.Data());
  
  TFile *f = TFile::Open(histoFName);
  f->ls();
  AliCDBEntry *entry = (AliCDBEntry*)f->Get("AliCDBEntry");
  AliZDCChMap *calibdata = dynamic_cast<AliZDCChMap*>  (entry->GetObject());
  calibdata->Print("");

}
