#include "VZERO/AliVZEROConst.h"

void MakeVZEROCalibEntryPbPbNew(Int_t run, Int_t outputRun, const char *inputCDB = "raw://", const char *outputCDB = "local://$ALICE_ROOT/OCDB/VZERO/PbPb", Float_t mipperadc = 1./3.){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  calibdaorg->GetMIPperADC(0);
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  Float_t *a = calibdaorg->GetPMGainsA();
  Float_t *b = calibdaorg->GetPMGainsB();

  for (Int_t i = 0; i < 64; ++i) {
    Float_t nPhPerMIP = (i < 32) ? 6950 : 33690; 
    Float_t hv = TMath::Exp((TMath::Log(kChargePerADC/(mipperadc*nPhPerMIP*calibdaorg->GetLightYields(i)*0.18*TMath::Qe()))-a[i])/b[i]);
    calibda->SetMeanHV(hv,i);
  }

  for (Int_t i = 0; i < 64; ++i) {
    Float_t mip = (i < 32) ? 6950 : 33690; 
    printf("%d   %.1f %.1f %.1f   %.1f %.1f %.1f  %.1f  %f\n",
	   i,
	   calibdaorg->GetMeanHV(i),
	   1./calibdaorg->GetMIPperADC(i),
	   mip*calibdaorg->GetLightYields(i)*0.18*TMath::Qe()*calibdaorg->GetGain(i)/0.6e-12,
	   calibda->GetMeanHV(i),
	   1./calibda->GetMIPperADC(i),
	   mip*calibda->GetLightYields(i)*0.18*TMath::Qe()*calibda->GetGain(i)/0.6e-12,
	   calibda->GetMIPperADC(i)/calibdaorg->GetMIPperADC(i),
	   calibda->GetDiscriThr(i));
  }
  for (Int_t i = 0; i < 64; ++i) {
    printf("Channel = %d    %d V\n",
	   i,TMath::Nint(calibda->GetMeanHV(i)));
  }

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Pb-Pb VZERO Calibration from RAW OCDB");
  AliCDBId id("VZERO/Calib/Data",outputRun,outputRun);

  man->SetDefaultStorage(outputCDB);
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
