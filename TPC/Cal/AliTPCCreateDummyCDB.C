#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <TRandom.h>
#include <TClonesArray.h>
#include <TFile.h>


#include <AliAlignObj.h>
#include <AliAlignObjAngles.h>
#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>


#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalDet.h"
#include "AliTPCParamSR.h"



#endif

//
// run number for the dummy file
Int_t gkDummyRun = 0;
char *gCDBpath   = "local://~/mycalib1";
AliCDBStorage* gStorLoc = 0;

//
Float_t gTSample    = 2.0000000e-07;
//
//
Float_t gMeanGain   = 1;
Float_t gSigmaGain  = 0;
//
Float_t gMeanTime0  = 0;
Float_t gSigmaTime0 = 0;
//
Float_t gMeanNoise  = 1;
Float_t gSigmaNoise = 0;
//
Float_t gMeanPRF    = 1;
Float_t gSigmaPRF   = 0;
//
Float_t gMeanPedestal  = 0;
Float_t gSigmaPedestal = 0;
//
// Missalignment
//
Float_t gSigmaDx    = 0;
Float_t gSigmaDy    = 0;
Float_t gSigmaDz    = 0;
Float_t gSigmaAngle = 0;




TObject* CreatePadObject(const char* shortName, const char* description, Float_t value, Float_t sigma=0)
{
  AliTPCCalPad *calPad = new AliTPCCalPad(shortName, description);
  for (UInt_t det=0; det<AliTPCROC::Instance()->GetNSectors(); det++){
    AliTPCCalROC *calROC = calPad->GetCalROC(det);
    for (UInt_t channel=0; channel<calROC->GetNchannels(); channel++){
      Float_t rvalue = gRandom->Gaus(value,sigma);
      calROC->SetValue(channel, rvalue);
    }
  }
  return calPad;
}


void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{
  AliCDBId id1(cdbPath, gkDummyRun, gkDummyRun); 
  gStorLoc->Put(object, id1, metaData); 
}
    

AliCDBMetaData* CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Marian Ivanov");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-06-00"); //root version
  md1->SetComment("The dummy values in this calibration file are for testing only");
  
  return md1;
}


void CDBAlignmentObjectCreation(const char *fileName, const char *arrayName, const char *detName){
  // make instance of storage
  AliCDBManager *CDB = AliCDBManager::Instance();
  AliCDBStorage* storLoc = CDB->GetStorage(gCDBpath);
  
  // create or get from a file the TClonesArray of alignment objects
  // for given detector, DET should be TPC, TRD ...
  TFile* f = TFile::Open(fileName,"READ");
  TClonesArray* tc = ((TClonesArray*) f->Get(arrayName));
  
  // AliCDBStorage::Put will make a file containing that TClonesArray in
  // the CDB directory structure
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetObjectClassName("TClonesArray");
  md->SetResponsible("Your name here");
  md->SetComment("Alignment objects storage example. Please give more details here.");
  AliCDBId id(Form("%s/Align/Data",detName),0,0); //you have to specify the run validity,
// although in the case of saving ideal objects makes not much sense
  storLoc->Put(tc,id,md);
  
  // In this way it's good enough because AliSimulation is then expecting
  // as default to find the files in $ALICE_ROOT/DET/Align/Data
  // Then just
  // AliSimulation sim
  // sim.Run()
  // Look at the AliInfo messages from AliSimulation::Run() to check that
  // the expected CDB objects have been retrieved.
}




void GenerateRndTPC(Float_t sigmatrx=0., Float_t sigmatry=0, Float_t sigmatrz=0, Float_t sigmarot = 0.){

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10000);
  TClonesArray &alobj = *array;
  
  TFile f("TPC_alignment.root","RECREATE");

  TRandom *rnd   = new TRandom(4357);
  AliAlignObjAngles o;
  Int_t j = 0;
  for (Int_t iLayer = AliAlignObj::kTPC1; iLayer <= AliAlignObj::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {

      Float_t dx = (rnd->Uniform()-0.5)*sigmatrx;
      Float_t dy = (rnd->Uniform()-0.5)*sigmatry;
      Float_t dz = (rnd->Uniform()-0.5)*sigmatrz;
      Float_t dpsi = (rnd->Uniform()-0.5)*sigmarot;
      Float_t dtheta = (rnd->Uniform()-0.5)*sigmarot;
      Float_t dphi = (rnd->Uniform()-0.5)*sigmarot;

      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *path = AliAlignObj::GetVolPath(volid);
      new(alobj[j]) AliAlignObjAngles(path, volid, dx, dy, dz, dpsi, dtheta, dphi);
      j++;
    }
  }
  f.cd();
  f.WriteObject(array,"TPCAlignObjs","kSingleKey");
  f.Close();
  array->Delete();
  CDBAlignmentObjectCreation("TPC_alignment.root","TPCAlignObjs","TPC");
}



void AliTPCCreateDummyCDB()
{
  cout << endl << "TPC :: Creating dummy CDB with event number " << gkDummyRun << endl;
  

  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage(gCDBpath);
  if (!gStorLoc)
    return;

  TObject* obj = 0;
  AliCDBMetaData* metaData = 0;
  
  //
  // Gain factor (relative) - normalized to 1 - spread 0
  //
  metaData = CreateMetaObject("AliTPCCalPad");  
  obj = CreatePadObject("PadGainFactor","TPC Gain Factor (local -pad- variations)", gMeanGain , gSigmaGain);
  StoreObject("TPC/Calib/PadGainFactor", obj, metaData);
  //
  // Time0 fluctuation   - normalized to 0  - spread 0.00 mus
  //
  metaData = CreateMetaObject("AliTPCCalPad");  
  obj = CreatePadObject("PadTime0","TPC Time 0  (local -pad- variations)", gMeanTime0 , gSigmaTime0);
  StoreObject("TPC/Calib/PadTime0", obj, metaData);
  //
  // Noise  fluctuation   - normalized to 1.0  - spread - 0.0 
  //
  metaData = CreateMetaObject("AliTPCCalPad");  
  obj = CreatePadObject("PadNoise","TPC Noise  (local -pad- variations)", gMeanNoise , gSigmaNoise);
  StoreObject("TPC/Calib/PadNoise", obj, metaData);
  //
  // PRF width fluctuation   - normalized to 0.  - spread - 0.0 
  //
  metaData = CreateMetaObject("AliTPCCalPad");  
  obj = CreatePadObject("PadPRF","TPC PRF  (local -pad- variations)", gMeanPRF , gSigmaPRF);
  StoreObject("TPC/Calib/PadPRF", obj, metaData);
  //
  // Pedestals
  //
  metaData = CreateMetaObject("AliTPCCalPad");  
  obj = CreatePadObject("PadPedestal","TPC pedestals  (local -pad- variations)", gMeanPedestal, gSigmaPedestal);
  StoreObject("TPC/Calib/Pedestals", obj, metaData);  
  //
  // Parameters 
  //
  metaData = CreateMetaObject("AliTPCParam");  
  AliTPCParam * param = new AliTPCParamSR;
  param->SetTSample(gTSample);
  param->Update();
  StoreObject("TPC/Calib/Parameters", param, metaData); 
  //
  //
  // generate random missalignemnt
  //
  GenerateRndTPC(gSigmaDx,gSigmaDy,gSigmaDz,gSigmaAngle);
}
