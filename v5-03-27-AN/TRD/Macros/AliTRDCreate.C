#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <TRandom.h>
#include <TSystem.h>
#include <TDatime.h>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>
#include <AliGeometry.h>
#include <AliPID.h>

#include "../TRD/AliTRDgeometry.h"

#include "../TRD/Cal/AliTRDCalROC.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalDet.h"

#endif

// run number for the dummy file
const Int_t gkDummyRun = 0;
AliCDBStorage* gStorLoc = 0;
Double_t t0minimum[540];
Double_t gainav[540];
Double_t vav[540];

TObject* CreatePadObjectT0Random(const char* shortName, const char* description, Float_t mean, Float_t sigma);
TObject* CreateDetT0Object(const char* shortName, const char* description);
Double_t funcpoly16(Double_t* x);
Double_t funcpoly12(Double_t* x);
Double_t funcpoly144(Double_t* x);
Double_t funcpoly16v(Double_t* x);
Double_t funcpoly12v(Double_t* x);
Double_t funcpoly144v(Double_t* x);
TObject* CreatePadObject(const char* shortName, const char* description, Float_t value);
TObject* CreatePadObjectRandom(const char* shortName, const char* description, Float_t mean, Float_t sigma);
TObject* CreatePadObjectRandomv(const char* shortName, const char* description, Float_t mean, Float_t sigma);
TObject* CreatePadObjectbridge(const char* shortName, const char* description);
TObject* CreatePadObjectbridgev(const char* shortName, const char* description);
TObject* CreateDetObject(const char* shortName, const char* description, Float_t value);
TObject* CreateDetObjectRandom(const char* shortName, const char* description, Float_t mean, Float_t sigma);
TObject* CreateDetObjectRandomg(const char* shortName, const char* description, Float_t mean, Float_t sigma);
TObject* CreateDetObjectRandomv(const char* shortName, const char* description, Float_t mean, Float_t sigma);
AliCDBMetaData* CreateMetaObject(const char* objectClassName);
void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
void AliTRDCreate(Bool_t residual = kFALSE);


//________________________________________________________________________________________________
TObject* CreatePadObjectT0Random(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  Double_t min = 0.0;
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      Double_t value[2400];
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	value[channel] = gRandom->Gaus(mean,sigma);
	if(channel == 0) {
	  min = value[0];
	  //printf("value[0] %f and min %f et gRandom %f\n",value[0],min, gRandom->Gaus(mean,sigma));
	}
	if(min > value[channel]) min = value[channel];
	//calROC->SetValue(channel, TMath::Abs(value));
      }
      t0minimum[det] = min;
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	//printf("min value for the det %d is %f and the channel is %f, we put %f\n",det,t0minimum[det],value[channel],(value[channel]-t0minimum[det]));
	calROC->SetValue(channel, (value[channel]-t0minimum[det]));
      }
    }
  
  return calPad;
}
//___________________________________________________________________________________________________
TObject* CreateDetT0Object(const char* shortName, const char* description)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    object->SetValue(det, t0minimum[det]);
  return object;
}
//__________________________________________________________________________________________________
Double_t funcpoly16(Double_t* x){

  //0.1

  Double_t par[3] = {-0.002678571429,0.040178571,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}
//_____________________________________________________________________________________________________
Double_t funcpoly12(Double_t* x){

  // 0.1

  Double_t par[3] = {-0.005,0.055,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}
//______________________________________________________________________________________________________
Double_t funcpoly144(Double_t* x){

  //0.1

  Double_t par[3] = {-0.00001956181536,0.00279339596,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}

//__________________________________________________________________________________________________
Double_t funcpoly16v(Double_t* x){

  //0.05 

  Double_t par[3] = {-0.0008928571429,0.013392857,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}
//_____________________________________________________________________________________________________
Double_t funcpoly12v(Double_t* x){

  //0.05

  Double_t par[3] = {-0.001666667,0.018333,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}
//___________________________________________________________________________________________________
Double_t funcpoly144v(Double_t* x){

  //0.05

  Double_t par[3] = {-0.000009780907668,0.001398669797,1};

  // sum landau + gaus with identical mean

  Double_t valLandau = par[0]*x[0]*x[0]+par[1]*x[0]+par[2];

  return valLandau;
}
//__________________________________________________________________________________________________
TObject* CreatePadObject(const char* shortName, const char* description, Float_t value)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	calROC->SetValue(channel, value);
      }
    }
  return calPad;
}
//___________________________________________________________________________________________________
TObject* CreatePadObjectRandom(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
    
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	Double_t value = gRandom->Gaus(mean,sigma);
	//cout << "value: " << value << endl;
	//if(value < 0) calROC->SetValue(channel, 0);
	calROC->SetValue(channel, value);
      }
    }
  
  return calPad;
}
//_______________________________________________________________________________________
TObject* CreatePadObjectRandomv(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
      vav[det]=0.0;
      Int_t nb = 0;
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	Double_t value = gRandom->Gaus(mean,sigma);
	vav[det] += value;
	nb++;
	//cout << "value: " << value << endl;
	//if(value < 0) calROC->SetValue(channel, 0);
	calROC->SetValue(channel, value);
      }
      vav[det] /= nb;
      for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel){
	Double_t value = calROC->GetValue(channel);
      	//cout << "value: " << value << endl;
	//if(value < 0) calROC->SetValue(channel, 0);
	calROC->SetValue(channel, value/vav[det]);
      }
    }
  
  return calPad;
}
//____________________________________________________________________________________________________
TObject* CreatePadObjectbridge(const char* shortName, const char* description)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  Double_t x[2];
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
  
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      Int_t rowMax = calROC->GetNrows();
      Int_t colMax = calROC->GetNcols();
      gainav[det] = 0.0;   

      // premier calcul
      for(Int_t row = 0; row < rowMax; ++row){
	for(Int_t col = 0; col < colMax; ++col){
	  x[0] = row;
	  x[1] = col;
	  if(rowMax == 12) {
	    Double_t value = funcpoly12(&x[0])*funcpoly144(&x[1]);
	    calROC->SetValue(col,row,value);
	    gainav[det] += value;
	  }
	  else {
	    Double_t value = funcpoly16(&x[0])*funcpoly144(&x[1]);
	    calROC->SetValue(col,row,value);
	    gainav[det] += value;
	  }
	}//col
      }//row
      gainav[det] /= (rowMax*colMax);

      //deuxieme calcul
      for(Int_t row = 0; row < rowMax; ++row){
	for(Int_t col = 0; col < colMax; ++col){
	  x[0] = row;
	  x[1] = col;
	  //Double_t value = funcpoly12(&x[0])*funcpoly144(&x[1]);
	  Double_t value = calROC->GetValue(col,row);
	  calROC->SetValue(col,row,value/gainav[det]);
	  //gainav[det] += value;
	}//col
      }//row

    }//det
  
  return calPad;
}
//__________________________________________________________________________________________
TObject* CreatePadObjectbridgev(const char* shortName, const char* description)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  Double_t x[2];
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    {
  
      AliTRDCalROC *calROC = calPad->GetCalROC(det);
      Int_t rowMax = calROC->GetNrows();
      Int_t colMax = calROC->GetNcols();
      vav[det] = 0.0;

      // premiere calcul
      for(Int_t row = 0; row < rowMax; ++row){
	for(Int_t col = 0; col < colMax; ++col){
	  x[0] = row;
	  x[1] = col;
	  if(rowMax == 12) {
	    Double_t value = funcpoly12v(&x[0])*funcpoly144v(&x[1]);
	    calROC->SetValue(col,row,value);
	    vav[det] += value;
	  }
	  else {
	    Double_t value = funcpoly16v(&x[0])*funcpoly144v(&x[1]);
	    calROC->SetValue(col,row,value);
	    vav[det] += value;
	  }
	}//col
      }//row
      vav[det] /= (rowMax*colMax);

      //deuxieme calcul
      for(Int_t row = 0; row < rowMax; ++row){
	for(Int_t col = 0; col < colMax; ++col){
	  x[0] = row;
	  x[1] = col;
	
	  //Double_t value = funcpoly12v(&x[0])*funcpoly144v(&x[1]);
	  Double_t value = calROC->GetValue(col,row);
	  calROC->SetValue(col,row,value/vav[det]);
	
	}//col
      }//row

    }//det
  
  return calPad;
}
//___________________________________________________________________________________________________
TObject* CreateDetObject(const char* shortName, const char* description, Float_t value)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    object->SetValue(det, value);
  return object;
}
//___________________________________________________________________________________________________
TObject* CreateDetObjectRandom(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det){
    Double_t value = gRandom->Gaus(mean,sigma);
    //cout << "value: " << value << endl;
    object->SetValue(det, value);
    //else cout << "value negative!" << endl;
  }
  return object;
}
//___________________________________________________________________________________________________
TObject* CreateDetObjectRandomg(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det){
    Double_t value = gRandom->Gaus(mean,sigma);
    //cout << "value: " << value << endl;
    object->SetValue(det, value/gainav[det]);
    //else cout << "value negative!" << endl;
  }
  return object;
}
//___________________________________________________________________________________________________
TObject* CreateDetObjectRandomv(const char* shortName, const char* description, Float_t mean, Float_t sigma)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det){
    Double_t value = gRandom->Gaus(mean,sigma);
    //cout << "value: " << value << endl;
    object->SetValue(det, value/vav[det]);
    //else cout << "value negative!" << endl;
  }
  return object;
}
//___________________________________________________________________________________________________
AliCDBMetaData* CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Raphaelle Bailhache");
  //md1->SetBeamPeriod(1);
  //md1->SetAliRootVersion("head"); //root version
  md1->SetComment("residual database TRD");
  
  return md1;
}
//___________________________________________________________________________________________________
void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{
  AliCDBId id1(cdbPath, gkDummyRun, 999999999); 
  id1.SetVersion(1);
  gStorLoc->Put(object, id1, metaData); 
}
//___________________________________________________________________________________________________
void AliTRDCreate(Bool_t residual)
{
 
  //************************random generator******************************************
  TDatime *datime = new TDatime();
  Int_t time = datime->GetTime();
  Int_t date = datime->GetDate();
  Int_t pid  = gSystem->GetPid();
  delete datime;
  Int_t seed = TMath::Abs(10000 * pid + time - date);
  gRandom->SetSeed(seed); 

  //*************************************************************************

  cout << endl << "TRD :: Creating dummy CDB with event number " << gkDummyRun << endl;
  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage("local://.");
  if (!gStorLoc)
    return;

  TObject* obj = 0;
  AliCDBMetaData* metaData = 0;

  if(!residual)
    {

      //Pad////////////////////////////////////////////////////////////////////

      metaData = CreateMetaObject("AliTRDCalPad");
      
      obj = CreatePadObjectT0Random("LocalT0","T0 (local variations)", 0, 0.2);
      StoreObject("TRD/Calib/LocalT0", obj, metaData);

      obj = CreatePadObjectbridge("LocalGainFactor","GainFactor (local variations)");
      StoreObject("TRD/Calib/LocalGainFactor", obj, metaData);

      obj = CreatePadObjectbridgev("LocalVdrift","TRD drift velocities (local variations)");
      StoreObject("TRD/Calib/LocalVdrift", obj, metaData);
      
      //Det//////////////////////////////////////////////////////////////////
      
      metaData = CreateMetaObject("AliTRDCalDet");
      
      obj = CreateDetObjectRandom("ChamberVdrift","TRD drift velocities (detector value)", 1.5, 0.08);
      StoreObject("TRD/Calib/ChamberVdrift", obj, metaData);

      obj = CreateDetT0Object("ChamberT0","T0 (detector value)");
      StoreObject("TRD/Calib/ChamberT0", obj, metaData);
      
      obj = CreateDetObjectRandom("ChamberGainFactor","GainFactor (detector value)", 1.0, 0.18);
      StoreObject("TRD/Calib/ChamberGainFactor", obj, metaData);

    }
  else
    {

      //Pad////////////////////////////////////////////////////////////////////
  
      metaData = CreateMetaObject("AliTRDCalPad");
      
      obj = CreatePadObjectRandomv("LocalVdrift","TRD drift velocities (local variations)", 1.5, 0.015);
      StoreObject("TRD/Calib/LocalVdrift", obj, metaData);
  
      obj = CreatePadObjectT0Random("LocalT0","T0 (local variations)", 0, 0.02);
      StoreObject("TRD/Calib/LocalT0", obj, metaData);

      obj = CreatePadObjectRandom("LocalGainFactor","GainFactor (local variations)", 1, 0.01);
      StoreObject("TRD/Calib/LocalGainFactor", obj, metaData);
  
      //Det//////////////////////////////////////////////////////////////////
      
      metaData = CreateMetaObject("AliTRDCalDet");
  
      obj = CreateDetT0Object("ChamberT0","T0 (detector value)");
      StoreObject("TRD/Calib/ChamberT0", obj, metaData);
      
    }  

}
