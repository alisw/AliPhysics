#include "ARVersion.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliESDVertex.h"
#include <TROOT.h>
#include "AliRun.h"
#include <TString.h>
#include "AliLog.h"
#endif

void UpdateCDBVertexDiamond(const char* objName, Double_t xmed, Double_t ymed, Double_t zmed, Double_t sigx, Double_t sigy, Double_t sigz) {
  // produce the mean vertex with the current AliRoot and store it in the
  // CDB
  // andrea.dainese@pd.infn.it
  
  AliCDBManager* man = AliCDBManager::Instance();
  if(!gSystem->Getenv("STORAGE")){
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }else{
      man->SetDefaultStorage(gSystem->Getenv("STORAGE"));
  }

  Int_t firstRun=0;
  if(gSystem->Getenv("FIRSTRUN")) firstRun=TString(gSystem->Getenv("FIRSTRUN")).Atoi();
  Int_t lastRun=AliCDBRunRange::Infinity();
  if(gSystem->Getenv("LASTRUN")) lastRun=TString(gSystem->Getenv("LASTRUN")).Atoi();
  TString objFullName="GRP/Calib/";

  if(objName=="MeanVertex" || objName=="MeanVertexSPD" || objName=="MeanVertexTPC"){
      objFullName += objName;
  }else{
     Printf("\"%s\" is not a valid calibration object. Exiting!",objName);
     return;
  }
  AliCDBId id(objFullName.Data(),firstRun,lastRun);
  AliCDBMetaData *metadata= new AliCDBMetaData();

  // Get root and AliRoot versions
  const char* rootv = gROOT->GetVersion();
  TString av(ALIROOT_VERSION);
  TString revnum(ALIROOT_REVISION);

  metadata->SetResponsible("prino@to.infn.it");
  metadata->SetComment("Default mean vertex position");
  metadata->SetAliRootVersion(av.Data());
  metadata->SetComment(Form("Default mean vertex produced with root version %s and AliRoot %s, revision number %s",rootv,av.Data(),revnum));
  

  Printf(Form("Storing in CDB the default mean vertex produced with root version %s and"
			  "AliRoot version %s, revision number %s", rootv, av.Data(), revnum));

  Double_t sigma[3],position[3];
  position[0]=xmed;
  position[1]=ymed;
  position[2]=zmed;
  sigma[0]=sigx;
  sigma[1]=sigy;
  sigma[2]=sigz;

  AliESDVertex *vertex = new AliESDVertex(position,sigma,"vtxmean");
  vertex->PrintStatus();

  man->Put(vertex,id,metadata);

}


void UpdateCDBVertexDiamondOld(Double_t xmed = 0., Double_t ymed = 0., Double_t sigx = 0.0060, Double_t sigy = 0.0060, Double_t sigz = 3.8) {
  // produce the mean vertex with the current AliRoot and store it in the
  // CDB
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  AliCDBId id("GRP/Calib/MeanVertex",0,AliCDBRunRange::Infinity());
  AliCDBId idTPC("GRP/Calib/MeanVertexTPC",0,AliCDBRunRange::Infinity());
  AliCDBId idSPD("GRP/Calib/MeanVertexSPD",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *metadata= new AliCDBMetaData();

  // Get root and AliRoot versions
  const char* rootv = gROOT->GetVersion();
  TString av(ALIROOT_VERSION);
  TString revnum(ALIROOT_REVISION);

  metadata->SetResponsible("prino@to.infn.it");
  metadata->SetComment("Default mean vertex position");
  metadata->SetAliRootVersion(av.Data());
  metadata->SetComment(Form("Default mean vertex produced with root version %s and AliRoot %s, revision number %s",rootv,av.Data(),revnum));
  

  Printf(Form("Storing in CDB the default mean vertex produced with root version %s and"
			  "AliRoot version %s, revision number %s", rootv, av.Data(), revnum));

  Double_t resolx=5./10000.; // this is error on the weighted mean (5 micron) 
  Double_t resoly=5./10000.; // this is error on the weighted mean (5 micron)
  Double_t sigma[3],position[3];
  position[0]=xmed;
  position[1]=ymed;
  position[2]=0.;
  sigma[0]=TMath::Sqrt(sigx*sigx+resolx*resolx);
  sigma[1]=TMath::Sqrt(sigy*sigy+resoly*resoly);
  sigma[2]=sigz;

  AliESDVertex *vertex = new AliESDVertex(position,sigma,"vtxmean");
  vertex->PrintStatus();

  man->Put(vertex,id,metadata);

  position[0]=xmed;
  position[1]=ymed;
  position[2]=0.;
  sigma[0]=TMath::Sqrt(sigx*sigx+resolx*resolx);
  sigma[1]=TMath::Sqrt(sigy*sigy+resoly*resoly);
  sigma[2]=sigz;

  AliESDVertex *vertexTPC = new AliESDVertex(position,sigma,"vtxmean");
  vertexTPC->PrintStatus();

  man->Put(vertexTPC,idTPC,metadata);

  position[0]=xmed;
  position[1]=ymed;
  position[2]=0.;
  sigma[0]=TMath::Sqrt(sigx*sigx+resolx*resolx);
  sigma[1]=TMath::Sqrt(sigy*sigy+resoly*resoly);
  sigma[2]=sigz;

  AliESDVertex *vertexSPD = new AliESDVertex(position,sigma,"vtxmean");
  vertexSPD->PrintStatus();

  man->Put(vertexSPD,idSPD,metadata);



}



