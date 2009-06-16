#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include "STEER/AliCDBEntry.h"
#include "STEER/AliCDBGrid.h"
#include "STEER/AliCDBId.h"
#include "STEER/AliCDBLocal.h"
#include "STEER/AliCDBManager.h"
#include "STEER/AliCDBMetaData.h"
#include "STEER/AliCDBPath.h"
#include "STEER/AliCDBRunRange.h"
#include "STEER/AliCDBStorage.h"
#include "ZDC/AliZDC.h"
#include "ZDC/AliZDCv3.h"
#include "ZDC/AliZDCRecoParam.h"
#include "ZDC/AliZDCRecoParampp.h"
#include "ZDC/AliZDCRecoParamPbPb.h"

#endif

void MakeZDCRecoParam(Int_t type=0){
//========================================================================
//
// Steering macro to create and store in OCDB
//       ZDC reconstruction parameters
//
// Contact: chiara.oppedisano@to.infn.it
//
//========================================================================

  AliCDBManager* cdb = AliCDBManager::Instance();
  //if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  if(type==0){
    AliZDCRecoParampp* zdcRecoParam = new AliZDCRecoParampp();
    // save in CDB storage
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Chiara Oppedisano");
    md->SetComment("p-p reconstruction parameters for ZDC");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetObjectClassName("AliZDCRecoParampp");
    AliCDBId id("ZDC/Calib/RecoParampp",0,AliCDBRunRange::Infinity());
    cdb->GetDefaultStorage()->Put(zdcRecoParam,id, md);
  }
  else if(type==1){ 
    //    
    TH1::AddDirectory(0);
    TH2::AddDirectory(0);
    //
    TFile * fileHistos = TFile::Open("$ALICE_ROOT/ZDC/GlauberMCHistos.root");
    fileHistos->cd();
    //
    TH2F *hZDCvsZEM = (TH2F*) fileHistos->Get("hZDCvsZEM");
    hZDCvsZEM->SetDirectory(0);
    //
    TH2F *hZDCCvsZEM = (TH2F*) fileHistos->Get("hZDCCvsZEM");
    hZDCCvsZEM->SetDirectory(0);
    //
    TH2F *hZDCAvsZEM = (TH2F*) fileHistos->Get("hZDCAvsZEM");
    hZDCAvsZEM->SetDirectory(0);
    //
    TH1D* hDist = (TH1D*) fileHistos->Get("hDist");
    hDist->SetDirectory(0);
    //
    TH1D* hbDist = (TH1D*) fileHistos->Get("hbDist");
    hbDist->SetDirectory(0);
    
    AliZDCRecoParamPbPb* zdcRecoParam = 
    	new AliZDCRecoParamPbPb(hZDCvsZEM, hZDCCvsZEM, 
		hZDCAvsZEM, hDist, hbDist, 0.1);
    //zdcRecoParam->Print("");
        
    // save in CDB storage
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Chiara Oppedisano");
    md->SetComment("A-A reconstruction parameters for ZDC");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetObjectClassName("AliZDCRecoParamPbPb");
    AliCDBId id("ZDC/Calib/RecoParamPbPb", 0, AliCDBRunRange::Infinity());
    cdb->GetDefaultStorage()->Put(zdcRecoParam, id, md);
    //
    fileHistos->Close();
  }
  else{
    printf("Event type not implemented\n");
    return;
  }

}
