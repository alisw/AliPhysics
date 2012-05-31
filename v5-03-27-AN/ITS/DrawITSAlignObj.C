#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1F.h>
#include <TString.h>
#include <Riostream.h>
#include "AliAlignObj.h"
#endif
void DrawITSAlignObj(Bool_t local=kFALSE) {  //
  // Draw distribution of alignment constants
  // Set TOCDB and STORAGE environment variables to run the macro
  // on an AliCDBEntry instead of on a file
  //

  const char* filenameOut="ITSfullModuleMisalignment.root";
  TClonesArray* ar = 0;
  
  AliCDBManager* cdb = AliCDBManager::Instance();

  // Activate CDB storage and load geometry from CDB
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    cdb->SetDefaultStorage(Storage.Data());
    cdb->SetRun(0);
    AliCDBEntry *entry = cdb->Get("GRP/Geometry/Data");
    if(!entry) Fatal("Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
    AliCDBEntry* eItsAlign = cdb->Get("ITS/Align/Data");
    ar = (TClonesArray*)eItsAlign->GetObject();
    if(!ar) {
      Fatal("Could not get the alignment-objects array from the CDB entry!");
      return;
    }
  }else{
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    cdb->SetRun(0);
    AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage
    const char* filename="ITSfullMisalignment.root";
    TFile* f = TFile::Open(filename);
    if(!f) {cerr<<"cannot open input file\n"; return;}
    ar = (TClonesArray*)f->Get("ITSAlignObjs");
    if(!ar) {
      Fatal("Could not get the alignment-objects array from the file %s!", filename);
      return;
    }
  }    
		  
  TCanvas  *cSPDsector = new TCanvas("cSPDsector","SPD sectors alignobjs",0,0,800,600);
  cSPDsector->Divide(3,2);		    
  TH1F* hxSPDsector = new TH1F("hxSPDsector","x in SPD sectors",100,-0.1,0.1);
  hxSPDsector->SetXTitle("[cm]");
  TH1F* hySPDsector = new TH1F("hySPDsector","y in SPD sectors",100,-0.1,0.1);
  hySPDsector->SetXTitle("[cm]");
  TH1F* hzSPDsector = new TH1F("hzSPDsector","z in SPD sectors",100,-0.1,0.1);
  hzSPDsector->SetXTitle("[cm]");
  TH1F* hrxSPDsector = new TH1F("hrxSPDsector","rx in SPD sectors",100,-0.5,0.5);
  hrxSPDsector->SetXTitle("[deg]");
  TH1F* hrySPDsector = new TH1F("hrySPDsector","ry in SPD sectors",100,-0.5,0.5);
  hrySPDsector->SetXTitle("[deg]");
  TH1F* hrzSPDsector = new TH1F("hrzSPDsector","rz in SPD sectors",100,-0.5,0.5);
  hrzSPDsector->SetXTitle("[deg]");

  TCanvas  *cSPDhalfstave = new TCanvas("cSPDhalfstave","SPD halfstaves alignobjs",0,0,800,600);
  cSPDhalfstave->Divide(3,2);		    
  TH1F* hxSPDhalfstave = new TH1F("hxSPDhalfstave","x in SPD halfstaves",100,-0.01,0.01);
  hxSPDhalfstave->SetXTitle("[cm]");
  TH1F* hySPDhalfstave = new TH1F("hySPDhalfstave","y in SPD halfstaves",100,-0.01,0.01);
  hySPDhalfstave->SetXTitle("[cm]");
  TH1F* hzSPDhalfstave = new TH1F("hzSPDhalfstave","z in SPD halfstaves",100,-0.01,0.01);
  hzSPDhalfstave->SetXTitle("[cm]");
  TH1F* hrxSPDhalfstave = new TH1F("hrxSPDhalfstave","rx in SPD halfstaves",100,-0.1,0.1);
  hrxSPDhalfstave->SetXTitle("[deg]");
  TH1F* hrySPDhalfstave = new TH1F("hrySPDhalfstave","ry in SPD halfstaves",100,-0.1,0.1);
  hrySPDhalfstave->SetXTitle("[deg]");
  TH1F* hrzSPDhalfstave = new TH1F("hrzSPDhalfstave","rz in SPD halfstaves",100,-0.5,0.5);
  hrzSPDhalfstave->SetXTitle("[deg]");

  TCanvas  *cSPDladder = new TCanvas("cSPDladder","SPD ladders alignobjs",0,0,800,600);
  cSPDladder->Divide(3,2);		    
  TH1F* hxSPDladder = new TH1F("hxSPDladder","x in SPD ladders",100,-0.01,0.01);
  hxSPDladder->SetXTitle("[cm]");
  TH1F* hySPDladder = new TH1F("hySPDladder","y in SPD ladders",100,-0.01,0.01);
  hySPDladder->SetXTitle("[cm]");
  TH1F* hzSPDladder = new TH1F("hzSPDladder","z in SPD ladders",100,-0.01,0.01);
  hzSPDladder->SetXTitle("[cm]");
  TH1F* hrxSPDladder = new TH1F("hrxSPDladder","rx in SPD ladders",100,-0.01,0.01);
  hrxSPDladder->SetXTitle("[deg]");
  TH1F* hrySPDladder = new TH1F("hrySPDladder","ry in SPD ladders",100,-0.01,0.01);
  hrySPDladder->SetXTitle("[deg]");
  TH1F* hrzSPDladder = new TH1F("hrzSPDladder","rz in SPD ladders",100,-0.01,0.01);
  hrzSPDladder->SetXTitle("[deg]");

  TCanvas  *cSDDladder = new TCanvas("cSDDladder","SDD ladders alignobjs",0,0,800,600);
  cSDDladder->Divide(3,2);		    
  TH1F* hxSDDladder = new TH1F("hxSDDladder","x in SDD ladders",100,-0.01,0.01);
  hxSDDladder->SetXTitle("[cm]");
  TH1F* hySDDladder = new TH1F("hySDDladder","y in SDD ladders",100,-0.01,0.01);
  hySDDladder->SetXTitle("[cm]");
  TH1F* hzSDDladder = new TH1F("hzSDDladder","z in SDD ladders",100,-0.01,0.01);
  hzSDDladder->SetXTitle("[cm]");
  TH1F* hrxSDDladder = new TH1F("hrxSDDladder","rx in SDD ladders",100,-0.01,0.01);
  hrxSDDladder->SetXTitle("[deg]");
  TH1F* hrySDDladder = new TH1F("hrySDDladder","ry in SDD ladders",100,-0.01,0.01);
  hrySDDladder->SetXTitle("[deg]");
  TH1F* hrzSDDladder = new TH1F("hrzSDDladder","rz in SDD ladders",100,-0.01,0.01);
  hrzSDDladder->SetXTitle("[deg]");

  TCanvas  *cSDDsensor = new TCanvas("cSDDsensor","SDD sensors alignobjs",0,0,800,600);
  cSDDsensor->Divide(3,2);		    
  TH1F* hxSDDsensor = new TH1F("hxSDDsensor","x in SDD sensors",100,-0.01,0.01);
  hxSDDsensor->SetXTitle("[cm]");
  TH1F* hySDDsensor = new TH1F("hySDDsensor","y in SDD sensors",100,-0.01,0.01);
  hySDDsensor->SetXTitle("[cm]");
  TH1F* hzSDDsensor = new TH1F("hzSDDsensor","z in SDD sensors",100,-0.01,0.01);
  hzSDDsensor->SetXTitle("[cm]");
  TH1F* hrxSDDsensor = new TH1F("hrxSDDsensor","rx in SDD sensors",100,-0.01,0.01);
  hrxSDDsensor->SetXTitle("[deg]");
  TH1F* hrySDDsensor = new TH1F("hrySDDsensor","ry in SDD sensors",100,-0.01,0.01);
  hrySDDsensor->SetXTitle("[deg]");
  TH1F* hrzSDDsensor = new TH1F("hrzSDDsensor","rz in SDD sensors",100,-0.01,0.01);
  hrzSDDsensor->SetXTitle("[deg]");


  TCanvas  *cSSDladder = new TCanvas("cSSDladder","SSD ladders alignobjs",0,0,800,600);
  cSSDladder->Divide(3,2);		    
  TH1F* hxSSDladder = new TH1F("hxSSDladder","x in SSD ladders",100,-0.01,0.01);
  hxSSDladder->SetXTitle("[cm]");
  TH1F* hySSDladder = new TH1F("hySSDladder","y in SSD ladders",100,-0.01,0.01);
  hySSDladder->SetXTitle("[cm]");
  TH1F* hzSSDladder = new TH1F("hzSSDladder","z in SSD ladders",100,-0.01,0.01);
  hzSSDladder->SetXTitle("[cm]");
  TH1F* hrxSSDladder = new TH1F("hrxSSDladder","rx in SSD ladders",100,-0.01,0.01);
  hrxSSDladder->SetXTitle("[deg]");
  TH1F* hrySSDladder = new TH1F("hrySSDladder","ry in SSD ladders",100,-0.01,0.01);
  hrySSDladder->SetXTitle("[deg]");
  TH1F* hrzSSDladder = new TH1F("hrzSSDladder","rz in SSD ladders",100,-0.01,0.01);
  hrzSSDladder->SetXTitle("[deg]");

  TCanvas  *cSSDsensor = new TCanvas("cSSDsensor","SSD sensors alignobjs",0,0,800,600);
  cSSDsensor->Divide(3,2);		    
  TH1F* hxSSDsensor = new TH1F("hxSSDsensor","x in SSD sensors",100,-0.01,0.01);
  hxSSDsensor->SetXTitle("[cm]");
  TH1F* hySSDsensor = new TH1F("hySSDsensor","y in SSD sensors",100,-0.01,0.01);
  hySSDsensor->SetXTitle("[cm]");
  TH1F* hzSSDsensor = new TH1F("hzSSDsensor","z in SSD sensors",100,-0.01,0.01);
  hzSSDsensor->SetXTitle("[cm]");
  TH1F* hrxSSDsensor = new TH1F("hrxSSDsensor","rx in SSD sensors",100,-0.01,0.01);
  hrxSSDsensor->SetXTitle("[deg]");
  TH1F* hrySSDsensor = new TH1F("hrySSDsensor","ry in SSD sensors",100,-0.01,0.01);
  hrySSDsensor->SetXTitle("[deg]");
  TH1F* hrzSSDsensor = new TH1F("hrzSSDsensor","rz in SSD sensors",100,-0.01,0.01);
  hrzSSDsensor->SetXTitle("[deg]");


  Int_t entries = ar->GetEntriesFast();
  Printf("number of alignment objects: %d",entries);
  //Bool_t overlaps;
  if(!AliGeomManager::GetGeometry()) return;
  AliGeomManager::ApplyAlignObjsToGeom(*ar); //,overlaps);

  AliAlignObj* a;
  Option_t *opt = NULL;
  Double_t tr[3];
  Double_t rot[3];
  TGeoHMatrix delta;


  TClonesArray *arrayOut = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobjOut = *arrayOut;
  Int_t j=0;

  for(Int_t i=0; i<entries; i++){
    a=(AliAlignObj*)ar->UncheckedAt(i);
    TString symName = a->GetSymName();
    UShort_t volUID = a->GetVolUID();
    //printf("VolId %d    %s\n",volUID,symName.Data());
    
    if(local) {
      a->GetLocalPars(tr,rot);
    } else {
      a->GetPars(tr,rot);
    }

    AliGeomManager::GetDeltaForBranch(*a,delta);
    //delta.Print();

    // create AliAlignObjParam with total delta
    if(i==0 || volUID!=0) {
      new(alobjOut[j]) AliAlignObjParams(symName.Data(),volUID,delta,kTRUE);
      j++;
    }


    // plots
    //
    if(!symName.Contains("SPD") && !symName.Contains("SDD") && !symName.Contains("SSD")) {
      a->Print(opt);
    }
    if(symName.Contains("SPD") && symName.Contains("Sector") && !symName.Contains("HalfStave")) {
      hxSPDsector->Fill(tr[0]);
      hySPDsector->Fill(tr[1]);
      hzSPDsector->Fill(tr[2]);
      hrxSPDsector->Fill(rot[0]);
      hrySPDsector->Fill(rot[1]);
      hrzSPDsector->Fill(rot[2]);
    }
    if(symName.Contains("SPD") && symName.Contains("Sector") && symName.Contains("HalfStave") && !symName.Contains("Ladder")) {
      hxSPDhalfstave->Fill(tr[0]);
      hySPDhalfstave->Fill(tr[1]);
      hzSPDhalfstave->Fill(tr[2]);
      hrxSPDhalfstave->Fill(rot[0]);
      hrySPDhalfstave->Fill(rot[1]);
      hrzSPDhalfstave->Fill(rot[2]);
    }
    if(symName.Contains("SPD") && symName.Contains("Sector") && symName.Contains("HalfStave") && symName.Contains("Ladder")) {
      hxSPDladder->Fill(tr[0]);
      hySPDladder->Fill(tr[1]);
      hzSPDladder->Fill(tr[2]);
      hrxSPDladder->Fill(rot[0]);
      hrySPDladder->Fill(rot[1]);
      hrzSPDladder->Fill(rot[2]);
    }
    if(symName.Contains("SDD") && symName.Contains("Ladder") && !symName.Contains("Sensor")) {
      hxSDDladder->Fill(tr[0]);
      hySDDladder->Fill(tr[1]);
      hzSDDladder->Fill(tr[2]);
      hrxSDDladder->Fill(rot[0]);
      hrySDDladder->Fill(rot[1]);
      hrzSDDladder->Fill(rot[2]);
    }
    if(symName.Contains("SDD") && symName.Contains("Ladder") && symName.Contains("Sensor")) {
      hxSDDsensor->Fill(tr[0]);
      hySDDsensor->Fill(tr[1]);
      hzSDDsensor->Fill(tr[2]);
      hrxSDDsensor->Fill(rot[0]);
      hrySDDsensor->Fill(rot[1]);
      hrzSDDsensor->Fill(rot[2]);
    }
    if(symName.Contains("SSD") && symName.Contains("Ladder") && !symName.Contains("Sensor")) {
      hxSSDladder->Fill(tr[0]);
      hySSDladder->Fill(tr[1]);
      hzSSDladder->Fill(tr[2]);
      hrxSSDladder->Fill(rot[0]);
      hrySSDladder->Fill(rot[1]);
      hrzSSDladder->Fill(rot[2]);
    }
    if(symName.Contains("SSD") && symName.Contains("Ladder") && symName.Contains("Sensor")) {
      hxSSDsensor->Fill(tr[0]);
      hySSDsensor->Fill(tr[1]);
      hzSSDsensor->Fill(tr[2]);
      hrxSSDsensor->Fill(rot[0]);
      hrySSDsensor->Fill(rot[1]);
      hrzSSDsensor->Fill(rot[2]);
    }


  }


  cSPDsector->cd(1);
  hxSPDsector->Draw();
  cSPDsector->cd(2);
  hySPDsector->Draw();
  cSPDsector->cd(3);
  hzSPDsector->Draw();
  cSPDsector->cd(4);
  hrxSPDsector->Draw();
  cSPDsector->cd(5);
  hrySPDsector->Draw();
  cSPDsector->cd(6);
  hrzSPDsector->Draw();

  cSPDhalfstave->cd(1);
  hxSPDhalfstave->Draw();
  cSPDhalfstave->cd(2);
  hySPDhalfstave->Draw();
  cSPDhalfstave->cd(3);
  hzSPDhalfstave->Draw();
  cSPDhalfstave->cd(4);
  hrxSPDhalfstave->Draw();
  cSPDhalfstave->cd(5);
  hrySPDhalfstave->Draw();
  cSPDhalfstave->cd(6);
  hrzSPDhalfstave->Draw();

  cSPDladder->cd(1);
  hxSPDladder->Draw();
  cSPDladder->cd(2);
  hySPDladder->Draw();
  cSPDladder->cd(3);
  hzSPDladder->Draw();
  cSPDladder->cd(4);
  hrxSPDladder->Draw();
  cSPDladder->cd(5);
  hrySPDladder->Draw();
  cSPDladder->cd(6);
  hrzSPDladder->Draw();


  cSDDladder->cd(1);
  hxSDDladder->Draw();
  cSDDladder->cd(2);
  hySDDladder->Draw();
  cSDDladder->cd(3);
  hzSDDladder->Draw();
  cSDDladder->cd(4);
  hrxSDDladder->Draw();
  cSDDladder->cd(5);
  hrySDDladder->Draw();
  cSDDladder->cd(6);
  hrzSDDladder->Draw();
 
  cSDDsensor->cd(1);
  hxSDDsensor->Draw();
  cSDDsensor->cd(2);
  hySDDsensor->Draw();
  cSDDsensor->cd(3);
  hzSDDsensor->Draw();
  cSDDsensor->cd(4);
  hrxSDDsensor->Draw();
  cSDDsensor->cd(5);
  hrySDDsensor->Draw();
  cSDDsensor->cd(6);
  hrzSDDsensor->Draw();

  cSSDladder->cd(1);
  hxSSDladder->Draw();
  cSSDladder->cd(2);
  hySSDladder->Draw();
  cSSDladder->cd(3);
  hzSSDladder->Draw();
  cSSDladder->cd(4);
  hrxSSDladder->Draw();
  cSSDladder->cd(5);
  hrySSDladder->Draw();
  cSSDladder->cd(6);
  hrzSSDladder->Draw();
 
  cSSDsensor->cd(1);
  hxSSDsensor->Draw();
  cSSDsensor->cd(2);
  hySSDsensor->Draw();
  cSSDsensor->cd(3);
  hzSSDsensor->Draw();
  cSSDsensor->cd(4);
  hrxSSDsensor->Draw();
  cSSDsensor->cd(5);
  hrySSDsensor->Draw();
  cSSDsensor->cd(6);
  hrzSSDsensor->Draw();


  TFile fout(filenameOut,"RECREATE");
  fout.cd();
  fout.WriteObject(arrayOut,"ITSAlignObjs","kSingleKey");
  fout.Close();
  
  return;
}

