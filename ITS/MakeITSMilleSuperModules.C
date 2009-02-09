#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TError.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliITSMisalignMaker.h"
#endif
void MakeITSMilleSuperModules(char *geomfile="geometry.root") {
//========================================================================
//
// Macro for building supermodules for ITS alignment
//
// Main author: M. Lunardon
//
//========================================================================

  const char* macroname = "MakeITSMilleSuperModules.C";

  if (!geomfile) { // look for default geometry
    // Activate CDB storage and load geometry from CDB
    AliCDBManager* cdb = AliCDBManager::Instance();
    if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    cdb->SetRun(0);
    
    AliCDBStorage* storage = NULL;
    
    TString compare("kTRUE");
    if(gSystem->Getenv("TOCDB") == compare.Data()){
      TString Storage = gSystem->Getenv("STORAGE");
      if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
	Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
	return;
      }
      storage = cdb->GetStorage(Storage.Data());
      if(!storage){
	Error(macroname,"Unable to open storage %s\n",Storage.Data());
	return;
      }
      AliCDBPath path("GRP","Geometry","Data");
      AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
      if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
      entry->SetOwner(0);
      TGeoManager* geom = (TGeoManager*) entry->GetObject();
      AliGeomManager::SetGeometry(geom);
    }else{
      AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage
    }
  }
  else 
    AliGeomManager::LoadGeometry(geomfile); //load geom from file
  //////////////////////////////////////////////////////////////////
  
  TClonesArray *array = new TClonesArray("AliAlignObjParams",2000);
  TClonesArray &alobj = *array;
  Int_t j = 0;
  
  // custom ITS supermodules are stored as AliAlignObjParams as follow:
  // symname, volid = SMsymname, SMvolid
  // matrix = TGeoHMatrix of the SM in the global c.s.
  // symname: the list of sensitive volumes, starting with "ITSMilleModuleList:" 
  //          keyword, is appended to the symname
  //          example:  "ITSMilleModuleList: N1 N2 N3-N4 ..."
  //                  where N is a ITS sensitive volume INDEX (0-2197)

  Char_t symname[4096];
  UShort_t volid;
  Char_t modlist[4096];
  Double_t t[3]; // global translation
  Double_t r[9]; // global rotation
  TGeoHMatrix m;

  // virtual volids start at layer 7
  volid=14336;

  // LAYERS: volids 14336 to 14341
  strcpy(modlist,"ITSMilleModuleList: 0-79");
  sprintf(symname,"%s  %s","ITS/SPD0",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 80-239");
  sprintf(symname,"%s  %s","ITS/SPD1",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 240-323");
  sprintf(symname,"%s  %s","ITS/SDD2",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 324-499");
  sprintf(symname,"%s  %s","ITS/SDD3",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 500-1247");
  sprintf(symname,"%s  %s","ITS/SSD4",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 1248-2197");
  sprintf(symname,"%s  %s","ITS/SSD5",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  // SPD BARREL: volid 14342
  strcpy(modlist,"ITSMilleModuleList: 0-239");
  sprintf(symname,"%s  %s","ITS/SPD/Barrel",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  
  // SPD HALF BARREL: volids 14343-14344
  strcpy(modlist,"ITSMilleModuleList: 0-39 80-159");
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel0",modlist); // up
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  strcpy(modlist,"ITSMilleModuleList: 40-79 160-239");
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel1",modlist); // down
  m.Clear(); // global frame
  //m.RotateZ(180.);
  //m.RotateY(180.); // just negY->posY
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;

  // SPD sectors: volids 14345 to 14354
  for (int is=0; is<10; is++) {
    // sect 0: 0-7  80-95
    // sect 1: 8-15 96-111
    // ...
    sprintf(modlist,"ITSMilleModuleList: %d-%d %d-%d",is*8,is*8+7,is*16+80,is*16+80+15);
    sprintf(symname,"ITS/SPD0/Sector%d",is);
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SPD sector\n");
      return;
    }
    sprintf(symname,"%s%d  %s","ITS/SPD0/Sector",is,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }

  // SPD staves: volids 14355 to 14414
  for (int is=0; is<60; is++) {
    // Stave0: 0-3
//     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d",is*4,is*4+3);
    strcpy(symname,AliITSgeomTGeo::GetSymName(is*4));
    char *clad=strstr(symname,"Ladder");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Ladder' in symname\n");
      return;
    }    
    //printf("Symname=%s\n",symname);
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SPD stave\n");
      return;
    }

    clad=strstr(symname,"HalfStave");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'HalfStave' in symname\n");
      return;
    }    
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }


  // SPD HALF staves: volids 14415 to 14534
  for (int is=0; is<120; is++) {
    // HalfStave0: 0-1
//     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d ",is*2,is*2+1);
    strcpy(symname,AliITSgeomTGeo::GetSymName(is*2));
    char *clad=strstr(symname,"Ladder");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Ladder' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SPD half-stave\n");
      return;
    }
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }

  // SPD HALF-HALF BARREL UP/DOWN FORWARD/BACKWARD (Y>0 and Z>0): volids 14535 to 14538
  strcpy(modlist,"ITSMilleModuleList: ");
  for (int ii=0; ii<40; ii++) {
    int ij=ii/2;
    if (!(ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  for (int ii=80; ii<160; ii++) {
    int ij=ii/2;
    if (!(ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel0fw",modlist); // up/fw
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  //-----
  strcpy(modlist,"ITSMilleModuleList: ");
  for (int ii=0; ii<40; ii++) {
    int ij=ii/2;
    if ((ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  for (int ii=80; ii<160; ii++) {
    int ij=ii/2;
    if ((ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel0bw",modlist); // up/bw
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  //-----
  strcpy(modlist,"ITSMilleModuleList: ");
  for (int ii=40; ii<80; ii++) {
    int ij=ii/2;
    if (!(ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  for (int ii=160; ii<240; ii++) {
    int ij=ii/2;
    if (!(ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel1fw",modlist); // down/fw
  m.Clear(); // global frame
  //m.RotateZ(180.);
  //m.RotateY(180.); // just negY->posY
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  //-----
  strcpy(modlist,"ITSMilleModuleList: ");
  for (int ii=40; ii<80; ii++) {
    int ij=ii/2;
    if ((ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  for (int ii=160; ii<240; ii++) {
    int ij=ii/2;
    if ((ij%2)) sprintf(modlist,"%s %d",modlist,ii); 
  }
  sprintf(symname,"%s  %s","ITS/SPD/HalfBarrel1bw",modlist); // down/bw
  m.Clear(); // global frame
  //m.RotateZ(180.);
  //m.RotateY(180.); // just negY->posY
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  //-----

  // SDD
  // start at volid=15000
  // - SDD2: 14 ladders da 6 elementi
  // - SDD3: 22 ladders da 8 elementi
  volid=15000;

  // SDD2 Ladders: volids 15000 to 15013
  for (int is=0; is<14; is++) {
    // Ladder: 0-5
    //     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d ",240+is*6,240+is*6+5);
    strcpy(symname,AliITSgeomTGeo::GetSymName(240+is*6));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SDD ladder\n");
      return;
    }
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }

  // SDD3 Ladders: volids 15014 to 15035
  for (int is=0; is<22; is++) {
    // Ladder: 0-7
    //     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d ",324+is*8,324+is*8+7);
    strcpy(symname,AliITSgeomTGeo::GetSymName(324+is*8));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SDD ladder\n");
      return;
    }
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }

  //----------

  // SSD
  // start at volid=16000
  // - SSD4: 34 ladders da 22 elementi
  // - SSD5: 38 ladders da 25 elementi
  volid=16000;

  // SSD4 Ladders: volids 16000 to 16033
  for (int is=0; is<34; is++) {
    // Ladder: 0-5
    //     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d ",500+is*22,500+is*22+21);
    strcpy(symname,AliITSgeomTGeo::GetSymName(500+is*22));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SSD ladder\n");
      return;
    }
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }

  // SSD5 Ladders: volids 16034 to 16071
  for (int is=0; is<38; is++) {
    // Ladder: 0-7
    //     // ...
    sprintf(modlist," ITSMilleModuleList: %d-%d ",1248+is*25,1248+is*25+24);
    strcpy(symname,AliITSgeomTGeo::GetSymName(1248+is*25));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SDD ladder\n");
      return;
    }
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }
  
  // SSD BARREL: volid 16072
  strcpy(modlist,"ITSMilleModuleList: 500-2197");
  sprintf(symname,"%s  %s","ITS/SSD/Barrel",modlist);
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  
  
  // SSD HALF BARREL: volids 16073-16074
  // SSD HALF BARREL UP: volid 16073
  strcpy(modlist,"ITSMilleModuleList: 500-873 1248-1722");
  sprintf(symname,"%s  %s","ITS/SSD/HalfBarrel0",modlist); // up
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  
  // SSD HALF BARREL DOWN: volid 16074
  strcpy(modlist,"ITSMilleModuleList: 874-1246 1723-2197");
  sprintf(symname,"%s  %s","ITS/SSD/HalfBarrel1",modlist); // down
  m.Clear(); // global frame
  new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
  *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
  j++; volid++;
  
  
  //----------
  // SSD HALF-LADDERS
  // - SSD4: 34x2=68 halfladders da 12+10 elementi
  // - SSD5: 38x2=76 halfladders da 12+13 elementi

  // SSD4 Ladders: volids 16075 to 16142
  for (int is=0; is<34; is++) {
    // Ladder: 0-33
    // first half
    sprintf(modlist," ITSMilleModuleList: %d-%d ",500+is*22,500+is*22+11);
    strcpy(symname,AliITSgeomTGeo::GetSymName(500+is*22));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SSD ladder\n");
      return;
    }
    strcat(symname,"Half0");
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;

    // second half
    sprintf(modlist," ITSMilleModuleList: %d-%d ",500+is*22+12,500+is*22+21);
    strcpy(symname,AliITSgeomTGeo::GetSymName(500+is*22));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SSD ladder\n");
      return;
    }
    strcat(symname,"Half1");
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;

  }

  // SSD5 Ladders: volids 16143 to 16218
  for (int is=0; is<38; is++) {
    //  first half
    sprintf(modlist," ITSMilleModuleList: %d-%d ",1248+is*25,1248+is*25+11);
    strcpy(symname,AliITSgeomTGeo::GetSymName(1248+is*25));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SDD ladder\n");
      return;
    }
    strcat(symname,"Half0");
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;

    // second half
    sprintf(modlist," ITSMilleModuleList: %d-%d ",1248+is*25+12,1248+is*25+24);
    strcpy(symname,AliITSgeomTGeo::GetSymName(1248+is*25));
    char *clad=strstr(symname,"Sensor");
    if (clad) *(clad-1) = NULL;
    else  {
      Error(macroname,"cannot find 'Sensor' in symname\n");
      return;
    }    
    if (AliGeomManager::GetMatrix(symname))
      m=(*AliGeomManager::GetMatrix(symname));
    else {
      Error(macroname,"cannot find matrix for SDD ladder\n");
      return;
    }
    strcat(symname,"Half1");
    strcat(symname,modlist);

    new(alobj[j]) AliAlignObjParams(symname, volid, m, kTRUE);
    *(strstr(symname,"ITSMilleModuleList"))=NULL ;printf("added module %s with volid %d \n",symname,volid);
    j++; volid++;
  }


  //////////////////////////////////////////////////////////////////
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSMilleSuperModules.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving ITS SuperModules as AliAlignObjs to the file %s", filename);
    f.cd();
    f.WriteObject(array,"ITSMilleSuperModules","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving ITS SuperModules in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Marcello Lunardon");
    md->SetComment("ITS super modules for ITS alignment with Millepede");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id, md);
  }
  
  array->Delete();
  return;
}
