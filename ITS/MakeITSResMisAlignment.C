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


TRandom3 rnd;
rnd.SetSeed(98723456);

void MakeITSResMisAlignment() {
//========================================================================
//
// Steering macro for ITS residual (realistic) misalignment
//
// Main author: L. Gaudichet
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================


  const char* macroname = "MakeITSResMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliCDBStorage* storage = NULL;

  if(TString(gSystem->Getenv("TOCDB")) == TString("kTRUE")){
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

  //   SETTINGS:
  //   - tranformations are defined by the "maximum transformation" 
  //     dx,dy,dz,dpsi,dtheta,dphi: then we take a Gaussian with sigma=dx/3
  //     and we cut it at +- 3*sigma = dx 
  //     (the option "unif" allows to sample from a uniform distr.)
  //   - units are cm and deg
  //   - transformations are defined in the local frame of the volume
  //     being misaligned
  //
  const Float_t kRadToDeg = 180./TMath::Pi();


  //=****************************************
  // misalignment of the whole ITS according to survey as reported by Werner Riegler (18/07/2008)
  // no smearing added (would clash with vertex constraint)
  //=****************************************
  Double_t its_dx     = -0.12;
  Double_t its_dy     = -0.07;
  Double_t its_dz     = 0.29;
  Double_t its_dpsi   = 0.;
  Double_t its_dtheta = 0.03;
  Double_t its_dphi   = 0.04;
  const char* unifits="fixed";

  //=****************************************
  // misalignment at the level of SPD sectors : source - A.Pepato
  //=****************************************
  Float_t spdsector_dx     = 0.0050/5.; //  50 micron (~tangetial, i.e. rphi) 
  Float_t spdsector_dy     = 0.0100/5.; // 100 micron (~radial)
  Float_t spdsector_dz     = 0.0100/5.; // 100 micron
  Float_t spdsector_dpsi   = 0.0100/30.*kRadToDeg/5.; // so as to have 100 micron difference at the two extremes
  Float_t spdsector_dtheta = 0.0100/30.*kRadToDeg/5.; // so as to have 100 micron difference at the two extremes
  Float_t spdsector_dphi   = 0.0050/1.5*kRadToDeg/5.; // so as to have 50 micron difference at the two extremes
  Bool_t unifspdsector=kFALSE;

  //=****************************************
  // misalignment at the level of SPD half-barrels : source - A.Pepato
  //=****************************************
  Float_t spdhalfbarrel_dx     = 0.000; // 200 micron  
  Float_t spdhalfbarrel_dy     = 0.000; // 200 micron 
  Float_t spdhalfbarrel_dz     = 0.000; // 200 micron
  Float_t spdhalfbarrel_dpsi   = 0.000/30.*kRadToDeg; // so as to have 100 micron difference at the two extremes
  Float_t spdhalfbarrel_dtheta = 0.000/30.*kRadToDeg; // so as to have 100 micron difference at the two extremes
  Float_t spdhalfbarrel_dphi   = 0.000/7.*kRadToDeg; // so as to have 100 micron difference at the two extremes

  //=****************************************
  // misalignment at the level of SPD barrel : source - A.Pepato
  //=****************************************
  Float_t spdbarrel_dx     = 0.000; // 1 mm (very pessimistic)  
  Float_t spdbarrel_dy     = 0.000; // 1 mm (very pessimistic)
  Float_t spdbarrel_dz     = 0.000; // 1 mm (very pessimistic)
  Float_t spdbarrel_dpsi   = 0.000/30.*kRadToDeg; // so as to have 500 micron difference at the two extremes
  Float_t spdbarrel_dtheta = 0.000/30.*kRadToDeg; // so as to have 500 micron difference at the two extremes
  Float_t spdbarrel_dphi   = 0.000/7.*kRadToDeg; // so as to have 500 micron difference at the two extremes
  

  //=****************************************
  // misalignment at the level of SDD and SSD layers: source
  //=****************************************
  Float_t sddlayer_dx     = 0.0000; 
  Float_t sddlayer_dy     = 0.0000; 
  Float_t sddlayer_dz     = 0.0000; 
  Float_t sddlayer_dpsi   = 0.0000;
  Float_t sddlayer_dtheta = 0.0000; 
  Float_t sddlayer_dphi   = 0.0000;

  Float_t ssdlayer_dx     = 0.0000;
  Float_t ssdlayer_dy     = 0.0000;
  Float_t ssdlayer_dz     = 0.0000;
  Float_t ssdlayer_dpsi   = 0.0000;
  Float_t ssdlayer_dtheta = 0.0000;
  Float_t ssdlayer_dphi   = 0.0000;
  

  //=****************************************
  // misalignment at the level of half-staves (SPD) : source - S.Moretto
  //                              ladders (SDD,SSD) : source -
  //=****************************************
  Float_t spdhalfstave_dx     = 0.0100/4.; // 100 micron  
  Float_t spdhalfstave_dy     = 0.0020/4.; // 20 micron 
  Float_t spdhalfstave_dz     = 0.0020/4.; // 20 micron
  Float_t spdhalfstave_dpsi   = 0.0020/7.*kRadToDeg/4.; // so as to have 20 micron difference at the two extremes
  Float_t spdhalfstave_dtheta = 0.0050/7.*kRadToDeg/4.; // so as to have 50 micron difference at the two extremes
  Float_t spdhalfstave_dphi   = 0.0050/0.7*kRadToDeg/4.; // so as to have 50 micron difference at the two extremes
  Bool_t unifspdhalfstave=kFALSE;

  Float_t sddladder_dx     = 0.0005; // 5 micron  
  Float_t sddladder_dy     = 0.0005; // 5 micron 
  Float_t sddladder_dz     = 0.0005; // 5 micron
  Float_t sddladder_dpsi   = 0.00; //  ?
  Float_t sddladder_dtheta = 0.00; //  ?
  Float_t sddladder_dphi   = 0.00; //  ?

  Float_t ssdladder_dx     = 0.0005; // 5 micron  
  Float_t ssdladder_dy     = 0.0005; // 5 micron 
  Float_t ssdladder_dz     = 0.0005; // 5 micron
  Float_t ssdladder_dpsi   = 0.00; //  ?
  Float_t ssdladder_dtheta = 0.00; //  ?
  Float_t ssdladder_dphi   = 0.00; //  ?


  //=****************************************
  // misalignment at the level of ladders (SPD) : source - R.Santoro
  //                              modules (SDD) : source - L.Gaudichet
  //                              modules (SSD) : source - 
  //=****************************************
  Float_t spdladder_dx     = 0.0010/5.; // 10 micron  
  Float_t spdladder_dy     = 0.0050/5.; // 50 micron 
  Float_t spdladder_dz     = 0.0010/5.; // 10 micron
  Float_t spdladder_dpsi   = 0.0001*kRadToDeg/5.; // 0.1 mrad
  Float_t spdladder_dtheta = 0.0001*kRadToDeg/5.; // 0.1 mrad
  Float_t spdladder_dphi   = 0.0001*kRadToDeg/5.; // 0.1 mrad

  Float_t sddmodule_dx     = 0.0045/5.; // 45 micron  
  Float_t sddmodule_dy     = 0.0045/5.; // 45 micron 
  Float_t sddmodule_dz     = 0.0105/5.; // 105 micron
  Float_t sddmodule_dpsi   = 0.00; // ?
  Float_t sddmodule_dtheta = 0.00; //  ?
  Float_t sddmodule_dphi   = 0.00; //  ?

  Float_t ssdmodule_dx     = 0.0050/5.; // 50 micron  
  Float_t ssdmodule_dy     = 0.0050/5.; // 50 micron 
  Float_t ssdmodule_dz     = 0.0050/5.; // 50 micron
  Float_t ssdmodule_dpsi   = 0.00; // ?
  Float_t ssdmodule_dtheta = 0.00; //  ?
  Float_t ssdmodule_dphi   = 0.00; //  ?
  //
  // END SETTINGS

  AliITSMisalignMaker alignMaker;

  //=****************************************
  // overall ITS misalignment :
  //=****************************************

  alignMaker.AddAlignObj("ITS",its_dx,its_dy,its_dz,its_dpsi,its_dtheta,its_dphi,unifits);


  //=****************************************
  // misalignment at the level of SPD barrel, half-barrels, and at the level
  // of SPD sectors
  //=****************************************

  Double_t vx,vy,vz,vpsi,vtheta,vphi;
  Double_t vxbarrel,vybarrel,vzbarrel,vpsibarrel,vthetabarrel,vphibarrel;

  //   barrel
  vxbarrel = AliMathBase::TruncatedGaus(0.,spdbarrel_dx/3,spdbarrel_dx);
  vxbarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dx/3,spdbarrel_dx);
  vybarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dy/3,spdbarrel_dy);
  vzbarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dz/3,spdbarrel_dz);
  vpsibarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dpsi/3,spdbarrel_dpsi);
  vthetabarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dtheta/3,spdbarrel_dtheta);
  vphibarrel = AliMathBase::TruncatedGaus(0,spdbarrel_dphi/3,spdbarrel_dphi);

  //  top half-barrel
  vx = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dx/3,spdhalfbarrel_dx);
  vy = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dy/3,spdhalfbarrel_dy);
  vz = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dz/3,spdhalfbarrel_dz);
  vpsi = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dpsi/3,spdhalfbarrel_dpsi);
  vtheta = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dtheta/3,spdhalfbarrel_dtheta);
  vphi = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dphi/3,spdhalfbarrel_dphi);

  vx += vxbarrel;
  vy += vybarrel;
  vz += vzbarrel;
  vpsi += vpsibarrel;
  vtheta += vthetabarrel;
  vphi += vphibarrel;

  alignMaker.AddSectorAlignObj(1,5,spdsector_dx,spdsector_dy,spdsector_dz,
			       spdsector_dpsi,spdsector_dtheta,spdsector_dphi,
			       vx,vy,vz,vpsi,vtheta,vphi,unifspdsector);

  //  bottom half-barrel
  vx = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dx/3,spdhalfbarrel_dx);
  vy = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dy/3,spdhalfbarrel_dy);
  vz = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dz/3,spdhalfbarrel_dz);
  vpsi = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dpsi/3,spdhalfbarrel_dpsi);
  vtheta = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dtheta/3,spdhalfbarrel_dtheta);
  vphi = AliMathBase::TruncatedGaus(0,spdhalfbarrel_dphi/3,spdhalfbarrel_dphi);

  vx += vxbarrel;
  vy += vybarrel;
  vz += vzbarrel;
  vpsi += vpsibarrel;
  vtheta += vthetabarrel;
  vphi += vphibarrel;

  alignMaker.AddSectorAlignObj(6,10,spdsector_dx,spdsector_dy,spdsector_dz,
			       spdsector_dpsi,spdsector_dtheta,spdsector_dphi,
			       vx,vy,vz,vpsi,vtheta,vphi,unifspdsector);


  //=****************************************
  // misalignment at the level of half-staves (SPD)/ladders (SDD,SSD) :
  //=****************************************

  alignMaker.AddAlignObj(0,-1,spdhalfstave_dx,spdhalfstave_dy,spdhalfstave_dz,spdhalfstave_dpsi,spdhalfstave_dtheta,spdhalfstave_dphi,0,0,0,0,0,0,unifspdhalfstave); // all SPD1 half-staves
  alignMaker.AddAlignObj(1,-1,spdhalfstave_dx,spdhalfstave_dy,spdhalfstave_dz,spdhalfstave_dpsi,spdhalfstave_dtheta,spdhalfstave_dphi,0,0,0,0,0,0,unifspdhalfstave); // all SPD2 half-staves

  alignMaker.AddAlignObj(2,-1,sddladder_dx,sddladder_dy,sddladder_dz,sddladder_dpsi,sddladder_dtheta,sddladder_dphi,0,0,0,0,0,0,kFALSE); // all SDD1 ladders
  alignMaker.AddAlignObj(3,-1,sddladder_dx,sddladder_dy,sddladder_dz,sddladder_dpsi,sddladder_dtheta,sddladder_dphi,0,0,0,0,0,0,kFALSE); // all SDD2 ladders

  alignMaker.AddAlignObj(4,-1,ssdladder_dx,ssdladder_dy,ssdladder_dz,ssdladder_dpsi,ssdladder_dtheta,ssdladder_dphi,0,0,0,0,0,0,kFALSE); // all SSD1 ladders
  alignMaker.AddAlignObj(5,-1,ssdladder_dx,ssdladder_dy,ssdladder_dz,ssdladder_dpsi,ssdladder_dtheta,ssdladder_dphi,0,0,0,0,0,0,kFALSE); // all SSD2 ladders

  //=****************************************
  // misalignment at the level of ladders (SPD)/modules (SDD,SSD) :
  //=****************************************

  alignMaker.AddAlignObj(0,spdladder_dx,spdladder_dy,spdladder_dz,spdladder_dpsi,spdladder_dtheta,spdladder_dphi,kFALSE); // all SPD1 ladders
  alignMaker.AddAlignObj(1,spdladder_dx,spdladder_dy,spdladder_dz,spdladder_dpsi,spdladder_dtheta,spdladder_dphi,kFALSE); // all SPD2 ladders

  alignMaker.AddAlignObj(2,sddmodule_dx,sddmodule_dy,sddmodule_dz,sddmodule_dpsi,sddmodule_dtheta,sddmodule_dphi,kFALSE); // all SDD1 modules
  alignMaker.AddAlignObj(3,sddmodule_dx,sddmodule_dy,sddmodule_dz,sddmodule_dpsi,sddmodule_dtheta,sddmodule_dphi,kFALSE); // all SDD2 modules

  alignMaker.AddAlignObj(4,ssdmodule_dx,ssdmodule_dy,ssdmodule_dz,ssdmodule_dpsi,ssdmodule_dtheta,ssdmodule_dphi,kFALSE); // all SSD1 modules
  alignMaker.AddAlignObj(5,ssdmodule_dx,ssdmodule_dy,ssdmodule_dz,ssdmodule_dpsi,ssdmodule_dtheta,ssdmodule_dphi,kFALSE); // all SSD2 modules


  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f.IsOpen()){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s",filename);
    f.cd();
    f.WriteObject(alignMaker.GetArray(),"ITSAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Andrea Dainese");
    md->SetComment("Alignment objects with residual ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(alignMaker.GetArray(),id,md);
  }


  return;
}

