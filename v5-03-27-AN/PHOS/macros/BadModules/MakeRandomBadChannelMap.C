#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "AliPHOSEmcBadChannelsMap.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#endif

void MakeRandomBadChannelMap(const TString modules="3", 
			     const Float_t fractionBad=0.02)
{
  // Create a random bad channel map for selected modules 
  // with a fraction of bad channel equal to fractionBad
  // Author Yuri Kharlov
  // 20.07.2008

  const Int_t nZ=56, nX=64;
  Int_t module;
  AliPHOSEmcBadChannelsMap badMap;
  TH2I *hBadChannels[5];
  TCanvas *cb[5];
  gStyle->SetOptStat(10);

  for (Int_t iModule=0; iModule<modules.Length(); iModule++) {
    const char chMod = modules[iModule];
    module = atoi(&chMod);
    if (module<1 || module>5) {
      printf("Wrong module number %d, can be from 1 to 5\n",module);
      exit(1);
    }
    printf("Creating bad channel map for module %d\n",module);

    hBadChannels[iModule] = new TH2I(Form("hBadChannels%d",module),
				     Form("Bad channels for module %d",module),
				     nZ,1,nZ, nX,1,nX);
    hBadChannels[iModule]->SetXTitle("z, cells");
    hBadChannels[iModule]->SetYTitle("x, cells");
    for(Int_t iZ=1; iZ<=nZ; iZ++){
      for(Int_t iX=1; iX<=nX; iX++){
	if (gRandom->Rndm() < fractionBad) {
 	  hBadChannels[iModule]->SetBinContent(iZ,iX,1);
 	  badMap.SetBadChannel(module,iZ,iX);
	}
      }
    }
    cb[iModule] = new TCanvas(Form("BadMap%d",module),
			      Form("Bad channels for module %d",module),
			      15+iModule*40,15+iModule*40,500,600) ;
    cb[iModule]->SetFillColor(kWhite);
    hBadChannels[iModule]->DrawClone("box") ;
  }
  

//   TFile * fout = new TFile("BadMap.root","recreate") ;
//   hBadMap->Write() ;
//   fout->Close() ;

  // put now result to local CDB
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://./");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment(Form("Random bad channel map with %.3f bad channels",fractionBad));
  AliCDBId id("PHOS/Calib/EmcBadChannels",0,999999);
  CDB->Put(&badMap,id, md);
 

}
