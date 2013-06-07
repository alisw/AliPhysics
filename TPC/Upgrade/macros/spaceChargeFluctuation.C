/*
  .x ~/rootlogon.C
  .L $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C+ 
//  GenerateMapRaw("/hera/alice/local/filtered/alice/data/2010/LHC10e/000129527/raw/merged_highpt_1.root","output.root", 50);
  GenerateMapRawIons( "/hera/alice/local/filtered/alice/data/2010/LHC10e/000129527/raw/merged_highpt_1.root","output.root", 25);
  
 */
#include "TMath.h"
#include "TRandom.h"
#include "TTreeStream.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "AliTPCParam.h"
#include "AliTPCcalibDB.h"
#include "AliTPCAltroMapping.h"
#include "AliAltroRawStream.h"
#include "AliSysInfo.h"
#include "AliTPCRawStreamV3.h"
#include "AliCDBManager.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliRawReaderRoot.h"
#include "AliRawReader.h"
#include "TH3.h"
#include "TH2.h"


TH1 * GenerateMapRawIons(const char *fileName="raw.root", const char *outputName="histo.root", Int_t maxEvents=25);


void spaceChargeFluctuation(Int_t mode=0){
  //
  //
  if (mode==0) GenerateMapRawIons();
  
}

void spaceChargeFluctuationMC(Int_t nframes){
  //
  // Toy MC to generate space charge fluctuation
  //
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("spaceChargeFluctuation.root","recreate");
  Double_t interactionRate=100000;
  Double_t driftTime=0.1;

  for (Int_t iframe=0; iframe<nframes; iframe++){
    Int_t nevents=gRandom->Poisson(interactionRate*driftTime);
    Int_t ntracksAll=0;
    TVectorD vecAllPhi128(128);
    TVectorD vecAllPhi36(36);
    for (Int_t ievent=0; ievent<nevents; ievent++){
      //
      //
      Int_t ntracks=gRandom->Exp(500)+gRandom->Exp(500); // to be taken form the distribution          
      ntracksAll+=ntracks; 
      for (Int_t isec=0; isec<128; isec++){
	vecAllPhi128(isec)+=ntracks*gRandom->Rndm()/128;
      }
      for (Int_t isec=0; isec<36; isec++){
	vecAllPhi36(isec)+=ntracks*gRandom->Rndm()/36;
      }
    }
    (*pcstream)<<"ntracks"<<
      "nevents="<<nevents<<                  // number of events withing time frame
      "ntracksAll="<<ntracksAll<<            // number of tracks within  time frame
      "vecAllPhi36="<<&vecAllPhi36<<         // number of tracks in phi bin (36 bins)    within time frame
      "vecAllPhi128="<<&vecAllPhi128<<       // number of tracks in phi bin (128 bins)   within time frame
      "\n";
  }
  delete pcstream;
}


// void spaceChargeFluctuationDraw(){
//   //
//   //
//   //
//   TFile *f = TFile::Open("spaceChargeFluctuation.root");
//   TTree * tree = (TTree*)f->Get("ntracks");
//   TCanvas * canvasFluc = new TCanvas("cavasFluc","canvasFluc");  
// }


TH1 * GenerateMapRaw(const char *fileName, const char *outputName, Int_t maxEvents){
  //
  // Generate 3D maps of electrons 
  // different threshold considered
  // Paramaters:
  //    fileName      - name of input raw file
  //    outputName    - name of output file with the space charge histograms 
  //    maxEvents     - grouping of the events
  // 
  //  
  TTreeSRedirector * pcstream  = new TTreeSRedirector(outputName, "recreate");

  const char *  ocdbpath =gSystem->Getenv("OCDB_PATH") ? gSystem->Getenv("OCDB_PATH"):"local://$ALICE_ROOT/OCDB/";  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,       AliMagF::kBeamTypepp, 2.76/2.));
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  AliTPCParam * param = AliTPCcalibDB::Instance()->GetParameters();
  //
  TStopwatch timer;
  timer.Start();
  TH3F * hisQ3D[3]={0};
  TH2F * hisQ2D[3]={0};
  TH1F * hisQ1D[3]={0};
  //
  for (Int_t ith=0; ith<3; ith++){
    char chname[100];
    snprintf(chname,100,"hisQ3D_Th%d",2*ith+2);
    hisQ3D[ith] = new TH3F(chname, chname,180, 0,2*TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,100,0,1100);
    snprintf(chname,100,"hisQ2D_Th%d",2*ith+2);
    hisQ2D[ith] = new TH2F(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    //
    snprintf(chname,100,"hisQ1D_Th%d",2*ith+2);
    hisQ1D[ith] = new TH1F(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    hisQ3D[ith]->SetDirectory(0);
    hisQ1D[ith]->SetDirectory(0);
    hisQ1D[ith]->SetDirectory(0);
  }
  //
  //  
  AliRawReader *reader = new AliRawReaderRoot(fileName);
  reader->Reset();
  AliAltroRawStream* stream = new AliAltroRawStream(reader);
  stream->SelectRawData("TPC");
  Int_t evtnr=0;
  Int_t chunkNr=0;
  //
  while (reader->NextEvent()) {
    if(evtnr>maxEvents) {
      chunkNr++;
      evtnr=0;
      pcstream->GetFile()->mkdir(Form("Chunk%d",chunkNr));
      pcstream->GetFile()->cd(Form("Chunk%d",chunkNr));
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Write(Form("His1D_%d",ith));
	hisQ2D[ith]->Write(Form("His2D_%d",ith));
	hisQ3D[ith]->Write(Form("His3D_%d",ith));
	(*pcstream)<<"histo"<<
	  Form("hist1D_%d.=",ith*2+2)<<hisQ1D[ith]<<
	  Form("hist2D_%d.=",ith*2+2)<<hisQ2D[ith]<<
	  Form("hist3D_%d.=",ith*2+2)<<hisQ2D[ith]<<
	  "\n";
	hisQ1D[ith]->Reset();
	hisQ2D[ith]->Reset();
	hisQ3D[ith]->Reset();
      }
    }
    cout<<"Chunk=\t"<<chunkNr<<"\tEvt=\t"<<evtnr<<endl;
    evtnr++;
    AliSysInfo::AddStamp(Form("Event%d",evtnr),evtnr);
    AliTPCRawStreamV3 input(reader,(AliAltroMapping**)mapping);
    //
    while (input.NextDDL()){
      Int_t sector = input.GetSector(); 
      while ( input.NextChannel() ) {
	Int_t    row    = input.GetRow();
	Int_t rowAbs    = row+ ((sector<36) ? 0: param->GetNRow(0));
	Int_t    pad    = input.GetPad();
	Int_t    nPads   = param->GetNPads(sector,row);
	Double_t localX  = param->GetPadRowRadii(sector,row); 
	Double_t localY  = (pad-nPads)*param->GetPadPitchWidth(sector);
	Double_t localPhi= TMath::ATan2(localY,localX);
	Double_t phi     = TMath::Pi()*((sector%18)+0.5)/9+localPhi;
	if (sector%36>=18) phi+=TMath::TwoPi();
	Double_t padLength=param->GetPadPitchLength(sector,row);
	//      Double_t padWidth=param->GetPadPitchWidth(sector);
	//    Double_t zLength=param->GetZLength();
	Double_t deltaPhi=hisQ3D[0]->GetXaxis()->GetBinWidth(0);
        Double_t deltaZ= hisQ3D[0]->GetZaxis()->GetBinWidth(0);
	//
	while ( input.NextBunch() ){
	  Int_t  startTbin    = (Int_t)input.GetStartTimeBin();
	  Int_t  bunchlength  = (Int_t)input.GetBunchLength();
	  const UShort_t *sig = input.GetSignals();	  
	  Int_t aboveTh[3]={0};
	  for (Int_t i=0; i<bunchlength; i++){ 
	    for (Int_t ith=0; ith<3; ith++){
	      if (sig[i]>(ith*2)+2) aboveTh[ith]++; 
	    }
	  }
	  for (Int_t ith=0; ith<3; ith++){
	    if (aboveTh[ith%3]>1){
	      for (Int_t i=0; i<bunchlength; i++){
		//
		// normalization
		//
		Double_t volume3D = padLength*(localX*deltaPhi)*deltaZ;
		if (sector%36<18) hisQ1D[ith]->Fill(rowAbs, sig[i]/volume3D);
		hisQ2D[ith]->Fill(phi,rowAbs, sig[i]/volume3D);
		hisQ3D[ith]->Fill(phi,rowAbs,startTbin+i,sig[i]/volume3D);
	      }
	    }
	  }
	}
      }
    }
  }
  timer.Print();
  delete pcstream;
  return 0;
}





TH1 * GenerateMapRawIons(const char *fileName, const char *outputName, Int_t maxEvents){
  //
  // Generate 3D maps of the space charge for the rad data maps
  // different threshold considered
  // Paramaters:
  //    fileName      - name of input raw file
  //    outputName    - name of output file with the space charge histograms 
  //    maxEvents     - grouping of the events
  // 
  //  
  TTreeSRedirector * pcstream  = new TTreeSRedirector(outputName, "recreate");
  const char *  ocdbpath =gSystem->Getenv("OCDB_PATH") ? gSystem->Getenv("OCDB_PATH"):"local://$ALICE_ROOT/OCDB/";  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,       AliMagF::kBeamTypepp, 2.76/2.));
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  AliTPCParam * param = AliTPCcalibDB::Instance()->GetParameters();
  TStopwatch timer;
  timer.Start();
  //   arrays of space charges - different elements corresponds to different threshold to accumulate charge
  TH3F * hisQ3D[3]={0};                // 3D maps space charge from drift volume
  TH2F * hisQ2D[3]={0};                
  TH1F * hisQ1D[3]={0};
  TH3F * hisQ3DROC[3]={0};             // 3D maps space charge from ROC
  TH2F * hisQ2DROC[3]={0};
  TH1F * hisQ1DROC[3]={0};
  Int_t nbinsRow=param->GetNRowLow()+param->GetNRowUp();
  Double_t *xbins = new Double_t[nbinsRow+1];
  xbins[0]=param->GetPadRowRadiiLow(0)-1;   //underflow bin
  for (Int_t ibin=0; ibin<param->GetNRowLow();ibin++) xbins[1+ibin]=param->GetPadRowRadiiLow(ibin);
  for (Int_t ibin=0; ibin<param->GetNRowUp();ibin++)  xbins[1+ibin+param->GetNRowLow()]=param->GetPadRowRadiiUp(ibin);
  

  //
  for (Int_t ith=0; ith<3; ith++){
    char chname[100];
    snprintf(chname,100,"hisQ3D_Th%d",2*ith+2);
    hisQ3D[ith] = new TH3F(chname, chname,180, 0,2*TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    snprintf(chname,100,"hisQ2D_Th%d",2*ith+2);
    hisQ2D[ith] = new TH2F(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ1D_Th%d",2*ith+2);
    hisQ1D[ith] = new TH1F(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    //
    snprintf(chname,100,"hisQ3DROC_Th%d",2*ith+2);
    hisQ3DROC[ith] = new TH3F(chname, chname,180, 0,2*TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    snprintf(chname,100,"hisQ2DROC_Th%d",2*ith+2);
    hisQ2DROC[ith] = new TH2F(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ1DROC_Th%d",2*ith+2);
    hisQ1DROC[ith] = new TH1F(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    hisQ1D[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ2D[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ3D[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ1DROC[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ2DROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ3DROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    //
    hisQ3D[ith]->SetDirectory(0);
    hisQ1D[ith]->SetDirectory(0);
    hisQ1D[ith]->SetDirectory(0);
    hisQ3DROC[ith]->SetDirectory(0);
    hisQ1DROC[ith]->SetDirectory(0);
    hisQ1DROC[ith]->SetDirectory(0);
  }
  //
  //  
  AliRawReader *reader = new AliRawReaderRoot(fileName);
  reader->Reset();
  AliAltroRawStream* stream = new AliAltroRawStream(reader);
  stream->SelectRawData("TPC");
  Int_t evtnr=0;
  Int_t chunkNr=0;
  // 

  while (reader->NextEvent()) {
    Double_t shiftZ= gRandom->Rndm()*250;
    //
    if(evtnr>=maxEvents) {
      chunkNr++;
      evtnr=0;
      pcstream->GetFile()->mkdir(Form("Chunk%d",chunkNr));
      pcstream->GetFile()->cd(Form("Chunk%d",chunkNr));
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Write(Form("His1DDrift_%d",ith));
	hisQ2D[ith]->Write(Form("His2DDrift_%d",ith));
	hisQ3D[ith]->Write(Form("His3DDrift_%d",ith));
	hisQ1DROC[ith]->Write(Form("His1DROC_%d",ith));
	hisQ2DROC[ith]->Write(Form("His2DROC_%d",ith));
	hisQ3DROC[ith]->Write(Form("His3DROC_%d",ith));
	(*pcstream)<<"histo"<<
	  Form("hist1D_%d.=",ith*2+2)<<hisQ1D[ith]<<
	  Form("hist2D_%d.=",ith*2+2)<<hisQ2D[ith]<<
	  Form("hist3D_%d.=",ith*2+2)<<hisQ2D[ith]<<
	  Form("hist1DROC_%d.=",ith*2+2)<<hisQ1DROC[ith]<<
	  Form("hist2DROC_%d.=",ith*2+2)<<hisQ2DROC[ith]<<
	  Form("hist3DROC_%d.=",ith*2+2)<<hisQ2DROC[ith];	  
      }
      (*pcstream)<<"histo"<<"\n";
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Reset();
	hisQ2D[ith]->Reset();
	hisQ3D[ith]->Reset();
	hisQ1DROC[ith]->Reset();
	hisQ2DROC[ith]->Reset();
	hisQ3DROC[ith]->Reset();
      }
    }
    cout<<"Chunk=\t"<<chunkNr<<"\tEvt=\t"<<evtnr<<endl;
    evtnr++;
    AliSysInfo::AddStamp(Form("Event%d",evtnr),evtnr);
    AliTPCRawStreamV3 input(reader,(AliAltroMapping**)mapping);
    //
    while (input.NextDDL()){
      Int_t sector = input.GetSector(); 
      while ( input.NextChannel() ) {
	Int_t    row    = input.GetRow();
	Int_t rowAbs    = row+ ((sector<36) ? 0: param->GetNRow(0));
	Int_t    pad    = input.GetPad();
	Int_t    nPads   = param->GetNPads(sector,row);
	Double_t localX  = param->GetPadRowRadii(sector,row); 
	Double_t localY  = (pad-nPads)*param->GetPadPitchWidth(sector);
	Double_t localPhi= TMath::ATan2(localY,localX);
	Double_t phi     = TMath::Pi()*((sector%18)+0.5)/9+localPhi;
	if (sector%36>=18) phi+=TMath::TwoPi();
	Double_t padLength=param->GetPadPitchLength(sector,row);
	//  Double_t padWidth=param->GetPadPitchWidth(sector);
	//  Double_t zLength=param->GetZLength();
	Double_t deltaPhi=hisQ3D[0]->GetXaxis()->GetBinWidth(0);
        Double_t deltaZ= hisQ3D[0]->GetZaxis()->GetBinWidth(0);
	//
	while ( input.NextBunch() ){
	  Int_t  startTbin    = (Int_t)input.GetStartTimeBin();
	  Int_t  bunchlength  = (Int_t)input.GetBunchLength();
	  const UShort_t *sig = input.GetSignals();	  
	  Int_t aboveTh[3]={0};
	  for (Int_t i=0; i<bunchlength; i++){ 
	    for (Int_t ith=0; ith<3; ith++){
	      if (sig[i]>(ith*2)+2) aboveTh[ith]++; 
	    }
	  }
	  for (Int_t ith=0; ith<3; ith++){
	    if (aboveTh[ith%3]>1){
	      for (Int_t i=0; i<bunchlength; i++){
		//
		// normalization
		//
		Double_t volume3D  = padLength*(localX*deltaPhi)*deltaZ;
		Double_t zIonDrift = (param->GetZLength()-startTbin*param->GetZWidth());
		zIonDrift+=shiftZ;
                if (TMath::Abs(zIonDrift)<param->GetZLength()){
		  if ((sector%36)>=18) zIonDrift*=-1;   // c side has opposite sign
		  if (sector%36<18) hisQ1D[ith]->Fill(localX, sig[i]/volume3D);
		  hisQ2D[ith]->Fill(phi,localX, sig[i]/volume3D);
		  hisQ3D[ith]->Fill(phi,localX,zIonDrift,sig[i]/volume3D);
		}
		//
		Double_t zIonROC = ((sector%36)<18)? shiftZ: -shiftZ;  // z position of the "ion disc" -  A side C side opposite sign
		if (sector%36<18) hisQ1DROC[ith]->Fill(localX, sig[i]/volume3D);
		hisQ2DROC[ith]->Fill(phi,localX, sig[i]/volume3D);
		hisQ3DROC[ith]->Fill(phi,localX,zIonROC,sig[i]/volume3D);
	      }
	    }
	  }
	}
      }
    }
  }
  timer.Print();
  delete pcstream;
  return 0;
}
