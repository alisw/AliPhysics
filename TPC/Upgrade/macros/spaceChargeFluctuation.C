/*
  .x $HOME/rootlogon.C
  .L $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C+ 
//  GenerateMapRaw("/hera/alice/local/filtered/alice/data/2010/LHC10e/000129527/raw/merged_highpt_1.root","output.root", 50);
  GenerateMapRawIons( kFALSE, "/hera/alice/local/filtered/alice/data/2010/LHC10e/000129527/raw/merged_highpt_1.root","output.root", 25);
  
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
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
#include "TLegend.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TStatToolkit.h"

#include "AliDCSSensor.h"
#include "AliCDBEntry.h"
#include "AliDCSSensorArray.h"

TH1 * GenerateMapRawIons(Int_t useGain,const char *fileName="raw.root", const char *outputName="histo.root", Int_t maxEvents=25);
void DoMerge();

void spaceChargeFluctuation(Int_t mode=0, Int_t arg0=0){
  //
  //
  if (mode==0) GenerateMapRawIons(arg0);  
  if (mode==1) DoMerge();  
}

void spaceChargeFluctuationMC(Int_t nframes){
  //
  // Toy MC to generate space charge fluctuation
  //
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("spaceChargeFluctuation.root","recreate");
  Double_t interactionRate=50000;
  Double_t driftTime=0.1;              
  const Double_t kB2C=-0.299792458e-3;
  Double_t eventMean=interactionRate*driftTime;
  Double_t trackMean=500;
  Double_t FPOT=1.0, EEND=3000;
  Double_t  EEXPO=0.8567;
  const Double_t XEXPO=-EEXPO+1, YEXPO=1/XEXPO;

  for (Int_t iframe=0; iframe<nframes; iframe++){
    printf("iframe=%d\n",iframe);
    Int_t nevents=gRandom->Poisson(interactionRate*driftTime);
    Int_t ntracksAll=0;
    TVectorD vecTracksPhi180(180);
    TVectorD vecTracksPhi36(36);
    TVectorD vecEPhi180(180);
    TVectorD vecEPhi36(36);
    Double_t dESum=0;
    for (Int_t ievent=0; ievent<nevents; ievent++){
      Int_t ntracks=gRandom->Exp(trackMean); // to be taken from the MB primary distribution      
      Float_t RAN = gRandom->Rndm();
      ntracks=TMath::Power((TMath::Power(FPOT,XEXPO)*(1-RAN)+TMath::Power(EEND,XEXPO)*RAN),YEXPO)/2.;
      ntracksAll+=ntracks; 
      for (Int_t itrack=0; itrack<ntracks; itrack++){
	Double_t phi  = gRandom->Rndm();
	vecTracksPhi180(Int_t(phi*180))+=1;
	vecTracksPhi36(Int_t(phi*36))  +=1;
	// simplified MC to get track length including loopers
	Double_t theta= gRandom->Rndm();
	Double_t pt   = gRandom->Exp(0.5)+0.05;
	Double_t crv  = TMath::Abs(5*kB2C/pt);   //GetC(b); // bz*kB2C/pt;
	Double_t deltaPhi=0;
	if (TMath::Abs(2*crv*(245-85)/2.) <1.) deltaPhi=TMath::ASin(crv*(245-85)/2.);
	else 
	  deltaPhi=TMath::Pi();
	Double_t dE=deltaPhi/crv;
	Double_t xloop=1;
	if (1./crv<250) {
	  xloop = TMath::Min(1./(TMath::Abs(theta)+0.0001),10.);
	  if (xloop<1) xloop=1;
	}
	dESum+=xloop*dE;
	if (itrack==0) (*pcstream)<<"track"<<
	  "pt="<<pt<<
	  "crv="<<crv<<
	  "theta="<<theta<<
	  "dE="<<dE<<
	  "xloop="<<xloop<<
	  "\n";
	
	vecEPhi180(Int_t(phi*180))     +=dE*xloop;
	vecEPhi36(Int_t(phi*36))       +=dE*xloop;
      }
      (*pcstream)<<"event"<<
	"ntracks="<<ntracks<<
	"nevents="<<nevents<<
	"\n";
    }
    (*pcstream)<<"ntracks"<<
      "eventMean="<<eventMean<<
      "trackMean="<<trackMean<<
      "nevents="<<nevents<<                       // number of events withing time frame
      "ntracksAll="<<ntracksAll<<                  // number of tracks within  time frame
      "dESum="<<dESum<<                            // sum of the energy loss
      "vecTracksPhi36.="<<&vecTracksPhi36<<         // number of tracks in phi bin (36 bins)    within time frame
      "vecTracksPhi180.="<<&vecTracksPhi180<<       // number of tracks in phi bin (180 bins)   within time frame
      "vecEPhi36.="<<&vecEPhi36<<         // number of tracks in phi bin (36 bins)    within time frame
      "vecEPhi180.="<<&vecEPhi180<<       // number of tracks in phi bin (180 bins)   within time frame
      "\n";
  }
  delete pcstream;
  //
  // Toy MC to simulate the space charge integral fluctuation
  //
  TFile * f = TFile::Open("spaceChargeFluctuation.root");
  TTree * treeStat = (TTree*)f->Get("ntracks");
  TTree * treedE = (TTree*)f->Get("track");
  TTree * treeEv = (TTree*)f->Get("event");
  Int_t nentries=treedE->Draw("dE*xloop","1","",1000000);
  Double_t meandE=TMath::Mean(nentries,treedE->GetV1());
  Double_t rmsdE=TMath::RMS(nentries,treedE->GetV1());
  treeStat->SetAlias("meandE",Form("(%f+0)",meandE));
  treeStat->SetAlias("rmsdE",Form("(%f+0)",rmsdE));
  nentries=treeEv->Draw("ntracks","1","",1000000);
  Double_t meanT=TMath::Mean(nentries,treeEv->GetV1());
  Double_t rmsT=TMath::RMS(nentries,treeEv->GetV1());
  treeStat->SetAlias("tracksMean",Form("(%f+0)",meanT));
  treeStat->SetAlias("tracksRMS",Form("(%f+0)",rmsT));
  //
  treeStat->Draw("nevents/eventMean>>hisEv(100,0.9,1.1)","");
  treeStat->Draw("ntracksAll/(eventMean*tracksMean)>>hisTrackAll(100,0.9,1)","","same");
  treeStat->Draw("vecTracksPhi180.fElements/(eventMean*tracksMean/180)>>hisTrackSector(100,0.9,1)","1/180","same");
  treeStat->Draw("vecEPhi180.fElements/(eventMean*tracksMean*meandE/180)>>hisdESector(100,0.9,1)","1/180","same");
  //
  //
  //
  treeStat->Draw("(ntracksAll/(eventMean*tracksMean)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean)");  //tracks All pull
  treeStat->Draw("(vecTracksPhi180.fElements/(eventMean*tracksMean/180.)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+180/(tracksMean*eventMean))");  //tracks spread
  treeStat->Draw("(vecEPhi180.fElements/(eventMean*tracksMean*meandE/180)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+180/(tracksMean*eventMean)+180*(rmsdE/meandE)**2/(eventMean*tracksMean))"); //dE spread

}
void FitMultiplicity(){
  //
  // Fit multiplicity distribution
  //
  TFile *fmult=TFile::Open("mult_dist_pbpb.root");
  TF1 f1("f1","[0]*(x+abs([2]))**(-abs([1]))",1,3000);
  TH1* his = (TH1*) fmult->Get("mult_dist_PbPb_normalizedbywidth");
  f1.SetParameters(his->GetEntries(),1,1);
  his->Fit(&f1,"","",2,3000);
  
  Double_t FPOT=1.0, EEND=3000, EEXPO= TMath::Abs(f1.GetParameter(1));
  EEXPO=0.8567;
  const Double_t XEXPO=-EEXPO+1, YEXPO=1/XEXPO;
  TH1F *hisr= new TH1F("aaa","aaa",4000,0,4000);
  for (Int_t i=0; i<400000; i++){
    Float_t RAN = gRandom->Rndm();
    hisr->Fill(TMath::Power((TMath::Power(FPOT,XEXPO)*(1-RAN)+TMath::Power(EEND,XEXPO)*RAN),YEXPO));
  }
  
}

// void spaceChargeFluctuationDraw(){
//   //
//   //
//   //
//   TFile *f = TFile::Open("spaceChargeFluctuation.root");
//   TTree * tree = (TTree*)f->Get("ntracks");
//   TCanvas * canvasFluc = new TCanvas("cavasFluc","canvasFluc");  
// }






TH1 * GenerateMapRawIons(Int_t useGainMap, const char *fileName, const char *outputName, Int_t maxEvents){
  //
  // Generate 3D maps of the space charge for the rad data maps
  // different threshold considered
  // Paramaters:
  //    useGainMap    - switch usage of the gain map
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
  AliTPCCalPad * gain = AliTPCcalibDB::Instance()->GetDedxGainFactor();
  AliTPCCalPad * noise = AliTPCcalibDB::Instance()->GetPadNoise();

  TStopwatch timer;
  timer.Start();
  //   arrays of space charges - different elements corresponds to different threshold to accumulate charge
  TH3D * hisQ3D[3]={0};                // 3D maps space charge from drift volume
  TH2D * hisQ2D[3]={0};                
  TH2D * hisQ2DRZ[3]={0};                
  TH1D * hisQ1D[3]={0};
  TH3D * hisQ3DROC[3]={0};             // 3D maps space charge from ROC
  TH2D * hisQ2DROC[3]={0};
  TH2D * hisQ2DRZROC[3]={0};                
  TH1D * hisQ1DROC[3]={0};
  Int_t nbinsRow=param->GetNRowLow()+param->GetNRowUp();
  Double_t *xbins = new Double_t[nbinsRow+1];
  xbins[0]=param->GetPadRowRadiiLow(0)-1;   //underflow bin
  for (Int_t ibin=0; ibin<param->GetNRowLow();ibin++) xbins[1+ibin]=param->GetPadRowRadiiLow(ibin);
  for (Int_t ibin=0; ibin<param->GetNRowUp();ibin++)  xbins[1+ibin+param->GetNRowLow()]=param->GetPadRowRadiiUp(ibin);
  

  //
  for (Int_t ith=0; ith<3; ith++){
    char chname[100];
    snprintf(chname,100,"hisQ3D_Th%d",2*ith+2);
    hisQ3D[ith] = new TH3D(chname, chname,180, 0,2*TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    snprintf(chname,100,"hisQ2D_Th%d",2*ith+2);
    hisQ2D[ith] = new TH2D(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ2DRZ_Th%d",2*ith+2);

    hisQ2DRZ[ith] = new TH2D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36),  125,-250,250);

    snprintf(chname,100,"hisQ1D_Th%d",2*ith+2);
    hisQ1D[ith] = new TH1D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    //
    snprintf(chname,100,"hisQ3DROC_Th%d",2*ith+2);
    hisQ3DROC[ith] = new TH3D(chname, chname,180, 0,2*TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    snprintf(chname,100,"hisQ2DROC_Th%d",2*ith+2);
    hisQ2DROC[ith] = new TH2D(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ2DRZROC_Th%d",2*ith+2);
    hisQ2DRZROC[ith] = new TH2D(chname,chname,param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36), 125,-250,250);
    snprintf(chname,100,"hisQ1DROC_Th%d",2*ith+2);
    hisQ1DROC[ith] = new TH1D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    hisQ1D[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ2D[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ2DRZ[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ3D[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ1DROC[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ2DROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ2DRZROC[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ3DROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    //
    hisQ3D[ith]->SetDirectory(0);
    hisQ2D[ith]->SetDirectory(0);
    hisQ2DRZ[ith]->SetDirectory(0);
    hisQ1D[ith]->SetDirectory(0);
    hisQ3DROC[ith]->SetDirectory(0);
    hisQ2DRZROC[ith]->SetDirectory(0);
    hisQ2DROC[ith]->SetDirectory(0);
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
    Double_t shiftZ= gRandom->Rndm()*250.;
    //
    if(evtnr>=maxEvents) {
      chunkNr++;
      pcstream->GetFile()->mkdir(Form("Chunk%d",chunkNr));
      pcstream->GetFile()->cd(Form("Chunk%d",chunkNr));
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Write(Form("His1DDrift_%d",ith));
	hisQ2D[ith]->Write(Form("His2DDrift_%d",ith));
	hisQ2DRZ[ith]->Write(Form("His2DDriftRZ_%d",ith));
	hisQ3D[ith]->Write(Form("His3DDrift_%d",ith));
	hisQ1DROC[ith]->Write(Form("His1DROC_%d",ith));
	hisQ2DROC[ith]->Write(Form("His2DROC_%d",ith));
	hisQ2DRZROC[ith]->Write(Form("His2DRZROC_%d",ith));
	hisQ3DROC[ith]->Write(Form("His3DROC_%d",ith));
	(*pcstream)<<"histo"<<
	  "events="<<evtnr<<
	  "useGain="<<useGainMap<<
	  Form("hist1D_%d.=",ith*2+2)<<hisQ1D[ith]<<
	  Form("hist2D_%d.=",ith*2+2)<<hisQ2D[ith]<<
	  Form("hist2DRZ_%d.=",ith*2+2)<<hisQ2DRZ[ith]<<
	  Form("hist3D_%d.=",ith*2+2)<<hisQ3D[ith]<<
	  Form("hist1DROC_%d.=",ith*2+2)<<hisQ1DROC[ith]<<
	  Form("hist2DROC_%d.=",ith*2+2)<<hisQ2DROC[ith]<<
	  Form("hist2DRZROC_%d.=",ith*2+2)<<hisQ2DRZROC[ith]<<
	  Form("hist3DROC_%d.=",ith*2+2)<<hisQ3DROC[ith];	  
      }
      (*pcstream)<<"histo"<<"\n";
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Reset();
	hisQ2D[ith]->Reset();
	hisQ2DRZ[ith]->Reset();
	hisQ3D[ith]->Reset();
	hisQ1DROC[ith]->Reset();
	hisQ2DROC[ith]->Reset();
	hisQ2DRZROC[ith]->Reset();
	hisQ3DROC[ith]->Reset();
      }
      evtnr=0;
    }
    cout<<"Chunk=\t"<<chunkNr<<"\tEvt=\t"<<evtnr<<endl;
    evtnr++;
    AliSysInfo::AddStamp(Form("Event%d",evtnr),evtnr);
    AliTPCRawStreamV3 input(reader,(AliAltroMapping**)mapping);
    //
    while (input.NextDDL()){
      Int_t sector = input.GetSector();  
      AliTPCCalROC * gainROC =gain->GetCalROC(sector);
      AliTPCCalROC * noiseROC =noise->GetCalROC(sector);
      while ( input.NextChannel() ) {
	Int_t    row    = input.GetRow();
	//	Int_t rowAbs    = row+ ((sector<36) ? 0: param->GetNRow(0));
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
	//	Double_t deltaPhi=hisQ3D[0]->GetXaxis()->GetBinWidth(0);
        //Double_t deltaZ= hisQ3D[0]->GetZaxis()->GetBinWidth(0);
	Double_t gainPad = gainROC->GetValue(row,pad); 
	Double_t noisePad = noiseROC->GetValue(row,pad); 
	//
	while ( input.NextBunch() ){
	  Int_t  startTbin    = (Int_t)input.GetStartTimeBin();
	  Int_t  bunchlength  = (Int_t)input.GetBunchLength();
	  const UShort_t *sig = input.GetSignals();	  
	  Int_t aboveTh[3]={0};
	  for (Int_t i=0; i<bunchlength; i++){ 
	    if (sig[i]<4*noisePad) continue;	    
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
		Double_t zIonDrift   =(param->GetZLength()-startTbin*param->GetZWidth());
		zIonDrift+=shiftZ;
		Double_t signal=sig[i];
		if (useGainMap) signal/=gainPad;
                if (TMath::Abs(zIonDrift)<param->GetZLength()){
		  if ((sector%36)>=18) zIonDrift*=-1;   // c side has opposite sign
		  if (sector%36<18) hisQ1D[ith]->Fill(localX, signal/padLength);
		  hisQ2D[ith]->Fill(phi,localX, signal/padLength);
		  hisQ2DRZ[ith]->Fill(localX, zIonDrift, signal/padLength);
		  hisQ3D[ith]->Fill(phi,localX,zIonDrift,signal/padLength);
		}
		//
		Double_t zIonROC = ((sector%36)<18)? shiftZ: -shiftZ;  // z position of the "ion disc" -  A side C side opposite sign
		if (sector%36<18) hisQ1DROC[ith]->Fill(localX, signal/padLength);
		hisQ2DROC[ith]->Fill(phi,localX, signal/padLength);
		hisQ2DRZROC[ith]->Fill(localX, zIonROC, signal/padLength);
		hisQ3DROC[ith]->Fill(phi,localX,zIonROC,signal/padLength);
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


void AnalyzeMaps1D(){
  //
  // Analyze space charge maps stored as shitograms in trees
  //
  TFile *  fhisto = new TFile("histo.root","update");
  TTree * tree = (TTree*)fhisto->Get("histo");
  if (!tree){
    TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom("histo.list","histo",0,1000,1);
    tree = chain->CopyTree("1");
    tree->Write("histo");
  }
  //
  //
  //
  TH1 *his1Th[3]={0,0,0};
  TF1 *fq1DStep= new TF1("fq1DStep","([0]+[1]*(x>134))/x**min(abs([2]),3)",85,245);  
  fq1DStep->SetParameters(1,-0.5,1);
  tree->Draw("hist1DROC_2.fArray:hist1D_2.fXaxis.fXbins.fArray>>his(40,85,245)","","prof");
  tree->GetHistogram()->Fit(fq1DStep);
  // normalize step between the IROC-OROC
  tree->SetAlias("normQ",Form("(1+%f*(hist1D_2.fXaxis.fXbins.fArray>136))",fq1DStep->GetParameter(1)/fq1DStep->GetParameter(0)));
  //
  {
    Int_t entries= tree->Draw("hist1DROC_2.fArray/(events*normQ)","1","goff");
    Double_t median=TMath::Median(entries,tree->GetV1());
    TCut cut10Median = Form("hist1DROC_2.fArray/(events*normQ)<%f",10*median);
    //
    tree->Draw("hist1DROC_2.fArray/(events*normQ):hist1D_2.fXaxis.fXbins.fArray>>his1Th0(40,86,245)",cut10Median+"","prof");
    his1Th[0] = tree->GetHistogram();
    tree->Draw("hist1DROC_4.fArray/(events*normQ):hist1D_2.fXaxis.fXbins.fArray>>his1Th1(40,86,245)",cut10Median+"","prof");
    his1Th[1] = tree->GetHistogram();
    tree->Draw("hist1DROC_6.fArray/(events*normQ):hist1D_2.fXaxis.fXbins.fArray>>his1Th2(40,86,245)",cut10Median+"","prof");
    his1Th[2]=tree->GetHistogram();
  }
  //
  TCanvas *canvasR = new TCanvas("canvasR","canvasR",600,500);
  canvasR->cd();
  for (Int_t i=0; i<3; i++){
    his1Th[i]->SetMarkerStyle(21);
    his1Th[i]->SetMarkerColor(i+2);
    fq1DStep->SetLineColor(i+2);
    his1Th[i]->Fit(fq1DStep,"","");
    his1Th[i]->GetXaxis()->SetTitle("r (cm)");
    his1Th[i]->GetYaxis()->SetTitle("#frac{N_{el}}{N_{ev}}(ADC/cm)");    
  }
  TLegend * legend  = new TLegend(0.11,0.11,0.7,0.39,"1D space Charge map (ROC part) (z,phi integrated)");
  for (Int_t i=0; i<3; i++){
    his1Th[i]->SetMinimum(0);fq1DStep->SetLineColor(i+2);
    his1Th[i]->Fit(fq1DStep,"qnr","qnr");
    if (i==0) his1Th[i]->Draw("");
    his1Th[i]->Draw("same");
    legend->AddEntry(his1Th[i],Form("Thr=%d Slope=%2.2f",2*i+2,fq1DStep->GetParameter(2)));
  }
  legend->Draw();
  legend->SaveAs("spaceCharge1d.png");
  //
  //
  //
}

void DoMerge(){
  //
  //
  //
  TFile *  fhisto = new TFile("histo.root","recreate");
  TTree * tree = 0;
  TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom("histo.list","histo",0,50,1);
  tree = chain->CopyTree("1");
  tree->Write("histo");
  delete fhisto;
}



void MakeFluctuationStudy3D(Int_t nhistos, Int_t nevents, Int_t niter){
  //
  //
  // 
  // Fix z binning before creating of official plot (PbPb  data)
  // 
  // 1. Make a summary integral   3D/2D/1D maps
  // 2. Create several maps with niter events  - Poisson flucturation in n
  // 3. Store results 3D maps in the tree (and also as histogram)  current and mean
  // 
  TTreeSRedirector * pcstream  = new TTreeSRedirector("histo.root", "update");
  TFile *  fhisto = pcstream->GetFile();
  
  TTree * tree = (TTree*)fhisto->Get("histo");
  tree->SetCacheSize(10000000000);

  TH1D * his1DROC=0,    * his1DROCSum=0,  * his1DROCN=0;
  TH1D * his1DDrift=0,  * his1DDriftSum=0, * his1DDriftN=0 ;
  TH2D * his2DROC=0,    * his2DROCSum=0,  * his2DROCN=0;
  TH2D * his2DRZROC=0,    * his2DRZROCSum=0,  * his2DRZROCN=0;
  TH2D * his2DDrift=0,  * his2DDriftSum=0, * his2DDriftN=0;  
  TH2D * his2DRZDrift=0,  * his2DRZDriftSum=0, * his2DRZDriftN=0;  
  TH3D * his3DROC=0,    * his3DROCSum=0,  * his3DROCN=0;
  TH3D * his3DDrift=0,  * his3DDriftSum=0, * his3DDriftN=0;
  //
  if (nhistos<0 || nhistos> tree->GetEntries()) nhistos = tree->GetEntries();
  Int_t  eventsPerChunk=0;
  tree->SetBranchAddress("hist1D_2.",&his1DDrift);
  tree->SetBranchAddress("hist1DROC_2.",&his1DROC);
  tree->SetBranchAddress("hist2D_2.",&his2DDrift);
  tree->SetBranchAddress("hist2DRZ_2.",&his2DRZDrift);
  tree->SetBranchAddress("hist2DROC_2.",&his2DROC);
  tree->SetBranchAddress("hist3D_2.",&his3DDrift);
  tree->SetBranchAddress("hist3DROC_2.",&his3DROC);
  tree->SetBranchAddress("hist2DRZROC_2.",&his2DRZROC);
  tree->SetBranchAddress("events",&eventsPerChunk);
  // 
  // 1. Make a summary integral   3D/2D/1D maps
  //
  Int_t neventsAll=0;
  for (Int_t i=0; i<nhistos; i++){
    tree->GetEntry(i);
    if (i%25==0) printf("%d\n",i);
    if (his1DROCSum==0)     his1DROCSum=new TH1D(*his1DROC);
    if (his1DDriftSum==0)   his1DDriftSum=new TH1D(*his1DDrift);
    if (his2DROCSum==0)     his2DROCSum=new TH2D(*his2DROC);
    if (his2DRZROCSum==0)     his2DRZROCSum=new TH2D(*his2DRZROC);
    if (his2DDriftSum==0)   his2DDriftSum=new TH2D(*his2DDrift);
    if (his2DRZDriftSum==0)   his2DRZDriftSum=new TH2D(*his2DRZDrift);
    if (his3DROCSum==0)     his3DROCSum=new TH3D(*his3DROC);
    if (his3DDriftSum==0)   his3DDriftSum=new TH3D(*his3DDrift);
    his1DROCSum->Add(his1DROC);
    his1DDriftSum->Add(his1DDrift);
    his2DROCSum->Add(his2DROC);
    his2DRZROCSum->Add(his2DRZROC);
    his2DDriftSum->Add(his2DDrift);
    his2DRZDriftSum->Add(his2DRZDrift);
    his3DROCSum->Add(his3DROC);
    his3DDriftSum->Add(his3DDrift);
    neventsAll+=eventsPerChunk;
  }
  //
  // 2. Create several maps with niter events  - Poisson flucturation in n
  //
  for (Int_t iter=0; iter<niter; iter++){
    printf("Itteration=\t%d\n",iter);
    Int_t ilast=gRandom->Rndm()*nhistos;
    Int_t nchunks=gRandom->Poisson(nevents)/eventsPerChunk;  // chunks with n typically 25 events
    for (Int_t i=0; i<nchunks; i++){
      //    ilast+=gRandom->Rndm()
      tree->GetEntry(gRandom->Rndm()*nhistos);
      if (i%10==0) printf("%d\t%d\n",iter, i);
      if (his1DROCN==0)     his1DROCN=new TH1D(*his1DROC);
      if (his1DDriftN==0)   his1DDriftN=new TH1D(*his1DDrift);
      if (his2DROCN==0)     his2DROCN=new TH2D(*his2DROC);
      if (his2DDriftN==0)   his2DDriftN=new TH2D(*his2DDrift);
      if (his2DRZROCN==0)     his2DRZROCN=new TH2D(*his2DRZROC);
      if (his2DRZDriftN==0)   his2DRZDriftN=new TH2D(*his2DRZDrift);
      if (his3DROCN==0)     his3DROCN=new TH3D(*his3DROC);
      if (his3DDriftN==0)   his3DDriftN=new TH3D(*his3DDrift);
      his1DROCN->Add(his1DROC);
      his1DDriftN->Add(his1DDrift);
      his2DROCN->Add(his2DROC);
      his2DRZDriftN->Add(his2DRZDrift);
      his2DRZROCN->Add(his2DRZROC);
      his2DDriftN->Add(his2DDrift);
      his3DROCN->Add(his3DROC);
      his3DDriftN->Add(his3DDrift);      
    } 
    //
    // 3. Store results 3D maps in the tree (and also as histogram)  current and mea
    //    
    Int_t eventsUsed=  nchunks*eventsPerChunk;
    (*pcstream)<<"fluctuation"<<
      "neventsAll="<<neventsAll<<   // total number of event to define mean
      "nmean="<<nevents<<         // mean number of events used
      "eventsUsed="<<eventsUsed<<         // number of chunks used for one fluct. study
      //
      // 1,2,3D histogram per group and total
      "his1DROCN.="<<his1DROCN<<
      "his1DROCSum.="<<his1DROCSum<<
      "his1DDriftN.="<<his1DDriftN<<
      "his1DDriftSum.="<<his1DDriftSum<<
      "his2DROCN.="<<his2DROCN<<
      "his2DROCSum.="<<his2DROCSum<<
      "his2DDriftN.="<<his2DDriftN<<
      "his2DDriftSum.="<<his2DDriftSum<<
      "his2DRZROCN.="<<his2DRZROCN<<
      "his2DRZROCSum.="<<his2DRZROCSum<<
      "his2DRZDriftN.="<<his2DRZDriftN<<
      "his2DRZDriftSum.="<<his2DRZDriftSum<<
      "his3DROCN.="<<his3DROCN<<
      "his3DROCSum.="<<his3DROCSum<<      
      "his3DDriftN.="<<his3DDriftN<<      
      "his3DDriftSum.="<<his3DDriftSum<<      
      "\n";      
    pcstream->GetFile()->mkdir(Form("Fluc%d",iter));
    pcstream->GetFile()->cd(Form("Fluc%d",iter));
    //
    his2DROCN->Write("his2DROCN");
    his2DRZROCN->Write("his2DRZROCN");
    //
    his2DROCSum->Write("his2DROCSum");        
    his2DRZROCSum->Write("his2DRZROCSum");
    //
    his2DDriftN->Write("his2DDriftN");
    his2DRZDriftN->Write("his2DRZDriftN");
    //
    his2DDriftSum->Write("his2DDriftSum");
    his2DRZDriftSum->Write("his2DRZDriftSum");
    //
    his3DROCN->Write("his3DROCN");
    his3DROCSum->Write("his3DROCSum");
    his3DDriftN->Write("his3DDriftN");
    his3DDriftSum->Write("his3DDriftSum");

    his1DROCN->Reset();
    his1DDriftN->Reset();
    his2DROCN->Reset();
    his2DRZDriftN->Reset();
    his2DRZROCN->Reset();
    his2DDriftN->Reset();
    his3DROCN->Reset();
    his3DDriftN->Reset();    
  }

  delete pcstream;


  //
  Double_t mean,rms;
  //
  mean=TMath::Median(his2DROCSum->GetNbinsX()*his2DROCSum->GetNbinsY(),his2DROCSum->GetArray());
  rms=TMath::Median(his2DROCSum->GetNbinsX()*his2DROCSum->GetNbinsY(),his2DROCSum->GetArray());
  his2DROCSum->SetMaximum(mean+4*rms);
  his2DROCSum->SetMinimum(mean-4*rms);

}


void DrawDCARPhiTrendTime(){
  //
  // Macros to draw the DCA correlation with the luminosity (estimated from the occupancy)
  //
  // A side and c side  0 differnt behaviour -
  // A side - space charge effect
  // C side - space charge effect+ FC charging: 
  //   Variables  to query from the QA/calibration DB - tree: 
  //   QA.TPC.CPass1.dcar_posA_0   -dca rphi in cm - offset
  //   Calib.TPC.occQA.Sum()       - luminosity is estimated using the mean occupancy per run
  //     
  TFile *fdb = TFile::Open("outAll.root");
  if (!fdb)  fdb = TFile::Open("http://www-alice.gsi.de/TPC/CPassMonitor/outAll.root"); 
  TTree * tree = (TTree*)fdb->Get("joinAll");
  tree->SetCacheSize(100000000);
  tree->SetMarkerStyle(25);
  
  //QA.TPC.CPass1.dcar_posA_0 QA.TPC.CPass1.dcar_posA_0_Err  QA.TPC.CPass1.meanMult  Calib.TPC.occQA.  DAQ.L3_magnetCurrent 
  
  TGraphErrors * grA = TStatToolkit::MakeGraphErrors(tree,"QA.TPC.CPass1.dcar_posA_0:Calib.TPC.occQA.Sum()*sign(DAQ.L3_magnetCurrent):2*QA.TPC.CPass1.dcar_posA_0_Err","run>190000&&QA.TPC.CPass1.status==1",25,2,0.5);
  TGraphErrors * grC = TStatToolkit::MakeGraphErrors(tree,"QA.TPC.CPass1.dcar_posC_0:Calib.TPC.occQA.Sum()*sign(DAQ.L3_magnetCurrent):2*QA.TPC.CPass1.dcar_posC_0_Err","run>190000&&QA.TPC.CPass1.status==1",25,4,0.5);
  Double_t mean,rms;
  TStatToolkit::EvaluateUni(grA->GetN(),grA->GetY(), mean,rms,grA->GetN()*0.8);
  grA->SetMinimum(mean-5*rms);
  grA->SetMaximum(mean+3*rms);
    
  
  grA->GetXaxis()->SetTitle("occ*sign(bz)");
  grA->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
  grA->Draw("ap");
  grC->Draw("p");
  TLegend* legend = new TLegend(0.11,0.11,0.5,0.3,"DCA_{rphi} as function of IR (2013)" );
  legend->AddEntry(grA,"A side","p");
  legend->AddEntry(grC,"C side","p");
  legend->Draw();
}



void DrawOpenGate(){
  //
  //  Make nice plot to demonstrate the space charge effect in run with the open gating grid
  //  For the moment the inmput is harwired - the CPass0 calibration data used
  //  Make nice drawing (with axis labels):
  //  To fix (longer term)
  //     the distortion map to be recalculated - using gaussian fit (currently we use mean)
  //     the histogram should be extended
  TFile f("/hera/alice/alien/alice/data/2013/LHC13g/000197470/cpass0/OCDB/root_archive.zip#meanITSVertex.root");
  TFile fref("/hera/alice/alien/alice/data/2013/LHC13g/000197584/cpass0/OCDB/root_archive.zip#meanITSVertex.root");
  //
  TTree * treeTOFdy=(TTree*)f.Get("TOFdy");
  TTree * treeTOFdyRef=(TTree*)fref.Get("TOFdy");
  treeTOFdy->AddFriend(treeTOFdyRef,"R");
  treeTOFdy->SetMarkerStyle(25);
  TTree * treeITSdy=(TTree*)f.Get("ITSdy");
  TTree * treeITSdyRef=(TTree*)fref.Get("ITSdy");
  treeITSdy->AddFriend(treeITSdyRef,"R");
  treeITSdy->SetMarkerStyle(25);
  TTree * treeVertexdy=(TTree*)f.Get("Vertexdy");
  TTree * treeVertexdyRef=(TTree*)fref.Get("Vertexdy");
  treeVertexdy->AddFriend(treeVertexdyRef,"R");
  treeVertexdy->SetMarkerStyle(25);

  //  treeITSdy->Draw("mean-R.mean:sector:abs(theta)","entries>50&&abs(snp)<0.1&&theta<0","colz")
  
  treeITSdy->Draw("mean-R.mean:sector:abs(theta)","entries>50&&abs(snp)<0.1&&theta>0","colz");
}


void DrawCurrent(const char * ocdb="/cvmfs/alice.gsi.de/alice/data/2013/OCDB/TPC/Calib/HighVoltage", Int_t run0=100000, Int_t run1=110000){
  //
  //
  //
  /*
    const char * ocdb="/cvmfs/alice.gsi.de/alice/data/2013/OCDB/TPC/Calib/HighVoltage";
    Int_t run0=197460;
    Int_t run1=197480;
  */
  const Int_t knpoints=100000;
  TVectorD vecTime(knpoints);
  TVectorD vecI(knpoints);
  Int_t npoints=0;
  for (Int_t irun=run0; irun<run1; irun++){
    TFile * f = TFile::Open(Form("%s/Run%d_%d_v1_s0.root",ocdb,irun,irun));
    if (!f) continue;
    AliCDBEntry *       entry = (AliCDBEntry *)f->Get("AliCDBEntry");    
    if (!entry) continue; 
    AliDCSSensorArray * array = (AliDCSSensorArray *)entry->GetObject();
    if (!array) continue;
    AliDCSSensor * sensor = array->GetSensor("TPC_VHV_D_I_MON");
    //sensor->Draw(Form("%d",irun));     
    TGraph *graph = sensor->GetGraph();
    for (Int_t ipoint=0; ipoint<graph->GetN(); ipoint++){
      vecTime[npoints]=sensor->GetStartTime()+graph->GetX()[ipoint]*3600;
      vecI[npoints]=graph->GetY()[ipoint];
      npoints++;
    }
  }
  TGraph * graph  = new TGraph(npoints, vecTime.GetMatrixArray(), vecI.GetMatrixArray());
  graph->Draw("alp");
  

}


void ConvertMaps(const char *inputfile="histo.root", const char *outputfile="SpaceCharge.root", Int_t event=0){
  //
  //
  //  const char *inputfile="histo.root";  const char *outputfile="SpaceCharge.root";  Int_t event=0
  //
  TFile *fhistos = TFile::Open(inputfile);
  TH2D *his2DROCN=0, *his2DROCSum=0, *his2DRZROCN=0, *his2DRZROCSum=0;
  TTree * treeHis = (TTree*)fhistos->Get("fluctuation");
  treeHis->SetBranchAddress("his2DROCN.",&his2DROCN);
  treeHis->SetBranchAddress("his2DROCSum.",&his2DROCSum);
  treeHis->SetBranchAddress("his2DRZROCN.",&his2DRZROCN);
  treeHis->SetBranchAddress("his2DRZROCSum.",&his2DRZROCSum);
  treeHis->GetEntry(event);
  //
  // Adjust the format for the AliTPCSpacecharge3D
  //
  TH3D *his2DROCNNew = new TH3D();
  TH3D *his2DROCSumNew = new TH3D();
  
  //
  TH2D *his2DRZROCNNorm = new TH2D(*his2DRZROCN);
  TH2D *his2DRZROCSumNorm = new TH2D(*his2DRZROCSum);

  for (Int_t i=1;i<his2DROCSum->GetXaxis()->GetNbins();i++){
    for (Int_t j=1;j<his2DROCSum->GetYaxis()->GetNbins();j++){
      if (his2DROCSum->GetXaxis()->GetBinCenter(j) < TMath::TwoPi()){
his2DROCSumNew->SetBinContent(i,j,0.499,his2DROCSum->GetBinContent(j,i));
      } else {
his2DROCSumNew->SetBinContent(i,j,-0.499,his2DROCSum->GetBinContent(j,i));
      }
    }
  }
  //
  // make normalization to C/cm**3 for RPHi
  //
  TH3D *his2DROCNNorm = new TH3D(*his2DROCNNew);
  TH3D *his2DROCSumNorm = new TH3D(*his2DROCSumNew);
  Double_t ePerADC = 500.;
  Double_t dV = 0;
  for (Int_t k=1;k<his2DROCSumNew->GetZaxis()->GetNbins();k++){
    for (Int_t i=1;i<his2DROCSumNew->GetXaxis()->GetNbins();i++){
      for (Int_t j=1;j<his2DROCSumNew->GetYaxis()->GetNbins();j++){
        dV = his2DROCSumNew->GetXaxis()->GetBinCenter(i)*1e-2 * his2DROCSumNew->GetYaxis()->GetBinWidth(j) * 2.5;
his2DROCSumNorm->SetBinContent(i,j,k,his2DROCSumNew->GetBinContent(i,j,k)*ePerADC*TMath::Qe()/dV); // Charge in C/m^3
      }
    }
  }
  //
  //  make normalization to C/cm**3 for RZ
  //
  for (Int_t i=1;i<his2DRZROCSum->GetXaxis()->GetNbins();i++){
    for (Int_t j=1;j<his2DRZROCSum->GetYaxis()->GetNbins();j++){
      dV = his2DRZROCSum->GetXaxis()->GetBinCenter(i)*1e-2 * his2DRZROCSum->GetYaxis()->GetBinWidth(j)*1e-2 * TMath::TwoPi()*his2DRZROCSum->GetXaxis()->GetBinCenter(i)*1e-2;
his2DRZROCSumNorm->SetBinContent(i,j,his2DRZROCSum->GetBinContent(i,j)*ePerADC*TMath::Qe()/dV); // Charge in C/m^3
    }
  }
  TFile *SCFile = new TFile(outputfile,"recreate");
  his2DROCSumNorm->Write("SpaceChargeInRPhi");
  his2DRZROCSumNorm->Write("SpaceChargeInRZ");
  SCFile->Close();

}


/*

  TFile f("histo.root")
  TH3* his0 = (TH3*)f.Get("Fluc0/his3DROCN")
  TH3* his1 = (TH3*)f.Get("Fluc1/his3DROCN")
  //
  //
  //
  AliTPCSpaceCharge3D *spaceCharge0 = new AliTPCSpaceCharge3D;
  spaceCharge0->SetInputSpaceCharge(his0,0.1);


  spaceCharge0->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  spaceCharge0->InitSpaceCharge3DDistortion();


  spaceCharge0->CreateHistoSCinZR(0.1,50,50)->Draw("surf1");

  spaceCharge0->CreateHistoDRPhiinZR(0,100,100)->Draw("colz");

  spaceCharge0->AddVisualCorrection(spaceCharge,1);


*/
