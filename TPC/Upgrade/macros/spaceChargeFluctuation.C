/*
  Macro to study the space charge fluctuations.
  3 functions using the ToyMC + analytical fomula to describe given MC results
 
  function to histogram space charge using the raw data ana anlyzing them
  To use given function - CPU conusming therefore batch farms used
  See  $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.sh macro to see example to run the code
  


  .x $HOME/NimStyle.C
  .x $HOME/rootlogon.C
  .L $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C+ 

 
  
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
#include "TStyle.h"
#include "AliTPCSpaceCharge3D.h"
#include "AliExternalTrackParam.h"
//
// constants
//
Double_t omegaTau=0.325;
//
// Function declaration
//
//   TOY MC
void spaceChargeFluctuationToyMC(Int_t nframes, Double_t interactionRate);
void spaceChargeFluctuationToyDraw();
void spaceChargeFluctuationToyDrawSummary();

//
// RAW data analysis
//
TH1 * GenerateMapRawIons(Int_t useGain,const char *fileName="raw.root", const char *outputName="histo.root", Int_t maxEvents=25);
void DoMerge();
void AnalyzeMaps1D();  // make nice plots
void MakeFluctuationStudy3D(Int_t nhistos, Int_t nevents, Int_t niter);
TH3D *  NormalizeHistoQ(TH3D * hisInput, Bool_t normEpsilon);
//
TH3D *  PermutationHistoZ(TH3D * hisInput, Double_t deltaZ);
TH3D *  PermutationHistoPhi(TH3D * hisInput, Double_t deltaPhi);
TH3D *  PermutationHistoLocalPhi(TH3D * hisInput, Int_t deltaPhi);
void MakeSpaceChargeFluctuationScan(Double_t scale, Int_t nfiles, Int_t sign);
void DrawFluctuationdeltaZ(Int_t stat=0, Double_t norm=10000);
void DrawFluctuationSector(Int_t stat=0, Double_t norm=10000);

void spaceChargeFluctuation(Int_t mode=0, Float_t arg0=0, Float_t arg1=0, Float_t arg2=0){
  //
  // function called from the shell script
  //
  gRandom->SetSeed(0);
  if (mode==0) GenerateMapRawIons(arg0);  
  if (mode==1) DoMerge();  
  if (mode==2) spaceChargeFluctuationToyMC(arg0,arg1);
  if (mode==3) MakeFluctuationStudy3D(10000, arg0, arg1);  
  if (mode==4) MakeSpaceChargeFluctuationScan(arg0,arg1,arg2); // param: scale, nfiles, sign Bz
  if (mode==5) {
    DrawFluctuationdeltaZ(arg0,arg1);
    DrawFluctuationSector(arg0,arg1);
  }
}


Double_t RndmdNchdY(Double_t s){
  //
  // dNch/deta - last 2 points inventeted (to find it somewhere ?)
  // 
  //  http://arxiv.org/pdf/1012.1657v2.pdf - table 1.  ALICE PbPb
  //  Scaled according s^0.15
  //  http://arxiv.org/pdf/1210.3615v2.pdf
  //  This we can cite. 
  //  Usage example::
  /*
    TH1F his550("his550","his550",1000,0,3000)
    for (Int_t i=0; i<300000; i++) his550->Fill(RndmdNchdY(5.5));
    his550->Draw();    
    TF1 f1("f1","[0]*x^(-(0.00001+abs([1])))",1,2000)
    f1.SetParameters(1,-1)
    his550->Fit("f1","","",10,3000);
    TH1F his276("his276","his276",1000,0,3000)
    for (Int_t i=0; i<300000; i++) his276->Fill(RndmdNchdY(2.76));
    his276->Draw();    

  */
  static TSpline3 * spline276=0;
  const Double_t sref=2.76; // reference s

  if (!spline276){
    // Refence multiplicities for 2.76 TeV
    // multplicity from archive except of the last  point was set to 0
    //
    const Double_t mult[20]={1601,  1294,   966,  649,   426,  261,  149,  76, 35,      0.001};
    const Double_t cent[20]={2.5,   7.5,    15,   25,    35,   45,   55,   65, 75,   100.};   
    TGraphErrors * gr = new TGraphErrors(10,cent,mult);
    spline276 = new TSpline3("spline276",gr);
  }
  Double_t norm = TMath::Power((s*s)/(sref*sref),0.15);
  spline276->Eval(gRandom->Rndm()*100.);
  return  spline276->Eval(gRandom->Rndm()*100.)*norm;
}





void pileUpToyMC(Int_t nframes){
  //
  //
  //
  /*
    Int)t nframes=1000;
   */
  TTreeSRedirector *pcstream = new TTreeSRedirector("pileup.root","recreate");
  Double_t central = 2350;
  Double_t pmean=5;
  TVectorD vectorT(nframes);
  //
  for (Int_t irate=1; irate<10; irate++){
    printf("rate\t%d\n",irate);
    for (Int_t iframe=0; iframe<nframes; iframe++){
      if (iframe%100000==0)printf("iframe=%d\n",iframe);
      Int_t ntracksAll=0;
      Int_t nevents=gRandom->Poisson(irate);
      Int_t ntracks=0; // to be taken from the MB primary distribution    
      Bool_t hasCentral=0;
      for (Int_t ievent=0; ievent<nevents; ievent++){
	ntracks=RndmdNchdY(5.5);
	ntracksAll+=ntracks; 
	if (ntracks>central) hasCentral = kTRUE;
      }    
      (*pcstream)<<"pileupFrame"<<
	"rate="<<irate<<
	"nevents="<<nevents<<
	"ntracks="<<ntracks<<
	"ntracksAll="<<ntracksAll<<
	"hasCentral"<<hasCentral<<
	"\n";
      vectorT[iframe]=ntracksAll;
    }
    Double_t mean   = TMath::Mean(nframes, vectorT.GetMatrixArray());
    Double_t rms    = TMath::RMS(nframes, vectorT.GetMatrixArray());
    Double_t median = TMath::Median(nframes, vectorT.GetMatrixArray());
    Double_t ord90  = TMath::KOrdStat(nframes,vectorT.GetMatrixArray() , Int_t(nframes*0.90));
    Double_t ord95  = TMath::KOrdStat(nframes,vectorT.GetMatrixArray() , Int_t(nframes*0.95));
    Double_t ord99  = TMath::KOrdStat(nframes,vectorT.GetMatrixArray() , Int_t(nframes*0.99));
    Double_t ord999  = TMath::KOrdStat(nframes,vectorT.GetMatrixArray() , Int_t(nframes*0.999));
    Double_t ord9999  = TMath::KOrdStat(nframes,vectorT.GetMatrixArray() , Int_t(nframes*0.9999));
    (*pcstream)<<"pileup"<<
      "rate="<<irate<<
      "mean="<<mean<<
      "rms="<<rms<<
      "median="<<median<<
      "ord90="<<ord90<<
      "ord95="<<ord95<<
      "ord99="<<ord99<<
      "ord999="<<ord999<<
      "ord9999="<<ord9999<<
      "\n";
  }
  delete pcstream;  
  // Draw
  pcstream = new TTreeSRedirector("pileup.root","update");
  TTree * treeStat = (TTree*)(pcstream->GetFile()->Get("pileup"));
  TTree * treeFrame = (TTree*)(pcstream->GetFile()->Get("pileupFrame"));
  Int_t  mentries =  treeStat->Draw("ord999","1","goff");
  Double_t maximum = TMath::MaxElement(mentries, treeStat->GetV1());
  const char * names[6]={"mean","median","ord90","ord95","ord99","ord999"};  
  const char * titles[6]={"Mean","Median","90 %","95 %","99 %","99.9 %"};  
  const Int_t mcolors[6]={1,2,3,4,6,7};  
  //
  //
  TF1 * f1 = new TF1("f1","[0]*x+[1]*sqrt(x)");
  Double_t par0=0;
  //
  TCanvas * canvasMult = new TCanvas("canvasCumul","canvasCumul");
  canvasMult->SetLeftMargin(0.13);
  TLegend * legend= new TLegend(0.14,0.6,0.45,0.89, "Effective dN_{ch}/d#eta");
  TGraphErrors *graphs[6]={0};  
  for (Int_t igr=0; igr<6; igr++){
    graphs[igr] = TStatToolkit::MakeGraphErrors(treeStat,Form("%s:rate",names[igr]),"1",21+(igr%5),mcolors[igr],0);
    graphs[igr]->SetMinimum(0);
    graphs[igr]->GetYaxis()->SetTitleOffset(1.3);
    graphs[igr]->SetMaximum(maximum*1.1);
    graphs[igr]->GetXaxis()->SetTitle("<N_{ev}>");
    graphs[igr]->GetYaxis()->SetTitle("dN_{ch}/d#eta");
    TF1 * f2 = new TF1("f2","[0]*x+[1]*sqrt(x)");
    f2->SetLineColor(mcolors[igr]);
    f2->SetLineWidth(0.5);
    if (igr>0) f2->FixParameter(0,par0);
    graphs[igr]->Fit(f2,"","");
    if (igr==0) par0=f2->GetParameter(0);
    if (igr==0) graphs[igr]->Draw("ap");
    graphs[igr]->Draw("p");
    legend->AddEntry(graphs[igr], titles[igr],"p");
  }
  legend->SetBorderSize(0);
  legend->Draw();

  canvasMult->SaveAs("effectiveMult.pdf");
  canvasMult->SaveAs("effectiveMult.png");
  gStyle->SetOptStat(0);
  TH2F * hisMult = new TH2F("ntracksNevent","ntracksnevents",9,1,10,100,0,2*maximum);
  {
    treeFrame->Draw("ntracksAll:rate>>ntracksNevent","","colz");
    hisMult->GetXaxis()->SetTitle("<N_{ev}>");
    hisMult->GetYaxis()->SetTitle("dN_{ch}/d#eta");
    hisMult->GetYaxis()->SetTitleOffset(1.3);
    hisMult->Draw("colz");
  }
  canvasMult->SaveAs("effectiveMultColz.pdf");
  canvasMult->SaveAs("effectiveMultColz.png");
  //
  //
  //
  TH2F * hisMult5 = new TH2F("ntracksNevent5","ntracksnEvents5",9,1,10,100,0,maximum);
  {
    treeFrame->Draw("ntracksAll:nevents>>ntracksNevent5","abs(rate-5)<0.5","colz");
    hisMult5->GetXaxis()->SetTitle("N_{ev}");
    hisMult5->GetYaxis()->SetTitle("dN_{ch}/d#eta");
    hisMult5->GetYaxis()->SetTitleOffset(1.3);
    hisMult5->Draw("colz");
  }
  canvasMult->SaveAs("effectiveMultF5.pdf");
  canvasMult->SaveAs("effectiveMultF5.png");

  {    
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    gStyle->SetOptTitle(1);
    TCanvas * canvasMultH = new TCanvas("canvasCumulH","canvasCumulH",700,700);
    canvasMultH->Divide(1,2);
    canvasMultH->cd(1);
    TH1F his550("his550","his550",1000,0,3000);
    TH1F his276("his276","his276",1000,0,3000);
    for (Int_t i=0; i<300000; i++) his550.Fill(RndmdNchdY(5.5));
    for (Int_t i=0; i<300000; i++) his276.Fill(RndmdNchdY(2.76));     
    TF1 f1("f1","[0]*x^(-(0.00001+abs([1])))",1,2000);
    f1.SetParameters(1,-1);
    his550.GetXaxis()->SetTitle("dN_{ch}/d#eta");
    his276.GetXaxis()->SetTitle("dN_{ch}/d#eta");
    his550.Fit("f1","","",10,3000);
    his276.Fit("f1","","",10,3000); 
    canvasMultH->cd(1)->SetLogx(1);
    canvasMultH->cd(1)->SetLogy(1);
    his550.Draw();    
    canvasMultH->cd(2)->SetLogx(1);
    canvasMultH->cd(2)->SetLogy(1);
    his276.Draw("");    
    canvasMultH->SaveAs("dNchdEta.pdf");
  }
  delete pcstream;
}

void spaceChargeFluctuationToyMC(Int_t nframes, Double_t interactionRate){
  //
  // Toy MC to generate space charge fluctuation, to estimate the fluctuation of the integral space charge in part of the
  // TPC
  // Parameters:
  //    nframes - number of frames to simulate 
  // 1. Make a toy simulation part for given setup
  // 2. Make a summary plots for given setups - see function spaceChargeFluctuationToyMCDraw()
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("spaceChargeFluctuation.root","recreate");
  Double_t driftTime=0.1;              
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
      "rate="<<interactionRate<<                  // interaction rate
      "eventMean="<<eventMean<<                   // mean number of events per frame
      "trackMean="<<trackMean<<                   // assumed mean of the tracks per event
      //       
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
  spaceChargeFluctuationToyDraw();
}


void spaceChargeFluctuationToyDraw(){
  //
  // Toy MC to simulate the space charge integral fluctuation
  // Draw function for given setup
  // for MC generation part see : void spaceChargeFluctuationToyMC
  TTreeSRedirector *pcstream = new TTreeSRedirector("spaceChargeFluctuation.root","update");
  TFile * f = pcstream->GetFile();
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
  nentries = treeStat->Draw("eventMean","","");
  Double_t meanEvents =TMath::Mean(nentries,treeStat->GetV1());  
  treeStat->SetMarkerStyle(21);
  treeStat->SetMarkerSize(0.4);
  //
  const Int_t kColors[6]={1,2,3,4,6,7};
  const Int_t kStyle[6]={20,21,24,25,24,25};
  const char  * htitles[6]={"Events","Tracks","Tracks #phi region (1/180)","Q #phi region (1/180)", "Tracks #phi region (1/36)","Q #phi region (1/36)"}; 
  const char  * hnames[6]={"Events","Tracks","TracksPhi180","QPhi180", "TracksPhi36","QPhi36"}; 

  TH1* hisFluc[6]={0};
  TH1* hisPull[6]={0};
  TVectorD *vecFitFluc[6]={0};
  TVectorD *vecFitFlucPull[6]={0};
  //
  // histograms
  //
  treeStat->Draw("nevents/eventMean>>hisEv(100,0.85,1.15)","");
  hisFluc[0]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("ntracksAll/(eventMean*tracksMean)>>hisTrackAll(100,0.85,1.1)","","same");
  hisFluc[1]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("vecTracksPhi180.fElements/(eventMean*tracksMean/180)>>hisTrackSector(100,0.85,1.1)","1/180","same");
  hisFluc[2]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("vecEPhi180.fElements/(eventMean*tracksMean*meandE/180)>>hisdESector(100,0.85,1.1)","1/180","same");
  hisFluc[3]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("vecTracksPhi36.fElements/(eventMean*tracksMean/36)>>hisTrackSector36(100,0.85,1.1)","1/36","same");
  hisFluc[4]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("vecEPhi36.fElements/(eventMean*tracksMean*meandE/36)>>hisdESector36(100,0.85,1.1)","1/36","same");
  hisFluc[5]=(TH1*)treeStat->GetHistogram()->Clone();
  //
  // pulls
  //
  treeStat->Draw("((nevents/eventMean)-1)/sqrt(1/eventMean)>>pullEvent(100,-6,6)","","err");  //tracks All pull 
  hisPull[0]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("(ntracksAll/(eventMean*tracksMean)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean)>>pullTrackAll(100,-6,6)","","err");  //tracks All pull 
  hisPull[1]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("(vecTracksPhi180.fElements/(eventMean*tracksMean/180.)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+180/(tracksMean*eventMean))>>pullTrack180(100,-6,6)","1/180","errsame");  //tracks spread
  hisPull[2]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("(vecEPhi180.fElements/(eventMean*tracksMean*meandE/180)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+180/(tracksMean*eventMean)+180*(rmsdE/meandE)**2/(eventMean*tracksMean))>>hisPulldE180(100,-6,6)","1/180","errsame"); //dE spread
  hisPull[3]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("(vecTracksPhi36.fElements/(eventMean*tracksMean/36.)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+36/(tracksMean*eventMean))>>pullTrack36(100,-6,6)","1/36","errsame");  //tracks spread
  hisPull[4]=(TH1*)treeStat->GetHistogram()->Clone();
  treeStat->Draw("(vecEPhi36.fElements/(eventMean*tracksMean*meandE/36)-1)/sqrt(1/eventMean+(tracksRMS/tracksMean)**2/eventMean+36/(tracksMean*eventMean)+36*(rmsdE/meandE)**2/(eventMean*tracksMean))>>hisPulldE36(100,-6,6)","1/36","errsame"); //dE spread
  hisPull[5]=(TH1*)treeStat->GetHistogram()->Clone();
  //
  
  for (Int_t ihis=0; ihis<6; ihis++) {
    vecFitFluc[ihis] = new TVectorD(3);
    vecFitFlucPull[ihis] =  new TVectorD(3);
    TF1 * fg = new TF1(Form("fg%d",ihis),"gaus");    
    fg->SetLineWidth(0.5); 
    fg->SetLineColor(kColors[ihis]); 
    hisFluc[ihis]->Fit(fg,"","");
    fg->GetParameters( vecFitFluc[ihis]->GetMatrixArray());
    hisPull[ihis]->Fit(fg,"","");
    fg->GetParameters( vecFitFlucPull[ihis]->GetMatrixArray());
    hisFluc[ihis]->SetName(Form("Fluctuation%s",hnames[ihis]));
    hisFluc[ihis]->SetTitle(Form("Fluctuation%s",htitles[ihis]));
    hisPull[ihis]->SetName(Form("Pull%s",hnames[ihis]));
    hisPull[ihis]->SetTitle(Form("Pull%s",htitles[ihis]));
  } 
  
  gStyle->SetOptStat(0);
  TCanvas * canvasQFluc= new TCanvas("SpaceChargeFluc","SpaceChargeFluc",600,700);
  canvasQFluc->Divide(1,2);
  canvasQFluc->cd(1);
  TLegend * legendFluc = new TLegend(0.11,0.55,0.45,0.89,"Relative fluctuation");
  TLegend * legendPull = new TLegend(0.11,0.55,0.45,0.89,"Fluctuation pulls");
  for (Int_t ihis=0; ihis<6; ihis++){
    hisFluc[ihis]->SetMarkerStyle(kStyle[ihis]);
    hisFluc[ihis]->SetMarkerColor(kColors[ihis]);
    hisFluc[ihis]->SetMarkerSize(0.8);
    if (ihis==0) hisFluc[ihis]->Draw("err"); 
    hisFluc[ihis]->Draw("errsame");    
    legendFluc->AddEntry(hisFluc[ihis],htitles[ihis]);
  }
  legendFluc->Draw();

  canvasQFluc->cd(2);
  for (Int_t ihis=0; ihis<6; ihis++){
    hisPull[ihis]->SetMarkerStyle(kStyle[ihis]);
    hisPull[ihis]->SetMarkerColor(kColors[ihis]);
    hisPull[ihis]->SetMarkerSize(0.8);
    if (ihis==0) hisPull[ihis]->Draw("err"); 
    hisPull[ihis]->Draw("errsame");    
    legendPull->AddEntry(hisPull[ihis],htitles[ihis]);
  }
  legendPull->Draw();
  //
  for (Int_t ihis=0; ihis<6; ihis++){
    hisFluc[ihis]->Write();
    hisPull[ihis]->Write();
  }
  (*pcstream)<<"summary"<<                             // summary information for given setup
    "meanEvents="<<meanEvents<<                        // mean number of events in the frame
    "meandE="<<meandE<<                                // mean "energy loss" of track
    "rmsdE="<<rmsdE<<                                  // rms 
    "meanT="<<meanT<<                                  // mean number of tracks per MB event
    "rmsT="<<rmsT<<                                    // rms of onumber of tracks
    //                                                 // fit of the relative fluctuation 
    "vflucE.="<<vecFitFluc[0]<<                        //         in events
    "vflucEP.="<<vecFitFlucPull[0]<<                   //         in events pull
    "vflucTr.="<<vecFitFluc[1]<<                       //         in tracks 
    "vflucTrP.="<<vecFitFlucPull[1]<<
    //
    "vflucTr180.="<<vecFitFluc[2]<<
    "vflucTr180P.="<<vecFitFlucPull[2]<<
    "vflucE180.="<<vecFitFluc[3]<<
    "vflucE180P.="<<vecFitFlucPull[3]<<
    //
    "vflucTr36.="<<vecFitFluc[4]<<
    "vflucTr36P.="<<vecFitFlucPull[4]<<
    "vflucE36.="<<vecFitFluc[5]<<
    "vflucE36P.="<<vecFitFlucPull[5]<<
    "\n"; 
  canvasQFluc->SaveAs("CanvasFluctuation.pdf");
  canvasQFluc->SaveAs("CanvasFluctuation.png");
  delete pcstream;

}

void spaceChargeFluctuationToyDrawSummary(){
  //
  // make a summary information plots using several runs with differnt mean IR setting
  // Input:
  //   space.list - list of root files produced by spaceChargeFluctuationToyDraw   
  // Output:
  //   canvas saved in current directory
  //
  TChain * chain = AliXRDPROOFtoolkit::MakeChain("space.list","summary",0,100);
  chain->SetMarkerStyle(21);
  const Int_t kColors[6]={1,2,3,4,6,7};
  const Int_t kStyle[6]={20,21,24,25,24,25};
  const char  * htitles[6]={"Events","Tracks","Tracks #phi region (1/180)","Q #phi region (1/180)", "Tracks #phi region (1/36)","Q #phi region (1/36)"}; 
  //  const char  * hnames[6]={"Events","Tracks","TracksPhi180","QPhi180", "TracksPhi36","QPhi36"}; 
  //
  Double_t meanT,rmsT=0; 
  Double_t meandE,rmsdE=0;
  Int_t entries = chain->Draw("meanT:rmsT:meandE:rmsdE","1","goff");
  meanT =TMath::Mean(entries, chain->GetV1());
  rmsT =TMath::Mean(entries, chain->GetV2());
  meandE =TMath::Mean(entries, chain->GetV3());
  rmsdE =TMath::Mean(entries, chain->GetV4());
  //
  //
  //
  TGraphErrors * graphs[6]={0};
  TF1 * functions[6]={0};

  graphs[5]=TStatToolkit::MakeGraphErrors(chain,"vflucE36.fElements[2]:meanEvents:0.025*vflucE36.fElements[2]","1",kStyle[5],kColors[5],1);
  graphs[4]=TStatToolkit::MakeGraphErrors(chain,"vflucTr36.fElements[2]:meanEvents:0.025*vflucTr36.fElements[2]","1",kStyle[4],kColors[4],1);
  graphs[3]=TStatToolkit::MakeGraphErrors(chain,"vflucE180.fElements[2]:meanEvents:0.025*vflucE180.fElements[2]","1",kStyle[3],kColors[3],1);
  graphs[2]=TStatToolkit::MakeGraphErrors(chain,"vflucTr180.fElements[2]:meanEvents:0.025*vflucTr180.fElements[2]","1",kStyle[2],kColors[2],1);
  graphs[1]=TStatToolkit::MakeGraphErrors(chain,"vflucTr.fElements[2]:meanEvents:0.025*vflucTr.fElements[2]","1",kStyle[1],kColors[1],1);
  graphs[0]=TStatToolkit::MakeGraphErrors(chain,"vflucE.fElements[2]:meanEvents:0.025*vflucE.fElements[2]","1",kStyle[0],kColors[0],1);

  functions[5]=new TF1("fe","sqrt(1+[0]**2+[1]/[2]+[1]*[3]**2/[2])/sqrt(x)",2000,200000);  
  functions[5]->SetParameters(rmsT/meanT,36.,meanT,rmsdE/meandE);
  functions[4]=new TF1("fe","sqrt(1+[0]**2+[1]/[2])/sqrt(x)",2000,200000);  
  functions[4]->SetParameters(rmsT/meanT,36.,meanT,0);
  functions[3]=new TF1("fe","sqrt(1+[0]**2+[1]/[2]+[1]*[3]**2/[2])/sqrt(x)",2000,200000);  
  functions[3]->SetParameters(rmsT/meanT,180.,meanT,rmsdE/meandE);
  functions[2]=new TF1("fe","sqrt(1+[0]**2+[1]/[2])/sqrt(x)",2000,200000);  
  functions[2]->SetParameters(rmsT/meanT,180.,meanT,0);
  functions[1]=new TF1("fe","sqrt(1+[0]**2)/sqrt(x)",2000,200000);  
  functions[1]->SetParameters(rmsT/meanT,0);
  functions[0]=new TF1("fe","sqrt(1)/sqrt(x)",2000,200000);  
  
  
  TCanvas *canvasF= new TCanvas("fluc","fluc",600,500);  
  //  TLegend *legend = new TLegend(0.5,0.65,0.89,0.89,"Relative fluctuation #sigma=#sqrt{1+#frac{#sigma_{T}^{2}}{#mu_{T}^{2}}}");
  TLegend *legendF = new TLegend(0.45,0.5,0.89,0.89,"Relative fluctuation of charge");
  for (Int_t ihis=0; ihis<4; ihis++){
    graphs[ihis]->SetMinimum(0.00);
    graphs[ihis]->SetMaximum(0.05);
    if (ihis==0) graphs[ihis]->Draw("ap");
    graphs[ihis]->GetXaxis()->SetTitle("events");
    graphs[ihis]->GetXaxis()->SetNdivisions(507);
    graphs[ihis]->GetYaxis()->SetTitle("#frac{#sigma}{#mu}");
    graphs[ihis]->Draw("p");    
    legendF->AddEntry(graphs[ihis],htitles[ihis],"p");
    if (functions[ihis]){
      functions[ihis]->SetLineColor(kColors[ihis]);
      functions[ihis]->SetLineWidth(0.5);
      functions[ihis]->Draw("same");
    }
  }
  legendF->Draw();
  canvasF->SaveAs("spaceChargeFlucScan.pdf");
  canvasF->SaveAs("spaceChargeFlucScan.png");

  TCanvas *canvasF36= new TCanvas("fluc36","fluc36",600,500);  
  //  TLegend *legend = new TLegend(0.5,0.65,0.89,0.89,"Relative fluctuation #sigma=#sqrt{1+#frac{#sigma_{T}^{2}}{#mu_{T}^{2}}}");
  TLegend *legendF36 = new TLegend(0.45,0.5,0.89,0.89,"Relative fluctuation of charge");
  for (Int_t ihis=0; ihis<6; ihis++){
    if (ihis==2 || ihis==3) continue;
    graphs[ihis]->SetMinimum(0.00);
    graphs[ihis]->SetMaximum(0.05);
    if (ihis==0) graphs[ihis]->Draw("ap");
    graphs[ihis]->GetXaxis()->SetTitle("events");
    graphs[ihis]->GetXaxis()->SetNdivisions(507);
    graphs[ihis]->GetYaxis()->SetTitle("#frac{#sigma}{#mu}");
    graphs[ihis]->Draw("p");    
    legendF36->AddEntry(graphs[ihis],htitles[ihis],"p");
    if (functions[ihis]){
      functions[ihis]->SetLineColor(kColors[ihis]);
      functions[ihis]->SetLineWidth(0.5);
      functions[ihis]->Draw("same");
    }
  }
  legendF36->Draw();
  canvasF36->SaveAs("spaceChargeFlucScan36.pdf");
  canvasF36->SaveAs("spaceChargeFlucScan36.png");


}



void FitMultiplicity(const char * fname="mult_dist_pbpb.root"){
  //
  // Fit multiplicity distribution using as a power law in limited range
  //  const char * fname="mult_dist_pbpb.root"
  TFile *fmult=TFile::Open(fname);
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
  gRandom->SetSeed(0);  //set initial seed to be independent for different jobs

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
  TH1D * hisQ1D[3]={0};
  TH1D * hisQ1DROC[3]={0};
  TH2D * hisQ2DRPhi[3]={0};                
  TH2D * hisQ2DRZ[3]={0};                
  TH2D * hisQ2DRPhiROC[3]={0};
  TH2D * hisQ2DRZROC[3]={0};                
  TH3D * hisQ3D[3]={0};                // 3D maps space charge from drift volume  
  TH3D * hisQ3DROC[3]={0};             // 3D maps space charge from ROC
  
  Int_t nbinsRow=param->GetNRowLow()+param->GetNRowUp();
  Double_t *xbins = new Double_t[nbinsRow+1];
  xbins[0]=param->GetPadRowRadiiLow(0)-1;   //underflow bin
  for (Int_t ibin=0; ibin<param->GetNRowLow();ibin++) xbins[1+ibin]=param->GetPadRowRadiiLow(ibin);
  for (Int_t ibin=0; ibin<param->GetNRowUp();ibin++)  xbins[1+ibin+param->GetNRowLow()]=param->GetPadRowRadiiUp(ibin);
  //
  for (Int_t ith=0; ith<3; ith++){
    char chname[100];
    // 1D
    snprintf(chname,100,"hisQ1D_Th%d",2*ith+2);
    hisQ1D[ith] = new TH1D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ1DROC_Th%d",2*ith+2);
    hisQ1DROC[ith] = new TH1D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    // 3D
    snprintf(chname,100,"hisQ3D_Th%d",2*ith+2);
    hisQ3D[ith] = new TH3D(chname, chname,360, 0,TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    snprintf(chname,100,"hisQ3DROC_Th%d",2*ith+2);
    hisQ3DROC[ith] = new TH3D(chname, chname,360, 0,TMath::TwoPi(),param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) ,125,-250,250);
    // 2D
    snprintf(chname,100,"hisQ2DRPhi_Th%d",2*ith+2);
    hisQ2DRPhi[ith] = new TH2D(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ2DRZ_Th%d",2*ith+2);
    hisQ2DRZ[ith] = new TH2D(chname,chname, param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36),  125,-250,250);
    //
    snprintf(chname,100,"hisQ2DRPhiROC_Th%d",2*ith+2);
    hisQ2DRPhiROC[ith] = new TH2D(chname,chname,180, 0,2*TMath::TwoPi(), param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36) );
    snprintf(chname,100,"hisQ2DRZROC_Th%d",2*ith+2);
    hisQ2DRZROC[ith] = new TH2D(chname,chname,param->GetNRow(0)+param->GetNRow(36) ,0, param->GetNRow(0)+param->GetNRow(36), 125,-250,250);
    //
    hisQ1D[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ1DROC[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ3D[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ3DROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    //
    hisQ2DRPhi[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ2DRZ[ith]->GetXaxis()->Set(nbinsRow,xbins);
    hisQ2DRPhiROC[ith]->GetYaxis()->Set(nbinsRow,xbins);
    hisQ2DRZROC[ith]->GetXaxis()->Set(nbinsRow,xbins);
    //
    hisQ1D[ith]->SetDirectory(0);
    hisQ1DROC[ith]->SetDirectory(0);
    hisQ3D[ith]->SetDirectory(0);
    hisQ3DROC[ith]->SetDirectory(0);
    //
    hisQ2DRPhi[ith]->SetDirectory(0);
    hisQ2DRZ[ith]->SetDirectory(0);
    hisQ2DRZROC[ith]->SetDirectory(0);
    hisQ2DRPhiROC[ith]->SetDirectory(0);
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
	hisQ2DRPhi[ith]->Write(Form("His2DRPhiDrift_%d",ith));
	hisQ2DRZ[ith]->Write(Form("His2DRZDrift_%d",ith));
	hisQ3D[ith]->Write(Form("His3DDrift_%d",ith));
	hisQ1DROC[ith]->Write(Form("His1DROC_%d",ith));
	hisQ2DRPhiROC[ith]->Write(Form("His2DRPhiROC_%d",ith));
	hisQ2DRZROC[ith]->Write(Form("His2DRZROC_%d",ith));
	hisQ3DROC[ith]->Write(Form("His3DROC_%d",ith));
	(*pcstream)<<"histo"<<
	  "events="<<evtnr<<
	  "useGain="<<useGainMap<<
	  Form("hist1D_%d.=",ith*2+2)<<hisQ1D[ith]<<
	  Form("hist2DRPhi_%d.=",ith*2+2)<<hisQ2DRPhi[ith]<<
	  Form("hist2DRZ_%d.=",ith*2+2)<<hisQ2DRZ[ith]<<
	  Form("hist3D_%d.=",ith*2+2)<<hisQ3D[ith]<<
	  Form("hist1DROC_%d.=",ith*2+2)<<hisQ1DROC[ith]<<
	  Form("hist2DRPhiROC_%d.=",ith*2+2)<<hisQ2DRPhiROC[ith]<<
	  Form("hist2DRZROC_%d.=",ith*2+2)<<hisQ2DRZROC[ith]<<
	  Form("hist3DROC_%d.=",ith*2+2)<<hisQ3DROC[ith];	  
      }
      (*pcstream)<<"histo"<<"\n";
      for (Int_t ith=0; ith<3; ith++){	
	hisQ1D[ith]->Reset();
	hisQ2DRPhi[ith]->Reset();
	hisQ2DRZ[ith]->Reset();
	hisQ3D[ith]->Reset();
	hisQ1DROC[ith]->Reset();
	hisQ2DRPhiROC[ith]->Reset();
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
	Int_t    pad    = input.GetPad();
	Int_t    nPads   = param->GetNPads(sector,row);
	Double_t localX  = param->GetPadRowRadii(sector,row); 
	Double_t localY  = (pad-nPads/2)*param->GetPadPitchWidth(sector);
	Double_t localPhi= TMath::ATan2(localY,localX);
	Double_t phi     = TMath::Pi()*((sector%18)+0.5)/9+localPhi;
	Double_t padLength=param->GetPadPitchLength(sector,row);
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
		Double_t shiftPhi = ((sector%36)<18) ? 0: TMath::TwoPi();
                if (TMath::Abs(zIonDrift)<param->GetZLength()){
		  if ((sector%36)>=18) zIonDrift*=-1;   // c side has opposite sign
		  if (sector%36<18) hisQ1D[ith]->Fill(localX, signal/padLength);
		  hisQ2DRPhi[ith]->Fill(phi+shiftPhi,localX, signal/padLength);
		  hisQ2DRZ[ith]->Fill(localX, zIonDrift, signal/padLength);
		  hisQ3D[ith]->Fill(phi,localX,zIonDrift,signal/padLength);
		}
		//
		Double_t zIonROC = ((sector%36)<18)? shiftZ: -shiftZ;  // z position of the "ion disc" -  A side C side opposite sign
		if (sector%36<18) hisQ1DROC[ith]->Fill(localX, signal/padLength);
		hisQ2DRPhiROC[ith]->Fill(phi+shiftPhi,localX, signal/padLength);
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


void DoMerge(){
  //
  // Merge results to the tree
  //
  TFile *  fhisto = new TFile("histo.root","recreate");
  TTree * tree = 0;
  TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom("histo.list","histo",0,100,1);
  chain->SetBranchStatus("hist3DROC_6*",kFALSE);
  chain->SetBranchStatus("hist3DROC_4*",kFALSE);
  tree = chain->CopyTree("1");
  tree->Write("histo");
  delete fhisto;
}




void AnalyzeMaps1D(){
  //
  // Analyze space charge maps stored as s hitograms in trees
  //
  TFile *  fhisto = new TFile("histo.root");
  TTree * tree = (TTree*)fhisto->Get("histo");
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
  canvasR->SaveAs("spaceCharge1d.png");
  canvasR->SaveAs("spaceCharge1d.eps");
  //
  //
  //
}
void MakeFluctuationStudy3D(Int_t nhistos, Int_t nevents, Int_t niter){
  //
  //
  // 
  // 
  // Input:
  //   nhistos - maximal number of histograms to be used for sum 
  //   nevents - number of events to make a fluctuation studies
  //   niter   - number of itterations
  // Algortihm: 
  // 1. Make a summary integral   3D/2D/1D maps
  // 2. Create several maps with niter events  - Poisson flucturation in n
  // 3. Store results 3D maps in the tree (and also as histogram)  current and mean
  //   

  TFile *  fhisto = TFile::Open("histo.root");
  TTree * tree = (TTree*)fhisto->Get("histo");
  tree->SetCacheSize(10000000000);

  TTreeSRedirector * pcstream  = new TTreeSRedirector("fluctuation.root", "update");
  

  TH1D * his1DROC=0,    * his1DROCSum=0,  * his1DROCN=0;
  TH1D * his1DDrift=0,  * his1DDriftSum=0, * his1DDriftN=0 ;
  TH2D * his2DRPhiROC=0,    * his2DRPhiROCSum=0,  * his2DRPhiROCN=0;
  TH2D * his2DRZROC=0,    * his2DRZROCSum=0,  * his2DRZROCN=0;
  TH2D * his2DRPhiDrift=0,  * his2DRPhiDriftSum=0, * his2DRPhiDriftN=0;  
  TH2D * his2DRZDrift=0,  * his2DRZDriftSum=0, * his2DRZDriftN=0;  
  TH3D * his3DROC=0,    * his3DROCSum=0,  * his3DROCN=0;
  TH3D * his3DDrift=0,  * his3DDriftSum=0, * his3DDriftN=0;
  //
  if (nhistos<0 || nhistos> tree->GetEntries()) nhistos = tree->GetEntries();
  Int_t  eventsPerChunk=0;
  tree->SetBranchAddress("hist1D_2.",&his1DDrift);
  tree->SetBranchAddress("hist1DROC_2.",&his1DROC);
  tree->SetBranchAddress("hist2DRPhi_2.",&his2DRPhiDrift);
  tree->SetBranchAddress("hist2DRZ_2.",&his2DRZDrift);
  tree->SetBranchAddress("hist2DRPhiROC_2.",&his2DRPhiROC);
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
    if (his2DRPhiROCSum==0)     his2DRPhiROCSum=new TH2D(*his2DRPhiROC);
    if (his2DRZROCSum==0)     his2DRZROCSum=new TH2D(*his2DRZROC);
    if (his2DRPhiDriftSum==0)   his2DRPhiDriftSum=new TH2D(*his2DRPhiDrift);
    if (his2DRZDriftSum==0)   his2DRZDriftSum=new TH2D(*his2DRZDrift);
    if (his3DROCSum==0)     his3DROCSum=new TH3D(*his3DROC);
    if (his3DDriftSum==0)   his3DDriftSum=new TH3D(*his3DDrift);
    his1DROCSum->Add(his1DROC);
    his1DDriftSum->Add(his1DDrift);
    his2DRPhiROCSum->Add(his2DRPhiROC);
    his2DRZROCSum->Add(his2DRZROC);
    his2DRPhiDriftSum->Add(his2DRPhiDrift);
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
    Int_t nchunks=gRandom->Poisson(nevents)/eventsPerChunk;  // chunks with n typically 25 events
    for (Int_t i=0; i<nchunks; i++){
      tree->GetEntry(gRandom->Rndm()*nhistos);
      if (i%10==0) printf("%d\t%d\n",iter, i);
      if (his1DROCN==0)     his1DROCN=new TH1D(*his1DROC);
      if (his1DDriftN==0)   his1DDriftN=new TH1D(*his1DDrift);
      if (his2DRPhiROCN==0)     his2DRPhiROCN=new TH2D(*his2DRPhiROC);
      if (his2DRPhiDriftN==0)   his2DRPhiDriftN=new TH2D(*his2DRPhiDrift);
      if (his2DRZROCN==0)     his2DRZROCN=new TH2D(*his2DRZROC);
      if (his2DRZDriftN==0)   his2DRZDriftN=new TH2D(*his2DRZDrift);
      if (his3DROCN==0)     his3DROCN=new TH3D(*his3DROC);
      if (his3DDriftN==0)   his3DDriftN=new TH3D(*his3DDrift);
      his1DROCN->Add(his1DROC);
      his1DDriftN->Add(his1DDrift);
      his2DRPhiROCN->Add(his2DRPhiROC);
      his2DRZDriftN->Add(his2DRZDrift);
      his2DRZROCN->Add(his2DRZROC);
      his2DRPhiDriftN->Add(his2DRPhiDrift);
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
      "his2DRPhiROCN.="<<his2DRPhiROCN<<
      "his2DRPhiROCSum.="<<his2DRPhiROCSum<<
      "his2DRPhiDriftN.="<<his2DRPhiDriftN<<
      "his2DRPhiDriftSum.="<<his2DRPhiDriftSum<<
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
    his2DRPhiROCN->Write("his2DRPhiROCN");
    his2DRZROCN->Write("his2DRZROCN");
    //
    his2DRPhiROCSum->Write("his2DRPhiROCSum");        
    his2DRZROCSum->Write("his2DRZROCSum");
    //
    his2DRPhiDriftN->Write("his2DRPhiDriftN");
    his2DRZDriftN->Write("his2DRZDriftN");
    //
    his2DRPhiDriftSum->Write("his2DRPhiDriftSum");
    his2DRZDriftSum->Write("his2DRZDriftSum");
    //
    his3DROCN->Write("his3DROCN");
    his3DROCSum->Write("his3DROCSum");
    his3DDriftN->Write("his3DDriftN");
    his3DDriftSum->Write("his3DDriftSum");

    his1DROCN->Reset();
    his1DDriftN->Reset();
    his2DRPhiROCN->Reset();
    his2DRZDriftN->Reset();
    his2DRZROCN->Reset();
    his2DRPhiDriftN->Reset();
    his3DROCN->Reset();
    his3DDriftN->Reset();    
  }

  delete pcstream;
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


void MakeSpaceChargeFluctuationScan(Double_t scale, Int_t nfilesMerge, Int_t sign){
  //
  //
  // Input:
  //   scale           - scaling of the space charge
  //   nfilesMerge     - amount of chunks to merge
  //                   - =0  all chunks used
  //                     <0  subset form full statistic
  //                     >0  subset from the  limited (1000 mean) stistic
  // Make fluctuation scan:
  //   1.) Shift of z disk - to show which granularity in time needed
  //   2.) Shift in sector - to show influence of the gass gain and epsilon
  //   3.) Smearing in phi - to define phi granularity needed
  //   4.) Rebin z
  //   5.) Rebin phi
  // For given SC setups the distortion on the space point and track level characterezed
  //    SpaceChargeFluc%d_%d.root
  //    SpaceChargeTrackFluc%d%d.root
  //
  

  //
  // Some constant definition
  //
  Int_t nitteration=100;    // number of itteration in the lookup
  Int_t fullNorm  =10000;  // normalization  fro the full statistic

  //
  // Init magnetic field and OCDB
  //
  
  Double_t bsign= sign;
  if (bsign>1) bsign=-1;
  const Int_t nTracks=2000;
  const char *ocdb="local://$ALICE_ROOT/OCDB/";
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(0);   
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", bsign, bsign, AliMagF::k5kG));   
  //
  

  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("SpaceChargeFluc%d_%d.root",nfilesMerge,sign),"recreate");
  TTreeSRedirector *pcstreamTrack = new TTreeSRedirector(Form("SpaceChargeTrackFluc%d_%d.root",nfilesMerge,sign),"recreate");
  TH1D *his1DROCN=0, *his1DROCSum=0; 
  TH2D *his2DRPhiROCN=0, *his2DRPhiROCSum=0, *his2DRZROCN=0, *his2DRZROCSum=0;
  TH3D *his3DROCN=0, *his3DROCSum=0; 
  TH1D *his1DROCNC=0, *his1DROCSumC=0; 
  TH2D *his2DRPhiROCNC=0, *his2DRPhiROCSumC=0, *his2DRZROCNC=0, *his2DRZROCSumC=0;
  TH3D *his3DROCNC=0, *his3DROCSumC=0; 
  TH1 * histos[8]={his1DROCN, his1DROCSum, his2DRPhiROCN, his2DRPhiROCSum, his2DRZROCN, his2DRZROCSum,  his3DROCN, his3DROCSum};
  Int_t neventsAll=0, neventsAllC=0;
  Int_t neventsChunk=0, neventsChunkC=0;
  const Double_t ePerADC = 500.; 
  const Double_t fgke0 = 8.854187817e-12;  
  //
  // 
  //
  const char *inputFile="fluctuation.root";  
  TObjArray * fileList = (gSystem->GetFromPipe("cat  fluctuation.list")).Tokenize("\n");
  if (fileList->GetEntries()==0) fileList->AddLast(new TObjString(inputFile));
  Int_t nfiles  = fileList->GetEntries();
  Int_t indexPer[1000];
  Double_t numbersPer[10000];
  for (Int_t i=0; i<nfiles; i++) numbersPer[i]=gRandom->Rndm();
  TMath::Sort(nfiles, numbersPer,indexPer);

  for (Int_t ifile=0; ifile<nfiles; ifile++){
    if (nfilesMerge>0 && ifile>=nfilesMerge) continue; // merge only limited amount if specified by argument
    TFile *fhistos = TFile::Open(fileList->At(indexPer[ifile])->GetName());
    if (!fhistos) continue;
    TTree * treeHis = (TTree*)fhistos->Get("fluctuation");
    if (!treeHis) { printf("file %s does not exist or tree does not exist\n",fileList->At(ifile)->GetName()); continue;}
    Int_t nchunks=treeHis->GetEntries();
    Int_t chunk=nchunks*gRandom->Rndm();
    treeHis->SetBranchAddress("his1DROCN.",&his1DROCNC);
    treeHis->SetBranchAddress("his1DROCSum.",&his1DROCSumC);
    treeHis->SetBranchAddress("his2DRPhiROCN.",&his2DRPhiROCNC);
    treeHis->SetBranchAddress("his2DRPhiROCSum.",&his2DRPhiROCSumC);
    treeHis->SetBranchAddress("his2DRZROCN.",&his2DRZROCNC);
    treeHis->SetBranchAddress("his2DRZROCSum.",&his2DRZROCSumC);
    treeHis->SetBranchAddress("his3DROCN.",&his3DROCNC);
    treeHis->SetBranchAddress("his3DROCSum.",&his3DROCSumC);
    treeHis->SetBranchAddress("neventsAll",&neventsAllC);
    treeHis->SetBranchAddress("eventsUsed",&neventsChunkC);
    treeHis->GetEntry(chunk);  
    neventsAll+=neventsAllC;
    neventsChunk+=neventsChunkC; 
    //
    TH1 * histosC[8]={ his1DROCNC, his1DROCSumC, his2DRPhiROCNC, his2DRPhiROCSumC, his2DRZROCNC, his2DRZROCSumC, his3DROCNC, his3DROCSumC};
    if (ifile==0) for (Int_t ihis=0; ihis<8; ihis++) histos[ihis] = (TH1*)(histosC[ihis]->Clone());
    if (ifile>0)  {
      for (Int_t ihis=0; ihis<8; ihis++) histos[ihis]->Add(histosC[ihis]);
    }
  }
  his1DROCN=(TH1D*)histos[0]; his1DROCSum=(TH1D*)histos[1];
  his2DRPhiROCN=(TH2D*)histos[2];  his2DRPhiROCSum=(TH2D*)histos[3];  his2DRZROCN=(TH2D*)histos[4];  his2DRZROCSum=(TH2D*)histos[5]; 
  his3DROCN=(TH3D*)histos[6];  his3DROCSum=(TH3D*)histos[7];
  //
  // Select input histogram
  //
  TH3D * hisInput= his3DROCSum;
  Int_t neventsCorr=0;                 // number of events used for the correction studies
  if (nfilesMerge>0){
    neventsCorr=neventsChunk;
    hisInput=his3DROCN;
  }else{
    neventsCorr=neventsAll;
    hisInput=his3DROCSum;
    hisInput->Scale(Double_t(fullNorm)/Double_t(neventsAll));
  }
  
  TObjArray *distortionArray = new TObjArray; 
  TObjArray *histoArray = new TObjArray; 
  //
  // Make a reference  - ideal distortion/correction
  //
  TH3D * his3DReference =  NormalizeHistoQ(hisInput,kFALSE); // q normalized to the Q/m^3
  his3DReference->Scale(scale*0.000001/fgke0); //scale back to the C/cm^3/epsilon0
  AliTPCSpaceCharge3D *spaceChargeRef = new AliTPCSpaceCharge3D;
  spaceChargeRef->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
  spaceChargeRef->SetInputSpaceCharge(his3DReference, his2DRPhiROCSum,his2DRPhiROCSum,1);
  spaceChargeRef->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
  spaceChargeRef->AddVisualCorrection(spaceChargeRef,1);
  spaceChargeRef->SetName("DistRef");
  his3DReference->SetName("hisDistRef");
  distortionArray->AddLast(spaceChargeRef);
  histoArray->AddLast(his3DReference);
  //
  // Draw histos
  TCanvas * canvasSC = new TCanvas("canvasSCDefault","canvasSCdefault",500,400);  
  canvasSC->SetRightMargin(0.12);
  gStyle->SetTitleOffset(0.8,"z");
  canvasSC->SetRightMargin(0.13);
  spaceChargeRef->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
  canvasSC->SaveAs(Form("canvasCreateHistoDRPhiinXY_Z10_%d_%d.pdf",nfilesMerge,sign));
  spaceChargeRef->CreateHistoDRinXY(10,250,250)->Draw("colz");
  canvasSC->SaveAs(Form("canvasCreateHistoDRinXY_Z10_%d_%d.pdf",nfilesMerge,sign));
  spaceChargeRef->CreateHistoSCinZR(0.05,250,250)->Draw("colz");
  canvasSC->SaveAs(Form("canvasCreateHistoSCinZR_Phi005_%d_%d.pdf",nfilesMerge,sign));
  spaceChargeRef->CreateHistoSCinXY(10.,250,250)->Draw("colz");
  canvasSC->SaveAs(Form("canvasCreateHistoSCinRPhi_Z10_%d_%d.pdf",nfilesMerge,sign));


  //
  // Make Z scan corrections
  // 
  if (1){
  for (Int_t  ihis=1; ihis<=9; ihis+=2){ 
    TH3 *his3DZ = PermutationHistoZ(his3DReference,16*(ihis));
    AliTPCSpaceCharge3D *spaceChargeZ = new AliTPCSpaceCharge3D;
    spaceChargeZ->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
    spaceChargeZ->SetInputSpaceCharge(his3DZ, his2DRPhiROCSum,his2DRPhiROCSum,1);
    spaceChargeZ->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
    spaceChargeZ->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
    spaceChargeZ->AddVisualCorrection(spaceChargeZ,100+ihis);    
    spaceChargeZ->SetName(Form("DistZ_%d", 16*(ihis)));
    his3DZ->SetName(Form("HisDistZ_%d", 16*(ihis)));
    distortionArray->AddLast(spaceChargeZ);
    histoArray->AddLast(his3DZ);
  }
  //
  // Make Sector scan corrections
  // 
  for (Int_t  ihis=1; ihis<=9; ihis+=2){ 
    TH3 *his3DSector = PermutationHistoPhi(his3DReference,TMath::Pi()*(ihis)/9.);
    AliTPCSpaceCharge3D *spaceChargeSector = new AliTPCSpaceCharge3D;
    spaceChargeSector->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
    spaceChargeSector->SetInputSpaceCharge(his3DSector, his2DRPhiROCSum,his2DRPhiROCSum,1);
    spaceChargeSector->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
    spaceChargeSector->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
    spaceChargeSector->AddVisualCorrection(spaceChargeSector,200+ihis);    
    spaceChargeSector->SetName(Form("DistSector_%d", ihis));
    his3DSector->SetName(Form("DistSector_%d", ihis));
    distortionArray->AddLast(spaceChargeSector);
    histoArray->AddLast(his3DSector);
  } 
  //
  // Make Local phi scan smear  corrections
  // 
  for (Int_t  ihis=1; ihis<=8; ihis++){ 
    TH3 *his3DSector = PermutationHistoLocalPhi(his3DReference,ihis);
    AliTPCSpaceCharge3D *spaceChargeSector = new AliTPCSpaceCharge3D;
    spaceChargeSector->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
    spaceChargeSector->SetInputSpaceCharge(his3DSector, his2DRPhiROCSum,his2DRPhiROCSum,1);
    spaceChargeSector->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
    spaceChargeSector->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
    spaceChargeSector->AddVisualCorrection(spaceChargeSector,300+ihis);    
    spaceChargeSector->SetName(Form("DistPhi_%d", ihis));
    his3DSector->SetName(Form("HisDistPhi_%d", ihis));
    distortionArray->AddLast(spaceChargeSector); 
    histoArray->AddLast(his3DSector);
  }
 //  // 
//   // Rebin Z
//   //
//   for (Int_t  ihis=2; ihis<=8;  ihis+=2){ 
//     TH3 *his3DSector = his3DReference->RebinZ(ihis,Form("RebinZ_%d",ihis));
//     AliTPCSpaceCharge3D *spaceChargeSector = new AliTPCSpaceCharge3D;
//     spaceChargeSector->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
//     spaceChargeSector->SetInputSpaceCharge(his3DSector, his2DRPhiROCSum,his2DRPhiROCSum,1);
//     spaceChargeSector->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
//     spaceChargeSector->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
//     spaceChargeSector->AddVisualCorrection(spaceChargeSector,300+ihis);    
//     spaceChargeSector->SetName(Form("RebinZ_%d", ihis));
//     his3DSector->SetName(Form("RebinZ_%d", ihis));
//     distortionArray->AddLast(spaceChargeSector); 
//     histoArray->AddLast(his3DSector);
//   }
//   //
//   // Rebin Phi
//   //
//   for (Int_t  ihis=2; ihis<=5; ihis++){ 
//     TH3 *his3DSector = his3DReference->RebinZ(ihis,Form("RebinPhi_%d",ihis));
//     AliTPCSpaceCharge3D *spaceChargeSector = new AliTPCSpaceCharge3D;
//     spaceChargeSector->SetOmegaTauT1T2(omegaTau*bsign,1,1); // Ne CO2
//     spaceChargeSector->SetInputSpaceCharge(his3DSector, his2DRPhiROCSum,his2DRPhiROCSum,1);
//     spaceChargeSector->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
//     spaceChargeSector->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
//     spaceChargeSector->AddVisualCorrection(spaceChargeSector,300+ihis);    
//     spaceChargeSector->SetName(Form("RebinZ_%d", ihis));
//     his3DSector->SetName(Form("RebinZ_%d", ihis));
//     distortionArray->AddLast(spaceChargeSector); 
//     histoArray->AddLast(his3DSector);
//   }
  }
  //
  // Space points scan
  //
  Int_t nx = his3DROCN->GetXaxis()->GetNbins();
  Int_t ny = his3DROCN->GetYaxis()->GetNbins();
  Int_t nz = his3DROCN->GetZaxis()->GetNbins();
  Int_t nbins=nx*ny*nz;
  TVectorF  vx(nbins), vy(nbins), vz(nbins), vq(nbins), vqall(nbins);
  //
  // charge in the ROC
  // for open gate data only fraction of ions enter to drift volume
  //
  const Int_t kbins=1000;
  Double_t deltaR[kbins], deltaZ[kbins],deltaRPhi[kbins], deltaQ[kbins];
  Int_t ndist = distortionArray->GetEntries();
  for (Int_t ix=1; ix<=nx; ix+=2){    // phi bin loop
    for (Int_t iy=1; iy<=ny; iy+=2){  // r bin loop
      Double_t phi= his3DROCN->GetXaxis()->GetBinCenter(ix);
      Double_t r  = his3DROCN->GetYaxis()->GetBinCenter(iy);
      Double_t x  = r*TMath::Cos(phi); 
      Double_t y  = r*TMath::Sin(phi); 
      //
      for (Int_t iz=1; iz<=nz; iz++){ // z bin loop
	Double_t z  = his3DROCN->GetZaxis()->GetBinCenter(iz);
	Double_t qN= his3DROCN->GetBinContent(ix,iy,iz);
	Double_t qSum= his3DROCSum->GetBinContent(ix,iy,iz);
	//	Double_t dV  in cm = dphi * r * dz  	in cm**3
	Double_t dV=   (his3DROCN->GetXaxis()->GetBinWidth(ix)*r)*his3DROCN->GetZaxis()->GetBinWidth(iz);
	Double_t norm= 1e6*ePerADC*TMath::Qe()/dV;  //normalization factor to the Q/m^3 inside of the ROC;	
	(*pcstream)<<"hisDump"<<
	  "neventsAll="<<neventsAll<<         // total number of events used for the Q reference
	  "nfiles="<<nfiles<<                 // number of files to define properties
	  "nfilesMerge="<<nfilesMerge<<       // number of files to define propertiesneventsCorr
	  "neventsCorr="<<neventsCorr<<       // number of events used to define the corection
	  "fullNorm="<<fullNorm<<             // in case full statistic used this is the normalization coeficient

	  "ix="<<ix<<     
	  "iy="<<iy<<
	  "iz="<<iz<<
	  // x,y,z
	  "x="<<x<<
	  "y="<<y<<
	  "z="<<z<<
	  // phi,r,z
	  "phi="<<phi<<
	  "r="<<r<<
	  "z="<<z<<
	  "norm="<<norm<<
	  "qN="<<qN<<
	  "qSum="<<qSum;
	for (Int_t idist=0; idist<ndist; idist++){
	  AliTPCCorrection * corr  = (AliTPCCorrection *)distortionArray->At(idist);
	  TH3 * his = (TH3*)histoArray->At(idist);
	  Double_t phi0= TMath::ATan2(y,x);
	  Int_t nsector=(z>=0) ? 0:18; 
	  Float_t distPoint[3]={x,y,z};
	  corr->CorrectPoint(distPoint, nsector);
	  Double_t r0=TMath::Sqrt(x*x+y*y);
	  Double_t r1=TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
	  Double_t phi1=TMath::ATan2(distPoint[1],distPoint[0]);
	  deltaR[idist]    = r1-r0;
	  deltaRPhi[idist] = (phi1-phi0)*r0;
	  deltaZ[idist]    = distPoint[2]-z;
	  deltaQ[idist]    = his->GetBinContent(ix,iy,iz);
	  //
	  (*pcstream)<<"hisDump"<<   //correct point - input point
	    Form("%sQ=",corr->GetName())<<deltaQ[idist]<<         
	    Form("%sDR=",corr->GetName())<<deltaR[idist]<<         
	    Form("%sDRPhi=",corr->GetName())<<deltaRPhi[idist]<<         
	    Form("%sDZ=",corr->GetName())<<deltaZ[idist];
	}
	(*pcstream)<<"hisDump"<<
	  "\n";
      }
    }
  }
  delete pcstream;
  //
  // generate track distortions
  //
  if (nTracks>0){
    for(Int_t nt=1; nt<=nTracks; nt++){
      gRandom->SetSeed(nt);
      TObjArray trackArray(10000);
      Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
      Double_t eta = gRandom->Uniform(-1, 1);
      Double_t pt = 1/(gRandom->Rndm()*5+0.00001); // momentum for f1
      Short_t psign=1;
      if(gRandom->Rndm() < 0.5){
	psign =1;
      }else{
	psign=-1;
      }      
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
      Double_t pxyz[3];
      pxyz[0]=pt*TMath::Cos(phi);
      pxyz[1]=pt*TMath::Sin(phi);
      pxyz[2]=pt*TMath::Tan(theta);
      Double_t vertex[3]={0,0,0};
      Double_t cv[21]={0};
      AliExternalTrackParam *t= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
      Double_t refX0=85.;
      Double_t refX1=1.;
      Int_t dir=-1;
      (*pcstreamTrack)<<"trackFit"<<
	"neventsAll="<<neventsAll<<         // total number of events used for the Q reference
	"nfiles="<<nfiles<<                 // number of files to define properties
	"nfilesMerge="<<nfilesMerge<<       // number of files to define propertiesneventsCorr
	"neventsCorr="<<neventsCorr<<       // number of events used to define the corection
	"fullNorm="<<fullNorm<<             // in case full statistic used this is the normalization coeficient
	"itrack="<<nt<<
	"input.="<<t;
      for (Int_t idist=0; idist<ndist; idist++){
	AliTPCCorrection * corr   = (AliTPCCorrection *)distortionArray->At(idist);
	AliExternalTrackParam *ot0= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
	AliExternalTrackParam *ot1= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
	AliExternalTrackParam *td0 =  corr->FitDistortedTrack(*ot0, refX0, dir,  0);
	AliExternalTrackParam *td1 =  corr->FitDistortedTrack(*ot1, refX1, dir,  0);
	trackArray.AddLast(td0);
	trackArray.AddLast(td1);
	trackArray.AddLast(ot0);
	trackArray.AddLast(ot1);
	char name0[100], name1[1000];
	char oname0[100], oname1[1000];
	snprintf(name0, 100, "T_%s_0.=",corr->GetName());
	snprintf(name1, 100, "T_%s_1.=",corr->GetName());
	snprintf(oname0, 100, "OT_%s_0.=",corr->GetName());
	snprintf(oname1, 100, "OT_%s_1.=",corr->GetName());
	(*pcstreamTrack)<<"trackFit"<<
	  name0<<td0<< 
	  name1<<td1<< 
	  oname0<<ot0<< 
	  oname1<<ot1; 
      }
      (*pcstreamTrack)<<"trackFit"<<"\n";
    }
  }
  delete pcstreamTrack;
  return;

}


void MakePlotPoisson3D(const char *inputfile="fluctuation.root", const char *outputfile="SpaceCharge.root", Int_t event=0){
  //
  // draw "standard" plot to show radial and theta dependence of the space charge distortion
  //
  //  const char *inputfile="fluctuation.root";  const char *outputfile="SpaceCharge.root";  Int_t event=0
  //
  TFile *fhistos = TFile::Open(inputfile);
  TH2D *his2DRPhiROCN=0, *his2DRPhiROCSum=0, *his2DRZROCN=0, *his2DRZROCSum=0;
  TH1D *his1DROCN=0, *his1DROCSum=0; 
  TH3D *his3DROCN=0, *his3DROCSum=0; 
  const Double_t ePerADC = 500.; 
  const Double_t fgke0 = 8.854187817e-12;  
  TTree * treeHis = (TTree*)fhistos->Get("fluctuation");
  treeHis->SetBranchAddress("his1DROCN.",&his1DROCN);
  treeHis->SetBranchAddress("his1DROCSum.",&his1DROCSum);
  treeHis->SetBranchAddress("his2DRPhiROCN.",&his2DRPhiROCN);
  treeHis->SetBranchAddress("his2DRPhiROCSum.",&his2DRPhiROCSum);
  treeHis->SetBranchAddress("his2DRZROCN.",&his2DRZROCN);
  treeHis->SetBranchAddress("his2DRZROCSum.",&his2DRZROCSum);
  treeHis->SetBranchAddress("his3DROCN.",&his3DROCN);
  treeHis->SetBranchAddress("his3DROCSum.",&his3DROCSum);
  treeHis->GetEntry(event);
  
  his3DROCSum->Scale(ePerADC*TMath::Qe()/fgke0); 

  AliTPCSpaceCharge3D *spaceChargeOrig = new AliTPCSpaceCharge3D;
  spaceChargeOrig->SetOmegaTauT1T2(0.0,1,1); // Ne CO2
  spaceChargeOrig->SetInputSpaceCharge(his3DROCSum, his2DRPhiROCSum,his2DRPhiROCSum,10*ePerADC*TMath::Qe());
  spaceChargeOrig->InitSpaceCharge3DPoisson(129, 129, 144,100);
  spaceChargeOrig->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
  spaceChargeOrig->AddVisualCorrection(spaceChargeOrig,1);
  //
  //AliTPCSpaceCharge3D *spaceChargeRef= spaceChargeOrig;



  //
  Int_t nfuns=5;
  Double_t dmax=0.75, dmin=-0.75;
  Double_t phiRange=18;
  TCanvas *canvasDistortionP3D = new TCanvas("canvasdistortionP3D","canvasdistortionP3D",1000,700);
  canvasDistortionP3D->SetGrid(1,1);
  canvasDistortionP3D->Divide(1,2);
  canvasDistortionP3D->cd(1)->SetGrid(1,1);    
  TLegend * legendR= new TLegend(0.11,0.11,0.45,0.35,"R scan (#Theta=0.1)");
  for (Int_t ifun1=0; ifun1<=nfuns; ifun1++){    
    Double_t rfun= 85.+ifun1*(245.-85.)/nfuns;
    TF1 *pf1 = new TF1("f1",Form("AliTPCCorrection::GetCorrSector(x,%f,0.1,1,1)",rfun),0,phiRange);
    pf1->SetMinimum(dmin);
    pf1->SetMaximum(dmax);
    pf1->SetNpx(360);
    pf1->SetLineColor(1+ifun1);
    pf1->SetLineWidth(2);    
    pf1->GetXaxis()->SetTitle("sector");
    pf1->GetXaxis()->SetNdivisions(530);
    pf1->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
    if (ifun1==0) pf1->Draw();
    pf1->Draw("same");
    legendR->AddEntry(pf1,Form("r=%1.0f",rfun));
  }
  legendR->Draw();
  //
  canvasDistortionP3D->cd(2)->SetGrid(1,1);
  TLegend * legendTheta= new TLegend(0.11,0.11,0.45,0.35,"#Theta scan (r=125 cm)");
  for (Int_t ifun1=0; ifun1<=nfuns; ifun1++){    
    Double_t tfun= 0.1+ifun1*(0.8)/nfuns;
    TF1 *pf1 = new TF1("f1",Form("AliTPCCorrection::GetCorrSector(x,125,%f,1,1)",tfun),0,phiRange);
    pf1->SetMinimum(dmin);
    pf1->SetMaximum(dmax);
    pf1->SetNpx(360);
    pf1->SetLineColor(1+ifun1);
    pf1->SetLineWidth(2);    
    pf1->GetXaxis()->SetTitle("sector");
    pf1->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
    pf1->GetXaxis()->SetNdivisions(530);
    if (ifun1==0) pf1->Draw();
    pf1->Draw("same");
    legendTheta->AddEntry(pf1,Form("#Theta=%1.2f",tfun));
  }
  legendTheta->Draw();

}

TH3D *  NormalizeHistoQ(TH3D * hisInput, Bool_t normEpsilon){
  //
  // Renormalize the histogram to the Q/m^3
  // Input:
  //   hisInput     - input 3D histogram
  //   normEpsilon  - flag - normalize to epsilon0
  //
  const Double_t ePerADC = 500.; 
  const Double_t fgkEpsilon0 = 8.854187817e-12;  
  TH3D * hisOutput= new TH3D(*hisInput);
  Int_t nx = hisInput->GetXaxis()->GetNbins();
  Int_t ny = hisInput->GetYaxis()->GetNbins();
  Int_t nz = hisInput->GetZaxis()->GetNbins();
  for (Int_t ix=1; ix<=nx; ix++){
    for (Int_t iy=1; iy<=ny; iy++){
      for (Int_t iz=1; iz<=nz; iz++){
	//	Double_t z = hisInput->GetZaxis()->GetBinCenter(iz);
	Double_t deltaRPhi = hisInput->GetXaxis()->GetBinWidth(ix)* hisInput->GetYaxis()->GetBinCenter(iy);
	Double_t deltaR= hisInput->GetYaxis()->GetBinWidth(iy);
	Double_t deltaZ= hisInput->GetYaxis()->GetBinWidth(iz);	
	Double_t volume= (deltaRPhi*deltaR*deltaZ)/1000000.;
	Double_t q   = hisInput->GetBinContent(ix,iy,iz)* ePerADC*TMath::Qe(); // Q in coulombs
	Double_t rho = q/volume;      // rpho - density in Q/m^3
	if (normEpsilon) rho/=fgkEpsilon0;
	hisOutput->SetBinContent(ix,iy,iz,rho);
      }
    }
  }
  return hisOutput;
}



TH3D *  PermutationHistoZ(TH3D * hisInput, Double_t deltaZ){
  //
  // Used to estimate the effect of the imperfection of the lookup tables as function of update frequency
  //
  // Permute/rotate the conten of the histogram in z direction
  // Reshufle/shift content -  Keeping the integral the same
  // Parameters:
  //    hisInput - input 3D histogram (phi,r,z)
  //    deltaZ   - deltaZ -shift of the space charge
  Double_t zmax=250;
  TH3D * hisOutput= new TH3D(*hisInput);
  Int_t nx = hisInput->GetXaxis()->GetNbins();
  Int_t ny = hisInput->GetYaxis()->GetNbins();
  Int_t nz = hisInput->GetZaxis()->GetNbins();
  //
  //
  for (Int_t ix=1; ix<=nx; ix++){
    for (Int_t iy=1; iy<=ny; iy++){
      for (Int_t iz=1; iz<=nz; iz++){
	Double_t zold = hisInput->GetZaxis()->GetBinCenter(iz);
	Double_t z=zold;
	if (z>0){
	  z+=deltaZ;
	  if (z<0) z+=zmax;
	  if (z>zmax) z-=zmax;
	}else{
	  z-=deltaZ;
	  if (z>0) z-=zmax;
	  if (z<-zmax) z+=zmax;	}
	Double_t kz= hisInput->GetZaxis()->FindBin(z);
	Double_t content = hisInput->GetBinContent(ix,iy,iz);
	hisOutput->SetBinContent(ix,iy,kz,content);
      }
    }
  }
  return hisOutput;
}




TH3D *  PermutationHistoPhi(TH3D * hisInput, Double_t deltaPhi){
  //
  // Used to estimate the effect of the imperfection of the lookup tables as function of update frequency
  //
  // Permute/rotate the conten of the histogram in phi
  // Reshufle/shift content -  Keeping the integral the same
  // Parameters:
  //    hisInput - input 3D histogram (phi,r,z)
  //    deltaPhi   - deltaPhi -shift of the space charge
  TH3D * hisOutput= new TH3D(*hisInput);
  Int_t nx = hisInput->GetXaxis()->GetNbins();
  Int_t ny = hisInput->GetYaxis()->GetNbins();
  Int_t nz = hisInput->GetZaxis()->GetNbins();
  //
  //
  for (Int_t iy=1; iy<=ny; iy++){
    for (Int_t iz=1; iz<=nz; iz++){
      for (Int_t ix=1; ix<=nx; ix++){
	Double_t phiOld = hisInput->GetXaxis()->GetBinCenter(ix);
	Double_t phi=phiOld;
	phi+=deltaPhi;
	if (phi<0) phi+=TMath::TwoPi();
	if (phi>TMath::TwoPi()) phi-=TMath::TwoPi();	
	Double_t kx= hisInput->GetXaxis()->FindBin(phi);
	Double_t content = hisInput->GetBinContent(ix,iy,iz);
	hisOutput->SetBinContent(kx,iy,iz,content);
      }
    }
  }
  return hisOutput;
}


TH3D *  PermutationHistoLocalPhi(TH3D * hisInput, Int_t deltaPhi){
  //
  // Used to estimate the effect of the imperfection of the lookup tables as function of update frequency
  // Use moving average of the content  instead of the content
  //
  // Parameters:
  //    hisInput - input 3D histogram (phi,r,z)
  //    deltaPhi   - moving average width
  TH3D * hisOutput= new TH3D(*hisInput);
  Int_t nx = hisInput->GetXaxis()->GetNbins();
  Int_t ny = hisInput->GetYaxis()->GetNbins();
  Int_t nz = hisInput->GetZaxis()->GetNbins();
  Int_t binSector=nx/18;
  //
  //
  for (Int_t iy=1; iy<=ny; iy++){
    for (Int_t iz=1; iz<=nz; iz++){
      for (Int_t ix=1; ix<=nx; ix++){
	Double_t sumRo=0,sumW=0;
	for (Int_t idx=-deltaPhi; idx<=deltaPhi; idx++){
	  Int_t index=ix+idx;
	  if (index<1) index+=nx+1;  // underflow and overflow bins
	  if (index>nx) index-=nx+1;
	  Double_t content = hisInput->GetBinContent(index,iy,iz);
	  sumRo+=content;
	  sumW++;
	}
	Double_t meanCont= sumRo/sumW;
	hisOutput->SetBinContent(ix,iy,iz,meanCont);	
	//printf("%d\t%f\n",ix,hisInput->GetBinContent(ix,iy,iz)/(hisInput->GetBinContent(ix,iy,iz)+meanCont));
      }	
    }
  }
  return hisOutput;
}



void ScanIterrationPrecision(TH3 * hisInput, Int_t offset){
  //
  //
  //
  for (Int_t iter=0; iter<=7; iter++){
    Int_t niter= 50.*TMath::Power(1.5,iter);
    AliTPCSpaceCharge3D *spaceChargeOrig = new AliTPCSpaceCharge3D;
    spaceChargeOrig->SetOmegaTauT1T2(0.0,1,1); // Ne CO2
    spaceChargeOrig->SetInputSpaceCharge(hisInput,0,0,1);
    spaceChargeOrig->InitSpaceCharge3DPoisson(129, 129, 144,niter);
    spaceChargeOrig->CreateHistoDRPhiinXY(10,250,250)->Draw("colz");
    spaceChargeOrig->AddVisualCorrection(spaceChargeOrig,offset+iter+1);
  }
}


void DrawFluctuationSector(Int_t stat, Double_t norm){
  //
  // Draw correction - correction  at rotated sector  
  // The same set of events used
  // Int_t stat=0; Double_t norm=10000;
  // 
  // Notes:
  //    1. (something wrong for the setup 2 pileups  -problem with data 24.07
  //
  //
  TFile *f0= TFile::Open(Form("SpaceChargeFluc%d.root",stat));
  TTree * tree0 = (TTree*)f0->Get("hisDump");
  tree0->SetCacheSize(1000000000);
  tree0->SetMarkerStyle(25);
  TObjArray * fitArray=new TObjArray(3);
  tree0->SetAlias("scNorm",Form("%f/neventsCorr",norm));
  //
  // Sector Scan
  //
  TH2 * hisSectorScan[5]={0};
  TH1 * hisSectorScanSigma[5]={0};
  for (Int_t ihis=0; ihis<5; ihis++){
    tree0->Draw(Form("(DistRefDR-DistSector_%dDR)*scNorm:r>>hisSec%d(50,84,245,100,-1,1)",ihis*2+1,ihis*2+1),"abs(z)<90","colzgoff");
    hisSectorScan[ihis]=(TH2*)tree0->GetHistogram();
    hisSectorScan[ihis]->FitSlicesY(0,0,-1,0,"QNR",fitArray);
    hisSectorScanSigma[ihis]=(TH1*)(fitArray->At(2)->Clone());
    hisSectorScanSigma[ihis]->SetMinimum(0);
    hisSectorScanSigma[ihis]->SetMaximum(0.2);
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TCanvas * canvasFlucSectorScan=new TCanvas("canvasFlucSectorScan","canvasFlucSectorScan",750,700);
  canvasFlucSectorScan->Divide(2,2,0,0);  
  gStyle->SetPadBorderMode(0);
  for (Int_t ihis=0; ihis<4; ihis++){
    canvasFlucSectorScan->cd(ihis+1)->SetLogz(1);
    hisSectorScan[ihis]->GetXaxis()->SetTitle("r (cm)");
    hisSectorScan[ihis]->GetYaxis()->SetTitle("#Delta_{R} (cm)");
    hisSectorScan[ihis]->Draw("colz");
    TLegend * legendSec=new TLegend(0.5,0.7,0.89,0.89);
    legendSec->AddEntry(hisSectorScan[ihis],Form("Sector #Delta %d",(ihis*2+1))); 
    legendSec->Draw();
  }
  canvasFlucSectorScan->SaveAs("canvasFlucSectorScan.pdf");
  canvasFlucSectorScan->SaveAs("canvasFlucSectorScan.png");
  //
  gStyle->SetOptTitle(0);
  TCanvas * canvasFlucSectorScanFit=new TCanvas("canvasFlucSectorScanFit","canvasFlucSectorScanFit",750,550);
  TLegend * legendSector = new TLegend(0.50,0.55,0.89,0.89,"Space charge: corr(sec)-corr(sec-#Delta_{sec})");
  for (Int_t ihis=0; ihis<5; ihis++){
    hisSectorScanSigma[ihis]->GetXaxis()->SetTitle("r (cm)");
    hisSectorScanSigma[ihis]->GetYaxis()->SetTitle("#sigma(#Delta_{R}) (cm)");
    hisSectorScanSigma[ihis]->SetMarkerStyle(21+ihis%5);
    hisSectorScanSigma[ihis]->SetMarkerColor(1+ihis%4);
    if (ihis==0) hisSectorScanSigma[ihis]->Draw("");
    hisSectorScanSigma[ihis]->Draw("same");
    legendSector->AddEntry(hisSectorScanSigma[ihis],Form("#Delta %d",(ihis*2+1)));
  }
  legendSector->Draw();
  canvasFlucSectorScanFit->SaveAs("canvasFlucSectorScanFit.pdf");
  canvasFlucSectorScanFit->SaveAs("canvasFlucSectorScanFit.png");
}



void DrawFluctuationdeltaZ(Int_t stat, Double_t norm){
  //
  // Draw correction - correction  shifted z  
  // The same set of events used
  //Int_t stat=0; Double_t norm=10000;
  Int_t deltaZ=16.;
  TFile *f0= TFile::Open(Form("SpaceChargeFluc%d.root",stat));
  TTree * tree0 = 0;
  if (f0) tree0 = (TTree*)f0->Get("hisDump");
  if (!tree0){
    tree0 = AliXRDPROOFtoolkit::MakeChainRandom("space.list","hisDump",0,10);
  }
  tree0->SetCacheSize(1000000000);
  tree0->SetMarkerStyle(25);
  TObjArray * fitArray=new TObjArray(3);  
  tree0->SetAlias("scNorm",Form("%f/neventsCorr",norm));
  //
  // DeltaZ Scan
  //
  TH2 * hisDeltaZScan[6]={0};
  TH1 * hisDeltaZScanSigma[6]={0};
  for (Int_t ihis=0; ihis<6; ihis++){
    tree0->Draw(Form("(DistRefDR-DistZ_%dDR)*scNorm:r>>hisZ%d(50,84,245,100,-1,1)",(ihis+1)*deltaZ,(ihis+1)*deltaZ),"abs(z/r)<1","colzgoff");
    hisDeltaZScan[ihis]=(TH2*)tree0->GetHistogram();
    hisDeltaZScan[ihis]->FitSlicesY(0,0,-1,0,"QNR",fitArray);
    hisDeltaZScanSigma[ihis]=(TH1*)(fitArray->At(2)->Clone());
    hisDeltaZScanSigma[ihis]->SetMinimum(0);
    hisDeltaZScanSigma[ihis]->SetMaximum(0.2);
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TCanvas * canvasFlucDeltaZScan=new TCanvas("canvasFlucDeltaZScan","canvasFlucDeltaZScan",700,700);
  canvasFlucDeltaZScan->Divide(3,2,0,0);  
  gStyle->SetPadBorderMode(0);
  for (Int_t ihis=0; ihis<6; ihis++){
    canvasFlucDeltaZScan->cd(ihis+1)->SetLogz(1);
    hisDeltaZScan[ihis]->GetXaxis()->SetTitle("r (cm)");
    hisDeltaZScan[ihis]->GetYaxis()->SetTitle("#Delta_{R} (cm)");
    hisDeltaZScan[ihis]->Draw("colz");
    TLegend * legendSec=new TLegend(0.5,0.7,0.89,0.89);
    legendSec->AddEntry(hisDeltaZScan[ihis],Form("DeltaZ #Delta %d",(ihis+1)*deltaZ)); 
    legendSec->Draw();
  }
  canvasFlucDeltaZScan->SaveAs(Form("canvasFlucDeltaZScan%d.pdf",stat));
  canvasFlucDeltaZScan->SaveAs(Form("canvasFlucDeltaZScan%d.png",stat));

  //
  gStyle->SetOptTitle(0);
  TCanvas * canvasFlucDeltaZScanFit=new TCanvas("canvasFlucDeltaZScanFit","canvasFlucDeltaZScanFit");
  TLegend * legendDeltaZ = new TLegend(0.50,0.55,0.89,0.89,"Space charge: corr(z_{ref})-corr(z_{ref}-#Delta_{z})");
  for (Int_t ihis=0; ihis<5; ihis++){
    hisDeltaZScanSigma[ihis]->GetXaxis()->SetTitle("r (cm)");
    hisDeltaZScanSigma[ihis]->GetYaxis()->SetTitle("#sigma(#Delta_{R}) (cm)");
    hisDeltaZScanSigma[ihis]->SetMarkerStyle(21+ihis%5);
    hisDeltaZScanSigma[ihis]->SetMarkerColor(1+ihis%4);
    if (ihis==0) hisDeltaZScanSigma[ihis]->Draw("");
    hisDeltaZScanSigma[ihis]->Draw("same");
    legendDeltaZ->AddEntry(hisDeltaZScanSigma[ihis],Form("#Delta %d (cm)",(ihis+1)*deltaZ));
  }
  legendDeltaZ->Draw();
  canvasFlucDeltaZScanFit->SaveAs(Form("canvasFlucDeltaZScanFit%d.pdf",stat));
  canvasFlucDeltaZScanFit->SaveAs(Form("canvasFlucDeltaZScanFit%d.png",stat));

}


void DrawDefault(Int_t stat){
  //
  // Draw correction - correction  shifted z  
  // The same set of events used
  //  Int_t stat=0
  TFile *f0= TFile::Open(Form("SpaceChargeFluc%d.root",stat));
  TTree * tree0 = (TTree*)f0->Get("hisDump");
  tree0->SetCacheSize(1000000000);
  tree0->SetMarkerStyle(25);
  tree0->SetMarkerSize(0.4);
  //  TObjArray * fitArray=new TObjArray(3);
  tree0->Draw("10000*DistRefDR/neventsCorr:r:z/r","abs(z/r)<0.9&&z>0&&rndm>0.8","colz");
}




void DrawTrackFluctuation(){
  //
  // Function to make a fluctuation figures for differnt multiplicities of pileup space charge
  // it is assumed that the text files  
  //
  //
  TObjArray arrayFit(3);
  const char *inputList;
  TH2F * hisCorrRef[5]={0};
  TH2F * hisCorrNo[5]={0};
  TH1  * hisCorrRefM[5], *hisCorrRefRMS[5];
  TH1  * hisCorrNoM[5], *hisCorrNoRMS[5];
  //
  // 1. Load chains for different statistic
  //  
  TCut cutOut="abs(T_DistRef_0.fX-OT_DistRef_0.fX)<0.1&&T_DistRef_0.fX>1&&abs(OT_DistRef_0.fP[4])<4";
  TCut cutOutF="abs(R.T_DistRef_0.fX-R.OT_DistRef_0.fX)<0.1&&R.T_DistRef_0.fX>1&&abs(R.OT_DistRef_0.fP[4])<4";
  TChain * chains[5]={0};
  TChain * chainR = AliXRDPROOFtoolkit::MakeChain("track0_1.list","trackFit",0,1000);
  chainR->SetCacheSize(1000000000);
  for (Int_t ichain=0; ichain<5; ichain++){
    chains[ichain] = AliXRDPROOFtoolkit::MakeChain(Form("track%d_1.list",2*(ichain+1)),"trackFit",0,1000);
    chains[ichain]->AddFriend(chainR,"R");
    chains[ichain]->SetCacheSize(1000000000);
    chains[ichain]->SetMarkerStyle(25);
    chains[ichain]->SetMarkerSize(0.5);
    chains[ichain]->SetAlias("meanNorm","(1+0.2*abs(neventsCorr/10000-1))"); // second order correction - renomalization of mean hardwired  
    // to be fitted?
  }
  //
  // 2. fill histograms if not available in file
  //    
  // 
  TFile *ftrackFluctuation = TFile::Open("trackFluctuation.root","update");
  for (Int_t ihis=0; ihis<5; ihis++){
    ftrackFluctuation->cd();
    hisCorrRef[ihis] = (TH2F*)(ftrackFluctuation->Get(Form("DeltaRPhiCorr%d",(ihis+1)*2000)));
    hisCorrNo[ihis]  = (TH2F*)(ftrackFluctuation->Get(Form("DeltaRPhi%d",(ihis+1)*2000)));
    if (hisCorrRef[ihis]==0) {
      chains[ihis]->Draw("(T_DistRef_0.fP[0]/meanNorm-neventsCorr*R.T_DistRef_0.fP[0]/10000):R.OT_DistRef_0.fP[4]>>his(10,-4,4,100,-0.25,0.25)",cutOut+cutOutF+"","colzgoff");
      hisCorrRef[ihis]=(TH2F*)(chains[ihis]->GetHistogram()->Clone());  
      hisCorrRef[ihis]->SetName(Form("DeltaRPhiCorr%d",(ihis+1)*2000));
      hisCorrRef[ihis]->SetTitle(Form("Corrected #Delta r#phi -  Pileup %d",(ihis+1)*2000));
      hisCorrRef[ihis]->GetXaxis()->SetTitle("1/p_{t} (1/GeV/c)");
      hisCorrRef[ihis]->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
      hisCorrRef[ihis]->Write();
      //
      chains[ihis]->Draw("(T_DistRef_0.fP[0]/meanNorm):R.OT_DistRef_0.fP[4]>>hisCorNo(10,-3,3,100,-4,4)",cutOut+cutOutF+"","colzgoff");
      hisCorrNo[ihis]=(TH2F*)(chains[ihis]->GetHistogram()->Clone());  
      hisCorrNo[ihis]->SetName(Form("DeltaRPhi%d",(ihis+1)*2000));
      hisCorrNo[ihis]->SetTitle(Form("Delta r#phi  = Pileup %d",(ihis+1)*2000));
      hisCorrNo[ihis]->GetXaxis()->SetTitle("1/p_{t} (1/GeV/c)");
      hisCorrNo[ihis]->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
      hisCorrNo[ihis]->Write();    
    }
  }
  ftrackFluctuation->Flush();
  //
  //
  //
  for (Int_t ihis=0; ihis<5; ihis++){
    hisCorrRef[ihis]->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    hisCorrRefM[ihis] = (TH1*)arrayFit.At(1)->Clone();
    hisCorrRefRMS[ihis] = (TH1*)arrayFit.At(2)->Clone();
    hisCorrRefM[ihis]->GetXaxis()->SetTitle("1/p_{t} (1/GeV/c)");
    hisCorrRefM[ihis]->GetYaxis()->SetTitle("#Delta_{r#phi} (cm)");
    hisCorrRefM[ihis]->SetMarkerStyle(20);
    hisCorrRefRMS[ihis]->SetMarkerStyle(21);
    hisCorrRefM[ihis]->SetMarkerColor(1);
    hisCorrRefRMS[ihis]->SetMarkerColor(2);
    hisCorrNo[ihis]->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    hisCorrNoM[ihis] = (TH1*)arrayFit.At(1)->Clone();
    hisCorrNoRMS[ihis] = (TH1*)arrayFit.At(2)->Clone();
  }

  //
  TCanvas *canvasMean = new TCanvas("canvasCorrectionMean","canvasCorrectionMean",900,1000);
  TCanvas *canvasMeanSummary = new TCanvas("canvasCorrectionMeanSummary","canvasCorrectionMeanSummary",700,600);

  canvasMean->Divide(3,5);
  gStyle->SetOptStat(0);
  for (Int_t ihis=0; ihis<5; ihis++){
    TLegend * legend = new TLegend(0.11,0.11,0.5,0.3,Form("Pile up %d",(ihis+1)*2000));
    canvasMean->cd(3*ihis+1);
    hisCorrNo[ihis]->Draw("colz"); 
    canvasMean->cd(3*ihis+2);
    hisCorrRef[ihis]->Draw("colz");    
    canvasMean->cd(3*ihis+3); 
    hisCorrRefM[ihis]->SetMaximum(0.25);
    hisCorrRefM[ihis]->SetMinimum(-0.25);
    hisCorrRefM[ihis]->Draw("");
    hisCorrRefRMS[ihis]->Draw("same");
    legend->AddEntry(hisCorrRefM[ihis],"Mean");
    legend->AddEntry(hisCorrRefRMS[ihis],"RMS");
    legend->Draw();
  }
  canvasMeanSummary->cd();
  TLegend * legendMeanSummary = new TLegend(0.5,0.6,0.89,0.89,"Space charge correction fluctuation in r#phi"); 
  for (Int_t ihis=4; ihis>=0; ihis--){    
    hisCorrRefRMS[ihis]->SetMarkerColor(1+ihis);
    hisCorrRefRMS[ihis]->SetMinimum(0);
    hisCorrRefRMS[ihis]->GetYaxis()->SetTitle("#sigma_{r#phi} (cm)");
    if (ihis==4) hisCorrRefRMS[ihis]->Draw("");
    hisCorrRefRMS[ihis]->Draw("same");
    legendMeanSummary->AddEntry(hisCorrRefRMS[ihis],Form("%d pile-up events",(ihis+1)*2000));
  }
  legendMeanSummary->Draw();

  canvasMean->SaveAs("canvasCorrectionMean.pdf"); 
  canvasMeanSummary->SaveAs("canvasCorrectionMeanSummary.pdf");
  //canvasMean->Write();
  //canvasMeanSummary->Write();
  ftrackFluctuation->Close();
}

void DrawTrackFluctuationZ(){
  //
  // Draw track fucutation dz
  //   
  const Int_t kColors[6]={1,2,3,4,6,7};
  const Int_t kStyle[6]={20,21,24,25,24,25};
  TObjArray arrayFit(3);
  TCut cutOut="abs(T_DistRef_0.fX-OT_DistRef_0.fX)<0.1&&T_DistRef_0.fX>1&&abs(OT_DistRef_0.fP[4])<4";
  TCut cutOutF="abs(R.T_DistRef_0.fX-R.OT_DistRef_0.fX)<0.1&&R.T_DistRef_0.fX>1&&abs(R.OT_DistRef_0.fP[4])<4";
  TChain * chains[5]={0};
  TChain * chainR = AliXRDPROOFtoolkit::MakeChain("track0_1.list","trackFit",0,1000);
  chainR->SetCacheSize(1000000000);
  for (Int_t ichain=0; ichain<5; ichain++){
    chains[ichain] = AliXRDPROOFtoolkit::MakeChain(Form("track%d_1.list",2*(ichain+1)),"trackFit",0,1000);
    chains[ichain]->AddFriend(chainR,"R");
    chains[ichain]->SetCacheSize(1000000000);
    chains[ichain]->SetMarkerStyle(25);
    chains[ichain]->SetMarkerSize(0.5);
  }
  //
  // 2.) Create 2D histo or read from files
  //
  TObjArray * arrayHisto = new TObjArray(25);
  TFile *ftrackFluctuationZ = TFile::Open("trackFluctuationZ.root","update");
  for (Int_t ihis=0; ihis<5; ihis++){
    ftrackFluctuationZ->cd();
    for (Int_t idz=0; idz<5; idz++){
      Int_t z= 16+idz*32;
      TH2 *his= (TH2*)ftrackFluctuationZ->Get(Form("TrackDz%d_PileUp%d",z, (ihis+1)*2000));
      if (!his){
	chains[ihis]->Draw(Form("T_DistZ_%d_0.fP[0]-T_DistRef_0.fP[0]:T_DistRef_0.fP[4]>>his(10,-4,4,100,-0.25,0.25)",z),cutOut+"","colz");
	his = (TH2*)(chains[ihis]->GetHistogram()->Clone());
	his->SetName(Form("TrackDz%d_PileUp%d",z, (ihis+1)*2000));
	his->Write();
      }
      arrayHisto->AddAtAndExpand(his,ihis*5+idz);
    }
  }
  ftrackFluctuationZ->Flush();

  //
  // 3.) Make fits 
  //
  TCanvas *canvasDz = new TCanvas("canvasDz","canvasDz",800,800);
  canvasDz->Divide(2,2,0,0);
  for (Int_t ihis=3; ihis>=0; ihis--){
    canvasDz->cd(ihis+1)->SetTicks(3);
    TLegend * legend  = new TLegend(0.31,0.51, 0.95,0.95,Form("Distortion due time/z delay (Pileup=%d)", (ihis+1)*2000)); 
    legend->SetBorderSize(0);
    for (Int_t idz=3; idz>=0; idz--){
      TH2 * his =  (TH2*)arrayHisto->At(ihis*5+idz);
      his->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
      TH1 * hisRMS = (TH1*)arrayFit.At(2)->Clone();
      hisRMS->SetMaximum(0.12);
      hisRMS->SetMinimum(0);
      hisRMS->GetXaxis()->SetTitle("1/p_{t} (GeV/c)");
      hisRMS->GetYaxis()->SetTitle("#sigma_{r#phi}(cm)");
      hisRMS->SetMarkerStyle(kStyle[idz]);
      hisRMS->SetMarkerColor(kColors[idz]);
      if (idz==3)     hisRMS->Draw();
      legend->AddEntry(hisRMS,Form("#Delta_{z}=%d (cm)",16+idz*32));
      hisRMS->Draw("same");
    }
    legend->Draw();
  }
  canvasDz->SaveAs("spaceChargeDeltaZScan.pdf");

}





void DrawTrackFluctuationFrame(){
  //
  // Function to make a fluctuation figures for differnt multiplicities of pileup space charge
  // it is assumed that the text files  
  //
  //
  TObjArray arrayFit(3);
  const char *inputList;
  TH2F * hisCorrRef[10]={0};
  TH2F * hisCorrNo[10]={0};
  TH1  * hisCorrRefM[10], *hisCorrRefRMS[10];
  TH1  * hisCorrNoM[10], *hisCorrNoRMS[10];
  //
  // 1. Load chains for different statistic
  //  
  TCut cutOut="abs(T_DistRef_0.fX-OT_DistRef_0.fX)<0.1&&T_DistRef_0.fX>1&&abs(OT_DistRef_0.fP[4])<4";
  TCut cutOutF="abs(R.T_DistRef_0.fX-R.OT_DistRef_0.fX)<0.1&&R.T_DistRef_0.fX>1&&abs(R.OT_DistRef_0.fP[4])<4";
  TCut cutFit="Entry$%4==0";  //use only subset of data for fit 

  TChain * chains[10]={0};
  TChain * chainR = AliXRDPROOFtoolkit::MakeChain("track0_1.list","trackFit",0,1000);
  chainR->SetCacheSize(1000000000);
  for (Int_t ichain=0; ichain<7; ichain++){
    chains[ichain] = AliXRDPROOFtoolkit::MakeChain(Form("track%d_1.list",2*(ichain+1)),"trackFit",0,1000);
    chains[ichain]->AddFriend(chainR,"R");
    chains[ichain]->SetCacheSize(1000000000);
    chains[ichain]->SetMarkerStyle(25);
    chains[ichain]->SetMarkerSize(0.5);
    chains[ichain]->SetAlias("meanNorm","(1+0.2*abs(neventsCorr/10000-1))"); // second order correction - renomalization of mean hardwired  
    chains[ichain]->SetAlias("dMean0","(neventsCorr*R.T_DistRef_0.fP[0]/10000)");
    chains[ichain]->SetAlias("dMeas0","T_DistRef_0.fP[0]");
    chains[ichain]->SetAlias("dMean1","(neventsCorr*R.T_DistRef_1.fP[0]/10000)");
    chains[ichain]->SetAlias("dMeas1","T_DistRef_1.fP[0]"); 
    for (Int_t ig=0; ig<10;ig++) chains[ichain]->SetAlias(Form("FR%d",ig),Form("(abs(Entry$-%d)<1000)",ig*2000+1000));
  }
  //
  // 2.  Get or Create histogram (do fit per frame)
  //   
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;  
  TString  fstringG="";              // global part
  fstringG+="dMean0++";  
  TVectorD vec0,vec1;
  TString  fstringF0="";              // global part
  for (Int_t ig=0; ig<10;ig++) fstringF0+=Form("FR%d++",ig);
  for (Int_t ig=0; ig<10;ig++) fstringF0+=Form("FR%d*dMean0++",ig);
  TString  fstringF1="";              // global part
  for (Int_t ig=0; ig<10;ig++) fstringF1+=Form("FR%d++",ig);
  for (Int_t ig=0; ig<10;ig++) fstringF1+=Form("FR%d*dMean0++",ig);
  for (Int_t ig=0; ig<10;ig++) fstringF1+=Form("FR%d*dMean0*abs(T_DistRef_0.fP[3])++",ig);
  for (Int_t ig=0; ig<10;ig++) fstringF1+=Form("FR%d*dMean0*(T_DistRef_0.fP[3]^2)++",ig);

  //
  //
  TH2F *hisA=0, *hisF0=0, *hisF1=0, *hisM=0;
  TObjArray * arrayHisto = new TObjArray(200);
  TFile *ftrackFit = TFile::Open("trackFluctuationFrame.root","update");
  for (Int_t ihis=0; ihis<7; ihis++){
    printf("\n\nProcessing frames\t%d\nnn",(ihis+1)*2000);
    hisM = (TH2F*)ftrackFit->Get(Form("hisMean_%d",(ihis+1)*2000));
    hisA = (TH2F*)ftrackFit->Get(Form("hisAll_%d",(ihis+1)*2000));
    hisF0 = (TH2F*)ftrackFit->Get(Form("hisFrame0_%d",(ihis+1)*2000));
    hisF1 = (TH2F*)ftrackFit->Get(Form("hisFrame1_%d",(ihis+1)*2000));
    if (!hisA){    
      ftrackFit->cd();
      TString * fitResultAll = TStatToolkit::FitPlane(chains[ihis],"dMeas0", fstringG.Data(),cutOut+cutOutF+cutFit, chi2,npoints,param,covar,-1,0, 40*2000, kFALSE);
      chains[ihis]->SetAlias("fitAll",fitResultAll->Data());  
      TString * fitResultF0 = TStatToolkit::FitPlane(chains[ihis],"dMeas0", fstringF0.Data(),cutOut+cutOutF+cutFit+"abs(dMeas0-fitAll)<0.3", chi2,npoints,vec0,covar,-1,0, 10*2000, kFALSE);
      chains[ihis]->SetAlias("fitF0",fitResultF0->Data());  
      TString * fitResultF1 = TStatToolkit::FitPlane(chains[ihis],"dMeas0", fstringF1.Data(),cutOut+cutOutF+cutFit+"abs(dMeas0-fitAll)<0.3", chi2,npoints,vec1,covar,-1,0, 10*2000, kFALSE);
      chains[ihis]->SetAlias("fitF1",fitResultF1->Data());  
      fitResultF0->Tokenize("++")->Print();
      chains[ihis]->Draw(Form("dMeas0-fitAll:T_DistRef_0.fP[4]>>hisAll_%d(20,-4,4,100,-0.25,0.25)",(ihis+1)*2000),cutOut+cutOutF,"colz",100000,0);   
      hisA = (TH2F*)chains[ihis]->GetHistogram();
      chains[ihis]->Draw(Form("dMeas0-fitF0:T_DistRef_0.fP[4]>>hisFrame0_%d(20,-4,4,100,-0.10,0.10)",(ihis+1)*2000),cutOut+cutOutF,"colz",20000,0);
      hisF0 = (TH2F*)chains[ihis]->GetHistogram();
      chains[ihis]->Draw(Form("dMeas0-fitF1:T_DistRef_0.fP[4]>>hisFrame1_%d(20,-4,4,100,-0.10,0.10)",(ihis+1)*2000),cutOut+cutOutF,"colz",20000,0);
      hisF1 = (TH2F*)chains[ihis]->GetHistogram();
      chains[ihis]->Draw(Form("dMeas0-dMean0:T_DistRef_0.fP[4]>>hisMean_%d(20,-4,4,100,-0.25,0.25)",(ihis+1)*2000),cutOut+cutOutF,"colz",100000,0);
      hisM = (TH2F*)chains[ihis]->GetHistogram();
      hisM->Write(); hisA->Write();hisF0->Write(); hisF1->Write();
      ftrackFit->Flush();
    }
  }

  for (Int_t ihis=0; ihis<7; ihis++){
    printf("\n\nProcessing frames\t%d\nnn",(ihis+1)*2000);
    hisM = (TH2F*)ftrackFit->Get(Form("hisMean_%d",(ihis+1)*2000));
    hisA = (TH2F*)ftrackFit->Get(Form("hisAll_%d",(ihis+1)*2000));
    hisF0 = (TH2F*)ftrackFit->Get(Form("hisFrame0_%d",(ihis+1)*2000));
    hisF1 = (TH2F*)ftrackFit->Get(Form("hisFrame1_%d",(ihis+1)*2000));
    arrayHisto->AddLast(hisA);
    arrayHisto->AddLast(hisF0);
    arrayHisto->AddLast(hisF1);
    arrayHisto->AddLast(hisM);
  }
  delete ftrackFit;
  //
  // 3. Draw figures
  //
  gStyle->SetOptStat(0);
  TCanvas *canvasFit   = new TCanvas("canvasFitFrame","canvasFitframe",900,700);
  canvasFit->Divide(3,2,0,0);
  for (Int_t ihis=1; ihis<7; ihis++){
    //
    canvasFit->cd(ihis);
    char hname[10000];
    snprintf(hname,1000,"hisAll_%d",(ihis+1)*2000);
    hisA = (TH2F*)arrayHisto->FindObject(hname);
    snprintf(hname,1000,"hisFrame0_%d",(ihis+1)*2000);
    hisF0 = (TH2F*)arrayHisto->FindObject(hname);
    snprintf(hname,1000,"hisFrame1_%d",(ihis+1)*2000);
    hisF1 = (TH2F*)arrayHisto->FindObject(hname);
    snprintf(hname,1000,"hisMean_%d",(ihis+1)*2000);
    hisM = (TH2F*)arrayHisto->FindObject(hname);
    //
    //
    hisM->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    TH1 * hisRA= (TH1*)arrayFit.At(2)->Clone();
    hisF0->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    TH1 * hisRF0= (TH1*)arrayFit.At(2)->Clone();    
    hisF1->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    TH1 * hisRF1= (TH1*)arrayFit.At(2)->Clone();    
    //
    hisRA->SetMarkerStyle(20);
    hisRF0->SetMarkerStyle(21);
    hisRF1->SetMarkerStyle(21);
    hisRA->SetMarkerColor(1);
    hisRF0->SetMarkerColor(4);
    hisRF1->SetMarkerColor(2);
    TF1 * f1a= new TF1("f1a","pol1");
    TF1 * f1f0= new TF1("f1a0","pol1");
    TF1 * f1f1= new TF1("f1a1","pol1");
    f1a->SetLineColor(1);
    f1f0->SetLineColor(4);
    f1f1->SetLineColor(2);
    hisRA->Fit(f1a);
    hisRF0->Fit(f1f0);
    hisRF1->Fit(f1f1);
    hisRF1->SetMinimum(0);
    hisRF1->SetMaximum(0.05);
    // hisRA->Draw();
    hisRF1->GetXaxis()->SetTitle("q/p_{T} (1/GeV)");
    hisRF1->GetYaxis()->SetTitle("#sigma_{r#phi} (cm)");
    hisRF1->Draw("");   
    hisRF0->Draw("same");   
    TLegend * legend = new TLegend(0.11,0.11,0.65,0.25, Form("Track residual r#phi distortion: N_{ion}=%d",(ihis+1)*2000));
    legend->AddEntry(hisRF0,"a_{0}+a_{1}#rho");
    legend->AddEntry(hisRF1,"a_{0}+(a_{1}+a_{2}z+a_{3}z^2)#rho");
    legend->SetBorderSize(0);
    legend->Draw();
  }
  //
  canvasFit->SaveAs("canvasFrameFitRPhiVersion0.pdf");
  canvasFit->SaveAs("canvasFrameFitRPhiVersion0.png");
  //
}
