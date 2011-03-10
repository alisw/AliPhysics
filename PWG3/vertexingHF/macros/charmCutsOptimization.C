#include <fstream>
#include <Riostream.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TKey.h>
#include <TObjectTable.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TText.h>

#include <AliMultiDimVector.h>
#include "AliHFMassFitter.h"
#include <AliSignificanceCalculator.h>

#include <fstream>

//global variables
Double_t nsigma=3;
Int_t decCh=0;
Int_t fitbtype=0;
Int_t rebin=2;
Double_t sigma=0.012;
Int_t pdg;
Double_t mass;

Double_t sigmaCut=0.035;//0.03;
Double_t errSgnCut=0.4;//0.35;
Double_t nSigmaMeanCut=4.;//3.;


ofstream outcheck;
ofstream outdetail;

Bool_t Data(TH1F* h,Double_t* rangefit,Bool_t writefit,Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Double_t& sigmafit, Int_t& status);
Bool_t BinCounting(TH1F* h, Double_t* rangefit, Bool_t writefit, Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Int_t& status);
Bool_t MC(TH1F* hs,TH1F* hb, Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Double_t& sigmaused, Int_t& status);

//decCh:
//- 0 = kDplustoKpipi
//- 1 = kD0toKpi
//- 2 = kDstartoKpipi
//- 3 = kDstoKKpi
//- 4 = kD0toKpipipi
//- 5 = kLambdactopKpi

//Note: writefit=kTRUE writes the root files with the fit performed but it also draw all the canvas, so if your computer is not powerfull enough I suggest to run it in batch mode (root -b)

Bool_t charmCutsOptimization(Bool_t isData=kTRUE,TString part="both"/*"A" anti-particle, "P" particle*/,TString centr="no",Bool_t writefit=kFALSE,Int_t minentries=50,Double_t *rangefit=0x0, Bool_t useBinCounting=kTRUE){
  outcheck.open("output.dat",ios::out);
  outdetail.open("discarddetails.dat",ios::out);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  //~/Lavoro/PbPb/tagli/SIGAOD33/mar02/cent3070/
  TString filename="AnalysisResults.root",dirname="PWG3_D2H_Significance",listname="coutputSig",mdvlistname="coutputmv";

  TString hnamemass="hMass_",hnamesgn="hSig_",hnamebkg="hBkg_";


  switch (decCh) {
  case 0:
    listname+="Dplus";
    mdvlistname+="Dplus";
    pdg=411;
    break;
  case 1:
    listname+="D0";
    mdvlistname+="D0";
    pdg=421;
    break;
  case 2:
    listname+="Dstar";
    mdvlistname+="Dstar";
    pdg=413;
    break;
  case 3:
    listname+="Ds";
    mdvlistname+="Ds";
    pdg=431;
    break;
  case 4:
    listname+="D04";
    mdvlistname+="D04";
    pdg=421;
    break;
  case 5:
    listname+="Lc";
    mdvlistname+="Lc";
    pdg=4122;
    break;
  default:
    cout<<decCh<<" is not allowed as decay channel "<<endl;
    return kFALSE;
  }
  mass=TDatabasePDG::Instance()->GetParticle(pdg)->Mass();

  if(part!="both"){
    listname.Append(part);
    mdvlistname.Append(part);
  }
  if(centr!="no"){
    listname.Append(centr);
    mdvlistname.Append(centr);
  }
  cout<<"Mass = "<<mass<<endl;
  
  Int_t countFitFail=0,countSgnfFail=0,countNoHist=0,countBkgOnly=0;

  outcheck<<"ptbin\tmdvGlobAddr\thistIndex\tSignif\tS\tB"<<endl;
  outdetail<<"ptbin\tmdvGlobAddr\thistIndex\trelErrS\t\tmean_F-mass (mass = "<<mass<<")"<<endl;
  TFile *fin=new TFile(filename.Data());
  if(!fin->IsOpen()){
    cout<<"File "<<filename.Data()<<" not found"<<endl;
    return kFALSE;
  }

  TDirectoryFile *dir=(TDirectoryFile*)fin->GetDirectory(dirname);
  if(!dir){
    cout<<"Directory "<<dirname<<" not found"<<endl;
    return kFALSE;
  }

  TList* histlist= (TList*)dir->Get(listname);
  if(!histlist) {
    cout<<listname<<" doesn't exist"<<endl;
    return kFALSE;
  }

  TList* listamdv= (TList*)dir->Get(mdvlistname);
  if(!listamdv) {
    cout<<mdvlistname<<" doesn't exist"<<endl;
    return kFALSE;
  }

  TH1F* hstat=(TH1F*)histlist->FindObject("fHistNEvents");
  TCanvas *cst=new TCanvas("hstat","Summary of statistics");
  if(hstat) {
    cst->cd();
    cst->SetGrid();
    hstat->Draw("htext0");
    cst->SaveAs("hstat.png");
  }else{
    cout<<"Warning! fHistNEvents not found in "<<listname.Data()<<endl;
  }

  Bool_t isMC=kFALSE;
  TH1F* htestIsMC=(TH1F*)histlist->FindObject("hSig_0");
  if(htestIsMC) isMC=kTRUE;

  Int_t nptbins=listamdv->GetEntries();
  Int_t nhist=(histlist->GetEntries()-1);//-1 because of fHistNevents
  if(isMC) nhist/=4; ///4 because hMass_, hSgn_,hBkg_,hRfl_
  Int_t count=0;
  Int_t *indexes= new Int_t[nhist];
  //initialize indexes[i] to -1
  for(Int_t i=0;i<nhist;i++){
    indexes[i]=-1;
  }

  TFile* fout=new TFile(Form("outputSignifMaxim.root"),"recreate");

  TH1F** hSigma=new TH1F*[nptbins];
  TH1F* hstatus=new TH1F("hstatus","Flag status",6,-0.5,5.5);
  hstatus->GetXaxis()->SetBinLabel(1,"fit fail");
  hstatus->GetXaxis()->SetBinLabel(2,"fit ok and good results");
  hstatus->GetXaxis()->SetBinLabel(3,"quality requirements not satisfied");
  hstatus->GetXaxis()->SetBinLabel(4,"only bkg fit ok");
  hstatus->GetXaxis()->SetBinLabel(5,"negative signif");
  hstatus->GetXaxis()->SetBinLabel(6,Form("< %d entries",minentries));

  //Check wheter histograms are filled
  if(isData){
    for(Int_t i=0;i<nhist;i++){
      TString name=Form("%s%d",hnamemass.Data(),i);
      TH1F* h=(TH1F*)histlist->FindObject(name.Data());

      if(!h){
	cout<<name<<" not found"<<endl;
	continue;
      }

      if(h->GetEntries()>minentries){
	//cout<<"Entries = "<<h->GetEntries()<<endl;
	if (h->Integral() > minentries){
	  cout<<i<<") Integral = "<<h->Integral()<<endl;
	  indexes[i]=i;
	  count++;
	}
      }
    }
  

    cout<<"There are "<<count<<" histogram with more than "<<minentries<<" entries"<<endl;
    if(count==0) {
      cout<<"No histogram to draw..."<<endl;
      return kFALSE;
    }
  }
  //create multidimvectors

  //for(Int_t i=0;i<1;i++){
  for(Int_t i=0;i<nptbins;i++){

    //multidimvectors for signal
    AliMultiDimVector *mdvS=(AliMultiDimVector*)listamdv->FindObject(Form("multiDimVectorPtBin%d",i));
    TString name=mdvS->GetName(),nameErr="err",setname="";
    
    setname=Form("S%s",name.Data());
    mdvS->SetName(setname.Data());
    outcheck<<"\n"<<mdvS->GetPtLimit(0)<<" < Pt <"<<mdvS->GetPtLimit(1)<<endl;

    AliMultiDimVector *mdvSerr=(AliMultiDimVector*)mdvS->Clone(setname.Data());
    setname=Form("%sS%s",nameErr.Data(),name.Data());
    mdvSerr->SetName(setname.Data());

    //multidimvectors for background
    setname=Form("B%s",name.Data());
    AliMultiDimVector *mdvB=(AliMultiDimVector*)mdvS->Clone(setname.Data());

    AliMultiDimVector *mdvBerr=(AliMultiDimVector*)mdvS->Clone(setname.Data());
    setname=Form("%sB%s",nameErr.Data(),name.Data());
    mdvBerr->SetName(setname.Data());

    //multidimvectors for significance
    setname=Form("Sgf%s",name.Data());
    AliMultiDimVector *mdvSgnf=(AliMultiDimVector*)mdvS->Clone(setname.Data());

    AliMultiDimVector *mdvSgnferr=(AliMultiDimVector*)mdvS->Clone(setname.Data());
    setname=Form("%sSgf%s",nameErr.Data(),name.Data());
    mdvSgnferr->SetName(setname.Data());

    hSigma[i]=new TH1F(Form("hSigmapt%d",i),Form("Sigma distribution pt bin %d (%.1f < pt < %.1f)",i,mdvSgnf->GetPtLimit(0),mdvSgnf->GetPtLimit(1)), 200,0.,0.05);

    Int_t nhistforptbin=mdvS->GetNTotCells();
    //Int_t nvarsopt=mdvS->GetNVariables();
 
    cout<<"nhistforptbin = "<<nhistforptbin<<endl;

    //loop on all histograms and do AliHFMassFitter
    //for(Int_t ih=0;ih<1;ih++){
    for(Int_t ih=0;ih<nhistforptbin;ih++){
      printf("Analyzing indexes[%d] = %d \n",ih+i*nhistforptbin,indexes[ih+i*nhistforptbin]);
      Int_t status=-1;
      if(isData && indexes[ih+i*nhistforptbin] == -1) {
	status=5;
	mdvSgnferr->SetElement(ih,0);
	mdvS->SetElement(ih,0);
	mdvSerr->SetElement(ih,0);
	mdvB->SetElement(ih,0);
	mdvBerr->SetElement(ih,0);

	continue;
      }
      outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin];
      TString name;
      TH1F* h=0x0;
      TH1F* g=0x0;
      Double_t signif=0, signal=0, background=0, errSignif=0, errSignal=0, errBackground=0,sigmafit=0;

      if(isData){
	name=Form("%s%d",hnamemass.Data(),indexes[ih+i*nhistforptbin]);
	h=(TH1F*)histlist->FindObject(name.Data());
	if(!h)continue;
	if(useBinCounting) {
	  if (h->GetEntries() >= minentries)
	    BinCounting(h,rangefit,writefit,signal,errSignal,background,errBackground,signif,errSignif,status);
	} else 
	  Data(h,rangefit,writefit,signal,errSignal,background,errBackground,signif,errSignif,sigmafit,status);
      }else{
	name=Form("%s%d",hnamesgn.Data(),ih+i*nhistforptbin);
	h=(TH1F*)histlist->FindObject(name.Data());
	if(!h){
	  cout<<name.Data()<<" not found"<<endl;
	  continue;
	}
	name=Form("%s%d",hnamebkg.Data(),ih+i*nhistforptbin);
	g=(TH1F*)histlist->FindObject(name.Data());
	if(!g){
	  cout<<name.Data()<<" not found"<<endl;
	  continue;
	}
	MC(h,g,signal,errSignal,background,errBackground,signif,errSignif,sigmafit,status);
      }
      hstatus->Fill(status);

      if(status==1){
	mdvSgnf->SetElement(ih,signif);
	mdvSgnferr->SetElement(ih,errSignif);
	mdvS->SetElement(ih,signal);
	mdvSerr->SetElement(ih,errSignal);
	mdvB->SetElement(ih,background);
	mdvBerr->SetElement(ih,errBackground);
	hSigma[i]->Fill(sigmafit);
      
      }else{
	mdvSgnf->SetElement(ih,0);
	mdvSgnferr->SetElement(ih,0);
	mdvS->SetElement(ih,0);
	mdvSerr->SetElement(ih,0);
	mdvB->SetElement(ih,0);
	mdvBerr->SetElement(ih,0);
	if(status==3){
	  mdvB->SetElement(ih,background);
	  mdvBerr->SetElement(ih,errBackground);
	}

      }

    }
  
  

    fout->cd();
    mdvS->Write();
    mdvB->Write();
    mdvSgnf->Write();

    mdvSerr->Write();
    mdvBerr->Write();
    mdvSgnferr->Write();
    hSigma[i]->Write();
 
  }
  

  TCanvas *cinfo=new TCanvas("cinfo","Status");
  cinfo->cd();
  cinfo->SetGrid();
  hstatus->Draw("htext0");

  fout->cd();
  hstatus->Write();

  fout->Close();

  outcheck<<"\nSummary:\n - Total number of histograms: "<<nhist<<"\n - "<<count<<" histograms with more than "<<minentries<<" entries; \n - Too few entries in histo "<<countNoHist<<" times;\n - Fit failed "<<countFitFail<<" times \n - no sense Signal/Background/Significance "<<countSgnfFail<<" times\n - only background "<<countBkgOnly<<" times"<<endl;
  outcheck.close();
  return kTRUE;
}

//this function fit the hMass histograms
//status = 0 -> fit fail
//         1 -> fit ok and good results
//         2 -> quality requirements not satisfied, try to fit with bkg only
//         3 -> only bkg fit ok
//         4 -> negative signif
//         5 -> not enough entries in the hisotgram
Bool_t Data(TH1F* h,Double_t* rangefit,Bool_t writefit, Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Double_t& sigmafit, Int_t& status){
  Int_t nbin=h->GetNbinsX();
  Double_t min=h->GetBinLowEdge(7);
  Double_t max=h->GetBinLowEdge(nbin-5)+h->GetBinWidth(nbin-5);

  min=h->GetBinLowEdge(1);
  max=h->GetBinLowEdge(nbin+1);

  if(rangefit) {
    min=rangefit[0];
    max=rangefit[1];
  }

  AliHFMassFitter fitter(h,min, max,rebin,fitbtype);
  fitter.SetInitialGaussianMean(mass);
  fitter.SetInitialGaussianSigma(sigma);

  //if(ih==0) fitter.InitNtuParam(Form("ntuPtbin%d",i));
  // fitter.SetHisto(h);
  // fitter.SetRangeFit(min,max);
  //fitter.SetRangeFit(1.68,2.05);

  //fitter.SetType(fitbtype,0);

  Bool_t ok=fitter.MassFitter(kFALSE);
  if(!ok){
    cout<<"FIT NOT OK!"<<endl;
    //countBkgOnly++;
    //outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t 0\t xxx"<<"\t bkgonly"<<endl;
    outcheck<<"\t 0\t xxx"<<"\t failed"<<endl;
    status=0;
    return kFALSE;
  }else{ //fit ok!

    if(writefit) fitter.WriteCanvas(h->GetName(),"./",nsigma);
    fitter.Signal(nsigma,sgn,errsgn);
    fitter.Background(nsigma,bkg,errbkg);
    Double_t meanfit=fitter.GetMean();
    sigmafit=fitter.GetSigma();
	  

    //if(ok==kTRUE && ( (sigmafit < 0.03) || (sigmafit < 0.04 && mdvS->GetPtLimit(0)>8.)) && sgn > 0 && bkg > 0){
    if(ok==kTRUE && ( (sigmafit < sigmaCut) ) && sgn > 0 && bkg > 0){
      Double_t errmeanfit=fitter.GetMassFunc()->GetParError(fitter.GetNFinalPars()-2);
      fitter.Significance(nsigma,sgnf,errsgnf);
      if(sgnf >0){
	      
	if(errsgn/sgn < errSgnCut && /*TMath::Abs(meanfit-mass)<0.015*/TMath::Abs(meanfit-mass)/errmeanfit < nSigmaMeanCut){
	  //outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t"<<signif<<" +- "<<errSignif<<"\t"<<sgn<<" +- "<<errsgn<<"\t"<<bkg<<" +- "<<errbkg<<endl;
	  outcheck<<"\t\t "<<sgnf<<" +- "<<errsgnf<<"\t"<<sgn<<" +- "<<errsgn<<"\t"<<bkg<<" +- "<<errbkg<<endl;
	  status=1;
	      
	}else{
	  status=2;
	  //outdetail<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t"<<errsgn/sgn<<"\t\t "<<(meanfit-mass)/errmeanfit<<endl;
	  outdetail<<"\t\t "<<errsgn/sgn<<"\t\t "<<(meanfit-mass)/errmeanfit<<endl;
	  ok=fitter.RefitWithBkgOnly(kFALSE);
	  if (ok){
	    status=3;
	    //countBkgOnly++;
	    Double_t bkg=0,errbkg=0.;
	    Double_t nsigmarange[2]={mass-nsigma*sigma,mass+nsigma*sigma};
	    fitter.Background(nsigmarange[0],nsigmarange[1],bkg,errbkg); 
	    //outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t 0\t "<<bkg <<"\t bkgonly"<<endl;
	    outcheck<<"\t\t  0\t "<<bkg <<"\t bkgonly"<<endl;
	  }else{
	    //outdetail<<i<<"\t\t "<<ih<<"\t\tnot able to refit with bkg obly"<<endl;
	    outdetail<<"\t\t \t\tnot able to refit with bkg obly"<<endl;
	    status=0;
	    return kFALSE;
	  }
	}//only bkg
      }//check signif>0
      else{ 
	status=4;
	//countSgnfFail++;
	cout<<"Setting to 0 (fitter results meaningless)"<<endl;
	outcheck<<"\t S || B || sgnf negative";

	return kFALSE;
      } 
    } //end fit ok!
  }
  outcheck<<endl;
  return kTRUE; 
}



//this function counts the entries in hSgn and hBgk
Bool_t MC(TH1F* hs,TH1F* hb, Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Double_t& sigmaused, Int_t& status){

  //do we want to use a fixed sigma or take the standard deviation of the signal histogram?
  sigmaused=hs->GetRMS();
  //sigmaused=sigma;

  Double_t nsigmarange[2]={mass-nsigma*sigmaused,mass+nsigma*sigmaused}; 
  cout<<"from "<<nsigmarange[0]<<" to "<<nsigmarange[1]<<endl;

  Int_t binnsigmarange[2]={hs->FindBin(nsigmarange[0]),hs->FindBin(nsigmarange[1])};//for bkg histo it's the same
  cout<<"bins "<<binnsigmarange[0]<<" e "<<binnsigmarange[1]<<endl;

  sgn=hs->Integral(binnsigmarange[0],binnsigmarange[1]);
  errsgn=TMath::Sqrt(sgn);
  bkg=hb->Integral(binnsigmarange[0],binnsigmarange[1]);
  errbkg=TMath::Sqrt(bkg);
  if(sgn+bkg>0.) sgnf=sgn/TMath::Sqrt(sgn+bkg);
  else {
    status=2;
    return kFALSE;
  }
  errsgnf=TMath::Sqrt(sgnf*sgnf/(sgn+bkg)/(sgn+bkg)*(1/4.*errsgn*errsgn+errbkg*errbkg)+sgnf*sgnf/sgn/sgn*errsgn*errsgn);
  status=1;
  return kTRUE; 

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// par[0], par[1] expo params, par[2], par[3] exclusion range
Bool_t reject = true;
Double_t ExpoBkgWoPeak(Double_t *x, Double_t *par){

  if( reject && x[0]>par[2] && x[0]<par[3] ){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + TMath::Exp(par[1]*x[0]) ;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this function fit the hMass histograms
//status = 0 -> fit fail
//         1 -> fit ok and good results
//         2 -> negative signif
Bool_t BinCounting(TH1F* h,Double_t* rangefit,Bool_t writefit, Double_t& sgn, Double_t& errsgn, Double_t& bkg, Double_t& errbkg, Double_t& sgnf, Double_t& errsgnf, Int_t& status){

  Int_t nbin=h->GetNbinsX();
  Double_t min=h->GetBinLowEdge(1);
  Double_t max=h->GetBinLowEdge(nbin+1);

  if(rangefit) {
    min=rangefit[0];
    max=rangefit[1];
  }

  // Bkg fit : exponential = A*exp(B*x) 
  reject = true;
  TF1 *fBkgFit = new TF1("fBkgFit",ExpoBkgWoPeak,min,max,4);
  fBkgFit->FixParameter(2,mass-nsigma*sigma);
  fBkgFit->FixParameter(3,mass+nsigma*sigma);
  TFitResultPtr r = h->Fit(fBkgFit,"RS+");
  Int_t ok = r;

  reject = false;
  TF1 *fBkgFct = new TF1("fBkgFct",ExpoBkgWoPeak,min,max,4);
  fBkgFct->SetLineStyle(2);
  for(Int_t i=0; i<4; i++) fBkgFct->SetParameter(i,fBkgFit->GetParameter(i));
  h->GetListOfFunctions()->Add(fBkgFct);
  TH1F * hBkgFct = (TH1F*)fBkgFct->GetHistogram();

  if(ok==-1){
    cout<<"FIT NOT OK!"<<endl;
    cout<<"\t 0\t xxx"<<"\t failed"<<endl;
    status=0;
    return kFALSE;
  } 
  else { //fit ok!
    status = 1;    
    Double_t binStartCount = h->FindBin(mass-nsigma*sigma);
    Double_t binEndCount = h->FindBin(mass+nsigma*sigma);
    Double_t counts=0., bkgcounts=0., errcounts=0., errbkgcounts=0.;
    for (Int_t ibin = binStartCount; ibin<=binEndCount; ibin++) {
      counts += h->GetBinContent( ibin );
      errcounts += counts ;
      Double_t center =  h->GetBinCenter(ibin);
      bkgcounts += hBkgFct->GetBinContent( hBkgFct->FindBin(center) );
      errbkgcounts += bkgcounts ;
    }
    bkg = bkgcounts;
    errbkg = TMath::Sqrt( errbkgcounts );
    sgn = counts - bkg ;
    if(sgn<0) status = 2; // significance < 0
    errsgn = TMath::Sqrt( counts + errbkg*errbkg );
    sgnf = sgn / TMath::Sqrt( sgn + bkg );
    errsgnf = TMath::Sqrt( sgnf*sgnf/(sgn+bkg)/(sgn+bkg)*(1/4.*errsgn*errsgn+errbkg*errbkg) + sgnf*sgnf/sgn/sgn*errsgn*errsgn );
    //    cout << " Signal "<<sgn<<" +- "<<errsgn<<", bkg "<<bkg<<" +- "<<errbkg<<", significance "<<sgnf<<" +- "<<errsgnf<<endl;

    if(writefit) {
      TString filename = Form("%sMassFit.root",h->GetName());
      TFile* outputcv = new TFile(filename.Data(),"recreate");      
      TCanvas* c = new TCanvas();
      c->SetName(Form("%s",h->GetName()));
      h->Draw();
      TPaveText *pavetext=new TPaveText(0.4,0.7,0.85,0.9,"NDC");     
      pavetext->SetBorderSize(0);
      pavetext->SetFillStyle(0);
      pavetext->AddText(Form("Signal = %4.2e #pm %4.2e",sgn,errsgn));
      pavetext->AddText(Form("Bkg = %4.2e #pm %4.2e",bkg,errbkg));
      pavetext->AddText(Form("Signif = %3.2f #pm %3.2f",sgnf,errsgnf));
      c->cd();
      pavetext->DrawClone();
      outputcv->cd();
      c->Write();
      outputcv->Close();
      delete outputcv;
      delete c;
    }

  }
  
  delete fBkgFit;
  delete fBkgFct;

  return kTRUE; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// which=0 plot significance
//      =1 plot signal
//      =2 plot background
//      =3 plot S/B
// maximize = kTRUE (default) if you want to fix the step of the variables not shown to the value that maximize the significance. Note that these values are saved in fixedvars.dat
// readfromfile = kTRUE (default is kFALSE) if you want to read the value fixed in a previous run of this function (e.g. significance or signal maximization)


void showMultiDimVector(Int_t n=2,Int_t which=0, Bool_t plotErrors=kFALSE,Bool_t readfromfile=kFALSE, Bool_t fixedrange=kFALSE, Bool_t fixedplane=kFALSE){

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetFrameBorderMode(0);

  TFile* fin=new TFile("outputSignifMaxim.root");
  if(!fin->IsOpen()){
    cout<<"outputSignifMaxim.root not found"<<endl;
    return;
  }

  if(n>2){
    cout<<"Error! cannot show "<<n+1<<" dimentions"<<endl;
    return;
  }

  TString name,title,namebis,shorttitle;
  switch (which){
  case 0:
    name="SgfmultiDimVectorPtBin";
    title="Significance";
    shorttitle="Sgnf";
    break;
  case 1:
    name="SmultiDimVectorPtBin";
    title="Signal";
    shorttitle="S";
    break;
  case 2:
    name="BmultiDimVectorPtBin";
    title="Background";
    shorttitle="B";
    break;
  case 3:
    name="SmultiDimVectorPtBin";
    namebis="BmultiDimVectorPtBin";
    title="Signal over Background ";
    shorttitle="SoB";
    break;
  // case 4:
  //   name="errBmultiDimVectorPtBin";
  //   title="Background (error)";
  //   break;
  }
 
  if(plotErrors && which!=3 && n==2){
    name.Prepend("err");
    title.Append(" Error") ;
    shorttitle.Append("Err");
  }

  Int_t nptbins=50;
  
  for(Int_t ip=0;ip<=nptbins;ip++){
    TString mdvname=Form("%s%d",name.Data(),ip);
    AliMultiDimVector* mdv=(AliMultiDimVector*)fin->Get(mdvname);
    if(!mdv){
      nptbins=ip;
      cout<<"Number of pt bins "<<ip<<endl;
      break;
    }
  }

  cout<<"Projecting "<<title.Data()<<" with respect to the maximization variable(s) [chose]"<<endl;
 
  Int_t variable[2]; //no more than 2D
  TString mdvname=Form("%s0",name.Data()), mdverrname="";//, mdvnamebis="", mdverrnamebis="";
  AliMultiDimVector* mdv=(AliMultiDimVector*)fin->Get(mdvname);
  AliMultiDimVector* mdverr=0x0;
  if(!mdv){
    cout<<mdvname.Data()<<" not found"<<endl;
    return;
  }
  
  Int_t nvarsopt=mdv->GetNVariables();
  //Int_t nfixed=nvarsopt-n;
 
  Int_t fixedvars[nvarsopt];
  Int_t allfixedvars[nvarsopt*nptbins];
 


  fstream writefixedvars;
  if(readfromfile) {
    //open file in read mode
    writefixedvars.open("fixedvars.dat",ios::in);
    Int_t longi=0;
    while(writefixedvars){
      writefixedvars>>allfixedvars[longi];
      longi++;
    }
  }
  else {
    //open file in write mode
    writefixedvars.open("fixedvars.dat",ios::out);
  }

  Bool_t freevars[nvarsopt];

  //ask variables for projection
  for(Int_t k=0;k<nvarsopt;k++){
    cout<<k<<" "<<mdv->GetAxisTitle(k)<<endl;
    freevars[k]=kTRUE;
  }
  cout<<"Choose "<<n<<" variable(s)"<<endl;
  for(Int_t j=0;j<n;j++){
    cout<<"var"<<j<<": ";
    cin>>variable[j];
    freevars[variable[j]]=kFALSE;
  }
  if(n==1) variable[1]=999;

  TCanvas* cvpj=new TCanvas(Form("proj%d",variable[0]),Form("%s wrt %s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data()));

  TMultiGraph* mg=new TMultiGraph(Form("proj%d",variable[0]),Form("%s wrt %s;%s;%s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),(mdv->GetAxisTitle(variable[0])).Data(),title.Data()));
  TLegend *leg=new TLegend(0.7,0.2,0.9,0.6,"Pt Bin");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  //set scale
  fstream writerange;
  Float_t axisrange[2*nptbins];
  if(fixedrange) {
    //open file in read mode
    writerange.open("rangehistos.dat",ios::in);
    Int_t longi=0;
    while(writerange){
      writerange>>axisrange[longi];
      longi++;
    }
  }
  else {
    //open file in write mode
    writerange.open("rangehistos.dat",ios::out);
  }

  for (Int_t i=0;i<nptbins;i++){   //loop on ptbins
    cout<<"\nPtBin = "<<i<<endl;

    //using AliSignificanceCalculator

    TString nameS,nameB,nameerrS,nameerrB;
    nameS.Form("SmultiDimVectorPtBin%d",i);
    nameerrS.Form("errSmultiDimVectorPtBin%d",i);
    nameB.Form("BmultiDimVectorPtBin%d",i);
    nameerrB.Form("errBmultiDimVectorPtBin%d",i);
 
    AliMultiDimVector* mdvS=(AliMultiDimVector*)fin->Get(nameS.Data());
    AliMultiDimVector* mdvB=(AliMultiDimVector*)fin->Get(nameB.Data());
    AliMultiDimVector* mdvBerr=(AliMultiDimVector*)fin->Get(nameerrS.Data());
    AliMultiDimVector* mdvSerr=(AliMultiDimVector*)fin->Get(nameerrB.Data());
    if(!(mdvS && mdvB && mdvSerr && mdvBerr)){
      cout<<"one of the multidimvector is not present"<<endl;
      return;
    }

    AliSignificanceCalculator *cal=new AliSignificanceCalculator(mdvS,mdvB,mdvSerr,mdvBerr,1.,1.);

    AliMultiDimVector* mvess=cal->GetSignificanceError();
    AliMultiDimVector* mvpur=cal->CalculatePurity();
    AliMultiDimVector* mvepur=cal->CalculatePurityError();

    Int_t ncuts=mdvS->GetNVariables();
    Int_t *maxInd=new Int_t[ncuts];
    Float_t *cutvalues=new Float_t[ncuts];
    //init
    // for(Int_t ind=0;ind<ncuts;ind++)maxInd[ind]=0;

    Float_t sigMax0=cal->GetMaxSignificance(maxInd,0); //look better into this!!

    for(Int_t ic=0;ic<ncuts;ic++){
      cutvalues[ic]=((AliMultiDimVector*)fin->Get(nameS.Data()))->GetCutValue(ic,maxInd[ic]);

      //setting step of fixed variables
      if(readfromfile){ //from file
	fixedvars[ic]=allfixedvars[i+ic];
      }

      if(!readfromfile) { //using the values which maximize the significance
       	fixedvars[ic]=maxInd[ic];
       	//write to output fixedvars.dat
       	writefixedvars<<fixedvars[ic]<<"\t";
      }
    }
    //output file: return after each pt bin
    if(!readfromfile) writefixedvars<<endl;

    printf("Maximum of significance for Ptbin %d found in bin:\n",i);
    for(Int_t ic=0;ic<ncuts;ic++)printf(" %d  ",maxInd[ic]);
    printf("\ncorresponding to cut:\n");
    for(Int_t ic=0;ic<ncuts;ic++)printf(" %f  ",cutvalues[ic]);
    
    printf("\nSignificance = %f +- %f\n",sigMax0,mvess->GetElement(maxInd,0));
    printf("Purity       = %f +- %f\n",mvpur->GetElement(maxInd,0),mvepur->GetElement(maxInd,i));
    
    if(which==3){
      //mdv=0x0;
      mdv=cal->CalculateSOverB();
      if(!mdv)cout<<mdv->GetName()<<" null"<<endl;
      //mdverr=0x0;
      mdverr=cal->CalculateSOverBError();
      if(!mdverr)cout<<mdverr->GetName()<<" null"<<endl;
    }else{

      //multidimvector
      mdvname=Form("%s%d",name.Data(),i);   
      mdv=(AliMultiDimVector*)fin->Get(mdvname);
      if(!mdv)cout<<mdvname.Data()<<" not found"<<endl;
      
      //multidimvector of errors
      mdverrname=Form("err%s%d",name.Data(),i);   
      mdverr=(AliMultiDimVector*)fin->Get(mdverrname);
      if(!mdverr)cout<<mdverrname.Data()<<" not found"<<endl;
    }
    printf("Global Address %d (%d)\n",(Int_t)mdv->GetGlobalAddressFromIndices(maxInd,0),(Int_t)mdv->GetNTotCells()*i+(Int_t)mdv->GetGlobalAddressFromIndices(maxInd,0));
    TString ptbinrange=Form("%.1f < p_{t} < %.1f GeV/c",mdv->GetPtLimit(0),mdv->GetPtLimit(1));

    Float_t maxval=0;

    if(n==2) {
      gStyle->SetPalette(1);
      Int_t steps[2];
      Int_t nstep[2]={mdv->GetNCutSteps(variable[0]),mdv->GetNCutSteps(variable[1])}; 
  
      TH2F* hproj=new TH2F(Form("hproj%d",i),Form("%s wrt %s vs %s (Ptbin%d %.1f<pt<%.1f);%s;%s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data(),i,mdv->GetPtLimit(0),mdv->GetPtLimit(1),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data()),nstep[0],mdv->GetMinLimit(variable[0]),mdv->GetMaxLimit(variable[0]),nstep[1],mdv->GetMinLimit(variable[1]),mdv->GetMaxLimit(variable[1]));
      if(fixedplane){
	hproj=mdv->Project(variable[0],variable[1],fixedvars,0);
	hproj->SetTitle(Form("%s wrt %s vs %s (Ptbin%d %.1f<pt<%.1f);%s;%s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data(),i,mdv->GetPtLimit(0),mdv->GetPtLimit(1),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data()));
      }else{
	for(Int_t ist1=0;ist1<nstep[0];ist1++){
	  steps[0]=ist1;
	  Int_t fillbin1=ist1+1;
	  if(mdv->GetCutValue(variable[0],0)>mdv->GetCutValue(variable[0],mdv->GetNCutSteps(variable[0])-1))fillbin1=nstep[0]-ist1;
	  for(Int_t ist2=0;ist2<nstep[1];ist2++){
	    steps[1]=ist2;
	    Int_t fillbin2=ist2+1;
	    if(mdv->GetCutValue(variable[1],0)>mdv->GetCutValue(variable[1],mdv->GetNCutSteps(variable[1])-1))fillbin2=nstep[1]-ist2;
	    Int_t* varmaxim=mdv->FindLocalMaximum(maxval,variable,steps,n,0);
	    hproj->SetBinContent(fillbin1,fillbin2,maxval);
	    delete varmaxim;
	  }
	}
      }
      if(fixedrange) {
	hproj->SetMinimum(axisrange[2*i]);
	hproj->SetMaximum(axisrange[2*i+1]);
      } else{
	writerange<<hproj->GetMinimum()<<"\t"<<hproj->GetMinimum()<<endl;
      }
      TCanvas* cvpj=new TCanvas(Form("proj%d%dpt%d",variable[0],variable[1],i),Form("%s wrt %s vs %s (Ptbin%d)",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data(),i));
      cvpj->cd();
      hproj->DrawClone("COLZtext");
      cvpj->SaveAs(Form("%s%s.png",shorttitle.Data(),cvpj->GetName()));
      delete hproj;
    }

    if(n==1){

      Int_t nbins=mdv->GetNCutSteps(variable[0]);
 
      Double_t *x=new Double_t[nbins];
      Double_t *y=new Double_t[nbins];
      Double_t *errx=new Double_t[nbins];
      Double_t *erry=new Double_t[nbins];

      for(Int_t k=0;k<nbins;k++){ //loop on the steps (that is the bins of the graph)
	//init
	x[k]=0;y[k]=0;
	errx[k]=0;erry[k]=0;
	
	x[k]=mdv->GetCutValue(variable[0],k);
	errx[k]=mdv->GetCutStep(variable[0])/2.;
	Int_t onevariable[1]={variable[0]};
	Int_t onestep[1]={k};
	ULong64_t gladd;

	Float_t maxval;
	Int_t* varmaxim=mdv->FindLocalMaximum(maxval,onevariable,onestep,n,0);
	y[k]=maxval;
	gladd=mdv->GetGlobalAddressFromIndices(varmaxim,0);
	cout<<gladd<<endl;
	delete varmaxim;
      
	erry[k]=mdverr->GetElement(gladd);
	
	cout<<mdv->GetAxisTitle(variable[0])<<" step "<<k<<" = "<<x[k]<<":"<<" y = "<<y[k]<<endl;
      }
      
      cout<<"----------------------------------------------------------"<<endl;
      TGraphErrors* gr=new TGraphErrors(nbins,x,y,errx,erry);
      gr->SetMarkerStyle(20+i);
      gr->SetMarkerColor(i+1);
      gr->SetLineColor(i+1);
      if(i>=9){
	gr->SetMarkerColor(i+2);
	gr->SetLineColor(i+2);
      }
      gr->SetMinimum(0);

      gr->SetName(Form("g1%d",i));
      mg->Add(gr,"P");
      leg->AddEntry(gr,ptbinrange.Data(),"p");
    }
  }
   
  if(n==1){
    cvpj->cd();
    mg->Draw("A");
    leg->Draw();
    cvpj->SaveAs(Form("%s%s.png",shorttitle.Data(),cvpj->GetName()));
  } else delete cvpj;
}

//draw sigma as a function of cuts

void DrawSigmas(TH2F* h2cuts){
  TFile *fin=0x0;
  TString fittype="ExpFit";
  Int_t ntot=5;
  if(fittype=="Pol2Fit") ntot=6;
  Int_t ihfirst=0,ihlast=1; //change this (must think on it and remember what I wanted to do!)
  for(Int_t ih=ihfirst;ih<ihlast;ih++){
    fin=new TFile(Form("h%d%s.root",ih,fittype.Data()));
    if(!fin) continue;
    TCanvas *cv=(TCanvas*)fin->Get(Form("cv1%s%d",fittype.Data(),ih));
    TF1* func=(TF1*)cv->FindObject("funcmass");
    Int_t sigma=func->GetParameter(ntot-1);
    //h2cuts->SetBinContent();
    //to be finished
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Some methods to get the index of the histogram corresponding to a set of cuts

// vct = the AliMultiDimVector correponding to ptbin: it can be from AnalysisResults.root or outputSignifMaxim.root
// ptbin = pt bin
// indices = array of the index of the cut variables (the dimension of the array must be equal to the number of variables maximized)

Int_t GetNHistFromIndices(AliMultiDimVector* vct,Int_t ptbin,Int_t* indices){
  cout<<"Calculating the index of the histogram corresponding to the following cut steps:"<<endl;
  for(Int_t i=0;i<vct->GetNVariables();i++){
    cout<<vct->GetAxisTitle(i)<<" --> "<<indices[i]<<endl;
  }
  cout<<"Pt bin "<<ptbin<<" from "<<vct->GetPtLimit(0)<<" to "<<vct->GetPtLimit(1)<<endl;
  cout<<"Info: total number of cells per multidim: "<<vct->GetNTotCells()<<endl;
  ULong64_t glindex=vct->GetGlobalAddressFromIndices(indices,0);
  cout<<"The histogram you want is:\t";
  return glindex+vct->GetNTotCells()*ptbin;
}

// vct = the AliMultiDimVector correponding to ptbin: it can be from AnalysisResults.root or outputSignifMaxim.root
// ptbin = pt bin
// values = array of the cut values (the dimension of the array must be equal to the number of variables maximized)

Int_t GetNHistFromValues(AliMultiDimVector* vct,Int_t ptbin,Float_t* values){
  cout<<"Calculating the index of the histogram corresponding to the following cut values:"<<endl;
  for(Int_t i=0;i<vct->GetNVariables();i++){
    cout<<vct->GetAxisTitle(i)<<" --> "<<values[i]<<endl;
  }
  cout<<"Pt bin "<<ptbin<<" from "<<vct->GetPtLimit(0)<<" to "<<vct->GetPtLimit(1)<<endl;
  cout<<"Info: total number of cells per multidim: "<<vct->GetNTotCells()<<endl;
  ULong64_t glindex=vct->GetGlobalAddressFromValues(values,0);
 
  cout<<"The histogram you want is:\t"<<glindex+vct->GetNTotCells()*ptbin<<endl;
  return glindex+vct->GetNTotCells()*ptbin;
}

// vct = the AliMultiDimVector correponding to ptbin: it can be from AnalysisResults.root or outputSignifMaxim.root
// ptbin = pt bin
// values = array of the cut values: the dimention can be <= number of variables maximized
// valsgiven = array of dimention = to the number of variables optimized. For each variable put kTRUE if the value is given (in values), kFALSE otherwise
// nhistinrange =  pass an integer which will contains the number of histogram returned (that is the dimention of the Int_t* returned)

//NB: Remember that the cut applied is the lower edge of the step where lower=looser

Int_t* GetRangeHistFromValues(AliMultiDimVector* vct,Int_t ptbin,Bool_t* valsgiven,Float_t* values,Int_t& nhistinrange){
  Int_t nvargiven=0;
  nhistinrange=1;

  Int_t nvar4opt=vct->GetNVariables();
  Float_t allvals[nvar4opt];
 
  for (Int_t i=0;i<nvar4opt;i++) {
    if(valsgiven[i]==kTRUE) {
      allvals[i]=values[nvargiven];
      nvargiven++;
    }
    else {
      nhistinrange+=vct->GetNCutSteps(i);
      allvals[i]=vct->GetCutValue(i,0);
      //allvals[i]=vct->GetCutValue(0,i); //ivar,icell
    }
  }
  cout<<nhistinrange<<" index will be returned"<<endl;
  Int_t *rangeofhistos=new Int_t[nhistinrange];

  if(nhistinrange==1){
    rangeofhistos[0]=GetNHistFromValues(vct,ptbin,allvals);
    cout<<"output"<<rangeofhistos[0]<<endl;
  }else{
    Int_t index[nvar4opt-nvargiven];
    Int_t k=0;
    for (Int_t i=0;i<nvar4opt;i++){
      if(valsgiven[i]==kFALSE) {
	//cout<<"kTRUE==>"<<i<<endl;
	index[k]=i;
	k++;
      }
    }

    for(Int_t i=0;i<nvar4opt-nvargiven;i++){ //loop on number of free variables
      cout<<"Info: incrementing "<<vct->GetAxisTitle(index[i])<<endl;
      for(Int_t j=0;j<vct->GetNCutSteps(i);j++){ //loop on steps of each free variable
	allvals[index[i]]=vct->GetCutValue(index[i],j);
	rangeofhistos[i*vct->GetNCutSteps(i)+j]=GetNHistFromValues(vct,ptbin,allvals);
      }
    }
  }
  return rangeofhistos;
}

// vct = the AliMultiDimVector correponding to ptbin: it can be from AnalysisResults.root or outputSignifMaxim.root
// ptbin = pt bin
// nhist = number of the histogram from which you want to have the cut values (returned)

Float_t* GetCutValuesFromNHist(AliMultiDimVector* vct,Int_t ptbin,Int_t nhist){
  ULong64_t totCells=vct->GetNTotCells();
  ULong64_t globadd=nhist-ptbin*totCells;
  const Int_t nvars=vct->GetNVariables();
  Float_t* cuts=new Float_t[nvars];
  Int_t ptinside;
  vct->GetCutValuesFromGlobalAddress(globadd,cuts,ptinside);
  return cuts;
}

// ptbin = pt bin
// values = array of the cut values: the dimention can be <= number of variables maximized
// valsgiven = array of dimention = to the number of variables optimized. For each variable put kTRUE if the value is given (in values), kFALSE otherwise
// 

void DrawPossibilities(Int_t ptbin,Bool_t* valsgiven,Float_t* values,TString path="./",Int_t decCh=1){
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptStat(0);

  Int_t nhists;
  TString filename="AnalysisResults.root";
  TString dirname="PWG3_D2H_Significance",listname="coutputSig",mdvlistname="coutputmv";
  TString centr="020";

  TFile *fin=new TFile(Form("%s%s",path.Data(),filename.Data()));
  if(!fin){
    cout<<path.Data()<<filename.Data()<<" not found"<<endl;
    return;
  }
  TDirectoryFile *dir=(TDirectoryFile*)fin->GetDirectory(dirname);
  if(!dir){
    cout<<"Directory "<<dirname<<" not found"<<endl;
    return;
  }
  switch (decCh) {
  case 0:
    listname+="Dplus";
    mdvlistname+="Dplus";
    
    break;
  case 1:
    listname+="D0";
    mdvlistname+="D0";
    
    break;
  case 2:
    listname+="Dstar";
    mdvlistname+="Dstar";
    
    break;
  case 3:
    listname+="Ds";
    mdvlistname+="Ds";
    
    break;
  case 4:
    listname+="D04";
    mdvlistname+="D04";
    
    break;
  case 5:
    listname+="Lc";
    mdvlistname+="Lc";
    
    break;
  default:
    cout<<decCh<<" is not allowed as decay channel "<<endl;
    return;
  }
  mdvlistname+=centr;
  listname+=centr;

  TList* listamdv= (TList*)dir->Get(mdvlistname);
  if(!listamdv) {
    cout<<mdvlistname<<" doesn't exist"<<endl;
    return;
  }

  AliMultiDimVector* vct=(AliMultiDimVector*)listamdv->FindObject(Form("multiDimVectorPtBin%d",ptbin));

  TFile* fin2;
  TString filehistname="";
  Int_t* indexes=GetRangeHistFromValues(vct,ptbin,valsgiven,values,nhists);
  for(Int_t i=0;i<nhists;i++){

    fin2=new TFile(Form("%sh%dExpMassFit.root",path.Data(),indexes[i]));
    if(!fin2){
      cout<<"File "<<indexes[i]<<" not found!"<<endl;
      continue;
    }else{
      TCanvas *cv=(TCanvas*)fin2->Get(Form("c1h%dExp",indexes[i]));
      if(!cv){
	cout<<"Canvas c1h"<<indexes[i]<<"Exp not found among";
	fin2->ls();
	continue;
      }else{
	cv->Draw();
	cv->SaveAs(Form("h%dExpMassFit.png",indexes[i]));
	fin2=0x0;
      }
    }
  }
}

void Merge2Bins(Int_t b1, Int_t b2,TString pathin="./",Int_t decCh=1,TString part="both"/*"A" anti-particle, "P" particle*/){

  if(b2!=b1+1){
    printf("The bins to be merget must be consecutive. Check! [b1 = %d, b2= %d]\n",b1,b2);
    return;
  }

  TFile *fin=new TFile(Form("%sAnalysisResults.root",pathin.Data()));
  if (!fin){
    cout<<"Input file not found"<<endl;
    return;
  }

  TString dirname="PWG3_D2H_Significance",listname="coutputSig",mdvlistname="coutputmv";

  switch (decCh) {
  case 0:
    listname+="Dplus";
    mdvlistname+="Dplus";
    break;
  case 1:
    listname+="D0";
    mdvlistname+="D0";
    break;
  case 2:
    listname+="Dstar";
    mdvlistname+="Dstar";
    break;
  case 3:
    listname+="Ds";
    mdvlistname+="Ds";
    break;
  case 4:
    listname+="D04";
    mdvlistname+="D04";
    break;
  case 5:
    listname+="Lc";
    mdvlistname+="Lc";
    break;
  default:
    cout<<decCh<<" is not allowed as decay channel "<<endl;
    return;
  }

  if(part!="both"){
    listname.Append(part);
    mdvlistname.Append(part);
  }

  TDirectoryFile *dir=(TDirectoryFile*)fin->GetDirectory(dirname);
  if(!dir){
    cout<<"Directory "<<dirname<<" not found"<<endl;
    return;
  }

  TList* histlist= (TList*)dir->Get(listname);
  if(!histlist) {
    cout<<listname<<" doesn't exist"<<endl;
    return;
  }

  TList* listamdv= (TList*)dir->Get(mdvlistname);
  if(!listamdv) {
    cout<<mdvlistname<<" doesn't exist"<<endl;
    return;
  }
  if (!gSystem->AccessPathName(Form("merged%d%d",b1,b2))) gSystem->Exec(Form("mkdir merged%d%d",b1,b2));
  gSystem->Exec(Form("cd merged%d%d",b1,b2));

  TFile* fout=new TFile("mergeAnalysisResults.root","recreate");

  fout->mkdir(dirname);
  TList* listmdvout=new TList();
  listmdvout->SetName(listamdv->GetName());
  listmdvout->SetOwner();
  //listmdvout->SetTitle(listamdv->GetTitle());
  TList* histlistout=new TList();
  histlistout->SetName(histlist->GetName());
  histlistout->SetOwner();
  //histlistout->SetTitle(histlist->GetTitle());

  AliMultiDimVector* mdvin1=(AliMultiDimVector*)listamdv->FindObject(Form("multiDimVectorPtBin%d",b1));
  AliMultiDimVector* mdvin2=(AliMultiDimVector*)listamdv->FindObject(Form("multiDimVectorPtBin%d",b2));

  Int_t ntotHperbin = mdvin1->GetNTotCells();
  if(mdvin2->GetNTotCells() != ntotHperbin) {
    cout<<"Error! Number of histos in pt bin "<<b1<<" = "<<ntotHperbin<<" != Number of histos in pt bin "<<b2<<" = "<<mdvin2->GetNTotCells()<<endl;
    return;
  }
  Int_t nvar1=mdvin1->GetNVariables();
  if(nvar1 != mdvin2->GetNVariables()){
    cout<<"Error! Mismatch in number of variables"<<endl;
    return;
  }
  
  Float_t newptbins[2]={mdvin1->GetPtLimit(0),mdvin2->GetPtLimit(1)};
  Float_t loosercuts[nvar1], tightercuts[nvar1];
  TString axistitles[nvar1];
  Int_t ncells[nvar1];

  for (Int_t ivar=0;ivar<nvar1;ivar++){
    loosercuts[ivar]=mdvin1->GetCutValue(ivar,0);
    if(loosercuts[ivar] - mdvin2->GetCutValue(ivar,0) < 1e-8) printf("Warning! The loose cut %s is different between the 2: %f and %f\n",mdvin1->GetAxisTitle(ivar).Data(),loosercuts[ivar],mdvin2->GetCutValue(ivar,0));
    tightercuts[ivar]=mdvin1->GetCutValue(ivar,mdvin1->GetNCutSteps(ivar)-1);
    if(tightercuts[ivar] - mdvin2->GetCutValue(ivar,mdvin1->GetNCutSteps(ivar)-1) < 1e-8) printf("Warning! The tight cut %s is different between the 2: %f and %f\n",mdvin1->GetAxisTitle(ivar).Data(),tightercuts[ivar],mdvin2->GetCutValue(ivar,mdvin2->GetNCutSteps(ivar)));
    axistitles[ivar]=mdvin1->GetAxisTitle(ivar);
    cout<<axistitles[ivar]<<"\t";
    ncells[ivar]=mdvin1->GetNCutSteps(ivar);
  }

  AliMultiDimVector* mdvout= new AliMultiDimVector(Form("multiDimVectorPtBin%d",b1),"MultiDimVector",1,newptbins,mdvin1->GetNVariables(),ncells,loosercuts,tightercuts,axistitles);
  cout<<"Info: writing mdv"<<endl;
  listmdvout->Add(mdvout);

  Bool_t isMC=kFALSE;
  TH1F* htestIsMC=(TH1F*)histlist->FindObject("hSgn_0");
  if(htestIsMC) isMC=kTRUE;
  Int_t nptbins=listamdv->GetEntries();
  Int_t nhist=(histlist->GetEntries()-1);//-1 because of fHistNevents
  if(isMC) nhist/=4; ///4 because hMass_, hSgn_,hBkg_,hRfl_

 cout<<"Merging bin from "<<mdvin1->GetPtLimit(0)<<" to "<<mdvin1->GetPtLimit(1)<<" and from "<<mdvin2->GetPtLimit(0)<<" to "<<mdvin2->GetPtLimit(1)<<endl;
 Int_t firsth1=b1*ntotHperbin,firsth2=b2*ntotHperbin; //firsth2 = (b1+1)*ntotHperbin
 Int_t lasth1=firsth1+ntotHperbin-1,lasth2=firsth2+ntotHperbin-1;
 cout<<"Histograms from "<<firsth1<<" to "<<lasth1<<" and "<<firsth2<<" to "<<lasth2<<endl;

 //add the others mdv to the list
 Int_t cnt=0;
 for(Int_t i=0;i<nptbins;i++){
   if(i!=b1 && i!=b2){
     AliMultiDimVector* vcttmp=(AliMultiDimVector*)listamdv->FindObject(Form("multiDimVectorPtBin%d",i));
     if(i>b2) {
       vcttmp->SetName(Form("multiDimVectorPtBin%d",b2+cnt));
       cnt++;
     }
     listmdvout->Add(vcttmp);
   }
 }

 histlistout->Add((TH1F*)histlist->FindObject("fHistNEvents"));

 Int_t ih2=firsth2;

 for(Int_t ih1=firsth1;ih1<lasth1;ih1++){
   TH1F* h1=(TH1F*)histlist->FindObject(Form("hMass_%d",ih1));
   if(!h1){
     cout<<"hMass_"<<ih1<<" not found!"<<endl;
     continue;
   }
   TH1F* h2=(TH1F*)histlist->FindObject(Form("hMass_%d",ih2));
   if(!h2){
     cout<<"hMass_"<<ih2<<" not found!"<<endl;
     continue;
   }
   //h1->SetName(Form("hMass_%d",cnt));
   h1->Add(h2);
   histlistout->Add(h1);
   ih2++;
   //cnt++;
   h1=0x0;
 }

 cnt=0;
 for(Int_t j=0;j<ntotHperbin*nptbins;j++){
   if(!(j>=firsth1 && j<lasth2)){
     TH1F* htmp=(TH1F*)histlist->FindObject(Form("hMass_%d",j));
     if(j>=lasth2){
       //cout<<lasth1<<" + "<<cnt<<endl;
       htmp->SetName(Form("hMass_%d",lasth1+cnt));
       cnt++;
     }
     histlistout->Add(htmp);
   }
 }

 fout->cd();
 ((TDirectoryFile*)fout->Get(dirname))->cd();
 listmdvout->Write(mdvlistname.Data(),TObject::kSingleKey);
 histlistout->Write(listname.Data(),TObject::kSingleKey);
 fout->Close();
}

void SubtractBkg(Int_t nhisto){

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptStat(0);

  TString fitType="Exp";
  TString filename=Form("h%d%sMassFit.root",nhisto,fitType.Data());

  TFile* fin=new TFile(filename.Data());
  if(!fin->IsOpen()){
    cout<<filename.Data()<<" not found, exit"<<endl;
    return;
  }

  TKey* key=((TKey*)((TList*)fin->GetListOfKeys())->At(fin->GetNkeys()-1));
  TCanvas* canvas=((TCanvas*)fin->Get(key->GetName()));
  if(!canvas){
    cout<<"Canvas not found"<<endl;
    return;
  }
  canvas->Draw();

  TH1F* hfit=(TH1F*)canvas->FindObject("fhistoInvMass");
  if(!hfit){
    canvas->ls();
    cout<<"Histogram not found"<<endl;
    return;
  }

  TF1* funcBkgRecalc=(TF1*)hfit->FindObject("funcbkgRecalc");
  if(!funcBkgRecalc){
    cout<<"Background fit function (final) not found"<<endl;
    return;
  }

  TF1* funcBkgFullRange=(TF1*)hfit->FindObject("funcbkgFullRange");
  if(!funcBkgFullRange){
    cout<<"Background fit function (side bands) not found"<<endl;
    return;
  }

  Int_t nbins=hfit->GetNbinsX();
  Double_t min=hfit->GetBinLowEdge(1), width=hfit->GetBinWidth(1);
  TH1F* hsubRecalc=(TH1F*)hfit->Clone("hsub");
  hsubRecalc->SetMarkerColor(kRed);
  hsubRecalc->SetLineColor(kRed);
  hsubRecalc->GetListOfFunctions()->Delete();
  TH1F* hsubFullRange=(TH1F*)hfit->Clone("hsub");
  hsubFullRange->SetMarkerColor(kGray+2);
  hsubFullRange->SetLineColor(kGray+2);
  hsubFullRange->GetListOfFunctions()->Delete();
  for(Int_t i=0;i<nbins;i++){
    //Double_t x=min+i*0.5*width;
    Double_t x1=min+i*width, x2=min+(i+1)*width;
    Double_t ycont=hfit->GetBinContent(i+1);
    Double_t y1=funcBkgRecalc->Integral(x1,x2) / width;//funcBkgRecalc->Eval(x);
    hsubRecalc->SetBinContent(i+1,ycont-y1);
    Double_t y2=funcBkgFullRange->Integral(x1,x2) / width;//funcBkgFullRange->Eval(x);
    hsubFullRange->SetBinContent(i+1,ycont-y2);
  }

  TCanvas* c=new TCanvas("c","subtraction");
  c->cd();
  hsubRecalc->DrawClone();
  hsubFullRange->DrawClone("sames");

  for(Int_t i=0;i<nbins;i++){
    if(hsubRecalc->GetBinContent(i+1)<0) hsubRecalc->SetBinContent(i+1,0);
    if(hsubFullRange->GetBinContent(i+1)<0) hsubFullRange->SetBinContent(i+1,0);
  }

  TCanvas *cvnewfits=new TCanvas("cvnewfits", "new Fits",1200,600);
  cvnewfits->Divide(2,1);

  AliHFMassFitter fitter1(hsubRecalc,min,min+nbins*width,1,1);
  fitter1.MassFitter(kFALSE);
  fitter1.DrawHere(cvnewfits->cd(1));

  AliHFMassFitter fitter2(hsubFullRange,min,min+nbins*width,1,1);
  fitter2.MassFitter(kFALSE);
  fitter2.DrawHere(cvnewfits->cd(2));

  canvas->SaveAs(Form("h%d%sMassFit.png",nhisto,fitType.Data()));
  c->SaveAs(Form("h%d%sSubtr.png",nhisto,fitType.Data()));
  cvnewfits->SaveAs(Form("h%d%sFitNew.png",nhisto,fitType.Data()));
}
