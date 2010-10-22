#include <fstream>
#include <Riostream.h>
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
#include <TObjectTable.h>
#include <TDatabasePDG.h>

#include <AliMultiDimVector.h>
#include "AliHFMassFitter.h"
#include <AliSignificanceCalculator.h>

#include <fstream>

//decCh:
//- 0 = kDplustoKpipi
//- 1 = kD0toKpi
//- 2 = kDstartoKpipi
//- 3 = kDstoKKpi
//- 4 = kD0toKpipipi
//- 5 = kLambdactopKpi

//Note: writefit=kTRUE writes the root files with the fit performed but it also draw all the canvas, so if your computer is not powerfull enough I suggest to run it in batch mode (root -b)

Bool_t charmCutsOptimization(Double_t nsigma=3,Int_t decCh=1,Int_t fitbtype=0,Int_t rebin=2,TString part="both"/*"A" anti-particle, "P" particle*/,Bool_t writefit=kTRUE,Double_t sigma=0.012,Int_t minentries=50,Double_t *rangefit=0x0,TString hname="hMass_"){

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);


  TString filename="AnalysisResults.root",dirname="PWG3_D2H_Significance",listname="coutputSig",mdvlistname="coutputmv";


  Int_t pdg;
  Double_t mass;

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

  cout<<"Mass = "<<mass<<endl;
  
  Int_t countFitFail=0,countSgnfFail=0,countNoHist=0,countBkgOnly=0;
  ofstream outcheck("output.dat");
  ofstream outdetail("discarddetails.dat");

  outcheck<<"ptbin\tmdvGlobAddr\thistIndex\tSignif\tS\tB"<<endl;
  outdetail<<"ptbin\tmdvGlobAddr\thistIndex\trelErrS\t\tmean_F-mass"<<endl;
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
    hstat->SaveAs("hstat.png");
  }else{
    cout<<"Warning! fHistNEvents not found in "<<listname.Data()<<endl;
  }

  Bool_t isMC=kFALSE;
  TH1F* htestIsMC=(TH1F*)histlist->FindObject("hSgn_0");
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

  //Check wheter histograms are filled
  for(Int_t i=0;i<nhist;i++){
    TString name=Form("%s%d",hname.Data(),i);
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

      if(indexes[ih+i*nhistforptbin] != -1){
	TString name=Form("%s%d",hname.Data(),indexes[ih+i*nhistforptbin]);
	TH1F* h=(TH1F*)histlist->FindObject(name.Data());

	Int_t nbin=((TH1F*)histlist->FindObject(name))->GetNbinsX();
	Double_t min=((TH1F*)histlist->FindObject(name))->GetBinLowEdge(7);
	Double_t max=((TH1F*)histlist->FindObject(name))->GetBinLowEdge(nbin-5)+((TH1F*)histlist->FindObject(name))->GetBinWidth(nbin-5);
	if(rangefit) {
	  min=rangefit[0];
	  max=rangefit[1];
	}

	AliHFMassFitter fitter(h,min, max,rebin,fitbtype);
	fitter.SetInitialGaussianMean(mass);
	fitter.SetInitialGaussianSigma(sigma);

	if(ih==0) fitter.InitNtuParam(Form("ntuPtbin%d",i));
	// fitter.SetHisto(h);
	// fitter.SetRangeFit(min,max);
	//fitter.SetRangeFit(1.68,2.05);

	//fitter.SetType(fitbtype,0);

	Bool_t ok=fitter.MassFitter(kFALSE);
        if(!ok){
	  cout<<"FIT NOT OK!"<<endl;
	  //   ok=fitter.RefitWithBkgOnly(kFALSE);
	  //   if (ok){ //onlybkg
	  countBkgOnly++;
	  //     Double_t bkg=0,errbkg=0.;
	  //     fitter.Background(nsigma,bkg,errbkg); 
	  outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t 0\t xxx"<<"\t bkgonly"<<endl;
	  //     mdvSgnf->SetElement(ih,0);
	  //     mdvSgnferr->SetElement(ih,0);
	  //     mdvS->SetElement(ih,0);
	  //     mdvSerr->SetElement(ih,0);
	  //     mdvB->SetElement(ih,bkg);
	  //     mdvBerr->SetElement(ih,errbkg);
	  //   }else{ //bkg fit failed
	  //     cout<<"Setting to 0"<<endl;
	  //     outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t fit failed also with only bkg"<<endl;
	  //     countFitFail++;
	  //     mdvSgnf->SetElement(ih,0);
	  //     mdvSgnferr->SetElement(ih,0);
	  //     mdvS->SetElement(ih,0);
	  //     mdvSerr->SetElement(ih,0);
	  //     mdvB->SetElement(ih,0);
	  //     mdvBerr->SetElement(ih,0);
	  //   }
	}else{ //fit ok!

	  if(writefit) fitter.WriteCanvas(Form("h%d",ih+i*nhistforptbin),"./",nsigma);
	  Double_t signif=0, signal=0, background=0, errSignif=0, errSignal=0, errBackground=0;
	  fitter.Signal(nsigma,signal,errSignal);
	  fitter.Background(nsigma,background,errBackground);
	  Double_t meanfit=fitter.GetMean();
	  Double_t sigmafit=fitter.GetSigma();
	  //outcheck<<"Sigma = "<<sigmafit<<endl;
	  hSigma[i]->Fill(sigmafit);

	  // if(sigmafit > 0.03){
	  //   //refit
	  //   fitter.Reset();
	  //   fitter.SetHisto(h);
	  //   fitter.SetRangeFit(1.8,1.93); //WARNING!! this is candidate dependant!! (change)
	  //   fitter.RebinMass(rebin);
	  //   fitter.SetInitialGaussianMean(mass);
	  //   fitter.SetInitialGaussianSigma(sigma);
	  //   ok=fitter.MassFitter(kFALSE);
	  //   if(ok){
	  //     meanfit=fitter.GetMean();
	  //     sigmafit=fitter.GetSigma();
	  //     fitter.Signal(nsigma,signal,errSignal);
	  //     fitter.Background(nsigma,background,errBackground);
	  //     if(writefit) fitter.WriteCanvas(Form("h%dbis",ih+i*nhistforptbin),"./",nsigma);
	  //     hSigma[i]->Fill(fitter.GetSigma());
	  //   }
	  // } //sigma check done and fit recalc

	  if(ok==kTRUE && sigmafit < 0.03 && signal > 0 && background > 0){
	    fitter.Significance(nsigma,signif,errSignif);
	    if(signif >0){
	      if(errSignal/signal < 0.35 && TMath::Abs(meanfit-mass)<0.01){
		outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t"<<signif<<" +- "<<errSignif<<"\t"<<signal<<" +- "<<errSignal<<"\t"<<background<<" +- "<<errBackground<<endl;
		mdvSgnf->SetElement(ih,signif);
		mdvSgnferr->SetElement(ih,errSignif);
		mdvS->SetElement(ih,signal);
		mdvSerr->SetElement(ih,errSignal);
		mdvB->SetElement(ih,background);
		mdvBerr->SetElement(ih,errBackground);
	      
	      }else{
		ok=fitter.RefitWithBkgOnly(kFALSE);
		if (ok){
		  countBkgOnly++;
		  Double_t bkg=0,errbkg=0.;
		  Double_t nsigmarange[2]={mass-nsigma*sigma,mass+nsigma*sigma};
		  fitter.Background(nsigmarange[0],nsigmarange[1],bkg,errbkg); 
		  outcheck<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t 0\t "<<bkg <<"\t bkgonly"<<endl;
		  outdetail<<i<<"\t\t "<<ih<<"\t\t"<<indexes[ih+i*nhistforptbin]<<"\t"<<errSignal/signal<<"\t\t "<<meanfit-mass<<endl;
		  mdvSgnf->SetElement(ih,0);
		  mdvSgnferr->SetElement(ih,0);
		  mdvS->SetElement(ih,0);
		  mdvSerr->SetElement(ih,0);
		  mdvB->SetElement(ih,bkg);
		  mdvBerr->SetElement(ih,errbkg);
		}
	      }//only bkg
	    }//check signif>0
	    else{ 
	      countSgnfFail++;
	      cout<<"Setting to 0 (fitter results meaningless)"<<endl;
	      outcheck<<"\t S || B || sgnf negative";
	      outcheck<<endl;
	      mdvSgnf->SetElement(ih,0);
	      mdvSgnferr->SetElement(ih,0);
	      mdvS->SetElement(ih,0);
	      mdvSerr->SetElement(ih,0);
	      mdvB->SetElement(ih,0);
	      mdvBerr->SetElement(ih,0);
	    } 
	  } //end fit ok!
	}
      }else{ //check on histo

	countNoHist++;
	cout<<"Setting to 0 (indexes = -1)"<<endl;
	outcheck<<"\t histo not accepted for fit";
	outcheck<<endl;
	mdvSgnf->SetElement(ih,0);
	mdvSgnferr->SetElement(ih,0);
	mdvS->SetElement(ih,0);
	mdvSerr->SetElement(ih,0);
	mdvB->SetElement(ih,0);
	mdvBerr->SetElement(ih,0);
	
      }

     
      cout<<mdvS->GetElement(ih)<<"\t"<<mdvB->GetElement(ih)<<endl;

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
 

  fout->Close();

  outcheck<<"\nSummary:\n - Total number of histograms: "<<nhist<<"\n - "<<count<<" histograms with more than "<<minentries<<" entries; \n - Too few entries in histo "<<countNoHist<<" times;\n - Fit failed "<<countFitFail<<" times \n - no sense Signal/Background/Significance "<<countSgnfFail<<" times\n - only background "<<countBkgOnly<<" times"<<endl;
  outcheck.close();
  return kTRUE;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// which=0 plot significance
//      =1 plot signal
//      =2 plot background
//      =3 plot S/B
// maximize = kTRUE (default) if you want to fix the step of the variables not shown to the value that maximize the significance. Note that these values are saved in fixedvars.dat
// readfromfile = kTRUE (default is kFALSE) if you want to read the value fixed in a previous run of this function (e.g. significance or signal maximization)


void showMultiDimVector(Int_t n=2,Int_t which=0, Bool_t plotErrors=kFALSE, Bool_t maximize=kTRUE,Bool_t readfromfile=kFALSE){

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetFrameBorderMode(0);

  if((maximize && readfromfile) || (!maximize && !readfromfile)){
    cout<<"Error! maximize & readfromfile cannot be both kTRUE or kFALSE"<<endl;
    return;
  }

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

  //ask variables for projection
  for(Int_t k=0;k<nvarsopt;k++){
    cout<<k<<" "<<mdv->GetAxisTitle(k)<<endl;
  }
  cout<<"Choose "<<n<<" variable(s)"<<endl;
  for(Int_t j=0;j<n;j++){
    cout<<"var"<<j<<": ";
    cin>>variable[j];	
  }
  if(n==1) variable[1]=999;

  TCanvas* cvpj=new TCanvas(Form("proj%d",variable[0]),Form("%s wrt %s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data()));

  TMultiGraph* mg=new TMultiGraph(Form("proj%d",variable[0]),Form("%s wrt %s;%s;%s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),(mdv->GetAxisTitle(variable[0])).Data(),title.Data()));
  TLegend *leg=new TLegend(0.7,0.2,0.9,0.6,"Pt Bin");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

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
    for(Int_t ind=0;ind<ncuts;ind++)maxInd[ind]=0;

    Float_t sigMax0=cal->GetMaxSignificance(maxInd,0);
    for(Int_t ic=0;ic<ncuts;ic++){
      cutvalues[ic]=((AliMultiDimVector*)fin->Get(nameS.Data()))->GetCutValue(ic,maxInd[ic]);

      //setting step of fixed variables
      if(readfromfile){ //from file
	fixedvars[ic]=allfixedvars[i+ic];
      }

      if(maximize) { //using the values which maximize the significance
	fixedvars[ic]=maxInd[ic];
	//write to output fixedvars.dat
	writefixedvars<<fixedvars[ic]<<"\t";
      }
    }
    //output file: return after each pt bin
    if(maximize) writefixedvars<<endl;

    printf("Maximum of significance for Ptbin %d found in bin:\n",i);
    for(Int_t ic=0;ic<ncuts;ic++){
      printf("   %d\n",maxInd[ic]);
      printf("corresponding to cut:\n");
      printf("   %f\n",cutvalues[ic]);
    }

    printf("Significance = %f +- %f\n",sigMax0,mvess->GetElement(maxInd,0));
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
    TString ptbinrange=Form("%.1f < p_{t} < %.1f GeV/c",mdv->GetPtLimit(0),mdv->GetPtLimit(1));

    if(n==2) {
      gStyle->SetPalette(1);
      TH2F* hproj=mdv->Project(variable[0],variable[1],fixedvars,0);
      hproj->SetTitle(Form("%s wrt %s vs %s (Ptbin%d %.1f<pt<%.1f);%s;%s",title.Data(),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data(),i,mdv->GetPtLimit(0),mdv->GetPtLimit(1),(mdv->GetAxisTitle(variable[0])).Data(),mdv->GetAxisTitle(variable[1]).Data()));

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

	fixedvars[variable[0]]=k; //variable[0] is the index of the variable on which we project. This variable must increase, the others have been set before

	x[k]=mdv->GetCutValue(variable[0],k);
	errx[k]=mdv->GetCutStep(variable[0])/2.;
	
	y[k]=mdv->GetElement(fixedvars,0);
	erry[k]=mdverr->GetElement(fixedvars,0);
	if(which==3){

	}
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
  nhistinrange=0;

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
    }
  }
  cout<<nhistinrange<<" index will be returned"<<endl;
  Int_t index[nvar4opt-nvargiven];
  Int_t k=0;
  for (Int_t i=0;i<nvar4opt;i++){
     if(valsgiven[i]==kFALSE) {
       //cout<<"kTRUE==>"<<i<<endl;
       index[k]=i;
       k++;
     }
  }

  Int_t *rangeofhistos=new Int_t[nhistinrange];
  for(Int_t i=0;i<nvar4opt-nvargiven;i++){ //loop on number of free variables
    cout<<"Info: incrementing "<<vct->GetAxisTitle(index[i])<<endl;
    for(Int_t j=0;j<vct->GetNCutSteps(i);j++){ //loop on steps of each free variable
      allvals[index[i]]=vct->GetCutValue(index[i],j);
      rangeofhistos[i*vct->GetNCutSteps(i)+j]=GetNHistFromValues(vct,ptbin,allvals);
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
