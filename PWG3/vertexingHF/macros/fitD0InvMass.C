#include <Riostream.h>
#include <TStyle.h>
#include <TFile.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <AliHFMassFitter.h>

#include <fstream>
#include <cmath>

//root[0] .x fitD0InvMass.C+

//input: nsigma os the number of sigma in which the significance (, signal ancg background) is calculated

//fit of histograms in InvMassD0.root. Must specify the list name where the histo are stored. Produce a root file with the fitted histos and a text file with signal, background and significance

void fitD0InvMass(Int_t rebin=0,TString listname="coutputmassD02",Int_t nsigma=3, TString pathin="./",TString pathout="./",Int_t btype=2){

  gStyle->SetOptFit(0111);
  gStyle->SetOptStat("nemrou");
  TString file = "D0InvMass.root";
  const Int_t npt=5;
  file.Prepend(pathin);
  cout<<"Opening "<<file<<endl;
  TFile *fin=new TFile(file.Data());
  if(!fin->IsOpen()){
    cout<<"File "<<file.Data()<<" not found, did you mean another? Please enter the name"<<endl;
    cin>>file;
    fin=new TFile(file.Data());
    if(!fin->IsOpen()){
      cout<<"File "<<file.Data()<<" not found, return"<<endl;
      return;
    }
  }

  TList *lista = (TList*) fin->Get(listname.Data());
  if(!lista) {
    cout<<listname<<" doesn't exist"<<endl;
    return;
  }

  TString fileout="table.dat";
  fileout.Prepend(pathout);
  ofstream out(fileout.Data());

  AliHFMassFitter *fitter=new AliHFMassFitter();

  cout<<"mean = "<<fitter->GetMean()<<"\tsigma= "<<fitter->GetSigma()<<endl;
  Bool_t init=kTRUE;

  out<<"Bakground type = "<<btype<<";\trebin = "<<rebin<<endl;
  Int_t i=0;
  for(Int_t ipt=0;ipt<npt;ipt++){ //npt

    out<<"************************************************************\n";
    out<<"* ptbin * entriesMass * entriesSgn * entriesBkg * entriesRfl *\n";

    //ALL w cuts
    TString nameall="histMass_";
    nameall+=(ipt+1);
    TH1F *hMass=(TH1F*)lista->FindObject(nameall);
    if(!hMass){
      cout<<nameall<<" not found"<<endl;
      return;
    }

    out<<"* "<< ipt+1 <<"  *  "<<hMass->GetEntries()<<"  *  ";
    //SIGNAL from MC
    TString namesgn="histSgn_";
    namesgn+=(ipt+1);
    TH1F *hSgn=(TH1F*)lista->FindObject(namesgn);
    if(!hSgn){
      cout<<namesgn<<" not found"<<endl;
      return;
    }
    
    out<<hSgn->GetEntries()<<"  *  ";
    //BACKGROUND from MC
    TString namebkg="histBkg_";
    namebkg+=(ipt+1);
    TH1F *hBkg=(TH1F*)lista->FindObject(namebkg);
    if(!hBkg){
      cout<<namebkg<<" not found"<<endl;
      return;
    }

    out<<hBkg->GetEntries()<<"  *  ";
    //REFLECTED SIGNAL from MC
    TString namerfl="histRfl_";
    namerfl+=(ipt+1);
    TH1F *hRfl=(TH1F*)lista->FindObject(namerfl);
    if(!hRfl){
      cout<<namerfl<<" not found"<<endl;
      return;
    }

    out<<hRfl->GetEntries()<<"  *\n";
    out<<"*************************************************\n";
    if(rebin!=0){
      hSgn->Rebin(rebin);
      hBkg->Rebin(rebin);
      hRfl->Rebin(rebin);
    }
    Double_t min, max;
    Int_t    nbin;
    Double_t sgn,errsgn,errsgn2=0,bkg,errbkg,errbkg2=0,sgnf,errsgnf,sgnfMC;
    Double_t meanfrfit,sigmafrfit, meanMC=hSgn->GetMean(),sigmaMC=hSgn->GetRMS();
    Double_t width=0;
    nbin=hMass->GetNbinsX();
    min=hMass->GetBinLowEdge(1);
    max=min+nbin*hMass->GetBinWidth(nbin);

    //FIT:

    TString namentu="ntupt3bin";
    if(init) fitter->InitNtuParam((char*)namentu.Data());
    // - all
    fitter->SetHisto(hMass);
    fitter->SetRangeFit(min, max);
    //fitter->SetRangeFit(1.83,1.89);
    fitter->SetType(btype,0);//(b,s)
    if(rebin!=0) fitter->RebinMass(rebin);
    Bool_t fitOK=fitter->MassFitter(kFALSE); //kFALSE = do not draw
    if(!fitOK) {
      out<<"Fit return kFALSE, skip "<<hMass->GetName()<<endl;
      fitter->Reset(); //delete histogram set
      continue;
    }
    width=fitter->GetHistoClone()->GetBinWidth(3);
    cout<<"\nChi^2 = "<<fitter->GetChiSquare()<<"\t Reduced Chi^2 = "<<fitter->GetReducedChiSquare()<<endl;
    meanfrfit=fitter->GetMean();
    sigmafrfit=fitter->GetSigma();
    out<<"mean = "<<meanfrfit<<" (MC "<<meanMC<<")"<<"\tsigma= "<<sigmafrfit<<" (MC "<<sigmaMC<<")"<<endl;
    //cout<<"old nbin = "<<hMass->GetNbinsX()<<"\tnew nbin = "<<fitter->GetHistoClone()->GetNbinsX()<<endl;
    if(meanfrfit<0 || sigmafrfit<0 || meanfrfit<min || meanfrfit>max) {
      cout<<"Fit failed, check"<<endl;
      out<<hMass->GetName();
      out<<": \nSgn not available"<<endl;
      out<<"Bkg  not available"<<endl;
      out<<"Sgnf not available"<<endl;
      out<<"*************************************************\n";
    }
    else{
      //cout<<"Writing..."<<endl;
      fitter->WriteHisto(pathout);
      hMass= fitter->GetHistoClone();
      //TH1F *hc=fitter->GetHistoClone();
      //TF1 *fbtest=hc->GetFunction("funcbkgRecalc"); //new version of ALiHFMassFitter
     
      fitter->FillNtuParam();
     
      Double_t limsx,limdx;
      limsx=meanfrfit-nsigma*sigmafrfit;
      limdx=meanfrfit+nsigma*sigmafrfit;
      // limsx=meanMC-nsigma*sigmaMC;
      //limdx=meanMC+nsigma*sigmaMC;


      //determine limit of nsigma in bins
      Int_t binsx,bindx;
      binsx=hMass->FindBin(limsx);
      if (limsx > hMass->GetBinCenter(binsx)) binsx++;
      bindx=hMass->FindBin(limdx);
      if (limdx < hMass->GetBinCenter(bindx)) bindx--;

      //reconvert bin in x
      Double_t sxr,dxr;
      sxr=hMass->GetBinLowEdge(binsx);
      dxr=hMass->GetBinLowEdge(bindx+1);

      fitter->Signal(sxr,dxr,sgn,errsgn);
      fitter->Background(sxr,dxr,bkg,errbkg);
      fitter->Significance(sxr,dxr,sgnf,errsgnf);
      

      Float_t inttot,intsgn,intsgnerr;

      Int_t np=-99;
      switch (btype){
      case 0: //expo
	np=2;
	break;
      case 1: //linear
	np=2;
	break;
      case 2: //pol2
	np=3;
	break;
      case 3: //no bkg
	np=1;
	break;
      }

      
      TF1 *fmass=hMass->GetFunction("funcmass");
      if (fmass){
	
	inttot=fmass->GetParameter(0);
	intsgn=fmass->GetParameter(np);
	intsgnerr=fmass->GetParError(np);
	//cout<<"i = "<<i<<"inttot = "<<inttot<<"\tintsgn = "<<intsgn<<"\tintsgnErr = "<<intsgnerr<<"\twidth = "<<width<<endl;
	Double_t errbkg2rel=errbkg/bkg*TMath::Sqrt((inttot-intsgn)/width/bkg);//intsgnerr/(inttot-intsgn)*TMath::Sqrt((inttot-intsgn)/width/bkg);
	Double_t errsgn2rel=errsgn/sgn*TMath::Sqrt(intsgn/width/sgn);
	errbkg2=errbkg2rel*bkg;
	errsgn2=errsgn2rel*sgn;
      }

           
      cout<<"bin sx = "<<binsx<<"\t xsx = "<<sxr<<endl;
      cout<<"bin dx = "<<bindx<<"\t xdx = "<<dxr<<endl;

      out<<hMass->GetName();
      Double_t sgnMC,bkgMC,rflMC;
      sgnMC=hSgn->Integral(binsx,bindx);
      bkgMC=hBkg->Integral(binsx,bindx);
      rflMC=hRfl->Integral(binsx,bindx);
      sgnfMC=sgnMC/TMath::Sqrt(sgnMC+bkgMC+rflMC);
      out<<": \nSgn  "<<sgn<<"+/-"<<errsgn<<" ("<<errsgn2<<")\tCompare with: "<<sgnMC<<endl;
      out<<"Bkg  "<<bkg<<"+/-"<<errbkg<<" ("<<errbkg2<<")\tCompare with: "<<bkgMC<<" + "<<rflMC<<" = "<<bkgMC+rflMC<<endl;
      out<<"Sgnf "<<sgnf<<"+/-"<<errsgnf<<"\tCompare with: "<<sgnfMC<<endl;
      out<<"sigmaS/S = "<<errsgn/sgn<<"\tCompare with 1/signif: "<<1./sgnfMC<<endl;
      out<<"nsigma considered for comparison = \ndx: "<<(dxr-meanfrfit)/sigmafrfit<<"\nsx: "<<(meanfrfit-sxr)/sigmafrfit<<endl;
      out<<"Mean = "<<meanfrfit<<"\tSigma = "<<sigmafrfit<<endl;
      out<<"*************************************************\n";
      i++;
    }
    
    sgn=0; bkg=0; sgnf=0;
    errsgn=0; errbkg=0; errsgnf=0;
    if (ipt == npt-1) fitter->WriteNtuple(pathout);

    out<<endl;
    fitter->Reset(); //delete histogram set

    init=kFALSE;

  }

  out.close();
}
