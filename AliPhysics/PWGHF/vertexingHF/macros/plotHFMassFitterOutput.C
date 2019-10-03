#include <Riostream.h>
#include <TH1F.h>
#include <TF1.h>
#include <TAttFill.h>
#include <TFile.h>
#include <TCanvas.h>
#include <AliHFMassFitter.h>

//root[0] .L plotHFMassFitterOutput.C+
//root[1] plotHFMassFitterOutput(2)

//input: pt bin you want to plot; path of the HFMassFitterOutput.root file

//three lines of fit:
//grey= background with parameters from the fit on the side-bands
//red= background with parameters from the total fit
//blue= total fit

void plotHFMassFitterOutput(Int_t ptbin, TString path="./"){
  TString stdfilename="HFMassFitterOutput.root";
  stdfilename.Prepend(path);

  cout<<"Opening "<<stdfilename<<endl;
  TFile *fin=new TFile(stdfilename.Data());
  if(!fin->IsOpen()){
    cout<<"File "<<stdfilename<<" not found"<<endl;
    return;
  }
  TString histoname="histMass_";
  histoname+=ptbin;

  TH1F *h=(TH1F*)fin->Get(histoname);
  if(!h) {
    cout<<histoname<<" is not here, sorry!"<<endl;
    histoname="fhistoInvMass";
    h=(TH1F*)fin->Get(histoname);
    if(!h){
      cout<<histoname<<" is not here, sorry!"<<endl;
      cout<<"Write the name of the histo required (e.g. histMassSum<factor>_<ptbin>): ";
      cin>> histoname;
      h=(TH1F*)fin->Get(histoname);
      if(!h){
	cout<<histoname<<" is not in "<<stdfilename<<" check it, please!"<<endl;
	  return;
      }
    } else {
      histoname+=ptbin;
      h->SetName(histoname);
    }
  }


    if(h){
    cout<<histoname<<" will be drawn!"<<endl;
    TCanvas *c1=new TCanvas("c1",histoname);
     //funcbkgfullrange_faint
    cout<<h<<endl;
    (h->GetFunction("funcbkgFullRange"))->SetLineColor(14);
    (h->GetFunction("funcbkgFullRange"))->SetLineStyle(4);

    Double_t xmin,xmax;
    (h->GetFunction("funcbkgFullRange"))->GetRange(xmin,xmax);
    Int_t nfreepar=h->GetFunction("funcbkgFullRange")->GetNumberFreeParameters();
    cout<<"nfreepar = "<<nfreepar<<endl;

    cout<<"Range = ("<<xmin<<", "<<xmax<<")"<<endl;
    cout<<"Bin Width "<<h->GetBinWidth(3)<<endl;
    cout<<"Initial parameters:\n";
    cout<<"par0= "<<h->GetFunction("funcbkgFullRange")->GetParameter(0)<<endl;
    cout<<"par1= "<<h->GetFunction("funcbkgFullRange")->GetParameter(1)<<endl;
    cout<<"par2= "<<h->GetFunction("funcbkgFullRange")->GetParameter(2)<<endl;
    cout<<"Formula= "<<h->GetFunction("funcbkgFullRange")->GetExpFormula()<<endl;

    TF1 *fmass=h->GetFunction("funcmass");
    cout<<"Parametri massa input:\n";
    cout<<"par0= "<<h->GetFunction("funcmass")->GetParameter(0)-h->GetFunction("funcmass")->GetParameter(nfreepar)<<endl;
    cout<<"par1= "<<h->GetFunction("funcmass")->GetParameter(1)<<endl;
    cout<<"par2= "<<h->GetFunction("funcmass")->GetParameter(2)<<endl;

 
    AliHFMassFitter *fitter=new AliHFMassFitter(h,xmin,xmax,1,2,0);
    fitter->SetSideBands(kFALSE);

    //funcbkgfullrange_deep
    
    TF1 *fbkgfullrange_deep=new TF1("fbkgfullrange_deep",fitter,&AliHFMassFitter::FitFunction4Bkg,xmin,xmax,nfreepar,"AliHFMassFitter","FitFunction4Bkg"); 
     
    
    fbkgfullrange_deep->SetParameter(0,(fmass->GetParameter(0)-fmass->GetParameter(nfreepar))); //lin exp pol2 no-bkg
    if(nfreepar>=2) fbkgfullrange_deep->SetParameter(1,fmass->GetParameter(1)); //lin exp pol2
    if(nfreepar==3) fbkgfullrange_deep->SetParameter(2,fmass->GetParameter(2)); //pol2

    fbkgfullrange_deep->SetMinimum(0);

    cout<<"Final parameters:\n";
    cout<<"par0= "<<fbkgfullrange_deep->GetParameter(0)<<"\tcompare with "<<fbkgfullrange_deep->Integral(xmin,xmax)<<endl;
    cout<<"par1= "<<fbkgfullrange_deep->GetParameter(1)<<endl;
    cout<<"par2= "<<fbkgfullrange_deep->GetParameter(2)<<endl;

    c1->cd();
    fbkgfullrange_deep->SetLineStyle(1);
    fbkgfullrange_deep->SetLineColor(2);
    h->Draw();
    fbkgfullrange_deep->Draw("sames");
    fmass->Draw("sames");
  
  } else {
    cout<<histoname<<" is not here, sorry!"<<endl;
    return;
  }
}
