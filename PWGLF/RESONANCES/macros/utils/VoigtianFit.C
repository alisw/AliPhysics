/* Author: Anders G. Knospe, The University of Texas at Austin
     Created: 28 January 2014

     This macro is a simple implementation of a Voigtian peak fit atop a background.  A Voigtian peak is a convolution of a Breit-Wigner peak and a Gaussian.  Here, it is implemented using the TMath::Voigt() function in ROOT.

     Input Parameters:
     h: input histogram containing a peak
     formula: character string containing the functional form of the background, as it would be implemented in a TF1.
       For example, a quadratic background can be implemented as "pol2(0)" or "[0]+[1]*x+[2]*x*x".
     par: array containing peak parameters {mass,resolution,width}
     fix: array containing three flags {fix mass,fix resolution,fix width}, which allow any of these peak parameters to be fixed during the fit
     range: array containing four values {R1,R2,R3,R4}
       R1 (R4) is the lower (upper) limit of the fit range
       R2 (R3) is the lower (upper) limit of the peak range.  The peak range is excluded from the fit initially so that the background can be estimated.  The peak range is included in all subsequent fits.
     file: file to which results may be saved (optional)
*/

int VoigtianFit(TH1D* h,char* formula,double* par,int* fix,double* range,TFile* file)
{
  if(!h){cerr<<"Error in VoigtianFit(): missing histogram"<<endl; return 1;}
  
  char name[500]; sprintf(name,h->GetName());
  cout<<"using VoigtianFit():\n  h="<<name<<"\n  background formula="<<formula<<endl;
  cout<<"  mass="<<par[0];
  if(fix[0]) cout<<" (fixed)"<<endl;
  else cout<<" (free)"<<endl;
  cout<<"  resolution="<<par[1];
  if(fix[1]) cout<<" (fixed)"<<endl;
  else cout<<" (free)"<<endl;
  cout<<"  width="<<par[2];
  if(fix[2]) cout<<" (fixed)"<<endl;
  else cout<<" (free)"<<endl;
  cout<<"  range="<<range[0]<<" "<<range[1]<<" "<<range[2]<<" "<<range[3]<<endl;
  if(file) cout<<"  file="<<file->GetName()<<endl;
  else cout<<"  output not saved"<<endl;

  int j;

  //create copy of histogram h with peak removed
  TH1D* a=(TH1D*) h->Clone(Form("%s_nopeak",name));
  for(j=h->GetXaxis()->FindBin(1.000001*range[1]);j<=h->GetXaxis()->FindBin(0.999999*range[2]);j++){
    a->SetBinContent(j,0.);
    a->SetBinError(j,0.);
  }

  //get initial estimate of background
  TF1* fb=new TF1(Form("%s_back",name),formula,range[0],range[3]);
  a->Fit(fb,"RQN");

  //define peak fit function
  int vp=fb->GetNpar();
  TF1* fp=new TF1(Form("%s_peak",name),Form("%s+[%i]*TMath::Voigt(x-[%i],[%i],[%i])",formula,vp,vp+1,vp+2,vp+3),range[0],range[3]);

  //set initial parameter values, only the peak height is free
  for(j=0;j<vp;j++){fp->SetParameter(j,fb->GetParameter(j)); fp->FixParameter(j,fb->GetParameter(j));}
  fp->SetParameter(vp,h->GetBinContent(h->GetXaxis()->FindBin(0.5*(range[2]-range[1])))-fb->Eval(0.5*(range[2]-range[1])));
  for(j=0;j<3;j++){fp->SetParameter(vp+j+1,par[j]); fp->FixParameter(vp+j+1,par[j]);}

  h->Fit(fp,"RQN");

  if(!fix[2]){//release width
    fp->ReleaseParameter(vp+3);
    fp->SetParError(vp+3,0.1*fp->GetParameter(vp+3));
    h->Fit(fp,"RQN");
  }

  if(!fix[1]){//release resolution 
    fp->ReleaseParameter(vp+2);
    fp->SetParError(vp+2,0.1*fp->GetParameter(vp+2));
    h->Fit(fp,"RQN");
  }

  if(!fix[0]){//release mass
    fp->ReleaseParameter(vp+1);
    fp->SetParError(vp+1,0.1*fp->GetParameter(vp+1));
    h->Fit(fp,"RQN");
  }

  //release background constant parameter
  fp->ReleaseParameter(0);
  fp->SetParError(0,fb->GetParError(0));
  h->Fit(fp,"RQN");

  //release other background parameters
  for(j=1;j<vp;j++){
    fp->ReleaseParameter(j);
    fp->SetParError(j,fb->GetParError(j));
  }
  h->Fit(fp,"RQN");

  //final fit
  cerr<<"doing final fit"<<endl;
  h->Fit(fp,"RNI");

  //save output
  if(file){
    file->cd();
    h->Write();
    a->Write();
    fb->Write();
    fp->Write();
  }

  return 0;
}
