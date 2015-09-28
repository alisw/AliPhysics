void RebinSpectrum(TH1** hi,TH1** hf,TF1* f,int n=1,int dofit=0,TH1* ht=0){
  /* Author: Anders G. Knospe, The University of Texas at Austin
     Created: 9 May 2014
     Last Update: 22 September 2015

     This macro rebins histograms, including cases with incompatible bins.

     Input Parameters:
     hi: an array of 1-dimensional histograms containing the original distribution
     hf: an array of 1-dimensional histograms with the new desired binning
     f: a TF1 fit function that is used to estimate how to divide the contents of bins that must be split
     n: the size of the arrays hi and hf
     dofit: if dofit=1 the histogram hi[0] will be fit with function f; if dofit==0 it is assumed that function f already describes the distribution (i.e., the fit was done outside this macro)
     ht: a 1-dimensional histogram that will be used for fitting if dofit==1; it should include the total uncertainties (statistical plus systematic, but you might want to exclude uncertainties correlated between pT bins)

     hi and hf are arrays so that systematic uncertainties can be accounted for.  hi[0] should contain the central values and statistical uncertainties.  The other histograms in hi[*] (if any) should contain the systematic uncertainties (more than one type of systematic uncertainty is allowed).  The output histogram hf[0] will contain the rebinned central values and statistical uncertainties.  The other output histograms hf[*] (if any) will contain the systematic uncertainties.  Going from hi[0] to hf[0], the statistical uncertainties are added in quadrature.  For the other histograms, the systematic uncertainties are added linearly.  It is assumed that all of the histograms in hi and ht have the same binning, and that all of the histograms in hf have the same binning.

     Explanation: The macro rebins from the input histogram to the output histogram.  Let A be a bin in the input histogram and let B be a bin in the output histogram.  If A lies entirely within B, the situation is simple: the contents of bin A will be added to bin B.  If A is split by one (or both) of the edges of B, a fit function is used to determine the fraction of the content of bin A that will be added to bin B.  The fitting may be done externally (set dofit==0), or this macro can be used to do the fit over a limited range (set dofit==1).  If dofit==1, any constraints applied to the fit function should be applied before it is passed to this macro.  Generally, the fit will be over A and the two bins adjacent to it.  Becuase the fit is done over only three bins, you cannot use fit functions with more than three free paramaters.  You should be careful about using functions with three free parameters, since you may end up fitting statistical noise.  For this reason, I use a simple exponential function in the example below.  You will need to judge for yourself whether a three-parameter function is really necessary to describe your distribution.  For the high-pT part of most spectra, an exponential should be sufficient.  Around the maximum of the pT distribution, it may be necessary to use a function with three free parameters.  If dofit==0 (i.e, you do the fit yourself), these warnings do no apply and you can use whatever complicated function you want.

     EXAMPLE: Here is some example code that may be useful:

     TFile* f=TFile::Open("measured_file.root");//open the file that contains your measured histogram

     TH1F* hi[3];//input histogram array
     int n=3;//size of arrays hi and hf
     hi[0]=(TH1F*) f->Get("measured_histogram_stat");//get your measured histogram with statistical uncertainties
     hi[1]=(TH1F*) f->Get("measured_histogram_sys1");//get your measured histogram with first type of systematic uncertainties
     hi[2]=(TH1F*) f->Get("measured_histogram_sys2");//get your measured histogram with second type of systematic uncertainties
     //expand the array as necessary to include as many types of systematic uncertainties as you need

     TH1F* ht=(TH1F*) hi[0]->Clone("ht");
     for(int j=1;j<=ht->GetNbinsX();j++) ht->SetBinError(j,sqrt(pow(hi[0]->GetBinError(j),2)+pow(hi[1]->GetBinError(j),2)));
     //ht is the input histogram with total (uncorrelated) uncertainties. It should have the central values of hi[0] and the uncertainties should be the sum of the statistical and systematic uncertainties in hi. (You may want to exclude sources of systematic uncertainty that are correlated between pT bins.)  This is the histogram that will be fit if dofit==1 (if dofit==0, you can set ht=0).

     //define your new binning
     int n=9;//the new number of bins here
     float bins[10]={0,1,2,3,};//new bin boundaries
     TH1F* hf[2];//output (rebinned) histogram array
     hf[0]=new TH1F("rebinned_stat","",n,bins);//will be filled with rebinned central values and statistical uncertainties (added in quadrature)
     hf[1]=new TH1F("rebinned_sys1","",n,bins);//will be filled with rebinned central values and first type of systematic uncertainties (added linearly)
     hf[2]=new TH1F("rebinned_sys2","",n,bins);//will be filled with rebinned central values and second type of systematic uncertainties (added linearly)
     //The number of histograms in hf should be the same as the number in hi.

     TF1* g=new TF1("fit","[0]*exp([1]*x)",0.,10.);//define your fit function
     g->SetParameters(1.,-1.);//set its parameters
     //You can do the fit in your own code, or let RebinSpectrum do it over a limited range as needed (depending on the value of dofit).  This macro can be used even if no fit is needed, but a placeholder fit function will still need to be defined.
     int dofit=1;//let RebinSpectrum do the fit

     gROOT->LoadMacro("*path/PWGLF/RESONANCES/macros/utils/RebinSpectrum.C");

     RebinSpectrum(hi,hf,f,n,dofit,ht);
     //hf now contains the rebinned histograms

     f->Close();
  */

  int j,k,l;
  double ai,bi,di,e,af,bf,df,v[100],u[100],d;

  for(j=0;j<n;j++) if(!hi[j]){cerr<<"Error in RebinSpectrum(): missing input histogram "<<j<<endl; return;}
  for(j=0;j<n;j++) if(!hf[j]){cerr<<"Error in RebinSpectrum(): missing output histogram "<<j<<endl; return;}
  if(!f){cerr<<"Error in RebinSpectrum(): missing fit function"<<endl; return;}
  if(n<1 || n>99){cerr<<"Error in RebinSpectrum(): invalid value for n "<<n<<endl; return;}
  if(dofit && !ht){cerr<<"Error in RebinSpectrum(): missing fit histogram"<<endl; return;}

  if(dofit) cerr<<"Info in RebinSpectrum(): will do fit inside macro RebinSpectrum.C"<<endl;
  else cerr<<"Info in RebinSpectrum(): using external fit"<<endl;

  for(j=0;j<n;j++) for(l=0;l<=hf[j]->GetNbinsX()+1;l++){
      //clear the output histograms
      hf[j]->SetBinContent(l,0.);
      hf[j]->SetBinError(l,0.);
    }

  for(l=1;l<=hf[0]->GetNbinsX();l++){
    af=hf[0]->GetXaxis()->GetBinLowEdge(l);
    bf=hf[0]->GetXaxis()->GetBinLowEdge(l+1);
    df=hf[0]->GetXaxis()->GetBinWidth(l);
    e=1.e-5*df;

    for(j=0;j<n;j++) v[j]=u[j]=0.;

    for(k=1;k<=hi[0]->GetNbinsX();k++){
      ai=hi[0]->GetXaxis()->GetBinLowEdge(k);
      bi=hi[0]->GetXaxis()->GetBinLowEdge(k+1);
      di=hi[0]->GetXaxis()->GetBinWidth(k);

      if(bi<=af+e || bf<=ai+e) continue;//bin k of hi completely outside bin l of hf
      else if(af<=ai+e && bi<=bf+e){
	//bin k of hi completely contained within bin l of hf
	d=1.;
      }else if(ai<=af+e || bf<=bi+e){
	//bin k of hi is split by the edge(s) of bin l of hf
	cerr<<"Info in RebinSpectrum(): splitting hi("<<ai<<","<<bi<<") hf("<<af<<","<<bf<<")"<<endl;
	if(dofit){
	  if(k==1 || ht->GetBinContent(k-1)<1.e-30){
	    //k is the first non-empty bin of ht, fit bin k and the two following bins
	    f->SetRange(ht->GetBinLowEdge(k),ht->GetBinLowEdge(k+3));
	  }else if(k==ht->GetNbinsX() || ht->GetBinContent(k+1)<1.e-30){
	    //k is the last non-empty bin of ht, fit bin k and the two preceeding bins
	    f->SetRange(ht->GetBinLowEdge(k-2),ht->GetBinLowEdge(k+1));
	  }else{
	    //k is neither the first nor the last non-empty bin of ht
	    f->SetRange(ht->GetBinLowEdge(k-1),ht->GetBinLowEdge(k+2));
	  }

	  ht->Fit(f,"NR");
	}

	if(ai<=af+e && bi<=bf+e){
	  //bin k of hi is split by the low edge of bin l of hf
	  d=f->Integral(af,bi)/f->Integral(ai,bi);
	}else if(af<=ai+e && bf<=bi+e){
	  //bin k of hi is split by the high edge of bin l of hf
	  d=f->Integral(ai,bf)/f->Integral(ai,bi);
	}else if(ai<=af+e && bf<=bi+e){
	  //bin k of hi completely contains bin l of hf
	  d=f->Integral(af,bf)/f->Integral(ai,bi);
	}else{
	  cerr<<"Error in RebinSpectrum(): undefined case: hi("<<ai<<","<<bi<<") hf("<<af<<","<<bf<<")"<<endl;
	  continue;
	}
      }

      for(j=0;j<n;j++){
	v[j]+=d*hi[j]->GetBinContent(k)*di;//add the content of bin k of hi[j] to the total
	if(!j) u[j]+=pow(d*hi[j]->GetBinError(k)*di,2);//add (in quadrature) the uncertainty of bin k of hi[0] to the total
	else u[j]+=d*hi[j]->GetBinError(k)*di;//add (linearly) the uncertainty of bin k of hi[j] to the total (for j>=1)
      }
    }

    for(j=0;j<n;j++){
      v[j]/=df;
      if(!j) u[j]=sqrt(u[j]);
      u[j]/=df;

      //fill the output histograms
      hf[j]->SetBinContent(l,v[j]);
      hf[j]->SetBinError(l,u[j]);
    }
  }

  return;
}
