#ifndef __FIT_UTILS_HH_
#define __FIT_UTILS_HH_

#include "TDatime.h"
#include "TRandom2.h"

//-------------------------------------------------------
// global variables
//-------------------------------------------------------

const char *g_fit_model_name[] = {"GP2", "GP6","G","NUM", "DG"};

// these are needed by minuit ... ugly, but that is the way it is ...
Int_t n_minuit;
Double_t x_minuit[100]; // 100 is a 'large' number, normally there are only 25 entries per scan
Double_t y_minuit[100]; 
Double_t ye_minuit[100];
Int_t nll_model = 1; // 0 = gauss; 1 = poisson
Double_t scale_minuit = 30.0;
Int_t fit_type_minuit = -1;

//-------------------------------------------------------
// Model of a gaussian 
//-------------------------------------------------------

Double_t fit_G(Double_t *x, Double_t *p)
{
	Double_t dx = x[0]-p[1];
	Double_t dx2 = dx*dx;
	Double_t Gx = TMath::Exp(-dx2/(2.0*p[2]*p[2]));
	return p[0]*Gx;
}

//-------------------------------------------------------
// Model of a double gaussian 
//-------------------------------------------------------

Double_t fit_DG(Double_t *x, Double_t *p)
{
	const Double_t sqrt2pi_inv = 1.0/TMath::Sqrt(TMath::TwoPi());
	Double_t dx  = x[0] - p[1]; //p[1] = center
	Double_t dx2 = dx*dx;
	Double_t N1 = sqrt2pi_inv/p[2];
	Double_t G1 = N1*TMath::Exp(-dx2/(2.0*p[2]*p[2]));
	Double_t N2 = sqrt2pi_inv/p[3];
	Double_t G2 = N2*TMath::Exp(-dx2/(2.0*p[3]*p[3])); //p[2], p[3] = width of each Gaussian
	return p[0] * (G1 + p[4]*G2); //kimc: p[0]/p[0]*p[4]: height, p[1]: common center, p[2]/p[3]: width
}

//-------------------------------------------------------
// Model of a gaussian and pol2 term
//-------------------------------------------------------

Double_t fit_GP2(Double_t *x, Double_t *p)
{
	Double_t dx = x[0]-p[1];
	Double_t dx2 = dx*dx;
	Double_t Gx = TMath::Exp(-dx2/(2.0*p[2]*p[2]));
	return p[0]*Gx*(1+p[3]*dx2);
}

//-------------------------------------------------------
// Model of a gaussian and up to  pol6 terms
//-------------------------------------------------------

Double_t fit_GP6(Double_t *x, Double_t *p)
{
	Double_t dx = x[0]-p[1];
	Double_t dx2 = dx*dx;
	Double_t dx4 = dx2*dx2;
	Double_t dx6 = dx4*dx2;
	Double_t Gx = TMath::Exp(-dx2/(2.0*p[2]*p[2]));
	return p[0]*Gx*(1+p[3]*dx2+p[4]*dx4+p[5]*dx6);
}

//-------------------------------------------------------
// return number of parameters of a given model
//-------------------------------------------------------

Int_t Get_number_par(Int_t fit_type)
{
	if      (fit_type == 0) return 4; // GP2
	else if (fit_type == 1) return 6; // GP6
	else if (fit_type == 2) return 3; // G
	else if (fit_type == 3) return 1; // NUM, this will be dummy parameters
	else if (fit_type == 4) return 5; // DG
	else return -1;
}

//-------------------------------------------------------
//Initialize the fit model
//-------------------------------------------------------

void Fit_model_init(Int_t fit_type, Double_t r_max, TF1 *model)
{
	const Double_t width_max = 0.5; // the maximum width of the gaussian term

	if (fit_type == 0) // GP2
	{
		model->SetParNames("R","#mu","#sigma","p2");
		model->SetParameters(r_max,0.0,0.03,4);
		model->SetParLimits(0,0.5*r_max,2.0*r_max);
		model->SetParLimits(2,0,width_max);
	}
	else if (fit_type == 1) // GP6
	{
		model->SetParNames("R","#mu","#sigma","p2","p4","p6");
		model->SetParameters(r_max,0.0,0.12,4,10,200);
		model->SetParLimits(0,0,2.0*r_max);
		model->SetParLimits(2,0.005,width_max);
	}
	else if (fit_type == 2) // G
	{
		model->SetParNames("R","#mu","#sigma");
		model->SetParameters(r_max,0.0,0.03);
		model->SetParLimits(0,0,2.0*r_max);
		model->SetParLimits(2,0.005,width_max);
	} 
	else if (fit_type == 4) // DG
	{
		model->SetParNames("R0", "#mu", "#sigma_{0}", "#sigma_{1}", "R1");
		model->SetParameters(r_max*0.9, 0.0, 0.1, 0.15, 0.1);
		model->SetParLimits(0, 0,    r_max); //h0
		model->SetParLimits(2, 0.01, width_max); //w0
		model->SetParLimits(3, 0.01, width_max); //w1
		model->SetParLimits(4, 0,    1); //h1
	}

	return;
}

Double_t Do_Numeric_Integration(
		Int_t n, Double_t *sep, Double_t *rate, Double_t *rate_err,
		Double_t *area, Double_t *rate_zero, Double_t *par, Double_t *par_err
		)
{
	// compute bin widths
	Double_t *widths = new Double_t [n];
	widths[0] = TMath::Abs(sep[0]-sep[1]);
	for(Int_t i=1;i<n-1;i++) widths[i] = 0.5*TMath::Abs(sep[i-1]-sep[i]) + 0.5*TMath::Abs(sep[i+1]-sep[i]);
	widths[n-1] = TMath::Abs(sep[n-1]-sep[n-2]);

	// find maximum rate
	rate_zero[0] = rate_zero[1] = 0;
	for(Int_t i=0;i<n;i++) {
		if (rate[i]>rate_zero[0]) {
			rate_zero[0] = rate[i];
			rate_zero[1] = rate_err[i];
		}
	}

	// compute area
	area[0]=area[1]=0.0;
	for(Int_t i=0;i<n;i++) {
		area[0] += (widths[i]*rate[i]);
		area[1] += (widths[i]*widths[i]*rate_err[i]*rate_err[i]);
	}
	area[1] = TMath::Sqrt(area[1]);

	// set dummy values
	Double_t chi2 = 1.0;
	par[0] = rate_zero[0];
	par_err[0] = rate_zero[1];

	// clean up
	delete [] widths;

	// return
	return chi2;
}

//-------------------------------------------------------
// approximation to the logarithm of the Gamma function
// https://en.wikipedia.org/wiki/Stirling%27s_approximation
//-------------------------------------------------------

Double_t ApproxLogGamma(Double_t x)
{
  Double_t a = 0.5*(TMath::Log(2.0*TMath::Pi())-TMath::Log(x));
  Double_t b = 12.0*x - (1.0/(10*x));
  Double_t c = TMath::Log(x+(1.0/b));
  Double_t d = x*(c-1);
  Double_t e = a+d;
  return e;
}


//-------------------------------------------------------
// negative log-likelihood function as required by minuit
//-------------------------------------------------------

void negative_log_likelihood(Int_t &npar, Double_t *gin, Double_t &nll, Double_t *par, Int_t iflag)
{
  Double_t sum = 0;
  for (Int_t i = 0; i < n_minuit; i++) {
    // get the model prediction for the current parameters
    Double_t f_i = 0;
    if (fit_type_minuit == 0) f_i = fit_GP2(&x_minuit[i], par);
    else if (fit_type_minuit == 1) f_i = fit_GP6(&x_minuit[i], par);
    else if (fit_type_minuit == 2) f_i = fit_G(&x_minuit[i], par);
    // compute the likelihood for this term
    if (nll_model == 0) { // gaussian model
      if(y_minuit[i]<1e-6 || ye_minuit[i] < 1e-6) continue; // avoid 'zero' rates
      Double_t a_i = 1./(TMath::Sqrt(2.0*TMath::Pi())*ye_minuit[i]);
      Double_t b_i =  (y_minuit[i]-f_i)/ye_minuit[i];
      Double_t g_i = a_i*TMath::Exp(-0.5*b_i*b_i);
      sum += TMath::Log(g_i);
    } else { // poisson model
      /*
      Double_t a_i = TMath::Power(f_i,y_minuit[i]);
      Double_t b_i = TMath::Gamma(y_minuit[i]+1);
      Double_t c_i = TMath::Exp(-f_i);
      sum += TMath::Log(a_i*c_i/b_i);
      */
      Double_t r = scale_minuit*y_minuit[i];
      Double_t a_i =r*TMath::Log(f_i);
      Double_t b_i = -ApproxLogGamma(r+1);
      Double_t c_i = -f_i;
      sum += (a_i + b_i + c_i);      
    }  
  } // end of loop over data points
  nll = -sum; // return negative log likelihood
}

 /*
void negative_log_likelihood(Int_t &npar, Double_t *gin, Double_t &nll, Double_t *par, Int_t iflag)
{
  Double_t sum = 0;
  for (Int_t i = 0; i < n_minuit; i++) {
    if(y_minuit[i]<1e-6 || ye_minuit[i] < 1e-6) continue; // avoid 'zero' rates
    Double_t nsig = y_minuit[i]/ye_minuit[i];
    Double_t f_i = fit_GP2(&x_minuit[i], par);
    if (nsig>4.0) { // gaussian model
      Double_t a_i = 1./(TMath::Sqrt(2.0*TMath::Pi())*ye_minuit[i]);
      Double_t b_i =  (y_minuit[i]-f_i)/ye_minuit[i];
      Double_t g_i = a_i*TMath::Exp(-0.5*b_i*b_i);
      sum += TMath::Log(g_i);
    } else { // poisson model
      Double_t a_i = TMath::Power(f_i,y_minuit[i]);
      Double_t b_i = TMath::Gamma(y_minuit[i]+1);
      Double_t c_i = TMath::Exp(-f_i);
      sum += TMath::Log(a_i*c_i/b_i);
    }  
  } // end of loop over data points
  nll = -sum; // return negative log likelihood
}
 */
//-------------------------------------------------------
// Main entry point for the fit of rate as a function of separation
// when using minuit
//-------------------------------------------------------

Double_t Fit_rate_separation_minuit(Int_t n, Double_t *sep, Double_t *rate, Double_t *rate_err, Int_t fit_type,
			     Double_t *area, Double_t *rate_zero, Double_t *par, Double_t *par_err,
			     Int_t scan, Int_t scan_type, Int_t bc)
// fit n points of rate(sep) using model given by fit_type
// output area and error, rate at zero and error, parameters and error
// returns chi2/dof ... returns -1 if fit does not converge
// (scan, scan_type and bc are only used to give a proper name to canvas if necessary)
{
  // if fit_type = 3, do numeric integration
  if(fit_type == 3)
  return Do_Numeric_Integration(n,sep,rate,rate_err,area,rate_zero,par,par_err);

  // compute some numbers needed later on
  // --> define the limits in the separation axis
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  if (sep[n-1]>sep[0]) { 
    sep_min = sep[0] - 0.5*TMath::Abs(sep[0]-sep[1]);
    sep_max = sep[n-1] + 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
  } else {
    sep_max = sep[0] + 0.5*TMath::Abs(sep[0]-sep[1]);
    sep_min = sep[n-1] - 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
  }

  // --> find maximum rate
  Double_t rate_max = 0;
  for(Int_t i=0;i<n;i++) {if (rate[i]>rate_max) rate_max = rate[i];}


  // set up the fit model
  TF1 *fit_model = NULL;
  Int_t nPar = 0; 
  if (fit_type == 0) {
    nPar = Get_number_par(fit_type);
    fit_model = new TF1("fit_model",fit_GP2, sep_min, sep_max, nPar);
  } else if (fit_type == 1) {
    nPar = Get_number_par(fit_type);
    fit_model = new TF1("fit_model",fit_GP6, sep_min, sep_max, nPar);
  } else if (fit_type == 2) {
    nPar = Get_number_par(fit_type);
    fit_model = new TF1("fit_model",fit_G, sep_min, sep_max, nPar);
  } else if (fit_type == 4) {
    nPar = Get_number_par(fit_type);
    fit_model = new TF1("fit_model",fit_DG, sep_min, sep_max, nPar);
  } else {
    cout << " Fit model " << fit_type << " not known " << endl;
    exit(-105);
  }

  
  // set up minuit
  // --> set up the model
  fit_type_minuit = fit_type;
  // --> fill in the data to be fit
  n_minuit = n;
  for(Int_t i=0;i<n_minuit;i++) {
    x_minuit[i] = sep[i];
    y_minuit[i] = rate[i];
    ye_minuit[i] = rate_err[i];
  }
  // define minuit
  TMinuit m(nPar);
  m.SetFCN(negative_log_likelihood);
  m.SetPrintLevel(-1); // -1 quiet, 0 normal, 1 verbose
  m.SetErrorDef(0.5); // 1 for chi2 fit, 0.5 for negative log-likelihood 
  // parameter no., name, start value, step size, range min., range max.
  // range min = range max = 0 -> no limits
  m.DefineParameter(0, "R", scale_minuit*rate_max, 0.01,
		    scale_minuit*0.5*rate_max, scale_minuit*2.0*rate_max);
  m.DefineParameter(1, "#mu", 0.0, 0.01, 0, 0);
  m.DefineParameter(2, "#sigma", 0.03, 0.01, 0.0,sep_max);
  if (fit_type == 0) m.DefineParameter(3, "p_{2}", 4.0, 0.01, 0,0);
  else if (fit_type == 1) {
    m.DefineParameter(3, "p_{2}", 4.0, 0.01, 0,0);
    m.DefineParameter(4, "p_{4}", 10.0, 0.01, 0,0);
    m.DefineParameter(5, "p_{6}", 20.0, 0.01, 0,0);      
  } 
  
  // now ready for minimization step
  m.Migrad();

  // if fit did not converge set output to zero
  Double_t chi2=-1;
  if(!m.fCstatu.Contains("CONVERGED")) {
     for(Int_t j=0;j<nPar;j++) {
      par[j] = 0.;
      par_err[j] = 0;
    }
    area[0]=area[1]=0;
    rate_zero[0]=rate_zero[1]=0;
    //return chi2;
  }

  // fit converged, get output
  for(Int_t i=0;i<nPar;i++) m.GetParameter(i, par[i], par_err[i]);
  TMatrixDSym cov(nPar); 
  m.mnemat( cov.GetMatrixArray(), nPar);
  
  // pass parameters to the model
  fit_model->SetParErrors(par_err);
  fit_model->SetParameters(par);

  // compute output
  // -- proxy of chi2
  chi2 = 0;
  for(Int_t i=0;i<n;i++) {
    Double_t rate_model = (fit_model->Eval(x_minuit[i]))/scale_minuit;
    Double_t pull = ((y_minuit[i]<1e-6 || ye_minuit[i] < 1e-6)? 0
		     :(y_minuit[i]-rate_model)/ye_minuit[i]);
    chi2 += (pull*pull);
  }
  Double_t ndf = n_minuit-nPar;
  chi2 = chi2/ndf;
  
  // -- area and rate
  Double_t epsilon = 0.1; //  Double_t epsilon = 0.001;
  area[0] = fit_model->Integral(sep_min, sep_max,epsilon);
  area[1] = fit_model->IntegralError(sep_min, sep_max,fit_model->GetParameters(),cov.GetMatrixArray(),epsilon);
  rate_zero[0] = fit_model->Eval(fit_model->GetParameter(1));
  rate_zero[1] = rate_zero[0]*(fit_model->GetParError(0))/(fit_model->GetParameter(0));    
  for(Int_t i=0;i<2;i++) {area[i]/=scale_minuit;rate_zero[i]/=scale_minuit;}
  
  // check if integral gives 'reasonable' value
  Double_t lim = 20.0;
  if ((area[0] > (lim*rate_max)) || (area[1] > (lim*rate_max)) ||
      (rate_zero[0] > (lim*rate_max)) || (rate_zero[1] > (lim*rate_max))) {
    for(Int_t j=0;j<nPar;j++) {
      par[j] = 0.;
      par_err[j] = 0;
    }
    area[0]=area[1]=0;
    rate_zero[0]=rate_zero[1]=0;
    return -2;
  }

  // plot if a particular bc is chosen
  if (bc > -1) {
    gStyle->SetOptFit(1);
    TGraphErrors *gr = new TGraphErrors(n,sep,rate,NULL,rate_err);
    char name[120];
    if (scan_type == 1) sprintf(name,"Scan_%d_x_bc_%d",scan,bc);
    if (scan_type == 2) sprintf(name,"Scan_%d_y_bc_%d",scan,bc);    
    TCanvas *c = new TCanvas(name,name,800,600);
    c->cd();
    gr->SetMarkerStyle(20);
    gr->Draw("ap");
    fit_model->SetParameter(0,par[0]/scale_minuit);
    fit_model->Draw("same");
    cout << name << " A= " << area[0] << "+/-" << area[1] << " R0= " << rate_zero[0] <<"+/-"<<rate_zero[1] << " chi2 " << chi2 << endl;
    if (!m.fCstatu.Contains("CONVERGED")) cout << "Fit did not converge for bc " << bc << " in " << name << endl;
  }
  
  // work done
  return chi2;
}

//-------------------------------------------------------
// Main entry point for the fit of rate as a function of separation
// using TGraphErrors for a chi2 fit
//-------------------------------------------------------

Double_t Fit_rate_separation(
		Int_t n, Double_t *sep, Double_t *rate, Double_t *rate_err, Int_t fit_type,
		Double_t *area, Double_t *rate_zero, Double_t *par, Double_t *par_err,
		const char* cName
		)
// fit n points of rate(sep) using model given by fit_type
// output area and error, rate at zero and error, parameters and error
// returns chi2/dof ... returns -1 if fit does not converge
{
	// if fit_type = 3, do numeric integration
	if (fit_type == 3) return Do_Numeric_Integration(n,sep,rate,rate_err,area,rate_zero,par,par_err);

	// define the limits in the separation axis
	Double_t sep_min = 0;
	Double_t sep_max = 0;
	if (sep[n-1]>sep[0])
	{ 
		sep_min = sep[0] - 0.5*TMath::Abs(sep[0]-sep[1]);
		sep_max = sep[n-1] + 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
	}
	else
	{
		sep_max = sep[0] + 0.5*TMath::Abs(sep[0]-sep[1]);
		sep_min = sep[n-1] - 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
	}

	// set up the model fit
	TF1 *fit_model = NULL;
	if      (fit_type == 0) { fit_model = new TF1("fit_model", fit_GP2, sep_min,sep_max,Get_number_par(fit_type)); }
	else if (fit_type == 1) { fit_model = new TF1("fit_model", fit_GP6, sep_min,sep_max,Get_number_par(fit_type)); }
	else if (fit_type == 2) { fit_model = new TF1("fit_model", fit_G,   sep_min,sep_max,Get_number_par(fit_type)); }
	else if (fit_type == 4) { fit_model = new TF1("fit_model", fit_DG,  sep_min,sep_max,Get_number_par(fit_type)); }
	else { 	cout << " Fit model " << fit_type << " not known " << endl; exit(-105); }

	// initialize the model
	// --> find maximum rate
	Double_t rate_max = 0;
	for (Int_t i=0;i<n;i++) { if (rate[i]>rate_max) rate_max = rate[i]; }

	// --> initialize
	Fit_model_init(fit_type, rate_max, fit_model);

	// define a TGraph to perform the fit
	TGraphErrors *gr = new TGraphErrors(n, sep, rate, NULL, rate_err);

	// fit and check
	TFitResultPtr r = gr->Fit("fit_model", "Q0RS");
	Double_t fitQuality = fit_model->GetChisquare()/((Double_t)fit_model->GetNDF()); //kimc
	string dName = cName; //Temporal copy, to prevent contamination (kimc)

	//-------------------------------------------

	const Int_t reFit_max = 20;
	const Double_t fitQuality_TL = 0.1; //Tolerance, lower (kimc)
	const Double_t fitQuality_TU = 10.; //Tolerance, upper (kimc)

    gStyle->SetOptFit(1);

	//Retry fit if it failed
    if ( !gMinuit->fCstatu.Contains("CONVERGED") ||
		 (fitQuality < fitQuality_TL) ||
		 (fitQuality > fitQuality_TU) )
	{
		//Temporary TH1 for pit parameters
		TH1F* H1Temp = new TH1F(Form("%s_temp", dName.c_str()), "", n, sep[0], sep[n-1]); H1Temp->Sumw2();
		for (int a=0; a<n; a++)
		{
			const int   xBin = H1Temp->GetXaxis()->FindBin(sep[a]);
			const float yVal = rate[a];
			const float yErr = rate_err[a];
			H1Temp->SetBinContent(xBin, yVal);
			H1Temp->SetBinError  (xBin, yErr);
		}
		const double tMax  = H1Temp->GetMaximum();
		const double tMean = H1Temp->GetMean();
		const double tRMS  = H1Temp->GetRMS();

		//Retry fit
		int reFit = 0;
		while (reFit < reFit_max)
		{
			cout <<Form("Retry fit for %s: %i...", dName.c_str(), reFit) <<endl;

			if (reFit > 0) //For GP2, GP6, and G, from 2nd iteration
			{
				TDatime DT;
				TRandom2 RD2(DT.GetTime() + reFit);
				double tPar[3] = {0};

				tPar[0] = RD2.Uniform(tMax*0.9, tMax*1.1);
				tPar[1] = RD2.Uniform(-2*fabs(tMean), 2*fabs(tMean));
				tPar[2] = RD2.Uniform(tRMS*0.5, tRMS*2);
				//cout <<Form(" fit seeds: max %4.3f, mean %4.3f, and sigma %4.3f\n", tPar[0],tPar[1],tPar[2]);

				//Release Gaussian parameters and set again
				for (int a=0; a<Get_number_par(fit_type); a++)
				{
					fit_model->ReleaseParameter(a);
					if (a<3) fit_model->SetParameter(a, tPar[a]);
					else     fit_model->SetParameter(a, 1.);
				}
			}

			r = gr->Fit("fit_model", "Q0RS");
			fitQuality = fit_model->GetChisquare()/((Double_t)fit_model->GetNDF()); //kimc
			if ( gMinuit->fCstatu.Contains("CONVERGED") &&
				 fitQuality > fitQuality_TL &&
				 fitQuality < fitQuality_TU ) break; //Converged: stop
			else reFit++;
		}//Refit

		if (reFit == reFit_max)
		{
			cout <<Form(" Fit failed: %s, chi2/NDF: %4.3f", dName.c_str(), fitQuality) <<endl;
			dName = Form("%s_FAIL", dName.c_str());
			fitQuality = -1;
		}
		else dName = Form("%s_reFit%i", dName.c_str(), reFit);
		H1Temp->Delete();
	}//Retry fit

	/*
    TCanvas *c = new TCanvas(dName.c_str(), dName.c_str(), 800, 600); c->cd();
    TH1 *h = (TH1*)gr->GetHistogram(); h->SetTitle(Form("%s;Separation (mm); Rate (Hz)", dName.c_str())); h->Draw();
    gr->SetMarkerStyle(20); gr->Draw("pe same");
    fit_model->Draw("same");
    c->Print(Form("../Fill-%d/QA_fits/%s.%s", g_vdm_Fill, dName.c_str(), FFormat));
    //c->Print(Form("../Fill-%d/QA_fits/%s.png", g_vdm_Fill, dName.c_str()));
	*/

	#if 1
	//Plot for public note, June 25
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	TCanvas* c1 = new TCanvas(dName.c_str(), dName.c_str(), 800, 600);
	c1->cd()->SetLogy();

	/*
	TString T1 = dName.c_str();
	T1.ReplaceAll("Fill", "Fill ");
	T1.ReplaceAll("_ODCBBD", ", ");
	T1.ReplaceAll("_OpticalIntensityCorrFBCT", "");
	T1.ReplaceAll("_VBAandVBC", "V0, ");
	T1.ReplaceAll("_TVX", " T0, ");
	T1.ReplaceAll("_F1", "");
	T1.ReplaceAll("_scan", "scan ");
	T1.ReplaceAll("_i", ", bunch index ");
	T1.ReplaceAll("_bc", " (ID ");
	T1.ReplaceAll("_x", "), horizontal");
	T1.ReplaceAll("_y", "), vertical");
	*/

	TH1F* H1 = (TH1F*)gr->GetHistogram();
	H1->GetXaxis()->SetTitleOffset(1.1);
	H1->GetXaxis()->SetRangeUser(-0.65, 0.65);
	H1->GetYaxis()->SetRangeUser(0.1, H1->GetMaximum()*2.0);
	//H1->SetTitle(Form("%s;Separation [mm];Rate [Hz]", T1.Data()));
	H1->SetTitle(";Separation [mm];Rate [Hz]");
	H1->DrawCopy();
    gr->SetMarkerStyle(20);
	gr->Draw("pe same");
	fit_model->SetLineStyle(2);
    fit_model->Draw("same");

    TLegend *L1 = new TLegend(0.4, 0.2, 0.6, 0.35);
    L1->SetMargin(0);
    L1->SetBorderSize(0);
    L1->SetTextAlign(13);
    L1->AddEntry((TObject*)0, "ALICE", "");
    L1->AddEntry((TObject*)0, "pp #sqrt{s} = 13 TeV", "");
	L1->Draw();

	c1->Print(Form("%s.eps", c1->GetName()));
	#endif

	//-------------------------------------------

	const Int_t    npar = Get_number_par(fit_type);
	const Double_t chi2 = fitQuality;
 
	if (chi2 == -1) //Failed fit - nullify results
	{
		for (Int_t j=0; j<npar; j++)
		{
			par[j] = 0.;
			par_err[j] = 0;
		}
		area[0] = area[1] = 0.;
		rate_zero[0] = rate_zero[1] = 0.;
		return chi2;
	}
	else
	{
		for (Int_t j=0; j<npar; j++)
		{
			par[j] = fit_model->GetParameter(j);
			par_err[j] = fit_model->GetParError(j);
		}
	}

	//Get R00 and its error
	Double_t epsilon = 0.1;
	area[0] = fit_model->Integral(sep_min, sep_max, epsilon);
	area[1] = fit_model->IntegralError(
			sep_min, sep_max, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray(), epsilon
			);
	rate_zero[0] = fit_model->Eval(fit_model->GetParameter(1));
	rate_zero[1] = rate_zero[0] * (fit_model->GetParError(0))/(fit_model->GetParameter(0));

	// check if integral gives 'reasonable' value
	const Double_t lim = 10.0;
	if ( (area[0] > (lim*rate_max)) ||
		 (area[1] > (lim*rate_max)) ||
		 (rate_zero[0] > (lim*rate_max)) ||
		 (rate_zero[1] > (lim*rate_max)) )
	{
		for (Int_t j=0; j<npar; j++)
		{
			par[j] = 0.;
			par_err[j] = 0;
		}
		area[0] = area[1] = 0;
		rate_zero[0] = rate_zero[1] = 0;
		return -2;
	}

	return chi2;
}

#endif
