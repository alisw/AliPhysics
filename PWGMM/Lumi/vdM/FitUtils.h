
const char *g_fit_model_name[] = {
  "GP2", "GP6","G","NUM"};

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
  if (fit_type == 0) return 4; // GP2
  if (fit_type == 1) return 6; // GP6
  if (fit_type == 2) return 3; // G
  if (fit_type == 3) return 1; // NUM, this will be dummy parameters    
  return -1;
}

//-------------------------------------------------------
//Initialize the fit model
//-------------------------------------------------------

void Fit_model_init(Int_t fit_type, Double_t r_max, TF1 *model)
{
  Double_t width_max = 0.5; // the maximum width of the gaussian term
  if (fit_type == 0) { // GP2
    model->SetParNames("R","#mu","#sigma","p2");
    model->SetParameters(r_max,0.0,0.03,4);
    model->SetParLimits(0,0.5*r_max,2.0*r_max);
    model->SetParLimits(2,0,width_max);
  } else if (fit_type == 1) { // GP6
    model->SetParNames("R","#mu","#sigma","p2","p4","p6");
    model->SetParameters(r_max,0.0,0.12,4,10,20);
    model->SetParLimits(0,0,2.0*r_max);
    model->SetParLimits(2,0,width_max);
  } else if (fit_type == 2) { // G
    model->SetParNames("R","#mu","#sigma");
    model->SetParameters(r_max,0.0,0.03);
    model->SetParLimits(0,0,2.0*r_max);
    model->SetParLimits(2,0,width_max);
  }
}

Double_t Do_Numeric_Integration(Int_t n, Double_t *sep, Double_t *rate, Double_t *rate_err, 
				Double_t *area, Double_t *rate_zero, Double_t *par, Double_t *par_err)
{
  // compute bin widths
  Double_t *widths = new Double_t [n];
  widths[0] = TMath::Abs(sep[0]-sep[1]);
  for(Int_t i=1;i<n-1;i++)
    widths[i] = 0.5*TMath::Abs(sep[i-1]-sep[i]) + 0.5*TMath::Abs(sep[i+1]-sep[i]);
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
// Main entry point for the fit of rate as a function of separation
//-------------------------------------------------------

Double_t Fit_rate_separation(Int_t n, Double_t *sep, Double_t *rate, Double_t *rate_err, Int_t fit_type,
			     Double_t *area, Double_t *rate_zero, Double_t *par, Double_t *par_err,
			     Int_t scan, Int_t scan_type, Int_t bc)
// fit n points of rate(sep) using model given by fit_type
// output area and error, rate at zero and error, parameters and error
// returns chi2/dof ... returns -1 if fit does not converge
// (scan, scan_type and bc are only used to give a proper name to canvas if necessary)
{
  // if fit_type = 3, do numeric integration
  return Do_Numeric_Integration(n,sep,rate,rate_err,area,rate_zero,par,par_err);
  
  // define the limits in the separation axis
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  if (sep[n-1]>sep[0]) { 
    sep_min = sep[0] - 0.5*TMath::Abs(sep[0]-sep[1]);
    sep_max = sep[n-1] + 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
  } else {
    sep_max = sep[0] + 0.5*TMath::Abs(sep[0]-sep[1]);
    sep_min = sep[n-1] - 0.5*TMath::Abs(sep[n-1]-sep[n-2]);
  }
  // set up the model fit
  TF1 *fit_model = NULL;
  if (fit_type == 0) {
    fit_model = new TF1("fit_model",fit_GP2, sep_min, sep_max,4);
  } else if (fit_type == 1) {
    fit_model = new TF1("fit_model",fit_GP6, sep_min, sep_max,6);
  } else if (fit_type == 2) {
    fit_model = new TF1("fit_model",fit_G, sep_min, sep_max,3);
  } else {
    cout << " Fit model " << fit_type << " not known " << endl;
    exit(-105);
  }
  // initialize the model
  // --> find maximum rate
  Double_t rate_max = 0;
  for(Int_t i=0;i<n;i++) {if (rate[i]>rate_max) rate_max = rate[i];}
  
  // --> initialize
  Fit_model_init(fit_type, rate_max, fit_model);
  
  // define a TGraph to perform the fit
  TGraphErrors *gr = new TGraphErrors(n,sep,rate,NULL,rate_err);
 
  // fit and check
  TFitResultPtr r = gr->Fit("fit_model","Q0RS");
  Double_t chi2=-1;
  Int_t npar = Get_number_par(fit_type);
  if (gMinuit->fCstatu.Contains("CONVERGED"))  chi2 = fit_model->GetChisquare()/((Double_t)fit_model->GetNDF());
  else {
    for(Int_t j=0;j<npar;j++) {
      par[j] = 0.;
      par_err[j] = 0;
    }
    area[0]=area[1]=0;
    rate_zero[0]=rate_zero[1]=0;
    return chi2;
  }
  for(Int_t j=0;j<npar;j++) {
    par[j] = fit_model->GetParameter(j);
    par_err[j] = fit_model->GetParError(j);
  }
  TMatrixDSym cov = r->GetCovarianceMatrix();
  Double_t epsilon = 0.001;
  area[0] = fit_model->Integral(sep_min, sep_max,epsilon);
  area[1] = fit_model->IntegralError(sep_min, sep_max,fit_model->GetParameters(),cov.GetMatrixArray(),epsilon);
  rate_zero[0] = fit_model->Eval(fit_model->GetParameter(1));
  rate_zero[1] = rate_zero[0]*(fit_model->GetParError(0))/(fit_model->GetParameter(0));    

  // check if integral gives 'reasonable' value
  Double_t lim = 10.0;
  if ((area[0] > (lim*rate_max)) || (area[1] > (lim*rate_max)) ||
      (rate_zero[0] > (lim*rate_max)) || (rate_zero[1] > (lim*rate_max))) {
    for(Int_t j=0;j<npar;j++) {
      par[j] = 0.;
      par_err[j] = 0;
    }
    area[0]=area[1]=0;
    rate_zero[0]=rate_zero[1]=0;
    return -2;
  }

  // plot if a particular bc is chosen
  if (bc > -1) {
    char name[120];
    if (scan_type == 1) sprintf(name,"Scan_%d_x_bc_%d",scan,bc);
    if (scan_type == 2) sprintf(name,"Scan_%d_y_bc_%d",scan,bc);    
    TCanvas *c = new TCanvas(name,name,800,600);
    c->cd();
    gr->SetMarkerStyle(20);
    gr->Draw("ap");
    fit_model->Draw("same");
    cout << name << " A= " << area[0] << "+/-" << area[1] << " R0= " << rate_zero[0] <<"+/-"<<rate_zero[1] << endl;
  }

  
  // work done
  return chi2;
}
