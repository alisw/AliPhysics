/** @file 
    @brief Implementation of an Resolution class */
#include "flow/AliFMDFlowResolution.h"
#include "flow/AliFMDFlowUtil.h"
#include "flow/AliFMDFlowBessel.h"
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
//#include <cmath>

//====================================================================
void 
AliFMDFlowResolution::Add(Double_t psiA, Double_t psiB) 
{ 
  Double_t diff    = NormalizeAngle(fOrder * (psiA - psiB));
  Double_t contrib = cos(diff);
  AliFMDFlowStat::Add(contrib);
}

//____________________________________________________________________
Double_t 
AliFMDFlowResolution::Correction(UShort_t k) const 
{ 
  Double_t e;
  return Correction(k, e); 
}

//____________________________________________________________________
Double_t 
AliFMDFlowResolution::Correction(UShort_t, Double_t& e2) const 
{ 
  e2 = fSqVar / fN;
  return sqrt(2) * sqrt(fabs(fAverage));
}

//____________________________________________________________________
void
AliFMDFlowResolution::Draw(Option_t* option) 
{
  TGraph* g = new TGraph(100);
  for (UShort_t i = 0; i < g->GetN(); i++) { 
    Double_t x = -1. + 2. / 100 * i;
    Double_t y = sqrt(2) * sqrt(fabs(x));
    g->SetPoint(i, x, y);
  }
  g->SetName("naive_res");
  g->SetTitle("Naive Resolution Function");
  g->GetHistogram()->SetXTitle("<cos(n(|#Psi_{A}-#Psi_{B}|))>");
  g->GetHistogram()->SetYTitle("<cos(n(#Psi-#Psi_{R}))>");
  g->GetHistogram()->SetDirectory(0);
  g->Draw(Form("lh %s", option));
}


//====================================================================
Double_t 
AliFMDFlowResolutionStar::Correction(UShort_t k, Double_t& e2) const 
{ 
  if (k > 4) return 0;
  Double_t delta = 0;
  Double_t chi   = Chi(fAverage, k, delta);
  Double_t dr    = 0;
  Double_t res   = Res(sqrt(2) * chi, k, dr);
  e2             = pow(dr * delta,2);
  return res;
}
//____________________________________________________________________
Double_t 
AliFMDFlowResolutionStar::Correction(UShort_t k) const 
{ 
  Double_t e;
  return Correction(k, e); 
}
//____________________________________________________________________
Double_t 
AliFMDFlowResolutionStar::Chi(Double_t res, UShort_t k, 
			      Double_t& delta) const 
{
  delta          = 1;
  Double_t chi   = 2;
  Double_t dr    = 0;
  for (UInt_t i = 0; i < 15; i++) { 
    if (Res(chi, k, dr) < res) chi += delta;
    else                       chi -= delta;
    delta /= 2;
  }
  return chi;
}
//____________________________________________________________________
Double_t 
AliFMDFlowResolutionStar::Res(Double_t chi, UShort_t k, Double_t& dr) const 
{ 
  // The resolution function is 
  // 
  //          sqrt(pi) x exp(-x/4) (f1(x^2/4) + f2(x^2/4))
  //   r(x) = --------------------------------------------
  //                          2 sqrt(2) 
  // 
  //        
  //        = c x (f1(y) - f2(y))
  //
  // where f1 is the modified Bessel function first kind I_{(k-1)/2}, 
  // and f2 is the modified Bessel function of the first kind
  // I_{(k+1)/2}, and 
  // 
  //          sqrt(pi) exp(-x^2/4) 
  //      c = --------------------,   y = x^2/4
  //              2 sqrt(2)
  // 
  // The derivative of the resolution function is 
  //
  //            c 
  //    r'(y) = - (4 sqrt(y) (f1'(y) - f2'(y)) - (4 y - 2)(f1(y) - f2(y)))
  //            2
  // 
  //            c                                    r(y)   
  //          = - (4 sqrt(y) (f1'(y) - f2'(y))) + --------- - sqrt(y) r(y)
  //		2             			  2 sqrt(y)       
  // 
  // Since dI_n(x)/dx = I_(n-1)(x) - n / x I_n(x), and substituting 
  // f3(y) = I_((k-3)/2)(y) 
  // 
  //            c  
  //    r'(y) = - ((4 - 2 k) f1(y) - (4 y + 2 k) f2(y) + 4 y f3(y))
  //            2   
  // 
  Double_t y  = chi * chi / 4;
  Double_t c  = sqrt(M_PI) * exp(-y) / (2 * sqrt(2));   
  
  // i[0] = I[(k-1)/2-1], i[1] = I[(k-1)/2], i[2] = I[(k-1)/2+1]
  Double_t i[3], di[3];
  AliFMDFlowBessel::Inu(Double_t(k-3)/2, Double_t(k+1)/2, y, i, di);
  
  Double_t r  = c * chi * (i[2] + i[1]);
  dr        = c / 2 * ((4 - 2*k)*i[1] - (4*y + 2*k)*i[2] + 4*y*i[0]);

  return r;  
}

//____________________________________________________________________
void
AliFMDFlowResolutionStar::Draw(Option_t* option) 
{
  TString opt(option);
  opt.ToLower();
  Bool_t chi = opt.Contains("chi");
  
  TH2* h = new TH2D(Form("star_%s_frame", (chi ? "chi" : "res")), 
		    Form("STAR %s Function for k=(1,2,4)", 
			 (chi ? "Chi" : "Resolution")),
		    100, 0, 1, 100, 0, (chi ? 3 : 1.5));
  h->SetXTitle("<cos(n(|#Psi_{A}-#Psi_{B}|))>");
  h->SetYTitle((chi ? "#chi" : "<cos(n(#Psi-#Psi_{R}))>"));
  h->SetStats(0);
  h->SetDirectory(0);
  h->Draw();

  for (UShort_t k = 1; k <= 4; k++) { 
    if (k == 3) continue;

    TGraphErrors* g = new TGraphErrors(100);
    for (UShort_t i = 0; i < g->GetN(); i++) { 
      Double_t e2 = 0;
      Double_t x  = 1. / 100 * i;
      Double_t y  = 0;
      Double_t c  = Chi(x, k, e2);
      if (chi) y  = c;
      else { 
	Double_t dr = 0;
	y           = Res(sqrt(2) * c, k, dr);
	e2          = pow(dr * e2,2);
      }
      g->SetPoint(i, x, y);
      g->SetPointError(i, 0, sqrt(e2));
    }
    g->SetLineColor(k);
    g->SetName(Form("star_%s_k%d", (chi ? "chi" : "res"), k));
    g->SetTitle(Form("STAR %s Function for k=%d", 
		     (chi ? "Chi" : "Resolution"), k));
    g->GetHistogram()->SetXTitle("<cos(n(|#Psi_{A}-#Psi_{B}|))>");
    g->GetHistogram()->SetYTitle((chi ? "#chi" : "<cos(n(#Psi-#Psi_{R}))>"));
    g->GetHistogram()->SetDirectory(0);
    g->Draw(Form("l %s same", option));
  }
}

//====================================================================
void 
AliFMDFlowResolutionTDR::Clear(Option_t*) 
{
  fN = 0;
  fLarge = 0;
}
//____________________________________________________________________
void 
AliFMDFlowResolutionTDR::Add(Double_t psi_a, Double_t psi_b)
{ 
  Double_t a = fabs(psi_a - psi_b);
  if (a >= .5 * M_PI) fLarge++;
  fN++;
}
//____________________________________________________________________
Double_t 
AliFMDFlowResolutionTDR::Correction(UShort_t k, Double_t& e2) const 
{ 
  // From nucl-ex/9711003v2 
  //
  //
  //   <cos n Delta phi> = 
  //
  //         sqrt(pi)
  //         -------- chi exp(-z) (I_((n-1)/2)(z) + I_((n+1)/2)(z))
  //            2
  // 
  // where z = chi^2 / 2
  //
  if (fLarge == 0) { 
    std::cerr << "TDR: K = 0" << std::endl;
    return -1;
  }
  if (fN == 0) { 
    std::cerr << "TDR: N = 0" << std::endl;
    return -1;
  }
  Double_t r     = Double_t(fLarge) / fN;
  Double_t echi2 = 0;
  Double_t y     = Chi2Over2(r, echi2);
  return Res(k, y, echi2, e2);
}
 
  

//____________________________________________________________________
Double_t 
AliFMDFlowResolutionTDR::Res(UShort_t k, Double_t y, Double_t echi2, 
			     Double_t& e2) const
{
  // y = chi^2 / 2
  Double_t chi   = sqrt(2 * y);
  Double_t c     = sqrt(M_PI) * exp(-y) / 2;

  // i[0] = I[(k-1)/2-1], i[1] = I[(k-1)/2], i[2] = I[(k-1)/2+1]
  Double_t i[3], di[3];
  AliFMDFlowBessel::Inu(Double_t(k-3)/2, Double_t(k+1)/2, y, i, di);
  
  Double_t r  = c * chi * (i[2] + i[1]);
  Double_t dr = c / 2 * ((4 - 2*k)*i[1] - (4*y + 2*k)*i[2] + 4*y*i[0]);
  e2        = dr * dr * echi2;
  return r;
}

//____________________________________________________________________
Double_t 
AliFMDFlowResolutionTDR::Correction(UShort_t k) const 
{ 
  Double_t e;
  return Correction(k, e); 
}
//____________________________________________________________________
Double_t 
AliFMDFlowResolutionTDR::Chi2Over2(Double_t r, Double_t& e2) const 
{
  // From nucl-ex/9711003v2 
  //
  // N equal to the total number of events, and k equal to the
  // number of events where |Psi_A - Psi_B| > pi/2, we get 
  //
  //   k   exp(-chi^2/2)                               k
  //   - = -------------     =>    chi^2 / 2 = - log 2 -
  //   N          2 		                       N
  //
  //                                               N
  //                         =>    chi^2 / 2 = log --   (*)
  //                                               k
  //            
  //                                                   2 k
  //                         =>    chi^2     = - 2 log ---
  //                                                    N
  //                                                            
  //                         =>    chi       = -/+ sgrt (- 2 log (2k/N))
  //
  // (*) this is what is returned.  
  //
  //  The error on chi is given by 
  //
  //              dchi                           1
  //    d^2chi = (------)^2 d^2(k/N) = - -------------------- d^2(k/N)
  //              d(k/N)                 4 (k/N)^2 log(2 k/N)
  // 
  // where 
  //                      k         k^2
  //               (1 - 2 -) dk^2 + --- dN^2
  //                      N	    N^2
  //   d^2(k/N) = --------------------------
  //                           N^2 
  // 
  // 
  // with dk = sqrt(k) and dN = sqrt(N), this reduces to 
  //
  //                     k      k^2        N k     k^2   k^2
  //              (1 - 2 -) k + ----- N    --- - 2 --- + ---
  //                     N      N^2         N       N     N
  //   d^2(k/N) = --------------------- = --------------------------
  //                           N^2                 N^2 
  //
  //              k N - k^2   k (N - k)   k (1 - k / N)
  //            = --------- = - ------- = - -----------
  //                 N^3      N   N^2     N      N
  //
  // Identifying r = k/N, we get 
  //
  //            
  //   d(k/N)   = sqrt(r (1 - r) / N) 
  //
  // which is the well-known result for Binomial errors.
  // Alternatively, one could compute these error using an Baysian
  // confidence level of say 68.3% (see for example 
  // http://home.fnal.gov/~paterno/images/effic.pdf). 
  // 
  // Our final expression for the error on chi then becomes 
  //
  //                      1	     
  //   d^2chi = - -------------------- r (1 - r) / N
  //		  4 r^2 log(2 r)
  //
  //                 1 - r            r - 1
  //          = - -------------- = ------------
  //              4 r N log(2 r)   4 k log(2 r)  
  if (r == 0) { 
    std::cerr << "TDR: Large/All = " << r << " <= 0!" << std::endl;
    return 0;
  }
  Double_t ratio = 1. / (2*r); // Double_t(fN) / (2 * fLarge);
  if (ratio <= 0) {
    std::cerr << "TDR: Large/All = " << ratio << " <= 0!" << std::endl;
    return -1;
  }
  Double_t chi2over2  = log(ratio);
  if (chi2over2 < 0) { 
    std::cerr << "TDR: log(" << ratio << ") = " << chi2over2 
	      << " < 0" << std::endl; 
    return -1;
  }
  if (fLarge != 0) e2 = (r - 1) / (4 * fLarge * log(2 * r));
  else             e2 = (r - 1) / (4 * r * log(2 * r));
  return chi2over2;
}

//____________________________________________________________________
void
AliFMDFlowResolutionTDR::Draw(Option_t* option) 
{
  TString opt(option);
  opt.ToLower();
  Bool_t chi = opt.Contains("chi");

  TH2* h = new TH2D(Form("tdr_%s_frame", (chi ? "chi" : "res")), 
		    Form("TDR %s Function for k=(1,2,4)", 
			 (chi ? "Chi" : "Resolution")),
		    100, 0, 1, 100, 0, (chi ? 3 : 1.5));
  h->SetXTitle("K/N");
  h->SetYTitle((chi ? "#chi" : "<cos(n(#Psi-#Psi_{R}))>"));
  h->SetStats(0);
  h->Draw();
  h->SetDirectory(0);

  for (UShort_t k = 1; k <= 4; k++) { 
    if (k == 3) continue;

    TGraphErrors* g = new TGraphErrors;
    Int_t i = 0;
    for (Double_t x = 0.02; x < 0.5; x += 0.01) { 
      Double_t e2 = 0;
      Double_t y  = 0;
      Double_t c  = Chi2Over2(x, e2);
      if (chi) y  = sqrt(2 * c);
      else { 
	Double_t dr = 0;
	y           = Res(k, c, e2, dr);
	e2          = dr;
      }
      g->SetPoint(i, x, y);
      g->SetPointError(i, 0, sqrt(fabs(e2)));
      i++;
    }
    g->SetLineColor(k);
    g->SetName(Form("tdr_%s_k%d", (chi ? "chi" : "res"), k));
    g->SetTitle(Form("TDR %s Function for k=%d", 
		     (chi ? "Chi" : "Resolution"), k));
    g->GetHistogram()->SetXTitle("K/N");
    g->GetHistogram()->SetYTitle((chi ? "#chi" : "<cos(n(#Psi-#Psi_{R}))>"));
    g->GetHistogram()->SetDirectory(0);
    g->Draw(Form("l %s %s", (k == 0 ? "h" : "same"), option));
    if (chi) break;
  }
}

//____________________________________________________________________
//
// EOF
//
