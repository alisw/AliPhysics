/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/

#include <stdlib.h>

#include "TH1.h"
#include "TRandom.h"
#include "TMath.h"

#include "AliTRDsim.h"
#include "AliTRDconst.h"

ClassImp(AliTRDsim)


const Float_t kD1 = kPeThick / kRaFoils;
const Float_t kD2 = kRaThick / (kRaFoils + 1);

//Root specials, to be removed

static TH1F *h100, *h101, *h102;

AliTRDsim::AliTRDsim()
{
  fNj=200;
  fIrst=1;
}

AliTRDsim::AliTRDsim(AliModule* mod, Int_t foil, Int_t gas)
{
  Float_t a1, z1, ro1, rad, abs;
  Float_t a2, z2, ro2;
  char * name[21];
  mod->AliGetMaterial(foil, name, a1, z1, ro1, rad, abs);
  mod->AliGetMaterial(gas, name, a2, z2, ro2, rad, abs);
  fOmega1 = 28.8*TMath::Sqrt(ro1*z1/a1);
  fOmega2 = 28.8*TMath::Sqrt(ro2*z2/a2);
}

void AliTRDsim::trd_sim()
{

  const Float_t amass[4] = { 5.11e-4,.13957,.4937,.10566 };
  const Double_t of[4] = { 20.9,24.4,14.27,26.9 };
  const Double_t og[4] = { .28,.7,.74,.74 };
  Int_t ifl = 0;
  Int_t ig = 1;
  Int_t nev = 1000;
  Double_t gamma = -10.;
  
  /* Local variables */
  static Float_t temp, pres;
  static Int_t i, j;
  static Float_t o, sigma[200];
  static Float_t trEn[10];
  static Double_t omega1, omega2;
  static Float_t am;
  static Int_t np;
  static Int_t ipa;
  
  /* ***********************************************************************
   */
  /*  TRD simulation - multimodule (regular rad.) */
  /*     after: M. CASTELLANO et al., */
  /*   COMP. PHYS. COMM. 51 (1988) 431 + COMP. PHYS. COMM. 61 (1990) 395 */
  
  /*   17.07.1998 - A.Andronic */
  /*   08.12.1998 - simplified version */
  
  ipa = 0;
  /* that's electron */
  am = amass[ipa];
  omega1 = of[ifl];
  /* plasma frequency: foil and gap */
  omega2 = og[ig - 1];
  if (gamma < -1e5) printf("*** Momentum steps !!! ***\n");
  if (gamma < 0. && gamma >= -1e5) {
    gamma = sqrt(gamma * gamma + am * am) / am;
    printf("*** Gamma (electron) = %f\n",gamma);
  }
  temp = 20.;
  pres = 1.;
  fBin = 100. /  fNj;
  /* binsize */
  fL = 1. - fBin / 2.;
  fU = fL + 100.;
  /*  setting the stage ................................... */
  for (j = 0; j < fNj; ++j) {
    /* getting the sigma values - for fixed energy values */
    o = fBin *  j + 1.;
    /* omega in keV */
    /* abs. in rad. (1 foi */
    sigma[j] = fsigmaRad(ifl, ig, o);
  }
  printf(" Working...\n");
  /*  sampling over some events ........................... */
  for (i = 0; i < nev; ++i) {
    xtr(gamma, omega1, omega2, ro1, ro2, sigma, np, trEn);
    /* TR: n, E */
    h101->Fill(np);
    /* sample nTR distr. */
    for (j = 0; j < np; ++j) {
      h102->Fill(trEn[j], 1. / fBin);
      /* sample the TR en. distr. */
    }
  }
  /* ------------------------------------------------------------------- */
  /*      else  !ns steps */
  /*      enddo  !imod */
  /* events */
  h100->Draw();
  h101->Draw();
  h102->Draw();
} /* trd_sim__ */

void AliTRDsim::xtr(Double_t gamma, Double_t omega1, Double_t omega2, Double_t ro1,
		    Double_t ro2,
		     Float_t *sigmaRad, Int_t &np, Float_t *trEn)
{
    /* Initialized data */

  static Double_t alfa = .0072973;
  static Double_t pi = 3.14159265;
  
  /* Local variables */
  static Double_t conv, a;
  static Int_t i, j;
  static Float_t o, w[200], omega;
  static Double_t tetan, stemp;
  static Float_t om;
  static Double_t sk;
  static Float_t wn[200];
  static Double_t cs1, cs2;
  static Double_t ro11, ro22, aux;
  static Float_t ntr;
  static Double_t sum;
  
  /************************************************************************
   ******/
  /*   TR: number and energy distr. */
  
  /* Function Body */
  sk = kD2 / kD1;
  /* -------------- starts with the TR spectrum ------------- */
  
  stemp = 0.;
  for (j = 0; j < fNj; ++j) {
    /* TR spectrum */
    omega = (fBin *  j + 1.) * 1e3;
    /* keV->eV */
    cs1 = omega1 / omega;
    cs2 = omega2 / omega;
    ro11 = omega * kD1 * 2.5 * (1. / (gamma * gamma) + cs1*cs1);
    ro22 = omega * kD1 * 2.5 * (1. / (gamma * gamma) + cs2*cs2);
    sum = 0.;
    for (i = 0; i < 10; ++i) {
/* 30 - it matters a bit */
      tetan = (pi * 2. *  (i+1) - (ro11 + sk * ro22)) / (sk + 1.);
      if (tetan < 0.) {
	tetan = 0.;
      }
      aux = 1. / (ro11 + tetan) - 1. / (ro22 + tetan);
      a = tetan * (aux * aux) * (1. - cos(ro11 + tetan));
      sum += a;
    }
    o = omega * .001;
    /* eV->keV */
    conv = 1. - exp(-kRaFoils * sigmaRad[j]);
    w[j] = alfa * 4. / (sigmaRad[j] * (sk + 1.)) * conv * sum;
    /* dW/domega */
    wn[j] = w[j] / o;
    /* dN/domega */
    stemp += wn[j];
    if (fIrst == 1) {
      h100->Fill(o, wn[j]);
      /* double precision not accepted */
    }
  }
  /* -------------- done with the spectrum ------------- */
  /* j (omega spectrum) */
  ntr = stemp * fBin;
  /* <nTR> (binsize corr.) */
  om = h100->GetMean();
  /* <Etr> */
  if (fIrst == 1) {
    /* prints the production */
    printf(" Produced TR - <n>, <E>: %5.2f  %6.2f  KeV\n",ntr,om);
    fIrst = 0;
  }
  /* prob. distr. */
  np = gRandom->Poisson(ntr);
  /* Np TR photons Poiss distr. from mean */
  for (j = 0; j < np; ++j) {
    /* TR energy (binsize corr.) */
    trEn[j] = hisran(wn, fNj, fL, fBin);
    /* their energy */
  }
}

Float_t AliTRDsim::fsigmaRad(Float_t ro1, Float_t ro2, Int_t ig, Float_t o)
{

  /* Local variables */
  static Float_t pres;
  static Double_t mumu;
  static Int_t j;
  static Double_t t;
  static Int_t i1, i2;
  static Double_t x1;
  static Double_t mu1, mu2, deo, omf[36], omg[36], muf[36], mug[36];
  
  static Bool_t first = kTRUE;
  
  /* cccccccccccccccccccccccccccccccccccccccccccc */
  /*    calculates sigma for radiator - one foil+one gap */
  
  if(first) {
    FILE* inp = fopen("po.x","r");
    for (j=0;j<36;++j) {
      fscanf(inp,"%lf %lf %lf",&omf[j],&muf[j],&mumu);
    }
    fclose(inp);
    inp = fopen("he.x","r");
    for (j=0;j<36;++j) {
      fscanf(inp,"%lf %lf %lf",&omg[j],&mug[j],&mumu);
    }
    fclose(inp);
    first=kFALSE;
  }
  /* first */
  x1 = o * .001;
  /* keV->MeV */
  if (x1 >= .001) {
    locate(omf, 36, x1, i1, deo);
    mu1 = muf[i1] - deo * (muf[i1] - muf[i1+1]) / (omf[i1+1] - omf[i1]);
    locate(omg, 36, x1, i2, deo);
    mu2 = mug[i2] - deo * (mug[i2] - mug[i2+1]) / (omg[i2+1] - omg[i2]);
    t = 273.16;
    /* gases at 0 C */
    return (mu1*ro1*kD1+mu2*293.16/t * ro2*kD2)/1e4;
    /* mu */
  } else {
    return 1e6;
  }
} 

Int_t AliTRDsim::locate(Double_t *xv, Int_t n, Double_t xval, 
		       Int_t &kl, Double_t &dx)
{
  /* -------------------------------------------------------------- */
  /*  locates a point (xval) in a 1-dim grid (xv(n)) --> iloc,dx,ier */
  /* -------------------------------------------------------------- */
  /* Function Body */
  if (xval >= xv[n-1]) return 1;
  if (xval < xv[0]) return -1;
  Int_t km,kh=n-1;
  kl=0;
  while(kh-kl>1) if(xval<xv[km=(kl+kh)/2]) kh=km; else kl=km;
  if(xval<xv[kl] || xval > xv[kl+1] || kl >= n-1) {
    printf("locate failed xv[%d] %f xval %f xv[%d] %f!!!\n",
	   kl,xv[kl],xval,kl+1,xv[kl+1]);
    exit(1);
  }
  dx=xval-xv[kl];
  return 0;
}

Float_t AliTRDsim::hisran(Float_t *y, Int_t n, Float_t xlo, Float_t xwid)
{
    /* Local variables */
    Float_t yinv, ytot=0;
    Int_t i;
    Float_t yr;

/*         SUBROUTINE TO GENERATE RANDOM NUMBERS */
/*         ACCORDING TO AN EMPIRICAL DISTRIBUTION */
/*         SUPPLIED BY THE USER IN THE FORM OF A HISTOGRAM */
/*         F. JAMES,    MAY, 1976 */

    if (y[n-1] != 1.) {

/*         INITIALIZE HISTOGRAM TO FORM CUMULATIVE DISTRIBUTION */

      ytot = 0.;
      for (i = 0; i < n; ++i) {
	if (y[i] < 0.) {
	  printf("hisran: found value y[%d] = %f\n",i,y[i]);
	  exit(1);
	}
	ytot += y[i];
	y[i] = ytot;
      }
      if (ytot <= 0.) {
	printf("hisran: total probability %f < 0\n",ytot);
	exit(1);
      }
      yinv = 1. / ytot;
      for (i = 0; i < n-1; ++i) {
	y[i] *= yinv;
      }
      y[n-1] = 1.;
    }
/*         NOW GENERATE RANDOM NUMBER BETWEEN 0 AND ONE */
    yr = gRandom->Rndm();
/*         AND TRANSFORM IT INTO THE CORRESPONDING X-VALUE */
    if(yr<=y[0]) return xlo + xwid * (yr / y[0]);
    else {
      Int_t km,kl=0,kh=n-1;
      while(kh-kl>1) if(yr<y[km=(kl+kh)/2]) kh=km; else kl=km;
      return xlo + xwid * (kl + (yr - y[kl]) / (y[kl + 1] - y[kl]));
    }
}

