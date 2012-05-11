#include <TMath.h>
#include <TH2.h>
#include <TProfile.h>

#include <iostream>

using namespace std;

const Double_t eMass  = 0.000511;               //electron mass
const Double_t piMass = 0.13957;               //pion mass
const Double_t kMass  = 0.493676999999999977;  //kaon mass
const Double_t pMass  = 0.938271999999999995;  //proton mass


Double_t fitf3G(Double_t* xx, Double_t* par);
Double_t FitFunc(Double_t* x, Double_t *par);
Double_t SigmaFunc(Double_t* x, Double_t *par);

TH2D* hDeDxVsP = 0;
TProfile* hMeanP = 0;
Double_t fixMIP     = 50.0;
Double_t fixPlateau = 75.0;

//______________________________________________________________________________
Double_t fitf3G(Double_t* xx, Double_t* par)
{
  //
  // Could speed up fit by forcing it to use <p>. In that way the parameters
  // could be amde statis cand only changed when going to a new p bin
  //
  Double_t p = xx[0];
  Double_t dedx = xx[1];

  const Int_t bin = hDeDxVsP->GetXaxis()->FindBin(p);
  
  if(hMeanP) {
    //    cout << "p before: " << p;
    p = hMeanP->GetBinContent(bin);
    //    cout << ", p after: " << p << endl;
  }
  const Int_t binStart = Int_t(par[0]);
  const Int_t binStop  = Int_t(par[1]);

  if(bin<binStart || bin>binStop) {
    
    cout << "Error: bin " << bin << " not inside inteval [" << binStart 
	 << "; " << binStop << "]" << endl;
    return 0;
  }

  const Int_t nParDeDx = Int_t(par[2]);
  
  Double_t* parDeDx = &par[3];

  Int_t offset = 4 + nParDeDx; // binStart + binStop + nParDeDx + optionDedx + nParDeDx parameters

  const Int_t nParSigma = Int_t(par[offset]);
  offset += 1; // nParSigma

  Double_t* parSigma = &par[offset];
  offset += 1 + nParSigma; // optionSigma + nParSigma parameters
  
  Double_t piMean = FitFunc(&p, parDeDx);
  Double_t pKeff = p*piMass/kMass; // corresponding p of a pion with same dE/dx
  Double_t kMean  = FitFunc(&pKeff, parDeDx);
  Double_t pPeff = p*piMass/pMass; // corresponding p of a pion with same dE/dx
  Double_t pMean  = FitFunc(&pPeff, parDeDx);

  const Double_t piSigma = SigmaFunc(&piMean, parSigma);
  const Double_t kSigma  = SigmaFunc(&kMean,  parSigma);
  const Double_t pSigma  = SigmaFunc(&pMean,  parSigma);


  const Int_t j = bin - binStart;
  const Double_t piYield   = par[j * 3 + offset + 0];
  const Double_t kYield    = par[j * 3 + offset + 1];
  const Double_t pYield    = par[j * 3 + offset + 2];
  
  return piYield* TMath::Gaus(dedx, piMean, piSigma, kTRUE)
    +    kYield * TMath::Gaus(dedx, kMean,  kSigma,  kTRUE)
    +    pYield * TMath::Gaus(dedx, pMean,  pSigma,  kTRUE);
}


//______________________________________________________________________________
Double_t FitFunc(Double_t* x, Double_t *par)
{
  static const Double_t bgMIP    = 0.5/piMass;
  static const Double_t beta2MIP = bgMIP*bgMIP / (1.0+bgMIP*bgMIP);
  //  static const Double_t betapowMIP = TMath;
  static const Double_t logMIP   = TMath::Log(1+bgMIP);
  
  Int_t option = Int_t(par[0]);
  Int_t specie = option;
  option = option%10;
  specie -= option;
  specie /= 10;


  Double_t bg = 0;
  switch (specie) {
    
  case 0: // pion
    bg = x[0]/piMass;
    break;
  case 1: // kaon
    bg = x[0]/kMass;
    break;
  case 2: // proton
    bg = x[0]/pMass;
    break;
  case 3: // electron
    bg = x[0]/eMass;
    break;
  default:
    cout << "Error in FitFunc: specie " << specie << " not supported!!!!!" << endl;
    return 0;
    break;
  }
    
  if(bg > 10000.0)
    bg = 10000.0;

  const Double_t beta2 = bg*bg / (1.0+bg*bg);
  
  switch (option) {
    
  case 1: // standard parametrisation
    {
      /*
	c0/beta^2 + c1 * log (1+x)
       */
      const Double_t c0 = par[1];
      const Double_t c1 = par[2];
      
      const Double_t value = c0/beta2 + c1*TMath::Log(1+bg);          
      return value;
    }
    break;
  case 2: // fix the dE/dx to 50 at 0.5 GeV/c
    {
      const Double_t c1 = par[1];
      const Double_t c0 = (fixMIP-par[1]*logMIP) * beta2MIP;
      
      const Double_t value = c0/beta2 + c1*TMath::Log(1+bg);          
      return value;
    }
    break;
  case 3: // fix the dE/dx to 50 at 0.5 GeV/c and the plateau to 75
    {
      /*
	a/beta^2 + b/c*log( (1+x)^c / (1 + d*(1+x)^c) )

	Assymptotic behavior:

	1) Small bg (and d small so that d*(1+x)^c << 1)

	a/beta^2 + b * log (1+x)
	
	So this is the same beavior as the standard expression. 

	2) Large bg where d*(1+x)^c >> 1
	a - b/c*log(d) = plateau
	-> d = exp(c*(a-plateau)/b)

       */
      const Double_t b = par[1];
      const Double_t a = (fixMIP-par[1]*logMIP) * beta2MIP;
      const Double_t c = par[2];
      const Double_t d = TMath::Exp(c*(a-fixPlateau)/b);

      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t powbg = TMath::Power(1.0+bg, c);

      const Double_t value = a/beta2 + b/c*TMath::Log(powbg/(1.0 + d*powbg));          
      return value;
    }
    break;
  case 4: // fix the dE/dx to 50 at 0.5 GeV/c and the plateau to 75
    {
      /*
	a/beta^2 + b/c*log( (1+x)^c / (1 + d*(1+x)^c) )

	Assymptotic behavior:

	1) Small bg (and d small so that d*(1+x)^c << 1)

	a/beta^2 + b * log (1+x)
	
	So this is the same beavior as the standard expression. 

	2) Large bg where d*(1+x)^c >> 1
	a - b/c*log(d) = plateau
	-> d = exp(c*(a-plateau)/b)

       */
      const Double_t a = par[1];
      const Double_t b = par[2];
      const Double_t c = par[3];
      const Double_t d = TMath::Exp(c*(a-fixPlateau)/b);

      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t powbg = TMath::Power(1.0+bg, c);

      const Double_t value = a/beta2 + b/c*TMath::Log(powbg/(1.0 + d*powbg));          
      return value;
    }
    break;
    // case 3: // fix the dE/dx to 50 at 0.5 GeV/c + powerlaw
    
    //   static const bgMIP    = 0.5/piMass;
    //   static const beta2MIP = bgMIP*bgMIP / (1.0+bgMIP*bgMIP);
    //   static const logMIP   = TMath::Log(1+bgMIP);
    
    //   const Double_t c1 = par[1];
    //   const Double_t c2 = par[2]; // expect it to be 0.75 from Bichsel (beta**-1.5 instead of -2)
    //   const Double_t c0 = (50.0-par[1]*logMIP) * beta2MIP;
    
    //   const Double_t value = TMathh::Power(c0/beta2, c2) + c1*TMath::Log(1+bg);          
    //   return value;
    
  default:
    break;
  }
  
  cout << "Error in FitFunc: option " << option << " not supported!!!!!" << endl;
  return 0;
}

//______________________________________________________________________________
Double_t SigmaFunc(Double_t* x, Double_t *par)
{
  Int_t option = Int_t(par[0]);
  
  switch (option) {
  case 1: // fixed sigma
    return par[1];
    break;
  case 2: // relative sigma
    return par[1]*x[0];
    break;
  case 3: // relative sigma + extrapolation
    return (par[1] + (x[0]-fixMIP)*par[2])*x[0];
    break;
  case 4: // relative sigma with dE/dx to some power close to 1
    return par[1]*TMath::Power(x[0], par[2]);
    break;
  case 5: // relative sigma with dE/dx to some power close to 1
    return TMath::Sqrt(par[1]*par[1]*x[0]*x[0] + par[2]*par[2]);
    break;
  default:
    break;
  }
  
  cout << "Error in SigmaFunc: option " << option << " not supported!!!!!" << endl;
  return 0;
}
