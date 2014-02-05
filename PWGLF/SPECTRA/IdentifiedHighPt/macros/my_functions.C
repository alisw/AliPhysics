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


class DeDxFitInfo : public TObject
{
 public:
  
  DeDxFitInfo();
  void Print(Option_t* option="") const;

  Double_t MIP;
  Double_t plateau;

  Int_t optionDeDx;
  Int_t nDeDxPar;
  Double_t parDeDx[8];

  Int_t optionSigma;
  Int_t nSigmaPar;
  Double_t parSigma[8];

  TString calibFileName;

  ClassDef(DeDxFitInfo, 1);    // Help class
};

//_____________________________________________________________________________
ClassImp(DeDxFitInfo)

DeDxFitInfo::DeDxFitInfo():
TObject(),
  MIP(0),
  plateau(0),
  optionDeDx(-1),
  nDeDxPar(-1),
  optionSigma(-1),
  nSigmaPar(-1),
  calibFileName("")
{
  // default constructor
  for(Int_t i = 0; i < 8; i++) {
    parDeDx[i]  = 0;
    parSigma[i] = 0;
  }
}

//_________________________________________________________
void DeDxFitInfo::Print(Option_t* option) const
{
  if(option) 
    cout << "Option: " << option << endl;

  cout << ClassName() << " : " << GetName() << endl  
       << "MIP: " << MIP << endl
       << "Plateau: " << plateau << endl
       << "OptionDeDx: " << optionDeDx << endl
       << "nDeDxPar: " << nDeDxPar << endl;
  for(Int_t i = 0; i < nDeDxPar; i++) {
    
    cout << "parDeDx[" << i << "] = " << parDeDx[i] << endl;
  }
  cout << "OptionSigma: " << optionSigma << endl
       << "nSigmaPar: " << nSigmaPar << endl;
  for(Int_t i = 0; i < nSigmaPar; i++) {
    
    cout << "parSigma[" << i << "] = " << parSigma[i] << endl;
  }
  
  if(calibFileName.IsNull()) {
    cout << "No eta calibration file." << endl; 
  } else {
    cout << "Eta calibration file: " << calibFileName.Data() << endl; 
  }
}

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
  
  

  Int_t option = TMath::Nint(par[0]);
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
  case 4: // just use bg
    bg = x[0];
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
  case 5: // fix the dE/dx to 50 at 0.5 GeV/c and the plateau to 75
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
      const Double_t e = par[4];

      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t powbg = TMath::Power(1.0+bg, c);

      const Double_t value = a/TMath::Power(beta2,e) + b/c*TMath::Log(powbg/(1.0 + d*powbg));          
      return value;
    }
    break;
  case 6:
    {
      /*
	a/beta^(e/2) + b/c*log( (1+x)^c / (1 + d*(1+x)^c) )

	Assymptotic behavior:

	1) Small bg (and d small so that d*(1+x)^c << 1)

	a/beta^(e/2) + b * log (1+x)
	
	So this is the same beavior as the standard expression. 

	2) Large bg where d*(1+x)^c >> 1
	a - b/c*log(d) = plateau
	-> d = exp(c*(a-plateau)/b)

	In this version we have 2 plateaus!!!!
	Plateau 1 = electrons!

       */
      if(specie==3)
	return fixPlateau;

      const Double_t a = par[1];
      const Double_t b = par[2];
      const Double_t c = par[3];
      const Double_t d = TMath::Exp(c*(a-par[5])/b);
      const Double_t e = par[4];
      
      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t powbg = TMath::Power(1.0+bg, c);

      const Double_t value = a/TMath::Power(beta2,e) + b/c*TMath::Log(powbg/(1.0 + d*powbg));          
      return value;
    }
    break;
  case 7: 
    {
      /*
	a/beta^(d/2) - b*log( c + 1.0/(1.0+x) )

	Assymptotic behavior:

	1) Small bg (and d small so that d*(1+x)^c << 1)

	a/beta^(d/2) - b * log (1+x)
	
	So this is the same beavior as the standard expression. 

	2) Large bg where c << 1
	-b*log(c) = plateau-a
	-> c = exp((a-plateau)/b)

       */
      const Double_t a = par[1];
      const Double_t b = par[2];
      const Double_t c = TMath::Exp((a-fixPlateau)/b);
      const Double_t d = par[3];
 
      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t value = a/TMath::Power(beta2,d) - b*TMath::Log(c + 1.0/(1.0+bg));          
      return value;
    }
  case 8: 
    {
      /*
	a/beta^(d/2) - b*log( c + 1.0/(1.0+x^e) )

	Assymptotic behavior:

	1) Small bg (and d small so that d*(1+x)^c << 1)

	a/beta^(d/2) - b * log (1+x^e)
	
	So this is the same beavior as the standard expression. 

	2) Large bg where c << 1
	-b*log(c) = plateau-a
	-> c = exp((a-plateau)/b)

       */
      const Double_t a = par[1];
      const Double_t b = par[2];
      const Double_t c = TMath::Exp((a-fixPlateau)/b);
      const Double_t d = par[3];
      const Double_t e = par[4];

      //      cout << bg << ": " << a << ", " << b << ", " << c << ", " << d << endl;

      const Double_t value = a/TMath::Power(beta2,d) - b*TMath::Log(c + 1.0/(1.0+bg) + e/(1.0+TMath::Power(bg, 2)));          
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
  case 6: // relative sigma with dE/dx to some power close to 1
    return par[1]*x[0]*TMath::Power(x[0]/50.0, par[2]);
    break;
  case 7: // 1/x^some power + constant
    return x[0]*(par[2]*TMath::Power(x[0], -par[3]) + par[1]);
    break;
  case 8: // for fitting relative sigma
    return par[2]*TMath::Power(x[0], -par[3]) + par[1];
    break;
  case 9: // for fitting relative sigma
    return par[1]+par[2]*x[0]+par[3]*x[0]*x[0];
    break;
  case 10: // for fitting relative sigma
    return par[1]+par[2]*x[0]+par[3]*x[0]*x[0]-par[4]/(x[0]*x[0]*x[0])-par[5]/(x[0]*x[0]);
    break;
  case 11: // for fitting relative sigma
    return x[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]-par[4]/(x[0]*x[0]*x[0])-par[5]/(x[0]*x[0]));
    break;
  case 12: // for fitting relative sigma
    return par[1]+par[2]*x[0]+par[3]*x[0]*x[0];
    break;
  case 13: // for fitting relative sigma
    return x[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]);
    break;
  case 14: // for fitting relative sigma
    return par[1]+par[2]*x[0];
    break;
  case 15: // for fitting relative sigma
    return x[0]*(par[1]+par[2]*x[0]);
    break;






  default:
    break;
  }
  
  cout << "Error in SigmaFunc: option " << option << " not supported!!!!!" << endl;
  return 0;
}
