#ifndef ALITRDSIM_H
#define ALITRDSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD simulation - multimodule (regular rad.)                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>

class TH1D;

class AliModule;

class AliTRDsim : public TObject {

 public:

  AliTRDsim();
  AliTRDsim(const AliTRDsim &s);
  AliTRDsim(AliModule *mod, Int_t foil, Int_t gap);
  virtual ~AliTRDsim();
  AliTRDsim &operator=(const AliTRDsim &s);

  virtual void     Copy(TObject &s) const;
  virtual void     Init();
  virtual Int_t    CreatePhotons(Int_t pdg, Float_t p
                               , Int_t &nPhoton, Float_t *ePhoton);
  virtual Int_t    TrPhotons(Float_t p, Float_t mass
                           , Int_t &nPhoton, Float_t *ePhoton);
  virtual Double_t Sigma(Double_t energykeV);
  virtual Double_t Interpolate(Double_t energyMeV
                             , Double_t *en, Double_t *mu, Int_t n);
  virtual Int_t    Locate(Double_t *xv, Int_t n, Double_t xval
                        , Int_t &kl, Double_t &dx);
  virtual Double_t Omega(Float_t rho, Float_t z, Float_t a)  { return (28.8 * TMath::Sqrt(rho * z / a)); };
  virtual Int_t    SelectNFoils(Float_t p);

          void     SetFoilThick(Float_t t)                   { fFoilThick = t;
                                                               SetSigma();                                  };
          void     SetGapThick(Float_t t)                    { fGapThick  = t;
                                                               SetSigma();                                  };
          void     SetFoilDens(Float_t d)                    { fFoilDens  = d; 
                                                               fFoilOmega = Omega(fFoilDens,fFoilZ,fFoilA);
                                                               SetSigma();                                  };
          void     SetFoilZ(Float_t z)                       { fFoilZ     = z; 
                                                               fFoilOmega = Omega(fFoilDens,fFoilZ,fFoilA); };
          void     SetFoilA(Float_t a)                       { fFoilA     = a; 
                                                               fFoilOmega = Omega(fFoilDens,fFoilZ,fFoilA); };
          void     SetGapDens(Float_t d)                     { fGapDens   = d;
                                                               fGapOmega  = Omega(fGapDens ,fGapZ ,fGapA );
                                                               SetSigma();                                  };
          void     SetGapZ(Float_t z)                        { fGapZ      = z;
                                                               fGapOmega  = Omega(fGapDens ,fGapZ ,fGapA ); };
          void     SetGapA(Float_t a)                        { fGapA      = a;
                                                               fGapOmega  = Omega(fGapDens ,fGapZ ,fGapA ); };
          void     SetTemp(Float_t t)                        { fTemp      = t; 
                                                               SetSigma();                                  };
          void     SetSigma();

  virtual Double_t GetMuPo(Double_t energyMeV);
  virtual Double_t GetMuCO(Double_t energyMeV);
  virtual Double_t GetMuXe(Double_t energyMeV);
  virtual Double_t GetMuBu(Double_t energyMeV);
  virtual Double_t GetMuMy(Double_t energyMeV);
  virtual Double_t GetMuN2(Double_t energyMeV);
  virtual Double_t GetMuO2(Double_t energyMeV);
  virtual Double_t GetMuHe(Double_t energyMeV);
  virtual Double_t GetMuAi(Double_t energyMeV);

          Float_t  GetFoilThick() const                      { return fFoilThick;     };
          Float_t  GetGapThick() const                       { return fGapThick;      };
          Float_t  GetFoilDens() const                       { return fFoilDens;      };
          Float_t  GetGapDens() const                        { return fGapDens;       };
          Double_t GetFoilOmega() const                      { return fFoilOmega;     };
          Double_t GetGapOmega() const                       { return fGapOmega;      };
          Float_t  GetTemp() const                           { return fTemp / 273.16; };
          TH1D    *GetSpectrum() const                       { return fSpectrum;      };

 protected:

          Int_t     fNFoilsDim;            //  Dimension of the NFoils array
          Int_t    *fNFoils;               //[fNFoilsDim] Number of foils in the radiator stack
          Double_t *fNFoilsUp;             //[fNFoilsDim] Upper momenta for a given number of foils
          Float_t   fFoilThick;            //  Thickness of the foils (cm)
          Float_t   fGapThick;             //  Thickness of the gaps between the foils (cm)

          Float_t   fFoilDens;             //  Density of the radiator foils (g/cm^3) 
          Float_t   fGapDens;              //  Density of the gas in the radiator gaps (g/cm^3)

          Double_t  fFoilOmega;            //  Plasma frequency of the radiator foils
          Double_t  fGapOmega;             //  Plasma frequency of the gas in the radiator gaps

          Float_t   fFoilZ;                //  Z of the foil material
          Float_t   fGapZ;                 //  Z of the gas in the gaps

          Float_t   fFoilA;                //  A of the foil material
          Float_t   fGapA;                 //  A of the gas in the gaps

          Float_t   fTemp;                 //  Temperature of the radiator gas (Kelvin)

          Int_t     fSpNBins;              //  Number of bins of the TR spectrum
          Float_t   fSpRange;              //  Range of the TR spectrum
          Float_t   fSpBinWidth;           //  Bin width of the TR spectrum
          Float_t   fSpLower;              //  Lower border of the TR spectrum
          Float_t   fSpUpper;              //  Upper border of the TR spectrum

          Double_t *fSigma;                //[fSpNBins] Array of sigma values

          TH1D     *fSpectrum;             //! TR photon energy spectrum

  ClassDef(AliTRDsim,2)                    //  Simulates TR photons

};
#endif
