#ifndef AliCocktailSmearing_H
#define AliCocktailSmearing_H

//########################################################################
//#                                                                      #
//#   Class to smear pt,phi,eta for HF cocktail (mother analysis task)   #
//#                                                                      #
//#  Authors:                                                            #
//#   Raphaelle Bailhache, Uni Frankfurt / Raphaelle.Bailhache@cern.ch   #
//#   Carsten Klein, Uni Frankfurt / Carsten.Klein@cern.ch               #
//#   Jerome Jung, Uni Frankfurt / s8286523@uni-frankfurt.de             #
//#   Sebastian Scheid, Uni Frankfurt / s.scheid@cern.ch                 #
//#                                                                      #
//########################################################################

// c++ includes

// ROOT includes

// forward declarations
class TObjArray;
class TLorentzVector;
class TVector3;
class TH1F;
class TH2F;
class TH3F;

class AliVParticle;

class AliCocktailSmearing {
  
public:
  AliCocktailSmearing(); ///< default constructor probably needed for AnalysisManager or such...
  virtual ~AliCocktailSmearing();
  
  void      	Print();
  void        ReadResoFile(TFile *fRes);
  void        ReadEffFile(TFile *fEff);
  void        SetSeed(UInt_t rndmseed);
  void        SetEffType(Int_t typeeff)                               { fUseEffType = typeeff; }
  
  // setters
  void        SetResolutionP(TObjArray *resArr, Bool_t b=kFALSE)      { fPResArr=resArr; fUseRelPResolution=b; }
  void        SetResolutionOpeningAngle(TObjArray *resArr)            { fOpeningAngleResArr=resArr; }
  //void      SetArrayCopy(TObjArray *array);
  
  // getters
  Bool_t      GetTypeEfficiency()     { return fUseEffType; }
  Double_t    GetEfficiency(const TLorentzVector& vec);
  
  Bool_t      GetIsRelPResolution()   { return fUseRelPResolution; }
  TObjArray*  GetArrayP()             { return fPResArr; }
  TObjArray*  GetArrayOpeningAngle()  { return fOpeningAngleResArr; }
  Bool_t      GetSmearingOton()       { return fOton;               }
  
protected:
  Double_t        GetSmearing(TObjArray *arr, Double_t x);
  void            SmearOpeningAngle(TLorentzVector &lv1, TLorentzVector &lv2);
  Double_t        eMass()  { return 0.5109989e-3; }
  
  TLorentzVector  Smear(AliVParticle * part);
  TLorentzVector  ApplySmearingOton(const TLorentzVector& vec, short ch);
  
  // variables
  Int_t      fUseEffType;                                // Use to get the type of efficiency (1D, 2D, 3D)
  TH1F       *fheffPt;                                    // histogram eff as a function of pt
  TH2F       *fheffPtEta;                                 // histogram eff as a function of pt,eta
  TH3F       *fheffPtEtaPhi;                              // histogram eff as a function of pt,eta,phi
  Bool_t      fUseRelPResolution;                         // Use to get the type of resolution
  TObjArray  *fPResArr;                                   // Array of resolution histos  old Run 1
  TObjArray  *fOpeningAngleResArr;                        // Array of resolution histos  old Run 1
  TObjArray  *fArrResoPt;                                 // Array of resolution histos  Run 2
  TObjArray  *fArrResoEta;                                // Array of resolution histos  Run 2
  TObjArray  *fArrResoPhi_Pos;                            // Array of resolution histos  Run 2
  TObjArray  *fArrResoPhi_Neg;                            // Array of resolution histos  Run 2
  Bool_t      fOton;                                      // Use to get the type of resolution
  TVector3    fZaxis;                                     // z axis

  AliCocktailSmearing(const AliCocktailSmearing &c); // not implemented
  AliCocktailSmearing& operator= (const AliCocktailSmearing &c); // not implemented
  
  ClassDef(AliCocktailSmearing,2)
};

#endif
