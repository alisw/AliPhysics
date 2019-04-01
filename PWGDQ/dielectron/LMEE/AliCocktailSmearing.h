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

class AliVParticle;

class AliCocktailSmearing {
  
public:
  AliCocktailSmearing(); ///< default constructor probably needed for AnalysisManager or such...
  virtual ~AliCocktailSmearing();
  
  void      	Print();
  void        ReadResoFile(TFile *fRes);
  void        SetSeed(UInt_t rndmseed);
  
  // setters
  void        SetResolutionP(TObjArray *resArr, Bool_t b=kFALSE)      { fPResArr=resArr; fUseRelPResolution=b; }
  void        SetResolutionOpeningAngle(TObjArray *resArr)            { fOpeningAngleResArr=resArr; }
  //void      SetArrayCopy(TObjArray *array);
  
  // getters
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
  
  ClassDef(AliCocktailSmearing,1)
};

#endif
