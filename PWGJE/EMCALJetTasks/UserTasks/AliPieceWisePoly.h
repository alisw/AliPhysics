#ifndef ALIPIECEWISEPOLY_H
#define ALIPIECEWISEPOLY_H

class TF1;

class AliPieceWisePoly {
public:
  AliPieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t xmin = 0, Double_t xmax = 1, Double_t* params = 0x0, Int_t smooth = 2);
  ~AliPieceWisePoly();
  double operator() (double* x, double* p = 0x0) {return Eval(x[0], p);};
  double Eval (double x, double* p = 0x0);
  void SetParam(Double_t* params);
  Double_t SumUp(Double_t constant, Int_t start, Int_t end, Int_t derivative, Int_t startn);
  Int_t GetNOfParam() const {return nOfFreeParams;};
  TF1* GetFunction() {return piecewisepolynom;};
  TF1* GetPartFunction(Int_t i);
  Int_t GetNParts() {return nParts;};
  Double_t* GetCuts() {return cuts;};
  Int_t* GetNParameters() {return polyParameters;};
  static TString ReadFSParameters(TString parameterFile, TF1** effFunctions);                             // Reads parameters from a fastSimParameter file into functions, additionally returns a string containing all parameters together
  static void ReadFSParametersFromString(TString parameterString, TF1** effFunctions);  // Reads parameters from a string
  static TF1* GetFSFunctionFromString(TString parameters, TString nameFunction);           // Gets function from parameter String
  
private:
  AliPieceWisePoly(const AliPieceWisePoly&)           ;  //Not implemented
  AliPieceWisePoly &operator=(const AliPieceWisePoly&);  //Not implemented
  
  Int_t doSmoothing;
  Int_t nParts;
  Double_t* cuts;
  Int_t nOfFreeParams;
  Int_t* polyParameters;
  TF1* piecewisepolynom;
};

#endif

