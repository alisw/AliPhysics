#ifndef AliCumulantsFunctions_H
#define AliCumulantsFunctions_H

//********************************* 
// functions and equations needed * 
// for calcualation of cumulants  *
// and final flow estimates       *   
// author: Ante Bilandzic         * 
// email: anteb@nikhef.nl         *
//********************************* 

#include "AliFlowCommonConstants.h"
#include "AliFlowCumuConstants.h"

class TH1;
class TProfile;
class TProfile2D;
class TProfile3D;

class TObjArray;
class TList;
class TFile;

class AliCumulantsFunctions{
 public:
  AliCumulantsFunctions();
  virtual ~AliCumulantsFunctions();
  AliCumulantsFunctions(TProfile2D *IntGen, TProfile3D *DiffGenRe, TProfile3D *DiffGenIm, TH1D *ifr, TH1D *dfr2, TH1D *dfr4, TH1D *dfr6, TH1D *dfr8, Double_t CvM);
  void Calculate();

 private:
  AliCumulantsFunctions(const AliCumulantsFunctions& fun);
  AliCumulantsFunctions& operator=(const AliCumulantsFunctions& fun);
  
  TProfile2D *fIG;   //
  TProfile3D *fDGRe; //
  TProfile3D *fDGIm; //
  TH1D *fifr;       //integrated flow final results 
  TH1D *fdfr2;      //differential flow final results //to be improved
  TH1D *fdfr4;      //differential flow final results //to be improved
  TH1D *fdfr6;      //differential flow final results //to be improved
  TH1D *fdfr8;      //differential flow final results //to be improved
  Double_t fAvMult; //avarage multiplicity
  
  ClassDef(AliCumulantsFunctions, 0);
};
#endif





