#ifndef ALITPCD_H
#define ALITPCD_H

// include files and class forward declarations


// Documentation
/**
  *   Object with the  TPC parameters 
  *   
  *
  *
  */
#include "TNamed.h"
#include "TTree.h"
class AliTPCPRF2D;
class AliTPCRF1D;
class AliTPCParam;
class TClonesArray;
class TTree;
class TDirectory;
R__EXTERN TDirectory *  gDirectory;
//class TBranch;

class AliTPCD : public TNamed{

public:
  AliTPCD(
	  Text_t *  name ="DIGIT",
	  AliTPCParam *param=0, 
	  AliTPCPRF2D *prf=0, 
	  AliTPCRF1D *prfz=0);
  ~AliTPCD();
  
public: 
 
  AliTPCParam & GetParam() {return *fParam;}
  //give us reference to the parameters
  AliTPCPRF2D &  GetPRF2D() {return *fPRF;}
  //give us reference to 2-dimensional pad response function object
  //this is responsible for respnse in x and y dimension
  AliTPCRF1D  &  GetRF() {return *fRF;}
  //give us reference to 1 dimensionl response
  //this is responsible for z-direction
  TClonesArray * GetArray() {return fDigits;}
  //return reference for digits array

  TTree * GetTree() { return fTreeD;}
  //return refeence to actual tree 
  Bool_t  SetTree(Int_t nevent=0, TDirectory *dir = gDirectory);
  //map tree from given directory
  Bool_t  MakeTree(Int_t nevent=0);
  //map tree from given directory
   void Fill();
 
protected:
  
public:
  AliTPCParam * fParam;     //geometry and gas parameters object 
  AliTPCPRF2D * fPRF;            //x and y pad response function object
  AliTPCRF1D * fRF;             //z (time) response function object
  TClonesArray *fDigits;    //aray of tpc Digits
  TTree   *fTreeD;        //tree 
private:
  AliTPCD *fpthis;  //pointer to object
  ClassDef(AliTPCD,2) 
};


#endif /* ALITPCD_H */
