#ifndef TTPCDIGITSH_H
#define TTPCDIGITSH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliTPCDigitsH             //
////////////////////////////////////////////////
 

// include files and class forward declarations

class TFile;   
class TH2F;
class AliH2F;
class TH1F;
class TPaveText;
class TPad;
class TCanvas;
class AliTPCParam;
class AliTPCD;
class TClonesArray;
class TTree;
#include "TObject.h"
#include "TUnixSystem.h"

class AliTPCDigitsH :public TObject{
public:   
  ~AliTPCDigitsH();
  AliTPCDigitsH(const AliTPCDigitsH &);
  AliTPCDigitsH();
  AliTPCDigitsH & operator = (const AliTPCDigitsH &);
  AliTPCD * GetDParam() {return fDParam;}
  void   SetDParam(AliTPCD * dig);
  AliTPCParam *&   GetParam();
  //get tpc parameters
  TFile * GetIn() {return fin;}  
  TFile * GetOut() {return fout;}
  void Anal();   //get two dimensional histogram time*pad*amplitude for given
  void SetSecRowTime(Int_t sec , Int_t row = 1, 
		   Int_t TimeN = 500 , 
		   Float_t TimeStart = 0, Float_t TimeStop = 500); 
  //set section row and timing histogram 
  //ranges and calculate histograms 
  void SetHisto( Int_t pad = 1 ); 
  virtual void Draw(Option_t * opt1 ="cont1"  , Option_t * opt2 = "error",
			Option_t * opt3 = "L" );
  
  void CloseFiles();  //close input and output file
  void DeleteHisto(const Text_t *namecycle);  //
 
  Bool_t  SetIO(const char *  inFile, const char* outFile);  
  Bool_t  SetTree(Int_t eventn =0 );
  void  SetParticles(Int_t sec = -1, Int_t row = -1 ,
  		   Int_t size1 = 30000,Int_t size2=300,
		   Bool_t all=1 );
  Bool_t SetEventN(Int_t EventN=0);  //set event number connected to analisys
  
  void SetbDelHisto(Bool_t bDelHisto) { fbDelHisto = bDelHisto;} //
  void SetThreshold(Int_t threshold = 10 ) { fThreshold = threshold ;} //
  /////////////  
  AliH2F &  GetHis1() { return  *fH2Digit;}
  AliH2F & GetHisBW() { return  *fH2DigitBW;}
  TH1F &  GetHis2() { return  *fH1Digit;}
  TH1F &  GetHis3() { return  *fH1Occu;}
  TPad & GetPad1() {return *fpad1;}   
  TPad & GetPad2() {return *fpad2;}   
  TPad & GetPad3() {return *fpad3;}  
  TCanvas &GetCanvas(){return *fcanvas;}
protected:

private:
  TFile *fin;   //input TTRE file with digits
  TFile *fout; //output file
  TPad * fpad1;
  TPad * fpad2;
  TPad * fpad3;
  TCanvas * fcanvas;
  AliTPCParam * fTPCParam;
  AliTPCD *fDParam;
  Bool_t fbIOState;  
  //state of Input and Output file : kTRUE if all is OK
  Bool_t fbDigState; 
  //state of Input digits tree kTRUE if ftree point to valid digit tree
  TClonesArray  *fDigits;  //pointer to digits object used inanal 
  TClonesArray  *fParticles; //pointer to particles array 
  TTree   *ftree;  //pointer to tree with digits
  Int_t fsec;  //actual sector of TPC  
  Int_t frow;  //actual row of TPC 
  Int_t fpad;  //actual pad
  Int_t fEventN; //event number connected to
  Int_t fTimeN;  //division number of the time
  Float_t fTimeStart,fTimeStop;  //start and stop time for histograms
  Int_t fOccuN;  //time division for occupancy calculating
  TPaveText *  fTitle;  //  
  Int_t fThreshold;   //treshold for the ocupancy definition
  Bool_t fbDelHisto;   //swith which tell if we want to delete produced histogram
  
  AliH2F *fH2Digit; //two dimensional histogram - sector row
  AliH2F *fH2DigitBW; //black and white histogram over treshold
  TH1F *fH1Digit; //one dimensional histogram sector row pad
  TH1F *fH1Occu; //occupancy graph  
  TH1F *fHParticles; //histogram of particles contributing to digits  
  TH1F *fHPartMultiplicity; 
  //histogram of particles multiplicity
  TH1F *fHAllP;  //histogram with accepted time bin for all  particles  
  TH1F *fHSecondaryP;
  //histogram with accepted time bin for secondary particles  
  
 ClassDef(AliTPCDigitsH,1)
};
#endif /* TTPCDIGITSH_H */
