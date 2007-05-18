#ifndef ALITREEDRAW_H
#define ALITREEDRAW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliTreeDraw                               //
//                                                                          //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////



#include <TObject.h>
#include <TObjArray.h>

class TH1;
class TH1F;
class TH2F;
class TTree;

class AliTreeDraw: public TObject{
public:
  AliTreeDraw();
  ~AliTreeDraw(){;}
  void SetTree(TTree *tree){fTree=tree;}
  const TH1 * GetRes() const{ return (TH1*)fRes;}
  const TH1 * GetMean() const{ return (TH1*)fMean;}
  const TObjArray* GetPoints() const {return fPoints;}
  void  ClearHisto();
  void  ClearPoints(){if (fPoints) fPoints->Clear();}
  //
  TH1F * DrawXY(const char * chx, const char *chy, const char* selection, 
		const char * quality,Int_t nbins, Float_t minx, Float_t maxx, 
		Float_t miny, Float_t maxy, Int_t nBinsRes=30);
  TH1F * DrawLogXY(const char * chx, const char *chy, const char* selection, 
		   const char * quality, Int_t nbins,Float_t minx, Float_t maxx, 
		   Float_t miny, Float_t maxy, Int_t nBinsRes=30); 
  TH1F * Eff(const char *variable, const char* selection, const char * quality, 
	     Int_t nbins,Float_t min, Float_t max); 
  TH1F * EffLog(const char *variable, const char* selection, const char * quality, 
	     Int_t nbins,Float_t min, Float_t max);
  //
  void   GetPoints3D(const char * label, const char * chpoints, const char* selection, TTree * tpoints, Int_t color=6, Float_t rmin=4.);

 private: 
  static void   AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle);
  static Double_t* CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax);
  static TH1F*  CreateEffHisto(TH1F* hGen, TH1F* hRec);
  static TH1F*  CreateResHisto(TH2F* hRes2, TH1F **phMean, 
				Bool_t drawBinFits = kTRUE,Bool_t overflowBinFits = kFALSE);

  TTree * fTree;    //the tree for visualization - NOT OWNER
  TH1F  * fRes;     //temporary histogram        - OWNER  
  TH1F  * fMean;    //temporary histogram        - OWNER
  TObjArray *fPoints;//                          - OWNER
  ClassDef(AliTreeDraw,0)
};

#endif
