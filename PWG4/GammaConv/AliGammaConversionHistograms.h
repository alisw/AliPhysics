#ifndef ALIGAMMACONVERSIONHISTOGRAMS_H
#define ALIGAMMACONVERSIONHISTOGRAMS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

#include "TString.h"
#include "Riostream.h"
#include <vector>

class TMap;
class TList;
class TH1F;
class TH2F;
class TH3F;

class AliGammaConversionHistograms{

 public: 
  
  AliGammaConversionHistograms();                                                                  //constructor
  AliGammaConversionHistograms(const AliGammaConversionHistograms & original);                     //copy constructor
  AliGammaConversionHistograms & operator = (const AliGammaConversionHistograms & original);       //assignment operator
  virtual ~AliGammaConversionHistograms();                                                         //virtual destructor
  

  //  TList * GetOutputContainer();
  void GetOutputContainer(TList *fOutputContainer);
  
  Int_t GetRBin(Double_t radius) const;
  Int_t GetPhiBin(Double_t phi) const;
  Int_t GetZBin(Double_t radius) const;
 
  void InitializeMappingValues(Int_t nPhiHistograms, Int_t nRHistograms, Int_t nBinsR, Double_t minRadius, Double_t maxRadius,Int_t nBinsPhi, Double_t minPhi, Double_t maxPhi);

  void AddMappingHistograms(Int_t nPhiHistograms, Int_t nRHistograms,Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle="", TString yAxisTitle="");

  /*
   * Adds a TH1F histogram to the histogram map and create a key for it 
   */
  void AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  /*
   * Adds a TH2F histogram to the histogram map and create a key for it 
   */
  void AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle="", TString yAxisTitle="", Int_t logAxis =-1);

  /*
   * Adds a TH3F histogram to the histogram map and create a key for it 
   */  
  void AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, Int_t nZBins, Double_t firstZ, Double_t lastZ, TString xAxisTitle="", TString yAxisTitle="", TString zAxisTitle="", Int_t logAxis =-1);

  /*
   * Create a logx binning suitable for dEdx plots
   */
  Bool_t BinLogAxis(const char* name, Int_t dim);


  /*
   *  Adds a TH1F Table to the table map and create a key for it
   */
  void AddTable(TString tableName,TString tableTitle,Int_t nXBins, const char * axesLabel[]);
    
 /*
  *  Adds a TH2F Table    
  */  
	
  void AddTable(TString tableName,TString tableTitle,Int_t nXBins, const char * axesXLabel[],Int_t nYBins, const char* axesYLabel[]);


  /*
  *  Adds a TH3F Table    
  */  
	
  void AddTable(TString tableName,TString tableTitle,Int_t nXBins, const char * axesXLabel[],Int_t nYBins, const char* axesYLabel[],Int_t nZBins, const char* axesZLabel[]);

  /*
   * Fills a TH1F histogram with the given name with the given value 
   */
  void FillHistogram(TString histogramName, Double_t xValue) const;

  /*
   * Fills a TH2F histogram with the given name with the given value 
   */
  void FillHistogram(TString histogramName, Double_t xValue, Double_t yValue) const;

  /*
   * Fills a TH3F histogram with the given name with the given value 
   */
  void FillHistogram(TString histogramName, Double_t xValue, Double_t yValue, Double_t zValue) const;

  /*
   * Fills a TH1F table with the given name with the given value
   */			
  void FillTable(TString tableName, Double_t xValue) const;

  /*  
   *  Fills a TH2F table with the given name with the given value
   */	
  void FillTable(TString tableName, Double_t xValue, Double_t yValue) const;

    /*  
   *  Fills a TH3F table with the given name with the given value
   */	
  void FillTable(TString tableName, Double_t xValue, Double_t yValue, Double_t zValue) const;

  /*
   *Returns a pointer to the histogram associated with name.
   */
   TObject* GetValue(const TString& name);


 private:
  TMap* fHistogramMap; // histogram map
  


  Int_t fNPhiIndex; //phi index
  Int_t fNRIndex; //r index
  Int_t fNZIndex; //z index
  //  Double_t fRBinLimits[8]; // Limits of the radius bins
  Double_t fRBinLimits[14]; // Limits of the radius bins
  Double_t fZBinLimits[13]; // Limits of the Z bins
  Double_t fMinRadius; //min radius cut
  Double_t fMaxRadius; //max radius cut
  Double_t fDeltaR; // delta r
  Double_t fMinPhi; //min phi
  Double_t fMaxPhi; // max phi
  Double_t fDeltaPhi;//delta phi

  TList * fMappingContainer; //mapping container
  TList * fBackgroundContainer; // background container
  TList * fDebugContainer; // debug container
  TList * fResolutionContainer; //resolution container
  TList * fMatchContainer; //match container
  TList * fESDContainer;//esd container
  TList * fMCContainer; // MC container
  TList * fTableContainer; // table container
  TList * fOtherContainer; // other container
  TList * f3DContainer; // 3D container
 
  ClassDef(AliGammaConversionHistograms,3)
};


#endif



