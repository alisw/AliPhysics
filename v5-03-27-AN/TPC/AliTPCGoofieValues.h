/*!\class AliTPCGoofieValues
   \brief TPC calibration class for Goofie values 

   Header: AliTPCGoofieValues.h,v 2.0.

   Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.

   See cxx source for full Copyright notice 

   TPC Calibration Class for GOOFIE values. Drift velocity, gas composition and the gain. 

   The class AliTPCGoofieValues allows the access to GoofieValues. 

   The only constructor is loading data from ASCI file. 

   The methods make Tgraphs and TSplines of the time dependence of the values. 

   One method allows save the the graphs and spline togather with tree of allvalues into file.        

*/


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for Goofie values                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TSystem.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include "TTree.h"
#include "TSpline.h"

class AliTPCGoofieValues : public TNamed{
public:
  AliTPCGoofieValues(); ///< default ctor
  //AliTPCGoofieValues(const char *fname); ///< ctor using log file, not implemented
  AliTPCGoofieValues(const char *fname); ///< ctor using an ASCII file
  virtual ~AliTPCGoofieValues();///< default dtor
  
  
  Long64_t GetLinesInFile();///< return lines in ASCII file 
  Double_t GetStartTime();///< return StartTime
  Double_t GetEndTime();///< return EndTime
  Double_t GetTimeOfRun();///< return TimeOfRun
  Double_t GetTempGrad(Double_t timeSec);///< return TempGrad
  
  Double_t EvalTempGrad(Double_t timeSec);///< evaluate temperature gradients for a certain time in seconds
  Double_t EvalAverageTemp(Double_t timeSec);///< evaluate average temperatures for a certain time in seconds
  Double_t EvalPress(Double_t timeSec);///< evaluate pressure for a certain time in seconds
  Double_t EvalVdrift(Double_t timeSec);///< evaluate  drift velocities for a certain time in seconds
  Double_t EvalVdriftcor(Double_t timeSec);///< evaluate drift velocities corrected for a certain time in seconds
  Double_t EvalGainF(Double_t timeSec);///< evaluate near gain for a certain time in seconds
  Double_t EvalGainN(Double_t timeSec);///< evaluate far gain for a certain time in seconds
  Double_t EvalCO2(Double_t timeSec);///< evaluate  CO2 content for a certain time in seconds
  Double_t EvalN2(Double_t timeSec);///< evaluate  N2 content for a certain time in seconds
  
  void FillAllGraphs();  ///< fill all the graphs after the tree
  void FillAllSplines();  ///< fill all the graphs after the splines
  
  void FillAverageTempGraph();///<graph of average temperatures 
  void FillTempGradGraph();///<graph of temperature gradients
  void FillPressGraph();///<graph of pressures
  void FillVdriftGraph();///<graph of drift velocities
  void FillVdriftcorGraph();///<graph of drift velocities corrected
  void FillGainFGraph();///<graph of near gain
  void FillGainNGraph();///<graph of far gain
  void FillCO2Graph();///<graph of CO2 content
  void FillN2Graph();///<graph of N2 content
  
  void FillAverageTempSpline();///< spline of average temperatures 
  void FillTempGradSpline();///<spline of  temperature gradients
  void FillPressSpline();///<spline of pressures
  void FillVdriftSpline();///<spline of drift velocities
  void FillVdriftcorSpline();///<spline of drift velocities corrected
  void FillGainFSpline();///<spline of near gain
  void FillGainNSpline();///<spline of far gain
  void FillCO2Spline();///<spline of CO2 content
  void FillN2Spline();///<spline of  N2 content
  TTree * GetTree(){return fGoofieValues;}
  void PrintTree(); ///< test: print tree values onto screen
  Bool_t    IsFolder() const {return kTRUE;}
protected:
  Long64_t fLinesInFile;// lines in ASCII file
  Double_t fStartTime;// StartTime
  Double_t fEndTime;// EndTime
  Double_t fTimeOfRun;//TimeOfRun
  Double_t fTempGrad;//TempGrad
  
  TGraph* fAverageTempGraph;//->graph of average temperatures 
  TGraph* fTempGradGraph;//->graph of temperature gradients
  TGraph* fPressGraph;//->graph of pressures
  TGraph* fVdriftGraph;//->graph of drift velocities
  TGraph* fVdriftcorGraph;//->graph of drift velocities corrected
  TGraph* fGainFGraph;//->graph of near gain
  TGraph* fGainNGraph;//->graph of far gain
  TGraph* fCO2Graph;//->graph of CO2 content
  TGraph* fN2Graph;//->graph of N2 content
  
  TSpline* fAverageTempSpline;// spline of average temperatures 
  TSpline* fTempGradSpline;//spline of  temperature gradients
  TSpline* fPressSpline;//spline of pressures
  TSpline* fVdriftSpline;//spline of drift velocities
  TSpline* fVdriftcorSpline;//spline of drift velocities corrected
  TSpline* fGainFSpline;//spline of near gain
  TSpline* fGainNSpline;//spline of far gain
  TSpline* fCO2Spline;//spline of CO2 content
  TSpline* fN2Spline;//spline of  N2 content
  
protected:
  TTree *fGoofieValues;   // tree with all Goofie values in branches

private:
  AliTPCGoofieValues(const AliTPCGoofieValues&); // Not implemented
  AliTPCGoofieValues& operator=(const AliTPCGoofieValues&); // Not implemented

  
  ClassDef(AliTPCGoofieValues,1)  //Basic ROOT object
};
