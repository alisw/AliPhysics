/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowAnalysisWithLYZEventPlane_H
#define AliFlowAnalysisWithLYZEventPlane_H

class AliFlowVector;
class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowLYZEventPlane;

class TString;
class TFile;
class TProfile;
class TH1F;
class TH1D;
class TList;
class Riostream;

// AliFlowAnalysisWithLYZEventPlane:
// Class to do flow analysis with the event plane from the LYZ method
// author: N. van der Kolk (kolk@nikhef.nl)


class AliFlowAnalysisWithLYZEventPlane {

 public:

  AliFlowAnalysisWithLYZEventPlane();                 //default constructor
  virtual ~AliFlowAnalysisWithLYZEventPlane();        //destructor
  
  virtual void   Init();
  virtual void   Make(AliFlowEventSimple* fEvent, AliFlowLYZEventPlane* fLYZEP);
  virtual void   Finish();

  // input files
  void	   SetSecondRunFileName(TString name)  { this->fSecondRunFileName = name ; }     // Sets input file name
  TString  GetSecondRunFileName() const	       { return this->fSecondRunFileName ; }     // Gets output file name
  void     SetSecondRunFile(TFile* file)       { this->fSecondRunFile = file ; }         // Sets first run file

  //output
  TList* GetHistList() const                   {return this->fHistList; }


 private:

  AliFlowAnalysisWithLYZEventPlane(const AliFlowAnalysisWithLYZEventPlane& aAnalysis);             // copy constructor
  AliFlowAnalysisWithLYZEventPlane& operator=(const AliFlowAnalysisWithLYZEventPlane& aAnalysis);  // assignment operator

  TFile*             fSecondRunFile ;         // pointer to file from second run
  TString            fSecondRunFileName;    

  //histograms
  TList*     fHistList;                       //list ro hold all histograms
  //input
  TProfile*  fSecondReDtheta;                 // input profile
  TProfile*  fSecondImDtheta;                 // input profile
  TProfile*  fFirstr0theta;                   // input profile
  TProfile*  fSecondVPt;                      // input profile
  //output
  TProfile*  fHistProFlow;                    //
  TProfile*  fHistProFlow2;                   //
  TProfile*  fHistProWr;                      //
  TProfile*  fHistProWrCorr;                  //
  TH1D*      fHistFlow;                       //
  TH1F*      fHistDeltaPhi;                   //
  TH1F*      fHistDeltaPhi2;                  //
  TH1F*      fHistDeltaPhihere;               //
  TH1F*      fHistPhiEP;                      //
  TH1F*      fHistPhiEPhere;                  //
  TH1F*      fHistPhiLYZ;                     //
  TH1F*      fHistPhiLYZ2;                    //
  
  AliFlowCommonHist* fCommonHists;            //
  AliFlowCommonHistResults* fCommonHistsRes;  //

  Int_t     fEventNumber;                     // event counter

  TVector2  *fQsum;                           // flow vector sum
  Double_t  fQ2sum;                           // flow vector sum squared
     

  ClassDef(AliFlowAnalysisWithLYZEventPlane, 0);          // lyz analysis 
};

 #endif

