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

// AliFlowAnalysisWithLYZEventPlane:
// Class to do flow analysis with the event plane from the LYZ method
// author: N. van der Kolk (kolk@nikhef.nl)


class AliFlowAnalysisWithLYZEventPlane {
 public:
  AliFlowAnalysisWithLYZEventPlane();
  virtual ~AliFlowAnalysisWithLYZEventPlane();
  
  virtual void   Init();
  virtual void   Make(AliFlowEventSimple* fEvent, AliFlowLYZEventPlane* fLYZEP);
  //  virtual void   Make(AliFlowEventSimple* anEvent);
  virtual void   Finish();

  // input files
  void	   SetFirstRunFileName(TString name) 	
    { this->fFirstRunFileName = name ; }      // Sets input file name
  TString  GetFirstRunFileName() const		
    { return this->fFirstRunFileName ; }      // Gets output file name
  void     SetFirstRunFile(TFile* file)         
    { this->fFirstRunFile = file ; }          // Sets first run file

  void	   SetSecondRunFileName(TString name) 	
    { this->fSecondRunFileName = name ; }     // Sets input file name
  TString  GetSecondRunFileName() const		
    { return this->fSecondRunFileName ; }     // Gets output file name
  void     SetSecondRunFile(TFile* file)         
    { this->fSecondRunFile = file ; }         // Sets first run file

  // Output 
  void	   SetOutFileName(TString name)    { this->fOutFileName = name ; } 
  // Sets output file name
  TString  GetOutFileName() const	   { return this->fOutFileName ; } 
  // Gets output file name
  TFile*   GetOutFile() const              { return this->fOutFile ; }     
  // Gets output file


 private:

  //  AliFlowAnalysisWithLYZEventPlane(const AliFlowAnalysisWithLYZEventPlane& lyz);    // AliFlowAnalysisWithLYZEventPlane object cannot be copied
  //  void operator=(const AliFlowAnalysisWithLYZEventPlane &lyz);   
 
  TFile*             fOutFile;                //! 
  TFile*             fFirstRunFile ;          //! pointer to file from first run
  TFile*             fSecondRunFile ;         //! pointer to file from second run
  TString            fFirstRunFileName;       //!
  TString            fSecondRunFileName;      //!
  TString            fOutFileName;            //!

  //histograms
  //input
  TProfile*  fSecondReDtheta;                 //!
  TProfile*  fSecondImDtheta;                 //!
  TProfile*  fFirstr0theta;                   //!
  TProfile*  fSecondVPt;                      //!
  //output
  TProfile*  fHistProFlow;                    //!
  TProfile*  fHistProFlow2;                   //!
  TProfile*  fHistProWr;                      //!
  TProfile*  fHistProWrCorr;                  //!
  TH1D*      fHistFlow;                       //!
  TH1F*      fHistDeltaPhi;                   //!
  TH1F*      fHistDeltaPhi2;                  //!
  TH1F*      fHistDeltaPhihere;               //!
  TH1F*      fHistPhiEP;                      //!
  TH1F*      fHistPhiEPhere;                  //!
  TH1F*      fHistPhiLYZ;                     //!
  TH1F*      fHistPhiLYZ2;                    //!
  TProfile*  fHistProR0theta;                 //!
  TProfile*  fHistProReDtheta;                //!
  TProfile*  fHistProImDtheta;                //!
  
  AliFlowCommonHist* fCommonHists;            //!
  AliFlowCommonHistResults* fCommonHistsRes;  //!

  Int_t     fEventNumber;  // event counter
  //  Double_t  fQX;           // local flow vector
  //  Double_t  fQY;           // local flow vector
  AliFlowVector  fQ;       // flow vector
  TVector2  fQsum;         // flow vector sum
  Double_t  fQ2sum;        // flow vector sum squared
  Double_t  fQtheta;       // flow vector projected on ref. angle theta
   
  //  AliFlowEventSimple*  fEvent ;               //!
  AliFlowTrackSimple*  fTrack ;               //!
  AliFlowLYZEventPlane* fLYZEP ;              //!

  ClassDef(AliFlowAnalysisWithLYZEventPlane, 0);          // lyz analysis 
};

 #endif
