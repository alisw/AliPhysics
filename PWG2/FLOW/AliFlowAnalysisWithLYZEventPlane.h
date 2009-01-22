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
  void           WriteHistograms(TString* outputFileName);

  void      SetEventNumber(Int_t n)      { this->fEventNumber = n; }
  Int_t     GetEventNumber() const       { return this->fEventNumber; }
  void      SetQ2sum(Double_t d)         { this->fQ2sum = d; }
  Double_t  GetQ2sum()                   { return this->fQ2sum; }

  //output
  TList*             GetHistList() const     {return this->fHistList; }
  AliFlowCommonHist* GetCommonHists() const  { return this->fCommonHists; }
  void               SetCommonHists(AliFlowCommonHist* aCommonHist)  
     { this->fCommonHists = aCommonHist; }
  AliFlowCommonHistResults* GetCommonHistsRes() const  
     { return this->fCommonHistsRes; }
  void               SetCommonHistsRes(AliFlowCommonHistResults* aCommonHistResult) 
     { this->fCommonHistsRes = aCommonHistResult; }

  // !!!!! make getters and setters for all histograms
  TProfile*  GetSecondReDtheta() {return this->fSecondReDtheta; } 
  void       SetSecondReDtheta(TProfile* aSecondReDtheta) 
    {this->fSecondReDtheta = aSecondReDtheta; }
  TProfile*  GetSecondImDtheta() {return this->fSecondImDtheta; }
  void       SetSecondImDtheta(TProfile* aSecondImDtheta)
    {this->fSecondImDtheta = aSecondImDtheta; }
  TProfile*  GetFirstr0theta()   {return this->fFirstr0theta; }
  void       SetFirstr0theta(TProfile* aFirstr0theta)
    {this->fFirstr0theta = aFirstr0theta; }
  TProfile*  GetHistProFlow()    {return this->fHistProFlow;}
  void       SetHistProFlow(TProfile* aHistProFlow)
    {this->fHistProFlow =aHistProFlow; }        
  TProfile*  GetHistProFlow2()   {return this->fHistProFlow2;} 
  void       SetHistProFlow2(TProfile* aHistProFlow2)
    {this->fHistProFlow2 = aHistProFlow2; }      
  TProfile*  GetHistProWr()      {return this->fHistProWr; }
  void       SetHistProWr(TProfile* aHistProWr)
    {this->fHistProWr = aHistProWr; }
  TProfile*  GetHistProWrCorr()  {return this->fHistProWrCorr; }
  void       SetHistProWrCorr(TProfile* aHistProWrCorr)
    {this->fHistProWrCorr = aHistProWrCorr; }
  TH1F*      GetHistQsumforChi() {return this->fHistQsumforChi; }
  void       SetHistQsumforChi(TH1F* aHistQsumforChi) 
    {this->fHistQsumforChi = aHistQsumforChi; }
  TH1F*      GetHistDeltaPhi()   {return this->fHistDeltaPhi; }  
  void       SetHistDeltaPhi(TH1F* aHistDeltaPhi)
    {this->fHistDeltaPhi = aHistDeltaPhi; }
  TH1F*      GetHistDeltaPhi2()  {return this->fHistDeltaPhi2; } 
  void       SetHistDeltaPhi2(TH1F* aHistDeltaPhi2)
    {this->fHistDeltaPhi2 = aHistDeltaPhi2; }
  TH1F*      GetHistDeltaPhihere() {return this->fHistDeltaPhihere; }
  void       SetHistDeltaPhihere(TH1F* aHistDeltaPhihere)
    {this->fHistDeltaPhihere = aHistDeltaPhihere; }
  TH1F*      GetHistPhiEP()      {return this->fHistPhiEP; }   
  void       SetHistPhiEP(TH1F* aHistPhiEP)
    {this->fHistPhiEP = aHistPhiEP; }  
  TH1F*      GetHistPhiEPhere()  {return this->fHistPhiEPhere; }
  void       SetHistPhiEPhere(TH1F* aHistPhiEPhere)
    {this->fHistPhiEPhere = aHistPhiEPhere; }       
  TH1F*      GetHistPhiLYZ()     {return this->fHistPhiLYZ; }  
  void       SetHistPhiLYZ(TH1F* aHistPhiLYZ)
    {this->fHistPhiLYZ = aHistPhiLYZ; }       
  TH1F*      GetHistPhiLYZ2()    {return this->fHistPhiLYZ2;}               
  void       SetHistPhiLYZ2(TH1F* aHistPhiLYZ2)
    {this->fHistPhiLYZ2 = aHistPhiLYZ2; }

  //input
  void       SetSecondRunList(TList* list) { this->fSecondRunList = list; }
  TList*     GetSecondRunList()            { return this->fSecondRunList; }

 private:

  AliFlowAnalysisWithLYZEventPlane(const AliFlowAnalysisWithLYZEventPlane& aAnalysis);             // copy constructor
  AliFlowAnalysisWithLYZEventPlane& operator=(const AliFlowAnalysisWithLYZEventPlane& aAnalysis);  // assignment operator

  //histograms
  TList*     fHistList;                       //list ro hold all histograms
  TList*     fSecondRunList;                  //list from Second LYZ run output
  //input
  TProfile*  fSecondReDtheta;                 // input profile
  TProfile*  fSecondImDtheta;                 // input profile
  TProfile*  fFirstr0theta;                   // input profile
  //output
  TProfile*  fHistProFlow;                    //
  TProfile*  fHistProFlow2;                   //
  TProfile*  fHistProWr;                      //
  TProfile*  fHistProWrCorr;                  //
  TH1F*      fHistQsumforChi;                 //
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
     

  ClassDef(AliFlowAnalysisWithLYZEventPlane, 1);          // lyz analysis 
};

 #endif

