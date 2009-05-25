/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWANALYSISWITHLYZEVENTPLANE_H
#define ALIFLOWANALYSISWITHLYZEVENTPLANE_H

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
  virtual void   GetOutputHistograms(TList *outputListHistos); //get pointers to all output histograms (called before Finish()) 
  virtual void   Finish();
  void           WriteHistograms(TString* outputFileName);
  void           WriteHistograms(TString outputFileName);

  void      SetEventNumber(Int_t n)      { this->fEventNumber = n; }
  Int_t     GetEventNumber() const       { return this->fEventNumber; }
  void      SetQ2sum(Double_t d)         { this->fQ2sum = d; }
  Double_t  GetQ2sum() const             { return this->fQ2sum; }

  //output
  TList*             GetHistList() const     {return this->fHistList; }
  AliFlowCommonHist* GetCommonHists() const  { return this->fCommonHists; }
  void               SetCommonHists(AliFlowCommonHist* const aCommonHist)  
     { this->fCommonHists = aCommonHist; }
  AliFlowCommonHistResults* GetCommonHistsRes() const  
     { return this->fCommonHistsRes; }
  void               SetCommonHistsRes(AliFlowCommonHistResults* const aCommonHistResult) 
     { this->fCommonHistsRes = aCommonHistResult; }

  // !!!!! make getters and setters for all histograms
  TProfile*  GetSecondReDtheta() const {return this->fSecondReDtheta; } 
  void       SetSecondReDtheta(TProfile* const aSecondReDtheta) 
    {this->fSecondReDtheta = aSecondReDtheta; }
  TProfile*  GetSecondImDtheta() const {return this->fSecondImDtheta; }
  void       SetSecondImDtheta(TProfile* const aSecondImDtheta)
    {this->fSecondImDtheta = aSecondImDtheta; }
  TProfile*  GetFirstr0theta() const   {return this->fFirstr0theta; }
  void       SetFirstr0theta(TProfile* const aFirstr0theta)
    {this->fFirstr0theta = aFirstr0theta; }
  
  TProfile*  GetHistProVetaRP() const   {return this->fHistProVetaRP;}
  void       SetHistProVetaRP(TProfile* const aHistProVetaRP)
    {this->fHistProVetaRP =aHistProVetaRP; }
  TProfile*  GetHistProVetaPOI() const  {return this->fHistProVetaPOI;} 
  void       SetHistProVetaPOI(TProfile* const aHistProVetaPOI)
    {this->fHistProVetaPOI = aHistProVetaPOI; } 
  TProfile*  GetHistProVPtRP() const    {return this->fHistProVPtRP;}
  void       SetHistProVPtRP(TProfile* const aHistProVPtRP)
    {this->fHistProVPtRP =aHistProVPtRP; } 
  TProfile*  GetHistProVPtPOI() const   {return this->fHistProVPtPOI;} 
  void       SetHistProVPtPOI(TProfile* const aHistProVPtPOI)
    {this->fHistProVPtPOI = aHistProVPtPOI; }      
  TProfile*  GetHistProWr() const       {return this->fHistProWr; }
  void       SetHistProWr(TProfile* const aHistProWr)
    {this->fHistProWr = aHistProWr; }
  TProfile*  GetHistProWrCorr() const {return this->fHistProWrCorr; }
  void       SetHistProWrCorr(TProfile* const aHistProWrCorr)
    {this->fHistProWrCorr = aHistProWrCorr; }
  TH1F*      GetHistQsumforChi() const {return this->fHistQsumforChi; }
  void       SetHistQsumforChi(TH1F* const aHistQsumforChi) 
    {this->fHistQsumforChi = aHistQsumforChi; }
  TH1F*      GetHistDeltaPhi()  const {return this->fHistDeltaPhi; }  
  void       SetHistDeltaPhi(TH1F* const aHistDeltaPhi)
    {this->fHistDeltaPhi = aHistDeltaPhi; }
  TH1F*      GetHistDeltaPhi2() const {return this->fHistDeltaPhi2; } 
  void       SetHistDeltaPhi2(TH1F* const aHistDeltaPhi2)
    {this->fHistDeltaPhi2 = aHistDeltaPhi2; }
  TH1F*      GetHistDeltaPhihere() const {return this->fHistDeltaPhihere; }
  void       SetHistDeltaPhihere(TH1F* const aHistDeltaPhihere)
    {this->fHistDeltaPhihere = aHistDeltaPhihere; }
  TH1F*      GetHistPhiEP() const     {return this->fHistPhiEP; }   
  void       SetHistPhiEP(TH1F* const aHistPhiEP)
    {this->fHistPhiEP = aHistPhiEP; }  
  TH1F*      GetHistPhiEPhere() const {return this->fHistPhiEPhere; }
  void       SetHistPhiEPhere(TH1F* const aHistPhiEPhere)
    {this->fHistPhiEPhere = aHistPhiEPhere; }       
  TH1F*      GetHistPhiLYZ()  const   {return this->fHistPhiLYZ; }  
  void       SetHistPhiLYZ(TH1F* const aHistPhiLYZ)
    {this->fHistPhiLYZ = aHistPhiLYZ; }       
  TH1F*      GetHistPhiLYZ2() const   {return this->fHistPhiLYZ2;}               
  void       SetHistPhiLYZ2(TH1F* const aHistPhiLYZ2)
    {this->fHistPhiLYZ2 = aHistPhiLYZ2; }

  //input
  void       SetSecondRunList(TList* const list) { this->fSecondRunList = list; }
  TList*     GetSecondRunList()  const           { return this->fSecondRunList; }

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
  TProfile*  fHistProVetaRP;                  //
  TProfile*  fHistProVetaPOI;                 //
  TProfile*  fHistProVPtRP;                   //
  TProfile*  fHistProVPtPOI;                  //
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

