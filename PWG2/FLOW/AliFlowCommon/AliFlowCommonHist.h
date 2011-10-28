/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:$ */

#ifndef AliFlowCommonHist_H
#define AliFlowCommonHist_H

// AliFlowCommonHist:
// Description: Class to organise common histograms for Flow Analysis
// authors: N.K A.B R.S

             
class AliFlowEventSimple;
class AliFlowTrackSimple;
class TH1F;
class TH2F;
class TH1D;
class TProfile;
class TCollection;
class TList;
class TBrowser;

class AliFlowCommonHist: public TNamed {

 public:

  AliFlowCommonHist();
  AliFlowCommonHist(const char *name,const char *title = "AliFlowCommonHist",Bool_t bookOnlyBasic = kFALSE);
  virtual ~AliFlowCommonHist();
  AliFlowCommonHist(const AliFlowCommonHist& aSetOfHists);

  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b);
  //merge function
  virtual Double_t  Merge(TCollection *aList);  
  //method to print stats
  void Print(Option_t* option = "") const;
 
  //fill method
  Bool_t FillControlHistograms(AliFlowEventSimple* anEvent,TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);
  
  //getters
  Double_t GetEntriesInPtBinRP(Int_t iBin);   //gets entries from fHistPtRP
  Double_t GetEntriesInPtBinPOI(Int_t iBin);  //gets entries from fHistPtPOI
  Double_t GetEntriesInEtaBinRP(Int_t iBin);  //gets entries from fHistEtaRP
  Double_t GetEntriesInEtaBinPOI(Int_t iBin); //gets entries from fHistEtaPOI
  Double_t GetMeanPt(Int_t iBin);             //gets the mean pt for this bin from fHistProMeanPtperBin   

  TH1F*     GetHistMultRP()          {return fHistMultRP; } ;  
  TH1F*     GetHistMultPOI()         {return fHistMultPOI; } ; 
  TH2F*     GetHistMultPOIvsRP()     {return fHistMultPOIvsRP; } ;
  TH1F*     GetHistPtRP()            {return fHistPtRP; } ;  
  TH1F*     GetHistPtPOI()           {return fHistPtPOI; } ;
  TH1F*     GetHistPtSub0()          {return fHistPtSub0; } ;
  TH1F*     GetHistPtSub1()          {return fHistPtSub1; } ;
  TH1F*     GetHistPhiRP()           {return fHistPhiRP; } ;  
  TH1F*     GetHistPhiPOI()          {return fHistPhiPOI; } ;  
  TH1F*     GetHistPhiSub0()         {return fHistPhiSub0; } ; 
  TH1F*     GetHistPhiSub1()         {return fHistPhiSub1; } ; 
  TH1F*     GetHistEtaRP()           {return fHistEtaRP; } ;  
  TH1F*     GetHistEtaPOI()          {return fHistEtaPOI;  } ;  
  TH1F*     GetHistEtaSub0()         {return fHistEtaSub0;  } ; 
  TH1F*     GetHistEtaSub1()         {return fHistEtaSub1;  } ; 
  TH2F*     GetHistPhiEtaRP()        {return fHistPhiEtaRP;  } ; 
  TH2F*     GetHistPhiEtaPOI()       {return fHistPhiEtaPOI;  } ; 
  TProfile* GetHistProMeanPtperBin() {return fHistProMeanPtperBin; } ;
  TH2F*     GetHistWeightvsPhi()     {return fHistWeightvsPhi; } ;
  TH1F*     GetHistQ()               {return fHistQ; } ;  
  TH1F*     GetHistAngleQ()          {return fHistAngleQ; }
  TH1F*     GetHistAngleQSub0()      {return fHistAngleQSub0; }
  TH1F*     GetHistAngleQSub1()      {return fHistAngleQSub1; }
  TProfile* GetHarmonic()            {return fHarmonic; } ; 
  TProfile* GetRefMultVsNoOfRPs()    {return fRefMultVsNoOfRPs; } ;
  TH1F*     GetHistRefMult()         {return fHistRefMult; } ; 
  TList*    GetHistList()            {return fHistList;} ;  

   
 private:

  AliFlowCommonHist& operator=(const AliFlowCommonHist& aSetOfHists);

  //define histograms here
  //control histograms
  Bool_t    fBookOnlyBasic;       // book and fill only control histos needed for all methods
  TH1F*     fHistMultRP;          // multiplicity for RP selection
  TH1F*     fHistMultPOI;         // multiplicity for POI selection
  TH2F*     fHistMultPOIvsRP;     // multiplicity for POI versus RP
  TH1F*     fHistPtRP;            // pt distribution for RP selection
  TH1F*     fHistPtPOI;           // pt distribution for POI selection
  TH1F*     fHistPtSub0;          // pt distribution for subevent 0
  TH1F*     fHistPtSub1;          // pt distribution for subevent 1
  TH1F*     fHistPhiRP;           // phi distribution for RP selection
  TH1F*     fHistPhiPOI;          // phi distribution for POI selection
  TH1F*     fHistPhiSub0;         // phi distribution for subevent 0
  TH1F*     fHistPhiSub1;         // phi distribution for subevent 1
  TH1F*     fHistEtaRP;           // eta distribution for RP selection
  TH1F*     fHistEtaPOI;          // eta distribution for POI selection
  TH1F*     fHistEtaSub0;         // eta distribution for subevent 0
  TH1F*     fHistEtaSub1;         // eta distribution for subevent 1
  TH2F*     fHistPhiEtaRP;        // eta vs phi for RP selection
  TH2F*     fHistPhiEtaPOI;       // eta vs phi for POI selection
  TProfile* fHistProMeanPtperBin; // mean pt for each pt bin (for POI selection)
  TH2F*     fHistWeightvsPhi;     // particle weight vs particle phi
  TH1F*     fHistQ;               // Qvector distribution
  TH1F*     fHistAngleQ;          // distribution of angle of Q vector
  TH1F*     fHistAngleQSub0;      // distribution of angle of subevent 0 Q vector
  TH1F*     fHistAngleQSub1;      // distribution of angle of subevent 1 Q vector
  TProfile* fHarmonic;            // harmonic 
  TProfile* fRefMultVsNoOfRPs;    // <reference multiplicity> versus # of RPs
  TH1F*     fHistRefMult;         // reference multiplicity distribution
  
  TList*    fHistList;            // list to hold all histograms  

  ClassDef(AliFlowCommonHist,3)   // macro for rootcint
};
#endif

