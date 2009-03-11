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
class TH1D;
class TProfile;
class TCollection;
class TList;
class TBrowser;

class AliFlowCommonHist: public TNamed {

 public:

  AliFlowCommonHist();
  AliFlowCommonHist(const char *name,const char *title = "AliFlowCommonHist");
  virtual ~AliFlowCommonHist();
  AliFlowCommonHist(const AliFlowCommonHist& aSetOfHists);

  Bool_t  IsFolder() const {return kTRUE;};
  //make fill methods here
  Bool_t FillControlHistograms(AliFlowEventSimple* anEvent);
  void Browse(TBrowser *b); 
  //make get methods here
  Double_t GetEntriesInPtBinRP(Int_t iBin);   //gets entries from fHistPtRP
  Double_t GetEntriesInPtBinPOI(Int_t iBin);  //gets entries from fHistPtPOI
  Double_t GetEntriesInEtaBinRP(Int_t iBin);  //gets entries from fHistEtaRP
  Double_t GetEntriesInEtaBinPOI(Int_t iBin); //gets entries from fHistEtaPOI
  Double_t GetMeanPt(Int_t iBin);             //gets the mean pt for this bin from fHistProMeanPtperBin   

  TH1F*     GetHistMultOrig()        {return fHistMultOrig;  } ;  
  TH1F*     GetHistMultRP()          {return fHistMultRP; } ;  
  TH1F*     GetHistMultPOI()         {return fHistMultPOI; } ;  
  TH1F*     GetHistPtRP()            {return fHistPtRP; } ;  
  TH1F*     GetHistPtPOI()           {return fHistPtPOI; } ;   
  TH1F*     GetHistPhiRP()           {return fHistPhiRP; } ;  
  TH1F*     GetHistPhiPOI()          {return fHistPhiPOI; } ;  
  TH1F*     GetHistEtaRP()           {return fHistEtaRP; } ;  
  TH1F*     GetHistEtaPOI()          {return fHistEtaPOI;  } ;   
  TProfile* GetHistProMeanPtperBin() {return fHistProMeanPtperBin; } ;
  TH1F*     GetHistQ()               {return fHistQ; } ;            
  TList*    GetHistList()            {return fHistList;} ;  

  virtual Double_t  Merge(TCollection *aList);  //merge function
  //method to print stats
  void Print(Option_t* option = "") const;

 
 private:

  AliFlowCommonHist& operator=(const AliFlowCommonHist& aSetOfHists);

  //define histograms here
  //control histograms
  TH1F*     fHistMultOrig;        //multiplicity before selection
  TH1F*     fHistMultRP;          //multiplicity for RP selection
  TH1F*     fHistMultPOI;         //multiplicity for POI selection
  TH1F*     fHistPtRP;            //pt distribution for RP selection
  TH1F*     fHistPtPOI;           //pt distribution for POI selection
  TH1F*     fHistPhiRP;           //phi distribution for RP selection
  TH1F*     fHistPhiPOI;          //phi distribution for POI selection
  TH1F*     fHistEtaRP;           //eta distribution for RP selection
  TH1F*     fHistEtaPOI;          //eta distribution for POI selection
  TProfile* fHistProMeanPtperBin; //mean pt for each pt bin (for POI selection)
  TH1F*     fHistQ;               //Qvector distribution

  TList*    fHistList;            //list to hold all histograms  

  ClassDef(AliFlowCommonHist,1)  // macro for rootcint
};
#endif

