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

class AliFlowCommonHist: public TObject {

 public:

  AliFlowCommonHist();
  AliFlowCommonHist(TString anInput);
  virtual ~AliFlowCommonHist();

  //make fill methods here
  Bool_t FillControlHistograms(AliFlowEventSimple* anEvent);
 
  //make get methods here
  Double_t GetEntriesInPtBin(Int_t iBin);   //gets entries from fHistPtDiff
  Double_t GetMeanPt(Int_t iBin);           //gets the mean pt for this bin from fHistProMeanPtperBin   

  TH1F*     GetfHistMultOrig()               {return fHistMultOrig;  } ;  
  TH1F*     GetfHistMultInt()                {return fHistMultInt; } ;  
  TH1F*     GetfHistMultDiff()               {return fHistMultDiff; } ;  
  TH1F*     GetfHistPtInt()                  {return fHistPtInt; } ;  
  TH1F*     GetfHistPtDiff()                 {return fHistPtDiff; } ;   
  TH1F*     GetfHistPhiInt()                 {return fHistPhiInt; } ;  
  TH1F*     GetfHistPhiDiff()                {return fHistPhiDiff; } ;  
  TH1F*     GetfHistEtaInt()                 {return fHistEtaInt; } ;  
  TH1F*     GetfHistEtaDiff()                {return fHistEtaDiff;  } ;   
  TProfile* GetfHistProMeanPtperBin()        {return fHistProMeanPtperBin; } ;
  TH1F*     GetfHistQ()                      {return fHistQ; } ;            
  TList*    GetHistList()                    {return fHistList;} ;  

  virtual Double_t  Merge(TCollection *aList);  //merge function
 
 private:

  AliFlowCommonHist(const AliFlowCommonHist& aSetOfHists);
  AliFlowCommonHist& operator=(const AliFlowCommonHist& aSetOfHists);

  //define histograms here
  //control histograms
  TH1F*     fHistMultOrig;        //multiplicity before selection
  TH1F*     fHistMultInt;         //multiplicity for integrated flow
  TH1F*     fHistMultDiff;        //multiplicity for differential flow
  TH1F*     fHistPtInt;           //pt distribution for integrated flow
  TH1F*     fHistPtDiff;          //pt distribution for differential flow
  TH1F*     fHistPhiInt;          //phi distribution for integrated flow
  TH1F*     fHistPhiDiff;         //phi distribution for differential flow
  TH1F*     fHistEtaInt;          //eta distribution for integrated flow
  TH1F*     fHistEtaDiff;         //eta distribution for differential flow
  TProfile* fHistProMeanPtperBin; //mean pt for each pt bin (for differential flow)
  TH1F*     fHistQ;               //Qvector distribution

  TList* fHistList;               //list to hold all histograms  

  ClassDef(AliFlowCommonHist,0);  // macro for rootcint
};
#endif

