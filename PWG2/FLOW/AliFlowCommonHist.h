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

class AliFlowCommonHist: public TObject {

 public:

  AliFlowCommonHist(TString input);
  virtual ~AliFlowCommonHist();

  //make fill methods here
  Bool_t FillControlHistograms(AliFlowEventSimple* Event);
 
  //make get methods here
  Double_t GetEntriesInPtBin(Int_t fBin);   //gets entries from fHistPtDiff
  Double_t GetMeanPt(Int_t fBin);           //gets the mean pt for this bin from fHistProMeanPtperBin   

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
   
  //  virtual Long64_t  Merge(TCollection *list);
 
 private:

  AliFlowCommonHist(const AliFlowCommonHist& aSetOfHists);
  AliFlowCommonHist& operator=(const AliFlowCommonHist& aSetOfHists);

  //define histograms here
  //control histograms
  TH1F*     fHistMultOrig;        
  TH1F*     fHistMultInt;        
  TH1F*     fHistMultDiff;       
  TH1F*     fHistPtInt;          
  TH1F*     fHistPtDiff;         
  TH1F*     fHistPhiInt;          
  TH1F*     fHistPhiDiff;         
  TH1F*     fHistEtaInt;          
  TH1F*     fHistEtaDiff;         
  TProfile* fHistProMeanPtperBin; 
  TH1F*     fHistQ;               
  
  ClassDef(AliFlowCommonHist,0);                 // macro for rootcint
};
#endif

