/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/*************************************
 *   AliFlowCommonHistResults:       *
 *   class to organize the common    *
 *   histograms for Flow Analysis    * 
 *                                   * 
 * authors: Naomi van der Kolk       *
 *           (kolk@nikhef.nl)        *  
 *          Raimond Snellings        *
 *           (snelling@nikhef.nl)    * 
 *          Ante Bilandzic           *
 *           (anteb@nikhef.nl)       * 
 * **********************************/
 
#ifndef ALIFLOWCOMMONHISTRESULTS_H
#define ALIFLOWCOMMONHISTRESULTS_H

class TH1F;
class TH1D;
class TCollection;
class TList;
class TBrowser;

class AliFlowCommonHistResults : public TNamed {

 public:
  AliFlowCommonHistResults();                                                           //default constructor
  AliFlowCommonHistResults(const char *name, const char *title = "AliFlowCommonHist");  //constructor
  virtual ~AliFlowCommonHistResults();                                                  //destructor

  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 

  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
  //                              !!!     to be removed    !!!
  //make fill methods here
  Bool_t FillIntegratedFlow(Double_t aV, Double_t anError);                //fill fHistIntFlow
  Bool_t FillDifferentialFlow(Int_t aBin, Double_t av, Double_t anError);  //fill fHistDiffFlow
  Bool_t FillChi(Double_t aChi);                                           //fill fHistChi
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 

  //make fill methods here
  Bool_t FillIntegratedFlowRP(Double_t aV, Double_t anError);                   //fill fHistIntFlowRP
  Bool_t FillChiRP(Double_t aChi);                                              //fill fHistChiRP
  Bool_t FillDifferentialFlowPtRP(Int_t aBin, Double_t av, Double_t anError);   //fill fHistDiffFlowPtRP
  Bool_t FillDifferentialFlowEtaRP(Int_t aBin, Double_t av, Double_t anError);  //fill fHistDiffFlowEtaRP 
  Bool_t FillIntegratedFlowPOI(Double_t aV, Double_t anError);                  //fill fHistIntFlowPOI
  Bool_t FillChiPOI(Double_t aChi);                                             //fill fHistChiPOI
  Bool_t FillDifferentialFlowPtPOI(Int_t aBin, Double_t av, Double_t anError);  //fill fHistDiffFlowPtPOI
  Bool_t FillDifferentialFlowEtaPOI(Int_t aBin, Double_t av, Double_t anError); //fill fHistDiffFlowEtaPOI   
             
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
  //                              !!!     to be removed    !!!          
  //make get methods here
  TH1D* GetHistDiffFlow()               {return fHistDiffFlow; } ; 
  TH1D* GetHistChi()                    {return fHistChi; } ;
  TH1D* GetHistIntFlow()                {return fHistIntFlow; } ;
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
    
  //RP = Reaction Plane (RP denotes particles used to determine the reaction plane)  
  TH1D *GetHistIntFlowRP() {return fHistIntFlowRP; } ; 
  TH1D *GetHistChiRP() {return fHistChiRP; } ;            
  TH1D *GetHistDiffFlowPtRP() {return fHistDiffFlowPtRP; } ; 
  TH1D *GetHistDiffFlowEtaRP() {return fHistDiffFlowEtaRP; } ; 
  //POI = Particles Of Interest (POI denotes particles of interest for the final physical results for int. and diff. flow)
  TH1D *GetHistIntFlowPOI() {return fHistIntFlowPOI; } ;
  TH1D *GetHistChiPOI() {return fHistChiPOI; } ; 
  TH1D *GetHistDiffFlowPtPOI() {return fHistDiffFlowPtPOI; } ; 
  TH1D *GetHistDiffFlowEtaPOI() {return fHistDiffFlowEtaPOI; } ; 
  
  TList*   GetHistList()                   {return fHistList;} ;  

  virtual Double_t  Merge(TCollection *aList);  //merge function
  void Print(Option_t* option = "") const;      //method to print stats

 private:

  AliFlowCommonHistResults(const AliFlowCommonHistResults& aSetOfResultHists);            //copy constructor
  AliFlowCommonHistResults& operator=(const AliFlowCommonHistResults& aSetOfResultHists); //assignment operator
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
  //                              !!!     to be removed    !!!
  TH1D*     fHistIntFlow;      //integrated flow
  TH1D*     fHistDiffFlow;     //differential flow
  TH1D*     fHistChi;          //resolution
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  //RP = Reaction Plane (RP denotes particles used to determine the reaction plane)  
  TH1D *fHistIntFlowRP;      //integrated flow
  TH1D *fHistChiRP;          //chi
  TH1D *fHistDiffFlowPtRP;   //differential flow (Pt)
  TH1D *fHistDiffFlowEtaRP;  //differential flow (Eta)
  //POI = Particles Of Interest (POI denotes particles of interest for the final physical results of int. and diff. flow)
  TH1D *fHistIntFlowPOI;     //integrated flow
  TH1D *fHistChiPOI;         //chi  
  TH1D *fHistDiffFlowPtPOI;  //differential flow (Pt)
  TH1D *fHistDiffFlowEtaPOI; //differential flow (Eta)     
  
  TList*    fHistList;         //list to hold all histograms

  ClassDef(AliFlowCommonHistResults,1)                 // macro for rootcint
};
 
#endif

