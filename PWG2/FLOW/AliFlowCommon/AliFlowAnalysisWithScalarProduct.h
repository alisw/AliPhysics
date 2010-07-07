/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHSCALARPRODUCT_H
#define ALIFLOWANALYSISWITHSCALARPRODUCT_H

#include "TString.h"

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TH1D;
class TH2D;
class TProfile;
class TList;
class TFile;
class Riostream;
class TDirectoryFile;

/////////////////////////////////////////////////////////////////////////////
// Description: Maker to analyze Flow from the Scalar Product method.
//               
// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl)              
/////////////////////////////////////////////////////////////////////////////
 
class AliFlowAnalysisWithScalarProduct {

 public:
 
   AliFlowAnalysisWithScalarProduct();            //default constructor
   virtual  ~AliFlowAnalysisWithScalarProduct();  //destructor
 
   void    Init();                                          //defines variables and histograms
   void    Make(AliFlowEventSimple* anEvent);               //calculates variables and fills histograms
   void    GetOutputHistograms(TList *outputListHistos);    //get pointers to all output histograms (called before Finish()) 
   void    Finish();                                        //saves histograms
   void    WriteHistograms(TString* outputFileName);        //writes histograms locally
   void    WriteHistograms(TString outputFileName);         //writes histograms locally
   void    WriteHistograms(TDirectoryFile *outputFileName); //writes histograms locally
    
   Double_t CalculateStatisticalError(Int_t bin, 
				      Double_t aStatErrorQaQb,
				      TProfile* aHistProUQ, 
				      TProfile* aHistProUQQaQb, 
				      TH1D** aHistSumOfWeights);

   void    SetDebug(Bool_t kt)   { this->fDebug = kt ; }
   Bool_t  GetDebug() const      { return this->fDebug ; }


   void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
   Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }

   //phi weights
   void    SetWeightsList(TList* const aWeightsList)  {this->fWeightsList = (TList*)aWeightsList->Clone();}
   TList*  GetWeightsList() const                     {return this->fWeightsList;}  
   void    SetUsePhiWeights(Bool_t const aPhiW)       {this->fUsePhiWeights = aPhiW;}
   Bool_t  GetUsePhiWeights() const                   {return this->fUsePhiWeights;}

   // Output 
   TList*    GetHistList() const    { return this->fHistList ; } // Gets output histogram list
   //histogram getters
   TProfile* GetHistProUQetaRP() const  {return this->fHistProUQetaRP;} 
   TProfile* GetHistProUQetaPOI() const {return this->fHistProUQetaPOI;}
   TProfile* GetHistProUQPtRP() const   {return this->fHistProUQPtRP;} 
   TProfile* GetHistProUQPtPOI() const  {return this->fHistProUQPtPOI;}
   TProfile* GetHistProQaQb() const     {return this->fHistProQaQb;}
   TProfile* GetHistProQaQbNorm() const {return this->fHistProQaQbNorm;}
   TH1D*     GetHistSumOfLinearWeights() const    {return this->fHistSumOfLinearWeights;}
   TH1D*     GetHistSumOfQuadraticWeights() const {return this->fHistSumOfQuadraticWeights;}

   TProfile* GetHistProUQQaQbPtRP() const   {return this->fHistProUQQaQbPtRP;}   
   TProfile* GetHistProUQQaQbEtaRP() const  {return this->fHistProUQQaQbEtaRP;}   
   TProfile* GetHistProUQQaQbPtPOI() const  {return this->fHistProUQQaQbPtPOI;}   
   TProfile* GetHistProUQQaQbEtaPOI() const {return this->fHistProUQQaQbEtaPOI;} 
   TH1D*     GetHistSumOfWeightsPtRP(Int_t i) const    {return this->fHistSumOfWeightsPtRP[i];}
   TH1D*     GetHistSumOfWeightsEtaRP(Int_t i) const   {return this->fHistSumOfWeightsEtaRP[i];}
   TH1D*     GetHistSumOfWeightsPtPOI(Int_t i) const   {return this->fHistSumOfWeightsPtPOI[i];}
   TH1D*     GetHistSumOfWeightsEtaPOI(Int_t i) const  {return this->fHistSumOfWeightsEtaPOI[i];}

   AliFlowCommonHist*        GetCommonHists() const    {return this->fCommonHists; }
   AliFlowCommonHistResults* GetCommonHistsRes() const {return this->fCommonHistsRes; }

   TH1D* GetHistQaQb() const       {return this->fHistQaQb;}
   TH1D* GetHistQaQbNorm() const   {return this->fHistQaQbNorm;}
   TH1D* GetHistQaQbCos() const    {return this->fHistQaQbCos;}
   TH1D* GetHistResolution() const {return this->fHistResolution;}
   TH1D* GetHistQaNorm() const     {return this->fHistQaNorm;}
   TH2D* GetHistQaNormvsMa() const {return this->fHistQaNormvsMa;}
   TH1D* GetHistQbNorm() const     {return this->fHistQbNorm;}
   TH2D* GetHistQbNormvsMb() const {return this->fHistQbNormvsMb;}
   TH2D* GetMavsMb() const         {return this->fHistMavsMb;}

   //histogram setters
   void SetHistProUQetaRP(TProfile* const aHistProUQetaRP)   
     {this->fHistProUQetaRP = aHistProUQetaRP;}
   void SetHistProUQetaPOI(TProfile* const aHistProUQetaPOI) 
     {this->fHistProUQetaPOI = aHistProUQetaPOI;}
   void SetHistProUQPtRP(TProfile* const aHistProUQPtRP)     
     {this->fHistProUQPtRP = aHistProUQPtRP;}
   void SetHistProUQPtPOI(TProfile* const aHistProUQPtPOI)   
     {this->fHistProUQPtPOI = aHistProUQPtPOI;}
   void SetHistProQaQb(TProfile* const aHistProQaQb)         
     {this->fHistProQaQb = aHistProQaQb;}
   void SetHistProQaQbNorm(TProfile* const aHistProQaQbNorm)         
     {this->fHistProQaQbNorm = aHistProQaQbNorm;}
   void SetHistSumOfLinearWeights(TH1D* const aHistSumOfLinearWeights) 
     {this->fHistSumOfLinearWeights = aHistSumOfLinearWeights;}
   void SetHistSumOfQuadraticWeights(TH1D* const aHistSumOfQuadraticWeights) 
     {this->fHistSumOfQuadraticWeights = aHistSumOfQuadraticWeights;}

   void SetHistProUQQaQbPtRP(TProfile* const aHistProUQQaQbPtRP)     
     {this->fHistProUQQaQbPtRP = aHistProUQQaQbPtRP;}   
   void SetHistProUQQaQbEtaRP(TProfile* const aHistProUQQaQbEtaRP)   
     {this->fHistProUQQaQbEtaRP = aHistProUQQaQbEtaRP;}   
   void SetHistProUQQaQbPtPOI(TProfile* const aHistProUQQaQbPtPOI)   
     {this->fHistProUQQaQbPtPOI = aHistProUQQaQbPtPOI;}   
   void SetHistProUQQaQbEtaPOI(TProfile* const aHistProUQQaQbEtaPOI) 
     {this->fHistProUQQaQbEtaPOI = aHistProUQQaQbEtaPOI;}
   void SetHistSumOfWeightsPtRP(TH1D* const aHistSumOfWeightsPtRP, Int_t const i) 
     {this->fHistSumOfWeightsPtRP[i] = aHistSumOfWeightsPtRP;}   
   void SetHistSumOfWeightsEtaRP(TH1D* const aHistSumOfWeightsEtaRP, Int_t const i) 
     {this->fHistSumOfWeightsEtaRP[i] = aHistSumOfWeightsEtaRP;}   
   void SetHistSumOfWeightsPtPOI(TH1D* const aHistSumOfWeightsPtPOI, Int_t const i) 
     {this->fHistSumOfWeightsPtPOI[i] = aHistSumOfWeightsPtPOI;}  
   void SetHistSumOfWeightsEtaPOI(TH1D* const aHistSumOfWeightsEtaPOI, Int_t const i) 
     {this->fHistSumOfWeightsEtaPOI[i] = aHistSumOfWeightsEtaPOI;}   

   void SetCommonHists(AliFlowCommonHist* const someCommonHists)              
     {this->fCommonHists = someCommonHists; }
   void SetCommonHistsRes(AliFlowCommonHistResults* const someCommonHistsRes) 
     {this->fCommonHistsRes = someCommonHistsRes; }

   void SetHistQaQb(TH1D* const aHistQaQb)
     {this->fHistQaQb = aHistQaQb;}
   void SetHistQaQbNorm(TH1D* const aHistQaQbNorm)
     {this->fHistQaQbNorm = aHistQaQbNorm;}
   void SetHistQaQbCos(TH1D* const aHistQaQbCos)
     {this->fHistQaQbCos = aHistQaQbCos;}
   void SetHistResolution(TH1D* const aHistResolution)
     {this->fHistResolution = aHistResolution;}
   void SetHistQaNorm(TH1D* const aHistQaNorm)
     {this->fHistQaNorm = aHistQaNorm;}
   void SetHistQaNormvsMa(TH2D* const aHistQaNormvsMa)
     {this->fHistQaNormvsMa = aHistQaNormvsMa;}
   void SetHistQbNorm(TH1D* const aHistQbNorm)
     {this->fHistQbNorm = aHistQbNorm;}
   void SetHistQbNormvsMb(TH2D* const aHistQbNormvsMb)
     {this->fHistQbNormvsMb = aHistQbNormvsMb;}
   void SetHistMavsMb(TH2D* const aHistMavsMb)
     {this->fHistMavsMb = aHistMavsMb;}
   
 private:
   AliFlowAnalysisWithScalarProduct(const AliFlowAnalysisWithScalarProduct& anAnalysis);            //copy constructor
   AliFlowAnalysisWithScalarProduct& operator=(const AliFlowAnalysisWithScalarProduct& anAnalysis); //assignment operator 
      
   Int_t      fEventNumber;      // event counter
   Bool_t     fDebug ;           // flag for analysis: more print statements

   Double_t   fRelDiffMsub;      // the relative difference the two subevent multiplicities can have

   TList*     fWeightsList;      // list holding input histograms with phi weights
   Bool_t     fUsePhiWeights;    // use phi weights
   TH1F*      fPhiWeightsSub0;   // histogram holding phi weights for subevent 0
   TH1F*      fPhiWeightsSub1;   // histogram holding phi weights for subevent 1

   TList*     fHistList;         //list to hold all output histograms  

   TProfile*  fHistProUQetaRP;   //uQ(eta) for RP
   TProfile*  fHistProUQetaPOI;  //uQ(eta) for POI
   TProfile*  fHistProUQPtRP;    //uQ(pt) for RP
   TProfile*  fHistProUQPtPOI;   //uQ(pt) for POI
   TProfile*  fHistProQaQb;      //average of QaQb 
   TProfile*  fHistProQaQbNorm;  //average of QaQb/MaMb
   TH1D*      fHistSumOfLinearWeights;     //holds sum of Ma*Mb
   TH1D*      fHistSumOfQuadraticWeights;  //holds sum of (Ma*Mb)^2
   
   TProfile*  fHistProUQQaQbPtRP;         //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbEtaRP;        //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbPtPOI;        //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbEtaPOI;       //holds weighted average of <QuQaQb>
   TH1D*      fHistSumOfWeightsPtRP[3];   //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsEtaRP[3];  //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsPtPOI[3];  //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsEtaPOI[3]; //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   
   AliFlowCommonHist*        fCommonHists;    //control histograms
   AliFlowCommonHistResults* fCommonHistsRes; //results histograms
   
   TH1D*      fHistQaQb;         // distribution of QaQb
   TH1D*      fHistQaQbNorm;     // distribution of QaQb/MaMb
   TH1D*      fHistQaQbCos;      // distribution of the angle between Qa and Qb (from Acos (va*vb))
   TH1D*      fHistResolution;   // distribution of cos(2(phi_a - phi_b))
   TH1D*      fHistQaNorm;       // distribution of Qa/Ma
   TH2D*      fHistQaNormvsMa;   // distribution of Qa/Ma vs Ma
   TH1D*      fHistQbNorm;       // distribution of Qb/Mb
   TH2D*      fHistQbNormvsMb;   // distribution of Qb/Mb vs Mb
   TH2D*      fHistMavsMb;       // Ma vs Mb
      
   ClassDef(AliFlowAnalysisWithScalarProduct,0)  // macro for rootcint
     };
 

#endif
