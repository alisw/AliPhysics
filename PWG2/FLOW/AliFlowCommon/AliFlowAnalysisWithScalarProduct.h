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
   void    Make(AliFlowEventSimple* anEvent);               //calls FillSP and FillmuQ
   void    FillSP(AliFlowEventSimple* anEvent);             //calculates variables and fills histograms for v2
   void    FillmuQ(AliFlowEventSimple* anEvent);            //calculates variables and fills histograms for uQ
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

   virtual void StoreFlags(); //store all booleans needed in Finish()
   virtual void AccessFlags(); //access all booleans needed in Finish()

   void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
   Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }

   //phi weights
   void    SetWeightsList(TList* const aWeightsList)  {this->fWeightsList = (TList*)aWeightsList->Clone();}
   TList*  GetWeightsList() const                     {return this->fWeightsList;}  
   void    SetUsePhiWeights(Bool_t const aPhiW)       {this->fUsePhiWeights = aPhiW;}
   Bool_t  GetUsePhiWeights() const                   {return this->fUsePhiWeights;}
   
   // correction for non-uniform acceptance:
   void   SetApplyCorrectionForNUA(Bool_t const acfNUA) {this->fApplyCorrectionForNUA = acfNUA;}
   Bool_t GetApplyCorrectionForNUA() const              {return this->fApplyCorrectionForNUA;}   
   // harmonic:     
   void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
   Int_t GetHarmonic() const {return this->fHarmonic;};
   // Total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb":  
   void SetTotalQvector(const char *tqv) {*this->fTotalQvector = tqv;};     
   
   // Output 
   TList*    GetHistList() const    { return this->fHistList ; } // Gets output histogram list
   //histogram getters
   TProfile* GetHistProFlags() const    {return this->fHistProFlags;};   
   TProfile* GetHistProUQetaRP() const  {return this->fHistProUQetaRP;} 
   TProfile* GetHistProUQetaPOI() const {return this->fHistProUQetaPOI;}
   TProfile* GetHistProUQetaAllEventsPOI() const {return this->fHistProUQetaAllEventsPOI;}
   TProfile* GetHistProUQPtRP() const   {return this->fHistProUQPtRP;} 
   TProfile* GetHistProUQPtPOI() const  {return this->fHistProUQPtPOI;}
   TProfile* GetHistProUQPtAllEventsPOI() const  {return this->fHistProUQPtAllEventsPOI;}
   TProfile* GetHistProQNorm() const    {return this->fHistProQNorm;}
   TProfile* GetHistProQaQb() const     {return this->fHistProQaQb;}
   TProfile* GetHistProQaQbNorm() const {return this->fHistProQaQbNorm;}
   TProfile* GetHistProQaQbReImNorm() const {return this->fHistProQaQbReImNorm;} 
   TProfile* GetHistProNonIsotropicTermsQ() const {return this->fHistProNonIsotropicTermsQ;} 
   TProfile* GetHistProNonIsotropicTermsU(Int_t rp, Int_t pe, Int_t sc) const {return this->fHistProNonIsotropicTermsU[rp][pe][sc];} 
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

   AliFlowCommonHist*        GetCommonHistsSP() const    {return this->fCommonHistsSP; }
   AliFlowCommonHistResults* GetCommonHistsResSP() const {return this->fCommonHistsResSP; }
   AliFlowCommonHist*        GetCommonHistsmuQ() const    {return this->fCommonHistsmuQ; }
   
   TH1D* GetHistQNorm() const      {return this->fHistQNorm;}
   TH1D* GetHistQaQb() const       {return this->fHistQaQb;}
   TH1D* GetHistQaQbNorm() const   {return this->fHistQaQbNorm;}
   TH2D* GetHistQNormvsQaQbNorm() const {return this->fHistQNormvsQaQbNorm;}
   TH1D* GetHistQaQbCos() const    {return this->fHistQaQbCos;}
   TH1D* GetHistResolution() const {return this->fHistResolution;}
   TH1D* GetHistQaNorm() const     {return this->fHistQaNorm;}
   TH2D* GetHistQaNormvsMa() const {return this->fHistQaNormvsMa;}
   TH1D* GetHistQbNorm() const     {return this->fHistQbNorm;}
   TH2D* GetHistQbNormvsMb() const {return this->fHistQbNormvsMb;}
   TH2D* GetMavsMb() const         {return this->fHistMavsMb;}

   //histogram setters   
   void SetHistProFlags(TProfile* const aHistProFlags) 
     {this->fHistProFlags = aHistProFlags;};  
   void SetHistProUQetaRP(TProfile* const aHistProUQetaRP) 
     {this->fHistProUQetaRP = aHistProUQetaRP;}
   void SetHistProUQetaPOI(TProfile* const aHistProUQetaPOI)
     {this->fHistProUQetaPOI = aHistProUQetaPOI;}
   void SetHistProUQetaAllEventsPOI(TProfile* const aHistProUQetaAllEventsPOI) 
     {this->fHistProUQetaAllEventsPOI = aHistProUQetaAllEventsPOI;}
   void SetHistProUQPtRP(TProfile* const aHistProUQPtRP)   
     {this->fHistProUQPtRP = aHistProUQPtRP;}
   void SetHistProUQPtPOI(TProfile* const aHistProUQPtPOI) 
     {this->fHistProUQPtPOI = aHistProUQPtPOI;}
   void SetHistProUQPtAllEventsPOI(TProfile* const aHistProUQPtAllEventsPOI)   
     {this->fHistProUQPtAllEventsPOI = aHistProUQPtAllEventsPOI;}
   void SetHistProQNorm(TProfile* const aHistProQNorm)
     {this->fHistProQNorm = aHistProQNorm;}
   void SetHistProQaQb(TProfile* const aHistProQaQb)       
     {this->fHistProQaQb = aHistProQaQb;}
   void SetHistProQaQbNorm(TProfile* const aHistProQaQbNorm)         
     {this->fHistProQaQbNorm = aHistProQaQbNorm;}     
   void SetHistProQaQbReImNorm(TProfile* const aHistProQaQbReImNorm)         
     {this->fHistProQaQbReImNorm = aHistProQaQbReImNorm;} 
   void SetHistProNonIsotropicTermsQ(TProfile* const aHistProNonIsotropicTermsQ)         
     {this->fHistProNonIsotropicTermsQ = aHistProNonIsotropicTermsQ;}           
   void SetHistProNonIsotropicTermsU(TProfile* const aHistProNonIsotropicTermsU, Int_t const rp, Int_t const pe, Int_t const sc)         
     {this->fHistProNonIsotropicTermsU[rp][pe][sc] = aHistProNonIsotropicTermsU;}          
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

   void SetCommonHistsSP(AliFlowCommonHist* const someCommonHists)              
     {this->fCommonHistsSP = someCommonHists; }
   void SetCommonHistsResSP(AliFlowCommonHistResults* const someCommonHistsRes) 
     {this->fCommonHistsResSP = someCommonHistsRes; }
   void SetCommonHistsmuQ(AliFlowCommonHist* const someCommonHists)              
     {this->fCommonHistsmuQ = someCommonHists; }
   

   void SetHistQNorm(TH1D* const aHistQNorm)
     {this->fHistQNorm = aHistQNorm;}
   void SetHistQaQb(TH1D* const aHistQaQb)
     {this->fHistQaQb = aHistQaQb;}
   void SetHistQaQbNorm(TH1D* const aHistQaQbNorm)
     {this->fHistQaQbNorm = aHistQaQbNorm;}
   void SetHistQNormvsQaQbNorm(TH2D* const aHistQNormvsQaQbNorm)
   {this->fHistQNormvsQaQbNorm = aHistQNormvsQaQbNorm;}
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
      
   Int_t      fEventNumber;           // event counter
   Bool_t     fDebug ;                // flag for analysis: more print statements
   Bool_t     fApplyCorrectionForNUA; // apply correction for non-uniform acceptance
   Int_t      fHarmonic;              // harmonic 
   TString   *fTotalQvector;          // total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"

   Double_t   fRelDiffMsub;      // the relative difference the two subevent multiplicities can have

   TList*     fWeightsList;      // list holding input histograms with phi weights
   Bool_t     fUsePhiWeights;    // use phi weights
   TH1F*      fPhiWeightsSub0;   // histogram holding phi weights for subevent 0
   TH1F*      fPhiWeightsSub1;   // histogram holding phi weights for subevent 1

   TProfile*  fHistProFlags;        // profile to hold all boolean flags needed in Finish()
   TProfile*  fHistProUQetaRP;      // uQ(eta) for RP (for events where both subevents are filled)
   TProfile*  fHistProUQetaPOI;     // uQ(eta) for POI (for events where both subevents are filled)
   TProfile*  fHistProUQetaAllEventsPOI; //uQ(eta) for POI (for events where 1 subevent may be empty)
   TProfile*  fHistProUQPtRP;       // uQ(pt) for RP (for events where both subevents are filled)
   TProfile*  fHistProUQPtPOI;      // uQ(pt) for POI (for events where both subevents are filled)
   TProfile*  fHistProUQPtAllEventsPOI;  // uQ(pt) for POI (for events where 1 subevent may be empty)
   TProfile*  fHistProQNorm;        // average of (Qa+Qb).Mod()
   TProfile*  fHistProQaQb;         // average of QaQb 
   TProfile*  fHistProQaQbNorm;     // average of QaQb/MaMb
   TProfile*  fHistProQaQbReImNorm; // average of Im[Qa/Ma], Re[Qa/Ma], Im[Qb/Mb], Re[Qb/Mb] 
   TProfile*  fHistProNonIsotropicTermsQ;          // 1st bin: sin, 2nd bin: cos 
   TProfile*  fHistProNonIsotropicTermsU[2][2][2]; // [RP/POI][pt/eta][sin/cos]  
   TH1D*      fHistSumOfLinearWeights;             // holds sum of Ma*Mb
   TH1D*      fHistSumOfQuadraticWeights;          // holds sum of (Ma*Mb)^2
   
   TProfile*  fHistProUQQaQbPtRP;         //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbEtaRP;        //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbPtPOI;        //holds weighted average of <QuQaQb>
   TProfile*  fHistProUQQaQbEtaPOI;       //holds weighted average of <QuQaQb>
   TH1D*      fHistSumOfWeightsPtRP[3];   //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsEtaRP[3];  //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsPtPOI[3];  //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
   TH1D*      fHistSumOfWeightsEtaPOI[3]; //holds sums of 0: Mq-1, 1: (Mq-1)^2, 2: (Mq-1)*Ma*Mb for each bin
    
   AliFlowCommonHist*        fCommonHistsSP;    // control histograms
   AliFlowCommonHistResults* fCommonHistsResSP; // results histograms
   AliFlowCommonHist*        fCommonHistsmuQ;    // control histograms
     
   TH1D*      fHistQNorm;        // distribution of (Qa+Qb)/(Ma+Mb)
   TH1D*      fHistQaQb;         // distribution of QaQb
   TH1D*      fHistQaQbNorm;     // distribution of QaQb/MaMb
   TH2D*      fHistQNormvsQaQbNorm; // distribution of (Qa+Qb)/(Ma+Mb) vs QaQb/MaMb
   TH1D*      fHistQaQbCos;      // distribution of the angle between Qa and Qb (from Acos (va*vb))
   TH1D*      fHistResolution;   // distribution of cos(2(phi_a - phi_b))
   TH1D*      fHistQaNorm;       // distribution of Qa/Ma
   TH2D*      fHistQaNormvsMa;   // distribution of Qa/Ma vs Ma
   TH1D*      fHistQbNorm;       // distribution of Qb/Mb
   TH2D*      fHistQbNormvsMb;   // distribution of Qb/Mb vs Mb
   TH2D*      fHistMavsMb;       // Ma vs Mb
      
   TList*     fHistList;         // list to hold all output histograms  

   ClassDef(AliFlowAnalysisWithScalarProduct,0)  // macro for rootcint
     };
 

#endif
