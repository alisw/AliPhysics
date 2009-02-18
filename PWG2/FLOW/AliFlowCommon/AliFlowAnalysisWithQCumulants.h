/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#ifndef ALIFLOWANALYSISWITHQCUMULANTS_H
#define ALIFLOWANALYSISWITHQCUMULANTS_H

#include "AliFlowCommonConstants.h"//needed as include

class TObjArray;
class TList;
class TFile;
class TGraph;

class TH1;
class TProfile;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowVector;

class AliFlowCommonHist;
class AliFlowCommonHistResults;

//================================================================================================================

class AliFlowAnalysisWithQCumulants{
 public:
  AliFlowAnalysisWithQCumulants();
  virtual ~AliFlowAnalysisWithQCumulants(); 
  
  virtual void Init();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();
  virtual void WriteHistograms(TString outputFileName);
 
//----------------------------------------------------------------------------------------------------------------
//                                            setters and getters                                                 
//----------------------------------------------------------------------------------------------------------------
  TList* GetHistList() const {return this->fHistList;} 
  
  void SetWeightsList(TList* wlist) {this->fWeightsList = wlist;}
  TList* GetWeightsList() const {return this->fWeightsList;}   
 
  void SetIntFlowResults(TH1D* const ifr) {this->fIntFlowResultsQC = ifr;};
  TH1D* GetIntFlowResults() const {return this->fIntFlowResultsQC;};
  
  void SetDiffFlowResults2nd(TH1D* const diff2nd) {this->fDiffFlowResults2ndOrderQC = diff2nd;};
  TH1D* GetDiffFlowResults2nd() const {return this->fDiffFlowResults2ndOrderQC;};
  
  void SetDiffFlowResults4th(TH1D* const diff4th) {this->fDiffFlowResults4thOrderQC = diff4th;};
  TH1D* GetDiffFlowResults4th() const {return this->fDiffFlowResults4thOrderQC;};
  
  void SetCovariances(TH1D* const cov) {this->fCovariances = cov;};
  TH1D* GetCovariances() const {return this->fCovariances;};
  
  void SetCommonHists2nd(AliFlowCommonHist* const ch2nd) {this->fCommonHists2nd = ch2nd;};
  AliFlowCommonHist* GetCommonHists2nd() const {return this->fCommonHists2nd;};
  
  void SetCommonHists4th(AliFlowCommonHist* const ch4th) {this->fCommonHists4th = ch4th;};
  AliFlowCommonHist* GetCommonHists4th() const {return this->fCommonHists4th;};
  
  void SetCommonHists6th(AliFlowCommonHist* const ch6th) {this->fCommonHists6th = ch6th;};
  AliFlowCommonHist* GetCommonHists6th() const {return this->fCommonHists6th;};
  
  void SetCommonHists8th(AliFlowCommonHist* const ch8th) {this->fCommonHists8th = ch8th;};
  AliFlowCommonHist* GetCommonHists8th() const {return this->fCommonHists8th;};
  
  void SetCommonHistsResults2nd(AliFlowCommonHistResults* const chr2nd) {this->fCommonHistsResults2nd = chr2nd;};
  AliFlowCommonHistResults* GetCommonHistsResults2nd() const {return this->fCommonHistsResults2nd;};
  
  void SetCommonHistsResults4th(AliFlowCommonHistResults* const chr4th) {this->fCommonHistsResults4th = chr4th;};
  AliFlowCommonHistResults* GetCommonHistsResults4th() const {return this->fCommonHistsResults4th;};
  
  void SetCommonHistsResults6th(AliFlowCommonHistResults* const chr6th) {this->fCommonHistsResults6th = chr6th;};
  AliFlowCommonHistResults* GetCommonHistsResults6th() const {return this->fCommonHistsResults6th;};
  
  void SetCommonHistsResults8th(AliFlowCommonHistResults* const chr8th) {this->fCommonHistsResults8th = chr8th;};
  AliFlowCommonHistResults* GetCommonHistsResults8th() const {return this->fCommonHistsResults8th;};
  
  void SetAverageMultiplicity(TProfile* const am) {this->fAvMultIntFlowQC = am;};
  TProfile* GetAverageMultiplicity() const {return this->fAvMultIntFlowQC;};
  
  void SetQvectorForEachEventX(TProfile* const qvfeex) {this->fQvectorForEachEventX = qvfeex;};
  TProfile* GetQvectorForEachEventX() const {return this->fQvectorForEachEventX;};

  void SetQvectorForEachEventY(TProfile* const qvfeey) {this->fQvectorForEachEventY = qvfeey;};
  TProfile* GetQvectorForEachEventY() const {return this->fQvectorForEachEventY;};
        
  void SetQCorrelations(TProfile* const QCorr) {this->fQCorrelations = QCorr;};
  TProfile* GetQCorrelations() const {return this->fQCorrelations;};
  
  void SetQProduct(TProfile* const qp) {this->fQProduct = qp;};
  TProfile* GetQProduct() const {return this->fQProduct;};
  
  void SetQVectorComponents(TProfile* const qvc) {this->fQvectorComponents = qvc;};
  TProfile* GetQVectorComponents() const {return this->fQvectorComponents;};
  
  void SetTwo1n1nPerPtBinRP(TProfile* const pb2PerPtBin1n1nRP) {this->f2PerPtBin1n1nRP = pb2PerPtBin1n1nRP;};
  TProfile* GetTwo1n1nPerPtBinRP() const {return this->f2PerPtBin1n1nRP;};
  
  void SetTwo2n2nPerPtBinRP(TProfile* const pb2PerPtBin2n2nRP) {this->f2PerPtBin2n2nRP = pb2PerPtBin2n2nRP;};
  TProfile* GetTwo2n2nPerPtBinRP() const {return this->f2PerPtBin2n2nRP;};
  
  void SetThree2n1n1nPerPtBinRP(TProfile* const pb3PerPtBin2n1n1nRP) {this->f3PerPtBin2n1n1nRP = pb3PerPtBin2n1n1nRP;};
  TProfile* GetThree2n1n1nPerPtBinRP() const {return this->f3PerPtBin2n1n1nRP;};
  
  void SetThree1n1n2nPerPtBinRP(TProfile* const pb3PerPtBin1n1n2nRP) {this->f3PerPtBin1n1n2nRP = pb3PerPtBin1n1n2nRP;};
  TProfile* GetThree1n1n2nPerPtBinRP() const {return this->f3PerPtBin1n1n2nRP;};
  
  void SetFour1n1n1n1nPerPtBinRP(TProfile* const pb4PerPtBin1n1n1n1nRP) {this->f4PerPtBin1n1n1n1nRP = pb4PerPtBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nPerPtBinRP() const {return this->f4PerPtBin1n1n1n1nRP;}; 

  void SetTwo1n1nPerEtaBinRP(TProfile* const pb2PerEtaBin1n1nRP) {this->f2PerEtaBin1n1nRP = pb2PerEtaBin1n1nRP;};
  TProfile* GetTwo1n1nPerEtaBinRP() const {return this->f2PerEtaBin1n1nRP;};
  
  void SetTwo2n2nPerEtaBinRP(TProfile* const pb2PerEtaBin2n2nRP) {this->f2PerEtaBin2n2nRP = pb2PerEtaBin2n2nRP;};
  TProfile* GetTwo2n2nPerEtaBinRP() const {return this->f2PerEtaBin2n2nRP;};
  
  void SetThree2n1n1nPerEtaBinRP(TProfile* const pb3PerEtaBin2n1n1nRP) {this->f3PerEtaBin2n1n1nRP = pb3PerEtaBin2n1n1nRP;};
  TProfile* GetThree2n1n1nPerEtaBinRP() const {return this->f3PerEtaBin2n1n1nRP;};
  
  void SetThree1n1n2nPerEtaBinRP(TProfile* const pb3PerEtaBin1n1n2nRP) {this->f3PerEtaBin1n1n2nRP = pb3PerEtaBin1n1n2nRP;};
  TProfile* GetThree1n1n2nPerEtaBinRP() const {return this->f3PerEtaBin1n1n2nRP;};
  
  void SetFour1n1n1n1nPerEtaBinRP(TProfile* const pb4PerEtaBin1n1n1n1nRP) {this->f4PerEtaBin1n1n1n1nRP = pb4PerEtaBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nPerEtaBinRP() const {return this->f4PerEtaBin1n1n1n1nRP;}; 
  
  void SetTwo1n1nPerPtBinPOI(TProfile* const pb2PerPtBin1n1nPOI) {this->f2PerPtBin1n1nPOI = pb2PerPtBin1n1nPOI;};
  TProfile* GetTwo1n1nPerPtBinPOI() const {return this->f2PerPtBin1n1nPOI;};
  
  void SetTwo2n2nPerPtBinPOI(TProfile* const pb2PerPtBin2n2nPOI) {this->f2PerPtBin2n2nPOI = pb2PerPtBin2n2nPOI;};
  TProfile* GetTwo2n2nPerPtBinPOI() const {return this->f2PerPtBin2n2nPOI;};
  
  void SetThree2n1n1nPerPtBinPOI(TProfile* const pb3PerPtBin2n1n1nPOI) {this->f3PerPtBin2n1n1nPOI = pb3PerPtBin2n1n1nPOI;};
  TProfile* GetThree2n1n1nPerPtBinPOI() const {return this->f3PerPtBin2n1n1nPOI;};
  
  void SetThree1n1n2nPerPtBinPOI(TProfile* const pb3PerPtBin1n1n2nPOI) {this->f3PerPtBin1n1n2nPOI = pb3PerPtBin1n1n2nPOI;};
  TProfile* GetThree1n1n2nPerPtBinPOI() const {return this->f3PerPtBin1n1n2nPOI;};
  
  void SetFour1n1n1n1nPerPtBinPOI(TProfile* const pb4PerPtBin1n1n1n1nPOI) {this->f4PerPtBin1n1n1n1nPOI = pb4PerPtBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nPerPtBinPOI() const {return this->f4PerPtBin1n1n1n1nPOI;}; 
  
  void SetTwo1n1nPerEtaBinPOI(TProfile* const pb2PerEtaBin1n1nPOI) {this->f2PerEtaBin1n1nPOI = pb2PerEtaBin1n1nPOI;};
  TProfile* GetTwo1n1nPerEtaBinPOI() const {return this->f2PerEtaBin1n1nPOI;};
  
  void SetTwo2n2nPerEtaBinPOI(TProfile* const pb2PerEtaBin2n2nPOI) {this->f2PerEtaBin2n2nPOI = pb2PerEtaBin2n2nPOI;};
  TProfile* GetTwo2n2nPerEtaBinPOI() const {return this->f2PerEtaBin2n2nPOI;};
  
  void SetThree2n1n1nPerEtaBinPOI(TProfile* const pb3PerEtaBin2n1n1nPOI) {this->f3PerEtaBin2n1n1nPOI = pb3PerEtaBin2n1n1nPOI;};
  TProfile* GetThree2n1n1nPerEtaBinPOI() const {return this->f3PerEtaBin2n1n1nPOI;};
  
  void SetThree1n1n2nPerEtaBinPOI(TProfile* const pb3PerEtaBin1n1n2nPOI) {this->f3PerEtaBin1n1n2nPOI = pb3PerEtaBin1n1n2nPOI;};
  TProfile* GetThree1n1n2nPerEtaBinPOI() const {return this->f3PerEtaBin1n1n2nPOI;};
  
  void SetFour1n1n1n1nPerEtaBinPOI(TProfile* const pb4PerEtaBin1n1n1n1nPOI) {this->f4PerEtaBin1n1n1n1nPOI = pb4PerEtaBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nPerEtaBinPOI() const {return this->f4PerEtaBin1n1n1n1nPOI;}; 
  
  void SetDirectCorrelations(TProfile* const dc) {this->fDirectCorrelations = dc;};
  TProfile* GetDirectCorrelations() const {return this->fDirectCorrelations;};
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
//----------------------------------------------------------------------------------------------------------------
 
 private:
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc);
  
  AliFlowTrackSimple* fTrack;                           //track
  TList*              fHistList;                        //list to hold all output histograms
  TList*              fDiffFlowList;                    //list to hold all histograms and profiles needed for differential flow
  TList*              fWeightsList;                     //list to hold all histograms with weights
    
  TProfile*           fAvMultIntFlowQC;                 //average selected multiplicity (for int. flow)
 
  TProfile*           fQvectorComponents;               //averages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, ...)
            
  TH1D*               fIntFlowResultsQC;                //integrated flow results from Q-cumulants
  TH1D*               fDiffFlowResults2ndOrderQC;       //differential flow results from 2nd order Q-cumulant
  TH1D*               fDiffFlowResults4thOrderQC;       //differential flow results from 4th order Q-cumulant
  TH1D*               fCovariances;                     //final results for covariances: 1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...
  
  TProfile*                  fQvectorForEachEventX;     //profile containing the x-components of Q-vectors from all events (to be removed)  
  TProfile*                  fQvectorForEachEventY;     //profile containing the y-components of Q-vectors from all events (to be removed)   
  TProfile*                  fQCorrelations;            //multi-particle correlations calculated from Q-vectors 
  TProfile*                  fQProduct;                 //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  
  TProfile*                  fDirectCorrelations;       //multi-particle correlations calculated with nested loop  

  //RP (Reaction Plane particles):  
  TProfile*                  fPtReq1nRP;                //real part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtImq1nRP;                //imaginary part of q-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtReq2nRP;                //real part of q-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fPtImq2nRP;                //imaginary part of q-vector evaluated in harmonic 2n for each pt-bin

  TProfile*                  f2PerPtBin1n1nRP;          //<<2'>>_{n|n} per pt-bin
  TProfile*                  f2PerPtBin2n2nRP;          //<<2'>>_{2n|2n} per pt-bin
  TProfile*                  f3PerPtBin2n1n1nRP;        //<<3'>>_{2n|n,n} per pt-bin
  TProfile*                  f3PerPtBin1n1n2nRP;        //<<3'>>_{n,n|2n} per pt-bin
  TProfile*                  f4PerPtBin1n1n1n1nRP;      //<<4'>>_{n,n|n,n} per pt-bin
  
  TProfile*                  fEtaReq1nRP;               //real part of q-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaImq1nRP;               //imaginary part of q-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaReq2nRP;               //real part of q-vector evaluated in harmonic 2n for each eta-bin
  TProfile*                  fEtaImq2nRP;               //imaginary part of q-vector evaluated in harmonic 2n for each eta-bin

  TProfile*                  f2PerEtaBin1n1nRP;         //<<2'>>_{n|n} per eta-bin
  TProfile*                  f2PerEtaBin2n2nRP;         //<<2'>>_{2n|2n} per eta-bin
  TProfile*                  f3PerEtaBin2n1n1nRP;       //<<3'>>_{2n|n,n} per eta-bin
  TProfile*                  f3PerEtaBin1n1n2nRP;       //<<3'>>_{n,n|2n} per eta-bin
  TProfile*                  f4PerEtaBin1n1n1n1nRP;     //<<4'>>_{n,n|n,n} per eta-bin  

  //POI (Particles Of Interest): 
  TProfile*                  fPtReq1nPrimePOI;          //real part of q'-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtImq1nPrimePOI;          //imaginary part of q'-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtReq2nPrimePOI;          //real part of q'-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fPtImq2nPrimePOI;          //imaginary part of q'-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fmPrimePerPtBin;           //number of particles selected both as POI and not as RP per each pt-bin
  
  TProfile*                  fPtReq1nPrimePrimePOI;     //real part of q''-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtImq1nPrimePrimePOI;     //imaginary part of q''-vector evaluated in harmonic n for each pt-bin
  TProfile*                  fPtReq2nPrimePrimePOI;     //real part of q''-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fPtImq2nPrimePrimePOI;     //imaginary part of q''-vector evaluated in harmonic 2n for each pt-bin
  TProfile*                  fmPrimePrimePerPtBin;      //number of particles selected both as RP and POI per each pt-bin
  
  TProfile*                  fEtaReq1nPrimePOI;         //real part of q'-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaImq1nPrimePOI;         //imaginary part of q'-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaReq2nPrimePOI;         //real part of q'-vector evaluated in harmonic 2n for each eta-bin
  TProfile*                  fEtaImq2nPrimePOI;         //imaginary part of q'-vector evaluated in harmonic 2n for each eta-bin
  TProfile*                  fmPrimePerEtaBin;          //number of particles selected both as POI and not as RP per each eta-bin
   
  TProfile*                  fEtaReq1nPrimePrimePOI;    //real part of q''-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaImq1nPrimePrimePOI;    //imaginary part of q''-vector evaluated in harmonic n for each eta-bin
  TProfile*                  fEtaReq2nPrimePrimePOI;    //real part of q''-vector evaluated in harmonic 2n for each eta-bin
  TProfile*                  fEtaImq2nPrimePrimePOI;    //imaginary part of q''-vector evaluated in harmonic 2n for each eta-bin
  TProfile*                  fmPrimePrimePerEtaBin;     //number of particles selected both as RP and POI in each eta-bin
 
  TProfile*                  f2PerPtBin1n1nPOI;         //<<2'>>_{n|n} per pt-bin
  TProfile*                  f2PerPtBin2n2nPOI;         //<<2'>>_{2n|2n} per pt-bin
  TProfile*                  f3PerPtBin2n1n1nPOI;       //<<3'>>_{2n|n,n} per pt-bin
  TProfile*                  f3PerPtBin1n1n2nPOI;       //<<3'>>_{n,n|2n} per pt-bin
  TProfile*                  f4PerPtBin1n1n1n1nPOI;     //<<4'>>_{n,n|n,n} per pt-bin

  TProfile*                  f2PerEtaBin1n1nPOI;        //<<2'>>_{n|n} per eta-bin
  TProfile*                  f2PerEtaBin2n2nPOI;        //<<2'>>_{2n|2n} per eta-bin
  TProfile*                  f3PerEtaBin2n1n1nPOI;      //<<3'>>_{2n|n,n} per eta-bin
  TProfile*                  f3PerEtaBin1n1n2nPOI;      //<<3'>>_{n,n|2n} per eta-bin
  TProfile*                  f4PerEtaBin1n1n1n1nPOI;    //<<4'>>_{n,n|n,n} per eta-bin  
 
  AliFlowCommonHist*         fCommonHists2nd;           //common control histograms (taking into account only the events with 2 and more particles) 
  AliFlowCommonHist*         fCommonHists4th;           //common control histograms (taking into account only the events with 4 and more particles) 
  AliFlowCommonHist*         fCommonHists6th;           //common control histograms (taking into account only the events with 6 and more particles) 
  AliFlowCommonHist*         fCommonHists8th;           //common control histograms (taking into account only the events with 8 and more particles) 
  
  AliFlowCommonHistResults*  fCommonHistsResults2nd;    //final results for 2nd order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults4th;    //final results for 4th order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults6th;    //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults*  fCommonHistsResults8th;    //final results for 8th order int. and diff. flow stored in the common histograms
      
  TH1D*                      f2pDistribution;           //distribution of <2>_{n|n}
  TH1D*                      f4pDistribution;           //distribution of <4>_{n,n|n,n}
  TH1D*                      f6pDistribution;           //distribution of <6>_{n,n,n|n,n,n} 
  TH1D*                      f8pDistribution;           //distribution of <8>_{n,n,n,n|n,n,n,n}
 
  Int_t                      fnBinsPt;                  //number of pt bins
  Double_t                   fPtMin;                    //minimum pt   
  Double_t                   fPtMax;                    //maximum pt    
  
  Int_t                      fnBinsEta;                 //number of eta bins
  Double_t                   fEtaMin;                   //minimum eta   
  Double_t                   fEtaMax;                   //maximum eta
  Int_t                      fEventCounter;             //counting the number of events    
   
  Bool_t                     fUsePhiWeights;            //phi weights
  Bool_t                     fUsePtWeights;             //v_2(pt) weights
  Bool_t                     fUseEtaWeights;            //v_2(eta) weights
                       
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





