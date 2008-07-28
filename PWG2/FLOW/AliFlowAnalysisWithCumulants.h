#ifndef AliFlowAnalysisWithCumulants_H
#define AliFlowAnalysisWithCumulants_H

//******************************* 
// flow analysis with cumulants *   
// author: Ante Bilandzic       * 
// email: anteb@nikhef.nl       *
//******************************* 

#include "AliFlowCommonConstants.h"
#include "AliFlowCumuConstants.h"

class TH1;
class TH3D;

class TObjArray;


class TProfile;
class TProfile2D;
class TProfile3D;


class TList;
class TFile;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

class AliFlowAnalysisWithCumulants {
 public:
  AliFlowAnalysisWithCumulants();
  virtual ~AliFlowAnalysisWithCumulants();
  
  virtual void CreateOutputObjects();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();

  // Output 
  
  // Output
  TList*   GetHistList() const             { return this->fHistList ; }     // Gets output histogram list //NEW

 private:
  AliFlowAnalysisWithCumulants(const AliFlowAnalysisWithCumulants& aAnalysis);
  AliFlowAnalysisWithCumulants& operator=(const AliFlowAnalysisWithCumulants& aAnalysis);

  AliFlowTrackSimple* fTrack;//track
  static const Int_t fgkQmax=AliFlowCumuConstants::kQmax;//needed for numerics
  static const Int_t fgkPmax=AliFlowCumuConstants::kPmax;//needed for numerics  
  static const Int_t fgkFlow=AliFlowCumuConstants::kFlow;//integrated flow coefficient to be calculated
  static const Int_t fgkMltpl=AliFlowCumuConstants::kMltpl;//the multiple in p=m*n (diff. flow) 
  static const Int_t fgknBins=100;//number of pt bins
 
  //TString  fHistFileName;      //! The output file name  // NEW: In the old version this was not commented out     
  //TFile* fHistFile; // histogram file for Cumulants // NEW: In the old version this was not commented out
  
  TList* fHistList; //list to hold all output histograms //NEW
  
  Double_t fAvM;//avarage SELECTED multiplicity

  Double_t fR0;//needed for numerics
  Double_t fPtMax;//maximum pt
  Double_t fPtMin;//minimum pt
  Double_t fBinWidth;//width of pt bin (in GeV)
      
  Double_t fAvQx;//<Q_x>
  Double_t fAvQy;//<Q_y>
  Double_t fAvQ2x;//<(Q_x)^2>
  Double_t fAvQ2y;//<(Q_y)^2>
   
  Int_t fNumberOfEvents;//number of events 
   
  
  TProfile*          fHistProAvM;       //avarage multiplicity
  
  TH1D*              fIntFlowResults;   //integrated flow final results
  
  TH1D*              fDiffFlowResults2; //differential flow final results //to be improved
  TH1D*              fDiffFlowResults4; //differential flow final results //to be improved
  TH1D*              fDiffFlowResults6; //differential flow final results //to be improved
  TH1D*              fDiffFlowResults8; //differential flow final results //to be improved
  
  TProfile2D*        fIntFlowGenFun;    // avarage of the generating function for integrated flow 
  TProfile3D*        fDiffFlowGenFunRe; // avarage of the generating function for differential flow (real part)
  TProfile3D*        fDiffFlowGenFunIm; // avarage of the generating function for differential flow (imaginary part)
 
  
  AliFlowCommonHist* fCommonHists;//control histograms
  //AliFlowCommonHistResults *fCommonHistsRes2, *fCommonHistsRes4, *fCommonHistsRes6, *fCommonHistsRes8;//final results//to be improved 
  
  Double_t fAvG[fgkPmax][fgkQmax];//avarage of the generating function used for integrated flow
  Int_t fBinEventEntries[fgknBins];//counts how many events have at least 1 particle in particular bin
  Int_t fBinNoOfParticles[fgknBins];//number of particles per bin
  Double_t fBinMeanPt[fgknBins];//mean pt per bin
  Double_t fBinEventDRe[fgknBins][fgkPmax][fgkQmax];//real part of the generating function used for differential flow
  Double_t fBinEventDIm[fgknBins][fgkPmax][fgkQmax];//imaginary part of the generating function used for differential flow
    
  ClassDef(AliFlowAnalysisWithCumulants, 0);
};
#endif





