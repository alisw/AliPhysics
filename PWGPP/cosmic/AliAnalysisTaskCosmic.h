#ifndef ALIANALYSISTASKCOSMIC_H
#define ALIANALYSISTASKCOSMIC_H 

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Analysis Task for the Quality Assurance of Cosmic Data
// Two track segments in the are matched in angle and charged. T
// The quality of the matching in is checked by comparing 
// the tarnsverse momenta and starting points of the track segments
//
// Author
// Andreas Morsch
// andreas.morsch@cern.ch



class TH1F;
class TH2F;

class TProfile;

class TList;
class TClonesArray;


#include "AliAnalysisTaskSE.h"

enum SelType {
   kPosC   = 0,
   kNegC   = 1,
   kPosZ   = 2,
   kNegZ   = 3,
   kGood   = 4,
   kBad    = 5
};   

class AliAnalysisTaskCosmic : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCosmic(const char *name = "AliAnalysisTaskCosmic");
  virtual ~AliAnalysisTaskCosmic() {}
  AliAnalysisTaskCosmic(const AliAnalysisTaskCosmic& task); 
  AliAnalysisTaskCosmic& operator=(const AliAnalysisTaskCosmic& task);
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 private:
  TList*          fHists;        // List of histograms
  TH1F*           fhPt[6];       // Pt distribution
  TH1F*           fhTheta[6];    // Eta distribution
  TH1F*           fhPhi[6];      // Phi distribution
  TH1F*           fhDPhi[6];     // DeltaPhi
  TH1F*           fhDTheta[6];   // DeltaTheta
  TH1F*           fhDZ[6];       // DeltaZ
  TH1F*           fhDX[6];       // DeltaX
  TH1F*           fhDY[6];       // DeltaY
  TH1F*           fhDPt[6];      // DeltaPt
  TH1F*           fhD1ovPt[6];   // Delta 1/Pt
  

  TH1F*           fhDPtovPt[6];  // DeltaPt/pt
  
      
  TH2F*           fhDZvsZ;       // dz vs z
  TH2F*           fhDZvsPhi;     // dz vs phi
    
  TH2F*           fhCh1Ch2;      // ch1 vs ch2
  TH2F*           fhPh1Ph2;      // phi1 vs phi2
  TH2F*           fhCl1Cl2G;     // #Clusters
  TH2F*           fhCl1Cl2B;     // #Clusters
    
  TProfile*       fpDPt[6];      // delta pt / <pt>
  TProfile*       fpDPtS[6];     // delta_pt / <error_pt>
  
  ClassDef(AliAnalysisTaskCosmic, 1); // example of analysis
};

#endif
