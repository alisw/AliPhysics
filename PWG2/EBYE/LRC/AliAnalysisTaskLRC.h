#ifndef ALIANALYSISTASKLRC_H
#define ALIANALYSISTASKLRC_H


// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This task is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data.  

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.0
// Version: 3.0.8  may 08


class TH1F;
class TH1D;
class AliESDEvent;
class TH2D;
class TProfile;
class TList;
#include "AliAnalysisTask.h"

class AliAnalysisTaskLRC : public AliAnalysisTask {

public:
 
 
  //Constructors
  
  AliAnalysisTaskLRC(const char *name = "AliAnalysisTaskLRC");
  virtual ~AliAnalysisTaskLRC() {}
  
  //AliAnalysisTask overloading
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // Setters 
  void SetForwardWindow(double StartETA,double EndETA);
  void SetBackwardWindow(double StartETA,double EndETA);
  void SetETAWindows(double _StartForwardETA,double _EndForwardETA,double _StartBakwardETA,double _EndBakwardETA);
  
  
  
  
private:
 
  AliESDEvent *fESD;    //ESD object
  
  // Total spectras (debugging)
  TH1F        *fHistPt; //Overal Pt spectrum
  TH1F        *fHistEta; //Overal Eta spectrum
  
 // Windows paramiters -----------------------------------
  
  double fStartForwardETA;  // Forward windos lover rapidity
  double fEndForwardETA;    // Forward window higer rapidity	
  double fStartBakwardETA;  // Bakward window lover rapidity
  double fEndBakwardETA;    // Bakward window higer rapidity

 
 
 //Otput List --------------------------------------------
  
  TList* fOutList;         // Output data container 
 
 // Output histogramms -----------------------------------

  TH2D* fHistNN;        // N-N 2D Profile
  TH2D* fHistPtN;	// Pt-N 2D Profile
  TH2D* fHistPtPt;	// Pt-Pt 2D Profile
  TH2D* fHistNberr;	// Nbackward error Profile
  TProfile* fProfdPtB;  // Used to store (in first bin) summary of PtB and its std diviation
  TProfile* fProfTestLRC; // Diognostic LRC Pt - N correlation

  // Supp. info for windows
  //Forward
  TH1D* fHistPtForward;   //Pt spectrum in Forward windows
  TH1D* fHistEtaForward;  //Eta spectrum in Forward windows
  TH1D* fHistNchForward;  //Nch spectrum in Forward windows
  
   //Bakward
  TH1D* fHistPtBakward;   //Pt spectrum in Bakward windows
  TH1D* fHistEtaBakward;  //Eta spectrum in Bakward windows
  TH1D* fHistNchBakward;  //Nch spectrum in Bakward windows
 
  
   
  AliAnalysisTaskLRC(const AliAnalysisTaskLRC&); // not implemented
  AliAnalysisTaskLRC& operator=(const AliAnalysisTaskLRC&); // not implemented
  
  ClassDef(AliAnalysisTaskLRC, 1); 
};

#endif
