#ifndef AliAnalysisTaskRLLYZNEWMETHOD_H
#define AliAnalysisTaskRLLYZNEWMETHOD_H

 
#include "AliAnalysisTaskRL.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TObjArray.h"
#include "TList.h"
#include "TProfile.h"

class AliFlowEvent;
class AliFlowTrack;
class AliFlowSelection;
class AliFlowMaker;
class AliESD;
class AliESDtrack;
 
class TFile;
class TTree;
class TObjArray;
class TH1F;
//class TProfile;


class AliAnalysisTaskRLLYZNewMethod : public AliAnalysisTaskRL {
 public:
  AliAnalysisTaskRLLYZNewMethod(const char *name);
  virtual ~AliAnalysisTaskRLLYZNewMethod();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  
 private:

  TFile*             fOutfile;         //! 
  TFile*             fFirstRunFile ;   //! pointer to file from first run
  TFile*             fSecondRunFile ;   //! pointer to file from second run
  AliESD*            fESD;             //! ESD object
  AliFlowEvent*      fFlowEvent;       //! flowevent object
  AliFlowTrack*      fFlowTrack;       //! 
  TObjArray*         fFlowTracks;      //! 
  AliFlowSelection*  fFlowSelect;      //! flowselection object
  AliFlowMaker*      fFlowMaker;       //! flowmaker object
  
  //histograms
  //input
  TProfile*  h1;    //!
  TProfile*  h2;    //!
  TH1F*      h3;    //!
  TProfile*  p1;    //!
  TProfile*  p2;    //!
  TProfile*  p3;    //!
  TProfile*  p4;    //!
  TProfile*  p5;    //!
  //output
  TProfile*  fHistProFlow;         //!
  TH1F*      fHistQtheta;          //!
  TProfile*  fHistProR0thetaHar2;  //!
  TProfile*  fHistProReDtheta;     //!
  TProfile*  fHistProImDtheta;     //!
  
   
  ClassDef(AliAnalysisTaskRLLYZNewMethod, 0);          // lyz analysis
};

 #endif
