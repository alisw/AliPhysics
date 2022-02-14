#ifndef ALIANALYSISTASKJETCOREPP_H
#define ALIANALYSISTASKJETCOREPP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// This task performs hadron-trigger recoil jet correlations 
// Output pT spectrum of jet given trigger pT 
// Author: filip krizek 16th March 2013
// *******************************************

class TH1F;
class TH1D;
class TH1I;
class TList;
class TProfile;
class TFile;
class TKey;
class AliESDEvent;
class AliGenPythiaEventHeader;
class AliMCEvent;    //FK//
class AliMCEventHandler; //FK//
class AliGenEventHeader; //FK//

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"

class AliJPtHardXection : public AliAnalysisTaskSE {
public:
   AliJPtHardXection();
   AliJPtHardXection(const char *name);
   AliJPtHardXection(const AliJPtHardXection& a); 
   AliJPtHardXection& operator=(const AliJPtHardXection& a); // not implemented
   virtual ~AliJPtHardXection();
   virtual void  LocalInit() {Init();}
   virtual void  Init();
   virtual void  UserCreateOutputObjects();
   virtual void  UserExec(Option_t *option);
   virtual void  Terminate(const Option_t*);
   virtual Bool_t Notify();

private:
   //private member functions
   TList *fOutputList;          //! output data container 
 
   TProfile*     fh1Xsec;   //! pythia cross section and trials
   TProfile*     fevweight;   //! pythia cross section and trials
   TH1F* 	 fimpactpar; //! impact parameter
   TH1F*         fh1Trials; //! trials are added
   TH1F*         fh1PtHard;  //! Pt har of the event...      
   TH1F*         fh1PtHardNoW;  //! Pt har of the event without weigt      
   TH1F*         fh1PtHardTrials;  //! Number of trials
   Float_t       fAvgTrials;       // Average number of trials

   ClassDef(AliJPtHardXection, 1);  //has to end with number larger than 0
};

#endif

