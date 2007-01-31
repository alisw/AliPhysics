#ifndef ALITOFQATASK_H
#define ALITOFQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the TOF performance in simulated data
// 
//*-- Silvia Arcelli 
//-Distributions to monitor the quality of the track matching, Time, PID
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include <TH1F.h> 
#include <TH2F.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 

class AliTOFQATask : public AliAnalysisTask {

public:
  AliTOFQATask(const char *name) ; //ctor
  AliTOFQATask(const AliTOFQATask & qatask); // copy constructor
  AliTOFQATask& operator=(const AliTOFQATask & qatask); // assignment operator
  virtual ~AliTOFQATask(); //dtor
  virtual void Exec(Option_t * opt = "") ;
  virtual void Init(Option_t * opt = "") ; 
  virtual void Terminate(Option_t * opt = "") ;

private:
  void BookHistos() ;
  void DrawHistos() ;
  void GetEfficiency() ;
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TH1F    * fhTOFMatch ; //Fraction of Matched tracks (kTIME)
  TH1F    * fhESDeffPhi ;   //Phi dst. of all ESD tracks (kTIME)
  TH1F    * fhESDeffTheta ; //Theta dst. of all ESD tracks (kTIME)
  TH1F    * fhESDeffMom ;   //Mom. dst. of all ESD tracks (kTIME)
  TH1F    * fhTOFeffPhi ;   //Phi dst.of ESD tracks matched with TOF
  TH1F    * fhTOFeffTheta ; //Theta dst.of ESD tracks matched with TOF
  TH1F    * fhTOFeffMom ;   //Mom. dst.of ESD tracks matched with TOF
  TH1F    * fhTOFeffPhiMT ; //Phi dst.of ESD tracks, good match with TOF(MC)
  TH1F    * fhTOFeffThetaMT;//Theta dst.of ESD tracks, good match with TOF(MC)
  TH1F    * fhTOFeffMomMT ; //Mom.dst.of ESD tracks, good match with TOF(MC)
  TH1F    * fhTOFsector ;   //sector distr.of matched TOF pads 
  TH1F    * fhTOFsectorMT ; //sector distr.of TOF pads, good match (MC)  
  TH1F    * fhTOFTime ; //Time Distribution of TOF Clusters matched to tracks 
  TH1F    * fhTOFDeltaTime; //TOF-exp.Time (pion mass hypothesis)
  TH1F    * fhTOFDeltaTimeMT; //TOF-exp.Time (pion mass), good match(MC)
  TH1F    * fhTOFIDSpecies; //ID-Sample Composition
  TH2F    * fhTOFMassVsMom; //Mass vs Momentum correlation
  TH1F    * fhTOFMass;      //reconstructed Mass from TOF
   
  ClassDef(AliTOFQATask, 0); //  TOF Quality Assurance analysis task 
};
#endif // ALITOFQATASK_H
