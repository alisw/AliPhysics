// $Id$
//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

#ifndef ALIANALYSISTASKD0TRIGGER_H
#define ALIANALYSISTASKD0TRIGGER_H


/** @file AliAnalysisTaskD0Trigger.h
    @author Gaute Ovrebekk
    @date   
    @brief An analysis task for the D0 Trigger
*/

class TH1F;
class TList;
class TObjArray;
class AliESDVertex;
class AliExternalTrackParam;
class AliAODVertex;

#include "AliAnalysisTaskSE.h"
#include <vector>

class AliAnalysisTaskD0Trigger : public AliAnalysisTaskSE {

 public: 
  AliAnalysisTaskD0Trigger();
  AliAnalysisTaskD0Trigger(const char *name,float cuts[7]);
  virtual ~AliAnalysisTaskD0Trigger() {}
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  virtual void  NotifyRun();
  
 private:
  /** copy constructor */
  AliAnalysisTaskD0Trigger(const AliAnalysisTaskD0Trigger&); 
  /** assignment operator */
  AliAnalysisTaskD0Trigger& operator=(const AliAnalysisTaskD0Trigger&); 

  TList *fOutputList; // list of output histograms

  void SingleTrackSelect(AliExternalTrackParam* t, AliESDVertex* pV);
  void RecD0(Int_t& nD0,AliESDVertex *pV,bool isHLT);

  //from AliD0toKpi
  Double_t InvMass(AliExternalTrackParam* d1, AliExternalTrackParam* d2);
  void CosThetaStar(AliExternalTrackParam* n, AliExternalTrackParam* p,Double_t &D0,Double_t &D0bar);
  Double_t PointingAngle(AliExternalTrackParam* n, AliExternalTrackParam* p, Double_t *pv, Double_t *sv);
  Double_t Pt(AliExternalTrackParam* d1, AliExternalTrackParam* d2);
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, const AliESDVertex *v, bool useKF);

  float fPtMin;                 // pt cut for decay, minimum [GeV/c]
  float fdca;                   // Distance between decay tracks [cm] ??
  float finvMass;               // Inv. mass half width [GeV]  
  float fcosThetaStar;          // Decay angle                           
  float fd0;                    // Distance from primary vertex for decay tracks [cm] 
  float fd0d0;                  // Product of d0 for the two decay tracks [cm^2]
  float fcosPoint;              // Pionting angle                                      

  Double_t fD0PDG;              // Mass of D0 from PDG              

  
  TH1F *fD0massHLT;             // D0 inv. mass plot from HLT             
  TH1F *fD0ptHLT;               // D0 pT plot from HLT      
  TH1F *fD0massOFF;             // D0 inv. mass plot from offline
  TH1F *fD0ptOFF;               // D0 pT plot from offline

  vector<AliExternalTrackParam*> fPos;   // vector for positive tracks      
  vector<AliExternalTrackParam*> fNeg;   // vector for negative tracks

  TObjArray *ftwoTrackArray;        // Array for the two decay products                      

  Int_t fTotalD0HLT;            // Conter for numbers of D0 from HLT   
  Int_t fTotalD0OFF;            // Conter for numbers of D0 from offline   
  Double_t fField;              // Magnetic Field                          

  Int_t fNevents;               // Counter for number of events          
  
  bool fuseKF;                  // Bool for switching to KF

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliAnalysisTaskD0Trigger, 1);

};

#endif
