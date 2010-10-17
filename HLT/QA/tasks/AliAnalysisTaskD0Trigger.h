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

  void SingleTrackSelect(AliExternalTrackParam*, AliESDVertex*);
  void RecD0(Int_t&,AliESDVertex *,bool);

  //from AliD0toKpi
  Double_t InvMass(AliExternalTrackParam* d1, AliExternalTrackParam* d2);
  void cosThetaStar(AliExternalTrackParam* n, AliExternalTrackParam* p,Double_t &D0,Double_t &D0bar);
  Double_t pointingAngle(AliExternalTrackParam* n, AliExternalTrackParam* p, Double_t *pv, Double_t *sv);
  Double_t Pt(AliExternalTrackParam* d1, AliExternalTrackParam* d2);
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, const AliESDVertex *v, bool useKF);

  /// pt cut for decay, minimum [GeV/c]
  float fPtMin;                                           
  /// Distance between decay tracks [cm] ??
  float fdca;                                             
  /// Inv. mass half width [GeV]
  float finvMass;                                         
  /// Decay angle
  float fcosThetaStar;                                    
  /// Distance from primary vertex for decay tracks [cm] 
  float fd0;                                              
  /// Product of d0 for the two decay tracks [cm^2]
  float fd0d0;                                            
  /// Pionting angle
  float fcosPoint;                                        

  Double_t mD0PDG;                                        

  /// D0 inv. mass plot
  TH1F *fD0massHLT;                                       
  TH1F *fD0ptHLT;                                         
  TH1F *fD0massOFF;                                       
  TH1F *fD0ptOFF;                                         

  vector<AliExternalTrackParam*> fPos;                    
  vector<AliExternalTrackParam*> fNeg;                    

  TObjArray *ftwoTrackArray;                              

  Int_t fTotalD0HLT;                                         
  Int_t fTotalD0OFF;                                     
  Double_t fField;                                        

  Int_t fNevents;                                         
  
  bool fuseKF;                                            

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliAnalysisTaskD0Trigger, 1);

};

#endif
