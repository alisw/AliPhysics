// author: I.Arsene, i.c.arsene@gsi.de
// Analysis task for J/psi hadron angular correlations
// Task to be updated 
//
#ifndef ALIANALYSISTASKJPSICORRELATION_H
#define ALIANALYSISTASKJPSICORRELATION_H

#include "TList.h"

#include "AliAnalysisTaskMultiDielectron.h"

class AliESDEvent;
class AliESDtrackCuts;

class AliAnalysisTaskJpsiCorrelation : public AliAnalysisTaskMultiDielectron {
  
public:
  AliAnalysisTaskJpsiCorrelation();
  AliAnalysisTaskJpsiCorrelation(const char *name);
  virtual ~AliAnalysisTaskJpsiCorrelation(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  void SetESDCuts(AliESDtrackCuts* cuts)     {fESDTrackCuts    = cuts;}
  
 private:
  TList fTreesList;           //! list with the trees 
  AliESDEvent     *fESD;      //! ESD object
  AliESDtrackCuts *fESDTrackCuts;    //! ESD cuts  	

  Int_t fIdxDielectron;       //  dielectron index
  Int_t fNjpsiPerEvent;       //  candidates per event
  Int_t    fSign;             //  charge
  Double_t fJpsiM;            //  candidate mass
  Double_t fJpsiPt;           //  candidate transverse momentum
  Double_t fJpsiPhi;          //  candidate azimuth
  Double_t fJpsiTheta;        //  candidate polar angle
  Double_t fJpsiY;            //  candidate rapidity
  Double_t fTrackPt;          //  hadron pt
  Double_t fTrackPhi;         //  hadron phi
  Double_t fTrackTheta;       //  hadron theta
  Double_t fTrackEta;         //  hadron eta

  Int_t fMultiDieleOutputs;   //  number of outputs

  AliAnalysisTaskJpsiCorrelation(const AliAnalysisTaskJpsiCorrelation &c);
  AliAnalysisTaskJpsiCorrelation& operator= (const AliAnalysisTaskJpsiCorrelation &c);
  
  ClassDef(AliAnalysisTaskJpsiCorrelation, 1); //Analysis Task for J/psi - hadron correlations
};
#endif
