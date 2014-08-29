// Jet v2 task using QA method, based on jet sample task (S.Aiola).
//
// Authors: Jason Mueller (CERN summer student 2014) & Alice Ohlson

#ifndef ALIANALYSISTASKEMCALJETVNQA_H
#define ALIANALYSISTASKEMCALJETVNQA_H

// $Id$

class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetv2QA : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetv2QA();
  AliAnalysisTaskEmcalJetv2QA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetv2QA();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetCentBins(Int_t n, Double_t* bins);
  void                        SetJetPtBins(Int_t n, Double_t* bins);

  Double_t		      GetJetv2(){return fJetv2;};
  void			      SetJetv2(Double_t jetv2){fJetv2 = jetv2;};
  Bool_t		      GetDoPtWeight(){return doPtWeight;};
  void			      SetDoPtWeight(Bool_t ptWeight){doPtWeight = ptWeight;};

 protected:
  Int_t			      nCentBins;
  Int_t                       nCentBins1;
  Double_t*		      centBins; //[nCentBins1]
  Int_t			      nJetPtBins;
  Int_t			      nJetPtBins1;
  Double_t*		      jetPtBins; //[nJetPtBins1]
  void                        ExecOnce();
  Bool_t                      FillHistograms();
  Bool_t                      Run();
  Double_t		      fJetv2;
  Bool_t		      doPtWeight;

  // General histograms
  TH1F                       *fHistTracksPt;            //!Track pt spectrum
  TH1F                       *fHistClustersPt;          //!Cluster pt spectrum
  TH1F                       *fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH1F                       *fHistLeadingJetPtCorr;    //!Leading jet pt spectrum, background subtracted
  TH2F                       *fHistJetsPhiEta;          //!Phi-Eta distribution of jets
  TH2F                       *fHistJetsPtArea;          //!Jet pt vs. area
  TH2F                       *fHistJetsPtLeadHad;       //!Jet pt vs. leading hadron
  TH2F                       *fHistJetsCorrPtArea;      //!Jet pt - bkg vs. area
  TH3F                       *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3F                       *fHistPtDEtaDPhiClusTrack; //!cluster pt, delta eta, delta phi to matched track
  TH1F                       *fDPhiJet;                 //!
  TH1F                       *fDPhiJetPythia;           //!
  TH1F                       *fDPhiEP;                  //!
  TH2D 			     *hGx;		        //!
  TH2D                       *hGy2;                     //!
  TH2D                       *hGxGy2;                   //!
  TH2D                       *hGy4;                     //!
  TH2D                       *hGx2;                     //!
  TH2D                       *hGx2Gy2;                  //!
  TH2D                       *hGxGy4;                   //!
  TH2D                       *hGy6;                     //!
  TH2D                       *hGx2Gy4;                  //!
  TH2D                       *hGxGy6;                   //!
  TH2D                       *hGy8;                     //!
  TH2D                       *hGy;                      //!
  TH2D                       *hN;                       //!
  TH2D                       *htv2std;                  //!
  TH2D                       *htjv2std;                 //!
  TH2D                       *htj2v2std;                //!
  TH2D                       *hV0jv2std;                //!
  TH3D                       *htdPsi;                   //!
  TH3D                       *htjdPsi;                  //!
  TH3D                       *htj2dPsi;                 //!
  TH3D                       *hV0jdPsi;                 //!
  TH2D		 	     *hAx;		        //!
  TH2D			     *hAxDijet;		        //!
  TH2D			     *hQx;		        //!
  TH2D			     *hQy;		        //!
  TH1F			     *hEventData;	        //!
  TH2F			     *hNTracks;		        //!
  TH2D			     *hNTracksCent;		//!
  TH2D                       *hGxTracks;                //!
  TH2D                       *hGyTracks;                //!
  TH2D                       *hGy2Tracks;               //!
  TH2D                       *hGxGy2Tracks;             //!
  TH2D                       *hGy4Tracks;               //!
  TH2D                       *htEPRes;                  //!
  TH2D                       *htjEPRes;                 //!
  TH2D                       *htj2EPRes;                //!

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters



 private:
  AliAnalysisTaskEmcalJetv2QA(const AliAnalysisTaskEmcalJetv2QA&);            // not implemented
  AliAnalysisTaskEmcalJetv2QA &operator=(const AliAnalysisTaskEmcalJetv2QA&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetv2QA, 3) // jet v2QA analysis task
};
#endif
