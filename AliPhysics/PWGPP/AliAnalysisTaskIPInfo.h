#ifndef ALIANALYSISTASKIPINFO_H
#define ALIANALYSISTASKIPINFO_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliESDfriend;
class AliESDEvent;
class AliESDVertex;
class AliIntSpotEstimator;

#include "AliAnalysisTask.h"

class AliAnalysisTaskIPInfo : public AliAnalysisTask 
{
 public:
  enum {kITSTPC,kTPC,kSPD,kNEst};
  //
  AliAnalysisTaskIPInfo(const char *name = "IPInfo");
			
  virtual ~AliAnalysisTaskIPInfo(); 
  //
  AliIntSpotEstimator* GetEstimator(Int_t i)                const {return i>=0&&i<kNEst ? fIPEst[i] : 0;}
  //
  void SetOptions(Int_t estID, Bool_t recoVtx=kFALSE, 
		  Double_t outcut=1e-4,	Int_t ntrIP=2,Int_t nPhiBins=12,Int_t nestb=1000,
		  Double_t estmin=-4e-2,Double_t estmax=6e-2,
		  Int_t ntrBins=10,Int_t ntMn=2,Int_t ntMx=32,
		  Int_t nPBins=14,Double_t pmn=0.2,Double_t pmx=3.,Bool_t fillNt=kFALSE);
  void SetIPCenIni(Int_t esdID, Double_t x=0,Double_t y=0,Double_t z=0);
  Int_t CreateSPDTracklets(TClonesArray& tracks);
  //
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //  
 protected:
  //
  // options for estimators creation
  Bool_t    fRecoVtx[kNEst];               //! request to refit the vertex for given estimator
  Int_t     fNTrMinIP[kNEst];              //! min tracks for IP estimator
  Int_t     fNPhiBins[kNEst];              //! n bins in phi for IP
  Int_t     fNEstb[kNEst];                 //! n of estimator bins
  Int_t     fNTrBins[kNEst];               //! n of vtx.mult. bins
  Int_t     fNPBins[kNEst];                //! n of track P bins
  Int_t     fNTrMin[kNEst];                //! min vtx multuplicity
  Int_t     fNTrMax[kNEst];                //! max vtx multuplicity
  Double_t  fOutCut[kNEst];                //! outliers cut level
  Double_t  fEstMin[kNEst];                //! lower estimator boundary
  Double_t  fEstMax[kNEst];                //! upper estimator boundary
  Double_t  fPMin[kNEst];                  //! lower P cut
  Double_t  fPMax[kNEst];                  //! upper P cut
  Double_t  fIPCenIni[kNEst][3];           //! initial estimate of IP Center
  Bool_t    fFillNt[kNEst];                //! request to fill ntuple
  //
  AliIntSpotEstimator* fIPEst[kNEst];     //! estimators
  AliESDEvent  *fESD;                     //! ESD object
  AliESDfriend *fESDfriend;               //! ESD friend object
  TList       *fOutput;                   //! list send on output slot 0
  TObjArray fTracks;                      //! temporary storage for extracted tracks
  static const Char_t* fEstNames[kNEst];  //  estimator names
  //
 private:    
  //
  AliAnalysisTaskIPInfo(const AliAnalysisTaskIPInfo&); // not implemented
  AliAnalysisTaskIPInfo& operator=(const AliAnalysisTaskIPInfo&); // not implemented
  AliESDVertex* ReconstructPrimaryVertexTPC() const;
  AliESDVertex* ReconstructPrimaryVertexITSTPC() const;
  
  ClassDef(AliAnalysisTaskIPInfo,1); // IP, vertexing and DCA resolutions analysis
};

#endif
