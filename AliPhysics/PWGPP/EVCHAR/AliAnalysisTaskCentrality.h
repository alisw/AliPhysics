#ifndef ALIANALYSISTASKCENTRALITY_H
#define ALIANALYSISTASKCENTRALITY_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskCentrality
//   author: Alberica Toia
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TH2F;
class TProfile;

class AliAnalysisTaskCentrality : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskCentrality();
  AliAnalysisTaskCentrality(const char *name);
  AliAnalysisTaskCentrality& operator= (const AliAnalysisTaskCentrality& ana);
  AliAnalysisTaskCentrality(const AliAnalysisTaskCentrality& c);
  virtual ~AliAnalysisTaskCentrality();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
  

 private:

  Int_t    fDebug;	   	//  Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  //
  TList   *fOutput;	   	//! list send on output slot 0
  TH1F *fhEzdc;            //!
  TH1F *fhEzem;            //!
  TH1F *fhNtracks;         //!
  TH1F *fhNtracklets;      //!
  TH1F *fhNclusters0;      //!
  TH1F *fhmultV0;          //!
  TH1F *fhmultFMD;         //!

  TProfile *fhEzemvsEzdc;            //!
  TProfile *fhNtracksvsEzdc;         //!
  TProfile *fhNtrackletsvsEzdc;      //!
  TProfile *fhNclusters0vsEzdc;      //!
  TProfile *fhmultV0vsEzdc;          //!
  TProfile *fhmultFMDvsEzdc;         //!
  TProfile *fhNtracksvsEzem;         //!
  TProfile *fhNtrackletsvsEzem;      //!
  TProfile *fhNclusters0vsEzem;      //!
  TProfile *fhmultV0vsEzem;          //!
  TProfile *fhmultFMDvsEzem;         //!
  TProfile *fhNtracksvsmultV0;       //!      
  TProfile *fhNtrackletsvsmultV0;    //!    
  TProfile *fhNclusters0vsmultV0;    //!
  TProfile *fhNtracksvsmultFMD;      //!
  TProfile *fhNtrackletsvsmultFMD;   //!
  TProfile *fhNclusters0vsmultFMD;   //!
  TProfile *fhmultV0vsmultFMD;       //!

  Int_t    fNev;		//  event counter
  Float_t  fBeamEnergy;		//  beam energy
  Int_t    fNmyTracksgen;       //  no. generated primary charged tracks 
  char     fTrigClass[100];	//  fired trigger classes
  //
  Double_t fxVertex;		//  X vertex from ITS
  Double_t fyVertex;		//  Y vertex from ITS
  Double_t fzVertex;		//  Z vertex from ITS
  Bool_t   fVertexer3d;		//  Is vertex from 3d vertexer?
  //
  Double_t fbMC;		// impact parameter from MC
  Int_t    fNpartTargMC;	// no. of participants for target nucleus from MC
  Int_t    fNpartProjMC;	// no. of participants for project nucleus from MC
  Int_t    fNNColl;             // Number of N-N collisions
  Int_t    fNNwColl;            // Number of N-Nwounded collisions
  Int_t    fNwNColl;            // Number of Nwounded-N collisons
  Int_t    fNwNwColl;           // Number of Nwounded-Nwounded collisions
  //
  Int_t    fNTracklets;		//  no. tracklets
  Int_t    fNSingleClusters;	//  no. single clusters
  Int_t    fNClusters[6];	//  no. clusters on 6 ITS layers
  Int_t    fNChips[2];	        //  no. chips on 2 SPD layers
  //
  Double_t fbZDC;		// impact parameter from ZDC reco
  Int_t    fNpartZDC;           // no. of participants from ZDC reco
  Double_t fbZDCA;		// impact parameter from ZDC reco side A
  Int_t    fNpartZDCA;          // no. of part. from ZDC reco side A
  Double_t fbZDCC;		// impact parameter from ZDC reco side C
  Int_t    fNpartZDCC;          // no. of part. from ZDC reco side C
  //
  UInt_t   fESDFlag;		//  ZDC ESD flags
  Float_t  fZNCEnergy;		//  ZNC Energy
  Float_t  fZPCEnergy;		//  ZPC Energy
  Float_t  fZNAEnergy;		//  ZNA Energy
  Float_t  fZPAEnergy;		//  ZPA Energy
  Float_t  fZEM1Energy;		//  ZEM1 Energy
  Float_t  fZEM2Energy;		//  ZEM2 Energy
  Float_t  fZNCtower[5];	//  ZNC 5 tower signals
  Float_t  fZPCtower[5];	//  ZPC 5 tower signals
  Float_t  fZNAtower[5];	//  ZNA 5 tower signals
  Float_t  fZPAtower[5];	//  ZPA 5 tower signals
  Float_t  fCentrZNC[2];	//  centroid over ZNC
  Float_t  fCentrZNA[2];	//  centroid over ZNA
  
  Int_t    fNTracks;		//  no. tracks
  Int_t    fNPmdTracks;		//  no. PMD tracks
  Float_t  fMultV0A;		//  multiplicity from V0 reco side A
  Float_t  fMultV0C;		//  multiplicity from V0 reco side C
  Float_t  fMultFMDA;      //  multiplicity from FMD on detector A
  Float_t  fMultFMDC;      //  multiplicity from FMD on detector C

  ClassDef(AliAnalysisTaskCentrality,1); 
};
#endif
