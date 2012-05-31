/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//==============================================================================
// AliHMPIDAnalysysTask - Class representing a basic analysis tool of HMPID data at  
// level of ESD.
// A set of histograms is created.
//==============================================================================

#ifndef AliHMPIDPERFORMANCETASK_H
#define AliHMPIDPERFORMANCETASK_H

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class TObject;
class TH1;
class TH1F;
class TParticle;
class TFile;
class AliESDtrack;
class AliESDEvent;
class AliESDtrackCuts;
class AliHMPIDCluster;
class AliESDVertex;
class AliCentrality;    

class AliHMPIDPerformanceTask : public AliAnalysisTaskSE {
 public:

  enum {kChamber = 7};
  AliHMPIDPerformanceTask();
  AliHMPIDPerformanceTask(const Char_t* name);
  AliHMPIDPerformanceTask& operator= (const AliHMPIDPerformanceTask& c);
  AliHMPIDPerformanceTask(const AliHMPIDPerformanceTask& c);
  virtual ~AliHMPIDPerformanceTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void   SetUseMC(Bool_t useMC) { fUseMC = useMC; }
  Bool_t Equal(Double_t x, Double_t y, Double_t tolerance);
  
 protected:
  
 private:     
  
  AliESDEvent  *fESD;                                    //! ESD object
  AliESDfriend *fESDfriend;                              //! ESD Friend
  AliESDtrackCuts *fEsdTrackCuts;                        //! ESD track cuts
  AliHMPIDCluster* fClu;                                 //! HMPID photon cluster
  AliMCEvent   *fMC;                                     //! MC event
  AliESDVertex *fEsdVtx;			         //! ESD Vertex
  AliCentrality *fCentrality;                            //! Centrality
  AliESDtrack       *fEsdTrack;                          //! ESD track
  AliESDfriendTrack *fEsdFriendTrack;                    //! ESD Friend track
  TObject       *fCalibObject;                           //! Calib object = clu from friend
  
  Bool_t       fUseMC;                                   // decide whether use or not the MC information
  TList        *fHmpHistList ;                           // list of histograms
  TH1F         *fHmpNevents;                             // number of events processed
  TH1I          *fHmpNevPerTrigClass;                    // Number of events per offline trigger class
  ULong64_t    fGlobalEventNumber;                       // Global event number
  Int_t        fvType;                                   // Tree input type 0; track 1; cluster
  Char_t*      fvFiredTriggerClasses;                    //fired trigger classes
  Int_t        fvRunNumber;                              // run number
  Int_t        fvBunchCrossNumber;                       // bunch crossing number
  Int_t        fvOrbitNumber;                            // orbit number
  Int_t        fvPeriodNumber;                           // period number
  Double_t     fvMagField;                               // magnetic field
  Double_t     fvVertexX;                                // vertex x coord
  Double_t     fvVertexY;                                // vertex y coord
  Double_t     fvVertexZ;                                // vertex z coord
  Int_t        fvVertexNContributors;                    // vertex number of contributors
  Double_t     fvPesd;                                   // esd track momentum
  Double_t     fvPhmpMag;                                // hmp track momentum magnitude
  Double_t     fvPhmp[3];                                // hmp track momentum 
  Double_t     fvCentrality;                             // centrality percentile
  Int_t        fvAcceptedTracks;                         // Number of accepted tracks by AliESDtrackcuts
  Int_t        fvRefMultTpc;                             // Reference multiplicity by AliESDtrackcuts
  Double_t     fvHmpChi2;                                // hmp chi2
  Int_t        fvHmpCluIndx ;                            // hmp cluster index => chamber =  cluidx%1000000/1000
  Float_t      fvHmpMipX;                                // hmp mip x
  Float_t      fvHmpMipY;                                // hmp mip y
  Int_t        fvHmpMipQ;                                // hmpid mip q
  Int_t        fvHmpMipNPhots;                           // hmp mip number of photons
  Double_t     fvHmpPid[5];                              // hmp pid
  Double_t     fvHmpSignal ;                             // hmp signal 
  Float_t      fvHmpTrkX;                                // hmp track x
  Float_t      fvHmpTrkY;                                // hmp track y
  Float_t      fvHmpTrkTheta;                            // hmp track theta
  Float_t      fvHmpTrkPhi;                              // hmp track phi
  Double_t     fvHmpCluQ;                                // hmp cluster q
  Double_t     fvHmpCluX;                                // hmp cluster x
  Double_t     fvHmpCluY;                                // hmp cluster y
  Int_t        fvHmpCluCh;                               // hmp cluster ch
  Int_t        fvHmpCluSize;                             // hmp cluster size
  Int_t        fvHmpCluBox;                              // hmp cluster box
  Int_t        fvHmpCluStatus;                           // hmp cluster status
  Int_t        fvEsdTrackAccepted;                       // esd track is accepted by track cuts or not
  Int_t        fvKinkIndex;                              // esd track kink index
  Double_t     fvTofSignal;                              // tof signal
  TTree        *fTree;                                   // tree with useful data for subsequent analysis



  ClassDef(AliHMPIDPerformanceTask,1);
};

#endif
