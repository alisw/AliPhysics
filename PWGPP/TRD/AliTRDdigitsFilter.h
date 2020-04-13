/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author:             *
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
//
// The task:
//
//
//  Author:
//
//

#ifndef ALITRDDIGITSFILTER_H
#define ALITRDDIGITSFILTER_H

#include "AliTRDdigitsTask.h"

#include <vector>
#include <list>

class TTreeStream;
class AliInputEventHandler;
class AliTRDdigitsManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliESD;
class AliESDtrack;
class AliAnalysisTask;
class AliESDInputHandler;
class AliESDv0KineCuts;
class AliAnalysisManager;
class AliCentrality;
class TSystem;
class TStyle;
class TROOT;
class Riostream;
class TChain;
class TArrayF;
class TFile;
class THnSparse;
class TH2;
class TF1;
class TH1;
class TObjArray;


class AliTRDdigitsFilter : public AliTRDdigitsTask {

public:

  AliTRDdigitsFilter(const char *name = "trd_digits_filter");
  virtual ~AliTRDdigitsFilter();

  virtual void   UserCreateOutputObjects();
  //virtual Bool_t UserNotify();
  virtual void   UserExec(Option_t *);
  virtual void   Process();
  //virtual void   Terminate(const Option_t*);

  //AliESDv0KineCuts* GetV0cuts() {return fV0cuts;}

  void AcceptTracks(TString label, EPID_t pid,
      Float_t minPt, Float_t maxPt, Float_t fraction);

  void AcceptEvents(TString label,
                    Float_t minCent, Float_t maxCent, Float_t fraction);

  void PrintSettings();

  struct TrackCrit {
    TString fLabel;
    EPID_t fPid;
    Float_t fMinPt;
    Float_t fMaxPt;
    Float_t fFraction;
  };

  std::vector<TrackCrit> fTrackCriteria; // criteria to accept tracks

  struct EventCrit {
    TString fLabel;
    Float_t fMinCent;
    Float_t fMaxCent;
    Float_t fFraction;
  };

  std::vector<EventCrit> fEventCriteria; // criteria to accept events


protected:

  //AliESDv0KineCuts *fV0cuts;        //! ESD V0 cuts

  //std::vector<EPID_t> fPidTags;     //! vector of PID info for all tracks


  //Bool_t ReadDigits();
  //void WriteDigits();

  //void SetupV0qa();
  //void FillV0PIDlist();
  //void ClearV0PIDlist();

  Bool_t PassTrackCuts(AliESDtrack *fESDTrack=0,Int_t thres=0);
  Bool_t PassTrackPIDCuts(AliESDtrack *fESDTrack=0);

private:
  //
  //
  //AliESDEvent *fESDevent;              //! ESD object

  //TObjArray *fOutputContainer;         //! output data container
  //AliESDtrackCuts *fESDtrackCuts;      //! basic cut variables for all non-V0 tracks
  //AliESDtrackCuts *fESDtrackCutsV0;    //! basic cut variables for all V0 tracks

  //TList   *fListQA;                    //! List with filter QA histograms

  //TFile* fDigitsInputFile;             //! Digits file for reading
  //TFile* fDigitsOutputFile;            //! Digits file for writing

  //Int_t fEventNoInFile;                //! Bookkeeping

  //AliTRDdigitsManager* fDigMan;        //! digits manager

  // Histograms
  THnSparse*  fhAcc;                        //! summary hist of all acc tracks
  TH1F*       fhEventCuts;                  //! statistics of event cuts
  TH1F*       fhPtTag;                      //! pT of PID-tagged tracks
  TH1F*       fhPtGood;                     //! pT spectrum after quality cuts
  TH1F*       fhPtAcc;                      //! pT spectrum of accepted track
  TH1F*       fhCent;                       //! Centrality of event
  TH1F*       fhCentAcc;                    //! Accepted events based on centrality
  //-------------------------------------------------------------------------

  //TH1F *fhPt[fgkNSpecies];            //! pT spectrum for different species

  AliTRDdigitsFilter(const AliTRDdigitsFilter&); // not implemented
  AliTRDdigitsFilter& operator=(const AliTRDdigitsFilter&); // not implemented

  ClassDef(AliTRDdigitsFilter, 1);
};
#endif
