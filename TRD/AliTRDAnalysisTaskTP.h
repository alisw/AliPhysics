#ifndef ALITRDANALYSISTASKTP_H
#define ALITRDANALYSISTASKTP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDAnalysisTaskTP.h 42548 2010-07-27 08:10:51Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Track point maker for the alignment of TRD                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

class TObjArray;
class AliVEvent;
class TTree;
class AliESDEvent;
class TH1D;
class TH2D;
class AliTrackPointArray;

class AliTRDAnalysisTaskTP:public AliAnalysisTaskSE {
public:
  AliTRDAnalysisTaskTP();
  AliTRDAnalysisTaskTP(const char *name);
  virtual ~AliTRDAnalysisTaskTP();

  //virtual void ConnectInputData(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

private:
  TObjArray          *fArrHists;         // histogram array 
  TObjArray          *fArrTTree;         // tree array
  TTree              *fTree;             //! alignment tree
  AliESDEvent        *fESD;              //! esd event
  TH2D               *fModpop;           //!
  TH1D               *fBug;              //!
  Int_t               fNevents;          //! number of events
  Int_t               fNtracks;          //! number of tracks
  Int_t               fNAcceptedTracks;  //! number of accepted tracks
  AliTrackPointArray *fArray;            //! pointer to the track points
  TFile              *fFile;             //! output file

  AliTRDAnalysisTaskTP(const AliTRDAnalysisTaskTP&);
  AliTRDAnalysisTaskTP& operator=(const AliTRDAnalysisTaskTP&);
  
  ClassDef(AliTRDAnalysisTaskTP,1)
};

#endif

