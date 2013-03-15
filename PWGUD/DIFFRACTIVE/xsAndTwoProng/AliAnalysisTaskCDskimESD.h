/*************************************************************************
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
//
// $Id$
//
// Task to skim ESD files, basic idea found at PWGGA/EMCALTasks as
// AliEsdSkimTask. This version is modified to produce a exact copy of the ESD
// stream containing only events which are of interest for the central
// diffractive analysis (PWGUD) and a control sample. In order to be able to
// process MC data as well, the usual concept using an output slot is for
// output data is violated and an output file is directly created without
// a corresponding output slot.
//
// Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALIANALYSISTASKCDSKIMESD_H
#define  ALIANALYSISTASKCDSKIMESD_H

#include "AliAnalysisTaskSE.h"

// $Id$

class TTree;
class TH1;
class TFile;
class TDirectory;

class AliAnalysisTaskCDskimESD : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskCDskimESD(const char *opt, Bool_t reduceGapEvents);
	AliAnalysisTaskCDskimESD();
	virtual ~AliAnalysisTaskCDskimESD();

	void UserExec(Option_t *opt);
	void UserCreateOutputObjects();

	void SetGapCondition(Int_t cond) { fRequestedGapCondition=cond; }
	void SetRefinedGapCondition(Int_t cond) {fRefinedGapCondition=cond; }

protected:
	// variables concerning the reconstructed information in the ESDs
	TTree*      fInputTree;             // input tree
	TTree*      fInputTreeCopy;         // input tree copy
	TTree*      fOutputTree;            // output tree

	// variables needed for the reproduction of the MC information
	Bool_t      fDoMC;                  // do we process MC information?
	TDirectory* fWrkDir;                // directory the analysis task used before
	                                    // it was manipulated to store the
	                                    // kinemactics truth
	TFile*      fKinematicsFile;        // output file containing the skimmed
	                                    // kinematics information
	TFile*      fTrackRefsFile;         // output file containing the skimmed
	                                    // TrackRefs information
	TTree*      fTEinput;               // input tree containing MC eventt headers
	TTree*      fTE;                    // tree containing the MC event headers
	TFile*      fgaliceFile;            // pointer to the galice.root output file

	// variables concerning the gaps used for the selection
	Bool_t      fDoGapCond;             // if true the gap condition for central
	                                    // diffractive analysis are stored
	Int_t       fCurrentGapCondition;   // gap condition
	Int_t       fRequestedGapCondition; // gap condition required for selection
	Double_t    fNoGapFraction;         // fraction of no-gap events to be stored
	Bool_t      fReduceGapEvents;       // do not store all V0 double-gap events
	                                    // necessary due to high double gap rates
	                                    // in MC productions (LHC10d{1,2,4,4a})
	Int_t       fRefinedGapCondition;   // gap condition with more stringent
	                                    // selection applied in case
	                                    // fReduceDgEvents is active
	Double_t    fReducedGapFraction;    // fraction of events stored altough they
	                                    // show only the requested and not the
	                                    // refined gap

	// selection monitoring
	TH1*        fStatsFlow;             // histogram indicating the number of
	                                    // input events
	TTree*      fSkimmingList;          // which events of the intial ones are
	                                    // skimmed?
	TString     fFileName;              // name of the input file
	Int_t       fEventNumberInFile;     // event index in the ESD file
	Int_t       fRunNumber;             // run number of the event
	UInt_t      fEventTime;             // event timestamp

private:
	//void CopygALICE(AliAnalysisManager* am); // copy objects in galice.root
	// not implemented:
	AliAnalysisTaskCDskimESD(const AliAnalysisTaskCDskimESD&);
	AliAnalysisTaskCDskimESD &operator=(const AliAnalysisTaskCDskimESD&);
	ClassDef(AliAnalysisTaskCDskimESD, 1); // Esd skimming task
};
#endif // ALIANALYSISTASKCDSKIMESD_H
