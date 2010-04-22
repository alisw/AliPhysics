/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//
// Implementation Class AliJCorranSetup
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE


#include  "AliJCORRANSetup.h" 
#include <iostream>
#include <fstream>

using namespace std;

//___________________________________________________________________

AliJCORRANSetup::AliJCORRANSetup():
  // GENERAL SELECTIONS:
  fSpecies(" "),
  fVtxRange(-999),
  // Trigger SELECTIONS:
  fOffLineTriggThreshold(-999),
  fOffLinePi0TriggThreshold(-999),
  fOffLineHadronTriggThreshold(-999),
  fOffLineAccHadronTriggThreshold(-999),
  fMinBiasScaleFact(-999),
  // CLUSTER SELECTIONS:
  fClusEnerMin(-999),
  fClusEnerMax(-999),
  fChi2Cut(-999),
  fProbPhotCut(-999),
  fClusNTowMin(-999),
  fClusNTowMax(-999),
  // TRACK SELECTIONS:
  fTrkPtMin(-999),
  fTrkPtMax(-999),
  // External file
  fWarnmapfile(" "),
  fComment(" ")
{
  //constructor
}

//___________________________________________________________________

AliJCORRANSetup::AliJCORRANSetup(const char *SelectionsFile):
  // GENERAL SELECTIONS:
  fSpecies(" "),
  fVtxRange(-999),
  // Trigger SELECTIONS:
  fOffLineTriggThreshold(-999),
  fOffLinePi0TriggThreshold(-999),
  fOffLineHadronTriggThreshold(-999),
  fOffLineAccHadronTriggThreshold(-999),
  fMinBiasScaleFact(-999),
  // CLUSTER SELECTIONS:
  fClusEnerMin(-999),
  fClusEnerMax(-999),
  fChi2Cut(-999),
  fProbPhotCut(-999),
  fClusNTowMin(-999),
  fClusNTowMax(-999),
  // TRACK SELECTIONS:
  fTrkPtMin(-999),
  fTrkPtMax(-999),
  // External file
  fWarnmapfile(" "),
  fComment(" ")
{
    //default constructor
    //
    //STEP 1: Get the input file that contains the correlation setup
    //
    ifstream selfile(SelectionsFile,ios::in);
    if(!selfile)
    {
	cerr << endl << "AliJCORRANSetup::AliJCORRANSetup" << endl
	    << "ERROR: Cannot open selections file "
	    << SelectionsFile << endl;
    }
    else
    {    
	selfile >> fComment
	    >> fSpecies      >> fComment
	    >> fVtxRange     >> fComment
	    >> fWarnmapfile  >> fComment
	    >> fComment
	    >> fOffLineTriggThreshold >> fComment
	    >> fOffLinePi0TriggThreshold >> fComment
	    >> fOffLineHadronTriggThreshold >> fComment
	    >> fOffLineAccHadronTriggThreshold >> fComment
	    >> fMinBiasScaleFact >> fComment
	    >> fClusEnerMin  >> fComment
	    >> fClusEnerMax  >> fComment
	    >> fChi2Cut      >> fComment
	    >> fProbPhotCut  >> fComment
	    >> fClusNTowMin  >> fComment
	    >> fClusNTowMax  >> fComment
	    >> fComment
	    >> fTrkPtMin     >> fComment
	    >> fTrkPtMax     >> fComment;
	cout<< "Selections Loaded from File:            " << SelectionsFile << endl;
    }

}

void AliJCORRANSetup::PrintOut(){
    //print object
    cout << "*******************************************************************"
	<< endl
	<< "          CORRELATION SELECTIONS:"        << endl 
	<< "General Selections:"                      << endl
	<< "Colliding fSpecies:                      " << fSpecies << endl
	<< "Vertex Range:                           " << fVtxRange << endl
        <<  "WarnMap File : "<< fWarnmapfile << endl
	<< "*******************************************************************"
	<< endl;
    cout<< "          EmCal Cluster Selections:"      << endl
	<<  "fOffLineTriggThreshold   " << fOffLineTriggThreshold << endl 
	<<  "fOffLinePi0TriggThreshold   " << fOffLinePi0TriggThreshold <<  endl
	<<  "fOffLineHadronTriggThreshold   " << fOffLineHadronTriggThreshold <<  endl
	<<  "fOffLineAccHadronTriggThreshold   " << fOffLineAccHadronTriggThreshold <<  endl
	<<  "fMinBiasScaleFact  " << fMinBiasScaleFact <<endl
	<< "EmCal Cluster Energy Range:             " << fClusEnerMin
	<< "->"                                       << fClusEnerMax << endl
	<< "Cluster Shape (Chi2 Cut):               " << fChi2Cut << endl;
    cout<< "Photon probability (ProbPhot Cut):      " << fProbPhotCut << endl;
    cout<< "Cluster Number Of Towers:               " << fClusNTowMin
	<< "->"                                       << fClusNTowMax << endl
	<<"********************************************************************"
	<< endl
	<< "          Track Selections:"              << endl
	<< "Track Transverse Momentum:              " << fTrkPtMin
	<< "->"                                       << fTrkPtMax << endl
	<< "***************************************************************"
	<< endl;
}

