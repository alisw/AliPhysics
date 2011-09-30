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


//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
// Finder Header Class 
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include "AliFastJetHeaderV1.h"

ClassImp(AliFastJetHeaderV1)

////////////////////////////////////////////////////////////////////////

AliFastJetHeaderV1::AliFastJetHeaderV1():
    AliJetHeader("AliFastJetHeaderV1"),
    fRparam(0.4),
    fRparamBkg(0.4),  
    fAlgorithm(fastjet::kt_algorithm),
    fBGAlgorithm(fastjet::kt_algorithm),
    fStrategy(fastjet::Best),
    fRecombScheme(fastjet::BIpt_scheme),
    fGhostEtaMax(2.0),
    fGhostArea(0.01),
    fActiveAreaRepeats(1),
    fAreaType(fastjet::active_area), 
    fPtMin(5.0),
    fMinCellEt(0.),
    fRapMax(0.),
    fRapMin(0.),
    fPhiMax(TMath::TwoPi()),
    fPhiMin(0),
    fBGMode(0),
    fUse4VectorArea(kFALSE)
{
  // Constructor
  
  Double_t rapmax = fGhostEtaMax - fRparam;
  Double_t rapmin = -fGhostEtaMax + fRparam;
  SetRapRange(rapmin, rapmax);
  
}

////////////////////////////////////////////////////////////////////////

void AliFastJetHeaderV1::PrintParameters() const
{
  // prints out parameters of jet algorithm

  cout << "FastJet algorithm  parameters:" << endl;
  
  cout << "-- Jet Definition --- " << endl;
  cout << "R " << fRparam << endl;
  cout << "Jet Algorithm " << fAlgorithm << endl; 
  cout << "Bkg Jet Algorithm " << fBGAlgorithm << endl;
  cout << "Strategy " << fStrategy << endl;  
  cout << "Recombination Scheme " << fRecombScheme << endl; 
  
  cout << "-- Ghosted Area Spec parameters --- " << endl;
  cout << "Ghost Eta Max " << fGhostEtaMax << endl;
  cout << "Ghost Area " << fGhostArea << endl;
  cout << "Active Area Repeats " << fActiveAreaRepeats << endl;
  
  cout << "-- Area Definition parameters --- " << endl;
  cout << "Area Type " << fAreaType << endl; 
  
  cout << "-- Cluster Sequence Area parameters --- " << endl;
  cout << "pt min " << fPtMin << endl; 
  
  cout << "-- Range Definition parameters --- " << endl;
  cout << " bkg rapidity range from  " << fRapMin << " to " << fRapMax << endl;
  cout << " bkg phi range from " << fPhiMin << " to " << fPhiMax << endl;

}
