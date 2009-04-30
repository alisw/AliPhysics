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
// SISCone (FastJet v2.3.4) finder algorithm interface
// Finder Header Class 
// Author: swensy.jangal@ires.in2p3.fr
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"

#include "AliSISConeJetHeader.h"

ClassImp(AliSISConeJetHeader)

////////////////////////////////////////////////////////////////////////

AliSISConeJetHeader::AliSISConeJetHeader():
    AliJetHeader("AliSISConeJetHeader"),
    fActiveAreaRepeats(1),
    fCaching(0),
    fConeRadius(0.4),
    fEffectiveRFact(1),
    fGhostArea(0.05),
    fGhostEtaMax(2.0),
    fMinJetPt(0),
    fNPassMax(0),
    fOverlapThreshold(0.5),
    fPhiMax(TMath::TwoPi()),
    fPhiMin(0),
    fRapMax(-0.9),
    fRapMin(0.9),
    fPtProtoJetMin(2),
    fSplitMergeScaleNumber(0),
    fSplitMergeStoppingScale(0)

{
  // Constructor
}

////////////////////////////////////////////////////////////////////////

void AliSISConeJetHeader::PrintParameters() const
{
  // prints out parameters of jet algorithm

  cout << "SISConeJet algorithm  parameters:"<<endl;

  cout<<"Cone Radius = "<<fConeRadius<<endl;
  cout<<"Overlap parameter = "<<fOverlapThreshold<<endl;
  cout<<"Maximum number of runs = "<<fNPassMax<<endl;
  cout<<"Pt min of protojets  = "<<fPtProtoJetMin<<endl;
  cout<<"Do we record cones of these events ? (0 = no, 1 = yes) = "<<fCaching<<endl;

  cout << "Background subtraction parameters :" <<endl;
  //cout<<"Kind of area used = "<<<<endl;
  cout<<"Eta max in which ghosts wil be generated = "<<fGhostEtaMax<<endl;
  cout<<"Ghost area = "<<fGhostArea<<endl;
  cout<<"Background will be studied in ["<<fRapMin<<","<<fRapMax<<"] in eta and ["<<fPhiMin<<","<<fPhiMax<<"] in phi"<<endl;
  //cout<<"Kind of recombination for split/merge procedure = "<<<<endl;
  //cout<<"Stopping scale for split/merge procedure = "<<<<endl;
  cout<<"Do we repeat active area calculus? (0 = no, 1 = yes) = "<<fActiveAreaRepeats<<endl;

  cout<<"Jets PtMin  = "<<fMinJetPt<<endl;
  
}
