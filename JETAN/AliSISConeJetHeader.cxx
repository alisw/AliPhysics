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

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include "AliSISConeJetHeader.h"

ClassImp(AliSISConeJetHeader)

////////////////////////////////////////////////////////////////////////

AliSISConeJetHeader::AliSISConeJetHeader():
    AliJetHeader("AliSISConeJetHeader"),
    fActiveAreaRepeats(1),
    fAreaTypeNumber(4),
    fBGAlgo(1),
    fNHardJets(2),
    fBGMode(1),
    fCaching(0),
    fConeRadius(0.7),
    fEffectiveRFact(1),
    fGhostEtaMax(4.0),
    fGhostArea(0.05),
    fGridScatter(1),
    fKtScatter(0.1),
    fMeanGhostKt(1e-100),
    fMinJetPt(2),
    fNPassMax(0),
    fOverlapThreshold(0.75),
    fPhiMax(TMath::TwoPi()),
    fPhiMin(0),
    fPtProtoJetMin(2),
    fRapMax(0.9),
    fRapMin(-0.9),
    fRRho(0.5),
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
  if (fAreaTypeNumber == 1) cout<<"Kind of area used = Active area"<<endl;
  if (fAreaTypeNumber == 2) cout<<"Kind of area used = Active area explicit ghosts"<<endl;
  if (fAreaTypeNumber == 3) cout<<"Kind of area used = One ghost passive area"<<endl;
  if (fAreaTypeNumber == 4) cout<<"Kind of area used = Passive area"<<endl;
  if (fAreaTypeNumber == 5) cout<<"Kind of area used = Voronoi"<<endl;
  if (fBGAlgo == 0) cout<<"Algorithm for rho calculus = kT"<<endl;
  if (fBGAlgo == 1) cout<<"Algorithm for rho calculus = Cambridge"<<endl;
  cout<<"Eta max in which ghosts wil be generated = "<<fGhostEtaMax<<endl;
  cout<<"Ghost area = "<<fGhostArea<<endl;
  cout<<"Background will be studied in ["<<fRapMin<<","<<fRapMax<<"] in eta and ["<<fPhiMin<<","<<fPhiMax<<"] in phi"<<endl;
  cout<<"Kind of recombination for split/merge procedure = SM_pttilde"<<endl;
  cout<<"Stopping scale for split/merge procedure = "<<fSplitMergeStoppingScale<<endl;
  cout<<"Do we repeat active area calculus? (0 = no, 1 = yes) = "<<fActiveAreaRepeats<<endl;
  cout<<"Fractional random fluctuations of the position of the ghosts on the y-phi grid = "<<fGridScatter<<endl;       
  cout<<"Fractional random fluctuations of the tranverse momentum of the ghosts on the y-phi grid = "<<fKtScatter<<endl;         
  cout<<"Average transverse momentum of the ghosts = "<<fMeanGhostKt<<endl;       

  cout<<"Jets PtMin  = "<<fMinJetPt<<endl;

}
