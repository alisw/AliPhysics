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


//////////////////////////////////////////////////////
//                                                  //
//     Class to specify cuts for track analysis     //
//     with AliTPCcalibTracks                       //
//                                                  //
//////////////////////////////////////////////////////

#include <iostream>
#include <TString.h>
#include <TChain.h>
#include <TList.h>
#include "AliTPCseed.h"
#include "AliVTrack.h"
#include "AliTPCcalibTracksCuts.h"

ClassImp(AliTPCcalibTracksCuts)


AliTPCcalibTracksCuts::AliTPCcalibTracksCuts():
  TNamed("calibTracksCuts", "calibTracksCuts"),
  fMinClusters(0),            // number of clusters
  fMinRatio(0),               // kMinRratio = 0.4
  fMax1pt(0),                 // kMax1pt = 0.5
  fEdgeYXCutNoise(0),         // kEdgeYXCutNoise = 0.13
  fEdgeThetaCutNoise(0)      // kEdgeThetaCutNoise = 0.018
{
   // 
   // default constructor
   // 
}


AliTPCcalibTracksCuts::AliTPCcalibTracksCuts(Int_t minClusters, Float_t minRatio, Float_t max1pt,
					     Float_t edgeXZCutNoise, Float_t edgeThetaCutNoise):
      TNamed("calibTracksCuts", "calibTracksCuts"),
      fMinClusters(minClusters),            // number of clusters
      fMinRatio(minRatio),                  // kMinRratio = 0.4
      fMax1pt(max1pt),                      // kMax1pt = 0.5
      fEdgeYXCutNoise(edgeXZCutNoise),      // kEdgeYXCutNoise = 0.13
      fEdgeThetaCutNoise(edgeThetaCutNoise)   // kEdgeThetaCutNoise = 0.018
{
   //
   // Constructor for AliTPCcalibTracksCuts
   // specify the cuts to be set on the processed tracks
   // default cuts are for comics
   //
}

AliTPCcalibTracksCuts::AliTPCcalibTracksCuts(AliTPCcalibTracksCuts *cuts):
  TNamed(cuts->GetName(), cuts->GetTitle()),
  fMinClusters(cuts->GetMinClusters()),             // number of clusters
  fMinRatio(cuts->GetMinRatio()),                   // kMinRratio = 0.4
  fMax1pt( cuts->GetMax1pt()),                      // kMax1pt = 0.5
  fEdgeYXCutNoise(cuts->GetEdgeYXCutNoise()),       // kEdgeYXCutNoise = 0.13
  fEdgeThetaCutNoise( cuts->GetEdgeThetaCutNoise())   // kEdgeThetaCutNoise = 0.018
{
  // 
  // copy constructor
  // 
}



AliTPCcalibTracksCuts::~AliTPCcalibTracksCuts(){
  //
  // Destructor
  //
  cout << "AliTPCcalibTracksCuts destructor called, nothing happend." << endl;
}




 AliTPCcalibTracksCuts  * AliTPCcalibTracksCuts::CreateCuts(char* ctype){
   // 
   // Create predefined cuts 
   // (creates AliTPCcalibTracksCuts object)
   // 
   // The following predefined sets of cuts can be selected:
   // laser:      20, 0.4, 0.5, 0.13, 0.018
   // cosmic:     20, 0.4, 0.5, 0.13, 0.018
   // lowflux:    20, 0.4, 5, 0.2, 0.0001
   // highflux:   20, 0.4, 5, 0.2, 0.0001
   // 
   
   TString cutType(ctype);
   cutType.ToUpper();
   AliTPCcalibTracksCuts *cuts = 0;
   if (cutType == "LASER")
     cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.13, 0.018);
   else if (cutType == "COSMIC")
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
   else if (cutType == "LOWFLUX")
     cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
   else if (cutType == "HIGHFLUX")
     cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
   else {
     cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
     cerr << "WARNING! unknown type '" << ctype << "', cuts set to default values for cosmics." << endl;
     cutType = "COSMIC";
   }
   cout << "Cuts were set to predefined set: " << cutType << endl;
   return cuts;
}



Int_t AliTPCcalibTracksCuts::AcceptTrack(const AliTPCseed * track) const {
  //
  // Function, that decides wheather a given track is accepted for 
  // the analysis or not. 
  // Returns 0 if a track is accepted or an integer different from 0 
  // to indicate the failed cut
  //
  
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
  if ( TMath::Abs(track->GetY() / track->GetX()) > fEdgeYXCutNoise )
    if ( TMath::Abs(track->GetTgl()) < fEdgeThetaCutNoise ) return 1;
  if (track->GetNumberOfClusters() < fMinClusters) return 2;
  Float_t ratio = track->GetNumberOfClusters() / (track->GetNFoundable() + 1.);
  if (ratio < fMinRatio) return 3;
  //   Float_t mpt = track->Get1Pt();       // Get1Pt() doesn't exist any more
  Float_t mpt = track->GetSigned1Pt();
  if (TMath::Abs(mpt) > fMax1pt) return 4;
  
  return 0;
}

Int_t AliTPCcalibTracksCuts::AcceptTrack(const AliVTrack * track) const {
  //
  // Function, that decides wheather a given track is accepted for 
  // the analysis or not. 
  // Returns 0 if a track is accepted or an integer different from 0 
  // to indicate the failed cut
  //
  
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
    AliExternalTrackParam trkprm;
    track->GetTrackParam(trkprm);
  if ( TMath::Abs(trkprm.GetY() / trkprm.GetX()) > fEdgeYXCutNoise )
    if ( TMath::Abs(trkprm.GetTgl()) < fEdgeThetaCutNoise ) return 1;
  if (track->GetTPCNcls() < fMinClusters) return 2;
  Float_t ratio = track->GetTPCNcls() / (track->GetTPCNclsF() + 1.);
  if (ratio < fMinRatio) return 3;
  //   Float_t mpt = track->Get1Pt();       // Get1Pt() doesn't exist any more
  Float_t mpt = trkprm.GetSigned1Pt();
  if (TMath::Abs(mpt) > fMax1pt) return 4;
  
  return 0;
}

void AliTPCcalibTracksCuts::Print(Option_t*) const {
  //
  // Print the cut contents
  //
  cout << "<AliTPCcalibTracksCuts>: The following cuts are specified: " << endl;
  cout << "fMinClusters: " << fMinClusters << endl;
  cout << "fMinRatio: " << fMinRatio << endl;
  cout << "fMax1pt: " << fMax1pt << endl;
  cout << "fEdgeYXCutNoise: " << fEdgeYXCutNoise << endl;
  cout << "fEdgeThetaCutNoise: " << fEdgeThetaCutNoise << endl;
}  // Prints out the specified cuts
