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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     Class to make and monitor alignment of the TPC                        //
// 
//         


#include "AliTracker.h"
#include "AliTPCcalibAlignment.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClusterParam.h"
#include "AliESDtrack.h"
#include "TTreeStream.h"
#include "AliTPCTracklet.h"

#include <iostream>
using namespace std;

/* TEST

  gSystem->Load("$ALICE_ROOT/TPC/TPCcalib/libTPCcalib.so");

  .L AliXRDPROOFtoolkit.cxx+
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("listcosmic.txt","esdTree",0,100,0)
//TChain * chain = tool.MakeChain("listpp.txt","esdTree",0,100,0)

  AliTPCcalibTracks::AddInfo(chain,"TPCClusterParam.root");
  AliTPCcalibTracks::AddCuts(chain,"LOWFLUX");

//  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 0., 10., 2); // NO FIELD
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 2); // FIELD
//  field->SetL3ConstField(0); //Using const. field in the barrel
  AliTracker::SetFieldMap(field,kFALSE);

   chain->SetBranchStatus("*",1);
  chain->Process("$ALICE_ROOT/TPC/TPCcalib/AliTPCSelectorTracks.cxx+"); 

  TFile f("AliTPCcalibAlignmentDebug.root")  
  TTree *t=f.Get("Alignment")

  t->Draw("tracklet1.fInner.GetY()-trackletR1.fInner.GetY()")
  t->Draw("tracklet1.fInner.GetSnp()-tracklet1.fOuter.GetSnp()")
  t->Draw("common1.GetSnp()-common2.GetSnp()")
  t->Draw("common1.GetY()-common2.GetY()")
  t->Draw("common1.GetZ()-common2.GetZ()")

  t->Draw("track1.GetY()-track2.GetY()","abs(track1.GetY()-track2.GetY())<0.2")
  t->Draw("track1.GetZ()-track2.GetZ()","abs(track1.GetZ()-track2.GetZ())<0.2")
  t->Draw("track1.GetY()-track2.GetY():track1.GetZ()-track2.GetZ()","abs(track1.GetY()-track2.GetY())<0.2&&abs(track1.GetZ()-track2.GetZ())<0.2")


 */

ClassImp(AliTPCcalibAlignment)

AliTPCcalibAlignment::AliTPCcalibAlignment():
  fDebugStream(0)
{
  //
  // Constructor
  //
   fDebugStream=new TTreeSRedirector("AliTPCcalibAlignmentDebug.root");
}

AliTPCcalibAlignment::~AliTPCcalibAlignment() {
  //
  // Destructor
  //
  delete fDebugStream;
}

void AliTPCcalibAlignment::Process(AliTPCseed *track) {
  //
  // Track processing
  //
  TObjArray tracklets=
    AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kKalman,kFALSE,20,2);
  TObjArray trackletsL=
    AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kLinear,kFALSE,20,2);
  TObjArray trackletsQ=
    AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kQuadratic,kFALSE,20,2);
  TObjArray trackletsR=
    AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kRiemann,kFALSE,20,2);
  tracklets.SetOwner();
  if (tracklets.GetEntries()==2) {
    AliTPCTracklet* t1=static_cast<AliTPCTracklet*>(tracklets[0]);
    AliTPCTracklet* t2=static_cast<AliTPCTracklet*>(tracklets[1]);
    AliExternalTrackParam *common1=0,*common2=0;
    if (AliTPCTracklet::PropagateToMeanX(*t1,*t2,common1,common2))
      (*fDebugStream)<<"Alignment"
		     <<"common1.="<<common1
		     <<"common2.="<<common2
		     <<"tracklet1.="<<tracklets[0]
		     <<"tracklet2.="<<tracklets[1]
		     <<"trackletL1.="<<trackletsL[0]
		     <<"trackletL2.="<<trackletsL[1]
		     <<"trackletQ1.="<<trackletsQ[0]
		     <<"trackletQ2.="<<trackletsQ[1]
		     <<"trackletR1.="<<trackletsR[0]
		     <<"trackletR2.="<<trackletsR[1]
		     <<"\n";
    delete common1;
    delete common2;
 }
}


