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

/* $Id$ */

//_________________________________________________________________________
//  Algorithm class for the reconstruction: clusterizer
//                                          track segment maker
//                                          particle identifier   
//                  
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TClonesArray.h"

// --- Standard library ---

#include <iomanip>

// --- AliRoot header files ---

#include "AliPHOSReconstructioner.h"
#include "AliPHOSClusterizer.h"

ClassImp(AliPHOSReconstructioner)

//____________________________________________________________________________
AliPHOSReconstructioner::AliPHOSReconstructioner(AliPHOSClusterizer * Clusterizer, 
						 AliPHOSTrackSegmentMaker * Tracker,
						 AliPHOSPID * Pid)
{
  // ctor
  
  fClusterizer        = Clusterizer ;
  fTrackSegmentMaker  = Tracker ;
  fPID                = Pid ; 
  fDebugReconstruction = kFALSE ;
} 


//____________________________________________________________________________
 void AliPHOSReconstructioner::Init(AliPHOSClusterizer * Clusterizer, 
						 AliPHOSTrackSegmentMaker * Tracker,
						 AliPHOSPID * Pid)
{
  fClusterizer        = Clusterizer ;
  fTrackSegmentMaker  = Tracker ;
  fPID                = Pid ; 
  fDebugReconstruction = kFALSE ;
} 

//____________________________________________________________________________
 void AliPHOSReconstructioner::Make(DigitsList * dl, RecPointsList * emccl, RecPointsList * ppsdl, 
				     TrackSegmentsList * trsl, RecParticlesList * rpl)
{
  // Launches the Reconstruction process in the sequence: Make the reconstructed poins (clusterize)
  //                                                      Make the track segments 
  //                                                      Make the reconstructed particles

  Int_t index ; 
  // Digit Debuging


  if  (fDebugReconstruction) 
    {
      cout << ">>>>>>>>>>>>>>>>>>>>>> DebugReconstruction  <<<<<<<<<<<<<<<<<<<<<<<<<<"  << endl ;
      cout << "DebugReconstruction>>> Digit list entries is " <<    dl->GetEntries() << endl ;
      AliPHOSDigit * digit;
      Bool_t calorimeter ;
      Float_t factor;
      cout << "DebugReconstruction>>>    Vol Id " << 
	  " Ene(MeV, KeV) "             <<                         
	  " Index "                    << 
	  " Nprim "                     << 
	  " Prim1 "                      << 
	  " Prim2 "                      << 
	  " Prim3 "                      <<  endl;  

      for (index = 0 ; index < dl->GetEntries() ; index++) {
	digit = (AliPHOSDigit * )  dl->At(index) ;
	calorimeter = fClusterizer->IsInEmc(digit);
	if (calorimeter) factor =1000. ; else factor=1000000.;
	cout << "DebugReconstruction>>>  " << 
	  setw(8)  <<  digit->GetId() << " "  <<
	  setw(3)  <<  calorimeter <<  
	  setw(10) <<  factor*fClusterizer->Calibrate(digit->GetAmp()) <<       "  "  <<                   
	  setw(6)  <<  digit->GetIndexInList() << "  "  << 
	  setw(5)  <<  digit->GetNprimary() <<"  "  << 
	  setw(5)  <<  digit->GetPrimary(1) <<"  "  << 
	  setw(5)  <<  digit->GetPrimary(2) <<"  "  << 
	  setw(5)  <<  digit->GetPrimary(3) << endl;  	 
      }

    }



  // Making Clusters
  if  (fDebugReconstruction)  cout << "DebugReconstruction>>> Start making reconstructed points (clusterizing)" << endl;
  fClusterizer->MakeClusters(dl, emccl, ppsdl);

  // mark the position of the RecPoints in the array
  AliPHOSEmcRecPoint * emcrp ; 
  for (index = 0 ; index < emccl->GetEntries() ; index++) {
    emcrp = (AliPHOSEmcRecPoint * )emccl->At(index) ; 
    emcrp->SetIndexInList(index) ; 
  }
  AliPHOSPpsdRecPoint * ppsdrp ; 
  for (index = 0 ; index < ppsdl->GetEntries() ; index++) {
    ppsdrp = (AliPHOSPpsdRecPoint * )ppsdl->At(index) ; 
    ppsdrp->SetIndexInList(index) ; 
  }

  if  (fDebugReconstruction) 
    {
      cout << "DebugReconstruction>>> Cluster emc list entries is " <<    emccl->GetEntries() << endl ;
      AliPHOSEmcRecPoint * recpoint;
      cout << "DebugReconstruction>>> Module "  << 
	  "Ene(MeV) "             <<                         
	  "Index "                    << 
	  "Multi "                    << 
	  "   X     "                      << 
	  "   Y     "                      << 
	  "   Z    "                      << 
	  " Lambda 1   "                     <<  
	  " Lambda 2   "                     <<
  	  "MaxEnergy(MeV) "                 <<
	"Nprim "                 <<
	"Prim1 "                 <<
	"Prim2 "                 <<
	"Prim3 "                 << endl;  
      for (index = 0 ; index < emccl->GetEntries() ; index++) {
	recpoint = (AliPHOSEmcRecPoint * )emccl->At(index) ; 
	TVector3  locpos;  recpoint->GetLocalPosition(locpos);
	Float_t lambda[2]; recpoint->GetElipsAxis(lambda);
	Int_t * primaries; 
	Int_t nprimaries;
	primaries = recpoint->GetPrimaries(nprimaries);
	cout << "DebugReconstruction>>>  " << 
	  setw(2) <<recpoint->GetPHOSMod() << " "  << 
	  setw(9) << 1000.*recpoint->GetTotalEnergy() <<       " "  <<                   
	  setw(6) <<  recpoint->GetIndexInList() << " "  << 
	  setw(5) <<  recpoint->GetMultiplicity() <<" "  << 
	  setw(8) <<  locpos.X() <<" "  << 
	  setw(8) <<  locpos.Y() <<" "  << 
	  setw(8) <<  locpos.Z() << " " <<
	  setw(10) << lambda[0] << "  " <<
	  setw(10) << lambda[1] << "  " <<
	  setw(9) << 1000*recpoint->GetMaximalEnergy() << "  " << 
	  setw(9) << nprimaries << "  "  <<
	  setw(4) << primaries[0] << "  "  <<
	  setw(4) << primaries[1] << "  "  <<
	  setw(4) << primaries[2] << "  "  << endl;	   	 
      }

      cout << "DebugReconstruction>>> Cluster ppsd list entries is " <<    ppsdl->GetEntries() << endl ;
      AliPHOSPpsdRecPoint * ppsdrecpoint;
      Text_t detector[4];
      cout << "DebugReconstruction>>> Module "  << 
	  "Det     "             <<      
	  "Ene(KeV) "             <<                         
	  "Index "                    << 
	  "Multi "                    << 
	  "   X     "                      << 
	  "   Y     "                      << 
	  "   Z         "                      << 
	"Nprim "                 <<
	"Prim1 "                 <<
	"Prim2 "                 <<
	"Prim3 "                 << endl;  
      for (index = 0 ; index < ppsdl->GetEntries() ; index++) {
	ppsdrecpoint = (AliPHOSPpsdRecPoint * ) ppsdl->At(index) ; 
	TVector3  locpos; ppsdrecpoint->GetLocalPosition(locpos);
	Int_t * primaries; 
	Int_t nprimaries;
        if (ppsdrecpoint->GetUp()) detector="CPV"; else detector="PC ";
	primaries = ppsdrecpoint->GetPrimaries(nprimaries);
	cout << "DebugReconstruction>>> " << 
	  setw(4) << ppsdrecpoint->GetPHOSMod() << "  "  << 
	  setw(4)  << detector << " "      <<
	  setw(9) << 1000000.*ppsdrecpoint->GetTotalEnergy() <<       " "  <<             
	  setw(6) <<  ppsdrecpoint->GetIndexInList() << " "  << 
	  setw(5) <<  ppsdrecpoint->GetMultiplicity() <<" "  << 
	  setw(8) <<  locpos.X() <<" "  << 
	  setw(8) <<  locpos.Y() <<" "  << 
	  setw(8) <<  locpos.Z() << " " <<
	  setw(9) << nprimaries << "  "  <<
	  setw(4) << primaries[0] << "  "  <<
	  setw(4) << primaries[1] << "  "  <<
	  setw(4) << primaries[2] << "  "  << endl;	   	 
      }
    }  


  if  (fDebugReconstruction)  cout << "DebugReconstruction>>>> Start making track segments(unfolding+tracksegments)" << endl;
  fTrackSegmentMaker->MakeTrackSegments(dl, emccl, ppsdl, trsl) ;   

  // mark the position of the TrackSegments in the array
  AliPHOSTrackSegment * trs ; 
  for (index = 0 ; index < trsl->GetEntries() ; index++) {
    trs = (AliPHOSTrackSegment * )trsl->At(index) ; 
    trs->SetIndexInList(index) ; 
  }
  if  (fDebugReconstruction){
    cout << "DebugReconstruction>>> Track segment list entries is " <<    trsl->GetEntries() << endl ;
    cout << "DebugReconstruction>>> Module "  << 
      "Ene(KeV) "             <<                         
      "Index "                    << 
      "   X      "                      << 
      "   Y      "                      << 
      "   Z       "                      <<
      " rX        "                      << 
      " rY        "                      << 
      " rZ      "                               << 
      "Nprim "                 <<
      "Prim1 "                 <<
      "Prim2 "                 <<
      "Prim3 "                 << endl;  

    for (index = 0 ; index < trsl->GetEntries() ; index++) {
      trs = (AliPHOSTrackSegment * )trsl->At(index) ; 
      TVector3 locpos; trs->GetPosition(locpos);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = trs->GetPrimariesEmc(nprimaries);
      cout << "DebugReconstruction>>> " << 
	setw(4) << trs->GetPHOSMod() << "  "  << 
    	setw(9) << 1000.*trs->GetEnergy() <<       " "  <<             
	setw(3) <<  trs->GetIndexInList() << " "  <<  
	setw(9) <<  locpos.X() <<" "  << 
	setw(9) <<  locpos.Y() <<" "  << 
	setw(9) <<  locpos.Z() << " " <<
	setw(10) <<  (trs->GetMomentumDirection()).X() << " " <<
	setw(10) <<  (trs->GetMomentumDirection()).Y() << " " <<
	setw(10) <<  (trs->GetMomentumDirection()).Z() << " " <<
	setw(4) << nprimaries << "  "  <<
	setw(4) << primaries[0] << "  "  <<
	setw(4) << primaries[1] << "  "  <<
	setw(4) << primaries[2] << "  "  << endl;
      }

  }
  if  (fDebugReconstruction)  cout << "DebugReconstruction>>>> Start making particles" << endl;
  
  fPID->MakeParticles(trsl, rpl) ; 
  
  // mark the position of the RecParticles in the array
  AliPHOSRecParticle * rp ; 
  for (index = 0 ; index < rpl->GetEntries() ; index++) {
    rp = (AliPHOSRecParticle * )rpl->At(index) ; 
    rp->SetIndexInList(index) ; 
  }
}
