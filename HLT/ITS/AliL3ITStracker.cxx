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

//-------------------------------------------------------------------------
//               Implementation of the HLT ITS tracker class
//    It reads AliITSclusterV2 clusters and HLT ESD tracks and creates 
//    AliITStrackV2 tracks. For details, see also RunHLTITS.C macro.
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

#include "AliESD.h"
#include "AliL3ITStrack.h"
#include "AliL3ITStracker.h"

ClassImp(AliL3ITStracker)

static Int_t CorrectForDeadZoneMaterial(AliITStrackV2 *t) {
  //--------------------------------------------------------------------
  // Correction for the material between the TPC and the ITS
  // (should it belong to the TPC code ?)
  //--------------------------------------------------------------------
  Double_t riw=80., diw=0.0053, x0iw=30; // TPC inner wall ? 
  Double_t rcd=61., dcd=0.0053, x0cd=30; // TPC "central drum" ?
  Double_t yr=12.8, dr=0.03; // rods ?
  Double_t zm=0.2, dm=0.40;  // membrane
  //Double_t rr=52., dr=0.19, x0r=24., yyr=7.77; //rails
  Double_t rs=50., ds=0.001; // something belonging to the ITS (screen ?)

  if (t->GetX() > riw) {
     if (!t->PropagateTo(riw,diw,x0iw)) return 1;
     if (TMath::Abs(t->GetY())>yr) t->CorrectForMaterial(dr);
     if (TMath::Abs(t->GetZ())<zm) t->CorrectForMaterial(dm);
     if (!t->PropagateTo(rcd,dcd,x0cd)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,dr,x0r); 
     if (!t->PropagateTo(rs,ds)) return 1;
  } else if (t->GetX() < rs) {
     if (!t->PropagateTo(rs,-ds)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,-dr,x0r); 
     if (!t->PropagateTo(rcd,-dcd,x0cd)) return 1;
     if (!t->PropagateTo(riw+0.001,-diw,x0iw)) return 1;
  } else {
    //    ::Error("CorrectForDeadZoneMaterial","track is already in the dead zone !");
    return 1;
  }
  
  return 0;
}

Int_t AliL3ITStracker::Clusters2Tracks(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions reconstructs HLT ITS tracks
  //--------------------------------------------------------------------
  TObjArray itsTracks(15000);

  {/* Read HLT ESD tracks */
    Int_t nentr;
    nentr=event->GetNumberOfHLTHoughTracks();
    Info("Clusters2Tracks", "Number of ESD HLT tracks: %d\n", nentr);
    while (nentr--) {

      AliESDHLTtrack *esd=event->GetHLTHoughTrack(nentr);
      if (esd->GetWeight() > 500) continue;

      AliL3ITStrack *t=0;
      try {
        t=new AliL3ITStrack(*esd,GetZ());
      } catch (const Char_t *msg) {
        Warning("Clusters2Tracks",msg);
        delete t;
        continue;
      }
      if (TMath::Abs(t->GetD())>5) {
	delete t;
	continue;
      }

      if (CorrectForDeadZoneMaterial(t)!=0) {
	Warning("Clusters2Tracks",
		"failed to correct for the material in the dead zone !\n");
	delete t;
	continue;
      }
      itsTracks.AddLast(t);
    }
  } /* End Read HLT ESD tracks */

  itsTracks.Sort();
  Int_t nentr=itsTracks.GetEntriesFast();
  Info("Clusters2Tracks", "Number of Selected for tracking HLT ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
       //       Info("Clusters2Tracks"," %d ",i);
       AliL3ITStrack *t=(AliL3ITStrack*)itsTracks.UncheckedAt(i);
       if (t==0) continue;              //this track has been already tracked
       Int_t tpcLabel=t->GetLabel(); //save the TPC track label
       //       Info("Clusters2Tracks","Pt:%f",1/t->Get1Pt());
       ResetTrackToFollow(*t);
       ResetBestTrack();

       for (FollowProlongation(); fI<kMaxLayer; fI++) {
	 while (TakeNextProlongation()) FollowProlongation();
       }

       if (fBestTrack.GetNumberOfClusters() == 0) continue;
       
       if (fConstraint[fPass]) {
	 ResetTrackToFollow(*t);
	 if (!RefitAt(3.7, &fTrackToFollow, &fBestTrack)) continue;
	 ResetBestTrack();
       }
       
       if (!fBestTrack.PropagateTo(3.,0.0028,65.19)) continue;
       if (!fBestTrack.PropagateToVertex()) continue;
       fBestTrack.SetLabel(tpcLabel);
       fBestTrack.CookdEdx();
       CookLabel(&fBestTrack,0.); //For comparison only
       // Specific to the AliL3ITStracker
       //       fBestTrack.UpdateESDtrack(AliESDtrack::kITSin);
       t->GetESDHLTtrack()->UpdateTrackParams(&fBestTrack);
       //
       UseClusters(&fBestTrack);
       delete itsTracks.RemoveAt(i);
       ntrk++;
     }
  }

  itsTracks.Delete();

  Info("Clusters2Tracks","Number of prolonged tracks: %d\n",ntrk);

  return 0;
}

Int_t AliL3ITStracker::PropagateBack(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions propagates reconstructed ITS tracks back
  //--------------------------------------------------------------------
  Int_t nentr=event->GetNumberOfHLTHoughTracks();
  Info("PropagateBack", "Number of HLT ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
     AliESDHLTtrack *esd=event->GetHLTHoughTrack(i);
     if (esd->GetWeight() > 500) continue;

     AliL3ITStrack *t=0;
     try {
        t=new AliL3ITStrack(*esd,GetZ());
     } catch (const Char_t *msg) {
        Warning("PropagateBack",msg);
        delete t;
        continue;
     }

     ResetTrackToFollow(*t);

     // propagete to vertex [SR, GSI 17.02.2003]
     // Start Time measurement [SR, GSI 17.02.2003], corrected by I.Belikov
     if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {
       if (fTrackToFollow.PropagateToVertex()) {
          fTrackToFollow.StartTimeIntegral();
       }
       fTrackToFollow.PropagateTo(3.,-0.0028,65.19);
     }

     fTrackToFollow.ResetCovariance(); fTrackToFollow.ResetClusters();
     if (RefitAt(49.,&fTrackToFollow,t)) {
        if (CorrectForDeadZoneMaterial(&fTrackToFollow)!=0) {
          Warning("PropagateBack",
                  "failed to correct for the material in the dead zone !\n");
          delete t;
          continue;
        }
        fTrackToFollow.SetLabel(t->GetLabel());
        //fTrackToFollow.CookdEdx();
        CookLabel(&fTrackToFollow,0.); //For comparison only
	//	cout<<" Backtrack "<<fTrackToFollow.GetLabel()<<" "<<fTrackToFollow.GetNumberOfClusters()<<" "<<fTrackToFollow.fFakeRatio<<endl;
	fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSout);
        //UseClusters(&fTrackToFollow);
        ntrk++;
     }
     delete t;
  }

  Info("PropagateBack","Number of back propagated HLT ITS tracks: %d\n",ntrk);

  return 0;
}

Int_t AliL3ITStracker::RefitInward(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions refits ITS tracks using the 
  // "inward propagated" TPC tracks
  //--------------------------------------------------------------------

  Int_t nentr=event->GetNumberOfHLTHoughTracks();
  Info("RefitInward", "The method is not yet implemented! %d",nentr);
  return 0;
}
