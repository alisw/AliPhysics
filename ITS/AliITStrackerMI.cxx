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
//               Implementation of the ITS tracker class
//    It reads AliITSclusterV2 clusters and creates AliITStrackV2 tracks
//                   and fills with them the ESD
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//                  Marian Ivanov, CERN, Marian.Ivanov@cern.ch 
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//     
//-------------------------------------------------------------------------
#include "AliITSrecoV2.h"
#include <TTree.h>
#include "AliITSgeom.h"
#include "AliTPCtrack.h"
#include "AliESD.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerMI.h"
#include "TMatrixD.h"
#include "AliHelix.h"



ClassImp(AliITStrackerMI)
ClassImp(AliITSRecV0Info)



AliITStrackerMI::AliITSlayer AliITStrackerMI::fgLayers[kMaxLayer]; // ITS layers

AliITStrackerMI::AliITStrackerMI(const AliITSgeom *geom) : AliTracker() {
  //--------------------------------------------------------------------
  //This is the AliITStrackerMI constructor
  //--------------------------------------------------------------------
  AliITSgeom *g=(AliITSgeom*)geom;

  Float_t x,y,z;
  Int_t i;
  for (i=1; i<kMaxLayer+1; i++) {
    Int_t nlad=g->GetNladders(i);
    Int_t ndet=g->GetNdetectors(i);

    g->GetTrans(i,1,1,x,y,z); 
    Double_t r=TMath::Sqrt(x*x + y*y);
    Double_t poff=TMath::ATan2(y,x);
    Double_t zoff=z;

    g->GetTrans(i,1,2,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    g->GetTrans(i,2,1,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    g->GetTrans(i,2,2,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    r*=0.25;

    new (fgLayers+i-1) AliITSlayer(r,poff,zoff,nlad,ndet);

    for (Int_t j=1; j<nlad+1; j++) {
      for (Int_t k=1; k<ndet+1; k++) { //Fill this layer with detectors
        Float_t x,y,zshift; g->GetTrans(i,j,k,x,y,zshift); 
        Double_t rot[9]; g->GetRotMatrix(i,j,k,rot);

        Double_t phi=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
        phi+=TMath::Pi()/2;
        if (i==1) phi+=TMath::Pi();
        Double_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);
        Double_t r=x*cp+y*sp;

        AliITSdetector &det=fgLayers[i-1].GetDetector((j-1)*ndet + k-1); 
        new(&det) AliITSdetector(r,phi); 
      } 
    }  

  }

  fI=kMaxLayer;

  fPass=0;
  fConstraint[0]=1; fConstraint[1]=0;

  Double_t xyz[]={kXV,kYV,kZV}, ers[]={kSigmaXV,kSigmaYV,kSigmaZV}; 
  SetVertex(xyz,ers);

  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=kLayersNotToSkip[i];
  fLastLayerToTrackTo=kLastLayerToTrackTo;

}

void AliITStrackerMI::SetLayersNotToSkip(Int_t *l) {
  //--------------------------------------------------------------------
  //This function set masks of the layers which must be not skipped
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=l[i];
}

Int_t AliITStrackerMI::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads ITS clusters
  //--------------------------------------------------------------------
  TBranch *branch=cTree->GetBranch("Clusters");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch !\n");
    return 1;
  }

  TClonesArray dummy("AliITSclusterV2",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  Int_t j=0;
  Int_t detector=0;
  for (Int_t i=0; i<kMaxLayer; i++) {
    Int_t ndet=fgLayers[i].GetNdetectors();
    Int_t jmax = j + fgLayers[i].GetNladders()*ndet;
    for (; j<jmax; j++) {           
      if (!cTree->GetEvent(j)) continue;
      Int_t ncl=clusters->GetEntriesFast();
      SignDeltas(clusters,GetZ());
      while (ncl--) {
        AliITSclusterV2 *c=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
	detector = c->GetDetectorIndex();
        fgLayers[i].InsertCluster(new AliITSclusterV2(*c));
      }
      clusters->Delete();
      //add dead zone virtual "cluster"      
      if (i<2){
	for (Float_t ydead = 0; ydead < 1.31 ; ydead+=(i+1.)*0.018){     
	  Int_t lab[4] = {0,0,0,detector};
	  Int_t info[3] = {0,0,0};
	  Float_t hit[5]={0,0,0.004/12.,0.001/12.,0};
	  if (i==0) hit[0] =ydead-0.4;
	  if (i==1) hit[0]=ydead-3.75; 
	  hit[1] =-0.04;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  hit[1]=-7.05;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  hit[1]=-7.15;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  hit[1] =0.06;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  hit[1]=7.05;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  hit[1]=7.25;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));       
	}
      }
      
    }
    //
    fgLayers[i].ResetRoad(); //road defined by the cluster density
    fgLayers[i].SortClusters();
  }

  return 0;
}

void AliITStrackerMI::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer; i++) fgLayers[i].ResetClusters();
}

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
  ::Error("CorrectForDeadZoneMaterial","track is already in the dead zone !");
    return 1;
  }
  
  return 0;
}

Int_t AliITStrackerMI::Clusters2Tracks(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------
  TObjArray itsTracks(15000);

  {/* Read ESD tracks */
    Int_t nentr=event->GetNumberOfTracks();
    Info("Clusters2Tracks", "Number of ESD tracks: %d\n", nentr);
    while (nentr--) {
      AliESDtrack *esd=event->GetTrack(nentr);

      if ((esd->GetStatus()&AliESDtrack::kTPCin)==0) continue;
      if (esd->GetStatus()&AliESDtrack::kTPCout) continue;
      if (esd->GetStatus()&AliESDtrack::kITSin) continue;

      AliITStrackV2 *t=0;
      try {
        t=new AliITStrackV2(*esd);
      } catch (const Char_t *msg) {
        Warning("Clusters2Tracks",msg);
        delete t;
        continue;
      }
      t->fD[0] = t->GetD(GetX(),GetY());
      t->fD[1] = t->GetZat(GetX())-GetZ(); 
      Double_t vdist = TMath::Sqrt(t->fD[0]*t->fD[0]+t->fD[1]*t->fD[1]);
      if (t->GetMass()<0.13) t->SetMass(0.13957); // MI look to the esd - mass hypothesys  !!!!!!!!!!!
      // write expected q
      t->fExpQ = TMath::Max(0.8*t->fESDtrack->GetTPCsignal(),30.);

      if (TMath::Abs(t->fD[0])>10) {
	delete t;
	continue;
      }

      if (TMath::Abs(vdist)>20) {
	delete t;
	continue;
      }
      if (TMath::Abs(1/t->Get1Pt())<0.120) {
	delete t;
	continue;
      }

      if (CorrectForDeadZoneMaterial(t)!=0) {
         Warning("Clusters2Tracks",
                 "failed to correct for the material in the dead zone !\n");
         delete t;
         continue;
      }
      t->fReconstructed = kFALSE;
      itsTracks.AddLast(t);
    }
  } /* End Read ESD tracks */

  itsTracks.Sort();
  Int_t nentr=itsTracks.GetEntriesFast();
  fTrackHypothesys.Expand(nentr);
  MakeCoeficients(nentr);
  Int_t ntrk=0;
  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
//       cerr<<fPass<<"    "<<i<<'\r';
       fCurrentEsdTrack = i;
       AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
       if (t==0) continue;              //this track has been already tracked
       if (t->fReconstructed&&(t->fNUsed<1.5)) continue;  //this track was  already  "succesfully" reconstructed
       if ( (TMath::Abs(t->GetD(GetX(),GetY()))  >3.) && fConstraint[fPass]) continue;
       if ( (TMath::Abs(t->GetZat(GetX())-GetZ())>3.) && fConstraint[fPass]) continue;

       Int_t tpcLabel=t->GetLabel(); //save the TPC track label       
       fI = 6;
       ResetTrackToFollow(*t);
       ResetBestTrack();
      
       FollowProlongationTree(t,i);


       SortTrackHypothesys(fCurrentEsdTrack,20,0);  //MI change
       //
       AliITStrackV2 * besttrack = GetBestHypothesys(fCurrentEsdTrack,t,15);
       if (!besttrack) continue;
       besttrack->SetLabel(tpcLabel);
       //       besttrack->CookdEdx();
       CookdEdx(besttrack);
       besttrack->fFakeRatio=1.;
       CookLabel(besttrack,0.); //For comparison only
       //       besttrack->UpdateESDtrack(AliESDtrack::kITSin);
       //
       
       UpdateESDtrack(besttrack,AliESDtrack::kITSin);
       
       if ( besttrack->GetNumberOfClusters()<5 && fConstraint[fPass]) {	 
	 continue;
       }
       if (besttrack->fChi2MIP[0]+besttrack->fNUsed>3.5) continue;
       if ( (TMath::Abs(besttrack->fD[0]*besttrack->fD[0]+besttrack->fD[1]*besttrack->fD[1])>0.1) && fConstraint[fPass])  continue;	 
       //delete itsTracks.RemoveAt(i);
       t->fReconstructed = kTRUE;
       ntrk++;                     
     }
     GetBestHypothesysMIP(itsTracks); 
  }

  //GetBestHypothesysMIP(itsTracks);
  //FindV0(event);

  itsTracks.Delete();
  //
  Int_t entries = fTrackHypothesys.GetEntriesFast();
  for (Int_t ientry=0;ientry<entries;ientry++){
    TObjArray * array =(TObjArray*)fTrackHypothesys.UncheckedAt(ientry);
    if (array) array->Delete();
    delete fTrackHypothesys.RemoveAt(ientry); 
  }

  fTrackHypothesys.Delete();
  Info("Clusters2Tracks","Number of prolonged tracks: %d\n",ntrk);

  return 0;
}



Int_t AliITStrackerMI::Clusters2Tracks(TTree *tpcTree, TTree *itsTree) {
  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------
  Int_t nentr=0; TObjArray itsTracks(15000);

   Warning("Clusters2Tracks(TTree *, TTree *)",
      "Will be removed soon !   Use Clusters2Tracks(AliESD *) instead.");

  {/* Read TPC tracks */ 
    AliTPCtrack *itrack=new AliTPCtrack; 
    TBranch *branch=tpcTree->GetBranch("tracks");
    if (!branch) {
       Error("Clusters2Tracks","Can't get the branch !");
       return 1;
    }
    tpcTree->SetBranchAddress("tracks",&itrack);
    nentr=(Int_t)tpcTree->GetEntries();

    Info("Clusters2Tracks","Number of TPC tracks: %d\n",nentr);

    for (Int_t i=0; i<nentr; i++) {
       tpcTree->GetEvent(i);
       AliITStrackV2 *t=0;
       try {
           t=new AliITStrackV2(*itrack);
       } catch (const Char_t *msg) {
           Warning("Clusters2Tracks",msg);
           delete t;
           continue;
       }
       if (TMath::Abs(t->GetD())>4) {
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
    delete itrack;
  }
  itsTracks.Sort();
  nentr=itsTracks.GetEntriesFast();


  AliITStrackV2 *otrack=&fBestTrack;
  TBranch *branch=itsTree->GetBranch("tracks");
  if (!branch) itsTree->Branch("tracks","AliITStrackV2",&otrack,32000,3);
  else branch->SetAddress(&otrack);

  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
       AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
       if (t==0) continue;           //this track has been already tracked
       Int_t tpcLabel=t->GetLabel(); //save the TPC track label

       ResetTrackToFollow(*t);
       ResetBestTrack();
       /*
       for (FollowProlongation(); fI<kMaxLayer; fI++) {
          while (TakeNextProlongation()) FollowProlongation();
       }
       */
       FollowProlongationTree(t,i);
       if (fBestTrack.GetNumberOfClusters() == 0) continue;

       if (fConstraint[fPass]) {
          ResetTrackToFollow(*t);
	  if (!RefitAt(3.7, &fTrackToFollow, &fBestTrack)) continue;
          ResetBestTrack();
       }

       fBestTrack.SetLabel(tpcLabel);
       //fBestTrack.CookdEdx();
       CookdEdx(&fBestTrack);

       CookLabel(&fBestTrack,0.); //For comparison only
       itsTree->Fill();
       //UseClusters(&fBestTrack);
       delete itsTracks.RemoveAt(i);
     }
  }

  nentr=(Int_t)itsTree->GetEntries();
  Info("Clusters2Tracks","Number of prolonged tracks: %d\n",nentr);

  itsTracks.Delete();

  return 0;
}

Int_t AliITStrackerMI::PropagateBack(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions propagates reconstructed ITS tracks back
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  Int_t nentr=event->GetNumberOfTracks();
  Info("PropagateBack", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
     AliESDtrack *esd=event->GetTrack(i);

     if ((esd->GetStatus()&AliESDtrack::kITSin)==0) continue;
     if (esd->GetStatus()&AliESDtrack::kITSout) continue;

     AliITStrackV2 *t=0;
     try {
        t=new AliITStrackV2(*esd);
     } catch (const Char_t *msg) {
        Warning("PropagateBack",msg);
        delete t;
        continue;
     }
     t->fExpQ = TMath::Max(0.8*t->fESDtrack->GetTPCsignal(),30.);

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
        fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSout);
        //UseClusters(&fTrackToFollow);
        ntrk++;
     }
     delete t;
  }

  Info("PropagateBack","Number of back propagated ITS tracks: %d\n",ntrk);

  return 0;
}

Int_t AliITStrackerMI::RefitInward(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions refits ITS tracks using the 
  // "inward propagated" TPC tracks
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  Int_t nentr=event->GetNumberOfTracks();
  Info("RefitInward", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);

    if ((esd->GetStatus()&AliESDtrack::kITSout) == 0) continue;
    if (esd->GetStatus()&AliESDtrack::kITSrefit) continue;
    if (esd->GetStatus()&AliESDtrack::kTPCout)
      if ((esd->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;

    AliITStrackV2 *t=0;
    try {
        t=new AliITStrackV2(*esd);
    } catch (const Char_t *msg) {
        Warning("RefitInward",msg);
        delete t;
        continue;
    }
    t->fExpQ = TMath::Max(0.8*t->fESDtrack->GetTPCsignal(),30.);
    if (CorrectForDeadZoneMaterial(t)!=0) {
       Warning("RefitInward",
               "failed to correct for the material in the dead zone !\n");
       delete t;
       continue;
    }

    ResetTrackToFollow(*t);
    fTrackToFollow.ResetClusters();

    if ((esd->GetStatus()&AliESDtrack::kTPCin)==0)
      fTrackToFollow.ResetCovariance();

    //Refitting...
    if (RefitAt(3.7, &fTrackToFollow, t)) {
       fTrackToFollow.SetLabel(t->GetLabel());
       //       fTrackToFollow.CookdEdx();
       CookdEdx(&fTrackToFollow);

       CookLabel(&fTrackToFollow,0.0); //For comparison only

       if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {//The beam pipe    
         Double_t a=fTrackToFollow.GetAlpha();
         Double_t cs=TMath::Cos(a),sn=TMath::Sin(a);
         Double_t xv= GetX()*cs + GetY()*sn;
         Double_t yv=-GetX()*sn + GetY()*cs;
         
         Double_t c=fTrackToFollow.GetC(), snp=fTrackToFollow.GetSnp();
         Double_t x=fTrackToFollow.GetX(), y=fTrackToFollow.GetY();
         Double_t tgfv=-(c*(x-xv)-snp)/(c*(y-yv) + TMath::Sqrt(1.-snp*snp));
         Double_t fv=TMath::ATan(tgfv);

         cs=TMath::Cos(fv); sn=TMath::Sin(fv);
         x = xv*cs + yv*sn;
         yv=-xv*sn + yv*cs; xv=x;

	 if (fTrackToFollow.Propagate(fv+a,xv)) {
            fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSrefit);
            //UseClusters(&fTrackToFollow);
            {
            AliITSclusterV2 c; c.SetY(yv); c.SetZ(GetZ());
            c.SetSigmaY2(GetSigmaY()*GetSigmaY());
            c.SetSigmaZ2(GetSigmaZ()*GetSigmaZ());
            Double_t chi2=fTrackToFollow.GetPredictedChi2(&c);
            //Double_t chi2=GetPredictedChi2MI(&fTrackToFollow,&c,fI);
            if (chi2<kMaxChi2)
	      if (fTrackToFollow.Update(&c,-chi2,0))
		//if (UpdateMI(&fTrackToFollow,&c,-chi2,0))
		fTrackToFollow.SetConstrainedESDtrack(chi2);            
            }
            ntrk++;
         }
       }
    }
    delete t;
  }

  Info("RefitInward","Number of refitted tracks: %d\n",ntrk);

  return 0;
}

AliCluster *AliITStrackerMI::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fgLayers[l].GetCluster(c);
}


void AliITStrackerMI::FollowProlongationTree(AliITStrackV2 * otrack, Int_t esdindex) 
{
  //--------------------------------------------------------------------
  // Follow prolongation tree
  //--------------------------------------------------------------------

  //setup tree of the prolongations
  //
  static AliITStrackV2 tracks[7][100];
  AliITStrackV2 *currenttrack;
  static AliITStrackV2 currenttrack1;
  static AliITStrackV2 currenttrack2;  
  static AliITStrackV2 backuptrack;
  Int_t ntracks[7];
  Int_t nindexes[7][100];
  Float_t normalizedchi2[100];
  for (Int_t ilayer=0;ilayer<6;ilayer++) ntracks[ilayer]=0;
  otrack->fNSkipped=0;
  new (&(tracks[6][0])) AliITStrackV2(*otrack);
  ntracks[6]=1;
  nindexes[6][0]=0;
  // 
  //
  // follow prolongations
  for (Int_t ilayer=5;ilayer>=0;ilayer--){
    //
    AliITSlayer &layer=fgLayers[ilayer]; 
    Double_t r=layer.GetR();
    ntracks[ilayer]=0;
    //
    //
   Int_t nskipped=0;
    Float_t nused =0;
    for (Int_t itrack =0;itrack<ntracks[ilayer+1];itrack++){
      //set current track
      if (ntracks[ilayer]>=100) break;  
      if (tracks[ilayer+1][nindexes[ilayer+1][itrack]].fNSkipped>0) nskipped++;
      if (tracks[ilayer+1][nindexes[ilayer+1][itrack]].fNUsed>2.) nused++;
      if (ntracks[ilayer]>15+ilayer){
	if (itrack>1&&tracks[ilayer+1][nindexes[ilayer+1][itrack]].fNSkipped>0 && nskipped>4+ilayer) continue;
	if (itrack>1&&tracks[ilayer+1][nindexes[ilayer+1][itrack]].fNUsed>2. && nused>3) continue;
      }

      new(&currenttrack1)  AliITStrackV2(tracks[ilayer+1][nindexes[ilayer+1][itrack]]);
      if (ilayer==3 || ilayer==1) {
	Double_t rs=0.5*(fgLayers[ilayer+1].GetR() + r);
	Double_t d=0.0034, x0=38.6;
	if (ilayer==1) {rs=9.; d=0.0097; x0=42;}
	if (!currenttrack1.PropagateTo(rs,d,x0)) {
	  continue;
	}
      }
      //
      //find intersection with layer
      Double_t x,y,z;  
      if (!currenttrack1.GetGlobalXYZat(r,x,y,z)) {
	continue;
      }
      Double_t phi=TMath::ATan2(y,x);
      Int_t idet=layer.FindDetectorIndex(phi,z);
      if (idet<0) {
	continue;
      }
      //propagate to the intersection
      const AliITSdetector &det=layer.GetDetector(idet);
      phi=det.GetPhi();
      new(&currenttrack2)  AliITStrackV2(currenttrack1);
      if (!currenttrack1.Propagate(phi,det.GetR())) {	
	continue;
      }
      currenttrack2.Propagate(phi,det.GetR());  //
      currenttrack1.SetDetectorIndex(idet);
      currenttrack2.SetDetectorIndex(idet);
      
      //
      //
      Double_t dz=7.5*TMath::Sqrt(currenttrack1.GetSigmaZ2() + 16.*kSigmaZ2[ilayer]);
      Double_t dy=7.5*TMath::Sqrt(currenttrack1.GetSigmaY2() + 16.*kSigmaY2[ilayer]);
      //
      Bool_t isBoundary=kFALSE;
      if (currenttrack1.GetY()-dy< det.GetYmin()+0.2) isBoundary = kTRUE;  
      if (currenttrack1.GetY()+dy> det.GetYmax()-0.2) isBoundary = kTRUE;
      if (currenttrack1.GetZ()-dz< det.GetZmin()+0.2) isBoundary = kTRUE;
      if (currenttrack1.GetZ()+dz> det.GetZmax()-0.2) isBoundary = kTRUE;
      
      if (isBoundary){ // track at boundary between detectors
	Float_t maxtgl = TMath::Abs(currenttrack1.GetTgl());
	if (maxtgl>1) maxtgl=1;
	dz = TMath::Sqrt(dz*dz+0.25*maxtgl*maxtgl);
	//
	Float_t maxsnp = TMath::Abs(currenttrack1.GetSnp());
	if (maxsnp>0.95) continue;
	//if (maxsnp>0.5) maxsnp=0.5;
	dy=TMath::Sqrt(dy*dy+0.25*maxsnp*maxsnp);
      }
      
      Double_t zmin=currenttrack1.GetZ() - dz; 
      Double_t zmax=currenttrack1.GetZ() + dz;
      Double_t ymin=currenttrack1.GetY() + r*phi - dy;
      Double_t ymax=currenttrack1.GetY() + r*phi + dy;
      layer.SelectClusters(zmin,zmax,ymin,ymax); 
      //
      // loop over all possible prolongations
      //
      Double_t msz=1./((currenttrack1.GetSigmaZ2() + 16.*kSigmaZ2[ilayer]));
      Double_t msy=1./((currenttrack1.GetSigmaY2() + 16.*kSigmaY2[ilayer]));
      if (fConstraint[fPass]){
	msy/=60; msz/=60.;
      }
      else{
	msy/=50; msz/=50.;
      }
      //
      const AliITSclusterV2 *c=0; Int_t ci=-1;
      Double_t chi2=12345.;
      Int_t deadzone=0;
      currenttrack = &currenttrack1;
      while ((c=layer.GetNextCluster(ci))!=0) { 
	if (ntracks[ilayer]>95) break; //space for skipped clusters  
	Bool_t change =kFALSE;  
	if (c->GetQ()==0 && (deadzone==1)) continue;
	Int_t idet=c->GetDetectorIndex();
	if (currenttrack->GetDetectorIndex()!=idet) {
	  const AliITSdetector &det=layer.GetDetector(idet);
	  Double_t y,z;
	  if (!currenttrack2.GetProlongationFast(det.GetPhi(),det.GetR(),y,z)) continue;
	  Float_t pz = (z - c->GetZ()) , py=(y - c->GetY());
	  if (pz*pz*msz+py*py*msy>1.) continue;
	  //
	  new (&backuptrack) AliITStrackV2(currenttrack2);
	  change = kTRUE;
	  currenttrack =&currenttrack2;
	  if (!currenttrack->Propagate(det.GetPhi(),det.GetR())) {
	    new (currenttrack) AliITStrackV2(backuptrack);
	    change = kFALSE;
	    continue;
	  }
	  currenttrack->SetDetectorIndex(idet);
	}
	else{
	  Float_t pz = (currenttrack->GetZ() - c->GetZ()) , py=(currenttrack->GetY() - c->GetY());
	  if (pz*pz*msz+py*py*msy>1.) continue;
	}

	chi2=GetPredictedChi2MI(currenttrack,c,ilayer); 
	if (chi2<kMaxChi2s[ilayer]){
	  if (c->GetQ()==0) deadzone=1;	    // take dead zone only once	  
	  if (ntracks[ilayer]>=100) continue;
	  AliITStrackV2 * updatetrack = new (&tracks[ilayer][ntracks[ilayer]]) AliITStrackV2(*currenttrack);
	  updatetrack->fClIndex[ilayer]=0;
	  if (change){
	    new (&currenttrack2) AliITStrackV2(backuptrack);
	  }
	  if (c->GetQ()!=0){
	    if (!UpdateMI(updatetrack,c,chi2,(ilayer<<28)+ci)) continue; 
	    updatetrack->SetSampledEdx(c->GetQ(),updatetrack->GetNumberOfClusters()-1); //b.b.
	  }
	  else {
	    updatetrack->fNDeadZone++;
	    updatetrack->fDeadZoneProbability=GetDeadZoneProbability(updatetrack->GetZ(),TMath::Sqrt(updatetrack->GetSigmaZ2()));
	  }
	  if (c->IsUsed()){
	    updatetrack->fNUsed++;
	  }
	  Double_t x0;
	  Double_t d=layer.GetThickness(updatetrack->GetY(),updatetrack->GetZ(),x0);
	  updatetrack->CorrectForMaterial(d,x0);	  
	  if (fConstraint[fPass]) {
	    updatetrack->fConstrain = fConstraint[fPass];
	    fI = ilayer;
	    Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
	    Double_t xyz[]={GetX(),GetY(),GetZ()};
	    Double_t ptfactor = 1;
	    Double_t ers[]={GetSigmaX()*ptfactor,GetSigmaY()*ptfactor,GetSigmaZ()};
	    Bool_t isPrim = kTRUE;
	    if (ilayer<4){
	      updatetrack->fD[0] = updatetrack->GetD(GetX(),GetY());
	      updatetrack->fD[1] = updatetrack->GetZat(GetX())-GetZ();
	      if ( TMath::Abs(updatetrack->fD[0]/(1.+ilayer))>0.4 ||  TMath::Abs(updatetrack->fD[1]/(1.+ilayer))>0.4) isPrim=kFALSE;
	    }
	    if (isPrim) updatetrack->Improve(d,xyz,ers);
	  } //apply vertex constrain	  	  
	  ntracks[ilayer]++;
	}  // create new hypothesy 
      } // loop over possible cluster prolongation      
      //      if (fConstraint[fPass]&&itrack<2&&currenttrack1.fNSkipped==0 && deadzone==0){	
      if (itrack<2&&currenttrack1.fNSkipped==0 && deadzone==0&&ntracks[ilayer]<100){	
	AliITStrackV2* vtrack = new (&tracks[ilayer][ntracks[ilayer]]) AliITStrackV2(currenttrack1);
	vtrack->fClIndex[ilayer]=0;
	fI = ilayer;
	Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
	Double_t xyz[]={GetX(),GetY(),GetZ()};
	Double_t ers[]={GetSigmaX(),GetSigmaY(),GetSigmaZ()};
	vtrack->Improve(d,xyz,ers);
	vtrack->fNSkipped++;
	ntracks[ilayer]++;
      }
      
    } //loop over track candidates
    //
    //
    Int_t accepted=0;
    
    Int_t golds=0;
    for (Int_t itrack=0;itrack<ntracks[ilayer];itrack++){
      normalizedchi2[itrack] = NormalizedChi2(&tracks[ilayer][itrack],ilayer); 
      if ( normalizedchi2[itrack]<3+0.5*ilayer) golds++;
      if (ilayer>4) accepted++;
      else{
	if ( fConstraint[fPass] && normalizedchi2[itrack]<kMaxNormChi2C[ilayer]+1) accepted++;
	if (!fConstraint[fPass] && normalizedchi2[itrack]<kMaxNormChi2NonC[ilayer]+1) accepted++;
      }
    }
    TMath::Sort(ntracks[ilayer],normalizedchi2,nindexes[ilayer],kFALSE);
    ntracks[ilayer] = TMath::Min(accepted,7+2*ilayer);
    if (ntracks[ilayer]<golds+2+ilayer) ntracks[ilayer]=TMath::Min(golds+2+ilayer,accepted);
    if (ntracks[ilayer]>90) ntracks[ilayer]=90; 
  } //loop over layers
  //printf("%d\t%d\t%d\t%d\t%d\t%d\n",ntracks[0],ntracks[1],ntracks[2],ntracks[3],ntracks[4],ntracks[5]);

  for (Int_t i=0;i<TMath::Min(20,ntracks[0]);i++) {
    AliITStrackV2 & track= tracks[0][nindexes[0][i]];
    if (!fConstraint[fPass]&&track.fNormChi2[0]>7.)continue;
    AddTrackHypothesys(new AliITStrackV2(track), esdindex);
  }
  for (Int_t i=0;i<TMath::Min(4,ntracks[1]);i++) {
    AliITStrackV2 & track= tracks[1][nindexes[1][i]];
    if (track.GetNumberOfClusters()<4) continue;
    if (!fConstraint[fPass]&&track.fNormChi2[1]>7.)continue;
    if (fConstraint[fPass]) track.fNSkipped+=1;
    if (!fConstraint[fPass]) {
      track.fD[0] = track.GetD(GetX(),GetY());   
      track.fNSkipped+=4./(4.+8.*TMath::Abs(track.fD[0]));
      if (track.fN+track.fNDeadZone+track.fNSkipped>6) {
	track.fNSkipped = 6-track.fN+track.fNDeadZone;
      }
    }
    AddTrackHypothesys(new AliITStrackV2(track), esdindex);
  }
  //}
  
  if (!fConstraint[fPass]){  
    for (Int_t i=0;i<TMath::Min(3,ntracks[2]);i++) {
      AliITStrackV2 & track= tracks[2][nindexes[2][i]];
      if (track.GetNumberOfClusters()<4) continue;
      if (!fConstraint[fPass]&&track.fNormChi2[2]>7.)continue;
      if (fConstraint[fPass]) track.fNSkipped+=2;      
      if (!fConstraint[fPass]){
	track.fD[0] = track.GetD(GetX(),GetY());
	track.fNSkipped+= 7./(7.+8.*TMath::Abs(track.fD[0]));
	if (track.fN+track.fNDeadZone+track.fNSkipped>6) {
	  track.fNSkipped = 6-track.fN+track.fNDeadZone;
	}
      }
      AddTrackHypothesys(new AliITStrackV2(track), esdindex);
    }
  }
}


AliITStrackerMI::AliITSlayer & AliITStrackerMI::GetLayer(Int_t layer) const
{
  //--------------------------------------------------------------------
  //
  //
  return fgLayers[layer];
}
AliITStrackerMI::AliITSlayer::AliITSlayer() {
  //--------------------------------------------------------------------
  //default AliITSlayer constructor
  //--------------------------------------------------------------------
  fN=0;
  fDetectors=0;
  fSkip = 0;
  fCurrentSlice=-1;
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) {
    fClusterWeight[i]=0;
    fClusterTracks[0][i]=-1;
    fClusterTracks[1][i]=-1;
    fClusterTracks[2][i]=-1;    
    fClusterTracks[3][i]=-1;    
  }
}

AliITStrackerMI::AliITSlayer::
AliITSlayer(Double_t r,Double_t p,Double_t z,Int_t nl,Int_t nd) {
  //--------------------------------------------------------------------
  //main AliITSlayer constructor
  //--------------------------------------------------------------------
  fR=r; fPhiOffset=p; fZOffset=z;
  fNladders=nl; fNdetectors=nd;
  fDetectors=new AliITSdetector[fNladders*fNdetectors];

  fN=0;
  fI=0;
  fSkip = 0;
  fRoad=2*fR*TMath::Sqrt(3.14/1.);//assuming that there's only one cluster
}

AliITStrackerMI::AliITSlayer::~AliITSlayer() {
  //--------------------------------------------------------------------
  // AliITSlayer destructor
  //--------------------------------------------------------------------
  delete[] fDetectors;
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) {
    fClusterWeight[i]=0;
    fClusterTracks[0][i]=-1;
    fClusterTracks[1][i]=-1;
    fClusterTracks[2][i]=-1;    
    fClusterTracks[3][i]=-1;    
  }
}

void AliITStrackerMI::AliITSlayer::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  for (Int_t i=0; i<kMaxClusterPerLayer;i++){
    fClusterWeight[i]=0;
    fClusterTracks[0][i]=-1;
    fClusterTracks[1][i]=-1;
    fClusterTracks[2][i]=-1;    
    fClusterTracks[3][i]=-1;  
  }
  
  fN=0;
  fI=0;
}

void AliITStrackerMI::AliITSlayer::ResetWeights() {
  //--------------------------------------------------------------------
  // This function reset weights of the clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) {
    fClusterWeight[i]=0;
    fClusterTracks[0][i]=-1;
    fClusterTracks[1][i]=-1;
    fClusterTracks[2][i]=-1;    
    fClusterTracks[3][i]=-1;  
  }
  for (Int_t i=0; i<fN;i++) {
    AliITSclusterV2 * cl = (AliITSclusterV2*)GetCluster(i);
    if (cl&&cl->IsUsed()) cl->Use();
  }

}

void AliITStrackerMI::AliITSlayer::ResetRoad() {
  //--------------------------------------------------------------------
  // This function calculates the road defined by the cluster density
  //--------------------------------------------------------------------
  Int_t n=0;
  for (Int_t i=0; i<fN; i++) {
     if (TMath::Abs(fClusters[i]->GetZ())<fR) n++;
  }
  //if (n>1) fRoad=2*fR*TMath::Sqrt(3.14/n);
  if (n>1) fRoad=2*fR*TMath::Sqrt(3.14/n);
}

Int_t AliITStrackerMI::AliITSlayer::InsertCluster(AliITSclusterV2 *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to this layer
  //--------------------------------------------------------------------
  if (fN==kMaxClusterPerLayer) {
    ::Error("InsertCluster","Too many clusters !\n");
    return 1;
  }
  fCurrentSlice=-1;
  if (fN==0) {fClusters[fN++]=c; return 0;}
  Int_t i=FindClusterIndex(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliITSclusterV2*));
  memmove(fY+i+1 ,fY+i,(fN-i)*sizeof(Float_t));
  memmove(fZ+i+1 ,fZ+i,(fN-i)*sizeof(Float_t));
  fClusters[i]=c; fN++;
  //
  AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
  Double_t y=fR*det.GetPhi() + c->GetY();
  if (y>2.*fR*TMath::Pi()) y -= 2.*fR*TMath::Pi();
  fY[i] = y;
  fZ[i] = c->GetZ(); 
  if (c->GetY()<det.GetYmin()) det.SetYmin(c->GetY());
  if (c->GetY()>det.GetYmax()) det.SetYmax(c->GetY());
  if (c->GetZ()<det.GetZmin()) det.SetZmin(c->GetZ());
  if (c->GetZ()>det.GetZmax()) det.SetZmax(c->GetZ());
			     
  return 0;
}

void  AliITStrackerMI::AliITSlayer::SortClusters()
{
  //
  //sort clusters
  //
  fYB[0]=10000000;
  fYB[1]=-10000000;
  for (Int_t i=0;i<fN;i++){
    if (fY[i]<fYB[0]) fYB[0]=fY[i];
    if (fY[i]>fYB[1]) fYB[1]=fY[i];
    fClusterIndex[i] = i;
  }
  //
  // fill slices
  fDy5 = (fYB[1]-fYB[0])/5.;
  fDy10 = (fYB[1]-fYB[0])/10.;
  fDy20 = (fYB[1]-fYB[0])/20.;
  for (Int_t i=0;i<6;i++)  fN5[i] =0;  
  for (Int_t i=0;i<11;i++) fN10[i]=0;  
  for (Int_t i=0;i<21;i++) fN20[i]=0;
  //  
  for (Int_t i=0;i<6;i++) {fBy5[i][0] =  fYB[0]+(i-0.75)*fDy5; fBy5[i][1] =  fYB[0]+(i+0.75)*fDy5;}
  for (Int_t i=0;i<11;i++) {fBy10[i][0] =  fYB[0]+(i-0.75)*fDy10; fBy10[i][1] =  fYB[0]+(i+0.75)*fDy10;} 
  for (Int_t i=0;i<21;i++) {fBy20[i][0] =  fYB[0]+(i-0.75)*fDy20; fBy20[i][1] =  fYB[0]+(i+0.75)*fDy20;}
  //
  //
  for (Int_t i=0;i<fN;i++) for (Int_t irot=-1;irot<=1;irot++){
    Float_t curY = fY[i]+irot*2.*TMath::Pi()*fR;
    if (curY<fYB[0]-fDy5) continue;
    if (curY>fYB[1]+fDy5) continue; 
    //
    // slice 10
    Float_t fslice = TMath::Nint(10.*(curY-fYB[0])/(fYB[1]-fYB[0]));
    Float_t ymiddle = fYB[0]+fslice*fDy10;
    for (Int_t di =-1;di<=1;di++){
      if (TMath::Abs(curY-(ymiddle+(float)di*fDy10))<0.75*fDy10){
	//
	Int_t slice = int(fslice+21.0001)-21+di;
	if (slice<0) continue;
	if (slice>10) continue;
	fClusters10[slice][fN10[slice]] = fClusters[i];
	fY10[slice][fN10[slice]] = curY;
	fZ10[slice][fN10[slice]] = fZ[i];
	fClusterIndex10[slice][fN10[slice]]=i;
	fN10[slice]++;
      }
    }
    //
    // slice 5
    fslice = TMath::Nint(5.*(curY-fYB[0])/(fYB[1]-fYB[0]));
    ymiddle = fYB[0]+fslice*fDy5;
    for (Int_t di =-1;di<=1;di++){
      if (TMath::Abs(curY-(ymiddle+(float)di*fDy5))<0.75*fDy5){
	//
	Int_t slice = int(fslice+21.0001)-21+di;
	if (slice<0) continue;
	if (slice>5) continue;
	fClusters5[slice][fN5[slice]] = fClusters[i];
	fY5[slice][fN5[slice]] = curY;
	fZ5[slice][fN5[slice]] = fZ[i];
	fClusterIndex5[slice][fN5[slice]]=i;
	fN5[slice]++;
      }
    }
    //
    // slice 20
    fslice = TMath::Nint(20.*(curY-fYB[0])/(fYB[1]-fYB[0]));
    ymiddle = fYB[0]+fslice*fDy20;
    for (Int_t di =-1;di<=1;di++){
      if (TMath::Abs(curY-(ymiddle+(float)di*fDy20))<0.75*fDy20){
	//
	Int_t slice = int(fslice+21.0001)-21+di;
	if (slice<0) continue;
	if (slice>20) continue;
	fClusters20[slice][fN20[slice]] = fClusters[i];
	fY20[slice][fN20[slice]] = curY;
	fZ20[slice][fN20[slice]] = fZ[i];
	fClusterIndex20[slice][fN20[slice]]=i;
	fN20[slice]++;
      }
    }
  }
}




Int_t AliITStrackerMI::AliITSlayer::FindClusterIndex(Float_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  Int_t ncl=0;
  const Float_t *zcl;  
  if (fCurrentSlice<0) {
    ncl = fN;
    zcl   = fZ;
  }
  else{
    ncl   = fNcs;
    zcl   = fZcs;;
  }
  
  if (ncl==0) return 0;
  Int_t b=0, e=ncl-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    //    if (z > fClusters[m]->GetZ()) b=m+1;
    if (z > zcl[m]) b=m+1;
    else e=m; 
  }
  return m;
}
/*
void AliITStrackerMI::AliITSlayer::
SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin, Double_t ymax) {
  //--------------------------------------------------------------------
  // This function sets the "window"
  //--------------------------------------------------------------------
  fI=FindClusterIndex(zmin); fZmax=zmax;
  fImax = TMath::Min(FindClusterIndex(zmax)+1,fN);
  Double_t circle=2*TMath::Pi()*fR;
  if (ymax>circle) { ymax-=circle; ymin-=circle; }
  fYmin=ymin; fYmax=ymax;
  fSkip = 0;
}
*/


void AliITStrackerMI::AliITSlayer::
SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin, Double_t ymax) {
  //--------------------------------------------------------------------
  // This function sets the "window"
  //--------------------------------------------------------------------
 
  Double_t circle=2*TMath::Pi()*fR;
  fYmin = ymin; fYmax =ymax;
  Float_t ymiddle = (fYmin+fYmax)*0.5;
  if (ymiddle<fYB[0]) {fYmin+=circle; fYmax+=circle;ymiddle+=circle;}
  else{
    if (ymiddle>fYB[1]) {fYmin-=circle; fYmax-=circle;ymiddle-=circle;}
  }
  //
  fCurrentSlice =-1;
  // defualt take all
  fClustersCs = fClusters;
  fClusterIndexCs = fClusterIndex;
  fYcs  = fY;
  fZcs  = fZ;
  fNcs  = fN;
  //
  //is in 20 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy20){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy20);
    if (slice<0) slice=0;
    if (slice>20) slice=20;
    Bool_t isOK = (fYmin>fBy20[slice][0]&&fYmax<fBy20[slice][1]);
    if (isOK) {
      fCurrentSlice=slice;
      fClustersCs = fClusters20[fCurrentSlice];
      fClusterIndexCs = fClusterIndex20[fCurrentSlice];
      fYcs  = fY20[fCurrentSlice];
      fZcs  = fZ20[fCurrentSlice];
      fNcs  = fN20[fCurrentSlice];
    }
  }  
  //
  //is in 10 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy10){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy10);
    if (slice<0) slice=0;
    if (slice>10) slice=10;
    Bool_t isOK = (fYmin>fBy10[slice][0]&&fYmax<fBy10[slice][1]);
    if (isOK) {
      fCurrentSlice=slice;
      fClustersCs = fClusters10[fCurrentSlice];
      fClusterIndexCs = fClusterIndex10[fCurrentSlice];
      fYcs  = fY10[fCurrentSlice];
      fZcs  = fZ10[fCurrentSlice];
      fNcs  = fN10[fCurrentSlice];
    }
  }  
  //
  //is in 5 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy5){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy5);
    if (slice<0) slice=0;
    if (slice>5) slice=5;
    Bool_t isOK = (fYmin>fBy5[slice][0]&&fYmax<fBy5[slice][1]);
    if ( isOK){
      fCurrentSlice=slice;
      fClustersCs = fClusters5[fCurrentSlice];
      fClusterIndexCs = fClusterIndex5[fCurrentSlice];
      fYcs  = fY5[fCurrentSlice];
      fZcs  = fZ5[fCurrentSlice];
      fNcs  = fN5[fCurrentSlice];
    }
  }  
  //  
  fI=FindClusterIndex(zmin); fZmax=zmax;
  fImax = TMath::Min(FindClusterIndex(zmax)+1,fNcs);
  fSkip = 0;
  fAccepted =0;
}

/*
const AliITSclusterV2 *AliITStrackerMI::AliITSlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------
  const AliITSclusterV2 *cluster=0;
  for (Int_t i=fI; i<fN; i++) {
    const AliITSclusterV2 *c=fClusters[i];
    if (c->GetZ() > fZmax) break;
    //    if (c->IsUsed()) continue;
    const AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + c->GetY();

    if (y>2.*fR*TMath::Pi()) y -= 2*fR*TMath::Pi();
    if (y>1.*fR*TMath::Pi() && fYmax<y) y -= 2*fR*TMath::Pi();

    if (y<fYmin) continue;
    if (y>fYmax) continue;
    cluster=c; ci=i;
    fI=i+1;
    break; 
  }
  return cluster;
}
*/

const AliITSclusterV2 *AliITStrackerMI::AliITSlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------

  if (fCurrentSlice<0){
    Double_t rpi2 = 2.*fR*TMath::Pi();
    for (Int_t i=fI; i<fImax; i++) {
      Double_t y = fY[i];
      if (fYmax<y) y -= rpi2;
      if (y<fYmin) continue;
      if (y>fYmax) continue;
      if (fClusters[i]->GetQ()==0&&fSkip==2) continue;
      ci=i;
      fI=i+1;
      return fClusters[i];
    }
  }
  else{
    for (Int_t i=fI; i<fImax; i++) {
      if (fYcs[i]<fYmin) continue;
      if (fYcs[i]>fYmax) continue;
      if (fClustersCs[i]->GetQ()==0&&fSkip==2) continue;
      ci=fClusterIndexCs[i];
      fI=i+1;
      return fClustersCs[i];
    }
  }
  return 0;
}



Int_t AliITStrackerMI::AliITSlayer::
FindDetectorIndex(Double_t phi, Double_t z) const {
  //--------------------------------------------------------------------
  //This function finds the detector crossed by the track
  //--------------------------------------------------------------------
  Double_t dphi=-(phi-fPhiOffset);
  if      (dphi <  0) dphi += 2*TMath::Pi();
  else if (dphi >= 2*TMath::Pi()) dphi -= 2*TMath::Pi();
  Int_t np=Int_t(dphi*fNladders*0.5/TMath::Pi()+0.5);
  if (np>=fNladders) np-=fNladders;
  if (np<0)          np+=fNladders;

  Double_t dz=fZOffset-z;
  Int_t nz=Int_t(dz*(fNdetectors-1)*0.5/fZOffset+0.5);
  if (nz>=fNdetectors) return -1;
  if (nz<0)            return -1;

  return np*fNdetectors + nz;
}

Double_t AliITStrackerMI::AliITSlayer::GetThickness(Double_t y,Double_t z,Double_t &x0)
const {
  //--------------------------------------------------------------------
  //This function returns the layer thickness at this point (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0085;
  x0=21.82;
  if (43<fR&&fR<45) { //SSD2
     Double_t dd=0.0034;
     d=dd;
     if (TMath::Abs(y-0.00)>3.40) d+=dd;
     if (TMath::Abs(y-1.90)<0.45) {d+=(0.013-0.0034);}
     if (TMath::Abs(y+1.90)<0.45) {d+=(0.013-0.0034);}
     for (Int_t i=0; i<12; i++) {
       if (TMath::Abs(z-3.9*(i+0.5))<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=0.0034; 
          break;
       }
       if (TMath::Abs(z+3.9*(i+0.5))<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=0.0034; 
          break;
       }         
       if (TMath::Abs(z-3.4-3.9*i)<0.50) {d+=(0.016-0.0034); break;}
       if (TMath::Abs(z+0.5+3.9*i)<0.50) {d+=(0.016-0.0034); break;}
     }
  } else 
  if (37<fR&&fR<41) { //SSD1
     Double_t dd=0.0034;
     d=dd;
     if (TMath::Abs(y-0.00)>3.40) d+=dd;
     if (TMath::Abs(y-1.90)<0.45) {d+=(0.013-0.0034);}
     if (TMath::Abs(y+1.90)<0.45) {d+=(0.013-0.0034);}
     for (Int_t i=0; i<11; i++) {
       if (TMath::Abs(z-3.9*i)<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=dd; 
          break;
       }
       if (TMath::Abs(z+3.9*i)<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=dd; 
          break;
       }         
       if (TMath::Abs(z-1.85-3.9*i)<0.50) {d+=(0.016-0.0034); break;}
       if (TMath::Abs(z+2.05+3.9*i)<0.50) {d+=(0.016-0.0034); break;}         
     }
  } else
  if (13<fR&&fR<26) { //SDD
     Double_t dd=0.0033;
     d=dd;
     if (TMath::Abs(y-0.00)>3.30) d+=dd;

     if (TMath::Abs(y-1.80)<0.55) {
        d+=0.016;
        for (Int_t j=0; j<20; j++) {
          if (TMath::Abs(z+0.7+1.47*j)<0.12) {d+=0.08; x0=9.; break;}
          if (TMath::Abs(z-0.7-1.47*j)<0.12) {d+=0.08; x0=9.; break;}
        } 
     }
     if (TMath::Abs(y+1.80)<0.55) {
        d+=0.016;
        for (Int_t j=0; j<20; j++) {
          if (TMath::Abs(z-0.7-1.47*j)<0.12) {d+=0.08; x0=9.; break;}
          if (TMath::Abs(z+0.7+1.47*j)<0.12) {d+=0.08; x0=9.; break;}
        } 
     }

     for (Int_t i=0; i<4; i++) {
       if (TMath::Abs(z-7.3*i)<0.60) {
          d+=dd;
          if (TMath::Abs(y-0.00)>3.30) d+=dd; 
          break;
       }
       if (TMath::Abs(z+7.3*i)<0.60) {
          d+=dd; 
          if (TMath::Abs(y-0.00)>3.30) d+=dd; 
          break;
       }
     }
  } else
  if (6<fR&&fR<8) {   //SPD2
     Double_t dd=0.0063; x0=21.5;
     d=dd;
     if (TMath::Abs(y-3.08)>0.5) d+=dd;
     //if (TMath::Abs(y-3.08)>0.45) d+=dd;
     if (TMath::Abs(y-3.03)<0.10) {d+=0.014;}
  } else
  if (3<fR&&fR<5) {   //SPD1
     Double_t dd=0.0063; x0=21.5;
     d=dd;
     if (TMath::Abs(y+0.21)>0.6) d+=dd;
     //if (TMath::Abs(y+0.21)>0.45) d+=dd;
     if (TMath::Abs(y+0.10)<0.10) {d+=0.014;}
  }

  return d;
}

Double_t AliITStrackerMI::GetEffectiveThickness(Double_t y,Double_t z) const
{
  //--------------------------------------------------------------------
  //Returns the thickness between the current layer and the vertex (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0028*3*3; //beam pipe
  Double_t x0=0;

  Double_t xn=fgLayers[fI].GetR();
  for (Int_t i=0; i<fI; i++) {
    Double_t xi=fgLayers[i].GetR();
    d+=fgLayers[i].GetThickness(y,z,x0)*xi*xi;
  }

  if (fI>1) {
    Double_t xi=9.;
    d+=0.0097*xi*xi;
  }

  if (fI>3) {
    Double_t xi=0.5*(fgLayers[3].GetR()+fgLayers[4].GetR());
    d+=0.0034*xi*xi;
  }

  return d/(xn*xn);
}

Int_t AliITStrackerMI::AliITSlayer::InRoad() const {
  //--------------------------------------------------------------------
  // This function returns number of clusters within the "window" 
  //--------------------------------------------------------------------
  Int_t ncl=0;
  for (Int_t i=fI; i<fN; i++) {
    const AliITSclusterV2 *c=fClusters[i];
    if (c->GetZ() > fZmax) break;
    if (c->IsUsed()) continue;
    const AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + c->GetY();

    if (y>2.*fR*TMath::Pi()) y -= 2*fR*TMath::Pi();
    if (y>1.*fR*TMath::Pi() && fYmax<y) y -= 2*fR*TMath::Pi();

    if (y<fYmin) continue;
    if (y>fYmax) continue;
    ncl++;
  }
  return ncl;
}

Bool_t 
AliITStrackerMI::RefitAt(Double_t xx,AliITStrackV2 *t,const AliITStrackV2 *c) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"
  //--------------------------------------------------------------------
  Int_t index[kMaxLayer];
  Int_t k;
  for (k=0; k<kMaxLayer; k++) index[k]=-1;
  Int_t nc=c->GetNumberOfClusters();
  for (k=0; k<nc; k++) { 
    Int_t idx=c->GetClusterIndex(k),nl=(idx&0xf0000000)>>28;
    index[nl]=idx; 
  }

  Int_t from, to, step;
  if (xx > t->GetX()) {
      from=0; to=kMaxLayer;
      step=+1;
  } else {
      from=kMaxLayer-1; to=-1;
      step=-1;
  }

  for (Int_t i=from; i != to; i += step) {
     AliITSlayer &layer=fgLayers[i];
     Double_t r=layer.GetR();
 
     {
     Double_t hI=i-0.5*step; 
     if (TMath::Abs(hI-1.5)<0.01 || TMath::Abs(hI-3.5)<0.01) {             
        Double_t rs=0.5*(fgLayers[i-step].GetR() + r);
        Double_t d=0.0034, x0=38.6; 
        if (TMath::Abs(hI-1.5)<0.01) {rs=9.; d=0.0097; x0=42;}
        if (!t->PropagateTo(rs,-step*d,x0)) {
          return kFALSE;
        }
     }
     }

     // remember old position [SR, GSI 18.02.2003]
     Double_t oldX=0., oldY=0., oldZ=0.;
     if (t->IsStartedTimeIntegral() && step==1) {
        t->GetGlobalXYZat(t->GetX(),oldX,oldY,oldZ);
     }
     //

     Double_t x,y,z;
     if (!t->GetGlobalXYZat(r,x,y,z)) { 
       return kFALSE;
     }
     Double_t phi=TMath::ATan2(y,x);
     Int_t idet=layer.FindDetectorIndex(phi,z);
     if (idet<0) { 
       return kFALSE;
     }
     const AliITSdetector &det=layer.GetDetector(idet);
     phi=det.GetPhi();
     if (!t->Propagate(phi,det.GetR())) {
       return kFALSE;
     }
     t->SetDetectorIndex(idet);

     const AliITSclusterV2 *cl=0;
     Double_t maxchi2=1000.*kMaxChi2;

     Int_t idx=index[i];
     if (idx>0) {
        const AliITSclusterV2 *c=(AliITSclusterV2 *)GetCluster(idx); 
        if (idet != c->GetDetectorIndex()) {
           idet=c->GetDetectorIndex();
           const AliITSdetector &det=layer.GetDetector(idet);
           if (!t->Propagate(det.GetPhi(),det.GetR())) {
             return kFALSE;
           }
           t->SetDetectorIndex(idet);
        }
        //Double_t chi2=t->GetPredictedChi2(c);
	Int_t layer = (idx & 0xf0000000) >> 28;;
        Double_t chi2=GetPredictedChi2MI(t,c,layer);
        if (chi2<maxchi2) { 
	  cl=c; 
	  maxchi2=chi2; 
	} else {
	  return kFALSE;
	}
     }
     /*
     if (cl==0)
     if (t->GetNumberOfClusters()>2) {
        Double_t dz=4*TMath::Sqrt(t->GetSigmaZ2()+kSigmaZ2[i]);
        Double_t dy=4*TMath::Sqrt(t->GetSigmaY2()+kSigmaY2[i]);
        Double_t zmin=t->GetZ() - dz;
        Double_t zmax=t->GetZ() + dz;
        Double_t ymin=t->GetY() + phi*r - dy;
        Double_t ymax=t->GetY() + phi*r + dy;
        layer.SelectClusters(zmin,zmax,ymin,ymax);

        const AliITSclusterV2 *c=0; Int_t ci=-1;
        while ((c=layer.GetNextCluster(ci))!=0) {
           if (idet != c->GetDetectorIndex()) continue;
           Double_t chi2=t->GetPredictedChi2(c);
           if (chi2<maxchi2) { cl=c; maxchi2=chi2; idx=ci; }
        }
     }
     */
     if (cl) {
       //if (!t->Update(cl,maxchi2,idx)) {
       if (!UpdateMI(t,cl,maxchi2,idx)) {
          return kFALSE;
       }
       t->SetSampledEdx(cl->GetQ(),t->GetNumberOfClusters()-1);
     }

     {
     Double_t x0;
     Double_t d=layer.GetThickness(t->GetY(),t->GetZ(),x0);
     t->CorrectForMaterial(-step*d,x0);
     }
                 
     // track time update [SR, GSI 17.02.2003]
     if (t->IsStartedTimeIntegral() && step==1) {
        Double_t newX, newY, newZ;
        t->GetGlobalXYZat(t->GetX(),newX,newY,newZ);
        Double_t dL2 = (oldX-newX)*(oldX-newX) + (oldY-newY)*(oldY-newY) + 
                       (oldZ-newZ)*(oldZ-newZ);
        t->AddTimeStep(TMath::Sqrt(dL2));
     }
     //

  }

  if (!t->PropagateTo(xx,0.,0.)) return kFALSE;
  return kTRUE;
}


Double_t AliITStrackerMI::GetNormalizedChi2(AliITStrackV2 * track, Int_t mode)
{
  //
  // calculate normalized chi2
  //  return NormalizedChi2(track,0);
  Float_t chi2 = 0;
  Float_t sum=0;
  Float_t *erry = GetErrY(fCurrentEsdTrack), *errz = GetErrZ(fCurrentEsdTrack);
  //  track->fdEdxMismatch=0;
  Float_t dedxmismatch =0;
  Float_t *ny = GetNy(fCurrentEsdTrack), *nz = GetNz(fCurrentEsdTrack); 
  if (mode<100){
    for (Int_t i = 0;i<6;i++){
      if (track->fClIndex[i]>0){
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else 
	  { cerry= track->fSigmaY[i]; cerrz = track->fSigmaZ[i];}
	cerry*=cerry;
	cerrz*=cerrz;	
	Float_t cchi2 = (track->fDy[i]*track->fDy[i]/cerry)+(track->fDz[i]*track->fDz[i]/cerrz);
	if (i>1){
	  Float_t ratio = track->fNormQ[i]/track->fExpQ;
	  if (ratio<0.5) {
	    cchi2+=(0.5-ratio)*10.;
	    //track->fdEdxMismatch+=(0.5-ratio)*10.;
	    dedxmismatch+=(0.5-ratio)*10.;	    
	  }
	}
	if (i<2 ||i>3){
	  AliITSclusterV2 * cl = (AliITSclusterV2*)GetCluster( track->fClIndex[i]);  
	  Double_t delta = cl->GetNy()+cl->GetNz()-ny[i]-nz[i];
	  if (delta>1) chi2 +=0.5*TMath::Min(delta/2,2.); 
	  if (i<2) chi2+=2*cl->GetDeltaProbability();
	}
	chi2+=cchi2;
	sum++;
      }
    }
    if (TMath::Abs(track->fdEdxMismatch-dedxmismatch)>0.0001){
      track->fdEdxMismatch = dedxmismatch;
    }
  }
  else{
    for (Int_t i = 0;i<4;i++){
      if (track->fClIndex[i]>0){
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else { cerry= track->fSigmaY[i]; cerrz = track->fSigmaZ[i];}
	cerry*=cerry;
	cerrz*=cerrz;
	chi2+= (track->fDy[i]*track->fDy[i]/cerry);
	chi2+= (track->fDz[i]*track->fDz[i]/cerrz);      
	sum++;
      }
    }
    for (Int_t i = 4;i<6;i++){
      if (track->fClIndex[i]>0){	
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else { cerry= track->fSigmaY[i]; cerrz = track->fSigmaZ[i];}
	cerry*=cerry;
	cerrz*=cerrz;	
	Float_t cerryb, cerrzb;
	if (ny[i+6]>0) {cerryb = erry[i+6]; cerrzb=errz[i+6];}
	else { cerryb= track->fSigmaY[i+6]; cerrzb = track->fSigmaZ[i+6];}
	cerryb*=cerryb;
	cerrzb*=cerrzb;
	chi2+= TMath::Min((track->fDy[i+6]*track->fDy[i+6]/cerryb),track->fDy[i]*track->fDy[i]/cerry);
	chi2+= TMath::Min((track->fDz[i+6]*track->fDz[i+6]/cerrzb),track->fDz[i]*track->fDz[i]/cerrz);      
	sum++;
      }
    }
  }
  if (track->fESDtrack->GetTPCsignal()>85){
    Float_t ratio = track->fdEdx/track->fESDtrack->GetTPCsignal();
    if (ratio<0.5) {
      chi2+=(0.5-ratio)*5.;      
    }
    if (ratio>2){
      chi2+=(ratio-2.0)*3; 
    }
  }
  //
  Double_t match = TMath::Sqrt(track->fChi22);
  if (track->fConstrain)  match/=track->GetNumberOfClusters();
  if (!track->fConstrain) match/=track->GetNumberOfClusters()-2.;
  if (match<0) match=0;
  Float_t deadzonefactor = (track->fNDeadZone>0) ? 3*(1.1-track->fDeadZoneProbability):0.;
  Double_t normchi2 = 2*track->fNSkipped+match+deadzonefactor+(1+(2*track->fNSkipped+deadzonefactor)/track->GetNumberOfClusters())*
    (chi2)/TMath::Max(double(sum-track->fNSkipped),
				1./(1.+track->fNSkipped));     
 
 return normchi2;
}


Double_t AliITStrackerMI::GetMatchingChi2(AliITStrackV2 * track1, AliITStrackV2 * track2)
{
  //
  // return matching chi2 between two tracks
  AliITStrackV2 track3(*track2);
  track3.Propagate(track1->GetAlpha(),track1->GetX());
  TMatrixD vec(5,1);
  vec(0,0)=track1->fP0-track3.fP0;
  vec(1,0)=track1->fP1-track3.fP1;
  vec(2,0)=track1->fP2-track3.fP2;
  vec(3,0)=track1->fP3-track3.fP3;
  vec(4,0)=track1->fP4-track3.fP4;
  //
  TMatrixD cov(5,5);
  cov(0,0) = track1->fC00+track3.fC00;
  cov(1,1) = track1->fC11+track3.fC11;
  cov(2,2) = track1->fC22+track3.fC22;
  cov(3,3) = track1->fC33+track3.fC33;
  cov(4,4) = track1->fC44+track3.fC44;
  
  cov(0,1)=cov(1,0) = track1->fC10+track3.fC10;
  cov(0,2)=cov(2,0) = track1->fC20+track3.fC20;
  cov(0,3)=cov(3,0) = track1->fC30+track3.fC30;
  cov(0,4)=cov(4,0) = track1->fC40+track3.fC40;
  //
  cov(1,2)=cov(2,1) = track1->fC21+track3.fC21;
  cov(1,3)=cov(3,1) = track1->fC31+track3.fC31;
  cov(1,4)=cov(4,1) = track1->fC41+track3.fC41;
  //
  cov(2,3)=cov(3,2) = track1->fC32+track3.fC32;
  cov(2,4)=cov(4,2) = track1->fC42+track3.fC42;
  //
  cov(3,4)=cov(4,3) = track1->fC43+track3.fC43;
  
  cov.Invert();
  TMatrixD vec2(cov,TMatrixD::kMult,vec);
  TMatrixD chi2(vec2,TMatrixD::kTransposeMult,vec);
  return chi2(0,0);
}

Double_t  AliITStrackerMI::GetDeadZoneProbability(Double_t zpos, Double_t zerr)
{
  //
  //  return probability that given point - characterized by z position and error  is in dead zone
  //
  Double_t probability =0;
  Double_t absz = TMath::Abs(zpos);
  Double_t nearestz = (absz<2)? 0.:7.1;
  if (TMath::Abs(absz-nearestz)>0.25+3*zerr) return 0;
  Double_t zmin=0, zmax=0;   
  if (zpos<-6.){
    zmin = -7.25; zmax = -6.95; 
  }
  if (zpos>6){
    zmin = 7.0; zmax =7.3;
  }
  if (absz<2){
    zmin = -0.75; zmax = 1.5;
  }
  probability = (TMath::Erf((zpos-zmin)/zerr) - TMath::Erf((zpos-zmax)/zerr))*0.5;
  return probability;
}


Double_t AliITStrackerMI::GetTruncatedChi2(AliITStrackV2 * track, Float_t fac)
{
  //
  // calculate normalized chi2
  Float_t chi2[6];
  Float_t *erry = GetErrY(fCurrentEsdTrack), *errz = GetErrZ(fCurrentEsdTrack);
  Float_t ncl = 0;
  for (Int_t i = 0;i<6;i++){
    if (TMath::Abs(track->fDy[i])>0){      
      chi2[i]= (track->fDy[i]/erry[i])*(track->fDy[i]/erry[i]);
      chi2[i]+= (track->fDz[i]/errz[i])*(track->fDz[i]/errz[i]);
      ncl++;
    }
    else{chi2[i]=10000;}
  }
  Int_t index[6];
  TMath::Sort(6,chi2,index,kFALSE);
  Float_t max = float(ncl)*fac-1.;
  Float_t sumchi=0, sumweight=0; 
  for (Int_t i=0;i<max+1;i++){
    Float_t weight = (i<max)?1.:(max+1.-i);
    sumchi+=weight*chi2[index[i]];
    sumweight+=weight;
  }
  Double_t normchi2 = sumchi/sumweight;
  return normchi2;
}


Double_t AliITStrackerMI::GetInterpolatedChi2(AliITStrackV2 * forwardtrack, AliITStrackV2 * backtrack)
{
  //
  // calculate normalized chi2
  //  if (forwardtrack->fNUsed>0.3*float(forwardtrack->GetNumberOfClusters())) return 10000;
  Int_t npoints = 0;
  Double_t res =0;
  for (Int_t i=0;i<6;i++){
    if ( (backtrack->fSigmaY[i]<0.000000001) || (forwardtrack->fSigmaY[i]<0.000000001)) continue;
    Double_t sy1 = forwardtrack->fSigmaY[i];
    Double_t sz1 = forwardtrack->fSigmaZ[i];
    Double_t sy2 = backtrack->fSigmaY[i];
    Double_t sz2 = backtrack->fSigmaZ[i];
    if (i<2){ sy2=1000.;sz2=1000;}
    //    
    Double_t dy0 = (forwardtrack->fDy[i]/(sy1*sy1) +backtrack->fDy[i]/(sy2*sy2))/(1./(sy1*sy1)+1./(sy2*sy2));
    Double_t dz0 = (forwardtrack->fDz[i]/(sz1*sz1) +backtrack->fDz[i]/(sz2*sz2))/(1./(sz1*sz1)+1./(sz2*sz2));
    // 
    Double_t nz0 = dz0*TMath::Sqrt((1./(sz1*sz1)+1./(sz2*sz2)));
    Double_t ny0 = dy0*TMath::Sqrt((1./(sy1*sy1)+1./(sy2*sy2)));
    //
    res+= nz0*nz0+ny0*ny0;
    npoints++;
  }
  if (npoints>1) return 
		   TMath::Max(TMath::Abs(0.3*forwardtrack->Get1Pt())-0.5,0.)+
		   //2*forwardtrack->fNUsed+
		   res/TMath::Max(double(npoints-forwardtrack->fNSkipped),
				  1./(1.+forwardtrack->fNSkipped));
  return 1000;
}
   




Float_t  *AliITStrackerMI::GetWeight(Int_t index) {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fgLayers[l].GetWeight(c);
}

void AliITStrackerMI::RegisterClusterTracks(AliITStrackV2* track,Int_t id)
{
  //---------------------------------------------
  // register track to the list
  for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){
    Int_t index = track->GetClusterIndex(icluster);
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].fN) continue;
    for (Int_t itrack=0;itrack<4;itrack++){
      if (fgLayers[l].fClusterTracks[itrack][c]<0){
	fgLayers[l].fClusterTracks[itrack][c]=id;
	break;
      }
    }
  }
}
void AliITStrackerMI::UnRegisterClusterTracks(AliITStrackV2* track, Int_t id)
{
  //---------------------------------------------
  // unregister track from the list
  for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){
    Int_t index = track->GetClusterIndex(icluster);
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].fN) continue;
    for (Int_t itrack=0;itrack<4;itrack++){
      if (fgLayers[l].fClusterTracks[itrack][c]==id){
	fgLayers[l].fClusterTracks[itrack][c]=-1;
      }
    }
  }
}
Float_t AliITStrackerMI::GetNumberOfSharedClusters(AliITStrackV2* track,Int_t id, Int_t list[6], AliITSclusterV2 *clist[6])
{
  //-------------------------------------------------------------
  //get number of shared clusters
  //-------------------------------------------------------------
  Float_t shared=0;
  for (Int_t i=0;i<6;i++) { list[i]=-1, clist[i]=0;}
  // mean  number of clusters
  Float_t *ny = GetNy(id), *nz = GetNz(id); 

 
  for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){
    Int_t index = track->GetClusterIndex(icluster);
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].fN) continue;
    if (ny[l]==0){
      printf("problem\n");
    }
    AliITSclusterV2 *cl = (AliITSclusterV2*)GetCluster(index);
    Float_t weight=1;
    //
    Float_t deltan = 0;
    if (l>3&&cl->GetNy()+cl->GetNz()>6) continue;
    if (l>2&&track->fNormQ[l]/track->fExpQ>3.5) continue;
    if (l<2 || l>3){      
      deltan = (cl->GetNy()+cl->GetNz()-ny[l]-nz[l]);
    }
    else{
      deltan = (cl->GetNz()-nz[l]);
    }
    if (deltan>2.0) continue;  // extended - highly probable shared cluster
    weight = 2./TMath::Max(3.+deltan,2.);
    //
    for (Int_t itrack=0;itrack<4;itrack++){
      if (fgLayers[l].fClusterTracks[itrack][c]>=0 && fgLayers[l].fClusterTracks[itrack][c]!=id){
	list[l]=index;
	clist[l] = (AliITSclusterV2*)GetCluster(index);
	shared+=weight; 
	break;
      }
    }
  }
  track->fNUsed=shared;
  return shared;
}

Int_t AliITStrackerMI::GetOverlapTrack(AliITStrackV2 *track, Int_t trackID, Int_t &shared, Int_t clusterlist[6],Int_t overlist[6])
{
  //
  // find first shared track 
  //
  // mean  number of clusters
  Float_t *ny = GetNy(trackID), *nz = GetNz(trackID); 
  //
  for (Int_t i=0;i<6;i++) overlist[i]=-1;
  Int_t sharedtrack=100000;
  Int_t tracks[24],trackindex=0;
  for (Int_t i=0;i<24;i++) {tracks[i]=-1;}
  //
  for (Int_t icluster=0;icluster<6;icluster++){
    if (clusterlist[icluster]<0) continue;
    Int_t index = clusterlist[icluster];
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (ny[l]==0){
      printf("problem\n");
    }
    if (c>fgLayers[l].fN) continue;
    //if (l>3) continue;
    AliITSclusterV2 *cl = (AliITSclusterV2*)GetCluster(index);
    //
    Float_t deltan = 0;
    if (l>3&&cl->GetNy()+cl->GetNz()>6) continue;
    if (l>2&&track->fNormQ[l]/track->fExpQ>3.5) continue;
    if (l<2 || l>3){      
      deltan = (cl->GetNy()+cl->GetNz()-ny[l]-nz[l]);
    }
    else{
      deltan = (cl->GetNz()-nz[l]);
    }
    if (deltan>2.0) continue;  // extended - highly probable shared cluster
    //
    for (Int_t itrack=3;itrack>=0;itrack--){
      if (fgLayers[l].fClusterTracks[itrack][c]<0) continue;
      if (fgLayers[l].fClusterTracks[itrack][c]!=trackID){
       tracks[trackindex]  = fgLayers[l].fClusterTracks[itrack][c];
       trackindex++;
      }
    }
  }
  if (trackindex==0) return -1;
  if (trackindex==1){    
    sharedtrack = tracks[0];
  }else{
    if (trackindex==2) sharedtrack =TMath::Min(tracks[0],tracks[1]);
    else{
      //
      Int_t track[24], cluster[24];
      for (Int_t i=0;i<trackindex;i++){ track[i]=-1; cluster[i]=0;}
      Int_t index =0;
      //
      for (Int_t i=0;i<trackindex;i++){
	if (tracks[i]<0) continue;
	track[index] = tracks[i];
	cluster[index]++;	
	for (Int_t j=i+1;j<trackindex;j++){
	  if (tracks[j]<0) continue;
	  if (tracks[j]==tracks[i]){
	    cluster[index]++;
	    tracks[j]=-1;
	  }
	}
	index++;
      }
      Int_t max=0;
      for (Int_t i=0;i<index;i++){
	if (cluster[index]>max) {
	  sharedtrack=track[index];
	  max=cluster[index];
	}
      }
    }
  }
  
  if (sharedtrack>=100000) return -1;
  //
  // make list of overlaps
  shared =0;
  for (Int_t icluster=0;icluster<6;icluster++){
    if (clusterlist[icluster]<0) continue;
    Int_t index = clusterlist[icluster];
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].fN) continue;
    AliITSclusterV2 *cl = (AliITSclusterV2*)GetCluster(index);
    if (l==0 || l==1){
      if (cl->GetNy()>2) continue;
      if (cl->GetNz()>2) continue;
    }
    if (l==4 || l==5){
      if (cl->GetNy()>3) continue;
      if (cl->GetNz()>3) continue;
    }
    //
    for (Int_t itrack=3;itrack>=0;itrack--){
      if (fgLayers[l].fClusterTracks[itrack][c]<0) continue;
      if (fgLayers[l].fClusterTracks[itrack][c]==sharedtrack){
	overlist[l]=index;
	shared++;      
      }
    }
  }
  return sharedtrack;
}


AliITStrackV2 *  AliITStrackerMI::GetBest2Tracks(Int_t trackID1, Int_t trackID2, Float_t th0, Float_t th1){
  //
  // try to find track hypothesys without conflicts
  // with minimal chi2;
  TClonesArray *arr1 = (TClonesArray*)fTrackHypothesys.At(trackID1);
  Int_t entries1 = arr1->GetEntriesFast();
  TClonesArray *arr2 = (TClonesArray*)fTrackHypothesys.At(trackID2);
  if (!arr2) return (AliITStrackV2*) arr1->UncheckedAt(0);
  Int_t entries2 = arr2->GetEntriesFast();
  if (entries2<=0) return (AliITStrackV2*) arr1->UncheckedAt(0);
  //
  AliITStrackV2 * track10=(AliITStrackV2*) arr1->UncheckedAt(0);
  AliITStrackV2 * track20=(AliITStrackV2*) arr2->UncheckedAt(0);
  if (TMath::Abs(1./track10->Get1Pt())>0.5+TMath::Abs(1/track20->Get1Pt())) return track10;

  for (Int_t itrack=0;itrack<entries1;itrack++){
    AliITStrackV2 * track=(AliITStrackV2*) arr1->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID1);
  }
  //
  for (Int_t itrack=0;itrack<entries2;itrack++){
    AliITStrackV2 * track=(AliITStrackV2*) arr2->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID2);
  }
  Int_t index1=0;
  Int_t index2=0;
  Float_t maxconflicts=6;
  Double_t maxchi2 =1000.;
  //
  // get weight of hypothesys - using best hypothesy
  Double_t w1,w2;
 
  Int_t list1[6],list2[6];
  AliITSclusterV2 *clist1[6], *clist2[6] ;
  RegisterClusterTracks(track10,trackID1);
  RegisterClusterTracks(track20,trackID2);
  Float_t conflict1 = GetNumberOfSharedClusters(track10,trackID1,list1,clist1);
  Float_t conflict2 = GetNumberOfSharedClusters(track20,trackID2,list2,clist2);
  UnRegisterClusterTracks(track10,trackID1);
  UnRegisterClusterTracks(track20,trackID2);
  //
  // normalized chi2
  Float_t chi21 =0,chi22=0,ncl1=0,ncl2=0;
  Float_t nerry[6],nerrz[6];
  Float_t *erry1=GetErrY(trackID1),*errz1 = GetErrZ(trackID1);
  Float_t *erry2=GetErrY(trackID2),*errz2 = GetErrZ(trackID2);
  for (Int_t i=0;i<6;i++){
     if ( (erry1[i]>0) && (erry2[i]>0)) {
       nerry[i] = TMath::Min(erry1[i],erry2[i]);
       nerrz[i] = TMath::Min(errz1[i],errz2[i]);
     }else{
       nerry[i] = TMath::Max(erry1[i],erry2[i]);
       nerrz[i] = TMath::Max(errz1[i],errz2[i]);
     }
     if (TMath::Abs(track10->fDy[i])>0.000000000000001){
       chi21 += track10->fDy[i]*track10->fDy[i]/(nerry[i]*nerry[i]);
       chi21 += track10->fDz[i]*track10->fDz[i]/(nerrz[i]*nerrz[i]);
       ncl1++;
     }
     if (TMath::Abs(track20->fDy[i])>0.000000000000001){
       chi22 += track20->fDy[i]*track20->fDy[i]/(nerry[i]*nerry[i]);
       chi22 += track20->fDz[i]*track20->fDz[i]/(nerrz[i]*nerrz[i]);
       ncl2++;
     }
  }
  chi21/=ncl1;
  chi22/=ncl2;
  //
  // 
  Float_t d1 = TMath::Sqrt(track10->fD[0]*track10->fD[0]+track10->fD[1]*track10->fD[1])+0.1;
  Float_t d2 = TMath::Sqrt(track20->fD[0]*track20->fD[0]+track20->fD[1]*track20->fD[1])+0.1;
  Float_t s1 = TMath::Sqrt(track10->GetSigmaY2()*track10->GetSigmaZ2());
  Float_t s2 = TMath::Sqrt(track20->GetSigmaY2()*track20->GetSigmaZ2());
  //
  w1 = (d2/(d1+d2)+ 2*s2/(s1+s2)+
	+s2/(s1+s2)*0.5*(chi22+2.)/(chi21+chi22+4.)
	+1.*TMath::Abs(1./track10->Get1Pt())/(TMath::Abs(1./track10->Get1Pt())+TMath::Abs(1./track20->Get1Pt()))
	);
  w2 = (d1/(d1+d2)+ 2*s1/(s1+s2)+
	s1/(s1+s2)*0.5*(chi21+2.)/(chi21+chi22+4.)
	+1.*TMath::Abs(1./track20->Get1Pt())/(TMath::Abs(1./track10->Get1Pt())+TMath::Abs(1./track20->Get1Pt()))
	);

  Double_t sumw = w1+w2;
  w1/=sumw;
  w2/=sumw;
  if (w1<w2*0.5) {
    w1 = (d2+0.5)/(d1+d2+1.);
    w2 = (d1+0.5)/(d1+d2+1.);
  }
  //  Float_t maxmax       = w1*track10->fChi2MIP[0]+w2*track20->fChi2MIP[0]+w1*conflict1+w2*conflict2+1.;
  //Float_t maxconflicts0 = w1*conflict1+w2*conflict2;
  //
  // get pair of "best" hypothesys
  //  
  Float_t * ny1 = GetNy(trackID1), * nz1 = GetNz(trackID1); 
  Float_t * ny2 = GetNy(trackID2), * nz2 = GetNz(trackID2); 

  for (Int_t itrack1=0;itrack1<entries1;itrack1++){
    AliITStrackV2 * track1=(AliITStrackV2*) arr1->UncheckedAt(itrack1);
    //if (track1->fFakeRatio>0) continue;
    RegisterClusterTracks(track1,trackID1);
    for (Int_t itrack2=0;itrack2<entries2;itrack2++){
      AliITStrackV2 * track2=(AliITStrackV2*) arr2->UncheckedAt(itrack2);

      //      Float_t current = w1*track1->fChi2MIP[0]+w2*track2->fChi2MIP[0];
      //if (track2->fFakeRatio>0) continue;
      Float_t nskipped=0;            
      RegisterClusterTracks(track2,trackID2);
      Int_t list1[6],list2[6];
      AliITSclusterV2 *clist1[6], *clist2[6] ;
      Float_t cconflict1 = GetNumberOfSharedClusters(track1,trackID1,list1,clist1);
      Float_t cconflict2 = GetNumberOfSharedClusters(track2,trackID2,list2,clist2);
      UnRegisterClusterTracks(track2,trackID2);
      //
      if (track1->fConstrain) nskipped+=w1*track1->fNSkipped;
      if (track2->fConstrain) nskipped+=w2*track2->fNSkipped;
      if (nskipped>0.5) continue;
      //
      //if ( w1*conflict1+w2*conflict2>maxconflicts0) continue;
      if (conflict1+1<cconflict1) continue;
      if (conflict2+1<cconflict2) continue;
      Float_t conflict=0;
      Float_t sumchi2=0;
      Float_t sum=0;
      for (Int_t i=0;i<6;i++){
	//
	Float_t c1 =0.; // conflict coeficients
	Float_t c2 =0.; 
	if (clist1[i]&&clist2[i]){
	  Float_t deltan = 0;
	  if (i<2 || i>3){      
	    deltan = (clist1[i]->GetNy()+clist1[i]->GetNz()-TMath::Max(ny1[i],ny2[i])-TMath::Max(nz1[i],nz2[i]));
	  }
	  else{
	    deltan = (clist1[i]->GetNz()-TMath::Max(nz1[i],nz2[i]));
	  }
	  c1  = 2./TMath::Max(3.+deltan,2.);	  
	  c2  = 2./TMath::Max(3.+deltan,2.);	  
	}
	else{
	  if (clist1[i]){
	    Float_t deltan = 0;
	    if (i<2 || i>3){      
	      deltan = (clist1[i]->GetNy()+clist1[i]->GetNz()-ny1[i]-nz1[i]);
	    }
	    else{
	      deltan = (clist1[i]->GetNz()-nz1[i]);
	    }
	    c1  = 2./TMath::Max(3.+deltan,2.);	  
	    c2  = 0;
	  }

	  if (clist2[i]){
	    Float_t deltan = 0;
	    if (i<2 || i>3){      
	      deltan = (clist2[i]->GetNy()+clist2[i]->GetNz()-ny2[i]-nz2[i]);
	    }
	    else{
	      deltan = (clist2[i]->GetNz()-nz2[i]);
	    }
	    c2  = 2./TMath::Max(3.+deltan,2.);	  
	    c1  = 0;
	  }	  
	}
	//
	Double_t chi21=0,chi22=0;
	if (TMath::Abs(track1->fDy[i])>0.) {
	  chi21 = (track1->fDy[i]/track1->fSigmaY[i])*(track1->fDy[i]/track1->fSigmaY[i])+
	    (track1->fDz[i]/track1->fSigmaZ[i])*(track1->fDz[i]/track1->fSigmaZ[i]);
	  //chi21 = (track1->fDy[i]*track1->fDy[i])/(nerry[i]*nerry[i])+
	  //  (track1->fDz[i]*track1->fDz[i])/(nerrz[i]*nerrz[i]);
	}else{
	  if (TMath::Abs(track1->fSigmaY[i]>0.)) c1=1;
	}
	//
	if (TMath::Abs(track2->fDy[i])>0.) {
	  chi22 = (track2->fDy[i]/track2->fSigmaY[i])*(track2->fDy[i]/track2->fSigmaY[i])+
	    (track2->fDz[i]/track2->fSigmaZ[i])*(track2->fDz[i]/track2->fSigmaZ[i]);
	  //chi22 = (track2->fDy[i]*track2->fDy[i])/(nerry[i]*nerry[i])+
	  //  (track2->fDz[i]*track2->fDz[i])/(nerrz[i]*nerrz[i]);
	}
	else{
	  if (TMath::Abs(track2->fSigmaY[i]>0.)) c2=1;
	}
	sumchi2+=w1*(1.+c1)*(1+c1)*(chi21+c1)+w2*(1.+c2)*(1+c2)*(chi22+c2);
	if (chi21>0) sum+=w1;
	if (chi22>0) sum+=w2;
	conflict+=(c1+c2);
      }
      Double_t norm = sum-w1*track1->fNSkipped-w2*track2->fNSkipped;
      if (norm<0) norm =1/(w1*track1->fNSkipped+w2*track2->fNSkipped);
      Double_t normchi2 = 2*conflict+sumchi2/sum;
      if ( normchi2 <maxchi2 ){	  
	index1 = itrack1;
	index2 = itrack2;
	maxconflicts = conflict;
	maxchi2 = normchi2;
      }      
    }
    UnRegisterClusterTracks(track1,trackID1);
  }
  //
  //  if (maxconflicts<4 && maxchi2<th0){   
  if (maxchi2<th0*2.){   
    Float_t orig = track10->fFakeRatio*track10->GetNumberOfClusters();
    AliITStrackV2* track1=(AliITStrackV2*) arr1->UncheckedAt(index1);
    track1->fChi2MIP[5] = maxconflicts;
    track1->fChi2MIP[6] = maxchi2;
    track1->fChi2MIP[7] = 0.01+orig-(track1->fFakeRatio*track1->GetNumberOfClusters());
    //    track1->UpdateESDtrack(AliESDtrack::kITSin);
    track1->fChi2MIP[8] = index1;
    fBestTrackIndex[trackID1] =index1;
    UpdateESDtrack(track1, AliESDtrack::kITSin);
  }  
  else if (track10->fChi2MIP[0]<th1){
    track10->fChi2MIP[5] = maxconflicts;
    track10->fChi2MIP[6] = maxchi2;    
    //    track10->UpdateESDtrack(AliESDtrack::kITSin);
    UpdateESDtrack(track10,AliESDtrack::kITSin);
  }   
  
  for (Int_t itrack=0;itrack<entries1;itrack++){
    AliITStrackV2 * track=(AliITStrackV2*) arr1->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID1);
  }
  //
  for (Int_t itrack=0;itrack<entries2;itrack++){
    AliITStrackV2 * track=(AliITStrackV2*) arr2->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID2);
  }

  if (track10->fConstrain&&track10->fChi2MIP[0]<kMaxChi2PerCluster[0]&&track10->fChi2MIP[1]<kMaxChi2PerCluster[1]
      &&track10->fChi2MIP[2]<kMaxChi2PerCluster[2]&&track10->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    //  if (track10->fChi2MIP[0]<kMaxChi2PerCluster[0]&&track10->fChi2MIP[1]<kMaxChi2PerCluster[1]
  //    &&track10->fChi2MIP[2]<kMaxChi2PerCluster[2]&&track10->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    RegisterClusterTracks(track10,trackID1);
  }
  if (track20->fConstrain&&track20->fChi2MIP[0]<kMaxChi2PerCluster[0]&&track20->fChi2MIP[1]<kMaxChi2PerCluster[1]
      &&track20->fChi2MIP[2]<kMaxChi2PerCluster[2]&&track20->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    //if (track20->fChi2MIP[0]<kMaxChi2PerCluster[0]&&track20->fChi2MIP[1]<kMaxChi2PerCluster[1]
    //  &&track20->fChi2MIP[2]<kMaxChi2PerCluster[2]&&track20->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    RegisterClusterTracks(track20,trackID2);  
  }
  return track10; 
 
}

void AliITStrackerMI::UseClusters(const AliKalmanTrack *t, Int_t from) const {
  //--------------------------------------------------------------------
  // This function marks clusters assigned to the track
  //--------------------------------------------------------------------
  AliTracker::UseClusters(t,from);

  AliITSclusterV2 *c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(0));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();
  c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(1));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();

}


void AliITStrackerMI::AddTrackHypothesys(AliITStrackV2 * track, Int_t esdindex)
{
  //------------------------------------------------------------------
  // add track to the list of hypothesys
  //------------------------------------------------------------------

  if (esdindex>=fTrackHypothesys.GetEntriesFast()) fTrackHypothesys.Expand(esdindex*2+10);
  //
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) {
    array = new TObjArray(10);
    fTrackHypothesys.AddAt(array,esdindex);
  }
  array->AddLast(track);
}

void AliITStrackerMI::SortTrackHypothesys(Int_t esdindex, Int_t maxcut, Int_t mode)
{
  //-------------------------------------------------------------------
  // compress array of track hypothesys
  // keep only maxsize best hypothesys
  //-------------------------------------------------------------------
  if (esdindex>fTrackHypothesys.GetEntriesFast()) return;
  if (! (fTrackHypothesys.At(esdindex)) ) return;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  Int_t entries = array->GetEntriesFast();
  //
  //- find preliminary besttrack as a reference
  Float_t minchi2=10000;
  Int_t maxn=0;
  AliITStrackV2 * besttrack=0;
  for (Int_t itrack=0;itrack<array->GetEntriesFast();itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track) continue;
    Float_t chi2 = NormalizedChi2(track,0);
    //
    Int_t tpcLabel=track->fESDtrack->GetTPCLabel();
    track->SetLabel(tpcLabel);
    CookdEdx(track);
    track->fFakeRatio=1.;
    CookLabel(track,0.); //For comparison only
    //
    //if (chi2<kMaxChi2PerCluster[0]&&track->fFakeRatio==0){
    if (chi2<kMaxChi2PerCluster[0]){
      if (track->GetNumberOfClusters()<maxn) continue;
      maxn = track->GetNumberOfClusters();
      if (chi2<minchi2){
	minchi2=chi2;
	besttrack=track;
      }
    }
    else{
      delete array->RemoveAt(itrack);
    }
  }
  if (!besttrack) return;
  //
  //
  //take errors of best track as a reference
  Float_t *erry = GetErrY(esdindex), *errz = GetErrZ(esdindex);
  Float_t *ny = GetNy(esdindex), *nz = GetNz(esdindex);
  for (Int_t i=0;i<6;i++) {
    if (besttrack->fClIndex[i]>0){
      erry[i] = besttrack->fSigmaY[i]; erry[i+6] = besttrack->fSigmaY[i+6];
      errz[i] = besttrack->fSigmaZ[i]; errz[i+6] = besttrack->fSigmaZ[i+6];
      ny[i]   = besttrack->fNy[i];
      nz[i]   = besttrack->fNz[i];
    }
  }
  //
  // calculate normalized chi2
  //
  Float_t * chi2        = new Float_t[entries];
  Int_t * index         = new Int_t[entries];  
  for (Int_t i=0;i<entries;i++) chi2[i] =10000;
  for (Int_t itrack=0;itrack<entries;itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (track){
      track->fChi2MIP[0] = GetNormalizedChi2(track, mode);            
      if (track->fChi2MIP[0]<kMaxChi2PerCluster[0]) 
	chi2[itrack] = track->fChi2MIP[0];
      else
	delete array->RemoveAt(itrack);	     
    }
  }
  //
  TMath::Sort(entries,chi2,index,kFALSE);
  besttrack = (AliITStrackV2*)array->At(index[0]);
  if (besttrack&&besttrack->fChi2MIP[0]<kMaxChi2PerCluster[0]){
    for (Int_t i=0;i<6;i++){
      if (besttrack->fClIndex[i]>0){
	erry[i] = besttrack->fSigmaY[i]; erry[i+6] = besttrack->fSigmaY[i+6];
	errz[i] = besttrack->fSigmaZ[i]; erry[i+6] = besttrack->fSigmaY[i+6];
	ny[i]   = besttrack->fNy[i];
	nz[i]   = besttrack->fNz[i];
      }
    }
  }
  //
  // calculate one more time with updated normalized errors
  for (Int_t i=0;i<entries;i++) chi2[i] =10000;  
  for (Int_t itrack=0;itrack<entries;itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (track){      
      track->fChi2MIP[0] = GetNormalizedChi2(track,mode);            
      if (track->fChi2MIP[0]<kMaxChi2PerCluster[0]) 
	chi2[itrack] = track->fChi2MIP[0]-0*(track->GetNumberOfClusters()+track->fNDeadZone); 
      else 
	delete array->RemoveAt(itrack);	
    }   
  }
  entries = array->GetEntriesFast();  
  //
  if (entries>0){
    TObjArray * newarray = new TObjArray();  
    TMath::Sort(entries,chi2,index,kFALSE);
    besttrack = (AliITStrackV2*)array->At(index[0]);
    if (besttrack){
      //
      for (Int_t i=0;i<6;i++){
	if (besttrack->fNz[i]>0&&besttrack->fNy[i]>0){
	  erry[i] = besttrack->fSigmaY[i]; erry[i+6] = besttrack->fSigmaY[i+6];
	  errz[i] = besttrack->fSigmaZ[i]; errz[i+6] = besttrack->fSigmaZ[i+6];
	  ny[i]   = besttrack->fNy[i];
	  nz[i]   = besttrack->fNz[i];
	}
      }
      besttrack->fChi2MIP[0] = GetNormalizedChi2(besttrack,mode);
      Float_t minchi2 = TMath::Min(besttrack->fChi2MIP[0]+5.+besttrack->fNUsed, double(kMaxChi2PerCluster[0]));
      Float_t minn = besttrack->GetNumberOfClusters()-3;
      Int_t accepted=0;
      for (Int_t i=0;i<entries;i++){
	AliITStrackV2 * track = (AliITStrackV2*)array->At(index[i]);	
	if (!track) continue;
	if (accepted>maxcut) break;
	track->fChi2MIP[0] = GetNormalizedChi2(track,mode);
	if (track->GetNumberOfClusters()<6 && (track->fChi2MIP[0]+track->fNUsed>minchi2)){
	  delete array->RemoveAt(index[i]);
	  continue;
	}
	if (track->fChi2MIP[0]+track->fNUsed<minchi2 && track->GetNumberOfClusters()>=minn){
	  accepted++;
	  //
	  newarray->AddLast(array->RemoveAt(index[i]));      
	  for (Int_t i=0;i<6;i++){
	    if (nz[i]==0){
	      erry[i] = track->fSigmaY[i]; erry[i+6] = track->fSigmaY[i+6];
	      errz[i] = track->fSigmaZ[i]; errz[i]   = track->fSigmaZ[i+6];
	      ny[i]   = track->fNy[i];
	      nz[i]   = track->fNz[i];
	    }
	  }
	}
	else{
	  delete array->RemoveAt(index[i]);
	}
      }
      array->Delete();
      delete fTrackHypothesys.RemoveAt(esdindex);
      fTrackHypothesys.AddAt(newarray,esdindex);
    }
    else{
      array->Delete();
      delete fTrackHypothesys.RemoveAt(esdindex);
    }
  }
  delete [] chi2;
  delete [] index;
}



AliITStrackV2 * AliITStrackerMI::GetBestHypothesys(Int_t esdindex, AliITStrackV2 * original, Int_t checkmax)
{
  //-------------------------------------------------------------
  // try to find best hypothesy
  // currently - minimal chi2 of track+backpropagated track+matching to the tpc track
  //-------------------------------------------------------------
  if (fTrackHypothesys.GetEntriesFast()<=esdindex) return 0;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) return 0;
  Int_t entries = array->GetEntriesFast();
  if (!entries) return 0;  
  Float_t minchi2 = 100000;
  AliITStrackV2 * besttrack=0;
  //
  AliITStrackV2 * backtrack    = new AliITStrackV2(*original);
  AliITStrackV2 * forwardtrack = new AliITStrackV2(*original);
  //
  for (Int_t i=0;i<entries;i++){    
    AliITStrackV2 * track = (AliITStrackV2*)array->At(i);    
    if (!track) continue;
    track->fChi2MIP[1] = 1000000;
    track->fChi2MIP[2] = 1000000;
    track->fChi2MIP[3] = 1000000;
    //
    // backtrack
    backtrack = new(backtrack) AliITStrackV2(*track); 
    backtrack->ResetCovariance();
    backtrack->ResetCovariance();
    backtrack->ResetClusters();
    Double_t x = original->GetX();
    if (!RefitAt(x,backtrack,track)) continue;
    track->fChi2MIP[1] = NormalizedChi2(backtrack,0);
    //for (Int_t i=2;i<6;i++){track->fDy[i]+=backtrack->fDy[i]; track->fDz[i]+=backtrack->fDz[i];}
    if (track->fChi2MIP[1]>kMaxChi2PerCluster[1]*6.)  continue;
    track->fChi22 = GetMatchingChi2(backtrack,original);

    if ((track->fConstrain) && track->fChi22>90.)  continue;
    if ((!track->fConstrain) && track->fChi22>30.)  continue;
    if ( track->fChi22/track->GetNumberOfClusters()>11.)  continue;


    if  (!(track->fConstrain)&&track->fChi2MIP[1]>kMaxChi2PerCluster[1])  continue;
    Bool_t isOK=kTRUE;
    /*
    for (Int_t i=0;i<6;i++){
      if (track->fClIndex[i]>0){
	Double_t dy1 = (track->fDy[i]/track->fSigmaY[i]);
	Double_t dz1 = (track->fDz[i]/track->fSigmaZ[i]);
	Double_t dy2 = (backtrack->fDy[i]/backtrack->fSigmaY[i]);
	Double_t dz2 = (backtrack->fDz[i]/backtrack->fSigmaZ[i]);
	if (TMath::Min(dy1*dy1+dz1*dz1,dy2*dy2+dz2*dz2)> kMaxChi2sR[i]) isOK =kFALSE;
	track->fDy[i+6] = backtrack->fDy[i];
	track->fDz[i+6] = backtrack->fDz[i];
	track->fSigmaY[i+6] = backtrack->fSigmaY[i];
	track->fSigmaZ[i+6] = backtrack->fSigmaZ[i];
      }
      else{
	if (i==5){
	  if (track->fClIndex[i-1]>0){
	    Double_t dy2 = (backtrack->fDy[i-1]/backtrack->fSigmaY[i-1]);
	    Double_t dz2 = (backtrack->fDz[i-1]/backtrack->fSigmaZ[i-1]);
	    if (dy2*dy2+dz2*dz2> kMaxChi2sR[i]) isOK =kFALSE;
	  }
	  else isOK = kFALSE;
	}
      }
    }
    */
    if(!isOK) continue;
    //
    //forward track - without constraint
    forwardtrack = new(forwardtrack) AliITStrackV2(*original);
    forwardtrack->ResetClusters();
    x = track->GetX();
    if (!RefitAt(x,forwardtrack,track))  continue;
    track->fChi2MIP[2] = NormalizedChi2(forwardtrack,0);    
    if  (track->fChi2MIP[2]>kMaxChi2PerCluster[2]*6.0)  continue;
    if  (!(track->fConstrain)&&track->fChi2MIP[2]>kMaxChi2PerCluster[2])  continue;
    
    track->fD[0] = forwardtrack->GetD(GetX(),GetY());
    track->fD[1] = forwardtrack->GetZat(GetX())-GetZ();
    forwardtrack->fD[0] = track->fD[0];
    forwardtrack->fD[1] = track->fD[1];    
    {
      Int_t list[6];
      AliITSclusterV2* clist[6];
      track->fChi2MIP[4] = GetNumberOfSharedClusters(track,esdindex,list,clist);      
      if ( (!track->fConstrain) && track->fChi2MIP[4]>1.0) continue;
    }
    
    track->fChi2MIP[3] = GetInterpolatedChi2(forwardtrack,backtrack);
    if  ( (track->fChi2MIP[3]>6.*kMaxChi2PerCluster[3])) continue;    
    if  ( (!track->fConstrain) && (track->fChi2MIP[3]>2*kMaxChi2PerCluster[3])) {
      track->fChi2MIP[3]=1000;
      continue; 
    }
    Double_t chi2 = track->fChi2MIP[0]+track->fNUsed;    
    //
    for (Int_t ichi=0;ichi<5;ichi++){
      forwardtrack->fChi2MIP[ichi] = track->fChi2MIP[ichi];
    }
    if (chi2 < minchi2){
      //besttrack = new AliITStrackV2(*forwardtrack);
      besttrack = track;
      besttrack->SetLabel(track->GetLabel());
      besttrack->fFakeRatio = track->fFakeRatio;
      minchi2   = chi2;
      original->fD[0] = forwardtrack->GetD(GetX(),GetY());
      original->fD[1] = forwardtrack->GetZat(GetX())-GetZ();
    }    
  }
  delete backtrack;
  delete forwardtrack;
  Int_t accepted=0;
  for (Int_t i=0;i<entries;i++){    
    AliITStrackV2 * track = (AliITStrackV2*)array->At(i);   
    if (!track) continue;
    if (accepted>checkmax || track->fChi2MIP[3]>kMaxChi2PerCluster[3]*6. || 
	(track->GetNumberOfClusters()<besttrack->GetNumberOfClusters()-1.)||
	track->fChi2MIP[0]>besttrack->fChi2MIP[0]+2.*besttrack->fNUsed+3.){
      delete array->RemoveAt(i);    
      continue;
    }
    else{
      accepted++;
    }
  }
  //
  array->Compress();
  SortTrackHypothesys(esdindex,checkmax,1);
  array = (TObjArray*) fTrackHypothesys.At(esdindex);
  besttrack = (AliITStrackV2*)array->At(0);  
  if (!besttrack)  return 0;
  besttrack->fChi2MIP[8]=0;
  fBestTrackIndex[esdindex]=0;
  entries = array->GetEntriesFast();
  AliITStrackV2 *longtrack =0;
  minchi2 =1000;
  Float_t minn=besttrack->GetNumberOfClusters()+besttrack->fNDeadZone;
  for (Int_t itrack=entries-1;itrack>0;itrack--){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track->fConstrain) continue;
    if (track->GetNumberOfClusters()+track->fNDeadZone<minn) continue;
    if (track->fChi2MIP[0]-besttrack->fChi2MIP[0]>0.0) continue;
    if (track->fChi2MIP[0]>4.) continue;
    minn = track->GetNumberOfClusters()+track->fNDeadZone;
    longtrack =track;
  }
  //if (longtrack) besttrack=longtrack;

  Int_t list[6];
  AliITSclusterV2 * clist[6];
  Float_t shared = GetNumberOfSharedClusters(besttrack,esdindex,list,clist);
  if (besttrack->fConstrain&&besttrack->fChi2MIP[0]<kMaxChi2PerCluster[0]&&besttrack->fChi2MIP[1]<kMaxChi2PerCluster[1]
      &&besttrack->fChi2MIP[2]<kMaxChi2PerCluster[2]&&besttrack->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    //if (besttrack->fChi2MIP[0]<kMaxChi2PerCluster[0]&&besttrack->fChi2MIP[1]<kMaxChi2PerCluster[1]
    //  &&besttrack->fChi2MIP[2]<kMaxChi2PerCluster[2]&&besttrack->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    RegisterClusterTracks(besttrack,esdindex);
  }
  //
  //
  if (shared>0.0){
    Int_t nshared;
    Int_t overlist[6];
    Int_t sharedtrack = GetOverlapTrack(besttrack, esdindex, nshared, list, overlist);
    if (sharedtrack>=0){
      //
      besttrack = GetBest2Tracks(esdindex,sharedtrack,10,5.5);     
      if (besttrack){
	shared = GetNumberOfSharedClusters(besttrack,esdindex,list,clist);
      }
      else return 0;
    }
  }  
  
  if (shared>2.5) return 0;
  if (shared>1.0) return besttrack;
  //
  // don't sign clusters if not gold track
  Double_t deltad = besttrack->GetD(GetX(),GetY());
  Double_t deltaz = besttrack->GetZat(GetX()) - GetZ();
  Double_t deltaprim = TMath::Sqrt(deltad*deltad+deltaz*deltaz);
  if (deltaprim>0.1 && (fConstraint[fPass])) return besttrack;
  if (TMath::Abs(deltad)>0.1)  return besttrack;
  if (TMath::Abs(deltaz)>0.1)  return besttrack;
  if (besttrack->fChi2MIP[0]>4.) return besttrack;
  if (besttrack->fChi2MIP[1]>4.) return besttrack;
  if (besttrack->fChi2MIP[2]>4.) return besttrack;
  if (besttrack->fChi2MIP[3]>4.) return besttrack;
  if (fConstraint[fPass]){
    //
    // sign clusters
    //
    Float_t *ny = GetNy(esdindex), *nz = GetNz(esdindex);
    for (Int_t i=0;i<6;i++){
      Int_t index = besttrack->fClIndex[i];
      if (index<=0) continue; 
      Int_t ilayer =  (index & 0xf0000000) >> 28;
      if (besttrack->fSigmaY[ilayer]<0.00000000001) continue;
      AliITSclusterV2 *c = (AliITSclusterV2*)GetCluster(index);     
      if (!c) continue;
      if (ilayer>3&&c->GetNy()+c->GetNz()>6) continue;
      if ( (c->GetNy()+c->GetNz() )> ny[i]+nz[i]+0.7) continue; //shared track
      if (  c->GetNz()> nz[i]+0.7) continue; //shared track
      if ( ilayer>2&& besttrack->fNormQ[ilayer]/besttrack->fExpQ>1.5) continue;
      //if (  c->GetNy()> ny[i]+0.7) continue; //shared track

      Bool_t cansign = kTRUE;
      for (Int_t itrack=0;itrack<entries; itrack++){
	AliITStrackV2 * track = (AliITStrackV2*)array->At(i);   
	if (!track) continue;
	if (track->fChi2MIP[0]>besttrack->fChi2MIP[0]+2.*shared+1.) break;
	if ( (track->fClIndex[ilayer]>0) && (track->fClIndex[ilayer]!=besttrack->fClIndex[ilayer])){
	  cansign = kFALSE;
	  break;
	}
      }
      if (cansign){
	if (TMath::Abs(besttrack->fDy[ilayer]/besttrack->fSigmaY[ilayer])>3.) continue;
	if (TMath::Abs(besttrack->fDz[ilayer]/besttrack->fSigmaZ[ilayer])>3.) continue;    
	if (!c->IsUsed()) c->Use();
      }
    }
  }
  return besttrack;
} 



void  AliITStrackerMI::GetBestHypothesysMIP(TObjArray &itsTracks)
{
  //
  // get "best" hypothesys
  //for (Int_t ilayer=0;ilayer<6;ilayer++) fgLayers[ilayer].ResetWeights();


  Int_t nentries = itsTracks.GetEntriesFast();
  for (Int_t i=0;i<nentries;i++){
    AliITStrackV2* track = (AliITStrackV2*)itsTracks.At(i);
    if (!track) continue;
    TObjArray * array = (TObjArray*) fTrackHypothesys.At(i);
    if (!array) continue;
    if (array->GetEntriesFast()<=0) continue;
    //
    AliITStrackV2* longtrack=0;
    Float_t minn=0;
    Float_t maxchi2=1000;
    for (Int_t j=0;j<array->GetEntriesFast();j++){
      AliITStrackV2* track = (AliITStrackV2*)array->At(j);
      if (!track) continue;
      if (track->GetNumberOfClusters()+track->fNDeadZone<minn) continue;
      if (track->GetNumberOfClusters()+track->fNDeadZone>minn) maxchi2 = track->fChi2MIP[0]; 
      if (track->fChi2MIP[0]>maxchi2) continue;
      minn= track->GetNumberOfClusters()+track->fNDeadZone;
      maxchi2 = track->fChi2MIP[0];
      longtrack=track;
      break;
    }    
    AliITStrackV2 * besttrack = (AliITStrackV2*)array->At(0);
    if (!longtrack) {longtrack = besttrack;}
    else besttrack= longtrack;
    if (besttrack){
      Int_t list[6];
      AliITSclusterV2 * clist[6];
      Float_t shared = GetNumberOfSharedClusters(longtrack,i,list,clist);
      //
      track->fNUsed = shared;      
      track->fNSkipped = besttrack->fNSkipped;
      track->fChi2MIP[0] = besttrack->fChi2MIP[0];
      if (shared>0){
	Int_t nshared;
	Int_t overlist[6]; 
	//
	Int_t sharedtrack = GetOverlapTrack(longtrack, i, nshared, list, overlist);
	//if (sharedtrack==-1) sharedtrack=0;
	if (sharedtrack>=0){       
	  besttrack = GetBest2Tracks(i,sharedtrack,10,5.5);			  
	}
      }   
      if (besttrack&&fConstraint[fPass]) 
      	UpdateESDtrack(besttrack,AliESDtrack::kITSin);
      //if (besttrack&&besttrack->fConstrain) 
      //	UpdateESDtrack(besttrack,AliESDtrack::kITSin);
      if (besttrack->fChi2MIP[0]+besttrack->fNUsed>1.5){
	if ( (TMath::Abs(besttrack->fD[0])>0.1) && fConstraint[fPass]) {
	  track->fReconstructed= kFALSE;
	}
	if ( (TMath::Abs(besttrack->fD[1])>0.1) && fConstraint[fPass]){
	  track->fReconstructed= kFALSE;
	}
      }       

    }    
  }
} 


void AliITStrackerMI::CookLabel(AliITStrackV2 *track,Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t tpcLabel=-1; 
     
  if ( track->fESDtrack)   tpcLabel =  TMath::Abs(track->fESDtrack->GetTPCLabel());

   track->fChi2MIP[9]=0;
   Int_t nwrong=0;
   for (Int_t i=0;i<track->GetNumberOfClusters();i++){
     Int_t cindex = track->GetClusterIndex(i);
     Int_t l=(cindex & 0xf0000000) >> 28;
     AliITSclusterV2 *cl = (AliITSclusterV2*)GetCluster(cindex);
     Int_t isWrong=1;
     for (Int_t ind=0;ind<3;ind++){
       if (tpcLabel>0)
	 if (cl->GetLabel(ind)==tpcLabel) isWrong=0;
     }
     track->fChi2MIP[9]+=isWrong*(2<<l);
     nwrong+=isWrong;
   }
   track->fFakeRatio = double(nwrong)/double(track->GetNumberOfClusters());
   if (tpcLabel>0){
     if (track->fFakeRatio>wrong) track->fLab = -tpcLabel;
     else
       track->fLab = tpcLabel;
   }
   
}



void AliITStrackerMI::CookdEdx(AliITStrackV2* track)
{
  //
  //
  //  Int_t list[6];
  //AliITSclusterV2 * clist[6];
  //  Int_t shared = GetNumberOfSharedClusters(track,index,list,clist);
  Float_t dedx[4];
  Int_t accepted=0;
  track->fChi2MIP[9]=0;
  for (Int_t i=0;i<track->GetNumberOfClusters();i++){
    Int_t cindex = track->GetClusterIndex(i);
    Int_t l=(cindex & 0xf0000000) >> 28;
    AliITSclusterV2 *cl = (AliITSclusterV2*)GetCluster(cindex);
    Int_t lab = TMath::Abs(track->fESDtrack->GetTPCLabel());
    Int_t isWrong=1;
    for (Int_t ind=0;ind<3;ind++){
      if (cl->GetLabel(ind)==lab) isWrong=0;
    }
    track->fChi2MIP[9]+=isWrong*(2<<l);
    if (l<2) continue;
    //if (l>3 && (cl->GetNy()>4) || (cl->GetNz()>4)) continue;  //shared track
    //if (l>3&& !(cl->GetType()==1||cl->GetType()==10)) continue;
    //if (l<4&& !(cl->GetType()==1)) continue;   
    dedx[accepted]= track->fdEdxSample[i];
    //dedx[accepted]= track->fNormQ[l];
    accepted++;
  }
  if (accepted<1) {
    track->SetdEdx(0);
    return;
  }
  Int_t indexes[4];
  TMath::Sort(accepted,dedx,indexes,kFALSE);
  Double_t low=0.;
  Double_t up=0.51;    
  Double_t nl=low*accepted, nu =up*accepted;  
  Float_t sumamp = 0;
  Float_t sumweight =0;
  for (Int_t i=0; i<accepted; i++) {
    Float_t weight =1;
    if (i<nl+0.1)     weight = TMath::Max(1.-(nl-i),0.);
    if (i>nu-1)     weight = TMath::Max(nu-i,0.);
    sumamp+= dedx[indexes[i]]*weight;
    sumweight+=weight;
  }
  track->SetdEdx(sumamp/sumweight);
}


void  AliITStrackerMI::MakeCoeficients(Int_t ntracks){
  //
  //
  fCoeficients = new Float_t[ntracks*48];
  for (Int_t i=0;i<ntracks*48;i++) fCoeficients[i]=-1.;
}


Double_t AliITStrackerMI::GetPredictedChi2MI(AliITStrackV2* track, const AliITSclusterV2 *cluster,Int_t layer) 
{
  //
  //
  //
  Float_t erry,errz;
  Float_t theta = track->GetTgl();
  Float_t phi   = track->GetSnp();
  phi = TMath::Sqrt(phi*phi/(1.-phi*phi));
  GetError(layer,cluster,theta,phi,track->fExpQ,erry,errz);
  Double_t chi2 = track->GetPredictedChi2MI(cluster->GetY(),cluster->GetZ(),erry,errz);
  Float_t ny,nz;
  GetNTeor(layer,cluster, theta,phi,ny,nz);  
  Double_t delta = cluster->GetNy()+cluster->GetNz()-nz-ny;
  if (delta>1){
    chi2+=0.5*TMath::Min(delta/2,2.);
    chi2+=2.*cluster->GetDeltaProbability();
  }
  //
  track->fNy[layer] =ny;
  track->fNz[layer] =nz;
  track->fSigmaY[layer] = erry;
  track->fSigmaZ[layer] = errz;
  //track->fNormQ[layer] = cluster->GetQ()/TMath::Sqrt(1+theta*theta+phi*phi);
  track->fNormQ[layer] = cluster->GetQ()/TMath::Sqrt((1.+ track->fP3*track->fP3)/(1.- track->fP2*track->fP2));
  return chi2;

}

Int_t    AliITStrackerMI::UpdateMI(AliITStrackV2* track, const AliITSclusterV2* cl,Double_t chi2,Int_t index) const 
{
  //
  //
  //
  Int_t layer = (index & 0xf0000000) >> 28;
  track->fClIndex[layer] = index;
  if ( (layer>1) &&track->fNormQ[layer]/track->fExpQ<0.5 ) {
    chi2+= (0.5-track->fNormQ[layer]/track->fExpQ)*10.;
    track->fdEdxMismatch+=(0.5-track->fNormQ[layer]/track->fExpQ)*10.;
  }
  return track->UpdateMI(cl->GetY(),cl->GetZ(),track->fSigmaY[layer],track->fSigmaZ[layer],chi2,index);
}

void AliITStrackerMI::GetNTeor(Int_t layer, const AliITSclusterV2* /*cl*/, Float_t theta, Float_t phi, Float_t &ny, Float_t &nz)
{
  //
  //get "mean shape"
  //
  if (layer==0){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.34;
    return;
  }
  if (layer==1){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.28;
    return;
  }
  
  if (layer>3){
    ny = 2.02+TMath::Abs(phi)*1.95;
    nz = 2.02+TMath::Abs(phi)*2.35;
    return;
  }
  ny  = 6.6-2.7*TMath::Abs(phi);
  nz  = 2.8-3.11*TMath::Abs(phi)+0.45*TMath::Abs(theta);
}



Int_t AliITStrackerMI::GetError(Int_t layer, const AliITSclusterV2*cl, Float_t theta, Float_t phi,Float_t expQ, Float_t &erry, Float_t &errz)
{
  //calculate cluster position error
  //
  Float_t nz,ny;
  GetNTeor(layer, cl,theta,phi,ny,nz);  
  erry   = TMath::Sqrt(cl->GetSigmaY2()); 
  errz   = TMath::Sqrt(cl->GetSigmaZ2()); 
  //
  // PIXELS
  if (layer<2){
    
    if (TMath::Abs(ny-cl->GetNy())>0.6)  {
      if (ny<cl->GetNy()){
	erry*=0.4+TMath::Abs(ny-cl->GetNy());
	errz*=0.4+TMath::Abs(ny-cl->GetNy());
      }else{
	erry*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
	errz*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
      }
    }
    if (TMath::Abs(nz-cl->GetNz())>1.)  {
      erry*=TMath::Abs(nz-cl->GetNz());
      errz*=TMath::Abs(nz-cl->GetNz());	      
    }
    erry*=0.85;
    errz*=0.85;
    erry= TMath::Min(erry,float(0.005));
    errz= TMath::Min(errz,float(0.03));
    return 10;
  }

//STRIPS
  if (layer>3){ 
    if (cl->GetNy()==100||cl->GetNz()==100){
      erry = 0.004;
      errz = 0.2;
      return 100;
    }
    if (cl->GetNy()+cl->GetNz()>12){
      erry = 0.06;
      errz = 0.57;
      return 100;
    }
    Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
    Float_t chargematch = TMath::Max(double(normq/expQ),2.);
    //
    if (cl->GetType()==1 || cl->GetType()==10 ){     							       
      if (chargematch<1.0 || (cl->GetNy()+cl->GetNz()<nz+ny+0.5)){
	errz = 0.043;
	erry = 0.00094;
	return 101;
      }
      if (cl->GetNy()+cl->GetNz()<nz+ny+1.2){
	errz = 0.06;
	erry =0.0013;
	return 102;
      }
      erry = 0.0027;
      errz = TMath::Min(0.028*(chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.15);
      return 103;
    }
    if (cl->GetType()==2 || cl->GetType()==11 ){ 
      erry = TMath::Min(0.0010*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.05);
      errz = TMath::Min(0.025*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.5);
      return 104;
    }
    
    if (cl->GetType()>100 ){     							       
      if ((chargematch+cl->GetNy()+cl->GetNz()-nz-ny<1.5)){
	errz = 0.05;
	erry = 0.00096;
	return 105;
      }
      if (cl->GetNy()+cl->GetNz()-nz-ny<1){
	errz = 0.10;
	erry = 0.0025;
	return 106;
      }

      errz = TMath::Min(0.05*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.4);
      erry = TMath::Min(0.003*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.05);
      return 107;
    }    
    Float_t diff = cl->GetNy()+cl->GetNz()-ny-nz;
    if (diff<1) diff=1;
    if (diff>4) diff=4;
        
    if (cl->GetType()==5||cl->GetType()==6||cl->GetType()==7||cl->GetType()==8){
      errz = 0.14*diff;
      erry = 0.003*diff;
      return 108;
    }  
    erry = 0.04*diff;
    errz = 0.06*diff;
    return 109;
  }
  //DRIFTS
  Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
  Float_t chargematch = normq/expQ;
  Float_t factorz=1;
  Int_t   cnz = cl->GetNz()%10;
  //charge match
  if (cl->GetType()==1){
    if (chargematch<1.25){
      erry =  0.0028*(1.+6./cl->GetQ());  // gold clusters
    }
    else{
      erry = 0.003*chargematch;
      if (cl->GetNz()==3) erry*=1.5;
    }
    if (chargematch<1.0){
      errz =  0.0011*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.002*(1+2*(chargematch-1.));
    }
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }
  if (cl->GetType()>1){
    if (chargematch<1){
      erry =  0.00385*(1.+6./cl->GetQ());  // gold clusters
      errz =  0.0016*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.0014*(1+3*(chargematch-1.));
      erry = 0.003*(1+3*(chargematch-1.));
    } 
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }

  if (TMath::Abs(cl->GetY())>2.5){
    factorz*=1+2*(TMath::Abs(cl->GetY())-2.5);
  }
  if (TMath::Abs(cl->GetY())<1){
    factorz*=1.+0.5*TMath::Abs(TMath::Abs(cl->GetY())-1.);
  }
  factorz= TMath::Min(factorz,float(4.));  
  errz*=factorz;

  erry= TMath::Min(erry,float(0.05));
  errz= TMath::Min(errz,float(0.05));  
  return 200;
}



void AliITStrackerMI::SignDeltas( TObjArray *ClusterArray, Float_t vz)
{
  //
  //  
  Int_t entries = ClusterArray->GetEntriesFast();
  if (entries<4) return;
  AliITSclusterV2* cluster = (AliITSclusterV2*)ClusterArray->At(0);
  Int_t layer = cluster->GetLayer();
  if (layer>1) return;
  Int_t index[10000];
  Int_t ncandidates=0;
  Float_t r = (layer>0)? 7:4;
  // 
  for (Int_t i=0;i<entries;i++){
    AliITSclusterV2* cl0 = (AliITSclusterV2*)ClusterArray->At(i);
    Float_t nz = 1+TMath::Abs((cl0->GetZ()-vz)/r);
    if (cl0->GetNy()+cl0->GetNz()<=5+2*layer+nz) continue;
    index[ncandidates] = i;  //candidate to belong to delta electron track
    ncandidates++;
    if (cl0->GetNy()+cl0->GetNz()>9+2*layer+nz) {
      cl0->SetDeltaProbability(1);
    }
  }
  //
  //  
  //
  for (Int_t i=0;i<ncandidates;i++){
    AliITSclusterV2* cl0 = (AliITSclusterV2*)ClusterArray->At(index[i]);
    if (cl0->GetDeltaProbability()>0.8) continue;
    // 
    Int_t ncl = 0;
    Float_t y[100],z[100],sumy,sumz,sumy2, sumyz, sumw;
    sumy=sumz=sumy2=sumyz=sumw=0.0;
    for (Int_t j=0;j<ncandidates;j++){
      if (i==j) continue;
      AliITSclusterV2* cl1 = (AliITSclusterV2*)ClusterArray->At(index[j]);
      //
      Float_t dz = cl0->GetZ()-cl1->GetZ();
      Float_t dy = cl0->GetY()-cl1->GetY();
      if (TMath::Sqrt(dz*dz+dy*dy)<0.2){
	Float_t weight = cl1->GetNy()+cl1->GetNz()-2;
	y[ncl] = cl1->GetY();
	z[ncl] = cl1->GetZ();
	sumy+= y[ncl]*weight;
	sumz+= z[ncl]*weight;
	sumy2+=y[ncl]*y[ncl]*weight;
	sumyz+=y[ncl]*z[ncl]*weight;
	sumw+=weight;
	ncl++;
      }
    }
    if (ncl<4) continue;
    Float_t det = sumw*sumy2  - sumy*sumy;
    Float_t delta=1000;
    if (TMath::Abs(det)>0.01){
      Float_t z0  = (sumy2*sumz - sumy*sumyz)/det;
      Float_t k   = (sumyz*sumw - sumy*sumz)/det;
      delta = TMath::Abs(cl0->GetZ()-(z0+k*cl0->GetY()));
    }
    else{
      Float_t z0  = sumyz/sumy;
      delta = TMath::Abs(cl0->GetZ()-z0);
    }
    if ( delta<0.05) {
      cl0->SetDeltaProbability(1-20.*delta);
    }   
  }
}


void AliITStrackerMI::UpdateESDtrack(AliITStrackV2* track, ULong_t flags) const
{
  //
  //
  track->UpdateESDtrack(flags);
  AliITStrackV2 * oldtrack = (AliITStrackV2*)(track->fESDtrack->GetITStrack());
  if (oldtrack) delete oldtrack; 
  track->fESDtrack->SetITStrack(new AliITStrackV2(*track));
}




void  AliITStrackerMI::FindV0(AliESD */*event*/)
{
  //
  // fast V0 finder
  //
  AliHelix helixes[30000];
  TObjArray trackarray(30000);
  Float_t dist[30000];
  AliITSRecV0Info vertex;
  //
  Int_t entries = fTrackHypothesys.GetEntriesFast();
  for (Int_t i=0;i<entries;i++){
    TObjArray * array = (TObjArray*)fTrackHypothesys.At(i);
    if (!array) continue;
    AliITStrackV2 * track = (AliITStrackV2*)array->At(fBestTrackIndex[i]);
    if (track){
      dist[i] = TMath::Sqrt(track->fD[0]*track->fD[0]+track->fD[1]*track->fD[1]);
      trackarray.AddAt(track,i);
      new (&helixes[i]) AliHelix(*track);
    }
  }
  for (Int_t itrack0=0;itrack0<entries;itrack0++){
    AliITStrackV2 * track0 = (AliITStrackV2*)trackarray.At(itrack0);
    if (!track0) continue;
    if (dist[itrack0]<0.2) continue;
    for (Int_t itrack1=itrack0+1;itrack1<entries;itrack1++){
      AliITStrackV2 * track1 = (AliITStrackV2*)trackarray.At(itrack1);
      if (!track1) continue;
      if (dist[itrack1]<0.2) continue;
      if (track1->fP4*track0->fP4>0) continue; //the same sign
      AliHelix *h1 = &helixes[itrack0];
      AliHelix *h2 = &helixes[itrack1];
      
      Double_t distance = TestV0(h1,h2,&vertex);
      if (distance>0.4) continue;
      if (vertex.fRr<0.3) continue;
      if (vertex.fRr>27) continue;
      Float_t v[3]={GetX(),GetY(),GetZ()};
      vertex.Update(v,track0->fMass, track1->fMass);
      
      if ((TMath::Abs(vertex.fInvMass-0.4976)<0.05 || TMath::Abs(vertex.fInvMass-0.0)<0.1||  TMath::Abs(vertex.fInvMass-1.1)<0.1)) {
	if (vertex.fPointAngle<0.85) continue;	
      }
      else{
	if (vertex.fPointAngle<0.98) continue;	
      }
      if (TMath::Abs(TMath::Abs(track0->fLab)-TMath::Abs(track1->fLab))<2){
	Float_t mindist = FindBestPair(itrack0,itrack1,&vertex);
	if (mindist>1) FindBestPair(itrack0,itrack1,&vertex);
      }
    }    
  }
}

Double_t  AliITStrackerMI::FindBestPair(Int_t esdtrack0, Int_t esdtrack1,AliITSRecV0Info *vertex)
{
  //
  // try to find best pair from the tree of track hyp.
  //
  TObjArray * array0 = (TObjArray*)fTrackHypothesys.At(esdtrack0);
  Int_t entries0 = array0->GetEntriesFast();
  TObjArray * array1 = (TObjArray*)fTrackHypothesys.At(esdtrack1);
  Int_t entries1 = array1->GetEntriesFast();  
  //
  //
  AliITSRecV0Info v0;
  AliITSRecV0Info vbest;
  Double_t criticalradius = vertex->fRr;
  Double_t mindist =1000;
  Int_t bestpair[2];
  //
  for (Int_t itrack0=0;itrack0<entries0;itrack0++){
    AliITStrackV2 * track0 = (AliITStrackV2*)array0->At(itrack0);
    if (!track0) continue;
    if (track0->fX<criticalradius-1) continue;
    if (track0->fX>criticalradius+5) continue;
    for (Int_t itrack1=0;itrack1<entries1;itrack1++){
      AliITStrackV2 * track1 = (AliITStrackV2*)array1->At(itrack1);
      if (!track1) continue;      
      if (track1->fX<criticalradius-1) continue;
      if (track1->fX>criticalradius+5) continue;

      AliHelix h0(*track0);
      AliHelix h1(*track1);
      Double_t distance = TestV0(&h0,&h1,&v0);
      if (v0.fRr>30) continue;
      if (v0.fRr<0.2) continue;
      Float_t v[3]={GetX(),GetY(),GetZ()};
      v0.Update(v,track0->fMass, track1->fMass);
      if (distance<mindist){
	mindist = distance;
	bestpair[0]=itrack0;
	bestpair[1]=itrack1;
	vbest = v0;	
      }      
    }
  }
  if (mindist<0.2){
    vbest.Dump();
    vertex->Dump();
  }
  return mindist;
}


void  AliITSRecV0Info::Update(Float_t vertex[3], Float_t mass1, Float_t mass2)
{
  //
  //Update V0 information
  //
  Float_t v[3] = {fXr[0]-vertex[0],fXr[1]-vertex[1],fXr[2]-vertex[2]};
  Float_t p[3] = {fPdr[0]+fPm[0], fPdr[1]+fPm[1],fPdr[2]+fPm[2]};

  Float_t vnorm2 = v[0]*v[0]+v[1]*v[1];
  Float_t vnorm3 = TMath::Sqrt(v[2]*v[2]+vnorm2);
  vnorm2 = TMath::Sqrt(vnorm2);
  Float_t pnorm2 = p[0]*p[0]+p[1]*p[1];
  Float_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
  pnorm2 = TMath::Sqrt(pnorm2);
  
  fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
  fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
  fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  
  Float_t e1    = TMath::Sqrt(mass1*mass1+
			      fPdr[0]*fPdr[0]+
			      fPdr[1]*fPdr[1]+
			      fPdr[2]*fPdr[2]);
  Float_t e2    = TMath::Sqrt(mass2*mass2+
			      fPm[0]*fPm[0]+
			      fPm[1]*fPm[1]+
			      fPm[2]*fPm[2]);
  
  fInvMass =  
    (fPm[0]+fPdr[0])*(fPm[0]+fPdr[0])+
    (fPm[1]+fPdr[1])*(fPm[1]+fPdr[1])+
    (fPm[2]+fPdr[2])*(fPm[2]+fPdr[2]);
  
  fInvMass = TMath::Sqrt((e1+e2)*(e1+e2)-fInvMass);


}



Double_t   AliITStrackerMI::TestV0(AliHelix *helix1, AliHelix *helix2, AliITSRecV0Info *vertex)
{
  //
  // test the helixes for the distnce calculate vertex characteristic
  //
  Float_t distance1,distance2;
  AliHelix & dhelix1 = *helix1;
  Double_t pp[3],xx[3];
  dhelix1.GetMomentum(0,pp,0);
  dhelix1.Evaluate(0,xx);      
  AliHelix &mhelix = *helix2;    
  //
  //find intersection linear
  //
  Double_t phase[2][2],radius[2];
  Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
  Double_t delta1=10000,delta2=10000;  
  
  if (points>0){
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
  }
  if (points==2){    
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
  }
  distance1 = TMath::Min(delta1,delta2);
  vertex->fDist1 = TMath::Sqrt(distance1);
  //
  //find intersection parabolic
  //
  points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
  delta1=10000,delta2=10000;  
  
  if (points>0){
    dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
  }
  if (points==2){    
    dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
  }
  
  distance2 = TMath::Min(delta1,delta2);
  vertex->fDist2 = TMath::Sqrt(distance2);      
  if (distance2<100){
    if (delta1<delta2){
      //get V0 info
      dhelix1.Evaluate(phase[0][0],vertex->fXr);
      dhelix1.GetMomentum(phase[0][0],vertex->fPdr);
      mhelix.GetMomentum(phase[0][1],vertex->fPm);
      dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],vertex->fAngle);
	vertex->fRr = TMath::Sqrt(radius[0]);
    }
    else{
      dhelix1.Evaluate(phase[1][0],vertex->fXr);
      dhelix1.GetMomentum(phase[1][0],vertex->fPdr);
      mhelix.GetMomentum(phase[1][1],vertex->fPm);
      dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],vertex->fAngle);
      vertex->fRr = TMath::Sqrt(radius[1]);
    }
  }
  //            
  //
  return  vertex->fDist2;
}

