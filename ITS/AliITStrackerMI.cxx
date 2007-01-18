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

//-------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    It reads AliITSRecPoint clusters and creates AliITStrackMI tracks
//                   and fills with them the ESD
//          Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch 
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//     
//-------------------------------------------------------------------------

#include <TMatrixD.h>
#include <TTree.h>
#include <TTreeStream.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliV0.h"
#include "AliHelix.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliITStrackerMI.h"
#include "AliTrackPointArray.h"
#include "AliAlignObj.h"

ClassImp(AliITStrackerMI)



AliITStrackerMI::AliITSlayer AliITStrackerMI::fgLayers[kMaxLayer]; // ITS layers
AliITStrackerMI::AliITStrackerMI():AliTracker(),
fI(0),
fBestTrack(),
fTrackToFollow(),
fTrackHypothesys(),
fBestHypothesys(),
fOriginal(),
fCurrentEsdTrack(),
fPass(0),
fAfterV0(kFALSE),
fLastLayerToTrackTo(0),
fCoeficients(0),
fEsd(0),
fDebugStreamer(0){
  //Default constructor
}


AliITStrackerMI::AliITStrackerMI(const AliITSgeom *geom) : AliTracker(),
fI(kMaxLayer),
fBestTrack(),
fTrackToFollow(),
fTrackHypothesys(),
fBestHypothesys(),
fOriginal(),
fCurrentEsdTrack(),
fPass(0),
fAfterV0(kFALSE),
fLastLayerToTrackTo(kLastLayerToTrackTo),
fCoeficients(0),
fEsd(0),
fDebugStreamer(0){
  //--------------------------------------------------------------------
  //This is the AliITStrackerMI constructor
  //--------------------------------------------------------------------
  fCoeficients = 0;
  fAfterV0     = kFALSE;
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
  for (Int_t i=0;i<100000;i++){
    fBestTrackIndex[i]=0;
  }
  //
  fDebugStreamer = new TTreeSRedirector("ITSdebug.root");

}

AliITStrackerMI::AliITStrackerMI(const AliITStrackerMI &tracker):AliTracker(tracker),
fI(tracker.fI),
fBestTrack(tracker.fBestTrack),
fTrackToFollow(tracker.fTrackToFollow),
fTrackHypothesys(tracker.fTrackHypothesys),
fBestHypothesys(tracker.fBestHypothesys),
fOriginal(tracker.fOriginal),
fCurrentEsdTrack(tracker.fCurrentEsdTrack),
fPass(tracker.fPass),
fAfterV0(tracker.fAfterV0),
fLastLayerToTrackTo(tracker.fLastLayerToTrackTo),
fCoeficients(tracker.fCoeficients),
fEsd(tracker.fEsd),
fDebugStreamer(tracker.fDebugStreamer){
  //Copy constructor
}

AliITStrackerMI & AliITStrackerMI::operator=(const AliITStrackerMI &tracker){
  //Assignment operator
  this->~AliITStrackerMI();
  new(this) AliITStrackerMI(tracker);
  return *this;
}


AliITStrackerMI::~AliITStrackerMI()
{
  //
  //destructor
  //
  if (fCoeficients) delete []fCoeficients;
  if (fDebugStreamer) {
    //fDebugStreamer->Close();
    delete fDebugStreamer;
  }
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
  TBranch *branch=cTree->GetBranch("ITSRecPoints");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch !\n");
    return 1;
  }

  TClonesArray dummy("AliITSRecPoint",10000), *clusters=&dummy;
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
        AliITSRecPoint *c=(AliITSRecPoint*)clusters->UncheckedAt(ncl);
	detector = c->GetDetectorIndex();
        fgLayers[i].InsertCluster(new AliITSRecPoint(*c));
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
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));
	  hit[1]=-7.05;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));
	  hit[1]=-7.15;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));
	  hit[1] =0.06;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));
	  hit[1]=7.05;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));
	  hit[1]=7.25;
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<2.) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab, hit, info));       
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

static Int_t CorrectForDeadZoneMaterial(AliITStrackMI *t) {
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
  fOriginal.Clear();
  fEsd = event;         // store pointer to the esd 
  {/* Read ESD tracks */
    Int_t nentr=event->GetNumberOfTracks();
    Info("Clusters2Tracks", "Number of ESD tracks: %d\n", nentr);
    while (nentr--) {
      AliESDtrack *esd=event->GetTrack(nentr);

      if ((esd->GetStatus()&AliESDtrack::kTPCin)==0) continue;
      if (esd->GetStatus()&AliESDtrack::kTPCout) continue;
      if (esd->GetStatus()&AliESDtrack::kITSin) continue;
      if (esd->GetKinkIndex(0)>0) continue;   //kink daughter
      AliITStrackMI *t=0;
      try {
        t=new AliITStrackMI(*esd);
      } catch (const Char_t *msg) {
        //Warning("Clusters2Tracks",msg);
        delete t;
        continue;
      }
      //t->fD[0] = t->GetD(GetX(),GetY());
      //t->fD[1] = t->GetZat(GetX())-GetZ();
      t->GetDZ(GetX(),GetY(),GetZ(),t->GetDP());              //I.B.
      Double_t vdist = TMath::Sqrt(t->GetD(0)*t->GetD(0)+t->GetD(1)*t->GetD(1));
      if (t->GetMass()<0.13) t->SetMass(0.13957); // MI look to the esd - mass hypothesys  !!!!!!!!!!!
      // write expected q
      t->SetExpQ(TMath::Max(0.8*t->GetESDtrack()->GetTPCsignal(),30.));

      if (esd->GetV0Index(0)>0 && t->GetD(0)<30){
	//track - can be  V0 according to TPC
      }
      else{	
	if (TMath::Abs(t->GetD(0))>10) {
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
	  //Warning("Clusters2Tracks",
	  //        "failed to correct for the material in the dead zone !\n");
	  delete t;
	  continue;
	}
      }
      t->SetReconstructed(kFALSE);
      itsTracks.AddLast(t);
      fOriginal.AddLast(t);
    }
  } /* End Read ESD tracks */

  itsTracks.Sort();
  fOriginal.Sort();
  Int_t nentr=itsTracks.GetEntriesFast();
  fTrackHypothesys.Expand(nentr);
  fBestHypothesys.Expand(nentr);
  MakeCoeficients(nentr);
  Int_t ntrk=0;
  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
//       cerr<<fPass<<"    "<<i<<'\r';
       fCurrentEsdTrack = i;
       AliITStrackMI *t=(AliITStrackMI*)itsTracks.UncheckedAt(i);
       if (t==0) continue;              //this track has been already tracked
       if (t->GetReconstructed()&&(t->GetNUsed()<1.5)) continue;  //this track was  already  "succesfully" reconstructed
       //if ( (TMath::Abs(t->GetD(GetX(),GetY()))  >3.) && fConstraint[fPass]) continue;
       //if ( (TMath::Abs(t->GetZat(GetX())-GetZ())>3.) && fConstraint[fPass]) continue;
       Float_t dz[2]; t->GetDZ(GetX(),GetY(),GetZ(),dz);              //I.B.
       if ( (TMath::Abs(dz[0])>3.) && fConstraint[fPass]) continue;
       if ( (TMath::Abs(dz[1])>3.) && fConstraint[fPass]) continue;

       Int_t tpcLabel=t->GetLabel(); //save the TPC track label       
       fI = 6;
       ResetTrackToFollow(*t);
       ResetBestTrack();
       FollowProlongationTree(t,i,fConstraint[fPass]);

       SortTrackHypothesys(fCurrentEsdTrack,20,0);  //MI change
       //
       AliITStrackMI * besttrack = GetBestHypothesys(fCurrentEsdTrack,t,15);
       if (!besttrack) continue;
       besttrack->SetLabel(tpcLabel);
       //       besttrack->CookdEdx();
       CookdEdx(besttrack);
       besttrack->SetFakeRatio(1.);
       CookLabel(besttrack,0.); //For comparison only
       UpdateESDtrack(besttrack,AliESDtrack::kITSin);

       /*       
       if ( besttrack->GetNumberOfClusters()<6 && fConstraint[fPass]) {	 
	 continue;
       }
       if (besttrack->fChi2MIP[0]+besttrack->fNUsed>3.5) continue;
       if ( (TMath::Abs(besttrack->fD[0]*besttrack->fD[0]+besttrack->fD[1]*besttrack->fD[1])>0.1) && fConstraint[fPass])  continue;	 
       //delete itsTracks.RemoveAt(i);
       */
       if (fConstraint[fPass]&&(!besttrack->IsGoldPrimary())) continue;  //to be tracked also without vertex constrain 


       t->SetReconstructed(kTRUE);
       ntrk++;                     
     }
     GetBestHypothesysMIP(itsTracks); 
  }

  //GetBestHypothesysMIP(itsTracks);
  UpdateTPCV0(event);
  FindV02(event);
  fAfterV0 = kTRUE;
  //GetBestHypothesysMIP(itsTracks);
  //
  itsTracks.Delete();
  //
  Int_t entries = fTrackHypothesys.GetEntriesFast();
  for (Int_t ientry=0;ientry<entries;ientry++){
    TObjArray * array =(TObjArray*)fTrackHypothesys.UncheckedAt(ientry);
    if (array) array->Delete();
    delete fTrackHypothesys.RemoveAt(ientry); 
  }

  fTrackHypothesys.Delete();
  fBestHypothesys.Delete();
  fOriginal.Clear();
  delete []fCoeficients;
  fCoeficients=0;
  Info("Clusters2Tracks","Number of prolonged tracks: %d\n",ntrk);
  
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

     AliITStrackMI *t=0;
     try {
        t=new AliITStrackMI(*esd);
     } catch (const Char_t *msg) {
       //Warning("PropagateBack",msg);
        delete t;
        continue;
     }
     t->SetExpQ(TMath::Max(0.8*t->GetESDtrack()->GetTPCsignal(),30.));

     ResetTrackToFollow(*t);

     // propagete to vertex [SR, GSI 17.02.2003]
     // Start Time measurement [SR, GSI 17.02.2003], corrected by I.Belikov
     if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {
       if (fTrackToFollow.PropagateToVertex(event->GetVertex())) {
          fTrackToFollow.StartTimeIntegral();
       }
       fTrackToFollow.PropagateTo(3.,-0.0028,65.19);
     }

     fTrackToFollow.ResetCovariance(10.); fTrackToFollow.ResetClusters();
     if (RefitAt(49.,&fTrackToFollow,t)) {
        if (CorrectForDeadZoneMaterial(&fTrackToFollow)!=0) {
          //Warning("PropagateBack",
          //        "failed to correct for the material in the dead zone !\n");
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
  RefitV02(event);
  Int_t nentr=event->GetNumberOfTracks();
  Info("RefitInward", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);

    if ((esd->GetStatus()&AliESDtrack::kITSout) == 0) continue;
    if (esd->GetStatus()&AliESDtrack::kITSrefit) continue;
    if (esd->GetStatus()&AliESDtrack::kTPCout)
      if ((esd->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;

    AliITStrackMI *t=0;
    try {
        t=new AliITStrackMI(*esd);
    } catch (const Char_t *msg) {
      //Warning("RefitInward",msg);
        delete t;
        continue;
    }
    t->SetExpQ(TMath::Max(0.8*t->GetESDtrack()->GetTPCsignal(),30.));
    if (CorrectForDeadZoneMaterial(t)!=0) {
      //Warning("RefitInward",
      //         "failed to correct for the material in the dead zone !\n");
       delete t;
       continue;
    }

    ResetTrackToFollow(*t);
    fTrackToFollow.ResetClusters();

    if ((esd->GetStatus()&AliESDtrack::kTPCin)==0)
      fTrackToFollow.ResetCovariance(10.);

    //Refitting...
    if (RefitAt(3.7, &fTrackToFollow, t,kTRUE)) {
       fTrackToFollow.SetLabel(t->GetLabel());
       //       fTrackToFollow.CookdEdx();
       CookdEdx(&fTrackToFollow);

       CookLabel(&fTrackToFollow,0.0); //For comparison only

       if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {//The beam pipe    
         AliESDtrack  *esdTrack =fTrackToFollow.GetESDtrack();
         esdTrack->UpdateTrackParams(&fTrackToFollow,AliESDtrack::kITSrefit);
         Float_t r[3]={0.,0.,0.};
         Double_t maxD=3.;
	 esdTrack->RelateToVertex(event->GetVertex(),GetBz(r),maxD);
         ntrk++;
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

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoPhysicalNode.h"

Bool_t AliITStrackerMI::GetTrackPoint(Int_t index, AliTrackPoint& p) const {
  //
  // Get track space point with index i
  //
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  AliITSRecPoint *cl = fgLayers[l].GetCluster(c);
  Int_t idet = cl->GetDetectorIndex();

  const char* name = AliAlignObj::SymName((AliAlignObj::ELayerID)
					  (l+AliAlignObj::kFirstLayer), idet);
  TGeoPNEntry *mapPN = gGeoManager->GetAlignableEntry(name);

  if (!mapPN) return kFALSE;
  TGeoPhysicalNode *node = mapPN->GetPhysicalNode();
  if (!node) {
    gGeoManager->MakeAlignablePN(name);
    node = mapPN->GetPhysicalNode();
  }
  if (!node) return kFALSE;
  TGeoHMatrix* matrix = node->GetMatrix();
  if (!matrix) return kFALSE;

  //
  // Calculate the global coordinates
  //
  Double_t localCoord[]  = {cl->GetDetLocalX(), 0, cl->GetDetLocalZ()};
  // LG  AliAlignObj makes life simple but keep in mind that alignable
  // LG  volume doesn't mean sensitive volume. There might be a shift
  // LG  between the 2, has it is here for the SPD :
  if (l<2) localCoord[1] = 0.01;
  // LG   !!!   Check for this when the new geometry comes   !!!
  Double_t globalCoord[3] = {0};
  matrix->LocalToMaster(localCoord, globalCoord);
  Float_t xyz[3]= {globalCoord[0], globalCoord[1], globalCoord[2]};

  //
  // Calculate the cov matrix
  //
  TGeoRotation rotMatrix(*matrix);
  TGeoRotation rotMatrixTr(rotMatrix.Inverse());
  Double_t sigmaArray[] = {cl->GetSigmaY2(),0,0, 0,0,0, 0,0,cl->GetSigmaZ2()};
  TGeoRotation sigmMatrix;
  sigmMatrix.SetMatrix( sigmaArray );
  sigmMatrix.MultiplyBy(&rotMatrixTr, kFALSE);
  sigmMatrix.MultiplyBy(&rotMatrix, kTRUE);
  const Double_t *globalSigma =  sigmMatrix.GetRotationMatrix();
  Float_t cov[6]= { globalSigma[0], globalSigma[1], globalSigma[2],
		                    globalSigma[4], globalSigma[5],
		                                    globalSigma[8] };

  p.SetXYZ(xyz[0],xyz[1],xyz[2],cov);
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer; 
  switch (l) {
  case 0:
    iLayer = AliAlignObj::kSPD1;
    break;
  case 1:
    iLayer = AliAlignObj::kSPD2;
    break;
  case 2:
    iLayer = AliAlignObj::kSDD1;
    break;
  case 3:
    iLayer = AliAlignObj::kSDD2;
    break;
  case 4:
    iLayer = AliAlignObj::kSSD1;
    break;
  case 5:
    iLayer = AliAlignObj::kSSD2;
    break;
  default:
    AliWarning(Form("Wrong layer index in ITS (%d) !",l));
    break;
  };
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,idet);
  p.SetVolumeID((UShort_t)volid);
  return kTRUE;
}

void AliITStrackerMI::FollowProlongationTree(AliITStrackMI * otrack, Int_t esdindex, Bool_t constrain) 
{
  //--------------------------------------------------------------------
  // Follow prolongation tree
  //--------------------------------------------------------------------
  //
  AliESDtrack * esd = otrack->GetESDtrack();
  if (esd->GetV0Index(0)>0){
    //
    // TEMPORARY SOLLUTION: map V0 indexes to point to proper track
    //                      mapping of esd track is different as its track in Containers
    //                      Need something more stable
    //                      Indexes are set back againg to the ESD track indexes in UpdateTPCV0
    for (Int_t i=0;i<3;i++){
      Int_t  index = esd->GetV0Index(i);
      if (index==0) break;
      AliESDv0 * vertex = fEsd->GetV0(index);
      if (vertex->GetStatus()<0) continue;     // rejected V0
      //
            if (esd->GetSign()>0) {
              vertex->SetIndex(0,esdindex);
            }
            else{
              vertex->SetIndex(1,esdindex);
            }
    }
  }
  TObjArray *bestarray = (TObjArray*)fBestHypothesys.At(esdindex);
  if (!bestarray){
    bestarray = new TObjArray(5);
    fBestHypothesys.AddAt(bestarray,esdindex);
  }

  //
  //setup tree of the prolongations
  //
  static AliITStrackMI tracks[7][100];
  AliITStrackMI *currenttrack;
  static AliITStrackMI currenttrack1;
  static AliITStrackMI currenttrack2;  
  static AliITStrackMI backuptrack;
  Int_t ntracks[7];
  Int_t nindexes[7][100];
  Float_t normalizedchi2[100];
  for (Int_t ilayer=0;ilayer<6;ilayer++) ntracks[ilayer]=0;
  otrack->SetNSkipped(0);
  new (&(tracks[6][0])) AliITStrackMI(*otrack);
  ntracks[6]=1;
  for (Int_t i=0;i<7;i++) nindexes[i][0]=0;
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
      if (tracks[ilayer+1][nindexes[ilayer+1][itrack]].GetNSkipped()>0) nskipped++;
      if (tracks[ilayer+1][nindexes[ilayer+1][itrack]].GetNUsed()>2.) nused++;
      if (ntracks[ilayer]>15+ilayer){
	if (itrack>1&&tracks[ilayer+1][nindexes[ilayer+1][itrack]].GetNSkipped()>0 && nskipped>4+ilayer) continue;
	if (itrack>1&&tracks[ilayer+1][nindexes[ilayer+1][itrack]].GetNUsed()>2. && nused>3) continue;
      }

      new(&currenttrack1)  AliITStrackMI(tracks[ilayer+1][nindexes[ilayer+1][itrack]]);
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
      new(&currenttrack2)  AliITStrackMI(currenttrack1);
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
      if (constrain){
	msy/=60; msz/=60.;
      }
      else{
	msy/=50; msz/=50.;
      }
      //
      const AliITSRecPoint *c=0; Int_t ci=-1;
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
	  new (&backuptrack) AliITStrackMI(currenttrack2);
	  change = kTRUE;
	  currenttrack =&currenttrack2;
	  if (!currenttrack->Propagate(det.GetPhi(),det.GetR())) {
	    new (currenttrack) AliITStrackMI(backuptrack);
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
	  AliITStrackMI * updatetrack = new (&tracks[ilayer][ntracks[ilayer]]) AliITStrackMI(*currenttrack);
	  updatetrack->SetClIndex(ilayer,0);
	  if (change){
	    new (&currenttrack2) AliITStrackMI(backuptrack);
	  }
	  if (c->GetQ()!=0){
	    if (!UpdateMI(updatetrack,c,chi2,(ilayer<<28)+ci)) continue; 
	    updatetrack->SetSampledEdx(c->GetQ(),updatetrack->GetNumberOfClusters()-1); //b.b.
	  }
	  else {
	    updatetrack->SetNDeadZone(updatetrack->GetNDeadZone()+1);
	    updatetrack->SetDeadZoneProbability(GetDeadZoneProbability(updatetrack->GetZ(),TMath::Sqrt(updatetrack->GetSigmaZ2())));
	  }
	  if (c->IsUsed()){
	    updatetrack->IncrementNUsed();
	  }
	  Double_t x0;
	  Double_t d=layer.GetThickness(updatetrack->GetY(),updatetrack->GetZ(),x0);
	  updatetrack->CorrectForMaterial(d,x0);	  
	  if (constrain) {
	    updatetrack->SetConstrain(constrain);
	    fI = ilayer;
	    Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
	    Double_t xyz[]={GetX(),GetY(),GetZ()};
	    Double_t ptfactor = 1;
	    Double_t ers[]={GetSigmaX()*ptfactor,GetSigmaY()*ptfactor,GetSigmaZ()};
	    Bool_t isPrim = kTRUE;
	    if (ilayer<4){
	      //updatetrack->fD[0] = updatetrack->GetD(GetX(),GetY());
	      //updatetrack->fD[1] = updatetrack->GetZat(GetX())-GetZ();
              updatetrack->GetDZ(GetX(),GetY(),GetZ(),updatetrack->GetDP()); //I.B.
	      if ( TMath::Abs(updatetrack->GetD(0)/(1.+ilayer))>0.4 ||  TMath::Abs(updatetrack->GetD(1)/(1.+ilayer))>0.4) isPrim=kFALSE;
	    }
	    if (isPrim) updatetrack->Improve(d,xyz,ers);
	  } //apply vertex constrain	  	  
	  ntracks[ilayer]++;
	}  // create new hypothesy 
      } // loop over possible cluster prolongation      
      //      if (constrain&&itrack<2&&currenttrack1.fNSkipped==0 && deadzone==0){	
      if (constrain&&itrack<2&&currenttrack1.GetNSkipped()==0 && deadzone==0&&ntracks[ilayer]<100){	
	AliITStrackMI* vtrack = new (&tracks[ilayer][ntracks[ilayer]]) AliITStrackMI(currenttrack1);
	vtrack->SetClIndex(ilayer,0);
	fI = ilayer;
	Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
	Double_t xyz[]={GetX(),GetY(),GetZ()};
	Double_t ers[]={GetSigmaX(),GetSigmaY(),GetSigmaZ()};
	vtrack->Improve(d,xyz,ers);
	vtrack->IncrementNSkipped();
	ntracks[ilayer]++;
      }

      if (constrain&&itrack<1&&TMath::Abs(currenttrack1.GetTgl())>1.1){  //big theta -- for low mult. runs
	AliITStrackMI* vtrack = new (&tracks[ilayer][ntracks[ilayer]]) AliITStrackMI(currenttrack1);
	vtrack->SetClIndex(ilayer,0);
	fI = ilayer;
	Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
	Double_t xyz[]={GetX(),GetY(),GetZ()};
	Double_t ers[]={GetSigmaX(),GetSigmaY(),GetSigmaZ()};
	vtrack->Improve(d,xyz,ers);
	vtrack->SetNDeadZone(vtrack->GetNDeadZone()+1);
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
	if ( constrain && normalizedchi2[itrack]<kMaxNormChi2C[ilayer]+1) accepted++;
	if (!constrain && normalizedchi2[itrack]<kMaxNormChi2NonC[ilayer]+1) accepted++;
      }
    }
    TMath::Sort(ntracks[ilayer],normalizedchi2,nindexes[ilayer],kFALSE);
    ntracks[ilayer] = TMath::Min(accepted,7+2*ilayer);
    if (ntracks[ilayer]<golds+2+ilayer) ntracks[ilayer]=TMath::Min(golds+2+ilayer,accepted);
    if (ntracks[ilayer]>90) ntracks[ilayer]=90; 
  } //loop over layers
  //printf("%d\t%d\t%d\t%d\t%d\t%d\n",ntracks[0],ntracks[1],ntracks[2],ntracks[3],ntracks[4],ntracks[5]);
  Int_t max = constrain? 20: 5;

  for (Int_t i=0;i<TMath::Min(max,ntracks[0]);i++) {
    AliITStrackMI & track= tracks[0][nindexes[0][i]];
    if (track.GetNumberOfClusters()<2) continue;
    if (!constrain&&track.GetNormChi2(0)>7.)continue;
    AddTrackHypothesys(new AliITStrackMI(track), esdindex);
  }
  for (Int_t i=0;i<TMath::Min(2,ntracks[1]);i++) {
    AliITStrackMI & track= tracks[1][nindexes[1][i]];
    if (track.GetNumberOfClusters()<4) continue;
    if (!constrain&&track.GetNormChi2(1)>7.)continue;
    if (constrain) track.IncrementNSkipped();
    if (!constrain) {
      track.SetD(0,track.GetD(GetX(),GetY()));   
      track.SetNSkipped(track.GetNSkipped()+4./(4.+8.*TMath::Abs(track.GetD(0))));
      if (track.GetNumberOfClusters()+track.GetNDeadZone()+track.GetNSkipped()>6) {
	track.SetNSkipped(6-track.GetNumberOfClusters()+track.GetNDeadZone());
      }
    }
    AddTrackHypothesys(new AliITStrackMI(track), esdindex);
  }
  //}
  
  if (!constrain){  
    for (Int_t i=0;i<TMath::Min(2,ntracks[2]);i++) {
      AliITStrackMI & track= tracks[2][nindexes[2][i]];
      if (track.GetNumberOfClusters()<3) continue;
      if (!constrain&&track.GetNormChi2(2)>7.)continue;
      if (constrain) track.SetNSkipped(track.GetNSkipped()+2);      
      if (!constrain){
	track.SetD(0,track.GetD(GetX(),GetY()));
	track.SetNSkipped(track.GetNSkipped()+7./(7.+8.*TMath::Abs(track.GetD(0))));
	if (track.GetNumberOfClusters()+track.GetNDeadZone()+track.GetNSkipped()>6) {
	  track.SetNSkipped(6-track.GetNumberOfClusters()+track.GetNDeadZone());
	}
      }
      AddTrackHypothesys(new AliITStrackMI(track), esdindex);
    }
  }
  
  if (!constrain){
    //
    // register best tracks - important for V0 finder
    //
    for (Int_t ilayer=0;ilayer<5;ilayer++){
      if (ntracks[ilayer]==0) continue;
      AliITStrackMI & track= tracks[ilayer][nindexes[ilayer][0]];
      if (track.GetNumberOfClusters()<1) continue;
      CookLabel(&track,0);
      bestarray->AddAt(new AliITStrackMI(track),ilayer);
    }
  }
  //
  // update TPC V0 information
  //
  if (otrack->GetESDtrack()->GetV0Index(0)>0){    
    Float_t fprimvertex[3]={GetX(),GetY(),GetZ()};
    for (Int_t i=0;i<3;i++){
      Int_t  index = otrack->GetESDtrack()->GetV0Index(i); 
      if (index==0) break;
      AliV0 * vertex = (AliV0*)fEsd->GetV0(index);
      if (vertex->GetStatus()<0) continue;     // rejected V0
      //
      if (otrack->GetSign()>0) {
	vertex->SetIndex(0,esdindex);
      }
      else{
	vertex->SetIndex(1,esdindex);
      }
      //find nearest layer with track info
      Double_t xrp[3]; vertex->GetXYZ(xrp[0],xrp[1],xrp[2]);  //I.B.
      Int_t nearestold  = GetNearestLayer(xrp);               //I.B.
      Int_t nearest     = nearestold; 
      for (Int_t ilayer =nearest;ilayer<8;ilayer++){
	if (ntracks[nearest]==0){
	  nearest = ilayer;
	}
      }
      //
      AliITStrackMI & track= tracks[nearest][nindexes[nearest][0]];
      if (nearestold<5&&nearest<5){
	Bool_t accept = track.GetNormChi2(nearest)<10; 
	if (accept){
	  if (track.GetSign()>0) {
	    vertex->SetParamP(track);
	    vertex->Update(fprimvertex);
	    //	    vertex->SetIndex(0,track.fESDtrack->GetID()); 
	    if (track.GetNumberOfClusters()>2) AddTrackHypothesys(new AliITStrackMI(track), esdindex);
	  }else{
	    vertex->SetParamN(track);
	    vertex->Update(fprimvertex);
	    //vertex->SetIndex(1,track.fESDtrack->GetID());
	    if (track.GetNumberOfClusters()>2) AddTrackHypothesys(new AliITStrackMI(track), esdindex);
	  }
	  vertex->SetStatus(vertex->GetStatus()+1);
	}else{
	  //  vertex->SetStatus(-2);  // reject V0  - not enough clusters
	}
      }
      // if (nearestold>3){
// 	Int_t indexlayer = (ntracks[0]>0)? 0:1;
// 	if (ntracks[indexlayer]>0){
// 	  AliITStrackMI & track= tracks[indexlayer][nindexes[indexlayer][0]];
// 	  if (track.GetNumberOfClusters()>4&&track.fNormChi2[indexlayer]<4){
// 	    vertex->SetStatus(-1);  // reject V0 - clusters before
// 	  }
// 	}
//      }
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

AliITStrackerMI::AliITSlayer::AliITSlayer():
fR(0),
fPhiOffset(0),
fNladders(0),
fZOffset(0),
fNdetectors(0),
fDetectors(0),
fN(0),
fDy5(0),
fDy10(0),
fDy20(0),
fClustersCs(0),
fClusterIndexCs(0),
fYcs(0),
fZcs(0),
fNcs(0),
fCurrentSlice(-1),
fZmax(0),
fYmin(0),
fYmax(0),
fI(0),
fImax(0),
fSkip(0),
fAccepted(0),
fRoad(0){
  //--------------------------------------------------------------------
  //default AliITSlayer constructor
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) {
    fClusterWeight[i]=0;
    fClusterTracks[0][i]=-1;
    fClusterTracks[1][i]=-1;
    fClusterTracks[2][i]=-1;    
    fClusterTracks[3][i]=-1;    
  }
}

AliITStrackerMI::AliITSlayer::
AliITSlayer(Double_t r,Double_t p,Double_t z,Int_t nl,Int_t nd):
fR(r),
fPhiOffset(p),
fNladders(nl),
fZOffset(z),
fNdetectors(nd),
fDetectors(0),
fN(0),
fDy5(0),
fDy10(0),
fDy20(0),
fClustersCs(0),
fClusterIndexCs(0),
fYcs(0),
fZcs(0),
fNcs(0),
fCurrentSlice(-1),
fZmax(0),
fYmin(0),
fYmax(0),
fI(0),
fImax(0),
fSkip(0),
fAccepted(0),
fRoad(0) {
  //--------------------------------------------------------------------
  //main AliITSlayer constructor
  //--------------------------------------------------------------------
  fDetectors=new AliITSdetector[fNladders*fNdetectors];
  fRoad=2*fR*TMath::Sqrt(3.14/1.);//assuming that there's only one cluster
}

AliITStrackerMI::AliITSlayer::AliITSlayer(const AliITSlayer& layer):
fR(layer.fR),
fPhiOffset(layer.fPhiOffset),
fNladders(layer.fNladders),
fZOffset(layer.fZOffset),
fNdetectors(layer.fNdetectors),
fDetectors(layer.fDetectors),
fN(layer.fN),
fDy5(layer.fDy5),
fDy10(layer.fDy10),
fDy20(layer.fDy20),
fClustersCs(layer.fClustersCs),
fClusterIndexCs(layer.fClusterIndexCs),
fYcs(layer.fYcs),
fZcs(layer.fZcs),
fNcs(layer.fNcs),
fCurrentSlice(layer.fCurrentSlice),
fZmax(layer.fZmax),
fYmin(layer.fYmin),
fYmax(layer.fYmax),
fI(layer.fI),
fImax(layer.fImax),
fSkip(layer.fSkip),
fAccepted(layer.fAccepted),
fRoad(layer.fRoad){
  //Copy constructor
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
    AliITSRecPoint * cl = (AliITSRecPoint*)GetCluster(i);
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


Int_t AliITStrackerMI::AliITSlayer::InsertCluster(AliITSRecPoint *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to this layer
  //--------------------------------------------------------------------
  if (fN==kMaxClusterPerLayer) {
    ::Error("InsertCluster","Too many clusters !\n");
    return 1;
  }
  fCurrentSlice=-1;
  fClusters[fN]=c;
  fN++;
  AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
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
  AliITSRecPoint **clusters = new AliITSRecPoint*[fN];
  Float_t *z                = new Float_t[fN];
  Int_t   * index           = new Int_t[fN];
  //
  for (Int_t i=0;i<fN;i++){
    z[i] = fClusters[i]->GetZ();
  }
  TMath::Sort(fN,z,index,kFALSE);
  for (Int_t i=0;i<fN;i++){
    clusters[i] = fClusters[index[i]];
  }
  //
  for (Int_t i=0;i<fN;i++){
    fClusters[i] = clusters[i];
    fZ[i]        = fClusters[i]->GetZ();
    AliITSdetector &det=GetDetector(fClusters[i]->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + fClusters[i]->GetY();
    if (y>2.*fR*TMath::Pi()) y -= 2.*fR*TMath::Pi();
    fY[i] = y;
  }
  delete[] index;
  delete[] z;
  delete[] clusters;
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
  for (Int_t i=0;i<fN;i++)
    for (Int_t irot=-1;irot<=1;irot++){
      Float_t curY = fY[i]+irot*TMath::TwoPi()*fR; 
      // slice 5
      for (Int_t slice=0; slice<6;slice++){
	if (fBy5[slice][0]<curY && curY<fBy5[slice][1]&&fN5[slice]<kMaxClusterPerLayer5){
	  fClusters5[slice][fN5[slice]] = fClusters[i];
	  fY5[slice][fN5[slice]] = curY;
	  fZ5[slice][fN5[slice]] = fZ[i];
	  fClusterIndex5[slice][fN5[slice]]=i;
	  fN5[slice]++;
	}
      }
      // slice 10
      for (Int_t slice=0; slice<11;slice++){
	if (fBy10[slice][0]<curY && curY<fBy10[slice][1]&&fN10[slice]<kMaxClusterPerLayer10){
	  fClusters10[slice][fN10[slice]] = fClusters[i];
	  fY10[slice][fN10[slice]] = curY;
	  fZ10[slice][fN10[slice]] = fZ[i];
	  fClusterIndex10[slice][fN10[slice]]=i;
	  fN10[slice]++;
	}
      }
      // slice 20
      for (Int_t slice=0; slice<21;slice++){
	if (fBy20[slice][0]<curY && curY<fBy20[slice][1]&&fN20[slice]<kMaxClusterPerLayer20){
	  fClusters20[slice][fN20[slice]] = fClusters[i];
	  fY20[slice][fN20[slice]] = curY;
	  fZ20[slice][fN20[slice]] = fZ[i];
	  fClusterIndex20[slice][fN20[slice]]=i;
	  fN20[slice]++;
	}
      }      
    }

  //
  // consistency check
  //
  for (Int_t i=0;i<fN-1;i++){
    if (fZ[i]>fZ[i+1]){
      printf("Bugg\n");
    }
  }
  //
  for (Int_t slice=0;slice<21;slice++)
  for (Int_t i=0;i<fN20[slice]-1;i++){
    if (fZ20[slice][i]>fZ20[slice][i+1]){
      printf("Bugg\n");
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


const AliITSRecPoint *AliITStrackerMI::AliITSlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------

  if (fCurrentSlice<0){
    Double_t rpi2 = 2.*fR*TMath::Pi();
    for (Int_t i=fI; i<fImax; i++) {
      Double_t y = fY[i];
      if (fYmax<y) y -= rpi2;
      if (fYmin>y) y += rpi2;
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
    const AliITSRecPoint *c=fClusters[i];
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

Bool_t AliITStrackerMI::RefitAt(Double_t xx,AliITStrackMI *t,
				const AliITStrackMI *c, Bool_t extra) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"
  // If "extra"==kTRUE, 
  //    the clusters from overlapped modules get attached to "t" 
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

     const AliITSRecPoint *cl=0;
     Double_t maxchi2=1000.*kMaxChi2;

     Int_t idx=index[i];
     if (idx>0) {
        const AliITSRecPoint *c=(AliITSRecPoint *)GetCluster(idx); 
	if (c){
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
     }

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
                 
     if (extra) { //search for extra clusters
        AliITStrackV2 tmp(*t);
        Double_t dz=4*TMath::Sqrt(tmp.GetSigmaZ2()+kSigmaZ2[i]);
        if (dz < 0.5*TMath::Abs(tmp.GetTgl())) dz=0.5*TMath::Abs(tmp.GetTgl());
        Double_t dy=4*TMath::Sqrt(t->GetSigmaY2()+kSigmaY2[i]);
        if (dy < 0.5*TMath::Abs(tmp.GetSnp())) dy=0.5*TMath::Abs(tmp.GetSnp());
        Double_t zmin=t->GetZ() - dz;
        Double_t zmax=t->GetZ() + dz;
        Double_t ymin=t->GetY() + phi*r - dy;
        Double_t ymax=t->GetY() + phi*r + dy;
        layer.SelectClusters(zmin,zmax,ymin,ymax);

        const AliITSRecPoint *c=0; Int_t ci=-1,cci=-1;
        Double_t maxchi2=1000.*kMaxChi2, tolerance=0.1;
        while ((c=layer.GetNextCluster(ci))!=0) {
           if (idet == c->GetDetectorIndex()) continue;

	   const AliITSdetector &det=layer.GetDetector(c->GetDetectorIndex());

	   if (!tmp.Propagate(det.GetPhi(),det.GetR())) continue;
           
	   if (TMath::Abs(tmp.GetZ() - c->GetZ()) > tolerance) continue;
           if (TMath::Abs(tmp.GetY() - c->GetY()) > tolerance) continue;

           Double_t chi2=tmp.GetPredictedChi2(c);
           if (chi2<maxchi2) { maxchi2=chi2; cci=ci; }
        }
        if (cci>=0) t->SetExtraCluster(i,(i<<28)+cci);
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

Bool_t 
AliITStrackerMI::RefitAt(Double_t xx,AliITStrackMI *t,const Int_t *clindex) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from array
  //--------------------------------------------------------------------
  Int_t index[kMaxLayer];
  Int_t k;
  for (k=0; k<kMaxLayer; k++) index[k]=-1;
  //
  for (k=0; k<kMaxLayer; k++) { 
    index[k]=clindex[k]; 
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
     if (step<0 && xx>r) break;  //
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

     const AliITSRecPoint *cl=0;
     Double_t maxchi2=1000.*kMaxChi2;

     Int_t idx=index[i];
     if (idx>0) {
        const AliITSRecPoint *c=(AliITSRecPoint *)GetCluster(idx); 
	if (c){
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

        const AliITSRecPoint *c=0; Int_t ci=-1;
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

Double_t AliITStrackerMI::GetNormalizedChi2(AliITStrackMI * track, Int_t mode)
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
      if (track->GetClIndex(i)>0){
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else 
	  { cerry= track->GetSigmaY(i); cerrz = track->GetSigmaZ(i);}
	cerry*=cerry;
	cerrz*=cerrz;	
	Float_t cchi2 = (track->GetDy(i)*track->GetDy(i)/cerry)+(track->GetDz(i)*track->GetDz(i)/cerrz);
	if (i>1){
	  Float_t ratio = track->GetNormQ(i)/track->GetExpQ();
	  if (ratio<0.5) {
	    cchi2+=(0.5-ratio)*10.;
	    //track->fdEdxMismatch+=(0.5-ratio)*10.;
	    dedxmismatch+=(0.5-ratio)*10.;	    
	  }
	}
	if (i<2 ||i>3){
	  AliITSRecPoint * cl = (AliITSRecPoint*)GetCluster( track->GetClIndex(i));  
	  Double_t delta = cl->GetNy()+cl->GetNz()-ny[i]-nz[i];
	  if (delta>1) chi2 +=0.5*TMath::Min(delta/2,2.); 
	  if (i<2) chi2+=2*cl->GetDeltaProbability();
	}
	chi2+=cchi2;
	sum++;
      }
    }
    if (TMath::Abs(track->GetdEdxMismatch()-dedxmismatch)>0.0001){
      track->SetdEdxMismatch(dedxmismatch);
    }
  }
  else{
    for (Int_t i = 0;i<4;i++){
      if (track->GetClIndex(i)>0){
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else { cerry= track->GetSigmaY(i); cerrz = track->GetSigmaZ(i);}
	cerry*=cerry;
	cerrz*=cerrz;
	chi2+= (track->GetDy(i)*track->GetDy(i)/cerry);
	chi2+= (track->GetDz(i)*track->GetDz(i)/cerrz);      
	sum++;
      }
    }
    for (Int_t i = 4;i<6;i++){
      if (track->GetClIndex(i)>0){	
	Float_t cerry, cerrz;
	if (ny[i]>0) {cerry = erry[i]; cerrz=errz[i];}
	else { cerry= track->GetSigmaY(i); cerrz = track->GetSigmaZ(i);}
	cerry*=cerry;
	cerrz*=cerrz;	
	Float_t cerryb, cerrzb;
	if (ny[i+6]>0) {cerryb = erry[i+6]; cerrzb=errz[i+6];}
	else { cerryb= track->GetSigmaY(i+6); cerrzb = track->GetSigmaZ(i+6);}
	cerryb*=cerryb;
	cerrzb*=cerrzb;
	chi2+= TMath::Min((track->GetDy(i+6)*track->GetDy(i+6)/cerryb),track->GetDy(i)*track->GetDy(i)/cerry);
	chi2+= TMath::Min((track->GetDz(i+6)*track->GetDz(i+6)/cerrzb),track->GetDz(i)*track->GetDz(i)/cerrz);      
	sum++;
      }
    }
  }
  if (track->GetESDtrack()->GetTPCsignal()>85){
    Float_t ratio = track->GetdEdx()/track->GetESDtrack()->GetTPCsignal();
    if (ratio<0.5) {
      chi2+=(0.5-ratio)*5.;      
    }
    if (ratio>2){
      chi2+=(ratio-2.0)*3; 
    }
  }
  //
  Double_t match = TMath::Sqrt(track->GetChi22());
  if (track->GetConstrain())  match/=track->GetNumberOfClusters();
  if (!track->GetConstrain()) match/=track->GetNumberOfClusters()-2.;
  if (match<0) match=0;
  Float_t deadzonefactor = (track->GetNDeadZone()>0) ? 3*(1.1-track->GetDeadZoneProbability()):0.;
  Double_t normchi2 = 2*track->GetNSkipped()+match+deadzonefactor+(1+(2*track->GetNSkipped()+deadzonefactor)/track->GetNumberOfClusters())*
    (chi2)/TMath::Max(double(sum-track->GetNSkipped()),
				1./(1.+track->GetNSkipped()));     
 
 return normchi2;
}


Double_t AliITStrackerMI::GetMatchingChi2(AliITStrackMI * track1, AliITStrackMI * track2)
{
  //
  // return matching chi2 between two tracks
  AliITStrackMI track3(*track2);
  track3.Propagate(track1->GetAlpha(),track1->GetX());
  TMatrixD vec(5,1);
  vec(0,0)=track1->GetY()   - track3.GetY();
  vec(1,0)=track1->GetZ()   - track3.GetZ();
  vec(2,0)=track1->GetSnp() - track3.GetSnp();
  vec(3,0)=track1->GetTgl() - track3.GetTgl();
  vec(4,0)=track1->Get1Pt() - track3.Get1Pt();
  //
  TMatrixD cov(5,5);
  cov(0,0) = track1->GetSigmaY2()+track3.GetSigmaY2();
  cov(1,1) = track1->GetSigmaZ2()+track3.GetSigmaZ2();
  cov(2,2) = track1->GetSigmaSnp2()+track3.GetSigmaSnp2();
  cov(3,3) = track1->GetSigmaTgl2()+track3.GetSigmaTgl2();
  cov(4,4) = track1->GetSigma1Pt2()+track3.GetSigma1Pt2();
  
  cov(0,1)=cov(1,0) = track1->GetSigmaZY()+track3.GetSigmaZY();
  cov(0,2)=cov(2,0) = track1->GetSigmaSnpY()+track3.GetSigmaSnpY();
  cov(0,3)=cov(3,0) = track1->GetSigmaTglY()+track3.GetSigmaTglY();
  cov(0,4)=cov(4,0) = track1->GetSigma1PtY()+track3.GetSigma1PtY();
  //
  cov(1,2)=cov(2,1) = track1->GetSigmaSnpZ()+track3.GetSigmaSnpZ();
  cov(1,3)=cov(3,1) = track1->GetSigmaTglZ()+track3.GetSigmaTglZ();
  cov(1,4)=cov(4,1) = track1->GetSigma1PtZ()+track3.GetSigma1PtZ();
  //
  cov(2,3)=cov(3,2) = track1->GetSigmaTglSnp()+track3.GetSigmaTglSnp();
  cov(2,4)=cov(4,2) = track1->GetSigma1PtSnp()+track3.GetSigma1PtSnp();
  //
  cov(3,4)=cov(4,3) = track1->GetSigma1PtTgl()+track3.GetSigma1PtTgl();
  
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


Double_t AliITStrackerMI::GetTruncatedChi2(AliITStrackMI * track, Float_t fac)
{
  //
  // calculate normalized chi2
  Float_t chi2[6];
  Float_t *erry = GetErrY(fCurrentEsdTrack), *errz = GetErrZ(fCurrentEsdTrack);
  Float_t ncl = 0;
  for (Int_t i = 0;i<6;i++){
    if (TMath::Abs(track->GetDy(i))>0){      
      chi2[i]= (track->GetDy(i)/erry[i])*(track->GetDy(i)/erry[i]);
      chi2[i]+= (track->GetDz(i)/errz[i])*(track->GetDz(i)/errz[i]);
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


Double_t AliITStrackerMI::GetInterpolatedChi2(AliITStrackMI * forwardtrack, AliITStrackMI * backtrack)
{
  //
  // calculate normalized chi2
  //  if (forwardtrack->fNUsed>0.3*float(forwardtrack->GetNumberOfClusters())) return 10000;
  Int_t npoints = 0;
  Double_t res =0;
  for (Int_t i=0;i<6;i++){
    if ( (backtrack->GetSigmaY(i)<0.000000001) || (forwardtrack->GetSigmaY(i)<0.000000001)) continue;
    Double_t sy1 = forwardtrack->GetSigmaY(i);
    Double_t sz1 = forwardtrack->GetSigmaZ(i);
    Double_t sy2 = backtrack->GetSigmaY(i);
    Double_t sz2 = backtrack->GetSigmaZ(i);
    if (i<2){ sy2=1000.;sz2=1000;}
    //    
    Double_t dy0 = (forwardtrack->GetDy(i)/(sy1*sy1) +backtrack->GetDy(i)/(sy2*sy2))/(1./(sy1*sy1)+1./(sy2*sy2));
    Double_t dz0 = (forwardtrack->GetDz(i)/(sz1*sz1) +backtrack->GetDz(i)/(sz2*sz2))/(1./(sz1*sz1)+1./(sz2*sz2));
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
		   res/TMath::Max(double(npoints-forwardtrack->GetNSkipped()),
				  1./(1.+forwardtrack->GetNSkipped()));
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

void AliITStrackerMI::RegisterClusterTracks(AliITStrackMI* track,Int_t id)
{
  //---------------------------------------------
  // register track to the list
  //
  if (track->GetESDtrack()->GetKinkIndex(0)!=0) return;  //don't register kink tracks
  //
  //
  for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){
    Int_t index = track->GetClusterIndex(icluster);
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].GetNumberOfClusters()) continue;
    for (Int_t itrack=0;itrack<4;itrack++){
      if (fgLayers[l].GetClusterTracks(itrack,c)<0){
	fgLayers[l].SetClusterTracks(itrack,c,id);
	break;
      }
    }
  }
}
void AliITStrackerMI::UnRegisterClusterTracks(AliITStrackMI* track, Int_t id)
{
  //---------------------------------------------
  // unregister track from the list
  for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){
    Int_t index = track->GetClusterIndex(icluster);
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    if (c>fgLayers[l].GetNumberOfClusters()) continue;
    for (Int_t itrack=0;itrack<4;itrack++){
      if (fgLayers[l].GetClusterTracks(itrack,c)==id){
	fgLayers[l].SetClusterTracks(itrack,c,-1);
      }
    }
  }
}
Float_t AliITStrackerMI::GetNumberOfSharedClusters(AliITStrackMI* track,Int_t id, Int_t list[6], AliITSRecPoint *clist[6])
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
    if (c>fgLayers[l].GetNumberOfClusters()) continue;
    if (ny[l]==0){
      printf("problem\n");
    }
    AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(index);
    Float_t weight=1;
    //
    Float_t deltan = 0;
    if (l>3&&cl->GetNy()+cl->GetNz()>6) continue;
    if (l>2&&track->GetNormQ(l)/track->GetExpQ()>3.5) continue;
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
      if (fgLayers[l].GetClusterTracks(itrack,c)>=0 && fgLayers[l].GetClusterTracks(itrack,c)!=id){
	list[l]=index;
	clist[l] = (AliITSRecPoint*)GetCluster(index);
	shared+=weight; 
	break;
      }
    }
  }
  track->SetNUsed(shared);
  return shared;
}

Int_t AliITStrackerMI::GetOverlapTrack(AliITStrackMI *track, Int_t trackID, Int_t &shared, Int_t clusterlist[6],Int_t overlist[6])
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
    if (c>fgLayers[l].GetNumberOfClusters()) continue;
    //if (l>3) continue;
    AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(index);
    //
    Float_t deltan = 0;
    if (l>3&&cl->GetNy()+cl->GetNz()>6) continue;
    if (l>2&&track->GetNormQ(l)/track->GetExpQ()>3.5) continue;
    if (l<2 || l>3){      
      deltan = (cl->GetNy()+cl->GetNz()-ny[l]-nz[l]);
    }
    else{
      deltan = (cl->GetNz()-nz[l]);
    }
    if (deltan>2.0) continue;  // extended - highly probable shared cluster
    //
    for (Int_t itrack=3;itrack>=0;itrack--){
      if (fgLayers[l].GetClusterTracks(itrack,c)<0) continue;
      if (fgLayers[l].GetClusterTracks(itrack,c)!=trackID){
       tracks[trackindex]  = fgLayers[l].GetClusterTracks(itrack,c);
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
    if (c>fgLayers[l].GetNumberOfClusters()) continue;
    AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(index);
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
      if (fgLayers[l].GetClusterTracks(itrack,c)<0) continue;
      if (fgLayers[l].GetClusterTracks(itrack,c)==sharedtrack){
	overlist[l]=index;
	shared++;      
      }
    }
  }
  return sharedtrack;
}


AliITStrackMI *  AliITStrackerMI::GetBest2Tracks(Int_t trackID1, Int_t trackID2, Float_t th0, Float_t th1){
  //
  // try to find track hypothesys without conflicts
  // with minimal chi2;
  TClonesArray *arr1 = (TClonesArray*)fTrackHypothesys.At(trackID1);
  Int_t entries1 = arr1->GetEntriesFast();
  TClonesArray *arr2 = (TClonesArray*)fTrackHypothesys.At(trackID2);
  if (!arr2) return (AliITStrackMI*) arr1->UncheckedAt(0);
  Int_t entries2 = arr2->GetEntriesFast();
  if (entries2<=0) return (AliITStrackMI*) arr1->UncheckedAt(0);
  //
  AliITStrackMI * track10=(AliITStrackMI*) arr1->UncheckedAt(0);
  AliITStrackMI * track20=(AliITStrackMI*) arr2->UncheckedAt(0);
  if (TMath::Abs(1./track10->Get1Pt())>0.5+TMath::Abs(1/track20->Get1Pt())) return track10;

  for (Int_t itrack=0;itrack<entries1;itrack++){
    AliITStrackMI * track=(AliITStrackMI*) arr1->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID1);
  }
  //
  for (Int_t itrack=0;itrack<entries2;itrack++){
    AliITStrackMI * track=(AliITStrackMI*) arr2->UncheckedAt(itrack);
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
  AliITSRecPoint *clist1[6], *clist2[6] ;
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
     if (TMath::Abs(track10->GetDy(i))>0.000000000000001){
       chi21 += track10->GetDy(i)*track10->GetDy(i)/(nerry[i]*nerry[i]);
       chi21 += track10->GetDz(i)*track10->GetDz(i)/(nerrz[i]*nerrz[i]);
       ncl1++;
     }
     if (TMath::Abs(track20->GetDy(i))>0.000000000000001){
       chi22 += track20->GetDy(i)*track20->GetDy(i)/(nerry[i]*nerry[i]);
       chi22 += track20->GetDz(i)*track20->GetDz(i)/(nerrz[i]*nerrz[i]);
       ncl2++;
     }
  }
  chi21/=ncl1;
  chi22/=ncl2;
  //
  // 
  Float_t d1 = TMath::Sqrt(track10->GetD(0)*track10->GetD(0)+track10->GetD(1)*track10->GetD(1))+0.1;
  Float_t d2 = TMath::Sqrt(track20->GetD(0)*track20->GetD(0)+track20->GetD(1)*track20->GetD(1))+0.1;
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
    AliITStrackMI * track1=(AliITStrackMI*) arr1->UncheckedAt(itrack1);
    //if (track1->fFakeRatio>0) continue;
    RegisterClusterTracks(track1,trackID1);
    for (Int_t itrack2=0;itrack2<entries2;itrack2++){
      AliITStrackMI * track2=(AliITStrackMI*) arr2->UncheckedAt(itrack2);

      //      Float_t current = w1*track1->fChi2MIP[0]+w2*track2->fChi2MIP[0];
      //if (track2->fFakeRatio>0) continue;
      Float_t nskipped=0;            
      RegisterClusterTracks(track2,trackID2);
      Int_t list1[6],list2[6];
      AliITSRecPoint *clist1[6], *clist2[6] ;
      Float_t cconflict1 = GetNumberOfSharedClusters(track1,trackID1,list1,clist1);
      Float_t cconflict2 = GetNumberOfSharedClusters(track2,trackID2,list2,clist2);
      UnRegisterClusterTracks(track2,trackID2);
      //
      if (track1->GetConstrain()) nskipped+=w1*track1->GetNSkipped();
      if (track2->GetConstrain()) nskipped+=w2*track2->GetNSkipped();
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
	if (TMath::Abs(track1->GetDy(i))>0.) {
	  chi21 = (track1->GetDy(i)/track1->GetSigmaY(i))*(track1->GetDy(i)/track1->GetSigmaY(i))+
	    (track1->GetDz(i)/track1->GetSigmaZ(i))*(track1->GetDz(i)/track1->GetSigmaZ(i));
	  //chi21 = (track1->fDy[i]*track1->fDy[i])/(nerry[i]*nerry[i])+
	  //  (track1->GetDz(i)*track1->GetDz(i))/(nerrz[i]*nerrz[i]);
	}else{
	  if (TMath::Abs(track1->GetSigmaY(i)>0.)) c1=1;
	}
	//
	if (TMath::Abs(track2->GetDy(i))>0.) {
	  chi22 = (track2->GetDy(i)/track2->GetSigmaY(i))*(track2->GetDy(i)/track2->GetSigmaY(i))+
	    (track2->GetDz(i)/track2->GetSigmaZ(i))*(track2->GetDz(i)/track2->GetSigmaZ(i));
	  //chi22 = (track2->fDy[i]*track2->fDy[i])/(nerry[i]*nerry[i])+
	  //  (track2->fDz[i]*track2->fDz[i])/(nerrz[i]*nerrz[i]);
	}
	else{
	  if (TMath::Abs(track2->GetSigmaY(i)>0.)) c2=1;
	}
	sumchi2+=w1*(1.+c1)*(1+c1)*(chi21+c1)+w2*(1.+c2)*(1+c2)*(chi22+c2);
	if (chi21>0) sum+=w1;
	if (chi22>0) sum+=w2;
	conflict+=(c1+c2);
      }
      Double_t norm = sum-w1*track1->GetNSkipped()-w2*track2->GetNSkipped();
      if (norm<0) norm =1/(w1*track1->GetNSkipped()+w2*track2->GetNSkipped());
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
    Float_t orig = track10->GetFakeRatio()*track10->GetNumberOfClusters();
    AliITStrackMI* track1=(AliITStrackMI*) arr1->UncheckedAt(index1);
    track1->SetChi2MIP(5,maxconflicts);
    track1->SetChi2MIP(6,maxchi2);
    track1->SetChi2MIP(7,0.01+orig-(track1->GetFakeRatio()*track1->GetNumberOfClusters()));
    //    track1->UpdateESDtrack(AliESDtrack::kITSin);
    track1->SetChi2MIP(8,index1);
    fBestTrackIndex[trackID1] =index1;
    UpdateESDtrack(track1, AliESDtrack::kITSin);
  }  
  else if (track10->GetChi2MIP(0)<th1){
    track10->SetChi2MIP(5,maxconflicts);
    track10->SetChi2MIP(6,maxchi2);    
    //    track10->UpdateESDtrack(AliESDtrack::kITSin);
    UpdateESDtrack(track10,AliESDtrack::kITSin);
  }   
  
  for (Int_t itrack=0;itrack<entries1;itrack++){
    AliITStrackMI * track=(AliITStrackMI*) arr1->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID1);
  }
  //
  for (Int_t itrack=0;itrack<entries2;itrack++){
    AliITStrackMI * track=(AliITStrackMI*) arr2->UncheckedAt(itrack);
    UnRegisterClusterTracks(track,trackID2);
  }

  if (track10->GetConstrain()&&track10->GetChi2MIP(0)<kMaxChi2PerCluster[0]&&track10->GetChi2MIP(1)<kMaxChi2PerCluster[1]
      &&track10->GetChi2MIP(2)<kMaxChi2PerCluster[2]&&track10->GetChi2MIP(3)<kMaxChi2PerCluster[3]){ 
    //  if (track10->fChi2MIP[0]<kMaxChi2PerCluster[0]&&track10->fChi2MIP[1]<kMaxChi2PerCluster[1]
  //    &&track10->fChi2MIP[2]<kMaxChi2PerCluster[2]&&track10->fChi2MIP[3]<kMaxChi2PerCluster[3]){ 
    RegisterClusterTracks(track10,trackID1);
  }
  if (track20->GetConstrain()&&track20->GetChi2MIP(0)<kMaxChi2PerCluster[0]&&track20->GetChi2MIP(1)<kMaxChi2PerCluster[1]
      &&track20->GetChi2MIP(2)<kMaxChi2PerCluster[2]&&track20->GetChi2MIP(3)<kMaxChi2PerCluster[3]){ 
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

  AliITSRecPoint *c=(AliITSRecPoint *)GetCluster(t->GetClusterIndex(0));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();
  c=(AliITSRecPoint *)GetCluster(t->GetClusterIndex(1));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();

}


void AliITStrackerMI::AddTrackHypothesys(AliITStrackMI * track, Int_t esdindex)
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
  AliITStrackMI * besttrack=0;
  for (Int_t itrack=0;itrack<array->GetEntriesFast();itrack++){
    AliITStrackMI * track = (AliITStrackMI*)array->At(itrack);
    if (!track) continue;
    Float_t chi2 = NormalizedChi2(track,0);
    //
    Int_t tpcLabel=track->GetESDtrack()->GetTPCLabel();
    track->SetLabel(tpcLabel);
    CookdEdx(track);
    track->SetFakeRatio(1.);
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
      if (track->GetConstrain() || track->GetNumberOfClusters()>5){  //keep best short tracks - without vertex constrain
	delete array->RemoveAt(itrack);
      }	 
    }
  }
  if (!besttrack) return;
  //
  //
  //take errors of best track as a reference
  Float_t *erry = GetErrY(esdindex), *errz = GetErrZ(esdindex);
  Float_t *ny = GetNy(esdindex), *nz = GetNz(esdindex);
  for (Int_t i=0;i<6;i++) {
    if (besttrack->GetClIndex(i)>0){
      erry[i] = besttrack->GetSigmaY(i); erry[i+6] = besttrack->GetSigmaY(i+6);
      errz[i] = besttrack->GetSigmaZ(i); errz[i+6] = besttrack->GetSigmaZ(i+6);
      ny[i]   = besttrack->GetNy(i);
      nz[i]   = besttrack->GetNz(i);
    }
  }
  //
  // calculate normalized chi2
  //
  Float_t * chi2        = new Float_t[entries];
  Int_t * index         = new Int_t[entries];  
  for (Int_t i=0;i<entries;i++) chi2[i] =10000;
  for (Int_t itrack=0;itrack<entries;itrack++){
    AliITStrackMI * track = (AliITStrackMI*)array->At(itrack);
    if (track){
      track->SetChi2MIP(0,GetNormalizedChi2(track, mode));            
      if (track->GetChi2MIP(0)<kMaxChi2PerCluster[0]) 
	chi2[itrack] = track->GetChi2MIP(0);
      else{
	if (track->GetConstrain() || track->GetNumberOfClusters()>5){  //keep best short tracks - without vertex constrain
	  delete array->RemoveAt(itrack);	     
	}
      }
    }
  }
  //
  TMath::Sort(entries,chi2,index,kFALSE);
  besttrack = (AliITStrackMI*)array->At(index[0]);
  if (besttrack&&besttrack->GetChi2MIP(0)<kMaxChi2PerCluster[0]){
    for (Int_t i=0;i<6;i++){
      if (besttrack->GetClIndex(i)>0){
	erry[i] = besttrack->GetSigmaY(i); erry[i+6] = besttrack->GetSigmaY(i+6);
	errz[i] = besttrack->GetSigmaZ(i); erry[i+6] = besttrack->GetSigmaY(i+6);
	ny[i]   = besttrack->GetNy(i);
	nz[i]   = besttrack->GetNz(i);
      }
    }
  }
  //
  // calculate one more time with updated normalized errors
  for (Int_t i=0;i<entries;i++) chi2[i] =10000;  
  for (Int_t itrack=0;itrack<entries;itrack++){
    AliITStrackMI * track = (AliITStrackMI*)array->At(itrack);
    if (track){      
      track->SetChi2MIP(0,GetNormalizedChi2(track,mode));            
      if (track->GetChi2MIP(0)<kMaxChi2PerCluster[0]) 
	chi2[itrack] = track->GetChi2MIP(0)-0*(track->GetNumberOfClusters()+track->GetNDeadZone()); 
      else
	{
	  if (track->GetConstrain() || track->GetNumberOfClusters()>5){  //keep best short tracks - without vertex constrain
	    delete array->RemoveAt(itrack);	
	  }
	}
    }   
  }
  entries = array->GetEntriesFast();  
  //
  //
  if (entries>0){
    TObjArray * newarray = new TObjArray();  
    TMath::Sort(entries,chi2,index,kFALSE);
    besttrack = (AliITStrackMI*)array->At(index[0]);
    if (besttrack){
      //
      for (Int_t i=0;i<6;i++){
	if (besttrack->GetNz(i)>0&&besttrack->GetNy(i)>0){
	  erry[i] = besttrack->GetSigmaY(i); erry[i+6] = besttrack->GetSigmaY(i+6);
	  errz[i] = besttrack->GetSigmaZ(i); errz[i+6] = besttrack->GetSigmaZ(i+6);
	  ny[i]   = besttrack->GetNy(i);
	  nz[i]   = besttrack->GetNz(i);
	}
      }
      besttrack->SetChi2MIP(0,GetNormalizedChi2(besttrack,mode));
      Float_t minchi2 = TMath::Min(besttrack->GetChi2MIP(0)+5.+besttrack->GetNUsed(), double(kMaxChi2PerCluster[0]));
      Float_t minn = besttrack->GetNumberOfClusters()-3;
      Int_t accepted=0;
      for (Int_t i=0;i<entries;i++){
	AliITStrackMI * track = (AliITStrackMI*)array->At(index[i]);	
	if (!track) continue;
	if (accepted>maxcut) break;
	track->SetChi2MIP(0,GetNormalizedChi2(track,mode));
	if (track->GetConstrain() || track->GetNumberOfClusters()>5){  //keep best short tracks - without vertex constrain
	  if (track->GetNumberOfClusters()<6 && (track->GetChi2MIP(0)+track->GetNUsed()>minchi2)){
	    delete array->RemoveAt(index[i]);
	    continue;
	  }
	}
	Bool_t shortbest = !track->GetConstrain() && track->GetNumberOfClusters()<6;
	if ((track->GetChi2MIP(0)+track->GetNUsed()<minchi2 && track->GetNumberOfClusters()>=minn) ||shortbest){
	  if (!shortbest) accepted++;
	  //
	  newarray->AddLast(array->RemoveAt(index[i]));      
	  for (Int_t i=0;i<6;i++){
	    if (nz[i]==0){
	      erry[i] = track->GetSigmaY(i); erry[i+6] = track->GetSigmaY(i+6);
	      errz[i] = track->GetSigmaZ(i); errz[i]   = track->GetSigmaZ(i+6);
	      ny[i]   = track->GetNy(i);
	      nz[i]   = track->GetNz(i);
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



AliITStrackMI * AliITStrackerMI::GetBestHypothesys(Int_t esdindex, AliITStrackMI * original, Int_t checkmax)
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
  AliITStrackMI * besttrack=0;
  //
  AliITStrackMI * backtrack    = new AliITStrackMI(*original);
  AliITStrackMI * forwardtrack = new AliITStrackMI(*original);
  Double_t xyzv[]={GetX(),GetY(),GetZ()};	
  Double_t ersv[]={GetSigmaX()/3.,GetSigmaY()/3.,GetSigmaZ()/3.};
  //
  for (Int_t i=0;i<entries;i++){    
    AliITStrackMI * track = (AliITStrackMI*)array->At(i);    
    if (!track) continue;
    Float_t sigmarfi,sigmaz;
    GetDCASigma(track,sigmarfi,sigmaz);
    track->SetDnorm(0,sigmarfi);
    track->SetDnorm(1,sigmaz);
    //
    track->SetChi2MIP(1,1000000);
    track->SetChi2MIP(2,1000000);
    track->SetChi2MIP(3,1000000);
    //
    // backtrack
    backtrack = new(backtrack) AliITStrackMI(*track); 
    if (track->GetConstrain()){
      if (!backtrack->PropagateTo(3.,0.0028,65.19)) continue;
      if (!backtrack->Improve(0,xyzv,ersv))         continue;      
      //if (!backtrack->PropagateTo(2.,0.0028,0))     continue; // This 
      //if (!backtrack->Improve(0,xyzv,ersv))         continue; // is
      //if (!backtrack->PropagateTo(1.,0.0028,0))     continue; // an over-kill
      //if (!backtrack->Improve(0,xyzv,ersv))         continue; //   (I.B.)
      //if (!backtrack->PropagateToVertex())          continue; //
      backtrack->ResetCovariance(10.);      
      //if (!backtrack->Improve(0,xyzv,ersv))         continue;            	        
    }else{
      backtrack->ResetCovariance(10.);
    }
    backtrack->ResetClusters();

    Double_t x = original->GetX();
    if (!RefitAt(x,backtrack,track)) continue;
    //
    track->SetChi2MIP(1,NormalizedChi2(backtrack,0));
    //for (Int_t i=2;i<6;i++){track->fDy[i]+=backtrack->fDy[i]; track->fDz[i]+=backtrack->fDz[i];}
    if (track->GetChi2MIP(1)>kMaxChi2PerCluster[1]*6.)  continue;
    track->SetChi22(GetMatchingChi2(backtrack,original));

    if ((track->GetConstrain()) && track->GetChi22()>90.)  continue;
    if ((!track->GetConstrain()) && track->GetChi22()>30.)  continue;
    if ( track->GetChi22()/track->GetNumberOfClusters()>11.)  continue;


    if  (!(track->GetConstrain())&&track->GetChi2MIP(1)>kMaxChi2PerCluster[1])  continue;
    Bool_t isOK=kTRUE;
    if(!isOK) continue;
    //
    //forward track - without constraint
    forwardtrack = new(forwardtrack) AliITStrackMI(*original);
    forwardtrack->ResetClusters();
    x = track->GetX();
    RefitAt(x,forwardtrack,track);
    track->SetChi2MIP(2,NormalizedChi2(forwardtrack,0));    
    if  (track->GetChi2MIP(2)>kMaxChi2PerCluster[2]*6.0)  continue;
    if  (!(track->GetConstrain())&&track->GetChi2MIP(2)>kMaxChi2PerCluster[2])  continue;
    
    //track->fD[0] = forwardtrack->GetD(GetX(),GetY());
    //track->fD[1] = forwardtrack->GetZat(GetX())-GetZ();
    forwardtrack->GetDZ(GetX(),GetY(),GetZ(),track->GetDP());   //I.B.
    forwardtrack->SetD(0,track->GetD(0));
    forwardtrack->SetD(1,track->GetD(1));    
    {
      Int_t list[6];
      AliITSRecPoint* clist[6];
      track->SetChi2MIP(4,GetNumberOfSharedClusters(track,esdindex,list,clist));      
      if ( (!track->GetConstrain()) && track->GetChi2MIP(4)>1.0) continue;
    }
    
    track->SetChi2MIP(3,GetInterpolatedChi2(forwardtrack,backtrack));
    if  ( (track->GetChi2MIP(3)>6.*kMaxChi2PerCluster[3])) continue;    
    if  ( (!track->GetConstrain()) && (track->GetChi2MIP(3)>2*kMaxChi2PerCluster[3])) {
      track->SetChi2MIP(3,1000);
      continue; 
    }
    Double_t chi2 = track->GetChi2MIP(0)+track->GetNUsed();    
    //
    for (Int_t ichi=0;ichi<5;ichi++){
      forwardtrack->SetChi2MIP(ichi, track->GetChi2MIP(ichi));
    }
    if (chi2 < minchi2){
      //besttrack = new AliITStrackMI(*forwardtrack);
      besttrack = track;
      besttrack->SetLabel(track->GetLabel());
      besttrack->SetFakeRatio(track->GetFakeRatio());
      minchi2   = chi2;
      //original->fD[0] = forwardtrack->GetD(GetX(),GetY());
      //original->fD[1] = forwardtrack->GetZat(GetX())-GetZ();
      forwardtrack->GetDZ(GetX(),GetY(),GetZ(),original->GetDP());    //I.B.
    }    
  }
  delete backtrack;
  delete forwardtrack;
  Int_t accepted=0;
  for (Int_t i=0;i<entries;i++){    
    AliITStrackMI * track = (AliITStrackMI*)array->At(i);   
    if (!track) continue;
    
    if (accepted>checkmax || track->GetChi2MIP(3)>kMaxChi2PerCluster[3]*6. || 
	(track->GetNumberOfClusters()<besttrack->GetNumberOfClusters()-1.)||
	track->GetChi2MIP(0)>besttrack->GetChi2MIP(0)+2.*besttrack->GetNUsed()+3.){
      if (track->GetConstrain() || track->GetNumberOfClusters()>5){  //keep best short tracks - without vertex constrain
	delete array->RemoveAt(i);    
	continue;
      }
    }
    else{
      accepted++;
    }
  }
  //
  array->Compress();
  SortTrackHypothesys(esdindex,checkmax,1);
  array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) return 0; // PH What can be the reason? Check SortTrackHypothesys
  besttrack = (AliITStrackMI*)array->At(0);  
  if (!besttrack)  return 0;
  besttrack->SetChi2MIP(8,0);
  fBestTrackIndex[esdindex]=0;
  entries = array->GetEntriesFast();
  AliITStrackMI *longtrack =0;
  minchi2 =1000;
  Float_t minn=besttrack->GetNumberOfClusters()+besttrack->GetNDeadZone();
  for (Int_t itrack=entries-1;itrack>0;itrack--){
    AliITStrackMI * track = (AliITStrackMI*)array->At(itrack);
    if (!track->GetConstrain()) continue;
    if (track->GetNumberOfClusters()+track->GetNDeadZone()<minn) continue;
    if (track->GetChi2MIP(0)-besttrack->GetChi2MIP(0)>0.0) continue;
    if (track->GetChi2MIP(0)>4.) continue;
    minn = track->GetNumberOfClusters()+track->GetNDeadZone();
    longtrack =track;
  }
  //if (longtrack) besttrack=longtrack;

  Int_t list[6];
  AliITSRecPoint * clist[6];
  Float_t shared = GetNumberOfSharedClusters(besttrack,esdindex,list,clist);
  if (besttrack->GetConstrain()&&besttrack->GetChi2MIP(0)<kMaxChi2PerCluster[0]&&besttrack->GetChi2MIP(1)<kMaxChi2PerCluster[1]
      &&besttrack->GetChi2MIP(2)<kMaxChi2PerCluster[2]&&besttrack->GetChi2MIP(3)<kMaxChi2PerCluster[3]){ 
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
  // Don't sign clusters if not gold track
  //
  if (!besttrack->IsGoldPrimary()) return besttrack;
  if (besttrack->GetESDtrack()->GetKinkIndex(0)!=0) return besttrack;   //track belong to kink
  //
  if (fConstraint[fPass]){
    //
    // sign clusters
    //
    Float_t *ny = GetNy(esdindex), *nz = GetNz(esdindex);
    for (Int_t i=0;i<6;i++){
      Int_t index = besttrack->GetClIndex(i);
      if (index<=0) continue; 
      Int_t ilayer =  (index & 0xf0000000) >> 28;
      if (besttrack->GetSigmaY(ilayer)<0.00000000001) continue;
      AliITSRecPoint *c = (AliITSRecPoint*)GetCluster(index);     
      if (!c) continue;
      if (ilayer>3&&c->GetNy()+c->GetNz()>6) continue;
      if ( (c->GetNy()+c->GetNz() )> ny[i]+nz[i]+0.7) continue; //shared track
      if (  c->GetNz()> nz[i]+0.7) continue; //shared track
      if ( ilayer>2&& besttrack->GetNormQ(ilayer)/besttrack->GetExpQ()>1.5) continue;
      //if (  c->GetNy()> ny[i]+0.7) continue; //shared track

      Bool_t cansign = kTRUE;
      for (Int_t itrack=0;itrack<entries; itrack++){
	AliITStrackMI * track = (AliITStrackMI*)array->At(i);   
	if (!track) continue;
	if (track->GetChi2MIP(0)>besttrack->GetChi2MIP(0)+2.*shared+1.) break;
	if ( (track->GetClIndex(ilayer)>0) && (track->GetClIndex(ilayer)!=besttrack->GetClIndex(ilayer))){
	  cansign = kFALSE;
	  break;
	}
      }
      if (cansign){
	if (TMath::Abs(besttrack->GetDy(ilayer)/besttrack->GetSigmaY(ilayer))>3.) continue;
	if (TMath::Abs(besttrack->GetDz(ilayer)/besttrack->GetSigmaZ(ilayer))>3.) continue;    
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
  //

  Int_t nentries = itsTracks.GetEntriesFast();
  for (Int_t i=0;i<nentries;i++){
    AliITStrackMI* track = (AliITStrackMI*)itsTracks.At(i);
    if (!track) continue;
    TObjArray * array = (TObjArray*) fTrackHypothesys.At(i);
    if (!array) continue;
    if (array->GetEntriesFast()<=0) continue;
    //
    AliITStrackMI* longtrack=0;
    Float_t minn=0;
    Float_t maxchi2=1000;
    for (Int_t j=0;j<array->GetEntriesFast();j++){
      AliITStrackMI* track = (AliITStrackMI*)array->At(j);
      if (!track) continue;
      if (track->GetGoldV0()) {
	longtrack = track;   //gold V0 track taken
	break;
      }
      if (track->GetNumberOfClusters()+track->GetNDeadZone()<minn) continue;
      Float_t chi2 = track->GetChi2MIP(0);
      if (fAfterV0){
	if (!track->GetGoldV0()&&track->GetConstrain()==kFALSE) chi2+=5;
      }
      if (track->GetNumberOfClusters()+track->GetNDeadZone()>minn) maxchi2 = track->GetChi2MIP(0);       
      //
      if (chi2 > maxchi2) continue;
      minn= track->GetNumberOfClusters()+track->GetNDeadZone();
      maxchi2 = chi2;
      longtrack=track;
    }    
    //
    //
    //
    AliITStrackMI * besttrack = (AliITStrackMI*)array->At(0);
    if (!longtrack) {longtrack = besttrack;}
    else besttrack= longtrack;
    //
    if (besttrack){
      Int_t list[6];
      AliITSRecPoint * clist[6];
      Float_t shared = GetNumberOfSharedClusters(longtrack,i,list,clist);
      //
      track->SetNUsed(shared);      
      track->SetNSkipped(besttrack->GetNSkipped());
      track->SetChi2MIP(0,besttrack->GetChi2MIP(0));
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
      if (besttrack&&fAfterV0){
	UpdateESDtrack(besttrack,AliESDtrack::kITSin);
      }
      if (besttrack&&fConstraint[fPass]) 
      	UpdateESDtrack(besttrack,AliESDtrack::kITSin);
      //if (besttrack&&besttrack->fConstrain) 
      //	UpdateESDtrack(besttrack,AliESDtrack::kITSin);
      if (besttrack->GetChi2MIP(0)+besttrack->GetNUsed()>1.5){
	if ( (TMath::Abs(besttrack->GetD(0))>0.1) && fConstraint[fPass]) {
	  track->SetReconstructed(kFALSE);
	}
	if ( (TMath::Abs(besttrack->GetD(1))>0.1) && fConstraint[fPass]){
	  track->SetReconstructed(kFALSE);
	}
      }       

    }    
  }
} 


void AliITStrackerMI::CookLabel(AliITStrackMI *track,Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t tpcLabel=-1; 
     
  if ( track->GetESDtrack())   tpcLabel =  TMath::Abs(track->GetESDtrack()->GetTPCLabel());

   track->SetChi2MIP(9,0);
   Int_t nwrong=0;
   for (Int_t i=0;i<track->GetNumberOfClusters();i++){
     Int_t cindex = track->GetClusterIndex(i);
     Int_t l=(cindex & 0xf0000000) >> 28;
     AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(cindex);
     Int_t isWrong=1;
     for (Int_t ind=0;ind<3;ind++){
       if (tpcLabel>0)
	 if (cl->GetLabel(ind)==tpcLabel) isWrong=0;
     }
     track->SetChi2MIP(9,track->GetChi2MIP(9)+isWrong*(2<<l));
     nwrong+=isWrong;
   }
   track->SetFakeRatio(double(nwrong)/double(track->GetNumberOfClusters()));
   if (tpcLabel>0){
     if (track->GetFakeRatio()>wrong) track->SetLabel(-tpcLabel);
     else
       track->SetLabel(tpcLabel);
   }
   
}



void AliITStrackerMI::CookdEdx(AliITStrackMI* track)
{
  //
  //
  //  Int_t list[6];
  //AliITSRecPoint * clist[6];
  //  Int_t shared = GetNumberOfSharedClusters(track,index,list,clist);
  Float_t dedx[4];
  Int_t accepted=0;
  track->SetChi2MIP(9,0);
  for (Int_t i=0;i<track->GetNumberOfClusters();i++){
    Int_t cindex = track->GetClusterIndex(i);
    Int_t l=(cindex & 0xf0000000) >> 28;
    AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(cindex);
    Int_t lab = TMath::Abs(track->GetESDtrack()->GetTPCLabel());
    Int_t isWrong=1;
    for (Int_t ind=0;ind<3;ind++){
      if (cl->GetLabel(ind)==lab) isWrong=0;
    }
    track->SetChi2MIP(9,track->GetChi2MIP(9)+isWrong*(2<<l));
    if (l<2) continue;
    //if (l>3 && (cl->GetNy()>4) || (cl->GetNz()>4)) continue;  //shared track
    //if (l>3&& !(cl->GetType()==1||cl->GetType()==10)) continue;
    //if (l<4&& !(cl->GetType()==1)) continue;   
    dedx[accepted]= track->GetSampledEdx(i);
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
  if (fCoeficients) delete []fCoeficients;
  fCoeficients = new Float_t[ntracks*48];
  for (Int_t i=0;i<ntracks*48;i++) fCoeficients[i]=-1.;
}


Double_t AliITStrackerMI::GetPredictedChi2MI(AliITStrackMI* track, const AliITSRecPoint *cluster,Int_t layer) 
{
  //
  //
  //
  Float_t erry,errz;
  Float_t theta = track->GetTgl();
  Float_t phi   = track->GetSnp();
  phi = TMath::Sqrt(phi*phi/(1.-phi*phi));
  GetError(layer,cluster,theta,phi,track->GetExpQ(),erry,errz);
  Double_t chi2 = track->GetPredictedChi2MI(cluster->GetY(),cluster->GetZ(),erry,errz);
  Float_t ny,nz;
  GetNTeor(layer,cluster, theta,phi,ny,nz);  
  Double_t delta = cluster->GetNy()+cluster->GetNz()-nz-ny;
  if (delta>1){
    chi2+=0.5*TMath::Min(delta/2,2.);
    chi2+=2.*cluster->GetDeltaProbability();
  }
  //
  track->SetNy(layer,ny);
  track->SetNz(layer,nz);
  track->SetSigmaY(layer,erry);
  track->SetSigmaZ(layer, errz);
  //track->fNormQ[layer] = cluster->GetQ()/TMath::Sqrt(1+theta*theta+phi*phi);
  track->SetNormQ(layer,cluster->GetQ()/TMath::Sqrt((1.+ track->GetTgl()*track->GetTgl())/(1.- track->GetSnp()*track->GetSnp())));
  return chi2;

}

Int_t    AliITStrackerMI::UpdateMI(AliITStrackMI* track, const AliITSRecPoint* cl,Double_t chi2,Int_t index) const 
{
  //
  //
  //
  Int_t layer = (index & 0xf0000000) >> 28;
  track->SetClIndex(layer, index);
  if ( (layer>1) &&track->GetNormQ(layer)/track->GetExpQ()<0.5 ) {
    chi2+= (0.5-track->GetNormQ(layer)/track->GetExpQ())*10.;
    track->SetdEdxMismatch(track->GetdEdxMismatch()+(0.5-track->GetNormQ(layer)/track->GetExpQ())*10.);
  }
  return track->UpdateMI(cl->GetY(),cl->GetZ(),track->GetSigmaY(layer),track->GetSigmaZ(layer),chi2,index);
}

void AliITStrackerMI::GetNTeor(Int_t layer, const AliITSRecPoint* /*cl*/, Float_t theta, Float_t phi, Float_t &ny, Float_t &nz)
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



Int_t AliITStrackerMI::GetError(Int_t layer, const AliITSRecPoint*cl, Float_t theta, Float_t phi,Float_t expQ, Float_t &erry, Float_t &errz)
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
    //factor 1.8 appears in new simulation
    //
    Float_t scale=1.8;
    if (cl->GetNy()==100||cl->GetNz()==100){
      erry = 0.004*scale;
      errz = 0.2*scale;
      return 100;
    }
    if (cl->GetNy()+cl->GetNz()>12){
      erry = 0.06*scale;
      errz = 0.57*scale;
      return 100;
    }
    Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
    Float_t chargematch = TMath::Max(double(normq/expQ),2.);
    //
    if (cl->GetType()==1 || cl->GetType()==10 ){     							       
      if (chargematch<1.0 || (cl->GetNy()+cl->GetNz()<nz+ny+0.5)){
	errz = 0.043*scale;
	erry = 0.00094*scale;
	return 101;
      }
      if (cl->GetNy()+cl->GetNz()<nz+ny+1.2){
	errz = 0.06*scale;
	erry =0.0013*scale;
	return 102;
      }
      erry = 0.0027*scale;
      errz = TMath::Min(0.028*(chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.15)*scale;
      return 103;
    }
    if (cl->GetType()==2 || cl->GetType()==11 ){ 
      erry = TMath::Min(0.0010*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.05)*scale;
      errz = TMath::Min(0.025*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.5)*scale;
      return 104;
    }
    
    if (cl->GetType()>100 ){     							       
      if ((chargematch+cl->GetNy()+cl->GetNz()-nz-ny<1.5)){
	errz = 0.05*scale;
	erry = 0.00096*scale;
	return 105;
      }
      if (cl->GetNy()+cl->GetNz()-nz-ny<1){
	errz = 0.10*scale;
	erry = 0.0025*scale;
	return 106;
      }

      errz = TMath::Min(0.05*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.4)*scale;
      erry = TMath::Min(0.003*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.05)*scale;
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



void   AliITStrackerMI::GetDCASigma(AliITStrackMI* track, Float_t & sigmarfi, Float_t &sigmaz)
{
  //
  //DCA sigmas parameterization
  //to be paramterized using external parameters in future 
  //
  // 
  sigmarfi = 0.004+1.4 *TMath::Abs(track->GetC())+332.*track->GetC()*track->GetC();
  sigmaz   = 0.011+4.37*TMath::Abs(track->GetC());
}


void AliITStrackerMI::SignDeltas( TObjArray *ClusterArray, Float_t vz)
{
  //
  //  
  Int_t entries = ClusterArray->GetEntriesFast();
  if (entries<4) return;
  AliITSRecPoint* cluster = (AliITSRecPoint*)ClusterArray->At(0);
  Int_t layer = cluster->GetLayer();
  if (layer>1) return;
  Int_t index[10000];
  Int_t ncandidates=0;
  Float_t r = (layer>0)? 7:4;
  // 
  for (Int_t i=0;i<entries;i++){
    AliITSRecPoint* cl0 = (AliITSRecPoint*)ClusterArray->At(i);
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
    AliITSRecPoint* cl0 = (AliITSRecPoint*)ClusterArray->At(index[i]);
    if (cl0->GetDeltaProbability()>0.8) continue;
    // 
    Int_t ncl = 0;
    Float_t y[100],z[100],sumy,sumz,sumy2, sumyz, sumw;
    sumy=sumz=sumy2=sumyz=sumw=0.0;
    for (Int_t j=0;j<ncandidates;j++){
      if (i==j) continue;
      AliITSRecPoint* cl1 = (AliITSRecPoint*)ClusterArray->At(index[j]);
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


void AliITStrackerMI::UpdateESDtrack(AliITStrackMI* track, ULong_t flags) const
{
  //
  //
  track->UpdateESDtrack(flags);
  AliITStrackMI * oldtrack = (AliITStrackMI*)(track->GetESDtrack()->GetITStrack());
  if (oldtrack) delete oldtrack; 
  track->GetESDtrack()->SetITStrack(new AliITStrackMI(*track));
  if (TMath::Abs(track->GetDnorm(1))<0.000000001){
    printf("Problem\n");
  }
}



Int_t AliITStrackerMI::GetNearestLayer(const Double_t *xr) const{
  //
  //Get nearest upper layer close to the point xr.
  // rough approximation 
  //
  const Float_t kRadiuses[6]={4,6.5,15.03,24.,38.5,43.7};
  Float_t radius = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
  Int_t res =6;
  for (Int_t i=0;i<6;i++){
    if (radius<kRadiuses[i]){
      res =i;
      break;
    }
  }
  return res;
}


void AliITStrackerMI::UpdateTPCV0(AliESD *event){
  //
  //try to update, or reject TPC  V0s
  //
  Int_t nv0s = event->GetNumberOfV0s();
  Int_t nitstracks = fTrackHypothesys.GetEntriesFast();

  for (Int_t i=0;i<nv0s;i++){
    AliESDv0 * vertex = event->GetV0(i);
    Int_t ip = vertex->GetIndex(0);
    Int_t im = vertex->GetIndex(1);
    //
    TObjArray * arrayp = (ip<nitstracks) ? (TObjArray*)fTrackHypothesys.At(ip):0;
    TObjArray * arraym = (im<nitstracks) ? (TObjArray*)fTrackHypothesys.At(im):0;
    AliITStrackMI * trackp = (arrayp!=0) ? (AliITStrackMI*)arrayp->At(0):0;
    AliITStrackMI * trackm = (arraym!=0) ? (AliITStrackMI*)arraym->At(0):0;
    //
    //
    if (trackp){
      if (trackp->GetNumberOfClusters()+trackp->GetNDeadZone()>5.5){
	if (trackp->GetConstrain()&&trackp->GetChi2MIP(0)<3) vertex->SetStatus(-100);
	if (!trackp->GetConstrain()&&trackp->GetChi2MIP(0)<2) vertex->SetStatus(-100); 
      }
    }

    if (trackm){
      if (trackm->GetNumberOfClusters()+trackm->GetNDeadZone()>5.5){
	if (trackm->GetConstrain()&&trackm->GetChi2MIP(0)<3) vertex->SetStatus(-100);
	if (!trackm->GetConstrain()&&trackm->GetChi2MIP(0)<2) vertex->SetStatus(-100); 
      }
    }
    if (vertex->GetStatus()==-100) continue;
    //
    Double_t xrp[3]; vertex->GetXYZ(xrp[0],xrp[1],xrp[2]);  //I.B.
    Int_t clayer = GetNearestLayer(xrp);                    //I.B.
    vertex->SetNBefore(clayer);        //
    vertex->SetChi2Before(9*clayer);   //
    vertex->SetNAfter(6-clayer);       //
    vertex->SetChi2After(0);           //
    //
    if (clayer >1 ){ // calculate chi2 before vertex
      Float_t chi2p = 0, chi2m=0;  
      //
      if (trackp){
	for (Int_t ilayer=0;ilayer<clayer;ilayer++){
	  if (trackp->GetClIndex(ilayer)>0){
	    chi2p+=trackp->GetDy(ilayer)*trackp->GetDy(ilayer)/(trackp->GetSigmaY(ilayer)*trackp->GetSigmaY(ilayer))+
	      trackp->GetDz(ilayer)*trackp->GetDz(ilayer)/(trackp->GetSigmaZ(ilayer)*trackp->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2p+=9;
	  }
	}
      }else{
	chi2p = 9*clayer;
      }
      //
      if (trackm){
	for (Int_t ilayer=0;ilayer<clayer;ilayer++){
	  if (trackm->GetClIndex(ilayer)>0){
	    chi2m+=trackm->GetDy(ilayer)*trackm->GetDy(ilayer)/(trackm->GetSigmaY(ilayer)*trackm->GetSigmaY(ilayer))+
	      trackm->GetDz(ilayer)*trackm->GetDz(ilayer)/(trackm->GetSigmaZ(ilayer)*trackm->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2m+=9;
	  }
	}
      }else{
	chi2m = 9*clayer;
      }
      vertex->SetChi2Before(TMath::Min(chi2p,chi2m));
      if (TMath::Min(chi2p,chi2m)/Float_t(clayer)<4) vertex->SetStatus(-10);  // track exist before vertex
    }
    
    if (clayer < 5 ){ // calculate chi2 after vertex
      Float_t chi2p = 0, chi2m=0;  
      //
      if (trackp&&TMath::Abs(trackp->GetTgl())<1.){
	for (Int_t ilayer=clayer;ilayer<6;ilayer++){
	  if (trackp->GetClIndex(ilayer)>0){
	    chi2p+=trackp->GetDy(ilayer)*trackp->GetDy(ilayer)/(trackp->GetSigmaY(ilayer)*trackp->GetSigmaY(ilayer))+
	      trackp->GetDz(ilayer)*trackp->GetDz(ilayer)/(trackp->GetSigmaZ(ilayer)*trackp->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2p+=9;
	  }
	}
      }else{
	chi2p = 0;
      }
      //
      if (trackm&&TMath::Abs(trackm->GetTgl())<1.){
	for (Int_t ilayer=clayer;ilayer<6;ilayer++){
	  if (trackm->GetClIndex(ilayer)>0){
	    chi2m+=trackm->GetDy(ilayer)*trackm->GetDy(ilayer)/(trackm->GetSigmaY(ilayer)*trackm->GetSigmaY(ilayer))+
	      trackm->GetDz(ilayer)*trackm->GetDz(ilayer)/(trackm->GetSigmaZ(ilayer)*trackm->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2m+=9;
	  }
	}
      }else{
	chi2m = 0;
      }
      vertex->SetChi2After(TMath::Max(chi2p,chi2m));
      if (TMath::Max(chi2m,chi2p)/Float_t(6-clayer)>9) vertex->SetStatus(-20);  // track not found in ITS
    }
  }
  //
}



void  AliITStrackerMI::FindV02(AliESD *event)
{
  //
  // V0 finder
  //
  //  Cuts on DCA -  R dependent
  //          max distance DCA between 2 tracks cut 
  //          maxDist = TMath::Min(kMaxDist,kMaxDist0+pvertex->GetRr()*kMaxDist);
  //
  const Float_t kMaxDist0      = 0.1;    
  const Float_t kMaxDist1      = 0.1;     
  const Float_t kMaxDist       = 1;
  const Float_t kMinPointAngle  = 0.85;
  const Float_t kMinPointAngle2 = 0.99;
  const Float_t kMinR           = 0.5;
  const Float_t kMaxR           = 220;
  //const Float_t kCausality0Cut   = 0.19;
  //const Float_t kLikelihood01Cut = 0.25;
  //const Float_t kPointAngleCut   = 0.9996;
  const Float_t kCausality0Cut   = 0.19;
  const Float_t kLikelihood01Cut = 0.45;
  const Float_t kLikelihood1Cut  = 0.5;
  const Float_t kCombinedCut     = 0.55;

  //
  //
  TTreeSRedirector &cstream = *fDebugStreamer;
  Int_t ntracks    = event->GetNumberOfTracks(); 
  Int_t nitstracks = fTrackHypothesys.GetEntriesFast();
  fOriginal.Expand(ntracks);
  fTrackHypothesys.Expand(ntracks);
  fBestHypothesys.Expand(ntracks);
  //
  AliHelix * helixes   = new AliHelix[ntracks+2];
  TObjArray trackarray(ntracks+2);     //array with tracks - with vertex constrain
  TObjArray trackarrayc(ntracks+2);    //array of "best    tracks" - without vertex constrain
  TObjArray trackarrayl(ntracks+2);    //array of "longest tracks" - without vertex constrain
  Bool_t * forbidden   = new Bool_t [ntracks+2];
  Int_t   *itsmap      = new Int_t  [ntracks+2];
  Float_t *dist        = new Float_t[ntracks+2];
  Float_t *normdist0   = new Float_t[ntracks+2];
  Float_t *normdist1   = new Float_t[ntracks+2];
  Float_t *normdist    = new Float_t[ntracks+2];
  Float_t *norm        = new Float_t[ntracks+2];
  Float_t *maxr        = new Float_t[ntracks+2];
  Float_t *minr        = new Float_t[ntracks+2];
  Float_t *minPointAngle= new Float_t[ntracks+2];
  //
  AliV0 *pvertex      = new AliV0;
  AliITStrackMI * dummy= new AliITStrackMI;
  dummy->SetLabel(0);
  AliITStrackMI  trackat0;    //temporary track for DCA calculation
  //
  Float_t primvertex[3]={GetX(),GetY(),GetZ()};
  //
  // make its -  esd map
  //
  for (Int_t itrack=0;itrack<ntracks+2;itrack++) {
    itsmap[itrack]        = -1;
    forbidden[itrack]     = kFALSE;
    maxr[itrack]          = kMaxR;
    minr[itrack]          = kMinR;
    minPointAngle[itrack] = kMinPointAngle;
  }
  for (Int_t itrack=0;itrack<nitstracks;itrack++){
    AliITStrackMI * original =   (AliITStrackMI*)(fOriginal.At(itrack));
    Int_t           esdindex =   original->GetESDtrack()->GetID();
    itsmap[esdindex]         =   itrack;
  }
  //
  // create its tracks from esd tracks if not done before
  //
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    if (itsmap[itrack]>=0) continue;
    AliITStrackMI * tpctrack = new AliITStrackMI(*(event->GetTrack(itrack)));
    //tpctrack->fD[0] = tpctrack->GetD(GetX(),GetY());
    //tpctrack->fD[1] = tpctrack->GetZat(GetX())-GetZ(); 
    tpctrack->GetDZ(GetX(),GetY(),GetZ(),tpctrack->GetDP());   //I.B.
    if (tpctrack->GetD(0)<20 && tpctrack->GetD(1)<20){
      // tracks which can reach inner part of ITS
      // propagate track to outer its volume - with correction for material
      CorrectForDeadZoneMaterial(tpctrack);  
    }
    itsmap[itrack] = nitstracks;
    fOriginal.AddAt(tpctrack,nitstracks);
    nitstracks++;
  }
  //
  // fill temporary arrays
  //
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    AliESDtrack *  esdtrack = event->GetTrack(itrack);
    Int_t          itsindex = itsmap[itrack];
    AliITStrackMI *original = (AliITStrackMI*)fOriginal.At(itsmap[itrack]);
    if (!original) continue;
    AliITStrackMI *bestConst  = 0;
    AliITStrackMI *bestLong   = 0;
    AliITStrackMI *best       = 0;    
    //
    //
    TObjArray * array    = (TObjArray*)  fTrackHypothesys.At(itsindex);
    Int_t       hentries = (array==0) ?  0 : array->GetEntriesFast();
    // Get best track with vertex constrain
    for (Int_t ih=0;ih<hentries;ih++){
      AliITStrackMI * trackh = (AliITStrackMI*)array->At(ih);
      if (!trackh->GetConstrain()) continue;
      if (!bestConst) bestConst = trackh;
      if (trackh->GetNumberOfClusters()>5.0){
	bestConst  = trackh;                         // full track -  with minimal chi2
	break;
      }
      if (trackh->GetNumberOfClusters()+trackh->GetNDeadZone()<=bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone())  continue;      
      bestConst = trackh;
      break;
    }
    // Get best long track without vertex constrain and best track without vertex constrain
    for (Int_t ih=0;ih<hentries;ih++){
      AliITStrackMI * trackh = (AliITStrackMI*)array->At(ih);
      if (trackh->GetConstrain()) continue;
      if (!best)     best     = trackh;
      if (!bestLong) bestLong = trackh;
      if (trackh->GetNumberOfClusters()>5.0){
	bestLong  = trackh;                         // full track -  with minimal chi2
	break;
      }
      if (trackh->GetNumberOfClusters()+trackh->GetNDeadZone()<=bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone())  continue;      
      bestLong = trackh;	
    }
    if (!best) {
      best     = original;
      bestLong = original;
    }
    //I.B. trackat0 = *bestLong;
    new (&trackat0) AliITStrackMI(*bestLong);
    Double_t xx,yy,zz,alpha; 
    bestLong->GetGlobalXYZat(bestLong->GetX(),xx,yy,zz);
    alpha = TMath::ATan2(yy,xx);    
    trackat0.Propagate(alpha,0);      
    // calculate normalized distances to the vertex 
    //
    Float_t ptfac  = (1.+100.*TMath::Abs(trackat0.GetC()));
    if ( bestLong->GetNumberOfClusters()>3 ){      
      dist[itsindex]      = trackat0.GetY();
      norm[itsindex]      = ptfac*TMath::Sqrt(trackat0.GetSigmaY2());
      normdist0[itsindex] = TMath::Abs(trackat0.GetY()/norm[itsindex]);
      normdist1[itsindex] = TMath::Abs((trackat0.GetZ()-primvertex[2])/(ptfac*TMath::Sqrt(trackat0.GetSigmaZ2())));
      normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
      if (!bestConst){
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<6) normdist[itsindex]*=2.;
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<5) normdist[itsindex]*=2.;
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<4) normdist[itsindex]*=2.;
      }else{
	if (bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()<6) normdist[itsindex]*=1.5;
	if (bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()<5) normdist[itsindex]*=1.5;
      }
    }
    else{      
      if (bestConst&&bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>4.5){
	dist[itsindex] = bestConst->GetD(0);
	norm[itsindex] = bestConst->GetDnorm(0);
	normdist0[itsindex] = TMath::Abs(bestConst->GetD(0)/norm[itsindex]);
	normdist1[itsindex] = TMath::Abs(bestConst->GetD(0)/norm[itsindex]);
	normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
      }else{
	dist[itsindex]      = trackat0.GetY();
	norm[itsindex]      = ptfac*TMath::Sqrt(trackat0.GetSigmaY2());
	normdist0[itsindex] = TMath::Abs(trackat0.GetY()/norm[itsindex]);
	normdist1[itsindex] = TMath::Abs((trackat0.GetZ()-primvertex[2])/(ptfac*TMath::Sqrt(trackat0.GetSigmaZ2())));
	normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
	if (TMath::Abs(trackat0.GetTgl())>1.05){
	  if (normdist[itsindex]<3) forbidden[itsindex]=kTRUE;
	  if (normdist[itsindex]>3) {
	    minr[itsindex] = TMath::Max(Float_t(40.),minr[itsindex]);
	  }
	}
      }
    }
    //
    //-----------------------------------------------------------
    //Forbid primary track candidates - 
    //
    //treetr->SetAlias("forbidden0","Tr0.fN<4&&Tr1.fN+Tr1.fNDeadZone>4.5");
    //treetr->SetAlias("forbidden1","ND<3&&Tr1.fN+Tr1.fNDeadZone>5.5");
    //treetr->SetAlias("forbidden2","ND<2&&Tr1.fClIndex[0]>0&&Tr1.fClIndex[0]>0");
    //treetr->SetAlias("forbidden3","ND<1&&Tr1.fClIndex[0]>0");
    //treetr->SetAlias("forbidden4","ND<4&&Tr1.fNormChi2[0]<2");
    //treetr->SetAlias("forbidden5","ND<5&&Tr1.fNormChi2[0]<1");
    //-----------------------------------------------------------
    if (bestConst){
      if (bestLong->GetNumberOfClusters()<4       && bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>4.5)               forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<3 && bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>5.5)               forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<2 && bestConst->GetClIndex(0)>0 && bestConst->GetClIndex(1)>0 ) forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<1 && bestConst->GetClIndex(0)>0)                              forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<4 && bestConst->GetNormChi2(0)<2)                             forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<5 && bestConst->GetNormChi2(0)<1)                             forbidden[itsindex]=kTRUE;      
      if (bestConst->GetNormChi2(0)<2.5) {
	minPointAngle[itsindex]= 0.9999;
	maxr[itsindex]         = 10;
      }
    }
    //
    //forbid daughter kink candidates
    //
    if (esdtrack->GetKinkIndex(0)>0) forbidden[itsindex] = kTRUE;
    Bool_t isElectron = kTRUE;
    Bool_t isProton   = kTRUE;
    Double_t pid[5];
    esdtrack->GetESDpid(pid);
    for (Int_t i=1;i<5;i++){
      if (pid[0]<pid[i]) isElectron= kFALSE;
      if (pid[4]<pid[i]) isProton= kFALSE;
    }
    if (isElectron){
      forbidden[itsindex]=kFALSE;	
      normdist[itsindex]*=-1;
    }
    if (isProton){
      if (normdist[itsindex]>2) forbidden[itsindex]=kFALSE;	
      normdist[itsindex]*=-1;
    }

    //
    // Causality cuts in TPC volume
    //
    if (esdtrack->GetTPCdensity(0,10) >0.6)  maxr[itsindex] = TMath::Min(Float_t(110),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(10,30)>0.6)  maxr[itsindex] = TMath::Min(Float_t(120),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(20,40)>0.6)  maxr[itsindex] = TMath::Min(Float_t(130),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(30,50)>0.6)  maxr[itsindex] = TMath::Min(Float_t(140),maxr[itsindex]);
    //
    if (esdtrack->GetTPCdensity(0,60)<0.4&&bestLong->GetNumberOfClusters()<3) minr[itsindex]=100;    
    //
    //
    if (kFALSE){
      cstream<<"Track"<<
	"Tr0.="<<best<<
	"Tr1.="<<((bestConst)? bestConst:dummy)<<
	"Tr2.="<<bestLong<<
	"Tr3.="<<&trackat0<<
	"Esd.="<<esdtrack<<
	"Dist="<<dist[itsindex]<<
	"ND0="<<normdist0[itsindex]<<
	"ND1="<<normdist1[itsindex]<<
	"ND="<<normdist[itsindex]<<
	"Pz="<<primvertex[2]<<
	"Forbid="<<forbidden[itsindex]<<
	"\n";
      //
    }
    trackarray.AddAt(best,itsindex);
    trackarrayc.AddAt(bestConst,itsindex);
    trackarrayl.AddAt(bestLong,itsindex);
    new (&helixes[itsindex]) AliHelix(*best);
  }
  //
  //
  //
  // first iterration of V0 finder
  //
  for (Int_t iesd0=0;iesd0<ntracks;iesd0++){
    Int_t itrack0 = itsmap[iesd0];
    if (forbidden[itrack0]) continue;
    AliITStrackMI * btrack0 = (AliITStrackMI*)trackarray.At(itrack0);
    if (!btrack0) continue;    
    if (btrack0->GetSign()>0) continue;
    AliITStrackMI *trackc0 = (AliITStrackMI*)trackarrayc.At(itrack0);
    //
    for (Int_t iesd1=0;iesd1<ntracks;iesd1++){
      Int_t itrack1 = itsmap[iesd1];
      if (forbidden[itrack1]) continue;

      AliITStrackMI * btrack1 = (AliITStrackMI*)trackarray.At(itrack1); 
      if (!btrack1) continue;
      if (btrack1->GetSign()<0) continue;
      Bool_t isGold = kFALSE;
      if (TMath::Abs(TMath::Abs(btrack0->GetLabel())-TMath::Abs(btrack1->GetLabel()))==1){
	isGold = kTRUE;
      }
      AliITStrackMI *trackc1 = (AliITStrackMI*)trackarrayc.At(itrack1);
      AliHelix &h1 = helixes[itrack0];
      AliHelix &h2 = helixes[itrack1];
      //
      // find linear distance
      Double_t rmin =0;
      //
      //
      //
      Double_t phase[2][2],radius[2];
      Int_t  points = h1.GetRPHIintersections(h2, phase, radius);
      if    (points==0)  continue;
      Double_t delta[2]={1000000,1000000};        
      rmin = radius[0];
      h1.ParabolicDCA(h2,phase[0][0],phase[0][1],radius[0],delta[0]);
      if (points==2){    
	if (radius[1]<rmin) rmin = radius[1];
	h1.ParabolicDCA(h2,phase[1][0],phase[1][1],radius[1],delta[1]);
      }
      rmin = TMath::Sqrt(rmin);
      Double_t distance = 0;
      Double_t radiusC  = 0;
      Int_t    iphase   = 0;
      if (points==1 || delta[0]<delta[1]){
	distance = TMath::Sqrt(delta[0]);
	radiusC  = TMath::Sqrt(radius[0]);
      }else{
	distance = TMath::Sqrt(delta[1]);
	radiusC  = TMath::Sqrt(radius[1]);
	iphase=1;
      }
      if (radiusC<TMath::Max(minr[itrack0],minr[itrack1]))    continue;
      if (radiusC>TMath::Min(maxr[itrack0],maxr[itrack1]))     continue; 
      Float_t maxDist  = TMath::Min(kMaxDist,Float_t(kMaxDist0+radiusC*kMaxDist1));      
      if (distance>maxDist) continue;
      Float_t pointAngle = h1.GetPointAngle(h2,phase[iphase],primvertex);
      if (pointAngle<TMath::Max(minPointAngle[itrack0],minPointAngle[itrack1])) continue;
      //
      //
      //      Double_t distance = TestV0(h1,h2,pvertex,rmin);
      //
      //       if (distance>maxDist)           continue;
      //       if (pvertex->GetRr()<kMinR)     continue;
      //       if (pvertex->GetRr()>kMaxR)     continue;
      AliITStrackMI * track0=btrack0;
      AliITStrackMI * track1=btrack1;
      //      if (pvertex->GetRr()<3.5){  
      if (radiusC<3.5){  
	//use longest tracks inside the pipe
	track0 = (AliITStrackMI*)trackarrayl.At(itrack0);
	track1 = (AliITStrackMI*)trackarrayl.At(itrack1);
      }      
      //
      //
      pvertex->SetParamN(*track0);
      pvertex->SetParamP(*track1);
      pvertex->Update(primvertex);
      pvertex->SetClusters(track0->ClIndex(),track1->ClIndex());  // register clusters

      if (pvertex->GetRr()<kMinR) continue;
      if (pvertex->GetRr()>kMaxR) continue;
      if (pvertex->GetV0CosineOfPointingAngle()<kMinPointAngle) continue;
//Bo:      if (pvertex->GetDist2()>maxDist) continue;
      if (pvertex->GetDcaV0Daughters()>maxDist) continue;
//Bo:        pvertex->SetLab(0,track0->GetLabel());
//Bo:        pvertex->SetLab(1,track1->GetLabel());
      pvertex->SetIndex(0,track0->GetESDtrack()->GetID());
      pvertex->SetIndex(1,track1->GetESDtrack()->GetID());
      //      
      AliITStrackMI * htrackc0 = trackc0 ? trackc0:dummy;      
      AliITStrackMI * htrackc1 = trackc1 ? trackc1:dummy;

      //
      //
      TObjArray * array0b     = (TObjArray*)fBestHypothesys.At(itrack0);
      if (!array0b&&pvertex->GetRr()<40 && TMath::Abs(track0->GetTgl())<1.1) 
	FollowProlongationTree((AliITStrackMI*)fOriginal.At(itrack0),itrack0, kFALSE);
      TObjArray * array1b    = (TObjArray*)fBestHypothesys.At(itrack1);
      if (!array1b&&pvertex->GetRr()<40 && TMath::Abs(track1->GetTgl())<1.1) 
	FollowProlongationTree((AliITStrackMI*)fOriginal.At(itrack1),itrack1, kFALSE);
      //
      AliITStrackMI * track0b = (AliITStrackMI*)fOriginal.At(itrack0);       
      AliITStrackMI * track1b = (AliITStrackMI*)fOriginal.At(itrack1);
      AliITStrackMI * track0l = (AliITStrackMI*)fOriginal.At(itrack0);       
      AliITStrackMI * track1l = (AliITStrackMI*)fOriginal.At(itrack1);
      
      Float_t minchi2before0=16;
      Float_t minchi2before1=16;
      Float_t minchi2after0 =16;
      Float_t minchi2after1 =16;
      Double_t xrp[3]; pvertex->GetXYZ(xrp[0],xrp[1],xrp[2]);  //I.B.
      Int_t maxLayer = GetNearestLayer(xrp);                   //I.B.
      
      if (array0b) for (Int_t i=0;i<5;i++){
	// best track after vertex
	AliITStrackMI * btrack = (AliITStrackMI*)array0b->At(i);
	if (!btrack) continue;
	if (btrack->GetNumberOfClusters()>track0l->GetNumberOfClusters()) track0l = btrack;     
	//	if (btrack->fX<pvertex->GetRr()-2.-0.5/(0.1+pvertex->GetAnglep()[2])) {
	if (btrack->GetX()<pvertex->GetRr()-2.) {
	  if ( (maxLayer>i+2|| (i==0)) && btrack->GetNumberOfClusters()==(6-i)&&i<3){
	    Float_t sumchi2= 0;
	    Float_t sumn   = 0;
	    if (maxLayer<3){   // take prim vertex as additional measurement
	      if (normdist[itrack0]>0 && htrackc0){
		sumchi2 += TMath::Min((3.-maxLayer)*normdist[itrack0]*normdist[itrack0],16.);
	      }else{
		sumchi2 +=  TMath::Min((3.-maxLayer)*(3*normdist[itrack0]*normdist[itrack0]+3.),16.);
	      }
	      sumn    +=  3-maxLayer;
	    }
	    for (Int_t ilayer=i;ilayer<maxLayer;ilayer++){
	      sumn+=1.;	      
	      if (!btrack->GetClIndex(ilayer)){
		sumchi2+=25;
		continue;
	      }else{
		Int_t c=( btrack->GetClIndex(ilayer) & 0x0fffffff);
		for (Int_t itrack=0;itrack<4;itrack++){
		  if (fgLayers[ilayer].GetClusterTracks(itrack,c)>=0 && fgLayers[ilayer].GetClusterTracks(itrack,c)!=itrack0){
		    sumchi2+=18.;  //shared cluster
		    break;
		  }
		}
		sumchi2+=btrack->GetDy(ilayer)*btrack->GetDy(ilayer)/(btrack->GetSigmaY(ilayer)*btrack->GetSigmaY(ilayer));
		sumchi2+=btrack->GetDz(ilayer)*btrack->GetDz(ilayer)/(btrack->GetSigmaZ(ilayer)*btrack->GetSigmaZ(ilayer));	       
	      }
	    }
	    sumchi2/=sumn;
	    if (sumchi2<minchi2before0) minchi2before0=sumchi2; 
	  }
	  continue;   //safety space - Geo manager will give exact layer
	}
	track0b       = btrack;
	minchi2after0 = btrack->GetNormChi2(i);
	break;
      }
      if (array1b) for (Int_t i=0;i<5;i++){
	// best track after vertex
	AliITStrackMI * btrack = (AliITStrackMI*)array1b->At(i);
	if (!btrack) continue;
	if (btrack->GetNumberOfClusters()>track1l->GetNumberOfClusters()) track1l = btrack;     
	//	if (btrack->fX<pvertex->GetRr()-2-0.5/(0.1+pvertex->GetAnglep()[2])){
	if (btrack->GetX()<pvertex->GetRr()-2){
	  if ((maxLayer>i+2 || (i==0))&&btrack->GetNumberOfClusters()==(6-i)&&(i<3)){
	    Float_t sumchi2= 0;
	    Float_t sumn   = 0;
	    if (maxLayer<3){   // take prim vertex as additional measurement
	      if (normdist[itrack1]>0 && htrackc1){
		sumchi2 +=  TMath::Min((3.-maxLayer)*normdist[itrack1]*normdist[itrack1],16.);
	      }else{
		sumchi2 += TMath::Min((3.-maxLayer)*(3*normdist[itrack1]*normdist[itrack1]+3.),16.);
	      }
	      sumn    +=  3-maxLayer;
	    }
	    for (Int_t ilayer=i;ilayer<maxLayer;ilayer++){
	      sumn+=1.;
	      if (!btrack->GetClIndex(ilayer)){
		sumchi2+=30;
		continue;
	      }else{
		Int_t c=( btrack->GetClIndex(ilayer) & 0x0fffffff);
		for (Int_t itrack=0;itrack<4;itrack++){
		  if (fgLayers[ilayer].GetClusterTracks(itrack,c)>=0 && fgLayers[ilayer].GetClusterTracks(itrack,c)!=itrack1){
		    sumchi2+=18.;  //shared cluster
		    break;
		  }
		}
		sumchi2+=btrack->GetDy(ilayer)*btrack->GetDy(ilayer)/(btrack->GetSigmaY(ilayer)*btrack->GetSigmaY(ilayer));
		sumchi2+=btrack->GetDz(ilayer)*btrack->GetDz(ilayer)/(btrack->GetSigmaZ(ilayer)*btrack->GetSigmaZ(ilayer));	       
	      }
	    }
	    sumchi2/=sumn;
	    if (sumchi2<minchi2before1) minchi2before1=sumchi2; 
	  }
	  continue;   //safety space - Geo manager will give exact layer	   
	}
	track1b = btrack;
	minchi2after1 = btrack->GetNormChi2(i);
	break;
      }
      //
      // position resolution - used for DCA cut
      Float_t sigmad = track0b->GetSigmaY2()+track0b->GetSigmaZ2()+track1b->GetSigmaY2()+track1b->GetSigmaZ2()+
	(track0b->GetX()-pvertex->GetRr())*(track0b->GetX()-pvertex->GetRr())*(track0b->GetSigmaSnp2()+track0b->GetSigmaTgl2())+
	(track1b->GetX()-pvertex->GetRr())*(track1b->GetX()-pvertex->GetRr())*(track1b->GetSigmaSnp2()+track1b->GetSigmaTgl2());
      sigmad =TMath::Sqrt(sigmad)+0.04;
      if (pvertex->GetRr()>50){
	Double_t cov0[15],cov1[15];
	track0b->GetESDtrack()->GetInnerExternalCovariance(cov0);
	track1b->GetESDtrack()->GetInnerExternalCovariance(cov1);
	sigmad = cov0[0]+cov0[2]+cov1[0]+cov1[2]+
	  (80.-pvertex->GetRr())*(80.-pvertex->GetRr())*(cov0[5]+cov0[9])+
	  (80.-pvertex->GetRr())*(80.-pvertex->GetRr())*(cov1[5]+cov1[9]);
	sigmad =TMath::Sqrt(sigmad)+0.05;
      }
      //       
      AliV0 vertex2;
      vertex2.SetParamN(*track0b);
      vertex2.SetParamP(*track1b);
      vertex2.Update(primvertex);
      //Bo:      if (vertex2.GetDist2()<=pvertex->GetDist2()&&(vertex2.GetV0CosineOfPointingAngle()>=pvertex->GetV0CosineOfPointingAngle())){
      if (vertex2.GetDcaV0Daughters()<=pvertex->GetDcaV0Daughters()&&(vertex2.GetV0CosineOfPointingAngle()>=pvertex->GetV0CosineOfPointingAngle())){
	pvertex->SetParamN(*track0b);
	pvertex->SetParamP(*track1b);
	pvertex->Update(primvertex);
	pvertex->SetClusters(track0b->ClIndex(),track1b->ClIndex());  // register clusters
	pvertex->SetIndex(0,track0->GetESDtrack()->GetID());
	pvertex->SetIndex(1,track1->GetESDtrack()->GetID());
      }
      pvertex->SetDistSigma(sigmad);
      //Bo:      pvertex->SetDistNorm(pvertex->GetDist2()/sigmad);       
      pvertex->SetNormDCAPrim(normdist[itrack0],normdist[itrack1]);
      //
      // define likelihhod and causalities
      //
      Float_t pa0=1, pa1=1, pb0=0.26, pb1=0.26;      
      if (maxLayer<1){
	Float_t fnorm0 = normdist[itrack0];
	if (fnorm0<0) fnorm0*=-3;
	Float_t fnorm1 = normdist[itrack1];
	if (fnorm1<0) fnorm1*=-3;
 	if (pvertex->GetAnglep()[2]>0.1 ||  (pvertex->GetRr()<10.5)&& pvertex->GetAnglep()[2]>0.05 || pvertex->GetRr()<3){
 	  pb0    =  TMath::Exp(-TMath::Min(fnorm0,Float_t(16.))/12.);
 	  pb1    =  TMath::Exp(-TMath::Min(fnorm1,Float_t(16.))/12.);
 	}
	pvertex->SetChi2Before(normdist[itrack0]);
	pvertex->SetChi2After(normdist[itrack1]);       
	pvertex->SetNAfter(0);
	pvertex->SetNBefore(0);
      }else{
	pvertex->SetChi2Before(minchi2before0);
	pvertex->SetChi2After(minchi2before1);
	 if (pvertex->GetAnglep()[2]>0.1 || ( pvertex->GetRr()<10.5 && pvertex->GetAnglep()[2]>0.05) || pvertex->GetRr()<3){
	   pb0    =  TMath::Exp(-TMath::Min(minchi2before0,Float_t(16))/12.);
	   pb1    =  TMath::Exp(-TMath::Min(minchi2before1,Float_t(16))/12.);
	 }
	 pvertex->SetNAfter(maxLayer);
	 pvertex->SetNBefore(maxLayer);      
      }
      if (pvertex->GetRr()<90){
	pa0  *= TMath::Min(track0->GetESDtrack()->GetTPCdensity(0,60),Float_t(1.));
	pa1  *= TMath::Min(track1->GetESDtrack()->GetTPCdensity(0,60),Float_t(1.));
      }
      if (pvertex->GetRr()<20){
	pa0  *= (0.2+TMath::Exp(-TMath::Min(minchi2after0,Float_t(16))/8.))/1.2;
	pa1  *= (0.2+TMath::Exp(-TMath::Min(minchi2after1,Float_t(16))/8.))/1.2;
      }
      //
      pvertex->SetCausality(pb0,pb1,pa0,pa1);
      //
      //  Likelihood calculations  - apply cuts
      //         
      Bool_t v0OK = kTRUE;
      Float_t p12 = pvertex->GetParamP()->GetParameter()[4]*pvertex->GetParamP()->GetParameter()[4];
      p12        += pvertex->GetParamN()->GetParameter()[4]*pvertex->GetParamN()->GetParameter()[4];
      p12         = TMath::Sqrt(p12);                             // "mean" momenta
      Float_t    sigmap0   = 0.0001+0.001/(0.1+pvertex->GetRr()); 
      Float_t    sigmap    = 0.5*sigmap0*(0.6+0.4*p12);           // "resolution: of point angle - as a function of radius and momenta

      Float_t causalityA  = (1.0-pvertex->GetCausalityP()[0])*(1.0-pvertex->GetCausalityP()[1]);
      Float_t causalityB  = TMath::Sqrt(TMath::Min(pvertex->GetCausalityP()[2],Float_t(0.7))*
					TMath::Min(pvertex->GetCausalityP()[3],Float_t(0.7)));
      //
      //Bo:      Float_t likelihood0 = (TMath::Exp(-pvertex->GetDistNorm())+0.1) *(pvertex->GetDist2()<0.5)*(pvertex->GetDistNorm()<5);
      Float_t lDistNorm = pvertex->GetDcaV0Daughters()/pvertex->GetDistSigma();
      Float_t likelihood0 = (TMath::Exp(-lDistNorm)+0.1) *(pvertex->GetDcaV0Daughters()<0.5)*(lDistNorm<5);

      Float_t likelihood1 = TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/sigmap)+
	0.4*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/(4.*sigmap))+
	0.4*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/(8.*sigmap))+
	0.1*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/0.01);
      //
      if (causalityA<kCausality0Cut)                                          v0OK = kFALSE;
      if (TMath::Sqrt(likelihood0*likelihood1)<kLikelihood01Cut)              v0OK = kFALSE;
      if (likelihood1<kLikelihood1Cut)                                        v0OK = kFALSE;
      if (TMath::Power(likelihood0*likelihood1*causalityB,0.33)<kCombinedCut) v0OK = kFALSE;
      
      //
      //
      if (kFALSE){	
	Bool_t gold = TMath::Abs(TMath::Abs(track0->GetLabel())-TMath::Abs(track1->GetLabel()))==1;
	cstream<<"It0"<<
	  "Tr0.="<<track0<<                       //best without constrain
	  "Tr1.="<<track1<<                       //best without constrain  
	  "Tr0B.="<<track0b<<                     //best without constrain  after vertex
	  "Tr1B.="<<track1b<<                     //best without constrain  after vertex 
	  "Tr0C.="<<htrackc0<<                    //best with constrain     if exist
	  "Tr1C.="<<htrackc1<<                    //best with constrain     if exist
	  "Tr0L.="<<track0l<<                     //longest best           
	  "Tr1L.="<<track1l<<                     //longest best
	  "Esd0.="<<track0->GetESDtrack()<<           // esd track0 params
	  "Esd1.="<<track1->GetESDtrack()<<           // esd track1 params
	  "V0.="<<pvertex<<                       //vertex properties
	  "V0b.="<<&vertex2<<                       //vertex properties at "best" track
	  "ND0="<<normdist[itrack0]<<             //normalize distance for track0
	  "ND1="<<normdist[itrack1]<<             //normalize distance for track1
	  "Gold="<<gold<<                         //
	  //	  "RejectBase="<<rejectBase<<             //rejection in First itteration
	  "OK="<<v0OK<<
	  "rmin="<<rmin<<
	  "sigmad="<<sigmad<<
	  "\n";
      }      
      //if (rejectBase) continue;
      //
      pvertex->SetStatus(0);
      //      if (rejectBase) {
      //	pvertex->SetStatus(-100);
      //}
      if (pvertex->GetV0CosineOfPointingAngle()>kMinPointAngle2) {
	//Bo:	  pvertex->SetESDindexes(track0->GetESDtrack()->GetID(),track1->GetESDtrack()->GetID());
	pvertex->SetIndex(0,track0->GetESDtrack()->GetID());//Bo: consistency 0 for neg
	pvertex->SetIndex(1,track1->GetESDtrack()->GetID());//Bo: consistency 1 for pos
	if (v0OK){
	  //	  AliV0vertex vertexjuri(*track0,*track1);
	  //	  vertexjuri.SetESDindexes(track0->fESDtrack->GetID(),track1->fESDtrack->GetID());
	  //	  event->AddV0(&vertexjuri);
	  pvertex->SetStatus(100);
	}
        pvertex->SetOnFlyStatus(kTRUE);
        pvertex->ChangeMassHypothesis(kK0Short);
	event->AddV0(pvertex);
      }
    }
  }
  //
  //
  // delete temporary arrays
  //  
  delete[] minPointAngle;
  delete[] maxr;
  delete[] minr;
  delete[] norm;
  delete[] normdist;
  delete[] normdist1;
  delete[] normdist0;
  delete[] dist;
  delete[] itsmap;
  delete[] helixes;
  delete   pvertex;
}


void AliITStrackerMI::RefitV02(AliESD *event)
{
  //
  //try to refit  V0s in the third path of the reconstruction
  //
  TTreeSRedirector &cstream = *fDebugStreamer;
  //
  Int_t  nv0s = event->GetNumberOfV0s();
  Float_t primvertex[3]={GetX(),GetY(),GetZ()};
  AliV0 v0temp;
  for (Int_t iv0 = 0; iv0<nv0s;iv0++){
    AliV0 * v0mi = (AliV0*)event->GetV0(iv0);
    if (!v0mi) continue;
    Int_t     itrack0   = v0mi->GetIndex(0);
    Int_t     itrack1   = v0mi->GetIndex(1);
    AliESDtrack *esd0   = event->GetTrack(itrack0);
    AliESDtrack *esd1   = event->GetTrack(itrack1);
    if (!esd0||!esd1) continue;
    AliITStrackMI tpc0(*esd0);  
    AliITStrackMI tpc1(*esd1);
    Double_t x,y,z; v0mi->GetXYZ(x,y,z); //I.B. 
    Double_t alpha =TMath::ATan2(y,x);   //I.B.
    if (v0mi->GetRr()>85){
      if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
	v0temp.SetParamN(tpc0);
	v0temp.SetParamP(tpc1);
	v0temp.Update(primvertex);
	if (kFALSE) cstream<<"Refit"<<
	  "V0.="<<v0mi<<
	  "V0refit.="<<&v0temp<<
	  "Tr0.="<<&tpc0<<
	  "Tr1.="<<&tpc1<<
	  "\n";
	//Bo:	if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	  v0mi->SetParamN(tpc0);
	  v0mi->SetParamP(tpc1);
	  v0mi->Update(primvertex);
	}
      }
      continue;
    }
    if (v0mi->GetRr()>35){
      CorrectForDeadZoneMaterial(&tpc0);
      CorrectForDeadZoneMaterial(&tpc1);
      if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
	v0temp.SetParamN(tpc0);
	v0temp.SetParamP(tpc1);
	v0temp.Update(primvertex);
	if (kFALSE) cstream<<"Refit"<<
	  "V0.="<<v0mi<<
	  "V0refit.="<<&v0temp<<
	  "Tr0.="<<&tpc0<<
	  "Tr1.="<<&tpc1<<
	  "\n";
	//Bo:	if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	  v0mi->SetParamN(tpc0);
	  v0mi->SetParamP(tpc1);
	  v0mi->Update(primvertex);
	}	
      }
      continue;
    }
    CorrectForDeadZoneMaterial(&tpc0);
    CorrectForDeadZoneMaterial(&tpc1);    
    //    if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
    if (RefitAt(v0mi->GetRr(),&tpc0, v0mi->GetClusters(0)) && RefitAt(v0mi->GetRr(),&tpc1, v0mi->GetClusters(1))){
      v0temp.SetParamN(tpc0);
      v0temp.SetParamP(tpc1);
      v0temp.Update(primvertex);
      if (kFALSE) cstream<<"Refit"<<
	"V0.="<<v0mi<<
	"V0refit.="<<&v0temp<<
	"Tr0.="<<&tpc0<<
	"Tr1.="<<&tpc1<<
	"\n";
      //Bo:      if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
      if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	v0mi->SetParamN(tpc0);
	v0mi->SetParamP(tpc1);
	v0mi->Update(primvertex);
      }	
    }    
  }
}







