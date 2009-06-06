/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id: AliITStrackerHLT.cxx 32466 2009-05-20 07:51:56Z hristov $ */

//-------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    It reads AliITSRecPoint clusters and creates AliHLTITSTrack tracks
//                   and fills with them the ESD
//          Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch 
//          Current support and development: 
//                     Andrea Dainese, andrea.dainese@lnl.infn.it
//     dE/dx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//     Params moved to AliITSRecoParam by: Andrea Dainese, INFN
//     Material budget from TGeo by: Ludovic Gaudichet & Andrea Dainese, INFN
//-------------------------------------------------------------------------

//#include <TMatrixD.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TRandom.h>
#include <TTreeStream.h>


#include "AliLog.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliV0.h"
#include "AliHelix.h"
#include "AliITSChannelStatus.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSReconstructor.h"
#include "AliITSClusterParam.h"
#include "AliITSsegmentation.h"
#include "AliITSCalibration.h"
#include "AliITSV0Finder.h"
#include "AliITStrackerHLT.h"
//#include "AliHLTTPCCATrackParam.h"
#include "AliHLTVertexer.h"



Bool_t AliITStrackerHLT::TransportToX( AliExternalTrackParam *t, double x ) const
{
  return t->PropagateTo( x, t->GetBz() );
}

Bool_t AliITStrackerHLT::TransportToPhiX( AliExternalTrackParam *t, double phi, double x ) const
{
  return t->Propagate( phi, x, t->GetBz() );
}


//------------------------------------------------------------------------
Int_t AliITStrackerHLT::UpdateMI(AliHLTITSTrack* track, const AliITSRecPoint* cl,Double_t /*chi2*/,Int_t index) const 
{
  //
  // Update ITS track
  //
  
  if (cl->GetQ()<=0) return 0;  // ingore the "virtual" clusters
 
  Int_t layer = (index & 0xf0000000) >> 28;
 
  // Take into account the mis-alignment (bring track to cluster plane)

  Double_t xTrOrig=track->GetX();
  if (!TransportToX( track, xTrOrig + cl->GetX() ) ) return 0;
  
  Double_t err2Y, err2Z;

  GetClusterErrors2( layer, cl, track, err2Y, err2Z );

  Double_t p[2]={ cl->GetY(), cl->GetZ()};
  Double_t cov[3]={err2Y, 0., err2Z};

  Int_t updated = 1;
  //if( layer!=2 && layer!=3 ) 
  updated = track->AliExternalTrackParam::Update(p,cov);

  int n = track->GetNumberOfClusters();
  track->SetClusterIndex(n,index);
  track->SetNumberOfClusters(n+1);      
  
  return updated;
}


AliHLTITSLayer AliITStrackerHLT::fgLayers[AliITSgeomTGeo::kNLayers]; // ITS layers

AliITStrackerHLT::AliITStrackerHLT():AliTracker(),
fEsd(0),
fUseTGeo(3),
fxOverX0Pipe(-1.),
fxTimesRhoPipe(-1.),
fDebugStreamer(0),
fITSChannelStatus(0),
fkDetTypeRec(0)
{
  //Default constructor
  Int_t i;
  for(i=0;i<4;i++) fSPDdetzcentre[i]=0.;
  for(i=0;i<2;i++) {fxOverX0Shield[i]=-1.;fxTimesRhoShield[i]=-1.;}
  for(i=0;i<6;i++) {fxOverX0Layer[i]=-1.;fxTimesRhoLayer[i]=-1.;}
}
//------------------------------------------------------------------------
AliITStrackerHLT::AliITStrackerHLT(const Char_t *geom) : AliTracker(),
fEsd(0),
fUseTGeo(3),
fxOverX0Pipe(-1.),
fxTimesRhoPipe(-1.),
fDebugStreamer(0),
fITSChannelStatus(0),
fkDetTypeRec(0)
{
  //--------------------------------------------------------------------
  //This is the AliITStrackerHLT constructor
  //--------------------------------------------------------------------
  if (geom) {
    AliWarning("\"geom\" is actually a dummy argument !");
  }

  for (Int_t i=1; i<AliITSgeomTGeo::GetNLayers()+1; i++) {
    Int_t nlad=AliITSgeomTGeo::GetNLadders(i);
    Int_t ndet=AliITSgeomTGeo::GetNDetectors(i);

    Double_t xyz[3], &x=xyz[0], &y=xyz[1], &z=xyz[2];
    AliITSgeomTGeo::GetOrigTranslation(i,1,1,xyz); 
    Double_t poff=TMath::ATan2(y,x);
    Double_t zoff=z;
    Double_t r=TMath::Sqrt(x*x + y*y);

    AliITSgeomTGeo::GetOrigTranslation(i,1,2,xyz);
    r += TMath::Sqrt(x*x + y*y);
    AliITSgeomTGeo::GetOrigTranslation(i,2,1,xyz);
    r += TMath::Sqrt(x*x + y*y);
    AliITSgeomTGeo::GetOrigTranslation(i,2,2,xyz);
    r += TMath::Sqrt(x*x + y*y);
    r*=0.25;

    new (fgLayers+i-1) AliHLTITSLayer(r,poff,zoff,nlad,ndet);

    for (Int_t j=1; j<nlad+1; j++) {
      for (Int_t k=1; k<ndet+1; k++) { //Fill this layer with detectors
        TGeoHMatrix m; AliITSgeomTGeo::GetOrigMatrix(i,j,k,m);
        const TGeoHMatrix *tm=AliITSgeomTGeo::GetTracking2LocalMatrix(i,j,k);
        m.Multiply(tm);
        Double_t txyz[3]={0.};
	xyz[0]=0.;xyz[1]=0.;xyz[2]=0.;
        m.LocalToMaster(txyz,xyz);
        r=TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
        Double_t phi=TMath::ATan2(xyz[1],xyz[0]);

        if (phi<0) phi+=TMath::TwoPi();
        else if (phi>=TMath::TwoPi()) phi-=TMath::TwoPi();

        AliHLTITSDetector &det=fgLayers[i-1].GetDetector((j-1)*ndet + k-1); 
        new(&det) AliHLTITSDetector(r,phi); 
	// compute the real radius (with misalignment)
        TGeoHMatrix mmisal(*(AliITSgeomTGeo::GetMatrix(i,j,k)));
        mmisal.Multiply(tm);
	xyz[0]=0.;xyz[1]=0.;xyz[2]=0.;
        mmisal.LocalToMaster(txyz,xyz);
        Double_t rmisal=TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	det.SetRmisal(rmisal);
	
      } // end loop on detectors
    } // end loop on ladders
  } // end loop on layers

  

  Double_t xyzVtx[]={AliITSReconstructor::GetRecoParam()->GetXVdef(),
		     AliITSReconstructor::GetRecoParam()->GetYVdef(),
		     AliITSReconstructor::GetRecoParam()->GetZVdef()}; 
  Double_t ersVtx[]={AliITSReconstructor::GetRecoParam()->GetSigmaXVdef(),
		     AliITSReconstructor::GetRecoParam()->GetSigmaYVdef(),
		     AliITSReconstructor::GetRecoParam()->GetSigmaZVdef()}; 
  SetVertex(xyzVtx,ersVtx);

  // store positions of centre of SPD modules (in z)
  Double_t tr[3];
  AliITSgeomTGeo::GetTranslation(1,1,1,tr);
  fSPDdetzcentre[0] = tr[2];
  AliITSgeomTGeo::GetTranslation(1,1,2,tr);
  fSPDdetzcentre[1] = tr[2];
  AliITSgeomTGeo::GetTranslation(1,1,3,tr);
  fSPDdetzcentre[2] = tr[2];
  AliITSgeomTGeo::GetTranslation(1,1,4,tr);
  fSPDdetzcentre[3] = tr[2];

  fUseTGeo = AliITSReconstructor::GetRecoParam()->GetUseTGeoInTracker();
  if(AliITSReconstructor::GetRecoParam()->GetExtendedEtaAcceptance() && fUseTGeo!=1 && fUseTGeo!=3) {
    AliWarning("fUseTGeo changed to 3 because fExtendedEtaAcceptance is kTRUE");
    fUseTGeo = 3;
  }

  for(Int_t i=0;i<2;i++) {fxOverX0Shield[i]=-1.;fxTimesRhoShield[i]=-1.;}
  for(Int_t i=0;i<6;i++) {fxOverX0Layer[i]=-1.;fxTimesRhoLayer[i]=-1.;}
  
  fDebugStreamer = new TTreeSRedirector("ITSdebug.root");

}
//------------------------------------------------------------------------
AliITStrackerHLT::AliITStrackerHLT(const AliITStrackerHLT &tracker):AliTracker(tracker),
fEsd(tracker.fEsd),
fUseTGeo(tracker.fUseTGeo),
fxOverX0Pipe(tracker.fxOverX0Pipe),
fxTimesRhoPipe(tracker.fxTimesRhoPipe),
fDebugStreamer(tracker.fDebugStreamer),
fITSChannelStatus(tracker.fITSChannelStatus),
fkDetTypeRec(tracker.fkDetTypeRec)
{
  //Copy constructor
  Int_t i;
  for(i=0;i<4;i++) {
    fSPDdetzcentre[i]=tracker.fSPDdetzcentre[i];
  }
  for(i=0;i<6;i++) {
    fxOverX0Layer[i]=tracker.fxOverX0Layer[i];
    fxTimesRhoLayer[i]=tracker.fxTimesRhoLayer[i];
  }
  for(i=0;i<2;i++) {
    fxOverX0Shield[i]=tracker.fxOverX0Shield[i];
    fxTimesRhoShield[i]=tracker.fxTimesRhoShield[i];
  }
}
//------------------------------------------------------------------------
AliITStrackerHLT & AliITStrackerHLT::operator=(const AliITStrackerHLT &tracker){
  //Assignment operator
  this->~AliITStrackerHLT();
  new(this) AliITStrackerHLT(tracker);
  return *this;
}
//------------------------------------------------------------------------
AliITStrackerHLT::~AliITStrackerHLT()
{
  //
  //destructor
  //
  if (fDebugStreamer) {
    //fDebugStreamer->Close();
    delete fDebugStreamer;
  }
  if(fITSChannelStatus) delete fITSChannelStatus;
}


//------------------------------------------------------------------------
void AliITStrackerHLT::ReadBadFromDetTypeRec() {
  //--------------------------------------------------------------------
  //This function read ITS bad detectors, chips, channels from AliITSDetTypeRec
  //i.e. from OCDB
  //--------------------------------------------------------------------

  if(!AliITSReconstructor::GetRecoParam()->GetUseBadZonesFromOCDB()) return;

  Info("ReadBadFromDetTypeRec","Reading info about bad ITS detectors and channels");

  if(!fkDetTypeRec) Error("ReadBadFromDetTypeRec","AliITSDetTypeRec nof found!\n");

  // ITS channels map
  if(fITSChannelStatus) delete fITSChannelStatus;
  fITSChannelStatus = new AliITSChannelStatus(fkDetTypeRec);

  // ITS detectors and chips
  Int_t i=0,j=0,k=0,ndet=0;
  for (i=1; i<AliITSgeomTGeo::GetNLayers()+1; i++) {
    Int_t nBadDetsPerLayer=0;
    ndet=AliITSgeomTGeo::GetNDetectors(i);    
    for (j=1; j<AliITSgeomTGeo::GetNLadders(i)+1; j++) {
      for (k=1; k<ndet+1; k++) {
        AliHLTITSDetector &det=fgLayers[i-1].GetDetector((j-1)*ndet + k-1);  
	det.ReadBadDetectorAndChips(i-1,(j-1)*ndet + k-1,fkDetTypeRec);
	if(det.IsBad()) {nBadDetsPerLayer++;}
      } // end loop on detectors
    } // end loop on ladders
    Info("ReadBadFromDetTypeRec",Form("Layer %d: %d bad out of %d",i-1,nBadDetsPerLayer,ndet*AliITSgeomTGeo::GetNLadders(i)));
  } // end loop on layers
  
  return;
}
//------------------------------------------------------------------------
Int_t AliITStrackerHLT::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads ITS clusters
  //--------------------------------------------------------------------
  TBranch *branch=cTree->GetBranch("ITSRecPoints");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch !\n");
    return 1;
  }

  static TClonesArray dummy("AliITSRecPoint",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  Int_t i=0,j=0,ndet=0;
  Int_t detector=0;
  for (i=0; i<AliITSgeomTGeo::GetNLayers(); i++) {
    ndet=fgLayers[i].GetNdetectors();
    Int_t jmax = j + fgLayers[i].GetNladders()*ndet;
    for (; j<jmax; j++) {           
      if (!cTree->GetEvent(j)) continue;
      Int_t ncl=clusters->GetEntriesFast();
      SignDeltas(clusters,GetZ());
 
      while (ncl--) {
        AliITSRecPoint *c=(AliITSRecPoint*)clusters->UncheckedAt(ncl);
        detector=c->GetDetectorIndex();

	if (!c->Misalign()) AliWarning("Can't misalign this cluster !");

        fgLayers[i].InsertCluster(new AliITSRecPoint(*c));
      }
      clusters->Delete();
      // add dead zone "virtual" cluster in SPD, if there is a cluster within 
      // zwindow cm from the dead zone      
      if (i<2 && AliITSReconstructor::GetRecoParam()->GetAddVirtualClustersInDeadZone()) {
	for (Float_t xdead = 0; xdead < AliITSRecoParam::GetSPDdetxlength(); xdead += (i+1.)*AliITSReconstructor::GetRecoParam()->GetXPassDeadZoneHits()) {
	  Int_t lab[4]   = {0,0,0,detector};
	  Int_t info[3]  = {0,0,i};
	  Float_t q      = 0.; // this identifies virtual clusters
	  Float_t hit[5] = {xdead,
			    0.,
			    AliITSReconstructor::GetRecoParam()->GetSigmaXDeadZoneHit2(),
			    AliITSReconstructor::GetRecoParam()->GetSigmaZDeadZoneHit2(),
			    q};
	  Bool_t local   = kTRUE;
	  Double_t zwindow = AliITSReconstructor::GetRecoParam()->GetZWindowDeadZone();
	  hit[1] = fSPDdetzcentre[0]+0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	  hit[1] = fSPDdetzcentre[1]-0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	  hit[1] = fSPDdetzcentre[1]+0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	  hit[1] = fSPDdetzcentre[2]-0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	  hit[1] = fSPDdetzcentre[2]+0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	  hit[1] = fSPDdetzcentre[3]-0.5*AliITSRecoParam::GetSPDdetzlength();
	  if (TMath::Abs(fgLayers[i].GetDetector(detector).GetZmax()-hit[1])<zwindow) 
	    fgLayers[i].InsertCluster(new AliITSRecPoint(lab,hit,info,local));
	}
      } // "virtual" clusters in SPD
      
    }
    //
    fgLayers[i].ResetRoad(); //road defined by the cluster density
    fgLayers[i].SortClusters();
  }

  dummy.Clear();

  return 0;
}
//------------------------------------------------------------------------
void AliITStrackerHLT::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<AliITSgeomTGeo::GetNLayers(); i++) fgLayers[i].ResetClusters();
}



//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CorrectForTPCtoITSDeadZoneMaterial(AliHLTITSTrack *t) {
  //--------------------------------------------------------------------
  // Correction for the material between the TPC and the ITS
  //--------------------------------------------------------------------
  if (t->GetX() > AliITSRecoParam::Getriw()) {   // inward direction 
      if (!t->PropagateToTGeo(AliITSRecoParam::Getriw(),1)) return 0;// TPC inner wall
      if (!t->PropagateToTGeo(AliITSRecoParam::Getrcd(),1)) return 0;// TPC central drum
      if (!t->PropagateToTGeo(AliITSRecoParam::Getrs(),1))  return 0;// ITS screen
  } else if (t->GetX() < AliITSRecoParam::Getrs()) {  // outward direction
      if (!t->PropagateToTGeo(AliITSRecoParam::Getrs(),1))        return 0;// ITS screen
      if (!t->PropagateToTGeo(AliITSRecoParam::Getrcd(),1))       return 0;// TPC central drum
      if (!t->PropagateToTGeo(AliITSRecoParam::Getriw()+0.001,1)) return 0;// TPC inner wall
  } else {
    printf("CorrectForTPCtoITSDeadZoneMaterial: Track is already in the dead zone !\n");
    return 0;
  }
  
  return 1;
}

#include "TStopwatch.h"

//------------------------------------------------------------------------
Int_t AliITStrackerHLT::Clusters2Tracks(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------
  std::cout<<"\n\n ITS starts...\n"<<std::endl;
  TStopwatch timer;
  
  fEsd = event;         // store pointer to the esd 
 
  {/* Read ESD tracks */
   
    Double_t pimass = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    for( int itr=0; itr<event->GetNumberOfTracks(); itr++ ){

      AliESDtrack *esdTrack = event->GetTrack(itr);

      if ((esdTrack->GetStatus()&AliESDtrack::kTPCin)==0) continue;
      if (esdTrack->GetStatus()&AliESDtrack::kTPCout) continue;
      if (esdTrack->GetStatus()&AliESDtrack::kITSin) continue;
      if (esdTrack->GetKinkIndex(0)>0) continue;   //kink daughter
      AliHLTITSTrack tMI(*esdTrack);
      AliHLTITSTrack *t = &tMI;
      
      if (tMI.GetMass()<0.9*pimass) t->SetMass(pimass); 
      t->SetExpQ(TMath::Max(0.8*t->GetESDtrack()->GetTPCsignal(),30.));

      if (!CorrectForTPCtoITSDeadZoneMaterial(t))  continue;

      Int_t tpcLabel=t->GetLabel(); //save the TPC track label       
      
      FollowProlongationTree(t);
      int nclu=0;
      for(Int_t i=0; i<6; i++) {
	if( t->GetClusterIndex(i)>=0 ) nclu++; 
      }
      cout<<"N assigned ITS clusters = "<<nclu<<std::endl;
      if( nclu>0 ){
	t->SetLabel(-1);//tpcLabel);
	t->SetFakeRatio(1.);
	CookLabel(t,0.); //For comparison only
	cout<<"label = "<<t->GetLabel()<<" / "<<tpcLabel<<endl;
	TransportToX(t, 0 );
	cout<<"\n fill track : parameters at "<<t->GetX()<<": "<< TMath::Sqrt(TMath::Abs(t->GetSigmaY2()))<<" "<< TMath::Sqrt(TMath::Abs(t->GetSigmaY2()))<<endl;
	//t->Print();
	UpdateESDtrack(t,AliESDtrack::kITSin);          
      }
    }
  } /* End Read ESD tracks */

  AliHLTVertexer vertexer;
  vertexer.SetESD( event );
  vertexer.FindPrimaryVertex();
  vertexer.FindV0s();


  timer.Stop();
  static double totalTime = 0;
  static int nEvnts = 0;
  totalTime+=timer.CpuTime();
  nEvnts++;
  std::cout<<"\n\n ITS tracker time = "<<totalTime/nEvnts<<" [s/ev]  for "<<nEvnts<<" events\n\n "<<std::endl;
  return 0;
}


//------------------------------------------------------------------------
Int_t AliITStrackerHLT::PropagateBack(AliESDEvent * /*event*/) {
  return 0;
}

//------------------------------------------------------------------------
Int_t AliITStrackerHLT::RefitInward(AliESDEvent * /*event*/ ) {
  return 0;
}
//------------------------------------------------------------------------
AliCluster *AliITStrackerHLT::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fgLayers[l].GetCluster(c);
}
//------------------------------------------------------------------------
Bool_t AliITStrackerHLT::GetTrackPoint(Int_t index, AliTrackPoint& p) const {
  //--------------------------------------------------------------------
  // Get track space point with index i
  //--------------------------------------------------------------------

  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  AliITSRecPoint *cl = fgLayers[l].GetCluster(c);
  Int_t idet = cl->GetDetectorIndex();

  Float_t xyz[3];
  Float_t cov[6];
  cl->GetGlobalXYZ(xyz);
  cl->GetGlobalCov(cov);
  p.SetXYZ(xyz, cov);
  p.SetCharge(cl->GetQ());
  p.SetDriftTime(cl->GetDriftTime());
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer; 
  switch (l) {
  case 0:
    iLayer = AliGeomManager::kSPD1;
    break;
  case 1:
    iLayer = AliGeomManager::kSPD2;
    break;
  case 2:
    iLayer = AliGeomManager::kSDD1;
    break;
  case 3:
    iLayer = AliGeomManager::kSDD2;
    break;
  case 4:
    iLayer = AliGeomManager::kSSD1;
    break;
  case 5:
    iLayer = AliGeomManager::kSSD2;
    break;
  default:
    AliWarning(Form("Wrong layer index in ITS (%d) !",l));
    break;
  };
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,idet);
  p.SetVolumeID((UShort_t)volid);
  return kTRUE;
}
//------------------------------------------------------------------------
Bool_t AliITStrackerHLT::GetTrackPointTrackingError(Int_t index, 
			AliTrackPoint& p, const AliESDtrack *t) {
  //--------------------------------------------------------------------
  // Get track space point with index i
  // (assign error estimated during the tracking)
  //--------------------------------------------------------------------

  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  const AliITSRecPoint *cl = fgLayers[l].GetCluster(c);
  Int_t idet = cl->GetDetectorIndex();

  const AliHLTITSDetector &det=fgLayers[l].GetDetector(idet);

  // tgphi and tglambda of the track in tracking frame with alpha=det.GetPhi
  Float_t detxy[2];
  detxy[0] = det.GetR()*TMath::Cos(det.GetPhi());
  detxy[1] = det.GetR()*TMath::Sin(det.GetPhi());
  Double_t alpha = t->GetAlpha();
  Double_t xdetintrackframe = detxy[0]*TMath::Cos(alpha)+detxy[1]*TMath::Sin(alpha);
  Float_t phi = TMath::ASin(t->GetSnpAt(xdetintrackframe,GetBz()));
  phi += alpha-det.GetPhi();
  Float_t tgphi = TMath::Tan(phi);

  Float_t tgl = t->GetTgl(); // tgl about const along track
  Float_t expQ = TMath::Max(0.8*t->GetTPCsignal(),30.);

  Float_t errlocalx,errlocalz;
  Bool_t addMisalErr=kFALSE;
  AliITSClusterParam::GetError(l,cl,tgl,tgphi,expQ,errlocalx,errlocalz,addMisalErr);

  Float_t xyz[3];
  Float_t cov[6];
  cl->GetGlobalXYZ(xyz);
  //  cl->GetGlobalCov(cov);
  Float_t pos[3] = {0.,0.,0.};
  AliCluster tmpcl((UShort_t)cl->GetVolumeId(),pos[0],pos[1],pos[2],errlocalx*errlocalx,errlocalz*errlocalz,0);
  tmpcl.GetGlobalCov(cov);

  p.SetXYZ(xyz, cov);
  p.SetCharge(cl->GetQ());
  p.SetDriftTime(cl->GetDriftTime());

  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer; 
  switch (l) {
  case 0:
    iLayer = AliGeomManager::kSPD1;
    break;
  case 1:
    iLayer = AliGeomManager::kSPD2;
    break;
  case 2:
    iLayer = AliGeomManager::kSDD1;
    break;
  case 3:
    iLayer = AliGeomManager::kSDD2;
    break;
  case 4:
    iLayer = AliGeomManager::kSSD1;
    break;
  case 5:
    iLayer = AliGeomManager::kSSD2;
    break;
  default:
    AliWarning(Form("Wrong layer index in ITS (%d) !",l));
    break;
  };
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,idet);

  p.SetVolumeID((UShort_t)volid);
  return kTRUE;
}


//------------------------------------------------------------------------
void AliITStrackerHLT::FollowProlongationTree(AliHLTITSTrack * track ) 
{
    //cout<<endl;
  for (Int_t ilayer=5; ilayer>=0; ilayer--) {
    //cout<<"\nLayer "<<ilayer<<endl;
    
    AliHLTITSLayer &layer=fgLayers[ilayer];
  
    //cout<<" shield material.. "<<endl;

    // material between SSD and SDD, SDD and SPD
    if (ilayer==3 && !CorrectForShieldMaterial(track,"SDD","inward")) continue;
    if (ilayer==1 && !CorrectForShieldMaterial(track,"SPD","inward")) continue;
    
    int idet;

    {
      // propagate to the layer radius
      
      Double_t r = layer.GetR(), phi,z;
      Double_t xToGo;
      //cout<<" propagate to layer R= "<<r<<" .."<<endl;
      if( !track->GetLocalXat(r,xToGo) ) continue;
      if( !TransportToX(track, xToGo) ) continue;

      // detector number
      
      if (!track->GetPhiZat(r,phi,z)) continue;
      idet=layer.FindDetectorIndex(phi,z);
      //cout<<" detector number = "<<idet<<endl;
   }


    //cout<<" correct for the layer material .. "<<endl;

    // correct for the layer material
    {
      Double_t trackGlobXYZ1[3];
      if (!track->GetXYZ(trackGlobXYZ1)) continue;
      CorrectForLayerMaterial(track,ilayer,trackGlobXYZ1,"inward");
    }

    // track outside layer acceptance in z
    
    if( idet<0 ) continue;
    
    // propagate to the intersection with the detector plane
     
    //cout<<" propagate to the intersection with the detector .. "<<endl;

    const AliHLTITSDetector &det=layer.GetDetector( idet );
    if (!TransportToPhiX( track, det.GetPhi(), det.GetR() ) ) continue;

    // DEFINITION OF SEARCH ROAD AND CLUSTERS SELECTION
    
    // road in global (rphi,z) [i.e. in tracking ref. system]
    
    Double_t zmin,zmax,ymin,ymax;
    
    //cout<<" ComputeRoad .. "<<endl;
    if (!ComputeRoad(track,ilayer,idet,zmin,zmax,ymin,ymax)) continue;
  
    //cout<<" Road: y["<<ymin<<","<<ymax<<"], z["<<zmin<<","<<zmax<<"] "<<endl;

    // select clusters in road
    
    //cout<<" SelectClusters .. "<<endl;
    layer.SelectClusters(zmin,zmax,ymin,ymax);     
    
    // Define criteria for track-cluster association
    
    Double_t msz = track->GetSigmaZ2() + 
      AliITSReconstructor::GetRecoParam()->GetNSigmaZLayerForRoadZ()*
      AliITSReconstructor::GetRecoParam()->GetNSigmaZLayerForRoadZ()*
      AliITSReconstructor::GetRecoParam()->GetSigmaZ2(ilayer);

    Double_t msy = track->GetSigmaY2() + 
      AliITSReconstructor::GetRecoParam()->GetNSigmaYLayerForRoadY()*
      AliITSReconstructor::GetRecoParam()->GetNSigmaYLayerForRoadY()*
      AliITSReconstructor::GetRecoParam()->GetSigmaY2(ilayer);

     msz *= AliITSReconstructor::GetRecoParam()->GetNSigma2RoadZNonC();
     msy *= AliITSReconstructor::GetRecoParam()->GetNSigma2RoadYNonC(); 
  
     msz = 1./msz; // 1/RoadZ^2
     msy = 1./msy; // 1/RoadY^2
     
     //
     // LOOP OVER ALL POSSIBLE TRACK PROLONGATIONS ON THIS LAYER
     //

     const AliITSRecPoint *cl=0; 
     Int_t clidx=-1;
     Double_t chi2trkcl=AliITSReconstructor::GetRecoParam()->GetMaxChi2(); // init with big value
     Bool_t deadzoneSPD=kFALSE;

     // check if the road contains a dead zone 
     Bool_t noClusters = !layer.GetNextCluster(clidx,kTRUE);
     
     Double_t dz=0.5*(zmax-zmin);
     Double_t dy=0.5*(ymax-ymin);
     
     Int_t dead = CheckDeadZone(track,ilayer,idet,dz,dy,noClusters); 

     // create a prolongation without clusters (check also if there are no clusters in the road)

     if (dead==1) { // dead zone at z=0,+-7cm in SPD
       deadzoneSPD=kTRUE;
     }
     
     clidx=-1;

     //cout<<" loop over clusters in the road .. "<<endl;

     // loop over clusters in the road     
     const AliITSRecPoint *bestCluster=0; 
     double bestChi2 = 1.e10;
     AliHLTITSTrack bestTrack( *track );
     int bestIdx = -1;
     while ((cl=layer.GetNextCluster(clidx))!=0) {        
       //cout<<" cluster: "<<cl->GetX()<<" "<<cl->GetY()<<" "<<cl->GetZ()<<endl;
       AliHLTITSTrack t(*track);
       if (cl->GetQ()==0 && deadzoneSPD==kTRUE) continue;
       
       Int_t idetc=cl->GetDetectorIndex();
       
       //cout<<" cluster detector: "<<idetc<<endl;

       if ( idet !=idetc ) { // new cluster's detector
	 const AliHLTITSDetector &detc=layer.GetDetector(idetc);
	 if (!TransportToPhiX( track, detc.GetPhi(),detc.GetR()) ) continue;
	 t = *track;
	 idet = idetc;
       }

       // take into account misalignment (bring track to real detector plane)

       if (!TransportToX( &t, t.GetX() + cl->GetX() ) ) continue;
       double chi2 = ( (t.GetZ()-cl->GetZ())*(t.GetZ()-cl->GetZ())*msz + 
		       (t.GetY()-cl->GetY())*(t.GetY()-cl->GetY())*msy   );
       //cout<<" chi2="<<chi2<<endl;
       if ( chi2 < bestChi2 ){
	 bestChi2 = chi2;
	 bestCluster = cl;
	 bestTrack = t;
	 bestIdx = clidx;
	 continue;
       }
     }

     //cout<<" best chi2= "<<bestChi2<<endl;

     if( bestCluster && bestChi2 <=1. ){

       //cout<<" cluster found "<<endl;
       *track = bestTrack;

 
      // calculate track-clusters chi2
       chi2trkcl = GetPredictedChi2MI(track,bestCluster,ilayer); 
       //cout<<" track-clusters chi2 = "<<chi2trkcl<<endl;
       //cout<<" max chi2 = "<<AliITSReconstructor::GetRecoParam()->GetMaxChi2s(ilayer)<<endl;

       // chi2 cut
       if (chi2trkcl < AliITSReconstructor::GetRecoParam()->GetMaxChi2s(ilayer)) {
	 if (bestCluster->GetQ()==0) deadzoneSPD=kTRUE; // only 1 prolongation with virtual cluster
	 //cout<<"set index.."<<endl;	 
	 //cout<<"set index ok"<<endl;
	 if (bestCluster->GetQ()!=0) { // real cluster	   
	   //cout<<" UpdateMI ... "<<endl;
	   if (!UpdateMI(track,bestCluster,chi2trkcl,(ilayer<<28)+bestIdx)) {
	     continue;
	   } 
	 }
       }          
     }
     //cout<<" goto next layer "<<endl;

  }
}


//------------------------------------------------------------------------
AliHLTITSLayer & AliITStrackerHLT::GetLayer(Int_t layer) const
{
  //--------------------------------------------------------------------
  //
  //
  return fgLayers[layer];
}













//------------------------------------------------------------------------
void AliITStrackerHLT::CookLabel(AliHLTITSTrack *track,Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t tpcLabel=-1; 
     
  if ( track->GetESDtrack())   tpcLabel =  TMath::Abs(track->GetESDtrack()->GetTPCLabel());

   Int_t nwrong=0;
   cout<<"cook label: nclu = "<<track->GetNumberOfClusters()<<endl;
   for (Int_t i=0;i<track->GetNumberOfClusters();i++){
     Int_t cindex = track->GetClusterIndex(i);
     //Int_t l=(cindex & 0xf0000000) >> 28;
     AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(cindex);
     Int_t isWrong=1;
     for (Int_t ind=0;ind<3;ind++){
       if (tpcLabel>0)
	 if (cl->GetLabel(ind)==tpcLabel) isWrong=0;
       AliDebug(2,Form("icl %d  ilab %d lab %d",i,ind,cl->GetLabel(ind)));
     }
      nwrong+=isWrong;
   }
   Int_t nclusters = track->GetNumberOfClusters();
   if (nclusters > 0) //PH Some tracks don't have any cluster
     track->SetFakeRatio(double(nwrong)/double(nclusters));
   cout<<"fake ratio = "<<track->GetFakeRatio()<<endl;
   if (tpcLabel>0){
     if (track->GetFakeRatio()>wrong) track->SetLabel(-tpcLabel);
     else
       track->SetLabel(tpcLabel);
   }
   AliDebug(2,Form(" nls %d wrong %d  label %d  tpcLabel %d\n",nclusters,nwrong,track->GetLabel(),tpcLabel));
   
}



void AliITStrackerHLT::GetClusterErrors2( Int_t layer, const AliITSRecPoint *cluster, AliHLTITSTrack* track, double &err2Y, double &err2Z ) const
{
  Float_t erry,errz;
  Float_t theta = track->GetTgl();
  Float_t phi   = track->GetSnp();
  phi = TMath::Abs(phi)*TMath::Sqrt(1./((1.-phi)*(1.+phi)));
  AliITSClusterParam::GetError(layer,cluster,theta,phi,track->GetExpQ(),erry,errz);
  err2Y = erry*erry;
  err2Z = errz*errz;
}


//------------------------------------------------------------------------
Double_t AliITStrackerHLT::GetPredictedChi2MI(AliHLTITSTrack* track, const AliITSRecPoint *cluster,Int_t layer) 
{
  //
  // Compute predicted chi2
  //

  AliHLTITSTrack t(*track);
  if (!t.Propagate( t.GetX()+cluster->GetX())) return 1000.;
 
  Double_t  err2Y,err2Z;  
  GetClusterErrors2( layer, cluster, &t, err2Y, err2Z );

  Double_t chi2 = t.GetPredictedChi2(cluster->GetY(),cluster->GetZ(),err2Y,err2Z);

  Float_t ny,nz; 
  Float_t theta = track->GetTgl();
  Float_t phi   = track->GetSnp();
  phi = TMath::Abs(phi)*TMath::Sqrt(1./((1.-phi)*(1.+phi)));
  AliITSClusterParam::GetNTeor(layer,cluster,theta,phi,ny,nz);  
  Double_t delta = cluster->GetNy()+cluster->GetNz()-nz-ny;
  if (delta>1){
    chi2+=0.5*TMath::Min(delta/2,2.);
    chi2+=2.*cluster->GetDeltaProbability();
  }
  return chi2;
}


//------------------------------------------------------------------------
void AliITStrackerHLT::SignDeltas(const TObjArray *clusterArray, Float_t vz)
{
  //
  // Clusters from delta electrons?
  //  
  Int_t entries = clusterArray->GetEntriesFast();
  if (entries<4) return;
  AliITSRecPoint* cluster = (AliITSRecPoint*)clusterArray->At(0);
  Int_t layer = cluster->GetLayer();
  if (layer>1) return;
  Int_t index[10000];
  Int_t ncandidates=0;
  Float_t r = (layer>0)? 7:4;
  // 
  for (Int_t i=0;i<entries;i++){
    AliITSRecPoint* cl0 = (AliITSRecPoint*)clusterArray->At(i);
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
    AliITSRecPoint* cl0 = (AliITSRecPoint*)clusterArray->At(index[i]);
    if (cl0->GetDeltaProbability()>0.8) continue;
    // 
    Int_t ncl = 0;
    Float_t y[100],z[100],sumy,sumz,sumy2, sumyz, sumw;
    sumy=sumz=sumy2=sumyz=sumw=0.0;
    for (Int_t j=0;j<ncandidates;j++){
      if (i==j) continue;
      AliITSRecPoint* cl1 = (AliITSRecPoint*)clusterArray->At(index[j]);
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
//------------------------------------------------------------------------
void AliITStrackerHLT::UpdateESDtrack(AliHLTITSTrack* track, ULong_t flags) const
{
  //
  // Update ESD track
  //
  track->UpdateESDtrack(flags);
  AliHLTITSTrack * oldtrack = (AliHLTITSTrack*)(track->GetESDtrack()->GetITStrack());
  if (oldtrack) delete oldtrack; 
  track->GetESDtrack()->SetITStrack(new AliHLTITSTrack(*track));
}




//------------------------------------------------------------------------
void AliITStrackerHLT::BuildMaterialLUT(TString material) {
  //--------------------------------------------------------------------
  // Fill a look-up table with mean material
  //--------------------------------------------------------------------

  Int_t n=1000;
  Double_t mparam[7];
  Double_t point1[3],point2[3];
  Double_t phi,cosphi,sinphi,z;
  // 0-5 layers, 6 pipe, 7-8 shields 
  Double_t rmin[9]={ 3.5, 5.5,13.0,22.0,35.0,41.0, 2.0, 8.0,25.0};
  Double_t rmax[9]={ 5.5, 8.0,17.0,26.0,41.0,47.0, 3.0,10.5,30.0};

  Int_t ifirst=0,ilast=0;  
  if(material.Contains("Pipe")) {
    ifirst=6; ilast=6;
  } else if(material.Contains("Shields")) {
    ifirst=7; ilast=8;
  } else if(material.Contains("Layers")) {
    ifirst=0; ilast=5;
  } else {
    Error("BuildMaterialLUT","Wrong layer name\n");
  }

  for(Int_t imat=ifirst; imat<=ilast; imat++) {
    Double_t param[5]={0.,0.,0.,0.,0.};
    for (Int_t i=0; i<n; i++) {
      phi = 2.*TMath::Pi()*gRandom->Rndm();
      cosphi = TMath::Cos(phi); sinphi = TMath::Sin(phi); 
      z = 14.*(-1.+2.*gRandom->Rndm()); // SPD barrel
      point1[0] = rmin[imat]*cosphi;
      point1[1] = rmin[imat]*sinphi;
      point1[2] = z;
      point2[0] = rmax[imat]*cosphi;
      point2[1] = rmax[imat]*sinphi;
      point2[2] = z;
      AliTracker::MeanMaterialBudget(point1,point2,mparam);
      for(Int_t j=0;j<5;j++) param[j]+=mparam[j];
    }
    for(Int_t j=0;j<5;j++) param[j]/=(Float_t)n;
    if(imat<=5) {
      fxOverX0Layer[imat] = param[1];
      fxTimesRhoLayer[imat] = param[0]*param[4];
    } else if(imat==6) {
      fxOverX0Pipe = param[1];
      fxTimesRhoPipe = param[0]*param[4];
    } else if(imat==7) {
      fxOverX0Shield[0] = param[1];
      fxTimesRhoShield[0] = param[0]*param[4];
    } else if(imat==8) {
      fxOverX0Shield[1] = param[1];
      fxTimesRhoShield[1] = param[0]*param[4];
    }
  }
  /*
  printf("%s\n",material.Data());
  printf("%f  %f\n",fxOverX0Pipe,fxTimesRhoPipe);
  printf("%f  %f\n",fxOverX0Shield[0],fxTimesRhoShield[0]);
  printf("%f  %f\n",fxOverX0Shield[1],fxTimesRhoShield[1]);
  printf("%f  %f\n",fxOverX0Layer[0],fxTimesRhoLayer[0]);
  printf("%f  %f\n",fxOverX0Layer[1],fxTimesRhoLayer[1]);
  printf("%f  %f\n",fxOverX0Layer[2],fxTimesRhoLayer[2]);
  printf("%f  %f\n",fxOverX0Layer[3],fxTimesRhoLayer[3]);
  printf("%f  %f\n",fxOverX0Layer[4],fxTimesRhoLayer[4]);
  printf("%f  %f\n",fxOverX0Layer[5],fxTimesRhoLayer[5]);
  */
  return;
}

//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CorrectForPipeMaterial(AliHLTITSTrack *t,
					      TString direction) {
  //-------------------------------------------------------------------
  // Propagate beyond beam pipe and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  // Define budget mode:
  // 0: material from AliITSRecoParam (hard coded)
  // 1: material from TGeo in one step (on the fly)
  // 2: material from lut
  // 3: material from TGeo in one step (same for all hypotheses)
  Int_t mode;
  switch(fUseTGeo) {
  case 0:
    mode=0; 
    break;    
  case 1:
    mode=1;
    break;    
  case 2:
    mode=2;
    break;
  case 3:
    mode=3; 
    break;
  case 4:
    mode=3;
    break;
  default:
    mode=0;
    break;
  }
  
  Float_t  dir = (direction.Contains("inward") ? 1. : -1.);
  Double_t rToGo=(dir>0 ? AliITSRecoParam::GetrInsidePipe() : AliITSRecoParam::GetrOutsidePipe());
  Double_t xToGo;
  if (!t->GetLocalXat(rToGo,xToGo)) return 0;

  Double_t xOverX0,x0,lengthTimesMeanDensity;

  switch(mode) {
  case 0:
    xOverX0 = AliITSRecoParam::GetdPipe();
    x0 = AliITSRecoParam::GetX0Be();
    lengthTimesMeanDensity = xOverX0*x0;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 1:
    if (!t->PropagateToTGeo(xToGo,1)) return 0;
    break;
  case 2:
    if(fxOverX0Pipe<0) BuildMaterialLUT("Pipe");  
    xOverX0 = fxOverX0Pipe;
    lengthTimesMeanDensity = fxTimesRhoPipe;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 3:
    double xOverX0PipeTrks, xTimesRhoPipeTrks;
    if (!t->PropagateToTGeo(xToGo,1,xOverX0,lengthTimesMeanDensity)) return 0;
    Double_t angle=TMath::Sqrt((1.+t->GetTgl()*t->GetTgl())/
			       ((1.-t->GetSnp())*(1.+t->GetSnp())));
    xOverX0PipeTrks = TMath::Abs(xOverX0)/angle;
    xTimesRhoPipeTrks = TMath::Abs(lengthTimesMeanDensity)/angle;
    xOverX0 = xOverX0PipeTrks;
    lengthTimesMeanDensity = xTimesRhoPipeTrks;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  }
  
  return 1;
}
//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CorrectForShieldMaterial(AliHLTITSTrack *t,
						TString shield,
						TString direction) {
  //-------------------------------------------------------------------
  // Propagate beyond SPD or SDD shield and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  // Define budget mode:
  // 0: material from AliITSRecoParam (hard coded)
  // 1: material from TGeo in steps of X cm (on the fly)
  //    X = AliITSRecoParam::GetStepSizeTGeo()
  // 2: material from lut
  // 3: material from TGeo in one step (same for all hypotheses)
  Int_t mode;
  switch(fUseTGeo) {
  case 0:
    mode=0; 
    break;    
  case 1:
    mode=1;
    break;    
  case 2:
    mode=2;
    break;
  case 3:
    mode=3;
    break;
  case 4:
    mode=3;
    break;
  default:
    mode=0;
    break;
  }


  Float_t  dir = (direction.Contains("inward") ? 1. : -1.);
  Double_t rToGo;
  Int_t    shieldindex=0;
  if (shield.Contains("SDD")) { // SDDouter
    rToGo=(dir>0 ? AliITSRecoParam::GetrInsideShield(1) : AliITSRecoParam::GetrOutsideShield(1));
    shieldindex=1;
  } else if (shield.Contains("SPD")) {        // SPDouter
    rToGo=(dir>0 ? AliITSRecoParam::GetrInsideShield(0) : AliITSRecoParam::GetrOutsideShield(0)); 
    shieldindex=0;
  } else {
    Error("CorrectForShieldMaterial"," Wrong shield name\n");
    return 0;
  }
  Double_t xToGo;
  if (!t->GetLocalXat(rToGo,xToGo)) return 0;

  Double_t xOverX0,x0,lengthTimesMeanDensity;
  Int_t nsteps=1;

  switch(mode) {
  case 0:
    xOverX0 = AliITSRecoParam::Getdshield(shieldindex);
    x0 = AliITSRecoParam::GetX0shield(shieldindex);
    lengthTimesMeanDensity = xOverX0*x0;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 1:
    nsteps= (Int_t)(TMath::Abs(t->GetX()-xToGo)/AliITSReconstructor::GetRecoParam()->GetStepSizeTGeo())+1;
    if (!t->PropagateToTGeo(xToGo,nsteps)) return 0; // cross the material and apply correction
    break;
  case 2:
    if(fxOverX0Shield[shieldindex]<0) BuildMaterialLUT("Shields");  
    xOverX0 = fxOverX0Shield[shieldindex];
    lengthTimesMeanDensity = fxTimesRhoShield[shieldindex];
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 3:    
    if (!t->PropagateToTGeo(xToGo,1,xOverX0,lengthTimesMeanDensity)) return 0;
    Double_t angle=TMath::Sqrt((1.+t->GetTgl()*t->GetTgl())/
			       ((1.-t->GetSnp())*(1.+t->GetSnp())));
    double xOverX0ShieldTrks = TMath::Abs(xOverX0)/angle;
    double xTimesRhoShieldTrks = TMath::Abs(lengthTimesMeanDensity)/angle;
    xOverX0 = xOverX0ShieldTrks;
    lengthTimesMeanDensity = xTimesRhoShieldTrks;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  }

  return 1;
}
//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CorrectForLayerMaterial(AliHLTITSTrack *t,
					       Int_t layerindex,
					       Double_t oldGlobXYZ[3],
					       TString direction) {
  //-------------------------------------------------------------------
  // Propagate beyond layer and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  // Define budget mode:
  // 0: material from AliITSRecoParam (hard coded)
  // 1: material from TGeo in stepsof X cm (on the fly)
  //    X = AliITSRecoParam::GetStepSizeTGeo()
  // 2: material from lut
  // 3: material from TGeo in one step (same for all hypotheses)
  Int_t mode;
  switch(fUseTGeo) {
  case 0:
    mode=0; 
    break;    
  case 1:
    mode=1;
    break;    
  case 2:
    mode=2;
    break;
  case 3:
    mode=3;
    break;
  case 4:
    mode=3;
    break;
  default:
    mode=0;
    break;
  }


  Float_t  dir = (direction.Contains("inward") ? 1. : -1.);

  Double_t r=fgLayers[layerindex].GetR();
  Double_t deltar=(layerindex<2 ? 0.10*r : 0.05*r);

  Double_t rToGo=TMath::Sqrt(t->GetX()*t->GetX()+t->GetY()*t->GetY())-deltar*dir;
  Double_t xToGo;
  if (!t->GetLocalXat(rToGo,xToGo)) return 0;  


  Double_t xOverX0=0.0,x0=0.0,lengthTimesMeanDensity=0.0;
  Int_t nsteps=1;

  // back before material (no correction)
  Double_t rOld,xOld;
  rOld=TMath::Sqrt(oldGlobXYZ[0]*oldGlobXYZ[0]+oldGlobXYZ[1]*oldGlobXYZ[1]);
  if (!t->GetLocalXat(rOld,xOld)) return 0;
  if (!TransportToX( t, xOld)) return 0;

  switch(mode) {
  case 0:
    xOverX0 = fgLayers[layerindex].GetThickness(t->GetY(),t->GetZ(),x0);
    lengthTimesMeanDensity = xOverX0*x0;
    lengthTimesMeanDensity *= dir;
    // Bring the track beyond the material
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 1:
    nsteps = (Int_t)(TMath::Abs(xOld-xToGo)/AliITSReconstructor::GetRecoParam()->GetStepSizeTGeo())+1;
    if (!t->PropagateToTGeo(xToGo,nsteps)) return 0; // cross the material and apply correction
    break;
  case 2:
    if(fxOverX0Layer[layerindex]<0) BuildMaterialLUT("Layers");  
    xOverX0 = fxOverX0Layer[layerindex];
    lengthTimesMeanDensity = fxTimesRhoLayer[layerindex];
    lengthTimesMeanDensity *= dir;
    // Bring the track beyond the material
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 3:
    nsteps = (Int_t)(TMath::Abs(xOld-xToGo)/AliITSReconstructor::GetRecoParam()->GetStepSizeTGeo())+1;
    if (!t->PropagateToTGeo(xToGo,nsteps,xOverX0,lengthTimesMeanDensity)) return 0;
    Double_t angle=TMath::Sqrt((1.+t->GetTgl()*t->GetTgl())/
			       ((1.-t->GetSnp())*(1.+t->GetSnp())));
    double xOverX0LayerTrks = TMath::Abs(xOverX0)/angle;
    double xTimesRhoLayerTrks = TMath::Abs(lengthTimesMeanDensity)/angle;
    xOverX0 = xOverX0LayerTrks;
    lengthTimesMeanDensity = xTimesRhoLayerTrks;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  }


  return 1;
}


//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CheckDeadZone(AliHLTITSTrack *track,
				     Int_t ilayer,Int_t idet,
				     Double_t dz,Double_t dy,
				     Bool_t noClusters) const {
  //-----------------------------------------------------------------
  // This method is used to decide whether to allow a prolongation 
  // without clusters, because there is a dead zone in the road.
  // In this case the return value is > 0:
  // return 1: dead zone at z=0,+-7cm in SPD
  // return 2: all road is "bad" (dead or noisy) from the OCDB
  // return 3: something "bad" (dead or noisy) from the OCDB
  //-----------------------------------------------------------------

  // check dead zones at z=0,+-7cm in the SPD
  if (ilayer<2 && !AliITSReconstructor::GetRecoParam()->GetAddVirtualClustersInDeadZone()) {
    Double_t zmindead[3]={fSPDdetzcentre[0] + 0.5*AliITSRecoParam::GetSPDdetzlength(),
			  fSPDdetzcentre[1] + 0.5*AliITSRecoParam::GetSPDdetzlength(),
			  fSPDdetzcentre[2] + 0.5*AliITSRecoParam::GetSPDdetzlength()};
    Double_t zmaxdead[3]={fSPDdetzcentre[1] - 0.5*AliITSRecoParam::GetSPDdetzlength(),
			  fSPDdetzcentre[2] - 0.5*AliITSRecoParam::GetSPDdetzlength(),
			  fSPDdetzcentre[3] - 0.5*AliITSRecoParam::GetSPDdetzlength()};
    for (Int_t i=0; i<3; i++)
      if (track->GetZ()-dz<zmaxdead[i] && track->GetZ()+dz>zmindead[i]) {
	AliDebug(2,Form("crack SPD %d",ilayer));
	return 1; 
      } 
  }

  // check bad zones from OCDB
  if (!AliITSReconstructor::GetRecoParam()->GetUseBadZonesFromOCDB()) return 0;

  if (idet<0) return 0;

  AliHLTITSDetector &det=fgLayers[ilayer].GetDetector(idet);  

  Int_t detType=-1;
  Float_t detSizeFactorX=0.0001,detSizeFactorZ=0.0001;
  if (ilayer==0 || ilayer==1) {        // ----------  SPD
    detType = 0;
  } else if (ilayer==2 || ilayer==3) { // ----------  SDD
    detType = 1;
    detSizeFactorX *= 2.;
  } else if (ilayer==4 || ilayer==5) { // ----------  SSD
    detType = 2;
  }
  AliITSsegmentation *segm = (AliITSsegmentation*)fkDetTypeRec->GetSegmentationModel(detType);
  if (detType==2) segm->SetLayer(ilayer+1);
  Float_t detSizeX = detSizeFactorX*segm->Dx(); 
  Float_t detSizeZ = detSizeFactorZ*segm->Dz(); 

  // check if the road overlaps with bad chips
  Float_t xloc,zloc;
  LocalModuleCoord(ilayer,idet,track,xloc,zloc);
  Float_t zlocmin = zloc-dz;
  Float_t zlocmax = zloc+dz;
  Float_t xlocmin = xloc-dy;
  Float_t xlocmax = xloc+dy;
  Int_t chipsInRoad[100];

  // check if road goes out of detector
  Bool_t touchNeighbourDet=kFALSE; 
  if (TMath::Abs(xlocmin)>0.5*detSizeX) {xlocmin=-0.5*detSizeX; touchNeighbourDet=kTRUE;} 
  if (TMath::Abs(xlocmax)>0.5*detSizeX) {xlocmax=+0.5*detSizeX; touchNeighbourDet=kTRUE;} 
  if (TMath::Abs(zlocmin)>0.5*detSizeZ) {zlocmin=-0.5*detSizeZ; touchNeighbourDet=kTRUE;} 
  if (TMath::Abs(zlocmax)>0.5*detSizeZ) {zlocmax=+0.5*detSizeZ; touchNeighbourDet=kTRUE;} 
  AliDebug(2,Form("layer %d det %d zmim zmax %f %f xmin xmax %f %f   %f %f",ilayer,idet,zlocmin,zlocmax,xlocmin,xlocmax,detSizeZ,detSizeX));

  // check if this detector is bad
  if (det.IsBad()) {
    AliDebug(2,Form("lay %d  bad detector %d",ilayer,idet));
    if(!touchNeighbourDet) {
      return 2; // all detectors in road are bad
    } else { 
      return 3; // at least one is bad
    }
  }

  Int_t nChipsInRoad = segm->GetChipsInLocalWindow(chipsInRoad,zlocmin,zlocmax,xlocmin,xlocmax);
  AliDebug(2,Form("lay %d nChipsInRoad %d",ilayer,nChipsInRoad));
  if (!nChipsInRoad) return 0;

  Bool_t anyBad=kFALSE,anyGood=kFALSE;
  for (Int_t iCh=0; iCh<nChipsInRoad; iCh++) {
    if (chipsInRoad[iCh]<0 || chipsInRoad[iCh]>det.GetNChips()-1) continue;
    AliDebug(2,Form("  chip %d bad %d",chipsInRoad[iCh],(Int_t)det.IsChipBad(chipsInRoad[iCh])));
    if (det.IsChipBad(chipsInRoad[iCh])) {
      anyBad=kTRUE;
    } else {
      anyGood=kTRUE;
    } 
  }

  if (!anyGood) {
    if(!touchNeighbourDet) {
      AliDebug(2,"all bad in road");
      return 2;  // all chips in road are bad
    } else {
      return 3; // at least a bad chip in road
    }
  }

  if (anyBad) {
    AliDebug(2,"at least a bad in road");
    return 3; // at least a bad chip in road
  } 


  if (!AliITSReconstructor::GetRecoParam()->GetUseSingleBadChannelsFromOCDB()
      || ilayer==4 || ilayer==5     // SSD
      || !noClusters) return 0;

  // There are no clusters in road: check if there is at least 
  // a bad SPD pixel or SDD anode 

  Int_t idetInITS=idet;
  for(Int_t l=0;l<ilayer;l++) idetInITS+=AliITSgeomTGeo::GetNLadders(l+1)*AliITSgeomTGeo::GetNDetectors(l+1);

  if (fITSChannelStatus->AnyBadInRoad(idetInITS,zlocmin,zlocmax,xlocmin,xlocmax)) {
    AliDebug(2,Form("Bad channel in det %d of layer %d\n",idet,ilayer));
    return 3;
  }
  //if (fITSChannelStatus->FractionOfBadInRoad(idet,zlocmin,zlocmax,xlocmin,xlocmax) > AliITSReconstructor::GetRecoParam()->GetMinFractionOfBadInRoad()) return 3;

  return 0;
}
//------------------------------------------------------------------------
Bool_t AliITStrackerHLT::LocalModuleCoord(Int_t ilayer,Int_t idet,
				       const AliHLTITSTrack *track,
				       Float_t &xloc,Float_t &zloc) const {
  //-----------------------------------------------------------------
  // Gives position of track in local module ref. frame
  //-----------------------------------------------------------------

  xloc=0.; 
  zloc=0.;

  if(idet<0) return kFALSE;

  Int_t ndet=AliITSgeomTGeo::GetNDetectors(ilayer+1); // layers from 1 to 6 

  Int_t lad = Int_t(idet/ndet) + 1;

  Int_t det = idet - (lad-1)*ndet + 1;

  Double_t xyzGlob[3],xyzLoc[3];

  AliHLTITSDetector &detector = fgLayers[ilayer].GetDetector(idet);
  // take into account the misalignment: xyz at real detector plane
  if(!track->GetXYZAt(detector.GetRmisal(),GetBz(),xyzGlob)) return kFALSE;

  if(!AliITSgeomTGeo::GlobalToLocal(ilayer+1,lad,det,xyzGlob,xyzLoc)) return kFALSE;

  xloc = (Float_t)xyzLoc[0];
  zloc = (Float_t)xyzLoc[2];

  return kTRUE;
}

//------------------------------------------------------------------------
Bool_t AliITStrackerHLT::ComputeRoad(AliHLTITSTrack* track,Int_t ilayer,Int_t idet,Double_t &zmin,Double_t &zmax,Double_t &ymin,Double_t &ymax) const {
  //--------------------------------------------------------------------
  // This function computes the rectangular road for this track
  //--------------------------------------------------------------------


  AliHLTITSDetector &det = fgLayers[ilayer].GetDetector(idet);
  // take into account the misalignment: propagate track to misaligned detector plane
  if (!TransportToPhiX( track, det.GetPhi(),det.GetRmisal() ) ) return kFALSE;

  Double_t dz=AliITSReconstructor::GetRecoParam()->GetNSigmaRoadZ()*
                    TMath::Sqrt(track->GetSigmaZ2() + 
		    AliITSReconstructor::GetRecoParam()->GetNSigmaZLayerForRoadZ()*
		    AliITSReconstructor::GetRecoParam()->GetNSigmaZLayerForRoadZ()*
		    AliITSReconstructor::GetRecoParam()->GetSigmaZ2(ilayer));
  Double_t dy=AliITSReconstructor::GetRecoParam()->GetNSigmaRoadY()*
                    TMath::Sqrt(track->GetSigmaY2() + 
		    AliITSReconstructor::GetRecoParam()->GetNSigmaYLayerForRoadY()*
		    AliITSReconstructor::GetRecoParam()->GetNSigmaYLayerForRoadY()*
		    AliITSReconstructor::GetRecoParam()->GetSigmaY2(ilayer));
      
  // track at boundary between detectors, enlarge road
  Double_t boundaryWidth=AliITSRecoParam::GetBoundaryWidth();
  if ( (track->GetY()-dy < det.GetYmin()+boundaryWidth) || 
       (track->GetY()+dy > det.GetYmax()-boundaryWidth) || 
       (track->GetZ()-dz < det.GetZmin()+boundaryWidth) ||
       (track->GetZ()+dz > det.GetZmax()-boundaryWidth) ) {
    Float_t tgl = TMath::Abs(track->GetTgl());
    if (tgl > 1.) tgl=1.;
    Double_t deltaXNeighbDets=AliITSRecoParam::GetDeltaXNeighbDets();
    dz = TMath::Sqrt(dz*dz+deltaXNeighbDets*deltaXNeighbDets*tgl*tgl);
    Float_t snp = TMath::Abs(track->GetSnp());
    if (snp > AliITSReconstructor::GetRecoParam()->GetMaxSnp()) return kFALSE;
    dy = TMath::Sqrt(dy*dy+deltaXNeighbDets*deltaXNeighbDets*snp*snp);
  } // boundary
  
  // add to the road a term (up to 2-3 mm) to deal with misalignments
  dy = TMath::Sqrt(dy*dy + AliITSReconstructor::GetRecoParam()->GetRoadMisal()*AliITSReconstructor::GetRecoParam()->GetRoadMisal());
  dz = TMath::Sqrt(dz*dz + AliITSReconstructor::GetRecoParam()->GetRoadMisal()*AliITSReconstructor::GetRecoParam()->GetRoadMisal());

  Double_t r = fgLayers[ilayer].GetR();
  zmin = track->GetZ() - dz; 
  zmax = track->GetZ() + dz;
  ymin = track->GetY() + r*det.GetPhi() - dy;
  ymax = track->GetY() + r*det.GetPhi() + dy;

  // bring track back to idead detector plane
  if (!TransportToPhiX( track, det.GetPhi(),det.GetR())) return kFALSE;

  return kTRUE;
}
