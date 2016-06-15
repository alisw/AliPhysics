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
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliV0.h"
#include "AliHelix.h"
#include "AliITSChannelStatus.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSReconstructor.h"
#include "AliITSClusterParam.h"
#include "AliITSsegmentation.h"
#include "AliITSCalibration.h"
#include "AliITSV0Finder.h"
#include "AliITStrackerHLT.h"
#include "TStopwatch.h"
//#include "AliHLTTPCCATrackParam.h"
//#include "AliHLTVertexer.h"
#include <vector>

using std::vector;
ClassImp(AliITStrackerHLT)

Bool_t AliITStrackerHLT::CheckTrack( const AliExternalTrackParam *t )
{
  // Check that the track parameters are reasonable in order to avoid floating point exeptions in AliExternalTrackParam routines
  
  bool ok = 1;
  for( int i=0; i<5; i++ ) ok = ok && finite(t->GetParameter()[i]);
  for( int i=0; i<15; i++ ) ok = ok && finite(t->GetCovariance()[i]);
  ok = ok && ( TMath::Abs(t->GetX()) < 500. );
  ok = ok && ( TMath::Abs(t->GetY()) < 500. );
  ok = ok && ( TMath::Abs(t->GetZ()) < 500. );
  ok = ok && ( TMath::Abs(t->GetSnp()) < .99 );
  ok = ok && ( TMath::Abs(t->GetSigned1Pt()) < 1./0.01 );
  return ok;
}

Bool_t AliITStrackerHLT::TransportToX( AliExternalTrackParam *t, double x ) const
{
  return CheckTrack(t) && t->PropagateTo( x, t->GetBz() );
}

Bool_t AliITStrackerHLT::TransportToPhiX( AliExternalTrackParam *t, double phi, double x ) const
{
  return CheckTrack(t) && t->Propagate( phi, x, t->GetBz() );
}



AliITStrackerHLT::AliITStrackerHLT()
  :AliTracker(),
   fRecoParam(0),
   fLayers(new AliHLTITSLayer[AliITSgeomTGeo::kNLayers]),
   fUseTGeo(2),
   fxOverX0Pipe(-1.),
   fxTimesRhoPipe(-1.), 
   fTracks(0),
   fITSOutTracks(0),
   fNTracks(0),
   fNITSOutTracks(0),
   fLoadTime(0),
   fRecoTime(0),
   fNEvents(0),
   fClusters(0),
   fNClusters(0)
{
  //Default constructor
  Int_t i;
  for(i=0;i<4;i++) fSPDdetzcentre[i]=0.;
  for(i=0;i<2;i++) {fxOverX0Shield[i]=-1.;fxTimesRhoShield[i]=-1.;}
  for(i=0;i<6;i++) {fxOverX0Layer[i]=-1.;fxTimesRhoLayer[i]=-1.;}
}
//------------------------------------------------------------------------
AliITStrackerHLT::AliITStrackerHLT(const Char_t *geom) 
: AliTracker(),
  fRecoParam(0),
  fLayers(new AliHLTITSLayer[AliITSgeomTGeo::kNLayers]),  
  fUseTGeo(2),
  fxOverX0Pipe(-1.),
  fxTimesRhoPipe(-1.),
  fTracks(0),
  fITSOutTracks(0),
  fNTracks(0),
  fNITSOutTracks(0),
  fLoadTime(0),
   fRecoTime(0),
  fNEvents(0),
  fClusters(0),
  fNClusters(0)
{
  //--------------------------------------------------------------------
  //This is the AliITStrackerHLT constructor
  //--------------------------------------------------------------------
  if (geom) {
    AliWarning("\"geom\" is actually a dummy argument !");
  }

  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }

  fRecoParam = AliITSReconstructor::GetRecoParam();
  if( !fRecoParam ){
    AliITSReconstructor *tmp = new AliITSReconstructor();
    tmp->Init();
    fRecoParam = AliITSRecoParam::GetLowFluxParam();
    tmp->AliReconstructor::SetRecoParam(fRecoParam);
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

    new (fLayers+i-1) AliHLTITSLayer(r,poff,zoff,nlad,ndet);

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

        AliHLTITSDetector &det=fLayers[i-1].GetDetector((j-1)*ndet + k-1); 
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

  
  Double_t xyzVtx[]={ fRecoParam->GetXVdef(),
		      fRecoParam->GetYVdef(),
		      fRecoParam->GetZVdef()}; 
  Double_t ersVtx[]={ fRecoParam->GetSigmaXVdef(),
		      fRecoParam->GetSigmaYVdef(),
		      fRecoParam->GetSigmaZVdef()}; 

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

  //fUseTGeo = fRecoParam->GetUseTGeoInTracker();
  //if(fRecoParam->GetExtendedEtaAcceptance() && fUseTGeo!=1 && fUseTGeo!=3) {
  //AliWarning("fUseTGeo changed to 3 because fExtendedEtaAcceptance is kTRUE");
  //fUseTGeo = 3;
  //}

  for(Int_t i=0;i<2;i++) {fxOverX0Shield[i]=-1.;fxTimesRhoShield[i]=-1.;}
  for(Int_t i=0;i<6;i++) {fxOverX0Layer[i]=-1.;fxTimesRhoLayer[i]=-1.;}
  
  Init();
}
//------------------------------------------------------------------------
AliITStrackerHLT::AliITStrackerHLT(const AliITStrackerHLT &tracker)
:AliTracker(tracker),
 fRecoParam( tracker.fRecoParam),
 fLayers(new AliHLTITSLayer[AliITSgeomTGeo::kNLayers]), 
 fUseTGeo(tracker.fUseTGeo),
 fxOverX0Pipe(tracker.fxOverX0Pipe),
 fxTimesRhoPipe(tracker.fxTimesRhoPipe), 
 fTracks(0),
 fITSOutTracks(0),
 fNTracks(0),
 fNITSOutTracks(0),
  fLoadTime(0),
   fRecoTime(0),
 fNEvents(0),
 fClusters(0),
 fNClusters(0)
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
  Init();
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
  delete[] fLayers;
  delete[] fTracks;
  delete[] fITSOutTracks;
  delete[] fClusters;
}

void AliITStrackerHLT::Init()
{
  BuildMaterialLUT("Layers");  
  BuildMaterialLUT("Pipe");  
  BuildMaterialLUT("Shields");  
}


void AliITStrackerHLT::StartLoadClusters( Int_t NOfClusters )
{
  // !
  delete[] fClusters;
  fClusters = new AliITSRecPoint[NOfClusters];
  fNClusters = 0;
}

void AliITStrackerHLT::LoadCluster( const AliITSRecPoint &cluster) 
{
  fClusters[fNClusters++] = cluster ;
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

  int nClustersTotal = 0;
  {
    Int_t j=0;
    for (int i=0; i<AliITSgeomTGeo::GetNLayers(); i++) {
      int ndet=fLayers[i].GetNdetectors();
      Int_t jmax = j + fLayers[i].GetNladders()*ndet;
      for (; j<jmax; j++) {           
	if (!cTree->GetEvent(j)) continue;
	nClustersTotal+=clusters->GetEntriesFast();      
	clusters->Delete();
      }
    }
  }
  StartLoadClusters(nClustersTotal);
  {
    Int_t j=0;
    for (int i=0; i<AliITSgeomTGeo::GetNLayers(); i++) {
      int ndet=fLayers[i].GetNdetectors();
      Int_t jmax = j + fLayers[i].GetNladders()*ndet;
      for (; j<jmax; j++) {           
	if (!cTree->GetEvent(j)) continue;
	Int_t ncl=clusters->GetEntriesFast(); 
	while (ncl--) {
	  LoadCluster( *( (AliITSRecPoint*)clusters->UncheckedAt(ncl)));
	}
	clusters->Delete();
      }
    }
  }

  dummy.Clear();

  return 0;
}

//------------------------------------------------------------------------
void AliITStrackerHLT::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<AliITSgeomTGeo::GetNLayers(); i++) fLayers[i].ResetClusters();
  delete[] fClusters;
  fClusters = 0;
  fNClusters=0;
}




void AliITStrackerHLT::Reconstruct( AliExternalTrackParam *tracksTPC, int *tracksTPCLab, int nTPCTracks )
{

  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  //--------------------------------------------------------------------

  fNEvents++;

  // Init clusters
 
  TStopwatch timerInit;

  for( int i=0; i<AliITSgeomTGeo::GetNLayers(); i++ ){ 
    fLayers[i].ResetClusters();
  }

  for( int icl=0; icl<fNClusters; icl++ ){   
    AliITSRecPoint &cl = fClusters[icl];
    if (!cl.Misalign()) AliWarning("Can't misalign this cluster !"); 
    fLayers[cl.GetLayer()].InsertCluster(&cl); 
  }

  for( int i=0; i<AliITSgeomTGeo::GetNLayers(); i++ ){ 
    fLayers[i].ResetRoad(); //road defined by the cluster density
    fLayers[i].SortClusters();
  }  
  timerInit.Stop();
  fLoadTime+=timerInit.RealTime();


  TStopwatch timer;

  Double_t pimass = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  delete[] fTracks;
  delete[] fITSOutTracks;
  fTracks = new AliHLTITSTrack[nTPCTracks];
  fITSOutTracks = new AliHLTITSTrack[nTPCTracks];
  fNTracks = 0;
  fNITSOutTracks = 0;
  for( int itr=0; itr<nTPCTracks; itr++ ){    

    AliHLTITSTrack tMI( tracksTPC[itr] );
    AliHLTITSTrack *t = &tMI;
    t->SetTPCtrackId( itr );
    t->SetMass(pimass); 
    t->SetExpQ(0);
    t->SetLabel(tracksTPCLab[itr]);

    //if (!CorrectForTPCtoITSDeadZoneMaterial(t))  continue;
      
    Int_t tpcLabel=t->GetLabel(); //save the TPC track label       
    
    FollowProlongationTree(t); 
    int nclu=0;
    for(Int_t i=0; i<6; i++) {
      if( t->GetClusterIndex(i)>=0 ) nclu++; 
    }
    //cout<<"N assigned ITS clusters = "<<nclu<<std::endl;
    t->SetLabel(-1);
    if( nclu>0 ){
      t->SetLabel(tpcLabel);
      t->SetFakeRatio(1.);
      CookLabel(t,.99); //For comparison only
      //cout<<"SG: label = "<<t->GetLabel()<<" / "<<tpcLabel<<endl;
    }

    CorrectForPipeMaterial(t);
   
    TransportToX(t, 0 );
    fTracks[fNTracks++] = *t;  
    //cout<<"SG: ITS: Bz = "<<t->GetBz()<<endl;

    if(  nclu>0 ){ // construct ITSOut track
      AliHLTITSTrack tOut(*t);
      if( FitOutward( &tOut ) ){
	fITSOutTracks[fNITSOutTracks++] = *t;  
      }
    }
  }

  timer.Stop();
  fRecoTime+=timer.RealTime();
}



//------------------------------------------------------------------------
Int_t AliITStrackerHLT::Clusters2Tracks(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------
  
  
  std::vector<AliExternalTrackParam> tracksTPC;
  std::vector<int> tracksTPCLab;
  tracksTPC.reserve(event->GetNumberOfTracks());

  for( int itr=0; itr<event->GetNumberOfTracks(); itr++ ){

    AliESDtrack *esdTrack = event->GetTrack(itr);
    //esdTrack->myITS = esdTrack->myTPC;
    if ((esdTrack->GetStatus()&AliESDtrack::kTPCin)==0) continue;
    //if (esdTrack->GetStatus()&AliESDtrack::kTPCout) continue;
    if (esdTrack->GetStatus()&AliESDtrack::kITSin) continue;
    if (esdTrack->GetKinkIndex(0)>0) continue;   //kink daughter
    
    AliHLTITSTrack t(*esdTrack);
    t.SetTPCtrackId( itr );
    tracksTPC.push_back( t );
    tracksTPCLab.push_back(esdTrack->GetLabel());
  }
  //for( int iter=0; iter<100; iter++){
  Reconstruct( &(tracksTPC[0]), &(tracksTPCLab[0]), tracksTPC.size() );
  //}

  for( int itr=0; itr<fNTracks; itr++ ){
    AliHLTITSTrack &t = fTracks[itr];    
    UpdateESDtrack(event->GetTrack(t.TPCtrackId()), &t, AliESDtrack::kITSin);          
    //event->GetTrack(t.TPCtrackId())->myITS = t;
  }
 

  //int hz = ( int ) ( (fRecoTime+fLoadTime) > 1.e-4 ? fNEvents / (fRecoTime+fLoadTime) : 0 );
  //int hz1 = ( int ) ( fRecoTime > 1.e-4 ? fNEvents / fRecoTime : 0 );
  //int hz2 = ( int ) ( fLoadTime > 1.e-4 ? fNEvents / fLoadTime : 0 );

  //std::cout<<"\n\nSG: ITS tracker time = "<<hz2<<" Hz load / "<<hz1<<" Hz reco ="
  //<<hz<<
  //" Hz ("
  //<<fLoadTime/fNEvents*1000<<"+"<<fRecoTime/fNEvents*1000.
  //<<" = "<<(fLoadTime + fRecoTime)/fNEvents*1000. 
  //<<" ms/ev), "<<fNEvents<<" events processed\n\n "<<std::endl;
  return 0;
}


AliCluster *AliITStrackerHLT::GetCluster(Int_t index) const 
{
  //       Return pointer to a given cluster
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fLayers[l].GetCluster(c);
}




//------------------------------------------------------------------------
void AliITStrackerHLT::FollowProlongationTree(AliHLTITSTrack * track ) 
{
  // FollowProlongationTree
  for (Int_t ilayer=5; ilayer>=0; ilayer--) {
   
    AliHLTITSLayer &layer=fLayers[ilayer];
  
    // material between SSD and SDD, SDD and SPD
    //if (ilayer==3 && !CorrectForShieldMaterial(track,1)) continue;
    //if (ilayer==1 && !CorrectForShieldMaterial(track,0)) continue;
    
    int idet;

    {            
      Double_t xloc, phi,z;
      if( !track->GetLocalXPhiZat( layer.GetR(), xloc, phi, z ) ) return;
      idet = layer.FindDetectorIndex(phi,z);
    }

    // track outside layer acceptance in z
   
    if( idet<0 ) continue;
    
    // propagate to the intersection with the detector plane     
    {
      const AliHLTITSDetector &det=layer.GetDetector( idet );
      if (!TransportToPhiX( track, det.GetPhi(), det.GetR() ) ) return;
      CorrectForLayerMaterial(track,ilayer);
    }

    // DEFINITION OF SEARCH ROAD AND CLUSTERS SELECTION
    
    // road in global (rphi,z) [i.e. in tracking ref. system]
    
    Double_t zmin,zmax,ymin,ymax;

    if (!ComputeRoad(track,ilayer,idet,zmin,zmax,ymin,ymax)) continue;
  
    layer.SelectClusters(zmin,zmax,ymin,ymax);     
    
    // Define criteria for track-cluster association
    
    Double_t msz = track->GetSigmaZ2() + 
      fRecoParam->GetNSigmaZLayerForRoadZ()*
      fRecoParam->GetNSigmaZLayerForRoadZ()*
      fRecoParam->GetSigmaZ2(ilayer);
    
    Double_t msy = track->GetSigmaY2() + 
      fRecoParam->GetNSigmaYLayerForRoadY()*
      fRecoParam->GetNSigmaYLayerForRoadY()*
      fRecoParam->GetSigmaY2(ilayer);
    
    msz *= fRecoParam->GetNSigma2RoadZNonC();
    msy *= fRecoParam->GetNSigma2RoadYNonC(); 
    
    msz = 1./msz; // 1/RoadZ^2
    msy = 1./msy; // 1/RoadY^2    
    
    const AliITSRecPoint *cl=0; 
    Int_t clidx=-1;     
    
    // loop over clusters in the road     
    
    const AliITSRecPoint *bestCluster=0; 
    double bestChi2 = 1.e10;
    AliHLTITSTrack bestTrack( *track );
    int bestIdx = -1;
     
    while( (cl=layer.GetNextCluster(clidx)) ){
      Int_t idetc=cl->GetDetectorIndex();
      if ( idet !=idetc ) { // new cluster's detector
	const AliHLTITSDetector &detc=layer.GetDetector(idetc);
	if (!TransportToPhiX( track, detc.GetPhi(),detc.GetR()) ) continue;
	idet = idetc;
      }  
      //double y,z;
      //if (! track->GetLocalYZat( layer.GetDetector(idetc).GetR() + cl->GetX(),y,z ) ) continue;
      double dz = track->GetZ() - cl->GetZ();
      double dy = track->GetY() - cl->GetY();
      double chi2 = dz*dz*msz + dy*dy*msy ;       
      if ( chi2 < bestChi2 ){
	bestChi2 = chi2;
	bestCluster = cl;
	bestTrack = *track;
	bestIdx = clidx;
	continue;
      }
    }
    
    if( !bestCluster || bestChi2 >2*10. ) continue;
    
    if (!TransportToX( &bestTrack, layer.GetDetector(bestCluster->GetDetectorIndex()).GetR() + bestCluster->GetX() ) ) continue;
    
    Double_t par[2]={ bestCluster->GetY(), bestCluster->GetZ()};
    Double_t cov[3]={ bestCluster->GetSigmaY2(), 0., bestCluster->GetSigmaZ2()};
    if( !bestTrack.AliExternalTrackParam::Update(par,cov) ) continue;
    
    *track = bestTrack;
    track->SetClusterIndex(track->GetNumberOfClusters(), (ilayer<<28)+bestIdx);
    track->SetNumberOfClusters(track->GetNumberOfClusters()+1);  
  }
}



Int_t AliITStrackerHLT::FitOutward(AliHLTITSTrack * track ) 
{
  // FitOutward
  track->ResetCovariance(100);

  for (Int_t iTrCl=track->GetNumberOfClusters()-1; iTrCl>=0; iTrCl--) {
    
    Int_t index = track->GetClusterIndex(iTrCl);
    Int_t ilayer=(index & 0xf0000000) >> 28;
    Int_t ic=(index & 0x0fffffff) >> 00;
    const AliHLTITSLayer &layer=fLayers[ilayer];
    AliITSRecPoint *cl = layer.GetCluster(ic); 
    int idet = cl->GetDetectorIndex();
    const AliHLTITSDetector &det=layer.GetDetector( idet );
 
    // material between SSD and SDD, SDD and SPD
    //if (ilayer==4 && !CorrectForShieldMaterial(track,1)) continue;
    //if (ilayer==2 && !CorrectForShieldMaterial(track,0)) continue;
    

    // propagate to the intersection with the detector plane     
    {
      if (!TransportToPhiX( track, det.GetPhi(), det.GetR()+ cl->GetX() ) ) return 0;
      CorrectForLayerMaterial(track,ilayer);
    }

    Double_t par[2]={ cl->GetY(), cl->GetZ()};
    Double_t cov[3]={ cl->GetSigmaY2(), 0., cl->GetSigmaZ2()};
    if( !track->AliExternalTrackParam::Update(par,cov) ) return 0;    
  }
  return 1;
}


//------------------------------------------------------------------------
AliHLTITSLayer & AliITStrackerHLT::GetLayer(Int_t layer) const
{
  //--------------------------------------------------------------------
  //
  //
  return fLayers[layer];
}



//------------------------------------------------------------------------
void AliITStrackerHLT::CookLabel(AliHLTITSTrack *track,Float_t wrong) const 
{
  // get MC label for the track

  Int_t mcLabel = -1;
  
  vector<int> labels;
  Int_t nClusters = track->GetNumberOfClusters();
  Int_t nClustersEff = 0;
  for (Int_t i=0; i<nClusters; i++){
    Int_t cindex = track->GetClusterIndex(i);
    //Int_t l=(cindex & 0xf0000000) >> 28;
    AliITSRecPoint *cl = (AliITSRecPoint*)GetCluster(cindex);    
    if ( cl->GetLabel(0) >= 0 ){ labels.push_back(cl->GetLabel(0)) ; nClustersEff++; }
    if ( cl->GetLabel(1) >= 0 ) labels.push_back(cl->GetLabel(1)) ;
    if ( cl->GetLabel(2) >= 0 ) labels.push_back(cl->GetLabel(2)) ;
  }
  std::sort( labels.begin(), labels.end() );

  labels.push_back( -1 ); // put -1 to the end	  
  int labelMax = -1, labelCur = -1, nLabelsMax = 0, nLabelsCurr = 0;

  for ( unsigned int iLab = 0; iLab < labels.size(); iLab++ ) {
    if ( labels[iLab] != labelCur ) {	      
      if ( labelCur >= 0 && nLabelsMax< nLabelsCurr ) {
	nLabelsMax = nLabelsCurr;
	labelMax = labelCur;
      }
      labelCur = labels[iLab];
      nLabelsCurr = 0;
    }
    nLabelsCurr++;
  }

  if( labelMax>=0 && nLabelsMax < wrong * nClustersEff ) labelMax = -labelMax;
  
  mcLabel = labelMax;
		
  track->SetLabel( mcLabel );
}









//------------------------------------------------------------------------
void AliITStrackerHLT::UpdateESDtrack(AliESDtrack *tESD, AliHLTITSTrack* track, ULong_t flags) const
{
  //
  // Update ESD track
  //
  tESD->UpdateTrackParams(track,flags);
  AliHLTITSTrack * oldtrack = (AliHLTITSTrack*)(tESD->GetITStrack());
  if (oldtrack) delete oldtrack; 
  tESD->SetITStrack(new AliHLTITSTrack(*track));
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


//------------------------------------------------------------------------
Int_t AliITStrackerHLT::CorrectForPipeMaterial(AliHLTITSTrack *t,
					       bool InwardDirection) {
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
  
  Float_t  dir = (InwardDirection ? 1. : -1.);
  Double_t rToGo= ( InwardDirection ? AliITSRecoParam::GetrInsidePipe() : AliITSRecoParam::GetrOutsidePipe());
  Double_t xToGo, phi,z;

  if (!t->GetLocalXPhiZat(rToGo,xToGo,phi,z)) return 0;

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
						 Int_t    shieldindex,
						bool InwardDirection) {
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


  Float_t  dir = (InwardDirection ? 1. : -1.);
  Double_t rToGo;

  if (shieldindex==1 ) { // SDDouter
    rToGo=(InwardDirection ? AliITSRecoParam::GetrInsideShield(1) : AliITSRecoParam::GetrOutsideShield(1));
  } else if (shieldindex==0 ) {        // SPDouter
    rToGo=(InwardDirection ? AliITSRecoParam::GetrInsideShield(0) : AliITSRecoParam::GetrOutsideShield(0)); 
  } else {
    Error("CorrectForShieldMaterial"," Wrong shield name\n");
    return 0;
  }
  Double_t xToGo, phi,z;

  if (!t->GetLocalXPhiZat(rToGo,xToGo,phi,z)) return 0;

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
    nsteps= (Int_t)(TMath::Abs(t->GetX()-xToGo)/fRecoParam->GetStepSizeTGeo())+1;
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
					       bool InwardDirection ){
  //-------------------------------------------------------------------
  // Propagate beyond layer and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  /*
    Double_t r=fLayers[layerindex].GetR();
    Double_t deltar=(layerindex<2 ? 0.10*r : 0.05*r);
    Double_t rToGo=TMath::Sqrt(t->GetX()*t->GetX()+t->GetY()*t->GetY());
    rToGo+= InwardDirection ?-deltar :deltar;
    Double_t xToGo;
    if (!t->GetLocalXat(rToGo,xToGo)) return 0;  
  */
  
  if(fxOverX0Layer[layerindex]<0) BuildMaterialLUT("Layers");

  Double_t lengthTimesMeanDensity = fxTimesRhoLayer[layerindex];
  if( !InwardDirection ) lengthTimesMeanDensity = -lengthTimesMeanDensity;

  return t->CorrectForMeanMaterial(fxOverX0Layer[layerindex],lengthTimesMeanDensity,kTRUE);
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

  AliHLTITSDetector &detector = fLayers[ilayer].GetDetector(idet);
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


  AliHLTITSDetector &det = fLayers[ilayer].GetDetector(idet);
  // take into account the misalignment: propagate track to misaligned detector plane

  double y,z,snp,cov[3];
  if( !track->GetYZAtPhiX( det.GetPhi(),det.GetRmisal(), y, z, snp, cov))return 0;


  Double_t dz=fRecoParam->GetNSigmaRoadZ()*
                    TMath::Sqrt(cov[2] + 
				fRecoParam->GetNSigmaZLayerForRoadZ()*
				fRecoParam->GetNSigmaZLayerForRoadZ()*
				fRecoParam->GetSigmaZ2(ilayer));
  Double_t dy=fRecoParam->GetNSigmaRoadY()*
                    TMath::Sqrt(cov[0] + 
				fRecoParam->GetNSigmaYLayerForRoadY()*
				fRecoParam->GetNSigmaYLayerForRoadY()*
				fRecoParam->GetSigmaY2(ilayer));
      
  // track at boundary between detectors, enlarge road
  Double_t boundaryWidth=AliITSRecoParam::GetBoundaryWidth();
  if ( (y-dy < det.GetYmin()+boundaryWidth) || 
       (y+dy > det.GetYmax()-boundaryWidth) || 
       (z-dz < det.GetZmin()+boundaryWidth) ||
       (z+dz > det.GetZmax()-boundaryWidth) ) {
    Float_t tgl = TMath::Abs(track->GetTgl());
    if (tgl > 1.) tgl=1.;
    Double_t deltaXNeighbDets=AliITSRecoParam::GetDeltaXNeighbDets();
    dz = TMath::Sqrt(dz*dz+deltaXNeighbDets*deltaXNeighbDets*tgl*tgl);
    if (snp > fRecoParam->GetMaxSnp()) return kFALSE;
    dy = TMath::Sqrt(dy*dy+deltaXNeighbDets*deltaXNeighbDets*snp*snp);
  } // boundary
  
  // add to the road a term (up to 2-3 mm) to deal with misalignments
  dy = TMath::Sqrt(dy*dy + fRecoParam->GetRoadMisal()*fRecoParam->GetRoadMisal());
  dz = TMath::Sqrt(dz*dz + fRecoParam->GetRoadMisal()*fRecoParam->GetRoadMisal());

  Double_t r = fLayers[ilayer].GetR();
  zmin = z - dz; 
  zmax = z + dz;
  ymin = y + r*det.GetPhi() - dy;
  ymax = y + r*det.GetPhi() + dy;

  return kTRUE;
}



Int_t AliITStrackerHLT::PropagateBack(AliESDEvent * /*event*/) 
{
  // dummy
  return 0;
}

Int_t AliITStrackerHLT::RefitInward(AliESDEvent * /*event*/ ) 
{
  // dummy
  return 0;
}


Bool_t AliITStrackerHLT::GetTrackPoint(Int_t /*index*/, AliTrackPoint& /*p*/) const 
{
  // dummy
  return 0;
}

Bool_t AliITStrackerHLT::GetTrackPointTrackingError(Int_t /*index*/, 
						    AliTrackPoint& /*p*/, const AliESDtrack */*t*/) 
{
  // dummy
  return 0;
}
