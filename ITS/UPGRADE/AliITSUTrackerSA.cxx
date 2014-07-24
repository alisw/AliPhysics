//--------------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    It reads AliITSUClusterPix clusters and and fills the ESD with tracks
//    
//    The algorithm implemented here takes inspiration from UniCA code of FIAS
//    group. 
//--------------------------------------------------------------------------------

#include <TBranch.h>
#include <TMath.h>
using TMath::Abs;
using TMath::Sort;
using TMath::Sqrt;
#include <TTree.h>
#include <algorithm>
using std::sort;

// Vc library
//#include "Vc/Vc"
//#include "AliITSUTrackerSAauxVc.h" // Structs and other stuff using Vc library
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliITSUClusterPix.h"
#include "AliITSUTrackerSA.h"
#include "AliITSUReconstructor.h"
#include "AliITSURecoDet.h"
#include "AliESDtrack.h"

#include <Riostream.h>

using std::cout;
using std::endl;
using std::flush;

//#include "AliITSUtrackSA.h"      // Some dedicated SA track class ?

ClassImp(AliITSUTrackerSA)

const Double_t AliITSUTrackerSA::fgkToler =  1e-6;// tolerance for layer on-surface check
const Double_t AliITSUTrackerSA::fgkChi2Cut =  600.f;

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA(AliITSUReconstructor* rec) :
fReconstructor(rec),
fITS(0),
fMatLUT(0),
fUseMatLUT(kFALSE),
fCurrMass(0.14),
//
fClustersTC(),
fChi2Cut( fgkChi2Cut ),
fPhiCut( 1  ),
fRPhiCut( 1 ),
fZCut( 1 )
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------
  if (rec) Init(rec);
}

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA(const AliITSUTrackerSA &t): AliTracker(t),
                                                               fReconstructor(t.fReconstructor),
                                                               fITS(t.fITS),
                                                               fMatLUT(t.fMatLUT),
                                                               fUseMatLUT(t.fUseMatLUT),
                                                               fCurrMass(t.fCurrMass),
  //
                                                               fClustersTC(),
                                                               fChi2Cut(fgkChi2Cut),
                                                               fPhiCut(),
                                                               fRPhiCut(),
                                                               fZCut()
{
  //--------------------------------------------------------------------
  // The copy constructor is protected
  //--------------------------------------------------------------------
}

//________________________________________________________________________________
AliITSUTrackerSA::~AliITSUTrackerSA()
{
  // d-tor
  delete fMatLUT;
}


//_________________________________________________________________________
void AliITSUTrackerSA::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  //
  fITS = rec->GetITSInterface();
  //
  // create material lookup table
  const int kNTest = 1000;
  const double kStepsPerCM=5;
  fMatLUT  = new AliITSUMatLUT(fITS->GetRMin(),fITS->GetRMax(),Nint(kStepsPerCM*(fITS->GetRMax()-fITS->GetRMin())));
  double zmn = 1e6;
  for (int ilr=fITS->GetNLayers();ilr--;) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if (zmn>Abs(lr->GetZMin())) zmn = Abs(lr->GetZMin());
    if (zmn>Abs(lr->GetZMax())) zmn = Abs(lr->GetZMax());
  }
  fMatLUT->FillData(kNTest,-zmn,zmn);
  //
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::Clusters2Tracks(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This is the main tracking function
  // The clusters must already be loaded
  //--------------------------------------------------------------------

  // Possibly, create the track "seeds" (combinatorial)

  // Possibly, increment the seeds with additional clusters (Kalman)

  // Possibly, (re)fit the found tracks

  // Three iterations:
  // - High momentum first;
  // - Low momentum with vertex constraint;
  // - Everything else.

  CellsCreation(0); 
  CellularAutomaton(event);
  // VertexFinding();
  // CellsCreation(1);
  // CellularAutomaton(event);

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::PropagateBack(AliESDEvent * event) {
  //--------------------------------------------------------------------
  // Here, we implement the Kalman smoother ?
  // The clusters must be already loaded
  //--------------------------------------------------------------------
  Int_t n=event->GetNumberOfTracks();
  Int_t ntrk=0;
  Int_t ngood=0;
  for (Int_t i=0; i<n; i++) {
    AliESDtrack *esdTrack=event->GetTrack(i);

    if ((esdTrack->GetStatus()&AliESDtrack::kITSin)==0) continue;

    AliITSUTrackCooked track(*esdTrack);

    track.ResetCovariance(10.); 

    int points[2*AliITSUAux::kMaxLayers];
    for (UInt_t k=0; k<2*AliITSUAux::kMaxLayers; k++) 
      points[k]=-1;
    Int_t nc=track.GetNumberOfClusters();
    for (Int_t k=0; k<nc; k++) {
      const int layer = (track.GetClusterIndex(k)&0xf0000000)>>28;
      const int idx = (track.GetClusterIndex(k)&0x0fffffff);
      points[layer<<1]=idx;
    }

    if (RefitTrack(&track,points,40,1)>=0) {

      CookLabel(&track, 0.); //For comparison only
      Int_t label=track.GetLabel();
      if (label>0) ngood++;

      esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSout);
      ntrk++;
    }
  }

  Info("PropagateBack","Back propagated tracks: %d",ntrk);
  if (ntrk)
    Info("PropagateBack","Good tracks/back propagated: %f",Float_t(ngood)/ntrk);

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::RefitInward(AliESDEvent * event) {
  //--------------------------------------------------------------------
  // Some final refit, after the outliers get removed by the smoother ?
  // The clusters must be loaded
  //--------------------------------------------------------------------
  Int_t n=event->GetNumberOfTracks();
  Int_t ntrk=0;
  Int_t ngood=0;
  for (Int_t i=0; i<n; i++) {
    AliESDtrack *esdTrack=event->GetTrack(i);

    if ((esdTrack->GetStatus()&AliESDtrack::kITSout)==0) continue;

    AliITSUTrackCooked track(*esdTrack);

    track.ResetCovariance(10.); 

    int points[2*AliITSUAux::kMaxLayers];
    for (UInt_t k=0; k<2*AliITSUAux::kMaxLayers; k++) 
      points[k]=-1;
    Int_t nc=track.GetNumberOfClusters();
    for (Int_t k=0; k<nc; k++) {
      const int layer = (track.GetClusterIndex(k)&0xf0000000)>>28;
      const int idx = (track.GetClusterIndex(k)&0x0fffffff);
      points[layer<<1]=idx;
    }

    if (RefitTrack(&track,points,1.8,1)>=0) { //2.1,1)>=0) {

      //if (!track.PropagateTo(1.8, 2.27e-3, 35.28*1.848)) continue;
      CookLabel(&track, 0.); //For comparison only
      Int_t label=track.GetLabel();
      if (label>0) ngood++;

      //cout << esdTrack->GetStatus() << " ";
      esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSrefit);
      //cout << esdTrack->GetStatus() << endl;
      ntrk++;
    } 
  }

  Info("RefitInward","Refitted tracks: %d",ntrk);
  if (ntrk)
    Info("RefitInward","Good tracks/refitted: %f",Float_t(ngood)/ntrk);
    
  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::LoadClusters(TTree *cluTree) {
  //--------------------------------------------------------------------
  // This function reads the ITSU clusters from the tree,
  // sort them, distribute over the internal tracker arrays, etc
  //--------------------------------------------------------------------
  fITS->LoadClusters(cluTree);
  fITS->ProcessClusters();
  //
  for(int iL=0; iL<7; ++iL) {
    fClustersTC[iL]=*fITS->GetLayerActive(iL)->GetClustersAddress();
    AliITSURecoLayer* lr = fITS->GetLayerActive(iL) ; // assign the layer which the cluster belongs to
    for(int iC=0;iC<fClustersTC[iL]->GetEntriesFast();++iC) {
      const AliITSUClusterPix *cl = (AliITSUClusterPix*)fClustersTC[iL]->At(iC);
      float pos[3];
      cl->GetGlobalXYZ(pos);
      float phi = TMath::PiOver2(); 
      if(Abs(pos[0])>1e-9) {
        phi=TMath::ATan2(pos[1]-GetY(),pos[0]-GetX());
        if(phi<0.f) phi+=TMath::TwoPi();
      } else if(pos[1]<0.f) phi *= 3.f ;
      AliITSURecoSens* sens = lr->GetSensorFromID(cl->GetVolumeId());
      const float alpha = sens->GetPhiTF();
      const float cov[3]={cl->GetSigmaZ2(),cl->GetSigmaYZ(),cl->GetSigmaY2()};
      
      fLayer[iL].AddPoint(pos,cov,phi,alpha);
      for ( int i=0 ; i<3; ++i ) {
        fLayer[iL].Points.back().Label[i] = (cl->GetLabel(i)<0) ? -1 : cl->GetLabel(i);
      }
    }
  
    fLayer[iL].Sort(); 

  }
  return 0;
}

//________________________________________________________________________________
void AliITSUTrackerSA::UnloadClusters() {
  //--------------------------------------------------------------------
  // This function unloads ITSU clusters from the RAM
  //--------------------------------------------------------------------
  for(int i=0;i<7;++i) 
    fLayer[i].Clear();
  for(int i=0;i<5;++i) 
    fCells[i].clear();

}

//________________________________________________________________________________
void AliITSUTrackerSA::CellularAutomaton(AliESDEvent *event) {

  // Here it's implemented the Cellular Automaton routine
  // Firstly the level of each doublet is set according to the level of
  // the neighbour doublets.
  // Doublet are considered to be neighbour if they share one point and the
  // phi and theta direction difference of the two is below a cut value.
  Int_t ntrk=0,ngood=0;

  for( int iL = 1; iL < 5; ++iL ) {
    for( size_t iC1 = 0; iC1 < fCells[iL-1].size(); ++iC1 ) {
      for( size_t iC2 = 0; iC2 < fCells[iL].size(); ++iC2 ) {
        if( fCells[iL-1][iC1].Points[1]==fCells[iL][iC2].Points[0] && 
            fCells[iL-1][iC1].Points[2]==fCells[iL][iC2].Points[1] && 
            fCells[iL-1][iC1].Level >= fCells[iL][iC2].Level - 1 ) {
          // The implementation of the curvature based matching has to be studied. 
          fCells[iL][iC2].Level = fCells[iL-1][iC1].Level+1;
          fCells[iL][iC2].Neighbours.push_back(iC1);
        }
      }
    }
  }
  
  for (int level = 5; level > 0; --level ) {
    vector<Road> roads; 
    roads.reserve(100); // should reserve() be based on number of clusters on outermost layer?
    for (int iCL=4; iCL >= level-1; --iCL ) {
      for (size_t iCell = 0; iCell < fCells[iCL].size(); ++iCell) {
        if ( fCells[iCL][iCell].Level != level ) continue;
        roads.push_back( Road(iCL,iCell) );
        for( size_t iN=0;iN<fCells[iCL][iCell].Neighbours.size(); ++iN ) {
          const int currD = iCL - 1;
          const int neigh = fCells[iCL][iCell].Neighbours[iN];
          if( iN > 0 ) roads.push_back(roads.back());
          CandidatesTreeTraversal(roads,neigh,currD);
        }
        fCells[iCL][iCell].Level = -1;
      }
    }

    // for(size_t iR=0; iR<roads.size(); ++iR) {
    //   cout << "ROAD " << iR << " | ";
    //   for(int i=0;i<5;++i) {
    //     if(roads[iR][i]<0) continue;
    //     else {
    //       if(roads[iR].Label==-1){
    //         roads[iR].Label = fCells[i][roads[iR][i]].Label;
    //         if(roads[iR].Label==-1) roads[iR].Label--;
    //       }
    //       if (fCells[i][roads[iR][i]].Label!=roads[iR].Label&&roads[iR].Label>-1) { 
    //         roads[iR].Label = -1;
    //         if(fCells[i][roads[iR][i]].Label==-1) roads[iR].Label--;
    //       }

    //       cout << fCells[i][roads[iR][i]].Label << " ";
    //     }
    //   }
    //   cout << " | " << roads[iR].Label << " | " << roads[iR].N << endl;
    // }
    vector<AliITSUTrackCooked> candidates;
    candidates.reserve(roads.size());

    for (size_t iR=0; iR<roads.size(); ++iR) { 
      if(roads[iR].N != level) {
        continue;
      }

      int points[2*AliITSUAux::kMaxLayers];
      for(unsigned int i=0;i<2*AliITSUAux::kMaxLayers;++i) points[i] = -1;
      for(int i=0;i<5;++i) {
        if(roads[iR].Elements[i]<0) continue;
        points[( i )<<1]=fLayer[ i ](fCells[i][roads[iR].Elements[i]].Points[0]);
        points[(i+1)<<1]=fLayer[i+1](fCells[i][roads[iR].Elements[i]].Points[1]);
        points[(i+2)<<1]=fLayer[i+2](fCells[i][roads[iR].Elements[i]].Points[2]);
      }

      candidates.push_back(AliITSUTrackCooked());
      
      InitTrackParams(candidates.back(),points);
      const double chi2 = RefitTrack( (AliExternalTrackParam*)&candidates.back(), points, 0. ,-1);

      if ( chi2 < 0. ) {
        // cout << "FAIL: " << chi2 << endl;
        // for(unsigned int i=0;i<2*AliITSUAux::kMaxLayers;++i) 
        //   cout << points[i] << " ";
        // cout << endl;
        candidates.back().SetChi2( 1e27 );
      } else candidates.back().SetChi2( chi2 );
      candidates.back().SetLabel(roads[iR].Label);
    }

    vector<int> index;
    index.reserve(candidates.size());
    for ( size_t i = 0; i < candidates.size(); ++i ) index.push_back(i);
    Comparison<AliITSUTrackCooked> comp(&candidates);
    sort(index.begin(),index.end(),comp);

    for ( size_t cand = 0; cand < candidates.size(); ++cand ) {
      const int ii = index[cand];

      if ( candidates[ii].GetChi2() < 0. ) continue;
      
      // cout << candidates[ii].GetChi2() << " " << candidates[ii].GetNumberOfClusters() << " | " << candidates[ii].GetLabel() << " | ";
      // for(int i=0;i<candidates[ii].GetNumberOfClusters();++i) {
      //   cout<< GetCluster(candidates[ii].GetClusterIndex(i))->GetLabel(0) << " ";
      // }
      // cout << endl;

      if( candidates[ii].GetChi2()/candidates[ii].GetNumberOfClusters() > fgkChi2Cut ) {      
        break;
      }
      bool goodTrack = true;
      for ( int point = 0; point < candidates[ii].GetNumberOfClusters(); ++point ) { 
        int layer = (candidates[ii].GetClusterIndex(point)&0xf0000000)>>28;
        int ind = (candidates[ii].GetClusterIndex(point)&0x0fffffff);

        if( (fLayer[ layer ].Points[ ind ].Used ) ) {
          goodTrack = false;
        }

      }
      if(!goodTrack) {
        continue;
      }
      for ( int point = 0; point < candidates[ii].GetNumberOfClusters(); ++point ) {
        int layer = (candidates[ii].GetClusterIndex(point)&0xf0000000)>>28;
        int ind = (candidates[ii].GetClusterIndex(point)&0x0fffffff);
        fLayer[ layer ].Points[ ind ].Used = true;
      }

      AliESDtrack outTrack;
      CookLabel((AliKalmanTrack*)&candidates[ii],0.f);
      ntrk++;
      if(candidates[ii].GetChi2()>0) ngood++;

      // cout << candidates[ii].GetChi2() << " " << candidates[ii].GetNumberOfClusters() << " | " << candidates[ii].GetLabel() << " | ";
      // for(int i=0;i<candidates[ii].GetNumberOfClusters();++i) {
      //   cout<< GetCluster(candidates[ii].GetClusterIndex(i))->GetLabel(0) << " ";
      // }
      // cout << endl;

      outTrack.UpdateTrackParams((AliKalmanTrack*)&candidates[ii],AliESDtrack::kITSin);
      outTrack.SetLabel(candidates[ii].GetLabel());
      event->AddTrack(&outTrack);
    }
  }
  Info("Clusters2Tracks","Reconstructed tracks: %d",ntrk);
  if (ntrk)
    Info("Clusters2Tracks","Good tracks/reconstructed: %f",Float_t(ngood)/ntrk);
}


//________________________________________________________________________________
void AliITSUTrackerSA::CellsCreation(const int &cutLevel) {
  // Make associations between two points on adjacent layers within an azimuthal window.
  // Under consideration:
  // - track parameter estimation using the primary vertex position
  // To do:
  // - last iteration

  float phiCut = 7.f;
  if( cutLevel==0 ) phiCut = fPhiCut;

  // Doublets creation
  vector<Cell> doublets[6];
  for( int iL = 0 ; iL < 6 ; ++iL ) {
    for ( int iC1 = 0 ; iC1 < fLayer[iL].N ; ++iC1 ) {
      for ( int iC2 = 0; iC2 < fLayer[iL+1].N ; ++iC2 ) {
        const float dPhi = Abs( fLayer[iL][iC1].Phi - fLayer[iL+1][iC2].Phi );
        if( dPhi < phiCut || Abs( dPhi - TMath::TwoPi() ) < phiCut) {
          doublets[iL].push_back(Cell(iC1,iC2));
          if(Abs(fLayer[iL][iC1].XYZ[0]-fLayer[iL+1][iC2].XYZ[0])<1e-32) {
            doublets[iL].back().Param[0] = 1e32;
          } else {
            doublets[iL].back().Param[0] = (fLayer[iL][iC1].XYZ[1]-fLayer[iL+1][iC2].XYZ[1])/(fLayer[iL][iC1].XYZ[0]-fLayer[iL+1][iC2].XYZ[0]);
          }
          const float r1  = Sqrt(fLayer[iL][iC1].XYZ[0] * fLayer[iL][iC1].XYZ[0] + fLayer[iL][iC1].XYZ[1] * fLayer[iL][iC1].XYZ[1]);
          const float r2  = Sqrt(fLayer[iL+1][iC2].XYZ[0] * fLayer[iL+1][iC2].XYZ[0] + fLayer[iL+1][iC2].XYZ[1] * fLayer[iL+1][iC2].XYZ[1]);
          doublets[iL].back().Param[1] = (fLayer[iL][iC1].XYZ[2]-fLayer[iL+1][iC2].XYZ[2])/(r1-r2);
          doublets[iL].back().Label=-1;
          for(int i=0;i<3;++i) {
            for(int j=0;j<3;++j) {
              if(fLayer[iL][iC1].Label[i]>-1&&fLayer[iL][iC1].Label[i]==fLayer[iL+1][iC2].Label[j])
                doublets[iL].back().Label = fLayer[iL][iC1].Label[i];
            }
          } 
        } else if( fLayer[iL+1][iC2].Phi - fLayer[iL][iC1].Phi > phiCut ) break;
      }

    }
  }

  // Triplets creation
  for( int iL = 5; iL > 0; --iL ) {
    fCells[iL-1].clear();
    for ( size_t iD2 = 0; iD2 < doublets[iL].size(); ++iD2 ) {
      for ( size_t iD1 = 0; iD1 < doublets[iL-1].size(); ++iD1 ) {
        const int id1 = doublets[iL-1][iD1].Points[1];
        const int id2 = doublets[iL][iD2].Points[0];
        if ( id1 == id2 ) {
          const int id3 = doublets[iL][iD2].Points[1];
          const float r3 = Sqrt( fLayer[iL+1][id3].XYZ[0] * fLayer[iL+1][id3].XYZ[0] + fLayer[iL+1][id3].XYZ[1] * fLayer[iL+1][id3].XYZ[1] );
          const float r2 = Sqrt( fLayer[iL][id2].XYZ[0] * fLayer[iL][id2].XYZ[0] + fLayer[iL][id2].XYZ[1] * fLayer[iL][id2].XYZ[1] );
          const float extrZ3 = doublets[iL-1][iD1].Param[1] * ( r3 - r2 ) + fLayer[iL][id2].XYZ[2] ;
          const int iii = doublets[iL-1][iD1].Points[0];
          if ( Abs ( extrZ3 - fLayer[iL+1][id3].XYZ[2] ) < fZCut ) {      
            fCells[iL-1].push_back(Cell(doublets[iL-1][iD1].Points[0],id2,id3));
            fCells[iL-1].back().Param[0] = Curvature(fLayer[iL+1][id3].XYZ[0],fLayer[iL+1][id3].XYZ[1],fLayer[iL][id2].XYZ[0],fLayer[iL][id2].XYZ[1],fLayer[iL-1][iii].XYZ[0],fLayer[iL-1][iii].XYZ[1]);
            fCells[iL-1].back().Param[1] = doublets[iL][iD2].Param[1];
            if(doublets[iL-1][iD1].Label==doublets[iL][iD2].Label&&doublets[iL][iD2].Label!=-1) 
              fCells[iL-1].back().Label=doublets[iL][iD2].Label;
            else
              fCells[iL-1].back().Label=-1;
          } 
        } 
      }
    }
  }

}

//______________________________________________________________________________
Bool_t AliITSUTrackerSA::InitTrackParams(AliITSUTrackCooked &track, int points[])
{
  // Set the initial guess on track kinematics for propagation.
  // Assume at least 3 points available
  int lrOcc[AliITSUAux::kMaxLayers], nCl=0;
  //
  // we will need endpoints and middle layer
  for (int i=fITS->GetNLayersActive()-1; i>=0; i--) {
    if (points[i<<0x1]>-1) {
      lrOcc[nCl++] = i;
      track.SetClusterIndex(i,points[i<<0x1]);
    }
  }

  if (nCl<3) {
    AliError(Form("Cannot estimate momentum of tracks with %d clusters",nCl));
    return kFALSE;
  }
  //
  const int lr0   = lrOcc[0];
  const int lr1   = lrOcc[nCl/2];
  const int lr2   = lrOcc[nCl-1];
  //
  const SpacePoint& cl0 = fLayer[lr0].Points[ points[lr0<<1] ];
  const SpacePoint& cl1 = fLayer[lr1].Points[ points[lr1<<1] ];
  const SpacePoint& cl2 = fLayer[lr2].Points[ points[lr2<<1] ];
  double cv = Curvature(cl0.XYZ[0],cl0.XYZ[1], cl1.XYZ[0],cl1.XYZ[1], cl2.XYZ[0],cl2.XYZ[1]);

  double tgl = (cl2.XYZ[2]-cl0.XYZ[2])/TMath::Sqrt((cl2.XYZ[0]-cl0.XYZ[0])*(cl2.XYZ[0]-cl0.XYZ[0])+(cl2.XYZ[1]-cl0.XYZ[1])*(cl2.XYZ[1]-cl0.XYZ[1]));
  //
  AliITSUClusterPix* clus = (AliITSUClusterPix*)fClustersTC[ lr0 ]->At( points[lr0<<1] );
  AliITSURecoLayer* lr = fITS->GetLayerActive(lr0);
  AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
  double x = sens->GetXTF() + clus->GetX();
  double alp = sens->GetPhiTF();
  //  printf("Alp: %f phi: %f\n",alp,phi);
  double par[5] = {clus->GetY(),clus->GetZ(),0,tgl,cv};
  double cov[15] = {
    5*5,
    0,  5*5,
    0,  0  , 0.7*0.7,
    0,  0,   0,       0.7*0.7,
    0,  0,   0,       0,      10
  };
  track.Set(x,alp,par,cov);

  return kTRUE;
}

//______________________________________________________________________________
void AliITSUTrackerSA::CandidatesTreeTraversal(vector<Road> &candidates, const int &iD, const int &doubl) {

  if ( doubl < 0 ) return;

  candidates.back().AddElement(doubl,iD);
  const int currentN = candidates.back().N;
  for ( size_t iN = 0; iN < fCells[doubl][iD].Neighbours.size(); ++iN ) {
    const int currD = doubl - 1 ;
    const int neigh = fCells[doubl][iD].Neighbours[iN];
    
    if ( iN > 0 ) {
      candidates.push_back(static_cast<Road>(candidates.back()));
      candidates.back().N = currentN;
    }

    CandidatesTreeTraversal(candidates,neigh,currD);
  }
  
  fCells[doubl][iD].Level = -1;

}

//______________________________________________________________________________
Double_t AliITSUTrackerSA::RefitTrack(AliExternalTrackParam* trc, 
              Int_t clInfo[2*AliITSUAux::kMaxLayers],
              Double_t rDest, Int_t stopCond)
{
  // refit track till radius rDest. 
  // if stopCond<0 : propagate till last cluster then stop
  // if stopCond==0: propagate till last cluster then try to go till limiting rDest, don't mind if fail
  // if stopCond>0 : rDest must be reached
  //
  // The clList should provide the indices of clusters at corresponding layer (as stored in the layer
  // TClonesArray, with convention (allowing for up to 2 clusters per layer due to the overlaps):
  // if there is a cluster on given layer I, then it should be stored at clInfo[2*I-1]
  // if there is an additional cluster on this layer, it goes to clInfo[2*I]
  // -1 means no cluster
  //
  double rCurr = Sqrt(trc->GetX()*trc->GetX() + trc->GetY()*trc->GetY());
  int dir,lrStart,lrStop;
  //
  dir = rCurr<rDest ? 1 : -1;
  lrStart = fITS->FindFirstLayerID(rCurr,dir);
  lrStop  = fITS->FindLastLayerID(rDest,dir); // lr id before which we have to stop
  //
  if (lrStop<0 || lrStart<0) AliFatal(Form("Failed to find start(%d) or last(%d) layers. "
             "Track from %.3f to %.3f",lrStart,lrStop,rCurr,rDest));
  //
  int nCl = 0;
  for (int i=2*fITS->GetNLayersActive();i--;) {if (clInfo[i]<0) continue; nCl++;}
  //
  AliExternalTrackParam tmpTr(*trc);
  double chi2 = 0;
  int iclLr[2],nclLr;
  int nclFit = 0;
  //
  int lrStop1 = lrStop+dir;
  for (int ilr=lrStart;ilr!=lrStop1;ilr+=dir) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if ( dir*(rCurr-lr->GetR(dir))>0) continue; // this layer is already passed
    int ilrA2,ilrA = lr->GetActiveID();
    // passive layer or active w/o hits will be traversed on the way to next cluster
    if (!lr->IsActive() || clInfo[ilrA2=(ilrA<<1)]<0) continue; 
    //
    // select the order in which possible 2 clusters (in case of the overlap) will be traversed and fitted
    nclLr=0;
    if (dir>0) { // clusters are stored in increasing radius order
      iclLr[nclLr++]=clInfo[ilrA2++];
      if (clInfo[ilrA2]>=0) iclLr[nclLr++]=clInfo[ilrA2];
    }
    else {
      if ( clInfo[ilrA2+1]>=0 ) iclLr[nclLr++]=clInfo[ilrA2+1];
      iclLr[nclLr++]=clInfo[ilrA2];
    }
    //
    Bool_t transportedToLayer = kFALSE;
    for (int icl=0;icl<nclLr;icl++) {
      AliITSUClusterPix* clus =  (AliITSUClusterPix*)lr->GetCluster(iclLr[icl]);
      AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
      if (!tmpTr.Rotate(sens->GetPhiTF())) return -1;
      //
      double xClus = sens->GetXTF()+clus->GetX();
      if (!transportedToLayer) {
  if (ilr!=lrStart && !TransportToLayerX(&tmpTr,lrStart,ilr,xClus)) return -1; // go to the entrance to the layer
  lrStart = ilr;
  transportedToLayer = kTRUE;
      }
      //
      if (!PropagateSeed(&tmpTr,xClus,fCurrMass)) return -1;
      //
      Double_t p[2]={clus->GetY(), clus->GetZ()};
      Double_t cov[3]={clus->GetSigmaY2(), clus->GetSigmaYZ(), clus->GetSigmaZ2()};
      double chi2cl = tmpTr.GetPredictedChi2(p,cov);
      chi2 += chi2cl;
      //
      if ( !tmpTr.Update(p,cov) ) return -1;
      if (++nclFit==nCl && stopCond<0) {
  *trc = tmpTr;
  return chi2; // it was requested to not propagate after last update
      }
    }
    //
  }
  // All clusters were succesfully fitted. Even if the track does not reach rDest, this is enough to validate it.
  // Still, try to go as close as possible to rDest.
  //
  if (lrStart!=lrStop) {
    if (!TransportToLayer(&tmpTr,lrStart,lrStop)) return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
    if (!GoToExitFromLayer(&tmpTr,fITS->GetLayer(lrStop),dir)) return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
  }
  // go to the destination radius. Note that here we don't select direction to avoid precision problems
  if (!tmpTr.GetXatLabR(rDest,rDest,GetBz(),0) || !PropagateSeed(&tmpTr,rDest,fCurrMass, 100, kFALSE)) {
    return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
  }
  *trc = tmpTr;

  return chi2;
}

//______________________________________________________________________________
Bool_t AliITSUTrackerSA::PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass, Double_t maxStep, Bool_t matCorr) 
{
  // propagate seed to given x applying material correction if requested
  const Double_t kEpsilon = 1e-5;
  Double_t xpos     = seed->GetX();
  Int_t dir         = (xpos<xToGo) ? 1:-1;
  Double_t xyz0[3],xyz1[3];
  //
  Bool_t updTime = dir>0 && seed->IsStartedTimeIntegral();
  if (matCorr || updTime) seed->GetXYZ(xyz1);   //starting global position
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t bz=GetBz();   // getting the local Bz
    if (!seed->PropagateTo(x,bz))  return kFALSE;
    double ds = 0;
    if (matCorr || updTime) {
      xyz0[0]=xyz1[0]; // global pos at the beginning of step
      xyz0[1]=xyz1[1];
      xyz0[2]=xyz1[2];
      seed->GetXYZ(xyz1);    //  // global pos at the end of step
      //
      if (matCorr) {
  Double_t xrho,xx0;
  ds = GetMaterialBudget(xyz0,xyz1,xx0,xrho); 
  if (dir>0) xrho = -xrho; // outward should be negative
  if (!seed->CorrectForMeanMaterial(xx0,xrho,mass)) return kFALSE;
      }
      else { // matCorr is not requested but time integral is
  double d0 = xyz1[0]-xyz0[0];
  double d1 = xyz1[1]-xyz0[1];
  double d2 = xyz1[2]-xyz0[2];  
  ds = TMath::Sqrt(d0*d0+d1*d1+d2*d2);
      }
    }
    if (updTime) seed->AddTimeStep(ds);
    //
    xpos = seed->GetX();
  }
  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliITSUTrackerSA::TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t rLim)
{
  // transport track from layerFrom to the entrance of layerTo or to rLim (if>0), wathever is closer
  //  
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  Bool_t limReached = kFALSE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) return kFALSE; // go till the end of current layer
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    if (rLim>0) {
      if (dir>0) {
  if (rLim<xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
      else {
  if (rLim>xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
    }
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
    if (limReached) break;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerSA::TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t xStop)
{
  // transport track from layerFrom to the entrance of layerTo but do not pass control parameter X
  //  
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) return kFALSE; // go till the end of current layer
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir); // R of the entrance to layer
    //
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    if ( (dir>0&&xToGo>xStop) || (dir<0&&xToGo<xStop) ) xToGo = xStop;
    //
#ifdef _ITSU_DEBUG_
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
#endif
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerSA::GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the exit from lr in direction dir, applying material corrections in steps specific for this layer
  // If check is requested, do this only provided the track has not exited the layer already
  double xToGo = lr->GetR(dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // on the surface or outside of the layer
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // on the surface or outside of the layer
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, lr->GetMaxStep())) return kFALSE;
  //
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerSA::GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the entrance of lr in direction dir, w/o applying material corrections.
  // If check is requested, do this only provided the track did not reach the layer already
  double xToGo = lr->GetR(-dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // already passed
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, 100, kFALSE)) return kFALSE;
  return kTRUE;
  //
}

//____________________________________________________
Double_t AliITSUTrackerSA::GetMaterialBudget(const double* pnt0,const double* pnt1, double& x2x0, double& rhol) const
{
  double par[7];
  if (fUseMatLUT && fMatLUT) {
    double d = fMatLUT->GetMatBudget(pnt0,pnt1,par);
    x2x0 = par[AliITSUMatLUT::kParX2X0];
    rhol = par[AliITSUMatLUT::kParRhoL];    
    return d;
  }
  else {
    MeanMaterialBudget(pnt0,pnt1,par);
    x2x0 = par[1];
    rhol = par[0]*par[4];    
    return par[4];
  }
}

//____________________________________________________________________
Double_t AliITSUTrackerSA::Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3) {

  //calculates the curvature of track
  Double_t den = (x3-x1)*(y2-y1)-(x2-x1)*(y3-y1);
  if(den*den<1e-32) return 0.;
  Double_t a = ((y3-y1)*(x2*x2+y2*y2-x1*x1-y1*y1)-(y2-y1)*(x3*x3+y3*y3-x1*x1-y1*y1))/den;
  if((y2-y1)*(y2-y1)<1e-32) return 0.;
  Double_t b = -(x2*x2-x1*x1+y2*y2-y1*y1+a*(x2-x1))/(y2-y1);
  Double_t c = -x1*x1-y1*y1-a*x1-b*y1;
  Double_t xc= -a/2.;

  if((a*a+b*b-4*c)<0) return 0.;
  Double_t rad = TMath::Sqrt(a*a+b*b-4*c)/2.;
  if(rad*rad<1e-32) return 1e16;

  if((x1>0 && y1>0 && x1<xc)) rad*=-1;
  if((x1<0 && y1>0 && x1<xc)) rad*=-1;
    //  if((x1<0 && y1<0 && x1<xc)) rad*=-1;
    // if((x1>0 && y1<0 && x1<xc)) rad*=-1;

  return 1/rad;

}

