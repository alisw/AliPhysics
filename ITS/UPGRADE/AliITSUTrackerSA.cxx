//-------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    It reads AliITSUClusterPix clusters and and fills the ESD with tracks
//-------------------------------------------------------------------------

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
const Double_t AliITSUTrackerSA::fgkChi2Cut =  10.f;

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA(AliITSUReconstructor* rec) :
fReconstructor(rec),
fITS(0),
fMatLUT(0),
fUseMatLUT(kFALSE),
fCurrMass(0.14),
//
fClusters(),
fClustersTC(),
fDoublets(),
fIndex(),
fNClusters(),
fNDoublets(),
fPhiCut(0.05),
fRPhiCut(0.03),
fZCut(0.01)
#ifdef __DEBUG__
,fCv(0x0)
,fMk(0x0)
,fLn(0x0)
,fTx(0x0)
#endif
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------
  for(Int_t i=0;i<7;++i) {
    fClusters[i].reserve(5000);
  }
  if (rec) Init(rec);
}

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA(const AliITSUTrackerSA &t):
AliTracker(t),
fReconstructor(t.fReconstructor),
fITS(t.fITS),
fMatLUT(t.fMatLUT),
fUseMatLUT(t.fUseMatLUT),
fCurrMass(t.fCurrMass),
//
fClusters(),
fClustersTC(),
fIndex(),
fNClusters(),
fNDoublets(),
fPhiCut(),
fRPhiCut(),
fZCut()
#ifdef __DEBUG__
,fCv(0x0)
,fMk(0x0)
,fLn(0x0)
,fTx(0x0)
#endif
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

  MakeDoublets();       // To be checked
  //MakeTriplets();       // Are triplets really necessary? MFT does not use them.
  CASelection(event);

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::PropagateBack(AliESDEvent * /*event*/) {
  //--------------------------------------------------------------------
  // Here, we implement the Kalman smoother ?
  // The clusters must be already loaded
  //--------------------------------------------------------------------

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::RefitInward(AliESDEvent * /*event*/) {
  //--------------------------------------------------------------------
  // Some final refit, after the outliers get removed by the smoother ?
  // The clusters must be loaded
  //--------------------------------------------------------------------

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
    TClonesArray *clCont=fClustersTC[iL];
    fNClusters[iL]=clCont->GetEntriesFast();
    Float_t phi[fNClusters[iL]];
    fIndex[iL] = new Int_t[fNClusters[iL]];
    AliITSURecoLayer* lr = fITS->GetLayerActive(iL) ; // assign the layer which the cluster belongs to
    for(int iC=0;iC<fNClusters[iL];++iC) {
      const AliITSUClusterPix *cl = (AliITSUClusterPix*)clCont->At(iC);
      float pos[3];
      cl->GetGlobalXYZ(pos);
      phi[iC] = pos[0]==0.f ? TMath::PiOver2() : TMath::ATan2(pos[1]-GetY(),pos[0]-GetX());
      AliITSURecoSens* sens = lr->GetSensorFromID(cl->GetVolumeId());
      //double x = sens->GetXTF() + clus->GetX();
      float angle= sens->GetPhiTF();
      #ifdef __DEBUG__
      int label=cl->GetLabel(0);
      //cout << "Guarda te che label" << label << endl;
      fClusters[iL].push_back(itsCluster(pos[0],pos[1],pos[2],cl->GetSigmaY2(),cl->GetSigmaZ2(),cl->GetSigmaYZ(),phi[iC],angle,label));
      #else
      fClusters[iL].push_back(itsCluster(pos[0],pos[1],pos[2],cl->GetSigmaY2(),cl->GetSigmaZ2(),cl->GetSigmaYZ(),phi[iC],angle));
      #endif
    }
    TMath::Sort(fNClusters[iL],phi,fIndex[iL],kFALSE);
  }
  #ifdef __DEBUG__
  //PrintInfo("clusters");
  DrawEvent("clusters");
  #endif
  return 0;
}

//________________________________________________________________________________
void AliITSUTrackerSA::UnloadClusters() {
  //--------------------------------------------------------------------
  // This function unloads ITSU clusters from the RAM
  //--------------------------------------------------------------------
  for(int i=0;i<7;++i) {
    fClusters[i].clear();
    fNClusters[i]=0;
    delete fIndex[i];
  }
  for(int i=0;i<6;++i) fDoublets[i].clear();
}

//________________________________________________________________________________
AliCluster *AliITSUTrackerSA::GetCluster(Int_t /*index*/) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  return 0;  // replace with an actual pointer
}

//________________________________________________________________________________
void AliITSUTrackerSA::CASelection(AliESDEvent *event) {
  // Here it's implemented the Cellular Automaton routine
  // Firstly the level of each doublet is set according to the level of
  // the neighbour doublets.
  // Doublet are considered to be neighbour if they share one point and the
  // phi and theta direction difference of the two is below a cut value.

  cout << "Begin of the CA selection" << endl;
  for( int iL = 1; iL < 6; ++iL ) {

    const itsCluster* clusters1 = &fClusters[iL-1][0];
    const itsCluster* clusters2 = &fClusters[iL][0];
    const itsCluster* clusters3 = &fClusters[iL+1][0];

    const nPlets* doublets1 = &fDoublets[iL-1][0];
    nPlets* doublets2 = &fDoublets[iL][0];

    for ( int iD2 = 0; iD2 < fNDoublets[iL]; ++iD2 ) {
      for ( int iD1 = 0; iD1 < fNDoublets[iL-1]; ++iD1 ) {
        const int id1 = doublets1[iD1].id1;
        const int id2 = doublets2[iD2].id0;
        if ( id1 == id2 ) {
          if ( doublets2[iD2].level <= ( doublets1[iD1].level + 1 ) ) {
            const int id3 = doublets2[iD2].id1;
            const float r3 = Sqrt( clusters3[id3].x * clusters3[id3].x + clusters3[id3].y * clusters3[id3].y );
            const float r2 = Sqrt( clusters2[id2].x * clusters2[id2].x + clusters2[id2].y * clusters2[id2].y );
            const float extrZ3 = doublets1[iD1].tanLambda * ( r3 - r2 ) + clusters2[id2].z ;
            //cout << extrZ3 << " " << clusters3[id3].z << " " << Abs ( extrZ3 - clusters3[id3].z ) << endl;
            if ( Abs ( extrZ3 - clusters3[id3].z ) < fZCut ) {
              //cout << "OK Z doublets: "<< iL-1 << "," << iD1 << "\t" << iL << "," <<iD2 << endl;
              const float det = (clusters1[id1].x - GetX())*(clusters2[id2].y - GetY()) - (clusters2[id2].x-GetX() )*(clusters1[id1].y - GetY()); // (GetX() - clusters2[id2].x)*(clusters1[id1].y - clusters2[id2].y) - (GetY() - clusters2[id2].y)*(clusters1[id1].x - clusters2[id2].x);
              //cout << det << endl;
              if ( Abs(det) <= 1e-12 ) {
                // linear extrapolation to next layer
                const float dsq = ( doublets1[iD1].tanPhi * (clusters3[id3].x + clusters2[id2].x) + clusters3[id3].y - clusters2[id2].y ) *
                ( doublets1[iD1].tanPhi * (clusters3[id3].x + clusters2[id2].x) + clusters3[id3].y - clusters2[id2].y ) / (1 + doublets1[iD1].tanPhi * doublets1[iD1].tanPhi );
                if ( dsq < fRPhiCut*fRPhiCut )  {
                  doublets2[iD2].level = doublets1[iD1].level+1;
                  doublets2[iD2].neighbours.push_back(iD1);
                }
              } else {
                const float r1sq = clusters1[id1].x * clusters1[id1].x + clusters1[id1].y * clusters1[id1].y ;
                const float rvsq = GetX() * GetX() + GetY() * GetY();
                const float deta = (rvsq - r1sq) * (clusters2[id2].y - GetY()) - (rvsq - r2*r2) * (clusters1[id1].y - GetY());
                const float detb = - (rvsq - r1sq) * (clusters2[id2].x - GetX()) + (rvsq - r2*r2) * (clusters1[id1].x - GetX()) ;
                const float a = deta/det ;
                const float b = detb/det ;
                const float c = -rvsq - a * GetX() - b * GetY();
                const float rc = Sqrt( a*a/4.f + b*b/4.f - c );
                const float d = Sqrt( (a/2.f + clusters3[id3].x) * (a/2.f + clusters3[id3].x) + (b/2.f + clusters3[id3].y) * (b/2.f + clusters3[id3].y) );
                //cout << d << " " << rc << " " << d - rc << endl;
                if ( Abs( d - rc ) < fRPhiCut ) {
                  doublets2[iD2].level = doublets1[iD1].level+1;
                  doublets2[iD2].neighbours.push_back(iD1);
                }
              }
            }
          }
        }
      }
    }
  }
  //#ifdef __DEBUG__
  //PrintInfo("doublets");
  //DrawEvent("doublets+level");
  //return;
  //#endif
  // Hic sunt leones: the following code could be optimised to be iterative. But now I don't have time.
  vector<trackC> tracks;
  for ( int level = 6; level >= 2 ; --level ) {
    cout << "level " << level << endl;
    vector<Road> roads;
    roads.clear();
    roads.reserve(100);
    for ( int doubl = 5; doubl >= level-1; --doubl ) {
      for ( int iD = 0; iD < fNDoublets[doubl]; ++iD ) {
        if ( fDoublets[doubl][iD].level == level ) {
          roads.push_back(Road());
          roads.back().AddElement(doubl,iD);
          cout << "\nseed " << iD << "("<<fDoublets[doubl][iD].id1<<","<<fDoublets[doubl][iD].id0;
          cout <<") in array " << doubl << " with level " << level;
          for ( unsigned int iN = 0; iN < fDoublets[doubl][iD].neighbours.size(); ++iN ) {
            const int currD = doubl - 1 ;
            const int neigh = fDoublets[doubl][iD].neighbours[iN];
            //if ( level != fDoublets[currD][neigh].level + 1 ) continue;
            if ( iN > 0 ) roads.push_back(static_cast<Road>(roads.back()));
            CandidatesTreeTraversal(roads,neigh,currD);
            cout << endl;// << roads.back() << endl;
          }
          fDoublets[doubl][iD].level = -1; // mark as used
        }
        //for ( int j = 0; j < 2*AliITSUAux::kMaxLayers; ++j ) cout << roads.back().fPoints[j] << " "; 
        //cout << endl;
      }
    }

    
    Double_t rDest = 0.;//fITS->GetRMax();
    DrawRoads(roads);
 
    vector<trackC> candidates;
    for ( size_t iR = 0; iR < roads.size(); ++iR ) {
      if ( roads[iR].fNElements!=level ) {
        cout << "JUMP" << endl;
        continue;
      } 
      candidates.push_back(trackC());
      candidates.back().fNPoints = level+1;
      for ( size_t j = 0; j < 6; ++j ) { 
        if ( roads[iR].fElements[j] == -1 ) continue;
        candidates.back().fPoints[j<<0x1] = fDoublets[j][roads[iR].fElements[j]].id0;
        candidates.back().fPoints[(j+1)<<0x1] = fDoublets[j][roads[iR].fElements[j]].id1;
        //cout << (j<<0x1) << " " << fDoublets[j][roads[iR].fElements[j]].id0 << "\t" << ((j+1)<<0x1) << " " << fDoublets[j][roads[iR].fElements[j]].id1 << endl;
      }
      //cout << endl;
      cout << "Candidate " << candidates.size() << ", number of points: " << level+1 << endl;
      //cout << candidates.back() << endl;
      InitTrackParams(candidates.back());
      candidates.back().fChi2 = RefitTrack( (AliExternalTrackParam*)&candidates.back(), candidates.back().fPoints, rDest ,-1);
      //cout << "Fit cnd: " << cand << " " << candidates[cand] << endl;
    }
  
    int index[candidates.size()];
    for ( size_t i = 0; i < candidates.size(); ++i ) index[i]=i;
    CompDesc comp(&candidates);
    sort(index,index+candidates.size(),comp);
    
    for ( size_t cand = 0; cand < candidates.size(); ++cand ) {
      const int ii = index[cand];
      //#ifdef __DEBUG__
      //cout << ii << " " << candidates[ii] << endl;
      //#endif
      if ( candidates[ii].fChi2 < 0. ) break;
      bool goodTrack = true;
      for ( unsigned int point = 0; point < 14; ++point ) { //-> here it's necessary to use a temporary array for the used clusters.
        if ( candidates[ii].fPoints[point] != -1 ) {
          if( !(fClusters[ point/2 ][ candidates[ii].fPoints[point] ].isUsed ) ) {
            fClusters[ point/2 ][ candidates[ii].fPoints[point] ].isUsed = true;
          } else {
            goodTrack = false;
          }
        }
      }
      //cout << endl;
      if ( goodTrack ) {
        tracks.push_back(candidates[ii]);
      }
    }
    //cout << "End of level " << level << endl;*/
  }

  bool joined[tracks.size()];
  MergeTracks(tracks,joined); 
  for ( unsigned int ii = 0; ii < tracks.size(); ++ii ) {
    if ( ! joined[ii] ) {
      AliESDtrack outTrack;
      outTrack.SetOuterParam((AliExternalTrackParam*)&tracks[ii],AliESDtrack::kITSpureSA);  
      event->AddTrack(&outTrack);
    }
  }
}

//________________________________________________________________________________
void AliITSUTrackerSA::MakeDoublets() {
  // Make associations between two points on adjacent layers within an azimuthal window.
  // Under consideration:
  // - track parameter estimation using the primary vertex position
  // To do:
  // - last iteration

  //cout << "Vertex of used by the tracker: " << GetX() << " " << GetY() << " " << GetZ() << endl;

  for( int iL = 0 ; iL < 6 ; ++iL ) {
    fNDoublets[iL] = 0;
    const itsCluster* clusters1 = &fClusters[iL][0];
    const itsCluster* clusters2 = &fClusters[iL+1][0];

    // 0 - 2Pi junction treatment (part I)
    for ( int iCC1 = 0 ; iCC1 < fNClusters[iL] ; ++iCC1 ) {
      bool flag = true;
      const int iC1 = fIndex[iL][iCC1];
      for ( int iCC2 = fNClusters[iL+1]-1; iCC2 >= 0 ; --iCC2 ) {
        const int iC2 = fIndex[iL+1][iCC2];
        if( (TMath::TwoPi() - (clusters2[iC2].phi-clusters1[iC1].phi) ) < fPhiCut ) {
          #ifdef __DEBUG__
          fDoublets[iL].push_back(nPlets(iC1,iC2,clusters1[iC1].pid,clusters2[iC2].pid));
          #else
          fDoublets[iL].push_back(nPlets(iC1,iC2));
          #endif
          fDoublets[iL][fNDoublets[iL]].tanPhi = (clusters1[iC1].y-clusters2[iC2].y)/(clusters1[iC1].x-clusters2[iC2].x);
          float r1  = Sqrt(clusters1[iC1].x * clusters1[iC1].x + clusters1[iC1].y * clusters1[iC1].y);
          //cout << clusters2[iC2].x * clusters2[iC2].x + clusters2[iC2].y * clusters2[iC2].y << flush << endl;
          float r2  = Sqrt(clusters2[iC2].x * clusters2[iC2].x + clusters2[iC2].y * clusters2[iC2].y);
          fDoublets[iL][fNDoublets[iL]].tanLambda = (clusters1[iC1].z-clusters2[iC2].z)/(r1-r2);
          ++fNDoublets[iL];
          flag = false;
        } else break;

      }
      if (flag) break;
    }


    // "Central" points
    for ( int iCC1 = 0 ; iCC1 < fNClusters[iL] ; ++iCC1 ) {
      const int iC1 = fIndex[iL][iCC1];
      for ( int iCC2 = 0; iCC2 < fNClusters[iL+1] ; ++iCC2 ) {
        const int iC2 = fIndex[iL+1][iCC2];
        if( Abs( clusters1[iC1].phi - clusters2[iC2].phi ) < fPhiCut ) {
          #ifdef __DEBUG__
          fDoublets[iL].push_back(nPlets(iC1,iC2,clusters1[iC1].pid,clusters2[iC2].pid));
          #else
          fDoublets[iL].push_back(nPlets(iC1,iC2));
          #endif
          fDoublets[iL][fNDoublets[iL]].tanPhi = (clusters1[iC1].y-clusters2[iC2].y)/(clusters1[iC1].x-clusters2[iC2].x);
          float r1  = Sqrt(clusters1[iC1].x * clusters1[iC1].x + clusters1[iC1].y * clusters1[iC1].y);
          float r2  = Sqrt(clusters2[iC2].x * clusters2[iC2].x + clusters2[iC2].y * clusters2[iC2].y);
          fDoublets[iL][fNDoublets[iL]].tanLambda = (clusters1[iC1].z-clusters2[iC2].z)/(r1-r2);
          ++fNDoublets[iL];
        } else if( clusters2[iC2].phi - clusters1[iC1].phi > fPhiCut ) break;

      }

    }

    // 0 - 2Pi junction treatment (part II)
    for ( int iCC1 = fNClusters[iL]-1; iCC1 > -1 ; --iCC1 ) {
      bool flag = true;
      const int iC1 = fIndex[iL][iCC1];
      for ( int iCC2 = 0; iCC2 < fNClusters[iL+1] ; ++iCC2 ) {
        const int iC2 = fIndex[iL+1][iCC2];
        if( (TMath::TwoPi() - (clusters1[iC1].phi-clusters2[iC2].phi) ) < fPhiCut ) {
          #ifdef __DEBUG__
          fDoublets[iL].push_back(nPlets(iC1,iC2,clusters1[iC1].pid,clusters2[iC2].pid));
          #else
          fDoublets[iL].push_back(nPlets(iC1,iC2));
          #endif
          fDoublets[iL][fNDoublets[iL]].tanPhi = (clusters1[iC1].y-clusters2[iC2].y)/(clusters1[iC1].x-clusters2[iC2].x);
          float r1  = Sqrt(clusters1[iC1].x * clusters1[iC1].x + clusters1[iC1].y * clusters1[iC1].y);
          float r2  = Sqrt(clusters2[iC2].x * clusters2[iC2].x + clusters2[iC2].y * clusters2[iC2].y);
          fDoublets[iL][fNDoublets[iL]].tanLambda = (clusters1[iC1].z-clusters2[iC2].z)/(r1-r2);
          ++fNDoublets[iL];
          flag = false;
        } else break;

      }

      if (flag) break;
    }

  }
  // #ifdef __DEBUG__
  // PrintInfo("doublets");
  // #endif
}
//______________________________________________________________________________
Bool_t AliITSUTrackerSA::InitTrackParams(trackC &track)
{
  // Set the initial guess on track kinematics for propagation.
  // Assume at least 3 points available
  int lrOcc[AliITSUAux::kMaxLayers], nCl=0;
  //
  // we will need endpoints and middle layer
  for (int i=fITS->GetNLayersActive();i--;) {
    if (track.fPoints[i<<0x1]>-1) {
      lrOcc[nCl++] = i;
      track.fInnermostLayer = i;
    }
  }
  track.fNPoints = nCl;
  track.fOutermostLayer = track.fInnermostLayer + nCl - 1; 
  if (nCl<3) {
    AliError(Form("Cannot estimate momentum of tracks with %d clusters",nCl));
    //cout << track << endl;
    return kFALSE;
  }
  //
  int lr0   = lrOcc[0];
  int lr1   = lrOcc[nCl/2];
  int lr2   = lrOcc[nCl-1];
  //
  //cout << (lr0<<0x1) << " " <<  (lr1<<0x1)<< " " <<  (lr2<<0x1) << endl;
  const itsCluster& cl0 = fClusters[lr0][ track.fPoints[lr0<<0x1] ];
  const itsCluster& cl1 = fClusters[lr1][ track.fPoints[lr1<<0x1] ];
  const itsCluster& cl2 = fClusters[lr2][ track.fPoints[lr2<<0x1] ];
  double cv = Curvature(cl0.x,cl0.y, cl1.x,cl1.y, cl2.x,cl2.y);
  double tgl = (cl2.z-cl0.z)/TMath::Sqrt((cl2.x-cl0.x)*(cl2.x-cl0.x)+(cl2.y-cl0.y)*(cl2.y-cl0.y));
  //  double phi = TMath::ATan2((cl2.y-cl1.y),(cl2.x-cl1.x));
  //
  AliITSUClusterPix* clus = (AliITSUClusterPix*)fClustersTC[ lr0 ]->At( track.fPoints[lr0<<0x1] );
  AliITSURecoLayer* lr = fITS->GetLayerActive(lr0);
  AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
  double x = sens->GetXTF() + clus->GetX();
  double alp = sens->GetPhiTF();
  //  printf("Alp: %f phi: %f\n",alp,phi);
  double par[5] = {clus->GetY(),clus->GetZ(),0,tgl,cv};
  double cov[15] = {
    5*5,
    0, 5*5,
    0, 0, 0.7*0.7,
    0,0,0,0.7*0.7,
    0,0,0,0,10
  };
  track.Set(x,alp,par,cov);
  //cout << track ;
  return kTRUE;
}

//______________________________________________________________________________
void AliITSUTrackerSA::CandidatesTreeTraversal(vector<Road> &candidates, const int &iD, const int &doubl) {

  if ( doubl < 0     ) {
    #ifdef __DEBUG__
    //cout << "ERROR IN CandidatesTreeTraversal" << endl;
    //cout << endl;
    #endif
    return;
  }
  candidates.back().AddElement(doubl,iD);
  cout << "\n\tTraversing the tree through " << iD << " in array " << doubl << " at level " << fDoublets[doubl][iD].level;
  //cout << doubl*2 << " ";
  for ( unsigned int iN = 0; iN < fDoublets[doubl][iD].neighbours.size(); ++iN ) {
    const int currD = doubl - 1 ;
    const int neigh = fDoublets[doubl][iD].neighbours[iN];
    //const int level = fDoublets[doubl][iD].level;
    //if ( level != fDoublets[currD][neigh].level + 1 ) continue;
    if ( iN > 0 ) {
      cout << " -> branching here!";
      candidates.push_back(static_cast<Road>(candidates.back()));
    }
    CandidatesTreeTraversal(candidates,neigh,currD);
  }
  
  fDoublets[doubl][iD].level = -1;

}

//______________________________________________________________________________
void AliITSUTrackerSA::MergeTracks(vector<trackC> &tracks, bool flag[]) {
  cout << "Merging tracks" << endl;
  for ( unsigned int iT = 0; iT < tracks.size(); ++iT ) { 
    flag[iT]=false; 
  }

  for ( unsigned int iT1 = 0; iT1 < tracks.size(); ++iT1 ) {
    if ( tracks[iT1].fNPoints > 4 || flag[iT1] ) continue;
    const int inL1 = tracks[iT1].fInnermostLayer;
    const int outL1 = tracks[iT1].fOutermostLayer;
    for ( unsigned int iT2 = iT1+1; iT2 < tracks.size(); ++iT2 ) {
      const int inL2 = tracks[iT2].fInnermostLayer;
      const int outL2 = tracks[iT2].fOutermostLayer;
      printf("%d: In/Out %d/%d\t%d: In/Out %d/%d",iT1,inL1,outL1,iT2,inL2,outL2);
      if ( outL1 < inL2 ) {
        cout << " <- check1" << flush;
        if ( TransportToLayer((AliExternalTrackParam*)&tracks[iT2],inL2,outL1) ) {
          cout << " <- transport";
          if ( tracks[iT2].Rotate( fClusters[outL1][ tracks[iT1].fPoints[ outL1<<0x1 ] ].phiM ) ) {
            const double chi2 = tracks[iT2].GetPredictedChi2((AliExternalTrackParam*)&tracks[iT1]);
            printf(" Merging candidates (%i,%i), chi2 = %f ",iT1,iT2,chi2);
            if ( chi2 < fgkChi2Cut ) {
              flag[iT2] = true;
              for ( int np = tracks[iT2].fInnermostLayer; np < tracks[iT2].fOutermostLayer; ++np ) {
                tracks[iT1].fPoints[np<<0x1] = tracks[iT2].fPoints[np<<0x1];
              }
              InitTrackParams(tracks[iT1]);
              RefitTrack((AliExternalTrackParam*)&tracks[iT1],tracks[iT1].fPoints,0,-1);
            } 
          } else { cout << " <- Failed rotation!!"; }
        } else { cout << " <- Failed transport!!"; }
      } else if ( inL1 > outL2 ) {
        cout << " <- check2" << flush;
        if ( TransportToLayer((AliExternalTrackParam*)&tracks[iT1],inL1,outL2) ) {
          cout << " <- transport";
          if ( tracks[iT1].Rotate( fClusters[outL2][ tracks[iT1].fPoints[ outL2<<0x1 ] ].phiM ) ) {
            const double chi2 = tracks[iT2].GetPredictedChi2((AliExternalTrackParam*)&tracks[iT1]);
            printf(" Merging candidates (%i,%i), chi2 = %f ",iT1,iT2,chi2);
            if ( chi2 < fgkChi2Cut ) {
              flag[iT2] = true;
              for ( int np = tracks[iT2].fInnermostLayer; np < tracks[iT2].fOutermostLayer; ++np ) {
                tracks[iT1].fPoints[np<<0x1] = tracks[iT2].fPoints[np<<0x1];
              }
              InitTrackParams(tracks[iT1]);
              RefitTrack((AliExternalTrackParam*)&tracks[iT1],tracks[iT1].fPoints,0,-1);
            }
          } else { cout << " <- Failed rotation!!"; } 
        } else { cout << " <- Failed transport!!"; }
      }
      cout << endl;
    }
  }
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
  //cout << "Start/End layers and direction " << lrStart << " " << lrStop << " " <<  dir << endl;
  //
  if (lrStop<0 || lrStart<0) AliFatal(Form("Failed to find start(%d) or last(%d) layers. "
  "Track from %.3f to %.3f",lrStart,lrStop,rCurr,rDest));
  //
  int nCl = 0;
  cout << "pid : ";
  for (int i=2*fITS->GetNLayersActive();i--;) {
    if (clInfo[i]<0) continue;
    cout << fClusters[i/2][clInfo[i]].pid;
    ++nCl;
  }
  cout << endl;
  //
  //cout << "#ptr: " << nCl << endl;
  AliExternalTrackParam tmpTr(*trc);
  double chi2 = 0;
  int iclLr[2],nclLr;
  int nclFit = 0;
  //
  int lrStop1 = lrStop+dir;
  for (int ilr=lrStart;ilr!=lrStop1;ilr+=dir) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if ( dir*(rCurr-lr->GetR(dir))>0) {
      cout << ilr << " passed!" << endl;
      continue;
    } // this layer is already passed
    int ilrA2,ilrA = lr->GetActiveID();
    // passive layer or active w/o hits will be traversed on the way to next cluster
    if (!lr->IsActive() || clInfo[ilrA2=(ilrA<<1)]<0) {
      //cout << ilr << " is inactive or without cluster for current candidates" << endl;
      continue;
    }
    //cout << "OK layer " << ilr << endl;
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
      if (!tmpTr.Rotate(sens->GetPhiTF())) { cout << "failed rotation" << endl; return -1; }
      //
      double xClus = sens->GetXTF()+clus->GetX();
      if (!transportedToLayer) {
        if (ilr!=lrStart && !TransportToLayerX(&tmpTr,lrStart,ilr,xClus)) {
          cout << "failed transport to the entrance" << endl; return -1; } // go to the entrance to the layer
          lrStart = ilr;
          transportedToLayer = kTRUE;
        }
        //
        if (!PropagateSeed(&tmpTr,xClus,fCurrMass)) { cout << "failed propagation of the seed X:" << xClus << endl; tmpTr.Print(); return -1; }
        //
        Double_t p[2]={clus->GetY(), clus->GetZ()};
        Double_t cov[3]={clus->GetSigmaY2(), clus->GetSigmaYZ(), clus->GetSigmaZ2()};
        double chi2cl = tmpTr.GetPredictedChi2(p,cov);
        chi2 += chi2cl;
        //
        if ( !tmpTr.Update(p,cov) ) {
          cout << "failed update of the covariance" << endl;
          return -1;
        }

        if (++nclFit==nCl && stopCond<0) {
          *trc = tmpTr;
          printf("Fit chi2: %f for %d clusters\n",chi2,nclFit);
          return chi2; // it was requested to not propagate after last update
        }
      }
      //
    }
    // All clusters were succesfully fitted. Even if the track does not reach rDest, this is enough to validate it.
    // Still, try to go as close as possible to rDest.
    //
    //  printf("Fit chi2: %f for %d clusters\n",chi2,nclFit);
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
      if (!seed->PropagateTo(x,bz)) {
        //cout << " failed PropagateTo " << endl;
        return kFALSE;
      }
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
    //printf("From %d, to %d, direction %d ",lFrom,lTo,dir);
    AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
    Bool_t checkFirst = kTRUE;
    Bool_t limReached = kFALSE;
    while(lFrom!=lTo) {
      if (lrFr) {
        if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) {
          cout << " Failed GoToExitFromLayer "; 
          return kFALSE; // go till the end of current layer
        }
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
        cout << " Failed GetXatLabR " << endl;
        return kFALSE;
      }
      if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
        //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
        //seed->Print("etp");
        cout << " Failed PropagateSeed " << endl;
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
  Double_t AliITSUTrackerSA::Curvature(Double_t x1,Double_t y1,Double_t
  x2,Double_t y2,Double_t x3,Double_t y3)
  {

    //calculates the curvature of track
    Double_t den = (x3-x1)*(y2-y1)-(x2-x1)*(y3-y1);
    if(den==0) return 0;
    Double_t a = ((y3-y1)*(x2*x2+y2*y2-x1*x1-y1*y1)-(y2-y1)*(x3*x3+y3*y3-x1*x1-y1*y1))/den;
    Double_t b = -(x2*x2-x1*x1+y2*y2-y1*y1+a*(x2-x1))/(y2-y1);
    Double_t c = -x1*x1-y1*y1-a*x1-b*y1;
    Double_t xc=-a/2.;

    if((a*a+b*b-4*c)<0) return 0;
    Double_t rad = TMath::Sqrt(a*a+b*b-4*c)/2.;
    if(rad==0) return 0;

    if((x1>0 && y1>0 && x1<xc)) rad*=-1;
    if((x1<0 && y1>0 && x1<xc)) rad*=-1;
    //  if((x1<0 && y1<0 && x1<xc)) rad*=-1;
    // if((x1>0 && y1<0 && x1<xc)) rad*=-1;

    return 1/rad;

  }

  #ifdef __DEBUG__
  //____________________________________________________
  void AliITSUTrackerSA::PrintInfo(TString what) {
    //
    if( what.Contains("clusters") ) {
      cout << "Dumping clusters info" << endl;
      for ( int i = 0; i < 7; ++i ) {
        cout << "**** Layer " << i << " ****" << endl;
        for ( int c = 0; c < fNClusters[i]; ++c ) {
          cout << "*** Cluster " << c << " ***" <<endl;
          cout << fClusters[i][fIndex[i][c]] << endl;
        }
      }
    }
    //
    if( what.Contains("doublets") ) {
      cout << "Dumping doublets info" << endl;
      for ( int i = 0; i < 6; ++i ) {
        cout << "**** Doublets array " << i << " ****" << endl;
        for ( int c = 0; c < fNDoublets[i]; ++c ) {
          cout << "*** Doublet " << c << " ***" <<endl;
          cout << fDoublets[i][c] << endl;
        }
      }
    }
  }

  //____________________________________________________
  void AliITSUTrackerSA::DrawEvent(TString what) {
    //

    const int size = 900;
    if (fCv == 0x0 ) {
      fCv = new TCanvas("cv_event","cv_event",size,size);
      fCv->Range(0,0,size,size);
    }
    if (fTx == 0x0 ) {
      fTx = new TText();
      fTx->SetTextSize(0.015);
      fTx->SetTextColor(kBlack);
      fTx->DrawText(0,0,"Bending plane");
    }

    const float cX = size/2.f; //pt.GetX();
    const float cY = size/2.f; //pt.GetY();
    const float sc = 8.5f;
    if( what.Contains("clusters") ) {
      if(fMk==0x0) fMk=new TMarker();
      fMk->SetMarkerStyle(20);
      fMk->SetMarkerColor(kRed);
      fMk->SetMarkerSize(0.3);
      for(int iL=0; iL<7; ++iL) {
        for(int iC=0; iC<fNClusters[iL];++iC) {
          fMk->DrawMarker(cX+fClusters[iL][iC].x*sc,cY+fClusters[iL][iC].y*sc);
          if( what.Contains("clusters+id") )  fTx->DrawText(cX+fClusters[iL][iC].x*sc,cY+fClusters[iL][iC].y*sc,Form("%i",fIndex[iL][iC])); 
        }
      }
    }
    //
    if( what.Contains("doublets") ) {
      if(fLn==0x0) fLn=new TLine();
      for(int iL=0; iL<6; ++iL) {
        for(int iD=0; iD<fNDoublets[iL]; ++iD) {
          const int id0 = fDoublets[iL][iD].id0;//fIndex[iL][fDoublets[iL][iD].id0];
          const int id1 = fDoublets[iL][iD].id1;//fIndex[iL+1][fDoublets[iL][iD].id1];
          fLn->DrawLine(cX+fClusters[iL][id0].x*sc,cY+fClusters[iL][id0].y*sc,cX+fClusters[iL+1][id1].x*sc,cY+fClusters[iL+1][id1].y*sc);
          if( what.Contains("doublets+level") ) fTx->DrawText(cX+(fClusters[iL][id0].x+fClusters[iL+1][id1].x)*sc/(2.f),cY+(fClusters[iL][id0].y+fClusters[iL+1][id1].y)*sc/(2.f),Form("%i",fDoublets[iL][iD].level));
        }
      }
    }
    //
  }

//____________________________________________________
void AliITSUTrackerSA::DrawRoads(vector<Road> &vec) {
  const int size = 900;
  const float cX = size/2.f; 
  const float cY = size/2.f; 
  const float sc = 8.5f;
  
  if (fCv == 0x0 ) {
    fCv = new TCanvas("cv_event","cv_event",size,size);
    fCv->Range(0,0,size,size);
  }
  if(fLn==0x0) {
    fLn=new TLine();
    fLn->SetLineColor(kRed);
  }
  for ( size_t iV = 0; iV < vec.size(); ++iV ) {
   cout << vec[iV] << endl;
    for( int iE=0; iE<6; ++iE ) {
      if (vec[iV].fElements[iE]!=-1) {
        const int iD = vec[iV].fElements[iE];
        const int id0 = fDoublets[iE][iD].id0;//fIndex[iL][fDoublets[iL][iD].id0];
        const int id1 = fDoublets[iE][iD].id1;//fIndex[iE+1][fDoublets[iE][iE].id1];
        fLn->DrawLine(cX+fClusters[iE][id0].x*sc,cY+fClusters[iE][id0].y*sc,cX+fClusters[iE+1][id1].x*sc,cY+fClusters[iE+1][id1].y*sc);
      }
    }
  }  
}


#endif
