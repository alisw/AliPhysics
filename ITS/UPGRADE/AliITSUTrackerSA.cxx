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

// Vc library
//#include "Vc/Vc"
//#include "AliITSUTrackerSAauxVc.h" // Structs and other stuff using Vc library  

#include "AliESDEvent.h"
#include "AliITSUClusterPix.h"
#include "AliITSUTrackerSA.h"

//#include "AliITSUtrackSA.h"      // Some dedicated SA track class ?  

ClassImp(AliITSUTrackerSA)

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA() : AliTracker(), 
  fClusters(),
  fClustersTC(),
  fDoublets(),
  fIndex(),
  fNClusters(),
  fNDoublets(),
  fPhiCut(0.05),
  fRPhiCut(0.01),
  fZCut(0.005)
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------
  for(Int_t i=0;i<7;++i) {
    fClusters[i].reserve(5000);
  }
}

//________________________________________________________________________________
AliITSUTrackerSA::AliITSUTrackerSA(const AliITSUTrackerSA &t): 
  AliTracker(t),
  fClusters(),
  fClustersTC(),
  fIndex(),
  fNClusters(),
  fNDoublets(),
  fPhiCut(),
  fRPhiCut(),
  fZCut()
{
  //--------------------------------------------------------------------
  // The copy constructor is protected
  //--------------------------------------------------------------------
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::Clusters2Tracks(AliESDEvent */*event*/) {
  //--------------------------------------------------------------------
  // This is the main tracking function
  // The clusters must already be loaded
  //--------------------------------------------------------------------

  // Possibly, create the track "seeds" (combinatorial)

  // Possibly, increment the seeds with additional clusters (Kalman)

  // Possibly, (re)fit the found tracks 

  // Tree iterations: 
  // - High momentum first;
  // - Low momentum with vertex constraint; 
  // - Everything else. 

  MakeDoublets();       // To be implemented
  //MakeTriplets();       // Are triplets really necessary? MFT does not use them. 
  CASelection();        // To be implemented
  GlobalFit();          // To be implemented
  ChiSquareSelection(); // To be implemented

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::PropagateBack(AliESDEvent */*event*/) {
  //--------------------------------------------------------------------
  // Here, we implement the Kalman smoother ?
  // The clusters must be already loaded
  //--------------------------------------------------------------------

  return 0;
}

//________________________________________________________________________________
Int_t AliITSUTrackerSA::RefitInward(AliESDEvent */*event*/) {
  //--------------------------------------------------------------------
  // Some final refit, after the outliers get removed by the smoother ?  
  // The clusters must be loaded
  //--------------------------------------------------------------------

  return 0;
}

Int_t AliITSUTrackerSA::LoadClusters(TTree *cluTree) {
  //--------------------------------------------------------------------
  // This function reads the ITSU clusters from the tree,
  // sort them, distribute over the internal tracker arrays, etc
  //--------------------------------------------------------------------
  
  for(Int_t i=0;i<7;++i) {
    TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",i));
    if (!br) return -1;
    br->SetAddress(&fClustersTC[i]);
  }    
  cluTree->GetEntry(0);
  
  for(int iL=0; iL<7; ++iL) {
    TClonesArray *clCont=&fClustersTC[iL];
    fNClusters[iL]=clCont->GetEntriesFast();
    Float_t phi[fNClusters[iL]];
    fIndex[iL] = new Int_t[fNClusters[iL]];
    for(int iC=0;iC<fNClusters[iL];++iC) {
      const AliITSUClusterPix *cl = (AliITSUClusterPix*)clCont->At(iC);
      float pos[3];
      cl->GetGlobalXYZ(pos);
      phi[iC] = pos[0]==0.f ? TMath::PiOver2() : TMath::ATan2(pos[1]-GetY(),pos[0]-GetX());
      float angle=0.f; // TO BE UNDERSTOOD: DO I STILL NEED THE STATION ANGLE IF I USE THE GLOBAL COVARIANCE MATRIX?
      fClusters[iL].push_back(itsCluster(pos[0],pos[1],pos[2],cl->GetSigmaY2(),cl->GetSigmaZ2(),cl->GetSigmaYZ(),phi[iC],angle));
    }
    TMath::Sort(fNClusters[iL],phi,fIndex[iL],kFALSE);
  }
  
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
void AliITSUTrackerSA::CASelection() {
  // Here it's implemented the Cellular Automaton routine
  // Firstly the level of each doublet is set according to the level of 
  // the neighbour doublets. 
  // Doublet are considered to be neighbour if they share one point and the 
  // phi and theta direction difference of the two is below a cut value.

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
	    const float r2 = Sqrt( clusters3[id2].x * clusters2[id2].x + clusters2[id2].y * clusters2[id2].y );
	    const float extrZ3 = doublets1[iD1].tanLambda * ( r3 - r2 ) + clusters3[id3].z ;
	    if ( Abs ( extrZ3 - clusters3[id3].z ) < fZCut ) {
	      const float det = (GetX() - clusters2[id2].x)*(GetY() - clusters1[id1].y) - (GetY() - clusters2[id2].y)*(GetX() - clusters1[id1].x);
	      if ( Abs(det) <= 1e-12 ) {
		// linear extrapolation to next layer
		const float dsq = ( doublets1[iD1].tanPhi * (clusters3[id3].x + clusters2[id2].x) + clusters3[id3].y - clusters2[id2].y ) * 
		  ( doublets1[iD1].tanPhi * (clusters3[id3].x + clusters2[id2].x) + clusters3[id3].y - clusters2[id2].y ) / (1 + doublets1[iD1].tanPhi * doublets1[iD1].tanPhi );
		if ( dsq < fRPhiCut*fRPhiCut )  {
		  doublets2[iD2].level += doublets1[iD1].level;
		  doublets2[iD2].neighbours.push_back(iD1);
		}
	      } else {
		const float r1sq = clusters1[id1].x * clusters1[id1].x + clusters1[id1].y * clusters1[id1].y ;
		const float rvsq = GetX() * GetX() + GetY() * GetY();
		const float deta = (r2*r2-rvsq) * (GetY() - clusters1[id1].y) - (r1sq-rvsq) * (GetY() - clusters2[id2].y);
		const float detb =  (r2*r2-rvsq) * (GetX() - clusters1[id1].x) - (r1sq-rvsq) * (GetX() - clusters2[id2].x) ;
		const float a = deta/det ;
		const float b = detb/det ;
		const float c = -rvsq - a * GetX() - b * GetY();
		const float rc = Sqrt( a*a/4.f + b*b/4.f - c );
		const float d = Sqrt( (a/2.f + clusters3[id3].x) * (a/2.f + clusters3[id3].x) + (b/2.f + clusters3[id3].y) * (b/2.f + clusters3[id3].y) );
		if ( ( d - rc ) < fRPhiCut ) {
		  doublets2[iD2].level += doublets1[iD1].level;
		  doublets2[iD2].neighbours.push_back(iD1);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for ( int level = 6; level >=3 ; --level ) {
    //vector<int> points[6];
    for ( int doubl = 5; doubl >= 0; --doubl ) {
      if( ( doubl + 1 - level ) < 0 ) break;
      for ( int iD = 0; iD < fNDoublets[doubl]; ++iD ) {
	if ( fDoublets[doubl][iD].level == level ) {
	  
	}
      }
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

  for( int iL = 0 ; iL < 6 ; ++iL ) {
    fNDoublets[iL] = 0; 
    const itsCluster* clusters1 = &fClusters[iL][0];
    const itsCluster* clusters2 = &fClusters[iL+1][0];

    // 0 - 2Pi junction treatment (part I)
    for ( int iC1 = 0 ; iC1 < fNClusters[iL] ; ++iC1 ) {
      bool flag = true;
      for ( int iC2 = fNClusters[iL]-1; iC2 > -1 ; --iC2 ) {
	
	if( (TMath::TwoPi() - (clusters2[iC2].phi-clusters1[iC1].phi) ) < fPhiCut ) {
	  fDoublets[iL].push_back(nPlets(iC1,iC2));
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

    
    // "Central" points 
    for ( int iC1 = 0 ; iC1 < fNClusters[iL] ; ++iC1 ) {

      for ( int iC2 = 0; iC2 < fNClusters[iL+1] ; ++iC2 ) {
	
	if( Abs( clusters1[iC1].phi - clusters2[iC2].phi ) < fPhiCut ) {
	  fDoublets[iL].push_back(nPlets(iC1,iC2));
	  fDoublets[iL][fNDoublets[iL]].tanPhi = (clusters1[iC1].y-clusters2[iC2].y)/(clusters1[iC1].x-clusters2[iC2].x);
	  float r1  = Sqrt(clusters1[iC1].x * clusters1[iC1].x + clusters1[iC1].y * clusters1[iC1].y);
	  float r2  = Sqrt(clusters2[iC2].x * clusters2[iC2].x + clusters2[iC2].y * clusters2[iC2].y);
	  fDoublets[iL][fNDoublets[iL]].tanLambda = (clusters1[iC1].z-clusters2[iC2].z)/(r1-r2);
	  ++fNDoublets[iL];
	} else if( clusters2[iC2].phi - clusters1[iC1].phi > fPhiCut ) break;
      
      }

    }

    // 0 - 2Pi junction treatment (part II)
    for ( int iC1 = fNClusters[iL]-1; iC1 > -1 ; --iC1 ) {
      bool flag = true;

      for ( int iC2 = 0; iC2 < fNClusters[iL+1] ; ++iC2 ) {

	if( (TMath::TwoPi() - (clusters1[iC1].phi-clusters2[iC2].phi) ) < fPhiCut ) { 
	  fDoublets[iL].push_back(nPlets(iC1,iC2));
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

}
