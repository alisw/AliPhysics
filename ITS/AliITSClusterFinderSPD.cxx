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


#include "AliITSClusterFinderSPD.h"
#include "AliRun.h"



ClassImp(AliITSClusterFinderSPD)

//----------------------------------------------------------
AliITSClusterFinderSPD::AliITSClusterFinderSPD
(AliITSsegmentation *seg, TClonesArray *digits, TClonesArray *recp)   
{
  // constructor
    fSegmentation=seg;
    fDigits=digits;
    fClusters=recp;
    fNclusters= fClusters->GetEntriesFast();
    SetDx();
    SetDz();
    SetMap();
    SetNCells();
}

//_____________________________________________________________________________
AliITSClusterFinderSPD::AliITSClusterFinderSPD()
{
  // constructor
  fSegmentation=0;
  fDigits=0;
  fClusters=0;
  fNclusters=0;
  SetDx();
  SetDz();
  SetMap();
  SetNCells();
  
}

//__________________________________________________________________________
AliITSClusterFinderSPD::AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fDz = source.fDz ;
  this->fDx = source.fDx ;
  this->fMinNCells = source.fMinNCells ;
  return;
}

//_________________________________________________________________________
AliITSClusterFinderSPD& 
  AliITSClusterFinderSPD::operator=(const AliITSClusterFinderSPD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fDz = source.fDz ;
  this->fDx = source.fDx ;
  this->fMinNCells = source.fMinNCells ;
  return *this;
}

//_____________________________________________________________________________
void AliITSClusterFinderSPD::SetMap()
{
  // set map
  if(!fMap) fMap=new AliITSMapA1(fSegmentation,fDigits);
  
}

//_____________________________________________________________________________

void AliITSClusterFinderSPD::Find1DClusters()
{
  // Find one dimensional clusters, i.e.
  // in r*phi(x) direction for each colunm in z direction
  
  AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
  
  // retrieve the parameters 
  Int_t fNofPixels = fSegmentation->Npz(); 
  Int_t fMaxNofSamples = fSegmentation->Npx();
  
  // read in digits -> do not apply threshold 
  // signal in fired pixels is always 1
  fMap->FillMap();
  
  Int_t nofFoundClusters = 0;
  
  Int_t k,it,m;
  for(k=0;k<fNofPixels;k++) {
    
    Int_t mmax = 10;  // a size of the window for the cluster finding
    
    for(it=0;it<fMaxNofSamples;it++) {
      
      Int_t lclx = 0;
      Int_t xstart = 0;
      Int_t xstop = 0;
      Int_t id = 0;
      Int_t ilcl =0;
      
      for(m=0;m<mmax;m++) {  // find the cluster inside the window
	id = it+m;
	if(id >= fMaxNofSamples) break;    // ! no possible for the fadc 
	
	if(fMap->GetHitIndex(k,id)) {   // start of the cluster
	  lclx += 1;
	  if(lclx == 1) xstart = id;
	  
	}
	
	if(lclx > 0 && !fMap->GetHitIndex(k,id)) {  
	  // end of cluster if a gap exists
	  xstop = id-1;
	  ilcl = 1;
	  break;
	}            
	
      }   //  end of m-loop
      
      if(lclx == 0 && ilcl == 0) it = id; // no cluster in the window,
      // continue the "it" loop
      
      if(id >= fMaxNofSamples && lclx == 0) break; // the x row finished
      
      if(id < fMaxNofSamples && ilcl == 0 && lclx > 0) {  
	// cluster end is outside of the window,
	mmax += 5;                 // increase mmax and repeat the cluster
	// finding
	it -= 1;
      }
      
      if(id >= fMaxNofSamples && lclx > 0) {  // the x row finished but
	xstop = fMaxNofSamples - 1;          // the end cluster exists
	ilcl = 1;
      } 
      
      // ---  Calculate z and x coordinates for one dimensional clusters
      
      if(ilcl == 1) {         // new cluster exists
	it = id;
	mmax = 10;
	    nofFoundClusters++;
	    Float_t clusterCharge = 0.;
            // get this from segmentation when this will be implemented
            Float_t zpitch = fSegmentation->Dpz(k+1); 
	    Float_t clusterZ, dummyX; 
            Int_t dummy=0;
            fSegmentation->GetCellCxz(dummy,k,dummyX,clusterZ);
            Float_t zstart = clusterZ - 0.5*zpitch;
            Float_t zstop = clusterZ + 0.5*zpitch;
	    Float_t clusterX = 0.;
            Int_t xstartfull = xstart;
            Int_t xstopfull = xstop;
	    Int_t clusterSizeX = lclx;
	    Int_t clusterSizeZ = 1;
	    
            Int_t its;
	    for(its=xstart; its<=xstop; its++) {
              Int_t firedpixel=0;
              if (fMap->GetHitIndex(k,its)) firedpixel=1; 
	      clusterCharge += firedpixel;
              clusterX +=its + 0.5;
	    }
            Float_t fRphiPitch = fSegmentation->Dpx(dummy);
	    clusterX /= (clusterSizeX/fRphiPitch); // center of gravity for x 
	    
	    
	    //printf("ClusterZ ClusterX %f %f \n",clusterZ, clusterX);
	    
            // Int_t nclusters = fClusters->GetEntriesFast();
	    //	    	    cout << nclusters << " clusters" << endl;
	    //            cout<< "Create point"<<endl;
	    
	    // Write the points (coordinates and some cluster information) to the
	    // AliITSRawClusterSPD object
	    
	    AliITSRawClusterSPD *clust = new AliITSRawClusterSPD(clusterZ,clusterX,clusterCharge,clusterSizeZ,clusterSizeX,xstart,xstop,xstartfull,xstopfull,zstart,zstop,k);
	    // fClusters->Add(point);
	    iTS->AddCluster(0,clust);
	    //    	    cout << "Cluster at Ladder: " << fLadder << ", Detector: " <<fDetector<<endl;
	    
	    // cout<<" end of cluster finding for Z pixel "<<endl;
	    
      }    // new cluster (ilcl=1)
    } // X direction loop (it)
  } // Z direction loop (k)
  return;
  
}

//_____________________________________________________________________________
void  AliITSClusterFinderSPD::GroupClusters()
{
  // Find two dimensional clusters, i.e. group one dimensional clusters
  // into two dimensional ones (go both in x and z directions).
  
  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;
  //printf("Group: fNclusters nofClusters %d %d\n",fNclusters, nofClusters);
  
  AliITSRawClusterSPD *clusterI;
  AliITSRawClusterSPD *clusterJ;
  
  Int_t *label = new Int_t [nofClusters];  // activate this for DEC machines
  Int_t i,j;
  for(i=0; i<nofClusters; i++) label[i] = 0;
  for(i=0; i<nofClusters; i++) {
    if(label[i] != 0) continue;
    for(j=i+1; j<nofClusters; j++) { 
      if(label[j] != 0) continue;
      clusterI = (AliITSRawClusterSPD*) fClusters->At(i);
      clusterJ = (AliITSRawClusterSPD*) fClusters->At(j);
      Bool_t pair = clusterI->Brother(clusterJ,fDz,fDx);
      if(pair) {     
	
	//    if((clusterI->XStop() == clusterJ->XStart()-1)||(clusterI->XStart()==clusterJ->XStop()+1)) cout<<"!! Diagonal cluster"<<endl;
	/*    
	      cout << "clusters " << i << "," << j << " before grouping" << endl;
	      clusterI->Print();
	      clusterJ->Print();
	*/    
	clusterI->Add(clusterJ);
	//        cout << "remove cluster " << j << endl;
	label[j] = 1;
	fClusters->RemoveAt(j);
	/*
	  cout << "cluster  " << i << " after grouping" << endl;
	  clusterI->Print();
	*/
      }  // pair
    } // J clusters  
    label[i] = 1;
  } // I clusters
  fClusters->Compress();
  //Int_t totalNofClusters = fClusters->GetEntriesFast();
  //cout << " Nomber of clusters at the group end ="<< totalNofClusters<<endl;
  
  delete [] nofClusters;
  delete nofClusters;
  
  return;
  
  
}
//_____________________________________________________________________________

void AliITSClusterFinderSPD::GetRecPoints()
{
  // get rec points
  AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
  
  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;
  //printf("GetRecP: fNclusters nofClusters %d %d\n",fNclusters, nofClusters);
  
  const Float_t kconv = 1.0e-4;
  const Float_t kRMSx = 12.0*kconv; // microns -> cm ITS TDR Table 1.3
  const Float_t kRMSz = 70.0*kconv; // microns -> cm ITS TDR Table 1.3

  Int_t i;
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSPD *clusterI = (AliITSRawClusterSPD*) fClusters->At(i);
    AliITSRecPoint rnew;
    rnew.SetX(clusterI->X()*kconv);
    rnew.SetZ(clusterI->Z()*kconv);
    rnew.SetQ(1.);
    rnew.SetdEdX(0.);
    rnew.SetSigmaX2(kRMSx*kRMSx);
    rnew.SetSigmaZ2(kRMSz*kRMSz);
    rnew.SetProbability(1.);
    iTS->AddRecPoint(rnew);
  } // I clusters
  
}
//_____________________________________________________________________________

void AliITSClusterFinderSPD::FindRawClusters()
{
  // find raw clusters
  Find1DClusters();
  GroupClusters();
  GetRecPoints();
}

