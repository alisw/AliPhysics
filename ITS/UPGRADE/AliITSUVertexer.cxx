#include <Riostream.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include "AliESDVertex.h"
#include "AliITSUClusterLines.h"
#include "AliITSUClusterPix.h"
#include "AliITSUVertexer.h"
#include "AliLog.h"
#include "AliStrLine.h"
#include "AliVertexerTracks.h"

using TMath::Abs;
using TMath::Sqrt;
using TMath::ATan2;
using TMath::TwoPi;
using TMath::BubbleLow;

using std::cout;
using std::endl;

//////////////////////////////////////////////////////////////////////
// This class is used to compute the position of all the primary    // 
// vertices in a single event using the upgraded ITS.               //
// Optimizations ongoing.                                           //
// Origin puccio@to.infn.it  Feb. 20 2014                           //
//////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________________
AliITSUVertexer::AliITSUVertexer(Double_t phicut, Double_t zcut, Double_t paircut,Double_t clcut, Int_t cclcut) : AliVertexer(),
														  fClusterContribCut(cclcut),
														  fClusterCut(clcut),
														  fClusterIndex(),
														  fClusterPhi(),
														  fClusters(),
														  fLines("AliStrLine",1000),
														  fLinesClusters("AliITSUClusterLines.h",1000),
														  fLinesPhi(0),
														  fNoClusters(0),
														  fNoLines(0),
														  fNoVertices(0),
														  fPairCut(paircut),
														  fPhiCut(phicut),
														  fZCut(zcut),
														  fUsedClusters(),
														  fUsedLines(),
														  fVertices(NULL) 
#ifdef MC_CHECK 
,fGoodLines(0),fGoodLinesPhi(0),fParticleId(0)
#endif
{
  // Standard I/O constructor
}

//_____________________________________________________________________________________________
AliITSUVertexer::~AliITSUVertexer() {
  // Destructor
  Reset();
}

//_____________________________________________________________________________________________
void AliITSUVertexer::FindVerticesForCurrentEvent() {
  // Try to find all the primary vertices in the current 
  fNoVertices=0;
  FindTracklets();
  if(fNoLines<2) { 
    //fVertices.push_back(AliESDVertex());
    return;// AliESDVertex();
  }
  
  //  fVertices.push_back(AliVertexerTracks::TrackletVertexFinder(&fLines,1));
  //fNoVertices=1;
  fUsedLines=new Short_t[fNoLines];
  for(UInt_t i=0;i<fNoLines;++i) fUsedLines[i]=-1;
  
  fNoClusters=0;
  for(UInt_t i1=0;i1<fNoLines;++i1) {
    if(fUsedLines[i1]!=-1) continue;
    AliStrLine* line1 = (AliStrLine*)fLines.At(i1);
    for(UInt_t i2=i1+1;i2<fNoLines;++i2) {
      if(fUsedLines[i2]!=-1) continue;
      AliStrLine* line2 = (AliStrLine*)fLines.At(i2);
      if(line1->GetDCA(line2)<=fPairCut) {
	//cout << fNoClusters <<" " << i1 << " " << i2 << " ";
	new(fLinesClusters[fNoClusters])AliITSUClusterLines(i1,line1,i2,line2);
	AliITSUClusterLines* current=(AliITSUClusterLines*)fLinesClusters.At(fNoClusters);
	Double_t p[3];
	current->GetVertex(p);
	if((p[0]*p[0]+p[1]*p[1])>=4) { // Beam pipe check
	  fLinesClusters.RemoveAt(fNoClusters);
	  fLinesClusters.Compress();
	  break;
	}
	fUsedLines[i1]=fNoClusters;
	fUsedLines[i2]=fNoClusters;
	for(UInt_t i3=0;i3<fNoLines;++i3) {
	  if(fUsedLines[i3]!=-1) continue;
	  AliStrLine *line3 = (AliStrLine*)fLines.At(i3);
	  //cout << p[0] << " " << p[1] << " " << p[2] << endl;
	  //line3->PrintStatus();
	  if(line3->GetDistFromPoint(p)<=fPairCut) {
	    //cout << i3 << " ";
	    current->Add(i3,line3);
	    fUsedLines[i3]=fNoClusters;
	    current->GetVertex(p);
	  }
	}
	++fNoClusters;
	//cout << endl;
	break;
      }
    }
  }
  
  fLinesClusters.Sort();

  for(UInt_t i0=0;i0<fNoClusters; ++i0) {
    Double_t p0[3],p1[3];
    AliITSUClusterLines *clu0 = (AliITSUClusterLines*)fLinesClusters.At(i0);
    clu0->GetVertex(p0);
    for(UInt_t i1=i0+1;i1<fNoClusters; ++i1) {
      AliITSUClusterLines *clu1 = (AliITSUClusterLines*)fLinesClusters.At(i1);
      clu1->GetVertex(p1);
      if (TMath::Abs(p0[2]-p1[2])<=fClusterCut) {
	Double_t distance=(p0[0]-p1[0])*(p0[0]-p1[0])+(p0[1]-p1[1])*(p0[1]-p1[1])+(p0[2]-p1[2])*(p0[2]-p1[2]);
	//Bool_t flag=kFALSE;
	if(distance<=fPairCut*fPairCut) {
	  UInt_t n=0;
	  Int_t *labels=clu1->GetLabels(n);
	  for(UInt_t icl=0; icl<n; ++icl) clu0->Add(labels[icl],(AliStrLine*)fLines.At(labels[icl]));
	  clu0->GetVertex(p0);
	  //flag=kTRUE;
	}
	fLinesClusters.RemoveAt(i1);
	fLinesClusters.Compress();
	fNoClusters--;
	i1--;
	//if(flag) i1=10;
      }
    }
  }

  fVertices=new AliESDVertex[fNoClusters];
  for(UInt_t i0=0; i0<fNoClusters; ++i0) {
    AliITSUClusterLines *clu0 = (AliITSUClusterLines*)fLinesClusters.At(i0);
    Int_t size=clu0->GetSize();
    if(size<fClusterContribCut&&fNoClusters>1) {
      fLinesClusters.RemoveAt(i0);
      fLinesClusters.Compress();
      fNoClusters--;
      continue;
    } 
    Double_t p0[3],cov[6];
    clu0->GetVertex(p0);
    clu0->GetCovMatrix(cov);
    if((p0[0]*p0[0]+p0[1]*p0[1])<1.98*1.98) {
      fVertices[fNoVertices++]=AliESDVertex(p0,cov,99999.,size);   
    }
  }
  
  return;// AliVertexerTracks::TrackletVertexFinder(&fLines,0);
}

//______________________________________________________________________
AliESDVertex* AliITSUVertexer::FindVertexForCurrentEvent(TTree *cluTree)
{
// Reconstruction of all the primary vertices in the event. It returns the vertex with the highest number of contributors.  
  Reset();
  for(Int_t i=0;i<3;++i) {
    TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",i));
    if (!br) return NULL;
    br->SetAddress(&fClusters[i]);
  }    
  cluTree->GetEntry(0);
  
  SortClusters();

  FindVerticesForCurrentEvent();
  if(fNoVertices<1) return NULL;
  return new AliESDVertex(fVertices[0]);
}  

//_____________________________________________________________________________________________
Int_t AliITSUVertexer::MatchPoints(UShort_t layer, Double_t anchor, Double_t *p0, Double_t *p1) {
  // Method for matching clusters with similar phi and z inside the range [zmin,zmax]
  AliDebug(3,Form("Matching points on layer %i...\n",layer));

  Double_t xyz[3][3];        // {{0: 'min', 1: 'mid', 2: 'max'},{0: 'x', 1: 'y', 2: 'z'}}
  AliITSUClusterPix* cl[3];  // {0: 'min', 1: 'mid', 2: 'max'}
  Double_t phi[3];           // {0: 'min', 1: 'mid', 2: 'max'}

  Int_t nocl=fClusters[layer]->GetEntriesFast();
  Int_t imin=0,imax=nocl!=0 ? nocl-1 : nocl,imid=(imax+imin)/2;

  Double_t a=0,b=0,c=0;
  if(layer==2) {
    a=(p0[0]-p1[0])*(p0[0]-p1[0])+(p0[1]-p1[1])*(p0[1]-p1[1]);
    b=(p1[0]-p0[0])*p0[0]+(p1[1]-p0[1])*p0[1];
    c=p0[0]*p0[0]+p0[1]*p0[1];
  }

  Int_t flag=-1;

  while (imax > imin) {
    while(fUsedClusters[layer][fClusterIndex[layer][imin]]==kTRUE&&imin<imax) ++imin;
    while(fUsedClusters[layer][fClusterIndex[layer][imax]]==kTRUE&&imax>0) --imax;
    if(imax<imin) return -1;    
    
    cl[0] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][imin]);
    cl[0]->GetGlobalXYZ(xyz[0]);
    phi[0] = fClusterPhi[layer][fClusterIndex[layer][imin]];
    AliDebug(4,Form("Cluster min: %i, %i, %f, %i\n",imin,cl[0]->GetLabel(0),phi[0]-anchor,fUsedClusters[layer][fClusterIndex[layer][imin]]));
    if((Abs(phi[0]-anchor)<=fPhiCut||Abs(Abs(phi[0]-anchor)-TMath::TwoPi())<=fPhiCut)) {
      AliDebug(4,Form("Ok min: %i \n",imin));
      if(layer==2) {
	if(fUsedClusters[layer][fClusterIndex[layer][imin]]) {
	  flag=1;
	  break;
	}
	c-=(xyz[0][0]*xyz[0][0]+xyz[0][1]*xyz[0][1]);
	Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
	Double_t zmin=z-fZCut,zmax=z+fZCut;
	if(zmin<=xyz[0][2]&&xyz[0][2]<=zmax) {
	  AliDebug(4,Form("Ok Z: %i \n",imin));
	  return fClusterIndex[layer][imin];
	}
	AliDebug(4,Form("No Z match: %i \n",imin));
	c=p0[0]*p0[0]+p0[1]*p0[1];
	flag=1;
	break;
      } else if(!fUsedClusters[layer][fClusterIndex[layer][imin]]) return fClusterIndex[layer][imin];
    }
    
    cl[2] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][imax]);
    cl[2]->GetGlobalXYZ(xyz[2]);
    phi[2] = fClusterPhi[layer][fClusterIndex[layer][imax]];
    AliDebug(4,Form("Cluster max: %i, %i, %f, %i\n",imax,cl[2]->GetLabel(0),phi[2]-anchor,fUsedClusters[layer][fClusterIndex[layer][imax]]));
    if((Abs(phi[2]-anchor)<=fPhiCut||Abs(Abs(phi[2]-anchor)-TMath::TwoPi())<=fPhiCut)) {
      AliDebug(4,Form("Ok max: %i \n",imax));
      if(layer==2) {
	if(fUsedClusters[layer][fClusterIndex[layer][imax]]) {
	  flag=3;
	  break;
	}
	c-=(xyz[2][0]*xyz[2][0]+xyz[2][1]*xyz[2][1]);
	Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
	Double_t zmin=z-fZCut,zmax=z+fZCut;
	if(zmin<=xyz[2][2]&&xyz[2][2]<=zmax) {
	  AliDebug(4,Form("Ok Z: %i \n",imax));
	  return fClusterIndex[layer][imax];
	}
	c=p0[0]*p0[0]+p0[1]*p0[1];
	AliDebug(4,Form("No Z match: %i \n",imax));
	flag=3;
	break;
      } else if(!fUsedClusters[layer][fClusterIndex[layer][imax]]) return fClusterIndex[layer][imax];
    }
    
    imid=(imax+imin)/2;
    if(imid==imin) return -1;
    Int_t step=1,sign=-1;
    while(fUsedClusters[layer][fClusterIndex[layer][imid]]) {
      imid=imid+step;
      sign*=-1;
      step=(step+sign)*(-1);
      if(imid>=imax||imid<=imin) return -1;
    }
    
    cl[1] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][imid]);
    cl[1]->GetGlobalXYZ(xyz[1]);
    phi[1] = fClusterPhi[layer][fClusterIndex[layer][imid]];
    AliDebug(4,Form("Cluster mid: %i, %i, %f, %i \n",imid,cl[1]->GetLabel(0),phi[1]-anchor,fUsedClusters[layer][fClusterIndex[layer][imid]]));
    //cout << imin << " " << imid << " " << imax << endl;
    //cout << fClusterIndex[layer][imid] << endl;
    if((Abs(phi[1]-anchor)<=fPhiCut||Abs(Abs(phi[1]-anchor)-TMath::TwoPi())<=fPhiCut)) {
      AliDebug(4,Form("Ok mid: %i \n",imid));
      if(layer==2) {
	c-=(xyz[1][0]*xyz[1][0]+xyz[1][1]*xyz[1][1]);
	Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
	Double_t zmin=z-fZCut,zmax=z+fZCut;
	if(zmin<=xyz[1][2]&&xyz[1][2]<=zmax) {
	  AliDebug(4,Form("Ok Z: %i \n",imid));
	  return fClusterIndex[layer][imid];
	}
	AliDebug(4,Form("No Z match: %i \n",imid));
	c=p0[0]*p0[0]+p0[1]*p0[1];
	flag=2;
	break;
      } else if(!fUsedClusters[layer][fClusterIndex[layer][imid]]) return fClusterIndex[layer][imid];
    }
    
    // determine which subarray to search
    if (phi[1] < anchor) {
      // change min index to search upper subarray
      AliDebug(4,"Case minor\n");
      imin = imid + 1;
      --imax;
    } else if (phi[1] > anchor) {
      AliDebug(4,"Case greater\n");
      // change max index to search lower subarray
      ++imin;
      imax = imid - 1;
    } else if(!fUsedClusters[layer][fClusterIndex[layer][imid]]) {
      return fClusterIndex[layer][imid];
    } else return -1;
  }

  if(flag>-1) {
    AliDebug(4,"Flag issued, starting forward backward check\n");
    Int_t start=imid;
    switch(flag) {
    case 1:
      start=imin;
      break;
    case 2:
      start=imid;
      break;
    case 3:
      start=imax;
      break;
    }
   
    Int_t curr=start-1;
    Bool_t lap=kFALSE;
    while(1){
      if(curr==-1&&!lap) {
	lap=kTRUE;
	curr=nocl-1;
      } else if(lap) break;
      cl[0] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][curr]);
      cl[0]->GetGlobalXYZ(xyz[0]);
      phi[0] = fClusterPhi[layer][fClusterIndex[layer][curr]];
      AliDebug(4,Form("Looking backward: %i, %i, %f, %i\n",curr,cl[0]->GetLabel(0),phi[0]-anchor,fUsedClusters[layer][fClusterIndex[layer][curr]]));
      if(!(Abs(phi[0]-anchor)<=fPhiCut||Abs(Abs(phi[0]-anchor)-TMath::TwoPi())<=fPhiCut)) break;
      if(fUsedClusters[layer][fClusterIndex[layer][curr]]) {
	curr--;
	continue;
      }
      
      AliDebug(4,Form("Ok backward: %i \n",curr));
      
      c-=(xyz[0][0]*xyz[0][0]+xyz[0][1]*xyz[0][1]);
      Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
      Double_t zmin=z-fZCut,zmax=z+fZCut;
      if(zmin<=xyz[0][2]&&xyz[0][2]<=zmax) {
	AliDebug(4,Form("Ok Z: %i \n",curr));
	return fClusterIndex[layer][curr];
      }
      AliDebug(4,Form("No Z match: %i \n",curr));
      c=p0[0]*p0[0]+p0[1]*p0[1];
      curr--;
    }
    
    lap=kFALSE;
    curr=start+1;
    while(1){
      if(curr==nocl&&!lap) {
	curr=0;
	lap=kTRUE;
      } else if(lap) break;
      cl[0] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][curr]);
      cl[0]->GetGlobalXYZ(xyz[0]);
      phi[0] = fClusterPhi[layer][fClusterIndex[layer][curr]];
      AliDebug(4,Form("Looking forward: %i, %i, %f, %i\n",curr,cl[0]->GetLabel(0),phi[0]-anchor,fUsedClusters[layer][fClusterIndex[layer][curr]]));
      if(!(Abs(phi[0]-anchor)<=fPhiCut||Abs(Abs(phi[0]-anchor)-TMath::TwoPi())<=fPhiCut)) break;
      if(fUsedClusters[layer][fClusterIndex[layer][curr]]) {
	curr++;
	continue;
      }
      AliDebug(4,Form("Ok forward: %i \n",curr));
      c-=(xyz[0][0]*xyz[0][0]+xyz[0][1]*xyz[0][1]);
      Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
      Double_t zmin=z-fZCut,zmax=z+fZCut;
      if(zmin<=xyz[0][2]&&xyz[0][2]<=zmax) {
	AliDebug(4,Form("Ok Z: %i \n",curr));
	return fClusterIndex[layer][curr];
      }
      AliDebug(4,Form("No Z match: %i \n",curr));
      c=p0[0]*p0[0]+p0[1]*p0[1];
      curr++;
    }
  }
  
  if(imax==imin&&imax!=0) {
    cl[0] = (AliITSUClusterPix*)fClusters[layer]->At(fClusterIndex[layer][imin]);
    cl[0]->GetGlobalXYZ(xyz[0]);
    phi[0] = fClusterPhi[layer][fClusterIndex[layer][imin]];
    AliDebug(4,Form("Cluster eq: %i, %i, %f, %i\n",imin,cl[0]->GetLabel(0),phi[0]-anchor,fUsedClusters[layer][fClusterIndex[layer][imin]]));
    if((Abs(phi[0]-anchor)<=fPhiCut||Abs(Abs(phi[0]-anchor)-TMath::TwoPi())<=fPhiCut)&&!fUsedClusters[layer][fClusterIndex[layer][imin]]) {
      AliDebug(4,Form("Ok eq: %i \n",imin));
      if(layer==2) {
        c-=(xyz[0][0]*xyz[0][0]+xyz[0][1]*xyz[0][1]);
        Double_t z=p0[2]+(p1[2]-p0[2])*(-b+Sqrt(b*b-a*c))/a;
        Double_t zmin=z-fZCut,zmax=z+fZCut;
        if(zmin<=xyz[0][2]&&xyz[0][2]<=zmax) {
	  AliDebug(4,Form("Ok Z eq: %i \n",imin));
	  return fClusterIndex[layer][imin];
	}	
	AliDebug(4,Form("No Z eq: %i \n",imin));
	c=p0[0]*p0[0]+p0[1]*p0[1];
      } else return fClusterIndex[layer][imin];
    }
  }
  // no match found :(
  return -1;

}

//_____________________________________________________________________________________________
void AliITSUVertexer::PrintStatus() const {
  // Prints all the cuts and important data members for the current status
  cout << "Cut on phi: " << fPhiCut << endl;
  cout << "Cut on z: " << fZCut << endl;
}

//_____________________________________________________________________________________________
void AliITSUVertexer::Reset() {
  // Resets the vertexer for a new event (or for its destruction) 
  AliDebug(2,"Resetting the vertexer...\n");
  fNoVertices=0;
  for(Int_t i=0;i<3;++i) {
    delete[] fUsedClusters[i];
    delete[] fClusterIndex[i];
    delete[] fClusterPhi[i];
  }
  if(fNoLines>2) {
    delete []fUsedLines;
  }

  delete[] fVertices;

  fLinesPhi=0;
  fLines.Clear();
  fLinesClusters.Clear();
  
  #ifdef MC_CHECK
  fGoodLines=0;
  fGoodLinesPhi=0;
  delete[] fParticleId;
  #endif

  //delete fVertices;
}

//_____________________________________________________________________________________________
void AliITSUVertexer::FindTracklets() {
  // It combines recpoints over the first three layers to define a list of tracklets
  // Here it uses the ordered list of points. Second and third layers are ordered using phi
  AliDebug(2,"Calling the trackleter...\n");

  UInt_t noPntL[3]; 
  for(UShort_t i=0;i<3;++i) {
    noPntL[i]=fClusters[i]->GetEntries();
    fUsedClusters[i]=new Bool_t[noPntL[i]];
    for(UInt_t ii=0;ii<noPntL[i];++ii) fUsedClusters[i][ii]=kFALSE;
  }
  
  #ifdef MC_CHECK
  fParticleId=new UInt_t[noPntL[0]];
  #endif

  fNoLines=0;
  UInt_t nolinesphi=0;
  Double_t p0[3],p1[3],pp2[3];
  for(UInt_t i0=0;i0<noPntL[0];++i0) {
    if(fUsedClusters[0][i0]) continue;
    AliITSUClusterPix* cluster0=(AliITSUClusterPix*)fClusters[0]->At(i0);
    cluster0->GetGlobalXYZ(p0);
    vector<Int_t> tmp;
    Int_t i1=0;
    Int_t label0=cluster0->GetLabel(0);  
    AliDebug(4,Form("Layer 0: %i\n",label0));
    while(i1>=0) {
      i1 = MatchPoints(1,fClusterPhi[0][i0]);
      if(i1<0) break;
      
      ++nolinesphi;
      AliITSUClusterPix* cluster1=(AliITSUClusterPix*)fClusters[1]->At(i1);
      cluster1->GetGlobalXYZ(p1);
      tmp.push_back(i1);
      fUsedClusters[1][i1]=kTRUE;

      Int_t i2 = MatchPoints(2,fClusterPhi[1][i1],p0,p1);
      if(i2<0) continue;
      
      #ifdef MC_CHECK
      CheckMC(i0,i1,i2);
      #endif

      AliITSUClusterPix* cluster2=(AliITSUClusterPix*)fClusters[2]->At(i2);
      cluster2->GetGlobalXYZ(pp2);

      fUsedClusters[0][i0]=kTRUE;
      fUsedClusters[2][i2]=kTRUE;
      
      // Errors to be checked... 
      
      Float_t cov0[6],cov1[6],cov2[6];
      cluster0->GetGlobalCov(cov0);
      cluster1->GetGlobalCov(cov1);
      cluster2->GetGlobalCov(cov2);
      
      //Error on tracklet direction near the vertex
      Double_t rad1=TMath::Sqrt(p0[0]*p0[0]+p0[1]*p0[1]);
      Double_t rad2=TMath::Sqrt(p1[0]*p1[0]+p1[1]*p1[1]);
      Double_t factor=(rad1+rad2)/(rad2-rad1);
      
      //Curvature error
      Double_t curvErr=0;
      Double_t bField=0.5;
      Double_t meanPtSelTrk=0.630;
      Double_t curvRadius=meanPtSelTrk/(0.3*bField)*100; //cm 
      Double_t dRad=TMath::Sqrt((p0[0]-p1[0])*(p0[0]-p1[0])+(p0[1]-p1[1])*(p0[1]-p1[1]));
      Double_t aux=dRad/2.+rad1;
      curvErr=TMath::Sqrt(curvRadius*curvRadius-dRad*dRad/4.)-TMath::Sqrt(curvRadius*curvRadius-aux*aux); //cm
      
      Double_t sq[3],wmat[9]={1,0,0,0,1,0,0,0,1};
      sq[0]=(cov0[0]+curvErr*curvErr/2.)*factor*factor;//cov1[0]+cov2[0]);
      sq[1]=(cov0[3]+curvErr*curvErr/2.)*factor*factor;//cov1[3]+cov2[3]);
      sq[2]=(cov0[5])*factor*factor;//cov1[5]+cov2[5]);	    

      // Multiple scattering
      Double_t meanPSelTrk=0.875;
      Double_t pOverMass=meanPSelTrk/0.140;
      Double_t beta2=pOverMass*pOverMass/(1+pOverMass*pOverMass);
      Double_t p2=meanPSelTrk*meanPSelTrk;
      Double_t rBP=1.98; // Beam pipe radius
      Double_t dBP=0.08/35.3; // 800 um of Be
      Double_t dL1=0.01; //approx. 1% of radiation length  
      Double_t theta2BP=14.1*14.1/(beta2*p2*1e6)*dBP;
      Double_t theta2L1=14.1*14.1/(beta2*p2*1e6)*dL1;
      Double_t rtantheta1=(rad2-rad1)*TMath::Tan(TMath::Sqrt(theta2L1));
      Double_t rtanthetaBP=(rad1-rBP)*TMath::Tan(TMath::Sqrt(theta2BP));
      for(Int_t ico=0; ico<3;ico++){    
	sq[ico]+=rtantheta1*rtantheta1*factor*factor/3.;
	sq[ico]+=rtanthetaBP*rtanthetaBP*factor*factor/3.;
      }
      
      if(sq[0]!=0) wmat[0]=1/sq[0];
      if(sq[1]!=0) wmat[4]=1/sq[1];
      if(sq[2]!=0) wmat[8]=1/sq[2];

      /*Int_t label0=cluster0->GetLabel(0);
      Int_t label1=cluster1->GetLabel(0);
      Int_t label2=cluster2->GetLabel(0);*/
      //cout << label0 << " " << label1 << " "<<label2 << endl;

      //printf("%f\t%f\t%f\n-\t%f\t%f\n-\t-\t%f\n",wmat[0],wmat[1],wmat[2],wmat[3],wmat[4],wmat[5]);
      //gSystem->Exec(Form("echo %i >> trkl",cluster0->GetLabel(0)+cluster1->GetLabel(0)));
      
      new(fLines[fNoLines++])AliStrLine(p0,sq,wmat,p1,kTRUE);
      break;
    }
    for(UInt_t itmp=0; itmp<tmp.size(); ++itmp) fUsedClusters[1][tmp.at(itmp)]=kFALSE;
    tmp.clear();
    //((AliStrLine*)fLines[fNoLines-1])->PrintStatus();
    //printf("(%f,%f,%f) and (%f,%f,%f)\n",p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
  }
  
  fLinesPhi=nolinesphi;
  //cout << Form("=======================================================\nNolinesphi: %i, nolines: %i\nGood Nolinesphi: %i, good nolines: %i\n",nolinesphi,fNoLines,fGoodLinesPhi,fGoodLines);
}
 
//___________________________________________________________________________
void AliITSUVertexer::SortClusters() { 
  // Reading of the clusters on the first three layer of the upgraded ITS
  for(Int_t i=0;i<3;++i) {
    TClonesArray *clr=fClusters[i];
    Int_t nocl=clr->GetEntriesFast();
    if(nocl==0) {
      fClusterPhi[i]=new Double_t[1];
      fClusterPhi[i][0]=-999999;
      fClusterIndex[i]=new Int_t[1];
      fClusterIndex[i][0]=0;
    } else {
      fClusterPhi[i]=new Double_t[nocl];
      fClusterIndex[i]=new Int_t[nocl];
      for(Int_t k=0;k<nocl;++k) {
	AliITSUClusterPix* cl=(AliITSUClusterPix*)clr->At(k);
	Double_t pt[3]; 
	cl->GetGlobalXYZ(pt);
	fClusterPhi[i][k]=ATan2(pt[1],pt[0]);
      }
      BubbleLow(nocl,fClusterPhi[i],fClusterIndex[i]);
    }
  }
}

#ifdef MC_CHECK
//___________________________________________________________________________
Bool_t AliITSUVertexer::CheckMC(UInt_t i0, UInt_t i1, UInt_t i2) {
  // Debugging function
  //cout << "Checking MC truth" << endl;
  Int_t label=0;
  Bool_t flag=kFALSE;
  AliITSUClusterPix* p0=(AliITSUClusterPix*)fClusters[0]->At(i0);
  AliITSUClusterPix* p1=(AliITSUClusterPix*)fClusters[1]->At(i1);
  for(Int_t i=0; i<3; ++i) {
    label=p0->GetLabel(i);
    for(Int_t j=0; j<3; ++j) {
      if(label==p1->GetLabel(j)&&label>0) {
        fGoodLinesPhi++;
	flag=kTRUE;			      
	break;
      }
    }
    if(flag) break;
  }
 
  if(!flag) return kFALSE;
  
  AliITSUClusterPix* p2=(AliITSUClusterPix*)fClusters[2]->At(i2);
  for(Int_t j=0; j<3; ++j) {
    if(label==p2->GetLabel(j)) {
      fParticleId[fGoodLines]=label;
      fGoodLines++;
      //cout << label << endl;
      return kTRUE;
    }
  }
  return kFALSE;
}
#endif
