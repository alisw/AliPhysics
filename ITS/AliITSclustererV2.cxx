//-------------------------------------------------------------------------
//            Implementation of the ITS clusterer V2 class
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
//uncomment the line below for running with V1 cluster finder classes 
//#define V1

#include "AliRun.h"

#include "AliITSclustererV2.h"
#include "AliITSclusterV2.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSSD.h"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include "AliITSgeom.h"
#include "AliITSdigitSPD.h"
#include "AliITSdigitSDD.h"
#include "AliITSdigitSSD.h"
#include "AliMC.h"

ClassImp(AliITSclustererV2)

extern AliRun *gAlice;

AliITSclustererV2::AliITSclustererV2():
fNModules(0),
fEvent(0),
fI(0),
fLastSPD1(0),
fNySPD(0),
fNzSPD(0),
fYpitchSPD(0),
fZ1pitchSPD(0),
fZ2pitchSPD(0),
fHwSPD(0),
fHlSPD(0),
fNySDD(0),
fNzSDD(0),
fYpitchSDD(0),
fZpitchSDD(0),
fHwSDD(0),
fHlSDD(0),
fYoffSDD(0),
fLastSSD1(0),
fYpitchSSD(0),
fHwSSD(0),
fHlSSD(0),
fTanP(0),
fTanN(0){
   //default constructor
 }
AliITSclustererV2::AliITSclustererV2(const AliITSgeom *geom):
fNModules(0),
fEvent(0),
fI(0),
fLastSPD1(0),
fNySPD(256),
fNzSPD(160),
fYpitchSPD(0.0050),
fZ1pitchSPD(0.0425),
fZ2pitchSPD(0.0625),
fHwSPD(0.64),
fHlSPD(3.48),
fNySDD(256),
fNzSDD(256),
fYpitchSDD(0.01825),
fZpitchSDD(0.02940),
fHwSDD(3.5085),
fHlSDD(3.7632),
fYoffSDD(0.0425),
fLastSSD1(0),
fYpitchSSD(0.0095),
fHwSSD(3.65),
fHlSSD(2.00),
fTanP(0.0275),
fTanN(0.0075) {
  //------------------------------------------------------------
  // Standard constructor
  //------------------------------------------------------------
  AliITSgeom *g=(AliITSgeom*)geom;

  Int_t mmax=geom->GetIndexMax();
  if (mmax>2200) {
     Fatal("AliITSclustererV2","Too many ITS subdetectors !"); 
  }
  Int_t m;
  for (m=0; m<mmax; m++) {
     Int_t lay,lad,det; g->GetModuleId(m,lay,lad,det);
     Float_t x,y,z;     g->GetTrans(lay,lad,det,x,y,z); 
     Double_t rot[9];   g->GetRotMatrix(lay,lad,det,rot);
     Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
     Double_t ca=TMath::Cos(alpha), sa=TMath::Sin(alpha);
     fYshift[m] = x*ca + y*sa;
     fZshift[m] = (Double_t)z;
     fNdet[m] = (lad-1)*g->GetNdetectors(lay) + (det-1);
     fNlayer[m] = lay-1;
  }
  fNModules = g->GetIndexMax();

  //SPD geometry  
  fLastSPD1=g->GetModuleIndex(2,1,1)-1;
  fYSPD[0]=0.5*fYpitchSPD;
  for (m=1; m<fNySPD; m++) fYSPD[m]=fYSPD[m-1]+fYpitchSPD; 
  fZSPD[0]=fZ1pitchSPD;
  for (m=1; m<fNzSPD; m++) {
    Double_t dz=fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=fZ2pitchSPD; 
    fZSPD[m]=fZSPD[m-1]+dz;
  }
  for (m=0; m<fNzSPD; m++) {
    Double_t dz=0.5*fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=0.5*fZ2pitchSPD; 
    fZSPD[m]-=dz;
  }

  //SSD geometry
  fLastSSD1=g->GetModuleIndex(6,1,1)-1;

}


Int_t AliITSclustererV2::Digits2Clusters(TTree *dTree, TTree *cTree) {
  //------------------------------------------------------------
  // This function creates ITS clusters
  //------------------------------------------------------------
  Int_t ncl=0;

  if (!dTree) {
    Error("Digits2Clusters","Can't get the tree with digits !");
    return 1;
  }

  TClonesArray *digitsSPD=new TClonesArray("AliITSdigitSPD",3000);
  dTree->SetBranchAddress("ITSDigitsSPD",&digitsSPD);
  TClonesArray *digitsSDD=new TClonesArray("AliITSdigitSDD",3000);
  dTree->SetBranchAddress("ITSDigitsSDD",&digitsSDD);
  TClonesArray *digitsSSD=new TClonesArray("AliITSdigitSSD",3000);
  dTree->SetBranchAddress("ITSDigitsSSD",&digitsSSD);

  TClonesArray *clusters=new TClonesArray("AliITSclusterV2",1000);
  TBranch *branch=cTree->GetBranch("Clusters");
  if (!branch) cTree->Branch("Clusters",&clusters);
  else branch->SetAddress(&clusters);

  Int_t mmax=(Int_t)dTree->GetEntries();
  if (mmax!=fNModules) {
    Error("Digits2Clusters","Number of entries != number of modules !");
    return 1;
  }

  for (fI=0; fI<mmax; fI++) {
    dTree->GetEvent(fI);

    if     (digitsSPD->GetEntriesFast()!=0) 
      FindClustersSPD(digitsSPD,clusters);
    else 
      if(digitsSDD->GetEntriesFast()!=0) 
	FindClustersSDD(digitsSDD,clusters);
      else if(digitsSSD->GetEntriesFast()!=0) 
	FindClustersSSD(digitsSSD,clusters);
    
    ncl+=clusters->GetEntriesFast();

    cTree->Fill();

    digitsSPD->Clear();
    digitsSDD->Clear();
    digitsSSD->Clear();
    clusters->Clear();
  }

  //cTree->Write();

  delete clusters;

  delete digitsSPD;
  delete digitsSDD;
  delete digitsSSD;

  //delete dTree;

  Info("Digits2Clusters","Number of found clusters : %d",ncl);

  return 0;
}

void AliITSclustererV2::Digits2Clusters(AliRawReader* rawReader) {
  //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader();
  if (!runLoader) {
    Error("Digits2Clusters", "no run loader found");
    return;
  }
  runLoader->LoadKinematics();
  AliLoader* itsLoader = runLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    Error("Digits2Clusters", "no loader for ITS found");
    return;
  }
  if (!itsLoader->TreeR()) itsLoader->MakeTree("R");
  TTree* cTree = itsLoader->TreeR();

  TClonesArray *array=new TClonesArray("AliITSclusterV2",1000);
  cTree->Branch("Clusters",&array);
  delete array;

  TClonesArray** clusters = new TClonesArray*[fNModules]; 
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    clusters[iModule] = NULL;
  }
  // one TClonesArray per module

  rawReader->Reset();
  AliITSRawStreamSPD inputSPD(rawReader);
  FindClustersSPD(&inputSPD, clusters);

  rawReader->Reset();
  AliITSRawStreamSDD inputSDD(rawReader);
  FindClustersSDD(&inputSDD, clusters);

  rawReader->Reset();
  AliITSRawStreamSSD inputSSD(rawReader);
  FindClustersSSD(&inputSSD, clusters);

  // write all clusters to the tree
  Int_t nClusters = 0;
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    array = clusters[iModule];
    if (!array) {
      Error("Digits2Clusters", "data for module %d missing!", iModule);
      array = new TClonesArray("AliITSclusterV2");
    }
    cTree->SetBranchAddress("Clusters", &array);
    cTree->Fill();
    nClusters += array->GetEntriesFast();
    delete array;
  }
  itsLoader->WriteRecPoints("OVERWRITE");

  delete[] clusters;

  Info("Digits2Clusters", "total number of found clusters in ITS: %d\n", 
       nClusters);
}

//**** Fast clusters *******************************
#include "TParticle.h"

//#include "AliITS.h"
#include "AliITSmodule.h"
#include "AliITSRecPoint.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITSRecPoint.h"


static void CheckLabels(Int_t lab[3]) {
  //------------------------------------------------------------
  // Tries to find mother's labels
  //------------------------------------------------------------
    Int_t label=lab[0];
    if (label>=0) {
       TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
       label=-3;
       while (part->P() < 0.005) {
          Int_t m=part->GetFirstMother();
          if (m<0) {
             Info("CheckLabels","Primary momentum: %f",part->P()); 
             break;
          }
          if (part->GetStatusCode()>0) {
             Info("CheckLabels","Primary momentum: %f",part->P()); 
             break;
          }
          label=m;
          part=(TParticle*)gAlice->GetMCApp()->Particle(label);
        }
        if      (lab[1]<0) lab[1]=label;
        else if (lab[2]<0) lab[2]=label;
        else ;//cerr<<"CheckLabels : No empty labels !\n";
    }
}

/*
static void CheckLabels(Int_t lab[3]) {
  //------------------------------------------------------------
  // Tries to find mother's labels
  //------------------------------------------------------------

  if(lab[0]<0 && lab[1]<0 && lab[2]<0) return; // In case of no labels just exit

  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  for (Int_t i=0;i<3;i++){
    Int_t label = lab[i];
    if (label>=0 && label<ntracks) {
      TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
      if (part->P() < 0.005) {
	Int_t m=part->GetFirstMother();
	if (m<0) {	
	  continue;
	}
	if (part->GetStatusCode()>0) {
	  continue;
	}
	lab[i]=m;       
      }
    }    
  }
  
}
*/
static void CheckLabels2(Int_t lab[10]) {
  //------------------------------------------------------------
  // Tries to find mother's labels
  //------------------------------------------------------------
  Int_t nlabels =0; 
  for (Int_t i=0;i<10;i++) if (lab[i]>=0) nlabels++;
  if(nlabels == 0) return; // In case of no labels just exit

  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();

  for (Int_t i=0;i<10;i++){
    Int_t label = lab[i];
    if (label>=0 && label<ntracks) {
      TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
      if (part->P() < 0.02) {
	  Int_t m=part->GetFirstMother();
	  if (m<0) {	
	    continue;
	  }
	  if (part->GetStatusCode()>0) {
	    continue;
	  }
	  lab[i]=m;       
      }
      else
	if (part->P() < 0.12 && nlabels>3) {
	  lab[i]=-2;
	  nlabels--;
	} 
    }
    else{
      if ( (label>ntracks||label <0) && nlabels>3) {
	lab[i]=-2;
	nlabels--;
      } 
    }
  }  
  if (nlabels>3){
    for (Int_t i=0;i<10;i++){
      if (nlabels>3){
	Int_t label = lab[i];
	if (label>=0 && label<ntracks) {
	  TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
	  if (part->P() < 0.1) {
	    lab[i]=-2;
	    nlabels--;
	  }
	}
      }
    }
  }

  //compress labels -- if multi-times the same
  Int_t lab2[10];
  for (Int_t i=0;i<10;i++) lab2[i]=-2;
  for (Int_t i=0;i<10  ;i++){
    if (lab[i]<0) continue;
    for (Int_t j=0;j<10 &&lab2[j]!=lab[i];j++){
      if (lab2[j]<0) {
	lab2[j]= lab[i];
	break;
      }
    }
  }
  for (Int_t j=0;j<10;j++) lab[j]=lab2[j];
  
}

static void AddLabel(Int_t lab[10], Int_t label) {

  if(label<0) return; // In case of no label just exit

  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  if (label>ntracks) return;
  for (Int_t i=0;i<10;i++){
    //    if (label<0) break;
    if (lab[i]==label) break;
    if (lab[i]<0) {
      lab[i]= label;
      break;
    }
  }
}

void AliITSclustererV2::RecPoints2Clusters
(const TClonesArray *points, Int_t idx, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Conversion AliITSRecPoint -> AliITSclusterV2 for the ITS 
  // subdetector indexed by idx 
  //------------------------------------------------------------
  TClonesArray &cl=*clusters;
  Int_t ncl=points->GetEntriesFast();
  for (Int_t i=0; i<ncl; i++) {
    AliITSRecPoint *p = (AliITSRecPoint *)points->UncheckedAt(i);
    Float_t lp[5];
    lp[0]=-(-p->GetDetLocalX()+fYshift[idx]); if (idx<=fLastSPD1) lp[0]*=-1; //SPD1
    lp[1]=  -p->GetZ()+fZshift[idx];
    lp[2]=p->GetSigmaDetLocX2();
    lp[3]=p->GetSigmaZ2();
    lp[4]=p->GetQ()*36./23333.;  //electrons -> ADC
    Int_t lab[4]; 
    lab[0]=p->GetLabel(0); lab[1]=p->GetLabel(1); lab[2]=p->GetLabel(2);
    lab[3]=fNdet[idx];
    CheckLabels(lab);
    Int_t dummy[3]={0,0,0};
    new (cl[i]) AliITSclusterV2(lab,lp, dummy);
  }  
} 

//***********************************

#ifndef V1

void AliITSclustererV2:: 
FindCluster(Int_t k,Int_t maxz,AliBin *bins,Int_t &n,Int_t *idx) {
  //------------------------------------------------------------
  // returns an array of indices of digits belonging to the cluster
  // (needed when the segmentation is not regular) 
  //------------------------------------------------------------
  if (n<200) idx[n++]=bins[k].GetIndex();
  bins[k].Use();

  if (bins[k-maxz].IsNotUsed()) FindCluster(k-maxz,maxz,bins,n,idx);
  if (bins[k-1   ].IsNotUsed()) FindCluster(k-1   ,maxz,bins,n,idx);
  if (bins[k+maxz].IsNotUsed()) FindCluster(k+maxz,maxz,bins,n,idx);
  if (bins[k+1   ].IsNotUsed()) FindCluster(k+1   ,maxz,bins,n,idx);
  /*
  if (bins[k-maxz-1].IsNotUsed()) FindCluster(k-maxz-1,maxz,bins,n,idx);
  if (bins[k-maxz+1].IsNotUsed()) FindCluster(k-maxz+1,maxz,bins,n,idx);
  if (bins[k+maxz-1].IsNotUsed()) FindCluster(k+maxz-1,maxz,bins,n,idx);
  if (bins[k+maxz+1].IsNotUsed()) FindCluster(k+maxz+1,maxz,bins,n,idx);
  */
}

void AliITSclustererV2::
FindClustersSPD(const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SPD cluster finder
  //------------------------------------------------------------
  Int_t kNzBins = fNzSPD + 2;
  const Int_t kMAXBIN=kNzBins*(fNySPD+2);

  Int_t ndigits=digits->GetEntriesFast();
  AliBin *bins=new AliBin[kMAXBIN];

  Int_t k;
  AliITSdigitSPD *d=0;
  for (k=0; k<ndigits; k++) {
     d=(AliITSdigitSPD*)digits->UncheckedAt(k);
     Int_t i=d->GetCoord2()+1;   //y
     Int_t j=d->GetCoord1()+1;
     bins[i*kNzBins+j].SetIndex(k);
     bins[i*kNzBins+j].SetMask(1);
  }
   
  Int_t n=0; TClonesArray &cl=*clusters;
  for (k=0; k<kMAXBIN; k++) {
     if (!bins[k].IsNotUsed()) continue;
     Int_t ni=0, idx[200];
     FindCluster(k,kNzBins,bins,ni,idx);
     if (ni==200) {
        Info("FindClustersSPD","Too big cluster !"); 
        continue;
     }
     Int_t milab[10];
     for (Int_t ilab=0;ilab<10;ilab++){
       milab[ilab]=-2;
     }

     d=(AliITSdigitSPD*)digits->UncheckedAt(idx[0]);
     Int_t ymin=d->GetCoord2(),ymax=ymin;
     Int_t zmin=d->GetCoord1(),zmax=zmin;

     for (Int_t l=0; l<ni; l++) {
        d=(AliITSdigitSPD*)digits->UncheckedAt(idx[l]);

        if (ymin > d->GetCoord2()) ymin=d->GetCoord2();
        if (ymax < d->GetCoord2()) ymax=d->GetCoord2();
        if (zmin > d->GetCoord1()) zmin=d->GetCoord1();
        if (zmax < d->GetCoord1()) zmax=d->GetCoord1();
	// MI addition - find all labels in cluster
	for (Int_t dlab=0;dlab<10;dlab++){
	  Int_t digitlab = (d->GetTracks())[dlab];
	  if (digitlab<0) continue;
	  AddLabel(milab,digitlab);	  
	}
	if (milab[9]>0) CheckLabels2(milab);
     }
     CheckLabels2(milab);
     //
     //Int_t idy = (fNlayer[fI]==0)? 2:3; 
     //for (Int_t iz=zmin; iz<=zmax;iz+=2)
     //Int_t idy = (ymax-ymin)/4.; // max 2 clusters
     Int_t idy = 0; // max 2 clusters
     if (fNlayer[fI]==0 &&idy<3) idy=3;
     if (fNlayer[fI]==1 &&idy<4) idy=4; 
     Int_t idz =3;
     for (Int_t iz=zmin; iz<=zmax;iz+=idz)
       for (Int_t iy=ymin; iy<=ymax;iy+=idy){
	 //
	 Int_t ndigits =0;
	 Float_t y=0.,z=0.,q=0.;	 
	 for (Int_t l=0; l<ni; l++) {
	   d=(AliITSdigitSPD*)digits->UncheckedAt(idx[l]);
	   if (zmax-zmin>=idz || ymax-ymin>=idy){
	     if (TMath::Abs( d->GetCoord2()-iy)>0.75*idy) continue;
	     if (TMath::Abs( d->GetCoord1()-iz)>0.75*idz) continue;
	   }
	   ndigits++;
	   Float_t qq=d->GetSignal();
	   y+=qq*fYSPD[d->GetCoord2()]; z+=qq*fZSPD[d->GetCoord1()]; q+=qq;   
	  
	 }     
	 if (ndigits==0) continue;
	 y/=q; z/=q;
	 y-=fHwSPD; z-=fHlSPD;
	 
	 Float_t lp[5];
	 lp[0]=-(-y+fYshift[fI]); if (fI<=fLastSPD1) lp[0]=-lp[0];
	 lp[1]=  -z+fZshift[fI];
	 // Float_t factor=TMath::Max(double(ni-3.),1.5);
	 Float_t factory=TMath::Max(ymax-ymin,1);
	 Float_t factorz=TMath::Max(zmax-zmin,1);
	 factory*= factory;
	 factorz*= factorz;	
	 //lp[2]= (fYpitchSPD*fYpitchSPD/12.)*factory;
	 //lp[3]= (fZ1pitchSPD*fZ1pitchSPD/12.)*factorz;
	 lp[2]= (fYpitchSPD*fYpitchSPD/12.);
	 lp[3]= (fZ1pitchSPD*fZ1pitchSPD/12.);
	 //lp[4]= q;
	 lp[4]= (zmax-zmin+1)*100 + (ymax-ymin+1);
	 
	 milab[3]=fNdet[fI];
	 d=(AliITSdigitSPD*)digits->UncheckedAt(idx[0]);
	 Int_t info[3] = {ymax-ymin+1,zmax-zmin+1,fNlayer[fI]};
	 new (cl[n]) AliITSclusterV2(milab,lp,info); n++; 	 
       }
  }
  
  delete [] bins;
}

void AliITSclustererV2::FindClustersSPD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SPD cluster finder for raw data
  //------------------------------------------------------------

  Int_t nClustersSPD = 0;
  Int_t kNzBins = fNzSPD + 2;
  Int_t kNyBins = fNySPD + 2;
  Int_t kMaxBin = kNzBins * kNyBins;
  AliBin *binsSPD = new AliBin[kMaxBin];
  AliBin *binsSPDInit = new AliBin[kMaxBin];
  AliBin *bins = NULL;

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsNewModule()) {
      Int_t iModule = input->GetPrevModuleID();

      // when all data from a module was read, search for clusters
      if (bins) { 
	clusters[iModule] = new TClonesArray("AliITSclusterV2");
	Int_t nClusters = 0;

	for (Int_t iBin = 0; iBin < kMaxBin; iBin++) {
	  if (bins[iBin].IsUsed()) continue;
	  Int_t nBins = 0;
	  Int_t idxBins[200];
	  FindCluster(iBin, kNzBins, bins, nBins, idxBins);
	  if (nBins == 200) {
	    Error("FindClustersSPD", "SPD: Too big cluster !\n"); 
	    continue;
	  }

	  Int_t label[4]; 
	  label[0] = -2;
	  label[1] = -2;
	  label[2] = -2;
//	  label[3] = iModule;
	  label[3] = fNdet[iModule];

	  Int_t ymin = (idxBins[0] / kNzBins) - 1;
	  Int_t ymax = ymin;
	  Int_t zmin = (idxBins[0] % kNzBins) - 1;
	  Int_t zmax = zmin;
	  for (Int_t idx = 0; idx < nBins; idx++) {
	    Int_t iy = (idxBins[idx] / kNzBins) - 1;
	    Int_t iz = (idxBins[idx] % kNzBins) - 1;
	    if (ymin > iy) ymin = iy;
	    if (ymax < iy) ymax = iy;
	    if (zmin > iz) zmin = iz;
	    if (zmax < iz) zmax = iz;
	  }

	  Int_t idy = 0; // max 2 clusters
	  if ((iModule <= fLastSPD1) &&idy<3) idy=3;
	  if ((iModule > fLastSPD1) &&idy<4) idy=4; 
	  Int_t idz =3;
	  for (Int_t iiz=zmin; iiz<=zmax;iiz+=idz)
	    for (Int_t iiy=ymin; iiy<=ymax;iiy+=idy){
	      //
	      Int_t ndigits =0;
	      Float_t y=0.,z=0.,q=0.;	 
	      for (Int_t idx = 0; idx < nBins; idx++) {
		Int_t iy = (idxBins[idx] / kNzBins) - 1;
		Int_t iz = (idxBins[idx] % kNzBins) - 1;
		if (zmax-zmin>=idz || ymax-ymin>=idy){
		  if (TMath::Abs(iy-iiy)>0.75*idy) continue;
		  if (TMath::Abs(iz-iiz)>0.75*idz) continue;
		}
		ndigits++;
		Float_t qBin = bins[idxBins[idx]].GetQ();
		y += qBin * fYSPD[iy]; 
		z += qBin * fZSPD[iz]; 
		q += qBin;   
	      }
	      if (ndigits==0) continue;
	      y /= q; 
	      z /= q;
	      y -= fHwSPD; 
	      z -= fHlSPD;

	      Float_t hit[5];  // y, z, sigma(y)^2, sigma(z)^2, charge
	      hit[0] = -(-y+fYshift[iModule]); 
	      if (iModule <= fLastSPD1) hit[0] = -hit[0];
	      hit[1] = -z+fZshift[iModule];
	      hit[2] = fYpitchSPD*fYpitchSPD/12.;
	      hit[3] = fZ1pitchSPD*fZ1pitchSPD/12.;
	      //	  hit[4] = q;
	      hit[4] = (zmax-zmin+1)*100 + (ymax-ymin+1);
	      //	  CheckLabels(label);
	      Int_t info[3]={ymax-ymin+1,zmax-zmin+1,fNlayer[iModule]};
	      new (clusters[iModule]->AddrAt(nClusters)) 
		AliITSclusterV2(label, hit,info); 
	      nClusters++;
	    }
	}

	nClustersSPD += nClusters;
	bins = NULL;
      }

      if (!next) break;
      bins = binsSPD;
      memcpy(binsSPD,binsSPDInit,sizeof(AliBin)*kMaxBin);
    }

    // fill the current digit into the bins array
    Int_t index = (input->GetCoord2()+1) * kNzBins + (input->GetCoord1()+1);
    bins[index].SetIndex(index);
    bins[index].SetMask(1);
    bins[index].SetQ(1);
  }

  delete [] binsSPDInit;
  delete [] binsSPD;

  Info("FindClustersSPD", "found clusters in ITS SPD: %d", nClustersSPD);
}


Bool_t AliITSclustererV2::IsMaximum(Int_t k,Int_t max,const AliBin *bins) {
  //------------------------------------------------------------
  //is this a local maximum ?
  //------------------------------------------------------------
  UShort_t q=bins[k].GetQ();
  if (q==1023) return kFALSE;
  if (bins[k-max].GetQ() > q) return kFALSE;
  if (bins[k-1  ].GetQ() > q) return kFALSE; 
  if (bins[k+max].GetQ() > q) return kFALSE; 
  if (bins[k+1  ].GetQ() > q) return kFALSE; 
  if (bins[k-max-1].GetQ() > q) return kFALSE;
  if (bins[k+max-1].GetQ() > q) return kFALSE; 
  if (bins[k+max+1].GetQ() > q) return kFALSE; 
  if (bins[k-max+1].GetQ() > q) return kFALSE;
  return kTRUE; 
}

void AliITSclustererV2::
FindPeaks(Int_t k,Int_t max,AliBin *b,Int_t *idx,UInt_t *msk,Int_t& n) {
  //------------------------------------------------------------
  //find local maxima
  //------------------------------------------------------------
  if (n<31)
  if (IsMaximum(k,max,b)) {
    idx[n]=k; msk[n]=(2<<n);
    n++;
  }
  b[k].SetMask(0);
  if (b[k-max].GetMask()&1) FindPeaks(k-max,max,b,idx,msk,n);
  if (b[k-1  ].GetMask()&1) FindPeaks(k-1  ,max,b,idx,msk,n);
  if (b[k+max].GetMask()&1) FindPeaks(k+max,max,b,idx,msk,n);
  if (b[k+1  ].GetMask()&1) FindPeaks(k+1  ,max,b,idx,msk,n);
}

void AliITSclustererV2::
MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m) {
  //------------------------------------------------------------
  //mark this peak
  //------------------------------------------------------------
  UShort_t q=bins[k].GetQ();

  bins[k].SetMask(bins[k].GetMask()|m); 

  if (bins[k-max].GetQ() <= q)
     if ((bins[k-max].GetMask()&m) == 0) MarkPeak(k-max,max,bins,m);
  if (bins[k-1  ].GetQ() <= q)
     if ((bins[k-1  ].GetMask()&m) == 0) MarkPeak(k-1  ,max,bins,m);
  if (bins[k+max].GetQ() <= q)
     if ((bins[k+max].GetMask()&m) == 0) MarkPeak(k+max,max,bins,m);
  if (bins[k+1  ].GetQ() <= q)
     if ((bins[k+1  ].GetMask()&m) == 0) MarkPeak(k+1  ,max,bins,m);
}

void AliITSclustererV2::
MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,AliITSclusterV2 &c) {
  //------------------------------------------------------------
  //make cluster using digits of this peak
  //------------------------------------------------------------
  Float_t q=(Float_t)bins[k].GetQ();
  Int_t i=k/max, j=k-i*max;

  c.SetQ(c.GetQ()+q);
  c.SetY(c.GetY()+i*q); 
  c.SetZ(c.GetZ()+j*q); 
  c.SetSigmaY2(c.GetSigmaY2()+i*i*q);
  c.SetSigmaZ2(c.GetSigmaZ2()+j*j*q);

  bins[k].SetMask(0xFFFFFFFE);
  
  if (bins[k-max].GetMask() == m) MakeCluster(k-max,max,bins,m,c);
  if (bins[k-1  ].GetMask() == m) MakeCluster(k-1  ,max,bins,m,c);
  if (bins[k+max].GetMask() == m) MakeCluster(k+max,max,bins,m,c);
  if (bins[k+1  ].GetMask() == m) MakeCluster(k+1  ,max,bins,m,c);
}

void AliITSclustererV2::
FindClustersSDD(AliBin* bins[2], Int_t nMaxBin, Int_t nzBins, 
		const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------
  Int_t ncl=0; TClonesArray &cl=*clusters;
  for (Int_t s=0; s<2; s++)
    for (Int_t i=0; i<nMaxBin; i++) {
      if (bins[s][i].IsUsed()) continue;
      Int_t idx[32]; UInt_t msk[32]; Int_t npeaks=0;
      FindPeaks(i, nzBins, bins[s], idx, msk, npeaks);

      if (npeaks>30) continue;
      if (npeaks==0) continue;

      Int_t k,l;
      for (k=0; k<npeaks-1; k++){//mark adjacent peaks
        if (idx[k] < 0) continue; //this peak is already removed
        for (l=k+1; l<npeaks; l++) {
           if (idx[l] < 0) continue; //this peak is already removed
           Int_t ki=idx[k]/nzBins, kj=idx[k] - ki*nzBins;
           Int_t li=idx[l]/nzBins, lj=idx[l] - li*nzBins;
           Int_t di=TMath::Abs(ki - li);
           Int_t dj=TMath::Abs(kj - lj);
           if (di>1 || dj>1) continue;
           if (bins[s][idx[k]].GetQ() > bins[s][idx[l]].GetQ()) {
              msk[l]=msk[k];
              idx[l]*=-1;
           } else {
              msk[k]=msk[l];
              idx[k]*=-1;
              break;
           } 
        }
      }

      for (k=0; k<npeaks; k++) {
        MarkPeak(TMath::Abs(idx[k]), nzBins, bins[s], msk[k]);
      }
        
      for (k=0; k<npeaks; k++) {
         if (idx[k] < 0) continue; //removed peak
         AliITSclusterV2 c;
         MakeCluster(idx[k], nzBins, bins[s], msk[k], c);
	 //mi change
	 Int_t milab[10];
	 for (Int_t ilab=0;ilab<10;ilab++){
	   milab[ilab]=-2;
	 }
	 Int_t maxi=0,mini=0,maxj=0,minj=0;
	 //AliBin *bmax=&bins[s][idx[k]];
	 //Float_t max = TMath::Max(TMath::Abs(bmax->GetQ())/5.,3.);
	 Float_t max=3;
	 for (Int_t di=-2; di<=2;di++)
	   for (Int_t dj=-3;dj<=3;dj++){
	     Int_t index = idx[k]+di+dj*nzBins;
	     if (index<0) continue;
	     if (index>=nMaxBin) continue;
	     AliBin *b=&bins[s][index];
	     if (TMath::Abs(b->GetQ())>max){
	       if (di>maxi) maxi=di;
	       if (di<mini) mini=di;
	       if (dj>maxj) maxj=dj;
	       if (dj<minj) minj=dj;
	       //
	       if(digits) {
		 if (TMath::Abs(di)<2&&TMath::Abs(dj)<2){
		   AliITSdigitSDD* d=(AliITSdigitSDD*)digits->UncheckedAt(b->GetIndex());
		   for (Int_t itrack=0;itrack<10;itrack++){
		     Int_t track = (d->GetTracks())[itrack];
		     if (track>=0) {
		       AddLabel(milab, track); 
		     }
		   }
		 }
	       }
	     }
	   }
	 
	 /* 
	    Float_t s2 = c.GetSigmaY2()/c.GetQ() - c.GetY()*c.GetY();
	    Float_t w=par->GetPadPitchWidth(sec);
	    c.SetSigmaY2(s2);
	    if (s2 != 0.) {
	    c.SetSigmaY2(c.GetSigmaY2()*0.108);
	    if (sec<par->GetNInnerSector()) c.SetSigmaY2(c.GetSigmaY2()*2.07);
	    }	 
	    s2 = c.GetSigmaZ2()/c.GetQ() - c.GetZ()*c.GetZ();
	    w=par->GetZWidth();
	    c.SetSigmaZ2(s2);
	    
	    if (s2 != 0.) {
	    c.SetSigmaZ2(c.GetSigmaZ2()*0.169);
	    if (sec<par->GetNInnerSector()) c.SetSigmaZ2(c.GetSigmaZ2()*1.77);
	    }
	 */

         c.SetSigmaY2(0.0030*0.0030);
         c.SetSigmaZ2(0.0020*0.0020);
         c.SetDetectorIndex(fNdet[fI]);

         Float_t y=c.GetY(),z=c.GetZ(), q=c.GetQ();
         y/=q; z/=q;
	 //
	 //Float_t s2 = c.GetSigmaY2()/c.GetQ() - y*y;
	 // c.SetSigmaY2(s2);
	 //s2 = c.GetSigmaZ2()/c.GetQ() - z*z;
         //c.SetSigmaZ2(s2);
	 //
         y=(y-0.5)*fYpitchSDD;
         y-=fHwSDD;
         y-=fYoffSDD;  //delay ?
         if (s) y=-y;

         z=(z-0.5)*fZpitchSDD;
         z-=fHlSDD;

         y=-(-y+fYshift[fI]);
         z=  -z+fZshift[fI];
         c.SetY(y);
         c.SetZ(z);
	 c.SetNy(maxj-minj+1);
	 c.SetNz(maxi-mini+1);
	 c.SetType(npeaks);
         c.SetQ(q/12.7);  //to be consistent with the SSD charges

         if (c.GetQ() < 20.) continue; //noise cluster
	 
	 if (digits) {	  
	   //	   AliBin *b=&bins[s][idx[k]];
	   //	   AliITSdigitSDD* d=(AliITSdigitSDD*)digits->UncheckedAt(b->GetIndex());
	   {
	     //Int_t lab[3];
	     //lab[0]=(d->GetTracks())[0];
	     //lab[1]=(d->GetTracks())[1];
	     //lab[2]=(d->GetTracks())[2];
	     //CheckLabels(lab);
	     CheckLabels2(milab); 
	     c.SetLabel(milab[0],0);
	     c.SetLabel(milab[1],1);
	     c.SetLabel(milab[2],2);
	     c.SetLayer(fNlayer[fI]);
	   }
	 }
         new (cl[ncl]) AliITSclusterV2(c); ncl++;
      }
    }
}

void AliITSclustererV2::
FindClustersSDD(const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------
  Int_t kNzBins = fNzSDD + 2;
  const Int_t kMAXBIN=kNzBins*(fNySDD+2);

  AliBin *bins[2];
  bins[0]=new AliBin[kMAXBIN];
  bins[1]=new AliBin[kMAXBIN];

  AliITSdigitSDD *d=0;
  Int_t i, ndigits=digits->GetEntriesFast();
  for (i=0; i<ndigits; i++) {
     d=(AliITSdigitSDD*)digits->UncheckedAt(i);
     Int_t y=d->GetCoord2()+1;   //y
     Int_t z=d->GetCoord1()+1;   //z
     Int_t q=d->GetSignal();
     if (q<3) continue;

     if (z <= fNzSDD) {
       bins[0][y*kNzBins+z].SetQ(q);
       bins[0][y*kNzBins+z].SetMask(1);
       bins[0][y*kNzBins+z].SetIndex(i);
     } else {
       z-=fNzSDD; 
       bins[1][y*kNzBins+z].SetQ(q);
       bins[1][y*kNzBins+z].SetMask(1);
       bins[1][y*kNzBins+z].SetIndex(i);
     }
  }
  
  FindClustersSDD(bins, kMAXBIN, kNzBins, digits, clusters);

  delete[] bins[0];
  delete[] bins[1];

}

void AliITSclustererV2::FindClustersSDD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SDD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSDD = 0;
  Int_t kNzBins = fNzSDD + 2;
  Int_t kMaxBin = kNzBins * (fNySDD+2);
  AliBin *binsSDDInit = new AliBin[kMaxBin];
  AliBin *binsSDD1 = new AliBin[kMaxBin];
  AliBin *binsSDD2 = new AliBin[kMaxBin];
  AliBin *bins[2] = {NULL, NULL};

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsNewModule()) {
      Int_t iModule = input->GetPrevModuleID();

      // when all data from a module was read, search for clusters
      if (bins[0]) { 
	clusters[iModule] = new TClonesArray("AliITSclusterV2");
	fI = iModule;
	FindClustersSDD(bins, kMaxBin, kNzBins, NULL, clusters[iModule]);
	Int_t nClusters = clusters[iModule]->GetEntriesFast();
	nClustersSDD += nClusters;
	bins[0] = bins[1] = NULL;
      }

      if (!next) break;
      bins[0]=binsSDD1;
      bins[1]=binsSDD2;
      memcpy(binsSDD1,binsSDDInit,sizeof(AliBin)*kMaxBin);
      memcpy(binsSDD2,binsSDDInit,sizeof(AliBin)*kMaxBin);
    }

    // fill the current digit into the bins array
    if(input->GetSignal()>=3) {
      Int_t iz = input->GetCoord1()+1;
      Int_t side = ((iz <= fNzSDD) ? 0 : 1);
      iz -= side*fNzSDD;
      Int_t index = (input->GetCoord2()+1) * kNzBins + iz;
      bins[side][index].SetQ(input->GetSignal());
      bins[side][index].SetMask(1);
      bins[side][index].SetIndex(index);
    }
  }

  delete[] binsSDD1;
  delete[] binsSDD2;
  delete[] binsSDDInit;

  Info("FindClustersSDD", "found clusters in ITS SDD: %d", nClustersSDD);
}



void AliITSclustererV2::
FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		Ali1Dcluster* pos, Int_t np,
		TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SSD cluster finder
  //------------------------------------------------------------
  TClonesArray &cl=*clusters;
  //
  Float_t tanp=fTanP, tann=fTanN;
  if (fI>fLastSSD1) {tann=fTanP; tanp=fTanN;}
  Int_t idet=fNdet[fI];
  Int_t ncl=0;
  //
  Int_t negativepair[30000];
  Int_t cnegative[3000];  
  Int_t cused1[3000];
  Int_t positivepair[30000];
  Int_t cpositive[3000];
  Int_t cused2[3000];
  for (Int_t i=0;i<3000;i++) {cnegative[i]=0; cused1[i]=0;}
  for (Int_t i=0;i<3000;i++) {cpositive[i]=0; cused2[i]=0;}
  for (Int_t i=0;i<30000;i++){negativepair[i]=positivepair[i]=0;}
  static Short_t pairs[1000][1000];
  memset(pairs,0,sizeof(Short_t)*1000000);
//   Short_t ** pairs = new Short_t*[1000];
//   for (Int_t i=0; i<1000; i++) {
//     pairs[i] = new Short_t[1000];
//     memset(pairs[i],0,sizeof(Short_t)*1000);
//   }  
  //
  // find available pairs
  //
  for (Int_t i=0; i<np; i++) {
    Float_t yp=pos[i].GetY()*fYpitchSSD; 
    if (pos[i].GetQ()<3) continue;
    for (Int_t j=0; j<nn; j++) {
      if (neg[j].GetQ()<3) continue;
      Float_t yn=neg[j].GetY()*fYpitchSSD;
      Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
      Float_t yt=yn + tann*zt;
      zt-=fHlSSD; yt-=fHwSSD;
      if (TMath::Abs(yt)<fHwSSD+0.01)
      if (TMath::Abs(zt)<fHlSSD+0.01*(neg[j].GetNd()+pos[i].GetNd())) {
	negativepair[i*10+cnegative[i]] =j;  //index
	positivepair[j*10+cpositive[j]] =i;
	cnegative[i]++;  //counters
	cpositive[j]++;	
	pairs[i][j]=100;
      }
    }
  }
  //
  for (Int_t i=0; i<np; i++) {
    Float_t yp=pos[i].GetY()*fYpitchSSD; 
    if (pos[i].GetQ()<3) continue;
    for (Int_t j=0; j<nn; j++) {
      if (neg[j].GetQ()<3) continue;
      if (cpositive[j]&&cnegative[i]) continue;
      Float_t yn=neg[j].GetY()*fYpitchSSD;
      Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
      Float_t yt=yn + tann*zt;
      zt-=fHlSSD; yt-=fHwSSD;
      if (TMath::Abs(yt)<fHwSSD+0.1)
      if (TMath::Abs(zt)<fHlSSD+0.15) {
	if (cnegative[i]==0) pos[i].SetNd(100);  // not available pair
	if (cpositive[j]==0) neg[j].SetNd(100);  // not available pair
	negativepair[i*10+cnegative[i]] =j;  //index
	positivepair[j*10+cpositive[j]] =i;
	cnegative[i]++;  //counters
	cpositive[j]++;	
	pairs[i][j]=100;
      }
    }
  }
  //
  Float_t lp[5];
  Int_t milab[10];
  Double_t ratio;
  
  //
  // sign gold tracks
  //
  for (Int_t ip=0;ip<np;ip++){
    Float_t ybest=1000,zbest=1000,qbest=0;
    //
    // select gold clusters
    if ( (cnegative[ip]==1) && cpositive[negativepair[10*ip]]==1){ 
      Float_t yp=pos[ip].GetY()*fYpitchSSD; 
      Int_t j = negativepair[10*ip];      
      ratio = (pos[ip].GetQ()-neg[j].GetQ())/(pos[ip].GetQ()+neg[j].GetQ());
      //
      Float_t yn=neg[j].GetY()*fYpitchSSD;
      Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
      Float_t yt=yn + tann*zt;
      zt-=fHlSSD; yt-=fHwSSD;
      ybest=yt; zbest=zt; 
      qbest=0.5*(pos[ip].GetQ()+neg[j].GetQ());
      lp[0]=-(-ybest+fYshift[fI]);
      lp[1]=  -zbest+fZshift[fI];
      lp[2]=0.0025*0.0025;  //SigmaY2
      lp[3]=0.110*0.110;  //SigmaZ2
      
      lp[4]=qbest;        //Q
      for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
      for (Int_t ilab=0;ilab<3;ilab++){
	milab[ilab] = pos[ip].GetLabel(ilab);
	milab[ilab+3] = neg[j].GetLabel(ilab);
      }
      //
      CheckLabels2(milab);
      milab[3]=(((ip<<10) + j)<<10) + idet; // pos|neg|det
      Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fI]};
      AliITSclusterV2 * cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
      ncl++;
      cl2->SetChargeRatio(ratio);    	
      cl2->SetType(1);
      pairs[ip][j]=1;
      if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	cl2->SetType(2);
	pairs[ip][j]=2;
      }
      cused1[ip]++;
      cused2[j]++;
    }
  }
    
  for (Int_t ip=0;ip<np;ip++){
    Float_t ybest=1000,zbest=1000,qbest=0;
    //
    //
    // select "silber" cluster
    if ( cnegative[ip]==1 && cpositive[negativepair[10*ip]]==2){
      Int_t in  = negativepair[10*ip];
      Int_t ip2 = positivepair[10*in];
      if (ip2==ip) ip2 =  positivepair[10*in+1];
      Float_t pcharge = pos[ip].GetQ()+pos[ip2].GetQ();
      if (TMath::Abs(pcharge-neg[in].GetQ())<10){
	//
	// add first pair
	if (pairs[ip][in]==100){  //
	  Float_t yp=pos[ip].GetY()*fYpitchSSD; 
	  Float_t yn=neg[in].GetY()*fYpitchSSD;
	  Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
	  Float_t yt=yn + tann*zt;
	  zt-=fHlSSD; yt-=fHwSSD;
	  ybest =yt;  zbest=zt; 
	  qbest =pos[ip].GetQ();
	  lp[0]=-(-ybest+fYshift[fI]);
	  lp[1]=  -zbest+fZshift[fI];
	  lp[2]=0.0025*0.0025;  //SigmaY2
	  lp[3]=0.110*0.110;  //SigmaZ2
	  
	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip].GetLabel(ilab);
	    milab[ilab+3] = neg[in].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip].GetQ()-neg[in].GetQ())/(pos[ip].GetQ()+neg[in].GetQ());
	  milab[3]=(((ip<<10) + in)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip].GetNd(),neg[in].GetNd(),fNlayer[fI]};
	  AliITSclusterV2 * cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
	  ncl++;
	  cl2->SetChargeRatio(ratio);    	
	  cl2->SetType(5);
	  pairs[ip][in] = 5;
	  if ((pos[ip].GetNd()+neg[in].GetNd())>6){ //multi cluster
	    cl2->SetType(6);
	    pairs[ip][in] = 6;
	  }
	}
	//
	// add second pair
	
	//	if (!(cused1[ip2] || cused2[in])){  //
	if (pairs[ip2][in]==100){
	  Float_t yp=pos[ip2].GetY()*fYpitchSSD;
	  Float_t yn=neg[in].GetY()*fYpitchSSD;
	  Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
	  Float_t yt=yn + tann*zt;
	  zt-=fHlSSD; yt-=fHwSSD;
	  ybest =yt;  zbest=zt; 
	  qbest =pos[ip2].GetQ();
	  lp[0]=-(-ybest+fYshift[fI]);
	  lp[1]=  -zbest+fZshift[fI];
	  lp[2]=0.0025*0.0025;  //SigmaY2
	  lp[3]=0.110*0.110;  //SigmaZ2
	  
	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip2].GetLabel(ilab);
	    milab[ilab+3] = neg[in].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip2].GetQ()-neg[in].GetQ())/(pos[ip2].GetQ()+neg[in].GetQ());
	  milab[3]=(((ip2<<10) + in)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip2].GetNd(),neg[in].GetNd(),fNlayer[fI]};
	  AliITSclusterV2 *cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
	  ncl++;
	  cl2->SetChargeRatio(ratio);    	
	  cl2->SetType(5);
	  pairs[ip2][in] =5;
	  if ((pos[ip2].GetNd()+neg[in].GetNd())>6){ //multi cluster
	    cl2->SetType(6);
	    pairs[ip2][in] =6;
	  }
	}	
	cused1[ip]++;
	cused1[ip2]++;
	cused2[in]++;
      }
    }    
  }
  
  //  
  for (Int_t jn=0;jn<nn;jn++){
    if (cused2[jn]) continue;
    Float_t ybest=1000,zbest=1000,qbest=0;
    // select "silber" cluster
    if ( cpositive[jn]==1 && cnegative[positivepair[10*jn]]==2){
      Int_t ip  = positivepair[10*jn];
      Int_t jn2 = negativepair[10*ip];
      if (jn2==jn) jn2 =  negativepair[10*ip+1];
      Float_t pcharge = neg[jn].GetQ()+neg[jn2].GetQ();
      //
      if (TMath::Abs(pcharge-pos[ip].GetQ())<10){
	//
	// add first pair
	//	if (!(cused1[ip]||cused2[jn])){
	if (pairs[ip][jn]==100){
	  Float_t yn=neg[jn].GetY()*fYpitchSSD; 
	  Float_t yp=pos[ip].GetY()*fYpitchSSD;
	  Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
	  Float_t yt=yn + tann*zt;
	  zt-=fHlSSD; yt-=fHwSSD;
	  ybest =yt;  zbest=zt; 
	  qbest =neg[jn].GetQ();
	  lp[0]=-(-ybest+fYshift[fI]);
	  lp[1]=  -zbest+fZshift[fI];
	  lp[2]=0.0025*0.0025;  //SigmaY2
	  lp[3]=0.110*0.110;  //SigmaZ2
	  
	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip].GetLabel(ilab);
	    milab[ilab+3] = neg[jn].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip].GetQ()-neg[jn].GetQ())/(pos[ip].GetQ()+neg[jn].GetQ());
	  milab[3]=(((ip<<10) + jn)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn].GetNd(),fNlayer[fI]};
	  AliITSclusterV2 * cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
	  ncl++;
	  cl2->SetChargeRatio(ratio);    	
	  cl2->SetType(7);
	  pairs[ip][jn] =7;
	  if ((pos[ip].GetNd()+neg[jn].GetNd())>6){ //multi cluster
	    cl2->SetType(8);
	    pairs[ip][jn]=8;
	  }
	}
	//
	// add second pair
	//	if (!(cused1[ip]||cused2[jn2])){
	if (pairs[ip][jn2]==100){
	  Float_t yn=neg[jn2].GetY()*fYpitchSSD; 
	  Double_t yp=pos[ip].GetY()*fYpitchSSD; 
	  Double_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
	  Double_t yt=yn + tann*zt;
	  zt-=fHlSSD; yt-=fHwSSD;
	  ybest =yt;  zbest=zt; 
	  qbest =neg[jn2].GetQ();
	  lp[0]=-(-ybest+fYshift[fI]);
	  lp[1]=  -zbest+fZshift[fI];
	  lp[2]=0.0025*0.0025;  //SigmaY2
	  lp[3]=0.110*0.110;  //SigmaZ2
	  
	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip].GetLabel(ilab);
	    milab[ilab+3] = neg[jn2].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip].GetQ()-neg[jn2].GetQ())/(pos[ip].GetQ()+neg[jn2].GetQ());
	  milab[3]=(((ip<<10) + jn2)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn2].GetNd(),fNlayer[fI]};
	  AliITSclusterV2* cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
	  ncl++;
	  cl2->SetChargeRatio(ratio);    	
	  pairs[ip][jn2]=7;
	  cl2->SetType(7);
	  if ((pos[ip].GetNd()+neg[jn2].GetNd())>6){ //multi cluster
	    cl2->SetType(8);
	    pairs[ip][jn2]=8;
	  }
	}
	cused1[ip]++;
	cused2[jn]++;
	cused2[jn2]++;
      }
    }    
  }
  
  for (Int_t ip=0;ip<np;ip++){
    Float_t ybest=1000,zbest=1000,qbest=0;
    //
    // 2x2 clusters
    //
    if ( (cnegative[ip]<5) && cpositive[negativepair[10*ip]]<5){ 
      Float_t minchargediff =4.;
      Int_t j=-1;
      for (Int_t di=0;di<cnegative[ip];di++){
	Int_t   jc = negativepair[ip*10+di];
	Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	if (TMath::Abs(chargedif)<minchargediff){
	  j =jc;
	  minchargediff = TMath::Abs(chargedif);
	}
      }
      if (j<0) continue;  // not proper cluster      
      Int_t count =0;
      for (Int_t di=0;di<cnegative[ip];di++){
	Int_t   jc = negativepair[ip*10+di];
	Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	if (TMath::Abs(chargedif)<minchargediff+3.) count++;
      }
      if (count>1) continue;  // more than one "proper" cluster for positive
      //
      count =0;
      for (Int_t dj=0;dj<cpositive[j];dj++){
	Int_t   ic  = positivepair[j*10+dj];
	Float_t chargedif = pos[ic].GetQ()-neg[j].GetQ();
	if (TMath::Abs(chargedif)<minchargediff+3.) count++;
      }
      if (count>1) continue;  // more than one "proper" cluster for negative
      
      Int_t jp = 0;
      
      count =0;
      for (Int_t dj=0;dj<cnegative[jp];dj++){
	Int_t   ic = positivepair[jp*10+dj];
	Float_t chargedif = pos[ic].GetQ()-neg[jp].GetQ();
	if (TMath::Abs(chargedif)<minchargediff+4.) count++;
      }
      if (count>1) continue;   
      if (pairs[ip][j]<100) continue;
      //
      //almost gold clusters
      Float_t yp=pos[ip].GetY()*fYpitchSSD; 
      Float_t yn=neg[j].GetY()*fYpitchSSD;
      Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
      Float_t yt=yn + tann*zt;
      zt-=fHlSSD; yt-=fHwSSD;
      ybest=yt; zbest=zt; 
      qbest=0.5*(pos[ip].GetQ()+neg[j].GetQ());
      lp[0]=-(-ybest+fYshift[fI]);
      lp[1]=  -zbest+fZshift[fI];
      lp[2]=0.0025*0.0025;  //SigmaY2
      lp[3]=0.110*0.110;  //SigmaZ2	
      lp[4]=qbest;        //Q
      for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
      for (Int_t ilab=0;ilab<3;ilab++){
	milab[ilab] = pos[ip].GetLabel(ilab);
	milab[ilab+3] = neg[j].GetLabel(ilab);
      }
      //
      CheckLabels2(milab);
      ratio = (pos[ip].GetQ()-neg[j].GetQ())/(pos[ip].GetQ()+neg[j].GetQ());
      milab[3]=(((ip<<10) + j)<<10) + idet; // pos|neg|det
      Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fI]};
      AliITSclusterV2 * cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info);
      ncl++;
      cl2->SetChargeRatio(ratio);    	
      cl2->SetType(10);
      pairs[ip][j]=10;
      if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	cl2->SetType(11);
	pairs[ip][j]=11;
      }
      cused1[ip]++;
      cused2[j]++;      
    }

  }
  
  //  
  for (Int_t i=0; i<np; i++) {
    Float_t ybest=1000,zbest=1000,qbest=0;
    Float_t yp=pos[i].GetY()*fYpitchSSD; 
    if (pos[i].GetQ()<3) continue;
    for (Int_t j=0; j<nn; j++) {
    //    for (Int_t di = 0;di<cpositive[i];di++){
    //  Int_t j = negativepair[10*i+di];
      if (neg[j].GetQ()<3) continue;
      if (cused2[j]||cused1[i]) continue;      
      if (pairs[i][j]>0 &&pairs[i][j]<100) continue;
      ratio = (pos[i].GetQ()-neg[j].GetQ())/(pos[i].GetQ()+neg[j].GetQ());      
      Float_t yn=neg[j].GetY()*fYpitchSSD;
      Float_t zt=(2*fHlSSD*tanp + yp - yn)/(tann+tanp);
      Float_t yt=yn + tann*zt;
      zt-=fHlSSD; yt-=fHwSSD;
      if (TMath::Abs(yt)<fHwSSD+0.01)
      if (TMath::Abs(zt)<fHlSSD+0.01*(neg[j].GetNd()+pos[i].GetNd())) {
        ybest=yt; zbest=zt; 
        qbest=0.5*(pos[i].GetQ()+neg[j].GetQ());
        lp[0]=-(-ybest+fYshift[fI]);
        lp[1]=  -zbest+fZshift[fI];
        lp[2]=0.0025*0.0025;  //SigmaY2
        lp[3]=0.110*0.110;  //SigmaZ2

        lp[4]=qbest;        //Q
	for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	for (Int_t ilab=0;ilab<3;ilab++){
	  milab[ilab] = pos[i].GetLabel(ilab);
	  milab[ilab+3] = neg[j].GetLabel(ilab);
	}
	//
	CheckLabels2(milab);
	milab[3]=(((i<<10) + j)<<10) + idet; // pos|neg|det
	Int_t info[3] = {pos[i].GetNd(),neg[j].GetNd(),fNlayer[fI]};
        AliITSclusterV2 * cl2 = new (cl[ncl]) AliITSclusterV2(milab,lp,info); 
	ncl++;
	cl2->SetChargeRatio(ratio);
	cl2->SetType(100+cpositive[j]+cnegative[i]);
	//cl2->SetType(0);
	/*
	  if (pairs[i][j]<100){
	  printf("problem:- %d\n", pairs[i][j]);
	  }
	  if (cnegative[i]<2&&cpositive[j]<2){
	  printf("problem:- %d\n", pairs[i][j]);
	  }
	*/
      }
    }
  }

//   for (Int_t i=0; i<1000; i++) delete [] pairs[i];
//   delete [] pairs;

}


void AliITSclustererV2::
FindClustersSSD(const TClonesArray *alldigits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SSD cluster finder
  //------------------------------------------------------------
  Int_t smaxall=alldigits->GetEntriesFast();
  if (smaxall==0) return;
  TObjArray *digits = new TObjArray;
  for (Int_t i=0;i<smaxall; i++){
    AliITSdigitSSD *d=(AliITSdigitSSD*)alldigits->UncheckedAt(i);
    if (d->GetSignal()<3) continue;
    digits->AddLast(d);
  }
  Int_t smax = digits->GetEntriesFast();
  if (smax==0) return;
  
  const Int_t MAX=1000;
  Int_t np=0, nn=0; 
  Ali1Dcluster pos[MAX], neg[MAX];
  Float_t y=0., q=0., qmax=0.; 
  Int_t lab[4]={-2,-2,-2,-2};
  
  AliITSdigitSSD *d=(AliITSdigitSSD*)digits->UncheckedAt(0);
  q += d->GetSignal();
  y += d->GetCoord2()*d->GetSignal();
  qmax=d->GetSignal();
  lab[0]=d->GetTrack(0); lab[1]=d->GetTrack(1); lab[2]=d->GetTrack(2);
  Int_t curr=d->GetCoord2();
  Int_t flag=d->GetCoord1();
  Int_t *n=&nn;
  Ali1Dcluster *c=neg;
  Int_t nd=1;
  Int_t milab[10];
  for (Int_t ilab=0;ilab<10;ilab++){
    milab[ilab]=-2;
  }
  milab[0]=d->GetTrack(0); milab[1]=d->GetTrack(1); milab[2]=d->GetTrack(2);

  for (Int_t s=1; s<smax; s++) {
      d=(AliITSdigitSSD*)digits->UncheckedAt(s);      
      Int_t strip=d->GetCoord2();
      if ((strip-curr) > 1 || flag!=d->GetCoord1()) {
         c[*n].SetY(y/q);
         c[*n].SetQ(q);
         c[*n].SetNd(nd);
	 CheckLabels2(milab);
         c[*n].SetLabels(milab);
         //Split suspiciously big cluster
	 /*
	 if (nd>10&&nd<16){
	   c[*n].SetY(y/q-0.3*nd);
	   c[*n].SetQ(0.5*q);
	   (*n)++;
	   if (*n==MAX) {
	     Error("FindClustersSSD","Too many 1D clusters !");
              return;
	   }
	   c[*n].SetY(y/q-0.0*nd);
	   c[*n].SetQ(0.5*q);
	   c[*n].SetNd(nd);
	   (*n)++;
	   if (*n==MAX) {
	     Error("FindClustersSSD","Too many 1D clusters !");
              return;
	   }
	   //
	   c[*n].SetY(y/q+0.3*nd);
	   c[*n].SetQ(0.5*q);
	   c[*n].SetNd(nd);
	   c[*n].SetLabels(milab);
	 }
	 else{
	 */
	 if (nd>4&&nd<25) {
	   c[*n].SetY(y/q-0.25*nd);
	   c[*n].SetQ(0.5*q);
	   (*n)++;
	   if (*n==MAX) {
	     Error("FindClustersSSD","Too many 1D clusters !");
	     return;
	   }
	   c[*n].SetY(y/q+0.25*nd);
	   c[*n].SetQ(0.5*q);
	   c[*n].SetNd(nd);
	   c[*n].SetLabels(milab);
	 }	 
         (*n)++;
         if (*n==MAX) {
          Error("FindClustersSSD","Too many 1D clusters !");
          return;
         }
         y=q=qmax=0.;
         nd=0;
         lab[0]=lab[1]=lab[2]=-2;
	 //
	 for (Int_t ilab=0;ilab<10;ilab++){
	   milab[ilab]=-2;
	 }
	 //
         if (flag!=d->GetCoord1()) { n=&np; c=pos; }
      }
      flag=d->GetCoord1();
      q += d->GetSignal();
      y += d->GetCoord2()*d->GetSignal();
      nd++;
      if (d->GetSignal()>qmax) {
         qmax=d->GetSignal();
         lab[0]=d->GetTrack(0); lab[1]=d->GetTrack(1); lab[2]=d->GetTrack(2);
      }
      for (Int_t ilab=0;ilab<10;ilab++) {
	if (d->GetTrack(ilab)>=0) AddLabel(milab, (d->GetTrack(ilab))); 
      }
      curr=strip;
  }
  c[*n].SetY(y/q);
  c[*n].SetQ(q);
  c[*n].SetNd(nd);
  c[*n].SetLabels(lab);
  //Split suspiciously big cluster
  if (nd>4 && nd<25) {
     c[*n].SetY(y/q-0.25*nd);
     c[*n].SetQ(0.5*q);
     (*n)++;
     if (*n==MAX) {
        Error("FindClustersSSD","Too many 1D clusters !");
        return;
     }
     c[*n].SetY(y/q+0.25*nd);
     c[*n].SetQ(0.5*q);
     c[*n].SetNd(nd);
     c[*n].SetLabels(lab);
  }
  (*n)++;
  if (*n==MAX) {
     Error("FindClustersSSD","Too many 1D clusters !");
     return;
  }

  FindClustersSSD(neg, nn, pos, np, clusters);
}

void AliITSclustererV2::FindClustersSSD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SSD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSSD = 0;
  const Int_t MAX = 1000;
  Ali1Dcluster clusters1D[2][MAX];
  Int_t nClusters[2] = {0, 0};
  Int_t lab[3]={-2,-2,-2};
  Float_t q = 0.;
  Float_t y = 0.;
  Int_t nDigits = 0;
  Int_t prevStrip = -1;
  Int_t prevFlag = -1;
  Int_t prevModule = -1;

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();

    if(input->GetSignal()<3 && next) continue;
    // check if a new cluster starts
    Int_t strip = input->GetCoord2();
    Int_t flag = input->GetCoord1();
    if ((!next || (input->GetModuleID() != prevModule)||
	 (strip-prevStrip > 1) || (flag != prevFlag)) &&
	(nDigits > 0)) {
      if (nClusters[prevFlag] == MAX) {
	Error("FindClustersSSD", "Too many 1D clusters !");
	return;
      }
      Ali1Dcluster& cluster = clusters1D[prevFlag][nClusters[prevFlag]++];
      cluster.SetY(y/q);
      cluster.SetQ(q);
      cluster.SetNd(nDigits);
      cluster.SetLabels(lab);

      //Split suspiciously big cluster
      if (nDigits > 4&&nDigits < 25) {
	cluster.SetY(y/q - 0.25*nDigits);
        cluster.SetQ(0.5*q);
	if (nClusters[prevFlag] == MAX) {
	  Error("FindClustersSSD", "Too many 1D clusters !");
	  return;
	}
	Ali1Dcluster& cluster2 = clusters1D[prevFlag][nClusters[prevFlag]++];
	cluster2.SetY(y/q + 0.25*nDigits);
	cluster2.SetQ(0.5*q);
	cluster2.SetNd(nDigits);
	cluster2.SetLabels(lab);
      }
      y = q = 0.;
      nDigits = 0;
    }

    if (!next || (input->GetModuleID() != prevModule)) {
      Int_t iModule = prevModule;

      // when all data from a module was read, search for clusters
      if (prevFlag >= 0) {
	clusters[iModule] = new TClonesArray("AliITSclusterV2");
	fI = iModule;
	FindClustersSSD(&clusters1D[0][0], nClusters[0], 
			&clusters1D[1][0], nClusters[1], clusters[iModule]);
	Int_t nClusters = clusters[iModule]->GetEntriesFast();
	nClustersSSD += nClusters;
      }

      if (!next) break;
      nClusters[0] = nClusters[1] = 0;
      y = q = 0.;
      nDigits = 0;
    }

    // add digit to current cluster
    q += input->GetSignal();
    y += strip * input->GetSignal();
    nDigits++;
    prevStrip = strip;
    prevFlag = flag;
    prevModule = input->GetModuleID();

  }

  Info("FindClustersSSD", "found clusters in ITS SSD: %d", nClustersSSD);
}

#else   //V1

#include "AliITSDetType.h"
#include "AliITS.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSClusterFinderSPD.h"

#include "AliITSresponseSDD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSClusterFinderSDD.h"

#include "AliITSsegmentationSSD.h"
#include "AliITSClusterFinderSSD.h"


void AliITSclustererV2::
FindClustersSPD(const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SPD cluster finding based on AliITSClusterFinderSPD
  //------------------------------------------------------------
  static AliITS *its=(AliITS*)gAlice->GetModule("ITS");
  static TClonesArray *points=its->RecPoints();
  static AliITSsegmentationSPD *seg=
         (AliITSsegmentationSPD *)its->DetType(0)->GetSegmentationModel();
  static AliITSClusterFinderSPD cf(seg, (TClonesArray*)digits, points);

  cf.FindRawClusters(fI);
  RecPoints2Clusters(points, fI, clusters);
  its->ResetRecPoints();

}

void AliITSclustererV2::
FindClustersSDD(const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SDD cluster finding based on AliITSClusterFinderSDD
  //------------------------------------------------------------
  static AliITS *its=(AliITS*)gAlice->GetModule("ITS");
  static TClonesArray *points=its->RecPoints();
  static AliITSresponseSDD *resp=
        (AliITSresponseSDD *)its->DetType(1)->GetResponseModel();
  static AliITSsegmentationSDD *seg=
         (AliITSsegmentationSDD *)its->DetType(1)->GetSegmentationModel();
  static AliITSClusterFinderSDD 
         cf(seg,resp,(TClonesArray*)digits,its->ClustersAddress(1));

  cf.FindRawClusters(fI);
  Int_t nc=points->GetEntriesFast();
  while (nc--) { //To be consistent with the SSD cluster charges
     AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(nc);
     p->SetQ(p->GetQ()/12.);
  }
  RecPoints2Clusters(points, fI, clusters);
  its->ResetClusters(1);
  its->ResetRecPoints();

}

void AliITSclustererV2::
FindClustersSSD(const TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SSD cluster finding based on AliITSClusterFinderSSD
  //------------------------------------------------------------
  static AliITS *its=(AliITS*)gAlice->GetModule("ITS");
  static TClonesArray *points=its->RecPoints();
  static AliITSsegmentationSSD *seg=
         (AliITSsegmentationSSD *)its->DetType(2)->GetSegmentationModel();
  static AliITSClusterFinderSSD cf(seg,(TClonesArray*)digits);

  cf.FindRawClusters(fI);
  RecPoints2Clusters(points, fI, clusters);
  its->ResetRecPoints();

}

#endif
