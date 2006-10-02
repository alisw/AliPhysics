/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
////////////////////////////////////////////////////////////////////////////
//            Implementation of the ITS clusterer V2 class                //
//                                                                        //
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch            //
//                                                                        //
///////////////////////////////////////////////////////////////////////////


#include "AliITSClusterFinderV2SSD.h"
#include "AliITSRecPoint.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSSD.h"
#include <TClonesArray.h>
#include "AliITSdigitSSD.h"

ClassImp(AliITSClusterFinderV2SSD)


AliITSClusterFinderV2SSD::AliITSClusterFinderV2SSD(AliITSDetTypeRec* dettyp):AliITSClusterFinderV2(dettyp),
fLastSSD1(0),
fYpitchSSD(0.0095),
fHwSSD(3.65),
fHlSSD(2.00),
fTanP(0.0275),
fTanN(0.0075){

  //Default constructor

  fLastSSD1=fDetTypeRec->GetITSgeom()->GetModuleIndex(6,1,1)-1;

}
 

void AliITSClusterFinderV2SSD::FindRawClusters(Int_t mod){

  //Find clusters V2
  SetModule(mod);
  FindClustersSSD(fDigits);

}

void AliITSClusterFinderV2SSD::FindClustersSSD(TClonesArray *alldigits) {
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
  
  const Int_t kMax=1000;
  Int_t np=0, nn=0; 
  Ali1Dcluster pos[kMax], neg[kMax];
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
	   if (*n==kMax) {
	     Error("FindClustersSSD","Too many 1D clusters !");
	     return;
	   }
	   c[*n].SetY(y/q+0.25*nd);
	   c[*n].SetQ(0.5*q);
	   c[*n].SetNd(nd);
	   c[*n].SetLabels(milab);
	 }	 
         (*n)++;
         if (*n==kMax) {
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
     if (*n==kMax) {
        Error("FindClustersSSD","Too many 1D clusters !");
        return;
     }
     c[*n].SetY(y/q+0.25*nd);
     c[*n].SetQ(0.5*q);
     c[*n].SetNd(nd);
     c[*n].SetLabels(lab);
  }
  (*n)++;
  if (*n==kMax) {
     Error("FindClustersSSD","Too many 1D clusters !");
     return;
  }

  FindClustersSSD(neg, nn, pos, np);
}


void AliITSClusterFinderV2SSD::RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters){

    //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStreamSSD inputSSD(rawReader);
  FindClustersSSD(&inputSSD,clusters);
  
}

void AliITSClusterFinderV2SSD::FindClustersSSD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SSD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSSD = 0;
  const Int_t kMax = 1000;
  Ali1Dcluster clusters1D[2][kMax];
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
      if (nClusters[prevFlag] == kMax) {
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
	if (nClusters[prevFlag] == kMax) {
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
	clusters[iModule] = new TClonesArray("AliITSRecPoint");
	fModule = iModule;
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

void AliITSClusterFinderV2SSD::
FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		Ali1Dcluster* pos, Int_t np,
		TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SSD cluster finder
  //------------------------------------------------------------
  TClonesArray &cl=*clusters;
  //
  Float_t tanp=fTanP, tann=fTanN;
  if (fModule>fLastSSD1) {tann=fTanP; tanp=fTanN;}
  Int_t idet=fNdet[fModule];
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
  for (Int_t i=0;i<30000;i++) {negativepair[i]=0; positivepair[i]=0;}
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
      lp[0]=-(-ybest+fYshift[fModule]);
      lp[1]=  -zbest+fZshift[fModule];
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
      Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fModule]};
      AliITSRecPoint * cl2;
      if(clusters){
	cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
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
      else{
	cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);	
	cl2->SetChargeRatio(ratio);    	
	cl2->SetType(1);
	pairs[ip][j]=1;
	if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	  cl2->SetType(2);
	  pairs[ip][j]=2;
	}
	cused1[ip]++;
	cused2[j]++;
	fDetTypeRec->AddRecPoint(*cl2);
      }
      ncl++;
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
	  lp[0]=-(-ybest+fYshift[fModule]);
	  lp[1]=  -zbest+fZshift[fModule];
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
	  Int_t info[3] = {pos[ip].GetNd(),neg[in].GetNd(),fNlayer[fModule]};

	  AliITSRecPoint * cl2;
	  if(clusters){
	    cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(5);
	    pairs[ip][in] = 5;
	    if ((pos[ip].GetNd()+neg[in].GetNd())>6){ //multi cluster
	      cl2->SetType(6);
	      pairs[ip][in] = 6;
	    }	    
	  }
	  else{
	    cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(5);
	    pairs[ip][in] = 5;
	    if ((pos[ip].GetNd()+neg[in].GetNd())>6){ //multi cluster
	      cl2->SetType(6);
	      pairs[ip][in] = 6;
	    }
	    
	    fDetTypeRec->AddRecPoint(*cl2);
	  }
	  ncl++;
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
	  lp[0]=-(-ybest+fYshift[fModule]);
	  lp[1]=  -zbest+fZshift[fModule];
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
	  Int_t info[3] = {pos[ip2].GetNd(),neg[in].GetNd(),fNlayer[fModule]};
	  
	  AliITSRecPoint * cl2;
	  if(clusters){
	    cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(5);
	    pairs[ip2][in] =5;
	    if ((pos[ip2].GetNd()+neg[in].GetNd())>6){ //multi cluster
	      cl2->SetType(6);
	      pairs[ip2][in] =6;
	    }
	  }
	  else{
	    cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(5);
	    pairs[ip2][in] =5;
	    if ((pos[ip2].GetNd()+neg[in].GetNd())>6){ //multi cluster
	      cl2->SetType(6);
	      pairs[ip2][in] =6;
	    }

	    fDetTypeRec->AddRecPoint(*cl2);
	  }
	  ncl++;
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
	  lp[0]=-(-ybest+fYshift[fModule]);
	  lp[1]=  -zbest+fZshift[fModule];
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
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn].GetNd(),fNlayer[fModule]};

	  AliITSRecPoint * cl2;
	  if(clusters){
	    cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(7);
	    pairs[ip][jn] =7;
	    if ((pos[ip].GetNd()+neg[jn].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      pairs[ip][jn]=8;
	    }

	  }
	  else{
	    cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(7);
	    pairs[ip][jn] =7;
	    if ((pos[ip].GetNd()+neg[jn].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      pairs[ip][jn]=8;
	    }

	    fDetTypeRec->AddRecPoint(*cl2);
	  }
	  ncl++;
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
	  lp[0]=-(-ybest+fYshift[fModule]);
	  lp[1]=  -zbest+fZshift[fModule];
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
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn2].GetNd(),fNlayer[fModule]};
	  AliITSRecPoint * cl2;
	  if(clusters){
	    cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    pairs[ip][jn2]=7;
	    cl2->SetType(7);
	    if ((pos[ip].GetNd()+neg[jn2].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      pairs[ip][jn2]=8;
	    }
	    
	  }
	  else{
	    cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	    cl2->SetChargeRatio(ratio);    	
	    pairs[ip][jn2]=7;
	    cl2->SetType(7);
	    if ((pos[ip].GetNd()+neg[jn2].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      pairs[ip][jn2]=8;
	    }
	    
	    fDetTypeRec->AddRecPoint(*cl2);
	  }

	  ncl++;
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
      lp[0]=-(-ybest+fYshift[fModule]);
      lp[1]=  -zbest+fZshift[fModule];
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
      Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fModule]};
      AliITSRecPoint * cl2;
      if(clusters){
	cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
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
      else{
	cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	cl2->SetChargeRatio(ratio);    	
	cl2->SetType(10);
	pairs[ip][j]=10;
	if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	  cl2->SetType(11);
	  pairs[ip][j]=11;
	}
	cused1[ip]++;
	cused2[j]++;      
	
	fDetTypeRec->AddRecPoint(*cl2);
      }      
      ncl++;
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
        lp[0]=-(-ybest+fYshift[fModule]);
        lp[1]=  -zbest+fZshift[fModule];
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
	Int_t info[3] = {pos[i].GetNd(),neg[j].GetNd(),fNlayer[fModule]};
	AliITSRecPoint * cl2;
	if(clusters){
	  cl2 = new (cl[ncl]) AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	  cl2->SetChargeRatio(ratio);
	  cl2->SetType(100+cpositive[j]+cnegative[i]);	  
	}
	else{
	  cl2 = new AliITSRecPoint(fModule,fDetTypeRec->GetITSgeom(),milab,lp,info);
	  cl2->SetChargeRatio(ratio);
	  cl2->SetType(100+cpositive[j]+cnegative[i]);
	  fDetTypeRec->AddRecPoint(*cl2);
	}
      	ncl++;
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



