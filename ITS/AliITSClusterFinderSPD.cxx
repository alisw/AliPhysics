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
////////////////////////////////////////////////////////////////////////////
//  Cluster finder                                                       ///
//  for Silicon pixels                                                    // 
//                                                                        //
////////////////////////////////////////////////////////////////////////////        

#include "AliITSClusterFinderSPD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRawClusterSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSdigitSPD.h"
#include "AliITSsegmentationSPD.h"
#include "AliLog.h"

//#define DEBUG

ClassImp(AliITSClusterFinderSPD)

//______________________________________________________________________
AliITSClusterFinderSPD::AliITSClusterFinderSPD():AliITSClusterFinder(),
fDz(0.0),
fDx(0.0),
fMinNCells(0){
    // constructor
}
//----------------------------------------------------------
AliITSClusterFinderSPD::AliITSClusterFinderSPD(AliITSDetTypeRec* dettyp):
AliITSClusterFinder(dettyp),
fDz(0.0),
fDx(0.0),
fMinNCells(0){
    // constructor

    SetDx();
    SetDz();
}
//----------------------------------------------------------
AliITSClusterFinderSPD::AliITSClusterFinderSPD(AliITSDetTypeRec* dettyp,
                                               TClonesArray *digits,
                                               TClonesArray *recp):
AliITSClusterFinder(dettyp,digits),
fDz(0.0),
fDx(0.0),
fMinNCells(0){
    // constructor

    SetDigits(digits);
    SetClusters(recp);
    SetDx();
    SetDz();
}
/*
//______________________________________________________________________
AliITSClusterFinderSPD::AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source) : AliITSClusterFinder(source) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Fatal("AliITSClusterFinderSPD","Copy constructor not allowed\n");
}
*/
//______________________________________________________________________
//AliITSClusterFinderSPD& AliITSClusterFinderSPD::operator=(const AliITSClusterFinderSPD& /* source */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  //Fatal("= operator","Assignment operator not allowed\n");
  //return *this;
//}
//______________________________________________________________________
void AliITSClusterFinderSPD::FindRawClusters(Int_t module){   
    // input of Cluster Finder
    Int_t   digitcount  = 0;
    Int_t   numberd     = 100000;
    Int_t   *digx       = new Int_t[numberd];
    Int_t   *digz       = new Int_t[numberd];
    Int_t   *digtr1     = new Int_t[numberd];
    Int_t   *digtr2     = new Int_t[numberd];
    Int_t   *digtr3     = new Int_t[numberd];
    Int_t   *digtr4     = new Int_t[numberd];
    //  output of Cluster Finder    
    Int_t   numberc     = 10000;
    Double_t *xcenterl   = new Double_t[numberc];
    Double_t *zcenterl   = new Double_t[numberc];
    Double_t *errxcenter = new Double_t[numberc];
    Double_t *errzcenter = new Double_t[numberc];
    Int_t   *tr1clus    = new Int_t[numberc];
    Int_t   *tr2clus    = new Int_t[numberc];
    Int_t   *tr3clus    = new Int_t[numberc];
    Int_t   nclus;

    SetModule(module);
    digitcount=0;
    Int_t ndigits = Digits()->GetEntriesFast();  
    if (!ndigits) return;

    AliITSdigitSPD *dig;
    Int_t ndig=0,i;
    /*
    AliDebug(4," ");
        scanf("%d",&ndig);
    */
    for(ndig=0; ndig<ndigits; ndig++) {
        dig= (AliITSdigitSPD*)GetDigit(ndig);
        digx[digitcount] = dig->GetCoord2()+1;  //starts at 1
        digz[digitcount] = dig->GetCoord1()+1;  //starts at 1
        digtr1[digitcount] = dig->GetTrack(0);
        digtr2[digitcount] = -3;
        digtr3[digitcount] = -3;
        AliDebug(5,Form("digtr1[%d]=%d fTracks[%d]=%d: ",
			digitcount,digtr1[digitcount],0,dig->GetTrack(0)));
        i=1;
        while(digtr1[digitcount]==dig->GetTrack(i) && i<dig->GetNTracks()) i++;
        AliDebug(5,Form(" fTracks[%d]=%d",i,dig->GetTrack(i)));
        if(i<dig->GetNTracks()){
            digtr2[digitcount] = dig->GetTrack(i);
            AliDebug(5,Form("digtr2[%d]=%d: ",digitcount,digtr2[digitcount]));
            while((digtr1[digitcount]==dig->GetTrack(i) || 
                   digtr2[digitcount]==dig->GetTrack(i))&&
                  i<=dig->GetNTracks()) i++;
            if(i<dig->GetNTracks()) digtr3[digitcount] = dig->GetTrack(i);
	    AliDebug(5,Form(" fTracks[%d]=%d digtr3[%d]=%d",
			    i,i<dig->GetNTracks()?dig->GetTrack(i):-1,digitcount,digtr3[digitcount]));
        } // end if
    //        if(GetDebug(4)) cout<<endl;
        digtr4[digitcount] = dig->GetSignal();
        digitcount++;
    } // end for ndig
    ClusterFinder(digitcount,digx,digz,digtr1,digtr2,digtr3,digtr4,
                  nclus,xcenterl,zcenterl,errxcenter,errzcenter,
                  tr1clus, tr2clus, tr3clus);
    DigitToPoint(nclus,xcenterl,zcenterl,errxcenter,errzcenter,
                 tr1clus, tr2clus, tr3clus);
    delete[] digx;
    delete[] digz;
    delete[] digtr1;
    delete[] digtr2;
    delete[] digtr3;
    delete[] digtr4;
    delete[] xcenterl;
    delete[] zcenterl;
    delete[] errxcenter;
    delete[] errzcenter;
    delete[] tr1clus;
    delete[] tr2clus;
    delete[] tr3clus;
}
//----------------------------------------------------------------------
void AliITSClusterFinderSPD::ClusterFinder(Int_t ndigits,Int_t digx[],
					   Int_t digz[],Int_t digtr1[],
					   Int_t digtr2[],Int_t digtr3[],
					   Int_t digtr4[],Int_t &nclus,
					   Double_t xcenter[],Double_t zcenter[],
					   Double_t errxcenter[],
					   Double_t errzcenter[],
					   Int_t tr1clus[],Int_t tr2clus[],
					   Int_t tr3clus[]){
    // Search for clusters of fired pixels (digits). Two digits are linked
    // inside a cluster if they are countiguous both in row or column
    // direction.  Diagonal digits are not linked.
    // xcenter, ycenter, zcenter are the coordinates of the center
    // of each found cluster, calculated from the averaging the corresponding
    // coordinate of the center of the linked digits. The coordinates are
    // given in the local reference sistem. 
    // errxcenter, errycenter, errzcenter are the errors associated to
    // the corresponding average.
    Int_t   if1, min, max, nd;
    Int_t   x1, z1, t1, t2, t3, t4;
    Int_t   ndx, ndz, ndxmin=0, ndxmax=0, ndzmin=0, ndzmax=0;
    Double_t dx, dz; 
    Int_t   i,k,ipos=0;
    Float_t xdum, zdum;      
    Int_t   kmax, sigmax;
    Double_t deltax, deltaz;
    Double_t ndig;
    Int_t   numberd = 10000;
    Int_t   *ifpad  = new Int_t[numberd];
    Int_t   *xpad   = new Int_t[numberd];
    Int_t   *zpad   = new Int_t[numberd];
    Int_t   *tr1pad = new Int_t[numberd];
    Int_t   *tr2pad = new Int_t[numberd];
    Int_t   *tr3pad = new Int_t[numberd];
    Int_t   *tr4pad = new Int_t[numberd];
    Int_t   *iclus  = new Int_t[numberd];

    nclus=1;
    for (i=0; i < ndigits ; i++){
        ifpad[i] = -1;
        iclus[i] = 0;
    } // end for i

    ifpad[0]=0;
    for (i=0; i < ndigits-1 ; i++) {
        if ( ifpad[i] == -1 ) { 
            nclus++;
            ipos++;
            ifpad[i]=nclus-1;
        } // end if ipad[i]
        for (Int_t j=i+1 ; j < ndigits ; j++)  {  
            if (ifpad[j]== -1 ) {
                dx = TMath::Abs(digx[i]-digx[j]);
                dz = TMath::Abs(digz[i]-digz[j]);
                // if ( ( dx+dz )==1 )  //clusters are not diagonal
                if(( dx+dz )==1 || (dx==1 && dz==1)){
                    //diagonal clusters allowed
                    ipos++;
                    ifpad[j] = ifpad[i];
                    
                    x1         = digx[j];
                    z1         = digz[j];
                    digx[j]    = digx[ipos];
                    digz[j]    = digz[ipos];
                    digx[ipos] = x1;
                    digz[ipos] = z1;
                    
                    t1 = digtr1[j];
                    t2 = digtr2[j];
                    t3 = digtr3[j];
                    t4 = digtr4[j];
                    digtr1[j] = digtr1[ipos];
                    digtr2[j] = digtr2[ipos];
                    digtr3[j] = digtr3[ipos];
                    digtr4[j] = digtr4[ipos];
                    digtr1[ipos] = t1;
                    digtr2[ipos] = t2;
                    digtr3[ipos] = t3;
                    digtr4[ipos] = t4;
                    
                    if1 = ifpad[j];
                    ifpad[j] = ifpad[ipos];
                    ifpad[ipos] = if1;
                } // end dx+dx...
            }// end if ifpad[j]== -1 
        } // end for j
    }//end loop on digits
    if ( ifpad[ndigits-1] == -1 ) {
        nclus++;
        ifpad[ndigits-1]=nclus-1;
    } // end if ifpad[ndigits-1] == -1
    
    for (i=0 ; i < ndigits ; i++) iclus[ifpad[i]]++;

    min=0;
    max=0;
    // loop on found clusters 
    for (i=0 ; i < nclus ; i++){
        min = max;
        max += iclus[i];
        deltax = GetSeg()->Dpx(0);
        if (iclus[i]!=1){
            //cluster with more than one digit
            nd=iclus[i];
            ndig=(Double_t) nd;
            Int_t count=0;
            for (k=min;k<min+nd;k++){
                xpad[count] = digx[k];	   
                zpad[count] = digz[k];

                tr1pad[count] = digtr1[k];	   
                tr2pad[count] = digtr2[k];	   
                tr3pad[count] = digtr3[k];	   
                tr4pad[count] = digtr4[k];	   
                
                count++; 
            } // end for k
            ndxmin = xpad[TMath::LocMin(nd,xpad)];
            ndxmax = xpad[TMath::LocMax(nd,xpad)];
            ndzmin = zpad[TMath::LocMin(nd,zpad)];
            ndzmax = zpad[TMath::LocMax(nd,zpad)];
            ndx = ndxmax - ndxmin+1;
            ndz = ndzmax - ndzmin+1;
            
            // calculate x and z coordinates of the center of the cluster
            GetSeg()->GetPadCxz(digx[min],digz[min]-1,xdum, zdum);
            
            if (ndx == 1) {
                xcenter[i] = xdum;
            }else{ 
                xcenter[i] = 0.;
                for (k=0;k<nd;k++) {
                    GetSeg()->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
                    xcenter[i] += (xdum / nd);
                } // end for k               
            } // end if ndx
            
            if (ndz == 1) {
                zcenter[i] = zdum;
            } else {
                zcenter[i] = 0.;
                for (k=0;k<nd;k++) {	      
                    GetSeg()->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
                    zcenter[i] += (zdum / nd);
                } // end for k
            } // end if ndz
            
            // error on points in x and z directions
            
            if (ndx == 1) {
                errxcenter[i] = deltax / TMath::Sqrt(12.);
            } else {
                errxcenter[i] = 0.;	 		
                for (k=0;k<nd;k++){ 
                    GetSeg()->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
                    errxcenter[i] += ((xdum-xcenter[i])*(xdum-xcenter[i]))/
                        (nd*(nd-1)); 
                } // end for k
                errxcenter[i] = TMath::Sqrt(errxcenter[i]);
            } // end if ndx	
            if (ndz == 1) {
                deltaz = GetSeg()->Dpz(digz[min]);	              
                errzcenter[i] = deltaz / TMath::Sqrt(12.);
            } else {
                errzcenter[i] = 0.;
                for (k=0;k<nd;k++){ 
                    GetSeg()->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
                    errzcenter[i] += ((zdum-zcenter[i])*(zdum-zcenter[i]))/
                        (nd*(nd-1));
                } // end for k
                errzcenter[i] = TMath::Sqrt(errzcenter[i]);
            } // end if ndz
            // take three track numbers for the cluster
            // choose the track numbers of the digit with higher signal 
            kmax = 0;
            sigmax = 0;
            for (k=0;k<nd;k++){
                if(tr4pad[k] > sigmax){
                    sigmax = tr4pad[k];
                    kmax   = k;
                } // end if tr4pad[k]
            } // end for k
            if(sigmax != 0) {
                tr1clus[i]= tr1pad[kmax];
                tr2clus[i]= tr2pad[kmax];
                tr3clus[i]= tr3pad[kmax];
            } else {
                tr1clus[i]= -2;
                tr2clus[i]= -2;
                tr3clus[i]= -2;
            } // end if sigmax
        } else {
            // cluster with single digit
            ndig= 1.;
            ndx = 1;
            ndz = 1;
            GetSeg()->GetPadCxz(digx[min],digz[min]-1,xdum,zdum);
            xcenter[i] = xdum;
            zcenter[i] = zdum;
            tr1clus[i]=digtr1[min];
            tr2clus[i]=digtr2[min];
            tr3clus[i]=digtr3[min];
            deltaz = GetSeg()->Dpz(digz[min]);
            errxcenter[i] = deltax / TMath::Sqrt(12.);
            errzcenter[i] = deltaz / TMath::Sqrt(12.);
        } // end if iclus[i]

        // store the cluster information to the AliITSRawCLusterSPD object
   

        //put the cluster center in local reference frame of the detector
        // and in microns
        xcenter[i] = xcenter[i] - GetSeg()->Dx()/2.; 
        zcenter[i] = zcenter[i] - GetSeg()->Dz()/2.;

        AliITSRawClusterSPD *clust = new AliITSRawClusterSPD(zcenter[i], //f
                                                             xcenter[i], //f
                                                             ndig, //f
                                                             ndz,ndx, //ii
                                                             ndxmin,ndxmax,//ii
                                                             (Double_t) ndzmin,
                                                             (Double_t) ndzmax,
                                                             0,GetModule());
	fDetTypeRec->AddCluster(0,clust);
        delete clust;
    }//end loop on clusters   
    delete[] ifpad;
    delete[] xpad ;
    delete[] zpad ;
    delete[] iclus;
    delete[] tr1pad;
    delete[] tr2pad;
    delete[] tr3pad;
    delete[] tr4pad;
}
//______________________________________________________----------------
void AliITSClusterFinderSPD::DigitToPoint(Int_t nclus,
					  Double_t *xcenter,Double_t *zcenter,
					  Double_t *errxcenter,
					  Double_t *errzcenter, 
					  Int_t *tr1clus, Int_t *tr2clus,
					  Int_t *tr3clus){
    // A point is associated to each cluster of SPD digits. The points
    // and their associated errors are stored in the file galiceSP.root.
    Double_t l[3],xg,zg;
    const Double_t kconv = 1.0e-4; // micron -> cm

    Int_t lay,lad,det;
    fDetTypeRec->GetITSgeom()->GetModuleId(fModule,lay,lad,det);
    Int_t ind=(lad-1)*fDetTypeRec->GetITSgeom()->GetNdetectors(lay)+(det-1);
    Int_t lyr=(lay-1);
    // get rec points
    for (Int_t i=0; i<nclus; i++){
        l[0] = kconv*xcenter[i];
        l[1] = kconv*GetSeg()->Dy()/2.;
        l[2] = kconv*zcenter[i];

        xg = l[0]; 
        zg = l[2]; 

        Double_t sigma2x = (kconv*errxcenter[i]) * (kconv*errxcenter[i]);
        Double_t sigma2z = (kconv*errzcenter[i]) * (kconv*errzcenter[i]);
        AliITSRecPoint rnew(fDetTypeRec->GetITSgeom());
        rnew.SetXZ(fModule,xg,zg);
        rnew.SetQ(1.);
        rnew.SetdEdX(0.);
        rnew.SetSigmaDetLocX2(sigma2x);
        rnew.SetSigmaZ2(sigma2z);
        rnew.SetLabel(tr1clus[i],0);
        rnew.SetLabel(tr2clus[i],1);
        rnew.SetLabel(tr3clus[i],2);
	rnew.SetDetectorIndex(ind);
	rnew.SetLayer(lyr);
	fDetTypeRec->AddRecPoint(rnew); 
    } // end for i
}
