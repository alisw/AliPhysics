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

// **************************************************************************
//  * The package was revised and changed by Boris Batiounia in the time     *
//  * period of March - June 2001                                            *
// **************************************************************************/
//
#include <Riostream.h>
#include <TArrayI.h>
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSdigitSSD.h"
#include "AliITSRawClusterSSD.h"
#include "AliITSRecPoint.h"
#include "AliITSMapA1.h"
#include "AliITSClusterFinderSSD.h"
#include "AliITSclusterSSD.h"
#include "AliITSpackageSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSgeom.h"

const Bool_t AliITSClusterFinderSSD::fgkSIDEP=kTRUE;
const Bool_t AliITSClusterFinderSSD::fgkSIDEN=kFALSE;

ClassImp(AliITSClusterFinderSSD)

//____________________________________________________________________
//
//  Constructor
//______________________________________________________________________
AliITSClusterFinderSSD::AliITSClusterFinderSSD():
AliITSClusterFinder(),
fITS(0),
fClusterP(0),
fNClusterP(0),
fClusterN(0),
fNClusterN(0),
fPackages(0),
fNPackages(0),
fDigitsIndexP(0),
fNDigitsP(0),
fDigitsIndexN(0),
fNDigitsN(0),
fPitch(0.0),
fTanP(0.0),
fTanN(0.0),
fPNsignalRatio(0.0),
fSFF(0),
fSFB(0){
    //Default constructor
}
//______________________________________________________________________
AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSsegmentation *seg,
                                               TClonesArray *digits):
AliITSClusterFinder(seg,0),
fITS(0),
fClusterP(0),
fNClusterP(0),
fClusterN(0),
fNClusterN(0),
fPackages(0),
fNPackages(0),
fDigitsIndexP(0),
fNDigitsP(0),
fDigitsIndexN(0),
fNDigitsN(0),
fPitch(0.0),
fTanP(0.0),
fTanN(0.0),
fPNsignalRatio(0.0),
fSFF(0),
fSFB(0){
    //Standard constructor

    SetDigits(digits);
    SetMap(new AliITSMapA1(GetSeg(),Digits()));
    fITS          = (AliITS*)gAlice->GetModule("ITS");
    fClusterP     = new TClonesArray ("AliITSclusterSSD",200);    
    fNClusterP    = 0;
    fClusterN     = new TClonesArray ("AliITSclusterSSD",200);   
    fNClusterN    = 0;
    fPackages     = new TClonesArray ("AliITSpackageSSD",200);    //packages  
    fNPackages    = 0;
    fDigitsIndexP = new TArrayI(300);
    fNDigitsP     = 0;
    fDigitsIndexN = new TArrayI(300);
    fNDigitsN     = 0;
    fPitch        = GetSeg()->Dpx(0);
    fPNsignalRatio= 7./8.;    // warning: hard-wired number
}
//______________________________________________________________________}
AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSsegmentation *seg,
                                               AliITSresponse *res):
AliITSClusterFinder(seg,res),
fITS(0),
fClusterP(0),
fNClusterP(0),
fClusterN(0),
fNClusterN(0),
fPackages(0),
fNPackages(0),
fDigitsIndexP(0),
fNDigitsP(0),
fDigitsIndexN(0),
fNDigitsN(0),
fPitch(0.0),
fTanP(0.0),
fTanN(0.0),
fPNsignalRatio(0.0),
fSFF(0),
fSFB(0){
    //Standard constructor

    fITS          = (AliITS*)gAlice->GetModule("ITS");
    fClusterP     = new TClonesArray ("AliITSclusterSSD",200);    
    fNClusterP    = 0;
    fClusterN     = new TClonesArray ("AliITSclusterSSD",200);   
    fNClusterN    = 0;
    fPackages     = new TClonesArray ("AliITSpackageSSD",200);    //packages  
    fNPackages    = 0;
    fDigitsIndexP = new TArrayI(300);
    fNDigitsP     = 0;
    fDigitsIndexN = new TArrayI(300);
    fNDigitsN     = 0;
    fPitch        = GetSeg()->Dpx(0);
    fPNsignalRatio= 7./8.;    // warning: hard-wired number
}
//______________________________________________________________________
AliITSClusterFinderSSD::~AliITSClusterFinderSSD(){
    // Default destructor

    fITS = 0;
    delete fClusterP;
    delete fClusterN;        
    delete fPackages;        
    delete fDigitsIndexP;        
    delete fDigitsIndexN; 
}
//______________________________________________________________________
void AliITSClusterFinderSSD::InitReconstruction(){
    // initialization of the cluster finder

    register Int_t i; //iterator

    for (i=0;i<fNClusterP;i++) fClusterP->RemoveAt(i);
    fNClusterP  =0;
    for (i=0;i<fNClusterN;i++) fClusterN->RemoveAt(i);
    fNClusterN=0;
    for (i=0;i<fNPackages;i++) fPackages->RemoveAt(i);
    fNPackages = 0;
    fNDigitsP  = 0;
    fNDigitsN  = 0;
    Float_t stereoP,stereoN;
    GetSeg()->Angles(stereoP,stereoN);
    CalcStepFactor(stereoP,stereoN);
    if(GetDebug(1)) cout<<"fSFF = "<<fSFF<<"  fSFB = "<<fSFB<<"\n";
}
//______________________________________________________________________
void AliITSClusterFinderSSD::FindRawClusters(Int_t module){
    // This function findes out all clusters belonging to one module
    // 1. Zeroes all space after previous module reconstruction
    // 2. Finds all neighbouring digits, create clusters
    // 3. If necesery, resolves for each group of neighbouring digits 
    //    how many clusters creates it.
    // 4. Colculate the x and z coordinate  
    Int_t lay, lad, detect;
    AliITSgeom *geom = fITS->GetITSgeom();

    SetModule(module);
    geom->GetModuleId(GetModule(),lay, lad, detect);
    if ( lay == 6 ) ((AliITSsegmentationSSD*)GetSeg())->SetLayer(6);
    if ( lay == 5 ) ((AliITSsegmentationSSD*)GetSeg())->SetLayer(5);

    InitReconstruction();  //ad. 1
    Map()->FillMap();
    FillDigitsIndex();
    SortDigits();
    FindNeighbouringDigits(); //ad. 2
    //SeparateOverlappedClusters();  //ad. 3
    ClustersToPackages();  //ad. 4
    Map()->ClearMap();
}
//______________________________________________________________________
void AliITSClusterFinderSSD::FindNeighbouringDigits(){
    //If there are any digits on this side, create 1st Cluster,
    // add to it this digit, and increment number of clusters
    register Int_t i;

    if ((fNDigitsP>0 )  && (fNDigitsN > 0 )) {
	Int_t currentstripNo;
	Int_t *dbuffer = new Int_t [300];   //buffer for strip numbers
	Int_t dnumber;    //curent number of digits in buffer
	TArrayI      &lDigitsIndexP = *fDigitsIndexP;
	TArrayI      &lDigitsIndexN = *fDigitsIndexN;
	TObjArray    &lDigits       = *(Digits());
	TClonesArray &lClusterP     = *fClusterP;
	TClonesArray &lClusterN     = *fClusterN;
	//process P side 
	dnumber = 1;
	dbuffer[0]=lDigitsIndexP[0];
	//If next digit is a neigh. of previous, adds to last clust. this digit
	for (i=1; i<fNDigitsP; i++) {
	    //reads new digit
	    currentstripNo = ((AliITSdigitSSD*)lDigits[lDigitsIndexP[i]])->
		                                            GetStripNumber(); 
	    //check if it is a neighbour of a previous one
	    if((((AliITSdigitSSD*)lDigits[lDigitsIndexP[i-1]])->
                                                            GetStripNumber()) 
	       ==  (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexP[i];
	    else{
		//create a new one side cluster
		new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,
							      Digits(),
							      fgkSIDEP); 
		dbuffer[0]=lDigitsIndexP[i];
		dnumber = 1;
	    } // end if else
	} // end loop over fNDigitsP
	new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,
						      Digits(),fgkSIDEP);
	//process N side 
	//for comments, see above
	dnumber = 1;
	dbuffer[0]=lDigitsIndexN[0];
	//If next digit is a neigh. of previous, adds to last clust. this digit
	for (i=1; i<fNDigitsN; i++) { 
	    currentstripNo = ((AliITSdigitSSD*)(lDigits[lDigitsIndexN[i]]))->
                                                            GetStripNumber();
	    if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexN[i-1]])->
                                                            GetStripNumber()) 
		 == (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexN[i];
	    else {
		new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,
                                                        Digits(),
                                                        fgkSIDEN);
		dbuffer[0]=lDigitsIndexN[i];
		dnumber = 1;
	    } // end if else
	} // end loop over fNDigitsN
	new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,
                                                   Digits(),fgkSIDEN);
	delete [] dbuffer;

    } // end condition on  NDigits 

    if (GetDebug(1)) cout<<"\n Found clusters: fNClusterP = "<<fNClusterP
		   <<"  fNClusterN ="<<fNClusterN<<"\n";
}
//______________________________________________________________________
void AliITSClusterFinderSSD::SeparateOverlappedClusters(){
    // overlapped clusters separation
    register Int_t i; //iterator
    Double_t  factor=0.75;            // How many percent must be lower signal 
                                     // on the middle one digit
                                     // from its neighbours
    Int_t    signal0;              //signal on the strip before the current one
    Int_t    signal1;              //signal on the current one signal
    Int_t    signal2;              //signal on the strip after the current one
    TArrayI *splitlist;              //  List of splits
    Int_t    numerofsplits=0;        // number of splits
    Int_t    initPsize = fNClusterP; //initial size of the arrays 
    Int_t    initNsize = fNClusterN; //we have to keep it because it will grow 
                                     // in this function and it doasn't make 
                                     // sense to pass through it again
    splitlist = new TArrayI(300);

    for (i=0;i<initPsize;i++){
	if (( ((AliITSclusterSSD*)(*fClusterP)[i])->
	      GetNumOfDigits())==1) continue;
	if (( ((AliITSclusterSSD*)(*fClusterP)[i])->
	      GetNumOfDigits())==2) continue;
        Int_t nj=(((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits()-1);
        for (Int_t j=1; j<nj; j++){
            signal1=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j);
            signal0=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j-1);
            signal2=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j+1);
            //if signal is less then factor*signal of its neighbours
            if (  (signal1<(factor*signal0)) && (signal1<(factor*signal2)) ){
		(*splitlist)[numerofsplits++]=j;
	    } // end if
	} // end loop over number of digits
	//split this cluster if necessary
	if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,fgkSIDEP);
	numerofsplits=0;
	//in signed places (splitlist)
    } // end loop over clusters on Pside

    for (i=0;i<initNsize;i++) {
	if (( ((AliITSclusterSSD*)(*fClusterN)[i])->
	      GetNumOfDigits())==1) continue;
	if (( ((AliITSclusterSSD*)(*fClusterN)[i])->
	      GetNumOfDigits())==2) continue;
	Int_t nj=(((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits()-1);
	for (Int_t j=1; j<nj; j++){
            signal1=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j);
            signal0=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j-1);
            signal2=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j+1);
            //if signal is less then factor*signal of its neighbours
            if (  (signal1<(factor*signal0)) && (signal1<(factor*signal2)) ) 
		(*splitlist)[numerofsplits++]=j;  
	} // end loop over number of digits 
	//split this cluster into more clusters
	if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,fgkSIDEN);
	numerofsplits=0;
	//in signed places (splitlist)
    } // end loop over clusters on Nside

    delete splitlist;
}
//______________________________________________________________________
void AliITSClusterFinderSSD::SplitCluster(TArrayI *list, Int_t nsplits,
					  Int_t index, Bool_t side){
    //This function splits one side cluster into more clusters
    //number of splits is defined by "nsplits"
    //Place of splits are defined in the TArray "list"
    // For further optimisation: Replace this function by two 
    // specialised ones (each for one side)
    // save one "if"
    //For comlete comments see AliITSclusterSSD::SplitCluster
    register Int_t i; //iterator
    AliITSclusterSSD* curentcluster;
    Int_t   *tmpdigits = new Int_t[100];
    Int_t    nn;

    // side true means P side
    if (side) {
	curentcluster =((AliITSclusterSSD*)((*fClusterP)[index])) ;
	for (i = nsplits; i>0 ;i--) {  
	    nn=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
	    new ((*fClusterP)[fNClusterP]) AliITSclusterSSD(nn,tmpdigits,
							    Digits(),side);
	    ( (AliITSclusterSSD*)((*fClusterP)[fNClusterP]) )->
                                                      SetLeftNeighbour(kTRUE);
	    //if left cluster had neighbour on the right before split 
	    //new should have it too
	    if ( curentcluster->GetRightNeighbour() ) 
		( (AliITSclusterSSD*)((*fClusterP)[fNClusterP]) )->
                                                     SetRightNeighbour(kTRUE);
	    else curentcluster->SetRightNeighbour(kTRUE); 
	    fNClusterP++;
	} // end loop over nplits
    } else {
	curentcluster =((AliITSclusterSSD*)((*fClusterN)[index]));
	for (i = nsplits; i>0 ;i--) {  
	    nn=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
	    new ((*fClusterN)[fNClusterN]) AliITSclusterSSD(nn,tmpdigits,
							    Digits(),side);
	    ((AliITSclusterSSD*)((*fClusterN)[fNClusterN]))->
                                                    SetRightNeighbour(kTRUE);
	    if (curentcluster->GetRightNeighbour())
		( (AliITSclusterSSD*)( (*fClusterN)[fNClusterN]) )->
                                                     SetRightNeighbour(kTRUE);
	    else curentcluster->SetRightNeighbour(kTRUE);      
	    fNClusterN++;
	} // end loop over nplits
    } // end if side
    delete []tmpdigits;
}
//______________________________________________________________________
Int_t AliITSClusterFinderSSD::SortDigitsP(Int_t start, Int_t end){
    // sort digits on the P side
    Int_t right;
    Int_t left;

    if (start != (end - 1) ){
	left=this->SortDigitsP(start,(start+end)/2);
	right=this->SortDigitsP((start+end)/2,end);  
	return (left || right);
    }else{ 
	left =  ((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexP)[start]]))->
                                                              GetStripNumber();
	right= ((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexP)[end]]))->
                                                              GetStripNumber();
	if( left > right ){
	    Int_t tmp = (*fDigitsIndexP)[start];
	    (*fDigitsIndexP)[start]=(*fDigitsIndexP)[end];
	    (*fDigitsIndexP)[end]=tmp;
	    return 1;
	}else return 0;
    } // end if
}
//______________________________________________________________________
Int_t AliITSClusterFinderSSD::SortDigitsN(Int_t start, Int_t end){
    // sort digits on the N side
    Int_t right;
    Int_t left;

    if (start != (end - 1)){
        left=this->SortDigitsN(start,(start+end)/2);
        right=this->SortDigitsN((start+end)/2,end);  
        return (left || right);
    }else{
        left =((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexN)[start]]))->
            GetStripNumber();
        right=((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexN)[end]]))->
            GetStripNumber();
        if ( left > right ){
            Int_t tmp = (*fDigitsIndexN)[start];
            (*fDigitsIndexN)[start]=(*fDigitsIndexN)[end];
            (*fDigitsIndexN)[end]=tmp;
            return 1;
        }else return 0;
    } // end if
}
//______________________________________________________________________
void AliITSClusterFinderSSD::FillDigitsIndex(){
    //Fill the indexes of the clusters belonging to a given ITS module
    Int_t pns=0, nns=0;
    Int_t tmp,bit,k;
    Int_t noentries;
    Int_t i;

    noentries = NDigits();

    Int_t* psidx = new Int_t [noentries*sizeof(Int_t)];
    Int_t* nsidx = new Int_t [noentries*sizeof(Int_t)]; 
    if (fDigitsIndexP==NULL) fDigitsIndexP = new TArrayI(noentries);
    if (fDigitsIndexN==NULL) fDigitsIndexN = new TArrayI(noentries);

    AliITSdigitSSD *dig;

    for ( i = 0 ; i< noentries; i++ ) {
        dig = (AliITSdigitSSD*)GetDigit(i);
        if(dig->IsSideP()) { 
            bit=1;
            tmp=dig->GetStripNumber();
            // I find this totally unnecessary - it's just a 
            // CPU consuming double check
            for( k=0;k<pns;k++){
                if (tmp==psidx[k]){
                    if (GetDebug(1)) cout<<"Such a digit exists \n";
                    bit=0;
                } // end if
            } // end for k
            // end comment 
            if(bit) {
                fDigitsIndexP->AddAt(i,fNDigitsP++);
                psidx[pns++]=tmp;
            } // end if bit
        } else {
            bit=1;
            tmp=dig->GetStripNumber();
            // same as above
            for( k=0;k<nns;k++){
                if (tmp==nsidx[k]){
                    if (GetDebug(1)) cout<<"Such a digit exists \n";
                    bit=0;
                } // end if
            } // for k
            // end comment
            if (bit) {
                fDigitsIndexN->AddAt(i,fNDigitsN++);
                nsidx[nns++] =tmp;
            } // end if bit
        } // end if
    } // end for i

    delete [] psidx;
    delete [] nsidx;

    if(GetDebug(1)) cout<<"Digits: P = "<<fNDigitsP<<" N = "<<fNDigitsN<<endl;
}
//______________________________________________________________________
void AliITSClusterFinderSSD::SortDigits(){
    // sort digits
    Int_t i;

    if(fNDigitsP>1) for (i=0;i<fNDigitsP-1;i++)
        if (SortDigitsP(0,(fNDigitsP-1-i))==0) break;
    if(fNDigitsN>1) for (i=0;i<fNDigitsN-1;i++)
        if(SortDigitsN(0,(fNDigitsN-1-i))==0) break;
}
//______________________________________________________________________
void AliITSClusterFinderSSD::FillClIndexArrays(Int_t* arrayP,Int_t *arrayN)
    const{
    // fill cluster index array
    register Int_t i;

    for (i=0; i<fNClusterP;i++) arrayP[i]=i;
    for (i=0; i<fNClusterN;i++) arrayN[i]=i;
}
//______________________________________________________________________
void AliITSClusterFinderSSD::SortClusters(Int_t* arrayP, Int_t *arrayN){
    // sort clusters
    Int_t i;

    if(fNClusterP>1) for (i=0;i<fNClusterP-1;i++)
        if (SortClustersP(0,(fNClusterP-1),arrayP)==0)  break;
    if(fNClusterN>1) for (i=0;i<fNClusterN-1;i++)
        if (SortClustersN(0,(fNClusterN-1),arrayN)==0)  break;
}
//______________________________________________________________________
Int_t AliITSClusterFinderSSD::SortClustersP(Int_t start, Int_t end,
					    Int_t *array){
    //Sort P side clusters
    Int_t right;
    Int_t left;

    if (start != (end - 1) ) {
        left=this->SortClustersP(start,(start+end)/2,array);
        right=this->SortClustersP((start+end)/2,end,array);  
        return (left || right);
    } else {
        left =((AliITSclusterSSD*)((*fClusterP)[array[start]]))->
            GetDigitStripNo(0);
        right=((AliITSclusterSSD*)((*fClusterP)[array[ end ]]))->
            GetDigitStripNo(0);
        if(left>right) {
            Int_t tmp = array[start];
            array[start]=array[end];
            array[end]=tmp;
            return 1;
        } else return 0;
    } // end if
}
//______________________________________________________________________
Int_t AliITSClusterFinderSSD::SortClustersN(Int_t start, Int_t end, 
                                            Int_t *array){
    //Sort N side clusters
    Int_t right;
    Int_t left;

    if (start != (end - 1) ) {
        left=this->SortClustersN(start,(start+end)/2,array);
        right=this->SortClustersN((start+end)/2,end,array);  
        return (left || right);
    } else {
        left =((AliITSclusterSSD*)((*fClusterN)[array[start]]))->
            GetDigitStripNo(0);
        right=((AliITSclusterSSD*)((*fClusterN)[array[ end ]]))->
            GetDigitStripNo(0);
        if( left > right) {
            Int_t tmp = array[start];
            array[start]=array[end];
            array[end]=tmp;
            return 1;
        } else return 0;
    } // end if
}
//______________________________________________________________________
void AliITSClusterFinderSSD::ClustersToPackages(){
    // fill packages   
    
    Int_t *oneSclP = new Int_t[fNClusterP];//I want to have sorted 1 S clusters
    Int_t *oneSclN = new Int_t[fNClusterN];//I can not sort it in TClonesArray
                                           //so, I create table of indexes and 
                                           //sort it
                                           //I do not use TArrayI on purpose
                                           //MB: well, that's not true that one
                                           //cannot sort objs in TClonesArray
    AliITSclusterSSD *currentP;
    AliITSclusterSSD *currentN;
    Int_t j1, j2;    

    //Fills in One Side Clusters Index Array
    FillClIndexArrays(oneSclP,oneSclN); 
    //Sorts filled Arrays
    //SortClusters(oneSclP,oneSclN);                   

    fNPackages=1;      
    new ((*fPackages)[0]) AliITSpackageSSD(fClusterP,fClusterN);

    //This part was includede by Boris Batiounia in March 2001.
    // Take all recpoint pairs (x coordinates) in both P and N sides  
    // to calculate z coordinates of the recpoints

    for (j1=0;j1<fNClusterP;j1++) {  
	currentP = GetPSideCluster(oneSclP[j1]);
	Double_t xP = currentP->GetPosition();
	Double_t signalP = currentP->GetTotalSignal();
	for (j2=0;j2<fNClusterN;j2++) {  
	    currentN = GetNSideCluster(oneSclN[j2]);
	    Double_t xN = currentN->GetPosition();
	    Double_t signalN = currentN->GetTotalSignal();
	    CreateNewRecPoint(xP,1,xN,1,signalP,signalN,currentP,currentN,
                           0.75);
	} // end for j2
    } // end for j1

    delete [] oneSclP;
    delete [] oneSclN;
}
//______________________________________________________________________
Bool_t AliITSClusterFinderSSD::CreateNewRecPoint(Double_t P,Double_t dP,
						 Double_t N, Double_t dN,
						 Double_t SigP,Double_t SigN, 
						 AliITSclusterSSD *clusterP,
						 AliITSclusterSSD *clusterN,
						 Stat_t prob){
    // create the recpoints
    const Double_t kADCtoKeV = 2.16; 
    // 50 ADC units -> 30000 e-h pairs; 1e-h pair -> 3.6e-3 KeV;
    // 1 ADC unit -> (30000/50)*3.6e-3 = 2.16 KeV 
    const Double_t kconv = 1.0e-4;
    const Double_t kRMSx = 20.0*kconv; 
    const Double_t kRMSz = 800.0*kconv;
    Int_t n=0;
    Int_t *tr;
    Int_t ntracks;

    if (GetCrossing(P,N)) {
        //GetCrossingError(dP,dN);
        dP = dN = prob = 0.0; // to remove unused variable warning.
        AliITSRawClusterSSD cnew;
        Int_t nstripsP=clusterP->GetNumOfDigits();
        Int_t nstripsN=clusterN->GetNumOfDigits();
        Double_t signal = 0;
        Double_t dedx = 0;
        if(SigP>SigN) {
            signal = SigP;
            dedx = SigP*kADCtoKeV;
        }else{
            signal = SigN;
            dedx = SigN*kADCtoKeV;
        } // end if SigP>SigN
        tr = (Int_t*) clusterP->GetTracks(n);
        ntracks = clusterP->GetNTracks();
        cnew.SetSignalP(SigP);
        cnew.SetSignalN(SigN);
        cnew.SetMultiplicity(nstripsP);
        cnew.SetMultN(nstripsN);
        cnew.SetQErr(TMath::Abs(SigP-SigN));
        cnew.SetNTrack(ntracks);
        fITS->AddCluster(2,&cnew);
        AliITSRecPoint rnew;
        rnew.SetX(P*kconv);
        rnew.SetZ(N*kconv);
        rnew.SetQ(signal);
        rnew.SetdEdX(dedx);
        rnew.SetSigmaX2( kRMSx* kRMSx); 
        rnew.SetSigmaZ2( kRMSz* kRMSz);
        rnew.fTracks[0]=tr[0];
        rnew.fTracks[1]=tr[1];
        rnew.fTracks[2]=tr[2];
        fITS->AddRecPoint(rnew);
        return kTRUE;
    } // end if
    return kFALSE;  
}
//______________________________________________________________________
void  AliITSClusterFinderSSD::CalcStepFactor(Double_t Psteo, Double_t Nsteo){
    // calculate the step factor for matching clusters
    // 95 is the pitch, 4000 - dimension along z ?
    Double_t dz=GetSeg()->Dz();

    fSFF = ( (Int_t)  (Psteo*dz/fPitch ) );// +1;
    fSFB = ( (Int_t)  (Nsteo*dz/fPitch ) );// +1;
}
//______________________________________________________________________
AliITSclusterSSD* AliITSClusterFinderSSD::GetPSideCluster(Int_t idx){
    // get P side clusters

    if((idx<0)||(idx>=fNClusterP)){
        Info("GetPSideCluster","0<index=%d<=%d out of range",idx,fNClusterP);
        return 0;
    }else{
        return (AliITSclusterSSD*)((*fClusterP)[idx]);
    } // end if
}
//______________________________________________________________________
AliITSclusterSSD* AliITSClusterFinderSSD::GetNSideCluster(Int_t idx){
    // get N side clusters

    if((idx<0)||(idx>=fNClusterN)){
        Info("GetNSideCluster","0<index=%d >= %d out of range",idx,fNClusterN);
        return 0;
    }else{
        return (AliITSclusterSSD*)((*fClusterN)[idx]);
    } // end if
}
//______________________________________________________________________
Bool_t AliITSClusterFinderSSD::GetCrossing (Double_t &P, Double_t &N){ 
    // get crossing
    // This function was rivised and changed by Boris Batiounia in March 2001
    Double_t dx = GetSeg()->Dx(); // detector size in x direction, microns
    Double_t dz = GetSeg()->Dz(); // detector size in z direction, microns
    Double_t xL; // x local coordinate
    Double_t zL; // z local coordinate
    Double_t x;  // x = xL + dx/2
    Double_t z;  // z = zL + dz/2
    Double_t xP; // x coordinate in the P side from the first P strip
    Double_t xN; // x coordinate in the N side from the first N strip
    Float_t stereoP,stereoN;

    GetSeg()->Angles(stereoP,stereoN);
    fTanP=TMath::Tan(stereoP);
    fTanN=TMath::Tan(stereoN);
    Double_t kP = fTanP; // Tangent of 0.0075 mrad
    Double_t kN = fTanN; // Tangent of 0.0275 mrad
    P *= fPitch;
    N *= fPitch; 

    xP = N;      // change the mistake for the P/N
    xN = P;      // coordinates correspondence in this function

    x = xP + kP*(dz*kN-xP+xN)/(kP+kN);
    z = (dz*kN-xP+xN)/(kP+kN); 
    xL = x - dx/2;
    zL = z - dz/2;
    P = xL;
    N = zL;  

    if(TMath::Abs(xL) > dx/2 || TMath::Abs(zL) > dz/2) return kFALSE;
    
    // Check that xL and zL are inside the detector for the 
    // correspondent xP and xN coordinates

    return kTRUE;   
}
//______________________________________________________________________
void AliITSClusterFinderSSD::GetCrossingError(Double_t& dP, Double_t& dN){
    // get crossing error
    Double_t dz, dx;

    dz = TMath::Abs(( dP + dN )*fPitch/(fTanP + fTanN) );
    dx = fPitch*(TMath::Abs(dP*(1 - fTanP/(fTanP + fTanN))) +
                 TMath::Abs(dN *fTanP/(fTanP + fTanN) ));
    dN = dz;
    dP = dx;
}
