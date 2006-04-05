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
//#include <Riostream.h>
#include <TArrayI.h>

//#include "AliRun.h"
//#include "AliITS.h"
#include "AliITSClusterFinderSSD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSMapA1.h"
#include "AliITSRawClusterSSD.h"
#include "AliITSRecPoint.h"
#include "AliITSdigitSSD.h"
#include "AliITSclusterSSD.h"
#include "AliITSpackageSSD.h"
#include "AliITSsegmentationSSD.h"
//#include "AliITSgeom.h"
#include "AliITSCalibrationSSD.h"
#include "AliLog.h"

const Bool_t AliITSClusterFinderSSD::fgkSIDEP=kTRUE;
const Bool_t AliITSClusterFinderSSD::fgkSIDEN=kFALSE;
const Int_t AliITSClusterFinderSSD::fgkNoiseThreshold=5;
const Int_t debug=0;

ClassImp(AliITSClusterFinderSSD)
  
  //____________________________________________________________________
  //
  //  Constructor
  //______________________________________________________________________
  AliITSClusterFinderSSD::AliITSClusterFinderSSD():
AliITSClusterFinder(),
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
AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSDetTypeRec* dettyp,
                                               TClonesArray *digits):
AliITSClusterFinder(dettyp,digits),
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
    fClusterP     = new TClonesArray ("AliITSclusterSSD",200);
    fNClusterP    = 0;
    fClusterN     = new TClonesArray ("AliITSclusterSSD",200);
    fNClusterN    = 0;
    fPackages     = new TClonesArray ("AliITSpackageSSD",200);    //packages
    fNPackages    = 0;
    fDigitsIndexP = new TArrayI(800);
    fNDigitsP     = 0;
    fDigitsIndexN = new TArrayI(800);
    fNDigitsN     = 0;
    fPitch        = GetSeg()->Dpx(0);
    fPNsignalRatio= 7./8.;    // warning: hard-wired number
}

//______________________________________________________________________
AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSDetTypeRec* dettyp):
AliITSClusterFinder(dettyp),
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

    fClusterP     = new TClonesArray ("AliITSclusterSSD",200);
    fNClusterP    = 0;
    fClusterN     = new TClonesArray ("AliITSclusterSSD",200);
    fNClusterN    = 0;
    fPackages     = new TClonesArray ("AliITSpackageSSD",200);    //packages
    fNPackages    = 0;
    fDigitsIndexP = new TArrayI(800);
    fNDigitsP     = 0;
    fDigitsIndexN = new TArrayI(800);
    fNDigitsN     = 0;
    fPitch        = GetSeg()->Dpx(0);
    fPNsignalRatio= 7./8.;    // warning: hard-wired number
}

//______________________________________________________________________
AliITSClusterFinderSSD::~AliITSClusterFinderSSD(){
  // Default destructor
  
  delete fClusterP;
  delete fClusterN;        
  delete fPackages;        
  delete fDigitsIndexP;        
  delete fDigitsIndexN; 
  delete fMap;
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
    //fSegmentation->Angles(stereoP,stereoN);
     GetSeg()->Angles(stereoP,stereoN);
   CalcStepFactor(stereoP,stereoN);
    if (debug) cout<<"fSFF = "<<fSFF<<"  fSFB = "<<fSFB<<"\n";
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
  
  //cout<<"clusterfinder: this is module "<<module<<endl;
  
  //    AliITSgeom *geom = fITS->GetITSgeom();
  
  if(!fDetTypeRec->GetITSgeom()) {
    Error("FindRawClusters","ITS geom is null!");
    return;
  }
  
  SetModule(module);
  fDetTypeRec->GetITSgeom()->GetModuleId(GetModule(),lay, lad, detect);
  //   geom->GetModuleId(module,lay, lad, detect);
  
  //cout<<"layer = "<<lay<<endl;
  
  if ( lay == 6 ) ((AliITSsegmentationSSD*)GetSeg())->SetLayer(6);
  if ( lay == 5 ) ((AliITSsegmentationSSD*)GetSeg())->SetLayer(5);
  
  
  InitReconstruction();  //ad. 1
  fMap->FillMap();
  
  FillDigitsIndex();
  if ( (fNDigitsP==0 )  || (fNDigitsN == 0 ))  return;
  
  SortDigits();
  FindNeighbouringDigits(module); //ad. 2

  //SeparateOverlappedClusters();  //ad. 3
  ClustersToPackages();  //ad. 4

  fMap->ClearMap();

}
//______________________________________________________________________
void AliITSClusterFinderSSD::FindNeighbouringDigits(Int_t module){
    //If there are any digits on this side, create 1st Cluster,
    // add to it this digit, and increment number of clusters

    register Int_t i,j;
    Int_t flag=0;

    //if ( (fNDigitsP==0 )  || (fNDigitsN == 0 ))  return;
    
    Int_t currentstripNo;
    Int_t *dbuffer = new Int_t [800];   //buffer for strip numbers
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
      currentstripNo = ((AliITSdigitSSD*)lDigits[lDigitsIndexP[i]])->GetStripNumber(); 
      
      if((((AliITSdigitSSD*)lDigits[lDigitsIndexP[i-1]])->GetStripNumber()) ==  (currentstripNo-1) ) 
	dbuffer[dnumber++]=lDigitsIndexP[i];
      
      else{
	  
	// check signal
	for(j=0;j<dnumber;j++) {
	  
	  if( (((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetSignal())
	      > fgkNoiseThreshold* ((AliITSCalibrationSSD*)GetResp(module))->
	      GetNoiseP( ((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetStripNumber()) ) 
	    flag+=1; 
	  
	}
	
	//if(flag==dnumber) {
	if(flag>0) {
	  //create a new one side cluster
	  new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,
							Digits(),
							fgkSIDEP); 
	  
	  flag=0;
	}
	
	flag=0;
	dbuffer[0]=lDigitsIndexP[i];
	dnumber = 1;
	
      } // end if else
      
    } // end loop over fNDigitsP
    
    // check signal
    for(j=0;j<dnumber;j++) {
      
      if( (((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetSignal())
	  > fgkNoiseThreshold*((AliITSCalibrationSSD*)GetResp(module))->
	  GetNoiseP( ((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetStripNumber()) ) 
	flag+=1; 
      
    }
    
    //if(flag==dnumber) {
    if(flag>0) {
      //create a new one side cluster
      new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,
						    Digits(),
						    fgkSIDEP); 
      
      flag=0;
    }
    
    flag=0;
    
    //process N side 
    //for comments, see above
    dnumber = 1;
    dbuffer[0]=lDigitsIndexN[0];
    //If next digit is a neigh. of previous, adds to last clust. this digit
    
    for (i=1; i<fNDigitsN; i++) { 
      currentstripNo = ((AliITSdigitSSD*)(lDigits[lDigitsIndexN[i]]))->GetStripNumber();
      
      if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexN[i-1]])->GetStripNumber()) == (currentstripNo-1) ) 
	dbuffer[dnumber++]=lDigitsIndexN[i];
      
      else {
	
	// check signal
	for(j=0;j<dnumber;j++) {
	  
	  if( (((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetSignal())
	      > fgkNoiseThreshold*((AliITSCalibrationSSD*)GetResp(module))->
	      GetNoiseN( ((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetStripNumber()) ) 
	    flag+=1; 
	  
	}
	
	//if(flag==dnumber) {
	if(flag>0) {
	  //create a new one side cluster
	  new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,
							Digits(),
							fgkSIDEN);
	  
	  flag=0;
	}
	
	flag=0;
	dbuffer[0]=lDigitsIndexN[i];
	dnumber = 1;
      } // end if else
    } // end loop over fNDigitsN
    
    // check signal
    for(j=0;j<dnumber;j++) {
      
      if( (((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetSignal())
	  > fgkNoiseThreshold*((AliITSCalibrationSSD*)GetResp(module))->
	  GetNoiseN( ((AliITSdigitSSD*)lDigits[dbuffer[j]])->GetStripNumber()) ) 
	flag+=1; 
      
    }
    
    //if(flag==dnumber) {
    if(flag>0) {
      //create a new one side cluster
      new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,
						    Digits(),
						    fgkSIDEN);
      
      flag=0;
    }
    
    flag=0;
    delete [] dbuffer;
    
    if (debug) cout<<"\n Found clusters: fNClusterP = "<<fNClusterP
		   <<"  fNClusterN ="<<fNClusterN<<"\n";
}
//______________________________________________________________________
void AliITSClusterFinderSSD::SeparateOverlappedClusters(){
    // overlapped clusters separation
    register Int_t i; //iterator
    Float_t  factor=0.75;            // How many percent must be lower signal 
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
    splitlist = new TArrayI(800);

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

    Int_t noentries;
    Int_t i;
    Int_t gain=0;
    Int_t signal=0;

    noentries = fDigits->GetEntriesFast();

    if (fDigitsIndexP==NULL) fDigitsIndexP = new TArrayI(noentries);
    if (fDigitsIndexN==NULL) fDigitsIndexN = new TArrayI(noentries);

    AliITSdigitSSD *dig;

    for ( i = 0 ; i< noentries; i++ ) {
      
      dig = (AliITSdigitSSD*)GetDigit(i);
      
      gain=(Int_t) ((AliITSCalibrationSSD*)GetResp(fModule))->GetGainP(dig->GetStripNumber());  
      signal=gain*dig->GetSignal();
      dig->SetSignal(signal);
      
      if(dig->IsSideP()) fDigitsIndexP->AddAt(i,fNDigitsP++);
      else fDigitsIndexN->AddAt(i,fNDigitsN++);
      
    } // end for i
    
    //    delete [] psidx;
    //delete [] nsidx;

    if (debug) cout<<"Digits :  P = "<<fNDigitsP<<"   N = "<<fNDigitsN<<endl;
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
void AliITSClusterFinderSSD::FillClIndexArrays(Int_t* arrayP, Int_t *arrayN) const{
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
	Float_t signalP = currentP->GetTotalSignal();
	for (j2=0;j2<fNClusterN;j2++) {  
	    currentN = GetNSideCluster(oneSclN[j2]);
	    Double_t xN = currentN->GetPosition();
	    Float_t signalN = currentN->GetTotalSignal();
	    CreateNewRecPoint(xP,1,xN,1,signalP,signalN,currentP,currentN,
			      0.75);
	} // end for j2
    } // end for j1

    delete [] oneSclP;
    delete [] oneSclN;
}
//______________________________________________________________________
Bool_t AliITSClusterFinderSSD::CreateNewRecPoint(Float_t P,Float_t dP,
						 Float_t N, Float_t dN,
						 Float_t SigP,Float_t SigN, 
						 AliITSclusterSSD *clusterP,
						 AliITSclusterSSD *clusterN,
						 Stat_t prob){
    // create the recpoints
    const Float_t kADCtoKeV = 2.16; 
    // 50 ADC units -> 30000 e-h pairs; 1e-h pair -> 3.6e-3 KeV;
    // 1 ADC unit -> (30000/50)*3.6e-3 = 2.16 KeV 
    const Float_t kconv = 1.0e-4;
    const Float_t kRMSx = 20.0*kconv; 
    const Float_t kRMSz = 800.0*kconv;
    Int_t n=0;
    Int_t *tr;
    Int_t ntracks;

    Int_t lay,lad,det;
    fDetTypeRec->GetITSgeom()->GetModuleId(fModule,lay,lad,det);
    Int_t ind=(lad-1)*fDetTypeRec->GetITSgeom()->GetNdetectors(lay)+(det-1);
    Int_t lyr=(lay-1);

    if (GetCrossing(P,N)) {
	//GetCrossingError(dP,dN);
	dP = dN = prob = 0.0; // to remove unused variable warning.
	AliITSRawClusterSSD cnew;
	Int_t nstripsP=clusterP->GetNumOfDigits();
	Int_t nstripsN=clusterN->GetNumOfDigits();
	Float_t signal = 0;
	Float_t dedx = 0;
	//	Float_t meannoiseP=clusterP->GetMeanNoise();
	//Float_t meannoiseN=clusterN->GetMeanNoise();
	//cout<<meannoiseP<<" "<<meannoiseN<<endl;
	if(SigP>SigN) {
	    signal = SigP;
	    dedx = SigP*kADCtoKeV;
	}else{
	    signal = SigN;
	    dedx = SigN*kADCtoKeV;
	} // end if SigP>SigN
     tr = (Int_t*) clusterP->GetTracks(n);
     ntracks = clusterP->GetNTracks();

     //cnew.SetDigitsClusterP(clusterP);
     //cnew.SetDigitsClusterN(clusterN);

     cnew.SetSignalP(SigP);
     cnew.SetSignalN(SigN);
     cnew.SetMultiplicity(nstripsP);
     cnew.SetMultN(nstripsN);
     cnew.SetQErr(TMath::Abs(SigP-SigN));
     cnew.SetNTrack(ntracks);
     //cnew.SetMeanNoiseP(meannoiseP);
     //cnew.SetMeanNoiseN(meannoiseN);
       fDetTypeRec->AddCluster(2,&cnew);
       //fITS->AddCluster(2,&cnew);
       //AliITSRecPoint rnew;
        AliITSRecPoint rnew(fDetTypeRec->GetITSgeom());
        rnew.SetXZ(fModule,P*kconv,N*kconv);
	//rnew.SetX(P*kconv);
	//rnew.SetZ(N*kconv);
     rnew.SetQ(signal);
     rnew.SetdEdX(dedx);
        rnew.SetSigmaDetLocX2( kRMSx* kRMSx);
	//     rnew.SetSigmaX2( kRMSx* kRMSx); 
    rnew.SetSigmaZ2( kRMSz* kRMSz);

         rnew.SetLabel(tr[0],0);
        rnew.SetLabel(tr[1],1);
        rnew.SetLabel(tr[2],2);
        rnew.SetDetectorIndex(ind);
        rnew.SetLayer(lyr);

    //rnew.fTracks[0]=tr[0];
    // rnew.fTracks[1]=tr[1];
    //rnew.fTracks[2]=tr[2];
     //rnew.SetMultP(nstripsP);
     //rnew.SetMultN(nstripsN);
        fDetTypeRec->AddRecPoint(rnew);
	//    fITS->AddRecPoint(rnew);
     return kTRUE;
    } // end if
    return kFALSE;  
}
//______________________________________________________________________
void  AliITSClusterFinderSSD::CalcStepFactor(Float_t Psteo, Float_t Nsteo){
    // calculate the step factor for matching clusters
    // 95 is the pitch, 4000 - dimension along z ?
    //Float_t dz=fSegmentation->Dz();
    Float_t dz=GetSeg()->Dz();

    fSFF = ( (Int_t)  (Psteo*dz/fPitch ) );// +1;
    fSFB = ( (Int_t)  (Nsteo*dz/fPitch ) );// +1;
}
//______________________________________________________________________
AliITSclusterSSD* AliITSClusterFinderSSD::GetPSideCluster(Int_t idx){
    // get P side clusters

    if((idx<0)||(idx>=fNClusterP)){
	printf("AliITSClusterFinderSSD::GetPSideCluster: index out of range\n");
	return 0;
    }else{
	return (AliITSclusterSSD*)((*fClusterP)[idx]);
    } // end if
}
//______________________________________________________________________
AliITSclusterSSD* AliITSClusterFinderSSD::GetNSideCluster(Int_t idx){
    // get N side clusters

    if((idx<0)||(idx>=fNClusterN)){
	printf("AliITSClusterFinderSSD::GetNSideCluster: index out of range\n");
	return 0;
    }else{
	return (AliITSclusterSSD*)((*fClusterN)[idx]);
    } // end if
}
//______________________________________________________________________
AliITSclusterSSD* AliITSClusterFinderSSD::GetCluster(Int_t idx, Bool_t side){
    // Get cluster

    return (side) ? GetPSideCluster(idx) : GetNSideCluster(idx);
}
//______________________________________________________________________
Bool_t AliITSClusterFinderSSD::GetCrossing (Float_t &P, Float_t &N){ 
    // get crossing
    // This function was rivised and changed by Boris Batiounia in March 2001
    Float_t dx = GetSeg()->Dx(); // detector size in x direction, microns
    Float_t dz = GetSeg()->Dz(); // detector size in z direction, microns
    //Float_t dx = fSegmentation->Dx(); // detector size in x direction, microns
    //Float_t dz = fSegmentation->Dz(); // detector size in z direction, microns
    //cout<<dx<<" "<<dz<<endl;
    Float_t xL; // x local coordinate
    Float_t zL; // z local coordinate
    Float_t x;  // x = xL + dx/2
    Float_t z;  // z = zL + dz/2
    Float_t xP; // x coordinate in the P side from the first P strip
    Float_t xN; // x coordinate in the N side from the first N strip
    Float_t stereoP,stereoN;

    //fSegmentation->Angles(stereoP,stereoN);
    GetSeg()->Angles(stereoP,stereoN);

    //cout<<stereoP<<" "<<stereoN<<" "<<P<<" "<<N<<endl;

    //cout<<" P="<<P<<", N="<<N<<endl;

    fTanP=TMath::Tan(stereoP);
    fTanN=TMath::Tan(stereoN);
    Float_t kP = fTanP; // Tangent of 0.0075 mrad
    Float_t kN = fTanN; // Tangent of 0.0275 mrad
    P *= fPitch;
    N *= fPitch; 

    xP = P;
    xN = N;
    //    xP = N;      // change the mistake for the P/N
    //xN = P;      // coordinates correspondence in this function

    z=(xN-xP+dz*kN)/(kP+kN);
    x=xP+kP*z;
    
    xL = x - dx/2;
    zL = z - dz/2;
    P = xL;
    N = zL;  
    //cout<<"P= "<<P<<" , N= "<<N<<" , dx= "<<dx<<endl;

    //cout<<"P="<<P<<", N="<<N<<endl;

    if(TMath::Abs(xL) > dx/2 || TMath::Abs(zL) > dz/2) return kFALSE;
    
    // Check that xL and zL are inside the detector for the 
    // correspondent xP and xN coordinates

    return kTRUE;   
}
//______________________________________________________________________
void AliITSClusterFinderSSD::GetCrossingError(Float_t& dP, Float_t& dN){
    // get crossing error
    Float_t dz, dx;

    dz = TMath::Abs(( dP + dN )*fPitch/(fTanP + fTanN) );
    dx = fPitch*(TMath::Abs(dP*(1 - fTanP/(fTanP + fTanN))) +
		 TMath::Abs(dN *fTanP/(fTanP + fTanN) ));
    dN = dz;
    dP = dx;
}
