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

#include <iostream.h>
#include <TArrayI.h>
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliITSRawCluster.h"
#include "AliITSRecPoint.h"
#include "AliITSMapA1.h"
#include "AliITSClusterFinderSSD.h"
#include "AliITSclusterSSD.h"
#include "AliITSpackageSSD.h"
#include "AliITSsegmentation.h"

const Bool_t AliITSClusterFinderSSD::fgkSIDEP=kTRUE;
const Bool_t AliITSClusterFinderSSD::fgkSIDEN=kFALSE;
const Int_t debug=0;

ClassImp(AliITSClusterFinderSSD)

//____________________________________________________________________
//
//  Constructor
//____________________________________________________________________
//                                     


AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSsegmentation *seg, TClonesArray *digits)   
{
//Standard constructor

    fSegmentation=seg;
    fDigits=digits;
    
    fMap = new AliITSMapA1(fSegmentation,fDigits);  
  
    fITS=(AliITS*)gAlice->GetModule("ITS");
    
    fClusterP =  new TClonesArray ("AliITSclusterSSD",200);    
    fNClusterP =0;   
    
    fClusterN=  new TClonesArray ("AliITSclusterSSD",200);   
    fNClusterN =0; 

    fPackages =  new TClonesArray ("AliITSpackageSSD",200);    //packages  
    fNPackages = 0;

        
    fDigitsIndexP      =  new TArrayI(300);
    fNDigitsP          =  0;
    
    fDigitsIndexN      =  new TArrayI(300);
    fNDigitsN          =  0;
    
    //SetAlpha1(1000);
    //SetAlpha2(1000);
    //SetAlpha3(1000);

    
    fPitch = fSegmentation->Dpx(0);
    Float_t StereoP,StereoN;
    fSegmentation->Angles(StereoP,StereoN);
    fTanP=TMath::Tan(StereoP);
    fTanN=TMath::Tan(StereoN);

    fPNsignalRatio=7./8.;         // warning: hard-wired number

}

//-------------------------------------------------------
AliITSClusterFinderSSD::~AliITSClusterFinderSSD() {
// Default destructor   

    delete fClusterP;
    delete fClusterN;        
    delete fPackages;        
    delete fDigitsIndexP;        
    delete fDigitsIndexN; 
    delete fMap;

}

//-------------------------------------------------------
void AliITSClusterFinderSSD::InitReconstruction()
{
// initialization of the cluster finder

  register Int_t i; //iterator

  for (i=0;i<fNClusterP;i++)
    {
      fClusterP->RemoveAt(i);
    }
  fNClusterP  =0;
  for (i=0;i<fNClusterN;i++)
    {
      fClusterN->RemoveAt(i);
    }
  fNClusterN=0;

  for (i=0;i<fNPackages;i++)
    {
      fPackages->RemoveAt(i);
    }

  fNPackages = 0;
  fNDigitsP=0;
  fNDigitsN=0;

  Float_t StereoP,StereoN;
  fSegmentation->Angles(StereoP,StereoN);

  CalcStepFactor (StereoP,StereoN);

  if (debug) cout<<"fSFF = "<<fSFF<<"  fSFB = "<<fSFB<<"\n";
}


//---------------------------------------------
void AliITSClusterFinderSSD::FindRawClusters() 
{
// This function findes out all clusters belonging to one module
// 1. Zeroes all space after previous module reconstruction
// 2. Finds all neighbouring digits
// 3. If necesery, resolves for each group of neighbouring digits 
//    how many clusters creates it.
// 4. Creates packages  
// 5. Creates clusters 

  InitReconstruction();  //ad. 1
  fMap->FillMap();
  FillDigitsIndex();
  SortDigits();
  FindNeighbouringDigits(); //ad. 2
  //SeparateOverlappedClusters();  //ad. 3
  ClustersToPackages();  //ad. 4
  AliITSRecPoint rnew;
  fMap->ClearMap();
}


//-------------------------------------------------
void AliITSClusterFinderSSD::FindNeighbouringDigits()
{
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
   //If next digit is a neighbour of previous, adds to last cluster this digit
   for (i=1; i<fNDigitsP; i++) {
     //reads new digit
     currentstripNo = ((AliITSdigitSSD*)lDigits[lDigitsIndexP[i]])->
                                                            GetStripNumber(); 
     //check if it is a neighbour of a previous one
     if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexP[i-1]])->GetStripNumber()) 
           ==  (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexP[i];
     else  { 
       //create a new one side cluster
       new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,Digits(),fgkSIDEP); 
       dbuffer[0]=lDigitsIndexP[i];
       dnumber = 1;
     }
   } // end loop over fNDigitsP
   new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,Digits(),fgkSIDEP);


   //process N side 
   //for comments, see above
   dnumber = 1;
   dbuffer[0]=lDigitsIndexN[0];
   //If next digit is a neighbour of previous, adds to last cluster this digit
   for (i=1; i<fNDigitsN; i++) { 
     currentstripNo = ((AliITSdigitSSD*)(lDigits[lDigitsIndexN[i]]))->
                                                            GetStripNumber();
     if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexN[i-1]])->GetStripNumber()) 
           == (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexN[i];
     else {
       new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,Digits(),fgkSIDEN);
       dbuffer[0]=lDigitsIndexN[i];
       dnumber = 1;
     }
   } // end loop over fNDigitsN
   new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,Digits(),fgkSIDEN);
   delete [] dbuffer; 
 
 } // end condition on  NDigits 

 if (debug) cout<<"\n Found clusters: fNClusterP = "<<fNClusterP<<"  fNClusterN ="<<fNClusterN<<"\n";

}  

//--------------------------------------------------------------

void AliITSClusterFinderSSD::SeparateOverlappedClusters()
{
// overlapped clusters separation

  register Int_t i; //iterator

  Float_t  factor=0.75;            // How many percent must be lower signal 
                                   // on the middle one digit
                                   // from its neighbours
  Int_t    signal0;                //signal on the strip before the current one
  Int_t    signal1;                //signal on the current one signal
  Int_t    signal2;                //signal on the strip after the current one
  TArrayI *splitlist;              //  List of splits
  Int_t    numerofsplits=0;        // number of splits
  Int_t    initPsize = fNClusterP; //initial size of the arrays 
  Int_t    initNsize = fNClusterN; //we have to keep it because it will grow 
                                   // in this function and it doasn't make 
                                   // sense to pass through it again

  splitlist = new TArrayI(300);

  for (i=0;i<initPsize;i++)
  {
    if (( ((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits())==1) continue;
    if (( ((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits())==2) continue;
        Int_t nj=(((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits()-1);
        for (Int_t j=1; j<nj; j++)
          {
            signal1=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j);
            signal0=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j-1);
            signal2=((AliITSclusterSSD*)(*fClusterP)[i])->GetDigitSignal(j+1);
            //if signal is less then factor*signal of its neighbours
            if (  (signal1<(factor*signal0)) && (signal1<(factor*signal2)) )
             {                                                               
               (*splitlist)[numerofsplits++]=j;  
	     }
	  } // end loop over number of digits
          //split this cluster if necessary
          if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,fgkSIDEP);
  	  numerofsplits=0;

	  //in signed places (splitlist)
  } // end loop over clusters on Pside

  for (i=0;i<initNsize;i++) {
    if (( ((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits())==1) continue;
    if (( ((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits())==2) continue;
       Int_t nj=(((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits()-1);
       for (Int_t j=1; j<nj; j++)
          {
            signal1=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j);
            signal0=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j-1);
            signal2=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j+1);
            //if signal is less then factor*signal of its neighbours
            if (  (signal1<(factor*signal0)) && (signal1<(factor*signal2)) ) 
               (*splitlist)[numerofsplits++]=j;  
          } // end loop over number of digits 
          //split this cluster into more clusters
          if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,fgkSIDEN);
	  numerofsplits=0;                                                              //in signed places (splitlist)
  } // end loop over clusters on Nside

  delete splitlist;
}

//-------------------------------------------------------
void AliITSClusterFinderSSD::SplitCluster(TArrayI *list, Int_t nsplits, Int_t index, Bool_t side)
{
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
  Int_t    NN;
  // side true means P side
  if (side) {
     curentcluster =((AliITSclusterSSD*)((*fClusterP)[index])) ;
     for (i = nsplits; i>0 ;i--) {  
         NN=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
         new ((*fClusterP)[fNClusterP]) AliITSclusterSSD(NN,tmpdigits,Digits(),side);
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
         NN=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
	 new ((*fClusterN)[fNClusterN]) AliITSclusterSSD(NN,tmpdigits,Digits(),side);
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


//-------------------------------------------------
Int_t AliITSClusterFinderSSD::SortDigitsP(Int_t start, Int_t end)
{
// sort digits on the P side

  
  Int_t right;
  Int_t left;
  if (start != (end - 1) ) 
    {
      left=this->SortDigitsP(start,(start+end)/2);
      right=this->SortDigitsP((start+end)/2,end);  
      return (left || right);
    }
  else
   { 
    left =  ((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexP)[start]]))->GetStripNumber();
    right= ((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexP)[end]]))->GetStripNumber();  
    if( left > right )
     {
       Int_t tmp = (*fDigitsIndexP)[start];
       (*fDigitsIndexP)[start]=(*fDigitsIndexP)[end];
       (*fDigitsIndexP)[end]=tmp;
       return 1;
     }
    else return 0;
   }
}


//--------------------------------------------------

Int_t AliITSClusterFinderSSD::SortDigitsN(Int_t start, Int_t end)
{
// sort digits on the N side

  Int_t right;
  Int_t left;
  if (start != (end - 1)) 
    {
      left=this->SortDigitsN(start,(start+end)/2);
      right=this->SortDigitsN((start+end)/2,end);  
      return (left || right);
    }
  else 
   {
    left =((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexN)[start]]))->GetStripNumber(); 
    right=((AliITSdigitSSD*)((*(Digits()))[(*fDigitsIndexN)[end]]))->GetStripNumber();
    if ( left > right )
      {
        Int_t tmp = (*fDigitsIndexN)[start];
        (*fDigitsIndexN)[start]=(*fDigitsIndexN)[end];
        (*fDigitsIndexN)[end]=tmp;
        return 1;
      }else return 0;
   }  
}


//------------------------------------------------
void AliITSClusterFinderSSD::FillDigitsIndex()
{
 //Fill the indexes of the clusters belonging to a given ITS module

 Int_t PNs=0, NNs=0;
 Int_t tmp,bit,k;
 Int_t N;
 Int_t i;
 
 N = fDigits->GetEntriesFast();

 Int_t* PSidx = new Int_t [N*sizeof(Int_t)];
 Int_t* NSidx = new Int_t [N*sizeof(Int_t)]; 
 if (fDigitsIndexP==NULL) fDigitsIndexP = new TArrayI(N);
 if (fDigitsIndexN==NULL) fDigitsIndexN = new TArrayI(N);
 
 AliITSdigitSSD *dig; 
 
 for ( i = 0 ; i< N; i++ ) {
      dig = (AliITSdigitSSD*)GetDigit(i);
      if(dig->IsSideP()) { 
           bit=1;
           tmp=dig->GetStripNumber();
	   // I find this totally unnecessary - it's just a 
	   // CPU consuming double check
           for( k=0;k<PNs;k++)
            {
             if (tmp==PSidx[k])
              {
                 if (debug) cout<<"Such a digit exists \n";
                 bit=0;
              }
           }   
	   // end comment 
        if(bit) {
            fDigitsIndexP->AddAt(i,fNDigitsP++);
            PSidx[PNs++]=tmp;
	}
      } else {
         bit=1;
         tmp=dig->GetStripNumber();
	 // same as above
         for( k=0;k<NNs;k++)
          {
           if (tmp==NSidx[k])
            {
             if (debug) cout<<"Such a digit exists \n";
             bit=0;
            }
          } 
	 // end comment
         if (bit) {
          fDigitsIndexN->AddAt(i,fNDigitsN++);
          NSidx[NNs++] =tmp;
	 }
      }
 }
   
 delete [] PSidx;
 delete [] NSidx;

 if (debug) cout<<"Digits :  P = "<<fNDigitsP<<"   N = "<<fNDigitsN<<endl;

}


//-------------------------------------------

void AliITSClusterFinderSSD::SortDigits()
{
// sort digits

  Int_t i;
  if(fNDigitsP>1) 
  for (i=0;i<fNDigitsP-1;i++)
  if (SortDigitsP(0,(fNDigitsP-1-i))==0) break;

  if(fNDigitsN>1) 
    for (i=0;i<fNDigitsN-1;i++)
  if(SortDigitsN(0,(fNDigitsN-1-i))==0) break;
}



//----------------------------------------------
void AliITSClusterFinderSSD::FillClIndexArrays(Int_t* arrayP, Int_t *arrayN)
{
// fill cluster index array

  register Int_t i;
  for (i=0; i<fNClusterP;i++)
    {
      arrayP[i]=i;
    }
  for (i=0; i<fNClusterN;i++)
    {
      arrayN[i]=i;
    }
}


//------------------------------------------------------
void AliITSClusterFinderSSD::SortClusters(Int_t* arrayP, Int_t *arrayN)
{
// sort clusters

  Int_t i;
  if(fNClusterP>1) 
    for (i=0;i<fNClusterP-1;i++)
      if (SortClustersP(0,(fNClusterP-1),arrayP)==0)  break;
    
    
  if(fNClusterN>1) 
    for (i=0;i<fNClusterN-1;i++)
      if (SortClustersN(0,(fNClusterN-1),arrayN)==0)  break;

}



//---------------------------------------------------
Int_t AliITSClusterFinderSSD::SortClustersP(Int_t start, Int_t end, Int_t *array)
{
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
  }

 
}



//-------------------------------------------------------
Int_t AliITSClusterFinderSSD::SortClustersN(Int_t start, Int_t end, Int_t *array)
{
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
    } 

}
    


//-------------------------------------------------------
void AliITSClusterFinderSSD::ClustersToPackages()
{  
// fill packages   
    
  Int_t *oneSclP = new Int_t[fNClusterP]; //I want to have sorted 1 S clusters
  Int_t *oneSclN = new Int_t[fNClusterN]; //I can not sort it in TClonesArray
                                          //so, I create table of indexes and 
                                          //sort it
                                          //I do not use TArrayI on purpose
                                          // MB: well, that's not true that one
                                          //cannot sort objects in TClonesArray
  
  //AliITSpackageSSD *currentpkg;
  AliITSclusterSSD *currentP;
  AliITSclusterSSD *currentN;
  Int_t j1, j2;    

//Fills in One Side Clusters Index Array
  FillClIndexArrays(oneSclP,oneSclN); 
//Sorts filled Arrays
  SortClusters(oneSclP,oneSclN);                   


  fNPackages=1;      
  new ((*fPackages)[0]) AliITSpackageSSD(fClusterP,fClusterN);
  //currentpkg = (AliITSpackageSSD*)((*fPackages)[0]);

  // Take all pairs of recpoints x coordinates in both sides  
  // to calculate z coordinates of the recpoints
  for (j1=0;j1<fNClusterP;j1++) {  
      currentP = GetPSideCluster(oneSclP[j1]);
      Double_t xP = currentP->GetPosition();
      //Int_t NxP = currentP->GetNumOfDigits();
      //cout<<"ClusterToP:New NxP,xP ="<<NxP<<","<<xP<<endl;
      Float_t signalP = currentP->GetTotalSignal();
    for (j2=0;j2<fNClusterN;j2++) {  
      currentN = GetNSideCluster(oneSclN[j2]);
      Double_t xN = currentN->GetPosition();
      //Int_t NxN = currentN->GetNumOfDigits();
      //cout<<"ClusterToP:New NxN,xN ="<<NxN<<","<<xN<<endl;
      Float_t signalN = currentN->GetTotalSignal();
      //cout<<"ClusterToP: signalP,signalN ="<<signalP<<","<<signalN<<endl;
      CreateNewRecPoint(xP,1,xN,1,signalP,signalN,currentP, currentN, 0.75);

    }
  }

   delete oneSclP;
   delete oneSclN;

}


//------------------------------------------------------
Bool_t AliITSClusterFinderSSD::    
CreateNewRecPoint(Float_t P, Float_t dP, Float_t N, Float_t dN,
                  Float_t SigP,Float_t SigN, 
                  AliITSclusterSSD *clusterP, AliITSclusterSSD *clusterN,
                  Stat_t prob)
{
// create the recpoints
  const Float_t kADCtoKeV = 2.16; 
  // 50 ADC units -> 30000 e-h pairs; 1e-h pair -> 3.6e-3 KeV;
  // 1 ADC unit -> (30000/50)*3.6e-3 = 2.16 KeV 
  const Float_t kconv = 1.0e-4; 

  const Float_t kRMSx = 20.0*kconv; 
  const Float_t kRMSz = 800.0*kconv; 

  Int_t stripP, stripN;
  Int_t n=0;
  Int_t *tr;
  if (GetCrossing(P,N)) {
  
     GetCrossingError(dP,dN);
     AliITSRawClusterSSD cnew;
     Int_t nstripsP=clusterP->GetNumOfDigits();
     Int_t nstripsN=clusterN->GetNumOfDigits();
     //Float_t signalP=clusterP->GetTotalSignal();
     //Float_t signalN=clusterN->GetTotalSignal();

     Float_t signal = 0;
     Float_t dedx = 0;
     if(SigP>SigN) {
       signal = SigP;
       dedx = SigP*kADCtoKeV;
     }else{
       signal = SigN;
       dedx = SigN*kADCtoKeV;
     }	       
     Float_t signalCut = TMath::Abs((SigP-SigN)/signal);

  //cnew.fSignalP=signalP;
  //cnew.fSignalN=signalN;
     cnew.fSignalP=SigP;
     cnew.fSignalN=SigN;
     cnew.fMultiplicity=nstripsP;
     cnew.fMultiplicityN=nstripsN;
     cnew.fQErr=signalCut;

     //fITS->AddCluster(2,&cnew);

     if(signalCut<0.18) fITS->AddCluster(2,&cnew);
     // the cut of the signals difference in P and N sides to
     // remove the "ghosts"

     fSegmentation->GetPadIxz(P,N,stripP,stripN);

     //cout<<"NewRec: P,N microns,stripP,stripN ="<<P<<","<<N<<","<<stripP<<","<<stripN<<endl;

     tr = (Int_t*) clusterP->GetTracks(n);
     //cout<<"!!!! ntr,tr0,tr1,tr2 ="<<n<<","<<tr[0]<<","<<tr[1]<<","<<tr[2]<<endl; 
     //cout<<"NewRec:AddRecP: xL,zL cm ="<<P*kconv<<","<<N*kconv<<endl;

     AliITSRecPoint rnew;
     rnew.SetX(P*kconv);
     rnew.SetZ(N*kconv);
     rnew.SetQ(signal);
     rnew.SetdEdX(dedx);
     rnew.SetSigmaX2( kRMSx* kRMSx); 
     rnew.SetSigmaZ2( kRMSz* kRMSz);
     //rnew.SetSigmaX2(1);   //This has to be verified
     //rnew.SetSigmaZ2(1);   //This has to be verified
     rnew.fTracks[0]=tr[0];
     rnew.fTracks[1]=tr[1];
     rnew.fTracks[2]=tr[2];

     //fITS->AddRecPoint(rnew);

     if(signalCut<0.18) fITS->AddRecPoint(rnew);
     // the cut of the signals difference in P and N sides to
     // remove the "ghosts"

     return kTRUE;
  }
     return kFALSE;  
}


//------------------------------------------------------
void  AliITSClusterFinderSSD::CalcStepFactor(Float_t Psteo, Float_t Nsteo )
{
// calculate the step factor for matching clusters


  // 95 is the pitch, 4000 - dimension along z ?


  Float_t dz=fSegmentation->Dz();

  fSFF = ( (Int_t)  (Psteo*dz/fPitch ) );// +1;
  fSFB = ( (Int_t)  (Nsteo*dz/fPitch ) );// +1;

}


//-----------------------------------------------------------
AliITSclusterSSD* AliITSClusterFinderSSD::GetPSideCluster(Int_t idx)
{
// get P side clusters




  if((idx<0)||(idx>=fNClusterP))
    {
      printf("AliITSClusterFinderSSD::GetPSideCluster  : index out of range\n");
      return 0;
    }
  else
    {
      return (AliITSclusterSSD*)((*fClusterP)[idx]);
    }
}

//-------------------------------------------------------
AliITSclusterSSD* AliITSClusterFinderSSD::GetNSideCluster(Int_t idx)
{
// get N side clusters

  if((idx<0)||(idx>=fNClusterN))
    {
      printf("AliITSClusterFinderSSD::GetNSideCluster  : index out of range\n");
      return 0;
    }
  else
    {
      return (AliITSclusterSSD*)((*fClusterN)[idx]);
    }
}

//--------------------------------------------------------
AliITSclusterSSD* AliITSClusterFinderSSD::GetCluster(Int_t idx, Bool_t side)
{
// Get cluster

  return (side) ? GetPSideCluster(idx) : GetNSideCluster(idx);
}

//_______________________________________________________________________

Bool_t AliITSClusterFinderSSD::GetCrossing (Float_t &P, Float_t &N) 
{ 
// get crossing 

  Float_t Dx = fSegmentation->Dx(); // detector size in x direction, microns
  Float_t Dz = fSegmentation->Dz(); // detector size in z direction, microns

  Float_t xL; // x local coordinate
  Float_t zL; // z local coordinate
  Float_t x;  // x = xL + Dx/2
  Float_t z;  // z = zL + Dz/2
  Float_t xP; // x coordinate in the P side from the first P strip
  Float_t xN; // x coordinate in the N side from the first N strip
  Float_t kP = fTanP; // Tangent of 0.0075 mrad
  Float_t kN = fTanN; // Tangent of 0.0275 mrad

  //cout<<"1 GetCros: strP,srtN,fPitch,Dx ="<<P<<","<<N<<","<<fPitch<<","<<fSegmentation->Dx()<<endl;


    P *= fPitch;
    N *= fPitch; 
    xP = N;      // change the mistake for the P/N
    xN = P;      // coordinates correspondence in this function
 
    //cout<<"2 GetCross: xP,xN,fTanP,fTanN ="<<xP<<","<<xN<<","<<fTanP<<","<<fTanN<<endl;     

    /*
    Float_t x = xP + 0.21428*(1100-xP+xN);
    Float_t z = 31428.571 - (xP-xN)/0.035;    
    */

    x = xP + kP*(Dz*kN-xP+xN)/(kP+kN);
    z = (Dz*kN-xP+xN)/(kP+kN); 
    xL = x - Dx/2;
    zL = z - Dz/2;

    //cout<<"3 GetCross: x,z,xL,zL ="<<x<<","<<z<<","<<xL<<","<<zL<<endl;     

    P = xL;
    N = zL;  

    if(TMath::Abs(xL) > Dx/2 || TMath::Abs(zL) > Dz/2) return kFALSE;
    // Check that xL and zL are inside the detector for the 
    // correspondent xP and xN coordinates

    return kTRUE;   
}


//_________________________________________________________________________

void AliITSClusterFinderSSD::GetCrossingError(Float_t& dP, Float_t& dN)
{
// get crossing error

  Float_t dz, dx;
  
  dz = TMath::Abs(( dP + dN )*fPitch/(fTanP + fTanN) );
  dx = fPitch*(TMath::Abs(dP*(1 - fTanP/(fTanP + fTanN))) +
               TMath::Abs(dN *fTanP/(fTanP + fTanN) ));
  
  dN = dz;
  dP = dx;
}







