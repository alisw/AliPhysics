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

/*

Adding rekonstruction facilities
Piotr Krzysztof Skowronski 
December 1999.
*/

/*
15 -18 V 2000
Eroor counting routines
Automatic combination routines improved (traps)

*/

#include "AliRun.h"
#include "AliITSClusterFinderSSD.h"
 

const Int_t debug=0;

ClassImp(AliITSClusterFinderSSD)

//____________________________________________________________________
//
//  Constructor
//____________________________________________________________________
//                                     


AliITSClusterFinderSSD::AliITSClusterFinderSSD(AliITSsegmentation *seg, TClonesArray *digits, TClonesArray *recp)   
{

    fSegmentation=seg;
    fDigits=digits;
    fRecPoints=recp;
    
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
    
    SetAlpha1(1000);
    SetAlpha2(1000);
    SetAlpha3(1000);

    
    fPitch = fSegmentation->Dpx(0);
    Float_t StereoP,StereoN;
    fSegmentation->Angles(StereoP,StereoN);
    fTanP=TMath::Tan(StereoP);
    fTanN=TMath::Tan(StereoN);

}

//-------------------------------------------------------
AliITSClusterFinderSSD::~AliITSClusterFinderSSD() {
   
    delete fClusterP;
    delete fClusterN;        
    delete fPackages;        
    delete fDigitsIndexP;        
    delete fDigitsIndexN; 
}

//-------------------------------------------------------
void AliITSClusterFinderSSD::InitReconstruction()
{

  register Int_t i; //iterator

  for(i=0;i<fNClusterP;i++)
    {
      fClusterP->RemoveAt(i);
    }
  fNClusterP  =0;
  for(i=0;i<fNClusterN;i++)
    {
      fClusterN->RemoveAt(i);
    }
  fNClusterN=0;

  for(i=0;i<fNPackages;i++)
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



//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

// This function findes out all clusters belonging to one module
// 1. Zeroes all space after previous module reconstruction
// 2. Finds all neighbouring digits
// 3. If necesery, resolves for each group of neighbouring digits 
//    how many clusters creates it.
// 4. Creates packages  
// 5. Creates clusters 

  InitReconstruction();  //ad. 1
  FillDigitsIndex();
  SortDigits();
  FindNeighbouringDigits(); //ad. 2
  SeparateOverlappedClusters();  //ad. 3
  ClustersToPackages();  //ad. 4
  ConsumeClusters();
  PackagesToPoints();   //ad. 5
  ReconstructNotConsumedClusters();
 
}


//-------------------------------------------------
void AliITSClusterFinderSSD::FindNeighbouringDigits()
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

 register Int_t i;

    //If there are any digits on this side, create 1st Cluster,
    // add to it this digit, and increment number of clusters

 if ((fNDigitsP>0 )  && (fNDigitsN > 0 )) {     

   Int_t currentstripNo;
   Int_t *dbuffer = new Int_t [300];   //buffer for strip numbers
   Int_t dnumber;    //curent number of digits in buffer
   TArrayI &lDigitsIndexP = *fDigitsIndexP;
   TArrayI &lDigitsIndexN = *fDigitsIndexN;
   TObjArray &lDigits=*fDigits;
   TClonesArray &lClusterP = *fClusterP;
   TClonesArray &lClusterN = *fClusterN;
  
   //process P side 
   dnumber = 1;
   dbuffer[0]=lDigitsIndexP[0];
   //If next digit is a neighbour of previous, adds to last cluster this digit
   for(i=1; i<fNDigitsP; i++) {
     //reads new digit
     currentstripNo = ((AliITSdigitSSD*)lDigits[lDigitsIndexP[i]])->
                                                            GetStripNumber(); 
     //check if it is a neighbour of a previous one
     if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexP[i-1]])->GetStripNumber()) 
           ==  (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexP[i];
     else  { 
       //create a new one side cluster
       new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,fDigits,SIDEP); 
       dbuffer[0]=lDigitsIndexP[i];
       dnumber = 1;
     }
   } // end loop over fNDigitsP
   new(lClusterP[fNClusterP++]) AliITSclusterSSD(dnumber,dbuffer,fDigits,SIDEP);


   //process N side 
   //for comments, see above
   dnumber = 1;
   dbuffer[0]=lDigitsIndexN[0];
   //If next digit is a neighbour of previous, adds to last cluster this digit
   for(i=1; i<fNDigitsN; i++) { 
     currentstripNo = ((AliITSdigitSSD*)(lDigits[lDigitsIndexN[i]]))->
                                                            GetStripNumber();
     if ( (((AliITSdigitSSD*)lDigits[lDigitsIndexN[i-1]])->GetStripNumber()) 
           == (currentstripNo-1) ) dbuffer[dnumber++]=lDigitsIndexN[i];
     else {
       new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,fDigits,SIDEN);
       dbuffer[0]=lDigitsIndexN[i];
       dnumber = 1;
     }
   } // end loop over fNDigitsN
   new(lClusterN[fNClusterN++]) AliITSclusterSSD(dnumber,dbuffer,fDigits,SIDEN);
   delete [] dbuffer; 
 
 } // end condition on  NDigits 

 if (debug) cout<<"\n Found clusters: fNClusterP = "<<fNClusterP<<"  fNClusterN ="<<fNClusterN<<"\n";

}  
/**********************************************************************/


void AliITSClusterFinderSSD::SeparateOverlappedClusters()
{
//************************************************
//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

  register Int_t i,j; //iterator

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

  for(i=0;i<initPsize;i++)
  {
    if (( ((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits())==1) continue;
    if (( ((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits())==2) continue;
        Int_t nj=(((AliITSclusterSSD*)(*fClusterP)[i])->GetNumOfDigits()-1);
        for(j=1; j<nj; j++)
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
          if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,SIDEP); 
	  numerofsplits=0;

	  //in signed places (splitlist)
  } // end loop over clusters on Pside

  for(i=0;i<initNsize;i++) {
    if (( ((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits())==1) continue;
    if (( ((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits())==2) continue;
       Int_t nj=(((AliITSclusterSSD*)(*fClusterN)[i])->GetNumOfDigits()-1);
       for(j=1; j<nj; j++)
          {
            signal1=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j);
            signal0=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j-1);
            signal2=((AliITSclusterSSD*)(*fClusterN)[i])->GetDigitSignal(j+1);
            //if signal is less then factor*signal of its neighbours
            if (  (signal1<(factor*signal0)) && (signal1<(factor*signal2)) ) 
               (*splitlist)[numerofsplits++]=j;  
          } // end loop over number of digits 
          //split this cluster into more clusters
          if(numerofsplits>0) SplitCluster(splitlist,numerofsplits,i,SIDEN);
	  numerofsplits=0;                                                              //in signed places (splitlist)
  } // end loop over clusters on Nside

  delete splitlist;
}

//-------------------------------------------------------
void AliITSClusterFinderSSD::SplitCluster(TArrayI *list, Int_t nsplits, Int_t index, Bool_t side)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

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
     for(i = nsplits; i>0 ;i--) {  
         NN=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
         new ((*fClusterP)[fNClusterP]) AliITSclusterSSD(NN,tmpdigits,fDigits,side);
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
     for(i = nsplits; i>0 ;i--) {  
         NN=curentcluster->SplitCluster((*list)[(i-1)],tmpdigits);
	 new ((*fClusterN)[fNClusterN]) AliITSclusterSSD(NN,tmpdigits,fDigits,side);
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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl
  
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
    left =  ((AliITSdigitSSD*)((*fDigits)[(*fDigitsIndexP)[start]]))->GetStripNumber();
    right= ((AliITSdigitSSD*)((*fDigits)[(*fDigitsIndexP)[end]]))->GetStripNumber();  
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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

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
    left =((AliITSdigitSSD*)((*fDigits)[(*fDigitsIndexN)[start]]))->GetStripNumber(); 
    right=((AliITSdigitSSD*)((*fDigits)[(*fDigitsIndexN)[end]]))->GetStripNumber();
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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

 //Fill the indexes of the clusters belonging to a given ITS module
 //Created by Piotr K. Skowronski, August 7 1999

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
 
 for(i = 0 ; i< N; i++ ) {
      dig=(AliITSdigitSSD*)fDigits->UncheckedAt(i);
      if(dig->IsSideP()) { 
           bit=1;
           tmp=dig->GetStripNumber();
	   // I find this totally unnecessary - it's just a 
	   // CPU consuming double check
           for(k=0;k<PNs;k++)
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
         for(k=0;k<NNs;k++)
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
   

 if (debug) cout<<"Digits :  P = "<<fNDigitsP<<"   N = "<<fNDigitsN<<endl;

}


//-------------------------------------------

void AliITSClusterFinderSSD::SortDigits()
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl



  Int_t i;
  if(fNDigitsP>1) 
  for(i=0;i<fNDigitsP-1;i++)
  if (SortDigitsP(0,(fNDigitsP-1-i))==0) break;

  if(fNDigitsN>1) 
    for(i=0;i<fNDigitsN-1;i++)
  if(SortDigitsN(0,(fNDigitsN-1-i))==0) break;
}



//----------------------------------------------
void AliITSClusterFinderSSD::FillClIndexArrays(Int_t* arrayP, Int_t *arrayN)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl


  register Int_t i;
  for(i=0; i<fNClusterP;i++)
    {
      arrayP[i]=i;
    }
  for(i=0; i<fNClusterN;i++)
    {
      arrayN[i]=i;
    }
}


//------------------------------------------------------
void AliITSClusterFinderSSD::SortClusters(Int_t* arrayP, Int_t *arrayN)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl


  Int_t i;
  if(fNClusterP>1) 
    for(i=0;i<fNClusterP-1;i++)
      if (SortClustersP(0,(fNClusterP-1),arrayP)==0)  break;
    
    
  if(fNClusterN>1) 
    for(i=0;i<fNClusterN-1;i++)
      if (SortClustersN(0,(fNClusterN-1),arrayN)==0)  break;

}



//---------------------------------------------------
Int_t AliITSClusterFinderSSD::SortClustersP(Int_t start, Int_t end, Int_t *array)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl



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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

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



//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl
  
    
  Int_t *oneSclP = new Int_t[fNClusterP]; //I want to have sorted 1 S clusters
  Int_t *oneSclN = new Int_t[fNClusterN]; //I can not sort it in TClonesArray
                                          //so, I create table of indexes and 
                                          //sort it
                                          //I do not use TArrayI on purpose
                                          // MB: well, that's not true that one
                                          //cannot sort objects in TClonesArray
  
  Int_t i,k;    //iterator
  AliITSpackageSSD *currentpkg;
  AliITSclusterSSD *currentP;
  AliITSclusterSSD *currentN;
  AliITSclusterSSD *fcurrN;

  Int_t     lastNclIndex=0;          //index of last N sidecluster
  Int_t     Nclpack=0;               //number of P clusters in current package
  Int_t     Pclpack=0;               //number of N clusters in current package
  Int_t     firstStripP;
  Int_t     lastStripP;
  Int_t     firstStripN;
  Int_t     lastStripN;
  Int_t     tmplastNclIndex;

//Fills in One Side Clusters Index Array
  FillClIndexArrays(oneSclP,oneSclN); 
//Sorts filled Arrays
  SortClusters(oneSclP,oneSclN);                   


  fNPackages=1;      
  new ((*fPackages)[0]) AliITSpackageSSD(fClusterP,fClusterN);
  currentpkg = (AliITSpackageSSD*)((*fPackages)[0]);
  
  //main loop on sorted P side clusters 
  for(i=0;i<fNClusterP;i++) {  
      //Take new P side cluster
      currentP = GetPSideCluster(oneSclP[i]);
      currentN = GetNSideCluster(oneSclN[lastNclIndex]);
      // take a new P cluster or not ?
      if(IsCrossing(currentP,currentN)) {
	  // don't take a new P cluster
          Pclpack++;
          currentP->AddCross(oneSclN[lastNclIndex]);
          currentN->AddCross(oneSclP[i]);
          currentpkg->AddPSideCluster(oneSclP[i]);
          if (Nclpack==0) {
              currentpkg->AddNSideCluster(oneSclN[lastNclIndex]);          
              Nclpack++;
	  }  
	  //check if previous N side clusters crosses with it too    
          for(k=1;k<fSFB+1;k++) {
	      //check if we are not out of array 
              if ((lastNclIndex-k)>-1) {
                  fcurrN = GetNSideCluster( oneSclN[lastNclIndex-k] );
                  if( IsCrossing(currentP,fcurrN) ) {
		      currentP->AddCross(oneSclN[lastNclIndex-k]);
                      fcurrN->AddCross(oneSclP[i]);      
		 } else break; //There is no sense to check rest of clusters
	      } else break;
	  }
          tmplastNclIndex=lastNclIndex;  
	  //Check if next N side clusters crosses with it too
          for(k=1;k<fSFF+1;k++) {
	      //check if we are not out of array 
              if ((tmplastNclIndex+k)<fNClusterN) {
                  fcurrN = GetNSideCluster( oneSclN[tmplastNclIndex+k] );
                  if(IsCrossing(currentP,fcurrN) ) {
                     lastNclIndex++;
                     fcurrN->AddCross(oneSclP[i]);
                     currentP->AddCross(oneSclN[tmplastNclIndex+k]);
                     currentpkg->AddNSideCluster(oneSclN[lastNclIndex]);
                     Nclpack++;
		  }else break;
	      } else break;
	  }
      
	  // end of package 
	  //if( IsCrossing )
      } else { 
         lastStripP  = currentP->GetLastDigitStripNo();
         lastStripN  = currentN->GetLastDigitStripNo();
         
         if((lastStripN-fSFF) < lastStripP ) {
	    //Take new PCluster
            if((Nclpack>0)&& (Pclpack>0)) {
              new ((*fPackages)[fNPackages]) AliITSpackageSSD(fClusterP,fClusterN);
              currentpkg = (AliITSpackageSSD*)((*fPackages)[fNPackages++]);
	    }
              
	    //so we have to take another cluster on side N and check it
	    //we have to check until taken cluster's last strip will 
	    //be > last of this cluster(P)
	    //so it means that there is no corresponding cluster on side N 
	    //to this cluster (P)
	    //so we have to continue main loop with new "lastNclIndex"
         
	    Nclpack=0;     
	    Pclpack=0;   
	    //We are not sure that next N cluter will cross with this P cluster
	    //There might one, or more, clusters that do not have any 
	    //corresponding cluster on the other side	   
	    for(;;) {     
	       lastNclIndex++;
	       //Check if we are not out of array
	       if (lastNclIndex<fNClusterN) {
		  currentN = GetNSideCluster(oneSclN[lastNclIndex]);
		  if ( IsCrossing(currentP, currentN) ){
		     Nclpack++;       
		     Pclpack++;
		     currentP->AddCross(oneSclN[lastNclIndex]);
		     currentN->AddCross(oneSclP[i]);
		     currentpkg->AddNSideCluster(oneSclN[lastNclIndex]);
		     currentpkg->AddPSideCluster(oneSclP[i]);          
		     //Check, if next N side clusters crosses with it too
		     tmplastNclIndex=lastNclIndex;
		     for(k=1;k<fSFF+1;k++) {
		       //Check if we are not out of array
		       if ((tmplastNclIndex+k)<fNClusterN) {
			  fcurrN = GetNSideCluster(oneSclN[tmplastNclIndex+k]);
			  if( IsCrossing(currentP, fcurrN) ) {
			     Nclpack++;
			     lastNclIndex++;
			     currentP->AddCross(oneSclN[tmplastNclIndex+k]);
			     fcurrN->AddCross(oneSclP[i]);
			     currentpkg->AddNSideCluster(oneSclN[tmplastNclIndex+k]);
			  }else break;
		       } else break; //we are out of array
		     } // end loop 
		     break;
		  } else {
		     firstStripP = currentP->GetFirstDigitStripNo();
		     firstStripN = currentN->GetFirstDigitStripNo();
		     if((firstStripN+fSFB)>=firstStripP) break;
		     else continue;
		  }  
	       } else goto EndOfFunction;
	    } // end for(;;)
	 }  else  //EndOfPackage
	   {
	     continue;     
	   } //else EndOfPackage
       
      }//else IsCrossing
	    
  }//main for

 EndOfFunction:
  if ((Nclpack<1) || (Pclpack<1)) fNPackages--;

   delete oneSclP;
   delete oneSclN;

}



//-----------------------------------------------
void AliITSClusterFinderSSD::PackagesToPoints()
{
 
 register Int_t i;
 AliITSpackageSSD *currentpkg;
 Int_t NumNcl;
 Int_t NumPcl;
 Int_t clusterIndex;
 Bool_t clusterSide;

 for(i=0;i<fNPackages;i++) {
    //get pointer to next package
    currentpkg = (AliITSpackageSSD*)((*fPackages)[i]); 
    NumNcl = currentpkg->GetNumOfClustersN();
    NumPcl = currentpkg->GetNumOfClustersP();
    if(debug) cout<<"New Package\nNumPcl ="<<NumPcl<<" NumNcl ="<<NumNcl<<"\n";

    while(((NumPcl>=2)&&(NumNcl>2))||((NumPcl>2)&&(NumNcl>=2))) {  
        //This is case of "big" pacakage
        //conditions see 3 lines above  
        //if in the big package exists cluster with one cross 
	//we can reconstruct this point without any geometrical ambiguities
	if(currentpkg->GetClusterWithOneCross(clusterIndex, clusterSide) )
        {
            ResolveClusterWithOneCross(currentpkg, clusterIndex, clusterSide);
        } else if (clusterIndex == -2) {
	    NumPcl = 0;
	    NumNcl = 0;
	    break;
        } else {
	    if ( (NumNcl==NumPcl) && (NumPcl<10)) {
               //if there is no cluster with one cross 
	       //we can resolve whole package at once 
	       //by finding best combination
	       //we can make combinations when NumNcl==NumPcl
	       if (ResolvePackageBestCombin(currentpkg)) { 
                   //sometimes creating combinations fail, 
		   //but it happens very rarely
		   NumPcl = 0;
		   NumNcl = 0;
		   break;
	       } else ResolveOneBestMatchingPoint(currentpkg);
	    } else {
	       ResolveOneBestMatchingPoint(currentpkg);
	    }
        }
	NumNcl = currentpkg->GetNumOfClustersN();
	NumPcl = currentpkg->GetNumOfClustersP();
    } // end while 
    if ((NumPcl==1)&&(NumNcl==1)) {
      ResolveSimplePackage(currentpkg);
       continue;
    }
    if (NumPcl==1) {
      ResolvePackageWithOnePSideCluster(currentpkg);
      continue;
    }
    if (NumNcl==1) {
      ResolvePackageWithOneNSideCluster(currentpkg);
      continue;
    }
    if ((NumPcl==2)&&(NumNcl==2)) {
      ResolveTwoForTwoPackage(currentpkg);
      continue; 
    }
    if ((NumPcl==0)&&(NumNcl>0)) { }
    if ((NumNcl==0)&&(NumPcl>0)) { }

  } // end loop over fNPackages


}


//----------------------------------------------------------

void AliITSClusterFinderSSD::
ResolveClusterWithOneCross(AliITSpackageSSD *currentpkg, Int_t clusterIndex, Bool_t clSide)
{

   if (clSide == SIDEP) ResolvePClusterWithOneCross(currentpkg,clusterIndex);
   else  ResolveNClusterWithOneCross(currentpkg,clusterIndex);

}


//---------------------------------------------------------
void AliITSClusterFinderSSD::
ResolvePClusterWithOneCross(AliITSpackageSSD *pkg, Int_t clusterIndex)
{

/*
There is cluster (side P) which crosses with only one cluster on the other 
side (N)

ie:
    \    /  \/   /
     \  /   /\  /
      \/   /  \/
      /\  /   /\
     /  \/   /  \    .....
    /   /\  /    \
   /   /  \/      \
  /   /   /\       \
 /   /   /  \       \
*/


  AliITSclusterSSD * clusterP;
  AliITSclusterSSD * clusterN;

  Float_t posClusterP;           //Cluster P position in strip coordinates
  Float_t posClusterN;           //Cluster N position in strip coordinates
  
  Float_t posErrorClusterP;
  Float_t posErrorClusterN;
  
  Float_t sigClusterP;
  Float_t sigClusterN;

  Float_t sigErrorClusterP;
  Float_t sigErrorClusterN;

  Float_t ps;
  Float_t ns;
  
  Float_t Chicomb;
  Int_t clusterIdx;

  if (debug) cout<<"ResolvePClusterWithOneCross\n";

  clusterP=pkg->GetPSideCluster(clusterIndex);
  posClusterP = GetClusterZ(clusterP);
  posErrorClusterP = clusterP->GetPositionError();
  sigClusterP = ps = clusterP->GetTotalSignal();
  sigErrorClusterP = clusterP->GetTotalSignalError(); 
  //carefully, it returns index in TClonesArray
  
  clusterIdx = clusterP->GetCross(0);
  clusterN=this->GetNSideCluster(clusterIdx);
  ns = clusterN->GetTotalSignal();
  posClusterN = GetClusterZ(clusterN);
  posErrorClusterN = clusterN->GetPositionError();
  pkg->DelCluster(clusterIndex,SIDEP);
  sigClusterN = ps/PNsignalRatio;
  // there is no sonse to check how signal ratio is far from perfect 
  // matching line if the if below it is true
  if (ns < sigClusterN) {
      sigClusterN=ns;
      if (debug) cout<<"n1 < p1/PNsignalRatio";
      if (debug) cout<<"Attempting to del cluster N "<<clusterIdx<<" ... \n";
      pkg->DelClusterOI(clusterIdx,SIDEN);
  } else {
      //Let's see how signal ratio is far from perfect matching line
      Chicomb = DistToPML(ps,ns);
      if (debug) cout<<"Chic "<<Chicomb<<"\n";
      if (Chicomb > falpha2) {
         //it is near, so we can risk throwing this cluster away too
	 if (debug) cout<<"Attempting to del cluster N "<<clusterIdx<<"...\n"; 
         pkg->DelClusterOI(clusterIdx,SIDEN);
      } else {
         clusterN->CutTotalSignal(sigClusterN);
	 if (debug) cout <<"Signal cut  |||||||||||||\n";
      }
  }
  sigErrorClusterN= clusterN->GetTotalSignalError(); 
  CreateNewRecPoint(posClusterP,posErrorClusterP,posClusterN,posErrorClusterN,
                    sigClusterP+sigClusterN, sigErrorClusterN+sigErrorClusterP,
      	            clusterP, clusterN, 0.75);

}


//------------------------------------------------------
void AliITSClusterFinderSSD::
ResolveNClusterWithOneCross(AliITSpackageSSD *pkg, Int_t clusterIndex)
{

  AliITSclusterSSD * clusterP;
  AliITSclusterSSD * clusterN;

  Float_t posClusterP;              //Cluster P position in strip coordinates
  Float_t posClusterN;              //Cluster N position in strip coordinates

  Float_t posErrorClusterP;
  Float_t posErrorClusterN;

  Float_t sigClusterP;
  Float_t sigClusterN;

  Float_t sigErrorClusterP;
  Float_t sigErrorClusterN;

  Float_t ps;
  Float_t ns;
  
  Float_t Chicomb;
  Int_t clusterIdx;
  
  if (debug) cout<<"ResolveNClusterWithOneCross\n"; 
  
  clusterN=pkg->GetNSideCluster(clusterIndex);
  posClusterN = GetClusterZ(clusterN);
  posErrorClusterN = clusterN->GetPositionError();
  sigClusterN = ns = clusterN->GetTotalSignal();
  sigErrorClusterN = clusterN->GetTotalSignalError(); 
  //carefully, it returns index in TClonesArray
  clusterIdx = clusterN->GetCross(0);
  clusterP=this->GetPSideCluster(clusterIdx);
  ps = clusterP->GetTotalSignal();
  posClusterP = GetClusterZ(clusterP);
  posErrorClusterP = clusterP->GetPositionError();
  pkg->DelCluster(clusterIndex,SIDEN);
  sigClusterP=ns*PNsignalRatio;
  // there is no sonse to check how signal ratio is far from perfect 
  // matching line if the if below it is true
  if (ps < sigClusterP) {
      sigClusterP = ps;
      if (debug) cout<<"ps < ns*PNsignalRatio";
      if (debug) cout<<"Attempting to del cluster P "<<clusterIdx<<" ... \n";
      pkg->DelClusterOI(clusterIdx,SIDEP);
  } else {
      //Let's see how signal ratio is far from perfect matching line
       Chicomb = DistToPML(ps,ns);
       if (debug) cout<<"Chic "<<Chicomb<<"\n";
       if (Chicomb > falpha2) {
          //it is near, so we can risk frowing this cluster away too
          if (debug) cout<<"Attempting to del cluster P "<<clusterIdx<<"...\n";
	  pkg->DelClusterOI(clusterIdx,SIDEP);
       } else {
          clusterN->CutTotalSignal(sigClusterP);
	  if (debug) cout <<"Signal cut  ------------\n";
       }
  }
  sigErrorClusterP= clusterP->GetTotalSignalError(); 
  CreateNewRecPoint( posClusterP,posErrorClusterP,posClusterN,posErrorClusterN,
                     sigClusterP+sigClusterN,sigErrorClusterN+sigErrorClusterP,
      	             clusterP, clusterN, 0.75);


}



//-------------------------------------------------------

Bool_t AliITSClusterFinderSSD::
ResolvePackageBestCombin(AliITSpackageSSD *pkg)
{
  if (debug) cout<<"NumNcl==NumPcl ("<<pkg->GetNumOfClustersN()
      <<"=="<<pkg->GetNumOfClustersP()<<"); Generating combinations ... \n";


  AliITSclusterSSD * clusterP;
  AliITSclusterSSD * clusterN;

  Int_t Ncombin; //Number of combinations
  Int_t itera;   //iterator
  Int_t sizet=1;   //size of array to allocate 
 
  Int_t NP = pkg->GetNumOfClustersP();
  for(itera =2; itera <= NP ;itera ++) {
     sizet=sizet*itera;
     if (sizet > 10000) {
       sizet=10000;
       break;
     }
  }

  Int_t** combin = new (Int_t*)[sizet]; //2D array to keep combinations in

  for(itera =0; itera <sizet;itera++) {
   combin[itera] = new Int_t[NP+1];
  }
     
  pkg->GetAllCombinations(combin,Ncombin,sizet);
 
  if (Ncombin==0) { 
    if (debug) cout<<"No combin Error";
    return kFALSE;
  }
//calculate which combination fits the best to perfect matching line 
  Int_t bc = GetBestComb(combin,Ncombin,NP,pkg); 
  if (debug) cout<<"  bc "<<bc <<"\n";
              
  for(itera =0; itera < NP; itera ++) {
      clusterP = pkg->GetPSideCluster(itera);
      //carefully here
      //becase AliITSclusterSSD::GetCross returns index in 
      //AliITSmoduleSSD.fClusterP, which is put to "combin"        
      clusterN = GetNSideCluster(combin[bc][itera]);
      CreateNewRecPoint(clusterP, clusterN, 0.75);
  } 

  for(itera =0; itera <sizet;itera++) {
     delete [](combin[itera]);
  }
  delete [] combin;  
  return kTRUE;

}


//----------------------------------------------------
void AliITSClusterFinderSSD::
ResolveOneBestMatchingPoint(AliITSpackageSSD *pkg)
{

 Int_t ni, pi;

 Int_t prvP, prvN;
 Int_t nextP, nextN=0;

  
 AliITSclusterSSD * clusterP;
 AliITSclusterSSD * clusterN;

 Bool_t split = kFALSE;

 if (debug) cout<<"ResolveOneBestMatchingPoint \n";

 GetBestMatchingPoint(pi, ni, pkg);
         
 clusterP = GetPSideCluster(pi);
 clusterN = GetNSideCluster(ni);
              
 CreateNewRecPoint(clusterP, clusterN, 0.75);
 
 if ((nextP=pkg->GetNextPIdx(pi))!=-1)
   if ((nextN=pkg->GetNextNIdx(ni))!=-1)
     if ((prvP=pkg->GetPrvPIdx(pi))!=-1)
       if ((prvN=pkg->GetPrvNIdx(ni))!=-1)
          if( !(GetPSideCluster(prvP)->IsCrossingWith(nextN)))
            if( !(GetPSideCluster(nextP)->IsCrossingWith(prvN)))
              {
                 split = kTRUE;
              }	  
 
 
 pkg->DelClusterOI(pi, SIDEP);
 pkg->DelClusterOI(ni, SIDEN);
 
 if (split) {
   if (debug) cout<<"spltting package ...\n";
   new ((*fPackages)[fNPackages]) AliITSpackageSSD(fClusterP,fClusterN);
   pkg->
     SplitPackage(nextP,nextN,(AliITSpackageSSD*)((*fPackages)[fNPackages++]));
 }

}


//--------------------------------------------------
void  AliITSClusterFinderSSD::ResolveSimplePackage(AliITSpackageSSD *pkg)
{

  AliITSclusterSSD * clusterP;
  AliITSclusterSSD * clusterN;

  clusterP         = pkg->GetPSideCluster(0);
        
  clusterN         = pkg->GetNSideCluster(0);

  CreateNewRecPoint(clusterP, clusterN, 1.0);


} 


//--------------------------------------------------
void AliITSClusterFinderSSD:: ResolvePackageWithOnePSideCluster(AliITSpackageSSD *pkg) {


/*
 \   \   \  /
  \   \   \/
   \   \  /\
    \   \/  \
     \  /\   \
      \/  \   \
      /\   \   \
     /  \   \   \  

*/


/************************/
/**************************/
/*XP SP itg jest tylko jeden nie musi byc tablica i w peetli nie trzeba po niej iterowcac*/
/***************************/


 Int_t k;
 AliITSclusterSSD * clusterP;
 AliITSclusterSSD * clusterN;
 
 Int_t NN=pkg->GetNumOfClustersN();
 Float_t sumsig = 0;
 
 Float_t   XP;
 Float_t   XPerr;
 
 Float_t * XN = new Float_t[NN];
 Float_t * XNerr = new Float_t[NN];
 
 Float_t * SP = new Float_t[NN];
 Float_t * SN = new Float_t[NN];
 
 Float_t * SPerr = new Float_t[NN];
 Float_t * SNerr = new Float_t[NN];
 
 Float_t p1;
 Float_t p1err;
 
 clusterP = pkg->GetPSideCluster(0);
 p1       = clusterP->GetTotalSignal();
   
 XP    = GetClusterZ(clusterP);
 XPerr = clusterP->GetPositionError();
 p1err = clusterP->GetTotalSignalError();
   
 for(k=0;k<NN;k++) {
    SN[k]    = pkg->GetNSideCluster(k)->GetTotalSignal();
    SNerr[k] = pkg->GetNSideCluster(k)->GetTotalSignalError();
    sumsig   += SN[k];
 }
 for(k=0;k<NN;k++) {
    clusterN = pkg->GetNSideCluster(k);
    SP[k]= p1*SN[k]/sumsig;
    SPerr[k] = p1err*SN[k]/sumsig;
    XN[k]=GetClusterZ(clusterN);
    XNerr[k]= clusterN->GetPositionError();
    CreateNewRecPoint(XP,XPerr, XN[k],XNerr[k], SP[k]+SN[k], SPerr[k]+SNerr[k], clusterP, clusterN, 1.0);

 }

}


//---------------------------------------------------------
void AliITSClusterFinderSSD::ResolvePackageWithOneNSideCluster(AliITSpackageSSD *pkg) {


/*
    \    /   /   /
     \  /   /   /
      \/   /   /
      /\  /   /
     /  \/   /
    /   /\  /  
   /   /  \/
  /   /   /\
 /   /   /  \  
*/


 Int_t k;
 AliITSclusterSSD * clusterP;
 AliITSclusterSSD * clusterN;
 
 Int_t NP=pkg->GetNumOfClustersP();
 Float_t sumsig = 0;
 
 Float_t * XP = new Float_t[NP];
 Float_t * XPerr = new Float_t[NP];
 
 Float_t   XN;
 Float_t   XNerr;
 
 Float_t * SP = new Float_t[NP];
 Float_t * SN = new Float_t[NP];
 
 Float_t * SPerr = new Float_t[NP];
 Float_t * SNerr = new Float_t[NP];
 
 Float_t n1;
 Float_t n1err;

 clusterN = pkg->GetNSideCluster(0);
 n1       = clusterN->GetTotalSignal();
   
 XN = GetClusterZ(clusterN);
 XNerr = clusterN->GetPositionError();
 
 n1err=clusterN->GetTotalSignalError();
   
 for(k=0;k<NP;k++) {
    SP[k] = pkg->GetPSideCluster(k)->GetTotalSignal();
    sumsig += SP[k];
    SPerr[k] = pkg->GetPSideCluster(k)->GetTotalSignalError();
 }

 for(k=0;k<NP;k++) {
    clusterP = pkg->GetPSideCluster(k);
    SN[k]= n1*SP[k]/sumsig;
    XP[k]=GetClusterZ(clusterP);
    XPerr[k]= clusterP->GetPositionError();
    SNerr[k] = n1err*SP[k]/sumsig;
    CreateNewRecPoint(XP[k],XPerr[k], XN,XNerr, SP[k]+SN[k], SPerr[k]+SNerr[k],clusterP, clusterN, 1.0);
  }
  
               
}

/**********************************************************************/
/*********      2  X   2  *********************************************/
/**********************************************************************/

void AliITSClusterFinderSSD::
ResolveTwoForTwoPackage(AliITSpackageSSD *pkg)
{

  AliITSclusterSSD *clusterP1 = pkg->GetPSideCluster(0);
  AliITSclusterSSD *clusterP2 = pkg->GetPSideCluster(1);
  AliITSclusterSSD *clusterN1 = pkg->GetNSideCluster(0);
  AliITSclusterSSD *clusterN2 = pkg->GetNSideCluster(1);

  AliITSclusterSSD *tmp;

  Float_t p1sig;
  Float_t p2sig;
  Float_t n1sig;
  Float_t n2sig;

  Float_t p1sigErr;
  Float_t p2sigErr;
  Float_t n1sigErr;
  Float_t n2sigErr;

  Float_t ZposP1;
  Float_t ZposP2;
  Float_t ZposN1;
  Float_t ZposN2;

  Float_t ZposErrP1;
  Float_t ZposErrP2;
  Float_t ZposErrN1;
  Float_t ZposErrN2;

  Float_t D12;

  Float_t Chicomb1;
  Float_t Chicomb2;

  if(clusterP1->GetDigitStripNo(0) > clusterP2->GetDigitStripNo(0)) {
     if (debug) cout<<"P strips flip\n";
     tmp       = clusterP1;
     clusterP1 = clusterP2;
     clusterP2 = tmp;
  }
  if(clusterN1->GetDigitStripNo(0) > clusterN2->GetDigitStripNo(0)) {
     if (debug) cout<<"N strips flip\n";
     tmp       = clusterN1;
     clusterN1 = clusterN2;
     clusterN2 = tmp;
  }
  p1sig = clusterP1->GetTotalSignal();
  p2sig = clusterP2->GetTotalSignal();
  n1sig = clusterN1->GetTotalSignal();
  n2sig = clusterN2->GetTotalSignal();


  p1sigErr = clusterP1->GetTotalSignalError();
  n1sigErr = clusterN1->GetTotalSignalError();
  p2sigErr = clusterP2->GetTotalSignalError();
  n2sigErr = clusterN2->GetTotalSignalError();
          
  ZposP1 = clusterP1->GetPosition();
  ZposN1 = clusterN1->GetPosition();
  ZposP2 = clusterP2->GetPosition();
  ZposN2 = clusterN2->GetPosition();

  ZposErrP1 = clusterP1->GetPositionError();
  ZposErrN1 = clusterN1->GetPositionError();
  ZposErrP2 = clusterP2->GetPositionError();
  ZposErrN2 = clusterN2->GetPositionError();

  //in this may be two types:
  //    1.When each cluster crosses with 2 clusters on the other side
  //            gives 2 ghosts and 2 real points
  //
  //    2.When two clusters (one an each side) crosses with only one on 
  //           the other side and two crosses (one on the each side) with 
  //           two on the other gives 2 or 3 points,
   
  if (debug) cout<<"2 for 2 ambiguity ...";
  /***************************/
  /*First sort of combination*/
  /***************************/
  
  if((clusterP1->GetCrossNo()==2) && (clusterN1->GetCrossNo()==2)) {  
     if (debug) cout<<"All clusters has 2 crosses\n";
     D12 = TMath::Sqrt((Float_t)((p1sig -p2sig)*(p1sig -p2sig) + (n1sig -n2sig)*(n1sig -n2sig)));
	    
     if(debug) cout<<"Distance between points in P(N) plane D12 = "<<D12<<"\n";
            
     /*********************************************/
     /*resolving ambiguities:                     */
     /*we count Chicomb's and we take combination */
     /*giving lower Chicomb as a real points      */
     /*Keep only better combinantion              */
     /*********************************************/
	    
     if (D12 > (falpha3*17.8768)) {
       if (debug) cout<<"decided to take only one pair \n";
       Chicomb1 = DistToPML(p1sig,n1sig) + DistToPML(p2sig,n2sig);
       Chicomb2 = DistToPML(p2sig,n1sig) + DistToPML(p1sig,n2sig);
       if (debug) {
	 cout<<" 00 11 combination : "<<Chicomb1<<" 01 10 combination : "<<Chicomb2<<" \n";
	 cout<<"p1 = "<<p1sig<<"  n1 = "<<n1sig<<"  p2 = "<<p2sig<<"  n2 = "<<n2sig<<"\n"; 
       }
       if (Chicomb1 < Chicomb2) {
	  if (debug) cout<<"00 11";
	  CreateNewRecPoint(ZposP1,ZposErrP1, ZposN1,ZposErrN1, p1sig+n1sig, p1sigErr+n1sigErr, clusterP1, clusterN1, 0.75);   
	  CreateNewRecPoint(ZposP2,ZposErrP2, ZposN2,ZposErrN2, p2sig+n2sig, p2sigErr+n2sigErr, clusterP2, clusterN2, 0.75); 

  
	  //second cominantion
       } else {
	  if (debug) cout<<"01 10";
	  CreateNewRecPoint(ZposP1,0, ZposN2,0, p1sig+n2sig, p1sigErr+n2sigErr, clusterP1, clusterN2, 0.75);   
	  CreateNewRecPoint(ZposP2,0, ZposN1,0, p2sig+n1sig, p2sigErr+n1sigErr, clusterP2, clusterN1, 0.75);    
       } //end second combinantion
       //if (D12 > falpha3*17.8768)
       //keep all combinations
     } else {
        if (debug) cout<<"We decide to take all points\n";
	CreateNewRecPoint(ZposP1,ZposErrP1, ZposN1,ZposErrN1, p1sig+n1sig, p1sigErr+n1sigErr, clusterP1, clusterN1, 0.5);   
	CreateNewRecPoint(ZposP2,ZposErrP2, ZposN2,ZposErrN2, p2sig+n2sig, p2sigErr+n2sigErr, clusterP2, clusterN2, 0.5);   
	CreateNewRecPoint(ZposP1,ZposErrP1, ZposN2,ZposErrN2, p1sig+n2sig, p1sigErr+n2sigErr, clusterP1, clusterN2, 0.5);   
	CreateNewRecPoint(ZposP2,ZposErrP2, ZposN1,ZposErrN1, p2sig+n1sig, p2sigErr+n1sigErr, clusterP2, clusterN1, 0.5);    
     }
  }
  // ad.2 Second type of combination
  else {
    Chicomb1 = DistToPML(p1sig,n1sig) + DistToPML(p2sig,n2sig);
    if (debug) cout<<"\nhere can be reconstructed 3 points: chicomb = "<<Chicomb1<<"\n"; 
  
    if (Chicomb1<falpha1) {
       if (debug) cout<<"\nWe decided to take 3rd point"; 
       if (clusterP1->GetCrossNo()==1) {
	  if (debug) cout<<"...  P1 has one cross\n"; 
	  n1sig = p1sig/PNsignalRatio;
	  p2sig = n2sig*PNsignalRatio;
              
	  clusterN1->CutTotalSignal(n1sig);
	  clusterP2->CutTotalSignal(p2sig);
       
	  CreateNewRecPoint(ZposP2,ZposErrP2, ZposN1,ZposErrN1, p2sig+n1sig, p2sigErr+n1sigErr, clusterP2, clusterN1, 0.5);    

	  n1sig = clusterN1->GetTotalSignal();
	  p2sig = clusterP2->GetTotalSignal();

       } else {
 	  if (debug) cout<<"...  N1 has one cross\n";
	      
	  n2sig=p2sig/PNsignalRatio;
	  p1sig=n1sig*PNsignalRatio;

	  clusterN2->CutTotalSignal(n2sig);
	  clusterP1->CutTotalSignal(p1sig);

	  CreateNewRecPoint(ZposP1,ZposErrP1, ZposN2,ZposErrN2, p1sig+n2sig, p1sigErr+n2sigErr, clusterP1, clusterN2, 0.5);   
	  
	  n2sig=clusterN2->GetTotalSignal();
	  p1sig=clusterP1->GetTotalSignal();
       }
    } else {
       if (debug) cout<<"\nWe decided  NOT to take 3rd point\n"; 
    }

    CreateNewRecPoint(ZposP1,ZposErrP1, ZposN1,ZposErrN1, p1sig+n1sig, p1sigErr+n1sigErr, clusterP1, clusterN1, 1.0);   
    CreateNewRecPoint(ZposP2,ZposErrP2, ZposN2,ZposErrN2, p2sig+n2sig, p2sigErr+n2sigErr, clusterP2, clusterN2, 1.0);   
  
  }

}



//------------------------------------------------------
Bool_t AliITSClusterFinderSSD::    
CreateNewRecPoint(Float_t P, Float_t dP, Float_t N, Float_t dN,
                  Float_t Sig,Float_t dSig, 
                  AliITSclusterSSD *clusterP, AliITSclusterSSD *clusterN,
                  Stat_t prob)
{

  const Float_t kdEdXtoQ = 2.778e+8;
  const Float_t kconv = 1.0e-4; 
  const Float_t kRMSx = 20.0*kconv; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 830.0*kconv; // microns->cm ITS TDR Table 1.3


  Float_t p=P;
  Float_t n=N;
  if (GetCrossing(P,N)) {
     GetCrossingError(dP,dN);
     AliITSRawClusterSSD cnew;
     Int_t nstripsP=clusterP->GetNumOfDigits();
     Int_t nstripsN=clusterN->GetNumOfDigits();
     cnew.fMultiplicity=nstripsP;
     cnew.fMultiplicityN=nstripsN;
     /*
     if (nstripsP > 100) {
        printf("multiplicity > 100 - increase the dimension of the arrays %d\n",nstripsP);
        nstripsP=100;
     }
     Int_t i;
     for(i=0;i<nstripsP;i++) {
       // check if 'clusterP->GetDigitStripNo(i)' returns the digit index
       cnew.fIndexMap[i] = clusterP->GetDigitStripNo(i); 
     } 
     if (nstripsN > 100) {
        printf("multiplicity > 100 - increase the dimension of the arrays %d\n",nstripsN);
        nstripsN=100;
     }
     for(i=0;i<nstripsN;i++) {
       // check if 'clusterN->GetDigitStripNo(i)' returns the digit index
       cnew.fIndexMapN[i] = clusterN->GetDigitStripNo(i); 
     }
     */
     cnew.fQErr=dSig;
     //cnew.fProbability=(float)prob; 
     fITS->AddCluster(2,&cnew);
     // add the rec point info
     AliITSRecPoint rnew;
     rnew.SetX(P*kconv);
     rnew.SetZ(N*kconv);
     rnew.SetQ(Sig);
     rnew.SetdEdX(Sig/kdEdXtoQ);
     rnew.SetSigmaX2( kRMSx* kRMSx); 
     rnew.SetSigmaZ2( kRMSz* kRMSz);
     rnew.SetProbability((float)prob);
     fITS->AddRecPoint(rnew);
     /*
     // it was
     fPointsM->AddLast( (TObject*) 
                   new AliITSRecPointSSD(P,dP,0,0,N,dN,Sig,dSig,fLayer,fLad,fDet,clusterP,clusterN,prob) );
     */
    
     if(debug) cout<<"\n"<<": ImpactPoint Created: X = "<<P<<" Z = "<<N<<";   P = "<<p<<"  N = "<<n<<"\n";
     return kTRUE;
  } else { 
     if (debug) {
       cout<<"\n"<<": ATENTION : stupid ImpactPoit Point X = "<<P<<" Z = "<<N<<";  P = "<<p<<"  N = "<<n<<"\n";
     }
     return kFALSE;
  }
}



/**********************************************************************/
Bool_t  AliITSClusterFinderSSD::
CreateNewRecPoint(AliITSclusterSSD *clusterP, AliITSclusterSSD *clusterN, Stat_t prob)
{
  Float_t posClusterP;  //Cluster P position in strip coordinates
  Float_t posClusterN;  //Cluster N position in strip coordinates
 
  Float_t posErrorClusterP;
  Float_t posErrorClusterN;

  Float_t sigClusterP;
  Float_t sigClusterN;

  Float_t sigErrorClusterP;
  Float_t sigErrorClusterN;

  posClusterP      = clusterP->GetPosition();
  posErrorClusterP = clusterP->GetPositionError();
  sigClusterP      = clusterP->GetTotalSignal();
  sigErrorClusterP = clusterP->GetTotalSignalError();

  posClusterN      = clusterN->GetPosition();
  posErrorClusterN = clusterN->GetPositionError();
  sigClusterN      = clusterN->GetTotalSignal();
  sigErrorClusterN = clusterN->GetTotalSignalError();

  return CreateNewRecPoint( posClusterP,posErrorClusterP,  posClusterN,posErrorClusterN, 
                     sigClusterP+sigClusterN, sigErrorClusterN+sigErrorClusterP,
      	             clusterP, clusterN, prob);

}

//--------------------------------------------------
Bool_t AliITSClusterFinderSSD::IsCrossing(AliITSclusterSSD* p, AliITSclusterSSD* n)
{
  Float_t x = p->GetPosition();
  Float_t y = n->GetPosition();
  return GetCrossing(x,y);
}


//----------------------------------------------
Bool_t AliITSClusterFinderSSD::Strip2Local( Float_t stripP, Float_t stripN, Float_t &Z,Float_t &X)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

/*
 Z = (stripN-stripP)*fFactorOne;  
 X = (stripN + stripP - fStripZeroOffset)*fFactorTwo;
*/

 Float_t P =stripP;
 Float_t N =stripN;

 GetCrossing(P,N);
 X=P;
 Z=N;
 if (debug) cout<<"P="<<stripP<<" N="<<stripN<<"   X = "<<X<<" Z = "<<Z<<"\n";
 if ((Z<2.1)&&(Z>-2.1)&&(X<3.65)&&(X>-3.65)) return true; 
    else return false;  

}


//-----------------------------------------------------------
Float_t  AliITSClusterFinderSSD::GetClusterZ(AliITSclusterSSD* clust)
{

  return clust->GetPosition();

}

//-------------------------------------------------------------
Int_t AliITSClusterFinderSSD::GetBestComb
(Int_t** comb,Int_t Ncomb, Int_t Ncl, AliITSpackageSSD * pkg)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

//returns index of best combination in "comb"
//comb : sets of combinations
//       in given combination on place "n" is index of 
//       Nside cluster coresponding to "n" P side cluster
//
//Ncomb : number of combinations == number of rows in "comb"
//
//Ncl   : number of clusters in each row  == number of columns in "comb"
//
//pkg   : package 


  Float_t currbval=-1;  //current best value, 
                        //starting value is set to number negative, 
                        //to be changed in first loop turn automatically
      
  Int_t  curridx=0;     //index of combination giving currently the best chi^2 
  Float_t chi;
  Float_t ps, ns;       //signal of P cluster and N cluster

  Int_t i,j;
  for(i=0;i<Ncomb;i++)
    {
      chi=0;
      
      for(j=0;j<Ncl;j++)
      {
        ps = pkg->GetPSideCluster(j)->GetTotalSignal();  //carrefully here, different functions
        ns = GetNSideCluster(comb[i][j])->GetTotalSignal(); 
        chi+=DistToPML(ps,ns);
      }
      
      if (currbval==-1) currbval = chi;
      if (chi<currbval)
        {
          curridx = i;
          currbval  =chi;
        }
    }


  return curridx;

}

//--------------------------------------------------
void AliITSClusterFinderSSD::GetBestMatchingPoint
(Int_t & ip, Int_t & in, AliITSpackageSSD* pkg )
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl
 
  if (debug) pkg->PrintClusters();
 
  Float_t ps, ns;  //signals on side p & n
  
  Int_t nc;       // number of crosses in given p side cluster
  Int_t n,p;
  AliITSclusterSSD *curPcl, *curNcl;  //current p side cluster
  Float_t  bestchi, chi;
  ip=p=pkg->GetPSideClusterIdx(0);
  in=n=pkg->GetNSideClusterIdx(0);
	
  bestchi=DistToPML(  pkg->GetPSideCluster(0)->GetTotalSignal(),
                  pkg->GetNSideCluster(0)->GetTotalSignal()  );
  Int_t i,j;
  for(i = 0; i< pkg->GetNumOfClustersP(); i++)
    {
      
      p =  pkg->GetPSideClusterIdx(i);   
      curPcl= GetPSideCluster(p);
      nc=curPcl->GetCrossNo();
      ps=curPcl->GetTotalSignal();
      
      for(j = 0; j< nc; j++)
        {
         n=curPcl->GetCross(j);
         curNcl= GetNSideCluster(n);
         ns=curNcl->GetTotalSignal();
         chi = DistToPML(ps, ns);
    
         if (chi <= bestchi)
          {
            bestchi = chi;
            in = n;
            ip = p;
          }
        }
    }
  
}


//------------------------------------------------------
void  AliITSClusterFinderSSD::CalcStepFactor(Float_t Psteo, Float_t Nsteo )
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl


  // 95 is the pitch, 4000 - dimension along z ?
  /*
  fSFF = ( (Int_t)  (Psteo*40000/95 ) );// +1;
  fSFB = ( (Int_t)  (Nsteo*40000/95 ) );// +1;
  */


  Float_t dz=fSegmentation->Dz();

  fSFF = ( (Int_t)  (Psteo*dz/fPitch ) );// +1;
  fSFB = ( (Int_t)  (Nsteo*dz/fPitch ) );// +1;

}


//-----------------------------------------------------------
AliITSclusterSSD* AliITSClusterFinderSSD::GetPSideCluster(Int_t idx)
{


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl



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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

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


//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl

  return (side) ? GetPSideCluster(idx) : GetNSideCluster(idx);
}


//--------------------------------------------------------
void AliITSClusterFinderSSD::ConsumeClusters()
{

 Int_t i;
 for(i=0;i<fNPackages;i++)
  {
    ((AliITSpackageSSD*)((*fPackages)[i]))->ConsumeClusters();
  }

}


//--------------------------------------------------------
void AliITSClusterFinderSSD::ReconstructNotConsumedClusters()
{
  Int_t i;
  AliITSclusterSSD *cluster;
  Float_t pos;
  Float_t sig;
  Float_t sigerr;
  Float_t x1,x2,z1,z2;

  Float_t Dx = fSegmentation->Dx();
  Float_t Dz = fSegmentation->Dz();
  
  const Float_t kdEdXtoQ = 2.778e+8; // GeV -> number of e-hole pairs
  const Float_t kconv = 1.0e-4; 
  const Float_t kRMSx = 20.0*kconv; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 830.0*kconv; // microns->cm ITS TDR Table 1.3
  
  for(i=0;i<fNClusterP;i++)
    {
      //printf("P side cluster: cluster number %d\n",i);
      cluster = GetPSideCluster(i);
      if (!cluster->IsConsumed())
        {
          if( (pos = cluster->GetPosition()) < fSFB)
            {
              sig = cluster->GetTotalSignal();
              
              sig += sig/PNsignalRatio;
              
              sigerr = cluster->GetTotalSignalError();
              x1 = -Dx/2 + pos *fPitch;
              z1 = -Dz/2;

              x2 = pos;
              z2 = 1;
              GetCrossing (x2,z2);
	      AliITSRawClusterSSD cnew;
	      Int_t nstripsP=cluster->GetNumOfDigits();
	      cnew.fMultiplicity=nstripsP;
	      cnew.fMultiplicityN=0;
	      /*
	      if (nstripsP > 100) {
		printf("multiplicity > 100 - increase the dimension of the arrays %d\n",nstripsP);
		nstripsP=100;
	      }
              Int_t k;
	      for(k=0;k<nstripsP;k++) {
		// check if 'clusterP->GetDigitStripNo(i)' returns 
		// the digit index and not smth else
		cnew.fIndexMap[k] = cluster->GetDigitStripNo(k); 
	      } 
	      */
	      cnew.fQErr=sigerr;
	      //cnew.fProbability=0.75; 
	      fITS->AddCluster(2,&cnew);
	      // add the rec point info
	      AliITSRecPoint rnew;
              rnew.SetX(kconv*(x1+x2)/2);
              rnew.SetZ(kconv*(z1+z2)/2);
              rnew.SetQ(sig);
              rnew.SetdEdX(sig/kdEdXtoQ);
	      rnew.SetSigmaX2( kRMSx* kRMSx);
	      rnew.SetSigmaZ2( kRMSz* kRMSz);
              rnew.SetProbability(0.75);
              fITS->AddRecPoint(rnew);
	      /*
              fPointsM->AddLast( (TObject*) 
                   new AliITSRecPointSSD((x1+x2)/2 ,0,0,0,(z1+z2)/2,0,sig,sigerr,fLayer,fLad,fDet,cluster, 0.75) );
	      */

              if(debug) cout<<"\n"<<": SINGLE SIDE ImpactPoint Created: X = "
                  <<(x1+x2)/2<<" Z = "<<(z1+z2)/2<<";   Pos = "<<pos;
            }
        }
    }

  for(i=0;i<fNClusterN;i++)
    {
      // printf("N side cluster: cluster number %d\n",i);
      cluster = GetNSideCluster(i);
      if (!cluster->IsConsumed())
        {
          if( (pos = cluster->GetPosition()) < fSFF)
            {
              sig = cluster->GetTotalSignal();
              
              sig += sig*PNsignalRatio;
              
              sigerr = cluster->GetTotalSignalError();
              
              x1 = -Dx/2 + pos *fPitch;
              z1 = Dz/2;

              x2 = 1;
              z2 = pos;
              
              GetCrossing (x2,z2);
	      AliITSRawClusterSSD cnew;
	      Int_t nstripsN=cluster->GetNumOfDigits();
	      cnew.fMultiplicity=0;
	      cnew.fMultiplicityN=nstripsN;
	      /*
	      if (nstripsN > 100) {
		printf("multiplicity > 100 - increase the dimension of the arrays %d\n",nstripsN);
		nstripsN=100;
	      }
              Int_t k;
	      for(k=0;k<nstripsN;k++) {
		// check if 'clusterP->GetDigitStripNo(i)' returns 
		// the digit index and not smth else
		cnew.fIndexMapN[k] = cluster->GetDigitStripNo(k); 
	      } 
	      */
	      cnew.fQErr=sigerr;
	      //cnew.fProbability=0.75; 
	      fITS->AddCluster(2,&cnew);
	      // add the rec point info
	      AliITSRecPoint rnew;
              rnew.SetX(kconv*(x1+x2)/2);
              rnew.SetZ(kconv*(z1+z2)/2);
              rnew.SetQ(sig);
              rnew.SetdEdX(sig/kdEdXtoQ);
	      rnew.SetSigmaX2( kRMSx* kRMSx);
	      rnew.SetSigmaZ2( kRMSz* kRMSz);
              rnew.SetProbability(0.75);
              fITS->AddRecPoint(rnew);
	      /*
              fPointsM->AddLast( (TObject*) 
                   new AliITSRecPointSSD((x1+x2)/2 ,0,0,0,(z1+z2)/2,0,sig,sigerr,fLayer,fLad,fDet,cluster, 0.75) );
	      */

              if(debug) cout<<"\n"<<": SINGLE SIDE ImpactPoint Created: X = "
                  <<(x1+x2)/2<<" Z = "<<(z1+z2)/2<<";   Pos = "<<pos;
 
            }
        }
    }


}

//_______________________________________________________________________

Bool_t AliITSClusterFinderSSD::GetCrossing (Float_t &P, Float_t &N) 
{ 

   Float_t Dx = fSegmentation->Dx();
   Float_t Dz = fSegmentation->Dz();

   //P==X
   //N==Z

    Float_t dx = 0.1;
    P *= fPitch;
    N *= fPitch; 
    
    //P = (kZ * fTan + N + P)/2.0;         // x coordinate
    //N = kZ - (P-N)/fTan;                 // z coordinate
    
    if(fTanP + fTanN == 0) return kFALSE;

    N = (N - P + fTanN * Dz) / (fTanP + fTanN);  // X coordinate
    P = P + fTanP * N;                           // Y coordinate

    P -= Dx/2;
    N -= Dz/2;

    //   N = -(N - P + kZ/2*(fTanP + fTanN))/(fTanP + fTanN);
    //  P = (fTanP*(N-kX/2-kZ*fTanN/2) + fTanN*(P-kX/2-kZ*fTanP/2) )/(fTanP + fTanN);
    
    if ( ( N < -(Dz/2+dx) ) || ( N > (Dz/2+dx) ) ) return kFALSE;
    if ( ( P < -(Dx/2+dx) ) || ( P > (Dx/2+dx) ) ) return kFALSE;

    return kTRUE;   
}


//_________________________________________________________________________

void AliITSClusterFinderSSD::GetCrossingError(Float_t& dP, Float_t& dN)
{
  Float_t dz, dx;
  
  dz = TMath::Abs(( dP + dN )*fPitch/(fTanP + fTanN) );
  dx = fPitch*(TMath::Abs(dP*(1 - fTanP/(fTanP + fTanN))) +
               TMath::Abs(dN *fTanP/(fTanP + fTanN) ));
  
  dN = dz;
  dP = dx;
}


