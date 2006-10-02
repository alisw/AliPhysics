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
/* $Id$ */

#include <Riostream.h>
#include <TClonesArray.h>
#include "AliITSpackageSSD.h"
#include "AliITSclusterSSD.h"

const Bool_t AliITSpackageSSD::fgkSIDEP=kTRUE;
const Bool_t AliITSpackageSSD::fgkSIDEN=kFALSE;

ClassImp(AliITSpackageSSD)
////////////////////////////////////////////////////////////////////////////
//Class describing set of AliITSoneSideClusterSSDs, which contact each other.
//Piotr Krzysztof Skowronski
//Warsaw University of Technology
//skowron@if.pw.edu.pl
//
//--------------------------------------------------------------------------
AliITSpackageSSD::AliITSpackageSSD():
fClustersN(0),
fClustersP(0),
fNclustersN(0),
fNclustersP(0),
fClusterNIndexes(0),
fClusterPIndexes(0){
  // constructor
}


/*******************************************************/

AliITSpackageSSD::AliITSpackageSSD
(TClonesArray *clustersP, TClonesArray *clustersN):
fClustersN(clustersN),
fClustersP(clustersP),
fNclustersN(0),
fNclustersP(0),
fClusterNIndexes(0),
fClusterPIndexes(0){
  // constructor
  fClustersP=clustersP;
  fClustersN=clustersN;
  	
  fNclustersN=0;
  fClusterNIndexes = new TArrayI(300); 
	
  fNclustersP=0;
  fClusterPIndexes = new TArrayI(300);		
}

/*******************************************************/


AliITSpackageSSD::AliITSpackageSSD(Int_t len, TClonesArray *clustersP, TClonesArray *clustersN):
fClustersN(clustersN),
fClustersP(clustersP),
fNclustersN(0),
fNclustersP(0),
fClusterNIndexes(0),
fClusterPIndexes(0){	
  // constructor
  fClustersP=clustersP;
  fClustersN=clustersN;

  fNclustersN=0;
  fClusterNIndexes = new TArrayI(len); 
	
  fNclustersP=0;
  fClusterPIndexes = new TArrayI(len);		
}


/*******************************************************/

AliITSpackageSSD::~AliITSpackageSSD()
{
  // destructor
  delete fClusterNIndexes;
  delete fClusterPIndexes;		
}

/*******************************************************/

AliITSpackageSSD::AliITSpackageSSD(const AliITSpackageSSD &package) : 
    TObject(package),
fClustersN(package.fClustersN),
fClustersP(package.fClustersP),
fNclustersN(package.fNclustersN),
fNclustersP(package.fNclustersP),
fClusterNIndexes(0),
fClusterPIndexes(0){
  // copy constractor
  Int_t i;  //iterator
 
  for ( i =0; i<fNclustersN;i++)
    {
      fClusterNIndexes[i]= package.fClusterNIndexes[i]; 
    }
  
  for ( i =0; i<fNclustersP;i++)
    {
      fClusterPIndexes[i]= package.fClusterPIndexes[i]; 
    }
  
  
}
/*******************************************************/

AliITSpackageSSD&  
AliITSpackageSSD::operator=( const AliITSpackageSSD & package)
{
  //  assignment operator
  //
  Int_t i;  //iterator
  
  if (this == &package) return *this;
  fClustersN = package.fClustersN;
  fClustersP = package.fClustersP;
  
  fNclustersN= package.fNclustersN;
  fNclustersP= package.fNclustersP;
  
  for ( i =0; i<fNclustersN;i++)
    {
      fClusterNIndexes[i]= package.fClusterNIndexes[i]; 
    }
  
  for ( i =0; i<fNclustersP;i++)
    {
      fClusterPIndexes[i]= package.fClusterPIndexes[i]; 
    }
  
  if (fgkDebug) cout << "Copying function was used\n<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>";
  
  return *this; 
  
}

/*******************************************************/

Int_t AliITSpackageSSD::GetNSideClusterIdx(Int_t index) const
{
  // get N-side cluster
  // 
  if ((index>-1)&&(index<fNclustersN))
    return (*fClusterNIndexes)[index];
  else 
    {
      cout << "AliITSpackageSSD::GetNSideClusterIdx  : Out of Range\n";
      return -1;
    }
}
/*******************************************************/


Int_t AliITSpackageSSD::GetPSideClusterIdx(Int_t index) const
{
  // get P-side cluster
  //
  if ((index>-1)&&(index<fNclustersP))
    return (*fClusterPIndexes)[index];
  else 
    {
      cout << "AliITSpackageSSD::GetPSideClusterIdx  : Out of Range\n";
      return -1;
    } 
}
/*******************************************************/
AliITSclusterSSD*  
AliITSpackageSSD::GetPSideCluster(Int_t index)
{
  // get Pside cluster from the TClonesArray of SSD clusters
  //
  return (AliITSclusterSSD*)((*fClustersP)[GetPSideClusterIdx(index)]);
}

/*******************************************************/

AliITSclusterSSD*  
AliITSpackageSSD::GetNSideCluster(Int_t index)
{
  // get Nside cluster from the TClonesArray of SSD clusters
  //
  return (AliITSclusterSSD*)((*fClustersN)[GetNSideClusterIdx(index)]);
}


/*******************************************************/

Bool_t AliITSpackageSSD::GetClusterWithOneCross
(Int_t & index, Bool_t& side)
{
  // select clusters with on cross 
  //
  if((fNclustersP==0)||(fNclustersN==0) )
    {
      printf("Empty package ((fNclustersP==0)||(fNclustersN==0))\n");
      index = -2;
      return kFALSE;
    } 
  Int_t ind;
  
  ind =(*fClusterPIndexes)[fNclustersP-1]; 
  if (  ( ((AliITSclusterSSD*)(*fClustersP)[ind]  )->GetCrossNo()  ) ==1 )
    {
      //index=ind;
      index =fNclustersP-1; 
      side=fgkSIDEP;
      return kTRUE;
    }
  
  ind =(*fClusterNIndexes)[fNclustersN-1]; 
  if (  (  ((AliITSclusterSSD*)(*fClustersN)[ind]  )->GetCrossNo() ) ==1  )
    {
      //index=ind;
      index = fNclustersN-1;
      side=fgkSIDEN;
      return kTRUE;
    }
  
  
  ind =(*fClusterPIndexes)[0];
  if (  ( ((AliITSclusterSSD*)(*fClustersP)[ind]  )->GetCrossNo() ) ==1 )
    {
      //index=ind;
      index = 0;
      side=fgkSIDEP;
      return kTRUE;
    }
  
  
  ind =(*fClusterNIndexes)[0];
  if (  ( ((AliITSclusterSSD*)(*fClustersN)[ind]  )->GetCrossNo()  ) ==1  )
    {
      //    index=ind;
      index = 0;  
      side=fgkSIDEN;
      return kTRUE;
    }
  
  
  //Add for to be shure 
  index = -1;
  return kFALSE;
  
}
/*******************************************************/

void AliITSpackageSSD::DelCluster(Int_t index, Bool_t side)
{
  // call DelPCluster or DelNCluster depending on side
  //
  if(side==fgkSIDEP) DelPCluster(index); else DelNCluster(index);
}
/*******************************************************/
void AliITSpackageSSD::DelPCluster(Int_t index)
{
  //it not deletes delete given cluster physically, 
  //but only complytely erase it from package
  //all clusters are deleted automatically when TClonesArray is deleted
  
  Int_t i;
  Int_t idx;
  Int_t clToDelIdx = GetPSideClusterIdx(index); //Index of cluster in TClonesArray 
  AliITSclusterSSD *clToDel = GetPSideCluster(index); //cluster to be deleted
  Int_t ncr = clToDel->GetCrossNo();
  
  for (i =0;i<ncr;i++)
    {
      idx = clToDel->GetCross(i);
      ((AliITSclusterSSD *)((*fClustersN)[idx])   )->DelCross(clToDelIdx);
    }
  
  
  for (i=index;i<fNclustersP-1;i++)
    {
      (*fClusterPIndexes)[i]=(*fClusterPIndexes)[i+1];
    }
  fNclustersP--; 
  if (fgkDebug) cout<<"Cluster P ("<<index<<") deleted\n";
  
  
  for (i=0;i<fNclustersN;i++)
    {
      if ( (GetNSideCluster(i)->GetCrossNo())==0) DelNCluster(i);
    }
}



/*******************************************************/
void AliITSpackageSSD::DelNCluster(Int_t index)
{
  //it not deletes delete given cluster physically, 
  //but only complytely erase it from package
  //all clusters are deleted automatically when TClonesArray is deleted
  
  Int_t i;
  Int_t idx;
  Int_t clToDelIdx = GetNSideClusterIdx(index); //Index of cluster in TClonesArray 
  AliITSclusterSSD *clToDel = GetNSideCluster(index); //cluster to be deleted
  Int_t ncr = clToDel->GetCrossNo();
  
  for (i =0;i<ncr;i++)
    {
      idx = clToDel->GetCross(i);
      ((AliITSclusterSSD *)((*fClustersP)[idx])   )->DelCross(clToDelIdx);
    }
  
  
  for (i=index;i<fNclustersN-1;i++)
    {
      (*fClusterNIndexes)[i]=(*fClusterNIndexes)[i+1];
    }
  fNclustersN--; 
  if (fgkDebug) cout<<"Cluster N ("<<index<<") deleted\n";
  
  for (i=0;i<fNclustersP;i++)
    {
      if ( (GetPSideCluster(i)->GetCrossNo())==0) DelPCluster(i);
    }
  
}


/*******************************************************/

void AliITSpackageSSD::DelPClusterOI(Int_t index)
{
  //This function looks like this, 
//because probably cut cluster is 
//on the beginning or on the end of package 
 Int_t i;
 if( ((*fClusterPIndexes)[0]) == index) 
  {
    DelPCluster(0);
    return;
  }
 else
  { 
   if( ((*fClusterPIndexes)[fNclustersP-1]) ==index)
    {
      DelPCluster(fNclustersP-1);
      return;
    }
   else
    {
     for (i=1;i<fNclustersP-1;i++)
       {
         if( ((*fClusterPIndexes)[i])==index)
	  {
            DelPCluster(i);
            return;
	  }
       }
    }
  }
 
 cout<<"AliITSpackageSSD - DelPClusterOI: index "<<index<<" not found\n";

}


/*******************************************************/

void AliITSpackageSSD::DelNClusterOI(Int_t index)
{
//This function looks like this, 
//because probably cluster to cut is 
//on the beginning or on the end of package 

 Int_t i;
 if( ((*fClusterNIndexes)[0])==index) 
  {
    DelNCluster(0);
    return;
  }
 else
  { 
   if( ((*fClusterNIndexes)[fNclustersN-1])==index)
    {
      DelNCluster(fNclustersN-1);
      return;
    }
   else
    {
     for (i=1;i<fNclustersN-1;i++)
       {
         if( ((*fClusterNIndexes)[i])==index)
         {
          DelNCluster(i);
          return;
         }
       }
    }
  }
 cout<<"AliITSpackageSSD - DelNClusterOI: index "<<index<<" not found\n";
}


/*******************************************************/


void AliITSpackageSSD::DelClusterOI(Int_t index, Bool_t side)
{
  // delete cluster
 if (side == fgkSIDEP)
  {    
    DelPClusterOI(index);
  } 
  else
  {
    DelNClusterOI(index);
  }

}


/**********************************************/


void  AliITSpackageSSD::GetAllCombinations(Int_t**array,Int_t &num,Int_t sizet)
{
  // get all combinations
  Int_t *takenNcl = new Int_t[fNclustersN];
  
  num=0;
  
  if (fgkDebug) PrintClusters();

  for (Int_t i=0;i<fNclustersP;i++)
   {
     takenNcl[i]=-1;
   }
  //see comment on the beginning of MakeCombin
  if (fgkDebug) cout<<"GetAllCombinations entered";
  MakeCombin (array,num,0,takenNcl,sizet);

  delete []takenNcl; 
}
/**********************************************/


void  AliITSpackageSSD::MakeCombin
   (Int_t**arr,Int_t& nu, Int_t np, Int_t *occup, Int_t sizet)

{
//ATTENTION: anybody watching this function
//AliITSclusterSSD::GetCrossNo() returns index of cluster in main array belonging to AliITSmodulesSSD
//however, we have pointer to that array (TClonesArray)
//we can not use 
//Get?SideCluster because it takes index from local look_up_table
 
 Int_t i,j;
 
 //this cluster
 AliITSclusterSSD *cl=GetPSideCluster(np);

 Int_t nc = cl->GetCrossNo();  //number of crosses for this cluster
 Int_t indcro;                 //index of given cluster on side N that 
                               // this cluster crosses with
   
 if (np == fNclustersP-1) {
   for (i=0;i<nc;i++) {
     indcro=cl->GetCross(i);
     if(IsFree(indcro,np,occup)) {
        occup[np]=indcro;
        for(j=0;j<fNclustersP;j++)  
        {
	   if (nu<sizet) arr[nu][j]=occup[j];
	   else {
	       continue;}
        }
      
        occup[np]=-1;
        if (nu<sizet-1) nu++;
     }
   }
  } else {
    for (i=0;i<nc;i++) {
       indcro=cl->GetCross(i);
       if(IsFree(indcro,np,occup)) {
      	  occup[np]=indcro;
	  if (nu<sizet) MakeCombin(arr,nu,(np+1),occup,sizet);
	  //else printf("MakeComb - exceeding array size!\n");
       }
    }
    occup[np]=-1;
  } 

}

/**********************************************/
Bool_t  AliITSpackageSSD::IsFree(Int_t idx, Int_t nn, const Int_t *lis) const
{
  // 
  for (Int_t i=0;i<nn;i++)
    {
      if (lis[i]==idx) return kFALSE;
    }
  return kTRUE;
}

/**********************************************/
void AliITSpackageSSD::PrintClusters()
{
  // print cluster info
Int_t i,j;
cout<<"SIDE P\n";
for (i=0;i<fNclustersP;i++)
 {
   cout<<i<<".  IO="<<GetPSideClusterIdx(i)<<" NC="<<GetPSideCluster(i)->GetCrossNo()<<"     C. IDXs : ";
   for (j=0;j<GetPSideCluster(i)->GetCrossNo();j++)
    {
      cout<<GetPSideCluster(i)->GetCross(j)<<" ";
    }
 //  if (GetPSideCluster(i)->GetSide()) cout<<" P\n";
 //  else cout<<"BAD SIDE ==N\n";
    cout<<"\n";
   
 }

cout <<"SIDE N\n";
for (i=0;i<fNclustersN;i++)
 {
   cout<<i<<".  IO="<<GetNSideClusterIdx(i)<<" NC="<<GetNSideCluster(i)->GetCrossNo()<<"     C. IDXs : ";
   for (j=0;j<GetNSideCluster(i)->GetCrossNo();j++)
    {
      cout<<GetNSideCluster(i)->GetCross(j)<<" ";
    }
 //  if (GetNSideCluster(i)->GetSide()) cout<<" N\n";
 //  else cout<<"BAD SIDE ==P\n";
    cout<<"\n";   
 }
 
}

/**********************************************/
void AliITSpackageSSD::ConsumeClusters()
{
  // consume cluster
  register Int_t i;
  
  for(i=0;i<fNclustersP;i++)
    {
      GetPSideCluster(i)->Consume();
    }
  
  for(i=0;i<fNclustersN;i++)
    {
      GetNSideCluster(i)->Consume();
    }
  
}

/**********************************************/

Int_t AliITSpackageSSD::GetNextPIdx(Int_t OI) const
{
 //Returns index of next P cluster OI in package; OI == Original Inedx (in TClonesArray)
 //if not egsist return -1;
 for (Int_t i =0; i<fNclustersP-1;i++)
  {
    if(GetPSideClusterIdx(i) == OI)
       return GetPSideClusterIdx(i+1);
  }
 return -1;
}

/**********************************************/
Int_t AliITSpackageSSD::GetPrvPIdx(Int_t OI) const
{
 //Returns index of previous P cluster  OI in package; OI == Original Inedx (in TClonesArray)
 //if not egsist return -1;
 
 for (Int_t i =1; i<fNclustersP;i++)
  {
    if(GetPSideClusterIdx(i) == OI)
       return GetPSideClusterIdx(i-1);
  }
  return -1;
}
/**********************************************/
Int_t AliITSpackageSSD::GetNextNIdx(Int_t OI) const
{
//Returns index of next N cluster OI in package; OI == Original Inedx (in TClonesArray)
 //if not egsist return -1;
 for (Int_t i =0; i<fNclustersN-1;i++)
  {
    if(GetNSideClusterIdx(i) == OI)
       return GetNSideClusterIdx(i+1);
  }
 return -1;

}
/**********************************************/
Int_t  AliITSpackageSSD::GetPrvNIdx(Int_t OI) const
{
 //Returns index of previous N cluster OI in package; OI == Original Inedx (in TClonesArray)
 //if not egsist return -1;
 
 for (Int_t i =1; i<fNclustersN;i++)
  {
    if(GetNSideClusterIdx(i) == OI)
       return GetNSideClusterIdx(i-1);
  }
  return -1;
 
}

void  AliITSpackageSSD::SplitPackage(Int_t pi, Int_t ni, AliITSpackageSSD* pkg)
{
  // split package of clusters
  Int_t p=-1, n=-1;
  Int_t i;
  for (i=0;i<fNclustersN;i++)
   {
     if((*fClusterNIndexes)[i]==ni) 
      {
        n = i;
	break;
      }
    } 
   
  for (i=0;i<fNclustersP;i++)
   {
     if((*fClusterPIndexes)[i]==pi) 
      {
        p = i;
	break;
      }
    }  
  if (fgkDebug) {
    cout<<" p = "<<p<<"  n = "<<n;
  }
  if ((p==-1)||(n==-1)) return;
  
  for (i=p;i<fNclustersP;i++)
   {
     pkg->AddPSideCluster(GetPSideClusterIdx(i));       
   }
  fNclustersP = p;
  
  for (i=n;i<fNclustersN;i++)
   {
     pkg->AddNSideCluster(GetNSideClusterIdx(i));       
   }
  fNclustersN = n;
 
  cout<<"  After split: fNclustersP = "<< fNclustersP<< "fNclustersN = "<< fNclustersN<<"\n";
}
