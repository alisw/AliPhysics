#include <iostream.h>

#include "AliITSdigit.h"
#include "AliITSclusterSSD.h"

ClassImp(AliITSclusterSSD)

AliITSclusterSSD::AliITSclusterSSD()
{
  // default constructor
	fSide        = kTRUE;
	fDigits      = 0;
	fNDigits     = 0;
	fDigitsIndex = 0;
        fNCrosses    = 0;
	fTotalSignal = -1;
	fNTrack      = -1;
	fLeftNeighbour  = kFALSE;
	fRightNeighbour = kFALSE;
	fCrossedClusterIndexes = new TArrayI(300);
	fConsumed=kFALSE;

}
/*************************************************************************/

AliITSclusterSSD::AliITSclusterSSD
  (Int_t ndigits, Int_t *DigitIndexes, 
   TObjArray *Digits, Bool_t side)
{
  // non-default constructor
  
	fNDigits = ndigits;
	fDigits = Digits;
	fSide = side;
	fDigitsIndex = new TArrayI(fNDigits,DigitIndexes );	
	fNCrosses    = 0;
	fCrossedClusterIndexes = new TArrayI(300);
	fLeftNeighbour  = kFALSE;
	fRightNeighbour = kFALSE;
	fTotalSignal =-1;
	fNTrack      = -1;
        fConsumed=kFALSE;
}
/*************************************************************************/
AliITSclusterSSD::~AliITSclusterSSD()
{
  // destructor
  delete fDigitsIndex;
  delete fCrossedClusterIndexes;
}
 
/*************************************************************************/
AliITSclusterSSD::AliITSclusterSSD(const AliITSclusterSSD &OneSCluster)
{
  // copy constructor
  if (this == &OneSCluster) return;
  fNDigits = OneSCluster.fNDigits;
  fSide=OneSCluster.fSide;
  fDigits=OneSCluster.fDigits;
  fDigitsIndex = new TArrayI(fNDigits);
  fLeftNeighbour  = OneSCluster.fLeftNeighbour;
  fRightNeighbour = OneSCluster.fRightNeighbour;
  fTotalSignal =-1;
  fNTrack      = -1;
  fNCrosses = OneSCluster.fNCrosses;
  fConsumed = OneSCluster.fConsumed;
  Int_t i;
  for (i = 0; i< fNCrosses ; i++)
   {
     fCrossedClusterIndexes[i] = OneSCluster.fCrossedClusterIndexes[i];
   }
  for (i = 0; i< fNDigits ; i++)
  {
    fDigitsIndex[i]=OneSCluster.fDigitsIndex[i];
  }
  return;
}	

/*************************************************************************/
AliITSclusterSSD& AliITSclusterSSD::operator=(const AliITSclusterSSD & OneSCluster)
{
  // assignment operator
  if (this == &OneSCluster) return *this;
  fNDigits = OneSCluster.fNDigits;
  fSide=OneSCluster.fSide;
  fDigits=OneSCluster.fDigits;
  fDigitsIndex = new TArrayI(fNDigits);
  fLeftNeighbour  = OneSCluster.fLeftNeighbour;
  fRightNeighbour = OneSCluster.fRightNeighbour;
  fTotalSignal =-1;
  fNTrack      = -1;
  fNCrosses = OneSCluster.fNCrosses;
  fConsumed = OneSCluster.fConsumed;
  Int_t i;
  for (i = 0; i< fNCrosses ; i++)
   {
     fCrossedClusterIndexes[i] = OneSCluster.fCrossedClusterIndexes[i];
   }
  for (i = 0; i< fNDigits ; i++)
  {
    fDigitsIndex[i]=OneSCluster.fDigitsIndex[i];
  }
  return *this;
}	

/*************************************************************************/
Int_t AliITSclusterSSD::SplitCluster(Int_t where, Int_t *outdigits)
{
//This methods generate data necessery to make new object of this class
//I choosen this way, because methods TClonesArray::Add* dont work
//so I have to use constraction: new (a[i]) Creator(params...);
//where 'a' is a TClonesArray 
//This method generate params - see AliITSmoduleSSD::SplitCluster;
		
	
  Int_t tmp = fNDigits;
  Int_t ind = 0;
    outdigits[ind++]=(*fDigitsIndex)[where];
                     //coping border strip (it is shared by this two clusters)
    for (Int_t i = (where+1); i < tmp; i++)
     {
       outdigits[ind++]=(*fDigitsIndex)[i];  //"moving" strips from this to the new one 
       (*fDigitsIndex)[i]=-1;   
       fNDigits--;   //deleting strips from this cluster
     } 
  return ind;		
}
		
/*******************************************************************/
Int_t AliITSclusterSSD::GetDigitStripNo(Int_t digit)
{
  // return strip no of a digit
	if (digit<0) return -1;
	return (digit>(fNDigits-1))?-1 :
            ((AliITSdigitSSD*)((*fDigits)[(*fDigitsIndex)[digit]]))->GetStripNumber();
}
/************************************************************/
Int_t AliITSclusterSSD::GetDigitSignal(Int_t digit)
{
  // returns digit signal
  Int_t index,signal;
  if (digit<0||digit>=fNDigits) return -1;
  index  = (*fDigitsIndex)[digit];
  signal = ((AliITSdigitSSD*)((*fDigits)[index]))->GetSignal();
  /*
  if(signal>1.e5) printf("GetDigitSignal: digit %d index %d signal %d\n",
			 digit,index, signal);
  */
  return  signal;
}

/***********************************************************/
void  AliITSclusterSSD::AddCross(Int_t clIndex)
{
  // add cluster cross to list of cluster crosses
  
	(*fCrossedClusterIndexes)[fNCrosses++] = clIndex;
}
/***********************************************************/

Int_t AliITSclusterSSD::GetCross(Int_t crIndex)
{
  // return crossing cluster
  
  return ((crIndex>-1)&&(crIndex<fNCrosses))?(*fCrossedClusterIndexes)[crIndex]:-1;
}
/***********************************************************/
Double_t AliITSclusterSSD::CentrOfGravity()
{
  // return center of gravity of the cluster
  Float_t ret=0;
  
  if (fLeftNeighbour) ret+=(GetDigitStripNo(0)*0.5*GetDigitSignal(0));
      else ret+=(GetDigitStripNo(0)*GetDigitSignal(0));
  
  if (fRightNeighbour) ret+=(GetDigitStripNo(fNDigits -1)*0.5*GetDigitSignal(fNDigits -1));
      else ret+=(GetDigitStripNo(fNDigits -1)*GetDigitSignal(fNDigits-1));
      
  for (Int_t i=1;i<fNDigits-1;i++)
    {
      ret +=GetDigitStripNo(i)*GetDigitSignal(i);
    }
  
  if (fTotalSignal<0) GetTotalSignal();
	
    return (ret/fTotalSignal);
}

/***********************************************************/
Float_t AliITSclusterSSD::GetTotalSignal()
{
  // return total signal
  
  if(fTotalSignal <0)
    {
      fTotalSignal=0;
      if (fNDigits ==1)  {
          fTotalSignal = (Float_t)GetDigitSignal(0);
	  //printf("1 digit: signal %d \n",GetDigitSignal(0)); 
	  return fTotalSignal;
      }

      if (fLeftNeighbour) fTotalSignal += (Float_t)(0.5*GetDigitSignal(0));
      else fTotalSignal += (Float_t) GetDigitSignal(0);
      //printf("GetTotalSignal :i  DigitSignal %d %d \n",0,GetDigitSignal(0)); 
		
      if (fRightNeighbour) fTotalSignal += (Float_t)(0.5*GetDigitSignal(fNDigits -1));
      else fTotalSignal += (Float_t)GetDigitSignal(fNDigits-1);
      //printf("GetTotalSignal :i  DigitSignal %d %d \n",fNDigits -1,GetDigitSignal(fNDigits -1)); 
				
      for (Int_t i = 1;i<fNDigits -1;i++) 
        {
      	  fTotalSignal += (Float_t)GetDigitSignal(i);
	  //printf("GetTotalSignal :i  DigitSignal %d %d \n",i,GetDigitSignal(i)); 
      	}
      //printf("GetTotalSignal: fNDigits %d fTotalSignal %.0f \n",fNDigits,fTotalSignal); 
    }
   return fTotalSignal;
}
/***********************************************************/

Float_t  AliITSclusterSSD::GetTotalSignalError()
{
  // return the error on the signal
  Float_t err =0;
  for (Int_t i =1; i<fNDigits -1; i++)
    {
      err+=0.1*GetDigitSignal(i);   
    } 
  if (GetLeftNeighbour())
   {
    err+=GetDigitSignal(0);
   }
  else
   {
    err+=0.1*GetDigitSignal(0);
   } 
  if (GetRightNeighbour())
   {
    err+=GetDigitSignal(fNDigits -1);
   }
  else
   {
    err+=0.1*GetDigitSignal(fNDigits -1);
   }
  return err; 
   
}

/***********************************************************/

void AliITSclusterSSD::DelCross(Int_t index)
{
  // remove cross clusters from the list of cross clusters
Int_t i,j; //iterators

for (i =0;i<fNCrosses;i++)
 {
  if ((*fCrossedClusterIndexes)[i] == index)
   {
     for (j=i;j<fNCrosses-1;j++)
      {
        (*fCrossedClusterIndexes)[j]=(*fCrossedClusterIndexes)[j+1];
      }
     fNCrosses--;
     return; 
   }
 }

}
/***********************************************************/

Int_t  *AliITSclusterSSD::GetTracks(Int_t &nt)
{
  // return the track number of the cluster
  Int_t *tidx=0;
  Int_t i, bit;
  nt=0;

  
   
   fNTrack =0;
   for (i=0;i<10;i++)
    {
     fTrack[i] = 0;
    }
   
   tidx=GetDigit(0)->GetTracks();
	 
   for (i = 0; i<3;i++)
   {
    fTrack[i]=tidx[i];
    if (fTrack[i] != 0) fNTrack++;
   }
   for (i = 1; i<fNDigits; i++)
   {
    tidx=GetDigit(i)->GetTracks();
    for (Int_t j = 0; j<3;j++)
    {
     bit = 1;
     if (tidx[j]==0) break;
     for (Int_t k = 0; k < fNTrack;k++)
      {
       if (tidx[j]==fTrack[k]) bit =0;
      }
     if (bit) fTrack[fNTrack++]=tidx[j];	
    }
   }
	
 
  if (fNTrack > 10)
   {
     cout<<"\n\n Error AliITSclusterSSD::GetTracks  OUT  "<<fNDigits<<"   "<<fNTrack<<"\n\n\n\n\n";
     
   }

 
 nt = fNTrack;
 if(!nt) return 0;
 return &(fTrack[0]);
}
/***********************************************************/

Double_t AliITSclusterSSD::GetPosition()
{
  // return position of the cluster
  Float_t ret;
  switch(fNDigits)
   {
     case  1:
       ret = GetDigitStripNo(0);
       break;
     case  2:
       ret = EtaAlgorithm();
       break;       
     default:
       ret = CentrOfGravity();   
   }
  return ret;
}

/***********************************************************/

Double_t AliITSclusterSSD::EtaAlgorithm()
{
  // algorithm for determing cluster position
  if (fNDigits != 2) return -1;

  Int_t strip1  = GetDigit(0)->GetStripNumber(); 
  Int_t strip2  = GetDigit(1)->GetStripNumber();
  Int_t signal1 = GetDigit(0)->GetSignal();
  Int_t signal2 = GetDigit(1)->GetSignal();

  Double_t eta;
  
 
  if (strip1<strip2)
   {
    eta = ((Double_t)signal2)/((Double_t)(signal1+signal2));
    if (eta<0.04) return strip1;
    if (eta>0.96) return strip2;
    return (strip1 + 0.43478261*eta + 0.2826087);   
   }
  else
  {
    eta = ((Double_t)signal1)/((Double_t)(signal1+signal2));
    if (eta<0.04) return strip2;
    if (eta>0.96) return strip1;
    return (strip2 + 0.43478261*eta + 0.2826087);   
  }

 
}

Double_t  AliITSclusterSSD::GetPositionError()
{
  // return the position error
 return (GetNumOfDigits()+1)/2;
}

Bool_t AliITSclusterSSD::IsCrossingWith(Int_t idx)
{
  // return the cluster to which he crosses
 for (Int_t i =0; i< fNCrosses;i++)
  {
    if (GetCross(i) == idx) return kTRUE;
  }
 return kFALSE;
}
