//Author:        Anders Strand Vestbo
//Last Modified: 14.12.2000

#include <iostream.h>

#include "AliL3Logging.h"
#include "AliL3ClustFinder.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3SpacePointData.h"

//
//Simple clusterfinder. Initial version by Tonko.
//

ClassImp(AliL3ClustFinder)

AliL3ClustFinder::AliL3ClustFinder()
{
 
  fDeconvTime=false;
  fDeconvPad=false;
  fXYErr = 0;
  fZErr = 0;
}

AliL3ClustFinder::AliL3ClustFinder(AliL3Transform *transform)
{

  fTransform = transform;

  fThreshold = 20;
  fDeconvTime = true;
  fDeconvPad = true;
}


AliL3ClustFinder::~AliL3ClustFinder()
{
  //destructor 
}


  void AliL3ClustFinder::SetTransformer( AliL3Transform *transform )
    {
    fTransform = transform;
    }


void AliL3ClustFinder::InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t nmaxpoints)
{
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow = firstrow;
  fLastRow = lastrow;
}


void AliL3ClustFinder::InitSlice(Int_t slice,Int_t patch,Int_t nmaxpoints)
{
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
}

void AliL3ClustFinder::SetOutputArray(AliL3SpacePointData *pt)
{
  
  fSpacePointData = pt;
}

void AliL3ClustFinder::Read(UInt_t ndigits,AliL3DigitRowData *ptr)
{
  fNDigitRowData = ndigits;
  fDigitRowData = ptr;

}

void AliL3ClustFinder::ProcessDigits()
{
  //Loop over rows, and call processrow
  
  
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fDigitRowData;
  
  for(Int_t i=fFirstRow; i<=fLastRow; i++)
    {
      fCurrentRow = i;
      ProcessRow(tempPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData) + tempPt->fNDigit*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
    }

  LOG(AliL3Log::kInformational,"AliL3ClustFinder","Cluster Finder")
    <<AliL3Log::kDec<<"ClusterFinder found "<<fNClusters<<" clusters"<<ENDLOG;
}

void AliL3ClustFinder::ProcessRow(AliL3DigitRowData* tempPt)
{
  
  Int_t cl_found,trackID[3];
  //UInt_t rows;
  
  //Reserve array for found clusters:
  resx aliresx[5000];
  bzero((char *)aliresx,5000*sizeof(resx));
  
  //local results
  resx *r;
//  UInt_t pres1[2500],pres2[2500];
  resx *pres1[2500];
  resx *pres2[2500];
  resx rr_local;
  resx *rr = &rr_local;
  
  //Set variables for each row:
  Int_t pres_cou1,pres_cou2;
//  UInt_t *r1,*r2;
  resx **r1,**r2;
  pres_cou1=pres_cou2=0;
  r1=pres2;
  r2=pres1;
  r=aliresx;
  
  //variables for each pad:
  UInt_t start;
  Int_t start_new=-1;
  Int_t go_on;
  Int_t cl_counter=0,rows;
  Int_t prev_start;
   
  UInt_t last_pad = 123456789;
  
  //Loop over digits in this row:
  for(start=0; start<tempPt->fNDigit; start++)
    {
      AliL3DigitData *bins = tempPt->fDigitData;
                  
      if(bins[start].fPad!=last_pad)
	{//new pad
	  
	  last_pad=bins[start].fPad;
	  
	  if(r2==pres2)
	    {
	      r1=pres2;
	      r2=pres1;
	    }
	  else
	    {
	      r1=pres1;
	      r2=pres2;
	    }
	  pres_cou1=pres_cou2;
	  pres_cou2=0;
	  
	  start_new=-1;

	  cl_counter=0;
	  prev_start=0;
	  
	  
	}//Looping over pads
      
      UInt_t av,charge;
      UInt_t mean;
      UInt_t flags;
      Int_t last_falling=0;
//      UInt_t *ri;
      resx **ri;
      UInt_t tmp_charge;
      
      //Set flag:
      flags=0;
      
      cl_counter++;
      if(cl_counter>MAX_C) {LOG(AliL3Log::kWarning,"AliL3ClustFinder::ProcessRow","Cluster Finder")
      <<AliL3Log::kDec<<"Too many clusters"<<ENDLOG; return;} //too many clusters
            
      
      if(fDeconvTime)
	{
	  
	  //this is a goto label:
	redo: ;
	//divide sequenz in two in case there is a change in slope:
	if(start_new>-1)
	  {
	    start=start_new;
	    start_new=-1;
	  }
	}
      //average and charge
      av=charge=0;
      
      //ok, create a new cluster:
      go_on=1;
      
      if(fDeconvTime)
	{
	  
	  last_falling=flags=0;
	}
      
      if(fDeconvPad)
	{
	  flags=0;
	}
      
      
      //block to introduce new variables:
      {
	UInt_t last_a;
	last_a=0;
	
	//Loop over this sequenz:
	while(1)
	  {
	    UInt_t aa;
	    //Int_t start_temp=start;
	    
	    //get the adc-value:
	    aa=bins[start].fCharge;
	    
	    
	    if(fDeconvTime)
	      {
		//Check if the last pixel in this sequenz is bigger or smaller than this:
		if(aa>last_a)
		  {
		    if(last_falling)
		      {
			start_new=start;
			break;
		      }
		  }
		else last_falling=1;
		last_a=aa;
		
	      }
	    //sum of total charge of this cluster on this pad
	    charge +=aa;
	    //this one is needed to determine the mean over all pads:
	    //av += start * (aa);
	    av+=bins[start].fTime*(aa);
	    
	    if((bins[start+1].fTime)!=(bins[start].fTime+1) || 
	       ((bins[start+1].fPad)!=(bins[start].fPad))) break;
	    
	    start++;
	  }//Loop over sequenz
      }//block of new variables
      
      if(fDeconvTime)
	{
	  if(start_new>0) flags = FLAG_DOUBLE_T;
	}
      //calculate mean in time direction for this sequenz:
      if(charge)
	{
	  //mean=cog for this sequenz:
	  mean=av/charge;
	}
      else
	{
	  charge=1;
	  mean=1;
	}
      //get the pointer to results on previous pad(s)
      ri=r1;
      
      //Center of gravity in pad direction
      //tmp_charge=charge*bins[start-1].fPad;
      tmp_charge=charge*bins[start].fPad;
 
      //compare with prevously stored results (previous pads)
      for(int k=0; k<pres_cou1; k++)
	{
	  //variable to store difference:
	  Int_t v;
	  
	  //structure that contains means (+more) of last pads(s)
	  resx *rr_tmp;
	  
	  //go get it
	  rr_tmp = (resx *) *ri;
	  
	  //difference between mean and mean on previous pads
	  v=mean-rr_tmp->mean;
	  
	  
	  if(v<-PARAM1) break;
	  
	  //goto next mean on previous pad
	  ri++;
	  
	  //is there a matching cluster on previous pad:
	  if(v<=PARAM1)
	    {
	      //ok, we got a new one, so we will add the new and the old
	      rr=&rr_local;
	      
	      //Add this function somewhere
	      memcpy7((UInt_t *)rr,(UInt_t *)rr_tmp);
	      
	      /* in case we found another sequenz-section on this pad
		 don't loop over all sequenzes on last pad but start
		 where we found a matching candidate for this sequenz-
		 section, adjust upper limit (=pres_cou1)*/
	      r1=ri;
	      pres_cou1-=1+k;
	      
	      if(fDeconvPad)
		{
		  if(charge>rr->scharge)
		    {
		      if(rr->falling)
			{
			  //previous state was falling
			  flags |= FLAG_DOUBLE_PAD;
			  rr_tmp->flags |= flags;
			  //and create a new one
			  break;
			}
		    }
		  else
		    {
		      rr->falling = 1;
		    }
		  
		}
	      //don't create a new one
	      go_on=0;
	      
	      //store the pointer to the matched in "2"
//	      r2[pres_cou2++] = (UInt_t) rr_tmp;
	      r2[pres_cou2++] =  rr_tmp;
	      
	      //calculate and fill the new means, charge etc. in rr_tmp
	      //charge on this pad to determine whether we cut in pad direction
	      rr->scharge = charge;
	      //unset FLAG_ONEPAD since we have at least 2 pads
	      rr->flags &= (~FLAG_ONEPAD);
	      //summ the charges for de/dx
	      rr->charge += charge;
	      //calculate mean in pad direction
	      rr->pad += tmp_charge;
	      //calculate mean in time direction
	      rr->t+=av;
	      //store the mean of this sequenz on this pad
	      rr->mean = mean;
	      
	      //Copy it back to storage
	      memcpy7((UInt_t *)rr_tmp,(UInt_t *)rr);
	      
	      //we found a matching one so stop looping over pads on previous pad
	      break;
	      
	    }//matching cluster on previous pads
	}//loop: means on previous pad
      
      //here we create from scratch new clusters
      if(go_on)
	{
	  //store this one beacuse it is the first
//	  r2[pres_cou2++]= (UInt_t )r;
	  r2[pres_cou2++]= r;
	  	  
	  mstore((UInt_t *)r,av,tmp_charge,charge,flags|FLAG_ONEPAD,mean,trackID);
	  r++;
	}

      if(fDeconvTime)
	{
	  
	  //maybe we still have to do something
	  if(start_new>=0) goto redo;
	}
      
    }//end of loop over digits
  
  
  //number of clusters found
  cl_found = r-aliresx;
  
  //tot_cl_found+=cl_found;
  
  //there was something on this padrow
  rows++;
   
  WriteClusters(cl_found,aliresx);
  
}

void AliL3ClustFinder::memcpy7(UInt_t *dst, UInt_t *src)
{
  Int_t i;
  for(i=0; i<7; i++)
    {
      *dst++ = *src++;
    }
  return;
}

void AliL3ClustFinder::mstore(UInt_t *r,UInt_t av,UInt_t pad,UInt_t ch,UInt_t flags,UInt_t mean,Int_t *trackID)
{
  resx *rr = (resx *) r;

  rr->mean = mean;
  rr->charge = ch;
  rr->flags =flags;
  rr->scharge = ch;
  rr->pad = pad;
  rr->t = av;
  rr->falling=0;  
  
  rr->trackID[0]=trackID[0];
  rr->trackID[1]=trackID[1]; 
  rr->trackID[2]=trackID[2];
  return;
}

void AliL3ClustFinder::WriteClusters(Int_t ncl,resx *r)
{
/* 
   if(fNClusters >= fMaxNClusters)
   {
   LOG(AliL3Log::kError,"AliL3ClustFinder::WriteClusters","Cluster Finder")
   <<AliL3Log::kDec<<"Too many clusters"<<ENDLOG;
   return;
   }
*/
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  if(!fXYErr || !fZErr)
    {
      LOG(AliL3Log::kError,"AliL3ClustFinder::WriteClusters","Cluster Errors")
	<<"No Errors"<<ENDLOG;
      return;
    }
  
  for(int j=0; j<ncl; j++)
    {
      if(r[j].flags & FLAG_ONEPAD) continue; //discard 1 pad clusters
      if(r[j].charge < fThreshold) continue; //noise cluster
      Float_t xyz[3];
      
      Float_t fpad=(Float_t)r[j].pad/(Float_t)r[j].charge;
      Float_t ftime=(Float_t)r[j].t/(Float_t)r[j].charge;
      
      if(fCurrentRow > 54) {thisrow = fCurrentRow-55; thissector = fCurrentSlice+36;}
      else {thisrow = fCurrentRow; thissector = fCurrentSlice;}
      fTransform->Raw2Local(xyz,thissector,thisrow,fpad,ftime);
      if(xyz[0]==0) LOG(AliL3Log::kError,"AliL3ClustFinder","Cluster Finder")
	<<AliL3Log::kDec<<"Zero cluster"<<ENDLOG;
      if(fNClusters >= fMaxNClusters)
	  {
	  LOG(AliL3Log::kError,"AliL3ClustFinder::WriteClusters","Cluster Finder")
	      <<AliL3Log::kDec<<"Too many clusters"<<ENDLOG;
	  return;
	  }  
      fSpacePointData[counter].fX = xyz[0];
      fSpacePointData[counter].fY = xyz[1];
      fSpacePointData[counter].fZ = xyz[2];
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fXYErr = fXYErr;
      fSpacePointData[counter].fZErr = fZErr;
      fSpacePointData[counter].fID = counter
                  +((fCurrentSlice&0x7f)<<25)+((fCurrentPatch&0x7)<<22);//uli

      
      fNClusters++;
      counter++;

    }

  
}
