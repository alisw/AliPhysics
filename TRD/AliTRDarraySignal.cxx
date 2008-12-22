/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

/* $Id: AliTRDarraySignal.cxx 25392 2008-04-23 19:40:29Z cblume $ */

/////////////////////////////////////////////////////////
//                                                     //
// Container Class for Signals                         //
//                                                     //
// Author:                                             //
//   Hermes Leon Vargas (hleon@ikf.uni-frankfurt.de)   //
//                                                     //
/////////////////////////////////////////////////////////

#include "AliTRDarraySignal.h"
//#include "AliLog.h"
#include "TArray.h" //for memset

ClassImp(AliTRDarraySignal)

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal()
                  :TObject()
                  ,fNdet(0)
                  ,fNrow(0)
                  ,fNcol(0)
                  ,fNtime(0)
                  ,fNdim(0)  
                  ,fSignal(0)
{

  //
  // AliTRDarraySignal default constructor
  //
	   
}

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal(Int_t nrow, Int_t ncol,Int_t ntime)
                  :TObject()
                  ,fNdet(0)
                  ,fNrow(0)
                  ,fNcol(0)
                  ,fNtime(0)
                  ,fNdim(0)
                  ,fSignal(0)
{
  //
  // AliTRDarraySignal constructor
  //
 
  Allocate(nrow,ncol,ntime);

}

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal(const AliTRDarraySignal &d)
                  :TObject()
		  ,fNdet(d.fNdet)
		  ,fNrow(d.fNrow)
		  ,fNcol(d.fNcol)
		  ,fNtime(d.fNtime)
		  ,fNdim(d.fNdim)
		  ,fSignal(0)
{
  //
  // AliTRDarraySignal copy constructor
  //

  fSignal = new Float_t[fNdim];
  memcpy(fSignal, d.fSignal, fNdim*sizeof(Float_t));

}

//_______________________________________________________________________
AliTRDarraySignal::~AliTRDarraySignal()
{
  //
  // AliTRDarraySignal destructor
  //

  if (fSignal)   
    {
      delete [] fSignal;
      fSignal=0;  
    }

}

//_______________________________________________________________________
AliTRDarraySignal &AliTRDarraySignal::operator=(const AliTRDarraySignal &d)
{
  //
  // Assignment operator
  //

  if (this==&d) 
    {
      return *this;
    }

  if (fSignal)
    {
      delete [] fSignal;
    }
  fNdet=d.fNdet;
  fNrow=d.fNrow;
  fNcol=d.fNcol;
  fNtime=d.fNtime;
  fNdim=d.fNdim;
  fSignal = new Float_t[fNdim];
  memcpy(fSignal,d.fSignal, fNdim*sizeof(Float_t));

  return *this;

}

//_______________________________________________________________________
void AliTRDarraySignal::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocates memory for an AliTRDarraySignal object with dimensions 
  // Row*Col*Time
  //

  fNrow=nrow;
  fNcol=ncol;
  fNtime=ntime;
  fNdim = nrow*ncol*ntime;
  if (fSignal)   
    {
      delete [] fSignal;
    }
  fSignal = new Float_t[fNdim];
  memset(fSignal,0,sizeof(Float_t)*fNdim);

}

//_______________________________________________________________________
Int_t AliTRDarraySignal::GetOverThreshold(Float_t threshold)
{
  //
  // Get the number of entries over the threshold 
  //

  Int_t counter=0;
  for(Int_t i=0; i<fNdim; i++)
    {
      if(fSignal[i]>threshold)
	{
	  counter++;
	}
    }
  return counter;

}

//_______________________________________________________________________
void AliTRDarraySignal::Compress(Float_t minval)
{
  //
  // Compress the array, setting values equal or 
  // below minval to zero (minval>=0)
  //

  Int_t counter=0;
  Int_t newDim=0;
  Int_t j;                 
  Int_t r=0;
  Int_t *longArr;
  longArr = new Int_t[fNdim];  
  Int_t k=0;

  //Initialize the array
  memset(longArr,0,sizeof(Int_t)*fNdim);

  for(Int_t i=0;i<fNdim; i++)
    {
      j=0;
      if(fSignal[i]<=minval) 
	{
	  for(k=i;k<fNdim;k++)
	    {
	      if(fSignal[k]<=minval)
		{
		  j=j+1;
		  longArr[r]=j;
		}
	      else
		{
		  break;
		}
	    } 
	  r=r+1;          
	}
      i=i+j;
    }

  //Calculate the size of the compressed array
  for(Int_t i=0; i<fNdim;i++)
    {
      if(longArr[i]!=0)   
	{
	  counter=counter+longArr[i]-1;
	}
    }
  newDim=fNdim-counter;   //New dimension

  //Fill the buffer of the compressed array
  Float_t* buffer;
  buffer = new Float_t[newDim];
  Int_t counterTwo=0;

  //Write the new array
  Int_t g=0;
  for(Int_t i=0; i<newDim; i++)
    {
      if(counterTwo<fNdim)
	{
	  if(fSignal[counterTwo]>minval)   
	    {
	      buffer[i]=fSignal[counterTwo];
	    }
	  if(fSignal[counterTwo]<=minval)   
	    {
	      buffer[i]=-(longArr[g]);
	      counterTwo=counterTwo+longArr[g]-1;
	      g++;
	    }  
	  counterTwo++;
	}
    }

  //Copy the buffer
  if(fSignal)
    {
      delete [] fSignal;
      fSignal=0;
    }
  fSignal = new Float_t[newDim];
  fNdim = newDim;
  for(Int_t i=0; i<newDim; i++)
    {
      fSignal[i] = buffer[i]; 
    }
  if(buffer)
    {
      delete [] buffer;
      buffer=0;
    } 
  if(longArr) 
    {
      delete [] longArr;
      longArr=0;
    }

}

//_______________________________________________________________________
void AliTRDarraySignal::Expand()
{
  //
  // Expand the array
  //

  //Check if the array has not been already expanded
  Int_t verif=0;
  for(Int_t i=0; i<fNdim; i++)
    {
      if(fSignal[i]<0)
	{
	  verif++;
	}
    }

  if(verif==0)
    {
      return;
    }

  Int_t *longArr; 
  longArr = new Int_t[fNdim];
  Int_t dimexp=0;
  memset(longArr,0,sizeof(Int_t)*fNdim);

  Int_t r2=0;
  for(Int_t i=0; i<fNdim;i++)
    {
      if(fSignal[i]<0)  
	{
	  longArr[r2]=(Int_t)(-fSignal[i]); 
	  r2++;
	}
    }

  //Calculate new dimensions
  for(Int_t i=0; i<fNdim;i++)
    {
      if(longArr[i]!=0)      
	{
	  dimexp=dimexp+longArr[i]-1;
	}
    }
  dimexp=dimexp+fNdim;   //Dimension of the expanded array

  //Write in the buffer the new array
  Float_t* bufferE;
  bufferE = new Float_t[dimexp];
  Int_t contaexp =0;    
  Int_t h=0;
  for(Int_t i=0; i<dimexp; i++)
    {
      if(fSignal[contaexp]>0)  
	{
	  bufferE[i]=fSignal[contaexp];
	}
      if(fSignal[contaexp]<0)  
 	{
	  for(Int_t j=0; j<longArr[h];j++)
	    {
	      bufferE[i+j]=0;
	    }
	  i=i+longArr[h]-1;
	  h++;
	}
      contaexp++;
    }

  //Copy the buffer
  if(fSignal)
    {
      delete [] fSignal;
      fSignal=0;
    }

  fSignal = new Float_t[dimexp];
  fNdim = dimexp;
  for(Int_t i=0; i<dimexp; i++)
    {
      fSignal[i] = bufferE[i]; 
    }

  if(bufferE) delete [] bufferE;
  if(longArr) delete [] longArr;

}
