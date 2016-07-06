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

#include <TArray.h>

#include "AliTRDarraySignal.h"
#include "AliTRDfeeParam.h"

ClassImp(AliTRDarraySignal)

Short_t *AliTRDarraySignal::fgLutPadNumbering = 0x0;

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal()
                  :TObject()
                  ,fNdet(0)
                  ,fNrow(0)
                  ,fNcol(0)
         	  ,fNumberOfChannels(0)
                  ,fNtime(0)
                  ,fNdim(0)  
                  ,fSignal(0)
{

  //
  // AliTRDarraySignal default constructor
  //

  CreateLut();
	   
}

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal(Int_t nrow, Int_t ncol,Int_t ntime)
                  :TObject()
                  ,fNdet(0)
                  ,fNrow(0)
                  ,fNcol(0)
	          ,fNumberOfChannels(0)
                  ,fNtime(0)
                  ,fNdim(0)
                  ,fSignal(0)
{
  //
  // AliTRDarraySignal constructor
  //

  CreateLut(); 
  Allocate(nrow,ncol,ntime);

}

//_______________________________________________________________________
AliTRDarraySignal::AliTRDarraySignal(const AliTRDarraySignal &d)
                  :TObject()
		  ,fNdet(d.fNdet)
		  ,fNrow(d.fNrow)
		  ,fNcol(d.fNcol)
                  ,fNumberOfChannels(d.fNumberOfChannels)
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

  delete [] fSignal;
  fSignal=0;  

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
  fNumberOfChannels = d.fNumberOfChannels;
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
  // Row*NumberOfNecessaryMCMs*ADCchannelsInMCM*Time
  // To be consistent with AliTRDarrayADC
  //

  fNrow=nrow;
  fNcol=ncol;
  fNtime=ntime;
  Int_t adcchannelspermcm = AliTRDfeeParam::GetNadcMcm(); 
  Int_t padspermcm = AliTRDfeeParam::GetNcolMcm(); 
  Int_t numberofmcms = fNcol/padspermcm; 
  fNumberOfChannels = numberofmcms*adcchannelspermcm;
  fNdim = nrow*fNumberOfChannels*ntime;
  if (fSignal)   
    {
      delete [] fSignal;
    }
  fSignal = new Float_t[fNdim];
  memset(fSignal,0,sizeof(Float_t)*fNdim);

}

//_______________________________________________________________________
Int_t AliTRDarraySignal::GetOverThreshold(Float_t threshold) const
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
  Int_t k=0;

  Int_t *longArr = new Int_t[fNdim];  

  if(longArr) 
    {

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
      Float_t* buffer = new Float_t[newDim];
      Int_t counterTwo=0;

      if(buffer)
        {

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

          delete [] buffer;
          buffer=0;

        } 

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

  if(fSignal)
    {

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

      Int_t dimexp=0;
      Int_t *longArr = new Int_t[fNdim];

      if(longArr)
	{

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
          Int_t contaexp =0;    
          Int_t h=0;
          Float_t* bufferE = new Float_t[dimexp];

          if(bufferE)
	    {

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
              delete [] fSignal;
              fSignal=0;
              fSignal = new Float_t[dimexp];
              fNdim = dimexp;
              for(Int_t i=0; i<dimexp; i++)
                {
                  fSignal[i] = bufferE[i]; 
                }

              delete [] bufferE;

	    }

          delete [] longArr;

        }

    }

}
//________________________________________________________________________________
void AliTRDarraySignal::Reset()
{
  //
  // Reset the array, the old contents are deleted
  // The array keeps the same dimensions as before
  //

  memset(fSignal,0,sizeof(Float_t)*fNdim);

}


//________________________________________________________________________________
void AliTRDarraySignal::CreateLut()
{
  //
  // Initializes the Look Up Table to relate
  // pad numbering and mcm channel numbering
  //

  if(fgLutPadNumbering)  return;
  
   fgLutPadNumbering = new Short_t[AliTRDfeeParam::GetNcol()];
   memset(fgLutPadNumbering,0,sizeof(Short_t)*AliTRDfeeParam::GetNcol());

  for(Int_t mcm=0; mcm<8; mcm++)
    {
      Int_t lowerlimit=0+mcm*18;
      Int_t upperlimit=18+mcm*18;
      Int_t shiftposition = 1+3*mcm;
      for(Int_t index=lowerlimit;index<upperlimit;index++)
	{
	  fgLutPadNumbering[index]=index+shiftposition;
	}
    }
}
