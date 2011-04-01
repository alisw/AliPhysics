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

/* $Id: AliTRDarrayDictionary.cxx 25392 2008-04-23 19:40:29Z cblume $ */

/////////////////////////////////////////////////////////
//                                                     //
// Container Class for Dictionary Info                 //
//                                                     //
// Author:                                             //
//   Hermes Leon Vargas (hleon@ikf.uni-frankfurt.de)   //
//                                                     //
/////////////////////////////////////////////////////////

#include <TArray.h>

#include "AliTRDarrayDictionary.h"
#include "AliTRDfeeParam.h"

ClassImp(AliTRDarrayDictionary)

Short_t *AliTRDarrayDictionary::fgLutPadNumbering = 0x0;

//________________________________________________________________________________
AliTRDarrayDictionary::AliTRDarrayDictionary()
                      :TObject()
                      ,fNdet(0)
                      ,fNrow(0)
                      ,fNcol(0)
	              ,fNumberOfChannels(0)
                      ,fNtime(0)
                      ,fNDdim(0)
                      ,fDictionary(0)
                      ,fFlag(kFALSE)
{
  //
  // AliTRDarrayDictionary default contructor
  //

  CreateLut();

}

//________________________________________________________________________________
AliTRDarrayDictionary::AliTRDarrayDictionary(Int_t nrow, Int_t ncol, Int_t ntime)
                      :TObject()
		      ,fNdet(0)
                      ,fNrow(0)
                      ,fNcol(0)
          	      ,fNumberOfChannels(0)
                      ,fNtime(0)
		      ,fNDdim(0)
		      ,fDictionary(0)
                      ,fFlag(kFALSE)

{
  //
  // AliTRDarrayDictionary contructor
  //

  CreateLut();
  Allocate(nrow,ncol,ntime);

}

//________________________________________________________________________________
AliTRDarrayDictionary::AliTRDarrayDictionary(const AliTRDarrayDictionary &a)
                      :TObject()
		      ,fNdet(a.fNdet)
                      ,fNrow(a.fNrow)
                      ,fNcol(a.fNcol)
	              ,fNumberOfChannels(a.fNumberOfChannels)
                      ,fNtime(a.fNtime)
		      ,fNDdim(a.fNDdim)
		      ,fDictionary(0)
                      ,fFlag(a.fFlag)
{
  //
  // AliTRDarrayDictionary copy constructor
  //

  fDictionary = new Int_t[fNDdim];
  for(Int_t i=0; i<fNDdim; i++)
    {
      fDictionary[i]=a.fDictionary[i];
    }

}

//________________________________________________________________________________
AliTRDarrayDictionary::~AliTRDarrayDictionary()
{
  //
  //   AliTRDarrayDictionary destructor
  //

  if(fDictionary)
    {
      delete [] fDictionary;
      fDictionary=0;
    }

}

//________________________________________________________________________________
AliTRDarrayDictionary &AliTRDarrayDictionary::operator=(const AliTRDarrayDictionary &a)
{
  //
  // Assignment operator
  //

  if(this==&a)
    {
      return *this;
    }

  if(fDictionary)
    {
      delete [] fDictionary;
    }
  fNdet=a.fNdet;
  fNDdim=a.fNDdim;
  fNrow=a.fNrow;
  fNcol=a.fNcol;
  fNumberOfChannels = a.fNumberOfChannels;
  fNtime=a.fNtime;
  fDictionary = new Int_t[fNDdim];
  for(Int_t i=0; i<fNDdim; i++)
    {
      fDictionary[i]=a.fDictionary[i];
    }
  fFlag=a.fFlag;
  return *this;

}

//________________________________________________________________________________
void AliTRDarrayDictionary::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocates memory for the dictionary array with dimensions
  // Row*NumberOfNecessaryMCMs*ADCchannelsInMCM*Time
  // To be consistent with AliTRDarrayADC
  // Object initialized to -1
  //

  fNrow=nrow;
  fNcol=ncol;
  fNtime=ntime;
  Int_t adcchannelspermcm = AliTRDfeeParam::GetNadcMcm(); 
  Int_t padspermcm = AliTRDfeeParam::GetNcolMcm(); 
  Int_t numberofmcms = fNcol/padspermcm;
  fNumberOfChannels = numberofmcms*adcchannelspermcm;
  fNDdim=nrow*fNumberOfChannels*ntime;
  if(fDictionary)
    {
      delete [] fDictionary;
      fDictionary=0;
    }
  fDictionary = new Int_t[fNDdim];
  memset(fDictionary,-1,sizeof(Int_t)*fNDdim);

}

//________________________________________________________________________________
void AliTRDarrayDictionary::Compress()
{
  //
  // Compress the array
  //


  //  AliDebug(1,"Compressing");
  Int_t counter=0;
  Int_t newDim=0;
  Int_t j;                 
  Int_t r=0;
  Int_t k=0;

  Int_t *longArr = new Int_t[fNDdim];  

  if(longArr) 
    {

      memset(longArr,0,sizeof(Int_t)*fNDdim);

      for(Int_t i=0;i<fNDdim; i++)
        {
          j=0;
          if(fDictionary[i]==-1)
   	    {
	      for(k=i;k<fNDdim;k++)
	        {
	          if(fDictionary[k]==-1)
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
      for(Int_t i=0; i<fNDdim;i++)
        {
          if(longArr[i]!=0)  
	    {
	      counter=counter+longArr[i]-1;
	    }
        }
      newDim=fNDdim-counter;   //Size of the compressed array

      //Fill the buffer of the compressed array
      Int_t* buffer = new Int_t[newDim];
      Int_t counterTwo=0;
      Int_t g=0;
      if(buffer)
        {
          for(Int_t i=0; i<newDim; i++)
            {
              if(counterTwo<fNDdim)
	        {
	          if(fDictionary[counterTwo]!=-1)
	            {
	              buffer[i]=fDictionary[counterTwo];
	            }
	          if(fDictionary[counterTwo]==-1)
	            {
	              buffer[i]=-(longArr[g]);
	              counterTwo=counterTwo+longArr[g]-1;
	              g++;
	            }  
	          counterTwo++;
	        }
            }

          //Copy the buffer
          if(fDictionary)
            {
              delete [] fDictionary;
              fDictionary=0;
            }
	  fDictionary = buffer;
          fNDdim = newDim;
        }
    
      delete [] longArr;
      longArr=0;
    }
  fFlag=kTRUE; // flag to signal compression

}

//________________________________________________________________________________
void AliTRDarrayDictionary::Expand()
{
  //  
  //  Expand the array
  //  

  Int_t dimexp=0;
  if(!IsCompressed())
    return;

  if(fDictionary&&fNDdim==1)
    { 
      dimexp = -fDictionary[0];	
      delete [] fDictionary;
      fDictionary=0;
      fDictionary = new Int_t[dimexp];
      fNDdim = dimexp;
      // Re-initialize the array
      memset(fDictionary,-1,sizeof(Int_t)*dimexp);
      return;
    }

  Int_t *longArr = new Int_t[fNDdim];
  if(longArr && fDictionary)
    {
      //Initialize the array
      memset(longArr,0,sizeof(Int_t)*fNDdim);

      Int_t r2=0;
      for(Int_t i=0; i<fNDdim;i++)
        {
          if((fDictionary[i]<0)&&(fDictionary[i]!=-1))  
	    {
	      longArr[r2]=-fDictionary[i]; 
	      r2++;
	    }
        }

      //Calculate new dimensions
      for(Int_t i=0; i<fNDdim;i++)
        {
          if(longArr[i]!=0)      
	    dimexp=dimexp+longArr[i]-1;
	  if(longArr[i]==0)
	    break;
        }
      dimexp=dimexp+fNDdim;  

      //Write in the buffer the new array
      Int_t contaexp =0;    
      Int_t h=0;
      Int_t* bufferE = new Int_t[dimexp];
      memset(bufferE,-1,sizeof(Int_t)*dimexp);

      if(bufferE)
	{

          for(Int_t i=0; i<dimexp; i++)
            {
              if(fDictionary[contaexp]>=-1)  
	        {
	          bufferE[i]=fDictionary[contaexp];
	        }
              if(fDictionary[contaexp]<-1)  
	        {
	          i=i+longArr[h]-1;
	          h++;
	        }
              contaexp++;
            }

          //Copy the buffer
          delete [] fDictionary;
          fDictionary=bufferE;
          fNDdim = dimexp;
	}
    }
  if (longArr)
    {
      delete [] longArr; 
    }

  fFlag=kFALSE; // flag to signal that is not compressed anymore

}
//________________________________________________________________________________
void AliTRDarrayDictionary::Reset()
{
  //
  // Reset the array, the old contents are deleted
  // and the data array elements are set to zero.
  //

  memset(fDictionary,0,sizeof(Int_t)*fNDdim);

}
//________________________________________________________________________________
Int_t AliTRDarrayDictionary::GetData(Int_t nrow, Int_t ncol, Int_t ntime) const
{
  //
  // Get the data using the pad numbering.
  // To access data using the mcm scheme use instead
  // the method GetDataByAdcCol
  //

  Int_t corrcolumn = fgLutPadNumbering[ncol];

  return fDictionary[(nrow*fNumberOfChannels+corrcolumn)*fNtime+ntime];

}
//________________________________________________________________________________
void AliTRDarrayDictionary::SetData(Int_t nrow, Int_t ncol, Int_t ntime, Int_t value)
{
  //
  // Set the data using the pad numbering.
  // To write data using the mcm scheme use instead
  // the method SetDataByAdcCol
  //

  Int_t colnumb = fgLutPadNumbering[ncol];

  fDictionary[(nrow*fNumberOfChannels+colnumb)*fNtime+ntime]=value;

}

//________________________________________________________________________________
void AliTRDarrayDictionary::CreateLut()
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
	  fgLutPadNumbering[index]= index+shiftposition;
	}
    }
}
