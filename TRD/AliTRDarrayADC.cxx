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

/* $Id: AliTRDarrayADC.cxx 25392 2008-04-23 19:40:29Z cblume $ */

////////////////////////////////////////////////////////
//                                                    //
// Container class for ADC values                     //
//                                                    //
// Author:                                            //
// Hermes Leon Vargas  (hleon@ikf.uni-frankfurt.de)   //
//                                                    // 
////////////////////////////////////////////////////////

#include "AliTRDarrayADC.h"
#include "Cal/AliTRDCalPadStatus.h"
#include "AliTRDfeeParam.h"
#include "AliTRDSignalIndex.h"
ClassImp(AliTRDarrayADC)

Short_t *AliTRDarrayADC::fgLutPadNumbering = 0x0;

//____________________________________________________________________________________
AliTRDarrayADC::AliTRDarrayADC()
               :TObject()
               ,fNdet(0)
               ,fNrow(0)
               ,fNcol(0)
	       ,fNumberOfChannels(0)
               ,fNtime(0) 
               ,fNAdim(0)
               ,fADC(0)
{
  //
  // AliTRDarrayADC default constructor
  //

  CreateLut();

}

//____________________________________________________________________________________
AliTRDarrayADC::AliTRDarrayADC(Int_t nrow, Int_t ncol, Int_t ntime)
               :TObject()
	       ,fNdet(0)               
               ,fNrow(0)
               ,fNcol(0)
	       ,fNumberOfChannels(0)
               ,fNtime(0) 
               ,fNAdim(0)
               ,fADC(0)
{
  //
  // AliTRDarrayADC constructor
  //

  CreateLut();
  Allocate(nrow,ncol,ntime);

}

//____________________________________________________________________________________
AliTRDarrayADC::AliTRDarrayADC(const AliTRDarrayADC &b)
               :TObject()
	       ,fNdet(b.fNdet)
               ,fNrow(b.fNrow)
               ,fNcol(b.fNcol)
	       ,fNumberOfChannels(b.fNumberOfChannels)
               ,fNtime(b.fNtime) 
               ,fNAdim(b.fNAdim)
               ,fADC(0)	 
{
  //
  // AliTRDarrayADC copy constructor
  //

  fADC =  new Short_t[fNAdim];
  memcpy(fADC,b.fADC, fNAdim*sizeof(Short_t));  

}

//____________________________________________________________________________________
AliTRDarrayADC::~AliTRDarrayADC()
{
  //
  // AliTRDarrayADC destructor
  //

  if(fADC)
    {
      delete [] fADC;
      fADC=0;
    }

}

//____________________________________________________________________________________
AliTRDarrayADC &AliTRDarrayADC::operator=(const AliTRDarrayADC &b)
{
  //
  // Assignment operator
  //

  if(this==&b)
    {
      return *this;
    }
  if(fADC)
    {
      delete [] fADC;
    }
  fNdet=b.fNdet;
  fNrow=b.fNrow;
  fNcol=b.fNcol;
  fNumberOfChannels = b.fNumberOfChannels;
  fNtime=b.fNtime;
  fNAdim=b.fNAdim;
  fADC = new Short_t[fNAdim];
  memcpy(fADC,b.fADC, fNAdim*sizeof(Short_t));

  return *this;

}

//____________________________________________________________________________________
void AliTRDarrayADC::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocate memory for an AliTRDarrayADC array with dimensions
  // Row*NumberOfNecessaryMCMs*ADCchannelsInMCM*Time
  //
  
  fNrow=nrow;
  fNcol=ncol;
  fNtime=ntime;
  Int_t adcchannelspermcm = AliTRDfeeParam::GetNadcMcm(); 
  Int_t padspermcm = AliTRDfeeParam::GetNcolMcm(); 
  Int_t numberofmcms = fNcol/padspermcm; 
  fNumberOfChannels = numberofmcms*adcchannelspermcm; 
  fNAdim=nrow*fNumberOfChannels*ntime;

  if(fADC)
    {
      delete [] fADC;
    }
  
  fADC = new Short_t[fNAdim];
  memset(fADC,0,sizeof(Short_t)*fNAdim);

}

//____________________________________________________________________________________
Short_t AliTRDarrayADC::GetDataBits(Int_t row, Int_t col, Int_t time) const
{
  //
  // Get the ADC value for a given position: row, col, time
  // Taking bit masking into account
  //
  // Adapted from code of the class AliTRDdataArrayDigits 
  //

  Short_t tempval = GetData(row,col,time);
  // Be aware of manipulations introduced by pad masking in the RawReader
  // Only output the manipulated Value
  CLRBIT(tempval, 10);
  CLRBIT(tempval, 11);
  CLRBIT(tempval, 12);
  return tempval;

}

//____________________________________________________________________________________
UChar_t AliTRDarrayADC::GetPadStatus(Int_t row, Int_t col, Int_t time) const
{
  // 
  // Returns the pad status stored in the pad signal
  //
  // Output is a UChar_t value
  // Status Codes:
  //               Noisy Masking:           2
  //               Bridged Left Masking     8
  //               Bridged Right Masking    8
  //               Not Connected Masking Digits
  //
  // Adapted from code of the class AliTRDdataArrayDigits
  //

  UChar_t padstatus = 0;
  Short_t signal = GetData(row,col,time);
  if(signal > 0 && TESTBIT(signal, 10)){
    if(TESTBIT(signal, 11))
      if(TESTBIT(signal, 12))
	padstatus = AliTRDCalPadStatus::kPadBridgedRight;
      else
	padstatus = AliTRDCalPadStatus::kNotConnected;
    else
      if(TESTBIT(signal, 12))
	padstatus = AliTRDCalPadStatus::kPadBridgedLeft;
      else
	padstatus = AliTRDCalPadStatus::kMasked;
  }

  return padstatus;

}

//____________________________________________________________________________________
void AliTRDarrayADC::SetPadStatus(Int_t row, Int_t col, Int_t time, UChar_t status)
{
  //
  // Setting the pad status into the signal using the Bits 10 to 14 
  // (currently used: 10 to 12)
  //
  // Input codes (Unsigned char):
  //               Noisy Masking:           2
  //               Bridged Left Masking     8
  //               Bridged Right Masking    8
  //               Not Connected Masking    32
  //
  // Status codes: Any masking:             Bit 10(1)
  //               Noisy masking:           Bit 11(0), Bit 12(0)
  //               No Connection masking:   Bit 11(1), Bit 12(0)
  //               Bridged Left masking:    Bit 11(0), Bit 12(1)
  //               Bridged Right masking:   Bit 11(1), Bit 12(1)
  // 
  // Adapted from code of the class AliTRDdataArrayDigits
  //

  Short_t signal = GetData(row,col,time);

  // Only set the Pad Status if the signal is > 0
  if(signal > 0)
    {
      switch(status)
	{
	case AliTRDCalPadStatus::kMasked:
	  SETBIT(signal, 10);
	  CLRBIT(signal, 11);
	  CLRBIT(signal, 12);
	  break;
	case AliTRDCalPadStatus::kNotConnected:
	  SETBIT(signal, 10);
	  SETBIT(signal, 11);
	  CLRBIT(signal, 12);
	  break;
	case AliTRDCalPadStatus::kPadBridgedLeft:
	  SETBIT(signal, 10);
	  CLRBIT(signal, 11);
	  SETBIT(signal, 12);
	  break;
	case AliTRDCalPadStatus::kPadBridgedRight:
	  SETBIT(signal, 10);
	  SETBIT(signal, 11);
	  SETBIT(signal, 12);
	default:
	  CLRBIT(signal, 10);
	  CLRBIT(signal, 11);
	  CLRBIT(signal, 12);
	}
      SetData(row, col, time, signal);
    }

}

//____________________________________________________________________________________
Bool_t AliTRDarrayADC::IsPadCorrupted(Int_t row, Int_t col, Int_t time)
{
  // 
  // Checks if the pad has any masking as corrupted (Bit 10 in signal set)
  // 
  // Adapted from code of the class AliTRDdataArrayDigits
  //

  Short_t signal = GetData(row,col,time);
  return (signal > 0 && TESTBIT(signal, 10)) ? kTRUE : kFALSE;

}

//____________________________________________________________________________________
void AliTRDarrayADC::Compress()
{
  //
  // Compress the array
  //

  Int_t counter=0;
  Int_t newDim=0;
  Int_t j;                  
  Int_t l;                  
  Int_t r=0;                
  Int_t s=0;                
  Int_t *longm;            
  longm = new Int_t[fNAdim];  
  Int_t *longz;            
  longz = new Int_t[fNAdim];
  Int_t k=0;
  memset(longz,0,sizeof(Int_t)*fNAdim);
  memset(longm,0,sizeof(Int_t)*fNAdim);

  for(Int_t i=0;i<fNAdim; i++)
    {
      j=0;
      if(fADC[i]==-1)
	{
	  for(k=i;k<fNAdim;k++)
	    {
	      if((fADC[k]==-1)&&(j<16000))   
		{
		  j=j+1;
		  longm[r]=j;                
		}
	      else
		{
		  break;
		}
	    }
	  r=r+1;            
	}
      l=16001;
      if(fADC[i]==0)
	{
	  for(k=i;k<fNAdim;k++)
	    {
	      if((fADC[k]==0)&&(l<32767))     
		{                             
		  l=l+1;
		  longz[s]=l;                
		}
	      else
		{
		  break;
		}
	    }
	  s=s+1;         
	}
      if(fADC[i]>0)
	{
	  i=i+1;
	}
      i=i+j+(l-16001-1); 
    }

  //Calculate the size of the compressed array
  for(Int_t i=0; i<fNAdim;i++)
    {
      if(longm[i]!=0)   
	{
	  counter=counter+longm[i]-1;
	}
      if(longz[i]!=0)  
	{
	  counter=counter+(longz[i]-16001)-1;
	}
    }
  newDim = fNAdim-counter;   //Dimension of the compressed array
  Short_t* buffer;
  buffer = new Short_t[newDim];
  Int_t counterTwo=0;

  //Fill the buffer of the compressed array
  Int_t g=0;
  Int_t h=0; 
  for(Int_t i=0; i<newDim; i++)
    {
      if(counterTwo<fNAdim)
	{
	  if(fADC[counterTwo]>0)
	    {
	      buffer[i]=fADC[counterTwo];
	    }
	  if(fADC[counterTwo]==-1)
	    {
	      buffer[i]=-(longm[g]);
	      counterTwo=counterTwo+longm[g]-1;
	      g++;
	    }  
	  if(fADC[counterTwo]==0)
	    {
	      buffer[i]=-(longz[h]); 
	      counterTwo=counterTwo+(longz[h]-16001)-1;
	      h++;
	    }  
	  counterTwo++;
	}
    }

  //Copy the buffer
  if(fADC)
    {
      delete [] fADC;
      fADC=0;
    }
  fADC = new Short_t[newDim];
  fNAdim = newDim;
  for(Int_t i=0; i<newDim; i++)
    {
      fADC[i] = buffer[i]; 
    }

  //Delete auxiliary arrays
  if(buffer)
    {
      delete [] buffer;
      buffer=0;
    } 
  if(longz) 
    {
      delete [] longz;
      longz=0;
    }
  if(longm) 
    {
      delete [] longm;
      longm=0;
    }

}

//____________________________________________________________________________________
void AliTRDarrayADC::Expand()
{
  //
  // Expand the array
  //

  //Check if the array has not been already expanded
  Int_t verif=0;
  for(Int_t i=0; i<fNAdim; i++)
    {
      if(fADC[i]<-1)
	{
	  verif++;
	}
    }
  
  if(verif==0)
    {
      //       AliDebug(1,"Nothing to expand");
      return;
    }

  Int_t *longz;
  longz = new Int_t[fNAdim];
  Int_t *longm;
  longm = new Int_t[fNAdim];
  Int_t dimexp=0;
  //Initialize arrays
  memset(longz,0,sizeof(Int_t)*fNAdim);
  memset(longm,0,sizeof(Int_t)*fNAdim);
  Int_t r2=0; 
  Int_t r3=0; 
  for(Int_t i=0; i<fNAdim;i++)
    {
      if((fADC[i]<0)&&(fADC[i]>=-16000))      
	{
	  longm[r2]=-fADC[i];
	  r2++;
	}
      if(fADC[i]<-16000)  
	{
	  longz[r3]=-fADC[i]-16001;  
	  r3++;
	}
    }
  //Calculate the new dimensions of the array
  for(Int_t i=0; i<fNAdim;i++)
    {
      if(longm[i]!=0)       
	{
	  dimexp=dimexp+longm[i]-1;
	}
      if(longz[i]!=0)      
	{
	  dimexp=dimexp+longz[i]-1;
	}
    }
  dimexp=dimexp+fNAdim;   

  //Write in the buffer the new array
  Short_t* bufferE;
  bufferE = new Short_t[dimexp];
  Int_t contaexp =0;     
  Int_t h=0;
  Int_t l=0;  
  for(Int_t i=0; i<dimexp; i++)
    {
      if(fADC[contaexp]>0)  
	{
	  bufferE[i]=fADC[contaexp];
	}

      if((fADC[contaexp]<0)&&(fADC[contaexp]>=-16000))  
	{
	  for(Int_t j=0; j<longm[h];j++)
	    {
	      bufferE[i+j]=-1;
	    }
	  i=i+longm[h]-1;
	  h++;
	}
      if(fADC[contaexp]<-16000)  
	{
	  for(Int_t j=0; j<longz[l];j++)
	    {
	      bufferE[i+j]=0;  
	    }
	  i=i+longz[l]-1;
	  l++;
	}
      contaexp++;
    }
  //Copy the buffer
  if(fADC)
    {
      delete [] fADC;
      fADC=0;
    }

  fADC = new Short_t[dimexp];
  fNAdim = dimexp;
  for(Int_t i=0; i<dimexp; i++)
    {
      fADC[i] = bufferE[i]; 
    }

  //Delete auxiliary arrays
  if(bufferE) delete [] bufferE;
  if(longm) delete [] longm;
  if(longz) delete [] longz;

}
//____________________________________________________________________________________
void AliTRDarrayADC::DeleteNegatives()
{

  //
  //This method modifies the digits array, changing the negative values (-1)
  //Produced during digitization into zero.
  //

  for(Int_t a=0; a<fNAdim; a++)
    {
      if(fADC[a]==-1)
	{
	  fADC[a]=0;
	}
    }
}
//________________________________________________________________________________
void AliTRDarrayADC::Reset()
{
  //
  // Reset the array, the old contents are deleted
  // The array keeps the same dimensions as before
  //
 
  memset(fADC,0,sizeof(Short_t)*fNAdim);
}
//________________________________________________________________________________
void AliTRDarrayADC::ConditionalReset(AliTRDSignalIndex* idx)
{
  //
  // Reset the array, the old contents are deleted
  // The array keeps the same dimensions as before
  //
 
  if(idx->GetNoOfIndexes()>25)
    memset(fADC,0,sizeof(Short_t)*fNAdim);
  else
    {
      Int_t row, col;
      while(idx->NextRCIndex(row, col)){
	for(Int_t time=0; time<fNtime; time++)
	  SetData(row, col, time, 0);
      }
    }

}

//________________________________________________________________________________
void AliTRDarrayADC::CreateLut()
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
