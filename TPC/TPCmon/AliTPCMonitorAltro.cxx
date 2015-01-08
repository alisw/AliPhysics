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
$Log$
Revision 1.2  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorAltro class
////
//// Class for decoding raw TPC data in the ALTRO format.
//// Data are transformed from 32 bit words to 10 bit or 40 bit respectively.
//// The whole payload is transformed at once and written to an array. 
//// The single channels are decoded starting from the trailer word which is
//// decoded by DecodeTrailer(Int_t pos)
////
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
/////////////////////////////////////////////////////////////////////////


#include "AliTPCMonitorAltro.h"
#include "AliLog.h" 
#include <Riostream.h>
using namespace std;
ClassImp(AliTPCMonitorAltro)  

//_____________________________________________________________________________________________
AliTPCMonitorAltro::AliTPCMonitorAltro(UInt_t* memory, Int_t size, Int_t fformat) :
  fverb(0),
  fmemory(memory),
  fsize(size),
  f40BitArray(0),
  f10BitArray(0),
  fdecoderPos(0),
  fallocate40BitArray(false),
  fallocate10BitArray(false),
  foffset(0),  
  fwrite10bit(0),
  fTrailerNWords(0), 
  fTrailerHwAddress(0),
  fTrailerDataPos(0),
  fTrailerBlockPos(0),
  fTrailerPos(0),
  fNextPos(0),
  ffilename()
{
  // Constructor: Set different CDH offsets for ROOT (0) and date format
  // Data offset can be changed via SetDataOffset(Int_t val)

  if(     fformat==0) foffset=7; // old CHD Format
  else if(fformat==1) foffset=8; // memory pointer from DATE (start CDH)
  else if(fformat==2) foffset=0; // memory pointer from ROOT (after CDH)
} 


//_____________________________________________________________________________________________
AliTPCMonitorAltro::AliTPCMonitorAltro(const AliTPCMonitorAltro &altro) :
  TNamed(altro),
  fverb(altro.fverb),
  fmemory(altro.fmemory),
  fsize(altro.fsize),
  f40BitArray(altro.f40BitArray),
  f10BitArray(altro.f10BitArray),
  fdecoderPos(altro.fdecoderPos),
  fallocate40BitArray(altro.fallocate40BitArray),
  fallocate10BitArray(altro.fallocate10BitArray),
  foffset(altro.foffset),
  fwrite10bit(altro.fwrite10bit),
  fTrailerNWords(altro.fTrailerNWords),
  fTrailerHwAddress(altro.fTrailerHwAddress),
  fTrailerDataPos(altro.fTrailerDataPos),
  fTrailerBlockPos(altro.fTrailerBlockPos),
  fTrailerPos(altro.fTrailerPos),
  fNextPos(altro.fNextPos),
  ffilename(altro.ffilename)
{
  // copy constructor

}

//_____________________________________________________________________________________________
AliTPCMonitorAltro &AliTPCMonitorAltro::operator =(const AliTPCMonitorAltro& altro)
{
  // assignement operator 
  
  if(this!=&altro){ 
    ((TNamed *)this)->operator=(altro);
    fverb=altro.fverb;
    fmemory=altro.fmemory;
    fsize=altro.fsize;
    f40BitArray=altro.f40BitArray;
    f10BitArray=altro.f10BitArray;
    fdecoderPos=altro.fdecoderPos;
    fallocate40BitArray=altro.fallocate40BitArray;
    fallocate10BitArray=altro.fallocate10BitArray;
    foffset=altro.foffset;
    fwrite10bit=altro.fwrite10bit;
    fTrailerNWords=altro.fTrailerNWords;
    fTrailerHwAddress=altro.fTrailerHwAddress;
    fTrailerDataPos=altro.fTrailerDataPos;
    fTrailerBlockPos=altro.fTrailerBlockPos;
    fTrailerPos=altro.fTrailerPos; 
    fNextPos=altro.fNextPos;
    ffilename = altro.ffilename;
  }
  return *this;
}


//_____________________________________________________________________________________________
AliTPCMonitorAltro::~AliTPCMonitorAltro() {
  // Destructor
  if(fallocate40BitArray == true) delete[] f40BitArray;
  if(fallocate10BitArray == true) delete[] f10BitArray;
 }

//_____________________________________________________________________________________________
void AliTPCMonitorAltro::Allocate40BitArray() 
{
  // create array for 40 bit decoded data
  fallocate40BitArray = true;
  f40BitArray = new long long[Get40BitArraySize()];
}

//_____________________________________________________________________________________________
void AliTPCMonitorAltro::Allocate10BitArray() 
{
  // Create array for 10 bit decoded data
  fallocate10BitArray = true;
  f10BitArray = new Short_t[Get10BitArraySize()];
} 

//_____________________________________________________________________________________________
long long *AliTPCMonitorAltro::Get40BitArray() 
{
  // Return pointer to array for 40 bit decoded data 
  return f40BitArray;
}

//_____________________________________________________________________________________________
Short_t *AliTPCMonitorAltro::Get10BitArray() 
{
  // Return pointer to array for 10 bit decoded data 
  return f10BitArray;
}

// //_____________________________________________________________________________________________
// Int_t AliTPCMonitorAltro::Get40BitArraySize() 
// {
//   // Return number of 40 bit words in payload 
//   return fmemory[fsize-1];
// }

// //_____________________________________________________________________________________________
// Int_t AliTPCMonitorAltro::Get10BitArraySize() 
// { 
//   // Return number of 10 bit words in payload 
//   return fmemory[fsize-1]*4;
// }

//_____________________________________________________________________________________________
void AliTPCMonitorAltro::Decodeto40Bit() 
{
  // Decode 32 bit words in fmemory to 40 bit words
  Long64_t blackbox = 0;
  Int_t rest;
  
  for(Int_t i = 0; i < Get40BitArraySize(); i++) 
    {
      rest = i%4;
      switch(rest) {
      case 0:
	blackbox = (fmemory[foffset])     + (((Long64_t)(fmemory[foffset+1]&fgk08BitOn))<<32);
	foffset +=1;
	break;
      case 1:
	blackbox = (fmemory[foffset]>>8 ) + (((Long64_t)(fmemory[foffset+1]&fgk16BitOn))<<24);
	foffset +=1;
	break;
      case 2:
	blackbox = (fmemory[foffset]>>16) + (((Long64_t)(fmemory[foffset+1]&fgk24BitOn))<<16);
	foffset +=1;
	break;
      case 3:
	blackbox = (fmemory[foffset]>>24) + (((Long64_t)(fmemory[foffset+1]         ))<< 8);
	foffset +=2;
	break;
      default:
	blackbox = 0;
	break;
      }
      f40BitArray[i] = blackbox;
    }
}

//_____________________________________________________________________________________________
void AliTPCMonitorAltro::Decodeto10Bit(Int_t equipment) 
{
  // Decode 32 bit words in fmemory to 10 bit words.
  // Write words to file if fwrite10bit ("Write 10 bit" in Monitor.C Gui ) is set

  Long64_t blackbox = 0;
  Int_t    rest     = 0;
  
  ofstream datout;
  if(fwrite10bit)
  {
    TString nameout=Form("%s_PayloadEquipmentId_%03i.txt",ffilename.Data(),equipment);
    datout.open(nameout.Data());
    AliInfo(Form("AliTPCMonitorAltro::decodeto10Bit :  Write Data to %s",nameout.Data()));
    if(foffset==0) datout <<  "Payload without CDH in 10bit words " << endl;
    else           datout  << "CDH in 32 bit hex words words (" << foffset << " lines )  followed by payload in 10 bit words " << endl;
    for(Int_t ih = 0; ih< foffset ; ih++)
    {
      datout << hex << fmemory[ih] << endl;
    }
  }
  
  for(Int_t ind = 0; ind < Get40BitArraySize(); ind++)
  {
    rest = ind%4;
    switch(rest)
    {
    case 0:
      blackbox = (fmemory[foffset])     + (((Long64_t)(fmemory[foffset+1]&fgk08BitOn))<<32);
      foffset +=1;
      break;
    case 1:
      blackbox = (fmemory[foffset]>>8 ) + (((Long64_t)(fmemory[foffset+1]&fgk16BitOn))<<24);
      foffset +=1;
      break;
    case 2:
      blackbox = (fmemory[foffset]>>16) + (((Long64_t)(fmemory[foffset+1]&fgk24BitOn))<<16);
      foffset +=1;
      break;
    case 3:
      blackbox = (fmemory[foffset]>>24) + (((Long64_t)(fmemory[foffset+1]         ))<< 8);
      foffset +=2;
      break;
    default:
      blackbox = 0;
      break;
    }
    f10BitArray[ind*4+0] = (Short_t)( blackbox & fgkmask10 )    ;
    f10BitArray[ind*4+1] = (Short_t)((blackbox & fgkmask20)>>10);
    f10BitArray[ind*4+2] = (Short_t)((blackbox & fgkmask30)>>20);
    f10BitArray[ind*4+3] = (Short_t)((blackbox & fgkmask40)>>30);
  }
  if(fwrite10bit)
  {
    for(Int_t ind = 0; ind < Get40BitArraySize(); ind++)
    {
      datout << dec <<  f10BitArray[ind*4+0] << "\n";
      datout << dec <<  f10BitArray[ind*4+1] << "\n";
      datout << dec <<  f10BitArray[ind*4+2] << "\n";
      datout << dec <<  f10BitArray[ind*4+3] << "\n";
    }
  }
  
  if(fwrite10bit) datout.close();
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorAltro::DecodeTrailer(Int_t pos)
{
  // Decode the trailer word starting at position pos in fmemory
  // Check if information leads to proper next trailer position
  if (GetAltroVersion()==0xaabb) return DecodeTrailerVbb(pos);
  fTrailerPos = pos;
  if(pos<=4) return 0;
  
  Long64_t  words         = 0;
  Long64_t  tail          = 0;
  Long64_t  trailer       = 0;
  Long64_t  carry         = 0;
  Long64_t  rest          = 0 ;
    
  for(Long64_t  iter = 0 ; iter<4;iter++)
    {
      carry =  f10BitArray[pos-iter] ; 
      carry = ( carry << (30- ((iter)*10)));
      trailer += carry ;
    }
  
  fTrailerHwAddress = (trailer & ((Long64_t )fgkTrailerMaskHardw)  );
  words             = (Long64_t )( (trailer & ((Long64_t )fgkTrailerMaskNWords))>>16);
  tail           = (Long64_t )( (trailer & ((Long64_t )fgkTrailerMaskTail   ))>>26 );
  
  if(words%4!=0) rest =  4-(words%4);
  fTrailerNWords      = words+rest ;
  fTrailerDataPos     = pos -4 -rest ;
  fTrailerBlockPos    = pos -4 ;
  fNextPos            = (pos -fTrailerNWords -4);
  
  if(       tail!=fgkTrailerTail        ) { AliError(Form("Could not read Trailer. \"Write 10bit\" for this event. Last Trailer line (2AA): %i. Supp.next Trailer line (2AA): %i ",pos,fNextPos));    return -1;      }
  else if(     fNextPos==-1           ) { /* was last channel  */   	                                                                                                                            return  0;      }
  else if(     fNextPos <0            ) { AliError("Next Trailer position < 0 ");                                                                                                                   return -1;      }
  else if((f10BitArray[fNextPos]!=682)) { AliError(Form("Could not find tail (2AA) at next supposed position %i",fNextPos));                                                                        return -1;      }
  else                                  {                                                                                                                                                           return fNextPos;}
  
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorAltro::DecodeTrailerVbb(Int_t pos)
{
  // Decode the trailer word starting at position pos in fmemory
  // Check if information leads to proper next trailer position
    fTrailerPos = pos;
  if(pos>Get10BitArraySize()-4) return 0;
  
  Long64_t  words         = 0;
  Long64_t  tail          = 0;
  Long64_t  trailer       = 0;
  Long64_t  carry         = 0;
  Long64_t  rest          = 0 ;
  
  for(Long64_t  iter = 0 ; iter<4;iter++)
  {
    carry =  f10BitArray[pos-iter] ;
    carry = ( carry << (30- ((iter)*10)));
    trailer += carry ;
  }
  
  fTrailerHwAddress = (trailer & ((Long64_t )fgkTrailerMaskHardw)  );
  words             = (Long64_t )( (trailer & ((Long64_t )fgkTrailerMaskNWords))>>16);
  tail           = (Long64_t )( (trailer & ((Long64_t )fgkTrailerMaskTail   ))>>26 );

    //Decide on read direction (sign) from the RCU trailer information
  Int_t sign=-1;
  if (GetAltroVersion()==0xaabb) sign=1;
  
  if(words%4!=0) rest =  4-(words%4);
  fTrailerNWords      = words+rest ;
  fTrailerDataPos     = pos +sign*(4+rest) ;
  fTrailerBlockPos    = pos +sign*4 ;
  fNextPos            = (pos +sign*(fTrailerNWords+4));
  
  if(       tail==fgkTrailerTailErr        ) {
    AliError(Form("Found Altro header with error marker (2AEE)[%0llx]: %i. Supp.next Trailer line (2AA)[%0x]: %i ",
                  tail,pos,f10BitArray[fNextPos],fNextPos));
    return -fNextPos;
  } else if ( tail!=fgkTrailerTail        ) {
    AliError(Form("Could not read Trailer. \"Write 10bit\" for this event. Last Trailer line (2AA)[%0llx]: %i. Supp.next Trailer line (2AA)[%0x]: %i ",
                  tail,pos,f10BitArray[fNextPos],fNextPos));    return -1;
  } else if(     fNextPos==Get10BitArraySize()+3           ) { /* was last channel  */
    return  0;
  }  else if(     fNextPos >=Get10BitArraySize()            ) {
    AliError(Form("Next Trailer position < 0 %i", fNextPos));
    return -1;
  } else if((f10BitArray[fNextPos]!=682)&&(f10BitArray[fNextPos]!=686)) {
    AliError(Form("Could not find tail (2AA) at next supposed position %i",fNextPos));
    return -1;
  } else {
    return fNextPos;
  }
  
}


