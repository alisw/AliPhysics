/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id:$*/

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode the SDD Raw Data from the CarlosRX format to  //
// a compressed format consisting in a word of 32 bit per cell   //
// The 32 bits for a data word are defined as follows:           //
//   31 control bit (0=data word, 1= control word)               //
//   30 -                                                        //
//   29  |                                                       //
//   28  |-> 4 bits to identify the Carlos (0-11) inside the DDL //
//   27 -                                                        //
//   26 detecot side (0= left, =right)                           //
//   25 -                                                        //
//   24  |                                                       //
//   23  |                                                       //
//   22  |                                                       //
//   21  |-> 8 bits to identify the anode number (0-255)         //
//   20  |                                                       //
//   19  |                                                       //
//   18 -                                                        //
//   17 -                                                        //
//   16  |                                                       //
//   15  |                                                       //
//   14  |                                                       //
//   13  |-> 8 bits to identify the time bin (0-255)             //
//   12  |                                                       //
//   11  |                                                       //
//   10 -                                                        //
//    9 -                                                        //
//    8  |                                                       //
//    7  |                                                       //
//    6  |                                                       //
//    5  |                                                       //
//    4  |-> 10 bit for the ADC counts                           //
//    3  |                                                       //
//    2  |                                                       //
//    1  |                                                       //
//    0 -                                                        //
//                                                               //
// Plus 2 types of control words:                                //
// - DDL identifier with the 4 more significant bits     = 1000  //
// - End of module data (needed by the Cluster Finder)   = 1111  //
//                                                               //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSCompressRawDataSDD.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"


ClassImp(AliITSCompressRawDataSDD)

AliITSCompressRawDataSDD::AliITSCompressRawDataSDD(TString filename):
TObject(),
fEventRange(kFALSE),
fFirstEvent(0),
fLastEvent(0)
{
  fNameFile=filename;
}
//______________________________________________________________________
void AliITSCompressRawDataSDD::Compress(){
  Int_t iev=0;
  if(fEventRange) iev=fFirstEvent;
  AliRawReader *rd;   
  if(fNameFile.Contains(".root")){
    rd=new AliRawReaderRoot(fNameFile.Data(),iev);
  }else{
    rd=new AliRawReaderDate(fNameFile.Data(),iev);
  }

  FILE *outtxt=fopen("data.txt","w");
  Int_t oldddl=-1;
  UInt_t word=0;
  do{
    rd->Reset();

    AliITSRawStreamSDD s(rd);
    while(s.Next()){
      if(rd->GetDDLID()!=oldddl){
	word=8<<28;
	word+=rd->GetDDLID();
	fprintf(outtxt,"%08X\n",word);
	oldddl=rd->GetDDLID();
      }
      if(s.IsCompletedModule()==kFALSE){
	word=s.GetCarlosId()<<27;
	word+=s.GetChannel()<<26;
	word+=s.GetCoord1()<<18;
	word+=s.GetCoord2()<<10;
	word+=s.GetSignal();
	fprintf(outtxt,"%08X\n",word);
      }
      if(s.IsCompletedModule()==kTRUE){
	word=15<<28;
	word+=s.GetCarlosId();
	fprintf(outtxt,"%08X\n",word);
      }
    }
    iev++;
    if(fEventRange && iev>fLastEvent) break;
  }while(rd->NextEvent());
  fclose(outtxt);
}
