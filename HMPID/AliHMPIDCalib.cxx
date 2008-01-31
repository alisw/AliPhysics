//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliHMPIDCalib.h" //class header
#include <fstream>
#include <TTree.h>



ClassImp(AliHMPIDCalib) 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::AliHMPIDCalib()
{
  //
  //constructor
  //
  Init();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::~AliHMPIDCalib()
{
  //
  //destructor
  //
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::Init()
{
  //
  //Init the q calc.
  //
   
    for(Int_t iDDL=0; iDDL< AliHMPIDCalib::kNDDL; iDDL++) 
      {
        faddl[iDDL]=kFALSE;

         for(Int_t ierr=0; ierr < AliHMPIDRawStream::kSumErr; ierr++)  fNumOfErr[iDDL][ierr]=0;               //reset errors for all DDLs

         for(Int_t row = 1; row <=AliHMPIDCalib::kNRows; row++){
          for(Int_t dil = 1; dil <=AliHMPIDCalib::kNDILOGICAdd; dil++){
            for(Int_t pad = 0; pad < AliHMPIDCalib::kNPadAdd; pad++){
                     fsq[iDDL][row][dil][pad]=0;
                    fsq2[iDDL][row][dil][pad]=0;                  
               }
          }
        }
      }
}//Init()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//void AliHMPIDCalib::FillPedestal(Int_t nDDL,Int_t row, Int_t dil,Int_t adr,Int_t q)
void AliHMPIDCalib::FillPedestal(Int_t pad,Int_t q)
{
  //
  //Called from the HMPIDda and fills the pedestal tree 
  //
  Int_t nDDL=0, row=0, dil=0, adr=0;
  //The decoding (abs. pad -> ddl,dil,...) is the same as in AliHMPIDDigit::Raw
  Int_t y2a[6]={5,3,1,0,2,4};

       nDDL=  2*AliHMPIDParam::A2C(pad)+AliHMPIDParam::A2P(pad)%2;             //DDL# 0..13
  Int_t tmp=  1+AliHMPIDParam::A2P(pad)/2*8+AliHMPIDParam::A2Y(pad)/6;         //temp variable
        row=   (AliHMPIDParam::A2P(pad)%2)? 25-tmp:tmp;                         //row r=1..24
        dil=  1+AliHMPIDParam::A2X(pad)/8;                                     //DILOGIC 
        adr=y2a[AliHMPIDParam::A2Y(pad)%6]+6*(AliHMPIDParam::A2X(pad)%8);    //ADDRESS 0..47 
        
  //  Printf("AbsPadNum: %d nDDL: %d tmp %d row: %d dil: %d adr: %d",pad,nDDL,tmp,row,dil,adr);
  //........... decoding done      
                                      
    if(q>0) { 
         fsq[nDDL][row][dil][adr]+=q;
        fsq2[nDDL][row][dil][adr]+=(q*q);
                       faddl[nDDL]=kTRUE;  
        }

        
        
        
}//FillPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillErrors(Int_t nDDL,Int_t nErrType, Int_t nErr)
{
  //
  //Fill decoding errors
  //
 // fhDdlDecodErrors[nDDL]->Fill(nErrType,nErr);      //select DDL, select bin and add the occurence
    fNumOfErr[nDDL][nErrType]+=nErr;
  
}//FillErrors()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::WriteErrors(ULong_t runNum, Int_t nDDL, Char_t* name, Int_t /*nEv*/)
{
  //
  //Write decoding errors to a txt file
  //
  
  if(faddl[nDDL]==kFALSE) return kFALSE;                   //if ddl is missing no error file is created
  ofstream outerr;  
  outerr.open(name);                                      //open error file
                                                             
                                                                outerr << Form("%d \n",(Int_t)runNum);
  for(Int_t  ierr=1; ierr < AliHMPIDRawStream::kSumErr; ierr++) outerr << Form("%2d\t",fNumOfErr[nDDL][ierr]);
                                                                outerr << Form("\n");  
  outerr.close();                                          //write error file
  return kTRUE;
    
}//FillErrors()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::CalcPedestal(ULong_t runNum, Int_t nDDL, Char_t* name, Int_t nEv)    
{
  //
  //Calculate pedestal for each pad  
  //
   
  if(faddl[nDDL]==kFALSE) return kFALSE;                   //if ddl is missing no ped file is created (and also for LDC selection). Check with Paolo what he checks for?!  
  Float_t mean=0,sigma=0;
  Float_t qs2m=0,qsm2=0;
  ofstream out;                                           //to write the pedestal text files
  Int_t inhard;
  //Printf("From AliHMPIDCalib::CalcPedestal: faddl[%i]= %i",nDDL,faddl[nDDL]);

  out.open(name);
  out << Form("%d \n",(Int_t)runNum);
  for(Int_t row = 1; row <= AliHMPIDCalib::kNRows; row++){
    for(Int_t dil = 1; dil <= AliHMPIDCalib::kNDILOGICAdd; dil++){
      for(Int_t pad = 0; pad < AliHMPIDCalib::kNPadAdd; pad++){
        
        mean = fsq[nDDL][row][dil][pad]/1.0/nEv;
        
        qs2m = fsq2[nDDL][row][dil][pad]/1.0/nEv;
        qsm2 = TMath::Power(fsq[nDDL][row][dil][pad]/1.0/nEv,2);
        sigma= TMath::Sqrt(qs2m-qsm2);
        
        inhard=((Int_t(mean))<<9)+Int_t(mean+3*sigma);
        out << Form("%2i %2i %2i %5.2f %5.2f %x\n",row,dil,pad,mean,sigma,inhard);
        }//adr
      }//dil
    }//row
    out.close();                                          //write pedestal file
  return kTRUE;
}//CaclPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
