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
#include "TTreePlayer.h"
#include <fstream>
#include <TTree.h>
//#include "AliHMPIDDigit.h"

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
void AliHMPIDCalib::FillPedestal(Int_t nDDL,Int_t row, Int_t dil,Int_t adr,Int_t q)
{
  //
  //Called from the HMPIDda and fills the pedestal tree 
  //
    if(q>0) { 
         fsq[nDDL][row][dil][adr]+=q;
        fsq2[nDDL][row][dil][adr]+=(q*q);
                       faddl[nDDL]=kTRUE;  
        }

}//FillPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::CalcPedestal(Int_t nDDL, Char_t* name, Int_t nEv)    
{
  //
  //Calculate pedestal for each pad  
  //
  Float_t mean=0,sigma=0;
  Float_t qs2m=0,qsm2=0;
  ofstream out;                                           //to write the pedestal text files
  Int_t inhard;
  if(faddl[nDDL]==kFALSE) return kFALSE;                   //if ddl is missing no ped file is created (and also for LDC selection). Check with Paolo what he checks for?!
  out.open(name);
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
