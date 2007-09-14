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
AliHMPIDCalib::AliHMPIDCalib():
   fPedTree(0x0),
   fa(-1),
   fd(-1),
   fr(-1),
   fq(-1),
   fPedHisto(0) 
{
  Init();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::~AliHMPIDCalib()
{
  //destructor
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::Init()
{
  //..
  for(Int_t iddl=0;iddl<11;iddl++) faddl[iddl]=kFALSE;
  //
  // Called to initialize the TTree in which the raw data words will be stored 
  //
  fPedHisto=new TH1F("fPedHisto","Temporary pedestal",4096,-0.5,4095.5);                                //init pedestal histo
  fPedTree = new TTree("PedTree","HMPID Pedestal Tree");                                                //init pedestal tree
  fPedTree->Branch("diladd",&fa,"diladd/I");
  fPedTree->Branch("dilnum",&fd,"dilnum/I");
  fPedTree->Branch("dilrow",&fr,"dilrow/I");
  fPedTree->Branch("qdc"   ,&fq,"qdc/F   ");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillPedestal(Int_t ddl,Int_t row, Int_t dil,Int_t adr,Int_t q)
{
  //
  //Called from the HMPIDda and fills the pedestal tree 
  //
  fq=q;
  fr=row;
  fd=dil;
  fa=adr;
  faddl[ddl]=kTRUE; 
  if(fq>-1) fPedTree->Fill();    
}//FillPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::CalcPedestal(Int_t ddl, Char_t* name)    
{
  //
  //Calculate pedestal for each pad  
  //
  
  ofstream out;                                           //to write the pedestal text files
  Double_t mean,sigma;
  Int_t inhard;
  if(faddl[ddl]==kFALSE) return kFALSE;                   //if ddl is missing no ped file is created (and also for LDC selection). Check with Paolo what he checks for?!
    out.open(name);
    for(Int_t row=1; row < 25; row++){
     for(Int_t dil=1; dil < 11; dil++){
      for(Int_t adr=0; adr < 48; adr++){
          fPedHisto->Reset();
          fPedTree->Draw("fq>>fPedHisto",Form("fr==%d&&fd==%d&&fa==%d",row,dil,adr));
          mean=fPedHisto->GetMean();
          sigma=fPedHisto->GetRMS();
          inhard=((Int_t(mean))<<9)+Int_t(mean+3*sigma);
	  out << Form("%2i %2i %2i %5.2f %5.2f %x\n",row,dil,adr,mean,sigma,inhard);
        }//adr
      }//dil
    }//row
    out.close();                                          //write pedestal file
  return kTRUE;
}//CaclPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
