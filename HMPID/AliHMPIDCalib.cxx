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
#include "AliHMPIDParam.h" //class header
#include "AliHMPIDRawStream.h" //class header
#include <fstream>
#include <TTree.h>



ClassImp(AliHMPIDCalib) 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::AliHMPIDCalib():
faddl(0x0),
fPadAdc(0x0),
fIsPad(0x0),
fFile(0x0),
fLdcId(0),
fTimeStamp(0),
fRunNum(0),
fSigCut(0),
fWritePads(0)
{
  //
  //constructor
  //
  faddl = new Bool_t[AliHMPIDRawStream::kNDDL];
  Int_t nPads =  (AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPcx+1)*(AliHMPIDParam::kMaxPcy+1);
  fPadAdc=new TH1I*[nPads];  
  fIsPad=new Bool_t[nPads];  
  for(Int_t np=0;np<nPads;np++) {fPadAdc[np]=0x0;   fIsPad[np]=kFALSE;}
  fWritePads=kFALSE;
  Init();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::~AliHMPIDCalib()
{
  //
  //destructor
  //
  if (faddl)     { delete [] faddl;   faddl = 0x0;   } 
  if (fPadAdc)   {delete  [] fPadAdc; fPadAdc=0x0;   }  
  if (fIsPad)   {delete  [] fIsPad; fIsPad=0x0;   }  
  if (fFile)     {delete     fFile;   fFile=0x0;     }  
  fLdcId=0;
  fTimeStamp=0;
  fRunNum=0;
  fSigCut=0;
  fWritePads=0;
}//dtor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::Init()
{
  //
  //Init the q calc.
  //Arguments: none
  //Return: none
  //
    fSigCut=3;  
  
    for(Int_t iDDL=0; iDDL< AliHMPIDRawStream::kNDDL; iDDL++) 
      {
         for(Int_t ierr=0; ierr <AliHMPIDRawStream::kSumErr ; ierr++) {
            fErr[iDDL][ierr]=0;
            }
        
        faddl[iDDL]=kFALSE;
           for(Int_t row = 1; row <=AliHMPIDRawStream::kNRows; row++){
              for(Int_t dil = 1; dil <=AliHMPIDRawStream::kNDILOGICAdd; dil++){
                for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
                       fsq[iDDL][row][dil][pad]=0;
                      fsq2[iDDL][row][dil][pad]=0;                  
                }//pad
            }//dil
          }//row
        }//DDL
}//Init()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::SetRunParams(ULong_t runNum,Int_t timeStamp, Int_t ldcId)
{
  //  
  //Set run parameters for the Pedestal and Error Files
  //Arguments: run number, time stamp and LDC Id
  //Returns: none
  //
  fRunNum=(Int_t)runNum;
  fTimeStamp=timeStamp;
  fLdcId=ldcId;
}//SetRunParams()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::SetSigCutFromFile(Char_t* name)
{
  //
  //Set Sigma Cut from the file on the LDC, if the input file is not present default value is set!
  //Arguments: the name of the SigmaCut file on the LDC
  //Returns: none
  //
  Int_t nSigCut=0;
  ifstream infile(name);
  if(!infile.is_open()) {fSigCut=3; return;}
  while(!infile.eof())
    {
    infile>>nSigCut;
  }
  infile.close();
  fSigCut=nSigCut; 
}//SetSigCutFromFile()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
void AliHMPIDCalib::InitHisto(Int_t q,Int_t histocnt,Char_t* name)
{
  //
  //Init the pad histos. For one DDL we have 11520 pads. ONLY if ENABLED!
  //Arguments: q-charge, the absolute number of the histogram (AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPcx+1)*(AliHMPIDParam::kMaxPcy+1) and the name of the histogram (unique) 
  //Returns: none
  //
  
 fFile->cd();
 Double_t lowbin,highbin=0;
// Printf("InitHisto: histocnt: %d Name: %s",histocnt,name);
 if(fIsPad[histocnt]==kTRUE) return;
 
 lowbin=q-40.5; highbin=q+40.5;  
 fPadAdc[histocnt]=new TH1I(name,name,81,lowbin,highbin);
 fPadAdc[histocnt]->Sumw2();
 fIsPad[histocnt]=kTRUE;
 
}//InitHisto()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillHisto(Int_t histocnt,Int_t q)
{
  //
  //Fill the ADC histograms for each pad
  //Arguments:  q-charge, the absolute number of the histogram (AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPcx+1)*(AliHMPIDParam::kMaxPcy+1)
  //Returns: none
  //
  fFile->cd();
  if(fIsPad[histocnt]==kFALSE) return;
  fPadAdc[histocnt]->Fill(q);
 
}//InitHisto()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::InitFile(Int_t ldcId)
{
  //
  //Initialize the ADC histo output file (one per LDC)
  //Arguments: LDC Id
  //Returns: none
  //
  fFile=new TFile(Form("HmpidPadsOnLdc%2d.root",ldcId),"RECREATE");
}//InitFile()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::CloseFile(Int_t /*ldcId*/)
{
  //
  //Close the ADC histo output file (one per LDC)
  //Arguments: LDC Id
  //Returns: none
  //
  fFile->cd();
  Int_t nPads = (AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPcx+1)*(AliHMPIDParam::kMaxPcy+1);
  for(Int_t np=0;np<nPads;np++) {if(fIsPad[np]==kTRUE) fPadAdc[np]->Write();} 
  fFile->Close();
}//CloseFile()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillPedestal(Int_t abspad,Int_t q)
{
  //
  //Called from the HMPIDda and fills the pedestal values
  //Arguments: absulote pad number as from AliHMPIDParam and q-charge
  //Returns: none
  //
  Int_t nDDL=0, row=0, dil=0, adr=0;
  //The decoding (abs. pad -> ddl,dil,...) is the same as in AliHMPIDDigit::Raw
  Int_t y2a[6]={5,3,1,0,2,4};

       nDDL=  2*AliHMPIDParam::A2C(abspad)+AliHMPIDParam::A2P(abspad)%2;              //DDL# 0..13
  Int_t tmp=  1+AliHMPIDParam::A2P(abspad)/2*8+AliHMPIDParam::A2Y(abspad)/6;          //temp variable
        row=   (AliHMPIDParam::A2P(abspad)%2)? 25-tmp:tmp;                            //row r=1..24
        dil=  1+AliHMPIDParam::A2X(abspad)/8;                                         //DILOGIC 
        adr=y2a[AliHMPIDParam::A2Y(abspad)%6]+6*(AliHMPIDParam::A2X(abspad)%8);       //ADDRESS 0..47 
  //........... decoding done      
        
                                   
    if(q>0) { 
         fsq[nDDL][row][dil][adr]+=q;
        fsq2[nDDL][row][dil][adr]+=(q*q);
                       faddl[nDDL]=kTRUE;  
        }

     Int_t histocnt=0;   
           histocnt=(nDDL)*11520+(row-1)*480+(dil-1)*48+adr;      //Histo counter for a single DDL  
     if(fWritePads==kTRUE)
     {
       InitHisto(q,histocnt,Form("hPad_Ch_%d_Pc_%d_Px_%d_Py_%d",
                                   AliHMPIDParam::A2C(abspad),AliHMPIDParam::A2P(abspad),AliHMPIDParam::A2X(abspad),AliHMPIDParam::A2Y(abspad))); 
       FillHisto(histocnt,q);
      }
            
}//FillPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillErrors(Int_t nDDL,Int_t eType, Int_t nErr)
{
  //
  //Fill decoding errors from AliHMPIDRawStream
  //Arguments: nDDL-DDL number, eType- error type as in AliHMPIDRawStream.h and the # of occurence for eType
  //Retutns: none
  //
    if(nErr<=0) return;
    if(eType < 0 || eType> AliHMPIDRawStream::kSumErr ) return;
    fErr[nDDL][eType]=fErr[nDDL][eType]+nErr;
  
}//FillErrors()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::WriteErrors(Int_t nDDL, Char_t* name, Int_t nEv)
{
  //
  //Write decoding errors to a txt file
  //Arguments: nDDL-DDL number, name of the error file and number of the read events
  //Retutns: kTRUE/kFALSE
  //
  
  if(faddl[nDDL]==kFALSE) return kFALSE;                                    //if ddl is missing no error file is created
  ofstream outerr;  outerr.open(name);                                      //open error file
  outerr << Form("%6s %2d\n","RunNum",(Int_t)fRunNum);                      //read run number
  outerr << Form("%6s %2d\n","LdcId" ,       fLdcId);                       //read LDC Id
  outerr << Form("%6s %2d\n","TimeSt",       fTimeStamp);                   //read time stamp
  outerr << Form("%6s %2d\n","NumEvt",       nEv);                          //read number of events processed

  for(Int_t  ierr=0; ierr <AliHMPIDRawStream::kSumErr; ierr++) outerr << Form("%2d\t",fErr[nDDL][ierr]); //write errors
                                                               outerr << Form("\n");                     //last break
  outerr.close();                                                                                        //write error file
  return kTRUE;
    
}//FillErrors()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::CalcPedestal(Int_t nDDL, Char_t* name, Int_t nEv)    
{
  //
  //Calculate pedestal for each pad  
  //Arguments: nDDL-DDL number, name of the pedestal file and number of the read events
  //Retutns: kTRUE/kFALSE
  //
  
   
  if(faddl[nDDL]==kFALSE) return kFALSE;                   //if ddl is missing no ped file is created (and also for LDC selection). Check with Paolo what he checks for?!  
  Float_t mean=0,sigma=0;
  Float_t qs2m=0,qsm2=0;
  ofstream out;                                           //to write the pedestal text files
  Int_t inhard;
  out.open(name);
  //out << Form("%2d %2d\n",(Int_t)runNum,nEv);
  out << Form("%6s %2d\n","RunNum",(Int_t)fRunNum);
  out << Form("%6s %2d\n","LdcId" ,       fLdcId);
  out << Form("%6s %2d\n","TimeSt",       fTimeStamp);
  out << Form("%6s %2d\n","NumEvt",       nEv);
  out << Form("%6s %2d\n","SigCut",       fSigCut);
  for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){
    for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){
      for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
        
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
