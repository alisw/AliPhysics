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
#include "AliHMPIDDigit.h" //class header
#include <fstream>
#include <TTree.h>
#include <TSystem.h>
#include <TTimeStamp.h>





ClassImp(AliHMPIDCalib) 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::AliHMPIDCalib():
faddl(0x0),
fsq(0x0),
fsq2(0x0),
fnpc(0x0),
fpedQ0(0x0),
fErr(0x0),
fPadAdc(0x0),
fIsPad(0x0),
fFile(0x0),
fLdcId(0),
fTimeStamp(0),
fRunNum(0),
fSigCut(0),
fnDDLInStream(0x0),
fnDDLOutStream(0x0),
fLargeHisto(kFALSE),
fSelectDDL(0),
fDeadMap(0x0),
fPedMeanMap(0x0),
fPedSigMap(0x0),
f1DPedMean(0x0),
f1DPedSigma(0x0),
fNumMaskedPads(0),
fNumDeadPads(0)
{
  //
  //constructor
  //
  
  faddl = new Bool_t[AliHMPIDRawStream::kNDDL];
  Int_t nPads =  (AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPcx+1)*(AliHMPIDParam::kMaxPcy+1);
 
  fpedQ0 = new Int_t***[AliHMPIDRawStream::kNDDL+1];
  fsq2   = new Float_t ***[AliHMPIDRawStream::kNDDL+1];
  fsq    = new Float_t ***[AliHMPIDRawStream::kNDDL+1];
  fnpc   = new Int_t ***[AliHMPIDRawStream::kNDDL+1];
    fErr = new Int_t*[AliHMPIDRawStream::kNDDL+1];
   
  fnDDLInStream  = new Int_t[AliHMPIDRawStream::kNDDL+1];
  fnDDLOutStream = new Int_t[AliHMPIDRawStream::kNDDL+1];

  
  for(Int_t iDDL=0;iDDL<AliHMPIDRawStream::kNDDL+1;iDDL++) {
    
      fErr[iDDL] = new Int_t[AliHMPIDRawStream::kSumErr+1];
    fpedQ0[iDDL] = new Int_t**[AliHMPIDRawStream::kNRows+1];
       fsq[iDDL] = new Float_t**[AliHMPIDRawStream::kNRows+1];
      fsq2[iDDL] = new Float_t**[AliHMPIDRawStream::kNRows+1];
      fnpc[iDDL] = new Int_t**[AliHMPIDRawStream::kNRows+1];
      
      for(Int_t iRow=0;iRow<AliHMPIDRawStream::kNRows+1;iRow++)  {
      
       fpedQ0[iDDL][iRow] = new Int_t*[AliHMPIDRawStream::kNDILOGICAdd+1];
          fsq[iDDL][iRow] = new Float_t*[AliHMPIDRawStream::kNDILOGICAdd+1];
         fsq2[iDDL][iRow] = new Float_t*[AliHMPIDRawStream::kNDILOGICAdd+1];
         fnpc[iDDL][iRow] = new Int_t*[AliHMPIDRawStream::kNDILOGICAdd+1];
      
        for(Int_t iDil=1;iDil<AliHMPIDRawStream::kNDILOGICAdd+1;iDil++){
      
         fpedQ0[iDDL][iRow][iDil] = new Int_t[AliHMPIDRawStream::kNPadAdd+1];
           fsq2[iDDL][iRow][iDil] = new Float_t[AliHMPIDRawStream::kNPadAdd+1];
            fsq[iDDL][iRow][iDil] = new Float_t[AliHMPIDRawStream::kNPadAdd+1];
           fnpc[iDDL][iRow][iDil] = new Int_t[AliHMPIDRawStream::kNPadAdd+1];
          }//iDil
      }//iRow
   }//iDDL
    
   for(Int_t iDDL=0;iDDL<AliHMPIDRawStream::kNDDL+1;iDDL++) {
        
     fnDDLInStream[iDDL]=-1;
     fnDDLOutStream[iDDL]=-1;
      
     for(Int_t iErr=0;iErr<AliHMPIDRawStream::kSumErr+1;iErr++)  {fErr[iDDL][iErr]=0;}
         
     for(Int_t iRow=0;iRow<AliHMPIDRawStream::kNRows+1;iRow++) {
        for(Int_t iDil=1;iDil<AliHMPIDRawStream::kNDILOGICAdd+1;iDil++) {
          for(Int_t iPad=1;iPad<AliHMPIDRawStream::kNPadAdd+1;iPad++) {
            fpedQ0[iDDL][iRow][iDil][iPad]=0;
               fsq[iDDL][iRow][iDil][iPad]=0;
              fsq2[iDDL][iRow][iDil][iPad]=0;
              fnpc[iDDL][iRow][iDil][iPad]=0;
        }//iPad
      }//iDil
     }//iRow
   }//iDDL
    
  fPadAdc=new TH1I*[nPads];  
  fIsPad=new Bool_t[nPads];  
  for(Int_t np=0;np<nPads;np++) {fPadAdc[np]=0x0;   fIsPad[np]=kFALSE;}

  Init();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDCalib::~AliHMPIDCalib()
{
  //
  //destructor
  //
  if (faddl)     { delete [] faddl;   faddl = 0x0;  } 
  if (fPadAdc)   { delete [] fPadAdc; fPadAdc=0x0;  }  
  if (fIsPad)    { delete [] fIsPad;  fIsPad=0x0;   }  
  if (fFile)     { delete    fFile;   fFile=0x0;    }  
  if (fDeadMap)  { delete    fDeadMap;fDeadMap=0x0; } 
  
  for(Int_t iErr=0;iErr<AliHMPIDRawStream::kSumErr+1;iErr++) { delete [] fErr[iErr];}  delete [] fErr;
  
  for(Int_t iDDL=0; iDDL< AliHMPIDRawStream::kNDDL; iDDL++) 
   for(Int_t iRow=0;iRow<AliHMPIDRawStream::kNRows+1;iRow++)         
     for(Int_t iDil=1;iDil<AliHMPIDRawStream::kNDILOGICAdd+1;iDil++)
      {
         delete [] fpedQ0[iDDL][iRow][iDil]; //del iPad
         delete []    fsq[iDDL][iRow][iDil]; //del iPad
         delete []   fsq2[iDDL][iRow][iDil]; //del iPad
         delete []   fnpc[iDDL][iRow][iDil]; //del iPad
       }
   for(Int_t iDDL=0; iDDL< AliHMPIDRawStream::kNDDL; iDDL++) 
     for(Int_t iRow=0;iRow<AliHMPIDRawStream::kNRows+1;iRow++)         
      {
        delete [] fpedQ0[iDDL][iRow];  //del iRow
          delete []  fsq[iDDL][iRow];  //del iRow
          delete [] fsq2[iDDL][iRow];  //del iRow
          delete [] fnpc[iDDL][iRow];  //del iRow
        }
       
   for(Int_t iDDL=0; iDDL< AliHMPIDRawStream::kNDDL; iDDL++) 
   {   
       delete [] fpedQ0[iDDL];        //del iRow
         delete [] fsq2[iDDL];        //del iRow
         delete []  fsq[iDDL];        //del iRow
         delete [] fnpc[iDDL];        //del iRow
     }
       
   delete [] fpedQ0;
   delete [] fsq2;
   delete [] fsq;
   delete [] fnpc;
    
  fpedQ0=0;    
    fsq2=0;
     fsq=0;
    fnpc=0;
    
  fLdcId=0;
  fTimeStamp=0;
  fRunNum=0;
  fSigCut=0;
  fLargeHisto=kFALSE;
  fSelectDDL=0;
  fNumMaskedPads=0;
  
  if (fPedMeanMap)   { delete [] fPedMeanMap; fPedMeanMap=0x0;  }  
  if (fPedSigMap)    { delete [] fPedSigMap;  fPedSigMap=0x0;   }  
  if (f1DPedMean)    { delete [] f1DPedMean;  f1DPedMean=0x0;   }
  if (f1DPedSigma)   { delete [] f1DPedSigma; f1DPedSigma=0x0;  }
  
}//dtor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::Init()
{
  //
  //Init the q calc.
  //Arguments: none
  //Return: none
  //
    
  fSigCut=3;  //the standard cut

 
  for(Int_t iDDL=0; iDDL< AliHMPIDRawStream::kNDDL; iDDL++) 
      {
         for(Int_t ierr=0; ierr <AliHMPIDRawStream::kSumErr ; ierr++) {
            fErr[iDDL][ierr]=0;
            }
        
        faddl[iDDL]=kFALSE;
      }//DDL
      
     Int_t      nbins[4]={14, 24, 10, 48};
     Double_t  binmin[4]={ 0,  1,  1,  0};
     Double_t  binmax[4]={14, 25, 11, 48};
     fDeadMap = new THnSparseD("fDeadMap","Dead Channel Map",4,nbins,binmin,binmax); 
     
   fPedMeanMap = new TH2F*[AliHMPIDParam::kMaxCh+1];
   fPedSigMap  = new TH2F*[AliHMPIDParam::kMaxCh+1];
   f1DPedMean  = new TH1F*[(AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPc+1)];
   f1DPedSigma = new TH1F*[(AliHMPIDParam::kMaxCh+1)*(AliHMPIDParam::kMaxPc+1)];
   
//   Int_t mapmin,mapmax;
   for(Int_t iCh=0;iCh<AliHMPIDParam::kMaxCh+1;iCh++)  
   {
    fPedMeanMap[iCh]=new TH2F(Form("fPedMeanMap%d",iCh),Form("fPedMeanMap%d;pad x;pad y;Mean pedestal (ADC)",iCh),160, 0,160,144,0,144);fPedMeanMap[iCh]->SetStats(kFALSE);
    fPedSigMap[iCh] =new TH2F(Form("fPedSigMap%d",iCh), Form("fPedSigMap%d;pad x;pad y;Sigma pedestal (ADC)",iCh), 160, 0,160,144,0,144);fPedSigMap[iCh]->SetStats(kFALSE);
   }
   for(Int_t iCh=0;iCh<=AliHMPIDParam::kMaxCh;iCh++)
   {
    for(Int_t iFee=0;iFee<6;iFee++)
     {
      f1DPedMean[6*iCh+iFee]  = new TH1F(Form("f1DPedMean_Ch%d_FEE_%d" , iCh,iFee),Form("Mean Pedestals, RICH %d, FEE %d;Channels;Mean pedestal (ADC)" ,iCh,iFee),3840,0,3840);f1DPedMean[6*iCh+iFee]->SetStats(kFALSE);
      f1DPedSigma[6*iCh+iFee] = new TH1F(Form("f1DPedSigma_Ch%d_FEE_%d" ,iCh,iFee),Form("Sigma Pedestal, RICH %d, FEE %d;Channels;Sigma pedestal (ADC)" ,iCh,iFee),3840,0,3840);f1DPedSigma[6*iCh+iFee]->SetStats(kFALSE);
     } 
   }
     
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
void AliHMPIDCalib::SetSigCutFromFile(TString hmpInFile)
{
  //
  //Set Sigma Cut from the file on the LDC, if the input file is not present default value is set!
  //Arguments: the name of the SigmaCut file on the LDC
  //Returns: none
  //
  Int_t nSigCut=0;
  ifstream infile(hmpInFile.Data());
  if(!infile.is_open()) {fSigCut=3; return;}
  while(!infile.eof())
    {
    infile>>nSigCut;
  }
  infile.close();
  if(nSigCut< 0 || nSigCut > 15 ) {Printf("WARNING: DAQ Sigma Cut from DAQ DB is out of bounds: %d, resetting it to 3!!!",nSigCut);nSigCut=3;}
  Printf("DAQ Sigma Cut from DAQ DB is: %d",nSigCut);
  fSigCut=nSigCut; 
}//SetSigCutFromFile()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::SetDeadChannelMapFromFile(TString hmpInFile)
{
  //
  //Set Dead Channel Map Cut from the file on the LDC, if the input file is not present default value is set!
  //Arguments: the name of the Dead Channel Map file on the LDC
  //Returns: none
  //
  Char_t header[256];
  Int_t ddl=0,row=0,dil=0,pad=0,ch=0,pc=0,chpadx=0,chpady=0,px=0,py=0,isitmasked=0;
  UInt_t dw=0;
  Double_t bin[4];
  ifstream infile(hmpInFile.Data());
  if(!infile.is_open()) {Printf("HMPID Dead Channel Map file cannot be opened!!!! No mask is applied!!!");return;}
  infile.getline(header,256);
  AliHMPIDDigit dig;
  while(!infile.eof())
    {
    infile>>ch>>chpadx>>chpady>>isitmasked;                     //read in masked coordinates; coordinates are in the module coordinate system
    pc=(chpadx/80)+2*(chpady/48);                               //get PC number
    px=chpadx-80*(chpadx/80);                                   //get pad X in PC coordinates
    py=chpady-48*(chpady/48);                                   //get pad Y in PC coordinates --- can we do it better??? -- just with one conversion???? clm
    if(!dig.Set(ch,pc,px,py,0) && isitmasked) {                 //in the AliHMPIDDigit:Set there is already a check if the coordinates makes sense
       dig.Raw(dw,ddl,row,dil,pad);
       bin[0]=ddl; bin[1]=row; bin[2]=dil; bin[3]=pad;
       fDeadMap->Fill(bin,1);
      }
    }
   infile.close();
  
}//SetDeadChannelMapFromFile()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCalib::FillPedestal(Int_t abspad,Int_t q)
{
  //
  //Called from the HMPIDda and fills the pedestal values
  //Arguments: absolute pad number as from AliHMPIDParam and q-charge
  //Returns: none
  //
  if(q<0) {
   AliError("Negative charge is read!!!!!!");
   return;
  }
  UInt_t w32;
  Int_t nDDL=0, row=0, dil=0, adr=0;
  //The decoding (abs. pad -> ddl,dil,...) is the same as in AliHMPIDDigit::Raw
  
  AliHMPIDDigit dig(abspad,q);
  dig.Raw(w32,nDDL,row,dil,adr);
  
  //........... decoding done      

     if(q>0) { 
        fsq[nDDL][row][dil][adr]+=q;
      fsq2[nDDL][row][dil][adr]+=q*q;
      fnpc[nDDL][row][dil][adr]++;                                                     //Count how many times the pad is good (can be different from the good DDL  count)
                       faddl[nDDL]=kTRUE; 
                     }
      else
      {
        fpedQ0[nDDL][row][dil][adr]++;                                                 //Count how many times a pad charge is zero
      }
      
     Int_t histocnt=0;   histocnt=(nDDL)*11520+(row-1)*480+(dil-1)*48+adr;             //Histo counter for a single DDL  
    
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
void AliHMPIDCalib::FillDDLCnt(Int_t iddl,Int_t inDDL, Int_t outDDL)
{
  //
  //Fill decoding DDL check from RawStream
  //Arguments: iddl - DDL under setting, inDDL- How many times the DDL is present in the raw stream, outDDL - How many time sthe DDL is succesfylly decoded
  //Retutns: none
  //
 
  if(inDDL==-1) return;
  if(fnDDLInStream[iddl]==-1) {fnDDLInStream[iddl]=0; fnDDLOutStream[iddl]=0;}
  fnDDLInStream[iddl]+=inDDL;
  fnDDLOutStream[iddl]+=outDDL;
 
  
}//FillDDLCnt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::WriteErrors(Int_t nDDL, Char_t* name, Int_t nEv)
{
  //
  //Write decoding errors to a txt file
  //Arguments: nDDL-DDL number, name of the error file and number of the read events
  //Retutns: kTRUE/kFALSE
  //  
  if(faddl[nDDL]==kFALSE) return kFALSE;                                                                 //if ddl is missing no error file is created
  ofstream outerr;  outerr.open(name);                                                                   //open error file
  outerr << Form("%8s %2d\n","RunNumber",(Int_t)fRunNum);                                                //read run number
  outerr << Form("%8s %2d\n","LdcId" ,          fLdcId);                                                 //read LDC Id
  outerr << Form("%8s %2d\n","TimeStamp",       fTimeStamp);                                             //read time stamp
  outerr << Form("%8s %2d\n","TotNumEvt",       nEv);                                                    //read number of total events processed
  outerr << Form("%8s %2d\n","TotDDLEvt",       fnDDLInStream[nDDL]);                                    //read number of bad events for DDL # nDDL processed
  outerr << Form("%8s %2d\n","NumBadEvt",       fnDDLInStream[nDDL]-fnDDLOutStream[nDDL]);               //read number of bad events for DDL # nDDL processed
  outerr << Form("%8s %2.2f\n","NBadE(%)",      (fnDDLInStream[nDDL]-fnDDLOutStream[nDDL])*100.0/nEv);   //read number of bad events (in %) for DDL # nDDL processed
  
  for(Int_t  ierr=0; ierr <AliHMPIDRawStream::kSumErr; ierr++) outerr << Form("%2d\t",fErr[nDDL][ierr]); //write errors
                                                               outerr << Form("\n");                     //last break
  /* write out pads with 0 charge read */
  for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){
    for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){
      for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
        if(fpedQ0[nDDL][row][dil][pad]>0) outerr<< Form("%2d %2d %2d %3d\n",row,dil,pad,fpedQ0[nDDL][row][dil][pad]);
      }
    }
  } 
                                                                                                                                                                                       
                                                               
  outerr.close();                                                                                        //write error file
  
  return kTRUE;
    
}//FillErrors()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCalib::CalcPedestal(Int_t nDDL, Char_t* name, Char_t *name2,Int_t nEv)    
{
  //
  //Calculate pedestal for each pad  
  //Arguments: nDDL-DDL number, name of the pedestal file and number of the read events
  //Retutns: kTRUE/kFALSE
  //
  
  if(faddl[nDDL]==kFALSE) return kFALSE;                   //if ddl is missing no ped file is created (and also for LDC selection). Check with Paolo what he checks for?!  

  Int_t feeOffset=196657;
  ofstream feeInput; feeInput.open(Form("%s",name2));      //write thr file for Fe2C
  
  Double_t mean=0,sigma=0, threshold=0;
  Double_t qs2m=0,qsm2=0;
  ofstream out;                                            //to write the pedestal text files
  Int_t inhard;
  Int_t nEvPerPad=0;
  Int_t pedbin=0;
  
  Int_t abspad,ch,pc,pcx,pcy,chX,chY,fee;
  Int_t binSp[4]={0};
  out.open(name);
  out << Form("%8s %2d\n","RunNumber",(Int_t)fRunNum);                                                //read run number
  out << Form("%8s %2d\n","LdcId" ,         fLdcId);                                                  //read LDC Id
  out << Form("%8s %2d\n","TimeStamp",      fTimeStamp);                                              //read time stamp
  out << Form("%8s %2d\n","TotNumEvt",      nEv);                                                     //read number of total events processed
  out << Form("%8s %2d\n","TotDDLEvt",      fnDDLInStream[nDDL]);                                     //read number of bad events for DDL # nDDL processed
  out << Form("%8s %2d\n","NumBadEvt",      fnDDLInStream[nDDL]-fnDDLOutStream[nDDL]);                //read number of bad events for DDL # nDDL processed
  out << Form("%8s %2f\n","NBadE(%)",       (fnDDLInStream[nDDL]-fnDDLOutStream[nDDL])*100.0/nEv);    //read number of bad events (in %) for DDL # nDDL processed
  out << Form("%8s %d\n","#SigCut",      fSigCut);                                                    //# of sigma cuts
      
  for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){
    feeInput << Form("0xabcdabcd \n");                                                                    //before each row we write a marker to separate the rows within a DDL                       
    
   
    for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){
      for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
        mean  = 50;sigma = 100;                                                   //init maen and sigma to a low value
        nEvPerPad=fnpc[nDDL][row][dil][pad];                                      //check how many times the pad was read out
        abspad=AliHMPIDRawStream::GetPad(nDDL,row,dil,pad);                       //get the absolute oad coordinate
        ch=AliHMPIDParam::A2C(abspad);                                            //get chamber number
        pc=AliHMPIDParam::A2P(abspad);                                            //get PC number
        pcx=AliHMPIDParam::A2X(abspad);                                           //get pad x in PC
        pcy=AliHMPIDParam::A2Y(abspad);                                           //get pad y in PC
        chX = (pc%2)*AliHMPIDParam::kPadPcX+pcx;                                  //get pad x in Ch   
        chY = (pc/2)*AliHMPIDParam::kPadPcY+pcy;                                  //get pad y in Ch
        binSp[0]=nDDL+1;binSp[1]=row;binSp[2]=dil;binSp[3]=pad+1;                 //set dead map coordinates for check
        
       if(nEvPerPad < 1 ) {                                                      //if the pad is bad then we assign 100  for the sigma and 50 for the mean
          mean  = AliHMPIDParam::kPadMeanZeroCharge;
          sigma = AliHMPIDParam::kPadSigmaZeroCharge;
          fNumDeadPads++;
        }
        else if(fDeadMap->GetBinContent(binSp)>0)                                 //check if channel is masked, if yes set maksed values
        {
          mean  = AliHMPIDParam::kPadMeanMasked;
          sigma = AliHMPIDParam::kPadSigmaMasked;
          fNumMaskedPads++;
        }
       else{            
         mean = fsq[nDDL][row][dil][pad]*1.0/nEvPerPad;
         qs2m = fsq2[nDDL][row][dil][pad]*1.0/nEvPerPad;
         qsm2 = TMath::Power(fsq[nDDL][row][dil][pad]*1.0/nEvPerPad,2); 
        sigma = TMath::Sqrt(TMath::Abs(qs2m-qsm2));
        }
        
        //The electronics takes the 32bit int as: first 9 bits for the pedestal and the second 9 bits for threshold
        threshold = mean+fSigCut*sigma;                                                                    
        if(mean > 511.0 || threshold > 511.0) {mean = AliHMPIDParam::kPadMeanMasked; threshold = AliHMPIDParam::kPadMeanMasked + 5.0 * AliHMPIDParam::kPadSigmaMasked; }
        //inhard=((Int_t(mean+fSigCut*sigma))<<9)+Int_t(mean);                                                 //right calculation, xchecked with Paolo 8/4/2008
        inhard=((Int_t(threshold))<<9)+Int_t(mean);                                                            //right calculation, xchecked with Paolo 8/4/2008
        
        out << Form("%2i %2i %2i %5.3f %5.3f %4.4x \n",row,dil,pad,mean,sigma,inhard);
        feeInput << Form("0x%4.4x\n",inhard);
        
        // fill histograms to be exported to AMORE    
        fPedMeanMap[ch]->SetTitle(Form("PedMeanMap%d RunNum: %d",ch,fRunNum));
        fPedSigMap[ch]->SetTitle(Form("PedSigmaMap%d RunNum: %d",ch,fRunNum));
        fPedMeanMap[ch]->Fill(chX,chY,mean);         
        fPedSigMap[ch]->Fill(chX,chY,sigma);         
        if(nDDL%2==0) pedbin = (24-row)*2*480+(10-dil)*48+pad;
        if(nDDL%2!=0) pedbin = (row*2-1)*480+(10-dil)*48+pad;
        pedbin = pedbin - 3840*(pedbin/3840);
        fee=AliHMPIDRawStream::GetFee(nDDL,row);
        f1DPedMean[6*(nDDL/2)+fee]->SetTitle(Form("PedMean_Ch%d_FEE_%d RunNum: %d",ch,fee,fRunNum));
        f1DPedSigma[6*(nDDL/2)+fee]->SetTitle(Form("PedSigma_Ch%d_FEE_%d RunNum: %d",ch,fee,fRunNum));
        f1DPedMean[6*(nDDL/2)+fee]->Fill(pedbin,mean);
        f1DPedSigma[6*(nDDL/2)+fee]->Fill(pedbin,sigma);

        }//adr==pad
        //we have to write up to 64 not 48 in the DILOGIC since they are daisy chained!
        //offset and format is defined for the Fe2C code
        for(Int_t idd=0;idd<16;idd++) feeInput << Form("0x%4.4x\n",idd+feeOffset);                 
      }//dil
      
      
    }//row
    out.close();                                          //write pedestal file
    feeInput.close();
 
  return kTRUE;
}//CaclPedestal()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
