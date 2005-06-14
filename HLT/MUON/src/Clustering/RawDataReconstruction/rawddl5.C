// ROOT includes
// #include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TParticle.h"
#include "TTree.h"
// //
// // // STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
// //
// // // MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
// // //
// //
// 
void rawddl5(Int_t evNumber1=0,Int_t evNumber2=0) 
{
  long int mappheader2,mapp2[123456];
  FILE *fb = fopen("lutraw51.dat","r");
  for(Int_t n=0;n<59136;n++)
    {
     fscanf(fb,"%ld",&mappheader2);
     fscanf(fb,"%ld",&mapp2[mappheader2]);
    }
  fclose(fb);
  long int mappheader4,mapp4[123456];
  FILE *fc = fopen("lutraw52.dat","r");
  for(Int_t n=0;n<59136;n++)
    {
     fscanf(fc,"%ld",&mappheader4);
     fscanf(fc,"%ld",&mapp4[mappheader4]);
    }
  fclose(fc);
	FILE *fp = fopen("ddl17.dat","w+");		 
	FILE *fq = fopen("ddl18.dat","w+");		 
	FILE *fr = fopen("ddl19.dat","w+");		 
	FILE *fs = fopen("ddl20.dat","w+");
   AliRunLoader * RunLoader = AliRunLoader::Open("galice.root","MUONFolder","REA
	   D");
  if (RunLoader ==0x0) {
                printf(">>> Error : Error Opening %s file \n","galice.root");
        return;
                      }
  // Loading MUON subsystem
     AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
        MUONLoader->LoadDigits("READ");

  // Creating MUON data container
     AliMUONData muondata(MUONLoader,"MUON","MUON");

     Int_t ievent, nevents;
     nevents = RunLoader->GetNumberOfEvents();
       printf(">>> No. of Event %d \n",nevents);
       AliMUONDigit * mDigit;
       
        Int_t dspheader[1000][4];
        Int_t dsphead[1000][100];
//   Start loop over events
  for (int ievent=evNumber1; ievent<= evNumber2; ievent++) // event start 
      {
	      printf("Event:%d\n",ievent+1);
    	  for(Int_t ii=0;ii<1000;ii++)
            {
          for(Int_t ij=0;ij<100;ij++)
             {
              dsphead[ii][ij]=0;
             }
            }
        for(Int_t ik=0;ik<1000;ik++)
            {
          for(Int_t il=0;il<4;il++)
             {
              dspheader[ik][il]=0;
             }
            }
       RunLoader->GetEvent(ievent);

       muondata.SetTreeAddress("D");

       Int_t ncathodes=2;
       for(Int_t icathode=0; icathode<ncathodes; icathode++) {
       muondata.GetCathode(icathode);

       Int_t ichamber;
       for( ichamber=8; ichamber<10; ichamber++) {
       Int_t idigit, ndigits;
       ndigits = (Int_t) muondata.Digits(ichamber)->GetEntriesFast(); //9 or 10
       for(idigit=0; idigit<ndigits; idigit++) {
       mDigit = static_cast<AliMUONDigit*>(muondata.Digits(ichamber)->At(idigit)
);
                   {                                              // pads start
               Int_t dsp = 0, xydata=0, manu=0, adc=0;
               Int_t index=0, counter=0;
                    Int_t iqpad = mDigit->Signal();              // charge per pad
                         if(iqpad>=4096)
	                           iqpad = 4095;
                    Int_t iPx = mDigit->PadX();               // pad number on X
                    Int_t iPy = mDigit->PadY();               // pad number on Y
									      
	Int_t modx = abs(iPx);
        if(icathode==0)
         {
 xydata = 2048*modx+iPy;
 manu=((mapp2[xydata]>>6)&1023);
 adc=(((mapp2[xydata]<<11)>>11)&63);
        dsp=iPy/80;
         }
              if(icathode == 1)
		{
 xydata = 256*modx+iPy;
 manu=((mapp4[xydata]>>6)&1023);
if(manu>0){
 adc=(((mapp4[xydata]<<11)>>11)&63);
if(manu<=625)
  dsp=0;
if(manu>625&&manu<=637)
  dsp=1;
if(manu>637&&manu<=656)
  dsp=2;
if(manu>656&&manu<=680)
  dsp=3;
if(manu>680&&manu<=716)
  dsp=4;
if(manu>716&&manu<=752)
  dsp=5;
if(manu>752&&manu<=788)
  dsp=6;
if(manu>788&&manu<=824)
  dsp=7;
if(manu>824&&manu<=860)
  dsp=8;
if(manu>860&&manu<=884)
  dsp=9;
if(manu>884&&manu<=903)
  dsp=10;
if(manu>903&&manu<=915)
  dsp=11;
if(manu>915&&manu<=924)
  dsp=12;
          }
                }
if(iPx<0)
  {
   if(ichamber==8)
     {
       dsp+=601;
       index=0;
     }
   if(ichamber==9)
     {
       dsp+=641;
       index=2;
     }
   dspheader[dsp][index]+=1;
   counter = dspheader[dsp][index];
   dsphead[dsp][counter]=4096*(64*(manu)+adc)+iqpad; 
  }
else
  {
   if(ichamber==8)
     {
       dsp+=621;
       index=1;
     }
   if(ichamber==9)
     {
       dsp+=661;
       index=3;
     }
    dspheader[dsp][index]+=1;
    counter = dspheader[dsp][index];
    dsphead[dsp][counter]=4096*(64*(manu)+adc)+iqpad; 
  }
               } //pad
	      }   //digit
	     }  // chamber
	    } //cathode
 for(Int_t dsp=601;dsp<680;dsp++)
    {
   for(Int_t aa=0;aa<4;aa++)
     {
     Int_t dspcount=dspheader[dsp][aa];
if(aa==0)
  {
     fwrite(&dsp,2,1,fp);
     fwrite(&dspcount,2,1,fp);
     for(Int_t ax=1;ax<=dspcount;ax++)
       fwrite(&dsphead[dsp][ax],4,1,fp);
  }
if(aa==1)
  {
     fwrite(&dsp,2,1,fq);
     fwrite(&dspcount,2,1,fq);
     for(Int_t ax=1;ax<=dspcount;ax++)
       fwrite(&dsphead[dsp][ax],4,1,fq);
  }
if(aa==2)
  {
     fwrite(&dsp,2,1,fr);
     fwrite(&dspcount,2,1,fr);
     for(Int_t ax=1;ax<=dspcount;ax++)
       fwrite(&dsphead[dsp][ax],4,1,fr);
  }
if(aa==3)
  {
     fwrite(&dsp,2,1,fs);
     fwrite(&dspcount,2,1,fs);
     for(Int_t ax=1;ax<=dspcount;ax++)
       fwrite(&dsphead[dsp][ax],4,1,fs);
  }
      }
     }
     } //event
fclose(fp);
fclose(fq);
fclose(fr);
fclose(fs);
}  // end of funtion	 
