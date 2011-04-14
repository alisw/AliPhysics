#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliITSDDLModuleMapSDD.h"
#include <TSystem.h>
#endif

/*  $Id$    */

// Macro to extract and print the configuration parameters of SDD modules 
// from JTAG hexadecimal files 
// Origin: F. Prino,   prino@to.infn.it

void DecodeJTAG(Int_t nDDL, Int_t nCarlos=0, TString type="PHY"){

  AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
  dmap->SetJun09Map();
  Int_t theMod=-1;
  if(nDDL>=0 && nDDL<23){
    theMod=dmap->GetModuleNumber(nDDL,nCarlos);
  }else if(nDDL>=240 && nDDL<500){
    theMod=nDDL;
    dmap->FindInDDLMap(theMod,nDDL,nCarlos);
  }else{
    printf("ERROR: wronng DDL/Module number %d\n",nDDL); 
    return;
  }
  TString crate;
  if(nDDL<8) crate="TOP";
  else if(nDDL<16) crate="MED";
  else if(nDDL<24) crate="BOT";

  Int_t numSlot=nDDL%8+1;
  Int_t car1,car2;
  if(nCarlos%2==0){
    car1=nCarlos;
    car2=nCarlos+1;
  }else{
    car2=nCarlos;
    car1=nCarlos-1;
  }
  TString filnam=Form("%s%d/%s_jtag_carlos%02d%02d.txt",crate.Data(),numSlot,
		      type.Data(),car1,car2);
  gSystem->Exec(Form("ls -l %s\n",filnam.Data()));
  FILE* fil=fopen(filnam.Data(),"r");
  UInt_t word;
  Int_t retcode;

  retcode=fscanf(fil,"%x\n",&word);
  Int_t modNum=word&0x0000000F;
  for(Int_t iMod=0; iMod<2; iMod++){
    retcode=fscanf(fil,"%x\n",&word);
    if(iMod==0 && word!=0xFFFFFF00){printf("ERROR - mod 0\n"); return;}
    if(iMod==1 && word!=0xFFFFFF11){printf("ERROR - mod 1\n"); return;}
    retcode=fscanf(fil,"%x\n",&word); // CarlosRX word
    if(iMod==0){
      Int_t maskMod=(word&0x000FFF00)>>8;  
      printf("Active modules: ");
      for(Int_t i=0;i<12;i++){ 
	if(maskMod&(1<<i)) printf(" 1 ");
	else printf(" 0 ");
      }
      printf("\n");
    }
    // CARLOS words 
    Bool_t doPrint=kFALSE;
    if(modNum+iMod==nCarlos)doPrint=kTRUE;
    if(doPrint) printf("*********** Module %d DDL %d (%s%d) Carlos %d **********\n",theMod,nDDL,crate.Data(),numSlot,modNum+iMod);
    retcode=fscanf(fil,"%x\n",&word);
    retcode=fscanf(fil,"%x\n",&word);
    Int_t th1=(word&0xFF000000)>>24;
    Int_t th0=(word&0x00FF0000)>>16;
    Int_t tl1=(word&0x0000FF00)>>8;
    Int_t tl0=(word&0x000000FF);
    if(doPrint) printf("High Thresholds = %d %d\n",th1,th0);
    if(doPrint) printf("Low Thresholds  = %d %d\n",tl1,tl0);
    retcode=fscanf(fil,"%x\n",&word);
    retcode=fscanf(fil,"%x\n",&word);
    for(Int_t iSide=0; iSide<2; iSide++){
      if(doPrint) printf("------------------ Hybrid %d -------------\n",iSide);
      // AMBRA WORDS
      for(Int_t iAmbra=0; iAmbra<4; iAmbra++){
	for(Int_t iw=0;iw<5;iw++){
	  retcode=fscanf(fil,"%x\n",&word);
	  if(doPrint && iw==1) {
	    printf("AMBRA %d --- SOP %d",iAmbra,word&0x000000FF);
	  }
	}
	if(word&0x0000FF00){
	  if(doPrint) printf(" --- Baselines present");
	  for(Int_t iwb=0;iwb<16;iwb++) retcode=fscanf(fil,"%x\n",&word);
	}else{
	  if(doPrint) printf(" --- No baseline equalization");
	}
	if(doPrint) printf("\n");
      }
      // PASCAL WORDS
      for(Int_t iPascal=0; iPascal<4; iPascal++){
	for(Int_t iw=0;iw<4;iw++){
	  retcode=fscanf(fil,"%x\n",&word);
	  if(doPrint && iw==2){
	    printf("PASCAL %d ",iPascal);
	    if(word&0xFF000000) printf(" --- ADC full");
	    else printf(" --- ADC 1/2");
	    if(word&0x00FF0000) printf(" --- AM  full");
	    else printf(" --- AM 1/2");
	    printf("\n");
	  }
	} 
      }
    }
  }
}
