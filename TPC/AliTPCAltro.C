#if !defined(__CINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include "AliTPCParamSR.h"
#include "AliTPCDigitsArray.h"
#include "AliSimDigits.h"
#include "AliTPCBuffer160.h"
#endif

int AliTPCAltro(char* FileName,Int_t eth=0){
  //eth is a threshold.
  //Digits stored into a file have an amplitude value greater than "eth"
  Int_t offset=1; //this should be equal to the threshold
  /*
    NB the amplitude values stored in the ALTRO file are shifted  by offset 
    because the range for each word goes from 0 to 1023, now due to zero suppression 
    values lower that the threshold never appear.
   */
  TFile *cf=TFile::Open(FileName);
  // old geometry (3.07)
  //AliTPCParamSR *param =(AliTPCParamSR *)cf->Get("75x40_100x60");
  // if new geometry comment out the line above and uncomment the one below
  AliTPCParamSR *param =(AliTPCParamSR *)cf->Get("75x40_100x60_150x60");
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
  digarr->Setup(param);
  
  char  cname[100];
  //old geometry
  //sprintf(cname,"TreeD_75x40_100x60_%d",eventn);
  // if new geometry comment out the line above and uncomment the one below
  Int_t eventn=0;
  sprintf(cname,"TreeD_75x40_100x60_150x60_%d",eventn);
  digarr->ConnectTree(cname);

  Int_t PSecNumber=-1;  //Previous Sector number
  Int_t PRowNumber=-1;  //Previous Row number  
  Int_t PPadNumber=-1;  //Previous Pad number
  Int_t PTimeBin=-1;    //Previous Time-Bin
  Int_t BunchLength=0;

  //AliTPCBuffer160 is used in write mode to generate AltroFormat.dat file
  AliTPCBuffer160 Buffer("AltroFormat.dat",1); 
  //number of entries in the tree
  Int_t nrows=Int_t(digarr->GetTree()->GetEntries());
  cout<<"Number of entries "<<nrows<<endl;
  ULong_t Count=0;
  Int_t nwords=0;
  Int_t numPackets=0;
  //ofstream ftxt("Data.txt");
  for (Int_t n=0; n<nrows; n++) {
    AliSimDigits *digrow=(AliSimDigits*)digarr->LoadEntry(n);
    Int_t sec,row; // sector and row number (in the TPC)
    param->AdjustSectorRow(digrow->GetID(),sec,row);   
    //cout<<"Sector:"<<sec<<" Row:"<<row<<endl;
    digrow->First();
    do{
      Short_t dig=digrow->CurrentDigit(); //adc
      Int_t time=digrow->CurrentRow(); //time
      Int_t pad =digrow->CurrentColumn(); // pad 
      if(dig>eth){
	Count++;
	//ftxt<<"Sec: "<<sec<<" Row: "<<row<<" Pad:"<<pad<<" Time: "<<time<<" ADC:"<<dig<<endl;
	//cout<<"Sec: "<<sec<<" Row: "<<row<<" Pad:"<<pad<<" Time: "<<time<<" ADC:"<<dig<<endl;
	if (PPadNumber==-1){
	  PSecNumber=sec;
	  PRowNumber=row;
	  PPadNumber=pad;
	  // PAmplitude=dig;
	  PTimeBin=time;
	  BunchLength=1;
	  Buffer.FillBuffer(dig-offset);
	  nwords++;
	}//end if
	else{
	  if ( (time==(PTimeBin+1)) &&
	       (PPadNumber==pad) &&
	       (PRowNumber==row) &&
	       (PSecNumber==sec)){
	    BunchLength++;
	  }//end if
	  else{
	    Buffer.FillBuffer(PTimeBin);
	    Buffer.FillBuffer(BunchLength+2);
	    nwords+=2;
	    if ((PPadNumber!=pad)||(PRowNumber!=row)||(PSecNumber!=sec)){
	      //Trailer is formatted and inserted!!
	      Buffer.WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
	      numPackets++;
	      nwords=0;
	    }//end if
	    
	    BunchLength=1;
	    PPadNumber=pad;
	    PRowNumber=row;
	    PSecNumber=sec;
	  }//end else
	  PTimeBin=time;
	  Buffer.FillBuffer(dig-offset);
	  nwords++;
	}//end else
      }//end if
    } while (digrow->Next());
  }//end for
  Buffer.FillBuffer(PTimeBin);
  Buffer.FillBuffer(BunchLength+2);
  nwords+=2;
  Buffer.WriteTrailer(nwords,PPadNumber,PRowNumber,PSecNumber);
  numPackets++;
  cout<<"There are "<<Count<<" Digits\n";
  cout<<"Packets "<<numPackets<<"\n";
  //ftxt.close();
  return 0;
}//end macro
