#if !defined(__CINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include "AliTPCParamSR.h"
#include "AliTPCDigitsArray.h"
#include "AliSimDigits.h"
#include "AliTPCBuffer160.h"
#endif

int AliTPCAltro(Int_t eth=0){
  //eth is a threshold.
  //Digits stored into a file have an amplitude value greater than "eth"
  Int_t offset=1; //this should be equal to the threshold
  /*
    NB the amplitude values stored in the ALTRO file are shifted  by offset 
    because the range for each word goes from 0 to 1023, now due to zero suppression 
    values lower that the threshold never appear.
   */

  Int_t PSecNumber=-1;  //Previous Sector number
  Int_t PRowNumber=-1;  //Previous Row number  
  Int_t PPadNumber=-1;  //Previous Pad number
  Int_t PTimeBin=-1;    //Previous Time-Bin
  Int_t BunchLength=0;
  //AliTPCBuffer160 is used in write mode to generate AltroFormat.dat file
  AliTPCBuffer160 Buffer("AltroFormat.dat",1); 
  ULong_t Count=0;
  Int_t nwords=0;
  Int_t numPackets=0;
  
  const char * inFile_new = "galice.root";
  AliRunLoader *rl = AliRunLoader::Open(inFile_new,"Event","read");

  Int_t nevents=rl->GetNumberOfEvents();
  cout<<"Number of Events:"<<nevents<<endl;
  Int_t choice=0;
  do{
    cout<<"Insert the event number: "; 
    cin>>choice;
  }while (choice<=0 || choice>nevents);
  rl->GetEvent(choice-1);
  AliLoader *tpcloader=rl->GetLoader("TPCLoader");
  tpcloader->LoadDigits();
  TTree *digitsTree=tpcloader->TreeD();

  AliSimDigits digrows, *dummy=&digrows;
  digitsTree->GetBranch("Segment")->SetAddress(&dummy);
  Stat_t nrows = digitsTree->GetEntries();
  cout<<"Number of entries (rows):"<<nrows<<endl;
  // get the TPC parameters
  rl->CdGAFile();
  AliTPCParamSR* param = AliTPC::LoadTPCParam(gFile);
  if (!param)
    cout<<"No TPC parameter"<<endl;
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
  digarr->Setup(param);
  digarr->ConnectTree(digitsTree);
  //ofstream ftxt("Data.txt");
  for (Int_t n=0; n<nrows; n++) {
    Int_t sec,row; // sector and row number (in the TPC)
    AliSimDigits *digrow=(AliSimDigits*)digarr->LoadEntry(n);
    param->AdjustSectorRow(digrow->GetID(),sec,row);   

    //cout<<"Sector:"<<sec<<" Row:"<<row<<endl;
    digrow->First();
    do{
      Short_t dig=digrow->CurrentDigit();  //adc
      Int_t time=digrow->CurrentRow();     //time
      Int_t pad =digrow->CurrentColumn();  // pad 
      //cout<<"dig:"<<dig<<" time:"<<time<<" pad:"<<pad<<endl;
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
