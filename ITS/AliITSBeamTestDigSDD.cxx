////////////////////////////////////////////////////
//  Class to define                               //
//  SDD beam test raw 2 dig conv.                 //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                // 
//                                                //
////////////////////////////////////////////////////
#include "AliITSdigitSDD.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSDDv2.h"
#include "AliITSRawStreamSDDv3.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSBeamTestDigSDD.h"
#include "AliITSBeamTest.h"

ClassImp(AliITSBeamTestDigSDD)

//_____________________________________________________________
  AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(): AliITSBeamTestDig()
{
  //
  // Constructor
  //
  SetBtPeriod();
  fSubEventAttributes=0;
  fThreshold=0;
  fStreamer=0;
}

//_____________________________________________________________
AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title)
{
  //
  // Constructor
  //
  SetBtPeriod();
  fSubEventAttributes=0;
  fThreshold=0;
  fStreamer=0;
}

//__________________________________________________________________
AliITSBeamTestDigSDD::~AliITSBeamTestDigSDD()
{
  //
  // Destructor
  //
  if(fSubEventAttributes) delete fSubEventAttributes;
  if(fStreamer) delete fStreamer;
 }


//_______________________________________________________________________
void AliITSBeamTestDigSDD::Exec(Option_t* /*opt*/)
{

  // Reads raw data and fills the tree of digits

  TBranch* branch = fTreeD->GetBranch("ITSDigitSDD");

  TClonesArray** digits = new TClonesArray*[fBt->GetNSDD()+fBt->GetNSPD()];
  Int_t* idig = new Int_t[fBt->GetNSDD()+fBt->GetNSPD()];
  for(Int_t idet=0;idet<(fBt->GetNSDD()+fBt->GetNSPD());idet++){
    digits[idet]=new TClonesArray("AliITSdigitSDD");
    idig[idet]=0;
  }


  switch(fBtPer){
  case kNov04:
    fStreamer = new AliITSRawStreamSDDv3(fReaderDate);
    break;
  case kAug04:
    fStreamer = new AliITSRawStreamSDDv2(fReaderDate);
    fReaderDate->RequireHeader(kFALSE);
    fReaderDate->ReadHeader();
   do{
      fSubEventAttributes = fReaderDate->GetSubEventAttributes();
    }while(fReaderDate->ReadHeader());
   fSDDEvType=GetEventType();
    if(fSDDEvType==1) fITSHeader->SetEventTypeSDD(kReal);
    if(fSDDEvType==2) fITSHeader->SetEventTypeSDD(kCalibration1);
    if(fSDDEvType==3) fITSHeader->SetEventTypeSDD(kCalibration2);
    fReaderDate->Reset();
    break;
  }


   fStreamer->SetLowCarlosThreshold(fThreshold,0);
   fStreamer->SetLowCarlosThreshold(fThreshold,1);

  //from raw data the signal is already decompressed..
  //set compressed fSignal of AliITSdigitSDD to -1000
  //set expanded fSignalExpanded of AliITSdigitSDD equal to fStreamer.GetSignal() 
  while(fStreamer->Next()){   

    Int_t ndet = fStreamer->GetChannel()+fBt->GetNSPD();
    Int_t anode = fStreamer->GetCoord1();

    /* if we are reading only one det, two wings
       if(fStreamer.GetChannel()==1) anode+=256; //wing 1 0-255, wing 2 256-511
    */

    /* bt august 2004 and november 2004: with only 1 carlos 
       channel 0 for one wing of one
       det, channel 1 for the wing of the second det*/

    const Int_t kdgt[3]={anode,fStreamer->GetCoord2(),-1000};
    const Int_t ktracks[10]={0,0,0,0,0,0,0,0,0,0};
    const Int_t khits[10]={0,0,0,0,0,0,0,0,0,0};
    const Float_t kcharges[10]={0,0,0,0,0,0,0,0,0,0};
 

    new ((*digits[ndet])[idig[ndet]]) AliITSdigitSDD(0,kdgt,ktracks,khits,kcharges,fStreamer->GetSignal());
    idig[ndet]++;

  }

  if(GetBtPeriod()==kNov04){
    Int_t jitter=fStreamer->ReadJitter();
    fITSHeader->SetJitterSDD(jitter);
  }
  for(Int_t n = fBt->GetNSPD();n<fBt->GetNSDD()+fBt->GetNSPD();n++){
    branch->SetAddress(&digits[n]);
    branch->Fill();
  }
      
  fTreeD->SetEntries(fBt->GetNSPD()+fBt->GetNSDD()+fBt->GetNSSD());
  fReaderDate->Reset();
  fTreeD->AutoSave();

  for(Int_t n=0;n<fBt->GetNSPD()+fBt->GetNSDD();n++){
    delete digits[n];
  }

  
  delete[] digits;
  delete[] idig;
  delete fStreamer;
}

  
//______________________________________
Int_t AliITSBeamTestDigSDD::GetEventType(){

  // defines the SDD event type:
  // 1: physics event (kReal)
  // 2: calibration 1 (kCalibration1, injector pulse)
  // 3: calibration 2 (kCalibration2, test pulse)
 
  fSDDEvType = 2;
  if(fSubEventAttributes[0]==0 && fSubEventAttributes[1]==0 && fSubEventAttributes[2]==0) fSDDEvType = 1;

  if(fSubEventAttributes[0]==2 && fSubEventAttributes[1]==0 && fSubEventAttributes[2]==0) fSDDEvType = 3;

  fSubEventAttributes = 0;
  return fSDDEvType;
}

//______________________________________________________________________
AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(const AliITSBeamTestDigSDD &bt):AliITSBeamTestDig(bt){
    // Copy constructor. 

  fSDDEvType = bt.fSDDEvType;
  fSubEventAttributes = bt.fSubEventAttributes;
  fBtPer = bt.fBtPer;
  fThreshold = bt.fThreshold;
  fStreamer = bt.fStreamer;
}
//______________________________________________________________________
AliITSBeamTestDigSDD& AliITSBeamTestDigSDD::operator=(AliITSBeamTestDigSDD &bt){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS beam test dig. It exits with an error.
    // Inputs:
    if(this==&bt) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSBeamTestDig");
    exit(1);
    return *this; //fake return
}




