////////////////////////////////////////////////////
//  Class to define                               //
//  SDD beam test raw 2 dig conv.                 //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                // 
//                                                //
////////////////////////////////////////////////////
#include "AliITSdigitSDD.h"
#include "AliRawReader.h"
#include "AliVMERawStream.h"
#include "AliITSRawStreamSDDv2.h"
#include "AliITSRawStreamSDDv3.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSBeamTestDigSDD.h"
#include "AliITSEventHeader.h"
#include "AliITSgeom.h"
#include <TBranch.h>
#include <TTree.h>

ClassImp(AliITSBeamTestDigSDD)

//_____________________________________________________________
  AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(): AliITSBeamTestDig(),
fSDDEvType(0),
fSubEventAttributes(0),
fBtPer(),
fThreshold(0),
fStreamer(0){
  //
  // Constructor
  //
  SetBtPeriod();
}

//_____________________________________________________________
AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title),
fSDDEvType(0),
fSubEventAttributes(0),
fBtPer(),
fThreshold(0),
fStreamer(0)
{
  //
  // Constructor
  //
  SetBtPeriod();
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

  TBranch* branch = fTreeD->GetBranch("ITSDigitsSDD");
  Int_t maxn=0;

  if(!fITSgeom){
    Error("Exec","fITSgeom is null!");
    return;
  }

  Int_t nsdd=0;
  Int_t nspd=0;
  Int_t nssd=0;
  for(Int_t nlay=1;nlay<=fITSgeom->GetNlayers();nlay++){
    for(Int_t nlad=1;nlad<=fITSgeom->GetNladders(nlay);nlad++){
      for(Int_t ndet=1;ndet<=fITSgeom->GetNdetectors(nlay);ndet++){
	Int_t index=fITSgeom->GetModuleIndex(nlay,nlad,ndet);
	if(fITSgeom->GetModuleTypeName(index)=="kSPD") nspd++;
	if(fITSgeom->GetModuleTypeName(index)=="kSDD") nsdd++;
	if(fITSgeom->GetModuleTypeName(index)=="kSSD") nssd++;
      }
    }
  }
  if(GetBtPeriod()==kAug04) maxn=nsdd;
  if(GetBtPeriod()==kNov04) maxn=nspd+nsdd+nssd;
  TClonesArray** digits = new TClonesArray*[maxn];
  Int_t* idig = new Int_t[maxn];
  for(Int_t idet=0;idet<maxn;idet++){
    digits[idet]=new TClonesArray("AliITSdigitSDD");
    idig[idet]=0;
  }

  switch(fBtPer){
  case kNov04:
    fStreamer = new AliITSRawStreamSDDv3(fReader);
    break;
  case kAug04:
    AliVMERawStream vmeStreamer(fReader);
    fReader->RequireHeader(kFALSE);
    while(fReader->ReadHeader()){
      fSubEventAttributes = fReader->GetSubEventAttributes();
    }
    
    fSDDEvType=GetEventType();
    if(fSDDEvType==1) fITSHeader->SetEventTypeSDD(kReal);
    if(fSDDEvType==2) fITSHeader->SetEventTypeSDD(kCalibration1);
    if(fSDDEvType==3) fITSHeader->SetEventTypeSDD(kCalibration2);
    fReader->Reset();
    fStreamer = new AliITSRawStreamSDDv2(fReader);
    break;
  }

 
  fStreamer->SetLowCarlosThreshold(fThreshold,0);
  fStreamer->SetLowCarlosThreshold(fThreshold,1);

  //from raw data the signal is already decompressed..
  //set compressed fSignal of AliITSdigitSDD to -1000
  //set expanded fSignalExpanded of AliITSdigitSDD equal to fStreamer.GetSignal() 
  while(fStreamer->Next()){   
    Int_t ndet =0;
    if(GetBtPeriod()==kNov04) ndet=fStreamer->GetChannel()+nspd;
    if(GetBtPeriod()==kAug04) ndet=fStreamer->GetChannel();
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
  for(Int_t n = 0;n<maxn;n++){
    branch->SetAddress(&digits[n]);
    branch->Fill();
  
  }
      
  fTreeD->SetEntries(maxn);
  fReader->Reset();
  fTreeD->AutoSave();

  for(Int_t n=0;n<maxn;n++){
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
 
  fSDDEvType=0;
  if(fSubEventAttributes[0]==0 && fSubEventAttributes[1]==0 && fSubEventAttributes[2]==0) fSDDEvType = 1;
  if(fSubEventAttributes[0]==1 && fSubEventAttributes[1]==0 && fSubEventAttributes[2]==0) fSDDEvType = 2;

  if(fSubEventAttributes[0]==2 && fSubEventAttributes[1]==0 && fSubEventAttributes[2]==0) fSDDEvType = 3;

  return fSDDEvType;
}

//______________________________________________________________________
AliITSBeamTestDigSDD::AliITSBeamTestDigSDD(const AliITSBeamTestDigSDD &bt):AliITSBeamTestDig(bt),
fSDDEvType(bt.fSDDEvType),
fSubEventAttributes(bt.fSubEventAttributes),
fBtPer(bt.fBtPer),
fThreshold(bt.fThreshold),
fStreamer(bt.fStreamer){
    // Copy constructor. 

}
//______________________________________________________________________
AliITSBeamTestDigSDD& AliITSBeamTestDigSDD::operator=(const AliITSBeamTestDigSDD &source){
    // Assignment operator. 
  this->~AliITSBeamTestDigSDD();
  new(this) AliITSBeamTestDigSDD(source);
  return *this;
}




