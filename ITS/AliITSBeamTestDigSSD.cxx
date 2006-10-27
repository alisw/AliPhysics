/////////////////////////////////////////
// Class for SSD raw2digits conv       //
//                                     //
// Author: Enrico Fragiacomo           //
// Date: October 2004                  //
////////////////////////////////////////

#include "AliITSgeom.h"
#include "AliITSdigitSSD.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSSDv1.h"
#include "AliITSBeamTestDigSSD.h"
#include <TBranch.h>
#include <TTree.h>

ClassImp(AliITSBeamTestDigSSD)

//_____________________________________________________________
AliITSBeamTestDigSSD::AliITSBeamTestDigSSD(): AliITSBeamTestDig(),
fFlagHeader(0) 
{
  //
  // Constructor
  //  
}

//_____________________________________________________________
AliITSBeamTestDigSSD::AliITSBeamTestDigSSD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title),
fFlagHeader(0){
  //
  // Constructor
  //
}

//__________________________________________________________________
AliITSBeamTestDigSSD::~AliITSBeamTestDigSSD() {
  //
  // Destructor
  //
}

//_______________________________________________________________________
void AliITSBeamTestDigSSD::Exec(Option_t* /*opt*/) {
  //Reads raw data for SSD, fill SSD digits tree, returns 1 if real data,
  //returns 2 if calibration (injector) data, returns 3 if calibration (test pul  //se) event

  TBranch* branch = fTreeD->GetBranch("ITSDigitsSSD");
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
  Int_t maxn=nspd+nsdd+nssd;

  TClonesArray** newdigits = new TClonesArray*[maxn];

  Int_t* idig = new Int_t[maxn];

  for (Int_t idet =0; idet <maxn;idet++) {
     newdigits[idet] = new TClonesArray("AliITSdigitSSD");
     idig[idet]=0;  
   }
  
  // this constructor sets the flag to select SSD data only 
  // the Next method below will then jump to SSD data for this event
  AliITSRawStreamSSDv1 str(fReader);

  // no selection of equipment 
  //fReader->SelectEquipment(-1);
  //fReader->SelectEquipment(17,102,102);

  while(str.Next()){   
    
    //if((str.GetADModule()!=2)&&(str.GetADModule()!=6)) continue;
  
    Int_t side = str.GetSideFlag();
    Int_t strip = str.GetStrip();
    Int_t signal = str.GetSignal();
    Int_t module = str.GetModuleID();
    Int_t module1= module-10+nsdd+nspd;
     if( (module<10) || (module>13) ) continue;
     if(module1>7) continue;
    const Int_t kdgt[3]={side,strip,signal};
     //  SSD modules 10, 11, 12 and 13
    new ( (*newdigits[module1])[idig[module1]] ) AliITSdigitSSD(kdgt);    
    idig[module1]++;
  } // end while
  
  for(Int_t n=0;n<maxn;n++){
    branch->SetAddress(&newdigits[n]);
    branch->Fill();  
  }
  
  fTreeD->SetEntries(maxn);
    
  fReader->Reset();
  
  fTreeD->AutoSave();
  
  for(Int_t n=0;n<maxn;n++){
    delete newdigits[n];
  }
  
  delete [] newdigits;
  delete idig;
}

  

