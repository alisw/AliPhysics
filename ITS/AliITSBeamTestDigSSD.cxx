/////////////////////////////////////////
// Class for SSD raw2digits conv       //
//                                     //
// Author: Enrico Fragiacomo           //
// Date: October 2004                  //
////////////////////////////////////////

#include "AliITSdigitSSD.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSSDv1.h"
#include "AliITSBeamTestDigSSD.h"
#include "AliITSBeamTest.h"

ClassImp(AliITSBeamTestDigSSD)

//_____________________________________________________________
AliITSBeamTestDigSSD::AliITSBeamTestDigSSD(): AliITSBeamTestDig() {
  //
  // Constructor
  //  
}

//_____________________________________________________________
AliITSBeamTestDigSSD::AliITSBeamTestDigSSD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title) {
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

  TBranch* branch = fTreeD->GetBranch("ITSDigitSSD");

  TClonesArray** newdigits = new TClonesArray*[fBt->GetNSSD()];

  Int_t* idig = new Int_t[fBt->GetNSSD()];

  for (Int_t idet =0; idet < fBt->GetNSSD();idet++) {
     newdigits[idet] = new TClonesArray("AliITSdigitSSD");
     idig[idet]=0;  
   }
  
  // this constructor sets the flag to select SSD data only 
  // the Next method below will then jump to SSD data for this event

  AliITSRawStreamSSDv1 str(fReaderDate);

  // no selection of equipment 
  //fReaderDate->SelectEquipment(-1);
  //fReaderDate->SelectEquipment(17,102,102);

  while(str.Next()){   
    
    //if((str.GetADModule()!=2)&&(str.GetADModule()!=6)) continue;
  
    Int_t side = str.GetSideFlag();
    Int_t strip = str.GetStrip();
    Int_t signal = str.GetSignal();
    Int_t module = str.GetModuleID();
    if( (module<10) || (module>13) ) continue;
    const Int_t kdgt[3]={side,strip,signal};
    
    //  SSD modules 10, 11, 12 and 13
    new ( (*newdigits[module-10])[idig[module-10]] ) AliITSdigitSSD(kdgt);    
    idig[module-10]++;
    
  } // end while
  
  for(Int_t n=0;n<fBt->GetNSSD();n++){
    branch->SetAddress(&newdigits[n]);
    branch->Fill();  
  }
  
  fTreeD->SetEntries(fBt->GetNSPD()+fBt->GetNSDD()+fBt->GetNSSD());
    
  fReaderDate->Reset();
  
  fTreeD->AutoSave();
  
  for(Int_t n=0;n<fBt->GetNSSD();n++){
    delete newdigits[n];
  }
  
}

  

