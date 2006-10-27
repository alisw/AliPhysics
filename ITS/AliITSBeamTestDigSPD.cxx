////////////////////////////////////////////////////
//  Class to define                               //
//  SPD beam test raw 2 dig conv.                 //
//                                                //
//  Origin: Jan Conrad  Jan.Conrad@cern.ch        //
//                                                //
////////////////////////////////////////////////////

#include "AliITSdigitSPD.h"
#include "AliRawReader.h"
#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSBeamTestDigSPD.h"
#include "AliITSgeom.h"
#include <TBranch.h>
#include <TTree.h>
#include "AliITSEventHeader.h"
ClassImp(AliITSBeamTestDigSPD)



//_____________________________________________________________
  AliITSBeamTestDigSPD::AliITSBeamTestDigSPD(): AliITSBeamTestDig(),
fFlagHeader(0){
  //
  // Constructor
  //
}

//_____________________________________________________________
AliITSBeamTestDigSPD::AliITSBeamTestDigSPD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title),
fFlagHeader(0){
  //
  // Constructor
  //
}

//__________________________________________________________________
AliITSBeamTestDigSPD::~AliITSBeamTestDigSPD()
{
  //
  // Destructor
  //

 }


//_______________________________________________________________________
void AliITSBeamTestDigSPD::Exec(Option_t* /*opt*/)
{
  //Reads raw data for SPD, fill SPD digits tree


  if(!fITSgeom){
    Error("Exec","fITSgeom is null!");
    return;
  }
  TBranch* branch = fTreeD->GetBranch("ITSDigitsSPD");

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

  for (Int_t idet =0; idet <maxn;idet++){
     newdigits[idet]=new TClonesArray("AliITSdigitSPD");
     idig[idet]=0;  
   }
  

  AliITSRawStreamSPD str(fReader);

  fReader->SelectEquipment(17,211,211);

  while(str.Next()){  

    const AliRawDataHeader* rdh = fReader->GetDataHeader();
    UChar_t blockAttributes = fReader->GetBlockAttributes();     
    UInt_t statusBits = fReader->GetStatusBits();     
    UInt_t orbitNumber = rdh->fEventID2; 			 
    UShort_t  bunchCross = rdh->fEventID1;      
      // UInt_t DataSize = rdh->fSize;				      
      //UChar_t L1TrigType = rdh->fL1TriggerType;			      
      //UInt_t MiniEvId = rdh->GetMiniEventID();			      
      // ULong64_t TriggerCL = rdh->GetTriggerClasses();     
      //ULong64_t ROI = rdh->GetROI();
      //      UChar_t Version =rdh->fVersion;                                 
      

    Int_t modID = str.GetModuleID();
    //    Int_t triggernumber = str.GetTriggerEventNumber();

    Int_t row = str.GetRow();
    Int_t col = str.GetColumn();

    const Int_t kdgt[3]={col,row,1};
    //    newdigits = new TClonesArray*[fBt->GetNSPD()];

    new ((*newdigits[modID])[idig[modID]]) AliITSdigitSPD(kdgt);
    
    idig[modID]++;
    
      fITSHeader->SetOrbitNumber(0,orbitNumber);
      fITSHeader->SetStatusBits(0,statusBits);
      fITSHeader->SetBlockAttributes(0,blockAttributes);
      fITSHeader->SetBunchCross(0,bunchCross);
      //fITSHeader->SetTriggerClass(0,TriggerCL);
      //fITSHeader->SetSubDet(0,
      //fITSHeader->SetMiniEvId(0,MiniEvId);
      //fITSHeader->SetVersion(0,Version);
      //fITSHeader->SetSubDet(0,Subdets);
      //fITSHeader->SetL1TriggerType(0,L1TrigType);

   // fITSHeader->SetOrbitNumberSPD(OrbitNumber);
         //printf("Bunch Crossing  = %x\n ",BunchCross);
     if ( blockAttributes != 0x3a ) {
       Info("Exec","Block Attribs  = %x\n ",blockAttributes);
     }  
    
     
    
    } // while(str.Next());
    

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

    delete[] newdigits;
    delete[] idig;

}

  
