////////////////////////////////////////////////////
//  Class to define                               //
//  SPD beam test raw 2 dig conv.                 //
//                                                //
//  Origin: Jan Conrad  Jan.Conrad@cern.ch        //
//                                                //
////////////////////////////////////////////////////


#include "AliITSdigitSPD.h"
#include "AliRawReaderDate.h"
#include "AliRawDataHeader.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSBeamTestDigSPD.h"
#include "AliITSBeamTest.h"


ClassImp(AliITSBeamTestDigSPD)



//_____________________________________________________________
  AliITSBeamTestDigSPD::AliITSBeamTestDigSPD(): AliITSBeamTestDig()
{
  //
  // Constructor
  //

  
}

//_____________________________________________________________
  AliITSBeamTestDigSPD::AliITSBeamTestDigSPD(const Text_t* name, const Text_t* title): AliITSBeamTestDig(name,title)
{
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


  TBranch* branch = fTreeD->GetBranch("ITSDigitSPD");
  TClonesArray** newdigits = new TClonesArray*[fBt->GetNSPD()];

   Int_t* idig = new Int_t[fBt->GetNSPD()];

  

  for (Int_t idet =0; idet < fBt->GetNSPD();idet++){
     newdigits[idet]=new TClonesArray("AliITSdigitSPD");
     idig[idet]=0;  
   }
  

 

  //cout <<"still here"<< endl;
  AliITSRawStreamSPD str(fReaderDate);

  fReaderDate->SelectEquipment(17,211,211);

  while(str.Next()){  

    const AliRawDataHeader* rdh = fReaderDate->GetDataHeader();
    UChar_t blockAttributes = fReaderDate->GetBlockAttributes();     
    UInt_t statusBits = fReaderDate->GetStatusBits();     
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



    //    if (modID == 0 || modID == 1) {
    // cout <<"Mod ID " << modID  <<" Row : "<< row << "Col : " << col << endl; 
     //}


    const Int_t kdgt[3]={row,col,1};
    
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
    

    for(Int_t n=0;n<fBt->GetNSPD();n++){
      branch->SetAddress(&newdigits[n]);
      branch->Fill();  
   }
    
    fTreeD->SetEntries(fBt->GetNSDD()+fBt->GetNSPD()+fBt->GetNSSD());
  
	
  
    fReaderDate->Reset();
    fTreeD->AutoSave();
   

    for(Int_t n=0;n<fBt->GetNSPD();n++){
    delete newdigits[n];
    }

    delete[] newdigits;
    delete[] idig;

}

  
