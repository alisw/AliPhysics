

#include <Riostream.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TTree.h> 
#include <TMath.h>

#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONDigitizerv1.h"
#include "AliMUONHit.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONPadHit.h"
#include "AliMUONTransientDigit.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

ClassImp(AliMUONDigitizerv1)

//___________________________________________
AliMUONDigitizerv1::AliMUONDigitizerv1() :
  AliDigitizer(),
  fHitMap(0),
  fTDList(0),
  fTDCounter(0),
  fDebug(0),
  fMask(0),
  fSignal(0)
{
// Default ctor - don't use it
  if (GetDebug()>2) 
    cerr<<"AliMUONDigitizerv1::AliMUONDigitizerv1"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
}

//___________________________________________
AliMUONDigitizerv1::AliMUONDigitizerv1(AliRunDigitizer* manager):
  AliDigitizer(manager),
  fHitMap(0),
  fTDList(0),
  fTDCounter(0),
  fDebug(0),
  fMask(0),
  fSignal(0)
{
// ctor which should be used
  if (GetDebug()>2) 
    cerr<<"AliMUONDigitizerv1::AliMUONDigitizerv1"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
}

//------------------------------------------------------------------------
AliMUONDigitizerv1::~AliMUONDigitizerv1()
{
// Destructor
}

//------------------------------------------------------------------------
void AliMUONDigitizerv1::AddTransientDigit(AliMUONTransientDigit * mTD)
{
  // Choosing the maping of the cathode plane of the chamber:
  Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
  fTDList->AddAtAndExpand(mTD, fTDCounter);
  fHitMap[iNchCpl]->SetHit( mTD->PadX(), mTD->PadY(), fTDCounter);
  fTDCounter++;
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizerv1::ExistTransientDigit(AliMUONTransientDigit * mTD)
{
  // Choosing the maping of the cathode plane of the chamber:
  Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
  return( fHitMap[iNchCpl]->TestHit(mTD->PadX(), mTD->PadY()) );
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizerv1::Init()
{

// Initialization 
  if (GetDebug()>2) Info("Init","AliMUONDigitizerv1::Init() starts");

  //Loaders (We assume input0 to be the output too)
  AliRunLoader * runloader;    // Input   loader
  AliLoader    * gime; 
 
  // Getting runloader
  runloader    = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  if (runloader == 0x0) {
    Error("Init","RunLoader is not in input file 0");
    return kFALSE;  // RunDigitizer is not working.  
  }
  // Getting MUONloader
  gime        = runloader->GetLoader("MUONLoader");
  gime->LoadHits("READ");
  gime->LoadDigits("RECREATE");
 
  return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONDigitizerv1::UpdateTransientDigit(Int_t track, AliMUONTransientDigit * mTD)
{
  // Choosing the maping of the cathode plane of the chamber:
  Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
  AliMUONTransientDigit *pdigit = 
    static_cast<AliMUONTransientDigit*>(fHitMap[iNchCpl]->GetHit(mTD->PadX(),mTD->PadY()));
  // update charge
  //
  Int_t iqpad  = mTD->Signal();        // charge per pad
  pdigit->AddSignal(iqpad);
  pdigit->AddPhysicsSignal(iqpad);		
  // update list of tracks
  //
  Int_t charge;   
  track=+ fMask;
  if (fSignal) charge = iqpad;
  //else         charge = kBgTag;
  else         charge = iqpad + fMask; 
	       
  pdigit->UpdateTrackList(track,charge);
}


//--------------------------------------------------------------------------
void AliMUONDigitizerv1::MakeTransientDigit(Int_t track, Int_t iHit, AliMUONHit * mHit)
{
  AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
  if (!pMUON) { 
    cerr<<"AliMUONDigitizerv1::Digitize Error:"
	<<" module MUON not found in the input file"<<endl;
  }
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit starts"<<endl;
  Int_t   ichamber = mHit->Chamber()-1;
  AliMUONChamber & chamber = pMUON->Chamber(ichamber);
  Float_t xhit = mHit->X();
  Float_t yhit = mHit->Y();
  Float_t zhit = mHit->Z();
  Float_t eloss= mHit->Eloss(); 
  Float_t tof  = mHit->Age();
  // Variables for chamber response from AliMUONChamber::DisIntegration
  Float_t newdigit[6][500];  // Pad information
  Int_t nnew=0;              // Number of touched Pads per hit
  Int_t digits[6];
  
  //
  // Calls the charge disintegration method of the current chamber 
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit calling AliMUONChamber::DisIngtegration starts"<<endl;
  chamber.DisIntegration(eloss, tof, xhit, yhit, zhit, nnew, newdigit);
  // Creating a new TransientDigits from hit
  for(Int_t iTD=0; iTD<nnew; iTD++) {
    digits[0] = Int_t(newdigit[1][iTD]);  // Padx of the Digit
    digits[1] = Int_t(newdigit[2][iTD]);  // Pady of the Digit
    digits[2] = Int_t(newdigit[5][iTD]);  // Cathode plane
    digits[3] = Int_t(newdigit[3][iTD]);  // Induced charge in the Pad
    if (fSignal) digits[4] = Int_t(newdigit[3][iTD]);
    else         digits[4] = 0;
    digits[5] = iHit+fMask;    // Hit number in the list
    if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit " <<
			"PadX "<< digits[0] << " " <<
			"PadY "<< digits[1] << " " <<
			"Plane " << digits[2] << " " <<
			"Charge " << digits[3] <<" " <<
			"Hit " << digits[5] << endl;
    // list of tracks
    Int_t charge;   
    track += fMask;
    if (fSignal) charge = digits[3];
    //else         charge = kBgTag;
    else         charge = digits[3] + fMask;
		 
    if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit Creating AliMUONTransientDigit"<<endl;
    AliMUONTransientDigit * mTD = new AliMUONTransientDigit(ichamber, digits);
    mTD->AddToTrackList(track,charge);
    if (!ExistTransientDigit(mTD)) {
      AddTransientDigit(mTD);
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit Adding TransientDigit"<<endl;
    }
    else {
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::MakeTransientDigit updating TransientDigit"<<endl;
      UpdateTransientDigit(track, mTD);
      delete mTD;
    }
  }
}
//-----------------------------------------------------------------------
void AliMUONDigitizerv1::Exec(Option_t* option)
{
  TString optionString = option;
  if (optionString.Data() == "deb") {
    Info("Digitize","Called with option deb ");
    fDebug = 3;
  }

  AliMUONChamber*   chamber;
  AliSegmentation*  c1Segmentation; //Cathode plane c1 of the chamber
  AliSegmentation*  c2Segmentation; //Cathode place c2 of the chamber

  if (GetDebug()>2) Info("Digitize","AliMUONDigitizerv1::Digitize() starts");
  fTDList = new TObjArray;

  //Loaders (We assume input0 to be the output too)
  AliRunLoader * runloader;    // Input   loader
  AliLoader    * gime; 
 
  // Getting runloader
  runloader    = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  if (runloader == 0x0) {
    Error("Digitize","RunLoader is not in input file 0");
    return;  // RunDigitizer is not working.  
  }
  // Getting MUONloader
  gime        = runloader->GetLoader("MUONLoader");
  if (gime->TreeH()==0x0) {
      if (GetDebug()>2) Info("Digitize","TreeH is not loaded yet. Loading...");
     gime->LoadHits("READ");
       if (GetDebug()>2) Info("Digitize","Now treeH is %#x. MUONLoader is %#x",gime->TreeH(),gime);
  }

  if (GetDebug()>2) Info("Digitize","Loaders ready");

  if (runloader->GetAliRun() == 0x0) runloader->LoadgAlice();
  gAlice = runloader->GetAliRun();

  // Getting Module MUON  
  AliMUON *pMUON  = (AliMUON *) gAlice->GetDetector("MUON");
  if (!pMUON) {
    Error("Digitize","Module MUON not found in the input file");
    return;
  }
  // Getting Muon data
  AliMUONData * muondata = pMUON->GetMUONData(); 
  muondata->SetLoader(gime);
  muondata->SetTreeAddress("H");

  // Loading Event
  Int_t currentevent = fManager->GetOutputEventNr();
  
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::Digitize() Event Number is "<<currentevent <<endl;
  if ( (currentevent<10)                                                 || 
       (Int_t(TMath::Log10(currentevent)) == TMath::Log10(currentevent) ) )
    cout <<"ALiMUONDigitizerv1::Digitize() Event Number is "<< currentevent <<endl;

  // Output file for digits same as hits
  //  AliRunLoader *  runloaderout  = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  //AliLoader * gimeout         = runloaderout->GetLoader("MUONLoader"); 
  // New branch per chamber for MUON digit in the tree of digits
  if (gime->TreeD() == 0x0) {
    gime->MakeDigitsContainer();
  }
  TTree* treeD = gime->TreeD();
  muondata->MakeBranch("D");
  muondata->SetTreeAddress("D");

  // Array of pointer of the AliMUONHitMapA1:
  //  two HitMaps per chamber, or one HitMap per cahtode plane
  fHitMap= new AliMUONHitMapA1* [2*AliMUONConstants::NCh()];

  //Loop over chambers for the definition AliMUONHitMap
  for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {
    chamber = &(pMUON->Chamber(i));
    c1Segmentation = chamber->SegmentationModel(1); // Cathode plane 1
    fHitMap[i] = new AliMUONHitMapA1(c1Segmentation, fTDList);
    c2Segmentation = chamber->SegmentationModel(2); // Cathode plane 2
    fHitMap[i+AliMUONConstants::NCh()] = new AliMUONHitMapA1(c2Segmentation, fTDList);
  }

// Loop over files to merge and to digitize
    fSignal = kTRUE;
    for (Int_t inputFile=0; inputFile<fManager->GetNinputs(); inputFile++) {
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::Digitize() Input File is "<<inputFile<<endl;


      // Connect MUON Hit branch
      if (inputFile > 0 ) {
	fSignal = kFALSE;
         runloader    = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
         if (runloader == 0x0) {
           cerr<<"AliMUONDigitizerv1::Digitize() RunLoader for inputFile "<<inputFile<< " not found !!! "<<endl;
         }
         gime  = runloader->GetLoader("MUONLoader");
         if (gime->TreeH() == 0x0) gime->LoadHits("READ");	
	 muondata->SetLoader(gime);
	 muondata->SetTreeAddress("H");
      }

      // Setting the address 
      TTree *treeH = gime->TreeH();
      if (treeH == 0x0) {
	Error("Digitize","Can not get TreeH from input %d",inputFile);
	Info("Digitize","Now treeH is %#x. MUONLoader is %#x",gime->TreeH(),gime);
	return;
      }
      if (GetDebug()>2) {
	cerr<<"AliMUONDigitizerv1::Exec inputFile is "<<inputFile<<" "<<endl;
	cerr<<"AliMUONDigitizerv1::Exec treeH" << treeH <<endl;
      }
	
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::Exec Setting tree addresses"<<endl;

      fMask = fManager->GetMask(inputFile);
      //
      // Loop over tracks
      Int_t itrack;
      Int_t ntracks = (Int_t) treeH->GetEntries();
      for (itrack = 0; itrack < ntracks; itrack++) {
	if (GetDebug()>2) cerr<<"AliMUONDigitizerv1::Exec itrack = "<<itrack<<endl;
	muondata->ResetHits();
	treeH->GetEvent(itrack);
	//
	//  Loop over hits
	Int_t ihit, ichamber;
	AliMUONHit* mHit;
	TClonesArray* hits = muondata->Hits();
	for(ihit = 0; ihit < hits->GetEntriesFast(); ihit++) {
	  mHit = static_cast<AliMUONHit*>(hits->At(ihit));
	  ichamber = mHit->Chamber()-1;  // chamber number
	  if (ichamber > AliMUONConstants::NCh()-1) {
	    cerr<<"AliMUONDigitizer: ERROR: "
		<<"fNch > AliMUONConstants::NCh()-1, fNch, NCh(): "
		<<ichamber<<", "<< AliMUONConstants::NCh()<<endl;
	    return;
	  }
	  chamber = &(pMUON->Chamber(ichamber));
	  //
	  //Dumping Hit content:
	  if (GetDebug()>2) {
	    cerr<<"AliMuonDigitizerv1::Exec ihit, ichamber, x, y, z, eloss " <<
	      ihit << " " << 
	      mHit->Chamber() << " " <<
	      mHit->X() << " " <<
	      mHit->Y() << " " <<
	      mHit->Z() << " " <<
	      mHit->Eloss() << " " << endl;
	  }
	  // 
	  // Inititializing Correlation
	  chamber->ChargeCorrelationInit();
	  if (ichamber < AliMUONConstants::NTrackingCh()) {
	    // Tracking Chamber
	    // Initialize hit position (cursor) in the segmentation model 
	    chamber->SigGenInit(mHit->X(), mHit->Y(), mHit->Z());
	  } else {
	    // Trigger Chamber 
	  }
	  MakeTransientDigit(itrack, ihit, mHit);
	} // hit loop
      } // track loop
    } // end file loop
    if (GetDebug()>2) cerr<<"AliMUONDigitizer::Exec End of hits, track and file loops"<<endl;

    // Loop on cathodes
    Int_t icat;
    for(icat=0; icat<2; icat++) {
      //
      // Filling Digit List
      Int_t tracks[kMAXTRACKS];
      Int_t charges[kMAXTRACKS];
      Int_t nentries = fTDList->GetEntriesFast();
      Int_t digits[6];
      for (Int_t nent = 0; nent < nentries; nent++) {
	AliMUONTransientDigit *address = (AliMUONTransientDigit*)fTDList->At(nent);
	if (address == 0) continue; 
	Int_t ich = address->Chamber();
	Int_t   q = address->Signal(); 
	chamber = &(pMUON->Chamber(ich));
	//
	//  Digit Response (noise, threshold, saturation, ...)
	AliMUONResponse * response = chamber->ResponseModel();
	q = response->DigitResponse(q,address);
	
	if (!q) continue;
	
	digits[0] = address->PadX();
	digits[1] = address->PadY();
	digits[2] = address->Cathode()-1;
	digits[3] = q;
	digits[4] = address->Physics();
	digits[5] = address->Hit();
	
	Int_t nptracks = address->GetNTracks();
	
	if (nptracks > kMAXTRACKS) {
	  if (GetDebug() >0) {
	    cerr<<"AliMUONDigitizer:Exec  nptracks > 10 "<<nptracks;
	    cerr<<"reset to max value "<<kMAXTRACKS<<endl;
	  }
	  nptracks = kMAXTRACKS;
	}
	if (nptracks > 2 && GetDebug() >2) {
	  cerr<<"AliMUONDigitizer::Exec  nptracks > 2 "<<nptracks<<endl;
	  //	printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,digits[0],digits[1],q);
	}
	for (Int_t tr = 0; tr < nptracks; tr++) {
	  tracks[tr]   = address->GetTrack(tr);
	  charges[tr]  = address->GetCharge(tr);
	}      //end loop over list of tracks for one pad
	// Sort list of tracks according to charge
	if (nptracks > 1) {
	  SortTracks(tracks,charges,nptracks);
	}
	if (nptracks < kMAXTRACKS ) {
	  for (Int_t i = nptracks; i < kMAXTRACKS; i++) {
	    tracks[i]  = -1;
	    charges[i] = 0;
	  }
	}
	
	// Add digits
	if (GetDebug()>3) cerr<<"AliMUONDigitzerv1::Exex TransientDigit to Digit"<<endl;
	if ( digits[2] == icat ) muondata->AddDigit(ich,tracks,charges,digits);
//	printf("test rm ich %d padX %d padY %d \n",ich, digits[0], digits[1]);
      }
      // Filling list of digits per chamber for a given cathode.
      muondata->Fill("D");
      muondata->ResetDigits();    
    } // end loop cathode
    fTDList->Delete();  
    
    for(Int_t ii = 0; ii < 2*AliMUONConstants::NCh(); ++ii) {
      if (fHitMap[ii]) {
	delete fHitMap[ii];
	fHitMap[ii] = 0;
      }
    }
    
    if (GetDebug()>2) 
      cerr<<"AliMUONDigitizer::Exec: writing the TreeD: "
	  <<treeD->GetName()<<endl;

    gime->WriteDigits("OVERWRITE");
    delete [] fHitMap;
    delete fTDList;
    muondata->ResetHits();
    gime->UnloadHits();
    gime->UnloadDigits();
}
//------------------------------------------------------------------------
void AliMUONDigitizerv1::SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr)
{
  //
  // Sort the list of tracks contributing to a given digit
  // Only the 3 most significant tracks are acctually sorted
  //
  
  //
  //  Loop over signals, only 3 times
  //
  
  Int_t qmax;
  Int_t jmax;
  Int_t idx[3] = {-2,-2,-2};
  Int_t jch[3] = {-2,-2,-2};
  Int_t jtr[3] = {-2,-2,-2};
  Int_t i,j,imax;
  
  if (ntr<3) imax=ntr;
  else imax=3;
  for(i=0;i<imax;i++){
    qmax=0;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == idx[i-1]) 
	 ||(i == 2 && (j == idx[i-1] || j == idx[i-2]))) continue;
      
      if(charges[j] > qmax) {
	qmax = charges[j];
	jmax=j;
      }       
    } 
    
    if(qmax > 0) {
      idx[i]=jmax;
      jch[i]=charges[jmax]; 
      jtr[i]=tracks[jmax];
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -2) {
      charges[i]=0;
      tracks[i]=0;
    } else {
      charges[i]=jch[i];
      tracks[i]=jtr[i];
    }
  }
}
