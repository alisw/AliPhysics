
//  Do the Digitization (Digit) from summable Digits (SDigit)
//  Allow the merging of signal file with background file(s).

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
#include "AliMUONDigitizerv2.h"
#include "AliMUONHit.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONPadHit.h"
#include "AliMUONTransientDigit.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

ClassImp(AliMUONDigitizerv2)

//___________________________________________
AliMUONDigitizerv2::AliMUONDigitizerv2() :
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
    cerr<<"AliMUONDigitizerv2::AliMUONDigitizerv2"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
}

//___________________________________________
AliMUONDigitizerv2::AliMUONDigitizerv2(AliRunDigitizer* manager):
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
    cerr<<"AliMUONDigitizerv2::AliMUONDigitizerv2"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
}

//------------------------------------------------------------------------
AliMUONDigitizerv2::~AliMUONDigitizerv2()
{
// Destructor
}

//------------------------------------------------------------------------
void AliMUONDigitizerv2::AddTransientDigit(AliMUONTransientDigit * mTD)
{
  // Choosing the maping of the cathode plane of the chamber:
  Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
  fTDList->AddAtAndExpand(mTD, fTDCounter);
  fHitMap[iNchCpl]->SetHit( mTD->PadX(), mTD->PadY(), fTDCounter);
  fTDCounter++;
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizerv2::ExistTransientDigit(AliMUONTransientDigit * mTD)
{
  // Choosing the maping of the cathode plane of the chamber:
  Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
  return( fHitMap[iNchCpl]->TestHit(mTD->PadX(), mTD->PadY()) );
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizerv2::Init()
{

// Initialization 
  if (GetDebug()>2) Info("Init","AliMUONDigitizerv2::Init() starts");

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
  gime->LoadSDigits("READ");
  gime->LoadDigits("RECREATE");
  
  return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONDigitizerv2::UpdateTransientDigit(Int_t /*track*/, AliMUONTransientDigit * mTD)
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
  // List of tracks and trackcharge
  Int_t itrack;
  for(itrack=0;itrack<kMAXTRACKS;itrack++) {
    pdigit->UpdateTrackList(mTD->Track(itrack), mTD->TrackCharge(itrack));
  }
}

//--------------------------------------------------------------------------
void AliMUONDigitizerv2::MakeTransientDigitFromSDigit(Int_t iChamber, AliMUONDigit * sDigit)
{
  Int_t digits[6];
  

  //
  // Creating a new TransientDigits from SDigit
  digits[0] = sDigit->PadX();  // Padx of the Digit
  digits[1] = sDigit->PadY();  // Pady of the Digit
  digits[2] = sDigit->Cathode()+1;  // Cathode plane
  digits[3] = sDigit->Signal();  // Induced charge in the Pad
  if (fSignal) digits[4] = sDigit->Signal();
  else         digits[4] = 0;
  digits[5] = sDigit->Hit();    // Hit number in the list
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::MakeTransientDigitFromSDigit " <<
		      "PadX "<< digits[0] << " " <<
		      "PadY "<< digits[1] << " " <<
		      "Plane " << digits[2] << " " <<
		      "Charge " << digits[3] <<" " <<
		      "Hit " << digits[5] << endl;   
		 
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::MakeTransientDigitFromSDigit Creating AliMUONTransientDigit"<<endl;
  AliMUONTransientDigit * mTD = new AliMUONTransientDigit(iChamber, digits);

  // List of tracks and trackcharge
  Int_t itrack;
  for(itrack=0;itrack<kMAXTRACKS;itrack++) {
    mTD->AddToTrackList(fMask+sDigit->Track(itrack), sDigit->TrackCharge(itrack));
  }

  if (!ExistTransientDigit(mTD)) {
    AddTransientDigit(mTD);
    if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::MakeTransientDigitFromSDigit Adding TransientDigit"<<endl;
  }
  else {
    if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::MakeTransientDigitFromSDigit updating TransientDigit"<<endl;
    UpdateTransientDigit(0, mTD);
    delete mTD;
  }
}

//-----------------------------------------------------------------------
void AliMUONDigitizerv2::Exec(Option_t* option)
{
  TString optionString = option;
  if (optionString.Data() == "deb") {
    Info("Digitize","Called with option deb ");
    fDebug = 3;
  }

  AliRun *aliRun;
  AliMUONChamber*   chamber;
  AliSegmentation*  c1Segmentation; //Cathode plane c1 of the chamber
  AliSegmentation*  c2Segmentation; //Cathode place c2 of the chamber

  if (GetDebug()>2) Info("Exec","AliMUONDigitizerv2::Exec() starts");
  fTDList = new TObjArray;

  //Loaders (We assume input0 to be the output too)
  AliRunLoader *runloader, *runloaderOut;    // Input /  Output loaders
  AliLoader    *gime, *gimeOut; 
 
  // Getting runloader
  runloader    = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  if (runloader == 0x0) {
    Error("Exec","RunLoader is not in input file 0");
    return;  // RunDigitizer is not working.  
  }
  // Getting MUONloader
  gime        = runloader->GetLoader("MUONLoader");
  if (gime->TreeS()==0x0) {
      if (GetDebug()>2) Info("Exec","TreeS is not loaded yet. Loading...");
     gime->LoadSDigits("READ");
       if (GetDebug()>2) Info("Exec","Now treeS is %#x. MUONLoader is %#x",gime->TreeS(),gime);
  }

  if (GetDebug()>2) Info("Exec","Loaders ready");

  if (runloader->GetAliRun() == 0x0) runloader->LoadgAlice();
  aliRun = runloader->GetAliRun();

  // Getting Module MUON  
  AliMUON *pMUON  = (AliMUON *) aliRun->GetDetector("MUON");
  if (!pMUON) {
    Error("Digitize","Module MUON not found in the input file");
    return;
  }
  // Getting Muon data
  AliMUONData * muondata = pMUON->GetMUONData(); 

  // Loading Event
  Int_t currentevent = fManager->GetOutputEventNr();
  
  if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::Exec() Event Number is "<<currentevent <<endl;
  if ( (currentevent<10)                                                 || 
       (Int_t(TMath::Log10(currentevent)) == TMath::Log10(currentevent) ) )
    cout <<"ALiMUONDigitizerv2::Exec() Event Number is "<< currentevent <<endl;

  // Output runloader 
  runloaderOut  = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  gimeOut = runloaderOut->GetLoader("MUONLoader");
  

  if (gimeOut->TreeD() == 0x0) {
    gimeOut->MakeDigitsContainer();
  }
  TTree* treeD = gimeOut->TreeD();
  muondata->SetLoader(gimeOut);
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
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::Exec() Input File is "<<inputFile<<endl;
      if (inputFile == 0) {
	 muondata->SetLoader(gime);
	 muondata->SetTreeAddress("S");	
      } else {
	// Connect MUON Hit branch	
	fSignal = kFALSE;
	runloader    = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
	if (runloader == 0x0) {
	  cerr<<"AliMUONDigitizerv2::Digitize() RunLoader for inputFile "<<inputFile<< " not found !!! "<<endl;
	}
	gime  = runloader->GetLoader("MUONLoader");
	if (gime->TreeS() == 0x0) gime->LoadSDigits("READ");	
	muondata->SetLoader(gime);
	muondata->SetTreeAddress("S");
      }
      
      // Setting the address 
      TTree *treeS = gime->TreeS();
      if (treeS == 0x0) {
	Error("Digitize","Can not get TreeS from input %d",inputFile);
	Info("Digitize","Now treeS is %#x. MUONLoader is %#x",gime->TreeS(),gime);
	return;
      }
      if (GetDebug()>2) {
	cerr<<"AliMUONDigitizerv2::Exec inputFile is "<<inputFile<<" "<<endl;
	cerr<<"AliMUONDigitizerv2::Exec treeS" << treeS <<endl;
      }
      
      if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::Exec Setting tree addresses"<<endl;
      
      fMask = fManager->GetMask(inputFile);
      //
      // Loop over SDigits
      Int_t ndig,k;
      AliMUONDigit *sDigit;
      TClonesArray * muonSDigits;
      for (Int_t ich = 0; ich < AliMUONConstants::NCh(); ich++) {  //  loop over chamber
	muondata->ResetSDigits();
	muondata->GetCathodeS(0);
	//TClonesArray *
	muonSDigits = muondata->SDigits(ich); 
	ndig = muonSDigits->GetEntriesFast();
	// 	printf("\n 1 Found %d Sdigits in %p chamber %d", ndig, muonSDigits,ich);
	Int_t n = 0;
	for (k = 0; k < ndig; k++) {
	  sDigit = (AliMUONDigit*) muonSDigits->UncheckedAt(k);
	  MakeTransientDigitFromSDigit(ich,sDigit);
	}
	muondata->ResetSDigits();
	muondata->GetCathodeS(1);
	//TClonesArray *
	muonSDigits = muondata->SDigits(ich); 
	ndig=muonSDigits->GetEntriesFast();
	// 	printf("\n 2 Found %d Sdigits in %p chamber %d", ndig, muonSDigits,ich);
	n = 0;
	for (k = 0; k < ndig; k++) {
	  sDigit = (AliMUONDigit*) muonSDigits->UncheckedAt(k);
	  MakeTransientDigitFromSDigit(ich,sDigit);
	}
	
      } // SDigits loop end loop over chamber
    } // end file loop
    if (GetDebug()>2) cerr<<"AliMUONDigitizerv2::Exec End of Sdigits, file loops"<<endl;
    
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
	if (GetDebug()>3) cerr<<"AliMUONDigitzerv2::Exec TransientDigit to Digit"<<endl;
	if ( digits[2] == icat ) muondata->AddDigit(ich,tracks,charges,digits);
// 	printf("test rm ich %d padX %d padY %d \n",ich, digits[0], digits[1]);
      }
      muondata->SetLoader(gimeOut);
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
    
    gimeOut->WriteDigits("OVERWRITE");
    delete [] fHitMap;
    delete fTDList;
    muondata->ResetSDigits();
    gime->UnloadSDigits();
    gimeOut->UnloadDigits();
}
//------------------------------------------------------------------------
void AliMUONDigitizerv2::SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr)
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
