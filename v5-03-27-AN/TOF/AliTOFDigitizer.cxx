/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________//
//                                                                         //
// This is a TTask that makes TOF-Digits out of TOF-SDigits.               //
// The simulation of the detector is performed at sdigits level:           //
// during digitization the unique task is the sum of all sdigits in the    //
// same pad.                                                               //
// Digits are written to TreeD in branch "TOF".                            //
//                                                                         //
// -- Author :  F. Pierella (Bologna University) pierella@bo.infn.it       //
//                                                                         //
//_________________________________________________________________________//

//#include "Riostream.h"

//#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
#include "TRandom.h"
#include "TObjArray.h"

#include "AliLoader.h"
#include "AliLog.h"
#include "AliDigitizationInput.h"
#include "AliRunLoader.h"
#include "AliRun.h"

#include "AliTOFcalib.h"
//#include "AliTOFChannelOnlineArray.h"
//#include "AliTOFChannelOnlineStatusArray.h"
//#include "AliTOFChannelOffline.h"
#include "AliTOFDigitizer.h"
#include "AliTOFdigit.h"
#include "AliTOFHitMap.h"
#include "AliTOFGeometry.h"
#include "AliTOFSDigit.h"
#include "AliTOF.h"

extern TRandom *gRandom;

extern AliRun *gAlice;


ClassImp(AliTOFDigitizer)

//___________________________________________
  AliTOFDigitizer::AliTOFDigitizer()  :
    AliDigitizer(),
    fDigits(new TClonesArray("AliTOFdigit",4000)),
    fSDigitsArray(new TClonesArray("AliTOFSDigit",1000)),
  fhitMap(0x0),
  fCalib(new AliTOFcalib())
{
  // Default ctor - don't use it
  InitDecalibration();
}

//___________________________________________
AliTOFDigitizer::AliTOFDigitizer(AliDigitizationInput* digInput): 
  AliDigitizer(digInput), 
  fDigits(new TClonesArray("AliTOFdigit",4000)),
  fSDigitsArray(new TClonesArray("AliTOFSDigit",1000)),
  fhitMap(0x0),
  fCalib(new AliTOFcalib())
{
  //ctor with RunDigitizer
  InitDecalibration();
}

//------------------------------------------------------------------------
AliTOFDigitizer::AliTOFDigitizer(const AliTOFDigitizer &source):
  AliDigitizer(source),
  fDigits(source.fDigits),
  fSDigitsArray(source.fSDigitsArray),
  fhitMap(source.fhitMap),
  fCalib(source.fCalib)
{
  // copy constructor
}

//------------------------------------------------------------------------
  AliTOFDigitizer& AliTOFDigitizer::operator=(const AliTOFDigitizer &source)
{
  // ass. op.
  
  if (this == &source)
    return *this;

  AliDigitizer::operator=(source);
  fDigits=source.fDigits;
  fSDigitsArray=source.fSDigitsArray;
  fhitMap=source.fhitMap;
  fCalib=source.fCalib;
  return *this;

}

//------------------------------------------------------------------------
AliTOFDigitizer::~AliTOFDigitizer()
{
  // Destructor
  delete fCalib;
  if (fDigits){
    fDigits->Delete();
    delete fDigits;
    fDigits=0x0;
  }
  if (fSDigitsArray){
    fSDigitsArray->Delete();
    delete fSDigitsArray;
    fSDigitsArray=0x0;
  }
}

//---------------------------------------------------------------------

void AliTOFDigitizer::Digitize(Option_t* /*option*/)
{
  //
  // Perform digitization and merging.
  // The algorithm is the following:
  // - a hitmap is created to check if a pad is already activated;
  // - an sdigits container is created to collect all sdigits from
  //   different files;
  // - sdigits are summed using the hitmap;
  // - the sdigits container is used to create the array of AliTOFdigit.
  //

  AliDebug(1, "");


  // get the ptr to TOF detector
  AliTOF * tof = (AliTOF *) gAlice->GetDetector("TOF") ;

  //Make branches

  const Int_t kSize = 20;
  char branchname[kSize];
  snprintf(branchname,kSize,"%s", tof->GetName ());
 
  AliRunLoader* outrl = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
  if (outrl == 0x0)
   {
     AliError("Can not find Run Loader in output folder.");
     return;
   }
   
  AliLoader* outgime = outrl->GetLoader("TOFLoader");
  if (outgime == 0x0)
   {
     AliError("Can not get TOF Loader from Output Run Loader.");
     return;
   }
  
  TTree* treeD = outgime->TreeD();
  if (treeD == 0x0)
   {
     outgime->MakeTree("D");
     treeD = outgime->TreeD();
   }
  //Make branch for digits (to be created in Init())
  tof->MakeBranchInTree(treeD,branchname,&fDigits,4000);

  // container for all summed sdigits (to be created in Init())
  //fSDigitsArray=new TClonesArray("AliTOFSDigit",1000);
  
  // create hit map (to be created in Init())
  fhitMap = new AliTOFHitMap(fSDigitsArray);
  
  // Loop over files to digitize

  for (Int_t inputFile=0; inputFile<fDigInput->GetNinputs();
       inputFile++) {
    ReadSDigit(inputFile);
   }

  // create digits
  CreateDigits();

  // free used memory for Hit Map in current event
  delete fhitMap;
  fSDigitsArray->Clear();

  treeD->Fill();

  AliDebug(2,"----------------------------------------");
  AliDebug(1,Form("%d digits have been created", fDigits->GetEntriesFast()));
  AliDebug(2,"----------------------------------------");

  outgime->WriteDigits("OVERWRITE");
  outgime->UnloadDigits();
  fDigits->Clear();

}

//---------------------------------------------------------------------

void AliTOFDigitizer::CreateDigits()
{
  // loop on sdigits container to fill the AliTOFdigit TClonesArray
  // start digitizing all the collected sdigits 

  Int_t ndump=0; // dump the first ndump created digits for each event

  // get the total number of collected sdigits
  Int_t ndig = fSDigitsArray->GetEntriesFast();

  Int_t  vol[5]={-1,-1,-1,-1,-1};  // location for a digit
  Int_t  digit[4] = {0,0,0,0};     // TOF digit variables
  Int_t tracknum[AliTOFSDigit::kMAXDIGITS]; // contributing tracks for the current slot
  for (Int_t aa=0; aa<AliTOFSDigit::kMAXDIGITS; aa++) tracknum[aa] = -1;

  for (Int_t k = 0; k < ndig; k++) {
    
    for (Int_t i=0; i<5; i++) vol[i] = -1;
    
    // Get the information for this digit
    AliTOFSDigit *tofsdigit = (AliTOFSDigit *) fSDigitsArray->UncheckedAt(k);
    
    Int_t nslot=tofsdigit->GetNDigits(); // get the number of slots
    // for current sdigit
    
    // TOF sdigit volumes (always the same for all slots)
    Int_t sector    = tofsdigit->GetSector(); // range [0-17]
    Int_t plate     = tofsdigit->GetPlate();  // range [0- 4]
    Int_t strip     = tofsdigit->GetStrip();  // range [0-14/18/19]
    Int_t padz      = tofsdigit->GetPadz();   // range [0- 1]
    Int_t padx      = tofsdigit->GetPadx();   // range [0-47]
    
    vol[0] = sector;
    vol[1] = plate;
    vol[2] = strip;
    vol[3] = padx;
    vol[4] = padz;
    
    //--------------------- QA section ----------------------
    // in the while, I perform QA
    Bool_t isSDigitBad = (sector<0 || sector>17 || plate<0 || plate >4 || padz<0 || padz>1 || padx<0 || padx>47);
    
    if (isSDigitBad)
      AliFatal(Form("strange sdigit found   %2d  %1d  %2d  %1d %2d", sector, plate, strip, padz, padx));
    //-------------------------------------------------------
    
    //------------------- Dump section ----------------------
    if (k<ndump) {
      AliInfo(Form("%2d-th digit: Sector %2d | Plate %1d | Strip %2d | PadZ %1d | PadX %2d ", k, sector, plate, strip, padz, padx));
      AliInfo("----------------------------------------------------");
    }
    // ------------------------------------------------------
    
    // start loop on number of slots for current sdigit
    for (Int_t islot = 0; islot < nslot; islot++) {
      for (Int_t aa=0; aa<4; aa++) digit[aa] = 0; // TOF digit variables
      for (Int_t aa=0; aa<AliTOFSDigit::kMAXDIGITS; aa++) tracknum[aa] = -1;
      
      Int_t tdc=tofsdigit->GetTdc(islot); digit[0]=tdc;
      Int_t adc=tofsdigit->GetAdc(islot); digit[1]=adc;

      //if (tdc>=8192) continue;//AdC

      tracknum[0]=tofsdigit->GetTrack(islot,0);
      tracknum[1]=tofsdigit->GetTrack(islot,1);
      tracknum[2]=tofsdigit->GetTrack(islot,2);
      
      // new with placement must be used
      // adding a TOF digit for each slot
      TClonesArray &aDigits = *fDigits;
      Int_t last=fDigits->GetEntriesFast();
      new (aDigits[last]) AliTOFdigit(tracknum, vol, digit);

    }
    
  } // end loop on sdigits - end digitizing all collected sdigits

  //Insert Decalibration 
  AliDebug(2,"in digitizer, create digits");
  DecalibrateTOFSignal();
}

//---------------------------------------------------------------------

void AliTOFDigitizer::ReadSDigit(Int_t inputFile )
{
  // Read sdigits for current event and inputFile; 
  // store them into the sdigits container
  // and update the hit map
  // SDigits from different files are assumed to
  // be created with the same simulation parameters.
  
  // creating the TClonesArray to store the digits
  static TClonesArray sdigitsClonesArray("AliTOFSDigit",  1000); 
  sdigitsClonesArray.Clear();

  // get the treeS from digInput
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inputFile));
  if (rl == 0x0)
   {
     AliError(Form("Can not find Run Loader in input %d folder.",inputFile));
     return;
   }

  AliLoader* gime = rl->GetLoader("TOFLoader");
  if (gime == 0x0)
   {
     AliError(Form("Can not get TOF Loader from Input %d Run Loader.",inputFile));
     return;
   }

  TTree* currentTreeS=gime->TreeS();
  if (currentTreeS == 0x0)
   {
     Int_t retval = gime->LoadSDigits();
     if (retval) 
      {
         AliError(Form("Error occured while loading S. Digits for Input %d",inputFile));
         return;
      }
     currentTreeS=gime->TreeS();
     if (currentTreeS == 0x0)
      {
         AliError(Form("Can not get S. Digits Tree for Input %d",inputFile));
         return;
      }
   } 
  // get the branch TOF inside the treeS
  TClonesArray * sdigitsDummyContainer=&sdigitsClonesArray;
  // check if the branch exist
  TBranch* tofBranch=currentTreeS->GetBranch("TOF");

  if(!tofBranch){
    AliFatal(Form("TOF branch not found for input %d",inputFile));
    return;
  }
  
  tofBranch->SetAddress(&sdigitsDummyContainer);           
  
  Int_t nEntries = (Int_t)tofBranch->GetEntries();                                

  // Loop through all entries in the tree
  Int_t nbytes = 0;
  
  Int_t  vol[5]; // location for a sdigit
  for (Int_t i=0; i<5; i++) vol[i] = -1;

  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
    
    // Import the tree
    nbytes += tofBranch->GetEvent(iEntry);
    
    // Get the number of sdigits
    Int_t ndig = sdigitsDummyContainer->GetEntriesFast();
    
    for (Int_t k=0; k<ndig; k++) {
      AliTOFSDigit *tofSdigit= (AliTOFSDigit*) sdigitsDummyContainer->UncheckedAt(k);
      
      for (Int_t i=0; i<5; i++) vol[i] = -1;

      // check the sdigit volume
      vol[0] = tofSdigit->GetSector();
      vol[1] = tofSdigit->GetPlate();
      vol[2] = tofSdigit->GetStrip();
      vol[3] = tofSdigit->GetPadx();
      vol[4] = tofSdigit->GetPadz();
      
      if (fhitMap->TestHit(vol) != kEmpty) {
	AliTOFSDigit *sdig = static_cast<AliTOFSDigit*>(fhitMap->GetHit(vol));
	sdig->Update(tofSdigit);

      } else {

	CollectSDigit(tofSdigit); // collect the current sdigit
	fhitMap->SetHit(vol);     // update the hitmap for location vol

      } // if (hitMap->TestHit(vol) != kEmpty)
      
    } // for (Int_t k=0; k<ndig; k++)

  } // end loop on entries

}


//_____________________________________________________________________________
void AliTOFDigitizer::CollectSDigit(const AliTOFSDigit * const sdigit)
{
  //
  // Add a TOF sdigit in container
  // new with placement must be used
  TClonesArray &aSDigitsArray = *fSDigitsArray;
  Int_t last=fSDigitsArray->GetEntriesFast();
  // make a copy of the current sdigit and
  // put it into tmp array
  new (aSDigitsArray[last]) AliTOFSDigit(*sdigit);
}

//_____________________________________________________________________________
void AliTOFDigitizer::InitDecalibration() const {
  //
  // Initialize TOF digits decalibration
  //

  fCalib->Init();
  /*
  fCalib->CreateCalArrays();
  fCalib->ReadSimHistoFromCDB("TOF/Calib", -1); // use AliCDBManager's number
  fCalib->ReadParOfflineFromCDB("TOF/Calib", -1); // use AliCDBManager's number
  */
}
//---------------------------------------------------------------------
void AliTOFDigitizer::DecalibrateTOFSignal() {
  //
  // Decalibrate TOF signals according to OCDB parameters
  //

  Double_t time=0., tot=0., corr=0.;
  Int_t deltaBC=0, l0l1=0, tdcBin=0;
  Int_t index = -1;
  Int_t detId[5] ={-1,-1,-1,-1,-1};
  UInt_t timestamp=0;

  Int_t ndigits = fDigits->GetEntriesFast();
  // Loop on TOF Digits
  for (Int_t i=0;i<ndigits;i++){
    AliTOFdigit * dig = (AliTOFdigit*)fDigits->At(i);
    detId[0] = dig->GetSector();
    detId[1] = dig->GetPlate();
    detId[2] = dig->GetStrip();
    detId[3] = dig->GetPadz();
    detId[4] = dig->GetPadx();
    dig->SetTdcND(dig->GetTdc()); // save the non decalibrated time

    index = AliTOFGeometry::GetIndex(detId); // The channel index    

    // Read Calibration parameters from the CDB
    // get digit info
    time = dig->GetTdc() * AliTOFGeometry::TdcBinWidth(); /* ps */
    tot = dig->GetToT() * AliTOFGeometry::ToTBinWidth() * 1.e-3; /* ns */
    deltaBC = 0;//dig->GetDeltaBC();
    l0l1 = 0;//dig->GetL0L1Latency();

    // get correction
    corr = fCalib->GetTimeCorrection(index, tot, deltaBC, l0l1, timestamp); /* ps */
    AliDebug(2, Form("calibrate index %d: time=%f (ps) tot=%f (ns) deltaBC=%d l0l1=%d timestamp=%d corr=%f (ps)",
		     index, time, tot, deltaBC, l0l1, timestamp, corr));

    // apply time correction
    time += corr;

    // convert in TDC bins and set digit
    //tdcBin = (Int_t)(time / AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
    tdcBin = TMath::Nint(time / AliTOFGeometry::TdcBinWidth()); //the corrected time (tdc counts)
    dig->SetTdc(tdcBin);

  }

  AliDebug(1,"Simulating miscalibrated digits");

  return;
}

//---------------------------------------------------------------------
/*
void AliTOFDigitizer::DecalibrateTOFSignal(){ // Old implementation

  // Read Calibration parameters from the CDB

  TObjArray * calOffline= fCalib->GetTOFCalArrayOffline();

  AliDebug(2,Form("Size of array for Offline Calibration = %i",calOffline->GetEntries()));

  // Initialize Quantities to Simulate ToT Spectra

  TH1F * hToT= fCalib->GetTOFSimToT();
  Int_t nbins = hToT->GetNbinsX();
  Float_t delta = hToT->GetBinWidth(1);
  Float_t maxch = hToT->GetBinLowEdge(nbins)+delta;
  Float_t minch = hToT->GetBinLowEdge(1);
  Float_t max=0,min=0; //maximum and minimum value of the distribution
  Int_t maxbin=0,minbin=0; //maximum and minimum bin of the distribution

  for (Int_t ii=nbins; ii>0; ii--){
    if (hToT->GetBinContent(ii)!= 0) {
      max = maxch - (nbins-ii-1)*delta;
      maxbin = ii; 
      break;}
  }
  for (Int_t j=1; j<nbins; j++){
    if (hToT->GetBinContent(j)!= 0) {
      min = minch + (j-1)*delta;
      minbin = j; 
      break;}
  }

  Float_t maxToT=max;
  Float_t minToT=min;
 
  Float_t maxToTDistr=hToT->GetMaximum();

  AliDebug (1, Form(" The minimum ToT = %f", minToT)); 
  AliDebug (1, Form(" The maximum ToT = %f", maxToT)); 
  AliDebug (1, Form(" The maximum peak in ToT = %f", maxToTDistr)); 
  
  // Loop on TOF Digits

  Bool_t isToTSimulated=kFALSE;
  Bool_t misCalibPars=kFALSE;
  if(hToT->GetEntries()>0)isToTSimulated=kTRUE;  
  Int_t ndigits = fDigits->GetEntriesFast();
  for (Int_t i=0;i<ndigits;i++){
    AliTOFdigit * dig = (AliTOFdigit*)fDigits->At(i);
    Int_t detId[5];
    detId[0] = dig->GetSector();
    detId[1] = dig->GetPlate();
    detId[2] = dig->GetStrip();
    detId[3] = dig->GetPadz();
    detId[4] = dig->GetPadx();
    dig->SetTdcND(dig->GetTdc()); // save the non decalibrated time
    if(isToTSimulated){  

      //A realistic ToT Spectrum was found in input, 
      //decalibrated TOF Digits likely to be simulated....
 
      Int_t index = AliTOFGeometry::GetIndex(detId); // The channel index    
      AliTOFChannelOffline *calChannelOffline = (AliTOFChannelOffline *)calOffline->At(index); //retrieve the info for time slewing 
      Double_t par[6];  // time slewing parameters
  
      //check whether we actually ask for miscalibration

      for (Int_t j = 0; j<6; j++){
	par[j]=(Double_t)calChannelOffline->GetSlewPar(j);
	if(par[j]!=0)misCalibPars=kTRUE;
      }
      AliDebug(2,Form(" Calib Pars = %f (0-th parameter for time slewing + time delay), %f, %f, %f, %f, %f ",par[0],par[1],par[2],par[3],par[4],par[5]));

      // Now generate Realistic ToT distribution from TestBeam Data. 
      // Tot is in ns, assuming a Matching Window of 10 ns.

      Float_t simToT = 0;
      Float_t trix = 0;
      Float_t triy = 0;
      Double_t timeCorr;
      Double_t tToT;
      while (simToT <= triy){
	trix = gRandom->Rndm(i);
	triy = gRandom->Rndm(i);
	trix = (maxToT-minToT)*trix + minToT; 
	triy = maxToTDistr*triy;
	Int_t binx=hToT->FindBin(trix);
	simToT=hToT->GetBinContent(binx);
      }
      // the generated ToT (ns)
      tToT= (Double_t) trix; // to apply slewing we start from ns..
      // transform TOF signal in ns
      AliDebug(2,Form(" The Initial Time (counts): %i: ",dig->GetTdc()));
      AliDebug(2,Form(" Time before miscalibration (ps) %e: ",dig->GetTdc()*(Double_t)AliTOFGeometry::TdcBinWidth()));
      // add slewing effect
      timeCorr=par[0] + tToT*(par[1] +tToT*(par[2] +tToT*(par[3] +tToT*(par[4] +tToT*par[5])))); 
      AliDebug(2,Form(" The Time slewing + delay (ns): %f: ",timeCorr));
      // add global time shift
      //convert to ps
      timeCorr*=1E3;
      Double_t timeMis = (Double_t)(dig->GetTdc())*(Double_t)AliTOFGeometry::TdcBinWidth();
      timeMis = timeMis+timeCorr;
      AliDebug(2,Form(" The Miscalibrated time (ps): %e: ",timeMis));

      // now update the digit info
 
      Int_t tdcCorr= (Int_t)(timeMis/AliTOFGeometry::TdcBinWidth());
      AliDebug(2,Form(" Final Time (counts): %i: ",tdcCorr));
      // Setting Decalibrated Time signal (TDC counts)    
      dig->SetTdc(tdcCorr);   
      // Setting realistic ToT signal (TDC counts) 
      tToT*=1E3; //back to ps  
      Int_t tot=(Int_t)(tToT/AliTOFGeometry::ToTBinWidth());//(factor 1E3 as input ToT is in ns)
      dig->SetToT(tot); 
      AliDebug(2,Form(" Final Time and ToT (counts): %d: , %d:",dig->GetTdc(),dig->GetToT()));
      if(tdcCorr<0){
	AliWarning (Form(" The bad Slewed Time(TDC counts)= %d ", tdcCorr)); 
	AliWarning(Form(" The bad ToT (TDC counts)= %d ", tot)); 
      }
    }
    else{
    // For Data with no Miscalibration, set ToT signal == Adc
      dig->SetToT((Int_t)(dig->GetAdc()/AliTOFGeometry::ToTBinWidth())); //remove the factor 10^3 just to have a reasonable ToT range for raw data simulation even in the case of non-realistic ToT distribution (n.b. fAdc is practically an arbitrary quantity, and ToT has no impact on the TOF reco for non-miscalibrated digits)
    }
  }

  if(!isToTSimulated)
    AliDebug(1,"Standard Production, no miscalibrated digits");
  else
    if(!misCalibPars)
      AliDebug(1,"Standard Production, no miscalibrated digits");
    else
      AliDebug(1,"Simulating miscalibrated digits");

  return;
}
*/
