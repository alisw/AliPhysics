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

#include <Riostream.h>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TRandom.h>

#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliRun.h"

#include "AliTOFCal.h"
#include "AliTOFcalib.h"
#include "AliTOFChannel.h"
#include "AliTOFDigitizer.h"
#include "AliTOFdigit.h"
#include "AliTOFHitMap.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFSDigit.h"
#include "AliTOF.h"

ClassImp(AliTOFDigitizer)

//___________________________________________
  AliTOFDigitizer::AliTOFDigitizer()  :
    AliDigitizer(),
    fGeom(0x0),
    fDigits(0x0),
    fSDigitsArray(0x0),
    fhitMap(0x0)
{
  // Default ctor - don't use it
}

//___________________________________________
AliTOFDigitizer::AliTOFDigitizer(AliRunDigitizer* manager): 
  AliDigitizer(manager), 
  fGeom(0x0),
  fDigits(0x0),
  fSDigitsArray(0x0),
  fhitMap(0x0)
{
  //ctor with RunDigitizer
}

//------------------------------------------------------------------------
AliTOFDigitizer::AliTOFDigitizer(const AliTOFDigitizer &source):
  AliDigitizer(source),
  fGeom(0x0), 
  fDigits(0),
  fSDigitsArray(0),
  fhitMap(0)
{
  // copy constructor
  this->fDigits=source.fDigits;
  this->fSDigitsArray=source.fSDigitsArray;
  this->fhitMap=source.fhitMap;
  this->fGeom=source.fGeom; 

}

//------------------------------------------------------------------------
  AliTOFDigitizer& AliTOFDigitizer::operator=(const AliTOFDigitizer &source)
{
  // ass. op.
  this->fDigits=source.fDigits;
  this->fSDigitsArray=source.fSDigitsArray;
  this->fhitMap=source.fhitMap;
  this->fGeom=source.fGeom; 
  return *this;

}

//------------------------------------------------------------------------
AliTOFDigitizer::~AliTOFDigitizer()
{
  // Destructor
}

//---------------------------------------------------------------------

void AliTOFDigitizer::Exec(Option_t* /*option*/)
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
  char branchname[20];
  sprintf (branchname, "%s", tof->GetName ());

  fDigits=new TClonesArray("AliTOFdigit",4000);
 
  AliRunLoader* outrl = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  if (outrl == 0x0)
   {
     AliError("Can not find Run Loader in output folder.");
     return;
   }
   
  outrl->CdGAFile();
  TFile *in=(TFile*)gFile;
  TDirectory *savedir=gDirectory;

  if (!in->IsOpen()) {
    AliWarning("Geometry file is not open default  TOF geometry will be used");
    fGeom = new AliTOFGeometryV5();
  }
  else {
    in->cd();
    fGeom = (AliTOFGeometry*)in->Get("TOFgeometry");
  }

  savedir->cd();

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
  fSDigitsArray=new TClonesArray("AliTOFSDigit",1000);
  
  // create hit map (to be created in Init())
  fhitMap = new AliTOFHitMap(fSDigitsArray, fGeom);
  
  // Loop over files to digitize

  for (Int_t inputFile=0; inputFile<fManager->GetNinputs();
       inputFile++) {
    ReadSDigit(inputFile);
   }

  // create digits
  CreateDigits();

  // free used memory for Hit Map in current event
  delete fhitMap;
  fSDigitsArray->Delete();
  delete fSDigitsArray;

  treeD->Fill();
 
  outgime->WriteDigits("OVERWRITE");
  outgime->UnloadDigits();
  fDigits->Delete();
  delete fDigits;

}

//---------------------------------------------------------------------

void AliTOFDigitizer::CreateDigits()
{
  // loop on sdigits container to fill the AliTOFdigit TClonesArray
  // start digitizing all the collected sdigits 

  Int_t ndump=0; // dump the first ndump created digits for each event

  // get the total number of collected sdigits
  Int_t ndig = fSDigitsArray->GetEntriesFast();

  for (Int_t k = 0; k < ndig; k++) {
    
    Int_t  vol[5];  // location for a digit
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
    
    if (isSDigitBad) {
      //AliFatal("strange sdigit found");
      AliFatal(Form("strange sdigit found   %3i  %2i  %2i  %3i    %3i", sector, plate, padz, padx, strip));
    }
    //-------------------------------------------------------
    
    //------------------- Dump section ----------------------
    if(k<ndump){
      cout << k << "-th | " << "Sector " << sector << " | Plate " << plate << " | Strip " << strip << " | PadZ " << padz << " | PadX " << padx << endl;
      cout << k << "-th sdigit" << endl;
      cout << "----------------------------------------------------"<< endl;
    }
    // ------------------------------------------------------
    
    // start loop on number of slots for current sdigit
    for (Int_t islot = 0; islot < nslot; islot++) {
      Float_t  digit[4] = {-1.,-1.,-1.,-1.};     // TOF digit variables
      Int_t tracknum[AliTOFSDigit::kMAXDIGITS];     // contributing tracks for the current slot
      
      Float_t tdc=tofsdigit->GetTdc(islot); digit[0]=tdc;
      Float_t adc=tofsdigit->GetAdc(islot); digit[1]=adc;
      
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

  AliTOFcalib * calib = new AliTOFcalib(fGeom);
  InitDecalibration(calib);
  DecalibrateTOFSignal(calib);
  delete calib;
}

//---------------------------------------------------------------------

void AliTOFDigitizer::ReadSDigit(Int_t inputFile )
{
  // Read sdigits for current event and inputFile; 
  // store them into the sdigits container
  // and update the hit map
  // SDigits from different files are assumed to
  // be created with the same simulation parameters.
  
  // get the treeS from manager
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
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
  TClonesArray * sdigitsDummyContainer= new TClonesArray("AliTOFSDigit",  1000); 

  // check if the branch exist
  TBranch* tofBranch=currentTreeS->GetBranch("TOF");

  if(!tofBranch){
    AliFatal(Form("TOF branch not found for input %d",inputFile));
  }
  
  tofBranch->SetAddress(&sdigitsDummyContainer);           
  
  Int_t nEntries = (Int_t)tofBranch->GetEntries();                                

  // Loop through all entries in the tree
  Int_t nbytes = 0;
  
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
    
    // Import the tree
    nbytes += tofBranch->GetEvent(iEntry);
    
    // Get the number of sdigits
    Int_t ndig = sdigitsDummyContainer->GetEntriesFast();
    
    for (Int_t k=0; k<ndig; k++) {
      AliTOFSDigit *tofSdigit= (AliTOFSDigit*) sdigitsDummyContainer->UncheckedAt(k);
      
      Int_t  vol[5]; // location for a sdigit
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
    sdigitsDummyContainer->Delete();

  } // end loop on entries

  delete sdigitsDummyContainer;

}


//_____________________________________________________________________________
void AliTOFDigitizer::CollectSDigit(AliTOFSDigit * sdigit)
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
void AliTOFDigitizer::InitDecalibration( AliTOFcalib *calib) const {
  // calib->ReadSimParFromCDB("TOF/Calib", 0); // original
  calib->ReadSimParFromCDB("TOF/Calib", -1); // use AliCDBManager's number
}
//---------------------------------------------------------------------
void AliTOFDigitizer::DecalibrateTOFSignal( AliTOFcalib *calib){

  // Read Calibration parameters from the CDB

  AliTOFCal * cal= calib->GetTOFCalSimArray();

  AliDebug(2,Form("Size of AliTOFCal = %i",cal->NPads()));
  for (Int_t ipad = 0 ; ipad<cal->NPads(); ipad++){
    AliTOFChannel *calChannel = cal->GetChannel(ipad);
    Float_t par[6];
    for (Int_t j = 0; j<6; j++){
      par[j]=calChannel->GetSlewPar(j);
    }
  }

  // Initialize Quantities to Simulate ToT Spectra


  TH1F * hToT= calib->GetTOFSimToT();
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
  

  // Loop on TOF Digits

  Bool_t dbEntry=kFALSE;
  Int_t ndigits = fDigits->GetEntriesFast();    
  for (Int_t i=0;i<ndigits;i++){
    AliTOFdigit * dig = (AliTOFdigit*)fDigits->At(i);
    Int_t detId[5];
    detId[0] = dig->GetSector();
    detId[1] = dig->GetPlate();
    detId[2] = dig->GetStrip();
    detId[3] = dig->GetPadz();
    detId[4] = dig->GetPadx();
    // For Data with no Miscalibration, set ToT signal == Adc
    dig->SetToT(dig->GetAdc());
    if(hToT->GetEntries()>0){  
      Float_t trix = 0;
      Float_t triy = 0;
      Float_t simToT = 0;
      while (simToT <= triy){
	trix = gRandom->Rndm(i);
	triy = gRandom->Rndm(i);
	trix = (maxToT-minToT)*trix + minToT; 
	triy = maxToTDistr*triy;
	Int_t binx=hToT->FindBin(trix);
	simToT=hToT->GetBinContent(binx);
      }
    // Setting realistic ToT signal (only for Miscalibrated Data)   
      dig->SetToT(trix);
    }
    Int_t index = calib->GetIndex(detId);     
    AliTOFChannel *calChannel = cal->GetChannel(index);
    // time slewing parameters
    Float_t par[6];
    for (Int_t j = 0; j<6; j++){
      par[j]=calChannel->GetSlewPar(j);
      if(par[j]!=0)dbEntry=kTRUE;
    }
    // the global time shift
    Float_t timedelay = calChannel->GetDelay();
    Float_t tToT= dig->GetToT();
    dig->SetTdcND(dig->GetTdc());
    Float_t tdc = ((dig->GetTdc())*AliTOFGeometry::TdcBinWidth()+32)*1.E-3; //tof signal in ns
    // add slewing effect
    Float_t timeoffset=par[0] + tToT*(par[1] +tToT*(par[2] +tToT*(par[3] +tToT*(par[4] +tToT*par[5])))); 
    Float_t timeSlewed = tdc+timeoffset;
    // add global time shift
    timeSlewed = timeSlewed + timedelay;

    // Setting Decalibrated Time signal    
    dig->SetTdc((timeSlewed*1E3-32)/AliTOFGeometry::TdcBinWidth());   
  }

  if(hToT->GetEntries()<=0 || !dbEntry){
    AliDebug(1,"Standard Production, no miscalibrated digits");   
  }else{
    AliDebug(1,"Miscalibrated digits");   
  }

  return;
}

