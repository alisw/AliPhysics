/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//_________________________________________________________________________
// Short description  
//
/*-- Author: Maxim Volkov (RRC KI)
             Dmitri Peressounko (RRC KI & SUBATECH)
             Yuri Kharlov (IHEP & SUBATECH)     */

//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"

// --- Standard library ---
//#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <netinet/in.h>

// --- AliRoot header files ---
#include "AliPHOSDigit.h"
#include "AliPHOSConTableDB.h"
#include "AliPHOSBeamTestEvent.h"
#include "AliPHOSRaw2Digits.h"
#include "AliPHOSv1.h"
#include "../EVGEN/AliGenBox.h"
#include "AliRun.h"

ClassImp(AliPHOSRaw2Digits)
  
  
//____________________________________________________________________________ 
  AliPHOSRaw2Digits::AliPHOSRaw2Digits():TTask() 
{
  fInName="";  
  fMK1 = 0x0123CDEF ;
  fMK2 = 0x80708070 ;
  fMK3 = 0x4321ABCD ;
  fMK4 = 0x80618061 ;
  fCKW = 0x4640E400 ;
  fDebug = kFALSE;             //  Debug flag
  fIsInitialized = kFALSE ;
  fTarget[0] = 0 ;
  fTarget[1] = 0 ;
  fTarget[2] = 0 ;
  fDigits = 0 ;
  fPHOSHeader =0 ;
  fEvent = 0 ;
  fctdb = 0;
}
//____________________________________________________________________________ 
  AliPHOSRaw2Digits::AliPHOSRaw2Digits(const char * filename):TTask("Default","") 
{
  fInName=filename;
  TString outname(fInName) ;
  outname.ToLower() ;
  outname.ReplaceAll(".fz",".root") ;
  outname.ReplaceAll(".gz","") ;
  SetTitle(outname) ;

  fMK1 = 0x0123CDEF ;
  fMK2 = 0x80708070 ;
  fMK3 = 0x4321ABCD ;
  fMK4 = 0x80618061 ;
  fCKW = 0x4640E400 ;
  fDebug = kFALSE;             //  Debug flag
  fIsInitialized = kFALSE ;
  fTarget[0] = 0 ;
  fTarget[1] = 0 ;
  fTarget[2] = 0 ;
  fDigits = 0 ;
  fPHOSHeader =0 ;
  fEvent = 0 ;
  fctdb = 0;
}

//____________________________________________________________________________ 
AliPHOSRaw2Digits::~AliPHOSRaw2Digits()
{
  if(fPHOSHeader)
    fPHOSHeader->Delete() ;
  if(fDigits){
    fDigits->Delete() ;
    delete fDigits ;
  }
  
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Exec(Option_t * option){
  //This is steering method performing all the conversion

  if(!fIsInitialized) //need initialization
    if(!Init())       //failed to initialize
      return ;

  ProcessRawFile() ;

  FinishRun() ;
} 
//____________________________________________________________________________ 
Bool_t AliPHOSRaw2Digits::Init(void){
  //Create PHOS geometry, sets magnetic field to zero, 
  //create Generator - to store target position, 
  //opens out file, creates TreeE and make initialization of contaniers


  if(fIsInitialized)
    return kTRUE;

  //Create PHOS
  new AliPHOSv1("PHOS","GPS2") ;

  //Set Magnetic field
  gAlice->SetField(0,2);  

  //Set positin of the virtex
  AliGenBox *gener = new AliGenBox(1);
  Float_t ox = fTarget[1]; 
  Float_t oy = fTarget[2]-460.; 
  Float_t oz = fTarget[0];
  gener->SetOrigin(ox, oy, oz);

  //  Create the output file
  TString outname("") ;
  if(strstr(GetTitle(),"root")){
    outname=GetTitle();
  }
  else{
    outname = fInName ;
    outname.ToLower() ;
    outname.ReplaceAll(".fz",".root") ;
  }
  TFile *rootfile = new TFile(outname,"recreate");
  rootfile->SetCompressionLevel(2);

  // Create the Root Trees
  gAlice->MakeTree("E") ;

  //Make container for digits
  fDigits = new TClonesArray("AliPHOSDigit",1000) ;

  //Fill now TreeE
  fPHOSHeader = new  AliPHOSBeamTestEvent() ;
  Int_t splitlevel = 0 ;
  Int_t bufferSize = 32000 ;    
  TBranch * headerBranch = gAlice->TreeE()->Branch("AliPHOSBeamTestEvent", 
						   "AliPHOSBeamTestEvent", 
						   &fPHOSHeader,bufferSize,splitlevel);
  headerBranch->SetName("AliPHOSBeamTestEvent") ;

  fIsInitialized = kTRUE ;
  return kTRUE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSRaw2Digits::ProcessRawFile(){

  //Method opens zebra file and reads successively bytes from it,
  //filling corresponding fields in header and digits.


  fStatus= -3 ;
  //First of all, open file and check if it is a zebra file

  char command[256];
  sprintf(command,"zcat %s",fInName.Data());
  FILE *dataFile = popen(command, "r");
  if (dataFile == NULL) {
    Warning("ProcessRawFile", " Cannot open file %s\n", fInName.Data() ) ;
    perror(fInName.Data()) ;
    fStatus = -1 ;
    return kFALSE ;
  }
  printf("Open pipe: %s\n",command);

  // Check if byte ordering is little-endian 
  UInt_t w = 0x12345678;
  Int_t swapo = memcmp(&w, "\x78\x56\x34\x12", sizeof(UInt_t)) == 0;
  if(fDebug)
    Info("ProcessRawFile", "swapo=%f\n", swapo ) ;
  
  
  UInt_t recBuf[300] ;

  // Read physical record control words 
  UInt_t nb = 8*sizeof(UInt_t);
  Int_t n = fread(recBuf, nb, 1, dataFile);
  if(static_cast<UInt_t>(n) != 1) {
    if (n < 0 )
      perror(fInName.Data());
    else
      Error("ProcessRawFile", "Could not read physical record control words" ) ;
    fStatus = -2 ;
    return kFALSE;
  }
  
  if(fDebug)
    Info("ProcessRawFile", "recbuf[0] = %d\n", recBuf[0] );
  
  // Check if it is a ZEBRA file and if the words are swapped 
  UInt_t swapi = 0 ;
  if (recBuf[0] != fMK1) {
    Swab4(recBuf, &w, 1);
    if (w != fMK1) {
      Error("ProcessRawFile", "Does not look like a ZEBRA file\n" ) ;
      pclose(dataFile) ;
      fStatus = -2 ;
      return kFALSE;
    }
    swapi=1 ;
  }
  
  if(fDebug){
    TString message ; 
    message  = "        w = %f\n" ; 
    message += "    swapi = %f\n" ; 
    Info("ProcessRawFile", message.Data(), w, swapi ) ; 
  }
  
  // Get number of words in physical record 
  UInt_t  nwphr ;
  if (swapi)
    Swab4(&recBuf[4],&nwphr,1);
  else 
    nwphr = recBuf[4];
  nwphr*=2; // 1998 -- Now we have 2 records 150 words each 


  //start loop over data  
  // Read rest of record 
  nb = (nwphr-8)*sizeof(UInt_t);
  n = fread(&recBuf[8], nb, 1, dataFile) ;
  if (static_cast<UInt_t>(n) != 1) {
    if (n < 0 ){
      perror(fInName.Data());
      fStatus = -2 ;
      return kFALSE;
    }
  }
  nb = nwphr *sizeof(UInt_t);

  UInt_t userVector[16] ;
  UInt_t zheader[12];    
  UShort_t pattern ;
  UShort_t scanning[32] ;
  UShort_t charge[12];
  UInt_t scaler[12]; 
  UShort_t tdc2228[32];
  
  //read untill the end of file
  fEvent=0 ;
  while(1){

    //    StartNewEvent() ;
    fDigits->Delete() ;
    gAlice->SetEvent(fEvent) ;
	  
    Int_t i ;
    for(i=0;i<16;i++)
      userVector[i]=*(recBuf+21+i);
    if(!swapi)
      Swab4(userVector, userVector, 16);     
    fPHOSHeader->SetUserVector(userVector) ;
    
    
    // ZEBRA event header
    for(i=0;i<12;i++)
      zheader[i]=*(recBuf+47+i);
    if(swapi)
      Swab4(zheader, zheader, 12);
    fPHOSHeader->SetHeader(zheader) ;
    
    // Swap input 
    if (swapi)
      Swab4(recBuf, recBuf, nwphr);
    
    /* Physical record control words */
    UInt_t * recptr = recBuf;  //Pointer to current position

    if(recptr[7] != 1) {
      Error("ProcessRawFile", "Cannot handle fast blocks" ) ; 
      fStatus = -2 ;
      return kFALSE;
    }    
    recptr += 8;
    
    // Logical record control words   
    UInt_t lrtyp = recptr[1];
    if (lrtyp != 3) {
      Error("ProcessRawFile", "Can not handle logical record type %d", lrtyp ) ;
      fStatus = -2 ;
      return kFALSE;
    }
    
    recptr += 2;
    if (recptr[0] != fCKW) {
      Error("ProcessRawFile", "Bad check word" ) ;
      fStatus = -2 ;
      return kFALSE;
    }
    
    UInt_t  nwuh = recptr[9];
    recptr += 10+nwuh;
    
    // Bank system words 
    UInt_t nd = recptr[8];	      /* Number of data words */
    recptr += 10;
    
    // Data words 
    UInt_t evtno = recptr[2];			/* Event number */
    
    if(fDebug)
       Info("ProcessRawFile", "evtno= %d", evtno);
    
    UInt_t nh = recptr[4];		     /* Number of header words in data bank */
    recptr += nh;
    
    // Unswap data from VME 
    if (swapi)
      Swab4(recptr, recptr, nd-nh-3);
    
    // Give buffer to monitor program 
    //  UInt_t esize = nd-nh-3;
    //    if (swapo)
    //       Swab2(recptr, recptr, esize); 
    // Two byte data are correct after this operation. 
    //But we're in trouble if the data array contains 4 byte data!
    
    // From now on deal with VME data (MSB first, or network byte order).
    
    
    // Fill the event with data from ADCs
    UChar_t *byteptr=(UChar_t*)recptr;
    
    // Trigger bit register  
    pattern=ntohs(*(UShort_t*)byteptr);
    fPHOSHeader->SetPattern(pattern) ;
    byteptr+=sizeof(UShort_t);
    
    // Either peak ADCs, 10 modulesX8=80 channels, 
    //or Kurchatov 64+2 channel ADC 
    //(the rest of the channels padded with 0xffff) 
    for(i=0;i<80;i++){
      Int_t peak = static_cast<Int_t>(ntohs(*(UShort_t*)byteptr));
      //make digit
      Int_t absID = fctdb->Raw2AbsId(i) ;
      if(absID > 0)
	new((*fDigits)[i])AliPHOSDigit(-1,absID,peak,0.,i) ;
      if(fDebug){
	if(peak>(UShort_t)1000)
	  Info("ProcessRawFile", "event= %d peak[%d] = %f", fEvent, i, peak);
      }
      byteptr+=sizeof(UShort_t);
    }
    
    // Scanning ADCs, 4 modulesX8=32 channels
    for(i=0;i<32;i++){
      scanning[i]=ntohs(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetScanning(scanning) ;
    
    // Charge ADCs, 1 moduleX12=12 channels
    for(i=0;i<12;i++){
      charge[i]=ntohs(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetCharge(charge) ;
    
    // Scalers, 1 moduleX12=12 (4 byte) channels
    for(i=0;i<12;i++){
      scaler[i]=ntohl(*(UInt_t*)byteptr);
      byteptr+=sizeof(UInt_t);
    }
    fPHOSHeader->SetScaler(scaler) ;
    
    // LeCroy TDC 2228A, 4 moduleX8=32 channels
    for(i=0;i<8;i++){
      tdc2228[i]=ntohs(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetTDC(tdc2228) ;

    WriteDigits() ;
    if(fDebug)
      Info("ProcessRawFile", "event= %d written", fEvent) ;
 
    // Read next record 
    UInt_t nb = nwphr *sizeof(UInt_t);
    n = fread( recBuf, nb,1,dataFile);
    if (n < 0 ){
      perror(fInName);
      fStatus = -2 ;
      return kFALSE;
    }
    if (static_cast<UInt_t>(n) != 1) {
      pclose(dataFile) ;
      fStatus = 1 ;
      return kTRUE ; //all read
    }
    fEvent++ ;
  }
  
  fStatus = 1 ;  
  return kTRUE ;  
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Swab4(void *from, void *to, size_t nwords){
  // The function swaps 4 bytes: byte#3<-->byte#0, byte#2<-->byte#1 
  register char *pf=static_cast<char*>(from) ;
  register char *pt=static_cast<char*>(to) ;
  register char c;
  while (nwords-- > 0 ) {
    c = pf[0];
    pt[0] = pf[3];
    pt[3] = c;
    c = pf[1];
    pt[1] = pf[2];
    pt[2] = c;
    pf += 4;
    pt += 4;
  }
}

//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Swab2(void *from, void *to, size_t nwords)
{ //The function swaps 2x2 bytes: byte#0<-->byte#1, byte#2<-->byte#3 
  register char *pf=static_cast<char*>(from) ;
  register char *pt=static_cast<char*>(to);
  register char c;   
  while (nwords-- > 0 ) {
    c = pf[0];
    pt[0] = pf[1];
    pt[1] = c;
    c = pf[2];
    pt[2] = pf[3];
    pt[3] = c;
    pf += 4;
    pt += 4;
  }
}

//____________________________________________________________________________ 
void AliPHOSRaw2Digits::FinishRun(){
  //Write geometry and header tree
  gAlice->Write(0,TObject::kOverwrite);
  gAlice->TreeE()->Write(0,TObject::kOverwrite);
  
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::WriteDigits(void){
  //In this method we create TreeD, write digits and Raw2Digits to it
  // and write Header to TreeE. Finally we write TreeD to root file 
  
  //Start from Digits
  fDigits->Sort() ;
  fDigits->Expand(fDigits->GetEntriesFast()) ;
  for(Int_t i=0;i<fDigits->GetEntriesFast(); i++)
    static_cast<AliPHOSDigit*>(fDigits->At(i))->SetIndexInList(i) ;

  char hname[30];
  sprintf(hname,"TreeD%d",fEvent);
  TTree * treeD = new TTree(hname,"Digits");
  //treeD->Write(0,TObject::kOverwrite);
  
  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = treeD->Branch("PHOS",&fDigits,bufferSize);
  digitsBranch->SetTitle("Default");
  
  // -- Create Digitizer branch
  Int_t splitlevel = 0 ;
  const AliPHOSRaw2Digits * d = this ;
  TBranch * digitizerBranch = treeD->Branch("AliPHOSRaw2Digits", 
					    "AliPHOSRaw2Digits", &d,bufferSize,splitlevel); 
  digitizerBranch->SetTitle("Default");
  
  digitsBranch->Fill() ;
  digitizerBranch->Fill() ; 
  treeD->Write(0,TObject::kOverwrite);
 
  delete treeD ;

  //Write header
  gAlice->TreeE()->Fill();
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Print(Option_t * option)const{
  
  TString message ;
  message  = "----------AliPHOSRaw2Digits---------- \n" ;
  message += "Input stream \n" ;
  message += "Current input  File: %s\n" ; 
  message += "Current output File: %s\n" ; 
  message += "Events processes in the last file %d\n" ; 
  message += "Input file status\n" ;
  switch (fStatus){
  case 0: message += "`Have not processed yet'\n" ;
    break ;
  case 1: message += "`Processed normally'\n" ;
    break ;
  case -1: message += "`File not found'\n" ;
    break ;
  case -2: message += "`Error in reading'\n" ;
    break ;
  case -3: message += "'Interupted'\n" ;
  default: ;
  }
  message += "Connection table: " ;
  if(fctdb)
    message += "%s %s \n" ; 
  else
    message += " no DB \n" ;

  Info("Print", message.Data(),  
       fInName.Data(), 
       GetTitle(), 
       fEvent, 
       fctdb->GetName(), fctdb->GetTitle() ) ; 
}
