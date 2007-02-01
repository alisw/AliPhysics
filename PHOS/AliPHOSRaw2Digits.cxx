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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.14  2007/01/23 10:27:37  alibrary
 * Adding include files where needed for latest ROOT
 *
 * Revision 1.13  2006/09/07 18:31:08  kharlov
 * Effective c++ corrections (T.Pocheptsov)
 *
 * Revision 1.12  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Class designed to convert raw data to aliroot format. We assume, that
// prototype is situated in the center of 3 PHOS module and convert prototype
// outpur to AliPHOSDigits. In addition, we fill branch of TreeE with 
// AliPHOSBeamTestEvent, contaning description of event(triggers etc).
// Note, that one byte per channel in raw data is transvormed to class 
// AliPHOSDigit, so finale zise increase ~ 100 times. So output can be split 
// into peases of reasonable size: each file can not contain more than 
// fMaxPerFile: if there are still events in raw file, then new directory 
// is created and header+digits files are written to it.
// 
// Use Case:
//   AliPHOSRaw2Digits * r = new AliPHOSRaw2Digits("path/input.file") ;
//                   //note, that it can be gzipped file!
//   //Set position of the target in the given run.
//   //Z axis along beam direction, from target to prototype (0-surface of prototype)
//   //X axis along columns of prototype (0-center of prototype)
//   //Y axis along raws of prototype    (0-center of prototype)
//   Double_t pos[3]={0,0,-452.} ;
//   r->SetTargetPosition(pos) ;
//   //Read/create connection Table:
//   TFile f("ConTableDB.root") ;
//   AliPHOSConTableDB * cdb = f.Get("AliPHOSConTableDB") ;
//   f.Close() ;
//   r->SetConTableDB(cdb) ;
//   r->ExecuteTask() ;
//
// As a result files galice.root and PHOS.Digits.root should be produced in 
// current dir, and, possibly, dirs 1,2,3... each with galice.root and PHOS.Digits.root,
// where the rest of data are written. 
//
/*-- Author: Maxim Volkov (RRC KI)
             Dmitri Peressounko (RRC KI & SUBATECH)
             Yuri Kharlov (IHEP & SUBATECH)     */

//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <Bytes.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

// --- Standard library ---

#include <stdio.h>

// --- AliRoot header files ---
#include "AliPHOSDigit.h"
#include "AliPHOSConTableDB.h"
#include "AliPHOSBeamTestEvent.h"
#include "AliPHOSRaw2Digits.h"
#include "AliRun.h"

ClassImp(AliPHOSRaw2Digits)
  
//____________________________________________________________________________ 
AliPHOSRaw2Digits::AliPHOSRaw2Digits() : 
  fDigits(0),
  fPHOSHeader(0),
  fctdb(0),
  fHeaderFile(0),
  fDigitsFile(0),
  fBeamEnergy(0.f),
  fMaxPerFile(20000),
  fEvent(0),
  fStatus(0),
  fInName(""),
  fDebug(kFALSE),
  fIsInitialized(kFALSE),
  fMK1(0x0123CDEF),
  fMK2(0x80708070),
  fMK3(0x4321ABCD),
  fMK4(0x80618061),
  fCKW(0x4640E400)
{
  //As one can easily see, this is constructor.
  fTarget[0] = 0 ;
  fTarget[1] = 0 ;
  fTarget[2] = 0 ;
}

//____________________________________________________________________________ 
AliPHOSRaw2Digits::AliPHOSRaw2Digits(const char * filename) : 
  TTask("Default",""),
  fDigits(0),
  fPHOSHeader(0),
  fctdb(0),
  fHeaderFile(0),
  fDigitsFile(0),
  fBeamEnergy(0.f),
  fMaxPerFile(20000),
  fEvent(0),
  fStatus(0),
  fInName(filename),
  fDebug(kFALSE),
  fIsInitialized(kFALSE),
  fMK1(0x0123CDEF),
  fMK2(0x80708070),
  fMK3(0x4321ABCD),
  fMK4(0x80618061),
  fCKW(0x4640E400)
{
  //this constructor should be normally used. Parameters: input file 
  TString outname(fInName) ;
  outname.ToLower() ;
  outname.ReplaceAll(".fz",".root") ;
  outname.ReplaceAll(".gz","") ;
  SetTitle(outname);

  fTarget[0] = 0 ;
  fTarget[1] = 0 ;
  fTarget[2] = 0 ;
}

//____________________________________________________________________________ 
AliPHOSRaw2Digits::AliPHOSRaw2Digits(AliPHOSRaw2Digits & r2d) :
  TTask(r2d.GetName(), r2d.GetTitle()),
  fDigits(r2d.fDigits),
  fPHOSHeader(r2d.fPHOSHeader),
  fctdb(new AliPHOSConTableDB(*r2d.fctdb)),
  fHeaderFile(new TFile(r2d.fHeaderFile->GetName(), "new" )),
  fDigitsFile(new TFile(r2d.fDigitsFile->GetName(), "new" )),
  fBeamEnergy(r2d.fBeamEnergy),
  fMaxPerFile(r2d.fMaxPerFile),
  fEvent(r2d.fEvent),
  fStatus(r2d.fStatus),
  fInName(r2d.fInName),
  fDebug(kFALSE),
  fIsInitialized(kFALSE),
  fMK1(r2d.fMK1),
  fMK2(r2d.fMK2),
  fMK3(r2d.fMK3),
  fMK4(r2d.fMK4),
  fCKW(r2d.fCKW)
{
  // cpy ctor. wrong. because dtor can delete fDigits twice (or n times you copy AliPHOSRaw2Digits)
  //because fHeaderFile and fDigitsFile will recreate existing files etc.
  fTarget[0] = r2d.fTarget[0] ;
  fTarget[1] = r2d.fTarget[1] ;
  fTarget[2] = r2d.fTarget[2] ;
}

//____________________________________________________________________________ 
AliPHOSRaw2Digits::~AliPHOSRaw2Digits()
{
//destructor
  if(fPHOSHeader)
    fPHOSHeader->Delete() ;
  if(fDigits){
    fDigits->Delete() ;
    delete fDigits ;
  }
  
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Exec(const Option_t *){
  //This is steering method performing all the conversion

  if(!fIsInitialized) //need initialization
    if(!Init())       //failed to initialize
      return ;

  ProcessRawFile() ;

} 
//____________________________________________________________________________ 
Bool_t AliPHOSRaw2Digits::Init(void){
  //Makes initialization of contaniers

  if(fIsInitialized)
    return kTRUE;

  //Make container for digits
  fDigits = new TClonesArray("AliPHOSDigit",1000) ;
  fPHOSHeader = new  AliPHOSBeamTestEvent() ;
  fIsInitialized = kTRUE ;
  return StartRootFiles() ;

}
//____________________________________________________________________________ 
Bool_t AliPHOSRaw2Digits::StartRootFiles(void ) const {
//   //Create PHOS geometry, sets magnetic field to zero, 
//   //create Generator - to store target position, 
//   //opens out file, creates TreeE 

//   //create gAlice if nececcary
//   if(!gAlice)
//     new AliRun("gAlice","The ALICE Off-line Simulation Framework") ;

//   //Create PHOS
//   if(!gAlice->GetModule("PHOS"))
//     new AliPHOSv1("PHOS","GPS2") ;

//   //Set Magnetic field
//   gAlice->SetField(0,2);  

//   //Set positin of the virtex
//   AliGenerator * gener = gAlice->Generator() ; 
//   if(!gener)    
//     gener = new AliGenBox(1);
//   Float_t ox = fTarget[1]; 
//   Float_t oy = fTarget[2]+460.; 
//   Float_t oz = fTarget[0];
//   gener->SetOrigin(ox, oy, oz);

//   //make directory 
//   Int_t nRootFile = (fEvent+1)/fMaxPerFile ;	
//   if(nRootFile){
//     char dname[20];
//     sprintf(dname,"%d",nRootFile) ;
//     if(gSystem->AccessPathName(dname)) //strange return: 0 if exists
//       if(gSystem->MakeDirectory(dname)!=0)
// 	Fatal("StartRootFiles","Can not make directory %s \n",dname) ;
    
//     if(!gSystem->ChangeDirectory(dname))
//       Fatal("StartRootFiles","Can not cd to %s\n",dname) ;
//   }

//   //  Create the output file
//   TString outname("") ;
//   if(strstr(GetTitle(),"root")){
//     outname=GetTitle();
//   }
//   else{
//     outname = fInName ;
//     outname.ToLower() ;
//     outname.ReplaceAll(".fz",".root") ;
//   }

//   fHeaderFile = new TFile(outname,"recreate");
//   fHeaderFile->SetCompressionLevel(2);
  
//   // Create the Root Trees
  
//   gime->MakeTree("E") ;
  
//   //Fill now TreeE
//   Int_t splitlevel = 0 ;
//   Int_t bufferSize = 32000 ;    
//   TBranch * headerBranch = gAlice->TreeE()->Branch("AliPHOSBeamTestEvent", 
// 						   "AliPHOSBeamTestEvent", 
// 						   &fPHOSHeader,bufferSize,splitlevel);
//   headerBranch->SetName("AliPHOSBeamTestEvent") ;

// //   if(fToSplit){
// //     fDigitsFile = new TFile("PHOS.Digits.root","recreate") ;
// //     fDigitsFile->SetCompressionLevel(2) ;
// //   }
   return kTRUE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSRaw2Digits::CloseRootFiles(void ){
  //cleans everething to start next root file
  if(fHeaderFile){
    printf("writing gAlice \n") ;
    fHeaderFile->cd() ;
    gAlice->Write(0,TObject::kOverwrite);
    gAlice->TreeE()->Write(0,TObject::kOverwrite);
  }

  delete gAlice ;
  
  if(fHeaderFile){
    fHeaderFile->Close() ;
    delete fHeaderFile ;
    fHeaderFile = 0;
  }   
	
  if(fDigitsFile){
    fDigitsFile->Close() ;
    delete fDigitsFile ;
    fDigitsFile = 0 ;
  }
  
  Int_t nRootFile = (fEvent-1)/fMaxPerFile ;	
  if(nRootFile){
    if(!gSystem->ChangeDirectory("../")){
     Fatal("CloseRootFile","Can not return to initial dir \n") ;      
     return kFALSE ;
    }
  }
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
    if((fEvent%fMaxPerFile == 0) && fEvent ){
      CloseRootFiles() ;
      StartRootFiles() ;
    }
    gAlice->SetEvent(fEvent%fMaxPerFile) ;

    //Set Beam Energy
    fPHOSHeader->SetBeamEnergy(fBeamEnergy) ;
	  
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
    pattern=net2host(*(UShort_t*)byteptr);
    fPHOSHeader->SetPattern(pattern) ;
    byteptr+=sizeof(UShort_t);
    
    // Either peak ADCs, 10 modulesX8=80 channels, 
    //or Kurchatov 64+2 channel ADC 
    //(the rest of the channels padded with 0xffff) 
    for(i=0;i<80;i++){
      Int_t peak = static_cast<Int_t>(net2host(*(UShort_t*)byteptr));
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
      scanning[i]=net2host(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetScanning(scanning) ;
    
    // Charge ADCs, 1 moduleX12=12 channels
    for(i=0;i<12;i++){
      charge[i]=net2host(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetCharge(charge) ;
    
    // Scalers, 1 moduleX12=12 (4 byte) channels
    for(i=0;i<12;i++){
      scaler[i]=net2host(*(UInt_t*)byteptr);
      byteptr+=sizeof(UInt_t);
    }
    fPHOSHeader->SetScaler(scaler) ;
    
    // LeCroy TDC 2228A, 4 moduleX8=32 channels
    for(i=0;i<8;i++){
      tdc2228[i]=net2host(*(UShort_t*)byteptr);
      byteptr+=sizeof(UShort_t);
    }
    fPHOSHeader->SetTDC(tdc2228) ;

    WriteDigits() ;
    if(fDebug)
      Info("ProcessRawFile", "event= %d written", fEvent) ;
 
    // Read next record 
    nb = nwphr *sizeof(UInt_t);
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
  CloseRootFiles() ;
  
  fStatus = 1 ;  
  return kTRUE ;  
}

//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Swab4(void *from, void *to, size_t nwords)const
{
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
void AliPHOSRaw2Digits::Swab2(void *from, void *to, size_t nwords)const
{ 
  //The function swaps 2x2 bytes: byte#0<-->byte#1, byte#2<-->byte#3 
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
void AliPHOSRaw2Digits::WriteDigits(void){
  //In this method we create TreeD, write digits and Raw2Digits to it
  // and write Header to TreeE. Finally we write TreeD to root file 
  
  //Start from Digits
  fDigits->Sort() ;
  fDigits->Expand(fDigits->GetEntriesFast()) ;
  for(Int_t i=0;i<fDigits->GetEntriesFast(); i++)
    static_cast<AliPHOSDigit*>(fDigits->At(i))->SetIndexInList(i) ;

  char hname[30];
  sprintf(hname,"TreeD%d",fEvent%fMaxPerFile);
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

  if(fDigitsFile)
    fDigitsFile->cd() ;
  digitsBranch->Fill() ;
  digitizerBranch->Fill() ; 
  treeD->Write(0,TObject::kOverwrite);
 
  delete treeD ;

  //Write header
  fHeaderFile->cd() ;
  gAlice->TreeE()->Fill();
}
//____________________________________________________________________________ 
void AliPHOSRaw2Digits::Print(const Option_t *)const{
  //prints current configuration and status.

  printf("----------AliPHOSRaw2Digits---------- \n") ;
  printf("Current input  File: %s\n",fInName.Data()) ; 
  printf("Current output File: %s\n", GetTitle()); 
  printf("Events processes in the last file %d\n",fEvent) ; 
  printf("Input file status\n") ;
  switch (fStatus){
  case 0: printf("`Have not processed yet'\n") ;
    break ;
  case 1: printf("`Processed normally'\n") ;
    break ;
  case -1: printf("`File not found'\n") ;
    break ;
  case -2: printf("`Error in reading'\n") ;
    break ;
  case -3: printf("'Interupted'\n") ;
  default: ;
  }
  printf("Connection table: " );
  if(fctdb)
    printf("%s %s \n",fctdb->GetName(), fctdb->GetTitle() ) ; 
  else
    printf(" no DB \n" );

}
