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

/*
  $Id$
*/

#include <stddef.h>
#include <iomanip>
#include <Riostream.h>
#include <fstream>
#include <AliRunLoader.h>
#include <AliLoader.h>
#include <AliITS.h>

#include "AliITSspdTestBeam.h"

ClassImp(AliITSspdTestBeam)

//----------------------------------------------------------------------
AliITSspdTestBeam::AliITSspdTestBeam(){
    // Default Constructor for the Task AliITSspdTestBeam.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     A default constructed AliITSspdTestBeam class

    fNBrst        = 0;
    fBrstSize     = 0;
    fBrst         = 0;
    fNData        = 0;
    fData         = 0;
    fHData        = 0;
    fTData        = 0;
    SetTerminationWord();
    fNEvents      = 0;
    fBuffSize     = 0;
    fBuff         = 0;
    fVersion      = 0;
    fITS          = 0;
    fNfiles       = 0;
    fMaxFiles     = 0;
    fFiles        = 0;
    fNeventsStart = 0;
    fNeventsEnd   = 0;
}
//----------------------------------------------------------------------
AliITSspdTestBeam::AliITSspdTestBeam(const Char_t *filename,const Char_t *opt,
                                     AliITS *its){
    // Standard Constructor for the Task AliITSspdTestBeam.
    // Inputs:
    //    const Char_t *filename   File where to read in the SPD test beam data
    //    const Char_t *opt        Option, 2002 => 2002 test beam data.
    // Outputs:
    //    none.
    // Return:
    //     A default constructed AliITSspdTestBeam class

    fNBrst        = 0;
    fBrstSize     = 0;
    fBrst         = 0;
    fNData        = 0;
    fData         = 0;
    fHData        = 0;
    fTData        = 0;
    SetTerminationWord();
    fNEvents      = 0;
    fBuffSize     = 0;
    fBuff         = 0;
    fVersion      = 0;
    fITS          = 0;
    fNfiles       = 0;
    fMaxFiles     = 0;
    fFiles        = 0;
    fNeventsStart = 0;
    fNeventsEnd   = 0;
    //
    fITS          = its;
    fNfiles       = 0;
    OpenInputFile(filename,0,-1);
    if(strcmp(opt,"2002")) cout << "2002 assumed" << endl;
}

//----------------------------------------------------------------------
AliITSspdTestBeam::~AliITSspdTestBeam(){
    // Destructor. Frees up any memory allocated or closes any files opened.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     nothing.

    DeletefBrst();
    DeletefNData();
    DeletefData();
    DeletefHData();
    DeletefTData();
    DeletefFiles();
    DeletefBrstSize();
    DeletefBuff();
    DeletefITS();
    DeletefNeventsStart();
    DeletefNeventsEnd();
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefBrst(){
    // Properly deletes fBrst object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;
    
    if(fBrst){
        for(i=0;i<fNBrst;i++) if(fBrst[i]){
            delete[] fBrst[i];
            fBrst[i] = 0;
        } // end for i/if
        delete[] fBrst;
        fBrst = 0;
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefNData(){
    // Properly deletes fNData object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,np=GetNumberOfPilots();

    if(fNData){
        for(i=0;i<np;i++) if(fNData[i]){
            delete[] fNData[i];
            fNData[i] = 0;
        } // end for i
        delete[] fNData;
        fNData = 0;
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefData(){
    // Properly deletes fData object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,j,np=GetNumberOfPilots();

    if(fData){
        for(i=0;i<np;i++) if(fData[i]){
            for(j=0;j<fNEvents;j++){
                delete[] fData[i][j];
                fData[i][j] = 0;
            } // end for j
            delete[] fData[i];
            fData[i] = 0;
        } // end for i
        delete[] fData;
        fData = 0;
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefHData(){
    // Properly deletes fHData object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,j,np=GetNumberOfPilots();

    if(fHData){
        for(i=0;i<np;i++) if(fHData[i]){
            for(j=0;j<fNEvents;j++){
                delete[] fHData[i][j];
                fHData[i][j] = 0;
            } // end for j
            delete[] fHData[i];
            fHData[i] = 0;
        } // end for i
        delete[] fHData;
        fHData = 0;
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefTData(){
    // Properly deletes fTData object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,j,np=GetNumberOfPilots();

    if(fTData){
        for(i=0;i<np;i++) if(fTData[i]){
            for(j=0;j<fNEvents;j++){
                delete[] fTData[i][j];
                fTData[i][j] = 0;
            } // end for j
            delete[] fTData[i];
            fTData[i] = 0;
        } // end for i
        delete[] fTData;
        fTData = 0;
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefFiles(){
    // Properly deletes fBrst object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;

    if(fFiles){
        for(i=0;i<fMaxFiles;i++){
            if(fFiles[i]!=0) delete fFiles[i];
        } // end for i
    } // end if
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::DeletefBrstSize(){
    // Properly deletes fBrstSize object.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,np=GetNumberOfPilots();

    if(fBrstSize){
        for(i=0;i<np;i++){
            if(fBrstSize[i]){
                delete[] fBrstSize[i];
                fBrstSize[i] = 0;
            } // end if
        } // end for i
        delete[] fBrstSize;
        fBrstSize = 0;
    } // end if
}
//----------------------------------------------------------------------
Int_t AliITSspdTestBeam::OpenInputFile(const Char_t *filename,Int_t start,Int_t end){
    // Opens input file for reading spd test beam data.
    // Inputs:
    //    const Char_t *filename    file name to read data from
    // Outputs:
    //    none.
    // Return:
    //    An error number. 0=success, -1=failure.
    Int_t stat = 0,i;

    if(fMaxFiles==0) {
        fMaxFiles     = 5;
        fFiles        = new ifstream*[fMaxFiles];
        for(i=0;i<fMaxFiles;i++) fFiles[i] = 0;
        fNeventsStart = new Int_t[fMaxFiles];
        fNeventsEnd   = new Int_t[fMaxFiles];
    } // end if
    if(fNfiles==fMaxFiles){// Need to expand array of open files.
        ifstream *tmp[fMaxFiles];
        Int_t st[fMaxFiles],en[fMaxFiles];
        for(i=0;i<fMaxFiles;i++) { // copy pointers into tmp
            tmp[i]    = fFiles[i];
            fFiles[i] = 0;
            st[i]     = fNeventsStart[i];
            en[i]     = fNeventsEnd[i];
        } // end for i
        delete fFiles;
        fMaxFiles    += 5;  // expand by 5.
        fFiles        = new ifstream*[fMaxFiles];
        fNeventsStart = new Int_t[fMaxFiles];
        fNeventsEnd   = new Int_t[fMaxFiles];
        for(i=0;i<fMaxFiles;i++) { // copy pointers back into fFiles
            fFiles[i]        = 0;  // and zero rest.
            fNeventsStart[i] = 0;
            fNeventsEnd[i]   = 0;
            if(i<fNfiles) {
                fFiles[i]        = tmp[i];
                tmp[i]           = 0;
                fNeventsStart[i] = st[i];
                fNeventsEnd[i]   = en[i];
            } // end if i<fNfiles
        } // end for i
        // the array of pointers tmp is deleted automatically.
    } // end if
    // Open file
    fFiles[fNfiles] = new ifstream(filename,ios::in|ios::binary);
    if(fFiles[fNfiles]==0){// file open error
        cout << "Error opening input file " << filename << endl;
        stat = -1;
        return stat;
    } // end if
    fNeventsStart[fNfiles] = start;
    fNeventsEnd[fNfiles]   = end;
    fNfiles++;
    return stat;
}
//----------------------------------------------------------------------
Int_t AliITSspdTestBeam::Read(Int_t i){
    // Read in one buffer's worth of the file.
    // Inputs:
    //    Int_t i  Which file from the array of files to be read in.
    // Outputs:
    //    none.
    // Return:
    //    size of file.
    Int_t filesize=0;

    fFiles[i]->seekg(0,ios::end);
    filesize  = fFiles[i]->tellg();
    fFiles[i]->seekg(0,ios::beg);
    if(fBuff) delete[] fBuff;
    fBuffSize = filesize;
    fBuff     = new UChar_t[fBuffSize];
    fFiles[i]->read((Char_t*)fBuff,fBuffSize);
    fFiles[i] ->close();
    return filesize;
}
//----------------------------------------------------------------------
Int_t AliITSspdTestBeam::Decode(){
    // Decode the fBuff read in.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    sets up the array fData fRH, and fRT
    Int_t stat=0;
    UInt_t *tr;
    union {
        UInt_t  *wd;
        UChar_t *bt;
    }u; // end union
    Int_t size;
    Int_t *ivnt,iburst,ip,np,i,j;
    AliITSspdTestBeamData  *d;
    AliITSspdTestBeamBurst *b;

    np = GetNumberOfPilots();
    cout << "Sizeof(Headder)=" << fRH.SizeOf() << endl;
    cout << "Sizeof(Tail)=" << fRT.SizeOf() << endl;
    ivnt = new Int_t[np];
    for(i=0;i<np;i++) ivnt[i] = 0;
    fRH.SetBuffer(fBuff);//Set Run Header
    size = fRT.SizeOf();
    size = fBuffSize - size;
    cout <<"fBuffSize-sizeof(AliITSspdTestBeamTail) "<< size << endl;
    fRT.SetBuffer(&(fBuff[size]));
    // Check termination
    size -= sizeof(UInt_t);
    tr   = (UInt_t*) &(fBuff[size]);
    if(!(*tr==fTermination)){
        cout << "Error Termination word not found at " << size << " tr=0x" <<
            hex << *tr << dec << endl;
        exit(-1);
    } // end if
    DeletefBrst();
    DeletefNData();
    DeletefData();
    DeletefHData();
    DeletefTData();
    DeletefBrstSize();
    fNEvents  = fRH.GetNumberOfEvents();
    if(fRT.GetNumberOfEvents()>fNEvents) fNEvents  = fRT.GetNumberOfEvents();
    fNBrst    = fNEvents/fRH.GetBurstSize();
    fBrst     = new AliITSspdTestBeamBurst*[fNBrst];
    fBrstSize = new Int_t*[np];
    fNData    = new Int_t*[np];
    fData     = new AliITSspdTestBeamData**[np];
    fHData    = new AliITSspdTestBeamData**[np];
    fTData    = new AliITSspdTestBeamData**[np];
    for(i=0;i<np;i++){
        fBrstSize[i] = new Int_t[fNBrst];
        fNData[i]    = new Int_t[fNEvents];
        fData[i]     = new AliITSspdTestBeamData*[fNEvents];
        fHData[i]    = new AliITSspdTestBeamData*[fNEvents];
        fTData[i]    = new AliITSspdTestBeamData*[fNEvents];
        for(j=0;j<fNEvents;j++){
            fNData[i][j] = 0;
            fData[i][j]  = 0;
            fHData[i][j] = 0;
            fTData[i][j] = 0;
        } // end for j
    } // end for i
    size      = fRH.SizeOf();
    u.bt      = &fBuff[size];
    //
    for(iburst=0;(*(u.wd) != fTermination)&&(u.wd<tr);iburst++){ 
        // loop over Bursts
        b   = new AliITSspdTestBeamBurst(u.bt);
        fBrst[iburst] = b;
        u.bt += b->SizeOf(); // increment wd byte wise
        for(ip=0;ip<np;ip++){  // loop over pilots
            // Get size of data stored for this pilot.
            AliITSTestBeamData::Swapit(4,u.bt);
            fBrstSize[ip][iburst] = (UInt_t) *(u.wd);
            u.bt += sizeof(UInt_t); // increment wd byte wise
            for(i=0;i<fBrstSize[ip][iburst];i++){ // loop over data
                AliITSTestBeamData::Swapit(4,u.bt);
                d = new AliITSspdTestBeamData(u.bt);
                switch (d->Mode()){
                case AliITSTestBeamData::kData :
                    fNData[ip][ivnt[ip]]++;
                    // set pointer to first data member
                    if(fData[ip][ivnt[ip]] == 0 ){
                        fData[ip][ivnt[ip]] = d;
                    } // end if
                    break;
                case AliITSTestBeamData::kHead :
                    fHData[ip][ivnt[ip]]   = d;
                    fTData[ip][ivnt[ip]]   = 0;
                    fNData[ip][ivnt[ip]]   = 0;
                    fData[ip][ivnt[ip]]    = 0;
                    break;
                case AliITSTestBeamData::kTail  : 
                case AliITSTestBeamData::kAbort :
                    fTData[ip][ivnt[ip]++] = d;
                    break;
                default:
                    cout << "Unknown Data Type: wd="<<hex<<*(u.wd)<<dec<<endl; 
                    break;
                } // end switch
                u.bt += d->SizeOf(); // increment wd byte wise
            } // end for i (next data word).
        } // end for loop over pilots (ip).
    } // end for loop over bursts
    delete[] ivnt;
    return stat;
}
//----------------------------------------------------------------------
Int_t AliITSspdTestBeam::DecodeModule(Int_t pilot,Int_t chip){
    // Determines the Module number based on the pilot and chip 
    // valules and the fVersion of the simulation.
    // Inputs:
    //    Int_t   pilot   Pilot number
    //    Int_t   chip    chip number
    // Outputs:
    //    none.
    // Return:
    //    The module number (see simulations geometry).

    switch (fVersion) {
    case 2002:
        if(pilot==0) return chip;
        if(pilot==1) return 2;
        if(pilot==2) return chip+3;
        break;
    default:
        if(pilot==0) return chip;
        if(pilot==1) return 2;
        if(pilot==2) return chip+3;
        break;
    } // end switch
    return -1;
}
//----------------------------------------------------------------------
Int_t AliITSspdTestBeam::DecodeColumn(Int_t pilot,Int_t chip,Int_t colm){
    // Determines the Column number based on the pilot, chip, and column 
    // valules and the fVersion of the simulation.
    // Inputs:
    //    Int_t   pilot   Pilot number
    //    Int_t   chip    chip number
    //    Int_t   colm    Column number
    // Outputs:
    //    none.
    // Return:
    //    The Column number (see simulations geometry).
    const Int_t colmperchip = 160/5; // Should be gotten from AliITSsegmentationSPD

    switch (fVersion) {
    case 2002:
        if(pilot==0) return colm;
        if(pilot==1) return colm+chip*colmperchip;
        if(pilot==2) return colm;
        break;
    default:
        if(pilot==0) return colm;
        if(pilot==1) return colm+chip*colmperchip;
        if(pilot==2) return colm;
        break;
    } // end switch
    return -1;
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::Digitize(Int_t evnt){
    // Write out ITS SPD Digits.
    // Inputs:
    //    Int_t   evnt   events number
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t p,i;
    Int_t module,chip,row,colm,digit[3];
    Int_t oldmodule=-1;

    fLoader->GetRunLoader()->SetEventNumber(evnt);
    fLoader->SetEvent();
    if(!(fLoader->TreeD())){
        fLoader->MakeTree("D");
    } // end if

    fITS->MakeBranch("D");
    //fITS->SetTreeAddress();
    fITS->SetTreeAddressD(fLoader->TreeD());

    for(p=0;p<GetNumberOfPilots();p++)
        for(i=0;i<fNData[p][evnt];i++){
            chip = fData[p][evnt]->Chip(i);
            module = DecodeModule(p,chip);
            row  = fData[p][evnt]->Row(i);
            colm = DecodeColumn(p,chip,fData[p][evnt]->Colm(i));
            digit[0] = row; digit[1] = colm, digit[2] = 1;
            fITS->AddRealDigit(0,digit);
            if(module!=oldmodule) { // New Module
                oldmodule= module;
                fLoader->TreeD()->Fill();
                fITS->ResetDigits();
            } // end if
    } // end for p
    fITS->ClearModules();
    fLoader->TreeD()->GetEntries();
    fLoader->TreeD()->AutoSave();
    fLoader->TreeD()->Reset();
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintHeadder(ostream *os){
    // Prints the Run Headder.
    // Inputs:
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    fRH.Print(os);
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintTrailer(ostream *os){
    // Prints the Run Trailer.
    // Inputs:
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    fRT.Print(os);
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintBurstInfo(Int_t i,ostream *os){
    // Prints the specified burts information
    // Inputs:
    //    Int_t    i     Specifies which burst to print.
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    fBrst[i]->Print(os);
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintEventData(Int_t evnt,ostream *os){
    // Prints the data for a specified event number.
    // Inputs:
    //    Int_t   evnt   events number
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t p,i;

    for(p=0;p<GetNumberOfPilots();p++)
        for(i=0;i<fNData[p][evnt];i++){
            *os << "Pilot=" << setw(3) << p << " ";
            fData[p][evnt]->Print(os,i);
            *os << endl;
    } // end for p
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintEventHead(Int_t evnt,ostream *os){
    // Prints the data Headder for a specified event number.
    // Inputs:
    //    Int_t   evnt   events number
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t p;

    for(p=0;p<GetNumberOfPilots();p++){
            *os << "Pilot=" << setw(3) << p << " " << *(fHData[p][evnt]);
    } // end for p
}
//----------------------------------------------------------------------
void AliITSspdTestBeam::PrintEventTail(Int_t evnt,ostream *os){
    // Prints the data Trailer for a specified event number.
    // Inputs:
    //    Int_t   evnt   events number
    //    ostream *os    Output stream where the data should go to
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t p;

    for(p=0;p<GetNumberOfPilots();p++){
            *os << "Pilot=" << setw(3) << p << " " << *(fTData[p][evnt]);
    } // end for p
}
//============================================================================
void AliITSTestBeamData::Swapit(Int_t i,UChar_t *a){
    // Swap the order of the bits
    // Inputs:
    //    Int_t   i  Size of UChar_t array a.
    //    UChar_t a  Array of bytes to be swapped.
    // Outputs:
    //    UChar_t a  Array with bytes swapped
    // Return:
    //    none.
    Int_t j;
    UChar_t c[i];

    for(j=0;j<i;j++) c[j] = a[i-j-1];
    for(j=0;j<i;j++) a[j] = c[j];
}
//============================================================================
void AliITSspdTestBeamHeader::Print(ostream *os){
    // print out the header information
    // Inputs:
    //    ostream *os  Pointer to the output stream.
    // Outputs:
    //    none.
    // Return:
    //    none.
/*
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
*/
    *os<<"Version: "<<GetVersion()<<" Written: "<<fUnion->fHead.fDate;
    *os<<" " << fUnion->fHead.fTime << endl;
    *os<<"Buffer Size [0], [1], [2]: " << GetBuffSize(0) << ",";
    *os<<GetBuffSize(1)<<"," <<GetBuffSize(2) << endl;
    *os<<"Test Pulse: " << GetTestPulse() << endl;
    *os<<"Trigger Mode General, [0], [1], [2]: " << GetTriggerMode();
    *os<<","<<GetTrigger(0)<<","<<GetTrigger(1)<< ",";
    *os<<GetTrigger(2) << endl;
    *os<<"Number of Events: " << GetNumberOfEvents() << " Burst Size: ";
    *os<<GetBurstSize() << endl;
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSspdTestBeamHeader &p){
    // Standard output operator. See Print
    // Inputs:
    //    ostream                 &os  the output stream.
    //    AliITSspdTestBeamHeader &p   the data to be printed out.
    // Outputs:
    //    none.
    // Return:
    //    ostream &os pointing now to the end of the present stream.

    p.Print(&os);
    return os;
}
//============================================================================
void AliITSspdTestBeamTail::Print(ostream *os){
    // print out the Tail information
    // Inputs:
    //    ostream *os  Pointer to the output stream.
    // Outputs:
    //    none.
    // Return:
    //    none.
/*
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
*/
    *os << "Number of Events: "<< GetNumberOfEvents() << " Written: "
        << fUnion->fTail.fDate;
    *os << " " << fUnion->fTail.fTime << endl;
    *os <<"Termination Flag: " << GetTermMode() << endl; 
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSspdTestBeamTail &p){
    // Standard output operator. See Print
    // Inputs:
    //    ostream                 &os  the output stream.
    //    AliITSspdTestBeamHeader &p   the data to be printed out.
    // Outputs:
    //    none.
    // Return:
    //    ostream &os pointing now to the end of the present stream.

    p.Print(&os);
    return os;
}
//============================================================================
void AliITSspdTestBeamBurst::Print(ostream *os){
    // print out the Burst information
    // Inputs:
    //    ostream *os  Pointer to the output stream.
    // Outputs:
    //    none.
    // Return:
    //    none.
/*
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
*/
    *os << "Burst Number: "<< GetEventNumber()<< " Transfers: " 
        << GetTransfers() << endl; 
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSspdTestBeamBurst &p){
    // Standard output operator. See Print
    // Inputs:
    //    ostream                 &os  the output stream.
    //    AliITSspdTestBeamHeader &p   the data to be printed out.
    // Outputs:
    //    none.
    // Return:
    //    ostream &os pointing now to the end of the present stream.

    p.Print(&os);
    return os;
}
//======================================================================
void AliITSspdTestBeamData::Print(ostream *os,Int_t i){
    // print out the the Test Beam Data information
    // Inputs:
    //    ostream *os  Pointer to the output stream.
    // Outputs:
    //    none.
    // Return:
    //    none.
/*
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
*/
  //*os << "Word=" << hex << (fUnion+i)->fBuf[0] << hex << (fUnion+i)->fBuf[1] 
  //               << hex << (fUnion+i)->fBuf[2] << hex << (fUnion+i)->fBuf[3]
    *os << "Word=" << hex << (fUnion+i)->fIBuff
        << dec;
    switch (this->Mode(i)){
    case AliITSTestBeamData::kData :
        *os << " kData chip=" << setw(3) << Chip(i); 
        *os << " Row="        << setw(3) << Row(i);
        *os << " Column="     << setw(3) << Colm(i);
        break;
    case AliITSTestBeamData::kHead :
        *os << " kHead Event Sync =" << EventCounter();
        break;
    case AliITSTestBeamData::kTail  :
        *os << " kTail Transmitted word count =" << TransWordCount();
        break;
    case AliITSTestBeamData::kAbort :
        *os << " kAbort Transmitted word count =" << TransWordCount();
        break;
    default:
        *os << " Unknown Data Type"; 
        break;
    } // end switch
    *os << endl;
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSspdTestBeamData &p){
    // Standard output operator. See Print
    // Inputs:
    //    ostream                 &os  the output stream.
    //    AliITSspdTestBeamHeader &p   the data to be printed out.
    // Outputs:
    //    none.
    // Return:
    //    ostream &os pointing now to the end of the present stream.

    p.Print(&os,0);
    return os;
}
