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

// A set of classes/routines that can read the SPD test beam data of 2002
// and create AliITSdigits. The posibility to use these routines to do the
// same for later and other detectors has yet to be demonstrated. At present
// there remains a bug in that the TreeE of event headders isn't created
// properly. See the macro AliITSspdTestBeam2Digits.C. The geometry from
// the class AliITSvSPD002 must be read in, one way or the other, so that
// the Geometry transoformation class AliITSgeom will prpoerly be inilized.

#include <stdlib.h>
#include <stddef.h>
#include <iomanip>
#include <Riostream.h>
#include <fstream>
#include <TArrayI.h>

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

    fRH           = 0;
    fRT           = 0;
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

    fRH           = 0;
    fRT           = 0;
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
    Int_t i,np;

    np = GetNumberOfPilots();
    fRH = 0;  // Just a Pointer into fBuff.
    fRT = 0;  // Just a Pointer into fBuff.
    if(fBrst){delete[] fBrst; fBrst = 0;}
    if(fNData)for(i=0;i<np;i++){
        if(fNData[i]) delete[] fNData[i];
    } // end if
    if(fNData) {delete[] fNData; fNData = 0;}
    if(fData)for(i=0;i<np;i++){
        if(fData[i]) delete[] fData[i];
    } // end if
    if(fData) delete[] fData;
    fData = 0;
    if(fHData)for(i=0;i<np;i++){
        if(fHData[i]) delete[] fHData[i];
    } // end if
    if(fHData) delete[] fHData;
    fHData = 0;
    if(fTData)for(i=0;i<np;i++){
        if(fTData[i]) delete[] fTData[i];
    } // end if
    if(fTData) delete[] fTData;
    fTData = 0;
    for(i=0;i<fMaxFiles;i++){
        if(fFiles[i]!=0) delete fFiles[i];
    } // end for i
    if(fBrstSize)for(i=0;i<np;i++){
        if(fBrstSize[i]) delete[] fBrstSize[i];
    } // end if
    if(fBrstSize) {delete[] fBrstSize; fBrstSize = 0;}
    delete[] fBuff;
    fITS = 0; //delete fITS;
    delete[] fFiles;
    delete[] fNeventsStart;
    delete[] fNeventsEnd;
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
        fNeventsStart = new Int_t[fMaxFiles];
        fNeventsEnd   = new Int_t[fMaxFiles];
    } // end if
    if(fNfiles==fMaxFiles){// Need to expand array of open files.
        ifstream **tmp = new ifstream*[fMaxFiles];
        TArrayI st(fMaxFiles);
	TArrayI  en(fMaxFiles);
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
	delete [] tmp;
    } // end if
    // Open file
#ifndef __DECCXX
    fFiles[fNfiles] = new ifstream(filename,ios::in|ios::binary);
#else
    fFiles[fNfiles] = new ifstream(filename,ios::in);
#endif
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
    ivnt = new Int_t[np];
    for(i=0;i<np;i++) ivnt[i] = 0;
    fRH  = (AliITSspdTestBeamHeader*) &(fBuff[0]); // Sets up the Run Header.
    fRT  = (AliITSspdTestBeamTail*)&(fBuff[fBuffSize-fRT->SizeOf()]);
    // Check termination
    tr   = (UInt_t*) &(fBuff[fBuffSize-(fRT->SizeOf())-sizeof(UInt_t)]);
    if(!(*tr==fTermination)){
        cout << "Error Termination word not found at "<<tr<<endl;
        exit(-1);
    } // end if
    if(fNData)for(i=0;i<np;i++){
        if(fNData[i]) delete[] fNData[i];
    } // end if
    if(fNData) {delete[] fNData; fNData = 0;}
    if(fData)for(i=0;i<np;i++){
        if(fData[i]) delete[] fData[i];
    } // end if
    if(fData) {delete[] fData; fData = 0;}
    if(fHData)for(i=0;i<np;i++){
        if(fHData[i]) delete[] fHData[i];
    } // end if
    if(fHData) {delete[] fHData; fHData = 0;}
    if(fTData)for(i=0;i<np;i++){
        if(fTData[i]) delete[] fTData[i];
    } // end if
    if(fTData) {delete[] fTData; fTData = 0;}
    if(fBrstSize)for(i=0;i<np;i++){
        if(fBrstSize[i]) delete[] fBrstSize[i];
    } // end if
    if(fBrstSize) {delete[] fBrstSize; fBrstSize = 0;}
    fNEvents  = fRH->GetNumberOfEvents();
    fNBrst    = fNEvents/fRH->GetBurstSize();
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
    if(fBrst){delete[] fBrst; fBrst = 0;}
    size      = fRH->SizeOf();
    u.bt      = &fBuff[size];
    //
    for(iburst=0;(*(u.wd) != fTermination)&&(u.wd<tr);iburst++){ // loop over Bursts
        b   = (AliITSspdTestBeamBurst *) u.wd;
        fBrst[iburst] = b;
        u.bt += b->SizeOf(); // increment wd byte wise
        for(ip=0;ip<np;ip++){  // loop over pilots
            // Get size of data stored for this pilot.
            fBrstSize[ip][iburst] = (UInt_t) u.wd;
            u.bt += sizeof(UInt_t); // increment wd byte wise
            for(i=0;i<fBrstSize[ip][iburst];i++){ // loop over data
                d = (AliITSspdTestBeamData *) u.wd;
                switch (d->Mode()){
                case AliITSTestBeamData::kData :
                    fNData[ip][ivnt[ip]]++;
                    // set pointer to first data member
                    if(fData[ip][ivnt[ip]] == 0 ) fData[ip][ivnt[ip]] = d;
                    break;
                case AliITSTestBeamData::kHead :
                    fNData[ip][ivnt[ip]]   = 0;
                    fData[ip][ivnt[ip]]    = 0;
                    fHData[ip][ivnt[ip]++] = d;
                    break;
                case AliITSTestBeamData::kTail  : 
                case AliITSTestBeamData::kAbort :
                    fTData[ip][ivnt[ip]] = d;
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
//============================================================================
void AliITSspdTestBeamHeader::Print(ostream *os)const{
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
    *os<<"Version: "<<fUnion.fHead.fVersion<<" Written: "<<fUnion.fHead.fDate;
    *os<<" " << fUnion.fHead.fTime << endl;
    *os<<"Buffer Size [0], [1], [2]: " << fUnion.fHead.fBuffSize[0] << ",";
    *os<<fUnion.fHead.fBuffSize[1]<<"," <<fUnion.fHead.fBuffSize[2] << endl;
    *os<<"Test Pulse: " << fUnion.fHead.fTestPulse << endl;
    *os<<"Trigger Mode General, [0], [1], [2]: " << fUnion.fHead.fTriggerMode;
    *os<<","<<fUnion.fHead.fTrigger[0]<<","<<fUnion.fHead.fTrigger[1]<< ",";
    *os<<fUnion.fHead.fTrigger[2] << endl;
    *os<<"Number of Events: " << fUnion.fHead.fNEvents << " Burst Size: ";
    *os<<fUnion.fHead.fBurstSize << endl;
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
void AliITSspdTestBeamTail::Print(ostream *os)const{
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
    *os << "Number of Events: "<< fUnion.fTail.fEvents << " Written: "
        << fUnion.fTail.fDate;
    *os << " " << fUnion.fTail.fTime << endl;
    *os <<"Termination Flag: " << fUnion.fTail.fTermMode << endl; 
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
void AliITSspdTestBeamBurst::Print(ostream *os)const{
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
    *os << "Burst Number: "<< fUnion.fBrst.fNumber << " Transfers: " 
        << fUnion.fBrst.fTransfers << endl; 
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
void AliITSspdTestBeamData::Print(ostream *os)const{
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
    *os << "Word=" << hex << fUnion.fBuf[0] << hex << fUnion.fBuf[1] 
                   << hex << fUnion.fBuf[2] << hex << fUnion.fBuf[3] << dec;
    switch (this->Mode()){
    case AliITSTestBeamData::kData :
        *os << " kData chip=" << setw(3) << fUnion.fDataD.fChip; 
        *os << " Row="        << setw(3) << fUnion.fDataD.fRow;
        *os << " Column="     << setw(3) << fUnion.fDataD.fColm;
        break;
    case AliITSTestBeamData::kHead :
        *os << " kHead Event Sync =" << fUnion.fDataH.fEventSync;
        break;
    case AliITSTestBeamData::kTail  :
        *os << " kTail Transmitted word count =" << fUnion.fDataT.fTrans;
        break;
    case AliITSTestBeamData::kAbort :
        *os << " kAbort Transmitted word count =" << fUnion.fDataA.fTrans;
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

    p.Print(&os);
    return os;
}
