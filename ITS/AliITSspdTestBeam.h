#include <Rtypes.h>
#ifndef ALIITSSPDTESTBEAM_H
#define ALIITSSPDTESTBEAM_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

#include <TTask.h>

//class ifstream;
class AliITS;
class AliITSspdTestBeamHeader;
class AliITSspdTestBeamTail;
class AliITSspdTestBeamBurst;
class AliITSspdTestBeamData;

class AliITSspdTestBeam : public TTask{
  public:
    AliITSspdTestBeam();
    AliITSspdTestBeam(const Char_t *filename,const Char_t *opt="2002",
                      AliITS *its=0);
    AliITSspdTestBeam(const AliITSspdTestBeam &s):TTask(s){if(this==&s) return;Error("Copy constructor","You are not allowed to make a copy of AliITSspdTestBeam");exit(1);} //Not to be used!
    AliITSspdTestBeam& operator=(AliITSspdTestBeam &s){if(this==&s) return *this;Error("operator=","You are not allowed to make a copy of AliITSspdTestBeam");exit(1);return *this;} //Not to be used!
    virtual ~AliITSspdTestBeam();
    //
    virtual Int_t OpenInputFile(const Char_t *filename,Int_t start=0,
                                Int_t end=-1);
    virtual Int_t Read(Int_t i=0);
    virtual Int_t Decode();
    virtual Int_t GetNumberOfPilots()const{return 3;}
  private:
    void SetTerminationWord(){fTermination=0xffffd9f0;}
    //
    AliITSspdTestBeamHeader  *fRH;    //! Run Header
    AliITSspdTestBeamTail    *fRT;    //! Run Trailer
    Int_t                     fNBrst; //! Number of burts (size of array).
    Int_t                   **fBrstSize; //! Size of each burst for each pilot
    AliITSspdTestBeamBurst  **fBrst;  //! Array of bursts.
    Int_t                   **fNData; //! array of the number of data points
    AliITSspdTestBeamData  ***fData;  //! Data
    AliITSspdTestBeamData  ***fHData; //! pointer to headers
    AliITSspdTestBeamData  ***fTData; //! pointer to Tail and Aborts
    UInt_t     fTermination;  //! Termination word
    Int_t      fNEvents;      // Number of events in file
    Int_t      fBuffSize;     // Read Buffere Size
    UChar_t   *fBuff;         // Read buffer
    AliITS    *fITS;          // Pointer to the ITS.
    Int_t      fNfiles;       // Number of input files to read from
    Int_t      fMaxFiles;     // The size of the pointer array fFiles.
    ifstream **fFiles;        //! Array of Pointer to the input streams
    Int_t     *fNeventsStart; // Starting event number for each file
    Int_t     *fNeventsEnd;   // Ending number of events for each file.

    ClassDef(AliITSspdTestBeam,1) // Task to read SPD test beam data
};
#endif
//======================================================================
#ifndef ALIITSTESTBEAMDATA_H
#define ALIITSTESTBEAMDATA_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

// Pure virtual class, no data
class AliITSTestBeamData{
  public:
    AliITSTestBeamData(){};
    virtual ~AliITSTestBeamData(){};
    virtual Int_t SizeOf()const{return 0;}
    enum {kData,kHead,kTail,kAbort,kFail};
  private:
};

#endif
//======================================================================
#ifndef ALIITSSPDTESTBEAMHEADER_H
#define ALIITSSPDTESTBEAMHEADER_H

#include <Riostream.h>

//class ostream;

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

class AliITSspdTestBeamHeader : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamHeader(){};
    virtual ~AliITSspdTestBeamHeader(){};
    virtual void Print(ostream *os)const;
    virtual Int_t GetNumberOfEvents()const{return fUnion.fHead.fNEvents;};
    virtual Int_t GetBurstSize(Int_t i=0)const{return fUnion.fHead.fBuffSize[i];};
    virtual Int_t SizeOf()const{return sizeof(AliITSspdTestBeamHeader);};
  private:
    union headder{
         struct {
            // Double_t must be 8 bytes long, Char_t 1 byte long,
            // Int_t 4 bytes long. 8*1+4*26+1*(2*12+3*20+3*44+3*8192)=24904 by
            Double_t fVersion; // Version number
            Char_t   fDate[12];// Date Field
            Char_t   fTime[12];// Time Field
            UInt_t   fBuffSize[3]; // Buffer Sizes
            UInt_t   fTestPulse;   // Test Pulse flag
            UInt_t   fTrigger[3];  // Array of triggers
            UInt_t   fTriggerMode; // Trigger mode flag
            UInt_t   fBurstSize;   // Burst Size
            UInt_t   fNEvents;     // Number of event to write
            UInt_t   fRes[16];            // not sure
            UChar_t  fMBDACS[3][20];      // not sure
            UChar_t  fA1DACS[3][44];      // not sure
            UChar_t  fA12Matrix[3][8192]; // not sure
        } fHead;
        UChar_t fBuf[24904]; // Equivalent char buffer
    } fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamHeader &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMTAIL_H
#define ALIITSSPDTESTBEAMTAIL_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

#include <Riostream.h>

//class ostream;

class AliITSspdTestBeamTail : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamTail(){};
    virtual ~AliITSspdTestBeamTail(){};
    virtual Int_t SizeOf()const{return sizeof(AliITSspdTestBeamTail);}
    virtual void Print(ostream *os)const;
  private:
    union tail{
        struct {
            // Char_t 1 byte long,
            // UInt_t 4 bytes long. size = 4+4+12*1+12*1 = 32 bytes
            UInt_t   fEvents;  // number of events written
            UInt_t   fTermMode;// Termination flag
            Char_t   fDate[12];// Date Field
            Char_t   fTime[12];// Time Field
        } fTail;
        UChar_t fBuf[32]; //Equivalent char buffer
    } fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamTail &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMBURST_H
#define ALIITSSPDTESTBEAMBURST_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

#include <Riostream.h>

//class ostream;

class AliITSspdTestBeamBurst : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamBurst(){};
    virtual ~AliITSspdTestBeamBurst(){};
    virtual Int_t SizeOf()const{return sizeof(AliITSspdTestBeamBurst);}
    virtual void Print(ostream *os)const;
  private:
    union tail{
        struct {
            // UInt_t 4 bytes long.
            UInt_t   fNumber;   // Burst Number
            UInt_t   fTransfers;// Number of Transfers
        } fBrst;
        ULong64_t fBuf;  // a strictly 64 bit long unsinged int
    } fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamBurst &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMDATA_H
#define ALIITSSPDTESTBEAMDATA_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

#include <Riostream.h>

//class ostream;

class AliITSspdTestBeamData : public AliITSTestBeamData {
  public:
    //
    AliITSspdTestBeamData(){};
    //AliITSspdTestBeamData(UInt_t *i){SetAddress(i);}
    virtual Int_t SizeOf()const{return sizeof(AliITSspdTestBeamData);}
    virtual ~AliITSspdTestBeamData(){};
    //
    virtual Bool_t IsHeader()const{return (fUnion.fDataH.fFlag==2||fUnion.fDataH.fFlag==3);}
    virtual Bool_t IsData()const{return (fUnion.fDataD.fFlag>3)&&(fUnion.fDataD.fFlag<8);}
    virtual Bool_t IsTrailer()const{return (fUnion.fDataT.fFlag==0);}
    virtual Bool_t IsAbort()const{return (fUnion.fDataA.fFlag==1);}
    virtual Int_t  Mode()const{if(IsData()) return kData;else if(IsHeader()) return kHead;else if(IsTrailer()) return kTail;else if(IsAbort()) return kAbort;else return kFail;}
    virtual Int_t  TransWordCount()const{if(IsTrailer()||IsAbort()) return fUnion.fDataT.fTrans;else return -1;}
    virtual Int_t  EventCounter()const{if(IsHeader()) return fUnion.fDataH.fEventSync;else return -1;}
    virtual Int_t Chip()const{if(IsData()) return fUnion.fDataD.fChip;else return -1;}
    virtual Int_t Row() const{if(IsData()) return fUnion.fDataD.fRow; else return -1;}
    virtual Int_t Colm()const{if(IsData()) return fUnion.fDataD.fColm;else return -1;}
    virtual void  Dataconst(Int_t &ch,Int_t &rw,Int_t &cl)const{ch=Chip();rw=Row();cl=Colm();}
    virtual void Print(ostream *os)const;
  private:
    union data{
        struct {
            unsigned fDate:12;    // Not Used
            unsigned fEventSync:6;// Event syncronization counter
            unsigned fFlag:3;     // =2or3 Header
            unsigned fPadding:12; // Not used
        } fDataH;
        struct {
            unsigned fColm:5;     // Pixel (Hit) Column Address
            unsigned fRow:8;      // Pixel Row Adress
            unsigned fChip:4;     // Pixel Chip Address
            unsigned fFlag:3;     // =4or5or6or7 Data
            unsigned fPadding:12; // Not used
        } fDataD;
        struct {
            unsigned fTrans:17;   // Transmitted word count
            unsigned fFlag:3;     // =0 Trialer
            unsigned fPadding:12; // Not used
        } fDataT;
        struct {
            unsigned fTrans:17;   // Transmitted word count
            unsigned fFlag:3;     // =1 Abort
            unsigned fPadding:12; // Not used
        } fDataA;
        UChar_t fBuf[4]; // Equivalent char buffer
    } fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamData &source);
#endif
