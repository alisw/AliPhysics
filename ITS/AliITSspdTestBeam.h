#include <Rtypes.h>
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
    virtual Int_t SizeOf(){return 0;}
    virtual Double_t Swap(Double_t a){union{Double_t b;UChar_t c[8];}d;d.b=a;
                                      Swapit(8,d.c);return d.b;};
    virtual UInt_t Swap(UInt_t a){union{UInt_t b;UChar_t c[4];}d;d.b=a;
                                  Swapit(4,d.c);return d.b;};
    static void Swapit(Int_t i,UChar_t *a);
    enum {kData,kHead,kTail,kAbort,kFail};
  private:
};

#endif
//======================================================================
#ifndef ALIITSSPDTESTBEAMHEADER_H
#define ALIITSSPDTESTBEAMHEADER_H

//class ostream;

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

class AliITSspdTestBeamHeader : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamHeader(){fUnion=0;};
    AliITSspdTestBeamHeader(UChar_t *f){SetBuffer(f);}
    virtual ~AliITSspdTestBeamHeader(){fUnion=0;};
    virtual void Print(ostream *os);
    virtual Int_t GetNumberOfEvents(){return Swap(fUnion->fHead.fNEvents);};
    virtual Int_t GetBuffSize(Int_t i=0){
        return Swap(fUnion->fHead.fBuffSize[i]);};
    virtual Int_t GetTestPulse(){return Swap(fUnion->fHead.fTestPulse);};
    virtual Int_t GetTrigger(Int_t i=0){
        return Swap(fUnion->fHead.fTrigger[i]);};
    virtual Int_t GetTriggerMode(){return Swap(fUnion->fHead.fTriggerMode);};
    virtual Int_t GetBurstSize(){return Swap(fUnion->fHead.fBurstSize);};
    virtual Double_t GetVersion(){return Swap(fUnion->fHead.fVersion);}
    virtual Int_t SizeOf(){return sizeof(union headder);};
    virtual void SetBuffer(UChar_t *f){
                            fUnion=(AliITSspdTestBeamHeader::headder*)f;}
    virtual UChar_t *GetBuffer(){return fUnion->fBuf;}
    virtual Char_t *GetDate(){return fUnion->fHead.fDate;}
    virtual Char_t *GetTime(){return fUnion->fHead.fTime;}
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
            UInt_t   fRes[16];            //
            UChar_t  fMBDACS[3][20];      //
            UChar_t  fA1DACS[3][44];      //
            UChar_t  fA12Matrix[3][8192]; //
        } fHead;
        UChar_t fBuf[24904];
    } *fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamHeader &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMTAIL_H
#define ALIITSSPDTESTBEAMTAIL_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
//class ostream;

class AliITSspdTestBeamTail : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamTail(){fUnion=0;};
    AliITSspdTestBeamTail(UChar_t *f){SetBuffer(f);}
    virtual ~AliITSspdTestBeamTail(){fUnion=0;};
    virtual Int_t GetNumberOfEvents(){return Swap(fUnion->fTail.fEvents);};
    virtual Int_t GetTermMode(){return Swap(fUnion->fTail.fTermMode);};
    virtual Int_t SizeOf(){return sizeof(union tail);}
    virtual void Print(ostream *os);
    virtual void SetBuffer(UChar_t *f){fUnion=(AliITSspdTestBeamTail::tail*)f;}
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
        UChar_t fBuf[32];
    } *fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamTail &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMBURST_H
#define ALIITSSPDTESTBEAMBURST_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
//class ostream;

class AliITSspdTestBeamBurst : public AliITSTestBeamData{
  public:
    AliITSspdTestBeamBurst(){fUnion=0;};
    AliITSspdTestBeamBurst(UChar_t *f){SetBuffer(f);}
    virtual ~AliITSspdTestBeamBurst(){fUnion=0;};
    virtual Int_t GetEventNumber(){return Swap(fUnion->fBrst.fNumber);};
    virtual Int_t GetTransfers(){return Swap(fUnion->fBrst.fTransfers);};
    virtual Int_t SizeOf(){return sizeof(union tail);}
    virtual void Print(ostream *os);
    virtual void SetBuffer(UChar_t *f){
                                  fUnion=(AliITSspdTestBeamBurst::tail*)f;}
  private:
    union tail{
        struct {
            // UInt_t 4 bytes long.
            UInt_t   fNumber;   // Burst Number
            UInt_t   fTransfers;// Number of Transfers
        } fBrst;
        ULong64_t fBuf;  // a strictly 64 bit long unsinged int
    } *fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamBurst &source);
#endif
//----------------------------------------------------------------------
#ifndef ALIITSSPDTESTBEAMDATA_H
#define ALIITSSPDTESTBEAMDATA_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
//class ostream;

class AliITSspdTestBeamData : public AliITSTestBeamData {
  public:
    //
    AliITSspdTestBeamData(){fUnion=0;};
    AliITSspdTestBeamData(UChar_t *f){SetBuffer(f);}
    virtual Int_t SizeOf(){return sizeof(union data);}
    virtual ~AliITSspdTestBeamData(){};
    //
    virtual Bool_t IsHeader(Int_t i=0){return ((fUnion+i)->fDataH.fFlag==2||
                                      (fUnion+i)->fDataH.fFlag==3);}
    virtual Bool_t IsData(Int_t i=0){return ((fUnion+i)->fDataD.fFlag>3)&&
                                   ((fUnion+i)->fDataD.fFlag<8);}
    virtual Bool_t IsTrailer(Int_t i=0){return ((fUnion+i)->fDataT.fFlag==0);}
    virtual Bool_t IsAbort(Int_t i=0){return ((fUnion+i)->fDataA.fFlag==1);}
    virtual Int_t  Mode(Int_t i=0){if(IsData(i)) return kData;
                          else if(IsHeader(i)) return kHead;
                          else if(IsTrailer(i)) return kTail;
                          else if(IsAbort(i)) return kAbort;
                          else return kFail;}
    virtual Int_t  TransWordCount(){if(IsTrailer()||IsAbort()) 
                                    return fUnion->fDataT.fTrans;
                                    else return -1;}
    virtual Int_t  EventCounter(){
        if(IsHeader()) return fUnion->fDataH.fEventSync;else return -1;}
    virtual Int_t Chip(Int_t i=0){if(IsData()) 
        return (fUnion+i)->fDataD.fChip; else return -1;}
    virtual Int_t Row(Int_t i=0) {if(IsData()) 
        return (fUnion+i)->fDataD.fRow; else return -1;}
    virtual Int_t Colm(Int_t i=0){if(IsData()) 
        return (fUnion+i)->fDataD.fColm;else return -1;}
    virtual void  Data(Int_t &ch,Int_t &rw,Int_t &cl,Int_t i=0){
        ch=Chip(i);rw=Row(i);cl=Colm(i);}
    virtual void Print(ostream *os,Int_t i=0);
    virtual void SetBuffer(UChar_t *f){fUnion=(AliITSspdTestBeamData::data*)f;}
  private:
    union data{
        struct { // Definingthe fDataH bit field
            unsigned fDate:11;    // Not Used
            unsigned fEventSync:6;// Event syncronization counter
            unsigned fFlag:3;     // =2or3 Header
            unsigned fPadding:12; // Not used
        } fDataH;
        struct { // Definingthe fDataD bit field
            unsigned fColm:5;     // Pixel (Hit) Column Address
            unsigned fRow:8;      // Pixel Row Adress
            unsigned fChip:4;     // Pixel Chip Address
            unsigned fFlag:3;     // =4or5or6or7 Data
            unsigned fPadding:12; // Not used
        } fDataD;
        struct { // Definingthe fDataT bit field
            unsigned fTrans:17;   // Transmitted word count
            unsigned fFlag:3;     // =0 Trialer
            unsigned fPadding:12; // Not used
        } fDataT;
        struct { // Definingthe fDataA bit field
            unsigned fTrans:17;   // Transmitted word count
            unsigned fFlag:3;     // =1 Abort
            unsigned fPadding:12; // Not used
        } fDataA;
        UChar_t fBuf[4];
        UInt_t  fIBuff;
    } *fUnion;
};
ostream &operator<<(ostream &os,AliITSspdTestBeamData &source);
#endif
//======================================================================
#ifndef ALIITSSPDTESTBEAM_H
#define ALIITSSPDTESTBEAM_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

#include <TTask.h>

//class ostream;
class AliITS;
class AliITSspdTestBeamHeader;
class AliITSspdTestBeamTail;
class AliITSspdTestBeamBurst;
class AliITSspdTestBeamData;
class AliLoader;

class AliITSspdTestBeam : public TTask{
  public:
    AliITSspdTestBeam();
    AliITSspdTestBeam(const Char_t *filename,const Char_t *opt="2002",
                      AliITS *its=0);
    virtual ~AliITSspdTestBeam();
    //
    virtual Int_t OpenInputFile(const Char_t *filename,Int_t start=0,
                                Int_t end=-1);
    virtual Int_t Read(Int_t i=0);
    virtual Int_t Decode();
    virtual Int_t DecodeModule(Int_t pilot,Int_t Chip);
    virtual Int_t DecodeColumn(Int_t pilot,Int_t Chip,Int_t Colm);
    virtual void Digitize(Int_t evnt);
    virtual void SetLoader(AliLoader *loader) {fLoader = loader;}
    virtual void SetITS(AliITS *its) {fITS = its;}
    virtual Int_t GetNumberOfPilots(){return 3;}
    virtual Int_t GetNumberOfEvents(){return fNEvents;}
    virtual void PrintHeadder(ostream *os = &cout);
    virtual void PrintTrailer(ostream *os = &cout);
    virtual void PrintBurstInfo(Int_t i,ostream *os=&cout);
    virtual void PrintEventData(Int_t i,ostream *os=&cout);
    virtual void PrintEventHead(Int_t i,ostream *os=&cout);
    virtual void PrintEventTail(Int_t i,ostream *os=&cout);
  private:
    void SetTerminationWord(){fTermination=0xf0d8ffff;}
    void DeletefBrst();
    void DeletefNData();
    void DeletefData();
    void DeletefHData();
    void DeletefTData();
    void DeletefFiles();
    void DeletefBrstSize();
    void DeletefBuff(){delete[] fBuff;};
    void DeletefITS(){fITS=0;}
    void DeletefNeventsStart(){delete[] fNeventsStart;}
    void DeletefNeventsEnd(){delete[] fNeventsEnd;}
    //
    AliITSspdTestBeamHeader   fRH;    //! Run Header
    AliITSspdTestBeamTail     fRT;    //! Run Trailer
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
    Int_t      fVersion;      // Test Beam Version number
    AliITS    *fITS;          // Pointer to the ITS.
    AliLoader *fLoader;       //! Pointer to AliITSLoader.
    Int_t      fNfiles;       // Number of input files to read from
    Int_t      fMaxFiles;     // The size of the pointer array fFiles.
    ifstream **fFiles;        //! Array of Pointer to the input streams
    Int_t     *fNeventsStart; // Starting event number for each file
    Int_t     *fNeventsEnd;   // Ending number of events for each file.
    ClassDef(AliITSspdTestBeam,1) // Task to read SPD test beam data
};
#endif
