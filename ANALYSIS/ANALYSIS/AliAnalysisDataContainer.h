#ifndef ALIANALYSISDATACONTAINER_H
#define ALIANALYSISDATACONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysDataContainer - Container of data of arbitrary type deriving
//      from TObject used for analysis. A container must be connected to the 
//      output data slot of a single analysis task (producer) , but also as 
//      input slot for possibly several other tasks (consumers). The connected 
//      slots must enforce the same data type as the container (or a derived type).
//      A container becomes the owner of the contained data once this was produced.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef ROOT_TStopwatch
#include "TStopwatch.h"
#endif

class TClass;
class TFile;
class TObjArray;
class TCollection;
class AliAnalysisTask;
class AliAnalysisDataWrapper;
class TBuffer;

using std::ostream;

class AliAnalysisDataContainer : public TNamed {

public:
enum ENotifyMessage {
   kDeleteData,
   kSaveData,
   kFileChange
};   
enum EAnalysisContainerFlags {
   kPostEventLoop = BIT(14),
   kSpecialOutput = BIT(15),
   kRegisterDataset = BIT(16),
   kExchangeData  = BIT(17),
   kTouchedFlag   = BIT(18)
};     
   AliAnalysisDataContainer();
   AliAnalysisDataContainer(const AliAnalysisDataContainer &cont);
   AliAnalysisDataContainer(const char *name, TClass *type);
   virtual ~AliAnalysisDataContainer();

   // Assignment
   AliAnalysisDataContainer &operator=(const AliAnalysisDataContainer &cont);
   // Getters
   TObject                  *GetData() const      {return fData;}
   const char               *GetFileName() const  {return fFileName.Data();}
   const char               *GetFolderName() const {return fFolderName.Data();}
   TFile                    *GetFile() const      {return fFile;}
   TClass                   *GetType() const;
   AliAnalysisTask          *GetProducer() const  {return fProducer;}
   TObjArray                *GetConsumers() const {return fConsumers;}
   virtual void              GetEntry(Long64_t ientry);
   // Setters
   void                      Reset()              {fData = 0; fDataReady = kFALSE; SetTouched(kFALSE);}
   void                      ResetDataReady()     {fDataReady = kFALSE;}
   virtual Bool_t            SetData(TObject *data, Option_t *option="");
   void                      SetDataOwned(Bool_t flag) {fOwnedData = flag;}
   void                      SetExchange(Bool_t flag) {TObject::SetBit(kExchangeData,flag);}
   void                      SetPostEventLoop(Bool_t flag=kTRUE) {TObject::SetBit(kPostEventLoop,flag);}
   void                      SetSpecialOutput(Bool_t flag=kTRUE) {TObject::SetBit(kSpecialOutput,flag);}
   void                      SetRegisterDataset(Bool_t flag=kTRUE) {TObject::SetBit(kRegisterDataset,flag);}
   void                      SetFileName(const char *filename);
   void                      SetFile(TFile *f) {fFile = f;}
   void                      SetProducer(AliAnalysisTask *prod, Int_t islot);
   void                      SetTouched(Bool_t flag=kTRUE)       {TObject::SetBit(kTouchedFlag,flag);}
   void                      AddConsumer(AliAnalysisTask *cons, Int_t islot);
   void                      DeleteData();
   // Wrapping
   AliAnalysisDataWrapper   *ExportData() const;
   void                      ImportData(AliAnalysisDataWrapper *pack);
   // Container status checking
   Bool_t                    IsDataReady() const  {return fDataReady;}
   Bool_t                    IsExchange() const      {return TObject::TestBit(kExchangeData);}
   Bool_t                    IsPostEventLoop() const {return TObject::TestBit(kPostEventLoop);}
   Bool_t                    IsSpecialOutput() const {return TObject::TestBit(kSpecialOutput);}
   Bool_t                    IsRegisterDataset() const {return TObject::TestBit(kRegisterDataset);}
   Bool_t                    IsTouched() const       {return TObject::TestBit(kTouchedFlag);}
   Bool_t                    IsOwnedData() const  {return fOwnedData;}
   Bool_t                    ClientsExecuted() const;
   Bool_t                    HasConsumers() const {return (fConsumers != 0);}
   Bool_t                    HasProducer() const  {return (fProducer != 0);}
   // Container merging
   virtual Long64_t          Merge(TCollection *list);
   // Send a notify signal to the container
   virtual void              NotifyChange(ENotifyMessage /*type*/) {;}
   // Print connected tasks/status
   void                      PrintContainer(Option_t *option="all", Int_t indent=0) const;

private:
   void                      SetType(TClass *type) {fType = type;}   

protected:
   Bool_t                    fDataReady;  // Flag that data is ready
   Bool_t                    fOwnedData;  // Flag data ownership
   TString                   fFileName;   // File storing the data
   TString                   fFolderName; // Folder name in the output file
   TFile                    *fFile;       //! Opened file
   TObject                  *fData;       // Contained data
   TClass                   *fType;       //! Type of contained data
   AliAnalysisTask          *fProducer;   // Analysis task to which the slot belongs
   TObjArray                *fConsumers;  // List of consumers of the data
   
   ClassDef(AliAnalysisDataContainer,2)  // Class describing a data container for analysis
};

//==============================================================================
//   AliAnalysysDataWrapper - A basic wrapper for exchanging via the network
// the data held by AliAnalysisDataContainer between the master and the client
// in PROOF case. 
//==============================================================================

class AliAnalysisDataWrapper : public TNamed {

public:

enum EAnalysisWrapperFlags {
   kDeleteData = BIT(14)
};     
   AliAnalysisDataWrapper() : TNamed(), fData(NULL) {}
   AliAnalysisDataWrapper(TObject *data);
   AliAnalysisDataWrapper(const AliAnalysisDataWrapper &other) 
                        : TNamed(other), fData(other.fData) {}
   virtual ~AliAnalysisDataWrapper();
   
   // Assignment
   AliAnalysisDataWrapper &operator=(const AliAnalysisDataWrapper &other);

   TObject                  *Data() const {return fData;}
   // Merging
   virtual Long64_t          Merge(TCollection *list);
   void                      SetDeleteData(Bool_t flag=kTRUE) {TObject::SetBit(kDeleteData,flag);}

protected:
   TObject                  *fData;       // Wrapped data

   ClassDef(AliAnalysisDataWrapper, 1) // Data wrapper class for exchange via the net
};

//==============================================================================
//   AliAnalysisFileDescriptor - A simple wrapper for data related to accessing
//      an analysis input file.
//==============================================================================
class AliAnalysisFileDescriptor : public TObject {

public:
   AliAnalysisFileDescriptor();
   AliAnalysisFileDescriptor(const TFile *file);
   AliAnalysisFileDescriptor(const AliAnalysisFileDescriptor &other);
   virtual ~AliAnalysisFileDescriptor();

   // Assignment
   AliAnalysisFileDescriptor &operator=(const AliAnalysisFileDescriptor &other);
   
   void         Done();
   
   const char  *GetGUID() const  {return fGUID;}
   Int_t        GetImage() const {return fImage;}
   const char  *GetLfn() const   {return fLfn;}
   Int_t        GetNreplicas() const {return fNreplicas;}
   Long64_t     GetOpenAt() const {return fOpenedAt;}
   Double_t     GetOpenTime() const {return fOpenTime;}
   const char  *GetPfn() const   {return fPfn;}
   Long64_t     GetReadBytes() const {return fReadBytes;}
   Long64_t     GetSize() const   {return fSize;}
   const char  *GetSE() const     {return fSE;}
   Double_t     GetThroughput() const {return fThroughput;}
   Double_t     GetProcessingTime() const {return fProcessingTime;}
   const char  *GetUrl() const    {return fUrl;}
   Bool_t       IsArchive() const {return fIsArchive;}
   virtual void Print(Option_t *option="") const;
   virtual void SavePrimitive(std::ostream &out, Option_t *option = "");
   void         SetOpenTime(Double_t time) {fOpenTime = time;}
   

protected:   
   TString      fLfn;       // logical file name
   TString      fGUID;      // GUID
   TString      fUrl;       // url for the opened copy
   TString      fPfn;       // physical file name
   TString      fSE;        // Storage element
   Bool_t       fIsArchive; // Archive flag
   Int_t        fImage;     // Image number
   Int_t        fNreplicas; // Number of replicas
   Long64_t     fStartBytes;// Total number of read bytes at start
   Long64_t     fReadBytes; // Number of bytes read
   Long64_t     fSize;      // Size of the file in bytes
   Long64_t     fOpenedAt;  // Absolute value for time when opened
   Double_t     fOpenTime;  // Time elapsed to open file
   Double_t     fProcessingTime; // Processing
   Double_t     fThroughput; // Throughput
   TStopwatch   fTimer;     //! Processing time
   
   ClassDef(AliAnalysisFileDescriptor,1)  // Class describing a a file processed in the analysis
};

#endif
