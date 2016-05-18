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
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysDataContainer - Container of data of arbitrary type deriving
//      from TObject used for analysis. A container must be connected to the 
//      output data slot of a single analysis task (producer) , but also as 
//      input slot for possibly several other tasks (consumers). The connected 
//      slots must enforce the same data type as the container (or a derived type).
//      A container becomes the owner of the contained data once this was produced.
//
// Containers should be defined by the analysis module using:
//
//   AliAnalysisModule::AddContainer(const char *name, TClass *type);
//
// A container should be connected to a producer:

//   AliAnalysisModule::ConnectOutput(AliAnalysisTask *task,
//                                    AliAnalysisDataContainer *cont)
// and to its consumers:
//
//   AliAnalysisModule::ConnectInput(AliAnalysisTask *task, Int_t islot,
//                                   AliAnalysisDataContainer *cont)
//
// The container will create an implicit connection between the producer task 
// and all consumers, which will become sub-tasks of the producer.
//
//==============================================================================

#include <Riostream.h>
#include <TMethodCall.h>

#include <TClass.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TROOT.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisTask.h"
#include "AliLog.h"
#include <TBuffer.h>

using std::endl;
using std::cout;
using std::ios;
using std::setiosflags;
using std::setprecision;
ClassImp(AliAnalysisDataContainer)

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer() : TNamed(),
                          fDataReady(kFALSE),
                          fOwnedData(kFALSE),
                          fFileName(),
                          fFolderName(),
                          fFile(NULL),
                          fData(NULL),
                          fType(NULL),
                          fProducer(NULL),
                          fConsumers(NULL)
{
// Dummy ctor.
}

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer(const char *name, TClass *type)
                         :TNamed(name,""),
                          fDataReady(kFALSE),
                          fOwnedData(kFALSE),
                          fFileName(),
                          fFolderName(),
                          fFile(NULL),
                          fData(NULL),
                          fType(type),
                          fProducer(NULL),
                          fConsumers(NULL)
{
// Default constructor.
   SetTitle(fType->GetName());
}

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer(const AliAnalysisDataContainer &cont)
                         :TNamed(cont),
                          fDataReady(cont.fDataReady),
                          fOwnedData(kFALSE),
                          fFileName(cont.fFileName),
                          fFolderName(cont.fFolderName),
                          fFile(NULL),
                          fData(cont.fData),
                          fType(NULL),
                          fProducer(cont.fProducer),
                          fConsumers(NULL)
{
// Copy ctor.
   GetType();
   if (cont.fConsumers) {
      fConsumers = new TObjArray(2);
      Int_t ncons = cont.fConsumers->GetEntriesFast();
      for (Int_t i=0; i<ncons; i++) fConsumers->Add(cont.fConsumers->At(i));
   }   
}

//______________________________________________________________________________
AliAnalysisDataContainer::~AliAnalysisDataContainer()
{
// Destructor. Deletes data ! (What happens if data is a container ???)
   if (fData && fOwnedData) delete fData;
   if (fConsumers) delete fConsumers;
}

//______________________________________________________________________________
AliAnalysisDataContainer &AliAnalysisDataContainer::operator=(const AliAnalysisDataContainer &cont)
{
// Assignment.
   if (&cont != this) {
      TNamed::operator=(cont);
      fDataReady = cont.fDataReady;
      fOwnedData = kFALSE;
      fFileName = cont.fFileName;
      fFolderName = cont.fFolderName;
      fFile = NULL;
      fData = cont.fData;
      GetType();
      fProducer = cont.fProducer;
      if (cont.fConsumers) {
         fConsumers = new TObjArray(2);
         Int_t ncons = cont.fConsumers->GetEntriesFast();
         for (Int_t i=0; i<ncons; i++) fConsumers->Add(cont.fConsumers->At(i));
      }   
   }   
   return *this;
}      

//______________________________________________________________________________
void AliAnalysisDataContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliAnalysisDataContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AliAnalysisDataContainer::Class(),this);
      this->SetDataOwned(kTRUE);
   } else {
      R__b.WriteClassBuffer(AliAnalysisDataContainer::Class(),this);
   }
}

//______________________________________________________________________________
void AliAnalysisDataContainer::AddConsumer(AliAnalysisTask *consumer, Int_t islot)
{
// Add a consumer for contained data;
   AliAnalysisDataSlot *slot = consumer->GetInputSlot(islot);
   if (!slot || !slot->GetType()) {
     cout<<"Consumer task "<< consumer->GetName()<<" does not have an input/type #"<<islot<<endl;
     //AliError(Form("Consumer task %s does not have an input #%i", consumer->GetName(),islot));
      return;
   }
   if (!slot->GetType()->InheritsFrom(GetType())) {
     cout<<"Data type "<<slot->GetTitle()<<" for input slot "<<islot<<" of task "<<consumer->GetName()<<" does not match container type "<<GetTitle()<<endl;  
     //AliError(Form("Data type %s for input slot %i of task %s does not match container type %s", slot->GetType()->GetName(),islot,consumer->GetName(),fType->GetName()));
      return;
   }   

   if (!fConsumers) fConsumers = new TObjArray(2);   
   fConsumers->Add(consumer);
   // Add the consumer task to the list of task of the producer
   if (fProducer && !fProducer->GetListOfTasks()->FindObject(consumer)) 
      fProducer->Add(consumer);
}      

//______________________________________________________________________________
Bool_t AliAnalysisDataContainer::ClientsExecuted() const
{
// Check if all client tasks have executed.
   TIter next(fConsumers);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) {
      if (!task->HasExecuted()) return kFALSE;
   }
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisDataContainer::DeleteData()
{
// Delete data if not needed anymore.
   if (!fDataReady || !ClientsExecuted()) {
     cout<<"Data not ready or not all clients of container "<<GetName()<<" executed. Data not deleted."<<endl;
     //AliWarning(Form("Data not ready or not all clients of container %s executed. Data not deleted.", GetName()));
      return;
   }
   if (!fOwnedData) {
     cout<<"Data not owned by container "<<GetName()<<". Not deleted."<<endl;
     //AliWarning(Form("Data not owned by container %s. Not deleted.", GetName()));
      return;
   }
   delete fData;
   fData = 0;
   fDataReady = kFALSE;
}   

//______________________________________________________________________________
TClass *AliAnalysisDataContainer::GetType() const
{
// Get class type for this slot.
   AliAnalysisDataContainer *cont = (AliAnalysisDataContainer*)this;
   if (!fType) cont->SetType(gROOT->GetClass(fTitle.Data()));
   if (!fType) printf("AliAnalysisDataContainer: Unknown class: %s\n", GetTitle());
   return fType;
}

//______________________________________________________________________________
void AliAnalysisDataContainer::GetEntry(Long64_t ientry)
{
// If data is ready and derives from TTree or from TBranch, this will get the
// requested entry in memory if not already loaded.
   if (!fDataReady || !GetType()) return;
   Bool_t istree = fType->InheritsFrom(TTree::Class());
   if (istree) {
      TTree *tree = (TTree*)fData;
      if (tree->GetReadEntry() != ientry) tree->GetEntry(ientry);
      return;
   }   
   Bool_t isbranch = fType->InheritsFrom(TBranch::Class());
   if (isbranch) {
      TBranch *branch = (TBranch*)fData;
      if (branch->GetReadEntry() != ientry) branch->GetEntry(ientry);
      return;
   }   
}   

//______________________________________________________________________________
Long64_t AliAnalysisDataContainer::Merge(TCollection *list)
{
// Merge a list of containers with this one. Containers in the list must have
// data of the same type.
   if (!list || !fData) return 0;
   AliInfo(Form("Merging %d containers %s\n", list->GetEntries(), GetName()));
   TMethodCall callEnv;
   if (fData->IsA())
      callEnv.InitWithPrototype(fData->IsA(), "Merge", "TCollection*");
   if (!callEnv.IsValid() && !list->IsEmpty()) {
      cout << "No merge interface for data stored by " << GetName() << ". Merging not possible !" << endl;
      return 1;
   }

   if (list->IsEmpty()) return 1;

   TIter next(list);
   AliAnalysisDataContainer *cont;
   // Make a list where to temporary store the data to be merged.
   TList *collectionData = new TList();
   Int_t count = 0; // object counter
   while ((cont=(AliAnalysisDataContainer*)next())) {
      TObject *data = cont->GetData();
      if (!data) continue;
      if (strcmp(cont->GetName(), GetName())) {
         cout << "Not merging containers with different names !" << endl;
         continue;
      }
      AliInfo(Form(" ... merging object %s\n", data->GetName()));
      collectionData->Add(data);
      count++;
   }
   callEnv.SetParam((Long_t) collectionData);
   callEnv.Execute(fData);
   delete collectionData;

   return count+1;
}

//______________________________________________________________________________
void AliAnalysisDataContainer::PrintContainer(Option_t *option, Int_t indent) const
{
// Print info about this container.
   TString ind;
   for (Int_t i=0; i<indent; i++) ind += " ";
   TString opt(option);
   opt.ToLower();
   TString ctype = "Container";
   if (IsExchange()) ctype = "Exchange container";
   Bool_t dep = (opt.Contains("dep"))?kTRUE:kFALSE;
   if (!dep) {
      if (IsPostEventLoop()) printf("%s%s: %s     DATA TYPE: %s POST_LOOP task\n", ind.Data(), ctype.Data(), GetName(), GetTitle());
      else printf("%s%s: %s     DATA TYPE: %s\n", ind.Data(), ctype.Data(), GetName(), GetTitle());
      if (!fProducer) 
//         printf("%s = Data producer: task %s\n",ind.Data(),fProducer->GetName());
//      else
         printf("%s= Not connected to a data producer\n",ind.Data());
      if (fConsumers && fConsumers->GetEntriesFast())
         printf("%s = Client tasks indented below:\n", ind.Data());
   }
   if (!IsExchange()) {
      if (!fFolderName.IsNull())
        printf("%s = Filename: %s  folder: %s\n", ind.Data(),fFileName.Data(), fFolderName.Data());
      else
        if (!fFileName.IsNull()) printf("%s = Filename: %s\n", ind.Data(),fFileName.Data());
   }     
   ((AliAnalysisDataContainer*)this)->SetTouched(kTRUE);  
   TIter next(fConsumers);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) task->PrintTask(option, indent+3);
}   

//______________________________________________________________________________
Bool_t AliAnalysisDataContainer::SetData(TObject *data, Option_t *)
{
// Set the data as READY only if it was published by the producer.
   // If there is no producer declared, this is a top level container.
   AliAnalysisTask *task;
   Bool_t init = kFALSE;
   Int_t i, nc;
   if (!fProducer) {
      if (data != fData) init = kTRUE;
      fData = data;
      fDataReady = kTRUE;
      if (fConsumers) {
         nc = fConsumers->GetEntriesFast();
         for (i=0; i<nc; i++) {
            task = (AliAnalysisTask*)fConsumers->At(i);
            task->CheckNotify(init);
         }
      }      
      return kTRUE;
   }
   // Check if it is the producer who published the data     
   if (fProducer->GetPublishedData()==data) {
      fData = data;
      fDataReady = kTRUE;
      if (fConsumers) {
         nc = fConsumers->GetEntriesFast();
         for (i=0; i<nc; i++) {
            task = (AliAnalysisTask*)fConsumers->At(i);
            task->CheckNotify();
         }
      }      
      return kTRUE;   
   } else {
     // Ignore data posting from other than the producer
//      cout<<"Data for container "<<GetName()<<" can be published only by producer task "<<fProducer->GetName()<<endl;
      //AliWarning(Form("Data for container %s can be published only by producer task %s", GetName(), fProducer->GetName()));   
      return kFALSE;           
   }              
}

//______________________________________________________________________________
void AliAnalysisDataContainer::SetFileName(const char *filename)
{
// The filename field can be actually composed by the actual file name followed
// by :dirname (optional):
// filename = file_name[:dirname]
// No slashes (/) allowed
  fFileName = filename;
  fFolderName = "";
  Int_t index = fFileName.Index(":");
  // Fill the folder name
  if (index >= 0) {
    fFolderName = fFileName(index+1, fFileName.Length()-index);
    fFileName.Remove(index);
  }  
  if (!fFileName.Length())
    Fatal("SetFileName", "Empty file name");   
  if (fFileName.Index("/")>=0)
    Fatal("SetFileName", "No slashes (/) allowed in the file name");   
}

//______________________________________________________________________________
void AliAnalysisDataContainer::SetProducer(AliAnalysisTask *prod, Int_t islot)
{
// Set the producer of data. The slot number is required for data type checking.
   if (fProducer) {
     cout<<"Data container "<<GetName()<<" already has a producer: "<<fProducer->GetName()<<endl;
     //AliWarning(Form("Data container %s already has a producer: %s",GetName(),fProducer->GetName()));
   } 
   if (fDataReady) {
     cout<<GetName()<<" container contains data - cannot change producer!"<<endl;
     //AliError(Form("%s container contains data - cannot change producer!", GetName()));
      return;
   }   
   AliAnalysisDataSlot *slot = prod->GetOutputSlot(islot);
   if (!slot) {
     cout<<"Producer task "<<prod->GetName()<<" does not have an output #"<<islot<<endl;
     //AliError(Form("Producer task %s does not have an output #%i", prod->GetName(),islot));
      return;
   }   
   if (!slot->GetType()->InheritsFrom(GetType())) {
     cout<<"Data type "<<slot->GetTitle()<<"for output slot "<<islot<<" of task "<<prod->GetName()<<" does not match container type "<<GetTitle()<<endl;
     //AliError(Form("Data type %s for output slot %i of task %s does not match container type %s", slot->GetType()->GetName(),islot,prod->GetName(),fType->GetName()));
      return;
   }   
   
   fProducer = prod;
   // Add all consumers as daughter tasks
   TIter next(fConsumers);
   AliAnalysisTask *cons;
   while ((cons=(AliAnalysisTask*)next())) {
      if (!prod->GetListOfTasks()->FindObject(cons)) prod->Add(cons);
   }   
}   

//______________________________________________________________________________
AliAnalysisDataWrapper *AliAnalysisDataContainer::ExportData() const
{
// Wraps data for sending it through the net.
   AliAnalysisDataWrapper *pack = 0;
   if (!fData) {
      Error("ExportData", "Container %s - No data to be wrapped !", GetName());
      return pack;
   } 
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (mgr && mgr->GetDebugLevel() > 1) printf("   ExportData: Wrapping data %s for container %s\n", fData->GetName(),GetName());
   pack = new AliAnalysisDataWrapper(fData);
   pack->SetName(fName.Data());
   return pack;
}

//______________________________________________________________________________
void AliAnalysisDataContainer::ImportData(AliAnalysisDataWrapper *pack)
{
// Unwraps data from a data wrapper.
   if (pack) {
      fData = pack->Data();
      if (!fData) {
         Error("ImportData", "No data was wrapped for container %s", GetName());
         return;
      }   
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      if (mgr && mgr->GetDebugLevel() > 1) printf("   ImportData: Unwrapping data %s for container %s\n", fData->GetName(),GetName());
      fDataReady = kTRUE;
      // Imported wrappers do not own data anymore (AG 13-11-07)
      pack->SetDeleteData(kFALSE);
   }   
}      
      
ClassImp (AliAnalysisDataWrapper)

//______________________________________________________________________________
AliAnalysisDataWrapper::AliAnalysisDataWrapper(TObject *data)
                       :TNamed(),
                        fData(data)
{
// Ctor.
   if (data) SetName(data->GetName());
}

//______________________________________________________________________________
AliAnalysisDataWrapper::~AliAnalysisDataWrapper()
{
// Dtor.
   if (fData && TObject::TestBit(kDeleteData)) delete fData;
}   

//______________________________________________________________________________
AliAnalysisDataWrapper &AliAnalysisDataWrapper::operator=(const AliAnalysisDataWrapper &other)
{
// Assignment.
   if (&other != this) {
      TNamed::operator=(other);
      fData = other.fData;
   }   
   return *this;
}

//______________________________________________________________________________
Long64_t AliAnalysisDataWrapper::Merge(TCollection *list)
{
// Merge a list of containers with this one. Containers in the list must have
// data of the same type.
   if (TH1::AddDirectoryStatus()) TH1::AddDirectory(kFALSE);
   if (!fData) return 0;
   if (!list || list->IsEmpty()) return 1;

   SetDeleteData();

   TMethodCall callEnv;
   if (fData->IsA())
      callEnv.InitWithPrototype(fData->IsA(), "Merge", "TCollection*");
   if (!callEnv.IsValid()) {
      cout << "No merge interface for data stored by " << GetName() << ". Merging not possible !" << endl;
      return 1;
   }

   TIter next1(list);
   AliAnalysisDataWrapper *cont;
   // Make a list where to temporary store the data to be merged.
   TList *collectionData = new TList();
   Int_t count = 0; // object counter
   // printf("Wrapper %s 0x%lx (data=%s) merged with:\n", GetName(), (ULong_t)this, fData->ClassName());
   while ((cont=(AliAnalysisDataWrapper*)next1())) {
      cont->SetDeleteData();
      TObject *data = cont->Data();
      if (!data) continue;
      // printf("   - %s 0x%lx (data=%s)\n", cont->GetName(), (ULong_t)cont, data->ClassName());
      collectionData->Add(data);
      count++;
   }
   callEnv.SetParam((Long_t) collectionData);
   callEnv.Execute(fData);
   delete collectionData;

   return count+1;
}

ClassImp(AliAnalysisFileDescriptor)

//______________________________________________________________________________
AliAnalysisFileDescriptor::AliAnalysisFileDescriptor()
                          :TObject(), fLfn(), fGUID(), fUrl(), fPfn(), fSE(),
                           fIsArchive(kFALSE), fImage(0), fNreplicas(0), 
                           fStartBytes(0), fReadBytes(0), fSize(0), fOpenedAt(0), 
                           fOpenTime(0.), fProcessingTime(0.), fThroughput(0.), fTimer()
{
// I/O constructor
}

//______________________________________________________________________________
AliAnalysisFileDescriptor::AliAnalysisFileDescriptor(const TFile *file)
                          :TObject(), fLfn(), fGUID(), fUrl(), fPfn(), fSE(),
                           fIsArchive(kFALSE), fImage(0), fNreplicas(0), 
                           fStartBytes(0), fReadBytes(0), fSize(0), fOpenedAt(0), 
                           fOpenTime(0.), fProcessingTime(0.), fThroughput(0.), fTimer()
{
// Normal constructor
   if (file->InheritsFrom("TAlienFile")) {
      fLfn =(const char*)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetLfn();", file));
      fGUID =(const char*)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetGUID();", file));
      fUrl =(const char*)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetUrl();", file));
      fPfn =(const char*)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetPfn();", file));
      fSE = (const char*)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetSE();", file));
      fImage = (Int_t)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetImage();", file));
      fNreplicas = (Int_t)gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetNreplicas();", file));
      fOpenedAt = gROOT->ProcessLine(Form("((TAlienFile*)%p)->GetOpenTime();", file));
      gROOT->ProcessLine(Form("((AliAnalysisFileDescriptor*)%p)->SetOpenTime(((TAlienFile*)%p)->GetElapsed());", this, file));
   } else {
      fLfn = file->GetName();
      fPfn = file->GetName();
      fUrl = file->GetName();
      fSE = "local";
      if (!fPfn.BeginsWith("/")) fPfn.Prepend(Form("%s/",gSystem->WorkingDirectory()));
      fOpenedAt = time(0);
   }
   fStartBytes = TFile::GetFileBytesRead();
   fIsArchive = file->IsArchive();
   fSize = file->GetSize();
}

//______________________________________________________________________________
AliAnalysisFileDescriptor::AliAnalysisFileDescriptor(const AliAnalysisFileDescriptor &other)
                          :TObject(other), fLfn(other.fLfn), fGUID(other.fGUID),
                           fUrl(other.fUrl), fPfn(other.fPfn), fSE(other.fSE),
                           fIsArchive(other.fIsArchive), fImage(other.fImage),
                           fNreplicas(other.fNreplicas), fStartBytes(other.fStartBytes), fReadBytes(other.fReadBytes),
                           fSize(other.fSize), fOpenedAt(other.fOpenedAt), fOpenTime(other.fOpenTime),
                           fProcessingTime(other.fProcessingTime), fThroughput(other.fThroughput), fTimer()
{
// CC
}

//______________________________________________________________________________
AliAnalysisFileDescriptor &AliAnalysisFileDescriptor::operator=(const AliAnalysisFileDescriptor &other)
{
// Assignment.
   if (&other == this) return *this;
   TObject::operator=(other);
   fLfn       = other.fLfn;
   fGUID      = other.fGUID;
   fUrl       = other.fUrl; 
   fPfn       = other.fPfn;
   fSE        = other.fSE;
   fIsArchive = other.fIsArchive;
   fImage     = other.fImage;
   fNreplicas = other.fNreplicas;
   fStartBytes = other.fStartBytes;;
   fReadBytes = other.fReadBytes;
   fSize      = other.fSize;
   fOpenedAt  = other.fOpenedAt;
   fOpenTime  = other.fOpenTime;
   fProcessingTime = other.fProcessingTime;
   fThroughput = other.fThroughput;
   return *this;
}

//______________________________________________________________________________
AliAnalysisFileDescriptor::~AliAnalysisFileDescriptor()
{
// Destructor
}

//______________________________________________________________________________
void AliAnalysisFileDescriptor::Done()
{
// Must be called at the end of processing, providing file->GetBytesRead() as argument.
   fTimer.Stop();
   const Double_t megabyte = 1048576.;
//   Long64_t stampnow = time(0);
   fReadBytes = TFile::GetFileBytesRead()-fStartBytes;
//   fProcessingTime = stampnow-fOpenedAt;
   fProcessingTime = fTimer.RealTime();
   Double_t readsize = fReadBytes/megabyte;
   fThroughput = readsize/fProcessingTime;
}   

//______________________________________________________________________________
void AliAnalysisFileDescriptor::Print(Option_t*) const
{
// Print info about the file descriptor
   const Double_t megabyte = 1048576.;
   printf("===== Logical file name: %s =====\n", fLfn.Data());
   printf("      Pfn: %s\n", fPfn.Data());
   printf("      url: %s\n", fUrl.Data());
   printf("      access time: %lld from SE: %s  image %d/%d\n", fOpenedAt, fSE.Data(), fImage, fNreplicas);
   printf("      open time: %g [sec]\n", fOpenTime);
   printf("      file size: %g [MB],  read size: %g [MB]\n", fSize/megabyte, fReadBytes/megabyte);
   printf("      processing time [sec]: %g\n", fProcessingTime);
   printf("      average throughput: %g [MB/sec]\n", fThroughput);
}

//______________________________________________________________________________
void AliAnalysisFileDescriptor::SavePrimitive(std::ostream &out, Option_t *)
{
// Stream info to file
   const Double_t megabyte = 1048576.;
   out << "#################################################################" << endl;
   out << "pfn          " << fPfn.Data() << endl;
   out << "url          " << fUrl.Data() << endl;
   out << "se           " << fSE.Data() << endl;
   out << "image        " << fImage << endl;
   out << "nreplicas    " << fNreplicas << endl;
   out << "openstamp    " << fOpenedAt << endl;
   std::ios_base::fmtflags original_flags = out.flags();
   out << setiosflags(std::ios::fixed) << std::setprecision(3);
   out << "opentime     " << fOpenTime << endl;
   out << "runtime      " << fProcessingTime << endl;
   out << "filesize     " << fSize/megabyte << endl;
   out << "readsize     " << fReadBytes/megabyte << endl;
   out << "throughput   " << fThroughput << endl;
   out.flags(original_flags);
}

