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

class TClass;
class TFile;
class TObjArray;
class AliAnalysisTask;
class AliESD;

class AliAnalysisDataContainer : public TNamed {

public:
enum ENotifyMessage {
   kDeleteData,
   kSaveData,
   kFileChange
};   
   AliAnalysisDataContainer();
   AliAnalysisDataContainer(const AliAnalysisDataContainer &cont);
   AliAnalysisDataContainer(const char *name, TClass *type);
   virtual ~AliAnalysisDataContainer();

   // Assignment
   AliAnalysisDataContainer &operator=(const AliAnalysisDataContainer &cont);
   // Getters
   TObject                  *GetData() const      {return fData;}
   TClass                   *GetType() const      {return fType;}
   AliAnalysisTask          *GetProducer() const  {return fProducer;}
   TObjArray                *GetConsumers() const {return fConsumers;}
   virtual void              GetEntry(Long64_t ientry);
   // Setters
   void                      ResetDataReady()     {fDataReady = kFALSE;}
   virtual Bool_t            SetData(TObject *data, Option_t *option="");
   void                      SetDataOwned(Bool_t flag) {fOwnedData = flag;}
   void                      OpenFile(const char *name, Option_t *option="RECREATE");
   void                      SetProducer(AliAnalysisTask *prod, Int_t islot);
   void                      AddConsumer(AliAnalysisTask *cons, Int_t islot);
   void                      DeleteData();
   // Container status checking
   Bool_t                    IsDataReady() const  {return fDataReady;}
   Bool_t                    IsOwnedData() const  {return fOwnedData;}
   Bool_t                    ClientsExecuted() const;
   Bool_t                    HasConsumers() const {return (fConsumers != 0);}
   Bool_t                    HasProducer() const  {return (fProducer != 0);}
   // Send a notify signal to the container
   virtual void              NotifyChange(ENotifyMessage /*type*/) {;}
   // Print connected tasks/status
   void                      PrintContainer(Option_t *option="all", Int_t indent=0) const;
   void                      WriteData();
   
protected:
   Bool_t                    fDataReady;  // Flag that data is ready
   Bool_t                    fOwnedData;  // Flag data ownership
   TFile                    *fFile;       // File storing the data
   TObject                  *fData;       // Contained data
   TClass                   *fType;       // Type of contained data
   AliAnalysisTask          *fProducer;   // Analysis task to which the slot belongs
   TObjArray                *fConsumers;  // List of consumers of the data
   
   ClassDef(AliAnalysisDataContainer,1)  // Class describing a data container for analysis
};
#endif
