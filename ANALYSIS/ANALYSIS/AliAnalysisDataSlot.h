#ifndef ALIANALYSISDATASLOT_H
#define ALIANALYSISDATASLOT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliAnalysisDataSlot
/// \brief AliAnalysysDataSlot
/// Class representing a data slot of an analysis task.
/// An analysis slot enforces a certain data type required by the Exec()
/// method of the analysis task. The slot must be connected to a data 
/// container with data of the same type.
///
/// The class should not be directly created by users - it is created by
/// each AliAnalysisTask when defining its input/output slots using:
///
///    AliAnalysisTask::SetInput(Int_t index, TClass *type);
///    AliAnalysisTask::SetOutput(TClass *type);
///
/// An existing data contaner (AliAnalysisDataContainer) can be connected to the
/// input/output slots of an analysis task using:
///
///   AliAnalysisModule::ConnectInput(AliAnalysisTask *task, Int_t islot,
///                                   AliAnalysisDataContainer *cont)
///   AliAnalysisModule::ConnectOutput(AliAnalysisTask *task,
///                                    AliAnalysisDataContainer *cont)
/// To connect a slot to a data container, the data types declared by both must
/// match.
/// \author Andrei Gheata
/// \date 31/05/2006



#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TClass;
class TTree;
class AliAnalysisDataContainer;
class AliAnalysisTask;

class AliAnalysisDataSlot : public TNamed {

public:
   AliAnalysisDataSlot() : TNamed(), fType(NULL), fParent(NULL), fContainer(NULL) {}
   AliAnalysisDataSlot(TClass *type, AliAnalysisTask *task);
   AliAnalysisDataSlot(const AliAnalysisDataSlot &slot);
   virtual ~AliAnalysisDataSlot() {}

   // Assignment
   AliAnalysisDataSlot &operator=(const AliAnalysisDataSlot &slot);
   // Connect some container to the slot
   Bool_t                    ConnectContainer(AliAnalysisDataContainer *cont);
   static Int_t              EnableBranch(const char *bname, TTree *tree);
   // Getters
   void                     *GetBranchAddress(const char *branch) const;
   Bool_t                    SetBranchAddress(const char *branch, void *address);
   TClass                   *GetType() const;
   AliAnalysisTask          *GetParent() const    {return fParent;}
   AliAnalysisDataContainer *GetContainer() const {return fContainer;}
   TObject                  *GetData() const;
   // Slot status checking
   Bool_t                    IsConnected() const  {return ((fContainer)?kTRUE:kFALSE);}
   Bool_t                    IsDataReady() const;

private:
   void                      SetType(TClass *type) {fType = type;}
   
protected:
   TClass                   *fType;       //!<! Type of the slot
   AliAnalysisTask          *fParent;     ///< Analysis task to which the slot belongs
   AliAnalysisDataContainer *fContainer;  ///< Container connected to the slot
   
   ClassDef(AliAnalysisDataSlot,1)  // Class describing an analysis data slot
};
#endif
