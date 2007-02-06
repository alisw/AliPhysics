#ifndef ALIANALYSISDATASLOT_H
#define ALIANALYSISDATASLOT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysDataSlot - Class representing a data slot of an analysis task.
//      An analysis slot enforces a certain data type required by the Exec()
//      method of the analysis task. The slot must be connected to a data 
//      container with data of the same type.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TClass;
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
   TClass                   *fType;       //! Type of the slot
   AliAnalysisTask          *fParent;     // Analysis task to which the slot belongs
   AliAnalysisDataContainer *fContainer;  // Container connected to the slot
   
   ClassDef(AliAnalysisDataSlot,1)  // Class describing an analysis data slot
};
#endif
