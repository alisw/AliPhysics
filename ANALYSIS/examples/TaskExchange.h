/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//______________________________________________________________________________
//
// Task that produces some generic data to be exchanged with some consumers
//
//______________________________________________________________________________

#ifndef TASKEXCHANGE_H
#define TASKEXCHANGE_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TList;
class TObjArray;

class TaskProducer : public AliAnalysisTaskSE {
 public:
    TaskProducer();
    TaskProducer(const char *name);
    virtual ~TaskProducer();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    
 private:
    TList           *fOutput;        //! Output list (do not stream)
    TObjArray       *fExchangedData; //! Arbitrary data (do not stream)
    
    TaskProducer(const TaskProducer&); // not implemented
    TaskProducer& operator=(const TaskProducer&); // not implemented
    
    ClassDef(TaskProducer, 1); // example of producer
};

//______________________________________________________________________________
//
// Task that uses generic data produced by other tasks
//
//______________________________________________________________________________

class TaskConsumer : public AliAnalysisTaskSE {
 public:
    TaskConsumer();
    TaskConsumer(const char *name);
    virtual ~TaskConsumer();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    
 private:
    TList           *fOutput;        //! Output list (do not stream)
    TNamed          *fImported1;     //! Arbitrary data (do not stream)
    TNamed          *fImported2;     //! Arbitrary data (do not stream)
    
    TaskConsumer(const TaskConsumer&); // not implemented
    TaskConsumer& operator=(const TaskConsumer&); // not implemented
    
    ClassDef(TaskConsumer, 1); // example of consumer task
};

#endif

