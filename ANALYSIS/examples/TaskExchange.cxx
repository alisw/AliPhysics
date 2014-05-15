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

//______________________________________________________________________________
// Producer task that exchanges data with other task using exchange containers
//______________________________________________________________________________
// 
#include "TaskExchange.h"

#include "AliAnalysisManager.h"
#include "TList.h"
#include "TObjArray.h"

ClassImp(TaskProducer)

//________________________________________________________________________
TaskProducer::TaskProducer() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutput(0),
    fExchangedData(0)
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
TaskProducer::TaskProducer(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutput(0),
    fExchangedData(0)
{
   // Constructor
   // Define input and output slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container

   // Output slot #1 used to publish the task private data
   DefineOutput(1, TList::Class());
    
   // Output slot #2 used to publish the exchanged data
   DefineOutput(2, TObjArray::Class()); // Type can be anything else
}

//________________________________________________________________________
TaskProducer::~TaskProducer()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fExchangedData;
}

//________________________________________________________________________
void TaskProducer::UserCreateOutputObjects()
{
   // Create histograms
   // Called once (on the worker node)
        
   fOutput = new TList();
   fOutput->SetOwner();  // IMPORTANT!
    
//    fOutput->Add(fHistPt);
//    fOutput->Add(fHistEta);
   // NEW HISTO added to fOutput here

   // Data to be exchanged can be initialized here, but updated in UserExec
   fExchangedData = new TObjArray();
   fExchangedData->SetOwner();
 
   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
// Do NOT post data for the exchange containers !
//   PostData(2, fExchangedData);
}

//________________________________________________________________________
void TaskProducer::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event

    //Here we fill the private output
    // ...
    PostData(1, fOutput);
    
    // And publish the exchanged data
    fExchangedData->Delete();
    TNamed *data = new TNamed(Form("<data by %s, ev. %d>", GetName(), fEntry), "");
    fExchangedData->Add(data);
    PostData(2, fExchangedData);
    Printf("    -- %s published: %s", GetName(), data->GetName());
}

//______________________________________________________________________________
// Consumer task reading data produced by other task
//______________________________________________________________________________
// 
ClassImp(TaskConsumer)

//________________________________________________________________________
TaskConsumer::TaskConsumer() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutput(0),
    fImported1(0),
    fImported2(0)
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
TaskConsumer::TaskConsumer(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutput(0),
    fImported1(0),
    fImported2(0)
{
   // Constructor
   // Define input and output slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container
   
   // Input slot #1 reads from an exchange container
   DefineInput(1, TObjArray::Class());
   // Input slot #2 reads from an exchange container
   DefineInput(2, TObjArray::Class());

   // Output slot #0 is reserved
   // Output slot #1 used to publish the task private data
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
TaskConsumer::~TaskConsumer()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   // Do not delete the imported data - we are not owners
}

//________________________________________________________________________
void TaskConsumer::UserCreateOutputObjects()
{
   // Create histograms
   // Called once (on the worker node)
        
   fOutput = new TList();
   fOutput->SetOwner();  // IMPORTANT!
    
//    fOutput->Add(fHistPt);
//    fOutput->Add(fHistEta);
   // NEW HISTO added to fOutput here

   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void TaskConsumer::UserExec(Option_t *) 
{
   // Main loop
   // Called for each event
   // This is how we get the actual exchange data
   TObjArray *exchange1 = (TObjArray*)GetInputData(1);
   if (!exchange1) {
      Error("UserExec", "Task %s could not read the exchanged data for slot 1", GetName());
      return;
   }
   fImported1 = (TNamed*)exchange1->At(0);
   TObjArray *exchange2 = (TObjArray*)GetInputData(1);
   if (!exchange2) {
      Error("UserExec", "Task %s could not read the exchanged data for slot 2", GetName());
      return;
   }
   fImported2 = (TNamed*)exchange2->At(0);
   Printf("    -- imported: %s and %s", fImported1->GetName(), fImported2->GetName());

   //Here we fill the private output
   // ...
   PostData(1, fOutput);    
}
