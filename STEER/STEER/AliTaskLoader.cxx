
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  class AliTaskLoader                                                                    //
//                                                                                         //
//  Container of all data needed for full                                                  //
//  description of each data type                                                          //
//  (Hits, Kine, ...)                                                                      //
//                                                                                         //
//  Each data loader has a basic standard setup of BaseLoaders                             //
//  which can be identuified by indexes (defined by EStdBasicLoaders)                      //
//  Data managed by these standard base loaders has fixed naming convention                //
//  e.g. - tree with hits is always named TreeH                                            //
//                     (defined in AliLoader::fgkDefaultHitsContainerName)                 //
//       - task DtectorName+Name defined                                                   //
//                                                                                         //
//  EStdBasicLoaders   idx     Object Type        Description                              //
//      kData           0    TTree or TObject     main data itself (hits,digits,...)       //
//      kTask           1        TTask            object producing main data               //
//      kQA             2        TTree                quality assurance tree               //
//      kQATask         3        TTask            task producing QA object                 //
//                                                                                         //
//                                                                                         //
//  User can define and add more basic loaders even Run Time.                              //
//  Caution: in order to save information about added base loader                          //
//  user must rewrite Run Loader to galice.file, overwriting old setup                     //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliTaskLoader.h"
#include "AliDataLoader.h"
#include "AliLog.h"

ClassImp(AliTaskLoader)

//______________________________________________________________________________
AliTaskLoader::AliTaskLoader(const TString& name, AliDataLoader* dl, 
			     TTask* parentaltask, Bool_t storeontop):
 AliBaseLoader(name,dl,storeontop),
 fParentalTask(parentaltask)
{
  //
  // Constructor
  //
}

//______________________________________________________________________________
void AliTaskLoader::Clean()
{
  //
  // Removes tasl from parental task
  // DO NOT DELETE OBJECT contrary to BaseLoader
  //
  AliDebug(1, Form("Clean %s %s",GetName(),GetDataLoader()->GetName()));
  TObject* obj = Get();
  if(obj)
    { 
      AliDebug(1, Form("cleaning %s.",GetName()));
      RemoveFromBoard(obj);
    }
}


//______________________________________________________________________________
void AliTaskLoader::RemoveFromBoard(TObject* obj)
{
  //
  // Removes the task "obj" from the board
  //
  GetParentalTask()->GetListOfTasks()->Remove(obj);
}

//______________________________________________________________________________
Int_t AliTaskLoader::AddToBoard(TObject* obj)
{
  //
  // Adds task "obj" to the board
  //
  TTask* task = dynamic_cast<TTask*>(obj);
  if (task == 0x0)
    {
      AliError("To TTask board can be added only tasks.");
      return 1;
    }
  GetParentalTask()->Add(task);
  return 0;
}

//______________________________________________________________________________
TObject* AliTaskLoader::Get() const
{
  //
  // Returns pointer to the current task
  //
  return (GetParentalTask()) ? GetParentalTask()->GetListOfTasks()->FindObject(GetName()) : 0x0;
}

//______________________________________________________________________________
TTask* AliTaskLoader::GetParentalTask() const
{
  //
  // Returns parental tasks for this task
  //
  return fParentalTask;
}



