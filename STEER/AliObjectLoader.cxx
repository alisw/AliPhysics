/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  class AliObjectLoader                                                                  //
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

#include "AliObjectLoader.h"
#include "AliDataLoader.h"
#include "AliLog.h"
#include <TFolder.h>

ClassImp(AliObjectLoader)

//______________________________________________________________________________
AliObjectLoader::AliObjectLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop):
  AliBaseLoader(name,dl,storeontop)
{
  //
  // Constructor
  //
}

//______________________________________________________________________________
TFolder* AliObjectLoader::GetFolder() const
{
  //
  // Returns pointer to the object folder
  //
  TFolder* df = GetDataLoader()->GetFolder();
  if (df == 0x0)
    {
      AliFatal("Data Folder is NULL");
    }
  return df;
}

//______________________________________________________________________________
void AliObjectLoader::RemoveFromBoard(TObject* obj)
{
  //
  // Removes "obj" from the board
  //
  GetFolder()->Remove(obj);
}

//______________________________________________________________________________
Int_t AliObjectLoader::AddToBoard(TObject* obj)
{
  //
  // Adds "obj" to the board
  //
  GetFolder()->Add(obj);
  return 0;
}

//______________________________________________________________________________
TObject* AliObjectLoader::Get() const
{
  //
  // Returns pointer to the object loader
  //
  return (GetFolder()) ? GetFolder()->FindObject(GetName()) : 0x0;
}



