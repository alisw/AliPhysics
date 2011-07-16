#ifndef ALIBASELOADER_H
#define ALIBASELOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//__________________________________________
////////////////////////////////////////////
//                                        //
//  class AliBaseLoader                   //
//                                        //
//  Loader responsible for one data type  //
//  i.e. Hits, Kine, etc.                 //
//  many objects type can be assciated    //
//  with one data type: storing object    //
//  (usually tree), task producing it,    //
//  Quality Assurance(QA), QA Task, and   //
//  others.                               //
//                                        //
//                                        //
////////////////////////////////////////////

class TObject;
class TFolder;
class AliDataLoader;

#include <TNamed.h>
#include <TDirectory.h>
#include <TString.h>

class AliBaseLoader: public TNamed
{
  public:
    AliBaseLoader();
    AliBaseLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop = kFALSE);
    
    virtual ~AliBaseLoader(){};
     
    virtual Int_t      Load(Option_t* opt="");
    virtual void       Unload();
    virtual Int_t      Reload();
    virtual Int_t      WriteData(Option_t* opt="");
    virtual void       Clean();
    virtual Int_t      Post();//Takes from file and sends to proper TFolder (Data Folder)
    virtual Int_t      Post(TObject* data);//Sends to proper TFolder (Data Folder)
    virtual TObject*   Get() const = 0; 
    Bool_t             IsLoaded()const{return fIsLoaded;}
    void               SetDataLoader(AliDataLoader* dl){fDataLoader = dl;}
    void               SetEventFolder(TFolder* /*ef*/){;}
    void               SetDoNotReload(Bool_t flag){fDoNotReload = flag;}
    Bool_t             DoNotReload() const {return fDoNotReload;}
    TDirectory*        GetDirectory() const;//returns pointer to directory where data are stored. 
    TObject*           GetFromDirectory(const char *name) const
      {return (GetDirectory())?GetDirectory()->Get(name):0x0;}    
   protected:
    
    virtual Int_t      AddToBoard(TObject* obj) = 0;//add to white board - board can be TTask or TFolder
    virtual void       RemoveFromBoard(TObject* obj) = 0;
    
    AliDataLoader*     GetDataLoader() const;

    Bool_t             fIsLoaded;    //!  flag indicating if data are loaded
    Bool_t             fStoreInTopOfFile;// if true, data are stored in top of file ->Indicates fDoNotReload == kTRUE

   private:

    AliBaseLoader(const AliBaseLoader&);   //Not implemented
    AliBaseLoader& operator=(const AliBaseLoader&); //Not implemented

    Bool_t             fDoNotReload; // if this flag is on object is not reloaded while GetEvent is called.
                                     //Specially important for tasks. Task loops over events while producing data, 
	                 //and has a base loader which writes it to file every processed event.
	                 //If this flag is not on, while taking next event, loader deletes task
	                 // and tries to get new one from file
    AliDataLoader*     fDataLoader;  //! pointer to Data Loader this Base Loader belongs to

 ClassDef(AliBaseLoader,1)    
};


#endif


