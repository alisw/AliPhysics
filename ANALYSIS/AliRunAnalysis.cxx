#include "AliRunAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TString.h>
#include <TObjString.h>
#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TObjArray.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliESDtrack.h>
#include <AliESD.h>

#include "AliEventCut.h"


ClassImp(AliRunAnalysis)
AliRunAnalysis::AliRunAnalysis():
 TTask("RunAnalysis","Alice Analysis Manager")	,
 fDirs(),
 fEventCut(0x0),
 fFileName("AliESDs.root"),
 fReadKinematics(kFALSE)
{
  //ctor
}
/*********************************************************/

AliRunAnalysis::~AliRunAnalysis()
{
  //dtor
  delete fDirs;
  delete fAnalysies;
  delete fEventCut;
}
/*********************************************************/

Int_t AliRunAnalysis::Run()
{
 //makes analysis
 
 Int_t currentdir = 0;
 Int_t ndirs;
 if (fDirs) //if array with directories is supplied by user
  {
    ndirs = fDirs->GetEntries(); //get the number if directories
  }
 else
  {
    ndirs = 0; //if the array is not supplied read only from current directory
  }

 /******************************/ 
 /*  Init Event                */ 
 /******************************/ 
 if (fAnalysies == 0x0)
  {
    Info("Run","No analysis present");
    return 0;
  }
  
 for (Int_t an = 0; an < fAnalysies->GetEntries(); an++)
  {
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies->At(an);
      analysis->Init();
  }

 do
  {
    TFile* file = OpenFile(currentdir);
    if (file == 0x0)
      {
        Error("Run","Cannot get File for dir no. %d",currentdir);
        currentdir++;
        continue;
      } 
    AliStack* stack = 0x0;
    AliRunLoader* rl = 0x0;
    if (fReadKinematics)
     {
       const TString& dirname = GetDirName(currentdir);
       if (dirname == "")
        {
         Error("Run","Can not get directory name");
         return 0x0;
        }
       TString filename = dirname +"/galice.root";

       rl = AliRunLoader::Open(filename);
       if ( rl == 0x0 )
        {
          Error("Run","Can't get Run Loader from dir %s",filename.Data());
          delete file;
          currentdir++;
          continue;
        }
       if( rl->LoadHeader() )
        {
          Error("Run","Error while loading Header from dir %s",filename.Data());
          delete file;
          delete rl;
          currentdir++;
          continue;
        }
       if( rl->LoadKinematics() )
        {
          Error("Run","Error while loading Kinematics from dir %s",filename.Data());
          delete file;
          delete rl;
          currentdir++;
          continue;
        }
     }
   
    file->cd();
    TIter keyiter(file->GetListOfKeys());
    TKey* key;
    while (( key = (TKey*)keyiter.Next() ))
     {
      if (key == 0x0)
       {
        if (GetDebug() > 2 )
          {
            Info("Run","No more keys.");
          }
        break;
       }

      TObject* esdobj = key->ReadObj();
      if (esdobj == 0x0)
       {
         if (GetDebug() > 2 )
           {
             Info("ReadNext","Key read NULL. Key Name is %s",key->GetName());
             key->Dump();
           }
         continue;
       }
      if (GetDebug() > 9 ) esdobj->Dump();
      AliESD* esd = dynamic_cast<AliESD*>(esdobj);

      if (esd == 0x0)
       {
         if (GetDebug() > 7 )
           {
             Info("ReadNext","It is not an ESD object");
           }
         delete esdobj;
         continue;
       }

      if (fReadKinematics)
       {
        TString esdname(esd->GetName());
        esdname.ReplaceAll("ESD","");
        Int_t nev = atoi(esdname);
        Info("Run","ESD name is %s, Event Number is %d",esd->GetName(),nev);
        if (rl->GetEvent(nev))
         {
           Error("Run","Error occured while RunLoader GetEvent %d",nev);
           delete esd;
           continue;
         }
        stack = rl->Stack();
        if (stack == 0x0) 
         {
           Error("Run","Can not get stack");
           delete esd;
           continue;
         }
       }
      /******************************/ 
      /*  Event Cut                 */ 
      /******************************/ 
      if (fEventCut)
       {
         if (fEventCut->Pass(esd))
          {
            if (GetDebug()) Info("Run","Event rejected by Event Cut");
            delete esd;
            continue; //Did not pass the 
          }
       }
      /******************************/ 
      /*  Process Event             */ 
      /******************************/ 
      for (Int_t an = 0; an < fAnalysies->GetEntries(); an++)
       {
           AliAnalysis* analysis = (AliAnalysis*)fAnalysies->At(an);
           analysis->ProcessEvent(esd,stack);
       }
      delete esd;
     }//end of loop over keys in file
     
   delete file;
   delete rl;
   currentdir++;
    
  }while (currentdir < ndirs);//end of loop over directories

 /******************************/ 
 /*  Finish Event              */ 
 /******************************/ 
 for (Int_t an = 0; an < fAnalysies->GetEntries(); an++)
  {
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies->At(an);
      analysis->Init();
  }

 return 0;   
}
/*********************************************************/

TFile* AliRunAnalysis::OpenFile(Int_t n)
{
//opens file with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename = dirname +"/"+ fFileName;
 TFile *ret = TFile::Open(filename.Data()); 

 if ( ret == 0x0)
  {
    Error("OpenFiles","Can't open file %s",filename.Data());
    return 0x0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open file  %s",filename.Data());
    return 0x0;
  }
 
 return ret;
}
/*********************************************************/

TString& AliRunAnalysis::GetDirName(Int_t entry)
{
//returns directory name of next one to read
  TString* retval;//return value
  if (fDirs ==  0x0)
   {
     retval = new TString(".");
     return *retval;
   }

  if ( (entry>fDirs->GetEntries()) || (entry<0))//if out of bounds return empty string
   {                                            //note that entry==0 is accepted even if array is empty (size=0)
     Error("GetDirName","Name out of bounds");
     retval = new TString();
     return *retval;
   }

  if (fDirs->GetEntries() == 0)
   { 
     retval = new TString(".");
     return *retval;
   }

  TClass *objclass = fDirs->At(entry)->IsA();
  TClass *stringclass = TObjString::Class();

  TObjString *dir = (TObjString*)objclass->DynamicCast(stringclass,fDirs->At(entry));

  if(dir == 0x0)
   {
     Error("GetDirName","Object in TObjArray is not a TObjString or its descendant");
     retval = new TString();
     return *retval;
   }
  if (gDebug > 0) Info("GetDirName","Returned ok %s",dir->String().Data());
  return dir->String();
}
/*********************************************************/

void  AliRunAnalysis::Add(AliAnalysis* a)
{
  //adds a to the list of analysis
  if (fAnalysies == 0x0) fAnalysies = new TObjArray();
  fAnalysies->Add(a);
}
