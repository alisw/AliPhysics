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

//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliReaderESDTree
//
// Reader for MUON ESD Tree (only for rec)
//
// finck@subatech.in2p3.fr
//
/////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TTree.h>
#include <TFile.h>

#include "AliAOD.h"
#include "AliESD.h"
#include "AliLog.h"
#include "AliReaderESDTree.h"
#include "AliRun.h"
#include "AliRunLoader.h"

ClassImp(AliReaderESDTree)

AliReaderESDTree::AliReaderESDTree(const Char_t* esdfilename, const Char_t* galfilename):
  AliReaderESD(esdfilename,galfilename),
  fTree(0x0)
{
//ctor
}

/********************************************************************/
AliReaderESDTree::~AliReaderESDTree()
{
//dtor 
 delete fTree;
}

/**********************************************************/
Int_t AliReaderESDTree::ReadNext()
{
//reads next event from fFile
//fRunLoader is for reading Kine
  
  AliDebug(1,"Entered");
    
  if (fEventSim == 0x0)  fEventSim = new AliAOD();
  if (fEventRec == 0x0)  fEventRec = new AliAOD();
  
  fEventSim->Reset();
  fEventRec->Reset();
        
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
    {
      if (fFile == 0x0)
        {
          fFile = OpenFile(fCurrentDir);//rl is opened here
          if (fFile == 0x0)
            {
              Error("ReadNext","Cannot get fFile for dir no. %d",fCurrentDir);
              fCurrentDir++;
              continue;
            }
          fCurrentEvent = 0;
        }

      static AliESD* esd = 0x0;
      fTree->SetBranchAddress("ESD", &esd);
      Int_t status = fTree->GetEvent(fCurrentEvent);

      if (!status)
        {
          AliDebug(2,Form("Cannot find event# %d in Tree", fCurrentEvent));
          fCurrentDir++;
          delete fTree;
          fTree = 0x0;
          delete fFile;//we have to assume there is no more ESD objects in the fFile
          fFile = 0x0;
          delete fRunLoader;
          fRunLoader = 0x0;
          continue;
        }

      ReadESD(esd);
      
      fCurrentEvent++;
      fNEventsRead++;
      return 0;//success -> read one event
    }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array  
   
  return 1; //no more directories to read
}

/**********************************************************/
TFile* AliReaderESDTree::OpenFile(Int_t n)
{
//opens fFile with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename = dirname +"/"+ fESDFileName;
 TFile *ret = TFile::Open(filename.Data()); 

 if (ret == 0x0)
  {
    Error("OpenFiles","Can't open fFile %s",filename.Data());
    return 0x0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open fFile  %s",filename.Data());
    return 0x0;
  }
 
 TString esdname = "esdTree";
 fTree = dynamic_cast<TTree*> (ret->Get(esdname));

 if (!fTree)
  {
    Error("OpenFiles","Can't open ESD Tree %s",esdname.Data());
    delete ret;
    return 0x0;

  }
 
 if (fReadSim )
  {
   fRunLoader = AliRunLoader::Open(dirname +"/"+ fGAlFileName);
   if (fRunLoader == 0x0)
    {
      Error("OpenFiles","Can't get RunLoader for directory %s",dirname.Data());
      delete fTree;
      delete ret;
      return 0x0;
    }
    
   fRunLoader->LoadHeader();
   if (fRunLoader->LoadKinematics())
    {
      Error("Next","Error occured while loading kinematics.");
      delete fRunLoader;
      delete fTree;
      delete ret;
      return 0x0;
    }
  }
   
 return ret;
}
