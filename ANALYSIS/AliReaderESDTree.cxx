#include "AliReaderESDTree.h"
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


#include <AliRun.h>
#include <AliRunLoader.h>

#include <AliESD.h>
#include "AliAOD.h"

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
  
  if (AliVAODParticle::GetDebug())
    Info("ReadNext","Entered");
    
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
	  if (AliVAODParticle::GetDebug() > 2 )
	    {
	      Info("ReadNext","Can not find event# %d in Tree", fCurrentEvent);
	    }
	  fCurrentDir++;
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
    return 0x0;

  }
 
 if (fReadSim )
  {
   fRunLoader = AliRunLoader::Open(dirname +"/"+ fGAlFileName);
   if (fRunLoader == 0x0)
    {
      Error("OpenFiles","Can't get RunLoader for directory %s",dirname.Data());
      delete ret;
      return 0x0;
    }
    
   fRunLoader->LoadHeader();
   if (fRunLoader->LoadKinematics())
    {
      Error("Next","Error occured while loading kinematics.");
      delete fRunLoader;
      delete ret;
      return 0x0;
    }
  }
   
 return ret;
}
