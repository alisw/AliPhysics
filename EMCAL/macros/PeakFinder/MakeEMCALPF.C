///
/// \file MakeEMCALPF.C
/// \ingroup EMCAL_SimRecDB
/// \brief Set OCDB for EMCAL Peak Finder
///
/// Create OCDB for Peak Finder vectods
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
///

#if !defined(__CINT__)
#include <TString.h>
#include <TStyle.h>
#include <TFile.h>

#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

#include "AliCaloPeakFinderVectors.h"
#endif

///
/// Main method
///
void MakeEMCALPF()
{
  const char* macroname = "MakeEMCALPF.C";
  
  TFile *f2 = new TFile("peakfindervectors2.root",  "read" ); 
  
  //AliCaloPeakFinderVectors *pfv =  (AliCaloPeakFinderVectors* )f2->GetKey( "AliCaloPeakFinderVectors"); 
  
  AliCaloPeakFinderVectors pfv =  *((AliCaloPeakFinderVectors* )f2->GetKey( "AliCaloPeakFinderVectors")); 
  
  f2->Close();
  
  {
    //TString Storage = "local://home/perthi/aliroot-current/OCDB/";
    TString storageName = "local://OCDB/";  
    
    if(!storageName.BeginsWith("local://") && !storageName.BeginsWith("alien://")) 
    {
      //Error(macroname ,"STORAGE variable set to %s is not valid. Exiting\n",storageName.Data());
      return;
    }
    
    //Info(macroname,"Saving PF objects in CDB storage %s",storageName.Data());
    
    AliCDBManager* cdb = AliCDBManager::Instance();
    
    AliCDBStorage* storage = cdb->GetStorage(storageName.Data());
    if(!storage)
    {
      //Error(macroname,"Unable to open storage %s\n",storageName.Data());
      return;
    }
    
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Per Thomas Hille");
    md->SetComment("Peak-Finder vectors for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Config/PeakFinder",0,AliCDBRunRange::Infinity());
    
    // if(pfv == 0)
    //   {
    // 	cout << "  ERROR !!!!!!!!!" << endl;
    //   }
    
    // else
    
    {
      storage->Put( &pfv,id,md);
    }
    
    //  delete md;
  }
  
  
  //
}

