////////////////////////////////////////////////////////////////////////////////
///
/// \file AliEMCALOCDBTenderConverter.cxx
/// \brief Script to convert OCDB files to tender used files
///
///  Script to convert OCDB files to tender used files
///
/// \author Jiri Kral, Jiri.Kral@cern.ch (Jyvaskyla)
///
////////////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__)
#include <TH2D.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliCaloCalibPedestal.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

///
/// Main method:
/// 
/// First it opens the OCDB file with bad map
/// Then it gets the bad map histograms and put it in a file.
///
/// \param runNum: reference run number to extract the parameters
/// \param outFileName: path and file name of output file with histograms
///
void AliEMCALOCDBTenderConverter( Int_t runNum, char *outFileName )
{
  Int_t i;
  char buf[100];
      
  // Create the OCDB manager
  AliCDBManager * man = AliCDBManager::Instance();
  
  // Point it to local storage
  // !!! careful, one must build an exact path of OCDB directories
  // and store the file in those
  // here "./OCDB/EMCAL/Calib/Pedestals/Run*.root) for masks
  AliCDBStorage * stor = man->GetStorage( "local://$ALICE_ROOT/OCDB");
  
  // Load the file data
  AliCaloCalibPedestal * ped = (AliCaloCalibPedestal*)(stor->Get("EMCAL/Calib/Pedestals", runNum)->GetObject());
  
  // Get the array of histos
  TObjArray map = ped->GetDeadMap();
  
  TFile * outFile = new TFile( outFileName, "RECREATE" );
  
  // Rename the histos and save
  for( i = 0; i < map.GetEntries(); i++ )
  {
    TH2D * histo = (TH2D*)(map[i]);
    printf("\n !!! EMCALBadChannelMap_Mod%d",i );
    sprintf( buf, "EMCALBadChannelMap_Mod%d", i );
    
    histo->SetName( buf );
    histo->SetTitle( buf );
    
    histo->Write();
  }
  
  // Cleanup
  delete outFile;
  delete ped;
  delete stor;
  delete man;
}