// Script to convert OCDB files to tender used files
// Author: Jiri Kral

#if !defined(__CINT__)
#include <TString.h>
#include <TH2.h>
#include <TF1.h>

#include "AliRun.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALGeoParams.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALOCDBTenderConverter( Int_t runNum, char *outFileName ){
	Int_t i;
	char buf[100];

	TH2D *histo;

	TFile *outFile;
	
	AliCDBManager *man;
	AliCDBStorage *stor;
	AliCaloCalibPedestal *ped;

	// created the OCDB manager
	man = AliCDBManager::Instance();

	// point it to local storage
	// !!! careful, one must build an exact path of OCDB directories
	// and store the file in those
	// here "./OCDB/EMCAL/Calib/Pedestals/Run*.root) for masks
	stor = man->GetStorage( "local://$ALICE_ROOT/OCDB");
	
	// load the file data
	ped = (AliCaloCalibPedestal*)(stor->Get("EMCAL/Calib/Pedestals", runNum)->GetObject());

	// get the array of histos
	TObjArray map = ped->GetDeadMap();

	outFile = new TFile( outFileName, "RECREATE" );

	// rename the histos and save
	for( i = 0; i < map.GetEntries(); i++ ){
		histo = (TH2D*)(map[i]);
                printf("\n !!! EMCALBadChannelMap_Mod%d",i );
		sprintf( buf, "EMCALBadChannelMap_Mod%d", i );

		histo->SetName( buf );
		histo->SetTitle( buf );

		histo->Write();
	}

	// cleanup
	delete outFile;
	delete ped;
	delete stor;
	delete man;
	
}