// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliDAQ.h"
#endif

/**
 * @file EnableHLTInGRP.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @ingroup alihlt_tutorial
 * @brief Patches GRP to process recorded data with HT in mode A.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q 'EnableHLTInGRP.C(\<runNumber\>,\<grpPath\>)'
 * </pre>
 * Where \<runNumber\> is the number of the run the GRP corresponds to. This
 * can be found from the file name of the GRP entry of interest in the CDB.
 * \<grpPath\> is the path to the local GRP/GRP/Data directory where the GRP
 * entry to be modified is stored. Setting a value for \<grpPath\> is optional
 * and is the current working directory by default.
 *
 * Data recorded in mode A will not have the HT bit set in the GRP's active
 * detector mask. Trying to process such data with the HT framework will fail.
 * To be able to process such runs the GRP must be downloaded and patched with
 * this macro. It will set the HT bit in the GRP and write a new version of
 * the GRP with and incremented version number.
 * 
 * @return true if the macro completed successfully and false otherwise.
 */
bool EnableHLTInGRP(Int_t runNumber,
		    const char* cdbURI = "local://$ALICE_ROOT/OCDB",
		    const char* grpPath = "local://./")
{
	AliCDBManager* cdb = AliCDBManager::Instance();
	cdb->SetDefaultStorage(cdbURI);
	cdb->SetSpecificStorage("GRP/GRP/Data", grpPath);
	cdb->SetRun(runNumber);
	AliCDBEntry* entry = (AliCDBEntry*)cdb->Get("GRP/GRP/Data")->Clone();
	if (entry == NULL)
	{
		cerr << "ERROR: Could not fetch the GRP entry for run " << runNumber << " and path " << grpPath << endl;
		return false;
	}
	AliGRPObject* grp = (AliGRPObject*)entry->GetObject();
	if (grp == NULL)
	{
		cerr << "ERROR: Could not fetch the GRP object for run " << runNumber << " and path " << grpPath << endl;
		return false;
	}
	grp->SetDetectorMask(grp->GetDetectorMask() | (0x1 << AliDAQ::kHLTId));
	entry->SetVersion(entry->GetId().GetVersion() + 1);
	cdb->Put(entry);
	return true;
}

/// Just print the usage if no parameters are given.
void EnableHLTInGRP()
{
	cout << "Usage: EnableHLTInGRP.C'(<runNumber> [, <grpPath>])'" << endl;
	cout << "Where <runNumber> is the number of the run the GRP corresponds to." << endl;
	cout << "This can be found from the file name of the GRP entry of interest in the CDB." << endl;
	cout << "<grpPath> is the path to the local GRP/GRP/Data directory where the GRP entry to be modified is stored." << endl;
	cout << "Setting a value for <grpPath> is optional and is the current working directory by default." << endl;
}
