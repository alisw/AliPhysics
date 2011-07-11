/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: DisplaydHLTData.C 37070 2009-11-20 13:53:08Z aszostak $

/**
 * \ingroup macros
 * \file CreateCDBFieldIntegrals.C
 * \brief Macro for creating and storing the magnetic field integrals in the CDB.
 *
 * This macro calculates the magnetic field integrals for the dipole magnet and
 * writes them into the CDB under the HLT/ConfigMUON/FieldIntegrals path.
 * These are then used by the online dHLT reconstruction components during
 * momentum estimation of the tracks found.
 * The simplest method to run this macro is with the following command:
 * \code
 *   > aliroot -b -q -l CreateCDBFieldIntegrals.C
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliHLTMUONConstants.h"
#include "AliMagF.h"
#include "TMap.h"
#include "TArrayD.h"
#include "TObjString.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

/**
 * Calculates the average magnetic field integral for the dipole magnet given the
 * L3 solenoid current and dipole magnet current.
 * \param [in] currentL3  The current in the L3 solenoid magnet.
 * \param [in] currentDip  The current in the dipole magnet.
 * \param [out] result  The magnetic field integral result.
 * \param [in] sqrts  The centre of mass energy for the beams.
 * \param [in] beamtype  The beam type as given the AliMagF, eg. p-p, A-A or none.
 * \note The sign of the current indicates the polarity setting for the magnet.
 * \returns  true if the magnetic field was calculated successfully and false otherwise.
 */
bool CalculateIntegral(
		Float_t currentL3, Float_t currentDip, Float_t& result,
		Float_t sqrts = 0, const char* beamtype = "none"
	)
{
	AliMagF* fld = AliMagF::CreateFieldMap(currentL3, currentDip, AliMagF::kConvLHC, kFALSE, sqrts, beamtype);
	if (fld == NULL) return false;
	
	TArrayD totals;
	for (int theta = 171; theta <= 179; ++theta)
	for (int phi = 0; phi < 360; phi += 10)
	{
		TVector3 dr;
		dr.SetMagThetaPhi(1, theta*TMath::Pi()/180., phi*TMath::Pi()/180.);
		Double_t dx = 1e-3 * dr.Mag();
		
		TVector3 r(0, 0, 0);
		int nsteps = 0;
		Double_t total = 0;
		while (r.Z() > -1700 && nsteps < 10000)
		{
			Double_t p[3] = {r.X(), r.Y(), r.Z()};
			Double_t b[3];
			fld->Field(p, b);
			total += b[0] * dx;
			++nsteps;
			r += dr;
		}
		
		totals.Set(totals.GetSize()+1);
		totals[totals.GetSize()-1] = total;
	}
	
	result = - TMath::Mean(totals.GetSize(), totals.GetArray());
	Float_t rms = TMath::RMS(totals.GetSize(), totals.GetArray());
	
	cout << "Integrated field value = " << result
		<< " +/- " << rms
		<< ", for L3 current = " << currentL3
		<< ", dipole current = " << currentDip << endl;
	return true;
}

/**
 * Calculates the magnetic field integrals for the dipole magnet and stores them
 * into the CDB.
 * \param cdbPath  The path to the local storage.
 * \param version  The version of the CDB entry.
 * \param firstRun = The first run number for which the CDB entry is valid.
 * \param lastRun = The last run number for which the CDB entry is valid.
 * \param sqrts  The centre of mass energy for the beams.
 * \param beamtype  The beam type as given the AliMagF, eg. p-p, A-A or none.
 */
void CreateCDBFieldIntegrals(
		const char* cdbPath = "local://$ALICE_ROOT/OCDB",
		Int_t version = 0,
		Int_t firstRun = 0,
		Int_t lastRun = AliCDBRunRange::Infinity(),
		Float_t sqrts = 0,
		const char* beamtype = "none"
	)
{
	// Get the CDB manager and storage.
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		cerr << "ERROR: Global CDB manager object does not exist." << endl;
		return;
	}
	AliCDBStorage* storage = cdbManager->GetStorage(cdbPath);
	if (storage == NULL)
	{
		cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
		return;
	}
	
	const int nL3Current = 5;
	const int nDipCurrent = 3;
	Float_t currentL3[nL3Current] = {-30e3, -12e3, 0, 12e3, 30e3};
	Float_t currentDip[nDipCurrent] = {-6e3, 0, 6e3};
	
	// Create and store the configuration parameters for the trigger reconstructor.
	TMap* params = new TMap;
	params->SetOwner(kTRUE);
	for (int i = 0; i < nL3Current; ++i)
	for (int j = 0; j < nDipCurrent; ++j)
	{
		if (currentL3[i]*currentDip[j] < 0) continue; // Skip invalid combinations.
		Float_t bfieldintegral;
		if (! CalculateIntegral(currentL3[i], currentDip[j], bfieldintegral, sqrts, beamtype)) continue;
		const char* paramName = Form("L3_current=%0.2e;Dipole_current=%0.2e", currentL3[i], currentDip[j]);
		const char* paramValue = Form("%8.8f", bfieldintegral);
		params->Add(new TObjString(paramName), new TObjString(paramValue));
	}
	
	const char* path = AliHLTMUONConstants::FieldIntegralsCDBPath();
	AliCDBId id(path, firstRun, lastRun, version);
	AliCDBMetaData* metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Magnetic field integral parameters for dimuon HLT.");
	storage->Put(params, id, metaData);
}
