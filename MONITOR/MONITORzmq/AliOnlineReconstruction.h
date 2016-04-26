// Author: Mihai Niculesu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliOnlineReconstruction_H
#define AliOnlineReconstruction_H

#include <AliReconstruction.h>
#include <AliCDBManager.h>
#include <AliGRPPreprocessor.h>

#include <TString.h>
#include <TEnv.h>

class AliOnlineReconstruction
{
public:
	AliOnlineReconstruction(int run);
	~AliOnlineReconstruction();
private:
	void StartOfRun();
	void EndOfRun();
	void FillRunsFromDatabase();
	int  RetrieveGRP(TString &gdc);
	void SetupReco();
	void ReconstructionLoop();

	int fRun;
	TString fDataSource;
	TEnv fSettings;

	AliReconstruction *fAliReco;
	AliCDBManager *fCDBmanager;

	AliOnlineReconstruction(const AliOnlineReconstruction&);
	AliOnlineReconstruction& operator=(const AliOnlineReconstruction&);
};

#endif
