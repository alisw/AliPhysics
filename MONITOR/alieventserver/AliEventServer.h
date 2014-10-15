// Author: Mihai Niculesu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEventServer_H
#define AliEventServer_H

#include <TObject.h>
#include <TString.h>

class AliEventServerReconstruction;
class AliDimIntNotifier;

class AliEventServer : public TQObject
{
public:
	AliEventServer();
	virtual ~AliEventServer();

	void StartOfRun(Int_t run);
	void EndOfRun(Int_t run);
private:
	void InitDIMListeners();
	void FillRunsFromDatabase();
	
	AliDimIntNotifier *fDimSORListener[5];
	AliDimIntNotifier *fDimEORListener[5];

	AliEventServerReconstruction* fRecoServer;

	AliEventServer(const AliEventServer&);
	AliEventServer& operator=(const AliEventServer&);
	ClassDef(AliEventServer, 0);
};

#endif
