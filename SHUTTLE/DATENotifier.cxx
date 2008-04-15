#include "DATENotifier.h"

#include "AliLog.h"

#include "AliShuttleTrigger.h"

//______________________________________________________________________
void DATENotifier::infoHandler()
{
// Info handler

	Int_t run = -1;
	if (getData())
		run = getInt();
	AliInfoGeneral("DATENotifier::infoHandler()", Form("DATE notification received for run %d...", run));
	fTrigger->Notify();
}

//______________________________________________________________________
void DATENotifier::errorHandler(int severity, int code, char *msg)
{
// Error handler

	AliInfoGeneral("DATENotifier::errorHandler()",
		Form("DIM Error: severity<%d>, code<%d> , message<%s>",
		severity, code, msg));
}
