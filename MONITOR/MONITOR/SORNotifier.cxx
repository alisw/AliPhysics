#include "SORNotifier.h"

#include "AliLog.h"
#include "AliOnlineRecoTrigger.h"

//______________________________________________________________________
void SORNotifier::infoHandler()
{
// Info handler

	Int_t run = -1;
	if (getData())
		run = getInt();
	AliInfoGeneral("SORNotifier::infoHandler()", Form("ECS SOR notification received for run %d...", run));
	fRun = run;
	fTrigger->Notify();
}

//______________________________________________________________________
void SORNotifier::errorHandler(int severity, int code, char *msg)
{
// Error handler

	AliInfoGeneral("SORNotifier::errorHandler()",
		Form("DIM Error: severity<%d>, code<%d> , message<%s>",
		severity, code, msg));
}
