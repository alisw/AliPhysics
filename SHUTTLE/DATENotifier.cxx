#include "DATENotifier.h"

#include "AliLog.h"

#include "AliShuttleTrigger.h"

//______________________________________________________________________
void DATENotifier::infoHandler()
{
// Info handler

	AliInfoGeneral("DATENotifier::infoHandler()",
			"DATE notification received ...");
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

//______________________________________________________________________
DATENotifier::DATENotifier(const DATENotifier& /*other*/):
DimInfo(), fTrigger(0) {
// copy constructor (not implemented)

}
