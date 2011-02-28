//
// Class AliRsnTOFCorrectionESD
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"

#include "AliRsnTOFCorrectionESD.h"

ClassImp(AliRsnTOFCorrectionESD)

//Bool_t         AliRsnTOFCorrectionESD::fgTOFcalibrateESD = kTRUE;
Bool_t         AliRsnTOFCorrectionESD::fgTOFcorrectTExp  = kTRUE;
Bool_t         AliRsnTOFCorrectionESD::fgTOFuseT0        = kTRUE;
Bool_t         AliRsnTOFCorrectionESD::fgTOFtuneMC       = kFALSE;
Double_t       AliRsnTOFCorrectionESD::fgTOFresolution   = 100.0;
AliTOFT0maker* AliRsnTOFCorrectionESD::fgTOFmaker        = 0x0;
AliTOFcalib*   AliRsnTOFCorrectionESD::fgTOFcalib        = 0x0;
Int_t          AliRsnTOFCorrectionESD::fgLastRun         = -1;

//_________________________________________________________________________________________________
AliRsnTOFCorrectionESD::AliRsnTOFCorrectionESD(Bool_t isMC, Double_t tofRes)
   : fOwnESDpid(kFALSE), fESDpid(0x0)
{
//
// Default constructor.
//

   fgTOFtuneMC     = isMC;
   fgTOFresolution = tofRes;
}

//_________________________________________________________________________________________________
AliRsnTOFCorrectionESD::AliRsnTOFCorrectionESD(const AliRsnTOFCorrectionESD& copy)
   : TObject(copy), fOwnESDpid(copy.fOwnESDpid), fESDpid(copy.fESDpid)
{
//
// Copy constructor
//
}

//_________________________________________________________________________________________________
AliRsnTOFCorrectionESD& AliRsnTOFCorrectionESD::operator=(const AliRsnTOFCorrectionESD& copy)
{
//
// Assignment operator
//

   fOwnESDpid = copy.fOwnESDpid;
   fESDpid    = copy.fESDpid;

   return (*this);
}

//_________________________________________________________________________________________________
void AliRsnTOFCorrectionESD::ProcessEvent(AliESDEvent *esd)
{
//
// Repeats the PID for current event.
// In order to avoid to repeat the initialization of calib object
// when this is not needed, this function uses the data-members
// to check if the run has changed.
//

   // compare run number with static data member
   Int_t run = esd->GetRunNumber();

   // initialize only if run number has changed
   if (run != fgLastRun) {
      AliInfo("============================================================================================");
      AliInfo(Form("*** CHANGING RUN NUMBER: PREVIOUS = %d --> CURRENT = %d ***", fgLastRun, run));
      AliInfo("============================================================================================");
      fgLastRun = run;

      AliCDBManager::Instance()->SetDefaultStorage("raw://");
      AliCDBManager::Instance()->SetRun(fgLastRun);

      if (fgTOFmaker) delete fgTOFmaker;
      if (fgTOFcalib) delete fgTOFcalib;

      fgTOFcalib = new AliTOFcalib();
      if (fgTOFtuneMC) {
         fgTOFcalib->SetRemoveMeanT0(kFALSE);
         fgTOFcalib->SetCalibrateTOFsignal(kFALSE);
      } else {
         fgTOFcalib->SetRemoveMeanT0(kTRUE);
         fgTOFcalib->SetCalibrateTOFsignal(kTRUE);
      }
      if (fgTOFcorrectTExp) fgTOFcalib->SetCorrectTExp(kTRUE);
      fgTOFcalib->Init();

      fgTOFmaker = new AliTOFT0maker(fESDpid, fgTOFcalib);
      fgTOFmaker->SetTimeResolution(fgTOFresolution);
   }

   // if the ESDpid object is not present, create it
   if (!fESDpid) {
      fESDpid = new AliESDpid;
      fOwnESDpid = kTRUE;
   }

   // repeat the calibration and PID computations
   /*if (fgTOFcalibrateESD)*/ fgTOFcalib->CalibrateESD(esd);
   if (fgTOFtuneMC) fgTOFmaker->TuneForMC(esd);
   if (fgTOFuseT0) {
      fgTOFmaker->ComputeT0TOF(esd);
      fgTOFmaker->ApplyT0TOF(esd);
      fESDpid->MakePID(esd, kFALSE, 0.);
   }
}
