#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TStopwatch.h>
  #include <TGrid.h>
  #include <TString.h>
  #include <TMath.h>

  #include "AliRawReader.h"
  #include "AliCaloRawStreamV3.h"
  #include "AliLog.h"
  #include "AliPHOSRawFitterv1.h"
  #include "AliCentralTrigger.h"
  #include "AliTriggerConfiguration.h"
  #include "AliTriggerClass.h"
  #include "AliCDBManager.h"
#endif

//-----------------------------------------------------------------------------
static Bool_t              firstEvent = kTRUE;
static Int_t               runNum;
static UInt_t              period;
static UInt_t              orbitID;
static UInt_t              bcID;
static AliRawReader       *reader;
TString GetTriggerClass(ULong64_t);

//-----------------------------------------------------------------------------
void ReadRawCaloTrigger(const TString fileName = "./", const TString calo="PHOS")
{

  if (fileName.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");

  reader = AliRawReader::Create(fileName);
  reader ->Reset();
  AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,calo);

  Int_t iev = 0;

  TStopwatch stopwatch;
  stopwatch.Start();
  Int_t    module,cellX,cellZ,caloFlag;

  while (reader->NextEvent()) {
    if (firstEvent) {
      firstEvent = kFALSE;
      runNum = reader->GetRunNumber();
      man = AliCDBManager::Instance();
      man ->SetRun(runNum);
    }
    ULong64_t triggerMask  = reader->GetClassMask();
    TString trclasses = GetTriggerClass(triggerMask);
    
    period  = reader->GetPeriod();
    orbitID = reader->GetOrbitID();
    bcID    = reader->GetBCID();
    iev++;
    if (!trclasses.Contains("CPHI1-B-NOPF-ALLNOTRD")) continue;

    AliInfoGeneral("",Form("Reading event %d of type %d, time %d, trig.class \"%s\"",
			   iev,reader->GetType(), reader->GetTimestamp(), trclasses.Data()));

    while (stream->NextDDL()) {
      while (stream->NextChannel()) {
	// Get geometry of the channel
	module   = stream->GetModule();
	cellX    = stream->GetCellX();
	cellZ    = stream->GetCellZ();
	caloFlag = stream->GetCaloFlag();

	Int_t nBunches = 0;
	while (stream->NextBunch()) {
	  const UShort_t *sig = stream->GetSignals();
	  Int_t sigLength = stream->GetBunchLength();
	  Int_t maxAmp  = 0;
	  for (Int_t i = 0; i < sigLength; i++) {
	    if (sig[i] >  maxAmp) {
	      maxAmp = sig[i];
	    }
	  }
	  printf("(m,x,z,c)=(%d,%2d,%2d,%d), bunch=%d, amplutide=%d\n",
		 module,cellX,cellZ,caloFlag,nBunches,maxAmp);
	  nBunches++;
	}
      }
    }
    stream->Reset();
  }
  stopwatch.Print();

}
//-----------------------------------------------------------------------------

TString GetTriggerClass(ULong64_t triggerMask)
{
  // Convert a trigger mask to a trigger class
  
  AliCentralTrigger aCTP;
  TString configstr("");
  TString trclasses;
  if (!aCTP.LoadConfiguration(configstr)) { // Load CTP config from OCDB
    AliInfoGeneral("","No trigger configuration found in OCDB! The trigger configuration information will not be used!");
    return trclasses;
  }
  aCTP.SetClassMask(triggerMask);
  AliTriggerConfiguration *config = aCTP.GetConfiguration();
  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
    AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
    if (trclass && trclass->GetMask()>0) {
      Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
      reader->LoadTriggerClass(trclass->GetName(),trindex);
      if (triggerMask & (1ull << trindex)) {
	trclasses += " ";
	trclasses += trclass->GetName();
	trclasses += " ";
      }
    }
  }
  return trclasses;
}
