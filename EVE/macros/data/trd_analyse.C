/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEvePointSet.h>

#include <AliPID.h>
#include <AliTRDhit.h>
#include <AliTRDarrayADC.h>
#include <AliTRDcluster.h>
#include <AliTRDtrackV1.h>
#include <AliTRDReconstructor.h>
#include <AliTRDrecoParam.h>
#include <AliTRDseedV1.h>
#include <AliEveTRDData.h>
#endif

//
// How to access the basic TRD data structures when a 
// pointer to an Eve opject is available. One can get a 
// valid pointer from the event display by accessing the 
// function ExportToCINT 
// 
// Usage:
// .L trd_analyse.C
// analyseXXX(ptr)
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
// 

//_______________________________________________________
void analyseHits(AliEveTRDHits *hits = 0x0)
{
// Simple print hits from a detector

  if(!hits) {
    Info("analyseHits", "Invalid hits set.");
    return;
  }

  AliTRDhit *h = 0x0;
  for(Int_t ih=0; ih<hits->GetN(); ih++){
    hits->PointSelected(ih);
  }
}

/* Obsolete code
 * 
 * 
//_______________________________________________________
void analyseDigits(AliEveTRDDigits *digits = 0x0)
{
// Simple print digits from a detector

  if(!digits) {
    Info("analyseDigits", "Invalid digits set.");
    return;
  }

  Int_t adc ;
  AliTRDdataArrayI *data = digits->GetData();
  Int_t nrows = data->GetNrow(),
        ncols = data->GetNcol(),
        ntbs  = data->GetNtime();
  data->Expand();
  printf("nrows[%d] ncols[%d] ntbs[%d]\n", nrows, ncols, ntbs);
  for (Int_t  row = 0;  row <  nrows;  row++)
  for (Int_t  col = 0;  col <  ncols;  col++)
    for (Int_t time = 0; time < ntbs; time++) {
      if((adc = data->GetDataUnchecked(row, col, time)) <= 1) continue;
      printf("r[%d] c[%d] t[%d] ADC[%d]\n", row, col, time, adc);
    }
  data->Compress(1);
}
*/

//_______________________________________________________
void analyseClusters(TEvePointSet *points = 0x0)
{
// print some info about clusters in one detector or 
// attached to tracks.
  if(!points) {
    Info("analyseClusters", "Invalid clusters set.");
    return;
  }

  AliTRDcluster *c = 0x0;
  for(Int_t ic=0; ic<points->Size(); ic++){
    if(!(c = (AliTRDcluster*)points->GetPointId(ic))) continue;
  
    c->Print();
  }
}

//_______________________________________________________
void analyseTracklet(AliEveTRDTracklet *line)
{
// print tracklet information
  if(!line) {
    Info("analyseTracklet", "Invalid tracklet.");
    return;
  }
  
  AliTRDseedV1 *tracklet = 0x0;
  tracklet = (AliTRDseedV1*)line->GetUserData();
  tracklet->Print();
}

//_______________________________________________________
void analyseTrack(AliEveTRDTrack *line)
{
// print tracklet information
  if(!line) {
    Info("analyseTrack", "Invalid track.");
    return;
  }
  
  AliTRDtrackV1 *track = 0x0;
  track = (AliTRDtrackV1*)line->GetUserData();
  
  AliTRDReconstructor *rec = new AliTRDReconstructor();
  rec->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  track->SetReconstructor(rec);

  rec->SetOption("!nn");
  track->CookPID();
  printf("PID LQ : "); for(int is=0; is<AliPID::kSPECIES; is++) printf("%s[%5.2f] ", AliPID::ParticleName(is), 1.E2*track->GetPID(is)); printf("\n");

  rec->SetOption("nn");
  track->CookPID();
  printf("PID NN : "); for(int is=0; is<AliPID::kSPECIES; is++) printf("%s[%5.2f] ", AliPID::ParticleName(is), 1.E2*track->GetPID(is)); printf("\n");
}

//_______________________________________________________
void trd_analyse()
{
  Info("trd_analyse", "******************************");
  Info("trd_analyse", "Example function which shows how to analyse TRD data exported from the EVE display.");
  Info("trd_analyse", "Usage : Load the macro. Select a TRD data object (digits, clusters, tracklet, track) and call the function \"ExportToCINT()\". Afterwards call the appropiate function (analyseXXX()) on the pointer.");
  Info("trd_analyse", "E.g. If \"tracklet\" is a pointer to a TRD track than one can call :"); 
  Info("trd_analyse", "analyseTrack(track)");

  //gROOT->LoadMacro("$ALICE_ROOT/EVE/macros/trd_analyse.C");
  return;
}
