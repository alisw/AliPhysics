#ifndef ALISTATS_H
#define ALISTATS_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliStats                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TDatime
#include <TDatime.h>
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif


// Forward class declarations
class TH1F;


class AliStats : public TObject {

public:
   AliStats(const char *filename = "", Int_t compmode = 0, Bool_t filter = kFALSE);
   AliStats(const AliStats &rhs);
   virtual ~AliStats();
   AliStats &operator=(const AliStats &rhs);

   void SetEvents(Int_t events) { fEvents = events; }
   void SetFirstId(Int_t run, Int_t event) { fRun = run; fFirstEvent = event; }
   void SetLastId(Int_t event) { fLastEvent = event; }
   void SetBeginTime() { fBegin.Set(); }
   void SetEndTime() { fEnd.Set(); }
   void SetFileSize(Double_t size) { fFileSize = size; }
   void SetCompressionFactor(Float_t comp) { fCompFactor = comp; }
   void Fill(Float_t time);

   Int_t       GetEvents() const { return fEvents; }
   Int_t       GetRun() const { return fRun; }
   Int_t       GetFirstEvent() const { return fFirstEvent; }
   Int_t       GetLastEvent() const { return fLastEvent; }
   TDatime    &GetBeginTime() { return fBegin; }
   TDatime    &GetEndTime() { return fEnd; }
   Double_t    GetFileSize() const { return fFileSize; }
   Int_t       GetCompressionMode() const { return fCompMode; }
   Float_t     GetCompressionFactor() const { return fCompFactor; }
   Bool_t      GetFilterState() const { return fFilter; }
   const char *GetFileName() const { return fFileName; }
   TH1F       *GetRTHist() const { return fRTHist; }

private:
   Int_t    fEvents;     // number of events in this file
   Int_t    fRun;   // run number of first event in file
   Int_t    fFirstEvent; // event number of first event in file
   Int_t    fLastEvent;  // event number of last event in file
   TDatime  fBegin;      // begin of filling time
   TDatime  fEnd;        // end of filling time
   TString  fFileName;   // name of file containing this data
   Double_t fFileSize;   // size of file
   Float_t  fCompFactor; // tree compression factor
   Int_t    fCompMode;   // compression mode
   Bool_t   fFilter;     // 3rd level filter on/off
   TH1F    *fRTHist;     // histogram of real-time to process chunck of data
   Float_t  fChunk;      //!chunk to be histogrammed

   ClassDef(AliStats,2)  // Statistics object
};

#endif
