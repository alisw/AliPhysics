#ifndef ALILHCREADER_H
#define ALILHCREADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class to read the file coming from DCS containing the information        //
//  from LHC.                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TString;
class TMap;
class TObjArray;

#include <TObject.h>

class AliLHCReader : public TObject{
 public:
	AliLHCReader(); // default ctor
	virtual ~AliLHCReader();
	TMap* ReadLHCDP(TString lhcDP);

	UInt_t GetStartTime() const {return fStartTime;}
	UInt_t GetEndTime() const {return fEndTime;}

	void SetStartTime(UInt_t startTime) {fStartTime=startTime;}
	void SetEndTime(UInt_t endTime) {fEndTime = endTime;}
	TObjArray* ReadSingleLHCDP(TString filename, TString alias);

 private:
	AliLHCReader& operator=(const AliLHCReader *reader); // assignment operator
	AliLHCReader(const AliLHCReader *reader); // copy ctor

	UInt_t fStartTime;   // start time of data taking
	UInt_t fEndTime;     // end time of data taking

	ClassDef(AliLHCReader,0)
		};
#endif
