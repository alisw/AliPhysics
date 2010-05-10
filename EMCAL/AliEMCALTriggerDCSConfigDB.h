#ifndef ALIEMCALTRIGGERDCSCONFIGDB_H
#define ALIEMCALTRIGGERDCSCONFIGDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 adapted from TRD: thanks!
*/

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliCDBEntry;

class AliEMCALTriggerDCSConfig;

class AliEMCALTriggerDCSConfigDB : public TObject 
{
public:
	
	static AliEMCALTriggerDCSConfigDB*  Instance();
	static void                         Terminate();

	void                                SetRun(Long64_t run);
	Long64_t                            GetRun() const { return fRun; }
	const AliEMCALTriggerDCSConfig*     GetTriggerDCSConfig();
	void                                GetSTUSegmentation(Int_t ss[], Int_t sp[]);
	Int_t                               GetTRUGTHRL0(Int_t iTRU);   
	
protected:

	// For caching see also implentation of GetCachedCDBObject in the .cxx file
	// For now, only one cached object but let the possibility to have more
	enum { kCDBCacheSize = 1    };   // Number of cached objects
	enum { kIDTriggerConfig = 0 };   // IDs of cached objects

	const TObject *GetCachedCDBObject(Int_t id);
  
	void           Invalidate();
    
	AliCDBEntry   *GetCDBEntry(const Char_t *cdbPath);
	const TObject *CacheCDBEntry(Int_t id, const Char_t *cdbPath);

	static AliEMCALTriggerDCSConfigDB* fgInstance;                 //  Instance of this class (singleton implementation)
	static Bool_t                      fgTerminated;               //  Defines if this class has already been terminated

	AliCDBEntry*                       fCDBEntries[kCDBCacheSize]; //  Cache for CDB entries
	TObject*                           fCDBCache[kCDBCacheSize];   //  Cache for calibration objects.

	Long64_t                           fRun;                       //  Run Number
  
 private:

	AliEMCALTriggerDCSConfigDB();                                  //  This is a singleton, constructor is private!  
	AliEMCALTriggerDCSConfigDB(const AliEMCALTriggerDCSConfigDB &c);   
	AliEMCALTriggerDCSConfigDB &operator=(const AliEMCALTriggerDCSConfigDB &c); 
	virtual ~AliEMCALTriggerDCSConfigDB();

	ClassDef(AliEMCALTriggerDCSConfigDB, 1)                         //  Provides central access to the CDB
};

#endif

