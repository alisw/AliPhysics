// $Id$
/**
 * @file makeGRPObject.C
 * @brief Creation of the GRP configuration object for HLT 
 *
 * <pre>
 * Usage: aliroot -b -q -l makeGRPObject.C'("arguments", "uri", \
 *                                          rangemin, rangemax)'
 * 
 * 
 * Usage: aliroot -b -q -l makeGRPObject.C'(l3Current, l3Polarity,         \
 *                                          dipoleCurrent, dipolePolarity, \
 *                                          "uri", rangemin, rangemax)'
 * </pre>
 *
 * Create the GRP OCDB object for the magnetic field settings.
 * The magnetic field in AliRoot is controlled by the AliMagF class and
 * the TGeoGlobalMagField container. AliMagF contains the field maps
 * for both the L3 and dipole magnets.
 * 
 * The following conventions for magnetic field polarity are known:
 * 1) kConvMap2005: used for the field mapping in 2005
 *    positive L3  current -> negative Bz
 *    positive Dip current -> positive Bx 
 * 2) kConvMapDCS2008: defined by the microswitches/cabling of power
 *                     converters as of 2008 - 1st half 2009
 *    positive L3  current -> positive Bz
 *    positive Dip current -> positive Bx
 * 3) kConvLHC : defined by LHC
 *    positive L3  current -> positive Bz
 *    positive Dip current -> negative Bx
 * 
 * For current data taking, only the last convention is relevant.
 *
 * Negative polarity corresponds to the value '1' in the GRP, while
 * positive to '0'.
 *
 * Parameters: <br>
 *      arguments      off, default, l3=current, l3=off,
 *                     dipole=current, dipole=off
 *      uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB
 *      rangemin (opt) default 0
 *      rangemax (opt) default 999999999
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void makeGRPObject(float l3Current,
		   int   l3Polarity,
		   float dipoleCurrent,
		   int   dipolePolarity,
		   const char* cdbUri=NULL,
		   int runmin=0,
		   int runmax=999999999,
		   bool bVerbose=true)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "can not get AliCDBManager" << end;
    exit;
  }
  TString storage;
  if (!man->IsDefaultStorageSet()) {
    if (cdbUri) {
      storage=cdbUri;
      if (storage.Contains("://")==0) {
	storage="local://"; storage+=cdbUri;
      }
    } else {
      storage="local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  // generate GRP object
  AliGRPObject* grpObj=new AliGRPObject;
  float cmsEnergy=14000;
  grpObj->SetBeamEnergy(cmsEnergy/0.120); // LHC convention
  grpObj->SetBeamType("p-p");
  grpObj->SetL3Current(l3Current,(AliGRPObject::Stats)0);
  grpObj->SetDipoleCurrent(dipoleCurrent,(AliGRPObject::Stats)0);  
  grpObj->SetL3Polarity(l3Polarity);  
  grpObj->SetDipolePolarity(dipolePolarity);
  grpObj->SetPolarityConventionLHC();                    // LHC convention +/+ current -> -/- field main components

  if (bVerbose) {
    AliGRPManager grpman;
    grpman.SetGRPEntry(grpObj);
    grpman.SetMagField();
  }

  // write object to OCDB
  AliCDBPath cdbPath("GRP/GRP/Data");
  AliCDBId cdbId(cdbPath, runmin, runmax);
  AliCDBMetaData cdbMetaData;
  cdbMetaData.SetResponsible("ALICE HLT");
  cdbMetaData.SetComment("Automatically produced GRP entry (AliHLTSimulation) for the magnetic field initialization of HLT components");
  man->Put(grpObj, cdbId, &cdbMetaData);
}

void makeGRPObject(const char* arguments="", 
		   const char* cdbUri=NULL,
		   int runmin=0,
		   int runmax=999999999)

{
  TString args=arguments;
  bool bHelp=args.IsNull();
  TObjArray* pTokens=args.Tokenize(" ");
  float l3Current=0;
  int l3Polarity=0;
  float dipoleCurrent=0;
  int dipolePolarity=0;
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      TString arg=((TObjString*)pTokens->At(i))->GetString();
      if (arg.CompareTo("off", TString::kIgnoreCase)==0) {
	  l3Current=0;
	  l3Polarity=0;
	  dipoleCurrent=0;
	  dipolePolarity=0;
      } else if (arg.CompareTo("default", TString::kIgnoreCase)==0) {
	  l3Current=30000;
	  l3Polarity=0;  // positive for positive field
	  dipoleCurrent=6000;
	  dipolePolarity=1; // negative for positive field
      } else if (arg.BeginsWith("l3=")) {
	arg.ReplaceAll("l3=", "");
	if (arg.CompareTo("off", TString::kIgnoreCase)==0) {
	  l3Current=0;
	  l3Polarity=0;
	} else if (arg.IsFloat()) {
	  l3Current=arg.Atof();
	  l3Polarity=l3Current<0 ? 1:0;
	  l3Current=TMath::Abs(l3Current);
	} else {
	  cerr << "Invalid parameter for key 'l3=', allowed is 'off' or current" << endl;
	}
	
      } else if (arg.BeginsWith("dipole=")) {
	arg.ReplaceAll("dipole=", "");
	if (arg.CompareTo("off", TString::kIgnoreCase)==0) {
	  dipoleCurrent=0;
	  dipolePolarity=0;
	} else if (arg.IsFloat()) {
	  dipoleCurrent=arg.Atof();
	  dipolePolarity=dipoleCurrent<0 ? 1:0;
	  dipoleCurrent=TMath::Abs(dipoleCurrent);
	} else {
	  cerr << "Invalid parameter for key 'dipole=', allowed is 'off' or current" << endl;
	}
      } else {
	if (arg.CompareTo("--help", TString::kIgnoreCase) &&
	    arg.CompareTo("-h", TString::kIgnoreCase)) {
	  cerr << "Unknown argument: " << arg << endl;
	}
	bHelp=true;
      }
    }
    delete pTokens;
  }

  if (!bHelp) {
    makeGRPObject(l3Current, l3Polarity, dipoleCurrent, dipolePolarity, cdbUri, runmin, runmax);
  } else {
    cout << "========================================================================" << endl;
    cout << "usage: aliroot -b -q -l makeGRPObject.C'(\"arguments\", \"uri\", rangemin, rangemax)'" << endl << endl;
    cout << "  arguments      off, default, l3=<current>, l3=off," << endl;
    cout << "                 dipole=<current>, dipole=off" << endl;
    cout << "  uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   " << endl;
    cout << "  rangemin (opt) default 0" << endl;
    cout << "  rangemax (opt) default 999999999" << endl << endl;
    cout << "alternative usage: aliroot -b -q -l makeGRPObject.C'(l3Current, l3Polarity, \\" << endl;
    cout << "                                         dipoleCurrent, dipolePolarity, \\" << endl;
    cout << "                                         \"uri\", rangemin, rangemax)'" << endl << endl;
    cout << "========================================================================" << endl;
  }
}
