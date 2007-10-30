/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliSurveyObj						   //
//  Retrieve and Convert survey data into ROOT Objects		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliSurveyObj.h"

//ROOT includes
//#include "TROOT.h"
//#include "Riostream.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TFile.h"
#include "TObjString.h"

//AliROOT includes
#include "AliLog.h"
#include "AliSurveyPoint.h"

ClassImp(AliSurveyObj)


const TString AliSurveyObj::fgkStorage = "alien://alice.cern.ch";  
const TString AliSurveyObj::fgkBaseFolder = "/alice/data/Reference";
const TString AliSurveyObj::fgkValidDetectors = "ACORDE,BABYFRAME,BACKFRAME,\
EMCAL,FMD,HMPID,ITS,L3 MAGNET,MUON,MUON ABSORBERS,MUON DIPOLE,PHOS,PMD,\
SPACEFRAME,SUPERSTRUCTURE,T0,TOF,TPC,TRD,V0,ZDC,GRP";
const TString AliSurveyObj::fgkGRPDetectors = "BABYFRAME,BACKFRAME,L3 MAGNET,\
SPACEFRAME,MUON DIPOLE,MUON ABSORBERS,GRP";
const TString AliSurveyObj::fgkMUONDetectors = "MUON,SUPERSTRUCTURE";

  
//_____________________________________________________________________________
AliSurveyObj::AliSurveyObj():
  TObject(),
  fTitle(""),
  fDate(""),
  fDetector(""),
  fURL("http://edms.cern.ch/"),
  fReportNr(-1),
  fVersion(-1),
  fObs(""),
  fCoordSys(""),
  fUnits(""),
  fNrColumns(-1),
  fColNames(""),
  fIsValid(kFALSE),
  fGridUser(""),
  fDataPoints(new TObjArray(1))
{
  // constructor
  fDataPoints->SetOwner(kTRUE);
}


//_____________________________________________________________________________
AliSurveyObj::~AliSurveyObj() {
  //destructor
  if (fDataPoints) {
    fDataPoints->Delete();
    fDataPoints = 0;
  }
}


//_____________________________________________________________________________
AliSurveyObj::AliSurveyObj(const AliSurveyObj& surveyObj):
  TObject(),
  fTitle(surveyObj.fTitle),
  fDate(surveyObj.fDate),
  fDetector(surveyObj.fDetector),
  fURL(surveyObj.fURL),
  fReportNr(surveyObj.fReportNr),
  fVersion(surveyObj.fVersion),
  fObs(surveyObj.fObs),
  fCoordSys(surveyObj.fCoordSys),
  fUnits(surveyObj.fUnits),
  fNrColumns(surveyObj.fNrColumns),
  fColNames(surveyObj.fColNames),
  fIsValid(surveyObj.fIsValid),
  fGridUser(surveyObj.fGridUser),
  fDataPoints(new TObjArray(1))
{
  // copy constructor
  TObject *curr = surveyObj.fDataPoints->First();
  while (curr != 0) {
    fDataPoints->Add(curr);
    curr = surveyObj.fDataPoints->After(curr);
  }
}


//_____________________________________________________________________________
AliSurveyObj& AliSurveyObj::operator=(const AliSurveyObj& surveyObj)
{
  // assignment operator
  if (this != &surveyObj) {
    fTitle = surveyObj.fTitle;
    fDate = surveyObj.fDate;
    fDetector = surveyObj.fDetector;
    fURL = surveyObj.fURL;
    fReportNr = surveyObj.fReportNr;
    fVersion = surveyObj.fVersion;
    fObs = surveyObj.fObs;
    fCoordSys = surveyObj.fCoordSys;
    fUnits = surveyObj.fUnits;
    fNrColumns = surveyObj.fNrColumns;
    fColNames = surveyObj.fColNames;
    fIsValid = surveyObj.fIsValid;
    fGridUser = surveyObj.fGridUser;
    TObject *curr = surveyObj.fDataPoints->First();
    while (curr != 0) {
      fDataPoints->Add(curr);
      curr = surveyObj.fDataPoints->After(curr);
    }
  }
  return *this;
}


//_____________________________________________________________________________
void AliSurveyObj::AddPoint(AliSurveyPoint* point) {
  // Adds a point to the TObjArray which containst the list of points
  fDataPoints->Add(point);
  return;
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::Connect(const char *gridUrl, const char *user) {
  // Connects to the grid

  // If the same "Grid" is alreay active, skip connection
  if (!gGrid || gridUrl != gGrid->GridUrl() ||
      (( user != "" ) && ( user != gGrid->GetUser() )) ) {
    // connection to the Grid
    AliInfo("\nConnecting to the Grid...");
    if (gGrid) {
      AliInfo(Form("gGrid = %x; GridUrl = %s; gGrid->GridUrl() = %s", 
		   gGrid, gridUrl, gGrid->GridUrl()));
      AliInfo(Form("User = %s; gGrid->GetUser() = %s",
		   user, gGrid->GetUser()));
    }
    TGrid::Connect(gridUrl,user);
  }
	
  if(!gGrid) {
    AliError("Connection failed!");
    return kFALSE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::OpenFile(TString openString) {
  // Opens the file and reads it to a buffer

  Printf("Open string: \n -> \"%s\"\n", openString.Data());

  if (openString.BeginsWith("alien://"))
    if (!Connect(fgkStorage.Data(), fGridUser.Data())) {
      AliError(Form("Error connecting to GRID"));
      return kFALSE;
    }

  TFile *file = TFile::Open(openString.Data(), "READ");
  if ( !file ) {
    AliError(Form("Error opening file \"%s\"", openString.Data()));
    return kFALSE;
  }

  Int_t size = file->GetSize();

  char *buf = new Char_t[size + 1];
  memset(buf, '\0', size + 1);

  file->Seek(0);  
  if ( file->ReadBuffer(buf, size) ) {
    AliError("Error reading file contents to buffer!");
    return kFALSE;
  }
  Printf("%d bytes read!\n", file->GetBytesRead());
  
  ParseBuffer(buf);

  file->Close();
  delete[] buf;
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::FillFromLocalFile(const Char_t* filename) {
  // Fills the object from a file in a local filesystem

  TString fullOpenString = "file://" + TString(filename) + "?filetype=raw";

  return OpenFile(fullOpenString);
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::IsValidDetector(TString detector) const {
  // Checks if the detector name is valid

  detector.ToUpper();

  TObjArray *dets = fgkValidDetectors.Tokenize(',');
  TObject *found = dets->FindObject(detector);
  dets->Delete();
  dets = 0;

  if (!found) return kFALSE;
  else return kTRUE;
}


//_____________________________________________________________________________
TString AliSurveyObj::RealFolderName(TString detector) const {
  // Returns the actual folder name for a given detector 
  // Some "detectors" don't have a folder of their own

  detector.ToUpper();
  TString folderName = detector;

  TObjArray *dets = fgkGRPDetectors.Tokenize(',');
  if (dets->FindObject(detector)) folderName = "GRP";
  dets->Delete();
  dets = 0;

  dets = fgkMUONDetectors.Tokenize(',');
  if (dets->FindObject(detector)) folderName = "MUON";
  dets->Delete();  
 


  dets = 0;

  return folderName;
}

//_____________________________________________________________________________
Bool_t AliSurveyObj::Fill(TString detector, Int_t reportNumber,
                          TString username) {
  // Fills the object from a file in the default storage location in AliEn.
  // The highest version available is selected.

  return Fill(detector, reportNumber, -1, username);
}

//_____________________________________________________________________________
Bool_t AliSurveyObj::Fill(TString detector, Int_t reportNumber,
                          Int_t reportVersion, TString username) {
  // Fills the object from a file in the default storage location in AliEn.
  // A specific version is selected.

  detector.ToUpper();
  
  // Check if <detector> is valid
  if (!IsValidDetector(detector)) {
    AliWarning(Form("Detector '%s' is not a valid detector/structure!", detector.Data()));
    return kFALSE;
  }

  // Some "detectors" don't have a folder of their own
  // TString detectorFolder = RealFolderName(detector);

  // Check if <year>, <reportNumber> and <reportVersion> are valid (roughly)
  if ((reportNumber < 1) || (reportVersion < -1) || (0 == reportVersion)) {
    AliError("Invalid parameter values for AliSurveyObj::Fill. (Report Number or Report Version)");
    return kFALSE;
  }

  // Check if the fGridUser is set, or specified
  if (username.Length() > 0) SetGridUser(username);
  else if (0 == fGridUser.Length()) {
    AliError("GRID username not specified and not previously set!");
    return kFALSE;
  }

  // Query AliEn for the available reports
  TGridResult *res = QueryReports(detector, -1, reportNumber, reportVersion);
  if (!res) AliError(Form("Error querying AliEn for detector '%s', \
                           report number '%d' and report version '%d'.",
			  detector.Data(), reportNumber, reportVersion));
  Int_t numberEntries = res->GetEntries();
  if (0 == numberEntries) {
    AliError(Form("No report found for detector '%s', report number '%d' and report version '%d'",
		  detector.Data(), reportNumber, reportVersion));
    return kFALSE;
  }

  TString fileNamePath = "";
  if (1 == numberEntries) fileNamePath = res->GetFileNamePath(0);
  else if (numberEntries > 1) {
    TString higherVerFNP = res->GetFileNamePath(0);
    Int_t lastYear = FileNamePathToReportYear(higherVerFNP);
    for (Int_t i = 1; i < numberEntries; ++i) {
      TString currFNP = res->GetFileNamePath(i);
      if (FileNamePathToReportVersion(currFNP) >
	  FileNamePathToReportVersion(higherVerFNP)) higherVerFNP = currFNP;
      if (lastYear != FileNamePathToReportYear(currFNP))
	AliWarning("Inconsistency detected, year differs for reports with the same report number! Please inform the responsible and check the report against the one in DCDB.");
    }
    fileNamePath = higherVerFNP;
  }

  TString fullOpenString = "alien://" + fileNamePath + "?filetype=raw";
  /*
  // Finally composes the full string
  TString fullOpenString = "alien://" + fgkBaseFolder + "/" + detectorFolder + "/RawSurvey/";
  fullOpenString += Form("%d/%d_v%d.txt?filetype=raw", year, reportNumber, reportVersion);
  */

  return OpenFile(fullOpenString);
}


//_____________________________________________________________________________
TString AliSurveyObj::FileNamePathToDetector(TString filename) const {
  // Get the report number from the complete path in the format:
  // /alice/data/Reference/HMPID/RawSurvey/2006/781282_v1.txt

  TString ret = "";

  if (filename.Length() > fgkBaseFolder.Length()) {
    ret = filename.Remove(0, fgkBaseFolder.Length());
    ret.Remove(TString::kLeading, '/');
    ret = ret(0, ret.First('/'));
    if (!IsValidDetector(ret)) ret = "";
  } 
  return ret;
}

///alice/cern.ch/user/r/rsilva/TRD/RawSurvey/2007/.816582_v2.txt/v1.0

//_____________________________________________________________________________
Int_t AliSurveyObj::FileNamePathToReportYear(TString filename) const {
  // Get the report year from the complete path in the format:
  // /alice/data/Reference/HMPID/RawSurvey/2006/781282_v1.txt

  TString ret = "";

  if (filename.Length() > fgkBaseFolder.Length()) {
    ret = filename.Remove(0, fgkBaseFolder.Length());
    ret.Remove(TString::kLeading, '/');
    Int_t beg = ret.First('/') + TString("RawSurvey/").Length() + 1;
    ret = ret(beg, ret.Last('/') - beg);
    return ret.Atoi();
  } 
  return -1;
}


//_____________________________________________________________________________
Int_t AliSurveyObj::FileNamePathToReportNumber(TString filename) const {
  // Get the report number from the complete path in the format:
  // /alice/data/Reference/HMPID/RawSurvey/2006/781282_v1.txt

  TString ret = "";

  if (filename.Length() > fgkBaseFolder.Length()) {
    ret = filename.Remove(0, fgkBaseFolder.Length());
    ret.Remove(TString::kLeading, '/');
    if ((ret.CountChar('/') > 3) || (ret.CountChar('.') > 1)) {
      AliWarning("Error getting the Report Number from the filename path!");
      return -1;
    }
    ret = ret(ret.Last('/') + 1 , ret.Last('_') - ret.Last('/') - 1);
    return ret.Atoi();
  } 
  AliWarning("Error getting the Report Number from the filename path!");
  return -1;
}


//_____________________________________________________________________________
Int_t AliSurveyObj::FileNamePathToReportVersion(TString filename) const {
  // Get the report version from the complete path in the format:
  // /alice/data/Reference/HMPID/RawSurvey/2006/781282_v1.txt

  TString ret = "";

  if (filename.Length() > fgkBaseFolder.Length()) {
    ret = filename.Remove(0, fgkBaseFolder.Length());
    ret.Remove(TString::kLeading, '/');
    if ((ret.CountChar('/') > 3) || (ret.CountChar('.') > 1)) {
      AliWarning("Error getting the Report Version from the filename path!");
      return -1;
    }
    ret = ret(ret.Last('_') + 1 + 1 , ret.Last('.') - ret.Last('_') - 1 - 1);
    return ret.Atoi();
  } 
  AliWarning("Error getting the Report Version from the filename path!");
  return -1;
}


//_____________________________________________________________________________
void AliSurveyObj::ListValidDetectors() {
  // List the valid detector names
  Printf("");
  Printf("Listing all valid detectors:\n");
  TObjArray *dets = fgkValidDetectors.Tokenize(',');
  for (int i = 0; i < dets->GetEntries(); ++i) 
    Printf("%s", ((TObjString *) dets->At(i))->GetString().Data());
  dets->Delete();
  dets = 0;
  Printf("");
  Printf("Some reports are stored in more general folders.");
  Printf("These reports can be opened using either name, the original or the");
  Printf("folder name. Example: 'SPACEFRAME' or 'GRP' are both valid when");
  Printf("opening a report for the Spaceframe.");
  Printf("");
  Printf("Detectors stored in 'MUON' folder:");
  dets = fgkMUONDetectors.Tokenize(',');
  for (int i = 0; i < dets->GetEntries(); ++i) 
    Printf("%s", ((TObjString *) dets->At(i))->GetString().Data());
  dets->Delete();
  dets = 0;
  Printf("");
  Printf("Detectors stored in 'GRP' folder:");
  dets = fgkGRPDetectors.Tokenize(',');
  for (int i = 0; i < dets->GetEntries(); ++i) 
    Printf("%s", ((TObjString *) dets->At(i))->GetString().Data());
  dets->Delete();
  dets = 0;
  return;
}


//_____________________________________________________________________________
TGridResult * AliSurveyObj::QueryReports(TString detector, Int_t year,
					 Int_t reportNumber,
					 Int_t reportVersion) {
  // Queries AliEn for existing reports matching the specified conditions
  TString lsArg = fgkBaseFolder;
  
  TString detectorFolder = "";
  if (detector.Length() > 0) {
    detector.ToUpper();
    // Check if <detector> is valid
    if (!IsValidDetector(detector)) {
      AliError(Form("Detector '%s' is not a valid detector/structure!",
		    detector.Data()));
      return 0;
    }
    // Some "detectors" don't have a folder of their own
    detectorFolder = "/" + RealFolderName(detector);
  } else detectorFolder = "/*";

  lsArg += detectorFolder + "/RawSurvey";

  TString yearFolder = "";
  if (year > 1950) yearFolder.Form("/%d", year);
  else yearFolder = "/*";

  TString reportFolder = "";
  if (reportNumber > 0) reportFolder.Form("/%d", reportNumber);
  else reportFolder = "/*";

  TString versionFolder = "";
  if (reportVersion > 0) versionFolder.Form("_v%d", reportVersion);
  else versionFolder = "_v*";

  lsArg += yearFolder + reportFolder + versionFolder + ".txt";

  Printf("");
  Printf(Form("Looking for:  %s \n", lsArg.Data()));
  
  // Check if fGridUser is set and Connect to AliEn
  if (0 == fGridUser.Length()) {
    AliError("To use this method it's necessary to call SetGridUser(...) in advance.");
    return 0;
  } else if (!Connect(fgkStorage.Data(), fGridUser.Data())) {
    AliError(Form("Error connecting to GRID"));
    return 0;
  }
  return gGrid->Ls(lsArg);
}


//_____________________________________________________________________________
Int_t AliSurveyObj::ListReports(TString detector, Int_t year, 
				Int_t reportNumber,
				Int_t reportVersion) {
  // Lists all available reports matching the specified conditions
  // Returns the number of reports found

  TGridResult *res = QueryReports(detector, year, reportNumber, reportVersion);
  
  if (0 == res) {
    AliError("Query failed.");
    return 0;
  }

  TString fn = "";
  Int_t numberEntries = res->GetEntries();

  if (numberEntries > 0) {
    Printf("");
    Printf(Form("%d reports found:", numberEntries));
    for (int i = 0; i < res->GetEntries(); ++i) {
      fn = res->GetFileNamePath(i);
      Printf(Form("Detector:%s\tYear:%d\tEDMS Report Number:%d\tVersion:%d",
		  FileNamePathToDetector(fn).Data(),
		  FileNamePathToReportYear(fn),
		  FileNamePathToReportNumber(fn),
		  FileNamePathToReportVersion(fn)));
    }
    delete res;
    return numberEntries;
  } else {
    AliInfo("No results found for the requested query.");
    delete res;
    return 0;
  }
}


//_____________________________________________________________________________
void AliSurveyObj::SetGridUser(TString username){
  // Set the username used to connect to the GRID
  fGridUser = username;
  return;
}


//_____________________________________________________________________________
TString &AliSurveyObj::Sanitize(TString str) {
  // Cleans up a line of new line and carriage return characters.
  // (Specially usefull for files created in Windows.)

  str.Remove(TString::kTrailing, '\r');
  str.Remove(TString::kTrailing, '\n');
  str.Remove(TString::kTrailing, '\r');
  if (!str.IsAscii()) {
    AliWarning("Warning: Non-ASCII characters!\n");
    str = "";
  }
  return str.Remove(TString::kBoth, ' ');
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::ParseBuffer(const Char_t* buf) {
  // Parses a character buffer assuming the format defined with the TS/SU
  // http://aliceinfo/Offline/Activities/Alignment/SurveyInformation.html

  // If the object is already filled clean it up
  if (fIsValid) Reset();

  // Copy the buffer to a TString to use Tokenize
  TString buffer = TString(buf);
  TObjArray *lines = buffer.Tokenize('\n');
  TObjArray *dataLine = NULL; // Used to Tokenize each point/line read
  TObjArray *colLine = NULL; // Used to Tokenize the column names

  // Some local variables declarations and initializations
  const Int_t kFieldCheck = 10;
  Bool_t check[kFieldCheck];
  TString tmpname = "";
  Float_t tmpx = 0.0, tmpy = 0.0, tmpz = 0.0;
  Float_t tmpprecX = 0.0, tmpprecY = 0.0, tmpprecZ = 0.0;
  Char_t tmptype = '\0';
  Bool_t tmptarg = kTRUE;
  AliSurveyPoint *dp = 0;
  for (Int_t i = 0; i < kFieldCheck; ++i) check[i] = kFALSE;

  Int_t nrLines = lines->GetEntries();
  Printf("Lines in file: %d\n", nrLines); 

  // The main cycle, the buffer is parsed a line at a time
  TString currLine = "", nextLine = "";
  for (Int_t i = 0; i < nrLines; ++i) {

    // Get the next line
    currLine = ((TObjString *)(lines->At(i)))->GetString().Data();
    nextLine = ((i + 1) < nrLines ? ((TObjString *)(lines->At(i + 1)))->GetString().Data() : "");
    currLine = Sanitize(currLine);
    nextLine = Sanitize(nextLine);
    // Printf("\n%d: \"\"%s\"\"\"\n", i + 1, currLine.Data());
    
    // Skip empty line
    if (0 == currLine.Length()) Printf("Info: Empty line skipped\n\n");

    // The line contains a keyword
    else if (currLine.BeginsWith(">") && !nextLine.BeginsWith(">")) {
      currLine.Remove(TString::kLeading, '>');
      currLine.Remove(TString::kTrailing, ':');
      currLine.Remove(TString::kBoth, ' ');
      nextLine.Remove(TString::kBoth, ' ');
      // Printf(" -> Field: \"%s\"\n", currLine.Data());
      
      if (currLine.BeginsWith("Title", TString::kIgnoreCase)) {
	// Report Title
	fTitle = nextLine;
	++i;
      } else if (currLine.BeginsWith("Date", TString::kIgnoreCase)) {
	// Report(measurement) Date
	fDate = nextLine;
	++i;
      } else if (currLine.BeginsWith("Subdetector", TString::kIgnoreCase)) {
	// Subdetector or structure
	fDetector = nextLine;
	++i;
      } else if (currLine.BeginsWith("Report URL", TString::kIgnoreCase)) {
	// Report URL in EDMS
	if (nextLine.BeginsWith("http://edms.cern.ch/document/", TString::kIgnoreCase) ||
	    nextLine.BeginsWith("https://edms.cern.ch/document/", TString::kIgnoreCase)) {
 	  fURL = nextLine;
	  nextLine.Remove(TString::kTrailing, '/');
	  nextLine = nextLine(nextLine.Last('/') + 1, nextLine.Length() - nextLine.Last('/') + 1);
	  
	  Int_t sscanftmp = 0;
	  if (1 != sscanf(nextLine.Data(), "%d", &sscanftmp)) {
	    AliError("Survey text file sintax error! (incorrectly formated Report URL)");
	    lines->Delete();
	    return kFALSE;
	  }
	  fReportNr = nextLine.Atoi();
	  //Printf(" $$ %d $$\n", fReportNr);
	  ++i;
	} else { 
	  // URL incorrectly formated
	  AliError("Survey text file sintax error! (incorrectly formated Report URL)");
	  return kFALSE;
	}
      } else if (currLine.BeginsWith("Version", TString::kIgnoreCase)) {
	// Report version
	if (!nextLine.IsDigit()) {
	  lines->Delete();
	  AliError("Survey text file sintax error! (incorrectly formated Report Version)");
	  return kFALSE;
	}
	fVersion = nextLine.Atoi();
	++i;
      } else if (currLine.BeginsWith("General Observations", TString::kIgnoreCase)) {
	// Observations
	fObs = "";
	// Can be more than 1 line. Loop until another keyword is found
	while (('>' != nextLine[0]) && (nextLine.Length() > 0) && (i < nrLines)) {	
	  fObs += (0 == fObs.Length()) ? nextLine : " / " + nextLine;
	  ++i;
	  nextLine = ((i + 1) < nrLines ? ((TObjString *)(lines->At(i + 1)))->GetString().Data() : "");
	  nextLine = Sanitize(nextLine);
	}
      } else if (currLine.BeginsWith("Coordinate System", TString::kIgnoreCase)) {
	// Coordinate System
	fCoordSys = nextLine;
	++i;
      } else if (currLine.BeginsWith("Units", TString::kIgnoreCase)) {
	// Measurement Unit
	fUnits = nextLine;
	++i;
      } else if (currLine.BeginsWith("Nr Columns", TString::kIgnoreCase)) {
	// Number of columns in the "Data" section
	if (!nextLine.IsDigit()) {
	  lines->Delete();
	  AliError("Survey text file sintax error! (incorrectly formated Number of Columns)");
	  return kFALSE;
	}
	fNrColumns = nextLine.Atoi();
	++i;
      } else if (currLine.BeginsWith("Column Names", TString::kIgnoreCase)) {
	// Column names separated by commas
	fColNames = nextLine;
	colLine = nextLine.Tokenize(',');
	if (colLine->GetEntries() != fNrColumns) {
	  AliError("Survey text file sintax error! (Declared number of Columns doesn't match number of column names)");
	  colLine->Delete();
	  lines->Delete();
	  return kFALSE;
	}
	++i;
      } else if (currLine.BeginsWith("Data", TString::kIgnoreCase)) {
	// Data section!
	while ((nextLine.Length() > 0) && ('>' != nextLine[0])) {

	  // Printf("Data LINE: \"%s\": %d\n", nextLine.Data(), nextLine.Length());

	  // What is the separator used between fields?
	  // The allowed are: comma (','), tab ('\t'), and space (' ')
          if (fNrColumns == nextLine.CountChar(',') + 1) dataLine = nextLine.Tokenize(',');
          else if (fNrColumns == nextLine.CountChar('\t') + 1) dataLine = nextLine.Tokenize('\t');
          else if (fNrColumns == nextLine.CountChar(' ') + 1) dataLine = nextLine.Tokenize(' ');
          else {
            // Error (No separator was found!)
            AliError("Survey text file syntax error! Error processing data line!");
            lines->Delete();
            return kFALSE;
          }
	  
	  if (dataLine->GetEntries() != fNrColumns) {
	    // The number of columns doesn't match the number specified in the header
	    AliError("Survey text file sintax error! (Number of entries in line is different from declared Number of Columns)");
	    dataLine->Delete();
	    lines->Delete();
	    return kFALSE;
	  }

	  // Reset the bits used to check if all the required fields are present
	  for (Int_t t = 0; t < kFieldCheck; ++t) check[t] = 0;

	  // Process the data line using the column names as index
	  for (Int_t j = 0; j < dataLine->GetEntries(); ++j) {
	    TString cn = ((TObjString *)(colLine->At(j)))->GetString();
	    TString value = ((TObjString *)(dataLine->At(j)))->GetString();
	    if (cn.BeginsWith("Point Name", TString::kIgnoreCase)) {
	      tmpname = value;
	      check[0] = kTRUE;
	    } else if (cn.BeginsWith("X", TString::kIgnoreCase)) {
	      tmpx = value.Atof();
	      check[1] = kTRUE;
	    } else if (cn.BeginsWith("Y", TString::kIgnoreCase)) {
	      tmpy = value.Atof();
	      check[2] = kTRUE;
	    } else if (cn.BeginsWith("Z", TString::kIgnoreCase)) {
	      tmpz = value.Atof();
	      check[3] = kTRUE;
	    } else if (cn.BeginsWith("Precision", TString::kIgnoreCase)) {
	      TString tmpCN = cn(0, cn.First('('));
	      Int_t precLength = TString("Precision").Length();
	      //Printf(" ====== %d ======= %d ====== \n", precLength, tmpCN.Length());
	      //Printf(" ====== %s ======= \n", tmpCN.Data());
	      if (precLength == tmpCN.Length()) {
		tmpprecX = tmpprecY = tmpprecZ = value.Atof();
		check[6] = kTRUE;
	      } else {
		TString axis = cn(precLength, tmpCN.Length() - precLength);
		if (axis.Contains('X', TString::kIgnoreCase)) {
		  tmpprecX = value.Atof();
		  check[7] = kTRUE;
		} else if (axis.Contains('Y', TString::kIgnoreCase)) {
		  tmpprecY = value.Atof();
		  check[8] = kTRUE;
		} else if (axis.Contains('Z', TString::kIgnoreCase)) {
		  tmpprecZ = value.Atof();
		  check[9] = kTRUE;
		} else {
		  AliError("Survey text file sintax error! (Precision column name invalid)");
		  dataLine->Delete();
		  lines->Delete();
		  return kFALSE;
		}
	      }
	    } else if (cn.BeginsWith("Point Type", TString::kIgnoreCase)) {
	      tmptype = value.Data()[0];
	      check[4] = kTRUE;
	    } else if (cn.BeginsWith("Target Used", TString::kIgnoreCase)) {
	      tmptarg = (value.Data()[0] == 'Y') ? kTRUE : kFALSE;
	      check[5] = kTRUE;
	    }

	    //Printf("--> %s\n", ((TObjString *)(dataLine->At(j)))->GetString().Data());
	  }

	  // Check if all the mandatory fields exist
	  Bool_t res = kTRUE, precInd = kTRUE;

	  // Target
	  if (kFALSE == check[5]) {
	    tmptarg = kTRUE;
	    check[5] = kTRUE;
	  }
	  
	  // Individual axis precisions
	  for (Int_t t = 7; t < 10; ++t) precInd &= check[t];
	  if ((kFALSE == check[6]) && (kTRUE == precInd)) check[6] = kTRUE;

	  for (Int_t t = 0; t < kFieldCheck - 3; ++t) {
	    //Printf("RES(%d): %d\n", t, check[t]);
	    res &= check[t];
	  }
	  if (kTRUE == res) {
	    dp = new AliSurveyPoint(tmpname, tmpx, tmpy, tmpz, tmpprecX, tmpprecY, tmpprecZ, tmptype, tmptarg);
	    dp->PrintPoint();
	    AddPoint(dp);
	  } else {
	    AliError("Parsing error processing data line!");
	    dataLine->Delete();
	    lines->Delete();
	    return kFALSE;
	  }
	  
	  dataLine->Delete();
	  dataLine = NULL;
	  ++i;
	  nextLine = ((i + 1) < nrLines ? ((TObjString *)(lines->At(i + 1)))->GetString().Data() : "");
	  nextLine = Sanitize(nextLine);
	}
      }
    } else {
      AliError("Survey text file sintax error!");
      lines->Delete();
      return kFALSE;
    }
  }
  lines->Delete();
  fIsValid = kTRUE;
  return kTRUE;
}

//_____________________________________________________________________________
void AliSurveyObj::Reset() {
  // Resets the AliSurveyObj to a clean object.
  // Used if the same object is filled more than once
  
  if (fDataPoints) {
    fDataPoints->Delete();
    fDataPoints = 0;
  }
  fTitle = "";
  fDate = "";
  fDetector = "";
  fURL = "http://edms.cern.ch/";
  fReportNr = -1;
  fVersion = -1;
  fObs = "";
  fCoordSys = "";
  fUnits = "";
  fNrColumns = -1;
  fColNames = "";
  fIsValid = kFALSE;
  fDataPoints = new TObjArray(1);
  fDataPoints->SetOwner(kTRUE);
  return;
}

