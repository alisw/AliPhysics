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
#include "TROOT.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "TFile.h"
#include "TObjString.h"

//AliROOT includes
#include "AliLog.h"

ClassImp(AliSurveyObj)
  
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
    for (Int_t i = 0; i < fDataPoints->GetEntries(); ++i) delete fDataPoints->At(i);
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
  TString storage = "alien://alice.cern.ch";

  Printf("TFile::Open string: \n -> \"%s\"\n", openString.Data());

  if (openString.BeginsWith("alien://"))
    if (!Connect(storage.Data(), fGridUser.Data())) {
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
Bool_t AliSurveyObj::Fill(TString detector, Int_t year, Int_t reportNumber,
                         Int_t reportVersion, TString username) {
  // Fills the object from a file in the default storage location in AliEn

  TString baseFolder = "/alice/data/Reference/";
  TString validDetectors = "ACORDE,BABYFRAME,BACKFRAME,EMCAL,FMD,HMPID,ITS,L3 MAGNET,MUON,MUON ABSORBERS,MUON DIPOLE,PHOS,PMD,SPACEFRAME,SUPERSTRUCTURE,T0,TOF,TPC,TRD,V0,ZDC";
  TString GRPDetectors = "BABYFRAME,BACKFRAME,L3 MAGNET,SPACEFRAME,MUON DIPOLE,MUON ABSORBERS";
  TString MUONDetectors = "MUON,SUPERSTRUCTURE";

  detector.ToUpper();
  
  // Check if <detector> is valid
  TObjArray *dets = validDetectors.Tokenize(',');
  if (!dets->FindObject(detector)) {
    AliError(Form("Detector '%s' is not a valid detector/structure!", detector.Data()));
    return kFALSE;
  }
  dets->Delete();
  dets = 0;

  // Some "detectors" don't have a folder of their own
  dets = GRPDetectors.Tokenize(',');
  if (dets->FindObject(detector)) detector = "GRP";
  dets->Delete();
  dets = 0;

  dets = MUONDetectors.Tokenize(',');
  if (dets->FindObject(detector)) detector = "MUON";
  dets->Delete();
  dets = 0;

  // Check if <year>, <reportNumber> and <reportVersion> are valid (roughly)
  if ((year < 1950) || (reportNumber < 1) || (reportVersion < 1)) {
    AliError("Invalid parameter values for AliSurveyObj::Fill. (Year, Report Number or Report Version)");
    return kFALSE;
  }

  // Finally composes the full string
  TString fullOpenString = "alien://" + baseFolder + detector + "/RawSurvey/";
  fullOpenString += Form("%d/%d_v%d.txt?filetype=raw", year, reportNumber, reportVersion);

  // Set the GRID username variable to connect to the GRID
  fGridUser = username;

  return OpenFile(fullOpenString);
}


//_____________________________________________________________________________
Bool_t AliSurveyObj::FillFromLocalFile(const Char_t* filename) {
  // Fills the object from a file in a local filesystem

  TString fullOpenString = "file://" + TString(filename) + "?filetype=raw";

  return OpenFile(fullOpenString);
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
  TString tmp_name = "";
  Float_t tmp_x = 0.0, tmp_y = 0.0, tmp_z = 0.0;
  Float_t tmp_precX = 0.0, tmp_precY = 0.0, tmp_precZ = 0.0;
  Char_t tmp_type = '\0';
  Bool_t tmp_targ = kTRUE;
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
	  
	  Int_t sscanf_tmp = 0;
	  if (1 != sscanf(nextLine.Data(), "%d", &sscanf_tmp)) {
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
	      tmp_name = value;
	      check[0] = kTRUE;
	    } else if (cn.BeginsWith("X", TString::kIgnoreCase)) {
	      tmp_x = value.Atof();
	      check[1] = kTRUE;
	    } else if (cn.BeginsWith("Y", TString::kIgnoreCase)) {
	      tmp_y = value.Atof();
	      check[2] = kTRUE;
	    } else if (cn.BeginsWith("Z", TString::kIgnoreCase)) {
	      tmp_z = value.Atof();
	      check[3] = kTRUE;
	    } else if (cn.BeginsWith("Precision", TString::kIgnoreCase)) {
	      TString tmpCN = cn(0, cn.First('('));
	      Int_t precLength = TString("Precision").Length();
	      //Printf(" ====== %d ======= %d ====== \n", precLength, tmpCN.Length());
	      //Printf(" ====== %s ======= \n", tmpCN.Data());
	      if (precLength == tmpCN.Length()) {
		tmp_precX = tmp_precY = tmp_precZ = value.Atof();
		check[6] = kTRUE;
	      } else {
		TString axis = cn(precLength, tmpCN.Length() - precLength);
		if (axis.Contains('X', TString::kIgnoreCase)) {
		  tmp_precX = value.Atof();
		  check[7] = kTRUE;
		} else if (axis.Contains('Y', TString::kIgnoreCase)) {
		  tmp_precY = value.Atof();
		  check[8] = kTRUE;
		} else if (axis.Contains('Z', TString::kIgnoreCase)) {
		  tmp_precZ = value.Atof();
		  check[9] = kTRUE;
		} else {
		  AliError("Survey text file sintax error! (Precision column name invalid)");
		  dataLine->Delete();
		  lines->Delete();
		  return kFALSE;
		}
	      }
	    } else if (cn.BeginsWith("Point Type", TString::kIgnoreCase)) {
	      tmp_type = value.Data()[0];
	      check[4] = kTRUE;
	    } else if (cn.BeginsWith("Target Used", TString::kIgnoreCase)) {
	      tmp_targ = (value.Data()[0] == 'Y') ? kTRUE : kFALSE;
	      check[5] = kTRUE;
	    }

	    //Printf("--> %s\n", ((TObjString *)(dataLine->At(j)))->GetString().Data());
	  }

	  // Check if all the mandatory fields exist
	  Bool_t res = kTRUE, precInd = kTRUE;

	  // Target
	  if (kFALSE == check[5]) {
	    tmp_targ = kTRUE;
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
	    dp = new AliSurveyPoint(tmp_name, tmp_x, tmp_y, tmp_z, tmp_precX, tmp_precY, tmp_precZ, tmp_type, tmp_targ);
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
    for (Int_t i = 0; i < fDataPoints->GetEntries(); ++i) delete fDataPoints->At(i);
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

