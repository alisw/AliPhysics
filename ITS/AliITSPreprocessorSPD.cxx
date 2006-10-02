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

/*
$Log$
Revision 1.3  2006/04/12 08:32:31  hristov
New SPD simulation (Massimo):
- speeding up of the diffusion code (Bjorn)
- creation  of CDB file with the dead channels, implementation
of the CDB reading, check of  the code (Henrik, Giuseppe, Domenico)
- final tuning of the diffusion model parameters (Romualdo)

Revision 1.2  2006/02/03 11:31:18  masera
Calibration framework improved (E. Crescio)

Revision 1.1  2005/10/11 12:31:50  masera
Preprocessor classes for SPD (Paul Nilsson)

*/

///////////////////////////////////////////////////////////////////////////
// AliITSPreprocessorSPD implementation by P. Nilsson 2005
// MAIN AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// The purpose of this class is to provide algorithms for identification
// of "bad channels" such as dead channels and noisy channels in the SPD
//
// Examples on how to use this class can be found in the macros:
//
// .findNoisyChannels.C (Locate and store noisy channels in the CDB)
// .readNoisyChannels.C (Read noisy channels from the CDB)
// .findDeadChannels.C  (Locate and store dead channels in the CDB)
// .readDeadChannels.C  (Read dead channels from the CDB)
//
// Modified by D. Elia, H. Tydesjo
// March 2006: Mixed up coordinates, bug fixed
//
///////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "AliITSPreprocessorSPD.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSPD.h" 
ClassImp(AliITSPreprocessorSPD)
//__________________________________________________________________________
AliITSPreprocessorSPD::AliITSPreprocessorSPD(void):
fITSLoader(0x0),
fRunLoader(0x0),
fThresholdRatio(5.),
fThreshold(5),
fMaximumNumberOfEvents(1000000),
fNumberOfModules(0),
fHighestModuleNumber(0),
fNumberOfColumns(0),
fNumberOfRows(0),
fSelectedAlgorithm(kOptimizedForRealData),
fGeometryMode(kALICEGeometry),
fNumberOfBadChannels(0),
  fIndex(0),
fInit(kFALSE),
fVMEMode(kFALSE),
fDigitsHistogram(0),
fBadChannelsObjArray(0),
fBadChannelsIntArray(0),
fBadChannelsIndexArray(0),
fBadChannelsContainer(0)
{
  // Default constructor for the SPD preprocessor
  //
  // Initialization has to be done by hand using Set* methods and the Open method
  //
  // Input :
  // Output:
  // Return: <empty/uninitialized AliITSPreprocessorSPD object>
}


//__________________________________________________________________________
AliITSPreprocessorSPD::AliITSPreprocessorSPD(const char *fileName, const char *mode,
					     const char *fileNameg, const Int_t maxNumberOfEvents):
  fITSLoader(0x0),
  fRunLoader(0x0),
  fThresholdRatio(0),
  fThreshold(0),
  fMaximumNumberOfEvents(1000000),
  fNumberOfModules(0),
  fHighestModuleNumber(0),
  fNumberOfColumns(0),
  fNumberOfRows(0),
  fSelectedAlgorithm(kOptimizedForRealData),
  fGeometryMode(kALICEGeometry),
  fNumberOfBadChannels(0),
  fIndex(0),
  fInit(kFALSE),
  fVMEMode(kFALSE),
  fDigitsHistogram(0),
  fBadChannelsObjArray(0),
  fBadChannelsIntArray(0),
  fBadChannelsIndexArray(0),
  fBadChannelsContainer(0)
{
  // Standard constructor for the SPD preprocessor
  //
  // Input : Name of digit file
  // Output: fInit
  // Return: Initialized SPD preprocessor object

  // Initialize
  AliITSPreprocessorSPD::Init();
  AliITSPreprocessorSPD::SetMaximumNumberOfEvents(maxNumberOfEvents);

  // Open and read the galice and digit files
  if (!AliITSPreprocessorSPD::Open(fileName, mode, fileNameg))
    {
      AliError("Failed to open file");
      fInit = kFALSE;
    }
}


//__________________________________________________________________________
AliITSPreprocessorSPD::AliITSPreprocessorSPD(const AliITSPreprocessorSPD &prep) :
TTask(prep),
  fITSLoader(prep.fITSLoader),
  fRunLoader(prep.fRunLoader),
  fThresholdRatio(prep.fThresholdRatio),
  fThreshold(prep.fThreshold),
  fMaximumNumberOfEvents(prep.fMaximumNumberOfEvents),
  fNumberOfModules(prep.fNumberOfModules),
  fHighestModuleNumber(prep.fHighestModuleNumber),
  fNumberOfColumns(prep.fNumberOfColumns),
  fNumberOfRows(prep.fNumberOfRows),
  fSelectedAlgorithm(prep.fSelectedAlgorithm),
  fGeometryMode(prep.fGeometryMode),
  fNumberOfBadChannels(prep.fNumberOfBadChannels),
  fIndex(prep.fIndex),
  fInit(prep.fInit),
  fVMEMode(prep.fVMEMode),
  fDigitsHistogram(prep.fDigitsHistogram),
  fBadChannelsObjArray(prep.fBadChannelsObjArray),
  fBadChannelsIntArray(prep.fBadChannelsIntArray),
  fBadChannelsIndexArray(prep.fBadChannelsIndexArray),
  fBadChannelsContainer(prep.fBadChannelsContainer)
{
  // Default copy constructor
  // Notice that only pointer addresses are copied!
  // Memory allocations of new objects are not done.


}


//__________________________________________________________________________
AliITSPreprocessorSPD& AliITSPreprocessorSPD::operator=(const AliITSPreprocessorSPD &prep)
{
  // Assignment operator

  AliError("Not implemented");

  if (this != &prep) { } // Do not delete this line

  return *this;
}


//__________________________________________________________________________
AliITSPreprocessorSPD::~AliITSPreprocessorSPD(void)
{
  // Destructor for the SPD preprocessor

  if (fDigitsHistogram)
    {
      // try fDigitsHistogram->Delete(); if the following crashes
      for (UInt_t module = 0; module < fNumberOfModules; module++)
	{
	  (*fDigitsHistogram)[module]->Delete();
	}
      delete fDigitsHistogram;
      fDigitsHistogram = 0;
    }

  if (fNumberOfBadChannels)
    {
      delete [] fNumberOfBadChannels;
      fNumberOfBadChannels = 0;
    }

  if (fBadChannelsIntArray)
    {
      delete [] fBadChannelsIntArray;
      fBadChannelsIntArray = 0;
    }

  if (fBadChannelsIndexArray)
    {
      delete [] fBadChannelsIndexArray;
      fBadChannelsIndexArray = 0;
    }

  if (fBadChannelsObjArray)
    {
      fBadChannelsObjArray->Delete();
      delete fBadChannelsObjArray;
      fBadChannelsObjArray = 0;
    }

  delete fRunLoader;
  fRunLoader = 0;

  delete fBadChannelsContainer;
  fBadChannelsContainer = 0;
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::Init(void)
{
  // Initialization of the SPD preprocessor
  //
  // Input : (void)
  // Output: Several logistical variables
  // Return: (void)

  // Default maximum number of events per histogram
  fMaximumNumberOfEvents = 1000000;

  // Default noisy channel removal algorithm
  fSelectedAlgorithm = kOptimizedForRealData;

  // Default noisy channel threshold and threshold ratio
  // (threshold for current bin content divided by the average neighboring bin contents)
  fThreshold = 10;
  fThresholdRatio = 5.;

  // Set default geometry mode (ALICE geometry). This also sets fNumberOfModules
  AliITSPreprocessorSPD::SetGeometryMode(kALICEGeometry);

  // The highest module number with found bad channels
  fHighestModuleNumber = 0;

  // Assume input is not from a VME file
  fVMEMode = kFALSE;

  // Initialization is complete
  fInit = kTRUE;
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::SetGeometryMode(UInt_t mode)
{
  // Set the geometry mode
  //
  // Input : Geometry mode (either kTestBeamGeometry or kALICEGeometry)
  // Output: 
  // Return: (void)

  fGeometryMode = mode;

  // In case of an input VME file, the number of modules has already been fixed.
  // Do not try to change it
  if (!fVMEMode)
    {
      if (mode == kALICEGeometry)
	{
	  fNumberOfModules = kNumberOfSPDModules;
	}
      else if (mode == kTestBeamGeometry)
	{
	  fNumberOfModules = kNumberOfTestBeamSPDModules;
	}
      else
	{
	  AliError("Unknown geometry mode, defaults to ALICE geometry");
	  fNumberOfModules = kNumberOfSPDModules;
	}
    }
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::SetFakeNoisyChannel(Int_t module, Int_t column, Int_t row)
{
  // Introduce a fake noisy channel in the hit histograms
  //
  // Input : Module, column and row numbers
  // Output: Updated hit histograms
  // Return: (void)

  if ((UInt_t)module < fNumberOfModules)
    {
      ((TH2F*)(*fDigitsHistogram)[module])->Fill(column, row, 1000000);
    }
  else
    {
      AliError("No such module number");
    }
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::Open(const char *fileName, const char *mode, const char *fileNameg)
{
  // Open digit file
  //
  // Input : Name and mode of ITS digit file, name of galice file
  // Output: Digits
  // Return: kTRUE if loaders succeed

  Bool_t status = kFALSE;
  Bool_t status0 = kFALSE;
  Bool_t status1 = kFALSE;
  Bool_t status2 = kFALSE;
  Bool_t status3 = kFALSE;

  // Only proceed if initialization has been done
  if (fInit)
    {
      TString m = (TString)mode;
      if (m == "daq" || m == "DAQ")
        {
          // Open the data file and get the run loader
          fRunLoader = AliRunLoader::Open(fileNameg);
          if (fRunLoader)
            {
              // Get the gAlice object
              status0 = AliITSPreprocessorSPD::GetgAlice();

              // Get the ITS digits
              if (status0) status1 = AliITSPreprocessorSPD::GetITSDigits(fileName);

              // Create the test beam object
              if (status1) status2 = AliITSPreprocessorSPD::CreateGeometryObj();

	      // Fill histograms with DAQ digits
	      if (status2) status3 = AliITSPreprocessorSPD::FillHistograms();

              status = status0 & status1 & status2 & status3;
            }
          else
            {
              AliError("Failed to get the run loader");
            }
        }
      else if (m == "vme" || m == "VME")
        {
	  // Open the VME file that contains histograms with fired channels as read by the stand-alone VME system
	  TFile *vmeFile = TFile::Open(fileName);

	  if (!vmeFile->IsOpen())
	    {
	      AliError("Could not open VME input file");
	    }
	  else
	    {
	      // Get the histograms from the VME file that contains all fired channels
	      status0 = AliITSPreprocessorSPD::GetVMEHistograms(vmeFile);

	      // Create the test beam object
	      if (status0) status1 = AliITSPreprocessorSPD::CreateGeometryObj();
	    }

	  // This boolean will be used to override any attempts of changing the number of modules by the user
	  // with the SetGeometryMode method. For VME files the number of modules will entirely be determined by
	  // the input VME root file, i.e. by the number of histograms in this file
	  fVMEMode = kTRUE;

	  status = status0 & status1;
        }
      else
        {
          AliError("Unknown filetype - assuming DAQ file");
        }

      // At this stage, the final number of modules will be known (safe to define arrays)
      // In case data is read from a VME root file, it will not be known before
      if (status)
	{
	  // Create the arrays for bookkeeping and storing the noisy channels
	  if (!fBadChannelsObjArray)
	    {
	      fBadChannelsObjArray = new TObjArray();
	    }
	  if (!fBadChannelsIndexArray)
	    {
	      // This array will contain the begin and end indices for each module, i.e. where to begin
	      // and stop reading the fBadChannelsObjArray for a certain module.
	      // The last position of the array will contain the end index of the last module
	      fBadChannelsIndexArray = new Int_t[fNumberOfModules + 1];
	      for (UInt_t module = 0; module < fNumberOfModules + 1; module++)
		{
		  fBadChannelsIndexArray[module] = 0;
		}
	    }
	}
    }
  else
    {
      AliError("SPD preprocessor not initialized. Can't load digits");
    }

  return status;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::GetgAlice(void)
{
  // Get the gAlice object
  //
  // Input : (void)
  // Output: gAlice object
  // Return: kTRUE if the gAlice object was found

  Bool_t status = kTRUE;

  // Load gAlice
  fRunLoader->LoadgAlice();
  if (!fRunLoader->GetAliRun())
    {
      AliError("gAlice not found on file. Aborting.");
      status = kFALSE;
    }

  return status;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::GetVMEHistograms(TFile *vmeFile)
{
  // Get pre-filled digit histograms from input VME file
  //
  // Input : pointer to VME file
  // Output: 
  // Return: kTRUE if histograms are found on file

  Bool_t status = kFALSE;

  // Get the file directory
  TDirectory *dir = (TDirectory *)vmeFile;

  // Get the number of keys (histograms in this case corresponding to modules/ladders)
  fNumberOfModules = dir->GetNkeys();
  if ((fNumberOfModules > 0) && (fNumberOfModules <= kNumberOfSPDModules))
    {
      status = kTRUE;

      // Create bad channel histograms
      fDigitsHistogram = new TObjArray(fNumberOfModules);

      // Create a key iterator
      TIter nextkey(dir->GetListOfKeys());
      TKey *key = 0;

      // Loop over all objects, read them in to memory one by one
      UInt_t module = 0;
      while ( (key = (TKey *)nextkey()) )
	{
	  (*fDigitsHistogram)[module++] = (TH2F *)key->ReadObj();

	  // For safety
	  if (module > fNumberOfModules)
	    {
	      AliError("Wrong number of keys in VME file");
	      status = kFALSE;
	      break;
	    }
	}
    }
  else
    {
      AliError("Wrong number of histograms in VME file");
    }
  return status;
}

//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::GetITSDigits(const char *fileName)
{
  // Get the ITS digits
  //
  // Input : name of digits file
  // Output: fITSLoader object, ITS digits
  // Return: kTRUE if loader succeed

  Bool_t status = kTRUE;

  // Load the gAlice and the header
  fRunLoader->LoadgAlice();
  fRunLoader->LoadHeader();

  // Get the ITS loader
  fITSLoader = (AliITSLoader*) fRunLoader->GetLoader("ITSLoader");
  if (!fITSLoader)
    {
      AliError("ITS loader not found");
      status = kFALSE;
    }
  else
    {
      // Open the digits file
      fITSLoader->SetDigitsFileName(fileName);
    }

  return status;
}


//__________________________________________________________________________
TClonesArray* AliITSPreprocessorSPD::CreateDigitsArray(void) const
{
  // Creation of the digits array
  //
  // Input : (void)
  // Output:
  // Return: Pointer to the SPD digits array

  // Create an array for 5 chips with 8192 channels each
  TClonesArray *digitsArray = new TClonesArray("AliITSdigitSPD", kNumberOfChannels);

  return digitsArray;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::CreateGeometryObj(void)
{
  // Creation of the geometry object
  //
  // This object is used to get the number of SPD half-staves
  //
  // Input : (void)
  // Output: fGeometryObj
  // Return: kTRUE if fGeometryObj has been created

  Bool_t status = true;

  // Create geometry object
  // fGeometryObj = new ...
  // if (!fGeometryObj) status = kFALSE;

  // Get SPD parameters
  // fNumberOfColumns = fGeometryObject->GetNumberOfColumns();
  // fNumberOfRows = fGeometryObject->GetNumberOfRows();

  fNumberOfColumns = kNumberOfColumns;
  fNumberOfRows = kNumberOfRows;

  return status;
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::CreateHistograms(void)
{
  // Create the noisy channel histograms
  //
  // Input : (void)
  // Output: Noisy channel histograms
  // Return: (void)

  TString moduleName;
  char n[4]; // For storing the module number (maximum module number is 240)

  // Create noisy channel histograms
  fDigitsHistogram = new TObjArray(fNumberOfModules);

  for (UInt_t i = 0; i < fNumberOfModules; i++)
    {
      // Histogram name
      moduleName = "module_";
      sprintf(n,"%d",i);
      moduleName.Append(n);

      (*fDigitsHistogram)[i] = new TH2F(moduleName,"Digits",
					fNumberOfColumns,-.5,(1.*fNumberOfColumns - .5),
					fNumberOfRows,-.5,(1.*fNumberOfRows - .5));
    }
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FillHistograms(void)
{
  // Fill the histograms with digits (hit positions of unclustered hits)
  //
  // (There is one digit histogram per SPD module, i.e. a half-stave, 10 chips)
  //
  // Input : No arguments (Empty digit histograms)
  // Output: Filled digit histograms
  // Return: kTRUE if histograms are filled with digits, kFALSE otherwise

  AliInfo("Filling noisy channel histograms");

  Bool_t status = kTRUE;
  AliITSdigitSPD *digitSPD = 0;
  UInt_t row;
  UInt_t column;
  TBranch *digitsSPDBranch;
  TTree *digitsTree;

  // Get the digits
  fITSLoader->LoadDigits("read");

  // Create noisy channel histograms
  AliITSPreprocessorSPD::CreateHistograms();

  // Create an empty SPD digits array
  TClonesArray *digitsArraySPD = AliITSPreprocessorSPD::CreateDigitsArray();

  // Get number of events
  UInt_t numberOfEvents = (fRunLoader->TreeE()) ? static_cast<UInt_t>(fRunLoader->TreeE()->GetEntries()) : 0;

  // Make sure we don't try to analyze more data than there actually is
  if (numberOfEvents < fMaximumNumberOfEvents)
    {
      fMaximumNumberOfEvents = numberOfEvents;
    }

  // Loop over all digits and put them in the corresponding histograms
  if (numberOfEvents > 0)
    {
      for (UInt_t event = 0; event < fMaximumNumberOfEvents; event++)
	{
	  // Get the current event
	  fRunLoader->GetEvent(event);

	  // Get the ITS digits tree
	  digitsTree = fITSLoader->TreeD();

	  // Disable all branches except the SPD branch to speed up the reading process
	  digitsTree->SetBranchStatus("ITSDigitsSPD",1);
	  digitsTree->SetBranchStatus("ITSDigitsSDD",0);
	  digitsTree->SetBranchStatus("ITSDigitsSSD",0);

	  // Reset the SPD digits array to make sure it is empty
	  digitsArraySPD->Clear();

	  // Get the SPD digits branch and set the address
	  digitsSPDBranch = digitsTree->GetBranch("ITSDigitsSPD");

	  digitsSPDBranch->SetAddress(&digitsArraySPD);

	  if (event%100 == 0) AliInfo(Form("Event #%d", event));

	  // Loop over all modules
	  UInt_t numberOfDigits = 0;
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      // Get event data for current module
	      digitsTree->GetEvent(module);

	      // Get the number of entries
	      numberOfDigits = digitsArraySPD->GetEntries();

	      // Loop over all digits
	      for (UInt_t digit = 0; digit < numberOfDigits; digit++)
		{
		  // Get the current digit
		  digitSPD = (AliITSdigitSPD*) digitsArraySPD->At(digit);
		  column = digitSPD->GetCoord1();
		  row = digitSPD->GetCoord2();

		  // Fill the digits histogram
		  ((TH2F*)(*fDigitsHistogram)[module])->Fill(column, row);
		  
		} // end digit loop
	    } // end module loop
	} // end event loop
    } // end numberOfEvents > 0

  else
    {
      status = kFALSE;
    }

  // Cleanup
  delete digitsArraySPD;
  digitsArraySPD = 0;

  return status;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FindDeadChannels(void)
{
  // Locate dead channels
  //
  // Input : Filled hit histograms
  // Output: TObjArray (fBadChannelsObjArray) with all identified dead channels
  // Return: kTRUE if dead channels have been found

  Bool_t status = kFALSE;

  // Proceed only if properly initialized
  if (fInit)
    {
      if (fVMEMode)
	{
	  // Initialize counters
	  fNumberOfBadChannels = new Int_t[fNumberOfModules];
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      fNumberOfBadChannels[module] = 0;
	    }

	  // Loop over all modules (intentional modularization - DCS will occationally want to
	  // look for noisy channels in certain modules, but not all)
	  fIndex = 0; // Global bad channels array counter (must be reset here)
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      status |= AliITSPreprocessorSPD::FindDeadChannelsInModule(module);
	    }
	}
      else
	{
	  AliError("Dead channels can only be found in data taken with stand-alone VME system");
	}
    }
  else
    {
      AliError("Not properly initialized");
    }

  return status;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FindDeadChannelsInModule(UInt_t module)
{
  // Locate dead channels
  // This method assumes that the preprocessor is operator in VME mode.
  // The algorithm is very simple. It assumes that the data was taken in
  // a mode where all working SPD pixels should respond as being hit.
  // fThreshold is used as the limit where everything below this value will
  // be considered as "dead".
  //
  // Input : Filled hit histograms
  // Output: TObjArray (fBadChannelsObjArray) with all identified noisy channels
  // Return: kTRUE if dead channels have been found

  // Store the index number for this module
  fBadChannelsIndexArray[module] = fIndex++;

  UInt_t xBin, numberOfXBins;
  UInt_t yBin, numberOfYBins;
  Float_t binContent;

  numberOfXBins = ((TH2F*)(*fDigitsHistogram)[module])->GetNbinsX();
  numberOfYBins = ((TH2F*)(*fDigitsHistogram)[module])->GetNbinsY();

  // Loop over all bins in this histogram
  for (xBin = 1; xBin <= numberOfXBins; xBin++)
    for (yBin = 1; yBin <= numberOfYBins; yBin++)
      {
	binContent = ((TH2F*)(*fDigitsHistogram)[module])->GetBinContent(xBin, yBin);

	// Store this channel/bin if outside accepted limits
	// A channel has to fire MORE times than the fThreshold value, or it will
	// be considered as "dead"
	if (binContent <= fThreshold)
	  {
	    // Store the dead channel in the array
	    // The channel object will be deleted in the destructor using the TObjArray Delete() method
	    // (The array will assume ownership of the channel object)
	    AliITSChannelSPD *channel = new AliITSChannelSPD((Int_t)(xBin - 1), (Int_t)(yBin - 1));

	    // Store the noisy channel in the array
	    fBadChannelsObjArray->Add(channel);

	    // Keep track of the number of bad channels in this module
	    fNumberOfBadChannels[module]++;
	    fIndex += 2;

	    // Keep track of the highest module number
	    if (module > fHighestModuleNumber) fHighestModuleNumber = module;

	    //AliInfo(Form("New dead pixel in (m,c,r) = (%d,%d,%d)", module, xBin - 1, yBin - 1));
	  }
      } // end bin loop

  return (fNumberOfBadChannels[module] > 0);
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FindNoisyChannels(void)
{
  // Locate noisy channels by searching through the digit histograms
  // (There is one digit histogram per SPD module, i.e. a half-stave, 5 chips)
  //
  // Input : Digits
  // Output: TObjArray (fBadChannelsObjArray) with all identified noisy channels
  // Return: kTRUE if noisy channels have been found

  Bool_t status = kFALSE;

  // Proceed only if properly initialized
  if (fInit)
    {
      // (For testing purposes, noisy channels can be inserted here)
      //SetFakeNoisyChannel(4,10,10);
      // Initialize counters
      fNumberOfBadChannels = new Int_t[fNumberOfModules];
      for (UInt_t module = 0; module < fNumberOfModules; module++)
	{
	  fNumberOfBadChannels[module] = 0;
	}

      // Scan through all the histogram bins and search for average filling deviations
      if ((fSelectedAlgorithm == kOptimizedForRealData) || (fSelectedAlgorithm == kOptimizedForRealDataRMS))
	{
	  // Noisy channel algorithm optimized for real data..........................
	  // Histograms can have any shape (both thresholds and quotients are used)
	  // This algorithm can be used to find noisy channels even if the data was
	  // taken with beam. All channels outside accepted limits (set by fThreshold
	  // and fThresholdRatio) will be marked as noisy

	  if (fSelectedAlgorithm == kOptimizedForRealData)
	    {
	      AliInfo("Searching for noisy channels (optimized for real data)");
	    }
	  else
	    {
	      AliInfo("Searching for noisy channels (optimized for real data, RMS version)");
	    }

	  // Loop over all modules (intentional modularization - DCS will occationally want to
	  // look for noisy channels in certain modules, but not all)
	  fIndex = 0; // Global bad channels array counter (must be reset here)
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      status |= AliITSPreprocessorSPD::FindNoisyChannelsInModuleAlgo0(module);
	    }
	} // end algorithm 0
      else
	{
	  // Noisy channel algorithm optimized for calibration data...........................
	  // Histograms will/should only contain noisy channels (only thresholds are used)
	  // Calibration histograms should have no background. The calibration data
	  // should have been taken without beam. All channels outside accepted limits
	  // (set by fThreshold) will be marked as noisy

	  AliInfo("Searching for noisy channels (optimized for calibration data)");

	  // Loop over all modules (intentional modularization - DCS will occationally want to
	  // look for noisy channels in certain modules, but not all)
	  fIndex = 0; // Global bad channels array counter (must be reset here)
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      status |= AliITSPreprocessorSPD::FindNoisyChannelsInModuleAlgo1(module);
	    }
	} // end algorithm 1
    } // end if fInit
  else
    {
      AliError("Not properly initialized");
    }

  return status;
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FindNoisyChannelsInModuleAlgo0(UInt_t module)
{
  // Locate the noisy channels in a module (optimized for real data)
  //
  // For each channel in these histograms, the algorithm checks the average
  // filling of the neighboring 3, 5 or 8 channels (depending on the location
  // of the channel in question; corner, border or inside), or compares with the
  // RMS of the neighbors. If the average is "much" lower, the channel will be
  // considered to be noisy. The default noisy-to-normal fraction is stored in the
  // fThresholdRatio varible. It can be set with the SetThresholdRatio method.
  // The channel also has to have fired more than a certain minimum, fThreshold.
  // It can be set with the SetThreshold method.
  //
  // To avoid difficulties with noisy channels that occur in pairs, the
  // neighboring channel with largest number of fillings will be removed from
  // the average calculation.
  //
  // NOTE: Since this method modifies the fBadChannelsObjArray and fBadChannelsIndexArray
  //       it is essential to initialize the fIndex counter before calling this module
  //       the first time. The bad channel data does not have to be ordered per module
  //       in the fBadChannelsObjArray, but the indices of where the data of a certain module
  //       starts has to be correct. A wrong fIndex can lead to segmentation violation
  //
  // Input : module number, filled digit histograms
  // Output: TObjArray (fBadChannelsObjArray) with all identified noisy channels,
  //         Int_t[] (fBadChannelsIndexArray) with fBadChannelsObjArray module indices,
  //         number of noisy channels in this module (global variable fNumberOfBadChannels[module])
  // Return: kTRUE if there are noisy channels in this module

  // Store the index number for this module
  fBadChannelsIndexArray[module] = fIndex++;

  UInt_t xBin, numberOfXBins;
  UInt_t yBin, numberOfYBins;
  UInt_t neighborXBin;
  UInt_t neighborYBin;
  UInt_t numberOfNeighboringBins;
  Float_t binContent;
  Float_t sumBinContent;
  Float_t neighborBinContent;
  Float_t maxBinContent;
  Float_t averageBinContent;
  Float_t ratio;

  numberOfXBins = (UInt_t)((TH2F*)(*fDigitsHistogram)[module])->GetNbinsX();
  numberOfYBins = (UInt_t)((TH2F*)(*fDigitsHistogram)[module])->GetNbinsY();

  // Loop over all bins in this histogram
  for (xBin = 1; xBin <= numberOfXBins; xBin++)
    for (yBin = 1; yBin <= numberOfYBins; yBin++)
      {
	numberOfNeighboringBins = 0;
	averageBinContent = 0.;
	sumBinContent = 0.;
	binContent = ((TH2F*)(*fDigitsHistogram)[module])->GetBinContent(xBin, yBin);
	maxBinContent = 0.;

	// Calculate the average pixel level on the surrounding pixels
	for (neighborXBin = xBin - 1; neighborXBin <= xBin + 1; neighborXBin++)
	  for (neighborYBin = yBin - 1; neighborYBin <= yBin + 1; neighborYBin++)
	    {
	      if ( !((neighborXBin == xBin) && (neighborYBin == yBin)) )
		{
		  // Only update the number of neighboring bins when we are not on the border
		  if ((neighborXBin > 0) && (neighborXBin <= numberOfXBins+1) &&
		      (neighborYBin > 0) && (neighborYBin <= numberOfYBins+1))
		    {
		      neighborBinContent =
			((TH2F*)(*fDigitsHistogram)[module])->GetBinContent(neighborXBin, neighborYBin);

		      if (fSelectedAlgorithm == kOptimizedForRealDataRMS)
			{
			  // RMS
			  sumBinContent += neighborBinContent*neighborBinContent;
			}
		      else
			{
			  // Geometrical mean
			  sumBinContent += neighborBinContent;
			}

		      if (neighborBinContent > maxBinContent) maxBinContent = neighborBinContent;

		      numberOfNeighboringBins++;
		    }
		}
	    }
	
	// Calculate the average. Remove the largest neighboring bin
	// (Correction for potential clusters of noisy channels)
	if (fSelectedAlgorithm == kOptimizedForRealDataRMS)
	  {
	    // Square the max bin content before removing it from the average calculation
	    maxBinContent *= maxBinContent;
	    // RMS
	    averageBinContent = TMath::Sqrt((sumBinContent - maxBinContent)/(Float_t)(numberOfNeighboringBins - 1));
	  }
	else
	  {
	    // Geometrical mean
	    averageBinContent = (sumBinContent - maxBinContent)/(Float_t)(numberOfNeighboringBins - 1);
	  }

	// Store this channel/bin if outside accepted limits
	// The threshold ratio is the threshold for the current bin content divided by the
	// average neighboring bin contents. The threshold bin content is the minimum number of
	// times a channel has to have fired to be called noisy
	ratio = (averageBinContent > 0) ? binContent/averageBinContent : 0.;
	if ( ((ratio >= fThresholdRatio) || (ratio == 0.)) && (binContent >= fThreshold) )
	  {
	    // Store the noisy channel in the array
	    // The channel object will be deleted in the destructor using the TObjArray Delete() method
	    // (The array will assume ownership of the channel object)
	    AliITSChannelSPD *channel = new AliITSChannelSPD((Int_t)(xBin - 1), (Int_t)(yBin - 1));

	    // Store the noisy channel in the array
	    fBadChannelsObjArray->Add(channel);

	    // Keep track of the number of bad channels in this module
	    fNumberOfBadChannels[module]++;
	    fIndex += 2;

	    // Keep track of the highest module number
	    if (module > fHighestModuleNumber) fHighestModuleNumber = module;

	    AliInfo(Form("New noisy pixel in (m,c,r) = (%d,%d,%d)", module, xBin - 1, yBin - 1));
	    AliInfo(Form("- Noisy pixel fired %d times, average neighborhood: %f",(Int_t)binContent,averageBinContent));
	  }
      } // end bin loop

  return (fNumberOfBadChannels[module] > 0);
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::FindNoisyChannelsInModuleAlgo1(UInt_t module)
{
  // Locate the noisy channels in a module (optimized for calibration data)
  //
  // This algorithm locates noisy channels by assuming original data was taken
  // in calibration mode. This should be done without beam and will thus only
  // contain data corresponding to background and noisy channels. The latter
  // should be clearly visible in this data so this algorithm simply assumes
  // that all histogram bins that are filled more than fThreshold times are
  // noisy.
  //
  // NOTE: Since this method modifies the fBadChannelsObjArray and fBadChannelsIndexArray
  //       it is essential to initialize the fIndex counter before calling this module
  //       the first time. The bad channel data does not have to be ordered per module
  //       in the fBadChannelsObjArray, but the indices of where the data of a certain module
  //       starts has to be correct. A wrong fIndex can lead to segmentation violation
  //
  // Input : module number, filled digit histograms
  // Output: TObjArray (fBadChannelsObjArray) with all identified noisy channels,
  //         Int_t[] (fBadChannelsIndexArray) with fBadChannelsObjArray module indices,
  //         number of noisy channels in this module (global variable fNumberOfBadChannels[module])
  // Return: kTRUE if there are noisy channels in this module

  // Store the index number for this module
  fBadChannelsIndexArray[module] = fIndex++;

  UInt_t xBin, numberOfXBins;
  UInt_t yBin, numberOfYBins;
  Float_t binContent;

  numberOfXBins = ((TH2F*)(*fDigitsHistogram)[module])->GetNbinsX();
  numberOfYBins = ((TH2F*)(*fDigitsHistogram)[module])->GetNbinsY();

  // Loop over all bins in this histogram
  for (xBin = 1; xBin <= numberOfXBins; xBin++)
    for (yBin = 1; yBin <= numberOfYBins; yBin++)
      {
	binContent = ((TH2F*)(*fDigitsHistogram)[module])->GetBinContent(xBin, yBin);

	// Store this channel/bin if outside accepted limits
	// The threshold bin content is the minimum number of times a channel has to have
	// fired to be called noisy
	if (binContent >= fThreshold)
	  {
	    // Store the noisy channel in the array
	    // The channel object will be deleted in the destructor using the TObjArray Delete() method
	    // (The array will assume ownership of the channel object)
	    AliITSChannelSPD *channel = new AliITSChannelSPD((Int_t)(xBin - 1), (Int_t)(yBin - 1));

	    // Store the noisy channel in the array
	    fBadChannelsObjArray->Add(channel);

	    // Keep track of the number of bad channels in this module
	    fNumberOfBadChannels[module]++;
	    fIndex += 2;

	    // Keep track of the highest module number
	    if (module > fHighestModuleNumber) fHighestModuleNumber = module;

	    AliInfo(Form("New noisy pixel in (m,c,r) = (%d,%d,%d)", module, xBin - 1, yBin - 1));
	    AliInfo(Form("- Noisy pixel fired %d times",(Int_t)binContent));
	  }
      } // end bin loop

  return (fNumberOfBadChannels[module] > 0);
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::PrintChannels(void)
{
  // Print all found bad channels to stdout
  //
  // Input : fBadChannelsObjArray
  // Output: (dump to stdout)
  // Return: (void)

  Int_t i = 0;
  Int_t j = 0;
  AliITSChannelSPD *channel = 0;

  // Print the bad channels stores in the array
  AliInfo("\nModule #\tColumn #\tRow #\n------------------------------------------------");
  for (UInt_t module = 0; module < fNumberOfModules; module++)
    {
      j = 0;
      while (j < fNumberOfBadChannels[module])
	{
	  channel = (AliITSChannelSPD *) fBadChannelsObjArray->At(i++);
	  std::cout << module << "\t\t" << channel->GetColumn() << "\t\t" << channel->GetRow() << std::endl;

	  // Go to next bad channel
	  j++;
	}
    }

  AliInfo(Form("%d bad channels were found", fBadChannelsObjArray->GetEntries()));
}


//__________________________________________________________________________
void AliITSPreprocessorSPD::MarkNoisyChannels(void)
{
  // WARNING: THIS METHOD DOESN'T WORK!!!
  //
  // Mark all identified noisy channels
  //
  // Input : List of noisy channels, original digits tree
  // Output: New digits tree containing SPD digits marked when noisy
  // Return: (void)
  //
  // The original digits tree (digitsTree) is cloned except for the SPD branch (ITSDigitSPD).
  // This branch is then redefined for each event and will contain all the original
  // information. All known noisy channels will be marked by using the TObject status bits
  // according to the following scheme. Dead channels are included for completeness. Note
  // that a dead channel will NEVER show up among digits..
  //
  // Interpretation of digit status bits (LSB):
  // Dead channel   Noisy channel   |  Integer
  // -----------------------------------------
  //      0               0         |     0
  //      0               1         |     1
  //      1               0         |     2
  //
  // meaning e.g. that a channel that is noisy will have the first bit set in its status bits

  // Do not continue unless we are processing DAQ data
  if (!fVMEMode)
    {
      AliInfo("Marking bad channels");

      // Create the storage container that will be used to access the bad channels
      if (!fBadChannelsContainer)
	{
	  // Add the bad channels array to the storage container
	  // (ownership is passed to the AliRunDataStorage object)
	  fBadChannelsContainer = new AliITSBadChannelsSPD();

	  // Convert the bad channels from TObjArray to Int_t[]
	  AliITSPreprocessorSPD::ConvertObjToIntArray();

	  // Store the arrays in the bad channels container object
	  const Int_t kBadChannelsArraySize =
	    2*fBadChannelsObjArray->GetEntries() + fNumberOfModules;
	  fBadChannelsContainer->Put(fBadChannelsIntArray, kBadChannelsArraySize,
				     fBadChannelsIndexArray, fNumberOfModules);
	}

      // Create the bad channels helper object
      // (will be used to find a bad channel within a TObjArray)
      AliITSBadChannelsAuxSPD *aux = new AliITSBadChannelsAuxSPD();

      AliITSdigitSPD *digitSPD = 0;
      UInt_t numberOfDigits;
      Int_t newDigit[3];
      Bool_t mark = kFALSE;

      TBranch *digitsBranch = 0;
      TTree *digitsTree;

      // Create an empty SPD digit array
      TObjArray *digitsArraySPD = new TObjArray();

      // Get the digits in update mode (we want to modify them if there are noisy channels)
      fITSLoader->UnloadDigits();
      fITSLoader->LoadDigits("update");

      // Get the number of events
      UInt_t numberOfEvents = (fRunLoader->TreeE()) ? static_cast<UInt_t>(fRunLoader->TreeE()->GetEntries()) : 0;

      // Loop over all events
      for (UInt_t event = 0; event < numberOfEvents; event++)
	{
	  if (event%100 == 0) AliInfo(Form("Event #%d", event));

	  // Get the current event
	  fRunLoader->GetEvent(event);

	  // Get the ITS digits tree
	  digitsTree = fITSLoader->TreeD();

	  // Get SPD branch that will contain all digits with marked noisy channels
	  digitsBranch = digitsTree->GetBranch("ITSDigitsSPD");
	  digitsBranch->SetAddress(&digitsArraySPD);

	  // Get the stored number of modules
	  UInt_t numberOfModules = (Int_t)digitsTree->GetEntries();
	  TObjArray **newDigitsArraySPD = new TObjArray*[numberOfModules];

	  Int_t *digitNumber = new Int_t[numberOfModules];
	  for (UInt_t m = 0; m < numberOfModules; m++)
	    {
	      newDigitsArraySPD[m] = new TObjArray();
	      digitNumber[m] = 0;
	    }

	  AliInfo(Form("ent = %d", (Int_t)digitsTree->GetEntries()));

	  // Reset the SPD digit arrays to make sure they are empty
	  digitsArraySPD->Clear();

	  // Get the SPD digits branch from the original digits tree and set the address
	  digitsBranch = digitsTree->GetBranch("ITSDigitsSPD");
	  digitsBranch->SetAddress(&digitsArraySPD);

	  // Loop over all modules
	  for (UInt_t module = 0; module < fNumberOfModules; module++)
	    {
	      // Get event data for current module
	      digitsTree->GetEvent(module);

	      // Get the hits in the current module
	      TObjArray *moduleObjArray = fBadChannelsContainer->CreateModuleObjArray(module);

	      // Get the number of entries
	      numberOfDigits = digitsArraySPD->GetEntries();

	      // Loop over all digits and all channels
	      for (UInt_t digit = 0; digit < numberOfDigits; digit++)
		{
		  // Get the current digit
		  digitSPD = (AliITSdigitSPD*) digitsArraySPD->At(digit);
		  newDigit[0] = digitSPD->GetCoord1(); // column
		  newDigit[1] = digitSPD->GetCoord2(); // row
		  newDigit[2] = digitSPD->GetSignal(); // signal

		  // Check if this channel is noisy
		  // (Compare with all stored channels in the bad channels array)
		  if (aux->Find(digitSPD, moduleObjArray))
		    {
		      // Set the mark flag and break the loop
		      mark = kTRUE;
		    }

		  // Store this digit in the SPD digits array using a placement new operation
		  new ((*newDigitsArraySPD[module])[digitNumber[module]]) AliITSdigitSPD(newDigit);

		  // Mark it if noisy and store in the noisy channel array
		  if (mark)
		    {
		      // Store this digit in the marked SPD digits array using a placement new operation
		      //new ((*badChannels[m])[numberOfBadChannels[m]]) AliITSChannelSPD(newBadChannel);
		      //new ((*newDigitsArraySPD[module])[digitNumber[module]]) AliITSdigitSPD(newDigit);

		      // Mark the original channel as noisy
		      ((*newDigitsArraySPD[module])[digitNumber[module]])->SetBit(kNoisyChannel);

		      mark = kFALSE;
		    }

		  digitNumber[module]++;

		} // end digit loop

	      // Cleanup
	      delete moduleObjArray;
	      moduleObjArray = 0;

	    } // end module loop

	  digitsBranch->Reset();
	  digitsBranch->ResetAddress();

	  // Cleanup
	  delete digitsArraySPD;
	  digitsArraySPD = 0;
	  digitsTree->Reset();

	  // WHY THIS RANGE?????????????????????????????????????????????????????????????????????
	  for (UInt_t n = 0; n < event; n++)
	    {
	      digitsTree->SetBranchAddress("ITSDigitsSPD", &newDigitsArraySPD[n]);
	      digitsTree->Fill();
	    }

	  digitsTree->AutoSave();

	  // Cleanup
	  for (UInt_t n = 0; n < event; n++)
	    {
	      delete newDigitsArraySPD[n];
	    }
	  delete [] newDigitsArraySPD;
	  newDigitsArraySPD = 0;
	  delete [] digitNumber;
	  digitNumber = 0;
	  delete digitsTree;
	  digitsTree = 0;

	} // end loop over all events

      // Unload the digits
      fITSLoader->UnloadDigits();

      // Cleanup
      delete aux;
      aux = 0;
    }
}

/*

//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::Store(AliCDBId &id, AliCDBMetaData *md)
{
  // Store the bad channels object in the calibration database
  // (See the corresponding run macro for further explanations)
  //
  // Input : fBadChannelsObjArray (now containing all found bad channels), object meta data
  // Output: Database file containing the bad channels
  // Return: kTRUE if successful

  Bool_t status = kFALSE;

  AliInfo("Storing bad channels");

  // Add the bad channels array to the storage container
  // (ownership is passed to the AliRunDataStorage object)
  fBadChannelsContainer = new AliITSBadChannelsSPD();

  // Convert the bad channels from TObjArray to Int_t[]
  AliITSPreprocessorSPD::ConvertObjToIntArray();

  // Store the arrays in the bad channels container object
  const Int_t kBadChannelsArraySize =
    2*fBadChannelsObjArray->GetEntries() + fNumberOfModules;
  fBadChannelsContainer->Put(fBadChannelsIntArray, kBadChannelsArraySize,
			     fBadChannelsIndexArray, fNumberOfModules);

  // Store the container
 if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
   //AliError("No storage set!");
   //	return status;
   AliCDBManager::Instance()->SetDefaultStorage("local://Calib");
 }
  
  if (AliCDBManager::Instance()->GetDefaultStorage()->Put(fBadChannelsContainer, id, md))
    {
      AliInfo("Bad channels object stored in database");
      status = kTRUE;
    }
  else
    {
      AliError("Failed to store object in database");
    }

  return status;
}

*/
//__________________________________________________________________________
void AliITSPreprocessorSPD::ConvertObjToIntArray()
{
  // Convert the bad channel TObjArray to an Int_t array
  //
  // Input : fBadChannelsObjArray (now containing all found bad channels)
  // Output: fBadChannelsIntArray
  // Return: (void)
  //
  // Data encoding:
  //                The TObjArray of this class (fBadChannelsObjArray) is converted to a sequential
  //                Int_t array (fBadChannelsIntArray) in this method. For each module, the first
  //                stored number is the number of bad channels in the current module. This is
  //                followed by all the columns and rows of the bad channels:
  //
  //                badChannelsArray =
  //                | N(m) | col0 | row0 | .. | colN(m) | N(m+1) | col0 | row0 | ...
  //                .   .......... module m .........   .   .... module m+1 ......
  //
  //                The bad channels index array (fBadChannelsIndexArray) contains the indices of
  //                the badChannelsArray, i.e. where the bad channels in certain module starts:
  //
  //                fBadChannelsObjArray =
  //                | i0 | i1 | .. | iM | (where M = the number of SPD modules)
  //
  //                e.g. i1 corresponds to the index of the badChannelsArray where N(1) is stored,
  //                i.e. the number of bad channels for module 1

  const Int_t kBadChannelsArraySize =
    2*fBadChannelsObjArray->GetEntries() + fNumberOfModules;
  fBadChannelsIntArray = new Int_t[kBadChannelsArraySize]; // Will be deleted in dtor
  AliITSChannelSPD *channel = 0;
  Int_t i = 0;
  Int_t j = 0;
  Int_t k = 0;

  // Loop over all modules
  for (UInt_t module = 0; module < fNumberOfModules; module++)
    {
      // Encode the number of bad channels of the current module
      fBadChannelsIntArray[k++] = fNumberOfBadChannels[module];

      // The columns and rows of the fBadChannelsObjArray will be stored sequentially
      // in the Int_t array
      j = 0;
      while (j < fNumberOfBadChannels[module])
	{
	  channel = (AliITSChannelSPD *) fBadChannelsObjArray->At(i++);
	  fBadChannelsIntArray[k++] = channel->GetColumn();
	  fBadChannelsIntArray[k++] = channel->GetRow();

	  // Go to next bad channel
	  j++;
	}
    }
}


//__________________________________________________________________________
Bool_t AliITSPreprocessorSPD::Store(AliCDBId& /*id*/, AliCDBMetaData* /*md*/, Int_t runNumber)
{
  // Store the bad channels object in the calibration database
  // (See the corresponding run macro for further explanations)
  //
  // Input : fBadChannelsObjArray (now containing all found bad channels), object meta data
  // Output: Database file containing the bad channels
  // Return: kTRUE if successful

  Bool_t status = kFALSE;

  AliInfo("Storing bad channels");
 
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliWarning("No storage set! Will use dummy one");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  }

  
  AliCDBEntry *entrySPD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSPD", runNumber);
  if(!entrySPD){
    AliWarning("Calibration object retrieval failed! Dummy calibration will be used.");
    AliCDBStorage *origStorage = AliCDBManager::Instance()->GetDefaultStorage();
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
	
    entrySPD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSPD", runNumber);
    AliCDBManager::Instance()->SetDefaultStorage(origStorage);
  }

  TObjArray *respSPD = (TObjArray *)entrySPD->GetObject();

  if ((! respSPD)) {
    AliWarning("Can not get calibration from calibration database !");
    return kFALSE;
  }
  

  Int_t i=0;
  AliITSChannelSPD *channel = 0;
  AliITSCalibrationSPD* res;
  for (Int_t module=0; module<respSPD->GetEntries(); module++) {
    res = (AliITSCalibrationSPD*) respSPD->At(module);
    Int_t j = 0;
    while (j < fNumberOfBadChannels[module]) {
      channel = (AliITSChannelSPD *) fBadChannelsObjArray->At(i++);
      res->AddDead(channel->GetColumn(),channel->GetRow());
      // Go to next bad channel
      j++;
    }
  }
    
  
  AliCDBManager::Instance()->Put(entrySPD);
  entrySPD->SetObject(NULL);
  entrySPD->SetOwner(kTRUE);

  delete entrySPD;
  status=kTRUE;
  return status;
}
