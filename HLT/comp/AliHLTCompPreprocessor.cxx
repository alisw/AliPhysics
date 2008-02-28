// $Id: AliHLTCompPreprocessor.cxx 23039 2007-12-13 20:53:02Z richterm $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Jenny Wagner                                          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/**
 * @file   AliHLTCompPreprocessor.cxx
 * @author Jenny Wagner, Matthias Richter
 * @brief  Implementation of the HLT preprocessor for the AliHLTComp library
 */

#include "AliHLTCompPreprocessor.h"
#include "AliCDBMetaData.h"
#include "AliPreprocessor.h"
#include "TObjString.h"
#include "TString.h"
#include "TList.h"
#include "TFile.h"

// necessary for huffman table to get information about the origin
#include "AliHLTCOMPHuffmanData.h"

ClassImp(AliHLTCompPreprocessor)

AliHLTCompPreprocessor::AliHLTCompPreprocessor() 
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTCompPreprocessor::fgkHuffmanFileId = "HuffmanData";

AliHLTCompPreprocessor::~AliHLTCompPreprocessor()
{
  // see header file for function documentation
}

void AliHLTCompPreprocessor::Initialize(Int_t /*run*/, UInt_t /*startTime*/, 
					UInt_t /*endTime*/)
{
  // see header file for function documentation
}


UInt_t AliHLTCompPreprocessor::Process(TMap* /*dcsAliasMap*/)
{
  // see header file for function documentation
  UInt_t retVal = 0;

  if (GetHuffmanTables() != 0) {
    // unable to fetch Huffman tables
    retVal = 1; 
    // but if set to 1, then also file from GetTempHisto won't be saved !!!
  }
	
  return retVal;
}

UInt_t AliHLTCompPreprocessor::GetHuffmanTables() 
{
  // see header file for function documentation

  UInt_t retVal = 0;

  // get all huffman tables stored at FXS	
  TList* HuffmanList = GetFileSources(AliPreprocessor::kHLT, fgkHuffmanFileId);
  // -> list of all DDL numbers that own a huffman table

  // if no huffman tables are produced, return 0 (successful end --> test if TPC, PHOS are used in current run!)
  // --> still open to be implemented!
  // if there is an error getting the sources return 1
  // else produce containers for each detector to be stored in the OCDB
  if (!HuffmanList) 
    {
      Log("No Huffman code tables for HLT");
      return 1;
    }
	  
  TList* TPCHuffmanList = new TList();
  TList* PHOSHuffmanList = new TList();
	
  // loop over all DDL numbers and put huffman tables into special containers
  // (one for each detector)
  for(Int_t ii=0; ii < HuffmanList->GetEntries(); ii++)
    {
      // get huffman table
      TObjString *objstr = (TObjString*) HuffmanList->At(ii);

      if(objstr == 0) // should not happen?! 
	{
	  char logging[1000];
	  sprintf(logging, "Error in Huffmanlist, no DDL at position %d", ii);
	  Log(logging);
	  continue;
	}

      TString fileName = GetFile(AliPreprocessor::kHLT, "HuffmanData", objstr->GetName());

      if (!(fileName.Length() > 0)) // error if local path/filename is not defined
	{
	  Log("Local file for current Huffman table is not properly defined.");
	  return 1;
	}

      TFile* currenthuffmanfile = new TFile(fileName, "READ");

      // if current huffman table file does not contain a table, return an error
      if ( currenthuffmanfile->Get("fHuffmanData") == NULL)
	{
	  char logging[1000];
	  sprintf(logging,"Local file %s does not contain a Huffman code table.", fileName.Data());
	  Log(logging);
	  //retVal = 1; // retVal must be zero to give other functions a chance to read their data
	  retVal = 0;
	}
	    
      TObject* huffmandata = (TObject*) currenthuffmanfile->Get("fHuffmanData");
      // class object not needed since container uses TObjects!
      //AliHLTCOMPHuffmanData* currenthuffmandata = (AliHLTCOMPHuffmanData*) huffmandata;

      // specifications necessary for sorting process
      TString detectororigin = "";//currenthuffmandata->GetOrigin();
      Int_t tablespec = 0;//currenthuffmandata->GetDataSpec();
	   
	    
      // plug them into a container:
      if(detectororigin == "PHOS") // belongs to PHOS table (one one!)
	{
	  PHOSHuffmanList->AddFirst(huffmandata);

	  if(PHOSHuffmanList->GetEntries() > 1)
	    {
	      Log("More than one table available for PHOS.");
	      // return: warning but go on...
	    }
	}
      else
	{
	  if(detectororigin == "TPC ") // belongs to TPC tables (six)
	    {
	      if(tablespec < 6)
		{
		  TPCHuffmanList->Add(huffmandata);

		  if(TPCHuffmanList->GetEntries() > 6)
		    {
		      Log("More than six tables available for TPC.");
		      // return warning but go on...
		    }
		}
	      else
		{
		  char logging[1000];
		  sprintf(logging, "Read data specification %d from Huffman table too large to belong to TPC.", tablespec);
		  Log(logging);
		  // retVal = 1; // retVal must be zero to give other functions a chance to read their data
		  retVal = 0;
		}
	    }
	  else // error! 
	    {
	      char logging[1000];
	      sprintf(logging, "Specified detector pattern %s does not define a valid detector.", detectororigin.Data());
	      Log(logging);
	      // retVal = 1; // retVal must be zero to give other functions a chance to read their data
	      retVal = 0;
	    }
	}

	    
    } // end loop over all DDLs

   
  // after loop all containers are filled and can be stored in OCDB
  AliCDBMetaData meta("Jenny Wagner");

  if (!(Store("CalibTPC", "HuffmanCodeTables", (TObject*) TPCHuffmanList, &meta, 0, kTRUE))) 
    {

      Log("Storing of TPCHuffmanList (Huffman code tables for TPC) to OCDB failed.");

      if (!(StoreReferenceData("CalibTPC", "HuffmanCodeTables", (TObject*) TPCHuffmanList, &meta)))
	{
	  Log("Storing of TPCHuffmanList (Huffman code tables for TPC) to reference storage failed.");

	  retVal = 1;
	}
    }
       
  if (!(Store("CalibPHOS", "HuffmanCodeTables", (TObject*) PHOSHuffmanList, &meta, 0, kTRUE))) 
    {

      Log("Storing of PHOSHuffmanList (Huffman code table for PHOS) to OCDB failed.");

      if (!(StoreReferenceData("CalibPHOS", "HuffmanCodeTables", (TObject*) PHOSHuffmanList, &meta)))
	{
	  Log("Storing of PHOSHuffmanList (Huffman code table for PHOS) to reference storage failed.");

	  retVal = 1;
	}
    }

  return retVal;
}
