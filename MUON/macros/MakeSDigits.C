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

/* $Id */

/// \ingroup macros
/// \file MakeSDigits.C
/// \brief Macro for generation SDigits from raw data for merging 
/// \author Indranil Das, HEP Division, SINP (indra.das@saha.ac.in, indra.ehep@gmail.com)
///
/// Usage : To run this code one should have $ALICE_ROOT/MUON/rootlogon.C 
///         in the current working directory. Then run it from command line as,
/// <pre>aliroot MakeSDigits.C\(\""reconstructed galicefile"\",\""rawdatafile"\",\""OCDB path"\",run-numer\)
///
///        where inputs are : 1. galice file of local reconstruction directory
///                           2. rootified rawdata file of local reconstruction directory.
///                           In case the raw data are ddl files, specify the path which has "raw0", "raw1".... events
///                           3. OCDB path
///                           4. run number
///                           5. switch to merge(true) or not merge(false) trigger digits 
///        and the output is : "MUON.SDigits.root" file with the same tree and event structure as produced in simulation directory
/// </pre>
/// Note:  
/// galice.root and raw.root cannot be a directory read from alien 
/// beacuse "muonLoader->WriteSDigits("OVERWRITE")" tries to write in reconstruction directory,
/// which will fail. For files in alien, make a local copy of the galice and raw.root file, and run the code
/// in that directory.
///
/// Merging Hints :    
/// Follow the merging procedure as specified in the \ref README_sim page, with 
/// <pre>MuonSim.MergeWith("recodir_that_has_MUON.SDigits.root_from_rawdata/galice.root",nofBackground)
/// </pre>
/// in the $ALICE_ROOT/MUON/runSimulation.C

/**********************************************************************
 Created on : 11/01/2010
 Purpose    : To create SDigits from RawData
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/


int MakeSDigits(const char* galiceFile="galice.root", const char* rawRootFile="./raw.root", 
		const char* ocdb = "local://$ALICE_ROOT/OCDB", int run=0, bool isMergeTrigger=true)
{
  //TGrid connect for alien ocdb
  TGrid::Connect("alien://");
  
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  man->SetRun(run);
  
  AliMpCDB::LoadDDLStore(true);
  
  AliRawReader *rawReader = AliRawReader::Create(rawRootFile);

  AliRunLoader* runLoader = AliRunLoader::Open(galiceFile);
  AliLoader* muonLoader = runLoader->GetDetectorLoader("MUON");
  
  for(Int_t iEvent=0;iEvent<runLoader->GetNumberOfEvents();iEvent++){
 
    cout<<"Running for Event : "<<iEvent<<endl;
   
    rawReader->NextEvent();
    runLoader->GetEvent(iEvent);
    
    muonLoader->LoadSDigits("update");
    muonLoader->CleanSDigits();
    if (!muonLoader->TreeS()) muonLoader->MakeSDigitsContainer();
    
    TTree* treeS = muonLoader->TreeS();
    
    AliMUONVDigitStore* sDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2S");
    sDigitStore->Connect(*treeS,true);
    AliMUONDigitMaker *digitMaker = new AliMUONDigitMaker(false);
    
    if(isMergeTrigger){
      AliMUONVTriggerStore* triggerStore = new AliMUONTriggerStoreV1;
      triggerStore->Connect(*treeS,true);
      digitMaker->SetMakeTriggerDigits(true);
      digitMaker->Raw2Digits(rawReader,sDigitStore,triggerStore);
    }else{
      digitMaker->Raw2Digits(rawReader,sDigitStore,0x0);
    }
    
    TIter next(sDigitStore->CreateIterator());
    AliMUONVDigit *mdigit;
    
    while ( (mdigit = reinterpret_cast<AliMUONVDigit *>(next())) ) {
      if(mdigit->DetElemId()<1100){
	mdigit->SetCharge(Float_t(mdigit->ADC()));
	mdigit->SetADC(0);
      }else{
	mdigit->SetCharge(1.0);
      }
      //mdigit->Print();
    }
    treeS->Fill();
   
    muonLoader->WriteSDigits("OVERWRITE");
    
    muonLoader->UnloadSDigits();
    
    if(isMergeTrigger)
      triggerStore->Delete();
    sDigitStore->Delete();
    digitMaker->Delete();
    
  }
  delete runLoader;

  return 0;
}