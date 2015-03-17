#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TGrid.h"
#include "TLeaf.h"
#include "TError.h"
#include "TFile.h"
#include "Riostream.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "TMath.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawVEquipment.h"
#include "AliRawDataHeader.h"
#include "AliRawData.h"
#include "TH1.h"
#include <map>
#include <string>
#include "AliDAQ.h"

using std::cout;
using std::endl;
namespace AAF {
  
namespace RAWMUON {

   std::map<int,int> ddlMap;

   void DisableBranches(TTree* tree)
  {
    TObjArray* list = tree->GetListOfLeaves();
    TIter next(list);
    TLeaf* leaf;
    
    TObjArray branchesToKeep;
    branchesToKeep.SetOwner(kTRUE);
    
    branchesToKeep.Add(new TObjString("rawevent"));
    branchesToKeep.Add(new TObjString("MUON"));
    branchesToKeep.Add(new TObjString("TRG"));
    branchesToKeep.Add(new TObjString("T0"));
    branchesToKeep.Add(new TObjString("VZERO"));
    branchesToKeep.Add(new TObjString("ZDC"));
    branchesToKeep.Add(new TObjString("ITSSPD"));
    
    TIter nit(&branchesToKeep);
    
    tree->SetBranchStatus("*",0);
    
    Bool_t on(kTRUE);
    
    while ( ( leaf = static_cast<TLeaf*>(next()) ) )
    {
      TString name(leaf->GetName());
      if (!name.BeginsWith("f") && name != "rawevent" )
      {
        TObjString* str;
        
        nit.Reset();
        
        on = kFALSE;
        
        while ( ( str = static_cast<TObjString*>(nit()) ) )
        {
          if ( name.BeginsWith(str->String()) )
          {
            on = kTRUE;
          }
        }
      }
      if ( on )
      {
        leaf->GetBranch()->SetStatus(1);
      }
    }
  }
  
  void DumpCDH(AliRawDataHeader *cdh)
  {
    // Dumps the CDH
    // ...
    cout << "        Size: " << cdh->fSize << endl;
    cout << "        Version: " << (Int_t)cdh->GetVersion() << endl;
    cout << "        Orbit: " << cdh->GetEventID2() << " Bunch-crossing: " << cdh->GetEventID1() << endl;
    cout << "        L1 trigger message: " << (UInt_t)cdh->GetL1TriggerMessage() << endl;
    cout << "        Participating sub-detectors: " << cdh->GetSubDetectors() << endl;
    cout << "        Block attributes: " << (Int_t)cdh->GetAttributes() << endl;
    cout << "        Status: " << cdh->GetStatus() << endl;
    cout << "        Mini event ID: " << cdh->GetMiniEventID() << endl;
    cout << "        Trigger classes: " << cdh->GetTriggerClasses() << endl;
    cout << "        ROI: " << cdh->GetROI() << endl;
  }

  bool CheckEvent(AliRawVEvent* rawEvent, bool debug)
  {
	  AliRawEventHeaderBase *rawEventHeader = rawEvent->GetHeader();

	  if (rawEventHeader->GetMagic() != 0xDA1E5AFE) {
	    Error("CheckEvent","Wrong magic number ( 0x%x != 0xDA1E5AFE )",rawEventHeader->GetMagic());
	    return false;
	  }

	  AliRawDataHeader* cdhRef = 0x0;

	  if (debug)
	  {
		  std::cout << "  *********** Event header ***********" << std::endl;
		  rawEventHeader->Print();
	  }

	  for(Int_t iSubEvent=0; iSubEvent < rawEvent->GetNSubEvents(); iSubEvent++)
	  {
	     AliRawVEvent *rawSubEvent = rawEvent->GetSubEvent(iSubEvent);
	     AliRawEventHeaderBase *rawSubEventHeader = rawSubEvent->GetHeader();

	     if  (debug)
	     {
	    	 std::cout << "    *********** Sub-event header (" << iSubEvent+1 << " / "
	    			 << rawEvent->GetNSubEvents() << " ) ***********" << std::endl;
	    	 rawSubEventHeader->Print("  ");
	     }

	     for(Int_t iEquipment=0; iEquipment < rawSubEvent->GetNEquipments(); iEquipment++)
	     {
	        AliRawVEquipment *rawEquip = rawSubEvent->GetEquipment(iEquipment);
	        AliRawEquipmentHeader *rawEquipHeader = rawEquip->GetEquipmentHeader();

	        if ( debug )
	        {
	        	cout << "      *********** Equipment event header (" << iEquipment+1 << " / " << rawSubEvent->GetNEquipments() << " ) ***********" << endl;
	        	rawEquipHeader->Print("    ");
	        }
	        AliRawData *rawData = rawEquip->GetRawData();
	        if ( !rawData ) continue; // protection needed if we disabled branches

	        ddlMap[rawEquipHeader->GetId()]++;
	     }
	  }

	  return true;
  }

  int CheckFile(const char* from, int maxevents)
  {
    TString sinputfile(from);
    
    if (sinputfile.BeginsWith("alien://"))
    {
      if (!gGrid)
      {
        TGrid::Connect("alien://");
      }
      if (!gGrid)
      {
        Error("CheckFile","Cannot connect to Grid !");
        return -1;
      }
    }
    
    TFile* rawFile = TFile::Open(from);
    
    if (!rawFile || !rawFile->IsOpen())
    {
      Error("CheckFile","Cannot open input file %s",from);
      return -3;
    }
    
    TTree *rawTree=(TTree *)rawFile->Get("RAW");
    if(!rawTree)
    {
      Error("CheckFile","Error getting RAW tree from file %s",from);
      return -2;
    }
    
    RAWMUON::DisableBranches(rawTree);

    AliRawVEvent *rawEvent=NULL;
    
    rawTree->SetBranchAddress("rawevent", &rawEvent);
    
    Int_t nEvents = rawTree->GetEntries();
    
    std::cout << "*******************************************" << std::endl;
    std::cout << "File: " << from << std::endl;
    std::cout << "GUID: " << rawFile->GetUUID().AsString() << std::endl;
    std::cout << "Total number of events: " << nEvents << std::endl;
    std::cout << "*******************************************" << std::endl;
    
    bool debug(false);

    if ( maxevents > 0 )
    {
    	nEvents = TMath::Min(maxevents,nEvents);
    }

    for(Int_t iEvent=0; iEvent < nEvents; ++iEvent)
    {
      rawTree->GetEntry(iEvent);
      if (debug)
      {
    	  	  std::cout << "  *********** Event " << iEvent << " *******" << std::endl;
      }
    	CheckEvent(rawEvent,debug);
//      DumpEvent(argv[0],rawEvent);
      rawEvent->Clear();
    }
    
    delete rawEvent;
    
    if ( debug )
    {
    	std::cout << "*******************************************" << std::endl;
    	std::cout << "EOF" << std::endl;
    	std::cout << "*******************************************" << std::endl;
    }
    
    delete rawTree;
    
    rawFile->Close();

    std::map<int,int>::const_iterator it;
    int dummy;

    for ( it = ddlMap.begin(); it != ddlMap.end(); ++it )
    {
    	std::cout << Form("DDL %5d (%8s) seen in %7d events",it->first,
    			 AliDAQ::DetectorNameFromDdlID(it->first,dummy),it->second) << std::endl;
    }
    
    return 0;
  }
  
}

  int FILTER_RAWMUON(const char* from, const char* to)
  {
    std::cout << "FILTER_RAWMUON(" << from << "," << to << ")" << std::endl;
    
    TString sinputfile(from);
    
    if (sinputfile.BeginsWith("alien://"))
    {
      if (!gGrid)
      {
        TGrid::Connect("alien://");
      }
      if (!gGrid)
      {
        Error("FILTER_RAWMUON","Cannot connect to Grid !");
        return -1;
      }
    }
    
    TFile* rawFile = TFile::Open(from);
    
    if (!rawFile || !rawFile->IsOpen())
    {
      Error("FILTER_RAWMUON","Cannot open input file %s",from);
      return -3;
    }
    
    TTree *rawTree=(TTree *)rawFile->Get("RAW");
    if(!rawTree)
    {
      Error("FILTER_RAWMUON","Error getting RAW tree from file %s",from);
      return -2;
    }
    
    TString snewfile(to);
    
    RAWMUON::DisableBranches(rawTree);
    
    TFile* newfile =  new TFile(snewfile.Data(),"recreate");
    
    TTree* newTree = rawTree->CloneTree();
    
    newTree->Print();
    
    newTree->AutoSave();
    
    newfile->Close();
    
    delete newfile;
    
    delete rawFile;
    
    return 0;
  }
  
}
