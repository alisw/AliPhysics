/// $Id$
///
/// Generate, for a given slat type, the list of translators and their buspatch
///
/// It's the equivalent of the St#_Nappes-crocus-v#.#.pdf files found in 
/// https://twiki.cern.ch/twiki/bin/view/ALICE/St345CrocusFlatCables
/// but sorted by slat type, and not by crocus.
///
/// \author Laurent Aphecetche, Subatech
///

#include "TString.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include "Riostream.h"
#include "TSystem.h"
#include <vector>
#include <utility>
#include "AliMpDetElement.h"
#include "AliMpDDLStore.h"
#include "AliMpCDB.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"
#include "AliMpBusPatch.h"
#include <algorithm>

void LoadMapping(Bool_t fromFile)
{
  if ( fromFile ) 
  {
    AliMpDataProcessor mp;
    {
      AliMpDataMap* datamap = mp.CreateDataMap("data");
      AliMpDataStreams dataStreams(datamap);
      AliMpDDLStore::ReadData(dataStreams);
    }
    {
      AliMpDataMap* datamap = mp.CreateDataMap("data_run");
      AliMpDataStreams dataStreams(datamap);
      AliMpManuStore::ReadData(dataStreams);
    }    
  }
  else
  {
    AliMpCDB::LoadAll2();
  }
}

void SlatTranslatorToBusPatches(const char* whichSlat="112200NR2")
{
  LoadMapping(kTRUE);
  
  TString file(gSystem->ExpandPathName("$ALICE_ROOT/MUON/mapping/data/station345/DetElemIdToSlatType.dat"));
  
  ifstream in(file);
  if (in.bad()) return;
  
  int detElemId;
  char s[80];
  char slatType[80];
  std::map<std::string,std::vector<int> > slats;
  
  while ( in.getline(s,80,'\n') )
  {
    if ( s[0] != '#' && strlen(s) > 2 )
    {
      sscanf(s,"%d %s",&detElemId,slatType);
      
      slats[slatType].push_back(detElemId);
    }
  }
  
  std::map<std::string,std::vector<int> >::const_iterator it;
  
  std::vector<int> v = slats[whichSlat];
  
  cout << "----------------------" << endl;
  cout << whichSlat << endl;
  cout << "----------------------" << endl;
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(v[0]);
  Int_t nofManus(0);
  Int_t nofChannels(0);
  
  for ( Int_t b = 0; b < de->GetNofBusPatches(); ++b ) 
  {
    Int_t busPatchId = de->GetBusPatchId(b);
    AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
    nofManus += bp->GetNofManus();
    cout << bp->GetTranslatorLabel() << " : " << bp->GetNofManus() << " manus" << endl;
    for ( Int_t im = 0; im < bp->GetNofManus(); ++im )
    {
      Int_t manuId = bp->GetManuId(im);
      nofChannels += de->NofChannelsInManu(manuId);
    }
  }

  cout << "Number of bus patches = " << de->GetNofBusPatches() << endl;
  cout << "Number of manus = " << nofManus << endl;
  cout << "Number of channels = " << nofChannels << endl;
  cout << "----------------------" << endl;

  std::sort(v.begin(),v.end());
  
  for ( size_t i = 0; i < v.size(); ++i ) 
  {
    Int_t detElemId = v[i];
    
    cout << Form("%04d ",detElemId) << endl;

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    for ( Int_t b = 0; b < de->GetNofBusPatches(); ++b ) 
    {
      Int_t busPatchId = de->GetBusPatchId(b);
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
      bp->Print();
    }
  }
  cout << endl;

  in.close();
}



