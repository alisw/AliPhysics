// $Id$
/**
 * @file makeGlobalHistoConfigObject.C
 * @brief Creation of Global Histo component configuration objects in OCDB
 *
 * <pre>
 * Usage:
 *  aliroot -b -q makeGlobalHistoConfigObject.C'("path", "key", "uri", runmin, runmax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param.
 * Many HLT components understand configuration strings containing
 * arguments and parameters just like the command line arguments.
 * This macro facilitates the creation of an appropriate object
 * from a parameter string.
 * As another approach the TObjString parameters are stored in a TMap
 * associated to a key. A TMap object is generated if 'key' is specified.
 *
 * Parameters: <br>
 * - path           path of the entry within the OCDB (e.g. HLT/ConfigHLT/GlobalHisto)
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   
 * - runmin (opt)   default 0
 * - runmax (opt)   default 999999999
 * The parameters of the component configuration should be defined inside function makeString().
 *
 * Note: The configuration procedure of an HLT component is not
 * restricted to that scheme. The implementation is up to the
 * developer and more complex objects are possible.
 *
 * @author Kalliopi.Kanaki@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void makeGlobalHistoConfigObject(const char* path, 
				 const char* key,
				 //const char* param,
				 const char* cdbUri,
				 int runmin=0,
				 int runmax=999999999,
				 int runNo=0)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "Cannot get AliCDBManager." << end;
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

  TMap* pMap=NULL;

  // load existing object and init TMap
  AliCDBEntry* pExisting=NULL;
  AliCDBStorage* pStorage=AliCDBManager::Instance()->GetDefaultStorage();
  if (key && pStorage->GetLatestVersion(path, runNo)>=0) {
    pExisting=pStorage->Get(path, runNo);
    if (pExisting->GetObject()->IsA() == TMap::Class()) {
      pMap=(TMap*)pExisting->GetObject()->Clone();
    }
  }  

  if (key && !pMap) pMap=new TMap;

  // here is the actual content of the configuration object
  TObject* obj=new TObjString(makeString()); //KK
  //TObject* obj=new TObjString(param); 
  if (pMap) {
    if (pMap->FindObject(key)) {
      pMap->Remove(new TObjString(key));
    }
    pMap->Add(new TObjString(key), obj);
    obj=pMap;
  }

  AliCDBPath cdbPath(path);
  AliCDBId cdbId(cdbPath, runmin, runmax);
  AliCDBMetaData* cdbMetaData=NULL;
  if (pExisting) cdbMetaData=pExisting->GetMetaData();
  else cdbMetaData=new AliCDBMetaData;
  man->Put(obj, cdbId, cdbMetaData);
}

// void makeGlobalHistoConfigObject(const char* path
//                                       ,//const char* param="",
// 				      ,const char* cdbUri=NULL
// 				      ,int runmin=0
// 				      ,int runmax=999999999)
// {
//   //makeGlobalHistoConfigObject(path, NULL, param, cdbUri, runmin, runmax);
//   makeGlobalHistoConfigObject(path, NULL, makeString(), cdbUri, runmin, runmax); //KK
// }
// 
// void makeComponentConfigurationObject(const char* path, 
// 				      int runNo,
// 				      const char* key,
// 				      const char* param)
// {
//   makeComponentConfigurationObject(path, key, param, NULL, 0, 999999999, runNo);
// }

void makeGlobalHistoConfigObject()
{
  cout << "===============================================================" << endl;
  cout << "usage: aliroot -b -q -l makeGlobalHistoConfigObject.C'(\"path\", NULL, \"uri\", rangemin, rangemax)'" << endl << endl;
  cout << "  path           path of the entry within the OCDB, e.g. HLT/ConfigHLT/GlobalHisto" << endl;
  cout << "  uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   " << endl;
  cout << "  rangemin (opt) default 0" << endl;
  cout << "  rangemax (opt) default 999999999" << endl;
  cout << "  The parameters of the component configuration should be defined inside function makeString()." << endl;
  cout << "===============================================================" << endl;
}

TString makeString(){
  
  TString s = "";

  s+="-max-track-count 8000 ";
  //s+="-max-track-count 8000 -fill-V0 -max-V0-count 200 ";
  
  s+="-histogram TrackPt(100,0,20) -size 1000 -expression Track_pt -title p_{T}_[GeV/c] -cut Track_TPCclus>0 ";
  s+="-histogram TrackPhi(180,0,360) -size 1000 -expression Track_phi -title #phi_(deg) -cut Track_TPCclus>0 ";
  s+="-histogram TrackMultiplicity(250,0,5000) -size 1000 -expression trackcount -title TrackMultiplicity ";
  s+="-histogram TrackEta(100,-2,2) -size 1000 -expression Track_eta -title #eta -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackTPCclus(200,0,200) -size 1000 -expression Track_TPCclus -title TPC_clusters/track -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackITSclus(7,0,7) -size 1000 -expression Track_ITSclus -title ITS_clusters/track ";
  s+="-histogram TrackTheta(90,0,180) -size 1000 -expression Track_theta -title #theta_(deg) -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackDCAr(100,-50,50) -size 1000 -expression Track_DCAr -title DCAr_[cm] -cut Track_TPCclus>0 ";
  s+="-histogram TrackCharge -size 1000 -expression Track_charge -title Polarity -cut Track_TPCclus>0 ";
 
  s+="-histogram VertexXY -size 1000 -expression vertexY:vertexX -title y:x_[cm] -cut nContributors>3 -opt colz ";
  s+="-histogram VertexX(50,-5,5)  -size 1000 -expression vertexX -title x_[cm] -cut nContributors>3 ";
  s+="-histogram VertexY(50,-5,5)  -size 1000 -expression vertexY -title y_[cm] -cut nContributors>3 ";
  s+="-histogram VertexZ(200,-50,50)  -size 1000 -expression vertexZ -title x_[cm] -cut nContributors>3 ";
  s+="-histogram VertexTrendX -size 1000 -expression vertexX:event -title vertexX_vs_event -cut nContributors>3 ";
  s+="-histogram VertexTrendY -size 1000 -expression vertexY:event -title vertexY_vs_event -cut nContributors>3 ";

  
  return s;
}
