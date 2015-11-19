// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TMap.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
#include <TRegexp.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>
#include <TEveTreeTools.h>

#include <AliAltroRawStreamV3.h>
#include <AliEveEventManager.h>
#include <AliEveTPCData.h>
#include <AliEveTPCSector2D.h>
#include <AliEveTPCSector3D.h>
#include <AliGeomManager.h>
#include <AliMagF.h>
#include <AliRawReader.h>

#include <AliGRPObject.h>
#include <AliCDBManager.h>
#include <AliCDBPath.h>
#include <AliCDBEntry.h>

#include <STEER/STEER/AliGRPManager.h>

#include <TPC/AliTPCcalibDB.h>
#include <TPC/AliTPCTransform.h>
#include <TPC/AliTPCclusterMI.h>
#include <TPC/AliTPCClustersRow.h>

#include <TPC/AliTPCParam.h>
#include <TPC/AliTPCRawStreamV3.h>
#include <TPC/AliTPCRecoParam.h>
#include <TPC/AliTPCReconstructor.h>

#include <HLT/BASE/AliHLTOUT.h>
#include <HLT/BASE/AliHLTSystem.h>
#include <HLT/BASE/AliHLTPluginBase.h>

#include <RAW/AliRawHLTManager.h>
#endif

// read TPC clusters compressed by HLT from RawReader and return output in OutputTree
TTree* readHLTClusters(AliRawReader *RawReader);

// add the HLT clusters from clustersTree to visualization
void renderHLTClusters(TTree* clustersTree);


// Macro to visualise rootified raw-data from TPC.
//
// Use tpc_raw(Int_t mode) in order to run it
// Needs that alieve_init() is already called
// mode = 1 - show only 2D sectors
// mode = 2 - show only 3D sectors
// mode = 3 - show both 2D and 3D sectors
void tpc_raw(Int_t mode = 3)
{
    printf("*** RAW TPC ***");
    
  gStyle->SetPalette(1, 0);

  AliEveEventManager::GetMaster()->AssertGeometry();
  AliEveEventManager::AssertMagField();
  
  AliRawReader *reader = AliEveEventManager::AssertRawReader();
  reader->Reset();
  
  AliTPCRawStreamV3 input(reader);
  reader->Select("TPC"); // ("TPC", firstRCU, lastRCU);

  AliEveTPCData *x = new AliEveTPCData;
  // x->SetLoadPedestal(5);
  x->SetLoadThreshold(5);
  x->SetAutoPedestal(kTRUE);

  x->LoadRaw(input, kTRUE, kTRUE);

  gEve->DisableRedraw();

  TEveElementList* sec2d = new TEveElementList("TPC 2D");
  gEve->AddElement(sec2d);

  TEveElementList* sec3d = new TEveElementList("TPC 3D");
  gEve->AddElement(sec3d);

  AliEveTPCSector2D *s;
  AliEveTPCSector3D *t;
  
  for (Int_t i=0; i<=35; ++i) {
    if (mode & 1) {
      s = new AliEveTPCSector2D(Form("2D sector %d",i));
      s->SetSectorID(i);
      s->SetAutoTrans(kTRUE); // place on proper 3D coordinates
      s->SetDataSource(x);
      s->SetFrameColor(36);
      sec2d->AddElement(s);
      s->IncRTS();
    }
    if (mode & 2) {
      t = new AliEveTPCSector3D(Form("3D sector %d",i));
      t->SetSectorID(i);
      t->SetAutoTrans(kTRUE);
      t->SetDataSource(x);
      sec3d->AddElement(t);
      t->IncRTS();
    }
  }
  
  // Display TPC clusters compressed by HLT
  TTree* hltClustersTree = readHLTClusters(reader); // read HLT compressed clusters from TPC from raw reader and output them in hltClustersTree
  if(hltClustersTree) renderHLTClusters(hltClustersTree);
  
  gEve->EnableRedraw();
  gEve->Redraw3D();
  
}

// read TPC Clusters compressed by HLT from a raw filename
TTree* readHLTClusters(AliRawReader *fRawReader)
{
/*************************
*  Raw IO Stuff
**************************/
  AliRawReader* rawReader = AliRawHLTManager::CreateRawReaderHLT(fRawReader, "TPC");
  rawReader->Select("TPC");

/*************************
*  HLT IO Stuff
**************************/
  AliHLTSystem* hltSystem=AliHLTPluginBase::GetInstance();

  AliHLTOUT* hltOut = AliHLTOUT::New(rawReader); 
  hltOut->Init();
  hltSystem->InitHLTOUT(hltOut);
  
/*************************
*  GRP Stuff
**************************/
  AliGRPManager manGRP;
  AliRecoParam* recoParam = InitRecoParam();
  
  AliEveEventManager* curEventMan = (AliEveEventManager*)gEve->GetCurrentEvent();
  
  const AliRunInfo* runInfo = manGRP.GetRunInfo();
  const THashTable* cosmicTriggers = manGRP.GetCosmicTriggers();
  const AliEventInfo* eventinfo = GetEventInfo();
  
  recoParam->SetEventSpecie(runInfo,*eventinfo,cosmicTriggers);
  
  AliTPCRecoParam *tpcRecoParam = (AliTPCRecoParam*)recoParam->GetDetRecoParam(1); // TPC has index 1
  tpcRecoParam->SetUseHLTClusters(3); // reconstruct only HLT clusters
	
/**************************
*  Reconstruction of Clusters
**************************/
  TTree* outputClustersTree = new TTree;
  
  AliTPCReconstructor *reconstructor = new AliTPCReconstructor();
  reconstructor->SetOption("useHLT");
  reconstructor->CreateTracker(); // this will set the option to reconstruct for only HLT clusters
  
  reconstructor->SetRecoParam(tpcRecoParam);
  reconstructor->SetRunInfo((AliRunInfo*)runInfo);
  reconstructor->SetEventInfo((AliEventInfo*)eventinfo);
  
  reconstructor->Reconstruct(rawReader, outputClustersTree);
  
  delete reconstructor;
  
  hltSystem->ReleaseHLTOUT(hltOut);
  
  return outputClustersTree;
}
 
void renderHLTClusters(TTree* clustersTree)
{
 
/**************************
*  Visualization of Clusters
**************************/
  const Int_t kMaxCl=100*160;
  Int_t fNColorBins = 5;
  
  TEvePointSet* clusters = new TEvePointSet(kMaxCl);
  clusters->SetOwnIds(kTRUE);

  TEvePointSetArray * cc = new TEvePointSetArray("TPC Clusters Colorized");
  cc->SetMainColor(kRed);
  cc->SetMarkerStyle(4);
  cc->SetMarkerSize(0.4);
  cc->InitBins("Cluster Charge", fNColorBins, 0., fNColorBins*60.);
  
  cc->GetBin(0)->SetMainColor(kGray);
  cc->GetBin(0)->SetMarkerSize(0.4);
  cc->GetBin(1)->SetMainColor(kBlue);
  cc->GetBin(1)->SetMarkerSize(0.42);
  cc->GetBin(2)->SetMainColor(kCyan);
  cc->GetBin(2)->SetMarkerSize(0.44);
  cc->GetBin(3)->SetMainColor(kGreen);
  cc->GetBin(3)->SetMarkerSize(0.46);
  cc->GetBin(4)->SetMainColor(kYellow);
  cc->GetBin(4)->SetMarkerSize(0.48);
  cc->GetBin(5)->SetMainColor(kRed);
  cc->GetBin(5)->SetMarkerSize(0.50);
  cc->GetBin(6)->SetMainColor(kMagenta);
  cc->GetBin(6)->SetMarkerSize(0.52);
  

// Loop over clusters
  Int_t nentries = clustersTree->GetEntriesFast();
	
  AliTPCClustersRow *clrow = new AliTPCClustersRow();
  clrow->SetClass("AliTPCclusterMI");
  //clrow->SetArray(kMaxCl);
  clustersTree->SetBranchAddress("Segment", &clrow);
  
  for (Int_t i=0; i<nentries; i++) {
    if (!clustersTree->GetEvent(i)) continue;
    
    TClonesArray *cl = clrow->GetArray();
    Int_t ncl = cl->GetEntriesFast();

    while (ncl--)
    {
      AliTPCclusterMI* clusterMI = (AliTPCclusterMI*) cl->At(ncl);

      AliCluster *c = (AliCluster*) cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
   
      cc->Fill(g[0], g[1], g[2], clusterMI->GetQ());
      clusters->SetNextPoint(g[0], g[1], g[2]);
      AliCluster *atp = new AliCluster(*clusterMI);
      clusters->SetPointId(atp);

    }
    cl->Clear();
  }

  delete clrow;
	  
  clusters->SetName("TPC Clusters");
  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters TPC"); // to be changed
  clusters->ApplyVizTag(viz_tag, "Clusters");
  
  cc->SetRnrSelf(kTRUE);

  gEve->AddElement(cc);

  return;
}

const AliEventInfo* GetEventInfo()
{
    // Fill the event info object
    
    AliCentralTrigger *aCTP = NULL;
    if (AliEveEventManager::AssertRawReader()) {
        fEventInfo.SetEventType(AliEveEventManager::AssertRawReader()->GetType());
        
        ULong64_t mask = AliEveEventManager::AssertRawReader()->GetClassMask();
        fEventInfo.SetTriggerMask(mask);
        UInt_t clmask = AliEveEventManager::AssertRawReader()->GetDetectorPattern()[0];
        fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(clmask));
        
        aCTP = new AliCentralTrigger();
        TString configstr("");
        if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
            AliError("No trigger configuration found in OCDB! The trigger configuration information will not be used!");
            delete aCTP;
            return 0;
        }
        aCTP->SetClassMask(mask);
        aCTP->SetClusterMask(clmask);
        
        if (AliEveEventManager::AssertRunLoader()) {
            AliCentralTrigger* rlCTP = AliEveEventManager::AssertRunLoader()->GetTrigger();
            if (rlCTP) {
                rlCTP->SetClassMask(mask);
                rlCTP->SetClusterMask(clmask);
            }
        }
    }
    else {
        fEventInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);
        
        if (AliEveEventManager::AssertRunLoader() && (!AliEveEventManager::AssertRunLoader()->LoadTrigger())) {
            aCTP = AliEveEventManager::AssertRunLoader()->GetTrigger();
            fEventInfo.SetTriggerMask(aCTP->GetClassMask());
            // get inputs from actp - just get
            AliESDHeader* esdheader = AliEveEventManager::AssertESD()->GetHeader();
            esdheader->SetL0TriggerInputs(aCTP->GetL0TriggerInputs());
            esdheader->SetL1TriggerInputs(aCTP->GetL1TriggerInputs());
            esdheader->SetL2TriggerInputs(aCTP->GetL2TriggerInputs());
            fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(aCTP->GetClusterMask()));
        }
        else {
            AliWarning("No trigger can be loaded! The trigger information will not be used!");
            return 0;
        }
    }
    
    AliTriggerConfiguration *config = aCTP->GetConfiguration();
    if (!config) {
        AliError("No trigger configuration has been found! The trigger configuration information will not be used!");
        if (AliEveEventManager::AssertRawReader()) delete aCTP;
        return 0;
    }
    
    TString declTriggerClasses;
    
    // Load trigger aliases and declare the trigger classes included in aliases
    AliCDBEntry * entry = AliCDBManager::Instance()->Get("GRP/CTP/Aliases");
    if (entry) {
        THashList * lst = dynamic_cast<THashList*>(entry->GetObject());
        if (lst) {
            lst->Sort(kSortDescending); // to avoid problems with substrings
            if (AliEveEventManager::AssertRawReader()) AliEveEventManager::AssertRawReader()->LoadTriggerAlias(lst);
            // Now declare all the triggers present in the aliases
            TIter iter(lst);
            TNamed *nmd = 0;
            while((nmd = dynamic_cast<TNamed*>(iter.Next()))){
                declTriggerClasses += " ";
                declTriggerClasses += nmd->GetName();
            }
        }
        else {
            AliError("Cannot cast the object with trigger aliases to THashList!");
        }
    }
    else {
        AliError("No OCDB entry for the trigger aliases!");
    }
    
    // Load trigger classes for this run
    UChar_t clustmask = 0;
    TString trclasses;
    ULong64_t trmask = fEventInfo.GetTriggerMask();
    const TObjArray& classesArray = config->GetClasses();
    Int_t nclasses = classesArray.GetEntriesFast();
    for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
        AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
        if (trclass && trclass->GetMask()>0) {
            Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
            if (AliEveEventManager::AssertESD()) AliEveEventManager::AssertESD()->SetTriggerClass(trclass->GetName(),trindex);
            if (AliEveEventManager::AssertRawReader()) AliEveEventManager::AssertRawReader()->LoadTriggerClass(trclass->GetName(),trindex);
            if (trmask & (1ull << trindex)) {
                trclasses += " ";
                trclasses += trclass->GetName();
                trclasses += " ";
                clustmask |= trclass->GetCluster()->GetClusterMask();
            }
        }
    }
    fEventInfo.SetTriggerClasses(trclasses);
    
    if (!aCTP->CheckTriggeredDetectors()) {
        if (AliEveEventManager::AssertRawReader()) delete aCTP;
        return 0;
    }
    
    if (AliEveEventManager::AssertRawReader()) delete aCTP;
    
    // everything went ok, return pointer
    return (&fEventInfo);
}

AliRecoParam* AliEveEventManager::InitRecoParam()
{
    // This is mostly a reap-off from reconstruction
    // The method accesses OCDB and retrieves all
    // the available reco-param objects from there.
    
    AliRecoParam *fgRecoParam = new AliRecoParam;
    const Int_t  kNDetectors = 14;
    
    static const TEveException kEH("AliEveEventManager::InitRecoParam");
    
    Bool_t isOK = kTRUE;
    
    if (fgRecoParam->GetDetRecoParamArray(kNDetectors)) {
        ::Info(kEH, "Using custom GRP reconstruction parameters");
    }
    else {
        ::Info(kEH, "Loading GRP reconstruction parameter objects");
        
        AliCDBPath path("GRP","Calib","RecoParam");
        AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
        if(!entry){
            ::Warning(kEH, "Couldn't find GRP RecoParam entry in OCDB");
            isOK = kFALSE;
        }
        else {
            TObject *recoParamObj = entry->GetObject();
            if (dynamic_cast<TObjArray*>(recoParamObj)) {
                // GRP has a normal TobjArray of AliDetectorRecoParam objects
                // Registering them in AliRecoParam
                fgRecoParam->AddDetRecoParamArray(kNDetectors,dynamic_cast<TObjArray*>(recoParamObj));
            }
            else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
                // GRP has only onse set of reco parameters
                // Registering it in AliRecoParam
                ::Info(kEH, "Single set of GRP reconstruction parameters found");
                dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
                fgRecoParam->AddDetRecoParam(kNDetectors,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
            }
            else {
                ::Error(kEH, "No valid GRP RecoParam object found in the OCDB");
                isOK = kFALSE;
            }
            entry->SetOwner(0);
        }
    }
    
    const char* fgkDetectorName[kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE" };
    
    
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
        
        if (fgRecoParam->GetDetRecoParamArray(iDet)) {
            ::Info(kEH, "Using custom reconstruction parameters for detector %s",fgkDetectorName[iDet]);
            continue;
        }
        
        ::Info(kEH, "Loading reconstruction parameter objects for detector %s",fgkDetectorName[iDet]);
        
        AliCDBPath path(fgkDetectorName[iDet],"Calib","RecoParam");
        AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
        if(!entry){
            ::Warning(kEH, "Couldn't find RecoParam entry in OCDB for detector %s",fgkDetectorName[iDet]);
            isOK = kFALSE;
        }
        else {
            TObject *recoParamObj = entry->GetObject();
            if (dynamic_cast<TObjArray*>(recoParamObj)) {
                // The detector has a normal TobjArray of AliDetectorRecoParam objects
                // Registering them in AliRecoParam
                fgRecoParam->AddDetRecoParamArray(iDet,dynamic_cast<TObjArray*>(recoParamObj));
            }
            else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
                // The detector has only onse set of reco parameters
                // Registering it in AliRecoParam
                ::Info(kEH, "Single set of reconstruction parameters found for detector %s",fgkDetectorName[iDet]);
                dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
                fgRecoParam->AddDetRecoParam(iDet,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
            }
            else {
                ::Error(kEH, "No valid RecoParam object found in the OCDB for detector %s",fgkDetectorName[iDet]);
                isOK = kFALSE;
            }
            entry->SetOwner(0);
            
        }
    }
    
    if(!isOK) {
        delete fgRecoParam;
        fgRecoParam = 0;
    }
    
    return fgRecoParam;
}


