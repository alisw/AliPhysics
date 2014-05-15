/**
 * @file   Raw2ESD.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 23 15:06:55 2014
 * 
 * @brief  Convert from raw data to ESD using the AliFMDReconstructor
 * 
 * 
 */
/** 
 * Default input file 
 * @ingroup FMD_script
 */
const char* df = "/data/alice/data/pp/LHC10c/raw/118561/physics_118561.root";

/** 
 * Convert from raw data to ESD using the AliFMDReconstructor.  This
 * illustrates the passes done in the official reconstruction.
 * 
 * @param file   Input raw data
 * @param nEv    Number of events to process (<=0 means all)
 * @param skip   Number of events to skip 
 * @param debug  Debug level
 *
 * @ingroup FMD_script
 */
void
Raw2ESD(const char* file=df, Int_t nEv=10, Int_t skip=300, Int_t debug=0)
{
  // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(118561);
  AliCDBManager::Instance()->SetDefaultStorageFromRun(118561);
  AliGeomManager::LoadGeometry("geometry.root");

  AliRawReader*        reader  = AliRawReader::Create(file);
  AliFMDReconstructor* reco    = new AliFMDReconstructor();
  reco->SetDiagnose(debug > 5);
  reco->Init();

  AliLog::SetModuleDebugLevel("FMD", debug);

  Int_t        event       = 0;
  TFile*       digitFile   = TFile::Open("reco_digits.rot", "RECREATE");

  TFile*       clusterFile   = TFile::Open("FMD.RecPoints.root", "RECREATE");
  TTree*       clusterTree = new TTree("cluster", "FMD digits");

  TFile*       esdFile     = TFile::Open("AliESDs.root", "RECREATE");
  TTree*       esdTree     = new TTree("esdTree", "ESD Treee");
  AliESDEvent* esd         = new AliESDEvent();

  esd->CreateStdContent();
  esd->WriteToTree(esdTree);
  while ((reader && reader->NextEvent())) {
    event++;
    // Check for skip events 
    if (skip > 0 && (event - skip) < 0) continue;
    // Check if we got enough events 
    if (nEv > 0 && (event-skip) > nEv) continue;

    digitFile->cd();
    TTree* digitTree  = new TTree("digit", "FMD digits");
    
    // Convert to digits first 
    reco->ConvertDigits(reader, digitTree);

    // Reconstruct to RecPoints 
    reco->Reconstruct(digitTree, clusterTree);

    // Set stuff on the ESD 
    esd->SetRunNumber(AliCDBManager::Instance()->GetRun());
    esd->SetEventNumberInFile(event-1);
    
    // Fill the ESD objet and write to file 
    reco->FillESD((TTree*)0, (TTree*)0, esd);
    esdTree->Fill();
    esd->Reset();

    digitFile->Write();
    delete digitTree;
  }
  esdFile->Write();
}

  
