//------------------------------------------------------------------------------
// Implementation of AliPerformanceMatch class. It checks the track matching 
// quality (residuals, pulls) bewteen tracking detectors TPC-ITS and TPC-TRD 
// with TPC as the reference detector. In addition, the ITS and TRD track 
// reconstruction efficiency with respect to TPC is calculated.
//
// The matchinig quality parameters are stored in the THnSparse histograms. 
// The analysis of these histograms to extract reference information is done in 
// the Analyse() function (post-processing).  
//  
// The histograms with reference information can be exported to the ROOT folder.
//
// Author: J.Otwinowski 17/10/2009 
// Changes by M.Knichel 22/10/2010 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceMatch * compObj = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCTRD");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderMatch" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Match.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"

#include "AliPerformanceMatch.h" 
#include "AliESDEvent.h" 
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliESDfriend.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliStack.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTRDtrackV1.h" 
#include "AliTreeDraw.h" 

using namespace std;

ClassImp(AliPerformanceMatch)

Bool_t AliPerformanceMatch::fgMergeTHnSparse = kFALSE;
Bool_t AliPerformanceMatch::fgUseMergeTHnSparse = kFALSE;

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch():
  AliPerformanceObject("AliPerformanceMatch"),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fFolderObj(0),
  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  
  fUseHLT(0)
{
  Init();
}

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch(Char_t* name="AliPerformanceMatch", Char_t* title="AliPerformanceMatch",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fFolderObj(0),
  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  fUseHLT(0)
{
  // named constructor	
  // 
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);

  Init();
}

//_____________________________________________________________________________
AliPerformanceMatch::~AliPerformanceMatch()
{
  // destructor
   
  if(fResolHisto) delete fResolHisto; fResolHisto=0;     
  if(fPullHisto)  delete fPullHisto;  fPullHisto=0;
  if(fTrackingEffHisto) delete fTrackingEffHisto; fTrackingEffHisto = 0x0;
  
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
  if(fFolderObj) delete fFolderObj; fFolderObj=0;
}

//_____________________________________________________________________________
void AliPerformanceMatch::Init(){

  //
  // Make performance histogrms
  //

  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 20.;

  Double_t *binsPt = 0;

  if (IsHptGenerator())  { 
        ptMax = 100.;
  } 
   binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);


  Double_t yMin = -0.02, yMax = 0.02;
  Double_t zMin = -12.0, zMax = 12.0;
  if(GetAnalysisMode() == 0 || GetAnalysisMode() == 1) { 
    yMin = -100.; yMax = 100.; 
    zMin = -100.; zMax = 100.; 
  }

  // res_y:res_z:res_phi,res_lambda:res_pt:y:z:eta:phi:pt:isRec
  Int_t binsResolHisto[9]={100,100,100,100,100,90,30,nPtBins,2};
  Double_t minResolHisto[9]={-1.,-1.,-0.03,-0.03,-0.2, 0., -1.5, ptMin,0};
  Double_t maxResolHisto[9]={ 1., 1., 0.03, 0.03, 0.2, 2.*TMath::Pi(), 1.5, ptMax,2};

  fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_pt:phi:eta:pt:isRec",9,binsResolHisto,minResolHisto,maxResolHisto);
  fResolHisto->SetBinEdges(7,binsPt);

  fResolHisto->GetAxis(0)->SetTitle("y-y_{ref} (cm)");
  fResolHisto->GetAxis(1)->SetTitle("z-z_{ref} (cm)");
  fResolHisto->GetAxis(2)->SetTitle("#phi-#phi_{ref} (rad)");
  fResolHisto->GetAxis(3)->SetTitle("#lambda-#lambda_{ref} (rad)");
  fResolHisto->GetAxis(4)->SetTitle("(p_{T}/p_{Tref}-1)");
  fResolHisto->GetAxis(5)->SetTitle("#phi_{ref} (rad)");
  fResolHisto->GetAxis(6)->SetTitle("#eta_{ref}");
  fResolHisto->GetAxis(7)->SetTitle("p_{Tref} (GeV/c)");
  fResolHisto->GetAxis(8)->SetTitle("isReconstructed");
  fResolHisto->Sumw2();

  //pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt:isRec
  Int_t binsPullHisto[9]={100,100,100,100,100,90,30,nPtBins,2};
  Double_t minPullHisto[9]={-5.,-5.,-5.,-5.,-5., 0.,-1.5, ptMin,0};
  Double_t maxPullHisto[9]={ 5., 5., 5., 5., 5., 2.*TMath::Pi(), 1.5, ptMax,2};
  fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:phi:eta:1pt:isRec",9,binsPullHisto,minPullHisto,maxPullHisto);

  fPullHisto->GetAxis(0)->SetTitle("(y-y_{ref})/#sigma");
  fPullHisto->GetAxis(1)->SetTitle("(z-z_{ref})/#sigma");
  fPullHisto->GetAxis(2)->SetTitle("(sin#phi-sin#phi_{ref})/#sigma");
  fPullHisto->GetAxis(3)->SetTitle("(tan#lambda-tan#lambda_{ref})/#sigma");
  fPullHisto->GetAxis(4)->SetTitle("(p_{Tref}/p_{T}-1)/#sigma");
  fPullHisto->GetAxis(5)->SetTitle("#phi_{ref} (rad)");
  fPullHisto->GetAxis(6)->SetTitle("#eta_{ref}");
  fPullHisto->GetAxis(7)->SetTitle("1/p_{Tref} (GeV/c)^{-1}");
  fPullHisto->GetAxis(8)->SetTitle("isReconstructed");
  fPullHisto->Sumw2();

  // -> has match:y:z:snp:tgl:phi:pt:ITSclusters
  Int_t binsTrackingEffHisto[5]    = { 2,   90,           100, 30,    7   };
  Double_t minTrackingEffHisto[5]  = {-0.5, 0.,            0, -1.5, -0.5 };
  Double_t maxTrackingEffHisto[5]  = { 1.5, 2*TMath::Pi(), 20, 1.5,  6.5 };
  
  fTrackingEffHisto = new THnSparseF("fTrackingEffHisto","has match:phi:pt:eta:ITSclusters",5,binsTrackingEffHisto,minTrackingEffHisto,maxTrackingEffHisto);
  fTrackingEffHisto->GetAxis(0)->SetTitle("IsMatching");
  fTrackingEffHisto->GetAxis(1)->SetTitle("phi (rad)");
  fTrackingEffHisto->GetAxis(2)->SetTitle("p_{T}");
  fTrackingEffHisto->GetAxis(3)->SetTitle("eta");
  fTrackingEffHisto->GetAxis(4)->SetTitle("number of ITS clusters");
  fTrackingEffHisto->Sumw2();

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderMatch","Analysis Matching Folder");
  
   // save merge status in object
  fMergeTHnSparseObj = fgMergeTHnSparse;

}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessITSTPC(Int_t iTrack, AliESDEvent *const esdEvent, AliStack* /*const stack*/, AliESDtrack *const esdTrack)
{
  //
  // addition to standard analysis - check if ITS stand-alone tracks have a match in the TPC
  // Origin: A. Kalwait
  // Modified: J. Otwinowski
  if(!esdEvent) return;
  if(!esdTrack) return;
  //  if(!esdFriendTrack) return;

  // ITS stand alone tracks with SPD layers 
  if(!(esdTrack->GetStatus() & AliESDtrack::kITSpureSA)) return;
  if(!(esdTrack->GetStatus() & AliESDtrack::kITSrefit)) return;
  if(esdTrack->GetNcls(0)<4) return;
   if(!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1)) return;
  
  const AliESDVertex* vtxESD = esdEvent->GetPrimaryVertexTracks();
  //   const AliESDVertex *vtxESD = esdEvent->GetPrimaryVertexTPC();
 AliESDtrack* tpcTrack2 = NULL;
  Bool_t hasMatch = kFALSE;
  for (Int_t jTrack = 0; jTrack < esdEvent->GetNumberOfTracks(); jTrack++) {
    // loop for all those tracks and check if the corresponding TPC track is found
    if (jTrack==iTrack) continue;
    AliESDtrack *trackTPC = esdEvent->GetTrack(jTrack);
    if (!trackTPC) continue;
    if (!trackTPC->GetTPCInnerParam()) continue;
    if(!(trackTPC->GetStatus() & AliESDtrack::kTPCrefit)) continue;
    
    // TPC nClust/track after first tracking pass
    // if(trackTPC->GetTPCNclsIter1()<fCutsRC->GetMinNClustersTPC()) continue;
    tpcTrack2 = AliESDtrackCuts::GetTPCOnlyTrack(esdEvent, jTrack);
    if(!tpcTrack2) continue;
    if(!tpcTrack2->RelateToVertex(vtxESD,esdEvent->GetMagneticField(),100.)) { delete tpcTrack2; continue; } 
    
     if(!fCutsRC->AcceptTrack(tpcTrack2)) { delete tpcTrack2; continue; }
    // check matching
    if (TMath::Abs(esdTrack->GetY() - tpcTrack2->GetY()) > 3) { delete tpcTrack2; continue; }
    if (TMath::Abs(esdTrack->GetSnp() - tpcTrack2->GetSnp()) > 0.2) { delete tpcTrack2; continue; }
    if (TMath::Abs(esdTrack->GetTgl() - tpcTrack2->GetTgl()) > 0.2) { delete tpcTrack2; continue; }
    
    hasMatch=kTRUE;
    break;
  }
  
  FillHistograms(tpcTrack2,esdTrack,hasMatch);     
  /*  if(tpcTrack2) { 
    delete tpcTrack2;
   
    }*/
}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessTPCITS(AliStack* /*const stack*/, AliESDtrack *const esdTrack)
{
  //
  // Match TPC and ITS min-bias tracks
  // at radius between detectors
  //
  if(!esdTrack) return;
  //   if(!esdFriendTrack) return;
   
  Bool_t isTPC = kFALSE;
  Bool_t isMatch = kFALSE;


  if(esdTrack->Charge()==0) return;
  if(!esdTrack->GetTPCInnerParam()) return;
  if(!(esdTrack->GetStatus()&AliESDtrack::kTPCrefit)) return;

  if(!fCutsRC->AcceptTrack(esdTrack)) return;

  isTPC = kTRUE;
  
  if( (esdTrack->GetStatus()&AliESDtrack::kITSrefit))
    isMatch = kTRUE;
  
  if(isTPC){
    Double_t vecTrackingEff[5] = { isMatch,esdTrack->Phi(), esdTrack->Pt(),esdTrack->Eta(),esdTrack->GetITSclusters(0) };
    fTrackingEffHisto->Fill(vecTrackingEff);
  }
}

//_____________________________________________________________________________
/*void AliPerformanceMatch::ProcessTPCTRD(AliStack* , AliESDtrack *const esdTrack, AliESDfriendTrack *const esdFriendTrack)
{
  return;
}*/

//_____________________________________________________________________________
void AliPerformanceMatch::FillHistograms(AliESDtrack *const refParam, AliESDtrack *const param, Bool_t isRec) 
{
  //
  // fill performance histograms 
  // (TPC always as reference)
  //

  if(!refParam) return;
  if(!param) return;
  if(!isRec) return;
  
  //
  // Deltas (dy,dz,dphi,dtheta,dpt)
  //
  Float_t delta[5] = {0};
  if(param && isRec) {
    delta[0] = param->GetY()-refParam->GetY();
    delta[1] = param->GetZ()-refParam->GetZ();
    delta[2] = TMath::ATan2(param->Py(),param->Px())-TMath::ATan2(refParam->Py(),refParam->Px());
    delta[3] = TMath::ATan2(param->Pz(),param->Pt())-TMath::ATan2(refParam->Pz(),refParam->Pt());
    if(refParam->Pt()) delta[4] = (param->Pt()-refParam->Pt())/refParam->Pt();
  }
  // 
  // Pulls (y,z,snp,tanl,1/pt)
  //
  Float_t sigma[5] = {0};
  Float_t pull[5] = {0};
  if(param && isRec) {
    sigma[0] = TMath::Sqrt(param->GetSigmaY2()+refParam->GetSigmaY2());   
    sigma[1] = TMath::Sqrt(param->GetSigmaZ2()+refParam->GetSigmaZ2());
    sigma[2] = TMath::Sqrt(param->GetSigmaSnp2()+refParam->GetSigmaSnp2());
    sigma[3] = TMath::Sqrt(param->GetSigmaTgl2()+refParam->GetSigmaTgl2());
    sigma[4] = TMath::Sqrt(param->GetSigma1Pt2()+refParam->GetSigma1Pt2());
    if(sigma[0]) pull[0] = delta[0] / sigma[0]; 
    if(sigma[1]) pull[1] = delta[1] / sigma[1]; 
    if(sigma[2]) pull[2] = (param->GetSnp()-refParam->GetSnp()) / sigma[2]; 
    if(sigma[3]) pull[3] = (param->GetTgl()-refParam->GetTgl()) / sigma[3]; 
    if(sigma[4]) pull[4] = (param->OneOverPt()-refParam->OneOverPt()) / sigma[4]; 
  }

  // Fill histograms
  Double_t vResolHisto[9] = {delta[0],delta[1],delta[2],delta[3],delta[4],refParam->Phi(),refParam->Eta(),refParam->Pt(),isRec};
  if(fabs(pull[4])<5)
    fResolHisto->Fill(vResolHisto);

  Double_t vPullHisto[9] = {pull[0],pull[1],pull[2],pull[3],pull[4],refParam->Phi(),refParam->Eta(),refParam->OneOverPt(),isRec};
  if(fabs(pull[4])<5)
    fPullHisto->Fill(vPullHisto);
}

//_____________________________________________________________________________
void AliPerformanceMatch::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend */*const esdFriend*/, const Bool_t bUseMC, const Bool_t /*bUseESDfriend*/)
{
  // Process comparison information 
  //
  if(!esdEvent) 
  {
    Error("Exec","esdEvent not available");
    return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  
  if(bUseMC)
  {
    if(!mcEvent) {
      Error("Exec","mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      Error("Exec","Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      Error("Exec","Stack not available");
      return;
    }
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      Error("Exec","Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);
  } 
  
  // use ESD friends
  /*  if(bUseESDfriend) {
    if(!esdFriend) {
      Error("Exec","esdFriend not available");
      return;
    }
    }*/

  // trigger
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = esdEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) return; 
  }

  // get TPC event vertex
  const AliESDVertex *vtxESD = esdEvent->GetPrimaryVertexTPC();
  if(vtxESD && (vtxESD->GetStatus()<=0)) return;

  //  Process events
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if(!track) continue;

    /*AliESDfriendTrack *friendTrack=0;
    if(bUseESDfriend) {
      friendTrack=esdFriend->GetTrack(iTrack);
      if(!friendTrack) continue;
      }*/

    if(GetAnalysisMode() == 0){
      if(!IsUseTOFBunchCrossing())
	ProcessTPCITS(stack,track);
      else
	if( track->GetTOFBunchCrossing(esdEvent->GetMagneticField())==0) 
	  ProcessTPCITS(stack,track);
    }
    /* else if(GetAnalysisMode() == 2) ProcessTPCTRD(stack,track,friendTrack);*/
    else if(GetAnalysisMode() == 1) {ProcessITSTPC(iTrack,esdEvent,stack,track);
    }
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }

}

//_____________________________________________________________________________
TH1F* AliPerformanceMatch::MakeResol(TH2F * his, Int_t integ, Bool_t type, Int_t cut){
  // Create resolution histograms
 
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoII(his,&hism,integ,kTRUE,cut);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliPerformanceMatch::Analyse() {
  // Analyse comparison information and store output histograms
  // in the folder "folderMatch"
  //
  TString selString;
  /*
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH1F *h2=0;
  TH2F *h2D=0;
  */
  TObjArray *aFolderObj = new TObjArray;

  // write results in the folder 
  // TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  // c->cd();

  // char name[256];
  // char title[256];

  if(GetAnalysisMode()==1 || GetAnalysisMode()==2) { 

  fResolHisto->GetAxis(8)->SetRangeUser(1.0,2.0); // only reconstructed
  fPullHisto->GetAxis(8)->SetRangeUser(1.0,2.0);  // only reconstructed
  for(Int_t i=0; i<5; i++) 
    {
      for(Int_t j=5; j<8; j++) 
	{
	  //if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
	  if(j!=6) fResolHisto->GetAxis(6)->SetRangeUser(0.0,1.49); // eta window
	  else fResolHisto->GetAxis(6)->SetRangeUser(-1.5,1.49);
	  fResolHisto->GetAxis(7)->SetRangeUser(0.01,10.); // pt threshold
	  
	  selString = "resol";
	  AddProjection(aFolderObj, "match", fResolHisto, i, j, &selString);
	  
	  
	  if(j!=6) fPullHisto->GetAxis(6)->SetRangeUser(0.0,1.49); // eta window
	  else  fPullHisto->GetAxis(6)->SetRangeUser(-1.5,1.49); // eta window
	  fPullHisto->GetAxis(7)->SetRangeUser(0.01,10.);  // pt threshold
	  selString = "pull";
	  AddProjection(aFolderObj, "match", fPullHisto, i, j, &selString);
	  
	}
    }
  }
  // 
  // TPC efficiency wrt ITS
  //
  if(GetAnalysisMode()==0) { 
    selString = "trackingeff";
    AddProjection(aFolderObj, "match", fTrackingEffHisto, 0, &selString);
    
    for(Int_t i=1; i<5; i++) 
      {
	// all ITS standalone tracks
	fTrackingEffHisto->GetAxis(0)->SetRange(1,fTrackingEffHisto->GetAxis(0)->GetNbins());
	selString = "trackingeff_all";
	AddProjection(aFolderObj, "match", fTrackingEffHisto, i, 3,&selString);
	
      // TPC tracks which has matching with TPC
	fTrackingEffHisto->GetAxis(0)->SetRange(2,2);
	selString = "trackingeff_tpc";
	AddProjection(aFolderObj, "match", fTrackingEffHisto, i, 3,&selString);
      }
  }

  printf("exportToFolder\n");
  fAnalysisFolder = ExportToFolder(aFolderObj);
  
  // delete only TObjArray
  if(fFolderObj) delete fFolderObj;
  fFolderObj = aFolderObj;  
  aFolderObj=0;  
  
  }


//_____________________________________________________________________________
TFolder* AliPerformanceMatch::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceMatch * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
       newFolder->Add(array->At(i));
       i++;
     }
  }
  
  return newFolder;
}
 
//_____________________________________________________________________________
TFolder* AliPerformanceMatch::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}

//_____________________________________________________________________________
Long64_t AliPerformanceMatch::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;
  
  Bool_t merge = ((fgUseMergeTHnSparse && fgMergeTHnSparse) || (!fgUseMergeTHnSparse && fMergeTHnSparseObj));

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  TObjArray* objArrayList = 0;
  objArrayList = new TObjArray();

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliPerformanceMatch* entry = dynamic_cast<AliPerformanceMatch*>(obj);
    if (entry == 0) continue; 
    if (merge) {
        if ((fResolHisto) && (entry->fResolHisto)) { fResolHisto->Add(entry->fResolHisto); }
        if ((fPullHisto) && (entry->fPullHisto)) { fPullHisto->Add(entry->fPullHisto); }
        if ((fTrackingEffHisto) && (entry->fTrackingEffHisto)) { fTrackingEffHisto->Add(entry->fTrackingEffHisto); }
    }
    // the analysisfolder is only merged if present
    if (entry->fFolderObj) { objArrayList->Add(entry->fFolderObj); }

    count++;
  }
  if (fFolderObj) { fFolderObj->Merge(objArrayList); } 
  // to signal that track histos were not merged: reset
  if (!merge) { fResolHisto->Reset(); fPullHisto->Reset(); fTrackingEffHisto->Reset(); }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;
return count;
}
