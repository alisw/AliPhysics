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

#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include "AliAODHeader.h"
//#include "AliTriggerAnalysis.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliGenPythiaEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"

#include "AliAnalysisTaskCounter.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCounter) ;
/// \endcond

//______________________________________________________________
/// Constructor.
//______________________________________________________________
AliAnalysisTaskCounter::AliAnalysisTaskCounter(const char *name) 
: AliAnalysisTaskSE(name), 
  fAcceptFastCluster(kTRUE),
  fZVertexCut(10.), 
  fTrackMultEtaCut(0.8),
  fAvgTrials(-1),
  fOutputContainer(0x0),
  fESDtrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
  //fTriggerAnalysis (new AliTriggerAnalysis),
  fCurrFileName(0), fCheckMCCrossSection(kFALSE),
  fUseAliCentrality(kFALSE), fCentralityClass("V0M"),
  fhNEvents(0),
  fhXVertex(0),    fhYVertex(0),    fhZVertex(0),
  fhXGoodVertex(0),fhYGoodVertex(0),fhZGoodVertex(0),
  fhCentrality(0), fhEventPlaneAngle(0),
  fh1Xsec(0),      fh1Trials(0)
{
  DefineOutput(1, TList::Class());
}

//______________________________________________
/// Default Constructor.
//_______________________________________________
AliAnalysisTaskCounter::AliAnalysisTaskCounter() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskCounter"),
    fAcceptFastCluster(kTRUE),
    fZVertexCut(10.),
    fTrackMultEtaCut(0.8),
    fAvgTrials(-1),
    fOutputContainer(0x0),
    fESDtrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
    //fTriggerAnalysis (new AliTriggerAnalysis),
    fCurrFileName(0), fCheckMCCrossSection(kFALSE),
    fUseAliCentrality(kFALSE), fCentralityClass("V0M"),
    fhNEvents(0),
    fhXVertex(0),    fhYVertex(0),    fhZVertex(0),
    fhXGoodVertex(0),fhYGoodVertex(0),fhZGoodVertex(0),
    fhCentrality(0), fhEventPlaneAngle(0),
    fh1Xsec(0),      fh1Trials(0)
{
}

//_______________________________________________
/// Destructor.
//_______________________________________________
AliAnalysisTaskCounter::~AliAnalysisTaskCounter()
{  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  
  if(fOutputContainer)
  {
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
  
  if(fESDtrackCuts)    delete fESDtrackCuts;
  //if(fTriggerAnalysis) delete fTriggerAnalysis;
  
}


//____________________________________________________
/// Init histogram pointers and add them to container.
//____________________________________________________
void AliAnalysisTaskCounter::UserCreateOutputObjects()
{  
  fOutputContainer = new TList();
  
  if(fCheckMCCrossSection)
  {
    fh1Xsec = new TH1F("hXsec","xsec from pyxsec.root",1,0,1);
    fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
    fOutputContainer->Add(fh1Xsec);
    
    fh1Trials = new TH1F("hTrials","trials root file",1,0,1);
    fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    fOutputContainer->Add(fh1Trials);
  }
  
  fhZVertex     = new TH1F("hZVertex", " Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZVertex);

  fhZGoodVertex     = new TH1F("hZGoodVertex", " Good Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZGoodVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZGoodVertex);
  
  fhXVertex     = new TH1F("hXVertex", " X vertex distribution"   , 200 , -2 , 2  ) ;
  fhXVertex->SetXTitle("v_{x} (cm)");
  fOutputContainer->Add(fhXVertex);
  
  fhXGoodVertex     = new TH1F("hXGoodVertex", " Good X vertex distribution"   , 200 , -2 , 2  ) ;
  fhXGoodVertex->SetXTitle("v_{x} (cm)");
  fOutputContainer->Add(fhXGoodVertex);
  
  fhYVertex     = new TH1F("hYVertex", " Y vertex distribution"   , 200 , -2 , 2  ) ;
  fhYVertex->SetXTitle("v_{y} (cm)");
  fOutputContainer->Add(fhYVertex);
  
  fhYGoodVertex     = new TH1F("hYGoodVertex", " Good Y vertex distribution"   , 200 , -2 , 2  ) ;
  fhYGoodVertex->SetXTitle("v_{y} (cm)");
  fOutputContainer->Add(fhYGoodVertex);
  
  fhCentrality   = new TH1F("hCentrality","Number of events in centrality bin, |vz|<10 cm, method <V0M> ",100,0.,100.) ;
  fhCentrality->SetXTitle("Centrality bin");
  fOutputContainer->Add(fhCentrality) ;  
  
  fhEventPlaneAngle=new TH1F("hEventPlaneAngle","Number of events in event plane, |vz|<10 cm, method <V0> ",100,0.,TMath::Pi()) ;
  fhEventPlaneAngle->SetXTitle("EP angle (rad)");
  fOutputContainer->Add(fhEventPlaneAngle) ;
  
  fhNEvents = new TH1I("hNEvents", "Number of analyzed events", 21, 0, 21) ;
  //fhNEvents->SetXTitle("Selection");
  fhNEvents->SetYTitle("# events");
  fhNEvents->GetXaxis()->SetBinLabel(1 ,"1  = PS");
  fhNEvents->GetXaxis()->SetBinLabel(2 ,"2  = 1  & ESD");
  fhNEvents->GetXaxis()->SetBinLabel(3 ,"3  = 2  & |Z|<10");
  fhNEvents->GetXaxis()->SetBinLabel(4 ,"4  = 2  & !track?");
  fhNEvents->GetXaxis()->SetBinLabel(5 ,"5  = 3  & 4");
  fhNEvents->GetXaxis()->SetBinLabel(6 ,"6  = 2  & V0AND");
  fhNEvents->GetXaxis()->SetBinLabel(7 ,"7  = 3  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(8 ,"8  = 4  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(9 ,"9  = 5  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(10,"10 = 2  & not pileup");
  fhNEvents->GetXaxis()->SetBinLabel(11,"11 = 2  & good vertex");
  fhNEvents->GetXaxis()->SetBinLabel(12,"12 = 3  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(13,"13 = 4  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(14,"14 = 6  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(15,"15 = 9  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(16,"16 = 10 & 11");
  fhNEvents->GetXaxis()->SetBinLabel(17,"17 = 6  & 10");
  fhNEvents->GetXaxis()->SetBinLabel(18,"18 = Reject EMCAL 1");
  fhNEvents->GetXaxis()->SetBinLabel(19,"19 = 18 & 3");
  fhNEvents->GetXaxis()->SetBinLabel(20,"20 = Reject EMCAL 2");
  fhNEvents->GetXaxis()->SetBinLabel(21,"21 = 20 & 3");

  fOutputContainer->Add(fhNEvents);

  fOutputContainer->SetOwner(kTRUE);
  
  PostData(1,fOutputContainer);
  
}

//_______________________________________________
/// Main event loop
/// It does the event counting, depending on different cuts
/// (see criteria in class description). 
/// It fills the event vertex, centrality and plane histograms 
//_______________________________________________
void AliAnalysisTaskCounter::UserExec(Option_t *) 
{
  //printf("___ Event __ %d __\n",(Int_t)Entry());
  
  Notify();
  
  fhNEvents->Fill(0.5);  
  
  AliVEvent * event = InputEvent();
  AliESDEvent * esdevent = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent * aodevent = dynamic_cast<AliAODEvent*> (event);
  
  // Init mag field for tracks in case of ESDs, needed, not clear why
  if (!TGeoGlobalMagField::Instance()->GetField() && esdevent) esdevent->InitMagneticField();

  TString triggerclasses = event->GetFiredTriggerClasses();

  //printf("Trigger class fired: %s \n",event->GetFiredTriggerClasses().Data());
  
  if (triggerclasses.Contains("FAST") && !triggerclasses.Contains("ALL") && !fAcceptFastCluster) 
  {
    //printf("Do not count events from fast cluster, trigger name %s\n",triggerclasses.Data());
    return;
  }

  fhNEvents->Fill(1.5);  
    
  //Initialize bools
  Bool_t bSelectVZ    = kFALSE;
  Bool_t bV0AND       = kFALSE; 
  Bool_t bPileup      = kFALSE;
  Bool_t bGoodV       = kFALSE;
  Bool_t bSelectTrack = kFALSE;   
  Int_t  trackMult    = 0;
  
  //---------------------------------
  //Get the primary vertex, cut on Z
  //---------------------------------
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v) ;
  fhXVertex->Fill(v[0]);
  fhYVertex->Fill(v[1]);
  fhZVertex->Fill(v[2]);
  
  if(TMath::Abs(v[2]) < fZVertexCut) 
  {
    bSelectVZ=kTRUE;
    fhNEvents->Fill(2.5);  
  }
  //else printf("Vertex out %f \n",v[2]);
  
  //--------------------------------------------------
  //Count tracks, cut on number of tracks in eta < 0.8
  //--------------------------------------------------
  Int_t nTracks   = event->GetNumberOfTracks() ;
  for (Int_t itrack =  0; itrack <  nTracks; itrack++)
  {////////////// track loop
    AliVTrack * track = (AliVTrack*)event->GetTrack(itrack) ; // retrieve track from esd
    
    //ESDs
    if(esdevent && !fESDtrackCuts->AcceptTrack((AliESDtrack*)track)) continue;
    
    //AODs
    if(aodevent && !((AliAODTrack*)track)->IsHybridGlobalConstrainedGlobal()) continue ;
    
    //Do not count tracks out of acceptance cut
    if(TMath::Abs(track->Eta())< fTrackMultEtaCut) trackMult++;
  }
  
  //printf("AliAnalysisTaskCounter::UserExec() - Track Mult %d \n",trackMult);
  
  //--------------------------------------------------
  // At least one track
  //--------------------------------------------------
  if (trackMult > 0) 
  {
    bSelectTrack = kTRUE; 
                  fhNEvents->Fill(3.5);
    if(bSelectVZ) fhNEvents->Fill(4.5);
  }
  
  //---------------------------------
  // V0AND
  //---------------------------------
  
  //if(esdevent) bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esdevent, AliTriggerAnalysis::kV0AND);
  AliVVZERO* v0 = fInputEvent->GetVZEROData();
  bV0AND = ((v0->GetV0ADecision()==1) && (v0->GetV0CDecision()==1));
  
  if(bV0AND)
  {
                                    fhNEvents->Fill(5.5);
    if (bSelectVZ)                  fhNEvents->Fill(6.5);
    if (bSelectTrack)               fhNEvents->Fill(7.5);
    if (bSelectVZ &&  bSelectTrack) fhNEvents->Fill(8.5);
  }
  
  //---------------------------------
  // Pileup
  //---------------------------------
  bPileup = event->IsPileupFromSPD(3, 0.8, 3., 2., 5.); //Default values, if not it does not compile
  //bPileup = event->IsPileupFromSPD();
  
  if (!bPileup)
  {
                fhNEvents->Fill(9.5);
    if(bV0AND)  fhNEvents->Fill(16.5);
  }
  
  //---------------------------------
  // Good vertex
  //---------------------------------
  bGoodV = CheckForPrimaryVertex();
  
  //Remove events with  vertex (0,0,0), bad vertex reconstruction
  if(TMath::Abs(v[0]) < 1.e-6 && 
     TMath::Abs(v[1]) < 1.e-6 && 
     TMath::Abs(v[2]) < 1.e-6) bGoodV = kFALSE;

  if(bGoodV) 
  {
    fhXGoodVertex->Fill(v[0]);
    fhYGoodVertex->Fill(v[1]);
    fhZGoodVertex->Fill(v[2]);
    
                     fhNEvents->Fill(10.5);
    if(bSelectVZ)    fhNEvents->Fill(11.5);
    if(bSelectTrack) fhNEvents->Fill(12.5);
    if(bV0AND)       fhNEvents->Fill(13.5);
    if(bSelectVZ && bSelectTrack && bV0AND)    
                     fhNEvents->Fill(14.5); 
    if(!bPileup)     fhNEvents->Fill(15.5); 

    if(TMath::Abs(v[2]) < 10.) 
    {
      if(fUseAliCentrality)
      {
        if(InputEvent()->GetCentrality() && fUseAliCentrality) 
          fhCentrality->Fill(InputEvent()->GetCentrality()->GetCentralityPercentile(fCentralityClass)) ;
      }
      else
      {
        AliMultSelection* multSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection") ;
        if(multSelection) fhCentrality->Fill(multSelection->GetMultiplicityPercentile(fCentralityClass, kTRUE));
      }
        
      if(InputEvent()->GetEventplane()) 
      {
        Float_t ep = InputEvent()->GetEventplane()->GetEventplane("V0", InputEvent());
      
        ep+=TMath::Pi()/2.; // put same range as for <Q> method, [0,pi]
        
        fhEventPlaneAngle->Fill(ep);
      }
    }
  
  }

  //printf("AliAnalysisTaskCounter::UserExec() : z vertex %d, good vertex %d, v0and %d, pile up %d, track mult %d\n ", bSelectVZ, bGoodV, bV0AND, bPileup, trackMult);
  
  // Events that could be rejected in EMCAL
  // LHC11a, SM4 and some SM3 events cut with this
  Bool_t bEMCALRejected = kFALSE;
  for (Int_t i = 0; i < InputEvent()->GetNumberOfCaloClusters(); i++)
  {
    AliVCluster *clus = InputEvent()->GetCaloCluster(i);
    if(clus->IsEMCAL())
    {    
      if ((clus->E() > 500 && clus->GetNCells() > 200 ) || clus->GetNCells() > 200)  
      {
        
        //printf("Counter: Reject event with cluster: E %f, ncells %d\n",clus->E(),clus->GetNCells());
        
                         fhNEvents->Fill(17.5); 
        if(bSelectVZ)    fhNEvents->Fill(18.5);
        bEMCALRejected = kTRUE;
        break;
      }
    }
  }
  
  //LHC11a, 3 last runs, cut with this
  if(!bEMCALRejected)
  {
    // Count number of cells in SM3 with energy larger than 0.1, cut on this number
    Int_t ncellsSM3 = 0;
    Int_t ncellsSM4 = 0;
    for(Int_t icell = 0; icell < event->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      if(event->GetEMCALCells()->GetAmplitude(icell) > 0.1 && event->GetEMCALCells()->GetCellNumber(icell)/(24*48)==3) ncellsSM3++;
      if(event->GetEMCALCells()->GetAmplitude(icell) > 0.1 && event->GetEMCALCells()->GetCellNumber(icell)/(24*48)==4) ncellsSM4++;
    }
    
    Int_t ncellcut = 21;
    if(triggerclasses.Contains("EMC")) ncellcut = 35;
    
    if( ncellsSM3 >= ncellcut || ncellsSM4 >= 100 )
    {
      //printf("Counter: reject event with ncells in SM3: ncells %d\n",ncells);

                       fhNEvents->Fill(19.5); 
      if(bSelectVZ)    fhNEvents->Fill(20.5);
    }
    
  }
  
  PostData(1,fOutputContainer);

}

//____________________________________________________
/// Check if the vertex was well reconstructed, 
/// copy of PCM.
//____________________________________________________
Bool_t AliAnalysisTaskCounter::CheckForPrimaryVertex()
{
  
  AliESDEvent * esdevent = dynamic_cast<AliESDEvent*> (InputEvent());
  AliAODEvent * aodevent = dynamic_cast<AliAODEvent*> (InputEvent());
  
  if(esdevent)
  {
    if(esdevent->GetPrimaryVertex()->GetNContributors() > 0)
    {
      return kTRUE;
    }
    
    if(esdevent->GetPrimaryVertex()->GetNContributors() < 1)
    {
      // SPD vertex
      if(esdevent->GetPrimaryVertexSPD()->GetNContributors() > 0)
      {
        return kTRUE;
        
      }
      if(esdevent->GetPrimaryVertexSPD()->GetNContributors() < 1)
      {
        return kFALSE;
      }
    }
  }
  else if(aodevent)
  {    
    if (aodevent->GetPrimaryVertex() != NULL)
    {
      if(aodevent->GetPrimaryVertex()->GetNContributors() > 0)
      {
        return kTRUE;
      }
    }
    
    if(aodevent->GetPrimaryVertexSPD() != NULL)
    {
      if(aodevent->GetPrimaryVertexSPD()->GetNContributors() > 0)
      {
        return kTRUE;
      }
      else
      {
        AliWarning(Form("Number of contributors from bad vertex type:: %s",aodevent->GetPrimaryVertex()->GetName()));
        return kFALSE;
      }
    }
  }
  else return kTRUE;
  
  return kFALSE;

}

//_____________________________________________
/// Put in the output some event summary histograms
//_____________________________________________
void AliAnalysisTaskCounter::FinishTaskOutput()
{  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputH = dynamic_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!inputH) return; 
  TH2F *histStat = dynamic_cast<TH2F*>(inputH->GetStatistics()); 
  TH2F *histBin0 = dynamic_cast<TH2F*>(inputH->GetStatistics("BIN0"));
  
  if(histStat)
    fOutputContainer->Add(histStat);
  else AliInfo("Stat histogram not available check, \n if ESDs, that AliPhysicsSelection was on, \n if AODs, if EventStat_temp.root exists");

  if(histBin0)
    fOutputContainer->Add(histBin0); 
  
}

//_____________________________________
/// Implemented Notify() to read the cross sections
/// and number of trials from pyxsec.root, values stored
/// in specific histograms.
//_____________________________________
Bool_t AliAnalysisTaskCounter::Notify()
{
  if(!fCheckMCCrossSection) return kTRUE;

  // Fetch the aod also from the input in,
  // have todo it in notify
  
  Float_t xsection = 0;
  Float_t trials   = 1;
  fAvgTrials = -1;
  
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if(!tree) return kFALSE;
  
  TFile *curfile = tree->GetCurrentFile();
  
  if(!curfile) return kFALSE;
  
  if(fCurrFileName == curfile->GetName()) return kFALSE;
  
  fCurrFileName = TString(curfile->GetName());
  
  if(!fh1Xsec||!fh1Trials)
  {
    AliInfo(Form("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__));
    return kFALSE;
  }

  
  Bool_t ok = PythiaInfoFromFile(fCurrFileName,xsection,trials);


    if(!ok || xsection==0){

    // Try to get the information for the header (AOD) if the pysec histo are not filled
       AliGenPythiaEventHeader *fPythiaHeader=0;
        
        AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
 
          if (aodMCH) {
            for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
            fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
              if (fPythiaHeader) break;
            }
	  }         
 
 if(!fPythiaHeader){ 
   AliError(Form("No pythia header found"));
   return kFALSE;
}

  Float_t pthard=0;

  pthard = fPythiaHeader->GetPtHard();
  xsection = fPythiaHeader->GetXsection();
  trials = fPythiaHeader->Trials();

  fh1Xsec->Fill("<#sigma>",xsection);
  fh1Trials->Fill("#sum{ntrials}",trials);
  
  AliInfo(Form("xs %f, trial %f, pt hard %f\n",xsection,trials, pthard));
 
  return kTRUE;
    }

    

  fh1Xsec->Fill("<#sigma>",xsection);
  
  // average trials per event
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  
  if(trials >= nEntries && nEntries > 0.) fAvgTrials = trials/nEntries;
  
  fh1Trials->Fill("#sum{ntrials}",trials);
  
  AliInfo(Form("xs %f, trial %f, avg trials %f\n",xsection,trials, fAvgTrials));
  
  AliDebug(1,Form("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName()));
  
  
 return kTRUE;

}

//_____________________________________________________________________________________________
/// This method gets and returns the 
/// \param file : either pyxsec.root (ESDs) or pysec_hists.root (AODs) files
/// \param xsec : cross section 
/// \param trials : number of event trials 
/// that should be located where the main data file is.
/// This is called in Notify and should provide the path to the AOD/ESD file
//_____________________________________________________________________________________________
Bool_t AliAnalysisTaskCounter::PythiaInfoFromFile(TString file,Float_t & xsec,Float_t & trials)
{
   xsec   = 0;
  trials = 1;
  
  if(file.Contains("root_archive.zip#"))
  {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else
  {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec)
  {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec)
    {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else
    {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if(!key)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else
  {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree)
    {
      fxsec->Close();
      return kFALSE;
    }
    
    UInt_t   ntrials  = 0;
    Double_t xsection = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }
   
  return kTRUE;
}

