/*************************************************************************
* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskBuildCorrTree.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliVVZERO.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "THashList.h"
#include "TMath.h"
#include <AliFilteredTrack.h>
#include <AliEventplane.h>
#include "AliAODMCParticle.h"
// #include<iostream>
#include <sstream>

ClassImp(AliAnalysisTaskBuildCorrTree)
//
//Task to create three particle correlations.
//Authors:
// Matthias Richter
// Paul Baetzing || pbatzing@cern.ch
//


AliAnalysisTaskBuildCorrTree::AliAnalysisTaskBuildCorrTree()
  : AliAnalysisTaskSE("AliAnalysisTaskCorrelation3p")
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption("")
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fVertex()
  , fCollisionType(AliAnalysisTaskBuildCorrTree::PbPb)
  , fperiod(AliAnalysisTaskBuildCorrTree::P11h)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fMcArray(NULL)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityEstimator("V0M")
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fRunNumberList(NULL)
  , fNruns(fNRunsP11h)
  , fRun(0)
  , fRunFillValue(0.0)
  , fEvent(0x0)
  , fTree(0x0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinPt(3.0)
  , fMaxPt(8.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
{
  // default constructor
  // 
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-5,-2.5,2.5,5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
}

AliAnalysisTaskBuildCorrTree::AliAnalysisTaskBuildCorrTree(const char *name, const char* opt)
  : AliAnalysisTaskSE(name)
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption(opt)
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fVertex()
  , fCollisionType(AliAnalysisTaskBuildCorrTree::PbPb)
  , fperiod(AliAnalysisTaskBuildCorrTree::P11h)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fMcArray(NULL)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityEstimator("V0M")
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fRunNumberList(NULL)
  , fNruns(fNRunsP11h)
  , fRun(0)
  , fRunFillValue(0.0)
  , fEvent(0x0)
  , fTree(0x0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinPt(3.0)
  , fMaxPt(8.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
  {
  // constructor with options
  //
  //
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-7.5,-2.5,2.5,7.5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
  }

int AliAnalysisTaskBuildCorrTree::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());   // reduced information tree
  return 0;
}

AliAnalysisTaskBuildCorrTree::~AliAnalysisTaskBuildCorrTree()
{
  // destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliAnalysisTaskBuildCorrTree::UserCreateOutputObjects()
{
  // create result objects and add to output list
  TH1::SetDefaultSumw2(kTRUE);//want the collection of weights on all histograms.
  fOutput = new THashList;
  fOutput->SetOwner();
    if(fCollisionType==pp)for(int i=0;i<fMBinEdges.GetSize();i++){
      fMBinEdges.AddAt(fMaxNumberOfTracksInPPConsidered*fMBinEdges.At(i)/fMBinEdges.At(fMBinEdges.GetSize()-1),i);
    }
  
  InitializeQAhistograms();
  OpenFile(2);
  fTree = new TTree("DstTree","Reduced ESD information");
  fEvent = new AliFilteredEvent("DstEvent");
  fTree->Branch("Event",&fEvent,16000,99);
  
  // all tasks must post data once for all outputs
  PostData(1, fOutput);
  PostData(2, fTree);
}

void AliAnalysisTaskBuildCorrTree::UserExec(Option_t* /*option*/)
{
  FillHistogram("Check",0.0);
  fEvent->ClearEvent();
  // process the event
  TObject* pInput=InputEvent();
  if (!pInput) {AliError("failed to get input");return;}
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  if(!pEvent){AliError(Form("input of wrong class type %s, expecting AliVEvent", pInput->ClassName()));return;}

  
  //Find out if it is AOD or ESD.  
  fisESD=pEvent->IsA()==AliESDEvent::Class();
  fisAOD=pEvent->IsA()==AliAODEvent::Class();
  if(fisESD)return;//ESD analysis not implemented
  fRun = pEvent->GetRunNumber();
  fEvent->SetRunNr(fRun);
  //Get the MC array:
  GetMCArray();
  
  
  TAxis* runnumberaxis= dynamic_cast<TH1D*>(fOutput->FindObject("EventsperRun"))->GetXaxis();
  if (runnumberaxis){double RunBin = runnumberaxis->FindBin(Form("%i",fRun));fRunFillValue = runnumberaxis->GetBinCenter(RunBin);}   
  
  //Get the runnumber and find which bin and fill value this corresponds to.
  GetCentralityAndVertex();
  if(!SelectEvent()) return;//events are rejected.
  // ep angle interval [todo, fill]
  AliEventplane *ep = pEvent->GetEventplane();
  if(ep){
    Double_t qx = 0,qy = 0;
    fEvent->SetEventPlane(ep->CalculateVZEROEventPlane(pEvent,10,2,qx,qy));
  }
  
  if(fCollisionType==AliAnalysisTaskBuildCorrTree::pp){ FillHistogram("centVsZVertex",fMultiplicity,fVertex[2]);fEvent->SetCentrality(fMultiplicity);fEvent->SetVertex(fVertex);}//only fill with selected events.
  if(fCollisionType==AliAnalysisTaskBuildCorrTree::PbPb){FillHistogram("centVsZVertex",fCentralityPercentile,fVertex[2]);fEvent->SetCentrality(fCentralityPercentile);fEvent->SetVertex(fVertex);}
  FillHistogram("EventsperRun", fRunFillValue);
  FillHistogram("NEventsVertex",fRunFillValue,fVertex[2]);
  FillHistogram("NEventsCent",fRunFillValue,fCentralityPercentile);
  FillHistogram("centVsNofTracks",fCentralityPercentile, GetTracks(pEvent));
  fTree->Fill();
  PostData(1, fOutput);  
  PostData(2, fTree);
}

void AliAnalysisTaskBuildCorrTree::FinishTaskOutput()
{
  // end of the processing
    TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject("trackCount")) ;
    if (hist) AliWarning(Form("FinishTaskOutput: %i events(s)" ,(int)hist->GetEntries()));
}

void AliAnalysisTaskBuildCorrTree::Terminate(Option_t *)
{
  // last action on the client
//   gObjectTable->Print();

}

Int_t AliAnalysisTaskBuildCorrTree::GetTracks(AliVEvent *pEvent)
{
  Int_t nofTracks = 0;
  nofTracks=pEvent->GetNumberOfTracks();
  FillHistogram("trackCount",nofTracks);
  for (int i=0; i<nofTracks; i++) {
    AliVParticle* t=pEvent->GetTrack(i);
    if (!t) continue;
    FillHistogram("TracksperRun",fRunFillValue);
    FillHistogram("trackUnselectedPt",t->Pt());
    FillHistogram("trackUnselectedPhi",t->Phi());
    FillHistogram("trackUnselectedTheta",t->Theta());
    if (!IsSelected(t)) continue;

    TClonesArray& tracks = *(fEvent->GetTracks());
    AliAODTrack *AODt = dynamic_cast<AliAODTrack*>(t);
    AliFilteredTrack *reducedParticle=new(tracks[fEvent->GetNtrks()]) AliFilteredTrack(*AODt);
    reducedParticle->SetMC(false);
    if(fCollisionType==AliAnalysisTaskBuildCorrTree::PbPb){
      FillHistogram("selectedTracksperRun",fRunFillValue);
      FillHistogram("NTracksVertexEta",fRunFillValue,fVertex[2],t->Eta());
      FillHistogram("NTracksCent",fRunFillValue,fCentralityPercentile);      
      FillHistogram("NTracksPhi",fRunFillValue,t->Phi());
      FillHistogram("NTrackspT",fRunFillValue,t->Pt());
    }
    
    FillHistogram("trackPt",t->Pt());
    FillHistogram("trackPhi",t->Phi());
    FillHistogram("trackTheta",t->Theta());
    fEvent->SetNtrks(fEvent->GetNtrks()+1);
  }
//   AliWarning("before if MCarray");
  if(fMcArray){
    int nofMCParticles = fMcArray->GetEntriesFast();
    FillHistogram("trackCount",nofMCParticles);
    for (int i=0;i<nofMCParticles;i++){
      AliVParticle* t =  (AliVParticle *) fMcArray->At(i);
      if (!t) continue;
//       AliWarning("before selection");
      if( t->Charge()!=0){//check if they are physical primary particles
	if (!IsSelected(t)) continue;
// 	AliWarning("MCTrack");
	FillHistogram("MCtrackPt",t->Pt());
	TClonesArray& tracks = *(fEvent->GetTracks());
	AliFilteredTrack *reducedParticle=new(tracks[fEvent->GetNtrks()]) AliFilteredTrack(*t);
	reducedParticle->SetMC(true);
	fEvent->SetNtrks(fEvent->GetNtrks()+1);
      }
    }
  }
  return nofTracks;
}




Bool_t AliAnalysisTaskBuildCorrTree::IsSelected(AliVParticle* p)
{
  //Performs selection cuts for tracks and triggers
  if (p->IsA()==AliESDtrack::Class() && IsSelectedTrackESD(p)) return IsSelectedTrack(p);
  if (p->IsA()==AliAODTrack::Class() && IsSelectedTrackAOD(p)) return IsSelectedTrack(p);
  if (p->IsA()==AliAODMCParticle::Class() && dynamic_cast<AliAODMCParticle*>(p)->IsPhysicalPrimary()) return IsSelectedTrack(p);
  return kFALSE;
}

Bool_t AliAnalysisTaskBuildCorrTree::IsSelectedTrack(AliVParticle* p)
{
  if (p->Pt()<=fMinPt) return kFALSE;
  if (fMaxPt>fMinPt && p->Pt()>fMaxPt) return kFALSE;
  float etatrigger=p->Eta();
  if (etatrigger<=-fAcceptancecut || etatrigger>=fAcceptancecut) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskBuildCorrTree::IsSelectedTrackAOD(AliVParticle* t)
{
  AliAODTrack *AODt = dynamic_cast<AliAODTrack*>(t);
  Bool_t isselected = kTRUE;
  Double_t DCAtang=-999.0;
  Double_t DCAlong=-999.0;
  GetDCA(DCAtang,DCAlong,AODt);
  if((AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2)))   FillHistogram("TrackDCAandonITS",DCAtang,DCAlong,1);
  else FillHistogram("TrackDCAandonITS",DCAtang,DCAlong,0);
  //Hybrid tracks give flat distributions
  if(AODt->IsHybridGlobalConstrainedGlobal()||AODt->TestFilterBit(BIT(4))||AODt->TestFilterBit(BIT(5))||AODt->TestFilterBit(BIT(6))) isselected = true;
  else isselected = false;
  if( (AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,1);
  if(!(AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,0);
  return isselected; 
}

Bool_t AliAnalysisTaskBuildCorrTree::IsSelectedTrackESD(AliVParticle* t)
{
  if(t)return kFALSE;//ESD is currently not supported
  else return kFALSE;
}

void AliAnalysisTaskBuildCorrTree::GetMCArray()
{//Gets the MCarray if we are in AOD.
    fMcArray = 0;
    AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (aodHandler){
      AliAODEvent *aod=aodHandler->GetEvent();
      if (aod ) {
	fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!fMcArray) AliError("Could not retrieve MC array!");
      }
      else AliError("Could not retrieve AOD event! MC is only supported in AOD.");
    }
}

void AliAnalysisTaskBuildCorrTree::GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt)
{
if(AODt->TestBit(AliAODTrack::kIsDCA)){
  DCAtang = AODt->DCA();
  DCAlong = AODt->ZAtDCA();
}
else{
  if(fVertexobj){
    Double_t fBzkg = dynamic_cast<AliAODEvent*>(InputEvent())->GetMagneticField();
    Double_t* dca = new Double_t[2];
    Double_t* dcacov = new Double_t[3];

    Double_t kVeryBigno = 1000000;
    if(AODt->PropagateToDCA(fVertexobj,fBzkg,kVeryBigno,dca,dcacov)){DCAtang=dca[0];DCAlong = dca[1];}
    else{DCAtang = -999;DCAlong=-999;}
    delete[] dca;
    delete[] dcacov;
    }
  }
} 


void AliAnalysisTaskBuildCorrTree::GetCentralityAndVertex()
{
  fVertexobj = InputEvent()->GetPrimaryVertex();
  if(fisAOD)
  {
    // Fill AliAODEvent interface specific information
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(dynamic_cast<AliAODEvent*>(InputEvent())->GetHeader());
    fCentrality =  dynamic_cast<AliCentrality*>(header->GetCentralityP());
    fMultiplicity = header->GetRefMultiplicity();
  }  
  if(fisESD)
  {
    fCentrality = dynamic_cast<AliESDEvent*>(InputEvent())->GetCentrality();
    fMultiplicity = 0.0 ;//not implemented atm
  }
  if(fCentrality)fCentralityPercentile = fCentrality->GetCentralityPercentile(fCentralityEstimator);

  if(fCentrality)FillHistogram("centrality",fCentralityPercentile);
  FillHistogram("multiplicity",fMultiplicity);
  //Get the primary Vertex
  if( fVertexobj ) {
    fVertex[0] = fVertexobj->GetX();
    fVertex[1] = fVertexobj->GetY();
    fVertex[2] = fVertexobj->GetZ();}
  else return;
}
Bool_t AliAnalysisTaskBuildCorrTree::SelectEvent()
{//This function provides the event cuts for this class.
  if(fCollisionType==pp){//With pp, the following cuts are applied:
    FillHistogram("Eventbeforeselection",fVertex[2],fMultiplicity,0);
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(TMath::Abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(fMultiplicity>fMaxNumberOfTracksInPPConsidered) return kFALSE;//Out of multiplicity bounds in pp, no histograms will be filled.
    if(InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.))return kFALSE;  //reject for pileup.
    FillHistogram("Eventafterselection",fVertex[2],fMultiplicity,0);
    
  }
  if(fCollisionType==PbPb){
    FillHistogram("Eventbeforeselection",fVertex[2],fMultiplicity,fCentralityPercentile);
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(TMath::Abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(!fCentrality){AliError("Centrality object not found.");return kFALSE;}//Centrality must be defined in the PbPb case.
    if(fCentrality->GetQuality()!=0)return kFALSE;//bad centrality.
    if(fCentralityPercentile<0) return kFALSE;//centrality is not defined
    if(fCentralityPercentile>fMaxMult) return kFALSE;//Out of centrality bounds in PbPb, will not fill any histogram.
    FillHistogram("Eventafterselection",fVertex[2],fMultiplicity,fCentralityPercentile);

  }
  return kTRUE;
}



void AliAnalysisTaskBuildCorrTree::InitializeQAhistograms()
{
  //Function that initializes the QA histograms 
  if (!fOutput) return;
  //QA histograms
  fOutput->Add(new TH1D("Check","Check",1,-0.5,0.5));

  fOutput->Add(new TH3D("TrackDCAandonITS","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers",50,-2,2,50,-5,5,2,-0.5,1.5));
  fOutput->Add(new TH3D("TrackDCAandonITSselected","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers for selected events",50,-2,2,50,-5,5,2,-0.5,1.5));
  fOutput->Add(new TH3D("Eventbeforeselection","Vertex vs Multiplicity vs Centrality before event selection.", 50,-15,15,50,0,12000,50,0,100));
  fOutput->Add(new TH3D("Eventafterselection","Vertex vs Multiplicity vs Centrality after event selection.", 50,-15,15,50,0,12000,50,0,100));
  fOutput->Add(new TH1D("trackCount", "trackCount", 1000,  0, 12000));
  fOutput->Add(new TH1D("trackUnselectedPt"   , "trackPt"   , 100,  0, 20));
  fOutput->Add(new TH1D("trackPt"   , "trackPt"   , 100,  0, 20));
  fOutput->Add(new TH1D("MCtrackPt"   , "trackPt"   , 100,  0, 20));
  fOutput->Add(new TH1D("trackUnselectedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackUnselectedTheta", "trackTheta",  180, 0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTheta", "trackTheta",  180, 0.0, TMath::Pi()));
//   fOutput->Add(new TH1D("Ntriggers","Number of triggers per event",50,-0.5,49.5));
  
  fOutput->Add(new TH1D("centrality", "Centrality",  100,  0, 100));
  fOutput->Add(new TH1D("multiplicity", "Multiplicity of tracklets",  100,  0, fMaxNumberOfTracksInPPConsidered));
  fOutput->Add(new TH2D("centVsNofTracks", "centVsNofTracks", 100, 0, 100, 100, 0, 12000));
  if(fCollisionType==PbPb)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, 100, 100, -10, 10));
  if(fCollisionType==pp)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, fMaxNumberOfTracksInPPConsidered, 100, -10, 10));

  

   //Initialize array for run numbers.
  Int_t runnumbersP10b[fNRunsP10b] = {117222, 117220, 117116, 117112, 117109, 117099, 117092, 117063, 117060, 117059, 117053, 117052, 117050, 117048, 116787, 116645, 116643, 116574, 116571, 116562, 116432, 116431, 116429,116403, 116402, 116372, 116360, 116358, 116288, 116102, 116081, 116079, 115521, 115414, 115406, 115401, 115399, 115393, 115369,115345, 115335, 115328,  115327, 115322, 115318, 115312, 115310,  115193, 115186, 115056, 114931, 114930, 114924, 114920, 114918, 114798, 114786};
  Int_t runnumbersP10c[fNRunsP10c] = {121040, 121039, 120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244,120079, 120076, 120073, 120072, 120069, 120067, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119163, 119161, 119159, 118561, 118560, 118558, 118556, 118518, 118512, 118507, 118506};
  Int_t runnumbersP10d[fNRunsP10d] = {126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168,126167, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843,  125842, 125633, 125632, 125630, 125628, 125296, 125295, 125186, 125156, 125140, 125139, 125134, 125133, 125101, 125100, 125097, 125085, 125083, 125023, 124751, 122375, 122374};
  Int_t runnumbersP10e[fNRunsP10e] = {130850, 130848, 130847, 130844, 130842, 130840, 130834, 130804, 130803, 130802, 130799, 130798, 130795, 130793, 130704, 130696, 130628, 130623, 130621, 130620, 130609, 130608, 130601, 130526, 130524, 130520, 130519, 130517, 130481, 130480, 130479, 130375, 130360, 130358, 130356, 130354, 130343, 130342, 130178, 130172, 130168, 130158, 130157, 130151, 130149, 129983, 129966, 129962, 129961, 129960, 129959, 129744, 129742, 129738, 129736, 129735, 129734, 129729, 129726, 129725, 129723, 129666, 129659, 129653, 129652, 129651, 129650, 129647, 129641, 129639, 129599, 129587, 129586, 129540, 129536, 129528, 129527, 129525, 129524, 129523, 129521, 129520, 129519, 129516, 129515, 129514, 129513, 129512, 129042, 128913, 128855, 128853, 128850, 128843, 128836, 128835, 128834, 128833, 128824, 128823, 128820, 128819, 128778, 128777, 128678, 128677, 128621, 128615, 128611, 128609, 128605, 128596, 128594, 128592, 128590, 128582, 128506, 128505, 128504, 128503, 128498, 128495, 128494, 128486, 128452, 128366}; 
  Int_t runnumbersP10h[fNRunsP10h] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};
  Int_t runnumbersP11a[fNRunsP11a] = {146860, 146859, 146858, 146856, 146824, 146817, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146402, 146369,146292, 146287, 146282, 146277, 146273, 146272, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146148, 146147, 146141, 146099, 146079, 146072, 146071, 146027, 146026, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145292, 145290, 145289, 145288};
  Int_t runnumbersP11h[fNRunsP11h] = {170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915};

  
  if (fperiod == AliAnalysisTaskBuildCorrTree::P10b){
    fNruns = fNRunsP10b;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10b[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::pp;
  }
  if (fperiod == AliAnalysisTaskBuildCorrTree::P10c){
    fNruns = fNRunsP10c;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10c[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::pp;
  }
  if (fperiod == AliAnalysisTaskBuildCorrTree::P10d){
    fNruns = fNRunsP10d;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10d[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::pp;
  }
  if (fperiod == AliAnalysisTaskBuildCorrTree::P10e){
    fNruns = fNRunsP10e;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10e[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::pp;
  }
  if (fperiod == AliAnalysisTaskBuildCorrTree::P11a){
    fNruns = fNRunsP11a;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11a[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::pp;
  }
  if (fperiod ==  AliAnalysisTaskBuildCorrTree::P10h){
    fNruns = fNRunsP10h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10h[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::PbPb;
  }
  if (fperiod == AliAnalysisTaskBuildCorrTree::P11h){
    fNruns = fNRunsP11h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11h[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskBuildCorrTree::PbPb;
  }
  
  
  //QA per run histograms:
  TH1D * eventsperrun 		= new TH1D("EventsperRun", "# Events per Run", fNruns, 0, 1);
  TH1D * TracksperRun 		= new TH1D("TracksperRun", "# tracks per Run", fNruns, 0,1);
  TH1D * selectedTracksperRun 	= new TH1D("selectedTracksperRun", "# selected tracks per Run", fNruns, 0,1);
  TH3D * NTracksVertexeta	= new TH3D("NTracksVertexEta","#selected tracks per run and vertex in eta bins",fNruns,0,1,100,-10.0,10.0,100,-3.0,3.0);
  TH2D * NEventsVertex		= new TH2D("NEventsVertex","Events per run and vertex",fNruns,0,1,100,-10.0,10.0);
  TH2D * NTracksCent		= new TH2D("NTracksCent","#selected tracks per run and vertex",fNruns,0,1,100,0.0,100.0);
  TH2D * NEventsCent		= new TH2D("NEventsCent","Events per run and vertex",fNruns,0,1,100,0.0,100.0);
  TH2D * NTracksPhi		= new TH2D("NTracksPhi","#selected tracks per run and phi", fNruns,0,1,180,0.,2*TMath::Pi());
  TH2D * NTracksPt 		= new TH2D("NTrackspT","#selected tracks per run and pT", fNruns, 0,1,100,0,20);
  for(int i=0; i<fNruns; i++){
    TString lable = Form("%i",fRunNumberList[i]);
    eventsperrun->GetXaxis()->SetBinLabel(i+1, lable);
    eventsperrun->GetXaxis()->LabelsOption("v");
    TracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
    TracksperRun->GetXaxis()->LabelsOption("v");
    selectedTracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
    selectedTracksperRun->GetXaxis()->LabelsOption("v");
    NTracksVertexeta->GetXaxis()->SetBinLabel(i+1,lable);
    NTracksVertexeta->GetXaxis()->LabelsOption("v"); 
    NEventsVertex->GetXaxis()->SetBinLabel(i+1,lable);
    NEventsVertex->GetXaxis()->LabelsOption("v"); 
    NTracksCent->GetXaxis()->SetBinLabel(i+1,lable);
    NTracksCent->GetXaxis()->LabelsOption("v"); 
    NEventsCent->GetXaxis()->SetBinLabel(i+1,lable);
    NEventsCent->GetXaxis()->LabelsOption("v"); 
    NTracksPhi->GetXaxis()->SetBinLabel(i+1,lable);
    NTracksPhi->GetXaxis()->LabelsOption("v"); 
    NTracksPt->GetXaxis()->SetBinLabel(i+1,lable);
    NTracksPt->GetXaxis()->LabelsOption("v"); 
  }
  fOutput->Add(eventsperrun);
  fOutput->Add(TracksperRun);
  fOutput->Add(selectedTracksperRun);
  fOutput->Add(NTracksVertexeta);
  fOutput->Add(NEventsVertex);
  fOutput->Add(NTracksCent);
  fOutput->Add(NEventsCent);
  fOutput->Add(NTracksPhi);
  fOutput->Add(NTracksPt);
}



void AliAnalysisTaskBuildCorrTree::SetMixingScheme(Int_t MaxNEventMix, Int_t MinNofTracksMix, TArrayD MBinEdges, TArrayD ZBinEdges)
{
  fMaxNEventMix= MaxNEventMix;
  fMinNofTracksMix = MinNofTracksMix;
  AliWarning(Form("%i", MBinEdges.GetSize()));
  for(int i=0; i<MBinEdges.GetSize()-1; ++i)
    if(MBinEdges.At(i) > MBinEdges.At(i+1)) AliFatal("edges are not sorted");
  for(int i=0; i<ZBinEdges.GetSize()-1; ++i)
    if(ZBinEdges.At(i) > ZBinEdges.At(i+1)) AliFatal("edges are not sorted");  
  fMBinEdges = MBinEdges;
  fZBinEdges = ZBinEdges;
  fMaxMult = fMBinEdges.At(fMBinEdges.GetSize()-1);
}


void AliAnalysisTaskBuildCorrTree::FillHistogram(const char* key, Double_t x)
{
  TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

void AliAnalysisTaskBuildCorrTree::FillHistogram(const char* key, Double_t x, Double_t y)
{
  TH2 * hist = dynamic_cast<TH2*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y) ;
  else AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}

void AliAnalysisTaskBuildCorrTree::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z)
{
  TH3 * hist = dynamic_cast<TH3*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y,z) ;
  else AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}




