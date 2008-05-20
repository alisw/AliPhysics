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

//----------------------------------------------------------------------------
//    Implementation of the BtoJPSItoEle reconstruction and analysis class
// Note: the two decay tracks are labelled: 0 (positive track)
//                                          1 (negative track)
// An example of usage can be found in the macro AliBtoJPSItoEleTest.C
//            Origin: G.E. Bruno    giuseppe.bruno@ba.infn.it            
//  based on Class for charm golden channel (D0->Kpi)
//----------------------------------------------------------------------------
#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TParticle.h>
#include "AliESDEvent.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliBtoJPSItoEle.h"
#include "AliBtoJPSItoEleAnalysis.h"
#include "AliLog.h"
#include "AliKFParticleBase.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumlab;
  Int_t mumpdg;
  Int_t gmumlab;
  Int_t gmumpdg;
  Int_t mumprongs;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} REFTRACK;

ClassImp(AliBtoJPSItoEleAnalysis)

//----------------------------------------------------------------------------
AliBtoJPSItoEleAnalysis::AliBtoJPSItoEleAnalysis():
fVertexOnTheFly(kFALSE),
fSim(kFALSE),
fOnlySignal(kFALSE),
fOnlyPrimaryJpsi(kFALSE),
fPID("TRDTPCparam"),
fPtCut(0.),
fd0Cut(0.), 
fMassCut(1000.),
fPidCut(0.),
//fKFPrimVertex(kFALSE),
//fKFTopConstr(kFALSE),
fKFSecondVertex(kFALSE)
{
  // Default constructor

  SetBCuts();
  SetVertex1();
}
//----------------------------------------------------------------------------
AliBtoJPSItoEleAnalysis::~AliBtoJPSItoEleAnalysis() {}
//----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::ApplySelection(const Char_t *inName,const Char_t *outName) const {
  // select candidates that pass fBCuts and write them to a new file

  TFile *inFile = TFile::Open(inName);

  TTree *treeBin=(TTree*)inFile->Get("TreeB");
  AliBtoJPSItoEleAnalysis *inAnalysis = (AliBtoJPSItoEleAnalysis*)inFile->Get("BtoJPSItoEleAnalysis");
  printf("+++\n+++  I N P U T   S T A T U S:\n+++\n");
  inAnalysis->PrintStatus();


  AliBtoJPSItoEle *d = 0; 
  treeBin->SetBranchAddress("BtoJPSItoEle",&d);
  Int_t entries = (Int_t)treeBin->GetEntries();

  printf("+++\n+++ Number of B candidates in input tree:  %d\n+++\n",entries);

  TTree *treeBout = new TTree("TreeB","Tree with selected B candidates");
  treeBout->Branch("BtoJPSItoEle","AliBtoJPSItoEle",&d,200000,0);


  Int_t okB=0;
  Int_t nSel = 0;

  for(Int_t i=0; i<entries; i++) {
    // get event from tree
    treeBin->GetEvent(i);

    if(fSim && fOnlySignal && !d->IsSignal()) continue; 

    // check if candidate passes selection (as B or Bbar)
    if(d->Select(fBCuts,okB)) {
      nSel++;
      treeBout->Fill();
    }

  }

  AliBtoJPSItoEleAnalysis *outAnalysis = (AliBtoJPSItoEleAnalysis*)inAnalysis->Clone("BtoJPSItoEleAnalysis");
  outAnalysis->SetBCuts(fBCuts);
  printf("------------------------------------------\n");
  printf("+++\n+++  O U T P U T   S T A T U S:\n+++\n");
  outAnalysis->PrintStatus();

  printf("+++\n+++ Number of B mesons in output tree:  %d\n+++\n",nSel);

  TFile* outFile = new TFile(outName,"recreate");
  treeBout->Write();
  outAnalysis->Write();
  outFile->Close();

  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEleAnalysis::CalculateTOFmass(Double_t mom,Double_t length,
					      Double_t time) const {
  // calculated the mass from momentum, track length from vertex to TOF
  // and time measured by the TOF
  if(length==0.) return -1000.;
  Double_t a = time*time/length/length;
  if(a > 1.) {
    a = TMath::Sqrt(a-1.);
  } else {
    a = -TMath::Sqrt(1.-a);
  }

  return mom*a;
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::FindCandidates(Int_t evFirst,Int_t evLast,
					const Char_t *outName) {
  // Find candidates and calculate parameters


  TString esdName="AliESDs.root";
  if(gSystem->AccessPathName(esdName.Data(),kFileExists)) {
    printf("AliBtoJPSItoEleAnalysis::FindCandidatesESD(): No ESDs file found!\n"); 
    return;
  }

  TString outName1=outName;
  TString outName2="nTotEvents.dat";

  Int_t    nTotEv=0,nBrec=0,nBrec1ev=0;
  Double_t dca;
  Double_t v2[3],mom[6],d0[2];
  Int_t    iTrkP,iTrkN,trkEntries;
  Int_t    nTrksP=0,nTrksN=0;
  Int_t    trkNum[2];
  Double_t tofmass[2];
  Int_t    okB=0;
  AliESDtrack *postrack = 0;
  AliESDtrack *negtrack = 0;

  // create the AliVertexerTracks object
  // (it will be used only if fVertexOnTheFly=kTrue)
  AliVertexerTracks *vertexer1 = new AliVertexerTracks();
  if(fVertexOnTheFly) {
    // open the mean vertex
    TFile *invtx = new TFile("AliESDVertexMean.root");
    AliESDVertex *initVertex = (AliESDVertex*)invtx->Get("vtxmean");
    invtx->Close();
    vertexer1->SetVtxStart(initVertex);
    delete invtx;
  }
  Int_t  skipped[2];
  Bool_t goodVtx1;
  
  // create tree for reconstructed decayes
  AliBtoJPSItoEle *ioBtoJPSItoEle=0;
  TTree *treeB = new TTree("TreeB","Tree with candidates");
  treeB->Branch("BtoJPSItoEle","AliBtoJPSItoEle",&ioBtoJPSItoEle,200000,0);

  // open file with tracks
  TFile *esdFile = TFile::Open(esdName.Data());
  AliESDEvent* event = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if(!tree) {
    Error("FindCandidatesESD", "no ESD tree found");
    return;
  }
  event->ReadFromTree(tree);

/*  if (fKFPrimVertex)
  AliRunLoader* runLoader = 0;
  {
    if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
    runLoader = AliRunLoader::Open(galName.Data());
    if (runLoader == 0x0) {
      cerr<<"Can not open session"<<endl;
      return;
    }
    cout << "Ok open galice.root" << endl;
    runLoader->LoadgAlice();

    gAlice = runLoader->GetAliRun();
    runLoader->LoadKinematics();
    runLoader->LoadHeader();
  } */

  // loop on events in file
  for(Int_t iEvent = evFirst; iEvent < tree->GetEntries(); iEvent++) {
    if(iEvent > evLast) break;
    tree->GetEvent(iEvent);
    Int_t ev = (Int_t)event->GetEventNumberInFile();
    printf("--- Finding B -> JPSI -> e+ e- in event %d\n",ev);
    // count the total number of events
    nTotEv++;

    trkEntries = (Int_t)event->GetNumberOfTracks();
    printf(" Number of tracks: %d\n",trkEntries);
    if(trkEntries<2) continue;

    // retrieve primary vertex from file
    if(!event->GetPrimaryVertex()) { 
      printf(" No vertex\n");
      continue;
    }
    event->GetPrimaryVertex()->GetXYZ(fV1);

    // call function which applies sigle-track selection and
    // separetes positives and negatives
    TObjArray trksP(trkEntries/2);
    Int_t    *trkEntryP   = new Int_t[trkEntries];
    TObjArray trksN(trkEntries/2);
    Int_t    *trkEntryN   = new Int_t[trkEntries];
    TTree *trkTree = new TTree();
    SelectTracks(event,trksP,trkEntryP,nTrksP,
	               trksN,trkEntryN,nTrksN);

    printf(" pos. tracks: %d    neg .tracks: %d\n",nTrksP,nTrksN);


    nBrec1ev = 0;

/*
//===================== PRIMARY VERTEX USING KF METHODS ==========================//

if (fKFPrimVertex)
{
  AliStack* stack = runLoader->Stack();

  class TESDTrackInfo
    {
    public:
    TESDTrackInfo(){}
    AliKFParticle fParticle; // KFParticle constructed from ESD track
    //Bool_t fPrimUsedFlag;    // flag says that the particle was used for primary vertex fit
    Bool_t fOK;              // is the track good enough
    Int_t mcPDG;             // Monte Carlo PDG code of the particle
    Int_t mcMotherID;        // Monte Carlo ID of its mother
    };

 // TESDTrackInfo ESDTrackInfo[trkEntries];
  TESDTrackInfo ESDTrackInfo[1000];
  for (Int_t iTr=0; iTr<trkEntries; iTr++)
    {
      TESDTrackInfo &info = ESDTrackInfo[iTr];
      info.fOK = 0;
      //info.fPrimUsedFlag = 0;
      info.mcPDG = -1;
      info.mcMotherID = -1;

      // track quality check

      AliESDtrack *pTrack = event->GetTrack(iTr);
      if( !pTrack  ) continue;
      if (pTrack->GetKinkIndex(0)>0) continue;
      if ( !( pTrack->GetStatus()&AliESDtrack::kITSrefit ) ) continue;
      //Int_t indi[12];
      //if( pTrack->GetITSclusters(indi) <5 ) continue;
      //Int_t PDG = ( pTrack->GetSigned1Pt() <0 ) ?321 :211;

      // take MC PDG

      Int_t mcID = TMath::Abs(pTrack->GetLabel());
      TParticle * part = stack->Particle(TMath::Abs(mcID));
      info.mcPDG = part->GetPdgCode();
      Int_t PDG = info.mcPDG;
      if( mcID>=0 ) info.mcMotherID = part->GetFirstMother();


      // Construct KFParticle for the track

      info.fParticle = AliKFParticle( *pTrack, PDG );
      info.fOK = 1;
    }

    // Find event primary vertex with KF methods

   AliKFVertex primVtx;
    {
   // const AliKFParticle * vSelected[trkEntries]; // Selected particles for vertex fit
   // Int_t vIndex[trkEntries];                    // Indices of selected particles
   // Bool_t vFlag[trkEntries];                    // Flags returned by the vertex finder
   const AliKFParticle * vSelected[1000]; // Selected particles for vertex fit
   Int_t vIndex[1000];                    // Indices of selected particles
   Bool_t vFlag[1000];                    // Flags returned by the vertex finder

    Int_t nSelected = 0;
    for( Int_t i = 0; i<trkEntries; i++){
      if(ESDTrackInfo[i].fOK ){
        vSelected[nSelected] = &(ESDTrackInfo[i].fParticle);
        vIndex[nSelected] = i;
        nSelected++;
      }
    }
    primVtx.ConstructPrimaryVertex( vSelected, nSelected, vFlag, 3. );
    //for( Int_t i = 0; i<nSelected; i++){
    //  if( vFlag[i] ) ESDTrackInfo[vIndex[i]].fPrimUsedFlag = 1;
    //}
    if( primVtx.GetNDF() <1 ) return; // Less then two tracks in primary vertex
    }

}*/

    // LOOP ON  POSITIVE  TRACKS
    for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
      if(iTrkP%1==0) printf("  Processing positive track number %d of %d\n",iTrkP,nTrksP);
	  
      // get track from track array
      postrack = (AliESDtrack*)trksP.UncheckedAt(iTrkP);
      trkNum[0] = trkEntryP[iTrkP];      

      // LOOP ON  NEGATIVE  TRACKS
      for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {

	// get track from tracks array
	negtrack = (AliESDtrack*)trksN.UncheckedAt(iTrkN);
	trkNum[1] = trkEntryN[iTrkN];      

	//
	// ----------- DCA MINIMIZATION ------------------
	//
	// find the DCA and propagate the tracks to the DCA
	Double_t b=event->GetMagneticField(); 
	AliESDtrack nt(*negtrack), pt(*postrack);
	dca = nt.PropagateToDCA(&pt,b);

	// define the AliESDv0 object
	AliESDv0 vertex2(nt,trkNum[0],pt,trkNum[1]);
	  
	// get position of the secondary vertex

        if (fKFSecondVertex){
          //Define the AliKFParticle Objects
          AliKFParticle trackP = AliKFParticle(pt,-11);
          AliKFParticle trackN = AliKFParticle(nt,11);
          //Construct the V0like mother
          AliKFParticle V0(trackP,trackN);
          //Get global position of the secondary vertex using KF methods
          v2[0] = V0.GetX();
          v2[1] = V0.GetY();
          v2[2] = V0.GetZ();
          mom[0] = trackP.GetPx(); mom[1] = trackP.GetPy(); mom[2] = trackP.GetPz();
          mom[3] = trackN.GetPx(); mom[4] = trackN.GetPy(); mom[5] = trackN.GetPz();
        }else{
          //Get position of the secondary vertex
          vertex2.GetXYZ(v2[0],v2[1],v2[2]);
          vertex2.GetPPxPyPz(mom[0],mom[1],mom[2]);
          vertex2.GetNPxPyPz(mom[3],mom[4],mom[5]);
        }
	goodVtx1 = kTRUE;
	
	// no vertexing if DeltaMass > fMassCut 
	if(fVertexOnTheFly) {
	  goodVtx1 = kFALSE;
	  if(SelectInvMass(mom)) {
	    // primary vertex from *other* tracks in the event
	    vertexer1->SetFieldkG(event->GetMagneticField());
	    skipped[0] = trkEntryP[iTrkP];
	    skipped[1] = trkEntryN[iTrkN];
	    vertexer1->SetSkipTracks(2,skipped);
	    AliESDVertex *vertex1onfly = 
	      (AliESDVertex*)vertexer1->FindPrimaryVertex(event); 
	    if(vertex1onfly->GetNContributors()>0) goodVtx1 = kTRUE;
	    vertex1onfly->GetXYZ(fV1);
	    //vertex1onfly->PrintStatus();
	    delete vertex1onfly;	
	  }
	}
       
/*      if (fKFSecondVertex&&fKFTopConstr&&fKFPrimVertex){
//====================== TOPOLOGICAL CONSTRAINT !!=====================================
        //Primary vertex constructed from ESD using KF methods!!!
        AliKFVertex primVtxCopy(*(event->GetPrimaryVertex()));

        //Subtract Daughters from primary vertex
        primVtxCopy -= trackP;
        primVtxCopy -= trackN;

        //Add V0 to the vertex in order to improve primary vertex resolution
        primVtxCopy += V0;

        //Set production vertex for V0
        V0.SetProductionVertex(primVtxCopy);

        //Recalculate primary vertex
        fV1[0] = primVtxCopy.GetX();
        fV1[1] = primVtxCopy.GetY();
        fV1[2] = primVtxCopy.GetZ();
//=====================================================================================
        }*/

	// impact parameters of the tracks w.r.t. the primary vertex
	d0[0] =  10000.*pt.GetD(fV1[0],fV1[1],b);
	d0[1] = -10000.*nt.GetD(fV1[0],fV1[1],b);

	// create the object AliBtoJPSItoEle
	AliBtoJPSItoEle theB(ev,trkNum,fV1,v2,dca,mom,d0);
	// select B's
	if(goodVtx1 && theB.Select(fBCuts,okB)) {
	  // get PID info from ESD
	  AliESDtrack *t0 = (AliESDtrack*)event->GetTrack(trkNum[0]);
	  Double_t esdpid0[5];
	  t0->GetESDpid(esdpid0);
	  if(t0->GetStatus()&AliESDtrack::kTOFpid) {
	    tofmass[0] = CalculateTOFmass(t0->GetP(),
					  t0->GetIntegratedLength(),
					  t0->GetTOFsignal());
	  } else {
	    tofmass[0] = -1000.;
	  }
	  AliESDtrack *t1 = (AliESDtrack*)event->GetTrack(trkNum[1]);
	  Double_t esdpid1[5];
	  t1->GetESDpid(esdpid1);
	  if(t1->GetStatus()&AliESDtrack::kTOFpid) {
	    tofmass[1] = CalculateTOFmass(t1->GetP(),
					  t1->GetIntegratedLength(),
					  t1->GetTOFsignal());
	  } else {
	    tofmass[1] = -1000.;
	  }

	  theB.SetPIDresponse(esdpid0,esdpid1);
	  theB.SetTOFmasses(tofmass);

	  // fill the tree
	  ioBtoJPSItoEle=&theB;
	  treeB->Fill();

	  nBrec++; nBrec1ev++;
	  ioBtoJPSItoEle=0; 
	}
	
	negtrack = 0;
      } // loop on negative tracks
      postrack = 0;
    }   // loop on positive tracks
    
    delete [] trkEntryP;
    delete [] trkEntryN;
    delete trkTree;

    printf(" Number of B candidates: %d\n",nBrec1ev);
  }    // loop on events in file


  printf("\n+++\n+++ Total number of events: %d\n+++\n",nTotEv);
  printf("\n+++\n+++ Total number of B candidates: %d\n+++\n",nBrec);

  delete vertexer1;

  esdFile->Close();

  // create a copy of this class to be written to output file
  AliBtoJPSItoEleAnalysis *copy = (AliBtoJPSItoEleAnalysis*)this->Clone("BtoJPSItoEleAnalysis");

  // add PDG codes to decay tracks in found candidates (in simulation mode)
  // and store tree in the output file
  if(!fSim) {
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeB->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  } else {
    printf(" Now adding information from simulation (PDG codes) ...\n");
    TTree *treeBsim = new TTree("TreeB","Tree with B candidates");
    SimulationInfo(treeB,treeBsim);
    delete treeB;
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeBsim->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  }

  // write to a file the total number of events
  FILE *outdat = fopen(outName2.Data(),"w");
  fprintf(outdat,"%d\n",nTotEv);
  fclose(outdat);

  return;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::PrintStatus() const {
  // Print parameters being used

  printf("Preselections:\n");
  printf("    fPtCut   = %f GeV\n",fPtCut);
  printf("    fd0Cut   = %f micron\n",fd0Cut);
  printf("    fMassCut = %f GeV\n",fMassCut);
  printf("    fPidCut  > %f \n",fPidCut);
  if(fVertexOnTheFly) printf("Primary vertex on the fly\n");
  if(fSim) { 
    printf("Simulation mode\n");
    if(fOnlySignal && !(fOnlyPrimaryJpsi)) printf("  Only signal goes to file\n");
    if(fOnlyPrimaryJpsi && !(fOnlySignal)) printf("  Only primary Jpsi go to file\n");
    if(fOnlyPrimaryJpsi && fOnlySignal) printf("  Both signal and primary Jpsi go to file\n");
  }
  printf("Cuts on candidates:\n");
  printf("    |M-MJPsi| [GeV]  < %f\n",fBCuts[0]);
  printf("    dca    [micron]  < %f\n",fBCuts[1]);
  printf("    cosThetaStar     < %f\n",fBCuts[2]);
  printf("    pTP     [GeV]    > %f\n",fBCuts[3]);
  printf("    pTN     [GeV]    > %f\n",fBCuts[4]);
  printf("    |d0P|  [micron]  < %f\n",fBCuts[5]);
  printf("    |d0N|  [micron]  < %f\n",fBCuts[6]);
  printf("    d0d0  [micron^2] < %f\n",fBCuts[7]);
  printf("    cosThetaPoint    > %f\n",fBCuts[8]);

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliBtoJPSItoEleAnalysis::SelectInvMass(const Double_t p[6]) const {
  // Apply preselection in the invariant mass of the pair

  Double_t mJPsi = 3.096916;
  Double_t mel   = 0.00510998902;

  Double_t energy[2];
  Double_t mom2[2],momTot2;

  mom2[0] = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  mom2[1] = p[3]*p[3] + p[4]*p[4] + p[5]*p[5];

  momTot2 = (p[0]+p[3])*(p[0]+p[3])+
            (p[1]+p[4])*(p[1]+p[4])+
            (p[2]+p[5])*(p[2]+p[5]);

  // J/Psi -> e+ e-
  energy[1] = TMath::Sqrt(mel*mel+mom2[1]);
  energy[0] = TMath::Sqrt(mel*mel+mom2[0]);

  Double_t minvJPsi = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);

  if(TMath::Abs(minvJPsi-mJPsi)  < fMassCut) return kTRUE;
  return kFALSE;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::SelectTracks(AliESDEvent *event,
        TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
        TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const {
  // Create two TObjArrays with positive and negative tracks and 
  // apply single-track preselection

  nTrksP=0,nTrksN=0;

  Int_t entr = event->GetNumberOfTracks();
 
  // transfer ITS tracks from ESD to arrays and to a tree
  for(Int_t i=0; i<entr; i++) {

    AliESDtrack *esdtrack = event->GetTrack(i);
    UInt_t status = esdtrack->GetStatus();

    if(!(status&AliESDtrack::kITSin)) continue;

    // single track selection
    if(!SingleTrkCuts(*esdtrack,event->GetMagneticField())) continue;

    if(esdtrack->GetSign()<0) { // negative track
      trksN.AddLast(esdtrack);
      trkEntryN[nTrksN] = i;
      nTrksN++;
    } else {                 // positive track
      trksP.AddLast(esdtrack);
      trkEntryP[nTrksP] = i;
      nTrksP++;
    }

  } // loop on ESD tracks

  return;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::SetBCuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8) {
  // Set the cuts for B selection
  fBCuts[0] = cut0;
  fBCuts[1] = cut1;
  fBCuts[2] = cut2;
  fBCuts[3] = cut3;
  fBCuts[4] = cut4;
  fBCuts[5] = cut5;
  fBCuts[6] = cut6;
  fBCuts[7] = cut7;
  fBCuts[8] = cut8;

  return;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::SetBCuts(const Double_t cuts[9]) {
  // Set the cuts for B selection

  for(Int_t i=0; i<9; i++) fBCuts[i] = cuts[i];

  return;
}
//-----------------------------------------------------------------------------
Bool_t 
AliBtoJPSItoEleAnalysis::SingleTrkCuts(const AliESDtrack& trk, Double_t b) const {
  // Check if track passes some kinematical cuts  
  // Magnetic field "b" (kG)

  if(TMath::Abs(1./trk.GetParameter()[4]) < fPtCut) 
    return kFALSE;
  if(TMath::Abs(10000.*trk.GetD(fV1[0],fV1[1],b)) < fd0Cut) 
    return kFALSE;
  //select only tracks with the "combined PID"
  UInt_t status = trk.GetStatus();
  if ((status&AliESDtrack::kESDpid)==0) return kTRUE;
  Double_t r[5];
  trk.GetESDpid(r);
  if(r[0] < fPidCut) return kFALSE;

  return kTRUE;
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::MakeTracksRefFile(AliRun *mygAlice,
					   Int_t evFirst,Int_t evLast) const {
  // Create a file with simulation info for the reconstructed tracks
  
  TFile *outFile = TFile::Open("BTracksRefFile.root","recreate");
  TFile *esdFile = TFile::Open("AliESDs.root");

  AliMC *mc = mygAlice->GetMCApp();
  
  Int_t      label;
  TParticle *part;  
  TParticle *mumpart;
  TParticle *gmumpart;
  REFTRACK   reftrk;
  
  AliESDEvent* event = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  event->ReadFromTree(tree);
  // loop on events in file
  for(Int_t iEvent=evFirst; iEvent<tree->GetEntries(); iEvent++) {
    if(iEvent>evLast) break;
    tree->GetEvent(iEvent);
    Int_t ev = (Int_t)event->GetEventNumberInFile();

    mygAlice->GetEvent(ev);

    Int_t nentr=(Int_t)event->GetNumberOfTracks();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",ev);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&reftrk,"lab/I:pdg:mumlab:mumpdg:gmumlab:gmumpdg:mumprongs:Vx/F:Vy:Vz:Px:Py:Pz");
//    reftree->Branch("rectracks",&reftrk,"lab/I:pdg:mumlab:mumpdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for(Int_t i=0; i<nentr; i++) {
      AliESDtrack *esdtrack = (AliESDtrack*)event->GetTrack(i);
      label = TMath::Abs(esdtrack->GetLabel());

      part = (TParticle*)mc->Particle(label); 
      reftrk.lab = label;
      reftrk.pdg = part->GetPdgCode();
      reftrk.mumlab = part->GetFirstMother();
      if(part->GetFirstMother()>=0) {
	mumpart = (TParticle*)mygAlice->GetMCApp()->Particle(part->GetFirstMother());
	reftrk.mumpdg = mumpart->GetPdgCode();
	reftrk.mumprongs = mumpart->GetNDaughters();
        reftrk.gmumlab = mumpart->GetFirstMother();
        if(mumpart->GetFirstMother()>=0) {
          gmumpart = (TParticle*)mygAlice->GetMCApp()->Particle(mumpart->GetFirstMother());
          reftrk.gmumpdg = gmumpart->GetPdgCode();
        }
      } else {
	reftrk.mumpdg=-1;
	reftrk.mumprongs=-1;
	reftrk.gmumpdg=-1;
        reftrk.gmumlab=part->GetFirstMother(); // If it hasn't any mother, then it has neither Gmother!
        // reftrk.gmumlab=-1; // If it hasn't any mother, then it has neither Gmother!
      }
      reftrk.Vx = part->Vx();
      reftrk.Vy = part->Vy();
      reftrk.Vz = part->Vz();
      reftrk.Px = part->Px();
      reftrk.Py = part->Py();
      reftrk.Pz = part->Pz();
      
      reftree->Fill();

    } // loop on tracks   

    outFile->cd();
    reftree->Write();

    delete reftree;
  } // loop on events

  esdFile->Close();
  outFile->Close();

  return;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEleAnalysis::SimulationInfo(TTree *treeBin,TTree *treeBout) const {
  // add pdg codes to candidate decay tracks (for sim)

  TString refFileName("BTracksRefFile.root");
  if(fSim && gSystem->AccessPathName(refFileName.Data(),kFileExists)) { 
    printf("AliBtoJPSItoEleAnalysis::SimulationInfo: no reference file found!\n"); 
    return;
  }
  TFile *refFile = TFile::Open(refFileName.Data());

  Char_t refTreeName[100];
  Int_t  event;
  Int_t  pdg[2],mumpdg[2],mumlab[2],gmumpdg[2],gmumlab[2];
  REFTRACK reftrk;

  // read-in reference tree for event 0 (the only event for Pb-Pb)
  sprintf(refTreeName,"Tree_Ref_%d",0);
  TTree *refTree0 = (TTree*)refFile->Get(refTreeName);
  refTree0->SetBranchAddress("rectracks",&reftrk);

  AliBtoJPSItoEle *theB = 0; 
  treeBin->SetBranchAddress("BtoJPSItoEle",&theB);
  treeBout->Branch("BtoJPSItoEle","AliBtoJPSItoEle",&theB,200000,0);

  Int_t entries = (Int_t)treeBin->GetEntries();

  for(Int_t i=0; i<entries; i++) {
    if(i%100==0) printf("  done %d candidates of %d\n",i,entries);    

    treeBin->GetEvent(i);
    event = theB->EventNo();

    if(event==0) { // always true for Pb-Pb (avoid to read-in tree every time)
      refTree0->GetEvent(theB->GetTrkNum(0));
      pdg[0]    = reftrk.pdg;
      mumpdg[0] = reftrk.mumpdg;
      mumlab[0] = reftrk.mumlab;
      gmumpdg[0] = reftrk.gmumpdg;
      gmumlab[0] = reftrk.gmumlab;
      refTree0->GetEvent(theB->GetTrkNum(1));
      pdg[1]    = reftrk.pdg;
      mumpdg[1] = reftrk.mumpdg;
      mumlab[1] = reftrk.mumlab;
      gmumpdg[1] = reftrk.gmumpdg;
      gmumlab[1] = reftrk.gmumlab;
    } else {
      sprintf(refTreeName,"Tree_Ref_%d",event);
      TTree *refTree = (TTree*)refFile->Get(refTreeName);
      refTree->SetBranchAddress("rectracks",&reftrk);
      refTree->GetEvent(theB->GetTrkNum(0));
      pdg[0]    = reftrk.pdg;
      mumpdg[0] = reftrk.mumpdg;
      mumlab[0] = reftrk.mumlab;
      gmumpdg[0] = reftrk.gmumpdg;
      gmumlab[0] = reftrk.gmumlab;
      refTree->GetEvent(theB->GetTrkNum(1));
      pdg[1]    = reftrk.pdg;
      mumpdg[1] = reftrk.mumpdg;
      mumlab[1] = reftrk.mumlab;
      gmumpdg[1] = reftrk.gmumpdg;
      gmumlab[1] = reftrk.gmumlab;
      delete refTree;
    }
    
    theB->SetPdgCodes(pdg);
    theB->SetMumPdgCodes(mumpdg);
    theB->SetGMumPdgCodes(gmumpdg);

    if(gmumpdg[0]==gmumpdg[1] &&              // Both GrandMothers are of the same sign
       (TMath::Abs(gmumpdg[0])==521   || TMath::Abs(gmumpdg[0])==511   ||  // GrandMother Bplus/Bminus or B0/B0bar
        TMath::Abs(gmumpdg[0])==523   || TMath::Abs(gmumpdg[0])==513   ||  // B0s/B0sbar
        TMath::Abs(gmumpdg[0])==515   || TMath::Abs(gmumpdg[0])==525   ||  // 
        TMath::Abs(gmumpdg[0])==531   || TMath::Abs(gmumpdg[0])==533   ||  // 
        TMath::Abs(gmumpdg[0])==535   || TMath::Abs(gmumpdg[0])==541   ||  // 
        TMath::Abs(gmumpdg[0])==543   || TMath::Abs(gmumpdg[0])==545   ||  // 
        TMath::Abs(gmumpdg[0])==10521 || TMath::Abs(gmumpdg[0])==10511 ||  //   all possible 
        TMath::Abs(gmumpdg[0])==10523 || TMath::Abs(gmumpdg[0])==10513 ||  //    B mesons
        TMath::Abs(gmumpdg[0])==20523 || TMath::Abs(gmumpdg[0])==20513 ||  // 
        TMath::Abs(gmumpdg[0])==10531 || TMath::Abs(gmumpdg[0])==10533 ||  // 
        TMath::Abs(gmumpdg[0])==20533 || TMath::Abs(gmumpdg[0])==10541 ||  // 
        TMath::Abs(gmumpdg[0])==20543 || TMath::Abs(gmumpdg[0])==10543 ||  // 
        TMath::Abs(gmumpdg[0])==4122  || TMath::Abs(gmumpdg[0])==4222  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4212  || TMath::Abs(gmumpdg[0])==4112  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4224  || TMath::Abs(gmumpdg[0])==4214  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4114  || TMath::Abs(gmumpdg[0])==4232  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4132  || TMath::Abs(gmumpdg[0])==4322  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4312  || TMath::Abs(gmumpdg[0])==4324  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4314  || TMath::Abs(gmumpdg[0])==4332  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4334  || TMath::Abs(gmumpdg[0])==4412  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4422  || TMath::Abs(gmumpdg[0])==4414  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4424  || TMath::Abs(gmumpdg[0])==4432  ||  // All possible B baryons
        TMath::Abs(gmumpdg[0])==4434  || TMath::Abs(gmumpdg[0])==4444      // All possible B baryons
       ) &&
       mumpdg[0]==443 && mumpdg[1]== 443 && 
       mumlab[0]==mumlab[1] && 
       reftrk.mumprongs==2 &&
       pdg[0]==-11 && pdg[1]==11  
     ) theB->SetSignal();

    else if (  // here consider the case of primary J/psi 
       mumpdg[0]==443 && mumpdg[1]== 443 &&
       pdg[0]==-11 && pdg[1]==11 &&
       mumlab[0]==mumlab[1] &&
       reftrk.mumprongs==2 &&
       ( gmumlab[0]<0                     ||   // really primary J/psi (without family. e.g. from Cocktail)
         TMath::Abs(gmumpdg[0])==100443   ||   // from Psi(2S)
         TMath::Abs(gmumpdg[0])==10441    ||   // from Csi_c0(1P)
         TMath::Abs(gmumpdg[0])==20443    ||   // from Csi_c1(1P)
         TMath::Abs(gmumpdg[0])==10443    ||   // from h_c(1P)
         TMath::Abs(gmumpdg[0])==445      ||   // from Csi_c2(1P)
         (gmumpdg[0]>=81 && gmumpdg[0]<=100)   // this is for MC internal use (e.g. J/psi from string)
       )
     ) theB->SetJpsiPrimary();

    // if(!fOnlySignal || theB->IsSignal()) treeBout->Fill();  
    
 // write it out 1) always if you have not asked for onlySignal or OnlyPrimaryJpsi (or both)
 //          or  2) if you have asked for Signal and it is Signal
 //          or  3) if you have asked for Primary Jpsi and it is a Primary Jpsi
    if ( (!fOnlySignal && !fOnlyPrimaryJpsi) || (fOnlySignal && theB->IsSignal()) 
        || (fOnlyPrimaryJpsi && theB->IsJpsiPrimary()) ) treeBout->Fill();
  }

  delete refTree0;

  refFile->Close();

  return;
}





