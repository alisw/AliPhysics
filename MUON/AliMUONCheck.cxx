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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONCheck
///
/// This class check the ESD tree, providing the matching with the trigger
/// response and designing useful plots (Pt, Y, ITS vertex, multiplicity).
///  Note that there is a special flag to turn on for pdc06 production. 
/// It also checks the TR tree giving hit densities on the two first and 
/// last planes of the spectrometer as well as the time of flight on these planes.
/// MUONkine() provides event stack and check if the event are generated with 
/// at least one muon and two muons (for PDC06).
/// DumpDigit() as a replacement of the function from MUONCheck.C macro.
///
/// \author Frederic Yermia, INFN Torino
//-----------------------------------------------------------------------------

#include "AliMUONCheck.h"
#include "AliMUONConstants.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMpCDB.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpCDB.h"
#include "AliMUONVDigitStore.h"

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliMagF.h"
#include "AliLog.h"

#include <TSystem.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>

/// \cond CLASSIMP
ClassImp(AliMUONCheck)
/// \endcond

//_____________________________________________________________________________
const TString& AliMUONCheck::GetDefaultOutFileName()
{
  /// Default output file name 
  static const TString kDefaultOutFileName = "output.txt";
  return kDefaultOutFileName;
}  

//_____________________________________________________________________________
AliMUONCheck::AliMUONCheck(const char* galiceFile, const char* esdFile,Int_t firstEvent, Int_t lastEvent,const char* outDir) 
: TObject(),
fFileName(galiceFile),
fFileNameSim(),
fesdFileName(esdFile),
fkOutDir(outDir),
fOutFileName(GetDefaultOutFileName()),
fFirstEvent(firstEvent),
fLastEvent(lastEvent)
{
  /// ctor  
}

//_____________________________________________________________________________
AliMUONCheck::AliMUONCheck(const char* galiceFile, const char* galiceFileSim,
                           const char* esdFile,Int_t firstEvent, Int_t lastEvent,
                           const char* outDir) 
: TObject(),
fFileName(galiceFile),
fFileNameSim(galiceFileSim),
fesdFileName(esdFile),
fkOutDir(outDir),
fOutFileName(GetDefaultOutFileName()),
fFirstEvent(firstEvent),
fLastEvent(lastEvent)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONCheck::~AliMUONCheck()
{
  /// Destructor
}

//_____________________________________________________________________________
void
AliMUONCheck::CheckESD(Bool_t pdc06TriggerResponse) 
{
  /// Check ESD files
  
  // Histograms
  TH1F * fhMUONVertex ; //! 
  TH1F * fhMUONMult   ; //!
  
  // create histograms 
  fhMUONVertex = new TH1F("hMUONVertex","ITS Vertex"                ,100, -25., 25.);
  fhMUONMult   = new TH1F("hMUONMult"  ,"Multiplicity of ESD tracks",10,  -0.5, 9.5);
  
  TH1F *hY = new TH1F("hY","Rapidity",100,-5.,-1.);
  TH1F *hPt = new TH1F("hPt","Pt",100, 0.,20.);
  
  // ------------->open the ESD file
  TFile* esdFile = TFile::Open(fesdFileName.Data());
  
  if (!esdFile || !esdFile->IsOpen()) 
  {
    AliError(Form("Error opening %s file \n",fesdFileName.Data()));
    return;
  }
  
  Int_t fSPLowpt=0     ; //!
  Int_t fSPHighpt=0    ; //!
  Int_t fSPAllpt=0     ; //!
  Int_t fSMLowpt=0     ; //!
  Int_t fSMHighpt =0   ; //!
  Int_t fSMAllpt=0     ; //!
  Int_t fSULowpt=0     ; //!
  Int_t fSUHighpt=0    ; //!
  Int_t fSUAllpt=0     ; //!
  Int_t fUSLowpt=0     ; //!
  Int_t fUSHighpt=0    ; //!
  Int_t fUSAllpt=0     ; //! 
  Int_t fLSLowpt=0     ; //!
  Int_t fLSHighpt=0    ; //! 
  Int_t fLSAllpt=0     ; //!
  
  Int_t fSLowpt=0      ; //!
  Int_t fSHighpt=0     ; //!
  
  Int_t fnTrackTrig=0  ; //!
  Int_t ftracktot=0    ; //!
  Int_t effMatch=0     ; //!
  
  TLorentzVector fV1;
  Float_t muonMass = 0.105658389;
  Double_t thetaX, thetaY, pYZ;
  Double_t fPxRec1, fPyRec1, fPzRec1, fE1;
  Int_t fZ1;
  
  AliESDEvent* fESD = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) 
  {
    Error("CheckESD", "no ESD tree found");
    AliError(Form("no ESD tree found"));
    return ;
  }
  fESD->ReadFromTree(tree);
  
  Int_t fnevents = tree->GetEntries();
  Int_t endOfLoop = fLastEvent+1;
  
  if ( fLastEvent == -1 ) endOfLoop = fnevents;
  Int_t ievent=0;
  Int_t nev=0;
  
  for (ievent = fFirstEvent; ievent < endOfLoop; ++ievent ) 
  {
    nev++;    
    if (tree->GetEvent(ievent) <= 0) 
    {
      Error("CheckESD", "no ESD object found for event %d", ievent);
      return ;
    }
    AliESDVertex* vertex = (AliESDVertex*) fESD->GetVertex();
    
    Double_t zVertex = 0. ;
    if (vertex) zVertex = vertex->GetZv();
    
    Int_t nTracks = (Int_t)fESD->GetNumberOfMuonTracks() ;
    ULong64_t trigword=fESD->GetTriggerMask();
    
    if(pdc06TriggerResponse)
    {
      if (trigword & 0x01) {
        fSPLowpt++;
      }
      
      if (trigword & 0x02){
        fSPHighpt++;
      }
      if (trigword & 0x04){
        fSPAllpt++;
      } 
      if (trigword & 0x08){
        fSMLowpt++;
      }  
      if (trigword & 0x010){
        fSMHighpt++;
      }
      if (trigword & 0x020){
        fSMAllpt++;
      } 
      if (trigword & 0x040){
        fSULowpt++;
      }  
      if (trigword & 0x080){
        fSUHighpt++;
      }   
      if (trigword & 0x100){
        fSUAllpt++;
      }  
      if (trigword & 0x200){
        fUSLowpt++;
      }
      
      if (trigword & 0x400){
        fUSHighpt++;
      }
      if (trigword & 0x800){
        fUSAllpt++;
      }
      if (trigword & 0x1000){
        fLSLowpt++;
      }
      
      if (trigword & 0x2000){
        fLSHighpt++;
      }
      
      if (trigword & 0x4000){
        fLSAllpt++;
      }
    }// if pdc06TriggerResponse
    else {
      if (trigword & 0x01) {
        fSLowpt++;
      }
      
      if (trigword & 0x02){
        fSHighpt++;
      }
      if (trigword & 0x04){
        fLSLowpt++;
      } 
      if (trigword & 0x08){
        fLSHighpt++;
      }  
      if (trigword & 0x010){
        fUSLowpt++;
      }
      if (trigword & 0x020){
        fUSHighpt++;
      }
    }
    
    Int_t tracktrig=0;
    
    for ( Int_t iTrack1 = 0; iTrack1<nTracks; ++iTrack1 ) 
    { //1st loop
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack1);
      
      // skip fake tracks (ghosts)
      if (!muonTrack->ContainTrackerData()) continue;
      
      ftracktot++;
      
      thetaX = muonTrack->GetThetaX();
      thetaY = muonTrack->GetThetaY();
      pYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      
      fPzRec1  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaY));
      fPxRec1  = fPzRec1 * TMath::Tan(thetaX);
      fPyRec1  = fPzRec1 * TMath::Tan(thetaY);
      fZ1 = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
      fE1 = TMath::Sqrt(muonMass * muonMass + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
      fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);
      // -----------> transverse momentum
      Float_t pt1 = fV1.Pt();
      // ----------->Rapidity
      Float_t y1 = fV1.Rapidity();
      
      if(muonTrack->GetMatchTrigger()) 
      {
        fnTrackTrig++;
        tracktrig++;
      }
      hY->Fill(y1);
      hPt->Fill(pt1);
    } // loop on track
    
    fhMUONVertex->Fill(zVertex) ;
    fhMUONMult->Fill(Float_t(nTracks)) ;
    
  } // loop over events
  
  AliInfo(Form("Terminate %s:", GetName())) ;
  
  effMatch=100*fnTrackTrig/ftracktot;
  
  if(pdc06TriggerResponse)
  {
    printf("=================================================================\n") ;
    printf("================  %s ESD SUMMARY    ==============\n", GetName()) ;
    printf("                                                   \n") ;
    printf("         Total number of processed events  %d      \n", nev) ;
    printf("\n")  ;
    printf("\n")  ;
    printf("Table 1:                                         \n") ;
    printf(" Global Trigger output       Low pt  High pt   All\n") ;
    printf(" number of Single Plus      :\t");
    printf("%i\t%i\t%i\t", fSPLowpt, fSPHighpt, fSPAllpt) ;
    printf("\n");
    printf(" number of Single Minus     :\t");
    printf("%i\t%i\t%i\t", fSMLowpt, fSMHighpt, fSMAllpt) ;
    printf("\n");
    printf(" number of Single Undefined :\t"); 
    printf("%i\t%i\t%i\t", fSULowpt, fSUHighpt, fSUAllpt) ;
    printf("\n");
    printf(" number of UnlikeSign pair  :\t"); 
    printf("%i\t%i\t%i\t", fUSLowpt, fUSHighpt, fUSAllpt) ;
    printf("\n");
    printf(" number of LikeSign pair    :\t");  
    printf("%i\t%i\t%i\t", fLSLowpt, fLSHighpt, fLSAllpt) ;
    printf("\n");
    printf("===================================================\n") ;
    printf("\n") ;
    printf("matching efficiency with the trigger for single tracks = %2d %% \n", effMatch);
    printf("================================================================\n") ;  printf("\n") ;
    
  }//if(pdc06TriggerResponse)
  
  gSystem->cd(fkOutDir);
  
  FILE *outtxt=fopen(fOutFileName.Data(),"a");
  
  if(pdc06TriggerResponse){
    fprintf(outtxt,"                                                   \n");
    fprintf(outtxt,"===================================================\n");
    fprintf(outtxt,"================      ESD SUMMARY    ==============\n");
    fprintf(outtxt,"                                                   \n");
    fprintf(outtxt,"    Total number of processed events  %d      \n", nev); 
    fprintf(outtxt,"\n");
    fprintf(outtxt,"\n");
    fprintf(outtxt,"Table 1:                                         \n");
    fprintf(outtxt," Global Trigger output       Low pt  High pt   All\n");
    fprintf(outtxt," number of Single Plus      :\t");
    fprintf(outtxt,"%i\t%i\t%i\t",fSPLowpt,fSPHighpt,fSPAllpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of Single Minus     :\t");
    fprintf(outtxt,"%i\t%i\t%i\t",fSMLowpt,fSMHighpt,fSMAllpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of Single Undefined :\t"); 
    fprintf(outtxt,"%i\t%i\t%i\t",fSULowpt,fSUHighpt,fSUAllpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of UnlikeSign pair  :\t"); 
    fprintf(outtxt,"%i\t%i\t%i\t",fUSLowpt,fUSHighpt,fUSAllpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of LikeSign pair    :\t");  
    fprintf(outtxt,"%i\t%i\t%i\t",fLSLowpt,fLSHighpt, fLSAllpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt,"===================================================\n");
    fprintf(outtxt,"\n");
    fprintf(outtxt,"matching efficiency with the trigger for single tracks = %2d %% \n", effMatch);
  }//if(pdc06TriggerResponse)
  
  else {
    
    fprintf(outtxt,"                                                   \n");
    fprintf(outtxt,"===================================================\n");
    fprintf(outtxt,"================      ESD SUMMARY    ==============\n");
    fprintf(outtxt,"                                                   \n");
    fprintf(outtxt,"    Total number of processed events  %d      \n", nev); 
    fprintf(outtxt,"\n");
    fprintf(outtxt,"\n");
    fprintf(outtxt,"Table 1:                                         \n");
    fprintf(outtxt," Global Trigger output       Low pt  High pt     \n");
    fprintf(outtxt," number of Single       :\t");
    fprintf(outtxt,"%i\t%i\t",fSLowpt,fSHighpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of UnlikeSign pair :\t"); 
    fprintf(outtxt,"%i\t%i\t",fUSLowpt,fUSHighpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt," number of LikeSign pair    :\t");  
    fprintf(outtxt,"%i\t%i\t",fLSLowpt,fLSHighpt);
    fprintf(outtxt,"\n");
    fprintf(outtxt,"===================================================\n");
    fprintf(outtxt,"\n");
    fprintf(outtxt,"matching efficiency with the trigger for single tracks = %2d %% \n", effMatch);
  }//else
  fclose(outtxt);
  
  TCanvas * c1 = new TCanvas("c1", "ESD", 400, 10, 600, 700) ;
  c1->Divide(1,2) ;
  c1->cd(1) ;
  fhMUONVertex->Draw() ;
  c1->cd(2) ;
  fhMUONMult->Draw() ;  
  c1->Print("VertexAndMul.eps") ; 
  TCanvas *c2 = new TCanvas("c2","ESD",400,10,600,700);
  c2->Divide(1,2);
  c2->cd(1);
  hY->Draw();
  c2->cd(2);
  hPt->Draw();
  c2->Print("YandPt.eps") ; 
}

//_____________________________________________________________________________
void
AliMUONCheck::CheckKine() 
{
  /// Check Stack 
  
  AliMUONMCDataInterface diSim(fFileNameSim.Data());
  if (!diSim.IsValid()) return;
  
  Int_t fnevents = diSim.NumberOfEvents();
  
  Int_t endOfLoop = fLastEvent+1;
  
  if ( fLastEvent == -1 ) endOfLoop = fnevents;
  
  Int_t nev=0;
  Int_t nmu=0;
  Int_t nonemu=0;
  Int_t ndimu=0;
  
  for ( Int_t ievent = fFirstEvent; ievent < endOfLoop; ++ievent ) 
  {
    Int_t nmu2=0;
    ++nev;  
    
    AliStack* stack = diSim.Stack(ievent);
    Int_t npa = stack->GetNprimary();
    Int_t npb = stack->GetNtrack(); 
    printf("Primary particles  %i   \n",npa); 
    printf("Sec particles  %i   \n",npb); 
    printf("=================================================================\n") ;
    printf("Primary particles listing:  \n"); 
    printf("=================================================================\n") ;
    for (Int_t i=0; i<npa; ++i) 
    {
      TParticle *p  = stack->Particle(i);
      p->Print("");
      Int_t pdg=p->GetPdgCode(); 
      
      if (abs(pdg) == 13) 
      {
        ++nmu2;
      }
    }
    printf("=================================================================\n") ;
    printf("=================================================================\n") ;
    
    printf("Secondaries particles listing:  \n"); 
    printf("=================================================================\n") ;
    for (Int_t i=npa; i<npb; ++i) 
    {
      stack->Particle(i)->Print("");
    }
    
    printf("=================================================================\n") ; 
    printf(">>> Event %d, Number of primary particles is %d \n",ievent, npa); 
    printf(">>> Event %d, Number of secondary articles is %d \n",ievent, npb-npa); 
    printf("=================================================================\n");
    if(nmu2>0)
    {
      printf(">>> Okay!!! Event %d with at least one muon on primary stack! \n",ievent); 
      ++nonemu;
    }
    
    if(nmu2==0)
    {
      printf(">>> Warning!!! Event %d without muon on primary stack! \n",ievent);     
      ++nmu;
    }
    
    if(nmu2>1)
    {
      printf(">>> Okay!!! Event %d with at least two muons on primary stack! \n",ievent); 
      ++ndimu; 
    }
    printf("=================================================================\n");  
    printf("                                                                  \n");
    printf("                                                                  \n") ;
  }//ievent
  
  printf("=================================================================\n") ;
  printf("               Total number of processed events  %d               \n", nev) ;
  printf("                                                                 \n") ;
  
  if(nmu>0)
  {
    printf("--->                       WARNING!!!                       <---\n"); 
    printf(" %i events without muon on primary stack \n",nmu); 
  }
  
  if(nmu==0)
  {
    printf("--->                          OKAY!!!                        <---\n"); 
    printf("  %i events generated with at least one muon on primary stack \n",nonemu);
  }
  
  if(ndimu>0)
  {
    printf("--->                          OKAY!!!                        <---\n"); 
    printf("  %i events generated with at least two muons on primary stack \n",ndimu); 
  }
  
  printf("                                                                 \n") ;
  printf("***                       Leaving MuonKine()                 *** \n");
  printf("**************************************************************** \n");
  
  gSystem->cd(fkOutDir);
  FILE *outtxt=fopen(fOutFileName.Data(),"a");
  fprintf(outtxt,"                                                   \n");
  fprintf(outtxt,"=================================================================\n");
  fprintf(outtxt,"================         MUONkine SUMMARY        ================\n");
  fprintf(outtxt,"\n");
  fprintf(outtxt,"=================================================================\n");
  fprintf(outtxt,"               Total number of processed events  %d              \n", nev) ;
  fprintf(outtxt,"                                                                 \n");
  
  if(nmu>0)
  {
    fprintf(outtxt,"                        ---> WARNING!!! <---                     \n"); 
    fprintf(outtxt,"  %i events without muon on primary stack \n",nmu); 
  }
  
  if(nmu==0)
  {
    fprintf(outtxt,"                         ---> OKAY!!! <---                       \n"); 
    fprintf(outtxt,"  %i events generated with at least one muon on primary stack \n",nonemu); 
  }
  
  if(ndimu>0)
  {
    fprintf(outtxt,"                         ---> OKAY!!! <---                       \n"); 
    fprintf(outtxt,"  %i events generated with at least two muons on primary stack \n",ndimu); 
  }
  
  fprintf(outtxt,"                                                                 \n") ;
  fprintf(outtxt,"***                       Leaving MuonKine()                 *** \n");
  fprintf(outtxt,"**************************************************************** \n");
  fclose(outtxt);
}

//_____________________________________________________________________________
void
AliMUONCheck::CheckTrackRef() 
{
  /// Check TrackRef files
  
  AliMUONMCDataInterface diSim(fFileNameSim.Data());
  if ( !diSim.IsValid() ) return;
  
  Int_t flag11=0,flag12=0,flag13=0,flag14=0;
  
  TH1F *tof01= new TH1F("tof01","TOF for first tracking plane",100,0.,100);
  tof01->SetXTitle("tof (ns)");
  TH1F *tof14= new TH1F("tof14","TOF for MT22",100,0.,100);
  tof14->SetXTitle("tof (ns)");
  
  TH1F   *hitDensity[4];
  hitDensity[0] =  new TH1F("TR_dhits01","",30,0,300);
  hitDensity[0]->SetFillColor(3);
  hitDensity[0]->SetXTitle("R (cm)");
  hitDensity[1] =  new TH1F("TR_dhits10","",30,0,300);
  hitDensity[1]->SetFillColor(3);
  hitDensity[1]->SetXTitle("R (cm)");
  hitDensity[2] =  new TH1F("TR_dhits11","",30,0,300);
  hitDensity[2]->SetFillColor(3);
  hitDensity[2]->SetXTitle("R (cm)");
  hitDensity[3] =  new TH1F("TR_dhits14","",30,0,300);
  hitDensity[3]->SetFillColor(3);
  hitDensity[3]->SetXTitle("R (cm)");
  
  Int_t fnevents = diSim.NumberOfEvents();
  
  Int_t endOfLoop = fLastEvent+1;
  
  if ( fLastEvent == -1 ) endOfLoop = fnevents;
  
  Int_t nev=0;
  Int_t ntot=fLastEvent+1-fFirstEvent;
  
  for ( Int_t ievent = fFirstEvent; ievent < endOfLoop; ++ievent ) 
  {
    Int_t  save=-99;
    ++nev;  
    
    Int_t nentries = diSim.NumberOfTrackRefs(ievent);
    
    for ( Int_t l=0; l<nentries; ++l )
    {
      TClonesArray* trackRefs = diSim.TrackRefs(ievent,l);
      if (!trackRefs) continue;
      
      Int_t nnn = trackRefs->GetEntriesFast();
      
      for ( Int_t k=0; k<nnn; ++k ) 
      {
        AliTrackReference *tref = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(k));
        Int_t label = tref->GetTrack();
        Float_t x     =    tref->X();        // x-pos of hit
        Float_t y     =    tref->Y();        // y-pos
        Float_t z     = tref->Z();
        
        Float_t r=TMath::Sqrt(x*x+y*y);
        Float_t time =    tref->GetTime();  
        
        Float_t wgt=1/(2*10*TMath::Pi()*r)/(ntot);
        
        if (save!=label){
          save=label;
          flag11=0;
          flag12=0;
          flag13=0;
          flag14=0;
        }
        
        if (save==label){
          
          //Ch 1, z=-526.16
          if (z<=-521&& z>=-531&&flag11==0){
            flag11=1;
            hitDensity[0]->Fill(r,wgt);
            tof01->Fill(1000000000*time,1);
          };
          
          //Ch 10, z=-1437.6
          if (z<=-1432&&z>=-1442&&flag12==0){
            flag12=1;
            hitDensity[1]->Fill(r,wgt);
          }
          
          //Ch 11, z=-1603.5
          if (z<=-1598&& z>=-1608&&flag13==0){
            flag13=1;
            hitDensity[2]->Fill(r,wgt);
          };
          
          //ch 14 z=-1720.5    
          if(z<=-1715&&z>=-1725&&flag14==0){
            flag14=1;
            hitDensity[3]->Fill(r,wgt);
            tof14->Fill(1000000000*time,1);
          }; 
          
        }//if save==label
        
      }//hits de tTR
      
    }//entree de tTR 
    
  }//evt loop
  
  gSystem->cd(fkOutDir);
  TCanvas *c6 = new TCanvas("c6","TOF",400,10,600,700);
  c6->Divide(1,2);
  c6->cd(1);
  
  tof01->Draw();
  c6->cd(2);
  tof14->Draw();
  c6->Print("tof_on_trigger.ps");
  
  TCanvas *c5 = new TCanvas("c5","TRef:Hits Density",400,10,600,700);
  c5->Divide(2,2);
  c5->cd(1);
  hitDensity[0]->Draw();
  c5->cd(2);
  hitDensity[1]->Draw();
  c5->cd(3);
  hitDensity[2]->Draw();
  c5->cd(4);
  hitDensity[3]->Draw();
  c5->Print("TR_Hit_densities.ps");
  printf("=================================================================\n") ;
  printf("================  %s Tref SUMMARY    ==============\n", GetName()) ;
  printf("                                                   \n") ;
  printf("         Total number of processed events  %d      \n", nev) ;
  printf("***                Leaving TRef()               *** \n");
  printf("*************************************************** \n");
}

//_____________________________________________________________________________
void 
AliMUONCheck::CheckOccupancy(Bool_t perDetEle) const
{
  /// Check occupancy for the first event selected
  
  Int_t dEoccupancyBending[14][26];
  Int_t dEoccupancyNonBending[14][26];
  Int_t cHoccupancyBending[14];
  Int_t cHoccupancyNonBending[14];
  Int_t totaloccupancyBending =0;
  Int_t totaloccupancyNonBending =0;
  
  Int_t dEchannelsBending[14][26];
  Int_t dEchannelsNonBending[14][26];
  Int_t cHchannelsBending[14];
  Int_t cHchannelsNonBending[14];
  Int_t totalchannelsBending =0;
  Int_t totalchannelsNonBending =0;
  
  Int_t nchambers = AliMUONConstants::NCh();
  
  AliMUONDataInterface di(fFileNameSim);
  
  AliMUONVDigitStore* digitStore = di.DigitStore(fFirstEvent);
  
  // Compute values
  for (Int_t ichamber=0; ichamber<nchambers; ++ichamber) 
  {
    cHchannelsBending[ichamber]=0;
    cHchannelsNonBending[ichamber]=0;
    cHoccupancyBending[ichamber]=0;
    cHoccupancyNonBending[ichamber]=0;
    
    for (Int_t idetele=0; idetele<26; idetele++) 
    {
      Int_t detele = 100*(ichamber +1)+idetele;
      
      if ( AliMpDEManager::IsValidDetElemId(detele) ) 
      {
        Int_t cathode(0);
        
        const AliMpVSegmentation* segbend = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detele, AliMp::kCath0);
        const AliMpVSegmentation* segnonbend = AliMpSegmentation::Instance()
          ->GetMpSegmentation(detele, AliMp::kCath1);
        
        if (AliMpDEManager::GetPlaneType(detele, AliMp::kCath0) != AliMp::kBendingPlane ) 
        {
          const AliMpVSegmentation* tmp = segbend;
          segbend    =  segnonbend;
          segnonbend =  tmp;
          cathode = 1;
        }  
        
        Int_t nchannels = segbend->NofPads();
        Int_t ndigits = digitStore->GetSize(detele,cathode);
	      dEchannelsBending[ichamber][idetele] = nchannels;
        dEoccupancyBending[ichamber][idetele] = ndigits;
	      cHchannelsBending[ichamber] += nchannels;
        cHoccupancyBending[ichamber] += ndigits;
	      totalchannelsBending += nchannels;
        totaloccupancyBending += ndigits;
        
        nchannels = segnonbend->NofPads();
        ndigits = digitStore->GetSize(detele,1-cathode);
        
	      dEchannelsNonBending[ichamber][idetele] = nchannels;
        dEoccupancyNonBending[ichamber][idetele] = ndigits;
	      cHchannelsNonBending[ichamber] += nchannels;
	      cHoccupancyNonBending[ichamber] += ndigits;
	      totalchannelsNonBending += nchannels;
	      totaloccupancyNonBending += ndigits;
	    }
      if (perDetEle) 
      {
        printf(">>> Detection element %4d has %5d channels in bending and %5d channels in nonbending \n",
               detele, dEchannelsBending[ichamber][idetele], dEchannelsNonBending[ichamber][idetele] ); 
      }
    }
    printf(">>> Chamber %2d has %6d channels in bending and %6d channels in nonbending \n",
           ichamber+1,  cHchannelsBending[ichamber], cHchannelsNonBending[ichamber]);
  }
  printf(">>Spectrometer has  %7d channels in bending and %7d channels in nonbending \n",
         totalchannelsBending, totalchannelsNonBending);
  
  // Output values
  
  for ( Int_t ichamber = 0; ichamber < nchambers; ++ichamber ) 
  {
    printf(">>> Chamber %2d  nChannels Bending %5d  nChannels NonBending %5d \n", 
         ichamber+1, 
         cHoccupancyBending[ichamber],
         cHoccupancyNonBending[ichamber]);  
    if ( cHchannelsBending[ichamber] != 0 && cHchannelsBending[ichamber] ) {               
      printf(">>> Chamber %2d  Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n", 
           ichamber+1, 
           100.*((Float_t) cHoccupancyBending[ichamber])/((Float_t) cHchannelsBending[ichamber]),
           100.*((Float_t) cHoccupancyNonBending[ichamber])/((Float_t) cHchannelsBending[ichamber]));
    }
  
    if ( perDetEle )
    {
      for(Int_t idetele=0; idetele<26; idetele++) 
      {
        Int_t detele = idetele + 100*(ichamber+1);
        if ( AliMpDEManager::IsValidDetElemId(detele) ) 
        {
          printf(">>> DetEle %4d nChannels Bending %5d  nChannels NonBending %5d \n", 
               idetele+100*(ichamber+1), 
               dEoccupancyBending[ichamber][idetele],
               dEoccupancyNonBending[ichamber][idetele]);  
               
          if ( dEchannelsBending[ichamber][idetele] != 0 && dEchannelsBending[ichamber][idetele] !=0 ) {     
            printf(">>> DetEle %4d Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n", 
                 idetele+100*(ichamber+1), 
                 100.*((Float_t) dEoccupancyBending[ichamber][idetele])/((Float_t) dEchannelsBending[ichamber][idetele]),
                 100.*((Float_t) dEoccupancyNonBending[ichamber][idetele])/((Float_t) dEchannelsBending[ichamber][idetele]));
          }       
        }
      }
    }
  }

  if ( totalchannelsBending != 0 && totalchannelsNonBending != 0 ) {
    printf(">>> Muon Spectrometer  Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n",  
         100.*((Float_t) totaloccupancyBending)/((Float_t) totalchannelsBending),
         100.*((Float_t) totaloccupancyNonBending)/((Float_t) totalchannelsNonBending));
  }       
}

//_____________________________________________________________________________
void AliMUONCheck::SetEventsToCheck(Int_t firstEvent, Int_t lastEvent)
{
  /// Set first and last event number to check
  
  fFirstEvent = firstEvent;
  fLastEvent = lastEvent;
}  
