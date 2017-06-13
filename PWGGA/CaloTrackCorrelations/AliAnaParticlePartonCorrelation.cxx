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

// --- ROOT system ---
#include "TH2F.h"
#include "TClass.h"

//---- ANALYSIS system ----
#include "AliAnaParticlePartonCorrelation.h" 
#include "AliMCEvent.h"  
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliVParticle.h"

/// \cond CLASSIMP
ClassImp(AliAnaParticlePartonCorrelation) ;
/// \endcond

//________________________________________________________________
/// Default Constructor. Initialize parameters.
//________________________________________________________________
AliAnaParticlePartonCorrelation::AliAnaParticlePartonCorrelation() :
AliAnaCaloTrackCorrBaseClass(),   
fhDeltaEtaNearParton(0), fhDeltaPhiNearParton(0), 
fhDeltaPtNearParton(0),  fhPtRatNearParton(0),
fhDeltaEtaAwayParton(0), fhDeltaPhiAwayParton(0), 
fhDeltaPtAwayParton(0),  fhPtRatAwayParton(0)
{
  InitParameters();
}

//________________________________________________________________
/// Create histograms to be saved in output file.
//________________________________________________________________
TList *  AliAnaParticlePartonCorrelation::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ParticlePartonHistos") ; 
  
  fhDeltaPhiNearParton  = new TH2F
  ("DeltaPhiNearParton","#phi_{particle} - #phi_{parton} vs p_{T particle}",
   200,0,120,200,0,6.4); 
  fhDeltaPhiNearParton->SetYTitle("#Delta #phi");
  fhDeltaPhiNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaPhiNearParton);
  
  fhDeltaEtaNearParton  = new TH2F
  ("DeltaEtaNearParton","#eta_{particle} - #eta_{parton} vs p_{T particle}",
   200,0,120,200,-2,2); 
  fhDeltaEtaNearParton->SetYTitle("#Delta #eta");
  fhDeltaEtaNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaEtaNearParton);
  
  fhDeltaPtNearParton  = new TH2F
  ("DeltaPtNearParton","#p_{T particle} - #p_{T parton} vs p_{T particle}",
   200,0,120,100,-10,10); 
  fhDeltaPtNearParton->SetYTitle("#Delta #p_{T}");
  fhDeltaPtNearParton->SetXTitle("p_{T particle} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtNearParton);
  
  fhPtRatNearParton  = new TH2F
  ("PtRatNearParton","#p_{T parton} / #p_{T particle} vs p_{T particle}",
   200,0,120,200,0,5); 
  fhPtRatNearParton->SetYTitle("ratio");
  fhPtRatNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhPtRatNearParton);
  
  fhDeltaPhiAwayParton  = new TH2F
  ("DeltaPhiAwayParton","#phi_{particle} - #phi_{parton} vs p_{T particle}",
   200,0,120,200,0,6.4); 
  fhDeltaPhiAwayParton->SetYTitle("#Delta #phi");
  fhDeltaPhiAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaPhiAwayParton);
  
  fhDeltaEtaAwayParton  = new TH2F
  ("DeltaEtaAwayParton","#eta_{particle} - #eta_{parton} vs p_{T particle}",
   200,0,120,200,-2,2); 
  fhDeltaEtaAwayParton->SetYTitle("#Delta #eta");
  fhDeltaEtaAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaEtaAwayParton);
  
  fhDeltaPtAwayParton  = new TH2F
  ("DeltaPtAwayParton","#p_{T particle} - #p_{T parton} vs p_{T particle}",
   200,0,120,100,-10,10); 
  fhDeltaPtAwayParton->SetYTitle("#Delta #p_{T}");
  fhDeltaPtAwayParton->SetXTitle("p_{T particle} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtAwayParton);
  
  fhPtRatAwayParton  = new TH2F
  ("PtRatAwayParton","#p_{T parton} / #p_{T particle} vs p_{T particle}",
   200,0,120,200,0,5); 
  fhPtRatAwayParton->SetYTitle("ratio");
  fhPtRatAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhPtRatAwayParton);
  
  return outputContainer;
}

//____________________________________________________
// Initialize the parameters of the analysis.
//____________________________________________________
void AliAnaParticlePartonCorrelation::InitParameters()
{
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("Partons");  
  AddToHistogramsName("AnaPartonCorr_");
}

//_____________________________________________________________________
// Print some relevant parameters set for the analysis.
//_____________________________________________________________________
void AliAnaParticlePartonCorrelation::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
} 

//__________________________________________________________
/// Particle-Parton Correlation Analysis, create AODs.
/// Add partons to the reference list of the trigger particle,
/// Partons are considered those in the first eight possitions in the stack
/// being 0, and 1 the 2 protons, and 6 and 7 the outgoing final partons.
//__________________________________________________________
void  AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD()
{
  if(!GetInputAODBranch())
    AliFatal(Form("No input particles in AOD with name branch < %s > ",GetInputAODName().Data()));
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s>",
             GetInputAODBranch()->GetClass()->GetName()));
  
  AliDebug(1,"Begin fill AODs");
  AliDebug(1,Form("In particle branch aod entries %d", GetInputAODBranch()->GetEntriesFast()));
  
  //Loop on stored AOD particles
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(!GetMC())
    {
      AliFatal("No Stack available, STOP");
      return; // coverity
    }
    
    if(GetMC()->GetNumberOfTracks() < 8)
    {
      AliWarning(Form("*** small number of particles, not a PYTHIA simulation? ***:  n tracks %d", GetMC()->GetNumberOfPrimaries()));
      continue ;
    }
    
    //Fill AOD reference only with partons
    
    //Array with reference to partons, initialize
    TObjArray * objarray  = NULL;
    Int_t nrefs = 0;
    
    AliVParticle * parton    = NULL ;
    for(Int_t ipr = 0;ipr < 8; ipr ++ )
    {
      parton = GetMC()->GetTrack(ipr) ;
      nrefs++;
      
      if(nrefs==1)
      {
        objarray = new TObjArray(0);
        objarray->SetName(GetAODObjArrayName());
        objarray->SetOwner(kFALSE);
      }
      objarray->Add(parton);
    }//parton loop
    
    if(objarray->GetEntriesFast() > 0) particle->AddObjArray(objarray);
    
  }//Aod branch loop
  
  AliDebug(1,"End fill AODs");
}

//_________________________________________________________________
// Particle-Parton Correlation Analysis, fill histograms.
//_________________________________________________________________
void  AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms()
{
  if(!GetInputAODBranch())
  {
    AliFatal(Form("No input particles in AOD with name branch < %s >",GetInputAODName().Data()));
    return; //coverity
  }
  
  AliDebug(1,"Begin parton correlation analysis, fill histograms");
  AliDebug(1,Form("In particle branch aod entries %d", GetInputAODBranch()->GetEntriesFast()));
  
  if(!GetMC())
  {
    AliFatal("No Stack available, STOP");
    return;// coverity
  }
  
  // Loop on stored AOD particles
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  AliVParticle *  mom = NULL ;
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    Float_t ptTrigg  = particle->Pt();
    Float_t phiTrigg = particle->Phi();
    Float_t etaTrigg = particle->Eta(); 
    Int_t imom = particle->GetLabel();
    Int_t iparent  = 2000;
    Int_t iawayparent = -1;
    
    TObjArray * objarray = particle->GetObjArray(GetAODObjArrayName());
    if(!(objarray) || (objarray->GetEntriesFast() < 7) )
    {
      AliFatal("Reference list with partons not filled, STOP analysis");
      return; // coverity
    }
    
    // Check and get indeces of mother and parton
    if      (imom < 8 ) 
      iparent = imom ;   //mother is already a parton
    else if (imom <  GetMC()->GetNumberOfTracks()) 
    {
      mom =  GetMC()->GetTrack(imom);
      if(mom)
      {
        iparent=mom->GetMother();
        //cout<<" iparent "<<iparent<<endl;
        while(iparent > 7 )
        {
          mom = GetMC()->GetTrack(iparent);
          if (mom)
          {
            imom = iparent ; //Mother label is of the inmediate parton daughter
            iparent = mom->GetMother();
          }
          else iparent = -1;
          //cout<<" while iparent "<<iparent<<endl;
        } 
      }
    }
    
    AliDebug(1,Form("N reference partons %d; labels:  mother %d, parent %d", objarray->GetEntriesFast(), imom, iparent));
    
    if(iparent < 0 || iparent > 8)
    {
      AliWarning(Form("Failed to find appropriate parton, index %d", iparent));
      continue ;
    }
    
    // Near parton is the parton that fragmented and created the mother
    AliVParticle * nearParton = (AliVParticle*) objarray->At(iparent);
    Float_t  ptNearParton    = nearParton->Pt();
    Float_t  phiNearParton   = nearParton->Phi() ;
    Float_t  etaNearParton   = nearParton->Eta() ;
    
    fhDeltaEtaNearParton->Fill(ptTrigg, etaTrigg-etaNearParton, GetEventWeight());
    fhDeltaPhiNearParton->Fill(ptTrigg, phiTrigg-phiNearParton, GetEventWeight());
    fhDeltaPtNearParton ->Fill(ptTrigg, ptTrigg-ptNearParton  , GetEventWeight());
    fhPtRatNearParton   ->Fill(ptTrigg, ptNearParton/ptTrigg   , GetEventWeight());
    
    if     (iparent == 7) iawayparent = 6;
    else if(iparent == 6) iawayparent = 7;
    else
    {
      AliWarning("Parent parton is not final state, skip");
      continue;
    }
    
    // Away parton is the other final parton.
    AliVParticle * awayParton = (AliVParticle*) objarray->At(iawayparent);
    Float_t  ptAwayParton    = awayParton->Pt();
    Float_t  phiAwayParton   = awayParton->Phi() ;
    Float_t  etaAwayParton   = awayParton->Eta() ;
    fhDeltaEtaAwayParton->Fill(ptTrigg, etaTrigg-etaAwayParton, GetEventWeight());
    fhDeltaPhiAwayParton->Fill(ptTrigg, phiTrigg-phiAwayParton, GetEventWeight());
    fhDeltaPtAwayParton ->Fill(ptTrigg, ptTrigg-ptAwayParton  , GetEventWeight());
    fhPtRatAwayParton   ->Fill(ptTrigg, ptAwayParton/ptTrigg  , GetEventWeight());
    
  }
  
  AliDebug(1,"End fill histograms");
} 
