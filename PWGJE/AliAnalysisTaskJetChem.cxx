/*************************************************************************
 *                                                                       *
 *                                                                       *
 *      Task for Jet Chemistry Analysis in PWG4 Jet Task Force Train     *
 *                                                                       *
 *                                                                       *
 *************************************************************************/

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

/* $Id: */

#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "THnSparse.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h" 
#include "AliAODInputHandler.h" 
#include "AliESDEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"

#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"


#include "AliPID.h" 
#include "AliExternalTrackParam.h"

#include "AliAnalysisTaskJetChem.h"

ClassImp(AliAnalysisTaskJetChem)

//____________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem()
   : AliAnalysisTaskFragmentationFunction()
   ,fK0Type(0)
   ,fFilterMaskK0(0)
   ,fListK0s(0)
   ,fV0QAK0(0)
   ,fFFHistosRecCutsK0Evt(0)      
   ,fFFHistosIMK0AllEvt(0)        
   ,fFFHistosIMK0Jet(0)           
   ,fFFHistosIMK0Cone(0)
   ,fFFHistosPhiCorrIMK0(0)
   ,fFFIMNBinsJetPt(0)    
   ,fFFIMJetPtMin(0) 
   ,fFFIMJetPtMax(0)
   ,fFFIMNBinsInvM(0) 
   ,fFFIMInvMMin(0)   
   ,fFFIMInvMMax(0)   
   ,fFFIMNBinsPt(0)      
   ,fFFIMPtMin(0)        
   ,fFFIMPtMax(0)        
   ,fFFIMNBinsXi(0)      
   ,fFFIMXiMin(0)        
   ,fFFIMXiMax(0)        
   ,fFFIMNBinsZ(0)       
   ,fFFIMZMin(0)         
   ,fFFIMZMax(0)
   ,fPhiCorrIMNBinsPt(0)
   ,fPhiCorrIMPtMin(0)
   ,fPhiCorrIMPtMax(0)
   ,fPhiCorrIMNBinsPhi(0)
   ,fPhiCorrIMPhiMin(0)
   ,fPhiCorrIMPhiMax(0)
   ,fPhiCorrIMNBinsInvM(0)
   ,fPhiCorrIMInvMMin(0)
   ,fPhiCorrIMInvMMax(0)
   ,fh1K0Mult(0)
   ,fh1dPhiJetK0(0)
{
   // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem(const char *name) 
  : AliAnalysisTaskFragmentationFunction(name)
  ,fK0Type(0)  
  ,fFilterMaskK0(0)
  ,fListK0s(0)
  ,fV0QAK0(0)
  ,fFFHistosRecCutsK0Evt(0)      
  ,fFFHistosIMK0AllEvt(0)        
  ,fFFHistosIMK0Jet(0)           
  ,fFFHistosIMK0Cone(0)
  ,fFFHistosPhiCorrIMK0(0)
  ,fFFIMNBinsJetPt(0)    
  ,fFFIMJetPtMin(0) 
  ,fFFIMJetPtMax(0)
  ,fFFIMNBinsInvM(0) 
  ,fFFIMInvMMin(0)   
  ,fFFIMInvMMax(0)   
  ,fFFIMNBinsPt(0)      
  ,fFFIMPtMin(0)        
  ,fFFIMPtMax(0)        
  ,fFFIMNBinsXi(0)      
  ,fFFIMXiMin(0)        
  ,fFFIMXiMax(0)        
  ,fFFIMNBinsZ(0)       
  ,fFFIMZMin(0)         
  ,fFFIMZMax(0)
  ,fPhiCorrIMNBinsPt(0)
  ,fPhiCorrIMPtMin(0)
  ,fPhiCorrIMPtMax(0)
  ,fPhiCorrIMNBinsPhi(0)
  ,fPhiCorrIMPhiMin(0)
  ,fPhiCorrIMPhiMax(0)
  ,fPhiCorrIMNBinsInvM(0)
  ,fPhiCorrIMInvMMin(0)
  ,fPhiCorrIMInvMMax(0)
  ,fh1K0Mult(0)
  ,fh1dPhiJetK0(0)     
{
  // constructor
  
  DefineOutput(1,TList::Class());  
}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem(const  AliAnalysisTaskJetChem &copy)
  : AliAnalysisTaskFragmentationFunction()
  ,fK0Type(copy.fK0Type)
  ,fFilterMaskK0(copy.fFilterMaskK0)
  ,fListK0s(copy.fListK0s)
  ,fV0QAK0(copy.fV0QAK0)
  ,fFFHistosRecCutsK0Evt(copy.fFFHistosRecCutsK0Evt)      
  ,fFFHistosIMK0AllEvt(copy.fFFHistosIMK0AllEvt)        
  ,fFFHistosIMK0Jet(copy.fFFHistosIMK0Jet)           
  ,fFFHistosIMK0Cone(copy.fFFHistosIMK0Cone)          
  ,fFFHistosPhiCorrIMK0(copy.fFFHistosPhiCorrIMK0)          
  ,fFFIMNBinsJetPt(copy.fFFIMNBinsJetPt)   
  ,fFFIMJetPtMin(copy.fFFIMJetPtMin)     
  ,fFFIMJetPtMax(copy.fFFIMJetPtMax)     
  ,fFFIMNBinsInvM(copy.fFFIMNBinsInvM)  
  ,fFFIMInvMMin(copy.fFFIMInvMMin)    
  ,fFFIMInvMMax(copy.fFFIMInvMMax)    
  ,fFFIMNBinsPt(copy.fFFIMNBinsPt)      
  ,fFFIMPtMin(copy.fFFIMPtMin)        
  ,fFFIMPtMax(copy.fFFIMPtMax)        
  ,fFFIMNBinsXi(copy.fFFIMNBinsXi)      
  ,fFFIMXiMin(copy.fFFIMXiMin)        
  ,fFFIMXiMax(copy.fFFIMXiMax)        
  ,fFFIMNBinsZ(copy.fFFIMNBinsZ)       
  ,fFFIMZMin(copy.fFFIMZMin)         
  ,fFFIMZMax(copy.fFFIMZMax) 
  ,fPhiCorrIMNBinsPt(copy.fPhiCorrIMNBinsPt)
  ,fPhiCorrIMPtMin(copy.fPhiCorrIMPtMin)
  ,fPhiCorrIMPtMax(copy.fPhiCorrIMPtMax)
  ,fPhiCorrIMNBinsPhi(copy.fPhiCorrIMNBinsPhi)
  ,fPhiCorrIMPhiMin(copy.fPhiCorrIMPhiMin)
  ,fPhiCorrIMPhiMax(copy.fPhiCorrIMPhiMax)
  ,fPhiCorrIMNBinsInvM(copy.fPhiCorrIMNBinsInvM)
  ,fPhiCorrIMInvMMin(copy.fPhiCorrIMInvMMin)
  ,fPhiCorrIMInvMMax(copy.fPhiCorrIMInvMMax)
  ,fh1K0Mult(copy.fh1K0Mult)
  ,fh1dPhiJetK0(copy.fh1dPhiJetK0)
{
  // copy constructor
  
}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem& AliAnalysisTaskJetChem::operator=(const AliAnalysisTaskJetChem& o)
{
  // assignment
  
  if(this!=&o){
    AliAnalysisTaskFragmentationFunction::operator=(o);

    fK0Type                         = o.fK0Type;
    fFilterMaskK0                   = o.fFilterMaskK0;
    fListK0s                        = o.fListK0s;
    fV0QAK0                         = o.fV0QAK0;
    fFFHistosRecCutsK0Evt           = o.fFFHistosRecCutsK0Evt;      
    fFFHistosIMK0AllEvt             = o.fFFHistosIMK0AllEvt;        
    fFFHistosIMK0Jet                = o.fFFHistosIMK0Jet;           
    fFFHistosIMK0Cone               = o.fFFHistosIMK0Cone;          
    fFFHistosPhiCorrIMK0            = o.fFFHistosPhiCorrIMK0;
    fFFIMNBinsJetPt                 = o.fFFIMNBinsJetPt;    
    fFFIMJetPtMin                   = o.fFFIMJetPtMin; 
    fFFIMJetPtMax                   = o.fFFIMJetPtMax;
    fFFIMNBinsPt                    = o.fFFIMNBinsPt;      
    fFFIMPtMin                      = o.fFFIMPtMin;        
    fFFIMPtMax                      = o.fFFIMPtMax;        
    fFFIMNBinsXi                    = o.fFFIMNBinsXi;      
    fFFIMXiMin                      = o.fFFIMXiMin;        
    fFFIMXiMax                      = o.fFFIMXiMax;        
    fFFIMNBinsZ                     = o.fFFIMNBinsZ;       
    fFFIMZMin                       = o.fFFIMZMin;         
    fFFIMZMax                       = o.fFFIMZMax;         
    fPhiCorrIMNBinsPt               = o.fPhiCorrIMNBinsPt;
    fPhiCorrIMPtMin                 = o.fPhiCorrIMPtMin;
    fPhiCorrIMPtMax                 = o.fPhiCorrIMPtMax;
    fPhiCorrIMNBinsPhi              = o.fPhiCorrIMNBinsPhi;
    fPhiCorrIMPhiMin                = o.fPhiCorrIMPhiMin;
    fPhiCorrIMPhiMax                = o.fPhiCorrIMPhiMax;
    fPhiCorrIMNBinsInvM             = o.fPhiCorrIMNBinsInvM;
    fPhiCorrIMInvMMin               = o.fPhiCorrIMInvMMin;
    fPhiCorrIMInvMMax               = o.fPhiCorrIMInvMMax;
    fh1K0Mult                       = o.fh1K0Mult;
    fh1dPhiJetK0                    = o.fh1dPhiJetK0;
  }
    
  return *this;
}

//_______________________________________________
AliAnalysisTaskJetChem::~AliAnalysisTaskJetChem()
{
  // destructor  

  if(fListK0s) delete fListK0s;
}

//________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AliFragFuncHistosInvMass(const char* name, 
									   Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
									   Int_t nInvMass, Float_t invMassMin, Float_t invMassMax,
									   Int_t nPt, Float_t ptMin, Float_t ptMax,
									   Int_t nXi, Float_t xiMin, Float_t xiMax,
									   Int_t nZ , Float_t zMin , Float_t zMax )
  : TObject()
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsInvMass(nInvMass)
  ,fInvMassMin(invMassMin)  
  ,fInvMassMax(invMassMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsXi(nXi) 
  ,fXiMin(xiMin)   
  ,fXiMax(xiMax)   
  ,fNBinsZ(nZ)  
  ,fZMin(zMin)    
  ,fZMax(zMax)    
  ,fh3TrackPt(0)
  ,fh3Xi(0)
  ,fh3Z(0)
  ,fh1JetPt(0)
  ,fNameFF(name)
{
  // default constructor

}

//______________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AliFragFuncHistosInvMass(const AliFragFuncHistosInvMass& copy)
  : TObject()
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsInvMass(copy.fNBinsInvMass)
  ,fInvMassMin(copy.fInvMassMin)  
  ,fInvMassMax(copy.fInvMassMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   
  ,fPtMax(copy.fPtMax)   
  ,fNBinsXi(copy.fNBinsXi) 
  ,fXiMin(copy.fXiMin)   
  ,fXiMax(copy.fXiMax)   
  ,fNBinsZ(copy.fNBinsZ)  
  ,fZMin(copy.fZMin)    
  ,fZMax(copy.fZMax)    
  ,fh3TrackPt(copy.fh3TrackPt)
  ,fh3Xi(copy.fh3Xi)
  ,fh3Z(copy.fh3Z)
  ,fh1JetPt(copy.fh1JetPt)
  ,fNameFF(copy.fNameFF)
{
  // copy constructor
}

//______________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass& AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::operator=(const AliAnalysisTaskJetChem::AliFragFuncHistosInvMass& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsJetPt   = o.fNBinsJetPt;
    fJetPtMin     = o.fJetPtMin;
    fJetPtMax     = o.fJetPtMax;
    fNBinsInvMass = o.fNBinsInvMass;
    fInvMassMin   = o.fInvMassMin;  
    fInvMassMax   = o.fInvMassMax;
    fNBinsPt      = o.fNBinsPt; 
    fPtMin        = o.fPtMin;   
    fPtMax        = o.fPtMax;   
    fNBinsXi      = o.fNBinsXi; 
    fXiMin        = o.fXiMin;   
    fXiMax        = o.fXiMax;   
    fNBinsZ       = o.fNBinsZ;  
    fZMin         = o.fZMin;    
    fZMax         = o.fZMax;    
    fh3TrackPt    = o.fh3TrackPt;
    fh3Xi         = o.fh3Xi;
    fh3Z          = o.fh3Z;
    fh1JetPt      = o.fh1JetPt;
    fNameFF       = o.fNameFF;
  }
    
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::~AliFragFuncHistosInvMass()
{ 
  // destructor 

  if(fh1JetPt)   delete fh1JetPt;
  if(fh3TrackPt) delete fh3TrackPt;
  if(fh3Xi)      delete fh3Xi;
  if(fh3Z)       delete fh3Z;
}

//_________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::DefineHistos()
{
  // book FF histos

  fh1JetPt   = new TH1F(Form("fh1FFJetPtIM%s", fNameFF.Data()),"",fNBinsJetPt,fJetPtMin,fJetPtMax);
  fh3TrackPt = new TH3F(Form("fh3FFTrackPtIM%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsPt, fPtMin, fPtMax);
  fh3Xi      = new TH3F(Form("fh3FFXiIM%s", fNameFF.Data()),"", fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsXi, fXiMin, fXiMax);
  fh3Z       = new TH3F(Form("fh3FFZIM%s", fNameFF.Data()),"", fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsZ, fZMin, fZMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{t} (GeV/c)", "entries"); 
  AliAnalysisTaskJetChem::SetProperties(fh3TrackPt,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","p_{t} (GeV/c)");
  AliAnalysisTaskJetChem::SetProperties(fh3Xi,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","#xi");
  AliAnalysisTaskJetChem::SetProperties(fh3Z,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","z");
}

//________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::FillFF(Float_t trackPt, Float_t invM, Float_t jetPt, Bool_t incrementJetPt)
{
  // fill FF
 
  if(incrementJetPt) fh1JetPt->Fill(jetPt);    
  fh3TrackPt->Fill(jetPt,invM,trackPt);
  
  Double_t z = 0.;
  if(jetPt>0) z = trackPt / jetPt;
  Double_t xi = 0;
  if(z>0) xi = TMath::Log(1/z);
  
  fh3Xi->Fill(jetPt,invM,xi);
  fh3Z->Fill(jetPt,invM,z);
}

//___________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  
  list->Add(fh3TrackPt);
  list->Add(fh3Xi);
  list->Add(fh3Z);
}

// ---


//_______________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AliFragFuncHistosPhiCorrInvMass(const char* name,
											 Int_t nPt,   Float_t ptMin,   Float_t ptMax,
											 Int_t nPhi,  Float_t phiMin,  Float_t phiMax,
											 Int_t nInvMass,  Float_t invMassMin,  Float_t invMassMax)
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fNBinsInvMass(nInvMass)
  ,fInvMassMin(invMassMin)  
  ,fInvMassMax(invMassMax)
  ,fh3PhiCorr(0) 
  ,fNamePhiCorr(name) 
{
  // default constructor
}

//____________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AliFragFuncHistosPhiCorrInvMass(const AliFragFuncHistosPhiCorrInvMass& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fNBinsInvMass(copy.fNBinsInvMass)
  ,fInvMassMin(copy.fInvMassMin)  
  ,fInvMassMax(copy.fInvMassMax)
  ,fh3PhiCorr(copy.fh3PhiCorr)
  ,fNamePhiCorr(copy.fNamePhiCorr)
{
  // copy constructor
}

//________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass& AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::operator=(const AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt      = o.fNBinsPt;
    fPtMin        = o.fPtMin;
    fPtMax        = o.fPtMax;
    fNBinsPhi     = o.fNBinsPhi;
    fPhiMin       = o.fPhiMin;
    fPhiMax       = o.fPhiMax;
    fNBinsInvMass = o.fNBinsInvMass;
    fInvMassMin   = o.fInvMassMin;  
    fInvMassMax   = o.fInvMassMax;
    
    fh3PhiCorr    = o.fh3PhiCorr;
    fNamePhiCorr  = o.fNamePhiCorr;
  }
  
  return *this;
}

//_________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::~AliFragFuncHistosPhiCorrInvMass()
{
  // destructor 
  
  if(fh3PhiCorr) delete fh3PhiCorr;
}

//__________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::DefineHistos()
{
  // book jet QA histos

  fh3PhiCorr  = new TH3F(Form("fh3PhiCorrIM%s", fNamePhiCorr.Data()), 
			 Form("%s: p_{t} - #phi - m_{inv} distribution",fNamePhiCorr.Data()), 
			 fNBinsPt, fPtMin, fPtMax, 
			 fNBinsPhi, fPhiMin, fPhiMax,
			 fNBinsInvMass, fInvMassMin, fInvMassMax);
  
  AliAnalysisTaskJetChem::SetProperties(fh3PhiCorr, "p_{t} (GeV/c)", "#phi", "m_{inv} (GeV/c^2)"); 
}

//___________________________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::FillPhiCorr(Float_t pt, Float_t phi, Float_t invM)
{
  // fill jet QA histos 

  fh3PhiCorr->Fill(pt, phi, invM);
}

//______________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AddToOutput(TList* list) const 
{
  // add histos to list

  list->Add(fh3PhiCorr);
}

//____________________________________________________
void AliAnalysisTaskJetChem::UserCreateOutputObjects()
{
  // create output objects

  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()");
 
  // create list of tracks and jets 
  
  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE);  

  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);

  fListK0s = new TList(); 
  fListK0s->SetOwner(kFALSE);

  //
  // Create histograms / output container
  //

  OpenFile(1);
  fCommonHistList = new TList();
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  	

  // histograms inherited from AliAnalysisTaskFragmentationFunction

  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1EvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fh1EvtSelection->GetXaxis()->SetBinLabel(2,"event selection: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(3,"event class: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(4,"vertex Ncontr: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(5,"vertex z: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(6,"vertex type: rejected");

  fh1EvtCent 	             = new TH1F("fh1EvtCent","centrality",100,0.,100.);
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 11,-.5, 10.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1Xsec                    = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Trials                  = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1PtHard                  = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fh1PtHardTrials            = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);

  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",10,-0.5,9.5);


  // histograms jetChem proper 

  fh1EvtMult 	             = new TH1F("fh1EvtMult","multiplicity",1200,0.,12000.);
  fh1K0Mult 	             = new TH1F("fh1K0Mult","K0 multiplicity",500,0.,500.);
  fh1dPhiJetK0               = new TH1F("fh1dPhiJetK0","",640,-1,5.4);



  fFFHistosRecCuts   	     = new AliFragFuncHistos("RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  fV0QAK0                   = new AliFragFuncQATrackHistos("V0QAK0",fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							    fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							    fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							    fQATrackHighPtThreshold);
  
  fFFHistosRecCutsK0Evt      = new AliFragFuncHistos("RecCutsK0Evt", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  
  fFFHistosIMK0AllEvt        = new AliFragFuncHistosInvMass("K0AllEvt", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
  
  fFFHistosIMK0Jet           = new AliFragFuncHistosInvMass("K0Jet", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
  
  
  fFFHistosIMK0Cone          = new AliFragFuncHistosInvMass("K0Cone", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
  
  fFFHistosPhiCorrIMK0       = new AliFragFuncHistosPhiCorrInvMass("K0",fPhiCorrIMNBinsPt, fPhiCorrIMPtMin, fPhiCorrIMPtMax, 
								   fPhiCorrIMNBinsPhi, fPhiCorrIMPhiMin, fPhiCorrIMPhiMax,  
								   fPhiCorrIMNBinsInvM , fPhiCorrIMInvMMin , fPhiCorrIMInvMMax);
  
  
  fV0QAK0->DefineHistos();
  fFFHistosRecCuts->DefineHistos();
  fFFHistosRecCutsK0Evt->DefineHistos();
  fFFHistosIMK0AllEvt->DefineHistos();
  fFFHistosIMK0Jet->DefineHistos();
  fFFHistosIMK0Cone->DefineHistos();
  fFFHistosPhiCorrIMK0->DefineHistos();

  const Int_t saveLevel = 5;
  if(saveLevel>0){
    
    fCommonHistList->Add(fh1EvtSelection);
    fCommonHistList->Add(fh1EvtCent);
    fCommonHistList->Add(fh1VertexNContributors);
    fCommonHistList->Add(fh1VertexZ);
    fCommonHistList->Add(fh1Xsec);
    fCommonHistList->Add(fh1Trials);
    fCommonHistList->Add(fh1PtHard);
    fCommonHistList->Add(fh1PtHardTrials);
    fCommonHistList->Add(fh1nRecJetsCuts);
    fCommonHistList->Add(fh1EvtMult);
    fCommonHistList->Add(fh1K0Mult);
    fCommonHistList->Add(fh1dPhiJetK0);

    fV0QAK0->AddToOutput(fCommonHistList);
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecCutsK0Evt->AddToOutput(fCommonHistList);

    fFFHistosIMK0AllEvt->AddToOutput(fCommonHistList);
    fFFHistosIMK0Jet->AddToOutput(fCommonHistList);
    fFFHistosIMK0Cone->AddToOutput(fCommonHistList);
    fFFHistosPhiCorrIMK0->AddToOutput(fCommonHistList);
  }
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fCommonHistList->GetEntries(); ++i){
    TH1 *h1 = dynamic_cast<TH1*>(fCommonHistList->At(i));
    if (h1) h1->Sumw2();
    else{
      THnSparse *hnSparse = dynamic_cast<THnSparse*>(fCommonHistList->At(i));
      if(hnSparse) hnSparse->Sumw2();
    }
  }
  
  TH1::AddDirectory(oldStatus);
}

//_______________________________________________
void AliAnalysisTaskJetChem::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);

  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(inputHandler->IsEventSelected() & AliVEvent::kMB){
    if(fDebug > 1)  Printf(" Trigger Selection: event ACCEPTED ... ");
    fh1EvtSelection->Fill(1.);
  } else {
    fh1EvtSelection->Fill(0.);
    if(inputHandler->InheritsFrom("AliESDInputHandler") && fUsePhysicsSelection){ // PhysicsSelection only with ESD input
      if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
      PostData(1, fCommonHistList);
      return;
    }
  }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD){
    if(fDebug>3) Printf("%s:%d ESDEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  fMCEvent = MCEvent();
  if(!fMCEvent){
    if(fDebug>3) Printf("%s:%d MCEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  // get AOD event from input/ouput
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
    fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
    if(fUseAODInputJets) fAODJets = fAOD;
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      fAODJets = fAOD;
      if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
    }
  }
  
  if(!fAODJets && !fUseAODInputJets){ // case we have AOD in input & output and want jets from output
    TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
      fAODJets = ((AliAODHandler*)outHandler)->GetAOD();
      if (fDebug > 1)  Printf("%s:%d jets from output AOD", (char*)__FILE__,__LINE__);
    }
  }
  
  if(fNonStdFile.Length()!=0){
    // case we have an AOD extension - fetch the jets from the extended output
    
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension not found for %s",fNonStdFile.Data());
    }
  }
  
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return;
  }
  if(!fAODJets){
    Printf("%s:%d AODEvent with jet branch not found", (char*)__FILE__,__LINE__);
    return;
  }
  

  // event selection  *****************************************
  
  Double_t centPercent = -1;
  if(fEventClass>0){
    Int_t cl = 0;
    if(handler && handler->InheritsFrom("AliAODInputHandler")){ 
      // since it is not supported by the helper task define own classes
      centPercent = fAOD->GetHeader()->GetCentrality();
      cl = 1;
      if(centPercent>10) cl = 2;
      if(centPercent>30) cl = 3;
      if(centPercent>50) cl = 4;
    }
    else {
      cl = AliAnalysisHelperJetTasks::EventClass();
      if(fESD) centPercent = fESD->GetCentrality()->GetCentralityPercentile("V0M"); // OB added
    }
    
    if(cl!=fEventClass){
      // event not in selected event class, reject event
      if (fDebug > 1) Printf("%s:%d event not in selected event class: event REJECTED ...",(char*)__FILE__,__LINE__);
      fh1EvtSelection->Fill(2.);
      PostData(1, fCommonHistList);
      return;
    }
  }
  
  // *** vertex cut ***
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(!nTracksPrim){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return;
  }
  
  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>fMaxVertexZ){
    if (fDebug > 1) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ()); 
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return; 
  }
  
  TString primVtxName(primVtx->GetName());
  
  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(5.);
    PostData(1, fCommonHistList);
    return;
  }
  
  if (fDebug > 1) Printf("%s:%d event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fh1EvtSelection->Fill(0.);
  fh1EvtCent->Fill(centPercent);
    

  //___ get MC information __________________________________________________________________

  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  
  if(fMCEvent){
     AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
     AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
     AliGenHijingEventHeader*  hijingGenHeader = 0x0;

     if(pythiaGenHeader){
	 if(fDebug>3) Printf("%s:%d pythiaGenHeader found", (char*)__FILE__,__LINE__);
	 nTrials = pythiaGenHeader->Trials();
	 ptHard  = pythiaGenHeader->GetPtHard();

	 fh1PtHard->Fill(ptHard);
	 fh1PtHardTrials->Fill(ptHard,nTrials);


     } else { // no pythia, hijing?

	 if(fDebug>3) Printf("%s:%d no pythiaGenHeader found", (char*)__FILE__,__LINE__);

         hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
         if(!hijingGenHeader){
            Printf("%s:%d no pythiaGenHeader or hjingGenHeader found", (char*)__FILE__,__LINE__);
         } else {
	    if(fDebug>3) Printf("%s:%d hijingGenHeader found", (char*)__FILE__,__LINE__);
	 }
     }

     fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
  }


  //____ fetch jets ______________________________________________________________

  Int_t nJCuts = GetListOfJets(fJetsRecCuts, kJetsRecAcceptance);
  Int_t nRecJetsCuts = 0;
  if(nJCuts>=0) nRecJetsCuts = fJetsRecCuts->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  if(nRecJetsCuts != nJCuts) Printf("%s:%d Mismatch selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  fh1nRecJetsCuts->Fill(nRecJetsCuts);
  

  //____ fetch particles __________________________________________________________
 
  Int_t nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());
  if(fTracksRecCuts->GetEntries() != nTCuts) 
    Printf("%s:%d Mismatch selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());
  fh1EvtMult->Fill(fTracksRecCuts->GetEntries());

  Int_t nK0s = GetListOfK0s(fListK0s,fK0Type);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  if(nK0s != fListK0s->GetEntries()) Printf("%s:%d Mismatch selected K0s: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  fh1K0Mult->Fill(fListK0s->GetEntries());

   // ___ V0 QA + K0 pt spectra all events _______________________________________________
  
  if(fListK0s->GetEntries()>0){
    
    for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 
      
      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
      if(!v0) continue;
    
      Float_t trackPt       = v0->Pt();
      Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
      Double_t invM         = v0->MassK0Short();
      Double_t jetPt        = fFFIMJetPtMin; // assign pro forma jet energy

      fV0QAK0->FillTrackQA(v0->Eta(), TVector2::Phi_0_2pi(v0->Phi()), v0->Pt()); 
      fFFHistosIMK0AllEvt->FillFF(trackPt, invM, jetPt, incrementJetPt);
    }
  }


  //____ fill FF histos  __________________________________________________________

  for(Int_t ij=0; ij<nRecJetsCuts; ++ij){

    AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));
    Double_t jetPt = jet->Pt();

    if(ij==0){ // leading jet
      
      TList* jettracklist = new TList();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
      
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt,  GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
      }

      if(GetFFMinNTracks()>0 && jettracklist->GetSize() <= GetFFMinNTracks()) isBadJet = kTRUE;
      if(isBadJet) continue; 

      
      for(Int_t it=0; it<jettracklist->GetSize(); ++it){

	AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));	
	if(!trackVP)continue;

	Float_t trackPt = trackVP->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	fFFHistosRecCuts->FillFF(trackPt, jetPt, incrementJetPt);
	if(nK0s>0) fFFHistosRecCutsK0Evt->FillFF(trackPt, jetPt, incrementJetPt);
      }
      
      delete jettracklist;

      // ---- K0s ---- 
      
      // fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );

      for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 

	AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
	if(!v0) continue;
 
	Double_t v0Mom[3];
	v0->PxPyPz(v0Mom);
	TVector3 v0MomVect(v0Mom);

	Double_t dPhiJetK0 = (jet->MomentumVector()->Vect()).DeltaPhi(v0MomVect);
	
	Float_t trackPt       = v0->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	Double_t invM         = v0->MassK0Short();
	
	fFFHistosIMK0Jet->FillFF(trackPt, invM, jetPt, incrementJetPt);
 	fFFHistosPhiCorrIMK0->FillPhiCorr(trackPt,TVector2::Phi_0_2pi(dPhiJetK0),invM);

	if(dPhiJetK0<fh1dPhiJetK0->GetXaxis()->GetXmin()) dPhiJetK0 += 2*TMath::Pi();
	fh1dPhiJetK0->Fill(dPhiJetK0);

      }

      if(fListK0s->GetSize() == 0){ // no K0: increment jet pt spectrum 

	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMK0Jet->FillFF(-1, -1, jetPt, incrementJetPt);
      }
      

      TList* jetConeK0list = new TList();
      Double_t sumPtK0     = 0.;
      Bool_t isBadJetK0    = kFALSE; // dummy, do not use

      GetJetTracksPointing(fListK0s, jetConeK0list, jet, GetFFRadius(), sumPtK0, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0);


      if(fDebug>2)Printf("%s:%d nK0s total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetConeK0list->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetConeK0list->GetSize(); ++it){ // loop K0s in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeK0list->At(it));
	if(!v0) continue;

	Double_t invM           = v0->MassK0Short();
	Float_t  trackPt        = v0->Pt();
	Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	//std::cout<<" trackPt "<<trackPt<<" invM "<<invM<<std::endl;
	fFFHistosIMK0Cone->FillFF(trackPt, invM, jetPt, incrementJetPt);
      }
     if(jetConeK0list->GetSize() == 0){ // no K0: increment jet pt spectrum 

	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMK0Cone->FillFF(-1, -1, jetPt, incrementJetPt);
      }

      delete jetConeK0list;
    }
  }
  
  fTracksRecCuts->Clear();
  fJetsRecCuts->Clear();
  fListK0s->Clear();

  //Post output data.
  PostData(1, fCommonHistList);    
}

// ____________________________________________________________________________________________
void AliAnalysisTaskJetChem::SetProperties(TH3F* h,const char* x, const char* y, const char* z)
{
  //Set properties of histos (x,y and z title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->SetZTitle(z);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->GetZaxis()->SetTitleColor(1);
}

// ____________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetChem::IsAccepteddEdx(const Double_t mom,const Double_t signal, AliPID::EParticleType n, const Double_t cutnSig) const{

  // apply TPC dE/dx cut similar as in AliTPCpidESD 
  // note: AliTPCidESD uses inner track param for momentum - not avaiable on AOD,  
  //       so we use global track momentum 
  // should use separate parametrisation for MC and data, but probably ALEPH param & 7% resolution used here anyway not the last word 
 
 
  const Double_t kBBMIP(50.);
  const Double_t kBBRes(0.07);
  //const Double_t kBBRange(5.);
  const Double_t kBBp1(0.76176e-1);
  const Double_t kBBp2(10.632);
  const Double_t kBBp3(0.13279e-4);
  const Double_t kBBp4(1.8631);
  const Double_t kBBp5(1.9479);

  Double_t mass=AliPID::ParticleMass(n); 
  Double_t betaGamma = mom/mass;

  const Float_t kmeanCorrection =0.1;
  Double_t bb = AliExternalTrackParam::BetheBlochAleph(betaGamma,kBBp1,kBBp2,kBBp3,kBBp4,kBBp5);
  Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
  Double_t bethe = bb * meanCorrection; // expected
  Double_t sigma = bethe * kBBRes;
        

  Double_t dedx = signal/kBBMIP; // measured

  Double_t nSig = (TMath::Abs(dedx - bethe))/sigma;
  
  if(nSig > cutnSig) return kFALSE; 

  return kTRUE;
}

//___________________________________________________________________
Bool_t AliAnalysisTaskJetChem::IsK0InvMass(const Double_t mass) const
{
  // K0 mass ? Use FF histo limits
  
  if(fFFIMInvMMin <= mass && mass <= fFFIMInvMMax) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Int_t AliAnalysisTaskJetChem::GetListOfK0s(TList *list, const Int_t type)
{
  // fill list of V0s selected according to type

  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }

  if(type==kTrackUndef) return 0;

  for(int i=0; i<fAOD->GetNumberOfV0s(); i++){ // loop over V0s
    
    AliAODv0* v0 = fAOD->GetV0(i);

    Bool_t isOnFly = v0->GetOnFlyStatus();
    
    if(!isOnFly &&  (type == kOnFly || type == kOnFlyPID || type == kOnFlydEdx || type == kOnFlyPrim)) continue; 
    if( isOnFly &&  (type == kOffl  || type == kOfflPID  || type == kOffldEdx  || type == kOfflPrim))  continue; 

    Double_t massK0 = v0->MassK0Short();

    if(!(IsK0InvMass(massK0))) continue; // moved invMass cut for HI - otherwise too slow

    if(type == kOnFlyPID || type == kOfflPID){
      
      AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0)); // slow 
      AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1)); // slow  
     
      // AOD pid - cuts strongly into signal

      AliAODTrack::AODTrkPID_t mpPIDNeg = trackNeg->GetMostProbablePID();
      AliAODTrack::AODTrkPID_t mpPIDPos = trackPos->GetMostProbablePID();
	
      if(!( (mpPIDNeg == AliAODTrack::kPion) && (mpPIDPos == AliAODTrack::kPion) ) ) continue;
    }
   
    if(type == kOnFlydEdx || type == kOffldEdx){

      AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0)); // slow 
      AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1)); // slow  

      AliAODPid*  aodPidPos = trackPos->GetDetPid();
      AliAODPid*  aodPidNeg = trackNeg->GetDetPid();

      Double_t  dEdxPos = aodPidPos->GetTPCsignal();
      Double_t  dEdxNeg = aodPidNeg->GetTPCsignal();

      Double_t momPos  = trackPos->P();
      Double_t momNeg  = trackNeg->P();

      Int_t cutnSigdEdx = 2;
      if(! (IsAccepteddEdx(momPos,dEdxPos,AliPID::kPion,cutnSigdEdx)) && (IsAccepteddEdx(momNeg,dEdxNeg,AliPID::kPion,cutnSigdEdx)) ) continue;
      
    }   
    
    if(type == kOnFlyPrim || type == kOfflPrim){
      
      AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0)); // slow 
      AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1)); // slow  
      
      //std::cout<<"  filer map trackPos "<<trackPos->GetFilterMap()<<" trackNeg "<<trackNeg->GetFilterMap()<<std::endl;
      
      if((fFilterMaskK0>0) && !(trackPos->TestFilterBit(fFilterMaskK0)))   continue;
      if((fFilterMaskK0>0) && !(trackNeg->TestFilterBit(fFilterMaskK0)))   continue;
    }
    
    list->Add(v0);
  }
  
  Int_t nK0s = list->GetSize();
  
  return nK0s;
}


