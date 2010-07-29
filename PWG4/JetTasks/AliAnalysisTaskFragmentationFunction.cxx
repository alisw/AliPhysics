/*************************************************************************
 *                                                                       *
 * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
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


#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include "AliAODInputHandler.h" 
#include "AliAODHandler.h" 
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliAnalysisTaskFragmentationFunction.h"


ClassImp(AliAnalysisTaskFragmentationFunction)

//____________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction()
   : AliAnalysisTaskSE()
   ,fESD(0)
   ,fAOD(0)
   ,fMCEvent(0)
   ,fBranchRecJets("jets")
   ,fBranchGenJets("")
   ,fTrackTypeGen(0)
   ,fJetTypeGen(0)
   ,fFilterMask(0)
   ,fTrackPtCut(0)
   ,fTrackEtaMin(0)
   ,fTrackEtaMax(0)
   ,fTrackPhiMin(0)
   ,fTrackPhiMax(0)
   ,fJetPtCut(0)
   ,fJetEtaMin(0)
   ,fJetEtaMax(0)
   ,fJetPhiMin(0)
   ,fJetPhiMax(0)
   ,fFFRadius(0)
   ,fDijetDeltaPhiCut(0)
   ,fDijetInvMassMin(0)
   ,fDijetInvMassMax(0)
   ,fDijetCDFcut(0)
   ,fDijetEMeanMin(0)
   ,fDijetEMeanMax(0)
   ,fDijetEFractionCut(0)
   ,fDijetEFraction(0)
   ,fTracksRec(0)
   ,fTracksRecCuts(0)
   ,fTracksGen(0)
   ,fJetsRec(0)
   ,fJetsRecCuts(0)
   ,fJetsGen(0)
   ,fQATrackHistosRec(0)
   ,fQATrackHistosRecCuts(0)
   ,fQATrackHistosGen(0)
   ,fQAJetHistosRec(0)
   ,fQAJetHistosRecCuts(0)
   ,fQAJetHistosRecCutsLeading(0)
   ,fQAJetHistosGen(0)
   ,fQAJetHistosGenLeading(0)
   ,fFFHistosRecCuts(0)
   ,fFFHistosRecLeading(0)
   ,fFFHistosRecLeadingTrack(0)
   ,fFFHistosGen(0)
   ,fFFHistosGenLeading(0)
   ,fFFHistosGenLeadingTrack(0)
   ,fQATrackHighPtThreshold(0)
   ,fFFNBinsJetPt(0)    
   ,fFFJetPtMin(0) 
   ,fFFJetPtMax(0)
   ,fFFNBinsPt(0)      
   ,fFFPtMin(0)        
   ,fFFPtMax(0)        
   ,fFFNBinsXi(0)      
   ,fFFXiMin(0)        
   ,fFFXiMax(0)        
   ,fFFNBinsZ(0)       
   ,fFFZMin(0)         
   ,fFFZMax(0)         
   ,fQAJetNBinsPt(0)   
   ,fQAJetPtMin(0)     
   ,fQAJetPtMax(0)     
   ,fQAJetNBinsEta(0)  
   ,fQAJetEtaMin(0)    
   ,fQAJetEtaMax(0)    
   ,fQAJetNBinsPhi(0)  
   ,fQAJetPhiMin(0)    
   ,fQAJetPhiMax(0)    
   ,fQATrackNBinsPt(0) 
   ,fQATrackPtMin(0)   
   ,fQATrackPtMax(0)   
   ,fQATrackNBinsEta(0)
   ,fQATrackEtaMin(0)  
   ,fQATrackEtaMax(0)  
   ,fQATrackNBinsPhi(0)
   ,fQATrackPhiMin(0)  
   ,fQATrackPhiMax(0)  
   ,fCommonHistList(0)
   ,fh1EvtSelection(0)
   ,fh1VertexNContributors(0)
   ,fh1VertexZ(0)
   ,fh1EvtMult(0)
   ,fh1nRecJetsCuts(0)
   ,fh1nGenJets(0)
{
   // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fMCEvent(0)
  ,fBranchRecJets("jets")
  ,fBranchGenJets("")
  ,fTrackTypeGen(0)
  ,fJetTypeGen(0)
  ,fFilterMask(0)
  ,fTrackPtCut(0)
  ,fTrackEtaMin(0)
  ,fTrackEtaMax(0)
  ,fTrackPhiMin(0)
  ,fTrackPhiMax(0)
  ,fJetPtCut(0)
  ,fJetEtaMin(0)
  ,fJetEtaMax(0)
  ,fJetPhiMin(0)
  ,fJetPhiMax(0)
  ,fFFRadius(0)
  ,fDijetDeltaPhiCut(0)
  ,fDijetInvMassMin(0)
  ,fDijetInvMassMax(0)
  ,fDijetCDFcut(0)
  ,fDijetEMeanMin(0)
  ,fDijetEMeanMax(0)
  ,fDijetEFractionCut(0)
  ,fDijetEFraction(0)
  ,fTracksRec(0)
  ,fTracksRecCuts(0)
  ,fTracksGen(0)
  ,fJetsRec(0)
  ,fJetsRecCuts(0)
  ,fJetsGen(0)
  ,fQATrackHistosRec(0)
  ,fQATrackHistosRecCuts(0)
  ,fQATrackHistosGen(0)
  ,fQAJetHistosRec(0)
  ,fQAJetHistosRecCuts(0)
  ,fQAJetHistosRecCutsLeading(0)
  ,fQAJetHistosGen(0)
  ,fQAJetHistosGenLeading(0)
  ,fFFHistosRecCuts(0)
  ,fFFHistosRecLeading(0)
  ,fFFHistosRecLeadingTrack(0)
  ,fFFHistosGen(0)
  ,fFFHistosGenLeading(0)
  ,fFFHistosGenLeadingTrack(0)
  ,fQATrackHighPtThreshold(0) 
  ,fFFNBinsJetPt(0)    
  ,fFFJetPtMin(0) 
  ,fFFJetPtMax(0)
  ,fFFNBinsPt(0)      
  ,fFFPtMin(0)        
  ,fFFPtMax(0)        
  ,fFFNBinsXi(0)      
  ,fFFXiMin(0)        
  ,fFFXiMax(0)        
  ,fFFNBinsZ(0)       
  ,fFFZMin(0)         
  ,fFFZMax(0)         
  ,fQAJetNBinsPt(0)   
  ,fQAJetPtMin(0)     
  ,fQAJetPtMax(0)     
  ,fQAJetNBinsEta(0)  
  ,fQAJetEtaMin(0)    
  ,fQAJetEtaMax(0)    
  ,fQAJetNBinsPhi(0)  
  ,fQAJetPhiMin(0)    
  ,fQAJetPhiMax(0)    
  ,fQATrackNBinsPt(0) 
  ,fQATrackPtMin(0)   
  ,fQATrackPtMax(0)   
  ,fQATrackNBinsEta(0)
  ,fQATrackEtaMin(0)  
  ,fQATrackEtaMax(0)  
  ,fQATrackNBinsPhi(0)
  ,fQATrackPhiMin(0)  
  ,fQATrackPhiMax(0)  
  ,fCommonHistList(0)
  ,fh1EvtSelection(0)
  ,fh1VertexNContributors(0)
  ,fh1VertexZ(0)
  ,fh1EvtMult(0)
  ,fh1nRecJetsCuts(0)
  ,fh1nGenJets(0)
{
  // constructor
  
  DefineOutput(1,TList::Class());
  

}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const  AliAnalysisTaskFragmentationFunction &copy)
  : AliAnalysisTaskSE()
  ,fESD(copy.fESD)
  ,fAOD(copy.fAOD)
  ,fMCEvent(copy.fMCEvent)
  ,fBranchRecJets(copy.fBranchRecJets)
  ,fBranchGenJets(copy.fBranchGenJets)
  ,fTrackTypeGen(copy.fTrackTypeGen)
  ,fJetTypeGen(copy.fJetTypeGen)
  ,fFilterMask(copy.fFilterMask)
  ,fTrackPtCut(copy.fTrackPtCut)
  ,fTrackEtaMin(copy.fTrackEtaMin)
  ,fTrackEtaMax(copy.fTrackEtaMax)
  ,fTrackPhiMin(copy.fTrackPhiMin)
  ,fTrackPhiMax(copy.fTrackPhiMax)
  ,fJetPtCut(copy.fJetPtCut)
  ,fJetEtaMin(copy.fJetEtaMin)
  ,fJetEtaMax(copy.fJetEtaMax)
  ,fJetPhiMin(copy.fJetPhiMin)
  ,fJetPhiMax(copy.fJetPhiMax)
  ,fFFRadius(copy.fFFRadius)
  ,fDijetDeltaPhiCut(copy.fDijetDeltaPhiCut)
  ,fDijetInvMassMin(copy.fDijetInvMassMin)
  ,fDijetInvMassMax(copy.fDijetInvMassMax)
  ,fDijetCDFcut(copy.fDijetCDFcut)
  ,fDijetEMeanMin(copy.fDijetEMeanMin)
  ,fDijetEMeanMax(copy.fDijetEMeanMax)
  ,fDijetEFractionCut(copy.fDijetEFractionCut)
  ,fDijetEFraction(copy.fDijetEFraction)
  ,fTracksRec(copy.fTracksRec)
  ,fTracksRecCuts(copy.fTracksRecCuts)
  ,fTracksGen(copy.fTracksGen)
  ,fJetsRec(copy.fJetsRec)
  ,fJetsRecCuts(copy.fJetsRecCuts)
  ,fJetsGen(copy.fJetsGen)
  ,fQATrackHistosRec(copy.fQATrackHistosRec)
  ,fQATrackHistosRecCuts(copy.fQATrackHistosRecCuts)
  ,fQATrackHistosGen(copy.fQATrackHistosGen)
  ,fQAJetHistosRec(copy.fQAJetHistosRec)
  ,fQAJetHistosRecCuts(copy.fQAJetHistosRecCuts)
  ,fQAJetHistosRecCutsLeading(copy.fQAJetHistosRecCutsLeading)
  ,fQAJetHistosGen(copy.fQAJetHistosGen)
  ,fQAJetHistosGenLeading(copy.fQAJetHistosGenLeading)
  ,fFFHistosRecCuts(copy.fFFHistosRecCuts)
  ,fFFHistosRecLeading(copy.fFFHistosRecLeading)
  ,fFFHistosRecLeadingTrack(copy.fFFHistosRecLeadingTrack)
  ,fFFHistosGen(copy.fFFHistosGen)
  ,fFFHistosGenLeading(copy.fFFHistosGenLeading)
  ,fFFHistosGenLeadingTrack(copy.fFFHistosGenLeadingTrack)
  ,fQATrackHighPtThreshold(copy.fQATrackHighPtThreshold) 
  ,fFFNBinsJetPt(copy.fFFNBinsJetPt)    
  ,fFFJetPtMin(copy.fFFJetPtMin) 
  ,fFFJetPtMax(copy.fFFJetPtMax)
  ,fFFNBinsPt(copy.fFFNBinsPt)      
  ,fFFPtMin(copy.fFFPtMin)        
  ,fFFPtMax(copy.fFFPtMax)        
  ,fFFNBinsXi(copy.fFFNBinsXi)      
  ,fFFXiMin(copy.fFFXiMin)        
  ,fFFXiMax(copy.fFFXiMax)        
  ,fFFNBinsZ(copy.fFFNBinsZ)       
  ,fFFZMin(copy.fFFZMin)         
  ,fFFZMax(copy.fFFZMax)         
  ,fQAJetNBinsPt(copy.fQAJetNBinsPt)   
  ,fQAJetPtMin(copy.fQAJetPtMin)     
  ,fQAJetPtMax(copy.fQAJetPtMax)     
  ,fQAJetNBinsEta(copy.fQAJetNBinsEta)  
  ,fQAJetEtaMin(copy.fQAJetEtaMin)    
  ,fQAJetEtaMax(copy.fQAJetEtaMax)    
  ,fQAJetNBinsPhi(copy.fQAJetNBinsPhi)  
  ,fQAJetPhiMin(copy.fQAJetPhiMin)    
  ,fQAJetPhiMax(copy.fQAJetPhiMax)    
  ,fQATrackNBinsPt(copy.fQATrackNBinsPt) 
  ,fQATrackPtMin(copy.fQATrackPtMin)   
  ,fQATrackPtMax(copy.fQATrackPtMax)   
  ,fQATrackNBinsEta(copy.fQATrackNBinsEta)
  ,fQATrackEtaMin(copy.fQATrackEtaMin)  
  ,fQATrackEtaMax(copy.fQATrackEtaMax)  
  ,fQATrackNBinsPhi(copy.fQATrackNBinsPhi)
  ,fQATrackPhiMin(copy.fQATrackPhiMin)  
  ,fQATrackPhiMax(copy.fQATrackPhiMax)  
  ,fCommonHistList(copy.fCommonHistList)
  ,fh1EvtSelection(copy.fh1EvtSelection)
  ,fh1VertexNContributors(copy.fh1VertexNContributors)
  ,fh1VertexZ(copy.fh1VertexZ)
  ,fh1EvtMult(copy.fh1EvtMult)
  ,fh1nRecJetsCuts(copy.fh1nRecJetsCuts)
  ,fh1nGenJets(copy.fh1nGenJets)
{
  // copy constructor

}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction& AliAnalysisTaskFragmentationFunction::operator=(const AliAnalysisTaskFragmentationFunction& o)
{
  // assignment
  
  if(this!=&o){

    AliAnalysisTaskSE::operator=(o);
    fESD                       = o.fESD;
    fAOD                       = o.fAOD;
    fMCEvent                   = o.fMCEvent;
    fBranchRecJets             = o.fBranchRecJets;
    fBranchGenJets             = o.fBranchGenJets;
    fTrackTypeGen              = o.fTrackTypeGen;
    fJetTypeGen                = o.fJetTypeGen;
    fFilterMask                = o.fFilterMask;
    fTrackPtCut                = o.fTrackPtCut;
    fTrackEtaMin               = o.fTrackEtaMin;
    fTrackEtaMax               = o.fTrackEtaMax;
    fTrackPhiMin               = o.fTrackPhiMin;
    fTrackPhiMax               = o.fTrackPhiMax;
    fJetPtCut                  = o.fJetPtCut;
    fJetEtaMin                 = o.fJetEtaMin;
    fJetEtaMax                 = o.fJetEtaMax;
    fJetPhiMin                 = o.fJetPhiMin;
    fJetPhiMax                 = o.fJetPhiMax;
    fFFRadius                  = o.fFFRadius;
    fDijetDeltaPhiCut          = o.fDijetDeltaPhiCut;
    fDijetInvMassMin           = o.fDijetInvMassMin;
    fDijetInvMassMax           = o.fDijetInvMassMax;
    fDijetCDFcut               = o.fDijetCDFcut;
    fDijetEMeanMin             = o.fDijetEMeanMin;
    fDijetEMeanMax             = o.fDijetEMeanMax;
    fDijetEFractionCut         = o.fDijetEFractionCut;
    fDijetEFraction            = o.fDijetEFraction;
    fTracksRec                 = o.fTracksRec;
    fTracksRecCuts             = o.fTracksRecCuts;
    fTracksGen                 = o.fTracksGen;
    fJetsRec                   = o.fJetsRec;
    fJetsRecCuts               = o.fJetsRecCuts;
    fJetsGen                   = o.fJetsGen;
    fQATrackHistosRec          = o.fQATrackHistosRec;
    fQATrackHistosRecCuts      = o.fQATrackHistosRecCuts;
    fQATrackHistosGen          = o.fQATrackHistosGen;
    fQAJetHistosRec            = o.fQAJetHistosRec;
    fQAJetHistosRecCuts        = o.fQAJetHistosRecCuts;
    fQAJetHistosRecCutsLeading = o.fQAJetHistosRecCutsLeading;
    fQAJetHistosGen            = o.fQAJetHistosGen;
    fQAJetHistosGenLeading     = o.fQAJetHistosGenLeading;
    fFFHistosRecCuts           = o.fFFHistosRecCuts;
    fFFHistosRecLeading        = o.fFFHistosRecLeading;
    fFFHistosRecLeadingTrack   = o.fFFHistosRecLeadingTrack;
    fFFHistosGen               = o.fFFHistosGen;
    fFFHistosGenLeading        = o.fFFHistosGenLeading;
    fFFHistosGenLeadingTrack   = o.fFFHistosGenLeadingTrack;
    fQATrackHighPtThreshold    = o.fQATrackHighPtThreshold; 
    fFFNBinsJetPt              = o.fFFNBinsJetPt;    
    fFFJetPtMin                = o.fFFJetPtMin; 
    fFFJetPtMax                = o.fFFJetPtMax;
    fFFNBinsPt                 = o.fFFNBinsPt;      
    fFFPtMin                   = o.fFFPtMin;        
    fFFPtMax                   = o.fFFPtMax;        
    fFFNBinsXi                 = o.fFFNBinsXi;      
    fFFXiMin                   = o.fFFXiMin;        
    fFFXiMax                   = o.fFFXiMax;        
    fFFNBinsZ                  = o.fFFNBinsZ;       
    fFFZMin                    = o.fFFZMin;         
    fFFZMax                    = o.fFFZMax;         
    fQAJetNBinsPt              = o.fQAJetNBinsPt;   
    fQAJetPtMin                = o.fQAJetPtMin;     
    fQAJetPtMax                = o.fQAJetPtMax;     
    fQAJetNBinsEta             = o.fQAJetNBinsEta;  
    fQAJetEtaMin               = o.fQAJetEtaMin;    
    fQAJetEtaMax               = o.fQAJetEtaMax;    
    fQAJetNBinsPhi             = o.fQAJetNBinsPhi;  
    fQAJetPhiMin               = o.fQAJetPhiMin;    
    fQAJetPhiMax               = o.fQAJetPhiMax;    
    fQATrackNBinsPt            = o.fQATrackNBinsPt; 
    fQATrackPtMin              = o.fQATrackPtMin;   
    fQATrackPtMax              = o.fQATrackPtMax;   
    fQATrackNBinsEta           = o.fQATrackNBinsEta;
    fQATrackEtaMin             = o.fQATrackEtaMin;  
    fQATrackEtaMax             = o.fQATrackEtaMax;  
    fQATrackNBinsPhi           = o.fQATrackNBinsPhi;
    fQATrackPhiMin             = o.fQATrackPhiMin;  
    fQATrackPhiMax             = o.fQATrackPhiMax;  
    fCommonHistList            = o.fCommonHistList;
    fh1EvtSelection            = o.fh1EvtSelection;
    fh1VertexNContributors     = o.fh1VertexNContributors;
    fh1VertexZ                 = o.fh1VertexZ;
    fh1EvtMult                 = o.fh1EvtMult;
    fh1nRecJetsCuts            = o.fh1nRecJetsCuts;
    fh1nGenJets                = o.fh1nGenJets; 
  }
    
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::~AliAnalysisTaskFragmentationFunction()
{
  // destructor
  
}



//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AliFragFuncHistos(const char* name, 
							 Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
							 Int_t nPt, Float_t ptMin, Float_t ptMax,
							 Int_t nXi, Float_t xiMin, Float_t xiMax,
							 Int_t nZ , Float_t zMin , Float_t zMax )
  : TObject()
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsXi(nXi) 
  ,fXiMin(xiMin)   
  ,fXiMax(xiMax)   
  ,fNBinsZ(nZ)  
  ,fZMin(zMin)    
  ,fZMax(zMax)    
  ,fh2TrackPt(0)
  ,fh2Xi(0)
  ,fh2Z(0)
  ,fh1JetPt(0)
  ,fName(name)
{
  // default constructor

}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AliFragFuncHistos(const AliFragFuncHistos& copy)
  : TObject()
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   
  ,fPtMax(copy.fPtMax)   
  ,fNBinsXi(copy.fNBinsXi) 
  ,fXiMin(copy.fXiMin)   
  ,fXiMax(copy.fXiMax)   
  ,fNBinsZ(copy.fNBinsZ)  
  ,fZMin(copy.fZMin)    
  ,fZMax(copy.fZMax)    
  ,fh2TrackPt(copy.fh2TrackPt)
  ,fh2Xi(copy.fh2Xi)
  ,fh2Z(copy.fh2Z)
  ,fh1JetPt(copy.fh1JetPt)
  ,fName(copy.fName)
{
  // copy constructor
}

//_______________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsJetPt = o.fNBinsJetPt;
    fJetPtMin   = o.fJetPtMin;
    fJetPtMax   = o.fJetPtMax;
    fNBinsPt    = o.fNBinsPt; 
    fPtMin      = o.fPtMin;   
    fPtMax      = o.fPtMax;   
    fNBinsXi    = o.fNBinsXi; 
    fXiMin      = o.fXiMin;   
    fXiMax      = o.fXiMax;   
    fNBinsZ     = o.fNBinsZ;  
    fZMin       = o.fZMin;    
    fZMax       = o.fZMax;    
    fh2TrackPt  = o.fh2TrackPt;
    fh2Xi       = o.fh2Xi;
    fh2Z        = o.fh2Z;
    fh1JetPt    = o.fh1JetPt;
    fName       = o.fName;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::~AliFragFuncHistos()
{
  // destructor 

  if(fh1JetPt)   delete fh1JetPt;
  if(fh2TrackPt) delete fh2TrackPt;
  if(fh2Xi)      delete fh2Xi;
  if(fh2Z)       delete fh2Z;
}

//_________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::DefineHistos()
{
  // book FF histos

  fh1JetPt   = new TH1F(Form("fh1FFJetPt%s", fName.Data()),"",fNBinsJetPt,fJetPtMin,fJetPtMax);
  fh2TrackPt = new TH2F(Form("fh2FFTrackPt%s",fName.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsPt, fPtMin, fPtMax);
  fh2Xi      = new TH2F(Form("fh2FFXi%s",fName.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsXi, fXiMin, fXiMax);
  fh2Z       = new TH2F(Form("fh2FFZ%s",fName.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsZ, fZMin, fZMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{T} [GeV/c]", "entries"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPt,"jet p_{T} [GeV/c]","p_{T} [GeV/c]","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi,"jet p_{T} [GeV/c]","#xi", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z,"jet p_{T} [GeV/c]","z","entries");
}

//_______________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::FillFF(Float_t trackPt, Float_t jetPt, Bool_t incrementJetPt)
{
  // fill FF
 
  if(incrementJetPt) fh1JetPt->Fill(jetPt);    
  fh2TrackPt->Fill(jetPt,trackPt);
  
  Double_t z = trackPt / jetPt;
  Double_t xi = 0;
  if(z>0) xi = TMath::Log(1/z);
  
  fh2Xi->Fill(jetPt,xi);
  fh2Z->Fill(jetPt,z);
}

//_________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  
  list->Add(fh2TrackPt);
  list->Add(fh2Xi);
  list->Add(fh2Z);
}


//_________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const char* name,
							       Int_t nPt,   Float_t ptMin,   Float_t ptMax,
							       Int_t nEta,  Float_t etaMin,  Float_t etaMax,
							       Int_t nPhi,  Float_t phiMin,  Float_t phiMax)
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsEta(nEta)
  ,fEtaMin(etaMin)
  ,fEtaMax(etaMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fh2EtaPhi(0)
  ,fh1Pt(0)
  ,fName(name)
{
  // default constructor
}

//____________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const AliFragFuncQAJetHistos& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsEta(copy.fNBinsEta)
  ,fEtaMin(copy.fEtaMin)
  ,fEtaMax(copy.fEtaMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fh2EtaPhi(copy.fh2EtaPhi)
  ,fh1Pt(copy.fh1Pt)
  ,fName(copy.fName)
{
  // copy constructor
}

//________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt  = o.fNBinsPt;
    fPtMin    = o.fPtMin;
    fPtMax    = o.fPtMax;
    fNBinsEta = o.fNBinsEta;
    fEtaMin   = o.fEtaMin;
    fEtaMax   = o.fEtaMax;
    fNBinsPhi = o.fNBinsPhi;
    fPhiMin   = o.fPhiMin;
    fPhiMax   = o.fPhiMax;
    fh2EtaPhi = o.fh2EtaPhi;
    fh1Pt     = o.fh1Pt;
    fName     = o.fName;
  }
  
  return *this;
}

//______________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::~AliFragFuncQAJetHistos()
{
  // destructor 
  
  if(fh2EtaPhi) delete fh2EtaPhi;
  if(fh1Pt)     delete fh1Pt;
}

//____________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::DefineHistos()
{
  // book jet QA histos

  fh2EtaPhi  = new TH2F(Form("fh2JetQAEtaPhi%s", fName.Data()), Form("%s: #eta - #phi distribution", fName.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt      = new TH1F(Form("fh1JetQAPt%s", fName.Data()), Form("%s: p_{T} distribution", fName.Data()), fNBinsPt, fPtMin, fPtMax);
	
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
}

//____________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::FillJetQA(Float_t eta, Float_t phi, Float_t pt)
{
  // fill jet QA histos 

  fh2EtaPhi->Fill( eta, phi);
  fh1Pt->Fill( pt );
}

//____________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AddToOutput(TList* list) const 
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh1Pt);
}


//___________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const char* name,
								   Int_t nPt, Float_t ptMin, Float_t ptMax,
								   Int_t nEta, Float_t etaMin, Float_t etaMax,
								   Int_t nPhi, Float_t phiMin, Float_t phiMax,
								   Float_t ptThresh) 
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsEta(nEta)
  ,fEtaMin(etaMin)
  ,fEtaMax(etaMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fHighPtThreshold(ptThresh)
  ,fh2EtaPhi(0)
  ,fh1Pt(0)
  ,fh2HighPtEtaPhi(0)
  ,fName(name)
{
  // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const AliFragFuncQATrackHistos& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsEta(copy.fNBinsEta)
  ,fEtaMin(copy.fEtaMin)
  ,fEtaMax(copy.fEtaMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fHighPtThreshold(copy.fHighPtThreshold)
  ,fh2EtaPhi(copy.fh2EtaPhi)
  ,fh1Pt(copy.fh1Pt)
  ,fh2HighPtEtaPhi(copy.fh2HighPtEtaPhi)
  ,fName(copy.fName)
{
  // copy constructor
}

// _____________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt         = o.fNBinsPt;
    fPtMin           = o.fPtMin;
    fPtMax           = o.fPtMax;
    fNBinsEta        = o.fNBinsEta;
    fEtaMin          = o.fEtaMin;
    fEtaMax          = o.fEtaMax;
    fNBinsPhi        = o.fNBinsPhi;
    fPhiMin          = o.fPhiMin;
    fPhiMax          = o.fPhiMax;
    fHighPtThreshold = o.fHighPtThreshold;
    fh2EtaPhi        = o.fh2EtaPhi;
    fh1Pt            = o.fh1Pt;
    fh2HighPtEtaPhi  = o.fh2HighPtEtaPhi;
    fName            = o.fName;
  }
  
  return *this;
}

//___________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::~AliFragFuncQATrackHistos()
{
  // destructor 
  
  if(fh2EtaPhi)       delete fh2EtaPhi;
  if(fh2HighPtEtaPhi) delete fh2HighPtEtaPhi;
  if(fh1Pt)           delete fh1Pt;
}

//______________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::DefineHistos()
{
  // book track QA histos

  fh2EtaPhi       = new TH2F(Form("fh2TrackQAEtaPhi%s", fName.Data()), Form("%s: #eta - #phi distribution", fName.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh2HighPtEtaPhi = new TH2F(Form("fh2TrackQAHighPtEtaPhi%s", fName.Data()), Form("%s: #eta - #phi distribution for high-p_{T}", fName.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt           = new TH1F(Form("fh1TrackQAPt%s", fName.Data()), Form("%s: p_{T} distribution", fName.Data()), fNBinsPt, fPtMin, fPtMax);
  
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2HighPtEtaPhi, "#eta", "#phi");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
}

//________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::FillTrackQA(Float_t eta, Float_t phi, Float_t pt)
{
  // fill track QA histos
    
  fh2EtaPhi->Fill( eta, phi);
  if(pt > fHighPtThreshold) fh2HighPtEtaPhi->Fill( eta, phi);
  fh1Pt->Fill( pt );	
}

//______________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh2HighPtEtaPhi);
  list->Add(fh1Pt);
}

//__________________________________________________________________
void AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()
{
  // create output objects
  
  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()");
 
  // create list of tracks and jets 
  
  fTracksRec = new TList();
  fTracksRec->SetOwner(kFALSE);  

  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE);  

  fTracksGen = new TList();
  fTracksGen->SetOwner(kFALSE);

  fJetsRec = new TList();
  fJetsRec->SetOwner(kFALSE);

  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);

  fJetsGen = new TList();
  fJetsGen->SetOwner(kFALSE);

  // fJetsKine = new TList();
  // fJetsKine->SetOwner(kTRUE); // delete AOD jets using mom from Kine Tree via TList::Clear()


  //
  // Create histograms / output container
  //

  OpenFile(1);
  fCommonHistList = new TList();
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  
  // Histograms	
  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 11,-.5, 10.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1EvtMult 	             = new TH1F("fh1EvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",100,0.,100.);
  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",10,-0.5,9.5);
  fh1nGenJets                = new TH1F("fh1nGenJets","generated jets per event",10,-0.5,9.5);


  fQATrackHistosRec          = new AliFragFuncQATrackHistos("Rec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
						 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, fQATrackHighPtThreshold);
  fQATrackHistosRecCuts      = new AliFragFuncQATrackHistos("RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
						 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, fQATrackHighPtThreshold);
  fQATrackHistosGen          = new AliFragFuncQATrackHistos("Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
						 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, fQATrackHighPtThreshold);
  

  fQAJetHistosRec            = new AliFragFuncQAJetHistos("Rec", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
					       fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
  fQAJetHistosRecCuts        = new AliFragFuncQAJetHistos("RecCuts", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
					       fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
  fQAJetHistosRecCutsLeading = new AliFragFuncQAJetHistos("RecCutsLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
					       fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
  fQAJetHistosGen            = new AliFragFuncQAJetHistos("Gen", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
					       fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
  fQAJetHistosGenLeading     = new AliFragFuncQAJetHistos("GenLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
					       fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
  

  fFFHistosRecCuts   	     = new AliFragFuncHistos("RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);
  fFFHistosRecLeading        = new AliFragFuncHistos("RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);
  fFFHistosRecLeadingTrack   = new AliFragFuncHistos("RecLeadingTrack", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);
  fFFHistosGen   	     = new AliFragFuncHistos("Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);
  fFFHistosGenLeading        = new AliFragFuncHistos("GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);
  fFFHistosGenLeadingTrack   = new AliFragFuncHistos("GenLeadingTrack", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, fFFNBinsPt, fFFPtMin, fFFPtMax, fFFNBinsXi, fFFXiMin, fFFXiMax,  
					    fFFNBinsZ , fFFZMin , fFFZMax);


  fQATrackHistosRec->DefineHistos();
  fQATrackHistosRecCuts->DefineHistos();
  fQATrackHistosGen->DefineHistos();

  fQAJetHistosRec->DefineHistos();
  fQAJetHistosRecCuts->DefineHistos();
  fQAJetHistosRecCutsLeading->DefineHistos();
  fQAJetHistosGen->DefineHistos();
  fQAJetHistosGenLeading->DefineHistos();

  fFFHistosRecCuts->DefineHistos();
  fFFHistosRecLeading->DefineHistos();
  fFFHistosRecLeadingTrack->DefineHistos();
  fFFHistosGen->DefineHistos();
  fFFHistosGenLeading->DefineHistos();
  fFFHistosGenLeadingTrack->DefineHistos();
  
  
  const Int_t saveLevel = 4;
  if(saveLevel>0){
    fCommonHistList->Add(fh1EvtSelection);
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecLeading->AddToOutput(fCommonHistList);
    fFFHistosRecLeadingTrack->AddToOutput(fCommonHistList);
    fFFHistosGen->AddToOutput(fCommonHistList);
    fFFHistosGenLeading->AddToOutput(fCommonHistList);
    fFFHistosGenLeadingTrack->AddToOutput(fCommonHistList);
  }
  if(saveLevel>1){
    fQATrackHistosRec->AddToOutput(fCommonHistList);
    fQATrackHistosRecCuts->AddToOutput(fCommonHistList);
    fQATrackHistosGen->AddToOutput(fCommonHistList);
    
    fQAJetHistosRec->AddToOutput(fCommonHistList);
    fQAJetHistosRecCuts->AddToOutput(fCommonHistList);
    fQAJetHistosRecCutsLeading->AddToOutput(fCommonHistList);
    fQAJetHistosGen->AddToOutput(fCommonHistList);
    fQAJetHistosGenLeading->AddToOutput(fCommonHistList);
    
    fCommonHistList->Add(fh1EvtMult);
    fCommonHistList->Add(fh1nRecJetsCuts);
    fCommonHistList->Add(fh1nGenJets);
  }
  if(saveLevel>2){
    fCommonHistList->Add(fh1VertexNContributors);
    fCommonHistList->Add(fh1VertexZ);    
  }
  if(saveLevel>3){
    
  }
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fCommonHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fCommonHistList->At(i));
    if (h1) h1->Sumw2();	
  }
  
  TH1::AddDirectory(oldStatus);
}

//_______________________________________________
void AliAnalysisTaskFragmentationFunction::Init()
{
  // Initialization
  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::Init()");

}

//_____________________________________________________________
void AliAnalysisTaskFragmentationFunction::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::UserExec()");
	

  if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);
  // Trigger selection
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(inputHandler->IsEventSelected()&AliVEvent::kMB){
    if(fDebug > 1)  Printf(" Trigger Selection: event ACCEPTED ... ");
    fh1EvtSelection->Fill(1.);
  } else {
    fh1EvtSelection->Fill(0.);
    if(inputHandler->InheritsFrom("AliESDInputHandler")){ // PhysicsSelection only with ESD input
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
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD  = ((AliAODHandler*)handler)->GetAOD();
      if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
    }
  }
  
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return;
  }
  
  
  // event selection (vertex) *****************************************
  
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(!nTracksPrim){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fh1EvtSelection->Fill(2.);
    PostData(1, fCommonHistList);
    return;
  }

  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>10){
    if (fDebug > 1) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ()); 
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return; 
  }

  TString primVtxName(primVtx->GetName());

  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return;
  }
  if (fDebug > 1) Printf("%s:%d primary vertex selection: event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fh1EvtSelection->Fill(5.);
  
  
  //___ fetch jets __________________________________________________________________________
  
  Int_t nJ = GetListOfJets(fJetsRec, kJetsRec);
  Int_t nRecJets = 0;
  if(nJ>=0) nRecJets = fJetsRec->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec jets: %d %d",(char*)__FILE__,__LINE__,nJ,nRecJets);
  if(nJ != nRecJets) Printf("%s:%d Mismatch Selected Rec Jets: %d %d",(char*)__FILE__,__LINE__,nJ,nRecJets);

  Int_t nJCuts = GetListOfJets(fJetsRecCuts, kJetsRecAcceptance);
  Int_t nRecJetsCuts = 0;
  if(nJCuts>=0) nRecJetsCuts = fJetsRecCuts->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  if(nRecJetsCuts != nJCuts) Printf("%s:%d Mismatch selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  fh1nRecJetsCuts->Fill(nRecJetsCuts);

  
  if(fJetTypeGen==kJetsKine || fJetTypeGen == kJetsKineAcceptance) fJetsGen->SetOwner(kTRUE); // kine aod jets allocated on heap, delete them with TList::Clear() 
  Int_t nJGen  = GetListOfJets(fJetsGen, fJetTypeGen);
  Int_t nGenJets = 0;
  if(nJGen>=0) nGenJets = fJetsGen->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
  if(nJGen != nGenJets) Printf("%s:%d Mismatch selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
  fh1nGenJets->Fill(nGenJets);


  //____ fetch particles __________________________________________________________
  
  Int_t nT = GetListOfTracks(fTracksRec, kTrackAOD);
  Int_t nRecPart = 0;
  if(nT>=0) nRecPart = fTracksRec->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,nRecPart);
  if(nRecPart != nT) Printf("%s:%d Mismatch selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,nRecPart);
  

  Int_t nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);
  Int_t nRecPartCuts = 0;
  if(nTCuts>=0) nRecPartCuts = fTracksRecCuts->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,nRecPartCuts);
  if(nRecPartCuts != nTCuts) Printf("%s:%d Mismatch selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,nRecPartCuts);
  fh1EvtMult->Fill(nRecPartCuts);


  Int_t nTGen = GetListOfTracks(fTracksGen,fTrackTypeGen);
  Int_t nGenPart = 0;
  if(nTGen>=0) nGenPart = fTracksGen->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
  if(nGenPart != nTGen) Printf("%s:%d Mismatch selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
  
  
  //____ analysis, fill histos ___________________________________________________
  
  // loop over tracks

  for(Int_t it=0; it<nRecPart; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRec->At(it));
    fQATrackHistosRec->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  for(Int_t it=0; it<nRecPartCuts; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(it));
    fQATrackHistosRecCuts->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  for(Int_t it=0; it<nGenPart; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksGen->At(it));
    fQATrackHistosGen->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  
  // loop over jets

  for(Int_t ij=0; ij<nRecJets; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRec->At(ij));
    fQAJetHistosRec->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
  }
  

  for(Int_t ij=0; ij<nRecJetsCuts; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRecCuts->At(ij));
    fQAJetHistosRecCuts->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());

    if(ij==0){ // leading jet
      
      fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
      
      TList* jettracklist = new TList();
      Double_t sumPt = 0.;
      Float_t leadTrackPt = 0.;
      
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet);
       } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt);
      }
      
      for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	Float_t trackPt = (dynamic_cast<AliVParticle*> (jettracklist->At(it)))->Pt();
	Float_t jetPt = jet->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	fFFHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt);
	
	if(it==0){ 
	  leadTrackPt = trackPt;
	  fFFHistosRecLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE);
	}
	fFFHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt);
      }
      
      delete jettracklist;
    }
  }
	

  for(Int_t ij=0; ij<nGenJets; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsGen->At(ij));
    fQAJetHistosGen->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
    
    if(ij==0){ // leading jet
    
      fQAJetHistosGenLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      
      TList* jettracklist = new TList();
      Double_t sumPt = 0.;
      Float_t leadTrackPt = 0.;
      
      if(GetFFRadius()<=0){
	GetJetTracksTrackrefs(jettracklist, jet);
      } else {
	GetJetTracksPointing(fTracksGen, jettracklist, jet, GetFFRadius(), sumPt);
      }
      
      for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	Float_t trackPt = (dynamic_cast<AliVParticle*>(jettracklist->At(it)))->Pt();
	Float_t jetPt = jet->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	fFFHistosGen->FillFF( trackPt, jetPt, incrementJetPt);
	
	if(it==0){ 
	  leadTrackPt = trackPt;
	  fFFHistosGenLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE);
	}
	fFFHistosGenLeading->FillFF( trackPt, leadTrackPt, incrementJetPt);
      }
      
      delete jettracklist;
    }
  }
  
  fTracksRec->Clear();
  fTracksRecCuts->Clear();
  fTracksGen->Clear();
  fJetsRec->Clear();
  fJetsRecCuts->Clear();
  fJetsGen->Clear();

  //Post output data.
  PostData(1, fCommonHistList);
  
}

//______________________________________________________________
void AliAnalysisTaskFragmentationFunction::Terminate(Option_t *) 
{
  // terminated

  if(fDebug > 1) printf("AliAnalysisTaskFragmentationFunction::Terminate() \n");
}  

//_________________________________________________________________________________
Int_t AliAnalysisTaskFragmentationFunction::GetListOfTracks(TList *list, Int_t type)
{
  // fill list of tracks selected according to type

  if(fDebug > 2) Printf("%s:%d Selecting tracks with %d", (char*)__FILE__,__LINE__,type);
  
  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }

  if(type==kTrackUndef) return -1;
  
  Int_t iCount = 0;
  if(type==kTrackAODCuts || type==kTrackAOD){

    // all rec. tracks, esd filter mask, eta range
    if(!fAOD) return -1;
    
    for(Int_t it=0; it<fAOD->GetNumberOfTracks(); ++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      
      if(type == kTrackAODCuts){
	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))   continue;
	if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax) continue;
	if(tr->Phi() < fTrackPhiMin || tr->Phi() > fTrackPhiMax) continue;
	if(tr->Pt()  < fTrackPtCut) continue;
      }
      list->Add(tr);
      iCount++;
    }
  }
  else if (type==kTrackKineAll || type==kTrackKineCharged || type==kTrackKineChargedAcceptance){
    // kine particles, all or rather charged
    if(!fMCEvent) return iCount;
    
    for(Int_t it=0; it<fMCEvent->GetNumberOfTracks(); ++it){
      AliMCParticle* part = (AliMCParticle*) fMCEvent->GetTrack(it);
      
      if(type == kTrackKineCharged || type == kTrackKineChargedAcceptance){
	if(part->Charge()==0) continue;
	
	if(type == kTrackKineChargedAcceptance && 
	   (       part->Eta() < fTrackEtaMin
		|| part->Eta() > fTrackEtaMax
		|| part->Phi() < fTrackPhiMin
		|| part->Phi() > fTrackPhiMax 
		|| part->Pt()  < fTrackPtCut)) continue;
      }
      
      list->Add(part);
      iCount++;
    }
  }
  else if (type==kTrackAODMCCharged || type==kTrackAODMCAll || type==kTrackAODMCChargedAcceptance) {
    // MC particles (from AOD), physical primaries, all or rather charged or rather charged within acceptance
    if(!fAOD) return -1;
    
    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    
    for(int it=0; it<tca->GetEntriesFast(); ++it){
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
      if(!part->IsPhysicalPrimary())continue;
      
      if (type==kTrackAODMCCharged || type==kTrackAODMCChargedAcceptance){
	if(part->Charge()==0) continue;
	if(type==kTrackAODMCChargedAcceptance && 
	   (     part->Eta() > fTrackEtaMax
	      || part->Eta() < fTrackEtaMin
	      || part->Phi() > fTrackPhiMax
	      || part->Phi() < fTrackPhiMin
	      || part->Pt()  < fTrackPtCut)) continue;
      }
      
      list->Add(part);
      iCount++;
    }
  }
  
  list->Sort();
  return iCount;
  
}
// _______________________________________________________________________________
Int_t AliAnalysisTaskFragmentationFunction::GetListOfJets(TList *list, Int_t type)
{
  // fill list of jets selected according to type
  
  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }

  if(type == kJetsRec || type == kJetsRecAcceptance){ // reconstructed jets

    if(fBranchRecJets.Length()==0){
      Printf("%s:%d no rec jet branch specified", (char*)__FILE__,__LINE__);
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    TClonesArray *aodRecJets = 0x0;
    if(fBranchRecJets.Length()) aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRecJets.Data()));
    if(!aodRecJets)             aodRecJets = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(fBranchRecJets.Data()));

    if(!aodRecJets){
      if(fBranchRecJets.Length()) Printf("%s:%d no reconstructed jet array with name %s found", (char*)__FILE__,__LINE__,fBranchRecJets.Data());

      if(fDebug>1)fAOD->Print();
      return 0;
    }

    Int_t nRecJets = 0;
    
    for(Int_t ij=0; ij<aodRecJets->GetEntries(); ++ij){

      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ij));
      if(!tmp) continue;
	
      if( tmp->Pt() < fJetPtCut ) continue;
      if( type == kJetsRecAcceptance &&
	  (    tmp->Eta() < fJetEtaMin
	    || tmp->Eta() > fJetEtaMax
	    || tmp->Phi() < fJetPhiMin
	    || tmp->Phi() > fJetPhiMax )) continue;
      
      list->Add(tmp);
	  
      nRecJets++;
    }

    list->Sort();
    return nRecJets;
  }
  else if(type == kJetsKine || type == kJetsKineAcceptance){
    
    // generated jets
    Int_t nGenJets = 0;
    
    if(!fMCEvent){
      if(fDebug>1) Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
      return 0;
    }
    
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(fMCEvent);
    if(!pythiaGenHeader){
      Printf("%s:%d no pythiaGenHeader found", (char*)__FILE__,__LINE__);
      return 0;
    }
    
    // fetch the pythia generated jets
    for(int ip=0; ip<pythiaGenHeader->NTriggerJets(); ++ip){
      
      Float_t p[4];
      AliAODJet *jet = new AliAODJet();
      pythiaGenHeader->TriggerJet(ip, p);
      jet->SetPxPyPzE(p[0], p[1], p[2], p[3]);

      if( type == kJetsKineAcceptance &&
          (    jet->Eta() < fJetEtaMin
            || jet->Eta() > fJetEtaMax
            || jet->Phi() < fJetPhiMin
            || jet->Phi() > fJetPhiMax )) continue;
      
      list->Add(jet);
      nGenJets++;
    }
    list->Sort();
    return nGenJets;
  }
  else if(type == kJetsGen || type == kJetsGenAcceptance ){

    if(fBranchGenJets.Length()==0){
      if(fDebug>1) Printf("%s:%d no gen jet branch specified", (char*)__FILE__,__LINE__);
      return 0;
    }
    
    TClonesArray *aodGenJets = 0x0;
    if(fBranchGenJets.Length()) aodGenJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchGenJets.Data()));
    if(!aodGenJets)             aodGenJets = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(fBranchGenJets.Data()));

    if(!aodGenJets){
      if(fDebug>0){
	if(fBranchGenJets.Length())         Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGenJets.Data());
      }
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    Int_t nGenJets = 0;
    
    for(Int_t ig=0; ig<aodGenJets->GetEntries(); ++ig){
	  
      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodGenJets->At(ig));
      if(!tmp) continue;
	  
      if( tmp->Pt() < fJetPtCut ) continue;
      if( type == kJetsGenAcceptance &&
	  (    tmp->Eta() < fJetEtaMin
	    || tmp->Eta() > fJetEtaMax
	    || tmp->Phi() < fJetPhiMin
	    || tmp->Phi() > fJetPhiMax )) continue;
      
      list->Add(tmp);
      
      nGenJets++;
    }
    list->Sort();
    return nGenJets;
  } 
  else{
    if(fDebug>0)Printf("%s:%d no such type %d",(char*)__FILE__,__LINE__,type);
    return 0;
  }
}


// __________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y)
{
  //Set properties of histos (x and y title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
}

// _________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(TH2* h,const char* x, const char* y, const char* z)
{
  //Set properties of histos (x,y and z title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->SetZTitle(z);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->GetZaxis()->SetTitleColor(1);
}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetJetTracksPointing(TList* inputlist, TList* outputlist, AliAODJet* jet, const Double_t radius,Double_t& sumPt)
{
  // fill list of tracks in cone around jet axis  

  sumPt = 0;

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));

    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){

      outputlist->Add(track);
      
      sumPt += track->Pt();
    }
  }
  
  outputlist->Sort();
}

// ___________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetJetTracksTrackrefs(TList* list, AliAODJet* jet)
{
  // list of jet tracks from trackrefs
  
  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  for (Int_t itrack=0; itrack<nTracks; itrack++) {
    
    AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(itrack));
    if(!track){
      AliError("expected ref track not found ");
      continue;
    }
        
    list->Add(track);
  }
  
  list->Sort();
}
