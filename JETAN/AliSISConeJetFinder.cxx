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
 

//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
//
// Author: swensy.jangal@ires.in2p3.fr 
//
// Last modification: Neutral cell energy included in the jet reconstruction
//
// Author: Magali.estienne@subatech.in2p3.fr 
//---------------------------------------------------------------------


#include <Riostream.h>
#include <TArrayF.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TClonesArray.h>

#include "AliHeader.h"
#include "AliJetKineReader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliJetUnitArray.h"
#include "AliSISConeJetFinder.h"
#include "AliSISConeJetHeader.h"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

// get info on how fastjet was configured
#include "fastjet/config.h"

#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<sstream>  // needed for internal io
#include<vector> 
#include <cmath> 

using namespace std;


ClassImp(AliSISConeJetFinder)


//____________________________________________________________________________

AliSISConeJetFinder::AliSISConeJetFinder():
AliJetFinder()
{
  // Constructor
}

//____________________________________________________________________________

AliSISConeJetFinder::~AliSISConeJetFinder()
{
  // destructor
}

//______________________________________________________________________________
void AliSISConeJetFinder::FindJets()
{

  // Pick up siscone header
  AliSISConeJetHeader *header      = (AliSISConeJetHeader*)fHeader;
  AliJetReaderHeader *readerheader = (AliJetReaderHeader*)fReader->GetReaderHeader();

  Int_t debug = header->GetDebug();                        // debug option
  Int_t fOpt  = fReader->GetReaderHeader()->GetDetector();


//******************************** SISCONE PLUGIN CONFIGURATION
// Here we look for SISCone parameters in the header and we define our plugin.  

  Double_t coneRadius                = header->GetConeRadius();       // cone radius
  Double_t overlapThreshold          = header->GetOverlapThreshold(); // overlap parameter
  Int_t    nPassMax                  = header->GetNPassMax();         // maximum number of passes
  Double_t ptProtoJetMin             = header->GetPtProtojetMin();    // pT min of protojets
  Double_t caching                   = header->GetCaching();          // do we record found cones for this set of data?

//  if (header->GetSplitMergeScale() == 0) fastjet::SISConePlugin::SplitMergeScale splitMergeScale = fastjet::SISConePlugin::SM_pttilde; // There's only one split merge scale
//  Double_t splitMergeStoppingScale = header->GetSplitMergeStoppingScale(); // Additional cut on pt_tilde of protojets

  fastjet::JetDefinition::Plugin * plugin;
  plugin = new fastjet::SISConePlugin(coneRadius, overlapThreshold, nPassMax, ptProtoJetMin, caching);

  vector<fastjet::PseudoJet> inputParticles; 

//********************************  READING OF INPUT PARTICLES
// Here we look for px, py pz and energy of each particle that we gather in a PseudoJet object, and we put all these PseudoJet in a vector of PseudoJets : input_particles. 

  Double_t pttrackcut = readerheader->GetPtCut(); 
  if (debug)cout<<"pT track cut for SISCone jet finder = "<<pttrackcut<<" GeV"<<endl;

  if(fOpt==0) 
  {
    TClonesArray *lvArray    = fReader->GetMomentumArray();

    // We check if lvArray is ok
    if(lvArray == 0)
    {
      cout << "Could not get the momentum array" << endl;
      delete plugin;
      return;
    }
	
    Int_t nIn = lvArray->GetEntries();
	
    if(nIn == 0)// nIn = Number of particles in the event
    {
      if (debug) cout << "entries = 0 ; Event empty !!!" << endl ;
      delete plugin;
      return;
    }
	
    Float_t px,py,pz,en;
	
    // Load input vectors
    for(Int_t i = 0; i < nIn; i++)
    { 
      TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
      px = lv->Px();
      py = lv->Py();
      pz = lv->Pz();
      en = lv->Energy();

      if (pttrackcut > TMath::Sqrt(px*px+py*py)) continue; //added by syssy
	    
      fastjet::PseudoJet inputPart(px,py,pz,en); 
      inputPart.set_user_index(i);
      inputParticles.push_back(inputPart); 
    }
  }
  else 
  {
    TClonesArray* fUnit = fReader->GetUnitArray();
    if(fUnit == 0) { cout << "Could not get the momentum array" << endl; 	    delete plugin; return; }
    Int_t         nIn = fUnit->GetEntries();
    if(nIn == 0) { if (debug) cout << "entries = 0 ; Event empty !!!" << endl ;  	    delete plugin; return; }
    // Information extracted from fUnitArray
    // load input vectors and calculate total energy in array
    Float_t pt,eta,phi,theta,px,py,pz,en;
    Int_t ipart = 0;
    for(Int_t i=0; i<nIn; i++) 
    {
      AliJetUnitArray *uArray = (AliJetUnitArray*)fUnit->At(i);
	  
      if(uArray->GetUnitEnergy()>0.)
      {
        // It is not necessary anymore to cut on particle pt
	pt    = uArray->GetUnitEnergy();
	eta   = uArray->GetUnitEta();
	phi   = uArray->GetUnitPhi();
	theta = EtaToTheta(eta);
	en    = (TMath::Abs(TMath::Sin(theta)) == 0) ? pt : pt/TMath::Abs(TMath::Sin(theta));
	px    = TMath::Cos(phi)*pt;
	py    = TMath::Sin(phi)*pt;
	pz    = en*TMath::TanH(eta);

        if (pttrackcut > TMath::Sqrt(px*px+py*py)) continue; 

	if(debug) cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << ", en: " << en << ", px: " << px << ", py: " << py << ", pz: " << pz << endl;
	    
	fastjet::PseudoJet inputPart(px,py,pz,en); // create PseudoJet object
	inputPart.set_user_index(ipart); //label the particle into Fastjet algortihm
	inputParticles.push_back(inputPart);  // back of the input_particles vector 
	ipart++;
      }
    } // End loop on UnitArray 
  }
  
//******************************** CHOICE OF JET AREA  
// Here we determine jets area for subtracting background later
// For more informations about jet areas see : The Catchment Area of Jets M. Cacciari, G. Salam and G. Soyez

  Double_t ghostEtamax        = header->GetGhostEtaMax();       // maximum eta in which a ghost can be generated
  Double_t ghostArea          = header->GetGhostArea();         // area of a ghost
  Int_t    activeAreaRepeats  = header->GetActiveAreaRepeats(); // do we repeat area calculation?
  Double_t gridScatter        = header->GetGridScatter();       // fractional random fluctuations of the position of the ghosts on the y-phi grid
  Double_t ktScatter          = header->GetKtScatter();         // fractional random fluctuations of the tranverse momentum of the ghosts on the y-phi grid
  Double_t meanGhostKt        = header->GetMeanGhostKt();       // average transverse momentum of the ghosts.

  Double_t areaTypeNumber  = header->GetAreaTypeNumber();       // the number determines jet area type 
  fastjet::AreaType areaType = fastjet::active_area;
  if (areaTypeNumber == 1)  areaType = fastjet::active_area;                 
  if (areaTypeNumber == 2)  areaType = fastjet::active_area_explicit_ghosts; 
  if (areaTypeNumber == 3)  areaType = fastjet::one_ghost_passive_area;      
  if (areaTypeNumber == 4)  areaType = fastjet::passive_area;                
  if (areaTypeNumber == 5)  areaType = fastjet::voronoi_area;                

  fastjet::AreaDefinition areaDef;

  if (areaTypeNumber < 5) 
  {
    fastjet::GhostedAreaSpec ghostSpec(ghostEtamax, activeAreaRepeats, ghostArea, gridScatter, ktScatter, meanGhostKt);
    areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
  }

  if (areaTypeNumber == 5)
  {
    Double_t effectiveRFact = header->GetEffectiveRFact();
    fastjet::VoronoiAreaSpec ghostSpec(effectiveRFact);
    areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
  }

//********************************
//********************************

  Bool_t bgMode = header->GetBGMode();// Here one choose to subtract BG or not

//********************************
//********************************

//********************************
//******************************** BG SUBTRACTION
  if (bgMode == 1)// BG subtraction
    {
    //******************************** JETS FINDING AND EXTRACTION
    fastjet::ClusterSequenceArea clustSeq(inputParticles, plugin, areaDef);

    // Here we extract inclusive jets with pt > ptmin, sorted by pt 
    Double_t ptMin = header->GetMinJetPt(); 
    vector<fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets(ptMin);
    vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets);
  
    //***************************** BACKGROUND SUBTRACTION

    // Set the rapidity-azimuth range within which to study background
    Double_t rapMin = header->GetRapMin();
    Double_t rapMax = header->GetRapMax();
    Double_t phiMin = header->GetPhiMin();
    Double_t phiMax = header->GetPhiMax();
    fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);

    // As SISCone gives too small areas to estimate background,
    // we have to choose between kT and Cambridge/Aachen algorithms
    // to estimate mean background rho
    fastjet::JetAlgorithm algorithm = fastjet::kt_algorithm;
    Int_t algo = header->GetBGAlgorithm();
    if (algo == 0) algorithm = fastjet::kt_algorithm;
    if (algo == 1) algorithm = fastjet::cambridge_algorithm;

    // Rho estimation :
    // Radius used to calculate rho, can be different from the one used to reconstruct jets
    Double_t RRho = header->GetRForRho(); 
    fastjet::JetDefinition jetDefinitionForRhoEstimation(algorithm, RRho);
    fastjet::ClusterSequenceArea csForRho(inputParticles, jetDefinitionForRhoEstimation, areaDef);
    vector<fastjet::PseudoJet> inclusiveJetsForRhoEstimation = sorted_by_pt(csForRho.inclusive_jets());

    // Number of hard jets not to be used to estimate rho
    Int_t NHardJets = header->GetNumberOfJetsToErase();
    if (debug) cout<<"Number of jets before subtraction : "<<inclusiveJetsForRhoEstimation.size()<<endl;
    if (debug) cout<<"Number of jets not to count into rho estimation : "<<NHardJets<<endl;
    Bool_t rhoOnJet = 0;
    if (inclusiveJetsForRhoEstimation.size() == 1 && inclusiveJetsForRhoEstimation[0].perp() > 5) rhoOnJet = 1;//if there's only one jet (not bg kind), rho = 0
    if (inclusiveJetsForRhoEstimation.size() <= NHardJets && inclusiveJetsForRhoEstimation.size() != 1  && inclusiveJetsForRhoEstimation.size() != 0)
    {
      //      inclusiveJetsForRhoEstimation.erase(inclusiveJetsForRhoEstimation.begin(),inclusiveJetsForRhoEstimation.begin()+1);
      inclusiveJetsForRhoEstimation.erase(inclusiveJetsForRhoEstimation.begin());
      cout<<"Only the 1st (hardest) jet of the event hasn't been taken into account for rho estimation"<<endl;
    }
    if (inclusiveJetsForRhoEstimation.size() > NHardJets) inclusiveJetsForRhoEstimation.erase(inclusiveJetsForRhoEstimation.begin(), inclusiveJetsForRhoEstimation.begin()+NHardJets);
    
    // Estimation of rho and fluctuations sigma

    //1st method
    Double_t rho      = 0;
    Double_t sigma    = 0;
    Double_t meanarea = 0;
    // 3rd argument   : 1 = use area 4 vector rather than simple area
    // last argument  : 0 = in case of explicit ghosts use
    if (inclusiveJetsForRhoEstimation.size() != 0) csForRho.get_median_rho_and_sigma(inclusiveJetsForRhoEstimation, range, 1, rho, sigma, meanarea, 1); // this gives the fluctuations also

    if (rhoOnJet) rho = 0;

    if(debug) cout<<"rho = "<<rho<<endl;

    // Vector of corrected jets
    vector<fastjet::PseudoJet> corrJets = sorted_by_pt(clustSeq.subtracted_jets(rho, ptMin));
  
    //***************************** JETS DISPLAY

    for (size_t j = 0; j < jets.size(); j++)
    {
//       // If the jet is only made of ghosts, continue.
//       if (clustSeq.is_pure_ghost(corrJets[j]) == 1) continue;

      // If the correction is > jet energy px = py = pz = e = 0
      if (corrJets[j].px() == 0 && corrJets[j].py() == 0 && corrJets[j].pz() == 0 && corrJets[j].E() == 0) continue;

      if(debug)
      {
	cout<<"********************************** Reconstructed jet(s) (not corrected)"<<endl;
	cout<<"Jet number "<<j+1<<" : "<<"Rapidity : "<<jets[j].rap()<<" Phi : "<<jets[j].phi()<<" pT : "<<jets[j].perp()<<" area : "<<clustSeq.area(jets[j])<<endl;
	cout<<"e = "<<jets[j].E()<<endl;
	
	cout<<"********************************** Corrected jet(s)"<<endl;
	cout<<"Jet number "<<j+1<<" : "<<"Rapidity : "<<corrJets[j].rap()<<" Phi : "<<corrJets[j].phi()<<" pT : "<<corrJets[j].perp()<<endl;
	cout<<"e = "<<corrJets[j].E()<<endl;
      }
      // Go to write AOD info
      AliAODJet aodjet (corrJets[j].px(), corrJets[j].py(), corrJets[j].pz(), corrJets[j].E());
      if(debug) aodjet.Print("");
      AddJet(aodjet);
    }
  }

//********************************
//******************************** NO BG SUBTRACTION

  if (bgMode == 0)// No BG subtraction
  {
    //******************************** JETS FINDING AND EXTRACTION
    fastjet::ClusterSequenceArea clustSeq(inputParticles, plugin, areaDef);
    // Here we extract inclusive jets with pt > ptmin, sorted by pt 
    Double_t ptMin = header->GetMinJetPt(); 
    vector<fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets(ptMin);
    vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets);

    //***************************** JETS DISPLAY

    for (size_t k = 0; k < jets.size(); k++)
    {
      if(debug)
      {
	cout<<"********************************** Reconstructed jet(s) (not corrected)"<<endl;
	cout<<"Jet number "<<k+1<<" : "<<"Rapidity : "<<jets[k].rap()<<" Phi : "<<jets[k].phi()<<" pT : "<<jets[k].perp()<<endl;
	cout<<"e = "<<jets[k].E()<<endl;
      }

      // Go to write AOD info
      Double_t area = clustSeq.area(jets[k]);
      AliAODJet aodjet (jets[k].px(), jets[k].py(), jets[k].pz(), jets[k].E());
      aodjet.SetEffArea(area,0);
      if(debug) aodjet.Print("");
      AddJet(aodjet);
    }
  }

  delete plugin;

}
 
//____________________________________________________________________________

Float_t  AliSISConeJetFinder::EtaToTheta(Float_t arg)
{
  //  return (180./TMath::Pi())*2.*atan(exp(-arg));
  return 2.*atan(exp(-arg));


}

//____________________________________________________________________________

void AliSISConeJetFinder::InitTask(TChain *tree)
{

  printf("SISCone jet finder initialization ******************");
  fReader->CreateTasks(tree);

}


//____________________________________________________________________________

void AliSISConeJetFinder::WriteJHeaderToFile() const
{
  fHeader->Write();
}

//____________________________________________________________________________

