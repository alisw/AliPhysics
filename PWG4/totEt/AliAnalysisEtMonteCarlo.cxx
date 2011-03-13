//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
	       //_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "TH2F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TList.h"
#include "AliESDCaloCluster.h"
#include "AliLog.h"

  using namespace std;

ClassImp(AliAnalysisEtMonteCarlo);


// ctor
AliAnalysisEtMonteCarlo::AliAnalysisEtMonteCarlo():AliAnalysisEt()
						  ,fImpactParameter(0)
						  ,fNcoll(0)
						  ,fNpart(0)
						  ,fHistPrimElectronEtaEET(0)  
						  ,fHistSecElectronEtaEET(0) 
						  ,fHistConvElectronEtaEET(0)  
						  ,fHistPrimGammaEtaEET(0)  
						  ,fHistPion0GammaEtaEET(0)  
						  ,fHistEtaGammaEtaEET(0)  
						  ,fHistOmega0GammaEtaEET(0)  
						  ,fHistSecGammaEtaEET(0)  

						  ,fHistPrimElectronEtaE(0) 
						  ,fHistSecElectronEtaE(0)
						  ,fHistConvElectronEtaE(0)
						  ,fHistPrimGammaEtaE(0) 
						  ,fHistPion0GammaEtaE(0)
						  ,fHistEtaGammaEtaE(0)
						  ,fHistOmega0GammaEtaE(0)
						  ,fHistSecGammaEtaE(0) 

						  ,fHistPrimElectronEtaERec(0) 
						  ,fHistSecElectronEtaERec(0)
						  ,fHistConvElectronEtaERec(0)
						  ,fHistPrimGammaEtaERec(0) 
						  ,fHistSecGammaEtaERec(0) 
						  ,fHistPion0GammaEtaERec(0)
						  ,fHistEtaGammaEtaERec(0)
						  ,fHistOmega0GammaEtaERec(0)

						  ,fHistAllERecEMC(0)	 
						  ,fHistGammaERecEMC(0)	
						  ,fHistElectronERecEMC(0)	

						  ,fHistElectronFirstMother(0)  
						  ,fHistElectronLastMother(0)  
						  ,fHistElectronFirstMotherEtaAcc(0)  
						  ,fHistElectronFirstMotherNPP(0)  
						  ,fHistElectronFirstMotherNPPAcc(0)  

						  ,fHistGammaFirstMother(0)  
						  ,fHistGammaLastMother(0)
						  ,fHistGammaFirstMotherEtaAcc(0)
						  ,fHistGammaFirstMotherNPP(0)
						  ,fHistGammaFirstMotherNPPAcc(0)
{
}

// dtor
AliAnalysisEtMonteCarlo::~AliAnalysisEtMonteCarlo() 
{
}

Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
  ResetEventValues();
	
  // Get us an mc event
  if(!ev){  
    AliFatal("ERROR: Event does not exist");
    return 0;
  }
  AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);
  if(!event){  
    AliFatal("ERROR: MC Event does not exist");
    return 0;
  }
	
  Double_t protonMass =fgProtonMass;
	
  // Hijing header
  AliGenEventHeader* genHeader = event->GenEventHeader();
  if(!genHeader){
    Printf("ERROR: Event generation header does not exist");   
    return 0;
  }
  AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  if (hijingGenHeader) {
    fImpactParameter = hijingGenHeader->ImpactParameter();
    fNcoll = hijingGenHeader->HardScatters(); // or should this be some combination of NN() NNw() NwN() NwNw() ?
    fNpart = hijingGenHeader->ProjectileParticipants() + hijingGenHeader->TargetParticipants(); 
    /*
      printf("Hijing: ImpactParameter %g ReactionPlaneAngle %g \n",
      hijingGenHeader->ImpactParameter(), hijingGenHeader->ReactionPlaneAngle());
      printf("HardScatters %d ProjecileParticipants %d TargetParticipants %d\n",
      hijingGenHeader->HardScatters(), hijingGenHeader->ProjectileParticipants(), hijingGenHeader->TargetParticipants()); 
      printf("ProjSpectatorsn %d ProjSpectatorsp %d TargSpectatorsn %d TargSpectatorsp %d\n",
      hijingGenHeader->ProjSpectatorsn(), hijingGenHeader->ProjSpectatorsp(), hijingGenHeader->TargSpectatorsn(), hijingGenHeader->TargSpectatorsp());
      printf("NN %d NNw %d NwN %d, NwNw %d\n",
      hijingGenHeader->NN(), hijingGenHeader->NNw(), hijingGenHeader->NwN(), hijingGenHeader->NwNw());
    */
  }
	
  /* // placeholder if we want to get some Pythia info later
     AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
     if (pythiaGenHeader) { // not Hijing; try with Pythia      
     printf("Pythia: ProcessType %d  GetPtHard %g \n",
     pythiaGenHeader->ProcessType(), pythiaGenHeader->GetPtHard());
     }
  */
	
  // Let's play with the stack!
  AliStack *stack = event->Stack();
	
  Int_t nPrim = stack->GetNtrack();
	
  for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {
		
      TParticle *part = stack->Particle(iPart);
      TParticle *partMom = 0;
      TParticle *partMomLast = 0;
		
      if (!part)
        {
	  Printf("ERROR: Could not get particle %d", iPart);
	  continue;
        }
		
      Int_t iPartMom = part->GetMother(0);
      Int_t iPartLastMom = part->GetMother(1);

      TParticlePDG *pdg = part->GetPDG(0);
      TParticlePDG *pdgMom = 0;
      TParticlePDG *pdgMomLast = 0;

      if (!pdg)
        {
	  //Printf("ERROR: Could not get particle PDG %d", iPart);
	  continue;
        }		
		
      if (iPartMom>0)
	{
	  partMom = stack->Particle(iPartMom);
	  pdgMom = partMom->GetPDG(0);
	}
			
      if (iPartLastMom>0)
	{
	  partMomLast = stack->Particle(iPartLastMom);
	  pdgMomLast = partMomLast->GetPDG(0);
	}
		
		
      Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle
		
      // Check if it is a primary particle
      //if (!stack->IsPhysicalPrimary(iPart)) continue;
		
      // it goes to the next particle in case it is not a physical primary or not a electron or gamma (want to check the contribution from secondary)
      if (!stack->IsPhysicalPrimary(iPart))
	{					
	  // this part is used only to check the origin of secondary electrons and gammas 
	  if (iPartMom>0)
	    {
	      //if (!stack->IsPhysicalPrimary(iPartMom)) continue;
				
	      if (pdgMom)
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    {
		      fHistElectronFirstMotherNPP->Fill(pdgMom->PdgCode());
		      // inside EMCal acceptance
		      if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
			fHistElectronFirstMotherNPPAcc->Fill(pdgMom->PdgCode());
		    }
		  else if (pdg->PdgCode() == fgGammaCode)
		    {		
		      fHistGammaFirstMotherNPP->Fill(pdgMom->PdgCode());
		      if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
			fHistGammaFirstMotherNPPAcc->Fill(pdgMom->PdgCode());
		    }
		}
	      else 
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    fHistElectronFirstMotherNPP->Fill(-598);
		  else if (pdg->PdgCode() == fgGammaCode)
		    fHistGammaFirstMotherNPP->Fill(-598);
					
		  continue;
		}
	    }
	  else
	    {
	      if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		fHistElectronFirstMotherNPP->Fill(-599);
	      else if (pdg->PdgCode() == fgGammaCode)
		fHistGammaFirstMotherNPP->Fill(-599);
				
	      continue;
	    }
			
	  if (!((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode) || (pdg->PdgCode() == fgGammaCode)))
	    continue;
	}
      // this part is used only to check the origin of physical primary electrons and gammas 
      else
	{					
	  if (iPartMom>0)
	    {
	      if (pdgMom)
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    fHistElectronFirstMother->Fill(pdgMom->PdgCode());
		  else if (pdg->PdgCode() == fgGammaCode)
		    fHistGammaFirstMother->Fill(pdgMom->PdgCode());
					
		  if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
		    {
		      if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
			fHistElectronFirstMotherEtaAcc->Fill(pdgMom->PdgCode());
		      else if (pdg->PdgCode() == fgGammaCode)
			fHistGammaFirstMotherEtaAcc->Fill(pdgMom->PdgCode());
		    }
		}
	      else 
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    fHistElectronFirstMother->Fill(-598);
		  else if (pdg->PdgCode() == fgGammaCode)
		    fHistGammaFirstMother->Fill(-598);
		}
	    }
	  else
	    {
	      if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		fHistElectronFirstMother->Fill(-599);
	      else if (pdg->PdgCode() == fgGammaCode)
		fHistGammaFirstMother->Fill(-599);
	    }
			
	  if (iPartLastMom>0)
	    {
	      if (pdgMomLast)
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    fHistElectronLastMother->Fill(pdgMomLast->PdgCode());
		  else if (pdg->PdgCode() == fgGammaCode)
		    fHistGammaLastMother->Fill(pdgMomLast->PdgCode());
		}
	      else
		{
		  if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		    fHistElectronLastMother->Fill(-598);
		  else if (pdg->PdgCode() == fgGammaCode)
		    fHistGammaLastMother->Fill(-598);
		}						
	    }
	  else
	    {
	      if ((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode))
		fHistElectronLastMother->Fill(-599);
	      else if (pdg->PdgCode() == fgGammaCode)
		fHistGammaLastMother->Fill(-599);
	    }
	}		

      //printf("MC: iPart %03d eta %4.3f phi %4.3f code %d charge %g \n", iPart, part->Eta(), part->Phi(), pdg->PdgCode(), pdg->Charge()); // tmp/debug printout
		
      // Check for reasonable (for now neutral and singly charged) charge on the particle
      //TODO:Maybe not only singly charged?
      if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;
		
      fMultiplicity++;
		
      if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
        {
	  // calculate E_T
	  if (
	      TMath::Abs(pdg->PdgCode()) == fgProtonCode ||
	      TMath::Abs(pdg->PdgCode()) == fgNeutronCode ||
	      TMath::Abs(pdg->PdgCode()) == fgLambdaCode ||
	      TMath::Abs(pdg->PdgCode()) == fgXiCode ||
	      TMath::Abs(pdg->PdgCode()) == fgXi0Code ||
	      TMath::Abs(pdg->PdgCode()) == fgOmegaCode
	      )
            {
	      if (pdg->PdgCode() > 0) { particleMassPart = - protonMass;}
	      if (pdg->PdgCode() < 0) { particleMassPart = protonMass;}
	    }
	  Double_t et = part->Energy() * TMath::Sin(part->Theta()) + particleMassPart;
			
	  // Fill up total E_T counters for each particle species
	  if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
	    {
	      fProtonEt += et;
	    }
	  if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
	    {
	      fPionEt += et;
	    }
	  if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
	    {
	      fChargedKaonEt += et;
	    }
	  if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
	    {
	      fMuonEt += et;
	    }
	  if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	    {
	      fElectronEt += et;
	    }
			
	  // some neutrals also
	  if(pdg->PdgCode() == fgNeutronCode)
	    {
	      fNeutronEt += et;
	    }
	  if(pdg->PdgCode() == fgAntiNeutronCode)
	    {
	      fAntiNeutronEt += et;
	    }
	  if(pdg->PdgCode() == fgGammaCode)
	    {
	      fGammaEt += et;
	    }
			
	  // Neutral particles
	  if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )
            {
	      fNeutralMultiplicity++;
	      fTotNeutralEt += et;
				
	      // inside EMCal acceptance
	      if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
		  fTotNeutralEtAcc += et;
		  fTotEtAcc += et;				

		  if(pdg->PdgCode() == fgGammaCode)
		    {
		      if (!stack->IsPhysicalPrimary(iPart))
			{
			  fHistSecGammaEtaE->Fill(part->Energy(),part->Eta());
			  fHistSecGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			}
		      else 
			{
			  if (pdgMom)
			    {
			      if (pdgMom->PdgCode() == fgPi0Code) 
				{
				  fHistPion0GammaEtaE->Fill(part->Energy(),part->Eta());
				  fHistPion0GammaEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			      else if (pdgMom->PdgCode() == fgOmega0Code) 
				{
				  fHistOmega0GammaEtaE->Fill(part->Energy(),part->Eta());
				  fHistOmega0GammaEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			      else if (pdgMom->PdgCode() == fgEtaCode)
				{
				  fHistEtaGammaEtaE->Fill(part->Energy(),part->Eta());
				  fHistEtaGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			      else
				{
				  fHistPrimGammaEtaE->Fill(part->Energy(),part->Eta());
				  fHistPrimGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			    }
			  else
			    {
			      fHistPrimGammaEtaE->Fill(part->Energy(),part->Eta());
			      fHistPrimGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			    }
			}
		    }
		}
            }
	  //Charged particles
	  else if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle())>1e-3 )
            {
	      fChargedMultiplicity++;
	      fTotChargedEt += et;
				
	      // inside EMCal acceptance
	      if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
		  fTotChargedEtAcc += et;
		  fTotEtAcc += et;
					
		  if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
		    {
		      fProtonEtAcc += et;
		    }
		  if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
		    {
		      fPionEtAcc += et;
		    }
		  if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
		    {
		      fChargedKaonEtAcc += et;
		    }
		  if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
		    {
		      fMuonEtAcc += et;
		    }
		  /*
		    if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
		    {
		    fElectronEtAcc += et;
		    }
		  */
		  if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
		    {
		      if (!stack->IsPhysicalPrimary(iPart))
			{
			  if (pdgMom)
			    {
			      if ((pdgMom->PdgCode() == fgGammaCode) && (stack->IsPhysicalPrimary(iPartMom)))
				{
				  fHistConvElectronEtaE->Fill(part->Energy(),part->Eta());
				  fHistConvElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			      else
				{
				  fHistSecElectronEtaE->Fill(part->Energy(),part->Eta());
				  fHistSecElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
				}
			    }
			  else
			    {
			      fHistSecElectronEtaE->Fill(part->Energy(),part->Eta());
			      fHistSecElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
			    }							
			}
		      else 
			{
			  fElectronEtAcc += et;
			  fHistPrimElectronEtaE->Fill(part->Energy(),part->Eta());
			  fHistPrimElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
			}
		    } // electron
		} // inside EMCal acceptance
				
	      //	  if (TrackHitsCalorimeter(part, event->GetMagneticField()))
	      if (TrackHitsCalorimeter(part)) // magnetic field info not filled?
		{
		  if (pdg->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
		  else fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
		}
	    }
	}
    }
    
  fTotEt = fTotChargedEt + fTotNeutralEt;
  fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;
    
  FillHistograms();
	
  return 0;    
}
//Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliMCEvent* mcEvent,AliESDEvent* realEvent)
Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{ // analyse MC and real event info
  //if(!mcEvent || !realEvent){
  if(!ev || !ev2){
    AliFatal("ERROR: Event does not exist");   
    return 0;
  }
  //AnalyseEvent(mcEvent);
  AnalyseEvent(ev);
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
  if(!mcEvent || !realEvent){  
    AliFatal("ERROR: mcEvent or realEvent does not exist");
    return 0;
  }

  AliStack *stack = mcEvent->Stack();
	
  // get all emcal clusters
  TRefArray* caloClusters = new TRefArray();
  realEvent->GetEMCALClusters( caloClusters );
	
  Int_t nCluster = caloClusters->GetEntries();
	
  // loop the clusters
  for (int iCluster = 0; iCluster < nCluster; iCluster++ ) 
    {
      AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
      Float_t caloE = caloCluster->E();

      UInt_t iPart = (UInt_t)TMath::Abs(caloCluster->GetLabel());
      TParticle *part  = stack->Particle(iPart);
      TParticle *partMom = 0;
		
      if (!part)
        {
	  Printf("No MC particle %d", iCluster);
	  continue;
        }
		
      // compare MC and Rec energies for all particles
      fHistAllERecEMC->Fill(part->Energy(),caloE);
		
      TParticlePDG *pdg = part->GetPDG(0);
      TParticlePDG *pdgMom = 0;
				
      if (!pdg)
        {
	  Printf("ERROR: Could not get particle PDG %d", iPart);
	  continue;
        }		

      Int_t iPartMom = part->GetMother(0);
		
      if (iPartMom>0)
	{
	  partMom = stack->Particle(iPartMom);
	  pdgMom = partMom->GetPDG(0);
	}			
		
      // Check if it is a primary particle
      //if (!stack->IsPhysicalPrimary(iPart)) continue;
      if (!stack->IsPhysicalPrimary(iPart)) // check whether particle is primary. we keep secondary electron and gamma for testing.
	{		
	  if (!((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode) || (pdg->PdgCode() == fgGammaCode)))
	    continue;
	} // end of primary particle check 
		
      if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;
		
      if(pdg->PdgCode() == fgGammaCode)
	{
	  // compare MC and Rec energies for gammas
	  fHistGammaERecEMC->Fill(part->Energy(),caloE);
			
	  if (!stack->IsPhysicalPrimary(iPart))
	    {
	      fHistSecGammaEtaERec->Fill(part->Energy(),part->Eta());
	    }
	  else
	    {
	      if (pdgMom)
		{
		  if (pdgMom->PdgCode() == fgPi0Code) 
		    {
		      fHistPion0GammaEtaERec->Fill(part->Energy(),part->Eta());
		    }
		  else if (partMom->GetPDG(0)->PdgCode() == fgOmega0Code) 
		    {
		      fHistOmega0GammaEtaERec->Fill(part->Energy(),part->Eta());
		    }
		  else if (partMom->GetPDG(0)->PdgCode() == fgEtaCode)
		    {
		      fHistEtaGammaEtaERec->Fill(part->Energy(),part->Eta());
		    }
		  else
		    {
		      fHistPrimGammaEtaERec->Fill(part->Energy(),part->Eta());
		    }
		}
	      else
		{
		  fHistPrimGammaEtaERec->Fill(part->Energy(),part->Eta());
		}
	    }
	} // gamma

      if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	{
	  // compare MC and Rec energies for electrons
	  fHistElectronERecEMC->Fill(part->Energy(),caloE);

	  if (!stack->IsPhysicalPrimary(iPart))
	    {
	      if (pdgMom)
		{
		  if ((pdgMom->PdgCode() == fgGammaCode) && (stack->IsPhysicalPrimary(iPartMom)))
		    {
		      fHistConvElectronEtaERec->Fill(part->Energy(),part->Eta());
		    }
		  else
		    {
		      fHistSecElectronEtaERec->Fill(part->Energy(),part->Eta());					
		    }					
		}
	      else
		{
		  fHistSecElectronEtaERec->Fill(part->Energy(),part->Eta());					
		}
	    }
	  else 
	    {
	      fHistPrimElectronEtaERec->Fill(part->Energy(),part->Eta());
	    }
	}
    } // end of loop over clusters	
  return 0;
}	

void AliAnalysisEtMonteCarlo::Init()
{ // init
  AliAnalysisEt::Init();
}

void AliAnalysisEtMonteCarlo::ResetEventValues()
{ // reset event values
  AliAnalysisEt::ResetEventValues();
	
  // collision geometry defaults for p+p:
  fImpactParameter = 0;
  fNcoll = 1;
  fNpart = 2;  
}

void AliAnalysisEtMonteCarlo::CreateHistograms()
{ // histogram related additions
  AliAnalysisEt::CreateHistograms();
  if (fTree) {
    fTree->Branch("fImpactParameter",&fImpactParameter,"fImpactParameter/D");
    fTree->Branch("fNcoll",&fNcoll,"fNcoll/I");
    fTree->Branch("fNpart",&fNpart,"fNpart/I");
  }
  fHistPrimElectronEtaEET = CreateEtaEHisto2D("fHistPrimElectronEtaEET","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistSecElectronEtaEET = CreateEtaEHisto2D("fHistSecElectronEtaEET","MC E_{T}, secondary (no conversion) electrons","E_{T}(GeV)");
  fHistConvElectronEtaEET = CreateEtaEHisto2D("fHistConvElectronEtaEET","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistPrimGammaEtaEET = CreateEtaEHisto2D("fHistPrimGammaEtaEET","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistPion0GammaEtaEET = CreateEtaEHisto2D("fHistPion0GammaEtaEET","MC E_{T}, #pi^{0}","E_{T}(GeV)");
  fHistOmega0GammaEtaEET = CreateEtaEHisto2D("fHistOmega0GammaEtaEET","MC E_{T}, #omega^{0}","E_{T}(GeV)");
  fHistEtaGammaEtaEET = CreateEtaEHisto2D("fHistEtaGammaEtaEET","MC E_{T}, #eta","E_{T}(GeV)");
  fHistSecGammaEtaEET = CreateEtaEHisto2D("fHistSecGammaEtaEET","MC E_{T}, secondary (no #pi^{0}, #eta or #omega) gammas","E_{T}(GeV)"); 

  fHistPrimElectronEtaE = CreateEtaEHisto2D("fHistPrimElectronEtaE","MC E_{T}, primary electrons","#");
  fHistSecElectronEtaE = CreateEtaEHisto2D("fHistSecElectronEtaE","MC E_{T}, secondary (no conversion) electrons","#");
  fHistConvElectronEtaE = CreateEtaEHisto2D("fHistConvElectronEtaE","MC E_{T}, electrons from conversion","#");
  fHistPrimGammaEtaE = CreateEtaEHisto2D("fHistPrimGammaEtaE","MC E_{T}, primary gammas","#"); 
  fHistPion0GammaEtaE = CreateEtaEHisto2D("fHistPion0GammaEtaE","MC E_{T}, #pi^{0}","#");
  fHistOmega0GammaEtaE = CreateEtaEHisto2D("fHistOmega0GammaEtaE","MC E_{T}, #omega^{0}","#");
  fHistEtaGammaEtaE = CreateEtaEHisto2D("fHistEtaGammaEtaE","MC E_{T}, #eta","#");
  fHistSecGammaEtaE = CreateEtaEHisto2D("fHistSecGammaEtaE","MC E_{T}, secondary (no #pi^{0}, #eta or #omega) gammas","#"); 

  fHistPrimElectronEtaERec = CreateEtaEHisto2D("fHistPrimElectronEtaERec","MC E_{T}, primary electrons","#");
  fHistSecElectronEtaERec = CreateEtaEHisto2D("fHistSecElectronEtaERec","MC E_{T}, secondary (no conversion) electrons","#");
  fHistConvElectronEtaERec = CreateEtaEHisto2D("fHistConvElectronEtaERec","MC E_{T}, electrons from conversion","#");
  fHistPrimGammaEtaERec = CreateEtaEHisto2D("fHistPrimGammaEtaERec","MC E_{T}, primary gammas","#"); 
  fHistPion0GammaEtaERec = CreateEtaEHisto2D("fHistPion0GammaEtaERec","MC E_{T}, #pi^{0}","#");
  fHistOmega0GammaEtaERec = CreateEtaEHisto2D("fHistOmega0GammaEtaERec","MC E_{T}, #omega","#");
  fHistEtaGammaEtaERec = CreateEtaEHisto2D("fHistEtaGammaEtaERec","MC E_{T}, #eta","#");
  fHistSecGammaEtaERec = CreateEtaEHisto2D("fHistSecGammaEtaERec","MC E_{T}, secondary (no #pi^{0}, #eta or #omega) gammas","#"); 

  fHistAllERecEMC = new TH2F("fHistAllERecEMC","E cluster Rec vs MC, all particles",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistAllERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistAllERecEMC->SetYTitle("E_{Rec}(GeV)");

  fHistElectronERecEMC = new TH2F("fHistElectronERecEMC","E cluster Rec vs MC, Electrons",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistElectronERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistElectronERecEMC->SetYTitle("E_{Rec}(GeV)");
	
  fHistGammaERecEMC = new TH2F("fHistGammaERecEMC","E cluster Rec vs MC, Gammas",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistGammaERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistGammaERecEMC->SetYTitle("E_{Rec}(GeV)");
	
  fHistElectronFirstMother = new TH1F("fHistElectronFirstMother","Electron First Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistElectronLastMother = new TH1F("fHistElectronLastMother","Electron Last Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistElectronFirstMotherEtaAcc = new TH1F("fHistElectronFirstMotherEtaAcc","Electron First Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistElectronFirstMotherNPP = new TH1F("fHistElectronFirstMotherNPP","Electron First Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistElectronFirstMotherNPPAcc = new TH1F("fHistElectronFirstMotherNPPAcc","Electron First Mother PDG Code Distribution",1201,-600.5,600.5);

  fHistGammaFirstMother = new TH1F("fHistGammaFirstMother","Gamma First Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistGammaLastMother = new TH1F("fHistGammaLastMother","Gamma Last Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistGammaFirstMotherEtaAcc = new TH1F("fHistGammaFirstMotherEtaAcc","Gamma First Mother PDG Code Distribution",1201,-600.5,600.5);	
  fHistGammaFirstMotherNPP = new TH1F("fHistGammaFirstMotherNPP","Gamma First Mother PDG Code Distribution",1201,-600.5,600.5);
  fHistGammaFirstMotherNPPAcc = new TH1F("fHistGammaFirstMotherNPP","Gamma First Mother PDG Code Distribution",1201,-600.5,600.5);	
}

void AliAnalysisEtMonteCarlo::FillOutputList(TList *list)
{//fill the output list
  AliAnalysisEt::FillOutputList(list);

  list->Add(fHistPrimElectronEtaEET);
  list->Add(fHistSecElectronEtaEET);
  list->Add(fHistConvElectronEtaEET);
  list->Add(fHistPrimGammaEtaEET); 
  list->Add(fHistPion0GammaEtaEET);
  list->Add(fHistOmega0GammaEtaEET);
  list->Add(fHistEtaGammaEtaEET);
  list->Add(fHistSecGammaEtaEET); 
	
  list->Add(fHistPrimElectronEtaE);
  list->Add(fHistSecElectronEtaE);
  list->Add(fHistConvElectronEtaE);
  list->Add(fHistPrimGammaEtaE); 
  list->Add(fHistPion0GammaEtaE);
  list->Add(fHistOmega0GammaEtaE);
  list->Add(fHistEtaGammaEtaE);
  list->Add(fHistSecGammaEtaE); 
	
  list->Add(fHistPrimElectronEtaERec);
  list->Add(fHistSecElectronEtaERec);
  list->Add(fHistConvElectronEtaERec);
  list->Add(fHistPrimGammaEtaERec); 
  list->Add(fHistPion0GammaEtaERec);
  list->Add(fHistOmega0GammaEtaERec);
  list->Add(fHistEtaGammaEtaERec);
  list->Add(fHistSecGammaEtaERec); 
	
  list->Add(fHistAllERecEMC);
  list->Add(fHistElectronERecEMC);
  list->Add(fHistGammaERecEMC);
	
  list->Add(fHistElectronFirstMother);
  list->Add(fHistElectronLastMother);
  list->Add(fHistElectronFirstMotherEtaAcc);
  list->Add(fHistElectronFirstMotherNPP);
  list->Add(fHistElectronFirstMotherNPPAcc);

  list->Add(fHistGammaFirstMother);
  list->Add(fHistGammaLastMother);
  list->Add(fHistGammaFirstMotherEtaAcc);
  list->Add(fHistGammaFirstMotherNPP);
  list->Add(fHistGammaFirstMotherNPPAcc);
}


bool AliAnalysisEtMonteCarlo::TrackHitsCalorimeter(TParticle* part, Double_t magField)
{
  //  printf(" TrackHitsCalorimeter - magField %f\n", magField);
  AliESDtrack *esdTrack = new AliESDtrack(part);
  // Printf("MC Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());
	
  Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);
	
  // if(prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());
	
  bool status = prop && 
    TMath::Abs(esdTrack->Eta()) < fEtaCutAcc && 
    esdTrack->Phi() > fPhiCutAccMin*TMath::Pi()/180. && 
    esdTrack->Phi() < fPhiCutAccMax*TMath::Pi()/180.;
  delete esdTrack;
	
  return status;
}

