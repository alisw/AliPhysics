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
#include <TMath.h>
#include <TList.h>
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TVector3.h"

//---- ANALYSIS system ----
#include "AliMCAnalysisUtils.h"
#include "AliCaloTrackReader.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMCAnalysisUtils) ;
/// \endcond

//________________________________________
/// Constructor
//________________________________________
AliMCAnalysisUtils::AliMCAnalysisUtils() : 
TObject(), 
fCurrentEvent(-1), 
fDebug(0),
fJetsList(new TList), 
fMCGenerator(kPythia),
fMCGeneratorString("PYTHIA"),
fDaughMom(),  fDaughMom2(),
fMotherMom(), fGMotherMom()
{}

//_______________________________________
/// Destructor.
//_______________________________________
AliMCAnalysisUtils::~AliMCAnalysisUtils() 
{  
  if (fJetsList)
  {
    fJetsList->Clear();
    delete fJetsList ;
  }     
}

//_____________________________________________________________________________________________
/// Check the first common ancestor of 2 clusters, given the most likely labels 
/// of the primaries generating such clusters.
//_____________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckCommonAncestor(Int_t index1, Int_t index2, 
                                              const AliCaloTrackReader* reader, 
                                              Int_t & ancPDG, Int_t & ancStatus, 
                                              TLorentzVector & momentum, TVector3 & prodVertex) 
{  
  Int_t label1[100];
  Int_t label2[100];
  label1[0]= index1;
  label2[0]= index2;
  Int_t counter1 = 0;
  Int_t counter2 = 0;
  
  if(label1[0]==label2[0])
  {
    //printf("AliMCAnalysisUtils::CheckCommonAncestor() - Already the same label: %d\n",label1[0]);
    counter1=1;
    counter2=1;
  }
  else
  {
    if(reader->ReadAODMCParticles())
    {      
      Int_t label=label1[0];
      if(label < 0) return -1;
      
      while(label > -1 && counter1 < 99)
      {
        counter1++;
        AliAODMCParticle * mom = (AliAODMCParticle *) (reader->GetMC())->GetTrack(label);
        if(mom)
        {
          label  = mom->GetMother() ;
          label1[counter1]=label;
        }
        //printf("\t counter %d, label %d\n", counter1,label);
      }
     
      //printf("Org label2=%d,\n",label2[0]);
      label=label2[0];
      if(label < 0) return -1;
      
      while(label > -1 && counter2 < 99)
      {
        counter2++;
        AliAODMCParticle * mom = (AliAODMCParticle *) (reader->GetMC())->GetTrack(label);
        if(mom)
        {
          label  = mom->GetMother() ;
          label2[counter2]=label;
        }
        //printf("\t counter %d, label %d\n", counter2,label);
      }
    }//AOD MC
    else
    { //Kine stack from ESDs
      AliMCEvent * mcevent = reader->GetMC();
      
      Int_t label=label1[0];
      while(label > -1 && counter1 < 99)
      {
        counter1++;
        TParticle * mom = mcevent->Particle(label);
        if(mom)
        {
          label  = mom->GetFirstMother() ;
          label1[counter1]=label;
        }
        //printf("\t counter %d, label %d\n", counter1,label);
      }
      
      //printf("Org label2=%d,\n",label2[0]);
      
      label=label2[0];
      while(label > -1 && counter2 < 99)
      {
        counter2++;
        TParticle * mom = mcevent->Particle(label);
        if(mom)
        {
          label  = mom->GetFirstMother() ;
          label2[counter2]=label;
        }
        //printf("\t counter %d, label %d\n", counter2,label);
      }
    }// Kine stack from ESDs
  }//First labels not the same
  
//  if((counter1==99 || counter2==99) && fDebug >=0)
//    printf("AliMCAnalysisUtils::CheckCommonAncestor() - Genealogy too large c1: %d, c2= %d\n", counter1, counter2);
  //printf("CheckAncestor:\n");
  
  Int_t commonparents = 0;
  Int_t ancLabel = -1;
  //printf("counters %d %d \n",counter1, counter2);
  for (Int_t c1 = 0; c1 < counter1; c1++)
  {
    for (Int_t c2 = 0; c2 < counter2; c2++)
    {
      if(label1[c1]==label2[c2] && label1[c1]>-1)
      {
        ancLabel = label1[c1];
        commonparents++;
        
        if(reader->ReadAODMCParticles())
        {
          AliAODMCParticle * mom = (AliAODMCParticle *) (reader->GetMC())->GetTrack(label1[c1]);
          
          if (mom)
          {
            ancPDG    = mom->GetPdgCode();
            ancStatus = mom->GetStatus();
            momentum.SetPxPyPzE(mom->Px(),mom->Py(),mom->Pz(),mom->E());
            prodVertex.SetXYZ(mom->Xv(),mom->Yv(),mom->Zv());
          }
        }
        else
        {
          TParticle * mom = reader->GetMC()->Particle(label1[c1]);
          
          if (mom)
          {
            ancPDG    = mom->GetPdgCode();
            ancStatus = mom->GetStatusCode();
            mom->Momentum(momentum);
            prodVertex.SetXYZ(mom->Vx(),mom->Vy(),mom->Vz());
          }
        }
        
        //First ancestor found, end the loops
        counter1=0;
        counter2=0;
      }//Ancestor found
    }//second cluster loop
  }//first cluster loop
  
  if(ancLabel < 0)
  {
    ancPDG    = -10000;
    ancStatus = -10000;
    momentum.SetXYZT(0,0,0,0);
    prodVertex.SetXYZ(-10,-10,-10);
  }
  
  return ancLabel;
}

//________________________________________________________________________________________
/// \return tag with primary particle at the origin of the cluster/track.
/// Here we just check if the event is AOD or ESD and then assign to the corresponding method.
/// See CheckOrigin() CheckOriginInStack() and CheckOriginInAOD()
//________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOrigin(const Int_t * label, Int_t nlabels,
                                      const AliCaloTrackReader* reader, Int_t calorimeter)
{
  if( nlabels <= 0 )
  {
    AliWarning("No MC labels available, please check!!!");
    return kMCBadLabel;
  }
  
  if ( !reader->GetMC() )
  {
    AliDebug(1,"MCEvent is not available, check analysis settings in configuration file, do nothing with MC!!");
    return -1;
  }
  
  Int_t tag = 0;
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();

  // Bad label
  if ( label[0] < 0 || label[0] >= nprimaries )
  {
    if(label[0] < 0)
      AliWarning(Form("*** bad label ***:  label %d", label[0]));
    
    if(label[0] >=  nprimaries)
      AliWarning(Form("*** large label ***:  label %d, n tracks %d", label[0], nprimaries));
    
    SetTagBit(tag,kMCUnknown);
    
    return tag; 
  } // Bad label
  
  TObjArray* arrayCluster = 0;
  if      ( calorimeter == AliCaloTrackReader::kEMCAL ) arrayCluster = reader->GetEMCALClusters();
  else if ( calorimeter == AliCaloTrackReader::kPHOS  ) arrayCluster = reader->GetPHOSClusters ();
    
  // Select where the information is, ESD-galice stack or AOD mcparticles branch

  if(reader->ReadStack())
  {
    tag = CheckOriginInStack(label, nlabels, reader->GetMC(), arrayCluster);
  }
  else if(reader->ReadAODMCParticles())
  {
    tag = CheckOriginInAOD  (label, nlabels, reader->GetMC(), arrayCluster);
  }
  
  return tag ;
}

//____________________________________________________________________________________________________
/// \return tag with primary particle at the origin of the cluster/track.
/// Here we just check if the event is AOD or ESD and then assign to the corresponding method 
/// and we have only one input MC label not multiple. 
/// See CheckOrigin() CheckOriginInStack() and CheckOriginInAOD()
//_____________________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOrigin(Int_t label, const AliCaloTrackReader* reader, Int_t calorimeter)
{    
  if( label < 0 )
  {
    AliWarning("No MC labels available, please check!!!");
    return kMCBadLabel;
  }
  
  Int_t labels[] = { label };
  
  return CheckOrigin(labels, 1,reader, calorimeter);  
}	

//__________________________________________________________________________________________
/// \return tag with primary particle(S) at the origin of the cluster/track.
/// Do this for ESDs, same things as in CheckOriginInAOD.
///
/// Generally speaking, label is the MC label of a reconstructed
/// entity (track, cluster, etc) for which we want to know something 
/// about its heritage, but one can also use it directly with stack 
/// particles not connected to reconstructed entities.
//__________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInStack(const Int_t *labels, Int_t nlabels,
                                             const AliMCEvent* mcevent, const TObjArray* arrayCluster)
{
  Int_t tag = 0;
  
  // Most significant particle contributing to the cluster
  Int_t label=labels[0]; 
    
  // MC particle of interest is the "mom" of the entity
  TParticle * mom = mcevent->Particle(label);
  Int_t iMom     = label;
  Int_t mPdgSign = mom->GetPdgCode();
  Int_t mPdg     = TMath::Abs(mPdgSign);
  Int_t mStatus  = mom->GetStatusCode() ;
  Int_t iParent  = mom->GetFirstMother() ;
  
  //if( label < 8 && fMCGenerator != kBoxLike ) AliDebug(1,Form("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d",iParent));
  
  // GrandParent of the entity
  TParticle * parent = NULL;
  Int_t pPdg = -1;
  Int_t pStatus =-1;
  if(iParent >= 0)
  {
    parent = mcevent->Particle(iParent);
    
    if(parent)
    {
      pPdg = TMath::Abs(parent->GetPdgCode());
      pStatus = parent->GetStatusCode();  
    }
  }
  else AliDebug(1,Form("Parent with label %d",iParent));
  
  AliDebug(2,"Cluster most contributing mother and its parent:");
  AliDebug(2,Form("\t Mother label %d, pdg %d, status %d",iMom, mPdg, mStatus));
  AliDebug(2,Form("\t Parent label %d, pdg %d, status %d",iParent, pPdg, pStatus));
  
  // Check if "mother" of entity is converted, if not, get the first non converted mother
  if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && mStatus == 0)
  {
    SetTagBit(tag,kMCConversion);
    
    //Check if the mother is photon or electron with status not stable
    while ((pPdg == 22 || pPdg == 11) && mStatus != 1)
    {
      //Mother
      iMom     = mom->GetFirstMother();
      
      if(iMom < 0) 
      {
        AliInfo(Form("pdg = %d, mother = %d, skip",pPdg,iMom));
        break;
      }
      
      mom      = mcevent->Particle(iMom);
      mPdgSign = mom->GetPdgCode();
      mPdg     = TMath::Abs(mPdgSign);
      mStatus  = mom->GetStatusCode() ;
      iParent  = mom->GetFirstMother() ;
      
      //if(label < 8 ) AliDebug(1,Form("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent));
      
      //GrandParent
      if(iParent >= 0)
      {
        parent = mcevent->Particle(iParent);
        if(parent)
        {
          pPdg = TMath::Abs(parent->GetPdgCode());
          pStatus = parent->GetStatusCode();  
        }
      }
      else {// in case of gun/box simulations
        pPdg    = 0;
        pStatus = 0;
        break;
      }
    }//while
    
    AliDebug(2,"Converted photon/electron:");
    AliDebug(2,Form("\t Mother label %d, pdg %d, status %d",iMom, mPdg, mStatus));
    AliDebug(2,Form("\t Parent label %d, pdg %d, status %d",iParent, pPdg, pStatus));
    
  } // mother and parent are electron or photon and have status 0
  else if((mPdg == 22 || mPdg == 11) && mStatus == 0)
  {
    //Still a conversion but only one electron/photon generated. Just from hadrons but not decays.
    if(pPdg == 2112 ||  pPdg == 211  ||  pPdg == 321 ||
       pPdg == 2212 ||  pPdg == 130  ||  pPdg == 13 )
    {
      SetTagBit(tag,kMCConversion);
      iMom     = mom->GetFirstMother();
      
      if(iMom < 0) 
      {
        AliInfo(Form("pdg = %d, mother = %d, skip",pPdg,iMom));
      }
      else
      {
        mom      = mcevent->Particle(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        
        AliDebug(2,"Converted hadron:");
        AliDebug(2,Form("\t Mother label %d, pdg %d, status %d",iMom, mPdg, mStatus));
      }
    } // hadron converted
    
    //Comment for the next lines, we do not check the parent of the hadron for the moment.
    //iParent =  mom->GetFirstMother() ;
    //if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
    
    //GrandParent
    //if(iParent >= 0){
    //	parent = stack->Particle(iParent);
    //	pPdg = TMath::Abs(parent->GetPdgCode());
    //}
  }  	  
  // conversion into electrons/photons checked  	  
  
  // first check for typical charged particles
  if     (mPdg     ==    13) SetTagBit(tag,kMCMuon);
  else if(mPdg     ==   211) SetTagBit(tag,kMCPion);
  else if(mPdg     ==   321) SetTagBit(tag,kMCKaon);
  else if(mPdgSign ==  2212) SetTagBit(tag,kMCProton);
  else if(mPdgSign == -2212) SetTagBit(tag,kMCAntiProton);
  else if(mPdgSign ==  2112) SetTagBit(tag,kMCNeutron);
  else if(mPdgSign == -2112) SetTagBit(tag,kMCAntiNeutron);
  
  //check for pi0 and eta (shouldn't happen unless their decays were turned off)
  else if(mPdg == 111)
  {
    
    SetTagBit(tag,kMCPi0Decay);
    
    AliDebug(2,"First mother is directly pi0, not decayed by generator");
    
    CheckOverlapped2GammaDecayESD(labels,nlabels, iMom, mcevent, tag); //set to kMCPi0 if 2 gammas in same cluster
  }
  else if(mPdg == 221)
  {
    SetTagBit(tag,kMCEtaDecay);
    
    AliDebug(2,"First mother is directly eta, not decayed by generator");
    
    CheckOverlapped2GammaDecayESD(labels,nlabels, iMom, mcevent, tag); //set to kMCEta if 2 gammas in same cluster
  }
  //Photons  
  else if(mPdg == 22)
  {
    SetTagBit(tag,kMCPhoton);
    
    if(pPdg == 111)
    {
      SetTagBit(tag,kMCPi0Decay);
      
      AliDebug(2,"PYTHIA pi0 decay photon,  parent pi0 with status 11");
      
      CheckOverlapped2GammaDecayESD(labels,nlabels, iParent, mcevent, tag); //set to kMCPi0 if 2 gammas in same cluster
                                                                            // In case it did not merge, check if the decay companion is lost
      if(!CheckTagBit(tag, kMCPi0) && !CheckTagBit(tag,kMCDecayPairInCalo) && !CheckTagBit(tag,kMCDecayPairLost))
        CheckLostDecayPairESD(arrayCluster,iMom, iParent, mcevent, tag);
      
      //printf("Bit set is Merged %d, Pair in calo %d, Lost %d\n",CheckTagBit(tag, kMCPi0),CheckTagBit(tag,kMCDecayPairInCalo),CheckTagBit(tag,kMCDecayPairLost));
    }
    else if (pPdg == 221)
    {
      SetTagBit(tag, kMCEtaDecay);
      
      AliDebug(2,"PYTHIA eta decay photon,  parent pi0 with status 11");
      
      CheckOverlapped2GammaDecayESD(labels,nlabels, iParent, mcevent, tag);//set to kMCEta if 2 gammas in same cluster
                                                                           // In case it did not merge, check if the decay companion is lost
      if(!CheckTagBit(tag, kMCEta) && !CheckTagBit(tag,kMCDecayPairInCalo) && !CheckTagBit(tag,kMCDecayPairLost))
        CheckLostDecayPairESD(arrayCluster,iMom, iParent, mcevent, tag);
    }
    else if(mStatus == 1)
    { //undecayed particle
      if(fMCGenerator == kPythia)
      {
        if(iParent < 8 && iParent > 5)
        {//outgoing partons
          if(pPdg == 22) SetTagBit(tag,kMCPrompt);
          else           SetTagBit(tag,kMCFragmentation);
        }//Outgoing partons 
        else  if(iParent <= 5)
        {
          SetTagBit(tag, kMCISR); //Initial state radiation
        }
        else  SetTagBit(tag,kMCUnknown);
      }//PYTHIA
      
      else if(fMCGenerator == kHerwig)
      {
        if(pStatus < 197)
        {//Not decay
          while(1)
          {
            if(parent)
            {
              if(parent->GetFirstMother()<=5) break;
              iParent = parent->GetFirstMother();
              parent  = mcevent->Particle(iParent);
              pStatus = parent->GetStatusCode();
              pPdg = TMath::Abs(parent->GetPdgCode());
            } else break;
          }//Look for the parton
          
          if(iParent < 8 && iParent > 5)
          {
            if(pPdg == 22) SetTagBit(tag,kMCPrompt);
            else           SetTagBit(tag,kMCFragmentation);
          }
          else SetTagBit(tag,kMCISR);//Initial state radiation
        }//Not decay
        else  SetTagBit(tag,kMCUnknown);
      }//HERWIG
    }
    else  SetTagBit(tag,kMCOtherDecay);
    
  }//Mother Photon
  
  // Electron check.  Where did that electron come from?
  else if(mPdg == 11)
  { //electron
    if(pPdg == 11 && parent)
    {
      Int_t iGrandma = parent->GetFirstMother();
      if(iGrandma >= 0)
      {
        TParticle* gma = mcevent->Particle(iGrandma); //get mother
        Int_t gPdg = TMath::Abs(gma->GetPdgCode());
        
        if      (gPdg == 23) { SetTagBit(tag,kMCZDecay); } //parent is Z-boson
        else if (gPdg == 24) { SetTagBit(tag,kMCWDecay); } //parent is W-boson
      }
    }
    
    SetTagBit(tag,kMCElectron);
    
    AliDebug(1,"Checking ancestors of electrons");
    
    if      (pPdg == 111) { SetTagBit(tag,kMCPi0Decay); SetTagBit(tag, kMCDecayDalitz) ; } //Pi0 Dalitz decay
    else if (pPdg == 221) { SetTagBit(tag,kMCEtaDecay); SetTagBit(tag, kMCDecayDalitz) ; } //Eta Dalitz decay
    else if((499 < pPdg && pPdg < 600)||(4999 < pPdg && pPdg < 6000)) { SetTagBit(tag,kMCEFromB); } //b-->e decay
    else if((399 < pPdg && pPdg < 500)||(3999 < pPdg && pPdg < 5000))
    { //check charm decay
      if(parent)
      {
        Int_t iGrandma = parent->GetFirstMother();
        if(iGrandma >= 0)
        {
          TParticle* gma = mcevent->Particle(iGrandma); //get mother of charm
          Int_t gPdg = TMath::Abs(gma->GetPdgCode());
          if((499 < gPdg && gPdg < 600)||(4999 < gPdg && gPdg < 6000)) SetTagBit(tag,kMCEFromCFromB); //b-->c-->e
          else SetTagBit(tag,kMCEFromC); //c-->e
        }
        else SetTagBit(tag,kMCEFromC); //c-->e
      }//parent
    }
    else
    {
      //if it is not from any of the above, where is it from?
      if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
      
      else SetTagBit(tag,kMCOtherDecay);
      
      //if(parent) AliDebug(1,Form("Status %d Electron from other origin: %s (pPdg = %d) %s (mpdg = %d)",mStatus,parent->GetName(),pPdg,mom->GetName(),mPdg));
    }
  }//electron check
   //Cluster was made by something else
  else
  {
    AliDebug(2,Form("\t Setting kMCUnknown for cluster from %s (pdg = %d, Parent pdg = %d)",
                    mom->GetName(),mPdg,pPdg));
    
    SetTagBit(tag,kMCUnknown);
  }
  
  return tag;
}

//__________________________________________________________________________________________
/// \return tag with primary particle(S) at the origin of the cluster/track.
/// Do this for AODs, same things as in CheckOriginInStack.
///
/// Generally speaking, label is the MC label of a reconstructed
/// entity (track, cluster, etc) for which we want to know something 
/// about its heritage, but one can also use it directly with stack 
/// particles not connected to reconstructed entities.
//__________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInAOD(const Int_t *labels, Int_t nlabels,
                                           const AliMCEvent* mcevent, const TObjArray* arrayCluster)
{  
  Int_t tag = 0;
  
  // Most significant particle contributing to the cluster
  Int_t label=labels[0];
    
  // Mother
  AliAODMCParticle * mom = (AliAODMCParticle *) mcevent->GetTrack(label);
  Int_t iMom     = label;
  Int_t mPdgSign = mom->GetPdgCode();
  Int_t mPdg     = TMath::Abs(mPdgSign);
  Int_t iParent  = mom->GetMother() ;
  
  //if(label < 8 && fMCGenerator != kBoxLike) AliDebug(1,Form("Mother is parton %d\n",iParent));
  
  //GrandParent
  AliAODMCParticle * parent = NULL ;
  Int_t pPdg = -1;
  if(iParent >= 0)
  {
    parent = (AliAODMCParticle *) mcevent->GetTrack(iParent);
    pPdg = TMath::Abs(parent->GetPdgCode());
  }
  else AliDebug(1,Form("Parent with label %d",iParent));
  
  AliDebug(2,"Cluster most contributing mother and its parent:");
  AliDebug(2,Form("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary()));
  AliDebug(2,Form("\t Parent label %d, pdg %d, Primary? %d, Physical Primary? %d",iParent, pPdg, parent?parent->IsPrimary():-1, parent?parent->IsPhysicalPrimary():-1));
  
  //Check if mother is converted, if not, get the first non converted mother
  if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && !mom->IsPrimary())
  {
    SetTagBit(tag,kMCConversion);
    
    // Check if the mother is photon or electron with status not stable
    while ((pPdg == 22 || pPdg == 11) && !mom->IsPhysicalPrimary())
    {
      // Mother
      iMom  = mom->GetMother();
      
      if(iMom < 0) 
      {
        AliInfo(Form("pdg = %d, mother = %d, skip",pPdg,iMom));
        break;
      }
      
      mom      = (AliAODMCParticle *) mcevent->GetTrack(iMom);
      mPdgSign = mom->GetPdgCode();
      mPdg     = TMath::Abs(mPdgSign);
      iParent  = mom->GetMother() ;
      //if(label < 8 ) AliDebug(1, Form("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent));
      
      // GrandParent
      if(iParent >= 0 && parent)
      {
        parent = (AliAODMCParticle *) mcevent->GetTrack(iParent);
        pPdg = TMath::Abs(parent->GetPdgCode());
      }
      // printf("\t While Mother label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary());
      // printf("\t While Parent label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iParent, pPdg, parent->IsPrimary(), parent->IsPhysicalPrimary()); 
      
    }//while	
    
    AliDebug(2,"AliMCAnalysisUtils::CheckOriginInAOD() - Converted photon/electron:");
    AliDebug(2,Form("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary()));
    AliDebug(2,Form("\t Parent label %d, pdg %d, Primary? %d, Physical Primary? %d",iParent, pPdg, parent?parent->IsPrimary():-1, parent?parent->IsPhysicalPrimary():-1));
    
  } // mother and parent are electron or photon and have status 0 and parent is photon or electron
  else if((mPdg == 22 || mPdg == 11) && !mom->IsPrimary())
  {
    // Still a conversion but only one electron/photon generated. Just from hadrons
    if(pPdg == 2112 ||  pPdg == 211 ||  pPdg == 321 ||  
       pPdg == 2212 ||  pPdg == 130 ||  pPdg == 13 )
    {
      SetTagBit(tag,kMCConversion);
      iMom     = mom->GetMother();
      
      if(iMom < 0) 
      {
        AliInfo(Form("pdg = %d, mother = %d, skip",pPdg,iMom));
      }
      else
      {
        mom      = (AliAODMCParticle *) mcevent->GetTrack(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        
        AliDebug(2,"AliMCAnalysisUtils::CheckOriginInAOD() - Converted hadron:");
        AliDebug(2,Form("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary()));
      }
    } // hadron converted
    
    //Comment for next lines, we do not check the parent of the hadron for the moment.
    //iParent =  mom->GetMother() ;
    //if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
    
    //GrandParent
    //if(iParent >= 0){
    //	parent = (AliAODMCParticle *) mcevent->GetTrack(iParent);
    //	pPdg = TMath::Abs(parent->GetPdgCode());
    //}
  }  
  
  //printf("Final mother mPDG %d\n",mPdg);
  
  // conversion into electrons/photons checked  
  
  //first check for typical charged particles
  if     (mPdg     ==    13) SetTagBit(tag,kMCMuon);
  else if(mPdg     ==   211) SetTagBit(tag,kMCPion);
  else if(mPdg     ==   321) SetTagBit(tag,kMCKaon);
  else if(mPdgSign ==  2212) SetTagBit(tag,kMCProton);
  else if(mPdgSign ==  2112) SetTagBit(tag,kMCNeutron);
  else if(mPdgSign == -2212) SetTagBit(tag,kMCAntiProton);
  else if(mPdgSign == -2112) SetTagBit(tag,kMCAntiNeutron);
  
  //check for pi0 and eta (shouldn't happen unless their decays were turned off)
  else if(mPdg == 111)
  {
    SetTagBit(tag,kMCPi0Decay);
    
    AliDebug(2,"First mother is directly pi0, not decayed by generator");
    
    CheckOverlapped2GammaDecayAOD(labels,nlabels, iMom, mcevent, tag); //set to kMCPi0 if 2 gammas in same cluster
  }
  else if(mPdg == 221)
  {
    SetTagBit(tag,kMCEtaDecay);   
    
    AliDebug(2,"First mother is directly eta, not decayed by generator");
    
    CheckOverlapped2GammaDecayAOD(labels,nlabels, iMom, mcevent, tag); //set to kMCEta if 2 gammas in same cluster
  }
  //Photons  
  else if(mPdg == 22)
  {
    SetTagBit(tag,kMCPhoton);
    
    if(pPdg == 111)
    {
      SetTagBit(tag,kMCPi0Decay);
      
      AliDebug(2,"Generator pi0 decay photon");
      
      CheckOverlapped2GammaDecayAOD(labels,nlabels, iParent, mcevent, tag); //set to kMCPi0 if 2 gammas in same cluster
                                                                            // In case it did not merge, check if the decay companion is lost
      if(!CheckTagBit(tag, kMCPi0) && !CheckTagBit(tag,kMCDecayPairInCalo) && !CheckTagBit(tag,kMCDecayPairLost))
      {
        CheckLostDecayPairAOD(arrayCluster,iMom, iParent, mcevent, tag);
      }
      
      //printf("Bit set is Merged %d, Pair in calo %d, Lost %d\n",CheckTagBit(tag, kMCPi0),CheckTagBit(tag,kMCDecayPairInCalo),CheckTagBit(tag,kMCDecayPairLost));
    }
    else if (pPdg == 221)
    {
      SetTagBit(tag, kMCEtaDecay);
      
      AliDebug(2,"Generator eta decay photon");
      
      CheckOverlapped2GammaDecayAOD(labels,nlabels, iParent, mcevent, tag); //set to kMCEta if 2 gammas in same cluster
                                                                            // In case it did not merge, check if the decay companion is lost
      if(!CheckTagBit(tag, kMCEta) && !CheckTagBit(tag,kMCDecayPairInCalo) && !CheckTagBit(tag,kMCDecayPairLost))
        CheckLostDecayPairAOD(arrayCluster,iMom, iParent, mcevent, tag);
    }
    else if( mom->IsPhysicalPrimary() && ( fMCGenerator == kPythia || fMCGenerator == kHerwig ) ) //undecayed particle
    {
      if(iParent < 8 && iParent > 5 )
      {//outgoing partons
        if(pPdg == 22) SetTagBit(tag,kMCPrompt);
        else SetTagBit(tag,kMCFragmentation);
      }//Outgoing partons
      else if( iParent <= 5 && ( fMCGenerator == kPythia || fMCGenerator == kHerwig ) )
      {
        SetTagBit(tag, kMCISR); //Initial state radiation
      }
      else SetTagBit(tag,kMCUnknown);
    }//Physical primary
    else SetTagBit(tag,kMCOtherDecay);
    
  } // Mother Photon
  
  // Electron check.  Where did that electron come from?
  else if(mPdg == 11)
  { //electron
    if(pPdg == 11 && parent)
    {
      Int_t iGrandma = parent->GetMother();
      if(iGrandma >= 0)
      {
        AliAODMCParticle* gma = (AliAODMCParticle*) mcevent->GetTrack(iGrandma);
        Int_t gPdg = TMath::Abs(gma->GetPdgCode());
        
        if      (gPdg == 23) { SetTagBit(tag,kMCZDecay); } //parent is Z-boson
        else if (gPdg == 24) { SetTagBit(tag,kMCWDecay); } //parent is W-boson
      }
    }
    
    SetTagBit(tag,kMCElectron);
    
    AliDebug(1,"Checking ancestors of electrons");
    
    if      (pPdg == 111) { SetTagBit(tag,kMCPi0Decay); SetTagBit(tag,kMCDecayDalitz);} //Pi0 Dalitz decay
    else if (pPdg == 221) { SetTagBit(tag,kMCEtaDecay); SetTagBit(tag,kMCDecayDalitz);} //Eta Dalitz decay
    else if((499 < pPdg && pPdg < 600)||(4999 < pPdg && pPdg < 6000)) { SetTagBit(tag,kMCEFromB);} //b-hadron decay
    else if((399 < pPdg && pPdg < 500)||(3999 < pPdg && pPdg < 5000))
    { //c-hadron decay check
      if(parent)
      {
        Int_t iGrandma = parent->GetMother();
        if(iGrandma >= 0)
        {
          AliAODMCParticle* gma = (AliAODMCParticle*) mcevent->GetTrack(iGrandma); //charm's mother
          Int_t gPdg = TMath::Abs(gma->GetPdgCode());
          if((499 < gPdg && gPdg < 600)||(4999 < gPdg && gPdg < 6000)) SetTagBit(tag,kMCEFromCFromB); //b-->c-->e decay
          else SetTagBit(tag,kMCEFromC); //c-hadron decay
        }
        else SetTagBit(tag,kMCEFromC); //c-hadron decay
      }//parent
    } else
    { //prompt or other decay
      
      TParticlePDG* foo  = TDatabasePDG::Instance()->GetParticle(pPdg);
      TParticlePDG* foo1 = TDatabasePDG::Instance()->GetParticle(mPdg);
      
      AliDebug(1,Form("Electron from other origin: %s (pPdg = %d) %s (mPdg = %d)",foo->GetName(), pPdg,foo1->GetName(),mPdg));
      
      if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
      else             SetTagBit(tag,kMCOtherDecay);
    }      
  }// electron check
   // cluster was made by something else
  else
  {
    AliDebug(1,Form("\t Setting kMCUnknown for cluster with pdg = %d, Parent pdg = %d",mPdg,pPdg));
    SetTagBit(tag,kMCUnknown);
  }
  
  return tag;
}

//_________________________________________________________________________________________
/// Check if cluster is formed from the contribution of 2 decay photons from pi0 or eta. 
/// Input in ESD stack.
//_________________________________________________________________________________________
void AliMCAnalysisUtils::CheckOverlapped2GammaDecayESD(const Int_t *labels, Int_t nlabels, 
                                                       Int_t mesonIndex, const AliMCEvent* mcevent,
                                                       Int_t &tag)
{  
  if(labels[0] < 0 || labels[0] > mcevent->GetNumberOfTracks() || nlabels <= 1)
  {
    AliDebug(2,Form("Exit : label[0] %d, n primaries %d, nlabels %d",labels[0],mcevent->GetNumberOfTracks(), nlabels));
    return;
  }
  
  TParticle * meson = mcevent->Particle(mesonIndex);
  Int_t mesonPdg    = meson->GetPdgCode();
  if(mesonPdg!=111 && mesonPdg!=221)
  {
    AliWarning(Form("Wrong pi0/eta PDG : %d",mesonPdg));
    return;
  }
  
  AliDebug(2,Form("%s, label %d",meson->GetName(), mesonIndex));
  
  //Check if meson decayed into 2 daughters or if both were kept.
  if(meson->GetNDaughters() != 2)
  {
    AliDebug(2,Form("Not overalapped. Number of daughters is %d, not 2",meson->GetNDaughters()));
    return;
  }
  
  //Get the daughters
  Int_t iPhoton0 = meson->GetDaughter(0);
  Int_t iPhoton1 = meson->GetDaughter(1);
  TParticle *photon0 = mcevent->Particle(iPhoton0);
  TParticle *photon1 = mcevent->Particle(iPhoton1);
  
  //Check if both daughters are photons
  if(photon0->GetPdgCode() != 22 || photon1->GetPdgCode()!=22)
  {
    AliDebug(2,Form("Not overalapped. PDG:  daughter 1 = %d, of daughter 2 = %d",photon0->GetPdgCode(),photon1->GetPdgCode()));
    return;
  }
  
  AliDebug(2,Form("Daughter labels : photon0 = %d, photon1 = %d",iPhoton0,iPhoton1));
  
  //Check if both photons contribute to the cluster
  Bool_t okPhoton0 = kFALSE;
  Bool_t okPhoton1 = kFALSE;
  
  AliDebug(3,"Labels loop:");
  
  Bool_t conversion = kFALSE;
  
  for(Int_t i = 0; i < nlabels; i++)
  {
    AliDebug(3,Form("\t  at begin:label %d/%d: %d, ok? photon1 %d, photon2 %d", i+1, nlabels, labels[i], okPhoton0, okPhoton1));
    
    //If we already found both, break the loop
    if(okPhoton0 && okPhoton1) break;
    
    Int_t index = labels[i];
    if      (iPhoton0 == index)
    {
      okPhoton0 = kTRUE;
      continue;
    }
    else if (iPhoton1 == index)
    {
      okPhoton1 = kTRUE;
      continue;
    }
    
    //Trace back the mother in case it was a conversion
    
    if(index >= mcevent->GetNumberOfTracks())
    {
      AliWarning(Form("Particle index %d larger than size of list %d!!",index,mcevent->GetNumberOfTracks()));
      continue;
    }
    
    TParticle * daught = mcevent->Particle(index);
    Int_t tmpindex = daught->GetFirstMother();		
    
    AliDebug(3,Form("\t Conversion? : mother %d",tmpindex));
    
    while(tmpindex>=0)
    {
      //MC particle of interest is the mother
      AliDebug(3,Form("\t \t parent index %d",tmpindex));
      daught   = mcevent->Particle(tmpindex);
      if      (iPhoton0 == tmpindex)
      {
        conversion = kTRUE;
        okPhoton0  = kTRUE;
        break;
      }
      else if (iPhoton1 == tmpindex)
      {
        conversion = kTRUE;
        okPhoton1  = kTRUE;
        break;
      }
      
      tmpindex = daught->GetFirstMother();
      
    }//While to check if pi0/eta daughter was one of these contributors to the cluster
    
    //if(i == 0 && (!okPhoton0 && !okPhoton1))
    //  AliDebug(1,Form("Something happens, first label should be from a photon decay!"));
    
  }//loop on list of labels
  
  //If both photons contribute tag as the corresponding meson.
  if(okPhoton0 && okPhoton1)
  {
    AliDebug(2,Form("%s OVERLAPPED DECAY", meson->GetName()));
    
    if(!CheckTagBit(tag,kMCConversion) && conversion) SetTagBit(tag,kMCConversion) ;
    
    if(mesonPdg == 111) SetTagBit(tag,kMCPi0);
    else                SetTagBit(tag,kMCEta);
  }
}	

//_________________________________________________________________________________________
/// Check if cluster is formed from the contribution of 2 decay photons from pi0 or eta. 
/// Input are AOD AliAODMCParticles.
//_________________________________________________________________________________________
void AliMCAnalysisUtils::CheckOverlapped2GammaDecayAOD(const Int_t *labels, Int_t nlabels, 
                                                       Int_t mesonIndex, const AliMCEvent* mcevent, 
                                                       Int_t & tag)
{  
  if(labels[0] < 0 || labels[0] > mcevent->GetNumberOfTracks() || nlabels <= 1)
  {
    AliDebug(2,Form("Exit : label[0] %d, n primaries %d, nlabels %d",labels[0],mcevent->GetNumberOfTracks(), nlabels));
    return;
  }
  
  AliAODMCParticle * meson = (AliAODMCParticle *) mcevent->GetTrack(mesonIndex);
  Int_t mesonPdg = meson->GetPdgCode();
  if(mesonPdg != 111 && mesonPdg != 221)
  {
    AliWarning(Form("Wrong pi0/eta PDG : %d",mesonPdg));
    return;
  }
  
  AliDebug(2,Form("pdg %d, label %d, ndaughters %d", mesonPdg, mesonIndex, meson->GetNDaughters()));
  
  //Get the daughters
  if(meson->GetNDaughters() != 2)
  {
    AliDebug(2,Form("Not overalapped. Number of daughters is %d, not 2",meson->GetNDaughters()));
    return;
  }
  
  Int_t iPhoton0 = meson->GetDaughter(0);
  Int_t iPhoton1 = meson->GetDaughter(1);
  //if((iPhoton0 == -1) || (iPhoton1 == -1)){
  //	if(fDebug > 2) 
  //		printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Exit : Not overlapped. At least a daughter do not exists : d1 %d, d2 %d \n", iPhoton0, iPhoton1);
  //	return;
  //}	
  AliAODMCParticle *photon0 = (AliAODMCParticle *) mcevent->GetTrack(iPhoton0);
  AliAODMCParticle *photon1 = (AliAODMCParticle *) mcevent->GetTrack(iPhoton1);
    
  //Check if both daughters are photons
  if(photon0->GetPdgCode() != 22 && photon1->GetPdgCode()!=22)
  {
    AliWarning(Form("Not overlapped. PDG:  daughter 1 = %d, of daughter 2 = %d",photon0->GetPdgCode(),photon1->GetPdgCode()));
    return;
  }
  
  AliDebug(2,Form("Daughter labels : photon0 = %d, photon1 = %d",iPhoton0,iPhoton1));
  
  //Check if both photons contribute to the cluster
  Bool_t okPhoton0 = kFALSE;
  Bool_t okPhoton1 = kFALSE;
  
  AliDebug(3,"Labels loop:");
  
  Bool_t conversion = kFALSE;
  
  for(Int_t i = 0; i < nlabels; i++)
  {
    AliDebug(3, Form("\t label %d/%d: %d, ok? %d, %d", i, nlabels, labels[i], okPhoton0, okPhoton1));
    
    if(labels[i]<0) continue;
    
    //If we already found both, break the loop
    if(okPhoton0 && okPhoton1) break;
    
    Int_t index = 	labels[i];
    if      (iPhoton0 == index)
    {
      okPhoton0 = kTRUE;
      continue;
    }
    else if (iPhoton1 == index)
    {
      okPhoton1 = kTRUE;
      continue;
    }
    
    //Trace back the mother in case it was a conversion
    
    if(index >= mcevent->GetNumberOfTracks())
    {
      AliWarning(Form("Particle index %d larger than size of list %d!!",index,mcevent->GetNumberOfTracks()));
      continue;
    }
    
    AliAODMCParticle * daught = (AliAODMCParticle*) mcevent->GetTrack(index);
    Int_t tmpindex = daught->GetMother();
    AliDebug(3,Form("Conversion? : mother %d",tmpindex));
    
    while(tmpindex>=0){
      
      //MC particle of interest is the mother
      AliDebug(3,Form("\t parent index %d",tmpindex));
      daught   = (AliAODMCParticle*) mcevent->GetTrack(tmpindex);
      //printf("tmpindex %d\n",tmpindex);
      if      (iPhoton0 == tmpindex)
      {
        conversion = kTRUE;
        okPhoton0  = kTRUE;
        break;
      }
      else if (iPhoton1 == tmpindex)
      {
        conversion = kTRUE;
        okPhoton1  = kTRUE;
        break;
      }
      
      tmpindex = daught->GetMother();
      
    }//While to check if pi0/eta daughter was one of these contributors to the cluster
    
    //if(i == 0 && (!okPhoton0 && !okPhoton1)) AliDebug(1,"Something happens, first label should be from a photon decay!");
    
  }//loop on list of labels
  
  //If both photons contribute tag as the corresponding meson.
  if(okPhoton0 && okPhoton1)
  {
    AliDebug(2,Form("%s OVERLAPPED DECAY",(TDatabasePDG::Instance()->GetParticle(mesonPdg))->GetName()));
    
    if(!CheckTagBit(tag,kMCConversion) && conversion)
    {
      AliDebug(2,"Second decay photon produced a conversion");
      SetTagBit(tag,kMCConversion) ;
    }
    
    if(mesonPdg == 111) SetTagBit(tag,kMCPi0);
    else                SetTagBit(tag,kMCEta);
  }	
}

//______________________________________________________________________________________________________
/// Check on ESDs if the current decay photon has the second photon companion lost.
//______________________________________________________________________________________________________
void    AliMCAnalysisUtils::CheckLostDecayPairESD(const TObjArray* arrayCluster, Int_t iMom, Int_t iParent,
                                                  const AliMCEvent* mcevent, Int_t & tag)
{  
  if(!arrayCluster || iMom < 0 || iParent < 0|| !mcevent) return;
  
  TParticle * parent= mcevent->Particle(iParent);
  
  if(parent->GetNDaughters()!=2)
  {
    SetTagBit(tag, kMCDecayPairLost);
    return ;
  }
  
  Int_t pairLabel = -1;
  if     ( iMom != parent->GetDaughter(0) ) pairLabel = parent->GetDaughter(0);
  else if( iMom != parent->GetDaughter(1) ) pairLabel = parent->GetDaughter(1);
  
  if(pairLabel<0)
  {
    SetTagBit(tag, kMCDecayPairLost);
    return ;
  }
  
  for(Int_t iclus = 0; iclus < arrayCluster->GetEntriesFast(); iclus++)
  {
    AliVCluster * cluster = (AliVCluster*) arrayCluster->At(iclus);
    for(UInt_t ilab = 0; ilab< cluster->GetNLabels(); ilab++)
    {
      Int_t label = cluster->GetLabels()[ilab];
      if ( label==pairLabel )
      {
        SetTagBit(tag, kMCDecayPairInCalo);
        return ;
      }
      else if ( label== iParent || label== iMom )
      {
        continue;
      }
      else // check the ancestry
      {
        TParticle * mother = mcevent->Particle(label);
        
        if ( !mother )
        {
          AliInfo(Form("MC Mother not available for label %d",label));
          continue;
        }
        
        Int_t momPDG = TMath::Abs(mother->GetPdgCode());
        if ( momPDG!=11 && momPDG!=22 ) continue;
        
        // Check if "mother" of entity is converted, if not, get the first non converted mother
        Int_t iParentClus = mother->GetFirstMother();
        if ( iParentClus < 0 ) continue;
        
        TParticle * parentClus = mcevent->Particle(iParentClus);
        if ( !parentClus ) continue;
        
        Int_t parentClusPDG    = TMath::Abs(parentClus->GetPdgCode());
        Int_t parentClusStatus = parentClus->GetStatusCode();
        
        if( parentClusPDG != 22 && parentClusPDG != 11 && parentClusStatus != 0 ) continue;
        
        //printf("Conversion\n");
        
        // Check if the mother is photon or electron with status not stable
        while ( (parentClusPDG == 22 || parentClusPDG == 11) && parentClusStatus != 1 )
        {
          //New Mother
          label            = iParentClus;
          momPDG           = parentClusPDG;
          
          iParentClus      = parentClus->GetFirstMother();
          if(iParentClus < 0) break;
          
          parentClus       = mcevent->Particle(iParentClus);
          if(!parentClus) break;
          
          parentClusPDG    = TMath::Abs(parentClus->GetPdgCode());
          parentClusStatus = parentClus->GetStatusCode() ;
        }//while
        
        if ( (momPDG == 22 || parentClusPDG ==22) && (label==pairLabel || iParentClus == pairLabel) )
        {
          SetTagBit(tag, kMCDecayPairInCalo);
          //printf("Conversion is paired\n");
          return ;
        }
        else continue;
        
      }
    }
  } // cluster loop
  
  SetTagBit(tag, kMCDecayPairLost);
}

//________________________________________________________________________________________________________
/// Check on AODs if the current decay photon has the second photon companion lost.
//________________________________________________________________________________________________________
void    AliMCAnalysisUtils::CheckLostDecayPairAOD(const TObjArray* arrayCluster, Int_t iMom, Int_t iParent,
                                                  const AliMCEvent* mcevent, Int_t & tag)
{
  if(!arrayCluster || iMom < 0 || iParent < 0|| !mcevent) return;
  
  AliAODMCParticle * parent = (AliAODMCParticle*) mcevent->GetTrack(iParent);
  
  //printf("*** Check label %d with parent %d\n",iMom, iParent);
  
  if(parent->GetNDaughters()!=2)
  {
    SetTagBit(tag, kMCDecayPairLost);
    //printf("\t ndaugh = %d\n",parent->GetNDaughters());
    return ;
  }
  
  Int_t pairLabel = -1;
  if     ( iMom != parent->GetDaughter(0) ) pairLabel = parent->GetDaughter(0);
  else if( iMom != parent->GetDaughter(1) ) pairLabel = parent->GetDaughter(1);
  
  if(pairLabel<0)
  {
    //printf("\t pair Label not found = %d\n",pairLabel);
    SetTagBit(tag, kMCDecayPairLost);
    return ;
  }
  
  //printf("\t *** find pair %d\n",pairLabel);
  
  for(Int_t iclus = 0; iclus < arrayCluster->GetEntriesFast(); iclus++)
  {
    AliVCluster * cluster = (AliVCluster*) arrayCluster->At(iclus);
    //printf("\t \t ** Cluster %d, nlabels %d\n",iclus,cluster->GetNLabels());
    for(UInt_t ilab = 0; ilab< cluster->GetNLabels(); ilab++)
    {
      Int_t label = cluster->GetLabels()[ilab];
      
      //printf("\t \t label %d\n",label);
      
      if ( label==pairLabel )
      {
        //printf("\t \t Pair found\n");
        SetTagBit(tag, kMCDecayPairInCalo);
        return ;
      }
      else if ( label== iParent || label== iMom )
      {
        //printf("\t \t skip\n");
        continue;
      }
      else // check the ancestry
      {
        AliAODMCParticle * mother = (AliAODMCParticle*) mcevent->GetTrack(label);
        
        if ( !mother )
        {
          AliInfo(Form("MC Mother not available for label %d",label));
          continue;
        }
        
        Int_t momPDG = TMath::Abs(mother->GetPdgCode());
        if ( momPDG!=11 && momPDG!=22 ) continue;
        
        // Check if "mother" of entity is converted, if not, get the first non converted mother
        Int_t iParentClus = mother->GetMother();
        if(iParentClus < 0) continue;
        
        AliAODMCParticle * parentClus =  (AliAODMCParticle*) mcevent->GetTrack(iParentClus);
        if(!parentClus) continue;
        
        Int_t parentClusPDG    = TMath::Abs(parentClus->GetPdgCode());
        Int_t parentClusStatus = parentClus->GetStatus();
        
        if ( parentClusPDG != 22 && parentClusPDG != 11 && parentClusStatus != 0 )
        {
          //printf("\t \t skip, not a conversion, parent: pdg %d, status %d\n",parentClusPDG,parentClusStatus);
          continue;
        }
        
        //printf("\t \t Conversion\n");
        
        // Check if the mother is photon or electron with status not stable
        while ((parentClusPDG == 22 || parentClusPDG == 11) && parentClusStatus != 1)
        {
          //New Mother
          label            = iParentClus;
          momPDG           = parentClusPDG;
          
          iParentClus      = parentClus->GetMother();
          if ( iParentClus < 0 ) break;
          
          parentClus       =  (AliAODMCParticle*) mcevent->GetTrack(iParentClus);
          if ( !parentClus ) break;
          
          parentClusPDG    = TMath::Abs(parentClus->GetPdgCode());
          parentClusStatus = parentClus->GetStatus() ;
        }//while
        
        if ( (momPDG == 22 || parentClusPDG ==22) && (label==pairLabel || iParentClus == pairLabel) )
        {
          SetTagBit(tag, kMCDecayPairInCalo);
          //printf("\t \t Conversion is paired: mom %d, parent %d\n",label,iParentClus);
          return ;
        }
        else
        {
          //printf("\t \t Skip, finally label %d, pdg %d, parent label %d, pdg %d, status %d\n",label,momPDG,iParentClus,parentClusPDG,parentClusStatus);
          continue;
        }
        
      }
    }
  } // cluster loop
  
  SetTagBit(tag, kMCDecayPairLost);
}

//_____________________________________________________________________
/// \return list of jets (TParticles) and index of most likely parton that originated it.
//_____________________________________________________________________
TList * AliMCAnalysisUtils::GetJets(const AliCaloTrackReader * reader)
{  
  AliMCEvent * mcevent = reader->GetMC();
  Int_t iEvent = reader->GetEventNumber();	
  AliGenEventHeader * geh = reader->GetGenEventHeader();
  
  if(fCurrentEvent!=iEvent)
  {
    fCurrentEvent = iEvent;
    fJetsList = new TList;
    Int_t nTriggerJets = 0;
    Float_t tmpjet[]={0,0,0,0};
		
    //printf("Event %d %d\n",fCurrentEvent,iEvent);
    //Get outgoing partons
    if(mcevent->GetNumberOfTracks() < 8) return fJetsList;
    TParticle * parton1 =  mcevent->Particle(6);
    TParticle * parton2 =  mcevent->Particle(7);
    
    AliDebug(2,Form("Parton 6 : %s, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f",
                    parton1->GetName(),parton1->Pt(),parton1->Energy(),parton1->Phi()*TMath::RadToDeg(),parton1->Eta()));
    AliDebug(2,Form("Parton 7 : %s, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f",
                    parton2->GetName(),parton2->Pt(),parton2->Energy(),parton2->Phi()*TMath::RadToDeg(),parton2->Eta()));
    
    // 		//Trace the jet from the mother parton
    // 		Float_t pt  = 0;
    // 		Float_t pt1 = 0;
    // 		Float_t pt2 = 0;
    // 		Float_t e   = 0;
    // 		Float_t e1  = 0;
    // 		Float_t e2  = 0;
    // 		TParticle * tmptmp = new TParticle;
    // 		for(Int_t i = 0; i< stack->GetNprimary(); i++){
    // 			tmptmp = stack->Particle(i);
		
    // 			if(tmptmp->GetStatusCode() == 1){
    // 				pt = tmptmp->Pt();
    // 				e =  tmptmp->Energy();			
    // 				Int_t imom = tmptmp->GetMother();
    // 				Int_t imom1 = 0;
    // 				//printf("1st imom %d\n",imom);
    // 				while(imom > 5){
    // 					imom1=imom;
    // 					tmptmp = stack->Particle(imom);
    // 					imom = tmptmp->GetMother();
    // 					//printf("imom %d	\n",imom);
    // 				}
    // 				//printf("Last imom %d %d\n",imom1, imom);
    // 				if(imom1 == 6) {
    // 					pt1+=pt;
    // 					e1+=e;				
    // 				}
    // 				else if (imom1 == 7){
    // 					pt2+=pt;
    // 					e2+=e;					}
    // 			}// status 1
    
    // 		}// for
		
    // 		printf("JET 1, pt %2.2f, e %2.2f; JET 2, pt %2.2f, e %2.2f \n",pt1,e1,pt2,e2);
		
		//Get the jet, different way for different generator
		//PYTHIA
    if(fMCGenerator == kPythia)
    {
      TParticle * jet =  0x0;
      AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) geh;
      nTriggerJets =  pygeh->NTriggerJets();
      AliDebug(2,Form("PythiaEventHeader: Njets: %d",nTriggerJets));
      
      for(Int_t i = 0; i< nTriggerJets; i++)
      {
        pygeh->TriggerJet(i, tmpjet);
        jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
        //Assign an outgoing parton as mother
        Float_t phidiff1 = TMath::Abs(jet->Phi()-parton1->Phi());		
        Float_t phidiff2 = TMath::Abs(jet->Phi()-parton2->Phi());
        if(phidiff1 > phidiff2) jet->SetFirstMother(7);
        else  jet->SetFirstMother(6);
        //jet->Print();
        AliDebug(1,Form("PYTHIA Jet %d: mother %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f",
                        i, jet->GetFirstMother(),jet->Pt(),jet->Energy(),jet->Phi()*TMath::RadToDeg(),jet->Eta()));
        fJetsList->Add(jet);			
      }
    }//Pythia triggered jets
    //HERWIG
    else if (fMCGenerator == kHerwig)
    {
      Int_t pdg = -1;		
      //Check parton 1
      TParticle * tmp = parton1;
      if(parton1->GetPdgCode()!=22)
      {
        while(pdg != 94){
          if(tmp->GetFirstDaughter()==-1) return fJetsList;
          tmp = mcevent->Particle(tmp->GetFirstDaughter());
          pdg = tmp->GetPdgCode();
        }//while
        
        //Add found jet to list
        TParticle *jet1 = new TParticle(*tmp);
        jet1->SetFirstMother(6);
        fJetsList->Add(jet1);
        //printf("jet 1:  first daughter %d, last daughter %d\n", tmp->GetFirstDaughter(), tmp->GetLastDaughter());
        //tmp = stack->Particle(tmp->GetFirstDaughter());
        //tmp->Print();
        //jet1->Print();
        AliDebug(1,Form("HERWIG Jet 1: mother %d, status %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f",
                        jet1->GetFirstMother(),jet1->GetStatusCode(),jet1->Pt(),jet1->Energy(),jet1->Phi()*TMath::RadToDeg(),jet1->Eta()));
      }//not photon
      
      //Check parton 2
      pdg = -1;
      tmp = parton2;
      if(parton2->GetPdgCode()!=22)
      {
        while(pdg != 94)
        {
          if(tmp->GetFirstDaughter()==-1) return fJetsList;
          tmp = mcevent->Particle(tmp->GetFirstDaughter());
          pdg = tmp->GetPdgCode();
        }//while
        
        //Add found jet to list
        TParticle *jet2 = new TParticle(*tmp);
        jet2->SetFirstMother(7);
        fJetsList->Add(jet2);
        //jet2->Print();
        AliDebug(2,Form("HERWIG Jet 2: mother %d, status %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f",
                        jet2->GetFirstMother(),jet2->GetStatusCode(),jet2->Pt(),jet2->Energy(),jet2->Phi()*TMath::RadToDeg(),jet2->Eta()));
        //Int_t first =  tmp->GetFirstDaughter();
        //Int_t last  =  tmp->GetLastDaughter();
        //printf("jet 2:  first daughter %d, last daughter %d, pdg %d\n",first, last, tmp->GetPdgCode());
				//	for(Int_t d = first ; d < last+1; d++){
        //						tmp = stack->Particle(d);
        //						if(i == tmp->GetMother())
        //							printf("Daughter n %d, Mother %d, name %s, status %d, pT %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
        //							d,tmp->GetMother(), tmp->GetName(), tmp->GetStatusCode(),tmp->Pt(),tmp->Energy(),tmp->Phi()*TMath::RadToDeg(),tmp->Eta());			   
        //			   }
        //tmp->Print();
      }//not photon
    }//Herwig generated jets
  }
  
  return fJetsList;
}

//______________________________________________________________________________
/// \return the kinematics of the particle that generated the signal, its pdg and its status and its label mother.
//______________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetDaughter(Int_t idaugh, Int_t label,
                                               const AliCaloTrackReader* reader,
                                               Int_t & pdg, Int_t & status, 
                                               Bool_t & ok, Int_t & daughlabel, TVector3 & prodVertex)
{
  fDaughMom.SetPxPyPzE(0,0,0,0);
  
  if(!reader->GetMC())
  {
    AliWarning("MCEvent is not available, check analysis settings in configuration file!!");
    
    ok = kFALSE;
    return fDaughMom;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if(label < 0 || label >= nprimaries)
  {
    ok = kFALSE;
    return fDaughMom;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    daughlabel = momP->GetDaughter(idaugh);
    
    if(daughlabel < 0 || daughlabel >= nprimaries)
    {
      ok = kFALSE;
      return fDaughMom;
    }
    
    TParticle * daughP = reader->GetMC()->Particle(daughlabel);
    daughP->Momentum(fDaughMom);
    pdg    = daughP->GetPdgCode();
    status = daughP->GetStatusCode();
    prodVertex.SetXYZ(daughP->Vx(),daughP->Vy(),daughP->Vz());
  }
  else if(reader->ReadAODMCParticles())
  {
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    daughlabel              = momP->GetDaughter(idaugh);
    
    if(daughlabel < 0 || daughlabel >= nprimaries)
    {
      ok = kFALSE;
      return fDaughMom;
    }
    
    AliAODMCParticle * daughP = (AliAODMCParticle *) reader->GetMC()->GetTrack(daughlabel);
    fDaughMom.SetPxPyPzE(daughP->Px(),daughP->Py(),daughP->Pz(),daughP->E());
    pdg    = daughP->GetPdgCode();
    status = daughP->GetStatus();
    prodVertex.SetXYZ(daughP->Xv(),daughP->Yv(),daughP->Zv());
  }
  
  ok = kTRUE;
  
  return fDaughMom;
}

//______________________________________________________________________________________________________
/// \return the kinematics of the particle that generated the signal.
//______________________________________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetMother(Int_t label, const AliCaloTrackReader* reader, Bool_t & ok)
{  
  Int_t pdg = -1; Int_t status = -1; Int_t momlabel = -1;
  
  return GetMother(label,reader,pdg,status, ok,momlabel);
}

//_________________________________________________________________________________________
// \return the kinematics of the particle that generated the signal.
//_________________________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetMother(Int_t label, const AliCaloTrackReader* reader,
                                             Int_t & pdg, Int_t & status, Bool_t & ok)
{  
  Int_t momlabel = -1;
  
  return GetMother(label,reader,pdg,status, ok,momlabel);
}

//______________________________________________________________________________________________________
/// \return the kinematics of the particle that generated the signal, its pdg and its status and its label mother.
//______________________________________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetMother(Int_t label, const AliCaloTrackReader* reader, 
                                             Int_t & pdg, Int_t & status, Bool_t & ok, Int_t & momlabel)
{  
  fMotherMom.SetPxPyPzE(0,0,0,0);
  
  if(!reader->GetMC()) 
  {
    AliWarning("MCEvent is not available, check analysis settings in configuration file, STOP!!");
    
    ok = kFALSE;
    return fMotherMom;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if(label < 0 || label >= nprimaries)
  {
    ok = kFALSE;
    return fMotherMom;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    momP->Momentum(fMotherMom);
    pdg      = momP->GetPdgCode();
    status   = momP->GetStatusCode();
    momlabel = momP->GetFirstMother();
  }
  else if(reader->ReadAODMCParticles())
  {
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    fMotherMom.SetPxPyPzE(momP->Px(),momP->Py(),momP->Pz(),momP->E());
    pdg      = momP->GetPdgCode();
    status   = momP->GetStatus();
    momlabel = momP->GetMother();
  }
  
  ok = kTRUE;
  
  return fMotherMom;
}

//___________________________________________________________________________________
/// \return the kinematics of the particle that generated the signal.
//___________________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetMotherWithPDG(Int_t label, Int_t pdg,
                                                    const AliCaloTrackReader* reader,
                                                    Bool_t & ok, Int_t & momlabel)
{  
  fGMotherMom.SetPxPyPzE(0,0,0,0);
  
  if ( !reader->GetMC() )
  {
    AliWarning("MCEvent is not available, check analysis settings in configuration file!!");
    
    ok = kFALSE;
    return fGMotherMom;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if ( label < 0 || label >= nprimaries )
  {  
    ok = kFALSE;
    return fGMotherMom;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    
    if(momP->GetPdgCode()==pdg)
    {
      AliDebug(2,"PDG of mother is already the one requested!");
      fGMotherMom.SetPxPyPzE(momP->Px(),momP->Py(),momP->Pz(),momP->Energy());
      
      ok=kTRUE;
      return fGMotherMom;
    }
    
    Int_t grandmomLabel = momP->GetFirstMother();
    Int_t grandmomPDG   = -1;
    TParticle * grandmomP = 0x0;
   
    while (grandmomLabel >=0 ) 
    {
      grandmomP   = reader->GetMC()->Particle(grandmomLabel);
      grandmomPDG = grandmomP->GetPdgCode();
      
      if(grandmomPDG==pdg)
      {
        momlabel = grandmomLabel;
        fGMotherMom.SetPxPyPzE(grandmomP->Px(),grandmomP->Py(),grandmomP->Pz(),grandmomP->Energy());
        break;
      }
      
      grandmomLabel =  grandmomP->GetFirstMother();
    }
    
    if(grandmomPDG!=pdg) AliInfo(Form("Mother with PDG %d, NOT found! \n",pdg));
  }
  else if(reader->ReadAODMCParticles())
  {    
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    
    if(momP->GetPdgCode()==pdg)
    {
      AliDebug(2,"PDG of mother is already the one requested!");
      fGMotherMom.SetPxPyPzE(momP->Px(),momP->Py(),momP->Pz(),momP->E());
      
      ok=kTRUE;
      return fGMotherMom;
    }
    
    Int_t grandmomLabel = momP->GetMother();
    Int_t grandmomPDG   = -1;
    AliAODMCParticle * grandmomP = 0x0;
    
    while (grandmomLabel >=0 ) 
    {
      grandmomP   = (AliAODMCParticle *) reader->GetMC()->GetTrack(grandmomLabel);
      grandmomPDG = grandmomP->GetPdgCode();
      if(grandmomPDG==pdg)
      {
        //printf("AliMCAnalysisUtils::GetMotherWithPDG(AOD) - mother with PDG %d FOUND! \n",pdg);
        momlabel = grandmomLabel;
        fGMotherMom.SetPxPyPzE(grandmomP->Px(),grandmomP->Py(),grandmomP->Pz(),grandmomP->E());
        break;
      }
      
      grandmomLabel =  grandmomP->GetMother();
    }
    
    if(grandmomPDG!=pdg) AliInfo(Form("Mother with PDG %d, NOT found!",pdg));
  }
  
  ok = kTRUE;
  
  return fGMotherMom;
}

//______________________________________________________________________________________________
/// \return the kinematics of the particle that generated the signal.
//______________________________________________________________________________________________
TLorentzVector AliMCAnalysisUtils::GetGrandMother(Int_t label, const AliCaloTrackReader* reader,
                                                  Int_t & pdg, Int_t & status, Bool_t & ok,
                                                  Int_t & grandMomLabel, Int_t & greatMomLabel)
{
  fGMotherMom.SetPxPyPzE(0,0,0,0);
  
  if ( !reader->GetMC() )
  {
    AliWarning("MCEvent is not available, check analysis settings in configuration file, STOP!!");
    
    ok = kFALSE;
    return fGMotherMom;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if ( label < 0 || label >= nprimaries )
  {  
    ok = kFALSE;
    return fGMotherMom;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    
    grandMomLabel = momP->GetFirstMother();
    
    TParticle * grandmomP = 0x0;
    
    if (grandMomLabel >=0 )
    {
      grandmomP   = reader->GetMC()->Particle(grandMomLabel);
      pdg    = grandmomP->GetPdgCode();
      status = grandmomP->GetStatusCode();
      
      fGMotherMom.SetPxPyPzE(grandmomP->Px(),grandmomP->Py(),grandmomP->Pz(),grandmomP->Energy());
      greatMomLabel =  grandmomP->GetFirstMother();
    }
  }
  else if(reader->ReadAODMCParticles())
  {
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    
    grandMomLabel = momP->GetMother();
    
    AliAODMCParticle * grandmomP = 0x0;
    
    if(grandMomLabel >=0 )
    {
      grandmomP   = (AliAODMCParticle *) reader->GetMC()->GetTrack(grandMomLabel);
      pdg    = grandmomP->GetPdgCode();
      status = grandmomP->GetStatus();
      
      fGMotherMom.SetPxPyPzE(grandmomP->Px(),grandmomP->Py(),grandmomP->Pz(),grandmomP->E());
      greatMomLabel =  grandmomP->GetMother();
      
    }
  }
  
  ok = kTRUE;
  
  return fGMotherMom;
}

//_______________________________________________________________________________________________________________
/// In case of an eta or pi0 decay into 2 photons, get the asymmetry in the energy of the photons.
//_______________________________________________________________________________________________________________
void AliMCAnalysisUtils::GetMCDecayAsymmetryAngleForPDG(Int_t label, Int_t pdg, const AliCaloTrackReader* reader,
                                                        Float_t & asym, Float_t & angle, Bool_t & ok)
{  
  if(!reader->GetMC())
  {
    AliWarning("MCEvent is not available, check analysis settings in configuration file, STOP!!");
    
    ok = kFALSE;
    return;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if ( label < 0 || label >= nprimaries )
  {  
    ok = kFALSE;
    return ;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    
    Int_t grandmomLabel = momP->GetFirstMother();
    Int_t grandmomPDG   = -1;
    TParticle * grandmomP = 0x0;
   
    while (grandmomLabel >=0 ) 
    {
      grandmomP   = reader->GetMC()->Particle(grandmomLabel);
      grandmomPDG = grandmomP->GetPdgCode();
      
      if(grandmomPDG==pdg) break;
      
      grandmomLabel =  grandmomP->GetFirstMother();
    }
    
    if(grandmomPDG==pdg && grandmomP->GetNDaughters()==2) 
    {
      TParticle * d1 = reader->GetMC()->Particle(grandmomP->GetDaughter(0));
      TParticle * d2 = reader->GetMC()->Particle(grandmomP->GetDaughter(1));
      
      if(d1->GetPdgCode() == 22 && d1->GetPdgCode() == 22)
      {
        asym = (d1->Energy()-d2->Energy())/grandmomP->Energy();
        d1->Momentum(fDaughMom );
        d2->Momentum(fDaughMom2);
        angle = fDaughMom.Angle(fDaughMom2.Vect());
      }
    }
    else 
    {
      ok = kFALSE;
      AliInfo(Form("Mother with PDG %d, not found!",pdg));
      return;
    }
  }
  else if(reader->ReadAODMCParticles())
  {    
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    
    Int_t grandmomLabel = momP->GetMother();
    Int_t grandmomPDG   = -1;
    AliAODMCParticle * grandmomP = 0x0;
   
    while (grandmomLabel >=0 ) 
    {
      grandmomP   = (AliAODMCParticle *) reader->GetMC()->GetTrack(grandmomLabel);
      grandmomPDG = grandmomP->GetPdgCode();
      
      if(grandmomPDG==pdg) break;
      
      grandmomLabel =  grandmomP->GetMother();
    }
    
    if(grandmomPDG==pdg && grandmomP->GetNDaughters()==2) 
    {
      AliAODMCParticle * d1 = (AliAODMCParticle *) reader->GetMC()->GetTrack(grandmomP->GetDaughter(0));
      AliAODMCParticle * d2 = (AliAODMCParticle *) reader->GetMC()->GetTrack(grandmomP->GetDaughter(1));
      
      if(d1->GetPdgCode() == 22 && d1->GetPdgCode() == 22)
      {
        asym = (d1->E()-d2->E())/grandmomP->E();
        fDaughMom .SetPxPyPzE(d1->Px(),d1->Py(),d1->Pz(),d1->E());
        fDaughMom2.SetPxPyPzE(d2->Px(),d2->Py(),d2->Pz(),d2->E());
        angle = fDaughMom.Angle(fDaughMom2.Vect());
      }
    }
    else 
    {
      ok = kFALSE;
      AliInfo(Form("Mother with PDG %d, not found! \n",pdg));
      return;
    }      
  }
  
  ok = kTRUE;
}

//_________________________________________________________________________________________________
/// \return the the number of daughters of a given MC particle.
//_________________________________________________________________________________________________
Int_t AliMCAnalysisUtils::GetNDaughters(Int_t label, const AliCaloTrackReader* reader, Bool_t & ok)
{  
  if ( !reader->GetMC() )
  {
    AliWarning("Stack is not available, check analysis settings in configuration file, STOP!!");
    
    ok=kFALSE;
    return -1;
  }
  
  Int_t nprimaries = reader->GetMC()->GetNumberOfTracks();
  if ( label < 0 || label >= nprimaries )
  {
    ok = kFALSE;
    return -1;
  }
  
  if(reader->ReadStack())
  {
    TParticle * momP = reader->GetMC()->Particle(label);
    
    ok=kTRUE;
    
    return momP->GetNDaughters();
  }
  else if(reader->ReadAODMCParticles())
  {    
    AliAODMCParticle * momP = (AliAODMCParticle *) reader->GetMC()->GetTrack(label);
    
    ok = kTRUE;
    
    return momP->GetNDaughters();
  }
  
  ok = kFALSE;
  
  return -1;
}

//_________________________________________________________________________________
/// Compare the primary depositing more energy with the rest,
/// if no photon/electron (conversion) or neutral meson as comon ancestor, 
/// consider it as other particle contributing.
/// Give as input the meson label in case it was a pi0 or eta merged cluster.
//_________________________________________________________________________________
Int_t AliMCAnalysisUtils::GetNOverlaps(const Int_t * label, UInt_t nlabels,
                                       Int_t mctag, Int_t mesonLabel,
                                       AliCaloTrackReader * reader, 
                                       Int_t *overpdg, Int_t *overlabel)
{  
  Int_t ancPDG = 0, ancStatus = -1;
  TVector3 prodVertex;
  Int_t ancLabel = 0;
  Int_t noverlaps = 0;
  Bool_t ok = kFALSE;
  
  for (UInt_t ilab = 1; ilab < nlabels; ilab++ )
  {
    ancLabel = CheckCommonAncestor(label[0],label[ilab],reader,ancPDG,ancStatus,fMotherMom,prodVertex);
    
    //printf("Overlaps, i %d: Main Label %d, second label %d, ancestor: Label %d, pdg %d - tag %d \n",
    //       ilab,label[0],label[ilab],ancLabel,ancPDG, mctag);
    
    Bool_t overlap = kFALSE;
    
    if     ( ancLabel < 0 )
    {
      overlap = kTRUE;
      //printf("\t \t \t No Label = %d\n",ancLabel);
    }
    else if ( ( ancPDG==111 || ancPDG==221 ) && 
              ( CheckTagBit(mctag,kMCPi0) ||  CheckTagBit(mctag,kMCEta) ) && 
              ( (mesonLabel != ancLabel) && mesonLabel >=0 ) ) // in case the label is not provided check it is larger than 0
    {
      //printf("\t \t  meson Label %d, ancestor Label %d\n",mesonLabel,ancLabel);
      overlap = kTRUE;
    }
    else if( ancPDG!=22 && TMath::Abs(ancPDG)!=11 && ancPDG != 111 && ancPDG != 221 )
    {
      //printf("\t \t \t Non EM PDG = %d\n",ancPDG);
      overlap = kTRUE ;
    }
    
    if( !overlap ) continue ;
    
    // We have at least one overlap
    
    //printf("Overlap!!!!!!!!!!!!!!\n");
    
    noverlaps++;
    
    // What is the origin of the overlap?
    Bool_t  mOK = 0,      gOK = 0;
    Int_t   mpdg = -999999,  gpdg = -1;
    Int_t   mstatus = -1, gstatus = -1;
    Int_t   gLabel = -1, ggLabel = -1;
    
    GetMother     (label[ilab],reader,mpdg,mstatus,mOK);
    fGMotherMom =
    GetGrandMother(label[ilab],reader,gpdg,gstatus,gOK, gLabel,ggLabel);
    
    //printf("\t Overlap!, mother pdg %d; grand mother pdg %d",mpdg,gpdg);
    
    if( ( mpdg == 22 || TMath::Abs(mpdg==11) ) &&
        ( gpdg == 22 || TMath::Abs(gpdg==11) ) &&
       gLabel >=0 )
    {
      Int_t labeltmp = gLabel;
      while( ( gpdg == 22 || TMath::Abs(gpdg==11) ) && gLabel >=0 )
      {
        mpdg=gpdg;
        fGMotherMom = GetGrandMother(labeltmp,reader,gpdg,gstatus,ok, gLabel,ggLabel);
        labeltmp=gLabel;
      }
    }
    overpdg  [noverlaps-1] = mpdg;
    overlabel[noverlaps-1] = label[ilab];
  }
  
  return noverlaps ;
}

//________________________________________________________
/// Print some relevant parameters set for the analysis.
//________________________________________________________
void AliMCAnalysisUtils::Print(const Option_t * opt) const
{  
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  
  printf("Debug level    = %d\n",fDebug);
  printf("MC Generator   = %s\n",fMCGeneratorString.Data());
  printf(" \n");
} 

//__________________________________________________
/// Print the assigned origins to this particle.
//__________________________________________________
void AliMCAnalysisUtils::PrintMCTag(Int_t tag) const
{  
  printf("AliMCAnalysisUtils::PrintMCTag() - tag %d \n    photon %d, conv %d, prompt %d, frag %d, isr %d, \n    pi0 decay %d, eta decay %d, other decay %d  pi0 %d,  eta %d \n    electron %d, muon %d,pion %d, proton %d, neutron %d, \n    kaon %d, a-proton %d, a-neutron %d, unk %d, bad %d\n",
         tag,
         CheckTagBit(tag,kMCPhoton),
         CheckTagBit(tag,kMCConversion),
         CheckTagBit(tag,kMCPrompt),
         CheckTagBit(tag,kMCFragmentation),
         CheckTagBit(tag,kMCISR),
         CheckTagBit(tag,kMCPi0Decay),
         CheckTagBit(tag,kMCEtaDecay),
         CheckTagBit(tag,kMCOtherDecay),
         CheckTagBit(tag,kMCPi0),
         CheckTagBit(tag,kMCEta),
         CheckTagBit(tag,kMCElectron),
         CheckTagBit(tag,kMCMuon), 
         CheckTagBit(tag,kMCPion),
         CheckTagBit(tag,kMCProton), 
         CheckTagBit(tag,kMCAntiNeutron),
         CheckTagBit(tag,kMCKaon), 
         CheckTagBit(tag,kMCAntiProton), 
         CheckTagBit(tag,kMCAntiNeutron),
         CheckTagBit(tag,kMCUnknown),
         CheckTagBit(tag,kMCBadLabel)
         );
} 

//__________________________________________________
/// Set the generator type.
//__________________________________________________
void AliMCAnalysisUtils::SetMCGenerator(Int_t mcgen)
{  
  fMCGenerator = mcgen ;
  if     (mcgen == kPythia) fMCGeneratorString = "PYTHIA";
  else if(mcgen == kHerwig) fMCGeneratorString = "HERWIG";
  else if(mcgen == kHijing) fMCGeneratorString = "HIJING";
  else
  {
    fMCGeneratorString = "";
    fMCGenerator       = kBoxLike ;
  }
}

//____________________________________________________
/// Set the generator type.
//____________________________________________________
void AliMCAnalysisUtils::SetMCGenerator(TString mcgen)
{  
  fMCGeneratorString = mcgen ;
  
  if     (mcgen == "PYTHIA") fMCGenerator = kPythia;
  else if(mcgen == "HERWIG") fMCGenerator = kHerwig;
  else if(mcgen == "HIJING") fMCGenerator = kHijing;
  else
  {
    fMCGenerator       = kBoxLike;
    fMCGeneratorString = "" ;
  }
}



