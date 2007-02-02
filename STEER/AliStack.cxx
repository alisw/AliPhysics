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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Particles stack class                                                    //
//  Implements the TMCVirtualStack of the Virtual Monte Carlo                //
//  Holds the particles transported during simulation                        //
//  Is used to compare results of reconstruction with simulation             //
//  Author A.Morsch                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

 
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include "AliLog.h"
#include "AliHit.h"
#include "AliModule.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliRunLoader.h"
#include "AliStack.h"

ClassImp(AliStack)

//_______________________________________________________________________
AliStack::AliStack():
  fParticles(0),
  fParticleMap(0),
  fParticleFileMap(0),
  fParticleBuffer(0),
  fCurrentTrack(0),
  fTreeK(0),
  fNtrack(0),
  fNprimary(0),
  fCurrent(-1),
  fCurrentPrimary(-1),
  fHgwmk(0),
  fLoadPoint(0),
  fEventFolderName(AliConfig::GetDefaultEventFolderName())
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliStack::AliStack(Int_t size, const char* evfoldname):
  fParticles(new TClonesArray("TParticle",1000)),
  fParticleMap(new TObjArray(size)),
  fParticleFileMap(0),
  fParticleBuffer(0),
  fCurrentTrack(0),
  fTreeK(0),
  fNtrack(0),
  fNprimary(0),
  fCurrent(-1),
  fCurrentPrimary(-1),
  fHgwmk(0),
  fLoadPoint(0),
  fEventFolderName(evfoldname)
{
  //
  //  Constructor
  //
}

//_______________________________________________________________________
AliStack::AliStack(const AliStack& st):
  TVirtualMCStack(st),
  fParticles(0),
  fParticleMap(0),
  fParticleFileMap(0),
  fParticleBuffer(0),
  fCurrentTrack(0),
  fTreeK(0),
  fNtrack(0),
  fNprimary(0),
  fCurrent(-1),
  fCurrentPrimary(-1),
  fHgwmk(0),
  fLoadPoint(0),
  fEventFolderName()
{
  //
  // Copy constructor
  //
  st.Copy(*this);
}

//_______________________________________________________________________
void AliStack::Copy(TObject&) const
{
  AliFatal("Not implemented!");
}

//_______________________________________________________________________
AliStack::~AliStack()
{
  //
  // Destructor
  //
  
  if (fParticles) {
    fParticles->Delete();
    delete fParticles;
  }
  delete fParticleMap;
  //PH???  delete fTreeK;
}

//
// public methods
//

//_____________________________________________________________________________
void AliStack::PushTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
                        Float_t *vpos, Float_t *polar, Float_t tof,
                        TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{ 
  //
  // Load a track on the stack
  //
  // done     0 if the track has to be transported
  //          1 if not
  // parent   identifier of the parent track. -1 for a primary
  // pdg    particle code
  // pmom     momentum GeV/c
  // vpos     position 
  // polar    polarisation 
  // tof      time of flight in seconds
  // mecha    production mechanism
  // ntr      on output the number of the track stored
  //

  //  const Float_t tlife=0;
  
  //
  // Here we get the static mass
  // For MC is ok, but a more sophisticated method could be necessary
  // if the calculated mass is required
  // also, this method is potentially dangerous if the mass
  // used in the MC is not the same of the PDG database
  //
    TParticlePDG* pmc =  TDatabasePDG::Instance()->GetParticle(pdg);
    if (pmc) {
	Float_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
	Float_t e=TMath::Sqrt(mass*mass+pmom[0]*pmom[0]+
			      pmom[1]*pmom[1]+pmom[2]*pmom[2]);
	
//    printf("Loading  mass %f ene %f No %d ip %d parent %d done %d pos %f %f %f mom %f %f %f kS %d m \n",
//	   mass,e,fNtrack,pdg,parent,done,vpos[0],vpos[1],vpos[2],pmom[0],pmom[1],pmom[2],kS);
  

	PushTrack(done, parent, pdg, pmom[0], pmom[1], pmom[2], e,
		 vpos[0], vpos[1], vpos[2], tof, polar[0], polar[1], polar[2],
		 mech, ntr, weight, is);
    } else {
	AliWarning(Form("Particle type %d not defined in PDG Database !", pdg));
	AliWarning("Particle skipped !");
    }
}

//_____________________________________________________________________________
void AliStack::PushTrack(Int_t done, Int_t parent, Int_t pdg,
  	              Double_t px, Double_t py, Double_t pz, Double_t e,
  		      Double_t vx, Double_t vy, Double_t vz, Double_t tof,
		      Double_t polx, Double_t poly, Double_t polz,
		      TMCProcess mech, Int_t &ntr, Double_t weight, Int_t is)
{ 
  //
  // Load a track on the stack
  //
  // done        0 if the track has to be transported
  //             1 if not
  // parent      identifier of the parent track. -1 for a primary
  // pdg         particle code
  // kS          generation status code
  // px, py, pz  momentum GeV/c
  // vx, vy, vz  position 
  // polar       polarisation 
  // tof         time of flight in seconds
  // mech        production mechanism
  // ntr         on output the number of the track stored
  //    
  // New method interface: 
  // arguments were changed to be in correspondence with TParticle
  // constructor.
  // Note: the energy is not calculated from the static mass but
  // it is passed by argument e.


  const Int_t kFirstDaughter=-1;
  const Int_t kLastDaughter=-1;


  TClonesArray &particles = *fParticles;
  TParticle* particle
    = new(particles[fLoadPoint++]) 
      TParticle(pdg, is, parent, -1, kFirstDaughter, kLastDaughter,
		px, py, pz, e, vx, vy, vz, tof);
   
  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(mech);

  if(!done) particle->SetBit(kDoneBit);

  //  Declare that the daughter information is valid
  particle->SetBit(kDaughtersBit);
  //  Add the particle to the stack
  
  
  fParticleMap->AddAtAndExpand(particle, fNtrack);//CHECK!!

  if(parent>=0) {
    particle = dynamic_cast<TParticle*>(fParticleMap->At(parent));
    if (particle) {
      particle->SetLastDaughter(fNtrack);
      if(particle->GetFirstDaughter()<0) particle->SetFirstDaughter(fNtrack);
    }
    else {
      AliError(Form("Parent %d does not exist",parent));
    }
  } 
  else { 
    //
    // This is a primary track. Set high water mark for this event
    fHgwmk = fNtrack;
    //
    // Set also number if primary tracks
    fNprimary = fHgwmk+1;
    fCurrentPrimary++;
  }
  ntr = fNtrack++;
}

//_____________________________________________________________________________
TParticle*  AliStack::PopNextTrack(Int_t& itrack)
{
  //
  // Returns next track from stack of particles
  //
  

  TParticle* track = GetNextParticle();

  if (track) {
    itrack = fCurrent;
    track->SetBit(kDoneBit);
  }
  else
    itrack = -1;
  
  fCurrentTrack = track;
  return track;
}

//_____________________________________________________________________________
TParticle*  AliStack::PopPrimaryForTracking(Int_t i)
{
  //
  // Returns i-th primary particle if it is flagged to be tracked,
  // 0 otherwise
  //
  
  TParticle* particle = Particle(i);
  
  if (!particle->TestBit(kDoneBit))
    return particle;
  else
    return 0;
}      

//_____________________________________________________________________________
void AliStack::PurifyKine()
{
  //
  // Compress kinematic tree keeping only flagged particles
  // and renaming the particle id's in all the hits
  //

  TObjArray &particles = *fParticleMap;
  int nkeep=fHgwmk+1, parent, i;
  TParticle *part, *father;
  TArrayI map(particles.GetLast()+1);

  // Save in Header total number of tracks before compression
  // If no tracks generated return now
  if(fHgwmk+1 == fNtrack) return;

  // First pass, invalid Daughter information
   for(i=0; i<fNtrack; i++) {
    // Preset map, to be removed later
    if(i<=fHgwmk) map[i]=i ; 
    else {
      map[i] = -99;
      if((part=dynamic_cast<TParticle*>(particles.At(i)))) {
//
//        Check of this track should be kept for physics reasons 
          if (KeepPhysics(part)) KeepTrack(i);
//
          part->ResetBit(kDaughtersBit);
          part->SetFirstDaughter(-1);
          part->SetLastDaughter(-1);
      }
    }
  }
  // Invalid daughter information for the parent of the first particle
  // generated. This may or may not be the current primary according to
  // whether decays have been recorded among the primaries
  part = dynamic_cast<TParticle*>(particles.At(fHgwmk+1));
  particles.At(part->GetFirstMother())->ResetBit(kDaughtersBit);
  // Second pass, build map between old and new numbering
  for(i=fHgwmk+1; i<fNtrack; i++) {
    if(particles.At(i)->TestBit(kKeepBit)) {
      
      // This particle has to be kept
      map[i]=nkeep;
      // If old and new are different, have to move the pointer
      if(i!=nkeep) particles[nkeep]=particles.At(i);
      part = dynamic_cast<TParticle*>(particles.At(nkeep));
      
      // as the parent is always *before*, it must be already
      // in place. This is what we are checking anyway!
      if((parent=part->GetFirstMother())>fHgwmk) 
	if(map[parent]==-99) Fatal("PurifyKine","map[%d] = -99!\n",parent);
	else part->SetFirstMother(map[parent]);

      nkeep++;
    }
  }
  
  // Fix daughters information
  for (i=fHgwmk+1; i<nkeep; i++) {
    part = dynamic_cast<TParticle*>(particles.At(i));
    parent = part->GetFirstMother();
    if(parent>=0) {
      father = dynamic_cast<TParticle*>(particles.At(parent));
      if(father->TestBit(kDaughtersBit)) {
      
	if(i<father->GetFirstDaughter()) father->SetFirstDaughter(i);
	if(i>father->GetLastDaughter())  father->SetLastDaughter(i);
      } else {
	// Initialise daughters info for first pass
	father->SetFirstDaughter(i);
	father->SetLastDaughter(i);
	father->SetBit(kDaughtersBit);
      }
    }
  }

  
  // Now loop on all registered hit lists
  TList* hitLists = gAlice->GetMCApp()->GetHitLists();
  TIter next(hitLists);
  TCollection *hitList;
  
  while((hitList = dynamic_cast<TCollection*>(next()))) {
    TIter nexthit(hitList);
    AliHit *hit;
    
    while((hit = dynamic_cast<AliHit*>(nexthit()))) {
	hit->SetTrack(map[hit->GetTrack()]);
    }
  }

  // 
  // This for detectors which have a special mapping mechanism
  // for hits, such as TPC and TRD
  //
  
   TObjArray* modules = gAlice->Modules();
   TIter nextmod(modules);
   AliModule *detector;
   while((detector = dynamic_cast<AliModule*>(nextmod()))) {
     detector->RemapTrackHitIDs(map.GetArray());
     detector->RemapTrackReferencesIDs(map.GetArray());
   }
   //
   gAlice->GetMCApp()->RemapTrackReferencesIDs(map.GetArray());

   // Now the output bit, from fHgwmk to nkeep we write everything and we erase
   if(nkeep>fParticleFileMap.GetSize()) fParticleFileMap.Set(Int_t (nkeep*1.5));

   for (i=fHgwmk+1; i<nkeep; ++i) {
     fParticleBuffer = dynamic_cast<TParticle*>(particles.At(i));
     fParticleFileMap[i]=static_cast<Int_t>(TreeK()->GetEntries());
     TreeK()->Fill();
     particles[i]=fParticleBuffer=0;
    }

   for (i=nkeep; i<fNtrack; ++i) particles[i]=0;

   Int_t toshrink = fNtrack-fHgwmk-1;
   fLoadPoint-=toshrink;

   
   for(i=fLoadPoint; i<fLoadPoint+toshrink; ++i) fParticles->RemoveAt(i);

   fNtrack=nkeep;
   fHgwmk=nkeep-1;
}

void AliStack::ReorderKine()
{
//
// In some transport code children might not come in a continuous sequence.
// In this case the stack  has  to  be reordered in order to establish the 
// mother daughter relation using index ranges.
//    
  if(fHgwmk+1 == fNtrack) return;

  //
  // Howmany secondaries have been produced ?
  Int_t nNew = fNtrack - fHgwmk - 1;
    
  if (nNew > 0) {
      Int_t i, j;
      TObjArray &particles = *fParticleMap;
      TArrayI map1(nNew);
      //
      // Copy pointers to temporary array
      TParticle** tmp = new TParticle*[nNew];
      
      for (i = 0; i < nNew; i++) {
	  if (particles.At(fHgwmk + 1 + i)) {
	      tmp[i] = (TParticle*) (particles.At(fHgwmk + 1 + i));
	  } else {
	      tmp[i] = 0x0;
	  }
	  map1[i] = -99;
      }
  
      
      //
      // Reset  LoadPoint 
      // 
      fLoadPoint = fHgwmk + 1;
      //
      // Re-Push particles into stack 
      // The outer loop is over parents, the inner over children.
      // -1 refers to the primary particle
      //
      for (i = -1; i < nNew-1; i++) {
	  Int_t ipa;
	  TParticle* parP;
	  if (i == -1) {
	      ipa  = tmp[0]->GetFirstMother();
	      parP =dynamic_cast<TParticle*>(particles.At(ipa));
	  } else {
	      ipa = (fHgwmk + 1 + i);
              // Skip deleted particles
	      if (!tmp[i])                          continue;
              // Skip particles without children
	      if (tmp[i]->GetFirstDaughter() == -1) continue;
	      parP = tmp[i];
	  }
          // Reset daughter information

	  Int_t idaumin = parP->GetFirstDaughter() - fHgwmk - 1;
	  Int_t idaumax = parP->GetLastDaughter()  - fHgwmk - 1;
	  parP->SetFirstDaughter(-1);
	  parP->SetLastDaughter(-1);
	  for (j = idaumin; j <= idaumax; j++) {
              // Skip deleted particles
	      if (!tmp[j])        continue;
              // Skip particles already handled
	      if (map1[j] != -99) continue;
	      Int_t jpa = tmp[j]->GetFirstMother();
              // Check if daughter of current parent
	      if (jpa == ipa) {
		  particles[fLoadPoint] = tmp[j];
		  // Re-establish daughter information
		  parP->SetLastDaughter(fLoadPoint);
		  if (parP->GetFirstDaughter() == -1) parP->SetFirstDaughter(fLoadPoint);
		  // Set Mother information
		  if (i != -1) {
		      tmp[j]->SetFirstMother(map1[i]);
		  } 
		  // Build the map
		  map1[j] = fLoadPoint;
		  // Increase load point
		  fLoadPoint++;
	      }
	  } // children
      } // parents

      delete[] tmp;

      //
      // Build map for remapping of hits
      // 
      TArrayI map(fNtrack);
      for (i = 0; i < fNtrack; i ++) {
	  if (i <= fHgwmk) {
	      map[i] = i;
	  } else{
	      map[i] = map1[i - fHgwmk -1];
	  }
      }
      
      // Now loop on all registered hit lists
      
      TList* hitLists = gAlice->GetMCApp()->GetHitLists();
      TIter next(hitLists);
      TCollection *hitList;
      
      while((hitList = dynamic_cast<TCollection*>(next()))) {
	  TIter nexthit(hitList);
	  AliHit *hit;
	  while((hit = dynamic_cast<AliHit*>(nexthit()))) {
	      hit->SetTrack(map[hit->GetTrack()]);
	  }
      }
  
  // 
  // This for detectors which have a special mapping mechanism
  // for hits, such as TPC and TRD
  //
  
      TObjArray* modules = gAlice->Modules();
      TIter nextmod(modules);
      AliModule *detector;
      while((detector = dynamic_cast<AliModule*>(nextmod()))) {
	  detector->RemapTrackHitIDs(map.GetArray());
	  detector->RemapTrackReferencesIDs(map.GetArray());
      }
      gAlice->GetMCApp()->RemapTrackReferencesIDs(map.GetArray());
  } // new particles poduced
}

Bool_t AliStack::KeepPhysics(TParticle* part)
{
    //
    // Some particles have to kept on the stack for reasons motivated
    // by physics analysis. Decision is put here.
    //
    Bool_t keep = kFALSE;
    //
    // Keep first-generation daughter from primaries with heavy flavor 
    //
    Int_t parent = part->GetFirstMother();
    if (parent >= 0 && parent <= fHgwmk) {
	TParticle* father = dynamic_cast<TParticle*>(Particles()->At(parent));
	Int_t kf = father->GetPdgCode();
	kf = TMath::Abs(kf);
	Int_t kfl = kf;
	// meson ?
	if  (kfl > 10) kfl/=100;
	// baryon
	if (kfl > 10) kfl/=10;
	if (kfl > 10) kfl/=10;
	if (kfl >= 4) {
	    keep = kTRUE;
	}
    }
    return keep;
}

//_____________________________________________________________________________
void AliStack::FinishEvent()
{
//
// Write out the kinematics that was not yet filled
//
  
// Update event header


  if (!TreeK()) {
//    Fatal("FinishEvent", "No kinematics tree is defined.");
//    Don't panic this is a probably a lego run
      return;
  }  
  
  CleanParents();
  if(TreeK()->GetEntries() ==0) {
    // set the fParticleFileMap size for the first time
    fParticleFileMap.Set(fHgwmk+1);
  }

  Bool_t allFilled = kFALSE;
  TObject *part;
  for(Int_t i=0; i<fHgwmk+1; ++i) 
    if((part=fParticleMap->At(i))) {
      fParticleBuffer = dynamic_cast<TParticle*>(part);
      fParticleFileMap[i]= static_cast<Int_t>(TreeK()->GetEntries());
      TreeK()->Fill();
      //PH      (*fParticleMap)[i]=fParticleBuffer=0;      
      fParticleBuffer=0;      
      fParticleMap->AddAt(0,i);      
      
      // When all primaries were filled no particle!=0
      // should be left => to be removed later.
      if (allFilled) AliWarning(Form("Why != 0 part # %d?\n",i));
    }
    else 
     {
      // // printf("Why = 0 part # %d?\n",i); => We know.
      // break;
      // we don't break now in order to be sure there is no
      // particle !=0 left.
      // To be removed later and replaced with break.
       if(!allFilled) allFilled = kTRUE;
     } 
} 
//_____________________________________________________________________________

void AliStack::FlagTrack(Int_t track)
{
  //
  // Flags a track and all its family tree to be kept
  //
  
  TParticle *particle;

  Int_t curr=track;
  while(1) {
    particle=dynamic_cast<TParticle*>(fParticleMap->At(curr));
    
    // If the particle is flagged the three from here upward is saved already
    if(particle->TestBit(kKeepBit)) return;
    
    // Save this particle
    particle->SetBit(kKeepBit);
    
    // Move to father if any
    if((curr=particle->GetFirstMother())==-1) return;
  }
}
 
//_____________________________________________________________________________
void AliStack::KeepTrack(Int_t track)
{ 
  //
  // Flags a track to be kept
  //
  
  fParticleMap->At(track)->SetBit(kKeepBit);
}

//_____________________________________________________________________________
void  AliStack::Reset(Int_t size) 
{
  //
  // Resets stack
  //
  
  fNtrack=0;
  fNprimary=0;
  fHgwmk=0;
  fLoadPoint=0;
  fCurrent = -1;
  fTreeK = 0x0;
  ResetArrays(size);
}

//_____________________________________________________________________________
void  AliStack::ResetArrays(Int_t size) 
{
  //
  // Resets stack arrays
  //

  if (fParticles) 
    fParticles->Clear();
  else
    fParticles = new TClonesArray("TParticle",1000);
  if (fParticleMap) {
    fParticleMap->Clear();
    if (size>0) fParticleMap->Expand(size);}
  else
    fParticleMap = new TObjArray(size);
}

//_____________________________________________________________________________
void AliStack::SetHighWaterMark(Int_t)
{
  //
  // Set high water mark for last track in event
  //

    
  fHgwmk = fNtrack-1;
  fCurrentPrimary=fHgwmk;
  // Set also number of primary tracks
  fNprimary = fHgwmk+1;
  fNtrack   = fHgwmk+1;      
}

//_____________________________________________________________________________
TParticle* AliStack::Particle(Int_t i)
{
  //
  // Return particle with specified ID
  
  //PH  if(!(*fParticleMap)[i]) {
  if(!fParticleMap->At(i)) {
    Int_t nentries = fParticles->GetEntriesFast();
    // algorithmic way of getting entry index
    // (primary particles are filled after secondaries)
    Int_t entry = TreeKEntry(i);
    // check whether algorithmic way and 
    // and the fParticleFileMap[i] give the same;
    // give the fatal error if not
    if (entry != fParticleFileMap[i]) {
      AliFatal(Form(
        "!! The algorithmic way and map are different: !!\n entry: %d map: %d",
	entry, fParticleFileMap[i])); 
    } 
      
    TreeK()->GetEntry(entry);
    new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
    fParticleMap->AddAt((*fParticles)[nentries],i);
  }
  //PH  return dynamic_cast<TParticle *>((*fParticleMap)[i]);
  return dynamic_cast<TParticle*>(fParticleMap->At(i));
}

//_____________________________________________________________________________
TParticle* AliStack::ParticleFromTreeK(Int_t id) const
{
// 
// return pointer to TParticle with label id
//
  Int_t entry;
  if ((entry = TreeKEntry(id)) < 0) return 0;
  if (fTreeK->GetEntry(entry)<=0) return 0;
  return fParticleBuffer;
}

//_____________________________________________________________________________
Int_t AliStack::TreeKEntry(Int_t id) const 
{
//
// return entry number in the TreeK for particle with label id
// return negative number if label>fNtrack
//
  Int_t entry;
  if (id<fNprimary)
    entry = id+fNtrack-fNprimary;
  else 
    entry = id-fNprimary;
  return entry;
}

//_____________________________________________________________________________
Int_t AliStack::GetCurrentParentTrackNumber() const
{
  //
  // Return number of the parent of the current track
  //
  
  TParticle* current = (TParticle*)fParticleMap->At(fCurrent);

  if (current) 
    return current->GetFirstMother();
  else {
    AliWarning("Current track not found in the stack");
    return -1;
  }  
}
 
//_____________________________________________________________________________
Int_t AliStack::GetPrimary(Int_t id)
{
  //
  // Return number of primary that has generated track
  //
  
  int current, parent;
  //
  parent=id;
  while (1) {
    current=parent;
    parent=Particle(current)->GetFirstMother();
    if(parent<0) return current;
  }
}
 
//_____________________________________________________________________________
void AliStack::DumpPart (Int_t i) const
{
  //
  // Dumps particle i in the stack
  //
  
  //PH  dynamic_cast<TParticle*>((*fParticleMap)[i])->Print();
  dynamic_cast<TParticle*>(fParticleMap->At(i))->Print();
}

//_____________________________________________________________________________
void AliStack::DumpPStack ()
{
  //
  // Dumps the particle stack
  //

  Int_t i;

  printf("\n\n=======================================================================\n");
  for (i=0;i<fNtrack;i++) 
    {
      TParticle* particle = Particle(i);
      if (particle) {
        printf("-> %d ",i); particle->Print();
        printf("--------------------------------------------------------------\n");
      }
      else 
        Warning("DumpPStack", "No particle with id %d.", i); 
    }	 

  printf("\n=======================================================================\n\n");
  
  // print  particle file map
  printf("\nParticle file map: \n");
  for (i=0; i<fNtrack; i++) 
      printf("   %d th entry: %d \n",i,fParticleFileMap[i]);
}


//_____________________________________________________________________________
void AliStack::DumpLoadedStack() const
{
  //
  // Dumps the particle in the stack
  // that are loaded in memory.
  //

  TObjArray &particles = *fParticleMap;
  printf(
	 "\n\n=======================================================================\n");
  for (Int_t i=0;i<fNtrack;i++) 
    {
      TParticle* particle = dynamic_cast<TParticle*>(particles[i]);
      if (particle) {
        printf("-> %d ",i); particle->Print();
        printf("--------------------------------------------------------------\n");
      }
      else { 	
        printf("-> %d  Particle not loaded.\n",i);
        printf("--------------------------------------------------------------\n");
      }	
    }
  printf(
	 "\n=======================================================================\n\n");
}

//
// protected methods
//

//_____________________________________________________________________________
void AliStack::CleanParents()
{
  //
  // Clean particles stack
  // Set parent/daughter relations
  //
  
  TObjArray &particles = *fParticleMap;
  TParticle *part;
  int i;
  for(i=0; i<fHgwmk+1; i++) {
    part = dynamic_cast<TParticle*>(particles.At(i));
    if(part) if(!part->TestBit(kDaughtersBit)) {
      part->SetFirstDaughter(-1);
      part->SetLastDaughter(-1);
    }
  }
}

//_____________________________________________________________________________
TParticle* AliStack::GetNextParticle()
{
  //
  // Return next particle from stack of particles
  //
  
  TParticle* particle = 0;
  
  // search secondaries
  //for(Int_t i=fNtrack-1; i>=0; i--) {
  for(Int_t i=fNtrack-1; i>fHgwmk; i--) {
      particle = dynamic_cast<TParticle*>(fParticleMap->At(i));
      if ((particle) && (!particle->TestBit(kDoneBit))) {
	  fCurrent=i;    
	  return particle;
      }   
  }    

  // take next primary if all secondaries were done
  while (fCurrentPrimary>=0) {
      fCurrent = fCurrentPrimary;    
      particle = dynamic_cast<TParticle*>(fParticleMap->At(fCurrentPrimary--));
      if ((particle) && (!particle->TestBit(kDoneBit))) {
	  return particle;
      } 
  }
  
  // nothing to be tracked
  fCurrent = -1;
 
  
  return particle;  
}
//__________________________________________________________________________________________

TTree* AliStack::TreeK()
{
//returns TreeK
  if (fTreeK)
   {
     return fTreeK;
   }
  else
   {
     AliRunLoader *rl = AliRunLoader::GetRunLoader(fEventFolderName);
     if (rl == 0x0)
      {
        AliFatal(Form("Can not get RunLoader from event folder named %s",fEventFolderName.Data()));
        return 0x0;//pro forma
      }
    fTreeK = rl->TreeK();
    if ( fTreeK )
     {
      ConnectTree();
     }
    else
     {
      //don't panic - could be Lego
      AliWarning(Form("Can not get TreeK from RL. Ev. Folder is %s",fEventFolderName.Data()));
     }
   }
  return fTreeK;//never reached
}
//__________________________________________________________________________________________

void AliStack::ConnectTree()
{
//
//  Creates branch for writing particles
//
  AliDebug(1, "Connecting TreeK");
  if (fTreeK == 0x0)
   {
    if (TreeK() == 0x0)
     {
      AliFatal("Parameter is NULL");//we don't like such a jokes
      return;
     }
    return;//in this case TreeK() calls back this method (ConnectTree) 
           //tree after setting fTreeK, the rest was already executed
           //it is safe to return now
   }

 //  Create a branch for particles   
  
  AliDebug(2, Form("Tree name is %s",fTreeK->GetName()));
   
  if (fTreeK->GetDirectory())
   {
     AliDebug(2, Form("and dir is %s",fTreeK->GetDirectory()->GetName()));
   }    
  else
    AliWarning("DIR IS NOT SET !!!");
  
  TBranch *branch=fTreeK->GetBranch(AliRunLoader::GetKineBranchName());
  if(branch == 0x0)
   {
    branch = fTreeK->Branch(AliRunLoader::GetKineBranchName(), "TParticle", &fParticleBuffer, 4000);
    AliDebug(2, "Creating Branch in Tree");
   }  
  else
   {
    AliDebug(2, "Branch Found in Tree");
    branch->SetAddress(&fParticleBuffer);
   }
  if (branch->GetDirectory())
   {
    AliDebug(1, Form("Branch Dir Name is %s",branch->GetDirectory()->GetName()));
   } 
  else
    AliWarning("Branch Dir is NOT SET");
}
//__________________________________________________________________________________________


void AliStack::BeginEvent()
{
// start a new event
 Reset();
}

//_____________________________________________________________________________
void AliStack::FinishRun()
{
// Clean TreeK information
}
//_____________________________________________________________________________

Bool_t AliStack::GetEvent()
{
//
// Get new event from TreeK

    // Reset/Create the particle stack
    fTreeK = 0x0;

    if (TreeK() == 0x0) //forces connecting
     {
      AliError("cannot find Kine Tree for current event");
      return kFALSE;
     }
      
    Int_t size = (Int_t)TreeK()->GetEntries();
    ResetArrays(size);
    return kTRUE;
}
//_____________________________________________________________________________

void AliStack::SetEventFolderName(const char* foldname)
{
 //Sets event folder name
 fEventFolderName = foldname;
}

Bool_t AliStack::IsStable(Int_t pdg) const
{
//
// Decide whether particle (pdg) is stable
//

    const Int_t kNstable = 14;
    Int_t i;

    Int_t pdgStable[kNstable] = {
	kGamma,             // Photon
	kElectron,          // Electron
	kMuonPlus,          // Muon 
	kPiPlus,            // Pion
	kKPlus,             // Kaon
	kProton,            // Proton 
	kNeutron,           // Neutron
	kLambda0,           // Lambda_0
	kSigmaMinus,        // Sigma Minus
	kSigma0,            // Sigma_0
	kSigmaPlus,         // Sigma Plus
	3312,               // Xsi Minus 
	3322,               // Xsi 
	3334,               // Omega
    };
    
    Bool_t isStable = kFALSE;
    for (i = 0; i < kNstable; i++) {
	if (pdg == TMath::Abs(pdgStable[i])) {
	    isStable = kTRUE;
	    break;
	}
    }

    return isStable;
}

Bool_t AliStack::IsPhysicalPrimary(Int_t index)
{
    //
    // Test if a particle is a physical primary according to the following definition:
    // Particles produced in the collision including products of strong and
    // electromagnetic decay and excluding feed-down from weak decays of strange
    // particles.
    //
    TParticle* p = Particle(index);
    Int_t ist = p->GetStatusCode();
    
    //
    // Initial state particle
    if (ist > 20) return kFALSE;
    
    Int_t pdg = TMath::Abs(p->GetPdgCode());
    
    if (!IsStable(pdg)) return kFALSE;
    
    if (index < GetNprimary()) {
//
// Particle produced by generator
	return kTRUE;
    } else {
//
// Particle produced during transport
//
// Check if this is a heavy flavor decay product
	Int_t imo =  p->GetFirstMother();
	TParticle* pm  = Particle(imo);
	Int_t mpdg = TMath::Abs(pm->GetPdgCode());
	Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));
	//
	// Light hadron
	if (mfl < 4) return kFALSE;
	
	//
	// Heavy flavor hadron produced by generator
	if (imo <  GetNprimary()) {
	    return kTRUE;
	}
	
	// To be sure that heavy flavor has not been produced in a secondary interaction
	// Loop back to the generated mother
	while (imo >=  GetNprimary()) {
	    imo = p->GetFirstMother();
	    pm  =  Particle(imo);
	}
	mpdg = TMath::Abs(pm->GetPdgCode());
	mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

	if (mfl < 4) {
	    return kFALSE;
	} else {
	    return kTRUE;
	} 
    } // produced by generator ?
} 

