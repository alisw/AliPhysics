// module identifier line...
// Author: Brian Thorsbro, 24/6-2014

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <set>
#include <time.h>

#include "THepMCParser.h"
#include "TObject.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "HepMC/IO_GenEvent.h"


using namespace std;


THepMCParser::THepMCParser(const char * infile) : fTree(0)
{
   HepMC::IO_BaseClass * events = new HepMC::IO_GenEvent(infile, std::ios::in);
   init(events);
   delete events;
   events = 0; // nullptr
}
THepMCParser::THepMCParser(HepMC::IO_BaseClass * events) : fTree(0)
{
   init(events);
}
void THepMCParser::init(HepMC::IO_BaseClass * events)
{
   int particlecount = 0;
   fTree = new TTree("treeEPOS","Tree EPOS");
   TClonesArray * array = new TClonesArray("TParticle");
   // array->BypassStreamer();
   fTree->Branch("Particles",&array); // more flags?
   THepMCParser::HeavyIonHeader_t heavyIonHeader;
   fTree->Branch("HeavyIonInfo", &heavyIonHeader, THepMCParser::fgHeavyIonHeaderBranchString);
   THepMCParser::PdfHeader_t pdfHeader;
   fTree->Branch("PdfInfo", &pdfHeader, THepMCParser::fgPdfHeaderBranchString);
   HepMC::GenEvent* evt =  0; // nullptr
   while ((evt = events->read_next_event())) {
      string errMsg1 = ParseGenEvent2TCloneArray(evt,array);
      string errMsg2 = ParseGenEvent2HeaderStructs(evt,heavyIonHeader,pdfHeader);
      if (errMsg1.length() == 0 && errMsg2.length() == 0) {
         fTree->Fill();
      } else {
         if (errMsg1.length() != 0) cerr << errMsg1 << endl;
         if (errMsg2.length() != 0) cerr << errMsg2 << endl;
      }
      particlecount += array->Capacity();
   }
//   array->Clear();
//   delete array;
//   array = 0; // nullptr
   cout << " parsed " << particlecount << " particles" << endl;
   if(array->Capacity() != array->GetEntries()) {
     cerr << ("Not all particles processed");
   }
}


TTree * THepMCParser::GetTTree()
{
   return fTree;
}
void THepMCParser::WriteTTreeToFile(const char *outfile)
{
   TFile * f = new TFile(outfile, "recreate");
   fTree->Write();
   delete f;
   f = 0; // nullptr
}



// Based on a validator written by Peter Hristov, CERN
bool THepMCParser::IsValidMotherDaughtersConsitency(bool useStdErr, bool requireSecondMotherBeforeDaughters)
{
   bool valid = true;
   TClonesArray * array = new TClonesArray("TParticle");
   TBranch* branch = fTree->GetBranch("Particles");
   branch->SetAddress(&array);
   Int_t count = branch->GetEntries();
   for (Int_t idx=0; idx<count; ++idx) {
      array->Clear();
      branch->GetEntry(idx); // "fill" the array
      Int_t nkeep = array->GetEntriesFast();
      for (Int_t i=0; i<nkeep; i++) {
         TParticle * part = (TParticle*)array->AddrAt(i);
         Int_t mum1 = part->GetFirstMother();
         Int_t mum2 = part->GetSecondMother();
         Int_t fd = part->GetFirstDaughter();
         Int_t ld = part->GetLastDaughter();
         if (mum1>-1 && i<mum1) {
            valid = false;
            if (useStdErr) cerr << "Problem: first_mother(" << mum1 << ") after daughter(" << i << ")" << endl;
         }
         if (mum2>-1 && i<mum2 && requireSecondMotherBeforeDaughters) {
            valid = false;
            if (useStdErr) cerr << "Problem: second_mother(" << mum2 << ") after daughter(" << i << ")" << endl;
         }
         if (fd > ld ) {
            valid = false;
            if (useStdErr) cerr << "Problem: first_daughter(" << fd << ") > last_daughter(" << ld << ")" << endl;
         }
         for (Int_t id=TMath::Max(fd,0); id<=ld; id++) {
            TParticle * daughter = (TParticle*)array->AddrAt(id);
            if (daughter->GetFirstMother() != i && daughter->GetSecondMother() != i) {
               valid = false;
               if (useStdErr) cerr << "Problem: mother("<< i << ") not registered as mother for either first_daughter("
                     << daughter->GetFirstMother() << ") or second_daughter("
                     << daughter->GetSecondMother() << ")" << endl;
            }
         }
      }
   }
   delete array;
   array = 0;
   return valid;
}

bool THepMCParser::IsValidParticleInvariantMass(bool useStdErr, bool includeStatusCode2Particles)
{
   bool valid = true;
   TClonesArray *array = new TClonesArray("TParticle");
   TBranch* branch = fTree->GetBranch("Particles");
   branch->SetAddress(&array);
   Int_t count = branch->GetEntries();
   for (Int_t idx=0; idx<count; ++idx) {
      array->Clear();
      branch->GetEntry(idx);
      Int_t nkeep = array->GetEntries();
      for (Int_t i=0; i<nkeep; i++) {
         TParticle * parton = (TParticle*)array->AddrAt(i);
         if (parton->GetStatusCode()==2 && !includeStatusCode2Particles) {
            continue;
         }
         TLorentzVector v;
         parton->Momentum(v);
         Double_t m1 = v.M(); // invariant mass from the particle in the HepMC event
         TParticlePDG *dbParton = parton->GetPDG();
         if (!dbParton) {
            if (useStdErr) cerr << "Warning: could not look up particle with PDG Code: " << parton->GetPdgCode() << endl << endl;
            continue;
         }
         Double_t m2 = dbParton->Mass();
         bool checkok;
         if (m2 == 0) {
            checkok = abs(m1) < 0.0001; // no such thing as negative mass...
         } else {
            checkok = abs(1 - m1/m2) < 0.01;
         }
         if (!checkok && useStdErr) {
            cerr << "Problem: " << GetParticleName(parton) << " HepMC:" << m1 << endl;
            cerr << ListReactionChain(array,i);
            cerr << endl;
         }
         if (!checkok)
            valid = false;
      }
   }
   delete array;
   array = 0;
   return valid;
}

bool THepMCParser::IsValidVertexInvariantMass(bool useStdErr, bool includeStatusCode2Particles)
{
   bool valid = true;
   TClonesArray * array = new TClonesArray("TParticle");
   TBranch* branch = fTree->GetBranch("Particles");
   branch->SetAddress(&array);
   Int_t count = branch->GetEntries();
   for (Int_t idx=0; idx<count; ++idx) {
      array->Clear();
      branch->GetEntry(idx); // "fill" the array
      TLorentzVector v_st1;
      TLorentzVector v_st4;
      Int_t nkeep = array->GetEntriesFast();
      for (Int_t i=0; i<nkeep; i++) {
         TParticle * parton = (TParticle*)array->AddrAt(i);
         TLorentzVector v_in;
         parton->Momentum(v_in);
         if (parton->GetStatusCode()==4) {
            v_st4 += v_in;
         } else if (parton->GetStatusCode()==1) {
            v_st1 += v_in;
         }
         if (!includeStatusCode2Particles) { // only check beam particles vs final particles
            continue;
         }
         Int_t fd = parton->GetFirstDaughter();
         Int_t ld = parton->GetLastDaughter();
         if (fd == -1) continue; // no daughters, continue loop
         Int_t mother2 = -1;
         TLorentzVector v_out;
         bool oneok = false; bool allok = true; bool agreemother2 = true;
         ostringstream daughterpdg;
         ostringstream motherpdg;
         for (Int_t id=TMath::Max(fd,0); id<=ld; id++) {
            TParticle * daughter = (TParticle*)array->AddrAt(id);
            if (fd==id) {
               daughter->Momentum(v_out);
               mother2 = daughter->GetSecondMother();
            } else {
               TLorentzVector d;
               daughter->Momentum(d);
               v_out += d;
               if (daughter->GetSecondMother() != mother2) agreemother2 = false;
            }
            if (daughter->GetFirstMother() == i) {
               oneok = true;
            } else {
               allok = false;
            }
            daughterpdg << " " << daughter->GetPdgCode();
         }
         motherpdg << " " << parton->GetPdgCode();
         if (mother2 > -1 && agreemother2) {
            TParticle * m2 = (TParticle*)array->AddrAt(mother2);
            TLorentzVector m2v;
            m2->Momentum(m2v);
            v_in += m2v;
            motherpdg << " " << m2->GetPdgCode();
         }
         if (oneok && allok && agreemother2) {
            bool checkok = abs(1 - v_in.M()/v_out.M()) < 0.1;
            if (!checkok) valid=false;
            if (!checkok && useStdErr) {
               //             cerr << "Problem: " << i << "[" << fd << "," << ld << "] PDG:" << motherpdg.str() << " ->" << daughterpdg.str() << "  Inv.mass: " << v_in.M() << " -> " << v_out.M() << endl;
               cerr << ListReactionChain(array,i);
               cerr << endl;
               //             cerr << v_in.Px() << " " << v_in.Py() << " " << v_in.Pz() << " " << v_in.E() << endl;
            } else if (useStdErr) {
               //             cerr << "OK: " << i << "[" << fd << "," << ld << "] PDG:" << motherpdg.str() << " ->" << daughterpdg.str() << "  Inv.mass: " << v_in.M() << " -> " << v_out.M() << endl;
            }
         }
      }
      bool checkok = abs(1 - v_st4.M()/v_st1.M()) < 0.001;
      if (!checkok) valid=false;
      if (!checkok && useStdErr) {
         cerr << " BEAM PARTICLES -> FINAL PARTICLES " << endl;
         cerr << " status 4 (" << v_st4.M() << ") -> status 1 (" << v_st1.M() << ")" << endl << endl;
      }
   }
   delete array;
   array = 0;
   return valid;
}

string THepMCParser::GetParticleName(TParticle * thisPart)
{
   TParticlePDG *dbPart = thisPart->GetPDG();
   ostringstream name;
   if (dbPart) {
      name << dbPart->GetName() << "(" << dbPart->PdgCode() << "," << dbPart->Mass() << ")";
   } else {
      name << thisPart->GetPdgCode() << "(NoDBinfo)";
   }
   return name.str();
}

string THepMCParser::ListReactionChain(TClonesArray * particles, Int_t particleId)
{
   ostringstream output;

   TParticle * part = (TParticle*)particles->AddrAt(particleId);
   Int_t m1id = part->GetFirstMother();
   Int_t m2id = part->GetSecondMother();
   if (m1id > 1) { // ignore the initial collision with beam particles
      ostringstream inStr;
      ostringstream outStr;
      TParticle * m1 = (TParticle*)particles->AddrAt(m1id);
      TLorentzVector v_in;
      m1->Momentum(v_in);
      inStr << GetParticleName(m1) << "[" << v_in.M() << "]";
      if (m2id > 1) {
         TParticle * m2 = (TParticle*)particles->AddrAt(m2id);
         TLorentzVector v_m2;
         m2->Momentum(v_m2);
         v_in += v_m2;
         inStr << " " << GetParticleName(m2) << "[" << v_m2.M() << "]";
      }
      Int_t fd = m1->GetFirstDaughter();
      Int_t ld = m1->GetLastDaughter();
      TLorentzVector v_out;
      part->Momentum(v_out);
      outStr << GetParticleName(part) << "[" << v_out.M() << "]";
      for (Int_t i=fd; i<=ld; ++i) {
         if (i!=particleId) {
            TParticle * d = (TParticle*)particles->AddrAt(i);
            TLorentzVector v_d;
            d->Momentum(v_d);
            v_out += v_d;
            outStr << " " << GetParticleName(d) << "[" << v_d.M() << "]";
         }
      }
      output << "Parent reaction, inv mass: " << v_in.M() << " -> " << v_out.M() << endl;
      output << " - partons: " << inStr.str() << " -> " << outStr.str() << endl;
   }
   Int_t fd = part->GetFirstDaughter();
   Int_t ld = part->GetLastDaughter();
   if (fd > -1) {
      ostringstream inStr;
      ostringstream outStr;
      TLorentzVector v_in;
      part->Momentum(v_in);
      inStr << GetParticleName(part) << "[" << v_in.M() << "]";

      TParticle * f = (TParticle*)particles->AddrAt(fd);
      m2id = f->GetSecondMother();
      if (m2id == particleId) {
         m2id = f->GetFirstMother();
      }
      if (m2id > -1) {
         TParticle * m2 = (TParticle*)particles->AddrAt(m2id);
         TLorentzVector v_m2;
         m2->Momentum(v_m2);
         v_in += v_m2;
         inStr << " " << GetParticleName(m2) << "[" << v_m2.M() << "]";
      }
      TLorentzVector v_out;
      f->Momentum(v_out);
      outStr << GetParticleName(f) << "[" << v_out.M() << "]";
      for (Int_t i=fd+1; i<=ld; ++i) {
         TParticle * d = (TParticle*)particles->AddrAt(i);
         TLorentzVector v_d;
         d->Momentum(v_d);
         v_out += v_d;
         outStr << " " << GetParticleName(d) << "[" << v_d.M() << "]";
      }
      output << "Child reaction, inv mass: " << v_in.M() << " -> " << v_out.M() << endl;
      output << " - partons: " << inStr.str() << " -> " << outStr.str() << endl;
   } else {
      output << "Child reaction" << endl << " - none" << endl;
   }

   return output.str();
}


string THepMCParser::ParseGenEvent2TCloneArray(HepMC::GenEvent * genEvent, TClonesArray * array, string momUnit, string lenUnit, bool requireSecondMotherBeforeDaughters)
{
   ostringstream errMsgStream;
   if (requireSecondMotherBeforeDaughters) {
      errMsgStream <<  " WARNING: requireSecondMotherBeforeDaughters not fully implemented yet!\n";
   }
   genEvent->use_units(momUnit, lenUnit);
   array->Clear();
   map<int,Int_t> partonMemory; // unordered_map in c++11 - but probably not much performance gain from that: log(n) vs log(1) where constant can be high
   Bool_t beamParticlesToTheSameVertex = kTRUE; // This is false if the beam particles do not end up in the same vertex. This happens for some HepMC files produced by agile-runmc

   
   // Check event with HepMC's internal validation algorithm
   if (!genEvent->is_valid()) {
      errMsgStream << "Error with event id: " << genEvent->event_number() << ", event is not valid!\n";
      return errMsgStream.str();
   }

   // Pull out the beam particles from the event
   const pair<HepMC::GenParticle *,HepMC::GenParticle *> beamparts = genEvent->beam_particles();

   // Four sanity checks:
   // - Beam particles exists and are not the same
   // - Both beam particles should have no production vertices, they come from the beams
   // - Both beam particles should have defined end vertices, as they both should contribute
   // - Both beam particles should have the exact same end vertex
   if (!beamparts.first || !beamparts.second || beamparts.first->barcode()==beamparts.second->barcode()) {
      errMsgStream << "Error with event id: " << genEvent->event_number() << ", beam particles doesn't exists or are the same\n";
      return errMsgStream.str();
   }
   if (beamparts.first->production_vertex() || beamparts.second->production_vertex()) {
      errMsgStream << "Error with event id: " << genEvent->event_number() << ", beam particles have production vertex/vertices...\n";
      return errMsgStream.str();
   }
   if (!beamparts.first->end_vertex() || !beamparts.second->end_vertex()) {
      errMsgStream << "Error with event id: " << genEvent->event_number() << ", beam particles have undefined end vertex/vertices...\n";
      return errMsgStream.str();
   }
   if (beamparts.first->end_vertex() != beamparts.second->end_vertex()) {
     
      errMsgStream << "Warning with event id: " << genEvent->event_number() << ", beam particles do not collide in the same end vertex.\n";// This was downgraded to a warning, because I noticed in Pythia 6 HepMC (generated via agile-runmc) files the barcode of the vertex was different, but the 4-vector of the vertices was the same.
      beamParticlesToTheSameVertex = kFALSE;
      if (beamparts.first->end_vertex()->position() != beamparts.second->end_vertex()->position()){
        errMsgStream << "Error with event id: " << genEvent->event_number() << ", beam particles do not collide in the same end vertex and the two vertices have different four vectors.\n"; 
        return errMsgStream.str();
      }
   }

   // Set the array to hold the number of particles in the event
   Int_t nParticlesInHepMC = genEvent->particles_size(); // FIXME
   array->Expand(genEvent->particles_size());

   // Create a TParticle for each beam particle
   new((*array)[0]) TParticle(
         beamparts.first->pdg_id(),
         beamparts.first->status(), // check if status has the same meaning
         -1, // no mother1
         -1, // no mother2
         -1, // first daughter not known yet
         -1, // last daughter not known yet
         beamparts.first->momentum().px(),
         beamparts.first->momentum().py(),
         beamparts.first->momentum().pz(),
         beamparts.first->momentum().e(),
         0, // no production vertex, so zero?
         0,
         0,
         0
   );
   partonMemory[beamparts.first->barcode()] = 0;
   new((*array)[1]) TParticle(
         beamparts.second->pdg_id(),
         beamparts.second->status(),
         -1, // no mother1
         -1, // no mother2
         -1, // first daughter not known yet
         -1, // last daughter not known yet
         beamparts.second->momentum().px(),
         beamparts.second->momentum().py(),
         beamparts.second->momentum().pz(),
         beamparts.second->momentum().e(),
         0, // no production vertex, so zero?
         0,
         0,
         0
   );
   partonMemory[beamparts.second->barcode()] = 1;
   
   Int_t arrayID = 2; // start counting IDs after the beam particles

   // If the beam particles end
   // in the same vertex, this has to be done only once, otherwise, it
   // has to be done twice (else we lose the "second beam part"
   // branch)

   for (Int_t ibeamPart =0 ; ibeamPart < 2; ibeamPart++) {
     //     std::cout << "1" << std::endl;
     
     HepMC::GenVertex * startVertex = 0;
     if (ibeamPart == 0) startVertex = beamparts.first  ->end_vertex();
     else                {
       startVertex = beamparts.second ->end_vertex();
       if (beamParticlesToTheSameVertex) break; // both beam particles end to the same vertex: no need to do this twice.
     }
     //     std::cout << "2" << std::endl;
     
     // Do first vertex as a special case for performance, since its often has the most daughters and both mothers are known
     Int_t firstDaughter = arrayID;
     for (HepMC::GenVertex::particles_out_const_iterator iter = startVertex->particles_out_const_begin();
          iter != startVertex->particles_out_const_end();
          ++iter) {
       new((*array)[arrayID]) TParticle(
                                        (*iter)->pdg_id(),
                                        (*iter)->status(),
                                        0, // beam particle 1
                                        1, // beam particle 2
                                        -1, // first daughter not known yet
                                        -1, // last daughter not known yet
                                        (*iter)->momentum().px(),
                                        (*iter)->momentum().py(),
                                        (*iter)->momentum().pz(),
                                        (*iter)->momentum().e(),
                                        beamparts.first->end_vertex()->position().x(),
                                        beamparts.first->end_vertex()->position().y(),
                                        beamparts.first->end_vertex()->position().z(),
                                        beamparts.first->end_vertex()->position().t()
                                        );
       partonMemory[(*iter)->barcode()] = arrayID;
       ++arrayID;
     }
     Int_t lastDaughter = arrayID-1;
     if (ibeamPart == 0){
       ((TParticle*)array->AddrAt(0))->SetFirstDaughter(firstDaughter); // beam particle 1
       ((TParticle*)array->AddrAt(0))->SetLastDaughter(lastDaughter);
       // If both beam particles end in the same vertex, we have to reference the daughters of the second beam particle here (this loop is only done once). Otherwise, we do it when ibeamPart is == 1
       if(beamParticlesToTheSameVertex) {
         ((TParticle*)array->AddrAt(1))->SetFirstDaughter(firstDaughter); // beam particle 2
         ((TParticle*)array->AddrAt(1))->SetLastDaughter(lastDaughter);
       }
     } else {
       ((TParticle*)array->AddrAt(1))->SetFirstDaughter(firstDaughter); // beam particle 2
       ((TParticle*)array->AddrAt(1))->SetLastDaughter(lastDaughter);

     }
     
     // Then we loop over all other vertices. 
     // Build vertex list by exploring tree and sorting such that
     // daughters comes after mothers
     list<HepMC::GenVertex*> vertexList;
     set<int> vertexSearchSet;
     ExploreVertex(startVertex,vertexList,vertexSearchSet,requireSecondMotherBeforeDaughters);

     // Analyze each vertex
     for (list<HepMC::GenVertex*>::iterator i = vertexList.begin(); i != vertexList.end(); ++i) {
       HepMC::GenVertex * vertex = (*i);
       //       std::cout << "Vertex: " << vertex->barcode() << std::endl;
       
       //      std::cout << "Processing vertex: " << vertex->barcode() << std::endl;
      
       // first establish mother-daughter relations (look at particles incoming in the vertex)
       HepMC::GenVertex::particles_in_const_iterator iterInParticles = vertex->particles_in_const_begin();
       if (iterInParticles == vertex->particles_in_const_end()) {
         return "Particle without a mother, and its not a beam particle!\n";
       }
       int motherA = partonMemory[(*iterInParticles)->barcode()];
       //      std::cout << "MotherA: "  << (*iterInParticles)->barcode() << ", " << motherA << std::endl;
      
       if (((TParticle*)array->AddrAt(motherA))->GetFirstDaughter() > -1 && motherA > 1) { // FIXME Temp hack: if motherA <2 it means it was not in the array (0 is the beam particle)
         errMsgStream << Form("Trying to assign new daughters to a particle that already has daughters defined! (motherA, barcode = %d, event = %d)\n",(*iterInParticles)->barcode(), genEvent->event_number());
         return errMsgStream.str();
         ((TParticle*)array->AddrAt(motherA))->Print();
       }
       ++iterInParticles;
       int motherB = -1;
       if (iterInParticles != vertex->particles_in_const_end()) {
         motherB = partonMemory[(*iterInParticles)->barcode()];
         //         std::cout << "MotherB: "  << (*iterInParticles)->barcode() << std::endl;
         if (((TParticle*)array->AddrAt(motherB))->GetFirstDaughter() > -1 && motherB > 1) { // FIXME Temp hack: if motherB < 2 it means it was not in the array (0 is the beam particle)
           errMsgStream << Form("Trying to assign new daughters to a particle that already has daughters defined! (motherB, barcode = %d, event = %d)\n",(*iterInParticles)->barcode(), genEvent->event_number());
           return errMsgStream.str();
         }
         ++iterInParticles;
         if (iterInParticles != vertex->particles_in_const_end() && (*iterInParticles)->pdg_id() != 93) { // If it's a string it could have many partons attached.. In that case, we ignore mother-daughter relationships
           //           std::cout << "Daughters" << std::endl;
           errMsgStream << Form ("Particle with more than two mothers! (Vertex Barcode: %d, PDG: %d)\n", (*iterInParticles)->barcode(), (*iterInParticles)->pdg_id()) ;
           //           return Form ("Particle with more than two mothers! (Vertex Barcode: %d, PDG: %d)", (*iterInParticles)->barcode(), (*iterInParticles)->pdg_id()) ;
            
         }
       }
       if (motherB > -1 && motherB < motherA) {
         int swap = motherA; motherA = motherB; motherB = swap;
       }

       // add the particles to the array, important that they are add in succession with respect to arrayID
       
       firstDaughter = arrayID;
       for (HepMC::GenVertex::particles_in_const_iterator iterOutParticles = vertex->particles_out_const_begin();
            iterOutParticles != vertex->particles_out_const_end();
            ++iterOutParticles) {
         //        std::cout << "Adding daughter: " << (*iterOutParticles)->barcode() << std::endl;
        
         new((*array)[arrayID]) TParticle(
                                          (*iterOutParticles)->pdg_id(),
                                          (*iterOutParticles)->status(),
                                          motherA, // mother 1
                                          motherB, // mother 2
                                          -1, // first daughter, if applicable, not known yet
                                          -1, // last daughter, if applicable, not known yet
                                          (*iterOutParticles)->momentum().px(),
                                          (*iterOutParticles)->momentum().py(),
                                          (*iterOutParticles)->momentum().pz(),
                                          (*iterOutParticles)->momentum().e(),
                                          vertex->position().x(),
                                          vertex->position().y(),
                                          vertex->position().z(),
                                          vertex->position().t()
                                          );
         partonMemory[(*iterOutParticles)->barcode()] = arrayID;
         ++arrayID;
       }
       lastDaughter = arrayID-1;
       if (lastDaughter < firstDaughter) {
         return "Vertex with no out particles, should not be possible!";
       }
       // update mother with daughter interval
       ((TParticle*)array->AddrAt(motherA))->SetFirstDaughter(firstDaughter);
       ((TParticle*)array->AddrAt(motherA))->SetLastDaughter(lastDaughter);
       if (motherB > -1) {
         ((TParticle*)array->AddrAt(motherB))->SetFirstDaughter(firstDaughter);
         ((TParticle*)array->AddrAt(motherB))->SetLastDaughter(lastDaughter);
       }
     }
   }

   // Printf("Particles in event: %d, added: %d, entries: %d/%d", nParticlesInHepMC, arrayID, array->GetEntries(), array->GetEntriesFast());
   std::cout << errMsgStream.str() << std::endl;
   
   return "";
}


void THepMCParser::ExploreVertex(HepMC::GenVertex * vertex, list<HepMC::GenVertex*> & vertexList, set<int> & vertexSearchSet, bool requireSecondMotherBeforeDaughters)
{
   // Prepare vertex list, and sort vertices so that mothers come before daughters (if required)

   // Loop over all particles exiting from our start vertex 
   for (HepMC::GenVertex::particles_out_const_iterator partOut = vertex->particles_out_const_begin();
         partOut != vertex->particles_out_const_end();
         ++partOut) {
      HepMC::GenVertex * testVertex = (*partOut)->end_vertex();
      if (testVertex) {
         bool foundVertex = vertexSearchSet.find((*partOut)->end_vertex()->barcode()) != vertexSearchSet.end(); // If this is true, the vertex was already found
         if (requireSecondMotherBeforeDaughters) {
            // redo this algorithem to move subtree instead of node....
            // its not completely error proof in its current implementation even though the error is extremely rare

            // Loop over all already found vertices if the vertex was already found, and remove the previous instance
            if (foundVertex) for (list<HepMC::GenVertex*>::iterator ivert = vertexList.begin(); ivert != vertexList.end(); ++ivert) {
               if ((*ivert)->barcode() == testVertex->barcode()) {
                  vertexList.erase(ivert);
                  cout << " it happened, the vertex parsing order had to be changed " << endl;
                  break;
               }
            } else {
              //otherwise, just add it to the list
               vertexSearchSet.insert((*partOut)->end_vertex()->barcode());
            }
            vertexList.push_back(testVertex);
            // Follow daughter recursively
            if (!foundVertex) ExploreVertex(testVertex,vertexList,vertexSearchSet,requireSecondMotherBeforeDaughters);

         } else {
            if (!foundVertex) {
              // If we don't care about having all mothers first, we just add it to list and set if not already found
               vertexSearchSet.insert((*partOut)->end_vertex()->barcode());
               vertexList.push_back(testVertex);
               ExploreVertex(testVertex,vertexList,vertexSearchSet,requireSecondMotherBeforeDaughters);
            }
         }
      }
   }
}



const char * THepMCParser::fgHeavyIonHeaderBranchString = "Ncoll_hard/I:Npart_proj:Npart_targ:Ncoll:spectator_neutrons:spectator_protons:N_Nwounded_collisions:Nwounded_N_collisions:Nwounded_Nwounded_collisions:impact_parameter/F:event_plane_angle:eccentricity:sigma_inel_NN";
const char * THepMCParser::fgPdfHeaderBranchString = "id1/I:id2:pdf_id1:pdf_id2:x1/D:x2:scalePDF:pdf1:pdf2";

string THepMCParser::ParseGenEvent2HeaderStructs(HepMC::GenEvent * genEvent, HeavyIonHeader_t & heavyIonHeader, PdfHeader_t & pdfHeader, bool fillZeroOnMissingHeavyIon, bool fillZeroOnMissingPdf)
{
   HepMC::HeavyIon * heavyIon = genEvent->heavy_ion();
   HepMC::PdfInfo * pdfInfo = genEvent->pdf_info();
   if ((!heavyIon && !fillZeroOnMissingHeavyIon) || (!pdfInfo && !fillZeroOnMissingPdf)) {
      return "HeavyIonInfo and/or PdfInfo not defined for this event, did you read it with IO_GenEvent?";
   }
   if (heavyIon) {
      heavyIonHeader.Ncoll_hard = heavyIon->Ncoll_hard();
      heavyIonHeader.Npart_proj = heavyIon->Npart_proj();
      heavyIonHeader.Npart_targ = heavyIon->Npart_targ();
      heavyIonHeader.Ncoll = heavyIon->Ncoll();
      heavyIonHeader.spectator_neutrons = heavyIon->spectator_neutrons();
      heavyIonHeader.spectator_protons = heavyIon->spectator_protons();
      heavyIonHeader.N_Nwounded_collisions = heavyIon->N_Nwounded_collisions();
      heavyIonHeader.Nwounded_N_collisions = heavyIon->Nwounded_N_collisions();
      heavyIonHeader.Nwounded_Nwounded_collisions = heavyIon->Nwounded_Nwounded_collisions();
      heavyIonHeader.impact_parameter = heavyIon->impact_parameter();
      heavyIonHeader.event_plane_angle = heavyIon->event_plane_angle();
      heavyIonHeader.eccentricity = heavyIon->eccentricity();
      heavyIonHeader.sigma_inel_NN = heavyIon->sigma_inel_NN();
   } else {
      heavyIonHeader.Ncoll_hard = 0;
      heavyIonHeader.Npart_proj = 0;
      heavyIonHeader.Npart_targ = 0;
      heavyIonHeader.Ncoll = 0;
      heavyIonHeader.spectator_neutrons = 0;
      heavyIonHeader.spectator_protons = 0;
      heavyIonHeader.N_Nwounded_collisions = 0;
      heavyIonHeader.Nwounded_N_collisions = 0;
      heavyIonHeader.Nwounded_Nwounded_collisions = 0;
      heavyIonHeader.impact_parameter = 0.0;
      heavyIonHeader.event_plane_angle = 0.0;
      heavyIonHeader.eccentricity = 0.0;
      heavyIonHeader.sigma_inel_NN = 0.0;
   }
   if (pdfInfo) {
      pdfHeader.id1 = pdfInfo->id1();
      pdfHeader.id2 = pdfInfo->id2();
      pdfHeader.pdf_id1 = pdfInfo->pdf_id1();
      pdfHeader.pdf_id2 = pdfInfo->pdf_id2();
      pdfHeader.x1 = pdfInfo->x1();
      pdfHeader.x2 = pdfInfo->x2();
      pdfHeader.scalePDF = pdfInfo->scalePDF();
      pdfHeader.pdf1 = pdfInfo->pdf1();
      pdfHeader.pdf2 = pdfInfo->pdf2();
   } else {
      pdfHeader.id1 = 0;
      pdfHeader.id2 = 0;
      pdfHeader.pdf_id1 = 0;
      pdfHeader.pdf_id2 = 0;
      pdfHeader.x1 = 0.0;
      pdfHeader.x2 = 0.0;
      pdfHeader.scalePDF = 0.0;
      pdfHeader.pdf1 = 0.0;
      pdfHeader.pdf2 = 0.0;
   }
   return "";
}







