// module identifier line...
// Author: Brian Thorsbro, 24/6-2014


#ifndef THEPMCPARSER_H
#define THEPMCPARSER_H

#include <string>
#include <list>
#include <set>
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

namespace HepMC {
  class IO_BaseClass;
  class GenVertex;
  class GenEvent;
}

class THepMCParser {

private:

   TTree * fTree;
   void init(HepMC::IO_BaseClass *);

   // 'name(pdgcode,invmass)' if db entry or 'pdgcode()' if no db entry
   static std::string GetParticleName(TParticle *);

   // Fold out the vertex tree in a list, such that daughters are always last
   static void ExploreVertex(HepMC::GenVertex *, std::list<HepMC::GenVertex*> &, std::set<int> &, bool);


public:

   // Header struct
   static const char * fgHeavyIonHeaderBranchString; // "Ncoll_hard/I,Npart_proj,Npart_targ,Ncoll,spectator_neutrons,spectator_protons,N_Nwounded_collisions,Nwounded_N_collisions,Nwounded_Nwounded_collisions,impact_parameter/F,event_plane_angle,eccentricity,sigma_inel_NN";
   struct HeavyIonHeader_t {
      Int_t   Ncoll_hard;                    // Number of hard scatterings
      Int_t   Npart_proj;                    // Number of projectile participants
      Int_t   Npart_targ;                    // Number of target participants
      Int_t   Ncoll;                         // Number of NN (nucleon-nucleon) collisions
      Int_t   spectator_neutrons;            // Number of spectator neutrons
      Int_t   spectator_protons;             // Number of spectator protons
      Int_t   N_Nwounded_collisions;         // Number of N-Nwounded collisions
      Int_t   Nwounded_N_collisions;         // Number of Nwounded-N collisons
      Int_t   Nwounded_Nwounded_collisions;  // Number of Nwounded-Nwounded collisions
      Float_t impact_parameter;              // Impact Parameter(in fm) of collision
      Float_t event_plane_angle;             // Azimuthal angle of event plane
      Float_t eccentricity;                  // eccentricity of participating nucleons in the transverse plane (as in phobos nucl-ex/0510031)
      Float_t sigma_inel_NN;                 // nucleon-nucleon inelastic (including diffractive) cross-section
   };
   static const char * fgPdfHeaderBranchString; // "id1/I,id2,pdf_id1,pdf_id2,x1/D,x2,scalePDF,pdf1,pdf2";
   struct PdfHeader_t {
      Int_t    id1;        // flavour code of first parton
      Int_t    id2;        // flavour code of second parton
      Int_t    pdf_id1;    // LHAPDF set id of first parton
      Int_t    pdf_id2;    // LHAPDF set id of second parton
      Double_t x1;         // fraction of beam momentum carried by first parton ("beam side")
      Double_t x2;         // fraction of beam momentum carried by second parton ("target side")
      Double_t scalePDF;   // Q-scale used in evaluation of PDF's   (in GeV)
      Double_t pdf1;       // PDF (id1, x1, Q) - x*f(x)
      Double_t pdf2;       // PDF (id2, x2, Q) - x*f(x)
   };

   // Default constructor/destructor stuff, don't inherit from this class unless you handle the tree pointer
   inline THepMCParser() : fTree(0) {;} // nullptr in c++11
   inline virtual ~THepMCParser() {;} // should be a memory management of the TTree...

   // The actual useful constructors, either take:
   //  - a file name for a file with HepMC data or
   //  - a HepMC event data structure
   THepMCParser(const char *);
   THepMCParser(HepMC::IO_BaseClass *);

   // Optional validators, set the argument to true for verbose output to STDERR
   // WARNING: including status code 2 may produce invalid flag when it is in fact valid, not recommended to use
   bool IsValidMotherDaughtersConsitency(bool useStdErr = false, bool requireSecondMotherBeforeDaughters = false);
   bool IsValidParticleInvariantMass(bool useStdErr = false, bool includeStatusCode2Particles = false);
   bool IsValidVertexInvariantMass(bool useStdErr = false, bool includeStatusCode2Particles = false);

   // Show the decay chain by point of view of some particle
   static std::string ListReactionChain(TClonesArray *, Int_t);

   // Access the TTree generated or write it to a file
   TTree * GetTTree();
   void WriteTTreeToFile(const char *);

   // ***** The constructors will parse all the events with this function for you *****
   // The work horse of the parser, takes one event and generates the corresponding TClonesArray of TParticles,
   // requireSecondMotherBeforeDaughters defaults to false
   //  - Success if length of return string is 0, the provided TClonesArray has been filled
   //  - Failure otherwise and the return string then contains the error message
   //
   // Note:
   // The TClonesArray must be initialized before hand
   // the array will be cleared and filled with TParticles
   // The capacity of the array will be set to the number of particles parsed
   // First mother will be before the daughters, but second mother may be after
   // The two first particles in the array will be the beam particles
   // The status code set on the particle is copied from HepMC, i.e.
   //  - 1 = final particle
   //  - 2 = transitory particle
   //  - 4 = beam particle
   // The function is static to enable customized wrappers around it
   static std::string ParseGenEvent2TCloneArray(HepMC::GenEvent *, TClonesArray *, bool requireSecondMotherBeforeDaughters = false);


   // Depending on the implementation of HepMC::IO_BaseClass there may be information on
   // heavy ions or parton distribution functions available. This function will pull them
   // out if available.
   // The caller must supply allocated structures which will then be filled.
   static std::string ParseGenEvent2HeaderStructs(HepMC::GenEvent *, HeavyIonHeader_t &, PdfHeader_t &, bool fillZeroOnMissingHeavyIon = true, bool fillZeroOnMissingPdf = true);


};

#endif
