#ifndef ALIGENHEPMCEVENTHEADER_H
#define ALIGENHEPMCEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenHepMCEventHeader : public AliGenEventHeader
{
public:
   AliGenHepMCEventHeader();
   AliGenHepMCEventHeader(const char* name);
   AliGenHepMCEventHeader(
         Int_t    Ncoll_hard,                    // Number of hard scatterings
         Int_t    Npart_proj,                    // Number of projectile participants
         Int_t    Npart_targ,                    // Number of target participants
         Int_t    Ncoll,                         // Number of NN (nucleon-nucleon) collisions
         Int_t    spectator_neutrons,            // Number of spectator neutrons
         Int_t    spectator_protons,             // Number of spectator protons
         Int_t    N_Nwounded_collisions,         // Number of N-Nwounded collisions
         Int_t    Nwounded_N_collisions,         // Number of Nwounded-N collisons
         Int_t    Nwounded_Nwounded_collisions,  // Number of Nwounded-Nwounded collisions
         Float_t  impact_parameter,              // Impact Parameter(in fm) of collision
         Float_t  event_plane_angle,             // Azimuthal angle of event plane
         Float_t  eccentricity,                  // eccentricity of participating nucleons in the transverse plane (as in phobos nucl-ex/0510031)
         Float_t  sigma_inel_NN,                 // nucleon-nucleon inelastic (including diffractive) cross-section
         Int_t    id1,        // flavour code of first parton
         Int_t    id2,        // flavour code of second parton
         Int_t    pdf_id1,    // LHAPDF set id of first parton
         Int_t    pdf_id2,    // LHAPDF set id of second parton
         Double_t x1,         // fraction of beam momentum carried by first parton ("beam side")
         Double_t x2,         // fraction of beam momentum carried by second parton ("target side")
         Double_t scalePDF,   // Q-scale used in evaluation of PDF's   (in GeV)
         Double_t pdf1,       // PDF (id1, x1, Q) - x*f(x)
         Double_t pdf2        // PDF (id2, x2, Q) - x*f(x)
   );
   virtual ~AliGenHepMCEventHeader() {}


   Int_t    Ncoll_hard() const {return fNcoll_hard;} // Number of hard scatterings
   Int_t    Npart_proj() const {return fNpart_proj;} // Number of projectile participants
   Int_t    Npart_targ() const {return fNpart_targ;} // Number of target participants
   Int_t    Ncoll() const {return fNcoll;} // Number of NN (nucleon-nucleon) collisions
   Int_t    spectator_neutrons() const {return fspectator_neutrons;} // Number of spectator neutrons
   Int_t    spectator_protons() const {return fspectator_protons;} // Number of spectator protons
   Int_t    N_Nwounded_collisions() const {return fN_Nwounded_collisions;} // Number of N-Nwounded collisions
   Int_t    Nwounded_N_collisions() const {return fNwounded_N_collisions;} // Number of Nwounded-N collisons
   Int_t    Nwounded_Nwounded_collisions() const {return fNwounded_Nwounded_collisions;} // Number of Nwounded-Nwounded collisions
   Float_t  impact_parameter() const {return fimpact_parameter;} // Impact Parameter(in fm) of collision
   Float_t  event_plane_angle() const {return fevent_plane_angle;} // Azimuthal angle of event plane
   Float_t  eccentricity() const {return feccentricity;} // eccentricity of participating nucleons in the transverse plane (as in phobos nucl-ex/0510031)
   Float_t  sigma_inel_NN() const {return fsigma_inel_NN;} // nucleon-nucleon inelastic (including diffractive) cross-section

   Int_t    id1() const {return fid1;} // flavour code of first parton
   Int_t    id2() const {return fid2;} // flavour code of second parton
   Int_t    pdf_id1() const {return fpdf_id1;} // LHAPDF set id of first parton
   Int_t    pdf_id2() const {return fpdf_id2;} // LHAPDF set id of second parton
   Double_t x1() const {return fx1;} // fraction of beam momentum carried by first parton ("beam side")
   Double_t x2() const {return fx2;} // fraction of beam momentum carried by second parton ("target side")
   Double_t scalePDF() const {return fscalePDF;} // Q-scale used in evaluation of PDF's   (in GeV)
   Double_t pdf1() const {return fpdf1;} // PDF (id1, x1, Q) - x*f(x)
   Double_t pdf2() const {return fpdf2;} // PDF (id2, x2, Q) - x*f(x)

   // convenience functions to check if the headers are containing information
   Bool_t   HeavyIonInfoValid();
   Bool_t   PDFValid();

protected:

   Int_t    fNcoll_hard;                    // Number of hard scatterings
   Int_t    fNpart_proj;                    // Number of projectile participants
   Int_t    fNpart_targ;                    // Number of target participants
   Int_t    fNcoll;                         // Number of NN (nucleon-nucleon) collisions
   Int_t    fspectator_neutrons;            // Number of spectator neutrons
   Int_t    fspectator_protons;             // Number of spectator protons
   Int_t    fN_Nwounded_collisions;         // Number of N-Nwounded collisions
   Int_t    fNwounded_N_collisions;         // Number of Nwounded-N collisons
   Int_t    fNwounded_Nwounded_collisions;  // Number of Nwounded-Nwounded collisions
   Float_t  fimpact_parameter;              // Impact Parameter(in fm) of collision
   Float_t  fevent_plane_angle;             // Azimuthal angle of event plane
   Float_t  feccentricity;                  // eccentricity of participating nucleons in the transverse plane (as in phobos nucl-ex/0510031)
   Float_t  fsigma_inel_NN;                 // nucleon-nucleon inelastic (including diffractive) cross-section

   Int_t    fid1;        // flavour code of first parton
   Int_t    fid2;        // flavour code of second parton
   Int_t    fpdf_id1;    // LHAPDF set id of first parton
   Int_t    fpdf_id2;    // LHAPDF set id of second parton
   Double_t fx1;         // fraction of beam momentum carried by first parton ("beam side")
   Double_t fx2;         // fraction of beam momentum carried by second parton ("target side")
   Double_t fscalePDF;   // Q-scale used in evaluation of PDF's   (in GeV)
   Double_t fpdf1;       // PDF (id1, x1, Q) - x*f(x)
   Double_t fpdf2;       // PDF (id2, x2, Q) - x*f(x)

   ClassDef(AliGenHepMCEventHeader, 2)  // Event header for HepMC event
};



#endif
