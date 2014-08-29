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

#include "AliGenHepMCEventHeader.h"
ClassImp(AliGenHepMCEventHeader)


AliGenHepMCEventHeader::AliGenHepMCEventHeader():
   fNcoll_hard(0),
   fNpart_proj(0),
   fNpart_targ(0),
   fNcoll(0),
   fspectator_neutrons(0),
   fspectator_protons(0),
   fN_Nwounded_collisions(0),
   fNwounded_N_collisions(0),
   fNwounded_Nwounded_collisions(0),
   fimpact_parameter(0.0),
   fevent_plane_angle(0.0),
   feccentricity(0.0),
   fsigma_inel_NN(0.0),
   fid1(0),
   fid2(0),
   fpdf_id1(0),
   fpdf_id2(0),
   fx1(0.0),
   fx2(0.0),
   fscalePDF(0.0),
   fpdf1(0.0),
   fpdf2(0.0)
{
   // Default Constructor
}

AliGenHepMCEventHeader::AliGenHepMCEventHeader(const char* name):
   AliGenEventHeader(name),
   fNcoll_hard(0),
   fNpart_proj(0),
   fNpart_targ(0),
   fNcoll(0),
   fspectator_neutrons(0),
   fspectator_protons(0),
   fN_Nwounded_collisions(0),
   fNwounded_N_collisions(0),
   fNwounded_Nwounded_collisions(0),
   fimpact_parameter(0.0),
   fevent_plane_angle(0.0),
   feccentricity(0.0),
   fsigma_inel_NN(0.0),
   fid1(0),
   fid2(0),
   fpdf_id1(0),
   fpdf_id2(0),
   fx1(0.0),
   fx2(0.0),
   fscalePDF(0.0),
   fpdf1(0.0),
   fpdf2(0.0)
{
   // Constructor
}

AliGenHepMCEventHeader::AliGenHepMCEventHeader(
      Int_t    var_Ncoll_hard,                    // Number of hard scatterings
      Int_t    var_Npart_proj,                    // Number of projectile participants
      Int_t    var_Npart_targ,                    // Number of target participants
      Int_t    var_Ncoll,                         // Number of NN (nucleon-nucleon) collisions
      Int_t    var_spectator_neutrons,            // Number of spectator neutrons
      Int_t    var_spectator_protons,             // Number of spectator protons
      Int_t    var_N_Nwounded_collisions,         // Number of N-Nwounded collisions
      Int_t    var_Nwounded_N_collisions,         // Number of Nwounded-N collisons
      Int_t    var_Nwounded_Nwounded_collisions,  // Number of Nwounded-Nwounded collisions
      Float_t  var_impact_parameter,              // Impact Parameter(in fm) of collision
      Float_t  var_event_plane_angle,             // Azimuthal angle of event plane
      Float_t  var_eccentricity,                  // eccentricity of participating nucleons in the transverse plane (as in phobos nucl-ex/0510031)
      Float_t  var_sigma_inel_NN,                 // nucleon-nucleon inelastic (including diffractive) cross-section
      Int_t    var_id1,        // flavour code of first parton
      Int_t    var_id2,        // flavour code of second parton
      Int_t    var_pdf_id1,    // LHAPDF set id of first parton
      Int_t    var_pdf_id2,    // LHAPDF set id of second parton
      Double_t var_x1,         // fraction of beam momentum carried by first parton ("beam side")
      Double_t var_x2,         // fraction of beam momentum carried by second parton ("target side")
      Double_t var_scalePDF,   // Q-scale used in evaluation of PDF's   (in GeV)
      Double_t var_pdf1,       // PDF (id1, x1, Q) - x*f(x)
      Double_t var_pdf2        // PDF (id2, x2, Q) - x*f(x)
):
   fNcoll_hard(var_Ncoll_hard),
   fNpart_proj(var_Npart_proj),
   fNpart_targ(var_Npart_targ),
   fNcoll(var_Ncoll),
   fspectator_neutrons(var_spectator_neutrons),
   fspectator_protons(var_spectator_protons),
   fN_Nwounded_collisions(var_N_Nwounded_collisions),
   fNwounded_N_collisions(var_Nwounded_N_collisions),
   fNwounded_Nwounded_collisions(var_Nwounded_Nwounded_collisions),
   fimpact_parameter(var_impact_parameter),
   fevent_plane_angle(var_event_plane_angle),
   feccentricity(var_eccentricity),
   fsigma_inel_NN(var_sigma_inel_NN),
   fid1(var_id1),
   fid2(var_id2),
   fpdf_id1(var_pdf_id1),
   fpdf_id2(var_pdf_id2),
   fx1(var_x1),
   fx2(var_x2),
   fscalePDF(var_scalePDF),
   fpdf1(var_pdf1),
   fpdf2(var_pdf2)
{
   // The Constructor
}

Bool_t AliGenHepMCEventHeader::HeavyIonInfoValid() {
   return fNcoll_hard != 0 ||
         fNpart_proj != 0 ||
         fNpart_targ != 0 ||
         fNcoll != 0 ||
         fspectator_neutrons != 0 ||
         fspectator_protons != 0 ||
         fN_Nwounded_collisions != 0 ||
         fNwounded_N_collisions != 0 ||
         fNwounded_Nwounded_collisions != 0 ||
         fimpact_parameter != 0.0 ||
         fevent_plane_angle != 0.0 ||
         feccentricity != 0.0 ||
         fsigma_inel_NN != 0.0;
}

Bool_t AliGenHepMCEventHeader::PDFValid() {
   return fid1 != 0 ||
         fid2 != 0 ||
         fpdf_id1 != 0 ||
         fpdf_id2 != 0 ||
         fx1 != 0.0 ||
         fx2 != 0.0 ||
         fscalePDF != 0.0 ||
         fpdf1 != 0.0 ||
         fpdf2 != 0.0;
}
