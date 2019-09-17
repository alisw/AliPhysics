// LHAPDF6.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the LHAPDF6 PDF plugin class.

#ifndef Pythia8_LHAPDF6_H
#define Pythia8_LHAPDF6_H

#include "Pythia8/PartonDistributions.h"
#include "LHAPDF/LHAPDF.h"

namespace Pythia8 {

//==========================================================================

// Global tracking of opened PDF sets.

//--------------------------------------------------------------------------

namespace LHAPDF6Interface {

//--------------------------------------------------------------------------

// Class to hold a PDF set, its information, and its uncertainty sets.

  class PdfSets {

  public:

    // Constructors.
    PdfSets() {;}
    PdfSets(string setName) : info(::LHAPDF::PDFSet(setName)),
      pdfs(vector< ::LHAPDF::PDF* >(info.size(), 0)) {;}

    // Access a PDF set.
    ::LHAPDF::PDF *operator[](unsigned int member) {
      if (!pdfs[member]) pdfs[member] = info.mkPDF(member);
      return pdfs[member];
    }

    // Get number of PDF sets.
    int size() {return pdfs.size();}

    // PDF sets and info.
    ::LHAPDF::PDFSet info;
    vector< ::LHAPDF::PDF* > pdfs;

  };

//--------------------------------------------------------------------------

// Class to globally track all open PDF sets.

  class PdfTracker {

  public:

    // Destructor, all PDFs from LHAPDF::mkPDF must be deleted.
    ~PdfTracker() {
      for (map<int, PdfSets>::iterator pdf = pdfs.begin();
           pdf != pdfs.end(); ++pdf)
        for (int iMem = 0; iMem < (int)pdf->second.size(); ++iMem)
          if (pdf->second.pdfs[iMem]) delete pdf->second.pdfs[iMem];
    }

    // Find and return a requested PDF set.
    PdfSets *find(string setName) {
      int id = ::LHAPDF::lookupLHAPDFID(setName, 0);
      if (id < 0) return 0;
      else if (pdfs.find(id) == pdfs.end()) pdfs[id] = PdfSets(setName);
      return &pdfs[id];
    }

  private:

    // Map to hold open PDF set information.
    map<int, PdfSets> pdfs;

  };

//--------------------------------------------------------------------------

// Define opened PDF sets global variable.

  PdfTracker pdfTracker;

}

//==========================================================================

// Provide interface to the LHAPDF6 library of parton densities.

class LHAPDF6 : public PDF {

public:

  // Constructor.
  LHAPDF6(int idBeamIn, string setName, int member, int, Info* infoPtr)
    : PDF(idBeamIn), pdf(0), extrapol(false)
    { init(setName, member, infoPtr); }

  // Allow extrapolation beyond boundaries (not implemented).
  void setExtrapolate(bool extrapolIn) {extrapol = extrapolIn;}

private:

  // The LHAPDF objects.
  LHAPDF6Interface::PdfSets *pdfs;
  ::LHAPDF::PDF *pdf;
  ::LHAPDF::Extrapolator *ext;
  bool extrapol;

  // Initialization of PDF set.
  void init(string setName, int member, Info* infoPtr);

  // Update parton densities.
  void xfUpdate(int id, double x, double Q2);

  // Check whether x and Q2 values fall inside the fit bounds.
  bool insideBounds(double x, double Q2) {
    return (x > pdf->xMin()  &&  x < pdf->xMax()
            && Q2 > pdf->q2Min() && Q2 < pdf->q2Max());}

  // Return the running alpha_s shipped with the LHAPDF set.
  double alphaS(double Q2) { return pdf->alphasQ2(Q2); }

  // Return quark masses used in the PDF fit.
  double muPDFSave, mdPDFSave, mcPDFSave, msPDFSave, mbPDFSave;
  double mQuarkPDF(int id) {
    if (abs(id) == 1) return mdPDFSave;
    if (abs(id) == 2) return muPDFSave;
    if (abs(id) == 3) return msPDFSave;
    if (abs(id) == 4) return mcPDFSave;
    if (abs(id) == 5) return mbPDFSave;
    return -1.;
 }

  // Calculate uncertainties using the LHAPDF prescription.
  void calcPDFEnvelope(int, double, double, int);
  void calcPDFEnvelope(pair<int,int>, pair<double,double>, double, int);
  PDFEnvelope pdfEnvelope;
  PDFEnvelope getPDFEnvelope() {return pdfEnvelope;}
  static const double PDFMINVALUE;

  int nMembersSave;
  int nMembers() { return nMembersSave; }

};

//--------------------------------------------------------------------------

// Constants.

const double LHAPDF6::PDFMINVALUE = 1e-10;

//--------------------------------------------------------------------------

// Initialize a parton density function from LHAPDF6.

void LHAPDF6::init(string setName, int member, Info *info) {
  isSet = false;

  // Initialize the LHAPDF sets.
  pdfs = LHAPDF6Interface::pdfTracker.find(setName);
  if (!pdfs) {
    info->errorMsg("Error in LHAPDF6::init: unknown PDF " + setName);
    return;
  } else if ((*pdfs).size() == 0) {
    info->errorMsg("Error in LHAPDF6::init: could not initialize PDF "
                   + setName);
    return;
  } else if (member >= (*pdfs).size()) {
    info->errorMsg("Error in LHAPDF6::init: " + setName
                   + " does not contain requested member");
    return;
  }
  pdf = (*pdfs)[member];
  isSet = true;

  // Store quark masses used in PDF fit.
  muPDFSave = pdf->info().get_entry_as<double>("MUp");
  mdPDFSave = pdf->info().get_entry_as<double>("MDown");
  mcPDFSave = pdf->info().get_entry_as<double>("MCharm");
  msPDFSave = pdf->info().get_entry_as<double>("MStrange");
  mbPDFSave = pdf->info().get_entry_as<double>("MBottom");

  nMembersSave  = pdf->info().get_entry_as<int>("NumMembers");

}

//--------------------------------------------------------------------------

// Give the parton distribution function set from LHAPDF6.

void LHAPDF6::xfUpdate(int, double x, double Q2) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  if (x < pdf->xMin() && !extrapol) x = pdf->xMin();
  if (x > pdf->xMax() )    x = pdf->xMax();
  if (Q2 < pdf->q2Min() ) Q2 = pdf->q2Min();
  if (Q2 > pdf->q2Max() ) Q2 = pdf->q2Max();

  // Update values.
  xg     = pdf->xfxQ2(21, x, Q2);
  xu     = pdf->xfxQ2(2,  x, Q2);
  xd     = pdf->xfxQ2(1,  x, Q2);
  xs     = pdf->xfxQ2(3,  x, Q2);
  xubar  = pdf->xfxQ2(-2, x, Q2);
  xdbar  = pdf->xfxQ2(-1, x, Q2);
  xsbar  = pdf->xfxQ2(-3, x, Q2);
  xc     = pdf->xfxQ2(4,  x, Q2);
  xb     = pdf->xfxQ2(5,  x, Q2);
  xgamma = pdf->xfxQ2(22, x, Q2);

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Calculate uncertainties using the LHAPDF prescription.

void LHAPDF6::calcPDFEnvelope(int idNow, double xNow, double Q2NowIn,
  int valSea) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNow < pdf->xMin() && !extrapol)
            ? pdf->xMin() : xNow;
  if (x1 > pdf->xMax() ) x1 = pdf->xMax();
  double Q2Now = (Q2NowIn < pdf->q2Min() )
               ? pdf->q2Min() : Q2NowIn;
  if (Q2Now > pdf->q2Max() ) Q2Now = pdf->q2Max();

  // Loop over the members.
  vector<double> xfCalc((*pdfs).size());
  for(int iMem = 0; iMem < (*pdfs).size(); ++iMem) {
    if (valSea==0 || (idNow != 1 && idNow != 2)) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(idNow, x1, Q2Now);
    } else if (valSea==1 && (idNow == 1 || idNow == 2 )) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(idNow, x1, Q2Now) -
        (*pdfs)[iMem]->xfxQ2(-idNow, x1, Q2Now);
    } else if (valSea==2 && (idNow == 1 || idNow == 2 )) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(-idNow, x1, Q2Now);
    }
  }

  // Calculate the uncertainty.
  ::LHAPDF::PDFUncertainty xfErr = (*pdfs).info.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;
}

//--------------------------------------------------------------------------

// Calculate uncertainties using the LHAPDF prescription.

void LHAPDF6::calcPDFEnvelope(pair<int,int> idNows, pair<double,double> xNows,
  double Q2NowIn, int valSea) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNows.first < pdf->xMin() && !extrapol)
            ? pdf->xMin() : xNows.first;
  if (x1 > pdf->xMax() ) x1 = pdf->xMax();
  double x2 = (xNows.second < pdf->xMin() && !extrapol)
            ? pdf->xMin() : xNows.second;
  if (x2 > pdf->xMax() ) x2 = pdf->xMax();
  double Q2Now = (Q2NowIn < pdf->q2Min() )
               ? pdf->q2Min() : Q2NowIn;
  if (Q2Now > pdf->q2Max() ) Q2Now = pdf->q2Max();

  // Loop over the members.
  vector<double> xfCalc((*pdfs).size());
  pdfEnvelope.pdfMemberVars.resize((*pdfs).size());
  for(int iMem = 0; iMem < (*pdfs).size(); ++iMem) {
    if        (valSea == 0 || (idNows.first != 1 && idNows.first != 2 ) ) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(idNows.first, x1, Q2Now);
    } else if (valSea == 1 && (idNows.first == 1 || idNows.first == 2)) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(idNows.first, x1, Q2Now)
        - (*pdfs)[iMem]->xfxQ2(-idNows.first, x1, Q2Now);
    } else if (valSea == 2 && (idNows.first == 1 || idNows.first == 2 )) {
      xfCalc[iMem] = (*pdfs)[iMem]->xfxQ2(-idNows.first, x1, Q2Now);
    }
    xfCalc[iMem] = max(0.0, xfCalc[iMem]);
    if        (valSea == 0 || (idNows.second != 1 && idNows.second != 2)) {
      xfCalc[iMem] /= max
        (PDFMINVALUE, (*pdfs)[iMem]->xfxQ2(idNows.second, x2, Q2Now));
    } else if (valSea == 1 && (idNows.second == 1 || idNows.second == 2 )) {
      xfCalc[iMem] /= max
        ((*pdfs)[iMem]->xfxQ2(idNows.second, x2, Q2Now) - (*pdfs)[iMem]->xfxQ2
         (-idNows.second, x2, Q2Now), PDFMINVALUE);
    } else if (valSea == 2 && (idNows.second == 1 || idNows.second == 2 )) {
      xfCalc[iMem] /= max
        ((*pdfs)[iMem]->xfxQ2(-idNows.second, x2, Q2Now), PDFMINVALUE);
    }
    pdfEnvelope.pdfMemberVars[iMem] = xfCalc[iMem];
  }

  // Calculate the uncertainty.
  ::LHAPDF::PDFUncertainty xfErr = (*pdfs).info.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;

}

//--------------------------------------------------------------------------

// Define external handles to the plugin for dynamic loading.

extern "C" LHAPDF6* newLHAPDF(int idBeamIn, string setName, int member,
                               Info* infoPtr) {
  return new LHAPDF6(idBeamIn, setName, member, 1, infoPtr);

}

extern "C" void deleteLHAPDF(LHAPDF6* pdf) {
  delete pdf;

}

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_LHAPDF6_H
