// main54.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program compares the internal and LHAPDF implementations of the
// NNPDF 2.3 QCD+QED sets, for results and for timing.
// Author: Juan Rojo.
// Warning: this example is constructed to work for LHAPDF5.
// There seem to be differences when instead comparing with LHAPDF6.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  cout<<"\n NNPDF2.3 QED LO phenomenology \n "<<endl;
  cout<<"\n Check access to NNPDF2.3 LO QED sets \n "<<endl;

  // Generator.
  Pythia pythia;

  // Access the PDFs.
  int idBeamIn = 2212;
  string xmlPath = "../xmldoc/";
  Info info;

  // Grid of studied points.
  string xpdf[] = {"x*g","x*d","x*u","x*s"};
  double xlha[] = {1e-5, 1e-1};
  double Q2[] = { 2.0, 10000.0 };
  string setName;
  string setName_lha;

  // For timing checks.
  int const nq = 200;
  int const nx = 200;

  // Loop over all internal PDF sets in Pythia8
  // and compare with their LHAPDF5 correspondents.
  for (int iFitIn = 3; iFitIn < 5; iFitIn++) {

    // Constructor for internal PDFs.
    NNPDF pdfs_nnpdf( idBeamIn, iFitIn, xmlPath, &info);

    // Constructor for LHAPDF.
    if (iFitIn == 3) setName = "LHAPDF5:NNPDF23_nlo_as_0119_qed.LHgrid";
    if (iFitIn == 4) setName = "LHAPDF5:NNPDF23_nnlo_as_0119_qed.LHgrid";
    LHAPDF pdfs_nnpdf_lha( idBeamIn, setName, &info);
    cout << "\n PDF set = " << setName << " \n" << endl;

    // Check quarks and gluons.
    for (int f = 0; f < 4; f++) {
      for (int iq = 0; iq < 2; iq++) {
        cout << "  " << xpdf[f] << ", Q2 = " << Q2[iq] << endl;
        cout << "   x \t     Pythia8\t   LHAPDF\t   diff(%) " << endl;
        for (int ix = 0; ix < 2; ix++) {
          double a = pdfs_nnpdf.xf( f, xlha[ix], Q2[iq]);
          double b = pdfs_nnpdf_lha.xf( f, xlha[ix], Q2[iq]);
          double diff = 1e2 * fabs((a-b)/b);
          cout << scientific << xlha[ix] << " " << a << " " << b
               << " " << diff << endl;
          if (diff > 1e-8) exit(-10);
        }
      }
    }

    // Check photon.
    cout << "\n Now checking the photon PDF \n" << endl;
    for (int iq = 0; iq < 2; iq++) {
      cout << "  " << "x*gamma" << ", Q2 = " << Q2[iq] << endl;
      cout << "   x \t     Pythia8\t LHAPDF\t diff(%)" << endl;
      for (int ix = 0; ix < 2; ix++) {
        double a = pdfs_nnpdf.xf( 22, xlha[ix], Q2[iq]);
        double b = pdfs_nnpdf_lha.xf( 22, xlha[ix], Q2[iq]);
        double diff = 1e2 * fabs((a-b)/b);
        cout << scientific << xlha[ix] << " " << a << " " << b
             << " " << diff << endl;
        if(diff > 1e-8) exit(-10);
      }
    }

    // Now check the timings for evolution.
    cout << "\n Checking timings " << endl;

    clock_t t1 = clock();
    for (int f = -4; f < 4; f++) {
      for (int iq = 0; iq < nq; iq++) {
        double qq2 = 2.0 * pow( 1e6 / 2.0, double(iq)/nq);
        for (int ix = 0; ix < nx; ix++) {
          double xx = 1e-6 * pow( 9e-1 / 1e-6, double(ix)/nx);
          pdfs_nnpdf.xf(f,xx,qq2);
        }
      }
    }
    clock_t t2 = clock();
    cout << " NNPDF internal timing = " << (t2-t1)/(double)CLOCKS_PER_SEC
         << endl;

    t1=clock();
    for (int f = -4; f < 4; f++) {
      for (int iq = 0; iq < nq; iq++) {
        double qq2 = 2.0 * pow(1e6 / 2.0, double(iq)/nq);
        for (int ix = 0; ix < nx; ix++) {
          double xx = 1e-6 * pow( 9e-1 / 1e-6, double(ix)/nx);
          pdfs_nnpdf_lha.xf(f,xx,qq2);
        }
      }
    }
    t2=clock();
    cout << " NNPDF LHAPDF   timing = " << (t2-t1)/(double)CLOCKS_PER_SEC
         << endl;

  } // End loop over NNPDF internal sets

  // Done.
  cout << "\n Checked that LHAPDF and Internal Pythia8 give identical"
       << " results\n" << endl;

  return 0;
}
