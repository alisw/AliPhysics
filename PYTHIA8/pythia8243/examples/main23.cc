// main23.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Examples how to write external derived classes
// that can be handed to Pythia for internal generation.
// 1) The MIXMAX random number generator; see include/Pythia8Plugins/MixMax.h.
// 2) A simple derived class for external beam momentum and vertex spread.
// 3) A simple derived class for external parton distributions.
// Warning: the parameters are not realistic.
// Two subruns are made, to compare runtime (and results)
// for the RanMar and MixMax random number generators.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/MixMax.h"

using namespace Pythia8;

//==========================================================================

// A derived class to set beam momentum and interaction vertex spread.

class MyBeamShape : public BeamShape {

public:

  // Constructor.
  MyBeamShape() {}

  // Initialize beam parameters.
  // In this particular example we will reuse the existing settings names
  // but with modified meaning, so init() in the base class can be kept.
  //virtual void init( Settings& settings, Rndm* rndmPtrIn);

  // Set the two beam momentum deviations and the beam vertex.
  virtual void pick();

};

//--------------------------------------------------------------------------

// Set the two beam momentum deviations and the beam vertex.
// Note that momenta are in units of GeV and vertices in mm,
// always with c = 1, so that e.g. time is in mm/c.

void MyBeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Set beam A transverse momentum deviation by a two-dimensional Gaussian.
  if (allowMomentumSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaPxA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxA  = sigmaPxA * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyA  = sigmaPyA * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevA * maxDevA);

    // Set beam A longitudinal momentum as a triangular shape.
    // Reuse sigmaPzA to represent maximum deviation in this case.
    if (sigmaPzA > 0.) {
      deltaPzA    = sigmaPzA * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) deltaPzA = -deltaPzA;
    }

    // Set beam B transverse momentum deviation by a two-dimensional Gaussian.
    do {
      totalDev = 0.;
      if (sigmaPxB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxB  = sigmaPxB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyB  = sigmaPyB * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevB * maxDevB);

    // Set beam B longitudinal momentum as a triangular shape.
    // Reuse sigmaPzB to represent maximum deviation in this case.
    if (sigmaPzB > 0.) {
      deltaPzB = sigmaPzB * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) deltaPzB = -deltaPzB;
    }
  }

  // Set beam vertex location by a two-dimensional Gaussian.
  if (allowVertexSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaVertexX > 0.) {
        gauss     = rndmPtr->gauss();
        vertexX   = sigmaVertexX * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexY > 0.) {
        gauss     = rndmPtr->gauss();
        vertexY   = sigmaVertexY * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevVertex * maxDevVertex);

    // Set beam B longitudinal momentum as a triangular shape.
    // This corresponds to two step-function beams colliding.
    // Reuse sigmaVertexZ to represent maximum deviation in this case.
    if (sigmaVertexZ > 0.) {
      vertexZ     = sigmaVertexZ * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) vertexZ = -vertexZ;

      // Set beam collision time flat between +-(sigmaVertexZ - |vertexZ|).
      // This corresponds to two step-function beams colliding (with v = c).
      vertexT = (2. * rndmPtr->flat() - 1.) * (sigmaVertexZ - abs(vertexZ));
    }

    // Add offset to beam vertex.
    vertexX      += offsetX;
    vertexY      += offsetY;
    vertexZ      += offsetZ;
    vertexT      += offsetT;
  }

}

//==========================================================================

// A simple scaling PDF. Not realistic; only to check that it works.

class Scaling : public PDF {

public:

  // Constructor.
  Scaling(int idBeamIn = 2212) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

};

//--------------------------------------------------------------------------

// No dependence on Q2, so leave out name for last argument.

void Scaling::xfUpdate(int, double x, double ) {

  // Valence quarks, carrying 60% of the momentum.
  double dv  = 4. * x * pow3(1. - x);
  double uv  = 2. * dv;

  // Gluons and sea quarks carrying the rest.
  double gl  = 2.  * pow5(1. - x);
  double sea = 0.4 * pow5(1. - x);

  // Update values
  xg    = gl;
  xu    = uv + 0.18 * sea;
  xd    = dv + 0.18 * sea;
  xubar = 0.18 * sea;
  xdbar = 0.18 * sea;
  xs    = 0.08 * sea;
  xc    = 0.04 * sea;
  xb    = 0.02 * sea;

  // Subdivision of valence and sea.
  xuVal = uv;
  xuSea = xubar;
  xdVal = dv;
  xdSea = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

int main() {

  // Number of events to generate. Max number of errors.
  int nEvent = 10000;
  int nAbort = 5;

  // Comparative statistics.
  double time[2];
  Hist nchRM("charged multiplicity RanMar", 100, -1., 399.);
  Hist nchMM("charged multiplicity MixMax", 100, -1., 399.);

  // Compare runtime (and results) for RanMar and MixMax random numbers.
  for (int iRNG = 0; iRNG < 2; ++iRNG) {
    clock_t start = clock();

    // Pythia generator.
    Pythia pythia;

    // Process selection.
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 20.");

    // LHC with acollinear beams in the x plane.
    // Use that default is pp with pz = +-7000 GeV, so this need not be set.
    pythia.readString("Beams:frameType = 3");
    pythia.readString("Beams:pxA = 1.");
    pythia.readString("Beams:pxB = 1.");

    // A class to generate beam parameters according to own parametrization.
    BeamShape* myBeamShape = new MyBeamShape();

    // Hand pointer to Pythia.
    // If you comment this out you get internal Gaussian-style implementation.
    pythia.setBeamShapePtr( myBeamShape);

    // Set up beam spread parameters - reused by MyBeamShape.
    pythia.readString("Beams:allowMomentumSpread = on");
    pythia.readString("Beams:sigmapxA = 0.1");
    pythia.readString("Beams:sigmapyA = 0.1");
    pythia.readString("Beams:sigmapzA = 5.");
    pythia.readString("Beams:sigmapxB = 0.1");
    pythia.readString("Beams:sigmapyB = 0.1");
    pythia.readString("Beams:sigmapzB = 5.");

    // Set up beam vertex parameters - reused by MyBeamShape.
    pythia.readString("Beams:allowVertexSpread = on");
    pythia.readString("Beams:sigmaVertexX = 0.3");
    pythia.readString("Beams:sigmaVertexY = 0.3");
    pythia.readString("Beams:sigmaVertexZ = 50.");
    // In MyBeamShape the time width is not an independent parameter.
    //pythia.readString("Beams:sigmaTime = 50.");

    // Optionally simplify generation and output.
    //pythia.readString("PartonLevel:all = off");
    pythia.readString("Next:numberShowEvent = 0");

    // A class to do random numbers externally. Hand pointer to Pythia.
    // You can provide four 32-bit unsigned integers to set random sequence.
    RndmEngine* rndm;
    if (iRNG == 1) {
      rndm = new MixMaxRndm( 0, 0, 0, 123);
      pythia.setRndmEnginePtr(rndm);
    }

    // Two classes to do the two PDFs externally. Hand pointers to Pythia.
    PDF* pdfAPtr = new Scaling(2212);
    PDF* pdfBPtr = new Scaling(2212);
    pythia.setPDFPtr( pdfAPtr, pdfBPtr);

    // Initialization.
    pythia.init();

    // Read out nominal energy.
    double eCMnom = pythia.info.eCM();

    // Local histograms.
    Hist eCM("center-of-mass energy deviation", 100, -20., 20.);
    Hist pXsum("net pX offset", 100, -1.0, 1.0);
    Hist pYsum("net pY offset", 100, -1.0, 1.0);
    Hist pZsum("net pZ offset", 100, -10., 10.);
    Hist pZind("individual abs(pZ) offset", 100, -10., 10.);
    Hist vtxX("vertex x position", 100, -1.0, 1.0);
    Hist vtxY("vertex y position", 100, -1.0, 1.0);
    Hist vtxZ("vertex z position", 100, -100., 100.);
    Hist vtxT("vertex time", 100, -100., 100.);
    Hist vtxZT("vertex |x| + |t|", 100, 0., 100.);

    // Begin event loop. Generate event.
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) {

        // List faulty events and quit if many failures.
        pythia.info.list();
        pythia.process.list();
        //pythia.event.list();
        if (++iAbort < nAbort) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      // Fill comparative histograms.
      int nch = 0;
      for (int i = 0; i < pythia.event.size(); ++i)
        if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nch;
      if (iRNG == 0) nchRM.fill( nch );
      else           nchMM.fill( nch );

      // Fill local histograms.
      double eCMnow = pythia.info.eCM();
      eCM.fill( eCMnow - eCMnom);
      pXsum.fill(  pythia.process[0].px() - 2. );
      pYsum.fill(  pythia.process[0].py() );
      pZsum.fill(  pythia.process[0].pz() );
      pZind.fill(  pythia.process[1].pz() - 7000. );
      pZind.fill( -pythia.process[2].pz() - 7000. );
      vtxX.fill(  pythia.process[0].xProd() );
      vtxY.fill(  pythia.process[0].yProd() );
      vtxZ.fill(  pythia.process[0].zProd() );
      vtxT.fill(  pythia.process[0].tProd() );
      double absSum = abs(pythia.process[0].zProd())
                    + abs(pythia.process[0].tProd());
      vtxZT.fill( absSum );

    // End of event loop. Statistics. Histograms.
    }
    pythia.stat();
    cout << eCM << pXsum << pYsum << pZsum << pZind
         << vtxX << vtxY << vtxZ << vtxT << vtxZT;

    // Remove created classes.
    delete myBeamShape;
    if (iRNG == 1) delete rndm;
    delete pdfAPtr;
    delete pdfBPtr;

    // Check time; end of loop over random number generators. Done.
    clock_t stop = clock();
    time[iRNG] = double(stop - start) / double(CLOCKS_PER_SEC);
  }
  cout << nchRM << nchMM;
  cout << "\n ========================================================="
       << "============" << endl;
  cout << fixed << setprecision(3) << "\n Time for run with RanMar is "
       << time[0] << " s and for run with MixMax " << time[1] << " s"<< endl;
  cout << "\n ========================================================="
       << "============" << endl;

  return 0;
}
