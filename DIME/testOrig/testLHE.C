// Read data from an LHE ascii file and create a mass plot of pi+pi- system
//
// mikael.mieskolainen@cern.ch

#include "Riostream.h"

void testLHE() {

   // Open LHE output ascii
   ifstream in;
   in.open("exerec1.dat");

   // LHE format variables
   Int_t nevents = 0;
   Int_t event = 0;
   Int_t npart = 0;

   Int_t pnr = 0;
   Int_t state = 0;
   Int_t pid = 0;

   Int_t md[4] = {0};
   Double_t v[4] = {0};
   Double_t mass = 0;
   Double_t empty = 0;

   // Mass histogram
   TH1* hM = new TH1D("hM","DIME #pi^{+}#pi^{-};M_{#pi^{+}#pi^{-}} #[]{GeV/#it{c}^{2}}", 100,  0.0, 4.0);

   // Event loop
   while (true) {

      // Read in event number and number of particles
      in >> event >> npart;

      if (!in.good()) break;
      ++nevents;

      TLorentzVector vSum(0,0,0,0);
      TLorentzVector pip(0,0,0,0);
      TLorentzVector pim(0,0,0,0);

      // Particle loop
      for (Int_t i = 0; i < npart; ++i) {
        in >> pnr >> state >> pid;
        in >> md[0] >> md[1] >> md[2] >> md[3];
        in >> v[0] >> v[1] >> v[2] >> v[3];
        in >> mass;
        in >> empty >> empty;

        // Create 4-vector
        printf("%0.3f %0.3f %0.3f %0.3f \n", v[0], v[1], v[2], v[3]);
        TLorentzVector vec(v[0], v[1], v[2], v[3]);

        if (pid == 211) { // Choose pions
           pip = vec;
        }
        else if (pid == -211) {
           pim = vec;
        }
      }
      // pi+pi- system
      vSum = pip + pim;

      // fill histogram with pi+pi- system invariant mass
      hM->Fill(vSum.M());

      printf("\n");
   }
   printf("Found %d events from the LHE output format file! \n\n", nevents);

   // Save plots as pdf
   hM->Draw(); c1->SaveAs("massLHE.pdf");

   // Close the text file
   in.close();
}
