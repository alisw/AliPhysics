void vertexmacro (const char *filename="galice.root",Int_t evNumber=0) {
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L vertexmacro.C   //this loads the macro in memory
//     Root > vertexmacro();     //by default process first event   
//     Root > vertexmacro("galice.root");   //process file galice.root
//     Root > vertexmacro("galice.root",0); // process first event 
//                                             from galice.root file
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } // end if gClassTable...

      
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file){ file = new TFile(filename);printf("Opening new file\n");}
   printf("Root input file is %s\n",filename);

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
//

   for(Int_t evnt=0;evnt<=evNumber;evnt++){

      // define some variables for later use.
      Int_t      nparticles = gAlice->GetEvent(evnt);
      if (nparticles <= 0) continue; // get next event

      FindVertexs(evnt,1.0,16.82);

   } // end for evnt
   return;
}
