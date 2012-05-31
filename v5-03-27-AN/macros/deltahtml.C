void deltahtml (Int_t delta=0) {
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
   gSystem->Load("libGeant3Dummy.so");   // a dummy version of Geant3
   gSystem->Load("PHOS/libPHOSdummy.so");        // the standard Alice classes 
   gSystem->Load("libgalice.so");        // the standard Alice classes 
   THtml html;
   html.MakeAll(delta);

}
