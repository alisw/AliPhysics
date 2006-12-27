//____________________________________________________________________
//
// $Id$
//
// Script that extracts the medium numbers corresponding to each
// registered detector.
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
/** @file    MediaTable.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Wed Dec 27 14:25:11 2006
    @brief   Script to make a file with a list of medium numbers for
    all active detectors.   
*/
/** Script to make a file with a list of medium numbers for all active
    detectors.  This should be executed in the same process as the one
    that makes the @c galice.root file.  For example
    @code 
    void 
    Simulate()
    {
      AliSimulation sim;
      AliLog::SetModuleDebugLevel("FMD", 1);
      sim.SetConfigFile("$(ALICE_ROOT)/FMD/Config.C");
      // sim.SetMakeSDigits("FMD");
      sim.SetMakeDigits("FMD ITS VZERO T0"); 
      sim.SetWriteRawData("FMD"); 
      gROOT->Macro("$(ALICE_ROOT)/FMD/scripts/MediaTable.C");
      TStopwatch w; 
      w.Start(); 
      sim.Run(1);  
      w.Stop(); 
      w.Print(); 
    }
    @endcode 
    It can also be executed after the fact, if one make sure to use
    the same configuration file as the one used for the simulation.
    For example 
    @verbatim 
    Root> gAlice->Init("$(ALICE_ROOT)/FMD/Config.C");
    Root> .x $ALICE_ROOT/FMD/scripts/MediaTable.C 
    Root> .q
    @endverbatim 
    The idea is to use the generated file to track where particles are
    produced in the simulation.  This is needed to do background
    studies in the FMD.  The file written contains a @c TObjArray of
    @c TObjString objects.  The @c TObjString objects has the name of
    the detector corresponding to the medium Id that the object is at
    in the array.  To read the file, do for example: 
    @code 
    TFile*      file   = TFile::Open("mediatable.root", "READ");
    TObjArray*  array  = (TObjArray*)(file->Get("mediatable"));
    Int_t       n      = array->GetEntries();
    TObjString* l      = 0;
    Int_t       min = 0;
    Int_t       j   = 0;
    for (j = 0; j < n; j++) {
      TObjString* s = static_cast<TObjString*>(array->At(j));
      if (!l || s->Compare(l) != 0) {
        if (l) 
	  std::cout << l->GetName() << "\t" << min 
                    << "\t" << j - 1 << std::endl;
        min = j;
        l = s;
      }
    }
    std::cout << l->GetName() << "\t" << min << "\t" << j << std::endl;
    file->Close();    
    @endcode 
    See also the script @c BackgroundStudy.C for more on how this is
    used in the FMD collaboration.   
*/
//____________________________________________________________________
void 
MediaTable()
{
  MediaMap m;
  TFile*     file   = TFile::Open("mediatable.root", "RECREATE");
  TObjArray* media  = gAlice->Modules();
  AliModule* module = 0;
  TIter next(media);
  TObjArray* array = new TObjArray(media->GetEntries());
  while ((module = static_cast<AliModule*>(next()))) {
    Int_t low  = module->LoMedium();
    Int_t high = module->HiMedium();
    for (Int_t j = low; j <=high; j++) {
      TObjString* o = new TObjString(module->GetName());
      std::cout << "Adding " << j << ":\t" << o->GetName() << std::endl;
      array->AddAtAndExpand(o, j);
    }
  }
  array->Write("mediatable", TObject::kSingleKey);
  file->Write();
  file->Close();
}

//
// EOF
//
