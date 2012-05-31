// $Id$
/*
 * Process and reconstruct the binary HLTOUT files as written to disk
 * by the HLTOUT formatter on the online HLT system.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-HLTOUT-bin.C'("pattern", runno, keepFiles)'
 *   aliroot -b -q rec-HLTOUT-bin.C'("")' # for help
 * </pre>
 *
 * In the online HLT, output of the HLTOUT formatter can be written to disk.
 * The resulting binary files contain the full DDL data. The macro
 * reconstructs the events from a specified sample of binary files.
 * 
 * @ingroup alihlt_programs
 * @author Matthias.Richter@ift.uib.no
 */
void rec_HLTOUT_bin(const char* input="", int runno=0, bool keepFiles=false)
{
  if (!input || input[0]==0) {
    cerr << "rec-HLTOUT-bin.C: Process and reconstruct binary HLTOUT files" << endl;
    cerr << "===============================================================" << endl;
    cerr << "usage: aliroot -b -q -l rec-HLTOUT-bin.C'(\"pattern\" " << endl;
    cerr << "                                          runNo," << endl;
    cerr << "                                          keepFiles)'" << endl << endl;
    cerr << "  Parameter:" << endl;
    cerr << "       pattern   input file pattern, e.g. \"HLT-Output*.bin\"" << endl;
    cerr << "       runNo     run no, default 0" << endl;
    cerr << "       keepFiles keep the intermediate files, default false" << endl;
    cerr << "===============================================================" << endl;
    return;
  }

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }

  TString inputPattern=input;
  if (!inputPattern.BeginsWith("/")) {
    inputPattern=gSystem->pwd();
    inputPattern+="/";
    inputPattern+=input;
  }
  TString workingDir=gSystem->TempDirectory();
  TUUID uuid;
  workingDir+="/"; workingDir+=uuid.AsString(); workingDir+="/";
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // prepare the input
  // in order to run the normal AliRoot reconstruction over the data, the files
  // need to be assembled in a format readable by the AliRawReaderFile
  // A shell script creates the rawx folders in the temporary working directory
  // and links one binary input file in each
  //
  TString command="c=0; for file in ";
  command+=inputPattern;
  command+=" ; do directory=";
  command+=workingDir;
  command+="raw$c; mkdir -p $directory; (cd $directory; ln -s $file HLT_7680.ddl; touch run";
  command+=runno;
  command+="); let c++; done";
 
  if (gSystem->Exec(command)!=0) {
    cerr << "failed to execute script" << endl << "================================================" << endl;
    cerr << command << endl;
    return;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput(workingDir);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetOption("HLT", "loglevel=0x7c");
  rec.Run();

  if (keepFiles) {
    cout << "event(s) reconstructed from " << workingDir << endl;
  } else {
    command="rm -r ";
    command+=workingDir;
    gSystem->Exec(command);
  }
}
