// $Id$
/**
 * @file compareDigitReaders.C
 * @brief Translate TPC DDL file into an ascii dump
 *
 * Usage:
 * <pre>
 *   aliroot -b -q -l raw2ascii.C'("infile", "outfile", slice, patch, "reader")' | tee compareDigitReaders.log
 * </pre>
 *
 * The macro translates a TPC DDL raw data file into an ascii dump using
 * a small chain with the AliHLTTPCDigitDumpComponent. The only mandatory
 * parameter is the input file.
 *   - outfile might be skipped -> default infile.ascii
 *   - optional parameter 'slice':
 *     can be skipped if infile obeys DDL naming rule
 *   - optional parameter 'patch':
 *     can be skipped if infile obeys DDL naming rule
 *   - optional parameter 'reader': default 'decoder'
 *     packed (AliRawStream), decoder (AliAltroDecoder)
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 */
void usage();

int raw2ascii(const char* infile, const char* outfile=NULL, int slice=-1, int part=-1, const char* reader="decoder")
{

  if (!infile) {
    usage();
    return 0;
  }

  AliHLTSystem gHLT(0x79);
  gHLT.LoadComponentLibraries("libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so");

  TString outFile;
  if (outfile) outFile=outfile;
  else outFile.Form("%s.ascii", infile);

  if (slice<0 || part<0) {
    TString number=infile;
    if (number.EndsWith(".ddl")) {
      number.ReplaceAll(".ddl", "");
      TObjArray* pTokens=number.Tokenize("_");
      if (pTokens) {
	if (pTokens->GetEntries()>1) {
	  number=((TObjString*)pTokens->At(pTokens->GetEntries()-1))->GetString();
	}
      }
      if (number.IsDigit()) {
	int ddlid=number.Atoi();
	if (ddlid<840) {
	  slice=(ddlid-768)/2;
	  part=ddlid%2; 
	} else {
	  slice=(ddlid-840)/4;
	  part=(ddlid%4)+2;
	}
      }
    }
    if (slice<0 || part<0) {
      cout << "error: can not derive slice and partition number from file name " << infile << endl;
      return -1;
    }
  }

  TString arg;
  arg.Form("-datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -datafile %s ", slice, slice, part, part, infile);
  AliHLTConfiguration pubconf("publisher", "FilePublisher", NULL , arg.Data());

  arg.Form("-datafile %s -digitreader %s -concatenate-blocks -idfmt= -unsorted", outFile.Data(), reader);
  AliHLTConfiguration dump("dump", "TPCDigitDump", "publisher", arg.Data());

  gHLT.BuildTaskList("dump");
  gHLT.Run();
}

int raw2ascii()
{
  usage();
  return 0;
}

void usage() {
  cout << "raw2ascii.C: translate TPC DDL file into an ascii dump" << endl << endl;
  cout << "    usage: aliroot -b -q -l raw2ascii.C'(\"infile\", \"outfile\", slice, patch, \"reader\")'" << endl;
  cout << "           outfile might be skipped -> default infile.ascii" << endl;
  cout << "           optional parameter 'slice':  can be skipped if infile obeys DDL naming rule" << endl;
  cout << "           optional parameter 'patch':  can be skipped if infile obeys DDL naming rule" << endl;
  cout << "           optional parameter 'reader': packed (AliRawStream), decoder (AliAltroDecoder)" << endl;
}

