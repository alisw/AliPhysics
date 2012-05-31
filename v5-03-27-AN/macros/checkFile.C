#if !defined( __CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <unistd.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#endif

Bool_t checkFile(const char* filename)
{
  // This macro returns FALSE in the following cases:
  // - the file doesn't exist
  // - the file is not readable
  // - the file is not a Root file
  // - the file is truncated
  // If the file is opened without any Error message,
  // the macro returns kTRUE

  // checkname is string with appended .check, i.e. Kinematics.root.check
  // The output of TFile::Open will be stored there
  TString checkname = filename;
  checkname += ".check";

  // remove checkname, if it already exists
  gSystem->Unlink(checkname.Data());

  // Copy and store the original stream desriptors
  Int_t origstdout = dup(fileno(stdout));
  Int_t origstderr = dup(fileno(stderr));

  // Redirect the output
  freopen(checkname.Data(),"w",stdout);
  dup2(fileno(stdout),fileno(stderr));

  // Try to open the file. The output goes to the file checkname
  TFile f(filename);

  // Restore the original stream descriptors
  dup2(origstdout,fileno(stdout));
  dup2(origstderr,fileno(stderr));

  // Prepare the grep command, executed by the shell, for example
  // grep Error Kinematics.root.check > /dev/null 2>&1
  TString command = "grep Error ";
  command += checkname;
  command += " > /dev/null 2>&1";
  // Execute the command. The return value is 0 if Error is found
  Bool_t result = gSystem->Exec(command.Data());

  // Remove the checkname
  gSystem->Unlink(checkname.Data());
  return result;
}
