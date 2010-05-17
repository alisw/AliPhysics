/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Author: Mihaela Gheata, 01/09/2008

//==============================================================================
//   AliAnalysisGrid - Base grid utility class. Provides interface for creating
// a personalized JDL, finding and creating a dataset.
//==============================================================================

#include "TSystem.h"
#include "AliAnalysisGrid.h"

ClassImp(AliAnalysisGrid)

//______________________________________________________________________________
AliAnalysisGrid::AliAnalysisGrid(const AliAnalysisGrid& other)
                :TNamed(other), fSpecialBits(0)
{
// Copy ctor.
}

//______________________________________________________________________________
AliAnalysisGrid &AliAnalysisGrid::operator=(const AliAnalysisGrid& other)
{
// Assignment.
   if (this != &other) {
      TNamed::operator=(other);
      fSpecialBits = other.fSpecialBits;
   }
   return *this;
}

//______________________________________________________________________________
Bool_t AliAnalysisGrid::CreateToken(const char *username)
{
// Check if a valid token exists - if not create one
   TString user = gSystem->Getenv("USER");
   if (!user.Length()) {
      printf("Error in AliAnalysisGrid::CreateToken: $USER environment empty");
      return kFALSE;
   }
   Int_t err_msg = gSystem->Exec("no_command > /dev/null 2>/dev/null");
   Int_t token_value = gSystem->Exec("bash alien-token-info > /dev/null 2>/dev/null");
   if (token_value == err_msg) {
      printf("Error in AliAnalysisGrid::CreateToken: You do not seem to have <alien-token-info> in your path.");
      return kFALSE;
   }
   
   Bool_t to_create_token = kFALSE;
   if (!token_value) {
      // Token still valid, check alien_API_USER environment
      TString token_user = gSystem->Getenv("alien_API_USER");
      if (token_user.Length()) {
         // Environment file sourced
         if (!username) return kTRUE;                // for default $USER
         if (token_user == username) return kTRUE;   // for <username>
         // A valid token existing, and environment sourced, but for a different user
         to_create_token = kTRUE;
      }
   } else {
      // Token not valid anymore for <username>. Call alien-token-init   
      to_create_token = kTRUE;
   }
   if (to_create_token) {   
      printf("______________________________________________________________________________________\n");
      printf("AliAnalysisGrid::CreateToken: Seems you need a token. Calling alien-token-init for you\n");
      printf("______________________________________________________________________________________\n");
      Int_t token_init = 0;
      if (username) token_init = gSystem->Exec(Form("alien-token-init %s", username));
      else          token_init = gSystem->Exec("alien-token-init");
      if (token_init == err_msg) {
         printf("   Woops - semms alien-token-init is not in your path...\n");
         return kFALSE;
      } else if (token_init != 0) {
         printf("   Woops - did not succeed...\n");
         return kFALSE;
      }
   }
   // We have a valid token - just source it
   printf("______________________________________________________________________________________\n");
   printf("AliAnalysisGrid::CreateToken: Your token needs to be sourced in the current shell\n");
   printf("   USE:    > source /tmp/gclient_env_%d\n", gSystem->GetUid(user));
   printf("______________________________________________________________________________________\n");
   return kFALSE;
}

//______________________________________________________________________________
AliAnalysisGrid::EPluginRunMode AliAnalysisGrid::GetRunMode() const
{
// Get the current run mode.
   if (TObject::TestBit(kTest)) return AliAnalysisGrid::kTest;
   if (TObject::TestBit(kOffline)) return AliAnalysisGrid::kOffline;
   if (TObject::TestBit(kSubmit)) return AliAnalysisGrid::kSubmit;
   if (TObject::TestBit(kMerge)) return AliAnalysisGrid::kMerge;
   return AliAnalysisGrid::kFull;
}
   
//______________________________________________________________________________
void AliAnalysisGrid::SetRunMode(const char *mode)
{
// Set the alien plugin run mode. All modes require presence of a valid token
// and sourcing the AliEn environment. Supported modes are:
// - full (default): Generates requested datasets, locally generates the JDL,
//                   saves existing analysis manager to the file analysis.root,
//                   generates analysis macro, execution and validation scripts,
//                   copies all these files to AliEn working space and submits 
//                   the job leaving user in an AliEn shell.
// - test          : Generates only 10 entries of the first requested dataset and
//                   copies this locally as wn.xml, generates all files from the
//                   full run mode except the JDL and executes the analysis locally.
//                   This mode can be used to test if the analysis may run in grid.
// - offline       : No dataset is produced, but all other files are locally generated.
//                   No file is copied in AliEn workspace. This mode can be used to
//                   customize the automatic JDL/analysis macro.
// - submit        : Datasets are generated in AliEn but the JDL and all the other
//                   files are supposed to exist in the local directory. The files
//                   are copied to AliEn and the job is submitted. This mode should
//                   be used in correlation with "offline mode" to submit customized
//                   analysis macro/jdl.
// - merge         : Only MergeOutputs() method called to merge the registered
//                   outputs of a job that finished.
   TString smode(mode);
   smode.ToLower();
   TObject::SetBit(kTest, kFALSE);
   TObject::SetBit(kOffline, kFALSE);
   TObject::SetBit(kSubmit, kFALSE);
   TObject::SetBit(kMerge, kFALSE);
   if (smode.Contains("test")) {
      TObject::SetBit(kTest, kTRUE);
      return;
   }
   if (smode.Contains("offline")) {
      TObject::SetBit(kOffline, kTRUE);
      SetUseCopy(kFALSE);
      SetCheckCopy(kFALSE);
      return;
   }
   if (smode.Contains("submit")) {
      TObject::SetBit(kSubmit, kTRUE);
      return;
   }
   if (smode.Contains("merge") || smode.Contains("terminate")) {
      TObject::SetBit(kMerge, kTRUE);
      return;
   }
   if (!smode.Contains("full")) {
      Warning("SetRunMode","Run mode \"%s\" not known. Supported modes: \"full\", \"test\", \"offline\", \"submit\" and \"merge\"", mode);
      Warning("SetRunMode","Run mode set to FULL");
   }   
}
