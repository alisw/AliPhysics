#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

void TestOrbits(Int_t Fill)
{
   // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // set branches
  Int_t orbit;
  Int_t c2orbit;
  Int_t aqflag;
  g_vdm_Tree->SetBranchAddress("orbits",&orbit);
  g_vdm_Tree->SetBranchAddress("c2orbit",&c2orbit);  
  g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);

  // loop
  Int_t n_eq = 0;
  Int_t n_eq_aq = 0;
  Int_t n_aq = 0;  
  Int_t n = g_vdm_Tree->GetEntries();
  for (Int_t i=0;i<n;i++) {
    g_vdm_Tree->GetEntry(i);
    if (i<100) cout << orbit << " " << c2orbit << endl;
    if (orbit == c2orbit) n_eq++;
    if (aqflag > 0) {
      n_aq++;
      if (orbit == c2orbit) n_eq_aq++;
    }
  }

  // print out results
  cout << " Entries " << n << endl;
  cout << " Equal " << n_eq << endl;
  cout << " AQ " << n_aq << endl;
  cout << " AQ equal " << n_eq << endl;  
}

void TestBCnumbering(Int_t Fill)
{
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // get bunch crossing information
  // -- number of bc
  Int_t nIBC = GetNumberInteractingBunchCrossings();

  // -- bunch indices
  Int_t *bunches = new Int_t [nIBC];
  GetBunchIndices(bunches);

  // -- bucket info
  // Int_t *BucketA = new Int_t [nIBC];
  // Int_t *BucketC = new Int_t [nIBC];
  Int_t BucketA[nIBC];
  Int_t BucketC[nIBC];  
  GetBucketInfo(BucketA,BucketC);

  // print indices
  for(Int_t j=0;j<nIBC;j++) {
    // my computation
    Int_t idx1 = ((Int_t) (BucketA[j]/10.0));
    Int_t idx2 = ((Int_t) (BucketC[j]/10.0));
    // pietros
    Int_t ib = bunches[j];
    Int_t i1=(ib+3220)%3564;
    Int_t i2=(ib+547)%3564;
    if ( i1 != idx1 || i2 != idx2)
      cout << j << " idx1 " << idx1 << " = " << i1
	   << " idx2 " << idx2 << " = " << i2 << endl;
  }
  // cleanup
  delete [] bunches;
  // delete [] BucketA;
  //delete [] BucketC;
}

void Test(Int_t Fill)
{
  // TestOrbits(Fill);
  TestBCnumbering(Fill);
}
