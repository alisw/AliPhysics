#include "AliFileUtilities.h"
#include <TSystem.h>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

// some simple unit tests for AliFileUtilities
// TODO: consider using a proper unit testing framework
// TODO: can we run these file tests in isolation
// TODO: we should use randomized file names (minizing risk to clash with true files)


void test1() {
  // assume we start from a clean directory
  gSystem->Exec("touch test_0.root");
  gSystem->Exec("touch test_1.root");
  gSystem->Exec("touch test_12.root");
  gSystem->Exec("mkdir test132134Ali");
  gSystem->Exec("touch test132134Ali/foo");
  assert(AliFileUtilities::CountLocalFiles("*.root") == 3);
  // the directory should also be counted
  assert(AliFileUtilities::CountLocalFiles("test*") == 4);
  AliFileUtilities::RemoveLocalFiles("test_?.root");
  assert(AliFileUtilities::CountLocalFiles("*.root") == 1);
  assert(AliFileUtilities::RemoveLocalFiles("test*") == 1);
  assert(AliFileUtilities::CountLocalFiles("test*") == 1);
  gSystem->Exec("rm -r -f test132134Ali");
}

void test2() {
  // test some more complicated stuff (with soft links)
  gSystem->Exec("touch test_0.root");
  gSystem->Exec("touch test_1.root");
  gSystem->Exec("touch test_12.root");
  gSystem->Exec("ln -s test_12.root test_2.root");
  assert(AliFileUtilities::CountLocalFiles("*.root") == 4);
  AliFileUtilities::RemoveLocalFiles("test_?.root");
  assert(AliFileUtilities::CountLocalFiles("*.root") == 1);
  gSystem->Exec("rm test_12.root");
}

void test3() {
  // check that globing works with directories
  // test some more complicated stuff (with soft links)
  gSystem->Exec("mkdir testDIR");
  gSystem->Exec("touch testDIR/test_0.root");
  gSystem->Exec("touch testDIR/test_1.root");
  gSystem->Exec("touch testDIR/test_12.root");
  assert(AliFileUtilities::CountLocalFiles("testDIR/*.root") == 3);
  assert(AliFileUtilities::RemoveLocalFiles("testDIR/*.root") == 3);
  assert(AliFileUtilities::CountLocalFiles("testDIR/*.root") == 0);
  assert(AliFileUtilities::CountLocalFiles("testDIR") == 1);
  // can remove empty dir
  assert(AliFileUtilities::RemoveLocalFiles("testDIR") == 1);
  assert(AliFileUtilities::CountLocalFiles("testDIR") == 0);
  // remaining cleanup
  gSystem->Exec("rm -r -f testDIR");
}

void test4() {
  gSystem->Exec("mkdir testDIR");
  gSystem->Exec("mkdir testDIR/testDIR2");
  gSystem->Exec("mkdir testDIRLINKED");
  gSystem->Exec("mkdir testDIR/testDIR2/foo");
  gSystem->Exec("touch testDIR/test_0.root");
  gSystem->Exec("touch testDIRLINKED/hello");
  gSystem->Exec("touch testDIR/test_1.root");
  gSystem->Exec("touch testDIR/test_12.root");
  gSystem->Exec("ln -s testDIRLINKED testDIR/BAR");
  assert(AliFileUtilities::Remove_All("testDIR") == 7);
  // cleanup
  gSystem->Exec("rm -r -f testDIRLINKED");
}

int main() {
  test1();
  test2();
  test3();
  test4();
  return 0;
}
