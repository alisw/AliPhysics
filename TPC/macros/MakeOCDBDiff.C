/// \file MakeOCDBDiff.C
///
/// ~~~{.cpp}
/// .L  $ALICE_ROOT/TPC/macros//MakeOCDBDiff.C
/// ~~~

Bool_t  MakeOCDBDiff(const char *ocdb1, const char *ocdb2){
  /// Compare by by byte the content of the OCDB entry
  /// Input parameters:
  ///   ocdb1 - path to the OCDB file1
  ///   ocdb2 - path to the OCDB file2
  /// Return value:
  ///   kTRUE - in case the content of the OCDB object (persistent part) is exactly the same
  ///   kFALSE - othewise

  /* 
     ocdb1="Run188720_192738_v2_s0.root"
     ocdb1="Run0_188719_v2_s0.root"
     ocdb2="Run192739_999999999_v2_s0.root" 
  */
  TFile * f1 = TFile::Open(ocdb1);
  TFile * f2 = TFile::Open(ocdb2);
  {if (!f1 || !f2){
    printf("Problem 0:  Files not accessible (%s,%s)\n",ocdb1,ocdb2);
    return kFALSE;
    }}
  AliCDBEntry * entry1 = (AliCDBEntry*)f1->Get("AliCDBEntry");
  AliCDBEntry * entry2 = (AliCDBEntry*)f2->Get("AliCDBEntry");
  {if (!entry1||!entry2){
    printf("Problem 1:  OCDB entry not available (%s,%s)\n",ocdb1,ocdb2);
    return kFALSE; 
    }}
  TObject* object1=entry1->GetObject();
  TObject* object2=entry2->GetObject();
  TMessage * file1 = new TMessage(TBuffer::kWrite);
  file1->WriteObject(object1);
  Int_t size1=file1->Length();  
  TMessage * file2 = new TMessage(TBuffer::kWrite);
  file2->WriteObject(object2);
  Int_t size2=file2->Length(); 
  {if (size1!=size2){
    printf("Problem 2:  OCDB entry of different size (%d,%d)",size1,size2);
    return kFALSE;
    }}
  Int_t countDiff=0;
  for (Int_t i=0; i<size1; i++)    if (file1->Buffer()[i]!=file2->Buffer()[i]) countDiff++;
  {if (countDiff>0){
    printf("Objects are different. %d different bytes\n", countDiff );
    return kFALSE;
    }}
  printf("OCDB Objects are the same\n");

  TClass *cl1= object1->Class();
  
  
  return kTRUE;
}
