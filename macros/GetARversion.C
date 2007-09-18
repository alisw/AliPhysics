const char* GetARversion(){
//void GetARversion(){
  // Returns AliRoot version extracted from what is found in the
  // $ALICE_ROOT/CVS/ directory
  // 
  TString vAli;
  const char* vFile = gSystem->ExpandPathName("$ALICE_ROOT/CVS/Tag");
  if(gSystem->AccessPathName(vFile)){
    vAli="HEAD";
  }else{
    TFile *fv= TFile::Open("$ALICE_ROOT/CVS/Tag?filetype=raw","READ");
    Int_t size = fv->GetSize();
    char *buf = new Char_t[size];
    memset(buf, '\0', size);
    fv->Seek(0);
    if ( fv->ReadBuffer(buf, size) ) {
      Warning("GetARversion.C","Error reading AliRoot version from file to buffer!");
      vAli="";
    }
    vAli = buf;
    if(vAli.Contains('\n')) vAli.Remove(vAli.First('\n'));
    if(vAli.Contains('v')) vAli.Remove(0,vAli.First('v'));
  }
  delete vFile;
  return vAli.Data();
}

