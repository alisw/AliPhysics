//$Id$

void runMakePileup(Char_t *path,Int_t startev,Int_t endev,Char_t *pileupdir="pileup",Int_t npiles=25,Int_t triggerevent=-1,Char_t *gfile=0)
{
  gSystem->Load("MakePileup_C.so");

  MakePileup(path,startev,endev,pileupdir,npiles,triggerevent,gfile);
}
