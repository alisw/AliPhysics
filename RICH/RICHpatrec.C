RICHpatrec (Int_t evNumber1=0,Int_t evNumber2=0)
{
  if(gAlice) delete gAlice;      
  AliRunLoader *al;
  AliRun *a;
  
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update")))
    Fatal("ReadAlice","galice.root broken");
  
  al->LoadgAlice();
  if(!gAlice) Fatal("ReadAlice","No gAlice in file");
  a=al->GetAliRun();
  
  AliRICHPatRec *detect = new AliRICHPatRec("RICH patrec algorithm","Reconstruction");
    

    for (int nev=0; nev<= evNumber2; nev++) {    // Event Loop
      Int_t nparticles = al->GetEvent(nev);
      cout <<endl<< "Processing event:" <<nev<<endl;
      cout << "Particles       :" <<nparticles<<endl;
      if (nev < evNumber1) continue;
      if (nparticles <= 0) return;
      detect->PatRec();
   } // event loop  
   
    printf("\nEnd of Macro*************************************\n");
}



