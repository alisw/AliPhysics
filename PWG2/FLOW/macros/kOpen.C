//////////////////////////////////////////////////////////////////////
// - old way -
//////////////////////////////////////////////////////////////////////

bool kOpen(int evtN=0)
{
 TString fileName = "galice.root" ;
 AliRunLoader* rl = AliRunLoader::Open(fileName.Data(),"MyEvent","read");
 rl->LoadgAlice();
 AliRun* gAlice = rl->GetAliRun();
 rl->LoadHeader();
 rl->LoadKinematics();
 Int_t fNumberOfEvents = rl->GetNumberOfEvents() ;
 cout << " Found :  " << fNumberOfEvents << "  event(s) ... " << endl ; 

 Int_t exitStatus = rl->GetEvent(evtN) ; if(exitStatus!=0) { return kFALSE ; }

 TTree* pKTree = (TTree*)rl->TreeK();	      // Particles TTree (KineTree)
 AliStack* pStack = gAlice->Stack();	      // Particles Stack (use "Label()" to get the number in the stack)

 // else if(rl)      // opens files one by one (unload and reload)
 // {
 //  rl->UnloadgAlice() ;
 //  rl->UnloadHeader() ;
 //  rl->UnloadKinematics() ;
 //  delete rl ; rl = 0 ;
 // }

 Int_t fNumberOfParticles = pKTree->GetEntries() ;
 Int_t nPart = pStack->GetNtrack() ;
 cout << " Event n. " << evtN << "  contains :  " << fNumberOfParticles << "  particles in the TTree  ( =  " << nPart << "  in the stack ) . " << endl ; 

 return kTRUE ; 
}

