#ifndef __MAKECINT__
 #ifndef __CINT__
  #include "alles.h"
  #include "TInterpreter.h"
  #include "TBranchClones.h"
  #include "TBranchElement.h"  
  #include "AliTPCTrackHits.h"
  #include "AliTRDtrackHits.h"
 #endif
#endif

void ConvertToNewIO(const char* name)
{
//  AliLoader::SetDebug();
  AliConfig* conf = AliConfig::Instance();
  TClass* detclass = AliDetector::Class(); 
  void* buff;
  Bool_t skipit = kFALSE;
  void * pbuf[100]; 
  TBranch* branches[100];
  Int_t nbranches = 0;
  
  AliRunLoader* rl = AliRunLoader::Open("galiceNewIO.root","Event","recreate");
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(100);
  
  AliRun* outAliRun;

  if (gAlice == 0x0)
   {
     outAliRun = new AliRun("OutgAlice","Output gAlice");
   }
  else
   {
    outAliRun = gAlice;
   } 
  gAlice = 0x0;
   
  outAliRun->SetRunLoader(rl);
  
  TFile* infile = TFile::Open(name);
  
  if (infile == 0x0)
   {
     ::Error("ConvertToNewIO.C","Can not open input file %s",name);
     return;
   }
  
  AliRun* inAliRun = (AliRun*)infile->Get("gAlice");
  
  if(inAliRun == 0x0)
   {
     ::Error("ConvertToNewIO.C","Can not find gAlice in input file");
     return;
   }
  gAlice = inAliRun;
//  inAliRun->GetEvent(0);
  
  TObjArray* modules = inAliRun->Modules();
  if (modules == 0x0)
   {
     ::Error("ConvertToNewIO.C","Can not get array with modules from AliRun");
     return;
   }
  TIter next(modules);
  AliModule* module;
  TObject * object;
  while ((object = next()))  
   {
     module = dynamic_cast<AliModule*>(object);
     if (module == 0x0) {
       ::Error("ConvertToNewIO.C","Can not cast module %#x ",object);
       object->Dump();
       continue;
     }

     outAliRun->AddModule(module);
     
     TClass* modclass = module->IsA();
     AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
     if (det) 
      {
        ::Info("ConvertToNewIO.C"," Adding %s to RL.",det->GetName());
        conf->Add(det,"Event");
        rl->AddLoader(det);
      }
   }

  TParticle* particleBuffer = new TParticle();
  /***************************************************/
  /****          Event to Event         **************/
  /***************************************************/

  TTree* treeE = (TTree*)infile->Get("TE");
  if (treeE == 0x0)
   {
     ::Error("ConvertToNewIO.C","Can not get TreeE from AliRun");
     return;
   }
  rl->MakeTree("E");
    
  AliHeader* inheader = new AliHeader();
  treeE->SetBranchAddress("Header",&inheader);
  
  Int_t nevents = (Int_t)treeE->GetEntries();
  for (Int_t i = 0; i<nevents; i++)
   {
     ::Info("ConvertToNewIO.C","Converting Event %d.",i);
     rl->SetEventNumber(i);
     treeE->GetEvent(i);
     /*****************************************/
     /*            H E A D E R                */
     /*****************************************/
     
     ::Info("ConvertToNewIO.C","Copying Header");
     AliHeader* outheader = rl->GetHeader();
     
     
     outheader->SetEvent(inheader->GetEvent());
     outheader->SetNvertex(inheader->GetNvertex());
     outheader->SetNprimary(inheader->GetNprimary());
     outheader->SetNtrack(inheader->GetNtrack());
     outheader->SetRun(inheader->GetRun());
     outheader->SetEventNrInRun(inheader->GetEventNrInRun());
     outheader->SetStack(inheader->Stack());
     outheader->SetGenEventHeader(inheader->GenEventHeader());
     rl->TreeE()->Fill();

     /*****************************************/
     /*       K I N E M A T I C S             */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying Kinematics.");
     TString treeKname("TreeK");
     treeKname+=i;
     TTree* treeK = (TTree*)infile->Get(treeKname);
     if (treeK)
      { 
       if (treeK->GetEntries() > 0)
       {   
        //I do this gimnastics to set directory correctly, without changing NewIO code 
        //only for this purpose
        rl->LoadKinematics("update");
        rl->MakeTree("K");
        rl->GetEventFolder()->Remove(rl->TreeK());
        delete rl->TreeK();
        
        treeK->SetBranchAddress("Particles",&particleBuffer);
        ::Info("ConvertToNewIO.C","Cloning TreeK ...");
        TTree* tk = treeK->CloneTree();
        ::Info("ConvertToNewIO.C","Cloning TreeK ... Done");
        tk->SetName(AliRunLoader::fgkKineContainerName);
        rl->GetEventFolder()->Add(tk);
        rl->WriteKinematics("OVERWRITE");
        rl->UnloadKinematics();
       }
       else
        {
          Info("ConvertToNewIO.C","Kinematics Tree is Empty");
        }
      }
     else
      {
        Error("ConvertToNewIO.C","Could not find Kinematics tree named %s",treeKname.Data());
      }
     delete treeK;
     /*****************************************/
     /*   T R A C K   R E F E R E N C E S     */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying Track Refs.");
     TString treeTRname("TreeTR");
     treeTRname+=i;
     TTree* treeTR = (TTree*)infile->Get(treeTRname);
     if (treeTR)
      { 
       if (treeTR->GetEntries() > 0)
       {   
        next.Reset();
	TObject * object;
        while ((object = next()))  
         {
	   module = dynamic_cast<AliModule*>(object);
	   if (module == 0x0) {
	     ::Error("ConvertToNewIO.C","Can not cast module %#x ",object);
	     object->Dump();
	     continue;
	   }
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           if (det) 
            {
              TClonesArray* trbuffer = new TClonesArray("AliTrackReference", 100);
              treeTR->SetBranchAddress(det->GetName(),&trbuffer);
            }
         }
        
        //I do this gimnastics to set directory correctly, without changing NewIO code 
        //only for this purpose
        rl->LoadTrackRefs("update");
        rl->MakeTrackRefsContainer();
        rl->GetEventFolder()->Remove(rl->TreeTR());
        delete rl->TreeTR();
        
        ::Info("ConvertToNewIO.C","Cloning TreeTR ...");
        TTree* tr = treeTR->CloneTree();
        ::Info("ConvertToNewIO.C","Cloning TreeTR ... Done");
        
        tr->SetName(AliRunLoader::fgkTrackRefsContainerName);
        rl->GetEventFolder()->Add(tr);
        rl->WriteTrackRefs("OVERWRITE");   
        rl->UnloadTrackRefs();
       }
       else
        {
          Info("ConvertToNewIO.C","Track References Tree is Empty");
        }
      }
     else
      {
        Error("ConvertToNewIO.C","Could not find Track Refs tree named %s",treeTRname.Data());
      }
     delete treeTR;
      
     /*****************************************/
     /*          H I T S                      */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying Hits.");
     TString treeHname("TreeH");
     treeHname+=i;
     TTree* treeH = (TTree*)infile->Get(treeHname);
     
     if (treeH)
      { 
       if (treeH->GetEntries() > 0)
        {   
         TObjArray* lob = treeH->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
        
         while ((module = (AliModule*)nextnewmodule()))
          {
            TClass* modclass = module->IsA();
            AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
            TClonesArray* ca = 0;
            if (det) 
             {
               AliLoader* loader = det->GetLoader();
               if (loader == 0x0)
                {
                  ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
                  continue;
                }

               TString mask(det->GetName());
               
               loader->LoadHits("update");
               loader->MakeTree("H");
               loaders->Add(loader);
               for(Int_t b=0; b<lob->GetEntries();b++)
                {
                  TBranch* branch = (TBranch*)lob->At(b);
                  TString bname(branch->GetName());//
                  if ( bname.BeginsWith(det->GetName()) )
                   {
                     ::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
                     ::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
                     TString contname(branch->GetClassName());
                  
                     Int_t splitlvl = branch->GetSplitLevel();
                    // if (splitlvl) splitlvl = 99;
                 
                     if( contname.CompareTo("TClonesArray") == 0)
                      { 
                        TBranchElement* belem = (TBranchElement*)branch;
                        ::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
                    
                        ca = new TClonesArray(belem->GetClonesName());
                        pbuf[nbranches] = ca;
	    
                        branch->SetAddress(&(pbuf[nbranches]));
                        ::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
 
                        ::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
                        branches[nbranches] = loader->TreeH()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
                        nbranches++;
                      }
                     else
                      {
                        TClass* bcl = gROOT->GetClass(branch->GetClassName());
                        if (bcl == 0x0)
                         {
                           ::Error("ConvertToNewIO.C","Can not get TClass object of class named %s",branch->GetClassName());
                           continue;
                         }
                        pbuf[nbranches] = bcl->New();
                        ::Info("ConvertToNewIO.C","Dumping buffer:");
                        ((TObject*)pbuf[nbranches])->Dump();
                        ::Info("ConvertToNewIO.C","Setting Adress:");
                        branch->SetAddress(&(pbuf[nbranches]));
                        ::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
                        branches[nbranches] =loader->TreeH()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
                        nbranches++;
                      }
                   }
                }//loop over branches
             }//if module is detector
          }//while loop over modules
         Int_t nentr = (Int_t)treeH->GetEntries();
         ::Info("ConvertToNewIO.C","Copying Hits . Number of entries %d  ... ",nentr);
               
//         ::Info("ConvertToNewIO.C","Getting event:");
//         nentr = 100;
         Int_t nl = loaders->GetEntries();
         for (Int_t e = 0; e < nentr; e++)
           {
             printf("%d\r",e);
             treeH->GetEntry(e);
             
             for (Int_t l = 0; l < nbranches; l++)
              {
//              printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
                branches[l]->SetAddress(&(pbuf[l]));    
              }
     
             for (Int_t l = 0; l < nl; l++)
              {
                AliLoader* loader = (AliLoader*)loaders->At(l);
//              printf("Filling %s\n",loader->GetName());
                loader->TreeH()->Fill();
              }
             #ifndef __MAKECINT__
              #ifndef __CINT__
               fflush(0);
              #endif
             #endif
            }
        printf("\n");
               
        ::Info("ConvertToNewIO.C","Copying Hits ...  Done");
        for (Int_t l = 0; l < nl; l++)
         {
           AliLoader* loader = (AliLoader*)loaders->At(l);
           loader->WriteHits("OVERWRITE");
           loader->UnloadHits();
         }
        delete loaders;
        for (Int_t l = 0; l < nbranches; l++)
         {
           delete (TObject*)pbuf[l];
         }
         nbranches = 0; 	    
       }
       else //tree has any entries
        {
          Info("ConvertToNewIO.C","Hits Tree is Empty");
        }
      }
     else //treeH
      {
       ::Warning("ConvertToNewIO.C","Could not get TreeH from in AliRun.");
      }
     delete treeH; 

     /*****************************************/
     /*          S  D i g i t s               */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying S Digits.");
     TString treeSname("TreeS");
     treeSname+=i;
     TTree* treeS = (TTree*)infile->Get(treeSname);
     if (treeS)
      { 
       if (treeS->GetEntries() > 0)
        {   
         TObjArray* lob = treeS->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
        
         while ((module = (AliModule*)nextnewmodule()))
          {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           TClonesArray* ca = 0;
           if (det) 
            {
              AliLoader* loader = det->GetLoader();
              if (loader == 0x0)
               {
                 ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
                 continue;
               }

              TString mask(det->GetName());
              
              loader->LoadSDigits("update");
              loader->MakeTree("S");
              loaders->Add(loader);
              for(Int_t b=0; b<lob->GetEntries();b++)
               {
                 TBranch* branch = (TBranch*)lob->At(b);
                 TString bname(branch->GetName());//
                 if ( bname.BeginsWith(det->GetName()) )
                  {
                    
                    ::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
                    ::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
                    TString contname(branch->GetClassName());
                 
                    Int_t splitlvl = branch->GetSplitLevel();
                   // if (splitlvl) splitlvl = 99;
                 
                    if ( contname.CompareTo("TClonesArray") == 0)
                     {
                       TBranchElement* belem = (TBranchElement*)branch;
                       ::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
                   
                       ca = new TClonesArray(belem->GetClonesName());
                       pbuf[nbranches] = ca;
	   
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);

                       ::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
                       branches[nbranches] = loader->TreeS()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                    else
                     {
                       TClass* bcl = gROOT->GetClass(branch->GetClassName());
                        if (bcl == 0x0)
                         {
                           ::Error("ConvertToNewIO.C","Can not get TClass object of class named %s",branch->GetClassName());
                           continue;
                         }
                       pbuf[nbranches] = bcl->New();
                       ::Info("ConvertToNewIO.C","Dumping buffer:");
                       ((TObject*)pbuf[nbranches])->Dump();
                       ::Info("ConvertToNewIO.C","Setting Adress:");
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
                       branches[nbranches] =loader->TreeS()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                  }
               }//loop over branches
            }//if module is detector
         }//while loop over modules
        Int_t nentr = (Int_t)treeS->GetEntries();
        ::Info("ConvertToNewIO.C","Copying SDigits. Number of entries %d  ... ",nentr);
              
//        ::Info("ConvertToNewIO.C","Getting event:");
//        nentr = 100;
        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            printf("%d\r",e);
            treeS->GetEntry(e);
            
            for (Int_t l = 0; l < nbranches; l++)
             {
//             printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
//               printf("Filling %s\n",loader->GetName());
               loader->TreeS()->Fill();
             }
            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
           }
        printf("\n");
              
        ::Info("ConvertToNewIO.C","Copying SDigits ...  Done");
        for (Int_t l = 0; l < nl; l++)
         {
           AliLoader* loader = (AliLoader*)loaders->At(l);
           loader->WriteSDigits("OVERWRITE");
           loader->UnloadSDigits();
         }
        delete loaders;
        for (Int_t l = 0; l < nbranches; l++)
         {
           delete (TObject*)pbuf[l];
         }
        nbranches = 0; 	    
       }
      else //tree has any entries
       {
         ::Info("ConvertToNewIO.C","S Digits Tree is Empty");
       }
      }
     else //treeS
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeS from AliRun.");
      }
     delete treeS; 

     /*****************************************/
     /*          D i g i t s                  */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying Digits.");
     TString treeDname("TreeD");
     treeDname+=i;
     TTree* treeD = (TTree*)infile->Get(treeDname);
     if (treeD)
      { 
       if (treeD->GetEntries() > 0)
        {   
       
         TObjArray* lob = treeD->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
         
         while ((module = (AliModule*)nextnewmodule()))
          {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           TClonesArray* ca = 0;
           if (det) 
            {
              AliLoader* loader = det->GetLoader();
              if (loader == 0x0)
               {
                 ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
                 continue;
               }

              TString mask(det->GetName());
              
              loader->LoadDigits("update");
              loader->MakeTree("D");
              loaders->Add(loader);
              for(Int_t b=0; b<lob->GetEntries();b++)
               {
                 TBranch* branch = (TBranch*)lob->At(b);
                 TString bname(branch->GetName());//
                 if ( bname.BeginsWith(det->GetName()) )
                  {
                    
                    ::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
                    ::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
                    TString contname(branch->GetClassName());
                 
                    Int_t splitlvl = branch->GetSplitLevel();
                   // if (splitlvl) splitlvl = 99;
                 
                    if ( contname.CompareTo("TClonesArray") == 0)
                     {
                       TBranchElement* belem = (TBranchElement*)branch;
                       ::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
                   
                       ca = new TClonesArray(belem->GetClonesName());
                       pbuf[nbranches] = ca;
	   
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);

                       ::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
                       branches[nbranches] = loader->TreeD()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                    else
                     {
                       TClass* bcl = gROOT->GetClass(branch->GetClassName());
                       if (bcl == 0x0)
                        {
                          ::Error("ConvertToNewIO.C","Can not get TClass object of class named %s",branch->GetClassName());
                          continue;
                        }
                       pbuf[nbranches] = bcl->New();
                       ::Info("ConvertToNewIO.C","Dumping buffer:");
                       ((TObject*)pbuf[nbranches])->Dump();
                       ::Info("ConvertToNewIO.C","Setting Adress:");
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
                       branches[nbranches] =loader->TreeD()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                  }
               }//loop over branches
            }//if module is detector
         }//while loop over modules
        Int_t nentr = (Int_t)treeD->GetEntries();
        ::Info("ConvertToNewIO.C","Copying Digits. Number of entries %d  ... ",nentr);
              
//        ::Info("ConvertToNewIO.C","Getting event:");
//        nentr = 100;
        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            printf("%d\r",e);
            treeD->GetEntry(e);
            
            for (Int_t l = 0; l < nbranches; l++)
             {
//             printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
//               printf("Filling %s\n",loader->GetName());
               loader->TreeD()->Fill();
             }
            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
           }
       printf("\n");
              
       ::Info("ConvertToNewIO.C","Copying Digits ...  Done");
       for (Int_t l = 0; l < nl; l++)
        {
          AliLoader* loader = (AliLoader*)loaders->At(l);
          loader->WriteDigits("OVERWRITE");
          loader->UnloadDigits();
        }
       delete loaders;
      }
      else //tree has any entries
       {
         Info("ConvertToNewIO.C","S Digits Tree is Empty");
       }
      }
     else //treeD
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeD from in AliRun.");
      }
     delete treeD; 
     for (Int_t l = 0; l < nbranches; l++)
      {
       delete (TObject*)pbuf[l];
      }
      nbranches = 0; 	    

     /*****************************************/
     /*          R e c  P o i n t s           */
     /*****************************************/
     ::Info("ConvertToNewIO.C","Copying RecPoints.");
     TString treeRname("TreeR");
     treeRname+=i;
     TTree* treeR = (TTree*)infile->Get(treeRname);
     if (treeR)
      { 
       if (treeR->GetEntries() > 0)
        {   
         TObjArray* lob = treeR->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
        
         while ((module = (AliModule*)nextnewmodule()))
          {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           TClonesArray* ca = 0;
           if (det) 
            {
              AliLoader* loader = det->GetLoader();
              if (loader == 0x0)
               {
                 ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
                 continue;
               }

              TString mask(det->GetName());
              
              loader->LoadRecPoints("update");
              loader->MakeTree("R");
              loaders->Add(loader);
              for(Int_t b=0; b<lob->GetEntries();b++)
               {
                 TBranch* branch = (TBranch*)lob->At(b);
                 TString bname(branch->GetName());//
                 if ( bname.BeginsWith(det->GetName()) )
                  {
                    
                    ::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
                    ::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
                    TString contname(branch->GetClassName());
                 
                    Int_t splitlvl = branch->GetSplitLevel();
                   // if (splitlvl) splitlvl = 99;
                 
                    if ( contname.CompareTo("TClonesArray") == 0)
                     {
                       TBranchElement* belem = (TBranchElement*)branch;
                       ::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
                   
                       ca = new TClonesArray(belem->GetClonesName());
                       pbuf[nbranches] = ca;
	   
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);

                       ::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
                       branches[nbranches] = loader->TreeR()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                    else
                     {
                       TClass* bcl = gROOT->GetClass(branch->GetClassName());
                       if (bcl == 0x0)
                        {
                          ::Error("ConvertToNewIO.C","Can not get TClass object of class named %s",branch->GetClassName());
                          continue;
                        }
                       pbuf[nbranches] = bcl->New();
                       ::Info("ConvertToNewIO.C","Dumping buffer:");
                       ((TObject*)pbuf[nbranches])->Dump();
                       ::Info("ConvertToNewIO.C","Setting Adress:");
                       branch->SetAddress(&(pbuf[nbranches]));
                       ::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
                       branches[nbranches] =loader->TreeR()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
                       nbranches++;
                     }
                  }
               }//loop over branches
            }//if module is detector
         }//while loop over modules
        Int_t nentr = (Int_t)treeR->GetEntries();
        ::Info("ConvertToNewIO.C","Copying RecPoints. Number of entries %d  ... ",nentr);
              
//        ::Info("ConvertToNewIO.C","Getting event:");
//        nentr = 100;
        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            printf("%d\r",e);
            treeR->GetEntry(e);
            
            for (Int_t l = 0; l < nbranches; l++)
             {
//             printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
//               printf("Filling %s\n",loader->GetName());
               loader->TreeR()->Fill();
             }
            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
           }
        printf("\n");
              
        ::Info("ConvertToNewIO.C","Copying RecPoints ...  Done");
        for (Int_t l = 0; l < nl; l++)
         {
           AliLoader* loader = (AliLoader*)loaders->At(l);
           loader->WriteRecPoints("OVERWRITE");
           loader->UnloadRecPoints();
         }
        delete loaders;
        for (Int_t l = 0; l < nbranches; l++)
         {
           delete (TObject*)pbuf[l];//delete branches buffers 
         }
        nbranches = 0; 	    
       }
      else //tree has any entries
       {
         Info("ConvertToNewIO.C","Rec Points Tree is Empty");
       }
      }
     else //treeR
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeR from in AliRun.");
      }
     delete treeR; 

    }//end of loop over events
    
  /***************************************************/
  /****             Run to Run           *************/
  /***************************************************/

  rl->WriteHeader("OVERWRITE");

  infile->cd();
  rl->CdGAFile();
  
  rl->WriteAliRun();
  rl->WriteRunLoader();

  ::Info("ConvertToNewIO.C","Done");
  
}

