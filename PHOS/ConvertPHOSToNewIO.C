#ifndef __MAKECINT__
 #ifndef __CINT__
#include <TROOT.h>
#include <TRint.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBS.h>
#include <TObjectTable.h>
#include <iostream.h>
#include <fstream.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TText.h>
#include <TTree.h>
#include <TBranch.h>
#include "AliModule.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliDetector.h"
#include "AliConfig.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include <TParticle.h>
#include "TBranchClones.h"
#include "TBranchElement.h"  

#include "PHOS/AliPHOSClusterizerv1.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "PHOS/AliPHOSSDigitizer.h"
#include "PHOS/AliPHOSTrackSegmentMakerv1.h"
#include "PHOS/AliPHOSPIDv1.h"

 #endif
#endif

void ConvertPHOSToNewIO(const char* name = "galice.root", TString branch ="Event" )
{
//  AliLoader::SetDebug();
  AliConfig* conf = AliConfig::Instance();
  TClass* detclass = AliDetector::Class(); 
  void * pbuf[100]; 
  TBranch* branches[100];
  Int_t nbranches = 0;
//  AliLoader::SetDebug(5); 
  AliRunLoader* rl = AliRunLoader::Open("galiceNewIO.root",branch,"recreate");
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);

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
  TString fSdig = "PHOS.SDigits."+branch+".root";
  TString fDig = "PHOS.Digits."+branch+".root";
  TString fRec = "PHOS.RecData."+branch+".root";

  TFile* insdfile  = TFile::Open(fSdig);
  TFile* indigfile = TFile::Open(fDig);
  TFile* inrecfile = TFile::Open(fRec);
  TFile* infile    = TFile::Open(name);

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
  
  TObjArray* modules = inAliRun->Modules();
  if (modules == 0x0)
   {
     ::Error("ConvertToNewIO.C","Can not get array with modules from AliRun");
     return;
   }
  TIter next(modules);
  AliModule* module;
  while ((module = (AliModule*)next()))  
   {
     outAliRun->AddModule(module);
     
     TClass* modclass = module->IsA();
     AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
     if (det) 
      {
        ////::Info("ConvertToNewIO.C"," Adding %s to RL.",det->GetName());
        conf->Add(det,branch);
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
     
     ////::Info("ConvertToNewIO.C","Copying Header");
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
     rl->SetNumberOfEventsPerFile(1000);
     /*****************************************/
     /*       K I N E M A T I C S             */
     /*****************************************/
     ////::Info("ConvertToNewIO.C","Copying Kinematics.");
     TString treeKname("TreeK");
     treeKname+=i;
     TTree* treeK = (TTree*)infile->Get(treeKname);
     if (treeK)
      { 
       if (treeK->GetEntries() > 0)
       {   
        //I do this gimnastics to set directory correctly, without changing NewIO code 
        //only for this purpose
        rl->SetNumberOfEventsPerFile(1000);
        rl->LoadKinematics("update");
        rl->MakeTree("K");
        rl->GetEventFolder()->Remove(rl->TreeK());
        delete rl->TreeK();
        
        treeK->SetBranchAddress("Particles",&particleBuffer);
        ////::Info("ConvertToNewIO.C","Cloning TreeK ...");
        TTree* tk = treeK->CloneTree();
        ////::Info("ConvertToNewIO.C","Cloning TreeK ... Done");
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
     ////::Info("ConvertToNewIO.C","Copying Track Refs.");
     TString treeTRname("TreeTR");
     treeTRname+=i;
     TTree* treeTR = (TTree*)infile->Get(treeTRname);
     if (treeTR)
      { 
       if (treeTR->GetEntries() > 0)
       {   
        next.Reset();
        while ((module = (AliModule*)next()))  
         {
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
	rl->SetNumberOfEventsPerFile(1000);
        rl->LoadTrackRefs("update");
        rl->MakeTrackRefsContainer();
        rl->GetEventFolder()->Remove(rl->TreeTR());
        delete rl->TreeTR();
        
        ////::Info("ConvertToNewIO.C","Cloning TreeTR ...");
        TTree* tr = treeTR->CloneTree();
        ////::Info("ConvertToNewIO.C","Cloning TreeTR ... Done");
        
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
     ////::Info("ConvertToNewIO.C","Copying Hits.");
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
            //TClonesArray* ca = 0;
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
                     ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
                     ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
                     TString contname(branch->GetClassName());
                  
                     Int_t splitlvl = branch->GetSplitLevel();
                    // if (splitlvl) splitlvl = 99;
                 
                     if( contname.CompareTo("TClonesArray") == 0)
                      { 
                        TBranchElement* belem = (TBranchElement*)branch;
                        //::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
                    
                        TClonesArray *ca = new TClonesArray(belem->GetClonesName());
                        pbuf[nbranches] = ca;
	    
                        branch->SetAddress(&(pbuf[nbranches]));
                        //::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
 
                        //::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
                        branches[nbranches] = loader->TreeH()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
                        nbranches++;
                      }
                     else
                      {
                        //::Info("ConvertToNewIO.C","Class Nmme is %s",branch->GetClassName());
                        TClass* bcl = gROOT->GetClass(branch->GetClassName());
                        if (bcl == 0x0)
                         {
                           ::Error("ConvertToNewIO.C","Can not get TClass object of class named %s",branch->GetClassName());
                           continue;
                         }
                        pbuf[nbranches] = bcl->New();
                        //::Info("ConvertToNewIO.C","Dumping buffer:");
                        //((TObject*)pbuf[nbranches])->Dump();
                        //::Info("ConvertToNewIO.C","Setting Adress:");
                        branch->SetAddress(&(pbuf[nbranches]));
                        //::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
                        branches[nbranches] =loader->TreeH()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
                        nbranches++;
                      }
                   }
                }//loop over branches
             }//if module is detector
          }//while loop over modules
         Int_t nentr = (Int_t)treeH->GetEntries();
         //::Info("ConvertToNewIO.C","Copying Hits . Number of entries %d  ... ",nentr);
               

         Int_t nl = loaders->GetEntries();
         for (Int_t e = 0; e < nentr; e++)
           {
             //printf("%d\n",e);
             treeH->GetEntry(e);
             
             for (Int_t l = 0; l < nbranches; l++)
              {
                //printf("Branch %d addr %#x\n",l,pbuf[l]);
                //printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
                branches[l]->SetAddress(&(pbuf[l]));    
              }
     
             for (Int_t l = 0; l < nl; l++)
              {
                AliLoader* loader = (AliLoader*)loaders->At(l);
                //printf("Filling %s\n",loader->GetName());
                loader->TreeH()->Fill();
              }
             #ifndef __MAKECINT__
              #ifndef __CINT__
               fflush(0);
              #endif
             #endif
            }
	 ////printf("\n");
               
        ////::Info("ConvertToNewIO.C","Copying Hits ...  Done");
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
     //::Info("ConvertToNewIO.C","Copying S Digits.\n\n\n");
     TString treeSname("TreeS");
     treeSname+=i;
     
     TTree* treeS = (TTree*)insdfile->Get(treeSname);
     if (treeS)
       { 
         TObjArray* lob = treeS->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
	 
         while ((module = (AliModule*)nextnewmodule()))
	   {
	     TClass* modclass = module->IsA();
	     AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
	     //TClonesArray* ca = 0;
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
			 
			 ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
			 ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
			 TString contname(branch->GetClassName());
			 
			 Int_t splitlvl = branch->GetSplitLevel();
			 // if (splitlvl) splitlvl = 99;
			 
			 if ( contname.CompareTo("TClonesArray") == 0)
			   {
			     TBranchElement* belem = (TBranchElement*)branch;
			     //::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
			     
			     TClonesArray * ca = new TClonesArray(belem->GetClonesName());
			     pbuf[nbranches] = ca;
			     
			     branch->SetAddress(&(pbuf[nbranches]));
			     //::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
			     
			     //::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
			     branches[nbranches] = loader->TreeS()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
			     nbranches++;
			   }
			 else
			   {
			     TClass* bcl = gROOT->GetClass(branch->GetClassName());
			     pbuf[nbranches] = bcl->New();
			     //::Info("ConvertToNewIO.C","Dumping buffer:");
			     //((TObject*)pbuf[nbranches])->Dump();
			     //::Info("ConvertToNewIO.C","Setting Adress:");
			     branch->SetAddress(&(pbuf[nbranches]));
			     //::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
			     branches[nbranches] =loader->TreeS()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
			     nbranches++;
			   }
		       }
		   }//loop over branches
            }//if module is detector
         }//while loop over modules
        TBranch* bbb = treeS->GetBranch("PHOS");
        Int_t nentr = (Int_t)bbb->GetEntries();
        ////::Info("ConvertToNewIO.C","Copying SDigits. Number of entries in PHSO branch %d  ... ",nentr);
             
        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            ////printf("%d\r",e);
//            treeS->GetEntry(e);
            bbb->GetEntry(e);

            for (Int_t l = 0; l < nbranches; l++)
             {
//             //printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
//               //printf("Filling %s\n",loader->GetName());
               loader->TreeS()->Fill();
             }
            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
           }
        ////printf("\n");
              
        ////::Info("ConvertToNewIO.C","Copying SDigits ...  Done");
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
     else //treeS
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeS from in AliRun.");
      }
     delete treeS; 

     /*****************************************/
     /*          D i g i t s                  */
     /*****************************************/
     //::Info("ConvertToNewIO.C","Copying Digits.\n\n\n");
     TString treeDname("TreeD");
     treeDname+=i;
     
     TTree* treeD = (TTree*)indigfile->Get(treeDname);
     if (treeD)
       { 
         TObjArray* lob = treeD->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
         
         while ((module = (AliModule*)nextnewmodule()))
	   {
	     TClass* modclass = module->IsA();
	     AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
	     //TClonesArray* ca = 0;
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
		     
		     ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
		     ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
		     TString contname(branch->GetClassName());
		     
		     Int_t splitlvl = branch->GetSplitLevel();
		     // if (splitlvl) splitlvl = 99;
                 
		     if ( contname.CompareTo("TClonesArray") == 0)
		       {
			 TBranchElement* belem = (TBranchElement*)branch;
			 ////::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
			 
			 TClonesArray * ca = new TClonesArray(belem->GetClonesName());
			 pbuf[nbranches] = ca;
			 
			 branch->SetAddress(&(pbuf[nbranches]));
			 ////::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
			 
			 ////::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
			 branches[nbranches] = loader->TreeD()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
			 nbranches++;
		       }
		     else
		       {
			 TClass* bcl = gROOT->GetClass(branch->GetClassName());
			 pbuf[nbranches] = bcl->New();
			 ////::Info("ConvertToNewIO.C","Dumping buffer:");
			 //((TObject*)pbuf[nbranches])->Dump();
			 ////::Info("ConvertToNewIO.C","Setting Adress:");
			 branch->SetAddress(&(pbuf[nbranches]));
			 ////::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
			 branches[nbranches] =loader->TreeD()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
			 nbranches++;
		       }
		   }
		   }//loop over branches
            }//if module is detector
         }//while loop over modules
	 TBranch* bbb = treeD->GetBranch("PHOS");
	 Int_t nentr = (Int_t)bbb->GetEntries();
	 ////::Info("ConvertToNewIO.C","Copying Digits. Number of entries %d  ... ",nentr);
              
        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            ////printf("%d\r",e);
            //treeD->GetEntry(e);
	    bbb->GetEntry(e);

            for (Int_t l = 0; l < nbranches; l++)
             {
//             //printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
//               //printf("Filling %s\n",loader->GetName());
               loader->TreeD()->Fill();
             }
            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
           }
	////printf("\n");
              
       ////::Info("ConvertToNewIO.C","Copying Digits ...  Done");
       for (Int_t l = 0; l < nl; l++)
        {
          AliLoader* loader = (AliLoader*)loaders->At(l);
          loader->WriteDigits("OVERWRITE");
          loader->UnloadDigits();
        }
       delete loaders;
       for (Int_t l = 0; l < nbranches; l++)
	 {
	   delete (TObject*)pbuf[l];
	 }
      nbranches = 0;
       }  
     else //treeD
       {
	 ::Warning("ConvertToNewIO.C","Could not get TreeD from in AliRun.");
       }
     delete treeD; 
     	    

     /*****************************************/
     /*          R e c  P o i n t s           */
     /*****************************************/
     ////::Info("ConvertToNewIO.C","Copying RecPoints.");
     TString treeRname("TreeR");
     treeRname+=i;

     TTree* treeR = (TTree*)inrecfile->Get(treeRname);
     if (treeR)
       {    
         TObjArray* lob = treeR->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
	 
         while ((module = (AliModule*)nextnewmodule()))
	   {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
          
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
		   if ( bname.Contains(det->GetName()) )
		     {
		       if(bname.Contains("Emc")||bname.Contains("Cpv"))
			 {
			   ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
			   ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
			   TString contname(branch->GetClassName());
			   
			   Int_t splitlvl = branch->GetSplitLevel();
			   // if (splitlvl) splitlvl = 99;
			   
			   if ( contname.CompareTo("TClonesArray") == 0)
			     {
			       TBranchElement* belem = (TBranchElement*)branch;
			       ////::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
			       
			       TClonesArray * ca = new TClonesArray(belem->GetClonesName());
			       pbuf[nbranches] = ca;
			       
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
			       
			       ////::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
			       branches[nbranches] = loader->TreeR()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			   else
			     {
			       TClass* bcl = gROOT->GetClass(branch->GetClassName());
			       pbuf[nbranches] = bcl->New();
			       ////::Info("ConvertToNewIO.C","Dumping buffer:");
			       //((TObject*)pbuf[nbranches])->Dump();
			       ////::Info("ConvertToNewIO.C","Setting Adress:");
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
			       branches[nbranches] =loader->TreeR()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			 }
		     }
		 }//loop over branches
	     }//if module is detector
	   }//while loop over modules
         TBranch* bbb = treeR->GetBranch("PHOSEmcRP");
	 Int_t nentr = (Int_t)bbb->GetEntries();
	 ////::Info("ConvertToNewIO.C","Copying RecPoints. Number of entries %d  ... ",nentr);
  
	
	 //Int_t nl = loaders->GetEntries();
	 //	  //printf(">>>>>>>>>>>>>>>>>>>>%d\n",nl);
	 bbb->SetAddress(&(pbuf[0]));

	 for (Int_t e = 0; e < nentr; e++)
	   {
	     ////printf("%d\r",e);
	    
	     bbb->GetEntry(e);
	     AliLoader* loader = (AliLoader*)loaders->At(0);
	     
	     TBranch* bbb = treeR->GetBranch("PHOSCpvRP");
	     Int_t nentr = (Int_t)bbb->GetEntries();
	     ////::Info("ConvertToNewIO.C","Copying RecPoints. Number of entries %d  ... ",nentr);
	     
	     bbb->SetAddress(&(pbuf[1]));
	    
	     for (Int_t e = 0; e < nentr; e++)
	       {
		 //		 //printf("%d\r",e);
		 bbb->GetEntry(e);
		 //		 AliLoader* loader = (AliLoader*)loaders->At(0);
	       }
	     ////printf("Filling %s\n",loader->GetName());
	     loader->TreeR()->Fill();

             #ifndef __MAKECINT__
              #ifndef __CINT__
              fflush(0);
             #endif
            #endif
	   }
	 ////printf("\n");


        ////::Info("ConvertToNewIO.C","Copying RecPoints ...  Done");
         Int_t nl = loaders->GetEntries();
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
     else //treeR
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeR from in AliRun.");
      }
     delete treeR; 
    
     /*****************************************/
     /*          Track Segments               */
     /*****************************************/
     ////::Info("ConvertToNewIO.C","Copying Tracks.");
     TString treeTname("TreeR");
     treeTname+=i;

     TTree* treeT = (TTree*)inrecfile->Get(treeTname);
     if (treeT)
       {    
         TObjArray* lob = treeT->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
	 
         while ((module = (AliModule*)nextnewmodule()))
	   {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           //TClonesArray* ca = 0;
           if (det) 
	     {
	       AliLoader* loader = det->GetLoader();
	       if (loader == 0x0)
		 {
		   ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
		   continue;
		 }
	       
	       TString mask(det->GetName());
	       
	       loader->LoadTracks("update");
	       loader->MakeTree("T");
	       loaders->Add(loader);
	       for(Int_t b=0; b<lob->GetEntries();b++)
		 {
		   TBranch* branch = (TBranch*)lob->At(b);
		   
		   TString bname(branch->GetName());//
		   if ( bname.Contains(det->GetName()) )
		     {
		       if(bname.Contains("TS"))
			 {
			   ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
			   ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
			   TString contname(branch->GetClassName());
			   
			   Int_t splitlvl = branch->GetSplitLevel();
			   // if (splitlvl) splitlvl = 99;
			   
			   if ( contname.CompareTo("TClonesArray") == 0)
			     {
			       TBranchElement* belem = (TBranchElement*)branch;
			       ////::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
			       
			       TClonesArray * ca = new TClonesArray(belem->GetClonesName());
			       pbuf[nbranches] = ca;
			       
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
			       
			       ////::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
			       branches[nbranches] = loader->TreeT()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			   else
			     {
			       TClass* bcl = gROOT->GetClass(branch->GetClassName());
			       pbuf[nbranches] = bcl->New();
			       ////::Info("ConvertToNewIO.C","Dumping buffer:");
			       //((TObject*)pbuf[nbranches])->Dump();
			       ////::Info("ConvertToNewIO.C","Setting Adress:");
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
			       branches[nbranches] =loader->TreeT()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			 }
		     }
		 }//loop over branches
	     }//if module is detector
	   }//while loop over modules
         TBranch* bbb = treeT->GetBranch("PHOSTS");
	 Int_t nentr = (Int_t)bbb->GetEntries();
	 ////::Info("ConvertToNewIO.C","Copying Tracks. Number of entries %d  ... ",nentr);
              

        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
	    // //printf("%d\r",e);
	    //bbb->GetEntry(e);
            //treeR->GetEntry(e);
            
            for (Int_t l = 0; l < nbranches; l++)
             {
	       //            //printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
	    bbb->GetEntry(e);
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
	       ////printf("Filling %s\n",loader->GetName());
	      loader->TreeT()->Fill();
             }
	    

            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
	  }
        ////printf("\n");
	
        ////::Info("ConvertToNewIO.C","Copying Tracks ...  Done");
        for (Int_t l = 0; l < nl; l++)
         {
           AliLoader* loader = (AliLoader*)loaders->At(l);
           loader->WriteTracks("OVERWRITE");
           loader->UnloadTracks();
         }
        delete loaders;
        for (Int_t l = 0; l < nbranches; l++)
         {
           delete (TObject*)pbuf[l];//delete branches buffers
         }
        nbranches = 0; 	    
      
      }
     else //treeT
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeR from in AliRun.");
      }
     delete treeT; 
  
    



     /*****************************************/
     /*          Rec Particles                */
     /*****************************************/
     ////::Info("ConvertToNewIO.C","Copying RecParticles.");
     TString treePname("TreeR");
     treePname+=i;

     TTree* treeP = (TTree*)inrecfile->Get(treeTname);
     if (treeP)
       {    
         TObjArray* lob = treeP->GetListOfBranches();
         TObjArray* loaders = new TObjArray();
         TIter nextnewmodule(outAliRun->Modules());
	 
         while ((module = (AliModule*)nextnewmodule()))
	   {
           TClass* modclass = module->IsA();
           AliDetector *det = (AliDetector*)(modclass->DynamicCast(detclass,module));
           //TClonesArray* ca = 0;
           if (det) 
	     {
	       AliLoader* loader = det->GetLoader();
	       if (loader == 0x0)
		 {
		   ::Error("ConvertToNewIO.C","Can not find loader from %s.",det->GetName());
		   continue;
		 }
	       
	       TString mask(det->GetName());
	       
	       loader->LoadRecParticles("update");
	       loader->MakeTree("P");
	       loaders->Add(loader);
	       for(Int_t b=0; b<lob->GetEntries();b++)
		 {
		   TBranch* branch = (TBranch*)lob->At(b);
		   
		   TString bname(branch->GetName());//
		   if ( bname.Contains(det->GetName()) )
		     {
		       if(bname.Contains("PHOSRP"))
			 {
			   ////::Info("ConvertToNewIO.C","Found branch %s.",branch->GetName());
			   ////::Info("ConvertToNewIO.C","Buffer Class Name %s.",branch->GetClassName());
			   TString contname(branch->GetClassName());
			   
			   Int_t splitlvl = branch->GetSplitLevel();
			   // if (splitlvl) splitlvl = 99;
			   
			   if ( contname.CompareTo("TClonesArray") == 0)
			     {
			       TBranchElement* belem = (TBranchElement*)branch;
			       ////::Info("ConvertToNewIO.C","Clones Class Name %s.",belem->GetClonesName());
			       
			       TClonesArray * ca = new TClonesArray(belem->GetClonesName());
			       pbuf[nbranches] = ca;
			       
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Buffer addrss %#x",pbuf[nbranches]);
			       
			       ////::Info("ConvertToNewIO.C","Creating branch for Clones SpliLvl = %d",splitlvl);
			       branches[nbranches] = loader->TreeP()->Branch(branch->GetName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			   else
			     {
			       TClass* bcl = gROOT->GetClass(branch->GetClassName());
			       pbuf[nbranches] = bcl->New();
			       ////::Info("ConvertToNewIO.C","Dumping buffer:");
			       //((TObject*)pbuf[nbranches])->Dump();
			       ////::Info("ConvertToNewIO.C","Setting Adress:");
			       branch->SetAddress(&(pbuf[nbranches]));
			       ////::Info("ConvertToNewIO.C","Creating branch SpliLvl = %d",splitlvl);
			       branches[nbranches] =loader->TreeP()->Branch(branch->GetName(),branch->GetClassName(),&(pbuf[nbranches]),4000,splitlvl);
			       nbranches++;
			     }
			 }
		     }
		 }//loop over branches
	     }//if module is detector
	   }//while loop over modules
         TBranch* bbb = treeP->GetBranch("PHOSRP");
	 Int_t nentr = (Int_t)bbb->GetEntries();
	 ////::Info("ConvertToNewIO.C","Copying RecParticles. Number of entries %d  ... ",nentr);
              

        Int_t nl = loaders->GetEntries();
        for (Int_t e = 0; e < nentr; e++)
          {
            ////printf("%d\r",e);
	    bbb->GetEntry(e);
            //treeR->GetEntry(e);
            
            for (Int_t l = 0; l < nbranches; l++)
             {
	       //            //printf("%s %#x \n", branches[l]->GetName(),pbuf[l]);
               branches[l]->SetAddress(&(pbuf[l]));    
             }
    
            for (Int_t l = 0; l < nl; l++)
             {
               AliLoader* loader = (AliLoader*)loaders->At(l);
	       ////printf("Filling %s\n",loader->GetName());
               loader->TreeP()->Fill();
             }
	    

            #ifndef __MAKECINT__
             #ifndef __CINT__
              fflush(0);
             #endif
            #endif
	  }
        ////printf("\n");
	
        ////::Info("ConvertToNewIO.C","Copying RecParticles ...  Done");
        for (Int_t l = 0; l < nl; l++)
         {
           AliLoader* loader = (AliLoader*)loaders->At(l);
           loader->WriteRecParticles("OVERWRITE");
           loader->UnloadRecParticles();
         }
        delete loaders;
        for (Int_t l = 0; l < nbranches; l++)
         {
           delete (TObject*)pbuf[l];//delete branches buffers 
         }
        nbranches = 0; 	    
      
      }
     else //treeP
      {
        ::Warning("ConvertToNewIO.C","Could not get TreeR from in AliRun.");
      }
     delete treeP; 



    }//end of loop over events
  /***************************************************/
  /****             Write Tasks          *************/
  /***************************************************/
  
  AliLoader * loader =  rl->GetLoader("PHOSLoader");
  /****             S Digits             *************/
  TTree * s = (TTree*)insdfile->Get("TreeS0");
  TBranch * bsd = s->GetBranch("AliPHOSSDigitizer");
  AliPHOSSDigitizer * sdig = 0;
  bsd->SetAddress(&sdig);
  bsd->GetEntry(0) ; 
  sdig->SetEventFolderName(sdig->GetName()) ; 
  sdig->SetName("PHOSSDigitizer");
  sdig->Print() ;  
  TFile fsd(loader->GetSDigitsFileName(), "update") ; 
  sdig->Write() ;
  fsd.Print() ; 
  fsd.Close() ; 
  delete s ;
  
  /****              Digits             *************/
  TTree * d = (TTree*)indigfile->Get("TreeD0");
  TBranch * bd = d->GetBranch("AliPHOSDigitizer");
  AliPHOSDigitizer * dig = 0 ;
  bd->SetAddress(&dig) ; 
  bd->GetEntry(0) ; 
  dig->SetEventFolderName(dig->GetName()) ; 
  dig->SetName("PHOSDigitizer");
  dig->Print() ;  
  //dig->Dump() ;  
  ////::Info("Digitizer","Print done");
  TFile fd(loader->GetDigitsFileName(), "update") ; 
  dig->Write() ;
  fd.Print() ; 
  fd.Close() ;
  delete d; 
 /****             Rec Data             *************/
  TTree * t = (TTree*)inrecfile->Get("TreeR0");
 /****             Clusterizer          *************/
  TBranch * bclu = t->GetBranch("AliPHOSClusterizer");
  AliPHOSClusterizerv1 * clu = 0;
  bclu->SetAddress(&clu);
  bclu->GetEntry(0) ;
  clu->SetEventFolderName(clu->GetName()) ; 
  clu->SetName("PHOSReconstructioner");
  clu->Print() ;  
  TFile fclu(loader->GetRecPointsFileName(), "update") ; 
  clu->Write() ;
  fclu.Print() ; 
  fclu.Close() ; 
  /****             TrackSegmentMaker      *************/
  TBranch * btr = t->GetBranch("AliPHOSTrackSegmentMaker");
  AliPHOSTrackSegmentMakerv1 * tra = 0 ;
  btr->SetAddress(&tra);
  btr->GetEntry(0) ;
  tra->SetEventFolderName(tra->GetName()) ; 
  tra->SetName("PHOSTracker");
  tra->Print() ;  
  TFile ftr(loader->GetTracksFileName(), "update") ; 
  tra->Write() ;
  ftr.Print() ; 
  ftr.Close() ; 

  /****             PID          *************/

  TBranch * bpid = t->GetBranch("AliPHOSPID");
  AliPHOSPIDv1 * pid = 0 ;
  bpid->SetAddress(&pid);
  bpid->GetEntry(0) ; 
  pid->SetEventFolderName(pid->GetName()) ; 
  pid->SetName("PHOSPIDTask");
  pid->Print() ;  
  TFile fpid(loader->GetRecParticlesFileName(), "update") ; 
  pid->Write() ;
  fpid.Print() ; 
  fpid.Close() ; 
  delete t;



  /***************************************************/
  /****             Run to Run           *************/
  /***************************************************/

  rl->WriteHeader("OVERWRITE");

  infile->cd();
  TGeometry* geo = inAliRun->GetGeometry();
  rl->CdGAFile();
  geo->Write();
  
  rl->WriteAliRun();
  rl->WriteRunLoader();

  ////::Info("ConvertToNewIO.C","Done");
  
}

