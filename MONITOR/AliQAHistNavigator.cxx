/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//
//  support class for the QA histogram viewer
//
//  origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include "AliQAHistNavigator.h"

ClassImp(AliQAHistNavigator)

//_________________________________________________________________________
AliQAHistNavigator::AliQAHistNavigator(Int_t run):
    fPFile( NULL ),
    fPCORRFile( NULL ),
    fPQAResultFile( NULL ),
    fRun( run ),
    fPCurrFile( NULL ),
    fPCurrDetector( NULL ),
    fPCurrLevel( NULL ),
    fPCurrItem( NULL ),
    fPListOfFiles( new AliQADirList() ),
    fLoopAllFiles(kTRUE),
    fLoopAllDetectors(kTRUE),
    fLoopAllLevels(kTRUE),
    fInitOK(kFALSE),
    fExpertMode(kFALSE),
    fExpertDirName("Expert"),
    fPEmptyList(new TList())
{
    fPEmptyList->AddLast(new AliQADirListItem(""));
    ReReadFiles();
}

//_________________________________________________________________________
AliQAHistNavigator::~AliQAHistNavigator()
{
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetHistogram(TH1*& hist)
{
    TString file = GetFileName();
    TString dir = GetDirName();
    TString histname = GetItemName();
    cout<<"GetHistogram: "<<file<<":"<<dir<<"/"<<histname<<endl;
    if (file==""||dir==""||histname=="") 
    {
        printf("GetItem: nothing to fetch...\n");
        return kFALSE;
    }
    if (!OpenCurrentDirectory()) return kFALSE;
    hist = dynamic_cast<TH1*>( gDirectory->FindKey(histname)->ReadObj() );
    if (!hist)
    {
        printf("GetItem: null pointer returned by gDirectory\n");
        return kFALSE;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Next()
{
    if (!fPCurrFile||!fPCurrDetector||!fPCurrLevel) return kFALSE;
    if (!(fPCurrItem=(AliQADirListItem*)GetItemList()->After(fPCurrItem)))
    {
        if (!(fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->After(fPCurrLevel)))
        {
            if (!(fPCurrDetector=(AliQADirList*)fPCurrFile->GetDirs()->After(fPCurrDetector)))
            {
                if (!(fPCurrFile=(AliQADirList*)fPListOfFiles->GetDirs()->After(fPCurrFile)))
                {
                    //we're at the end of everything
                    if (fLoopAllFiles)
                    {
                        //rewind to the beginning
                        fPCurrFile = (AliQADirList*)fPListOfFiles->GetDirs()->First();
                        fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->First();
                        fPCurrLevel = (AliQADirList*) fPCurrDetector->GetDirs()->First();
                        fPCurrItem = (AliQADirListItem*)GetItemList()->First();
                        OpenCurrentFile();
                        OpenCurrentDirectory();
                        printf("----------------back at the beginning!\n");
                    } else return kFALSE; //no rewind, we finish
                } else //if there is a next file
                {
                    fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->First();
                    fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->First();
                    fPCurrItem=(AliQADirListItem*)GetItemList()->First();
                    OpenCurrentFile();
                    OpenCurrentDirectory();
                }
            } else //if there is a next detector
            {
                fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->First();
                fPCurrItem=(AliQADirListItem*)GetItemList()->First();
                OpenCurrentDirectory();
            }
        } else //if there is a next level
        {
            fPCurrItem=(AliQADirListItem*)GetItemList()->First();
            OpenCurrentDirectory();
        }
    } else //if there is a next histgram
    {
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Prev()
{
    if (!fPCurrLevel||!fPCurrDetector||!fPCurrFile) return kFALSE;
    if (!(fPCurrItem=(AliQADirListItem*)GetItemList()->Before(fPCurrItem)))
    {
        if (!(fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->Before(fPCurrLevel)))
        {
            if (!(fPCurrDetector=(AliQADirList*)fPCurrFile->GetDirs()->Before(fPCurrDetector)))
            {
                if (!(fPCurrFile=(AliQADirList*)fPListOfFiles->GetDirs()->Before(fPCurrFile)))
                {
                    //we're at the end of everything
                    if (fLoopAllFiles)
                    {
                        //rewind to the beginning
                        fPCurrFile = (AliQADirList*)fPListOfFiles->GetDirs()->Last();
                        fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->Last();
                        fPCurrLevel = (AliQADirList*) fPCurrDetector->GetDirs()->Last();
                        fPCurrItem = (AliQADirListItem*)GetItemList()->Last();
                        OpenCurrentFile();
                        OpenCurrentDirectory();
                        printf("----------------back at the end!\n");
                    } else return kFALSE; //no rewind, we finish
                } else //if there is a next file
                {
                    fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->Last();
                    fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->Last();
                    fPCurrItem=(AliQADirListItem*)GetItemList()->Last();
                    OpenCurrentFile();
                    OpenCurrentDirectory();
                }
            } else //if there is a next detector
            {
                fPCurrLevel=(AliQADirList*)fPCurrDetector->GetDirs()->Last();
                fPCurrItem=(AliQADirListItem*)GetItemList()->Last();
                OpenCurrentDirectory();
            }
        } else //if there is a next level
        {
            fPCurrItem=(AliQADirListItem*)GetItemList()->Last();
            OpenCurrentDirectory();
        }
    } else //if there is a next histgram
    {
    }
    return kTRUE;
}

//_________________________________________________________________________
void AliQAHistNavigator::SetExpertMode(Bool_t mode)
{
    //sets the expert mode
    Bool_t oldmode = fExpertMode;
    fExpertMode = mode;
    if (fExpertMode!=oldmode) fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::OpenCurrentFile()
{
    //open current file
    if (fPFile) fPFile->Close();
    if (!(fPFile->Open(GetFileName(),"READ")))
    {
        return kFALSE;
        cout<<"There is no file: "<<GetFileName()<<endl;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::OpenCurrentDirectory()
{
    //Open current directory
    if (!gDirectory->cd(GetDirName())) return kFALSE;
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetFile( TString file )
{
    //Set a new file to read from
    AliQADirList* tmp = (AliQADirList*)fPListOfFiles->GetDirs()->FindObject ( file.Data() );
    if (!tmp) return kFALSE;
    fPCurrFile = tmp;
    OpenCurrentFile();
    fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->First();
    fPCurrLevel = (AliQADirList*)fPCurrDetector->GetDirs()->First();
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetFile( Int_t file )
{
    //Set a new file to read from
    printf("AliQAHistNavigator::SetFile(%i)\n",file);
    AliQADirList* tmp = (AliQADirList*)fPListOfFiles->GetDirs()->At(file);
    if (!tmp) return kFALSE;
    fPCurrFile = tmp;
    OpenCurrentFile();
    fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->First();
    fPCurrLevel = (AliQADirList*)fPCurrDetector->GetDirs()->First();
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetDetector( TString det )
{
    //Set a new detector
    AliQADirList* tmp = (AliQADirList*)fPCurrFile->GetDirs()->FindObject( det.Data() );
    if (!tmp) return kFALSE;
    fPCurrDetector = tmp;
    fPCurrLevel = (AliQADirList*)fPCurrDetector->GetDirs()->First();
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetDetector( Int_t det )
{
    //Set a new detector
    printf("AliQAHistNavigator::SetDetector(%i)\n",det);
    AliQADirList* tmp = (AliQADirList*)fPCurrFile->GetDirs()->At( det );
    if (!tmp) return kFALSE;
    fPCurrDetector = tmp;
    fPCurrLevel = (AliQADirList*)fPCurrDetector->GetDirs()->First();
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetLevel( TString level )
{
    //Set a new level
    AliQADirList* tmp = (AliQADirList*)fPCurrDetector->GetDirs()->FindObject( level.Data() );
    if (!tmp) return kFALSE;
    fPCurrLevel = tmp;
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetLevel( Int_t level )
{
    //Set a new level
    AliQADirList* tmp = (AliQADirList*)fPCurrDetector->GetDirs()->At( level );
    if (!tmp) return kFALSE;
    fPCurrLevel = tmp;
    fPCurrItem = (AliQADirListItem*)GetItemList()->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetItem( TString hist )
{
    //set the new item
    AliQADirListItem* tmp = (AliQADirListItem*)GetItemList()->FindObject( hist.Data() );
    if (!tmp) return kFALSE;
    fPCurrItem = tmp;
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetItem( Int_t hist )
{
    //set the new item
    AliQADirListItem* tmp = (AliQADirListItem*)GetItemList()->At( hist );
    if (!tmp) return kFALSE;
    fPCurrItem = tmp;
    return kTRUE;
}

//_________________________________________________________________________
TList* AliQAHistNavigator::GetItemList()
{
    //returns the current list of histograms, if none, returns empty list
    TList* itemlist=NULL;
    if (fExpertMode)
    {
        AliQADirList* expertlist = (AliQADirList*)fPCurrLevel->GetDirs()->FindObject(fExpertDirName);
        if (expertlist) itemlist = expertlist->GetItems();
        else
        {
            //this is a bit of a hack, but it will always return something sensible
            AliQADirListItem* it = (AliQADirListItem*)fPEmptyList->First();
            it->SetParent(fPCurrLevel);
            itemlist = fPEmptyList;
        }
    } else
    {
        itemlist = fPCurrLevel->GetItems();
    }
    return itemlist;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetDetectorName()
{
    //Get name of current detector
    if (!fPCurrDetector) return "";
    TString name = fPCurrDetector->GetName();
    return name;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetLevelName()
{
    //Get name of current leve
    if (!fPCurrLevel) return "";
    TString name = fPCurrLevel->GetName();
    return name;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetFileName()
{
    //Get name of current file
    if (!fPCurrFile) return "";
    TString file = fPCurrFile->GetName();
    return file;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetDirName()
{
    //get the name of dir containing current item
    if (!fPCurrItem) return "";
    AliQADirList* d=(AliQADirList*)fPCurrItem->GetParent();
    TString path;
    do
    {
        path = d->GetName() + path;
        path = "/" + path;
        d=d->GetParent();
    }
    while (d->GetParent());
    return path;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetPath(AliQADirListItem* item)
{
    //Get the full path to teh directory containing item
    AliQADirList* d=item->GetParent();
    TString path = "";//item->GetString();
    do
    {
        TString sep = (d->GetParent()) ? "/" : ":/" ;
        path = d->GetName() + sep + path;
    }
    while (d=d->GetParent());
    return path;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetItemName()
{
    //Get name of current item
    if (!fPCurrItem) return "";
    return fPCurrItem->GetString();
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetListOfFiles()
{
    //scan directory for QA files
    delete fPListOfFiles;
    fPListOfFiles = new AliQADirList();

    TString CORRFileName;
    TString QAResultFileName;

    TString macdir(".");
    gSystem->ExpandPathName(macdir);

    void* dirhandle = gSystem->OpenDirectory(macdir.Data());
    if(dirhandle != 0)
    {
        const char* filename;
        TString runstr = "";
        TString revstr = "";
        runstr += fRun;
        TString reg;
        reg+= ".*QA\\.";
        reg+= (fRun==0) ? "[0-9].*" : revstr.Data();
        reg+= "\\.root$";
        TPRegexp reHist(reg);
        TPRegexp reMerged("Merged.QA.Data.[0-9]*.root");
        TPRegexp reCORR("CORR.QA.[0-9]*.root");
        TPRegexp reQA("QA.root");
        std::list<string> names;
        while((filename = gSystem->GetDirEntry(dirhandle)) != 0)
        {
            if (reCORR.Match(filename))
            {
                CORRFileName = filename;
                continue;
            }
            if (reQA.Match(filename))
            {
                QAResultFileName = filename;
                continue;
            }
            if (reMerged.Match(filename))
            {
                names.clear();
                names.push_back(filename);
                break;
            }
            if (reHist.Match(filename))
            {
                names.push_back(filename);
            }
        }
        if (!fPCORRFile) fPCORRFile = new TFile(CORRFileName,"READ");
        if (!fPQAResultFile) fPQAResultFile = new TFile(QAResultFileName,"READ");
        if (names.empty())
        {
            printf("GetListOfFiles: no files matching...\n");
            return kFALSE;
        }
        names.sort();
        char fullName[1000];
        for (std::list<string>::iterator si=names.begin(); si!=names.end(); ++si)
        {
          sprintf(fullName,"%s", si->c_str());
          AliQADirList* f = new AliQADirList();
          f->SetName(fullName);
          fPListOfFiles->GetDirs()->AddLast(f);
        }
    }
    else
    {
        gSystem->FreeDirectory(dirhandle);
        return kFALSE;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::CloneDirStructure()
{
    //scan all files
    if (!GetListOfFiles()) return kFALSE;
    if (fPListOfFiles->GetDirs()->GetEntries()==0) return kFALSE;
    TIter fileiter(fPListOfFiles->GetDirs());
    AliQADirList* f;
    while ((f = (AliQADirList*)fileiter.Next()))
    {
        TString filename = f->GetName();
        TFile file(filename);
        if (!Crawl(f)) continue;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Crawl(AliQADirList* dir)
{
    TString oldkeyname;
    TString keyname;
    //scans the current directory and puts result in dir
    TString pwd = gDirectory->GetPath();
    //returns false if dir in file empty
    TList* keys = gDirectory->GetListOfKeys();
    if (!keys) return kFALSE;
    if (keys->GetEntries()==0) return kFALSE;
    TIter keyiter(keys);
    TKey* key;
    while ((key = dynamic_cast <TKey* > (keyiter.Next()) ))
    {
        keyname = key->GetName();
        if (keyname==oldkeyname) continue; //do not copy cycles
        oldkeyname = keyname;
        TString classname=key->GetClassName();
        if (!classname) return kFALSE;
        if (classname=="TDirectoryFile")
        {
            gDirectory->cd(key->GetName());
            AliQADirList* newdir = new AliQADirList();
            if (!Crawl(newdir))
            {
                gDirectory->cd(pwd);
                continue;
            }
            gDirectory->cd(pwd);

            newdir->SetName(keyname);
            newdir->SetParent(dir);
            dir->GetDirs()->AddLast(newdir);
        }
        else
        {
            AliQADirListItem* item = new AliQADirListItem(keyname);
            item->SetParent(dir);
            dir->GetItems()->AddLast(item);
        }
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::ReReadFiles()
{
    //close, open and rescan the files
    if (!CloneDirStructure())
    {
        printf("AliQAHistNavigator::AliQAHistNavigator(): error reading files\n");
        return kFALSE;
    }
    fPCurrFile = (AliQADirList*)fPListOfFiles->GetDirs()->First();
    if (fPCurrFile)
    {
        fPCurrDetector = (AliQADirList*)fPCurrFile->GetDirs()->First();
        if (fPCurrDetector)
        {
            fPCurrLevel = (AliQADirList*) fPCurrDetector->GetDirs()->First();
            if (fPCurrLevel)
            {
                fPCurrItem = (AliQADirListItem*) GetItemList()->First();
                if (fPCurrItem)
                {
                    fInitOK = kTRUE;
                    OpenCurrentFile();
                    OpenCurrentDirectory();
                }
            }
        }
    }
    return kTRUE;
}

//_________________________________________________________________________
ClassImp(AliQADirList)
//_________________________________________________________________________
AliQADirList::AliQADirList():
    TNamed(),
    fPParent(NULL),
    fPItems(new TList()),
    fPDirs(new TList())
{
    //ctor
}

//_________________________________________________________________________
AliQADirList::~AliQADirList()
{
    //dtor
    if (fPParent) delete fPParent;
    delete fPItems;
    delete fPDirs;
}

//_________________________________________________________________________
ClassImp(AliQADirListItem)
//_________________________________________________________________________
AliQADirListItem::AliQADirListItem(const char* s):
    TObjString(s),
    fPParent(NULL)
{
    //ctor
}

//_________________________________________________________________________
AliQADirListItem::~AliQADirListItem()
{
    //dtor
    if (fPParent) delete fPParent;
}

