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
AliQAHistNavigator::AliQAHistNavigator(Int_t run, Int_t rev):
    fPFile( NULL ),
    fRun( run ),
    fCyc( rev ),
    fPCurrFile( NULL ),
    fPCurrDetector( NULL ),
    fPCurrLevel( NULL ),
    fPCurrHistName( NULL ),
    fPListOfFiles( new TList() ),
    fLoopAllFiles(kTRUE),
    fLoopAllDetectors(kTRUE),
    fLoopAllLevels(kTRUE)
{
    if (CloneDirStructure())
    {
        fPCurrFile = (TList*)fPListOfFiles->First();
        if (fPCurrFile)
        {
            fPCurrDetector = (TList*)fPCurrFile->First();
            if (fPCurrDetector)
            {
                fPCurrLevel = (TList*) fPCurrDetector->First();
                if (fPCurrLevel)
                {
                    fPCurrHistName = (TObjString*) fPCurrLevel->First();
                    if (fPCurrHistName)
                    {
                        OpenCurrentFile();
                        OpenCurrentDirectory();
                    }
                }
            }
        }
    } else
    {
        printf("AliQAHistNavigator::AliQAHistNavigator(): error reading files\n");
    }
    
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetNextHistogram(TH1*& hist)
{
    //moves to the next histogram from the list and tries to get it.
    if (!Next()) 
    {
        hist = NULL;
        return kFALSE;
    }
    if (GetHistName()=="")
        if (!Next())
        {
            hist = NULL;
            return kTRUE;
        }
    return GetHistogram(hist);
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetHistogram(TH1*& hist)
{
    gDirectory->GetObject(GetFileName()+":"+GetDirName()+"/"+GetHistName(),hist);
    if (!hist)
    {
        printf("GetHistogram: null pointer returned by gDirectory\n");
        return kFALSE;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetPrevHistogram(TH1*& hist)
{
    //moves to the prev histogram from the list and tries to get it.
    if (!Prev()) 
    {
        hist = NULL;
        return kFALSE;
    }
    if (GetHistName()=="")
        if (!Prev())
        {
            hist = NULL;
            return kTRUE;
        }
    return GetHistogram(hist);
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Next()
{
    if (!fPCurrHistName) return kFALSE;
    if (!(fPCurrHistName=(TObjString*)fPCurrLevel->After(fPCurrHistName)))
    {
        if (!(fPCurrLevel=(TList*)fPCurrDetector->After(fPCurrLevel)))
        {
            if (!(fPCurrDetector=(TList*)fPCurrFile->After(fPCurrDetector)))
            {
                if (!(fPCurrFile=(TList*)fPListOfFiles->After(fPCurrFile)))
                {
                    //we're at the end of everything
                    if (fLoopAllFiles)
                    {
                        //rewind to the beginning
                        fPCurrFile = (TList*)fPListOfFiles->First();
                        fPCurrDetector = (TList*)fPCurrFile->First();
                        fPCurrLevel = (TList*) fPCurrDetector->First();
                        fPCurrHistName = (TObjString*) fPCurrLevel->First();
                        OpenCurrentFile();
                        OpenCurrentDirectory();
                        printf("----------------back at the beginning!\n");
                    } else return kFALSE; //no rewind, we finish
                } else //if there is a next file
                {
                    fPCurrDetector = (TList*)fPCurrFile->First();
                    fPCurrLevel=(TList*)fPCurrDetector->First();
                    fPCurrHistName=(TObjString*)fPCurrLevel->First();
                    cout<<GetFileName()<<":"<<GetDetectorName()<<"/"<<GetLevelName()<<"/"<<GetHistName()<<endl;
                    OpenCurrentFile();
                    OpenCurrentDirectory();
                }
            } else //if there is a next detector
            {
                fPCurrLevel=(TList*)fPCurrDetector->First();
                fPCurrHistName=(TObjString*)fPCurrLevel->First();
                cout<<GetDetectorName()<<"/"<<GetLevelName()<<"/"<<GetHistName()<<endl;
                OpenCurrentDirectory();
            }
        } else //if there is a next level
        {
            fPCurrHistName=(TObjString*)fPCurrLevel->First();
            cout<<GetLevelName()<<"/"<<GetHistName()<<endl;
            OpenCurrentDirectory();
            cout<<GetHistName()<<endl;
        }
    } else //if there is a next histgram
    {
        cout<<GetHistName()<<endl;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Prev()
{
    if (!fPCurrHistName) return kFALSE;
    if (!(fPCurrHistName=(TObjString*)fPCurrLevel->Before(fPCurrHistName)))
    {
        if (!(fPCurrLevel=(TList*)fPCurrDetector->Before(fPCurrLevel)))
        {
            if (!(fPCurrDetector=(TList*)fPCurrFile->Before(fPCurrDetector)))
            {
                if (!(fPCurrFile=(TList*)fPListOfFiles->Before(fPCurrFile)))
                {
                    //we're at the end of everything
                    if (fLoopAllFiles)
                    {
                        //rewind to the beginning
                        fPCurrFile = (TList*)fPListOfFiles->Last();
                        fPCurrDetector = (TList*)fPCurrFile->Last();
                        fPCurrLevel = (TList*) fPCurrDetector->Last();
                        fPCurrHistName = (TObjString*) fPCurrLevel->Last();
                        OpenCurrentFile();
                        OpenCurrentDirectory();
                        printf("----------------back at the end!\n");
                    } else return kFALSE; //no rewind, we finish
                } else //if there is a next file
                {
                    fPCurrDetector = (TList*)fPCurrFile->Last();
                    fPCurrLevel=(TList*)fPCurrDetector->Last();
                    fPCurrHistName=(TObjString*)fPCurrLevel->Last();
                    cout<<GetFileName()<<":"<<GetDetectorName()<<"/"<<GetLevelName()<<"/"<<GetHistName()<<endl;
                    OpenCurrentFile();
                    OpenCurrentDirectory();
                }
            } else //if there is a next detector
            {
                fPCurrLevel=(TList*)fPCurrDetector->Last();
                fPCurrHistName=(TObjString*)fPCurrLevel->Last();
                cout<<GetDetectorName()<<"/"<<GetLevelName()<<"/"<<GetHistName()<<endl;
                OpenCurrentDirectory();
            }
        } else //if there is a next level
        {
            fPCurrHistName=(TObjString*)fPCurrLevel->Last();
            cout<<GetLevelName()<<"/"<<GetHistName()<<endl;
            OpenCurrentDirectory();
            cout<<GetHistName()<<endl;
        }
    } else //if there is a next histgram
    {
        cout<<GetHistName()<<endl;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::OpenCurrentFile()
{
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
    if (!gDirectory->cd(GetDirName())) return kFALSE;
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetFile( TString file )
{
    TList* tmp = (TList*)fPListOfFiles->FindObject ( file.Data() );
    if (!tmp) return kFALSE;
    fPCurrFile = tmp;
    OpenCurrentFile();
    fPCurrDetector = (TList*)fPCurrFile->First();
    fPCurrLevel = (TList*)fPCurrDetector->First();
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetFile( Int_t file )
{
    printf("AliQAHistNavigator::SetFile(%i)\n",file);
    TList* tmp = (TList*)fPListOfFiles->At(file);
    if (!tmp) return kFALSE;
    fPCurrFile = tmp;
    OpenCurrentFile();
    fPCurrDetector = (TList*)fPCurrFile->First();
    fPCurrLevel = (TList*)fPCurrDetector->First();
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetDetector( TString det )
{
    TList* tmp = (TList*)fPCurrFile->FindObject( det.Data() );
    if (!tmp) return kFALSE;
    fPCurrDetector = tmp;
    fPCurrLevel = (TList*)fPCurrDetector->First();
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetDetector( Int_t det )
{
    printf("AliQAHistNavigator::SetDetector(%i)\n",det);
    TList* tmp = (TList*)fPCurrFile->At( det );
    if (!tmp) return kFALSE;
    fPCurrDetector = tmp;
    fPCurrLevel = (TList*)fPCurrDetector->First();
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetLevel( TString level )
{
    TList* tmp = (TList*)fPCurrDetector->FindObject( level.Data() );
    if (!tmp) return kFALSE;
    fPCurrLevel = tmp;
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetLevel( Int_t level )
{
    TList* tmp = (TList*)fPCurrDetector->At( level );
    if (!tmp) return kFALSE;
    fPCurrLevel = tmp;
    fPCurrHistName = (TObjString*)fPCurrLevel->First();
    OpenCurrentDirectory();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetHist( TString hist )
{
    TObjString* tmp = (TObjString*)fPCurrLevel->FindObject( hist.Data() );
    if (!tmp) return kFALSE;
    fPCurrHistName = tmp;
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::SetHist( Int_t hist )
{
    TObjString* tmp = (TObjString*)fPCurrLevel->At( hist );
    if (!tmp) return kFALSE;
    fPCurrHistName = tmp;
    return kTRUE;
}

//_________________________________________________________________________
void AliQAHistNavigator::PrintDebugInfo()
{
    if (!fPCurrHistName) {cout<<"no more histograms"<<endl;return;};
    if (!fPCurrLevel) {cout<<"no more levels"<<endl;return;};
    if (!fPCurrDetector) {cout<<"no more detectors"<<endl;return;};
    cout<<"AliQAHistNavigator state information:"<<endl;
    cout<<"hist: "<<GetHistName()<<endl;
    cout<<"dir:  "<<GetDirName()<<endl;
    cout<<"$PWD: ";gDirectory->pwd();cout<<endl;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetDetectorName()
{
    if (!fPCurrDetector) return "";
    TString name = fPCurrDetector->GetName();
    return name;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetLevelName()
{
    if (!fPCurrLevel) return "";
    TString name = fPCurrLevel->GetName();
    return name;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetFileName()
{
    if (!fPCurrFile) return "";
    TString file = fPCurrFile->GetName();
    return file;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetDirName()
{
    TString dir = "/"+ GetDetectorName()+"/"+GetLevelName();
    return dir;
}

//_________________________________________________________________________
TString AliQAHistNavigator::GetHistName()
{
    if (!fPCurrHistName) return "";
    return fPCurrHistName->GetString();
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::DumpList( TString filename )
{
    ofstream fout(filename);
    if (fout.bad()) return kFALSE;
    TString lastlevel="";
    TString lastdet="";
    if (!Next()) return kTRUE;
    do
    {
        if (GetLevelName()!=lastlevel)
        {
            fout<<GetDetectorName()<<"/"<<GetLevelName()<<":"<<endl;
            lastlevel=GetLevelName();
        }
        fout << GetHistName() << endl;
    } while (Next());
    fout.close();
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::ReadList( TString filename )
{
    TString line = filename;
    //TString level="";
    //TString lastlevel="";
    //TString det="";
    //TString lastdet="";
    //TRegexp detRegexp("^[a-zA-Z0-9]*");
    //TRegexp typRegexp("[a-zA-Z0-9]*");

    //ifstream fin(filename);
    //if (fin.bad()) return kFALSE;
    //delete fPListOfFiles;
    //fPListOfFiles = new TList();
    //TList* pDetector=new TList();
    //TList* pLevel=new TList();
    //
    //line.ReadLine(fin,kFALSE);
    //if (line=="") return kFALSE;
    //cout<<"read line: "<<line<<endl;
    //det = line(detRegexp);
    //level = line(typRegexp, line.Index("/")+1);
    //pDetector->SetName(det);
    //lastdet = det;
    //pLevel->SetName(level);
    //lastlevel = level;
    //pLevel->AddLast(new TObjString("DO NOT REMOVE THIS LINE"));
    //
    //while (!fin.eof())
    //{
    //    line.ReadLine(fin,kFALSE);
    //    cout<<"read line: "<<line<<endl;
    //    if (line.EndsWith(":"))
    //    {
    //        det = line(detRegexp);
    //        level = line(typRegexp, line.Index("/")+1);
    //        if (det!=lastdet)
    //        {
    //            pDetector->AddLast(pLevel);
    //            fPListOfFiles->AddLast(pDetector);
    //            pDetector = new TList();
    //            pDetector->SetName(det);
    //            cout<<"new detector: "<<det<<endl;
    //            lastdet = det;
    //            pLevel = new TList();
    //            pLevel->SetName(level);
    //            cout<<"new level: "<<level<<endl;
    //            lastlevel = level;
    //            continue;
    //        }
    //        if (level!=lastlevel)
    //        {
    //            pDetector->AddLast(pLevel);
    //            pLevel = new TList();
    //            pLevel->SetName(level);
    //            cout<<"new level: "<<level<<endl;
    //            lastlevel = level;
    //            continue;
    //        }
    //    }
    //    if (line.BeginsWith("//")) continue;
    //    if (line.BeginsWith("#")) continue;
    //    pLevel->AddLast(new TObjString(line));        
    //    cout<<"added line: "<<line<<endl;
    //}
    //
    //fPCurrDetector = (TList*)fPListOfFiles->First();
    //fPCurrLevel = (TList*) fPCurrDetector->First();
    //fPCurrHistName = (TObjString*) fPCurrLevel->First();
    //OpenCurrentFile();
    //OpenCurrentDirectory();
    //fPListOfFiles->Print();
    return kFALSE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::GetListOfFiles()
{
    delete fPListOfFiles;
    fPListOfFiles = new TList();

    TString macdir(".");
    gSystem->ExpandPathName(macdir);

    void* dirhandle = gSystem->OpenDirectory(macdir.Data());
    if(dirhandle != 0)
    {
        const char* filename;
        TString runstr = "";
        TString revstr = "";
        runstr += fRun;
        revstr += fCyc;
        TString reg;
        reg+= ".*QA\\.";
        reg+= (fRun==0) ? "[0-9].*" : runstr.Data();
        reg+= "\\.";
        reg+= (fRun==0) ? "[0-9].*" : revstr.Data();
        reg+= "\\.root$";
        cout<<reg<<endl;
        TPRegexp re(reg);
        std::list<string> names;
        while((filename = gSystem->GetDirEntry(dirhandle)) != 0)
        {
            if(re.Match(filename))
            {
                names.push_back(filename);
            }
        }
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
          TList* f = new TList();
          f->SetName(fullName);
          fPListOfFiles->AddLast(f);
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
    if (!GetListOfFiles()) return kFALSE;
    if (fPListOfFiles->GetEntries()==0) return kFALSE;
    TIter fileiter(fPListOfFiles);
    TList* f;
    while ((f = (TList*)fileiter.Next()))
    {
        TString filename = f->GetName();
        cout<<filename<<endl;
        TFile file(filename);
        if (!Crawl(f)) continue;
    }
    return kTRUE;
}

//_________________________________________________________________________
Bool_t AliQAHistNavigator::Crawl(TList* dir)
{
    TString pwd = gDirectory->GetPath();
    //returns false if dir in file empty
    TList* keys = gDirectory->GetListOfKeys();
    if (!keys) return kFALSE;
    if (keys->GetEntries()==0) return kFALSE;
    TIter keyiter(keys);
    TKey* key;
    while ((key = dynamic_cast <TKey* > (keyiter.Next()) ))
    {
        TString classname=key->GetClassName();
        if (!classname) return kFALSE;
        if (classname=="TDirectoryFile")
        {
            gDirectory->cd(key->GetName());
            gDirectory->pwd();
            TList* newdir = new TList();
            if (!Crawl(newdir))
            {
                gDirectory->cd(pwd);
                continue;
            }
            gDirectory->cd(pwd);

            newdir->SetName(key->GetName());
            dir->AddLast(newdir);
        }
        else
        {
            cout<<key->GetName()<<endl;
            dir->AddLast(new TObjString(key->GetName()));
        }
    }
    return kTRUE;
}
