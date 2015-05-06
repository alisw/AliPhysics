#include "AliStorageDatabase.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliOnlineReconstructionUtil.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"

#include <iostream>
#include <fstream>

#include <TSQLRow.h>
#include <TSQLResult.h>
#include <TThread.h>
#include <TSystem.h>
#include <TFile.h>
#include <TEnv.h>

using namespace std;

AliStorageDatabase::AliStorageDatabase() :
fHost(""),
fPort(""),
fDatabase(""),
fUID(""),
fPassword(""),
fTable(""),
fServer(0),
fStoragePath("")
{
    TThread::Lock();
    ifstream configFile (GetConfigFilePath());
    if (configFile.is_open())
    {
        string line;
        int from,to;
        while(configFile.good())
        {
            getline(configFile,line);
            from = line.find("\"")+1;
            to = line.find_last_of("\"");
            if(line.find("HOST=")==0){
                fHost=line.substr(from,to-from);
            }
            else if(line.find("PORT=")==0){
                fPort=line.substr(from,to-from);
            }
            else if(line.find("DATABASE=")==0){
                fDatabase=line.substr(from,to-from);
            }
            else if(line.find("USER=")==0){
                fUID=line.substr(from,to-from);
            }
            else if(line.find("PASS=")==0){
                fPassword=line.substr(from,to-from);
            }
            else if(line.find("TABLE=")==0){
                fTable=line.substr(from,to-from);
            }
            else if(line.find("STORAGE_PATH=")==0){
                fStoragePath=line.substr(from,to-from);
            }
            
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout << "DATABASE -- Unable to open file" <<endl;}
    
    cout<<"DATABASE -- connecting to server:"<<Form("mysql://%s:%s/%s",fHost.c_str(),fPort.c_str(),fDatabase.c_str())<<fUID.c_str()<<fPassword.c_str()<<endl;
    fServer = TSQLServer::Connect(Form("mysql://%s:%s/%s",fHost.c_str(),fPort.c_str(),fDatabase.c_str()),fUID.c_str(),fPassword.c_str());
    
    TThread::UnLock();
    cout<<"Connected"<<endl;
}

AliStorageDatabase::~AliStorageDatabase(){
    if (fServer) {delete fServer;}
}

void AliStorageDatabase::InsertEvent(int runNumber,
                                     int eventNumber,
                                     char *system,
                                     int multiplicity,
                                     char *filePath,
                                     ULong64_t triggerMask,
                                     ULong64_t triggerMaskNext50)
{
    TSQLResult *res = fServer->Query(Form("select * FROM %s WHERE run_number = %d AND event_number = %d AND permanent = 1;",fTable.c_str(),runNumber,eventNumber));
    TSQLRow *row = res->Next();
    
    cout<<"DATABASE -- insterting:"<<Form("REPLACE INTO %s (run_number,event_number,system,multiplicity,permanent,file_path,trigger_mask,trigger_mask_next) VALUES (%d,%d,'%s',%d,0,'%s',%llu,%llu);",fTable.c_str(),runNumber,eventNumber,system,multiplicity,filePath,triggerMask,triggerMaskNext50)<<endl;
    
    if(!row)
    {
        res = fServer->Query(Form("REPLACE INTO %s (run_number,event_number,system,multiplicity,permanent,file_path,trigger_mask,trigger_mask_next) VALUES (%d,%d,'%s',%d,0,'%s',%llu,%llu);",fTable.c_str(),runNumber,eventNumber,system,multiplicity,filePath,triggerMask,triggerMaskNext50));
    }
    
    delete row;
    delete res;
}

bool AliStorageDatabase::MarkEvent(struct eventStruct event)
{
    TSQLResult *res;
    res = fServer->Query(Form("UPDATE %s SET permanent = 1 WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber));
    if(!res)
    {
        cout<<"DATABASE -- couldn't update permanent flag"<<endl;
        delete res;
        return 0;
    }
    else
    {
        cout<<"DATABASE -- permanent flag updated"<<endl;
        
        res = fServer->Query(Form("UPDATE %s SET file_path = '%s' WHERE run_number = %d AND event_number = %d;",fTable.c_str(),Form("%s/permEvents.root",fStoragePath.c_str()), event.runNumber,event.eventNumber));
        if(!res)
        {
            cout<<"DATABASE -- couldn't update file's path. Unsetting permanent flag"<<endl;
            res = fServer->Query(Form("UPDATE %s SET permanent = 0 WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber));
            delete res;
            return 0;
        }
        else
        {
            cout<<"DATABASE -- event marked"<<endl;
            delete res;
            return 1;
        }
    }
}

bool AliStorageDatabase::UpdateEventPath(struct eventStruct event,const char *newPath)
{
    TSQLResult* res;
    res = fServer->Query(Form("UPDATE %s SET file_path = '%s' WHERE run_number = %d AND event_number = %d;",fTable.c_str(),newPath,event.runNumber,event.eventNumber));
    if(!res)
    {
        cout<<"DATABASE -- couldn't update file's path"<<endl;
        delete res;
        return 0;
    }
    else
    {
        cout<<"DATABASE -- path updated for event:"<<event.eventNumber<<endl;
        delete res;
        return 1;
    }
}


vector<serverListStruct> AliStorageDatabase::GetList(struct listRequestStruct list)
{
    cout<<"LIST:"<< list.runNumber[0]<<"\t"<<list.runNumber[1]<<"\t"<<list.eventNumber[0]<<"\t"<<list.eventNumber[1]<<"\t"<< list.multiplicity[0]<<"\t"<< list.multiplicity[1]<<"\t"<<list.marked[0]<<"\t"<< list.marked[1]<<"\t"<<list.system[0]<<"\t"<<list.system[1]<<"\t"<<list.triggerClass<<endl;
    
    ULong64_t triggerMask;
    ULong64_t triggerMaskNext50;
    
    vector<serverListStruct> eventsVector;
    
    if(strcmp(list.triggerClass,"No trigger selection")!=0 && strcmp(list.triggerClass,"")!=0)
    {
        
        // craete CDB manager:
        TEnv settings;
        settings.ReadFile(AliOnlineReconstructionUtil::GetPathToServerConf(), kEnvUser);
        const char *cdbPath = settings.GetValue("cdb.defaultStorage", "");
        AliCDBManager *man = AliCDBManager::Instance();
        man->SetDefaultStorage(cdbPath);
        
        // get list of runs stored in Storage's database:
        vector<int> runs = GetListOfRuns();
        
        // get trigger classes for all runs:
        AliCDBEntry *cdbEntry;
        AliTriggerConfiguration *cfg;
        TObjArray trarr;
        AliCDBPath path("GRP/CTP/Config");
        
        AliTriggerClass* trgclass;
        
        for(int i=0;i<runs.size();i++)
        {
            if(runs[i] > list.runNumber[0] && runs[i] < list.runNumber[1])
            {
                man->SetRun(runs[i]);
                cdbEntry = man->Get(path);
                cfg = (AliTriggerConfiguration*)cdbEntry->GetObject();
                trarr = cfg->GetClasses();
                
                triggerMask = 0;
                triggerMaskNext50 = 0;
                
                for (int j=0;j<trarr.GetEntriesFast();j++)
                {
                    trgclass = (AliTriggerClass*)trarr.At(j);
                    if(strcmp(trgclass->GetName(),list.triggerClass)==0)
                    {
                        triggerMask = trgclass->GetMask();
                        triggerMaskNext50 = trgclass->GetMaskNext50();
                    }
                }
                
                TThread::Lock();
                TSQLResult *result = NULL;
                
                result =  fServer->Query(Form("SELECT * FROM %s WHERE run_number = %d AND event_number >= %d AND event_number <= %d AND multiplicity >= %d AND multiplicity <= %d AND (permanent = %d OR permanent = %d) AND (system = '%s' OR system = '%s') AND ((trigger_mask & %llu) > 0 OR (trigger_mask_next & %llu) > 0) ORDER BY run_number,event_number;",
                                              fTable.c_str(),
                                              runs[i],
                                              list.eventNumber[0],
                                              list.eventNumber[1],
                                              list.multiplicity[0],
                                              list.multiplicity[1],
                                              list.marked[0],
                                              list.marked[1],
                                              list.system[0],
                                              list.system[1],
                                              triggerMask,
                                              triggerMaskNext50));
                
                cout<<"Query:"<<Form("SELECT * FROM %s WHERE run_number = %d AND event_number >= %d AND event_number <= %d AND multiplicity >= %d AND multiplicity <= %d AND (permanent = %d OR permanent = %d) AND (system = '%s' OR system = '%s') AND ((trigger_mask & %llu) > 0 OR (trigger_mask_next & %llu) > 0) ORDER BY run_number,event_number;",
                                     fTable.c_str(),
                                     runs[i],
                                     list.eventNumber[0],
                                     list.eventNumber[1],
                                     list.multiplicity[0],
                                     list.multiplicity[1],
                                     list.marked[0],
                                     list.marked[1],
                                     list.system[0],
                                     list.system[1],
                                     triggerMask,
                                     triggerMaskNext50)<<endl;
                
                TThread::UnLock();
                TSQLRow *row;
                
                while((row = result->Next()))
                {
                    serverListStruct resultList;
                    
                    resultList.runNumber = atoi(row->GetField(0));
                    resultList.eventNumber = atoi(row->GetField(1));
                    strcpy(resultList.system, row->GetField(2));
                    resultList.multiplicity = atoi(row->GetField(3));
                    resultList.marked = atoi(row->GetField(4));
                    strcpy(resultList.triggerClass,list.triggerClass);
                    
                    eventsVector.push_back(resultList);
                    delete row;
                }
                delete result;
                
            }
        }
    }
    else
    {
        TThread::Lock();
        TSQLResult *result = NULL;
        
        result =  fServer->Query(Form("SELECT * FROM %s WHERE run_number >= %d AND run_number <= %d AND event_number >= %d AND event_number <= %d AND multiplicity >= %d AND multiplicity <= %d AND (permanent = %d OR permanent = %d) AND (system = '%s' OR system = '%s') ORDER BY run_number,event_number;",
                                      fTable.c_str(),
                                      list.runNumber[0],
                                      list.runNumber[1],
                                      list.eventNumber[0],
                                      list.eventNumber[1],
                                      list.multiplicity[0],
                                      list.multiplicity[1],
                                      list.marked[0],
                                      list.marked[1],
                                      list.system[0],
                                      list.system[1]));
        
        cout<<"Query:"<<Form("SELECT * FROM %s WHERE run_number >= %d AND run_number <= %d AND event_number >= %d AND event_number <= %d AND multiplicity >= %d AND multiplicity <= %d AND (permanent = %d OR permanent = %d) AND (system = '%s' OR system = '%s') ORDER BY run_number,event_number;",
                             fTable.c_str(),
                             list.runNumber[0],
                             list.runNumber[1],
                             list.eventNumber[0],
                             list.eventNumber[1],
                             list.multiplicity[0],
                             list.multiplicity[1],
                             list.marked[0],
                             list.marked[1],
                             list.system[0],
                             list.system[1])<<endl;
        
        TThread::UnLock();
        TSQLRow *row;
        
        while((row = result->Next()))
        {
            serverListStruct resultList;
            
            resultList.runNumber = atoi(row->GetField(0));
            resultList.eventNumber = atoi(row->GetField(1));
            strcpy(resultList.system, row->GetField(2));
            resultList.multiplicity = atoi(row->GetField(3));
            resultList.marked = atoi(row->GetField(4));
            strcpy(resultList.triggerClass,list.triggerClass);
            
            eventsVector.push_back(resultList);
            delete row;
        }
        delete result;
        
    }
    
    return eventsVector;
}

AliESDEvent* AliStorageDatabase::GetEvent(struct eventStruct event)
{
    cout<<"database - get event:"<<event.runNumber<<"\t"<<event.eventNumber<<endl;
    string pathToFile = GetFilePath(event);
    
    cout<<"DATABASE -- path to file:"<<pathToFile<<endl;
    
    if(!strcmp(pathToFile.c_str(),""))
    {
        cout<<"DATABASE -- no such file in database"<<endl;
        return NULL;
    }
    
    TFile *tmpFile = new TFile(pathToFile.c_str(),"read");
    if(!tmpFile)
    {
        cout<<"DATABASE -- couldn't open temp file"<<endl;
        return NULL;
    }
    tmpFile->cd(Form("run%d",event.runNumber));
    AliESDEvent *data = (AliESDEvent*)gDirectory->Get(Form("event%d;1",event.eventNumber));
    tmpFile->Close();
    delete tmpFile;
    
    if(data)
    {
        cout<<"DATABASE -- read event:"<<data->GetEventNumberInFile()<<endl;
    }
    else
    {
        cout<<"DATABASE -- event is corrupted"<<endl;
    }
    
    return data;
}

void AliStorageDatabase::RemoveEvent(struct eventStruct event)
{
    TSQLResult* res;
    res = fServer->Query(Form("DELETE FROM %s WHERE run_number = %d AND event_number = %d",fTable.c_str(),event.runNumber,event.eventNumber));
    delete res;
}

void AliStorageDatabase::RemoveEventsWithPath(string path)
{
    TSQLResult *res = fServer->Query(Form("DELETE FROM %s WHERE file_path = \"%s\";",fTable.c_str(),path.c_str()));
    delete res;
}

string AliStorageDatabase::GetFilePath(struct eventStruct event)
{
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber));
    TSQLRow *row;
    row = result->Next();
    if(row)
    {
        string path(row->GetField(5));
        delete row;
        return path;
    }
    else
    {
        return "";
    }
}

vector<int> AliStorageDatabase::GetListOfRuns()
{
    vector<int> resultingVector;
    
    TSQLResult *result = fServer->Query(Form("SELECT run_number FROM %s ORDER BY run_number;",fTable.c_str()));
    TSQLRow *row;
    
    int currentRun,prevRun=-1;
    
    while((row=result->Next()))
    {
        currentRun = atoi(row->GetField(0));
        if(currentRun!=prevRun)
        {
            resultingVector.push_back(currentRun);
            prevRun=currentRun;
        }
    }
    delete row;
    return resultingVector;
}

AliESDEvent* AliStorageDatabase::GetNextEvent(struct eventStruct event)
{
    cout<<"Database:"<<event.runNumber<<"\t"<<event.eventNumber<<endl;
    
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));
    
    TSQLRow *row;
    bool isCurrentEvent=false;
    struct eventStruct nextEvent;
    
    while((row = result->Next()))
    {
        if(isCurrentEvent)
        {
            nextEvent.runNumber = atoi(row->GetField(0));
            nextEvent.eventNumber = atoi(row->GetField(1));
            return GetEvent(nextEvent);
        }
        
        //if current event found
        if(atoi(row->GetField(0))==event.runNumber && atoi(row->GetField(1))==event.eventNumber)
        {
            isCurrentEvent=true;
        }
        
        delete row;
    }
    
    return NULL;
}

AliESDEvent* AliStorageDatabase::GetPrevEvent(struct eventStruct event)
{
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number DESC;",fTable.c_str()));
    
    TSQLRow *row;
    bool isCurrentEvent=false;
    struct eventStruct nextEvent;
    
    while((row = result->Next()))
    {
        if(isCurrentEvent)
        {
            nextEvent.runNumber = atoi(row->GetField(0));
            nextEvent.eventNumber = atoi(row->GetField(1));
            return GetEvent(nextEvent);
        }
        
        //if current event found
        if(atoi(row->GetField(0))==event.runNumber && atoi(row->GetField(1))==event.eventNumber)
        {
            isCurrentEvent=true;
        }
        
        delete row;
    }
    delete result;
    return NULL;
}

struct eventStruct AliStorageDatabase::GetOldestEvent()
{
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s WHERE permanent = 0 ORDER BY run_number,event_number;",fTable.c_str()));
    
    TSQLRow *row;
    struct eventStruct oldestEvent = {0,0};
    
    if((row = result->Next()))
    {
        oldestEvent.runNumber = atoi(row->GetField(0));
        oldestEvent.eventNumber = atoi(row->GetField(1));
        delete row;
    }
    else
    {
        cout<<"DATABASE -- NO OLDEST EVENT FOUND. Storage may be corrupted."<<endl;
    }
    return oldestEvent;
}

AliESDEvent* AliStorageDatabase::GetLastEvent()
{
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));
    
    TSQLRow *row;
    struct eventStruct lastEvent = {0,0};
    
    while((row = result->Next()))
    {
        lastEvent.runNumber = atoi(row->GetField(0));
        lastEvent.eventNumber = atoi(row->GetField(1));
        delete row;
    }
    cout<<"Last event is:"<<lastEvent.eventNumber<<endl;
    return GetEvent(lastEvent);
    
}

AliESDEvent* AliStorageDatabase::GetFirstEvent()
{
    cout<<"Database - first"<<endl;
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number DESC;",fTable.c_str()));
    
    TSQLRow *row;
    struct eventStruct firstEvent = {0,0};
    
    while((row = result->Next()))
    {
        firstEvent.runNumber = atoi(row->GetField(0));
        firstEvent.eventNumber = atoi(row->GetField(1));
        delete row;
    }
    cout<<"First event is:"<<firstEvent.eventNumber<<endl;
    return GetEvent(firstEvent);
    
}
