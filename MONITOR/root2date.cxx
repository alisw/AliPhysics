// Author: Filimon Roukoutakis 02/08/2006
//         Cvetan Cheshkov

/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#include "root2date.h"

int Root2Date(AliRawVEvent *gdcRootEvent, unsigned char *gdcDateEvent, char *ddlDir) {

 unsigned char *p=gdcDateEvent;
 int ldcCounter, equipmentCounter, chunkSize;
 AliRawVEquipment *aliEquipment=NULL;
 AliRawEquipmentHeader *aliEquipmentHeader=NULL;
 AliRawEventHeaderBase *aliHeader=NULL;
 AliRawVEvent *ldcRootEvent=NULL;
 
 aliHeader=gdcRootEvent->GetHeader();

 char runNbFileName[256];
 sprintf(runNbFileName,"%s/run%u",ddlDir,aliHeader->Get("RunNb"));
 ofstream runNbFile(runNbFileName);
 runNbFile.close();

 memcpy(p, aliHeader->HeaderBaseBegin(), chunkSize=aliHeader->HeaderBaseSize());
 p+=chunkSize;
 memcpy(p, aliHeader->HeaderBegin(), chunkSize=aliHeader->HeaderSize()); // Write DATE GDC header
 p+=chunkSize;
 memcpy(p, aliHeader->GetExtendedData(), chunkSize=aliHeader->GetExtendedDataSize());
 p+=chunkSize;
 for(ldcCounter=0; ldcCounter<gdcRootEvent->GetNSubEvents(); ldcCounter++) {
  ldcRootEvent=gdcRootEvent->GetSubEvent(ldcCounter);
  aliHeader=ldcRootEvent->GetHeader();
  memcpy(p, aliHeader->HeaderBaseBegin(), chunkSize=aliHeader->HeaderBaseSize());
  p+=chunkSize;
  memcpy(p, aliHeader->HeaderBegin(), chunkSize=aliHeader->HeaderSize()); // Write DATE LDC header
  p+=chunkSize;
  memcpy(p, aliHeader->GetExtendedData(), chunkSize=aliHeader->GetExtendedDataSize());
  p+=chunkSize;
  for(equipmentCounter=0; equipmentCounter<ldcRootEvent->GetNEquipments(); equipmentCounter++) {
   aliEquipment=ldcRootEvent->GetEquipment(equipmentCounter);
   aliEquipmentHeader=aliEquipment->GetEquipmentHeader();
   if(aliEquipmentHeader->GetEquipmentSize()) {
    memcpy(p, aliEquipmentHeader->HeaderBegin(), chunkSize=aliEquipmentHeader->HeaderSize()); // Write DATE Equipment header
    p+=chunkSize;
   }
   memcpy(p, aliEquipment->GetRawData()->GetBuffer(), chunkSize=aliEquipment->GetRawData()->GetSize()); // Write Equipment payload (including CDH)
   // Write ddl files if requested by the user
   if (ddlDir && aliEquipmentHeader->GetEquipmentSize()) {
     Int_t ddlIndex;
     Int_t detId = AliDAQ::DetectorIDFromDdlID(aliEquipmentHeader->GetId(),ddlIndex);
     char ddlFileName[256];
     sprintf(ddlFileName,"%s/%s",ddlDir,AliDAQ::DdlFileName(detId,ddlIndex));
     FILE *ddlFile;
     if((ddlFile=fopen(ddlFileName, "wb"))) {
       fwrite(p, chunkSize, 1, ddlFile);
     }
     fclose(ddlFile);
   }
   p+=chunkSize;
  }
 }
 
 return(p-gdcDateEvent);
 
}
