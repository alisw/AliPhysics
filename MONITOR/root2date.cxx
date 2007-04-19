// Author: Filimon Roukoutakis 02/08/2006
//         Cvetan Cheshkov

/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#include "root2date.h"

int Root2Date(AliRawEvent *gdcRootEvent, unsigned char *gdcDateEvent) {

 unsigned char *p=gdcDateEvent;
 int ldcCounter, equipmentCounter, chunkSize;
 AliRawEquipment *aliEquipment=NULL;
 AliRawEquipmentHeader *aliEquipmentHeader=NULL;
 AliRawEventHeaderBase *aliHeader=NULL;
 AliRawEvent *ldcRootEvent=NULL;
 
 aliHeader=gdcRootEvent->GetHeader();
 memcpy(p, aliHeader->HeaderBaseBegin(), chunkSize=aliHeader->HeaderBaseSize());
 p+=chunkSize;
 memcpy(p, aliHeader->HeaderBegin(), chunkSize=aliHeader->HeaderSize()); // Write DATE GDC header
 p+=chunkSize;
 for(ldcCounter=0; ldcCounter<gdcRootEvent->GetNSubEvents(); ldcCounter++) {
  ldcRootEvent=gdcRootEvent->GetSubEvent(ldcCounter);
  aliHeader=ldcRootEvent->GetHeader();
  memcpy(p, aliHeader->HeaderBaseBegin(), chunkSize=aliHeader->HeaderBaseSize());
  p+=chunkSize;
  memcpy(p, aliHeader->HeaderBegin(), chunkSize=aliHeader->HeaderSize()); // Write DATE LDC header
  p+=chunkSize;
  for(equipmentCounter=0; equipmentCounter<ldcRootEvent->GetNEquipments(); equipmentCounter++) {
   aliEquipment=ldcRootEvent->GetEquipment(equipmentCounter);
   aliEquipmentHeader=aliEquipment->GetEquipmentHeader();
   if(aliEquipmentHeader->GetEquipmentSize()) {
    memcpy(p, aliEquipmentHeader->HeaderBegin(), chunkSize=aliEquipmentHeader->HeaderSize()); // Write DATE Equipment header
    p+=chunkSize;
   }
   memcpy(p, aliEquipment->GetRawData()->GetBuffer(), chunkSize=aliEquipment->GetRawData()->GetSize()); // Write Equipment payload (including CDH)
   p+=chunkSize;
  }
 }
 
 return(p-gdcDateEvent);
 
}
