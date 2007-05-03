#ifndef AliFemtoCutMonitorHandler_hh
#define AliFemtoCutMonitorHandler_hh


#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoEvent.h"
#include "Infrastructure/AliFemtoTrack.h"
#include "Infrastructure/AliFemtoV0.h"
#include "Infrastructure/AliFemtoKink.h"
#include "Infrastructure/AliFemtoPair.h" //Gael 12/04/02
#include "Infrastructure/AliFemtoParticleCollection.h" // Gael 19/06/02
#include "Infrastructure/AliFemtoCutMonitorCollection.h"
#include "Base/AliFemtoCutMonitor.h"

class AliFemtoCutMonitorHandler{
  
 public:
  
  AliFemtoCutMonitorHandler();
  AliFemtoCutMonitorHandler(const AliFemtoCutMonitorHandler& aHan);
  virtual ~AliFemtoCutMonitorHandler();
  AliFemtoCutMonitorHandler& operator=(const AliFemtoCutMonitorHandler& aHan);

  AliFemtoCutMonitorCollection* PassMonitorColl(); 
  AliFemtoCutMonitorCollection* FailMonitorColl(); 
  AliFemtoCutMonitor* PassMonitor(int n); 
  AliFemtoCutMonitor* FailMonitor(int n); 
  void AddCutMonitor(AliFemtoCutMonitor* cutMoni1, AliFemtoCutMonitor* cutMoni2); 
  void AddCutMonitor(AliFemtoCutMonitor* cutMoni); 
  void AddCutMonitorPass(AliFemtoCutMonitor* cutMoni); 
  void AddCutMonitorFail(AliFemtoCutMonitor* cutMoni); 
  void FillCutMonitor(const AliFemtoEvent* event, bool pass); 
  void FillCutMonitor(const AliFemtoTrack* track, bool pass); 
  void FillCutMonitor(const AliFemtoV0* v0, bool pass); 
  void FillCutMonitor(const AliFemtoKink* kink, bool pass);
  void FillCutMonitor(const AliFemtoPair* pair, bool pass);//Gael 11/04/02
  void FillCutMonitor(const AliFemtoParticleCollection* partColl);// Gael 19/06/02
  void FillCutMonitor(const AliFemtoEvent* event, const AliFemtoParticleCollection* partColl);// Gael 19/06/02
  void Finish();
  
 private:
  bool fCollectionsEmpty;
  AliFemtoCutMonitorCollection* fPassColl; 
  AliFemtoCutMonitorCollection* fFailColl; 
#ifdef __ROOT__  
  ClassDef(AliFemtoCutMonitorHandler, 0)
#endif  
  
};

inline AliFemtoCutMonitorCollection* AliFemtoCutMonitorHandler::PassMonitorColl() { return fPassColl;}
inline AliFemtoCutMonitorCollection* AliFemtoCutMonitorHandler::FailMonitorColl() { return fFailColl;}

#endif
