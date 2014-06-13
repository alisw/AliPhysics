// Main authors: Mihai Niculescu 2014

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliThreadedSocket_H
#define AliThreadedSocket_H

#include <TQObject.h>

class TThread;
class AliNetMessage;

namespace zmq {
	class context_t;
}

class AliThreadedSocket : public TQObject
{
public:
	enum EOpenMode{READ, WRITE};

	AliThreadedSocket(zmq::context_t *context, EOpenMode mode);
	virtual ~AliThreadedSocket();

	Bool_t Start();
	Bool_t Stop();
	Bool_t Kill();
	
	zmq::context_t* GetContext() const;
	TThread* GetThread() const;
		
	void Started(); // *SIGNAL*
	void Stopped(); // *SIGNAL*

	void Continue();
	
protected:
  AliThreadedSocket(const AliThreadedSocket&);            // Not implemented
  AliThreadedSocket& operator=(const AliThreadedSocket&); // Not implemented

	// reimplement these in a derived class
	static void* RunThrdRead(void* arg);
	static void* RunThrdWrite(void* arg);

	zmq::context_t* fContext;
	TThread* fThread;
	EOpenMode fOpenMode;


  ClassDef(AliThreadedSocket, 0);  

};
#endif
