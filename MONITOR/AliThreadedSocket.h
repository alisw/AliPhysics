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
	enum EOpenMode{kREAD, kWRITE};

	AliThreadedSocket(zmq::context_t *context, EOpenMode mode);
	virtual ~AliThreadedSocket();

	Bool_t Start();
	Bool_t Stop();
	Bool_t Kill();
	void Wait();
	
	zmq::context_t* GetContext() const;
	TThread* GetThread() const;
	EOpenMode GetMode() const { return fOpenMode; }
		
	void Started(); // *SIGNAL*
	void Stopped(); // *SIGNAL*

	void Continue();
	
protected:
  AliThreadedSocket(const AliThreadedSocket&);            // Not implemented
  AliThreadedSocket& operator=(const AliThreadedSocket&); // Not implemented

	// Reimplement these in a derived Class
	virtual void RunThrdRead(); // function to run in a thread when in Read mode
	virtual void RunThrdWrite(); // function to run in a thread when in Write mode

	TThread* fThread;
	zmq::context_t* fContext;
	EOpenMode fOpenMode;

private:
	static void* Dispatch(void* arg)
	{
		AliThreadedSocket* th = static_cast<AliThreadedSocket*>(arg);
		
		if(th->GetMode()==kREAD)
			th->RunThrdRead();
		else
			th->RunThrdWrite();
			
			return NULL;
	}

  ClassDef(AliThreadedSocket, 0);  
};
#endif
