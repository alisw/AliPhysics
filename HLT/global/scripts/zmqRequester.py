#!/usr/bin/python
import zmq
import sys

request=b"ALIESD"
endpoint="tcp://localhost:60200"
if len(sys.argv)>1:
  request=sys.argv[1];
if len(sys.argv)>2:
  endpoint=sys.argv[2];

#  Prepare our context and sockets
context = zmq.Context()
socket = context.socket(zmq.REQ)
socket.connect(endpoint)

#send request
socket.send(request)
print("sent request: "+request)

#receive reply (multipart)
msg = socket.recv_multipart();
print "topic: "+str(msg[0])
print "messagesize: "+str(sys.getsizeof(msg[1]))
