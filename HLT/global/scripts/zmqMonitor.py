#!/usr/bin/python
import zmq,sys

endpoint = "tcp://localhost:5556"
subscription = ""
if len(sys.argv)>1:
    subscription=sys.argv[1]
if len(sys.argv)>2:
    endpoint=sys.argv[2]

#  Prepare our context and sockets
context = zmq.Context()
socket = context.socket(zmq.SUB)
print "connect to: "+endpoint
socket.connect(str(endpoint))
print "subscribe to: "+subscription
socket.setsockopt(zmq.SUBSCRIBE, subscription)

while True:
  msg = socket.recv_multipart();
  print "topic: "+str(msg[0])
  print "messagesize: "+str(sys.getsizeof(msg[1]))
