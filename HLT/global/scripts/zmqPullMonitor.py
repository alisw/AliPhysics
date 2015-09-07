#!/usr/bin/python
import zmq,sys

endpoint = "tcp://*:60202"
subscription = ""
if len(sys.argv)>1:
    subscription=sys.argv[1]
if len(sys.argv)>2:
    endpoint=sys.argv[2]

#  Prepare our context and sockets
context = zmq.Context()
socket = context.socket(zmq.PULL)
print "connect to: "+endpoint
socket.bind(str(endpoint))
print "subscribe to: "+subscription
#socket.setsockopt(zmq.SUBSCRIBE, subscription)

while True:
  msg = socket.recv_multipart();
  print "_________________________"
  i=0;
  for message in msg:
    if i==0:
      print "topic: "+str(message)
      i=1
    elif i==1:
      print "  messagesize: "+str(sys.getsizeof(message))
      i=0
