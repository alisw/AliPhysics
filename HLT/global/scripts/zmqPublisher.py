#!/usr/bin/python
import zmq,time,sys

port="60201"
topic="ALIESDV0HLT"
if (len(sys.argv)>1):
  topic=sys.argv[1]
if (len(sys.argv)>2):
  port=sys.argv[2]

#  Prepare our context and sockets
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind('tcp://*:'+port)
while True:
  print "publishing on "+port
  socket.send(topic,zmq.SNDMORE)
  socket.send("payload part",0)
  time.sleep(1)

