#! /usr/bin/env python

import os
import sys
import jinja2

def render_matchbox(matchbox_template_path, context):
  path, filename = os.path.split(matchbox_template_path)
  return jinja2.Environment(
    loader=jinja2.FileSystemLoader(path or './')
  ).get_template(filename).render(context)

def generate_matchbox_from_template(templatename):
  matchboxname = templatename
  matchboxname = matchboxname.replace(".tmpl", ".in")
  # decode configuration parameters from environment variables set by AliDPG
  configparams = {}
  for k, v in os.environ.items():
    if "CONFIG" in k:
      configparams[k] = v

  with open(matchboxname, 'w') as matchboxwriter:
    rendered = render_matchbox(templatename, configparams)
    #print(rendered)
    matchboxwriter.write(rendered)

if __name__ == "__main__":
  generate_matchbox_from_template(sys.argv[1])  
