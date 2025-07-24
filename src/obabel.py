#!/usr/bin/env python3

import cherrypy
import simplejson
import json
from openbabel import pybel
import tempfile

#from main import *

#InChI=1S/C21H30O2/c1-5-6-7-8-15-12-18(22)20-16-11-14(2)9-10-17(16)21(3,4)23-19(20)13-15/h11-13,16-17,22H,5-10H2,1-4H3/t16-,17-/m1/s1

class SmileRenderer(object):

	@cherrypy.expose
	def svg(self, data=None):
		#data = cherrypy.request.json
		#accept = cherrypy.request.headers['Accept'] 
		#input_format = cherrypy.request.headers['Content-Type']

		#if input_format != 'chemical/x-daylight-smiles':
		#	raise cherrypy.HTTPError(status=415)
		if data is None:
			data = cherrypy.request.params.get("data")
		if data[0:6] == "InChI=":
			input_format = "inchi"
		else:
			input_format = "smi"
		mol = pybel.readstring(input_format, data)
		svg = mol.write('svg')
		cherrypy.response.headers['Content-Type'] = "image/svg+xml"
		return svg.encode('utf8')

	@cherrypy.expose
	def png(self, smiles=None):
		mol = pybel.readstring('smi', smiles)
		cherrypy.response.headers['Content-Type'] = "image/png"

		with tempfile.NamedTemporaryFile(prefix="render-molecule-", suffix="png") as fp:
			mol.write("_png2", fp.name, overwrite=True)
			with open(fp.name, mode="rb") as f:
				return(f.read())

	@cherrypy.expose
	def index(self):
		return "Render Molecule API\n"

if __name__ == '__main__':
	conf = {
			'/': {
				#'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
				#'tools.sessions.on': True,
				#'tools.response_headers.on': True,
				#'tools.response_headers.headers': [('Content-Type','application/json')],
			}
		}
	cherrypy.config.update({
		"server.socket_port": 8555,
	})
	cherrypy.quickstart( SmileRenderer(), '/', conf )
