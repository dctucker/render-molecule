#!/usr/bin/env python3

from config import Config
import cherrypy
import simplejson
import json
from rdkit import Chem
from rdkit.Chem import Draw, rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
import PIL.ImageOps
import PIL.ImageEnhance
import numpy as np

import tempfile

#from main import *

#InChI=1S/C21H30O2/c1-5-6-7-8-15-12-18(22)20-16-11-14(2)9-10-17(16)21(3,4)23-19(20)13-15/h11-13,16-17,22H,5-10H2,1-4H3/t16-,17-/m1/s1

def set_options(drawer):
	opts = drawer.drawOptions()
	opts.additionalAtomLabelPadding = 0.125
	opts.setBackgroundColour(config.getcolor('bg'))
	opts.setAtomPalette({n: config.getcolor(n) for n in [-1, 0, 1, 6, 7, 8, 9, 15, 16, 17, 35, 53, 201]})
	#opts.addStereoAnnotation = True

def get_molecule(data=None, h=False):
	mol = Chem.MolFromSmiles(data)
	if mol is None:
		raise cherrypy.HTTPError(status=400, message='Syntax error')

	mol = Chem.Mol(mol)

	Chem.rdDepictor.SetPreferCoordGen(False)
	if h:
		mol = Chem.AddHs(mol)

	if not mol.GetNumConformers():
		Chem.rdDepictor.Compute2DCoords(mol, useRingTemplates=True)

	Chem.rdDepictor.StraightenDepiction(mol)
	#Chem.rdMolTransforms.CanonicalizeMol(mol)
	#Chem.rdDepictor.NormalizeDepiction(mol, -1, 0)

	#Draw.MolsToGridImage([mol, moltmp])
	return mol

class SmileRenderer(object):

	@cherrypy.expose
	def svg(self, data=None, h=False, size=None):
		if size is None:
			size = (-1, -1)

		mol = get_molecule(data, h)

		drawer = rdMolDraw2D.MolDraw2DSVG(*size)
		set_options(drawer)
		drawer.DrawMolecule(mol)
		drawer.FinishDrawing()
		svg = drawer.GetDrawingText()

		cherrypy.response.headers['Content-Type'] = 'image/svg+xml'
		return svg.encode('utf8')

	@cherrypy.expose
	def png(self, data=None, h=False, size=None):
		default_size = config.getsize('default')
		max_size = config.getsize('max')
		if size is None or size[0] > max_size[0] or size[1] > max_size[1]:
			size = default_size

		mol = get_molecule(data, h)

		drawer = rdMolDraw2D.MolDraw2DCairo(*size)
		set_options(drawer)
		drawer.DrawMolecule(mol)
		drawer.FinishDrawing()
		png = drawer.GetDrawingText()

		buffer = io.BytesIO()
		buffer.write(png)
		buffer.seek(0)
		png = buffer.read()
		cherrypy.response.headers['Content-Type'] = 'image/png'
		return png

	@cherrypy.expose
	def render(self, name=None, format=None, h=None, size=None):
		if h is None or len(h) == 0 or h == '0' or h == 'False' or h == 'false':
			h = False
		else:
			h = True

		if size is not None and len(size) > 0:
			size = tuple(int(w) for w in size.split('x'))
		else:
			size = None

		data = name
		if format is None:
			output_format = config.get('server', 'format')
		else:
			output_format = format.strip('.')

		if output_format not in ('svg', 'png'):
			raise cherrypy.HTTPError(status=400, message='Unknown output format')

		if output_format == 'png':
			return self.png(data, h, size)
		elif output_format == 'svg':
			return self.svg(data, h, size)

	@cherrypy.expose
	def index(self):
		return "Render Molecule API\n"

if __name__ == '__main__':
	global config
	config = Config()
	debug = config.getboolean('server', 'debug')
	conf = {
			'/': {
				#'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
				#'tools.sessions.on': True,
				#'tools.response_headers.on': True,
				#'tools.response_headers.headers': [('Content-Type', 'application/json')],
				'request.show_tracebacks': debug,
			}
		}
	cherrypy.config.update({
		'server.socket_host': config.get('server', 'host'),
		'server.socket_port': config.getint('server', 'port'),
	})
	cherrypy.quickstart( SmileRenderer(), config.get('server', 'path'), conf )
