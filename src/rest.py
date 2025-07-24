#!/usr/bin/env python3

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
	opts.setBackgroundColour((0.06274509804, 0.07450980392, 0.09411764706))
	opts.setAtomPalette({
		-1 : (0.8, 0.8, 0.8),
		0  : (0.9, 0.9, 0.9),
		1  : (0.625, 0.625, 0.625),
		6  : (0.9, 0.9, 0.9),
		7  : (0.5, 0.6, 1.0),
		8  : (1.0, 0.3, 0.3),
		9  : (0.1, 0.9, 0.9),
		15 : (1.0, 0.5, 0.0),
		16 : (1.0, 0.8, 0.3),
		17 : (0.0, 0.802, 0.0),
		35 : (0.71, 0.4, 0.07),
		53 : (0.89, 0.004, 1),
		201: (0.68, 0.85, 0.90),
	})
	#opts.addStereoAnnotation = True

def get_molecule(data=None, h=False):
	mol = Chem.MolFromSmiles(data)
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
	def svg(self, data=None, h=False, size=(-1,-1)):
		mol = get_molecule(data, h)

		drawer = rdMolDraw2D.MolDraw2DSVG(*size)
		set_options(drawer)
		drawer.DrawMolecule(mol)
		drawer.FinishDrawing()
		svg = drawer.GetDrawingText()

		cherrypy.response.headers['Content-Type'] = "image/svg+xml"
		return svg.encode('utf8')

	@cherrypy.expose
	def png(self, data=None, h=False, size=(300,300)):
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
		cherrypy.response.headers['Content-Type'] = "image/png"
		return png

	@cherrypy.expose
	def render(self, name=None, format=None, h=None, size=None):
		if h is None or len(h) == 0 or h == "0" or h == "False" or h == "false":
			h = False
		else:
			h = True

		if size is None or size == "":
			size = (-1,-1)
		else:
			size = (int(w) for w in size.split("x"))

		if format is None:
			data, input_format = name.split(".")
		else:
			data = name
			input_format = format.strip('.')

		if input_format == "png":
			if size[0] > 1000 or size[1] > 1000:
				size = (300,300)
			return self.png(data, h, size)
		elif input_format == "svg":
			return self.svg(data, h, size)

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
