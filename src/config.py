import os
import configparser

def hex2tuple(color: str) -> list[float,float,float]:
	color = color.strip('#')
	if len(color) == 3:
		return tuple((int(s, 16) * 0x11) / 255.0 for s in color)
	elif len(color) == 6:
		return tuple(int(color[i:i+2], 16) / 255.0 for i in range(0, len(color), 2))
	return None

class Config(configparser.ConfigParser):
	def __init__(self):
		super().__init__(defaults={
			'host': '127.0.0.1',
			'port': '8555',
			'format': 'svg',
			'path': '/',
			'debug': False,
		}, default_section='server')
		self.read_dict({
			'size': {
				'default': '300x300',
				'max': '1000x1000',
			},
			'colors': {
				'bg': '#101218',
				'fg': '#cccccc',
				0   : '#e6e6e6',
				1   : '#999999',
				6   : '#e6e6e6',
				7   : '#7f99ff',
				8   : '#ff4e4e',
				9   : '#1ae6e6',
				15  : '#ff7f00',
				16  : '#ffcc4d',
				17  : '#00cd00',
				35  : '#b66612',
				53  : '#e301ff',
				201 : '#aed9e5',
			}
		})
		self.read([
			os.getenv('MOLECULE_CFG', 'molecule.cfg'),
			os.path.expanduser('~/.config/molecule.cfg'),
			'/etc/molecule.cfg',
		])

	def getsize(self, key: str) -> tuple[int,int]:
		return tuple(int(d) for d in self.get('size', key).split('x'))

	def getcolor(self, key: str) -> tuple[float,float,float]:
		return hex2tuple(self.get('colors', str(key), fallback=self.get('colors','fg')))

