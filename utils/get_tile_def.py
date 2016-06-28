#!/usr/bin/python

import xml.etree.ElementTree as ET
from pyresample import utils
import pickle
import urllib2

length=109800

#https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_20150622T000000_21000101T000000_ZZ_0001
FNAME="S2A_OPER_GIP_TILPAR_20150622T000000_21000101T000000_ZZ_0001"

TILPAR_URL="https://sentinel.esa.int/documents/247904/1955685/"+FNAME

FNAME=FNAME+".kml"

tiles = urllib2.urlopen(TILPAR_URL)
with open(FNAME,'wb') as output:
  output.write(tiles.read())

tiles.close()

tree = ET.parse(FNAME)
root = tree.getroot()

s2tiles={}

for pm in root.iter('{http://www.opengis.net/kml/2.2}Placemark'):
	tilename=None
	epsg=None
	utm_ul_x=None
	utm_ul_y=None

	for name in pm.iter('{http://www.opengis.net/kml/2.2}name'):
		tilename=name.text
	for simple in pm.iter('{http://www.opengis.net/kml/2.2}SimpleData'):
		if (simple.attrib['name']=='epsg'):
			epsg=simple.text
		if(simple.attrib['name']=='utm_ul_x'):
			utm_ul_x=simple.text
		if(simple.attrib['name']=='utm_ul_y'):
			utm_ul_y=simple.text

	extent=(float(utm_ul_x),float(utm_ul_y)-length,float(utm_ul_x)+length,float(utm_ul_y))

	s2tiles[tilename]=[epsg,extent]

f=open('s2tiles.pickle','w')
pickle.dump(s2tiles,f)
f.close()

