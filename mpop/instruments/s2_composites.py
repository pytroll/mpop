from mpop.imageo.geo_image import GeoImage

def s2_truecolor(self):
	
	self.check_channels('B02','B03','B04')
    
	ch1 = self['B04'].data
	ch2 = self['B03'].data
	ch3 = self['B02'].data

	img = GeoImage((ch1, ch2, ch3),
				self.area,
                                 self.time_slot,
                                 fill_value=None,
                                 mode="RGB")

	img.enhance(stretch="linear")
	#img.enhance(stretch="histogram")
	img.enhance(gamma=2.0)
    
	return img

s2_truecolor.prerequisites = set(['B02', 'B03','B04'])
msi=[s2_truecolor]
