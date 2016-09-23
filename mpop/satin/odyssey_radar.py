import Image
import glob
import os
from ConfigParser import ConfigParser
import numpy as np
import numpy.ma as ma
from mpop import CONFIG_PATH

import pyresample
import logging

import h5py

LOG = logging.getLogger(__name__)

ODIM_H5_FIELD_NAMES = {
   'TH': 'total_power',      # uncorrected reflectivity, horizontal
   'TV': 'total_power',      # uncorrected reflectivity, vertical
   'DBZH': 'reflectivity',    # corrected reflectivity, horizontal
   'DBZV': 'reflectivity',    # corrected reflectivity, vertical
   'ZDR': 'differential_reflectivity',    # differential reflectivity
   'RHOHV': 'cross_correlation_ratio',
   'LDR': 'linear_polarization_ratio',
   'PHIDP': 'differential_phase',
   'KDP': 'specific_differential_phase',
   'SQI': 'normalized_coherent_power',
   'SNR': 'signal_to_noise_ratio',
   'VRAD': 'velocity',
   'WRAD': 'spectrum_width',
   'QIND': 'quality_index',
   'RATE': 'precip',         # precip
   'ACRR': 'accu_precip',      # 1 hour ACCU
}


def load(satscene, *args, **kwargs):
   """Loads the *channels* into the satellite *scene*.
   """
   #
   # Dataset information
   #
   # Read config file content
   conf = ConfigParser()
   conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

   values = {"orbit": satscene.orbit,
          "satname": satscene.satname,
          "number": satscene.number,
          "instrument": satscene.instrument_name,
          "satellite": satscene.fullname
          }

   # projection info
   projectionName = conf.get("radar-level2", "projection")
   projection = pyresample.utils.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)
   satscene.area = projection
   
   for chn_name in satscene.channels_to_load:
      filename = os.path.join(
         satscene.time_slot.strftime(conf.get("radar-level2", "dir", raw=True)) % values,
         satscene.time_slot.strftime(conf.get(chn_name,  "filename", raw=True)) % values )

      # Load data from the h5 file
      LOG.debug("filename: "+filename)
      filenames=glob.glob(str(filename))

      if len(filenames) == 0:
         LOG.debug("no input file found: "+filename)
         print "no input file found:"+filename
         quit()
      else:
         filename = glob.glob(str(filename))[0]
      
      # open the file
      hfile = h5py.File(filename, 'r')
      odim_object = hfile['what'].attrs['object']
      if odim_object != 'COMP':
         raise NotImplementedError('object: %s not implemented.' % (odim_object))
      else:
         # File structure
         
         #>>> hfile.keys()
         #[u'dataset1', u'dataset2', u'how', u'what', u'where']


         #>>> for f in hfile['what'].attrs.keys():
         #...  print "hfile['what'].attrs['",f,"']=", hfile['what'].attrs[f]
         #
         #hfile['what'].attrs[' object ']= COMP
         #hfile['what'].attrs[' version ']= H5rad 2.0
         #hfile['what'].attrs[' date ']= 20151201
         #hfile['what'].attrs[' time ']= 060000
         #hfile['what'].attrs[' source ']= ORG:247

         #>>> for f in hfile['where'].attrs.keys():
         #...  print "hfile['where'].attrs['",f,"']=", hfile['where'].attrs[f]
         #
         #hfile['where'].attrs[' projdef ']= +proj=laea +lat_0=55.0 +lon_0=10.0 +x_0=1950000.0 +y_0=-2100000.0 +units=m +ellps=WGS84
         #hfile['where'].attrs[' xsize ']= 1900
         #hfile['where'].attrs[' ysize ']= 2200
         #hfile['where'].attrs[' xscale ']= 2000.0
         #hfile['where'].attrs[' yscale ']= 2000.0
         #hfile['where'].attrs[' LL_lon ']= -10.4345768386
         #hfile['where'].attrs[' LL_lat ']= 31.7462153193
         #hfile['where'].attrs[' UL_lon ']= -39.5357864125
         #hfile['where'].attrs[' UL_lat ']= 67.0228327583
         #hfile['where'].attrs[' UR_lon ']= 57.8119647501
         #hfile['where'].attrs[' UR_lat ']= 67.6210371028
         #hfile['where'].attrs[' LR_lon ']= 29.4210386356
         #hfile['where'].attrs[' LR_lat ']= 31.9876502779

         # hfile['how'].attrs['nodes'] 
         # list of radar in composite

         #>>> for f in hfile['dataset1']['what'].attrs.keys():
         #...  print "hfile['dataset1'][what].attrs['",f,"']=", hfile['dataset1']['what'].attrs[f]
         #
         #hfile['dataset1'][what].attrs[' product ']= COMP
         #hfile['dataset1'][what].attrs[' startdate ']= 20151201
         #hfile['dataset1'][what].attrs[' starttime ']= 055000
         #hfile['dataset1'][what].attrs[' enddate ']= 20151201
         #hfile['dataset1'][what].attrs[' endtime ']= 060500
         #hfile['dataset1'][what].attrs[' quantity ']= RATE
         #hfile['dataset1'][what].attrs[' gain ']= 1.0
         #hfile['dataset1'][what].attrs[' offset ']= 0.0
         #hfile['dataset1'][what].attrs[' nodata ']= -9999000.0
         #hfile['dataset1'][what].attrs[' undetect ']= -8888000.0
         #>>> for f in hfile['dataset2']['what'].attrs.keys():
         #...  print "hfile['dataset2'][what].attrs['",f,"']=", hfile['dataset2']['what'].attrs[f]
         #
         #hfile['dataset2'][what].attrs[' product ']= COMP
         #hfile['dataset2'][what].attrs[' startdate ']= 20151201
         #hfile['dataset2'][what].attrs[' starttime ']= 055000
         #hfile['dataset2'][what].attrs[' enddate ']= 20151201
         #hfile['dataset2'][what].attrs[' endtime ']= 060500
         #hfile['dataset2'][what].attrs[' quantity ']= QIND
         #hfile['dataset2'][what].attrs[' gain ']= 1.0
         #hfile['dataset2'][what].attrs[' offset ']= 0.0
         #hfile['dataset2'][what].attrs[' nodata ']= -9999000.0
         #hfile['dataset2'][what].attrs[' undetect ']= -8888000.0

         _xsize = hfile['where'].attrs['xsize']
         _ysize = hfile['where'].attrs['ysize']
         #nbins= _xsize * _ysize

         #projection = hfile['where'].attrs['projdef']
         
         datasets = [k for k in hfile if k.startswith('dataset')]
         datasets.sort()
         nsweeps = len(datasets)
         
         try:
            ds1_what = hfile[datasets[0]]['what'].attrs
         except KeyError:
            # if no how group exists mock it with an empty dictionary
            ds1_what = {}
         
         _type = ''
         if 'product' in ds1_what:
            LOG.debug("product: "+ds1_what['product'])
            if ds1_what['product'] == 'COMP':
               if 'quantity' in ds1_what:
                  _type = ds1_what['quantity']
                  LOG.debug("product_type: "+_type)

                  #for chn_name in satscene.channels_to_load:
                  #   if chn_name == _type:

                  raw_data = hfile[datasets[0]]['data1']['data'][:]
                  raw_data = raw_data.reshape(_ysize,_xsize)
         
                  # flag no data
                  if 'nodata' in ds1_what:
                     nodata = ds1_what['nodata']
                     data = np.ma.masked_equal(raw_data, nodata)
                  else:
                     data = np.ma.masked_array(raw_data)
         
                  mask = np.ma.masked_array( raw_data == nodata )
                  mask = np.ma.masked_equal( mask, False)
            
                  # flag undetect data 
                  if 'undetect' in ds1_what:
                     undetect = ds1_what['undetect']
                     data[data == undetect] = np.ma.masked
                        
                  #from trollimage.image import Image as trollimage
                  #img = trollimage(mask, mode="L", fill_value=[1,1,1]) # [0,0,0] [1,1,1]
                  #from trollimage.colormap import rainbow
                  #img.colorize(rainbow)
                  #img.show()
                  #quit()

                  # gain/offset adjustment
                  if 'offset' in ds1_what:
                     offset = ds1_what['offset']
                  else:
                     offset = 0.0
                     
                  if 'gain' in ds1_what:
                     gain = ds1_what['gain']
                  else:
                     gain = 1.0

                  data *= gain + offset
                  
                  satscene[chn_name] = data
                  satscene[chn_name+'-MASK'] = mask

                  LOG.debug(" *** channel:"+chn_name)
                  
                  if _type == 'DBZH':
                     units = 'dBZ'
                  
                  if _type == 'RATE':
                     units = 'mm/h'

                  if _type == 'ACRR':
                     units = 'mm'
                     
                  satscene[chn_name].info["units"] = units
                  LOG.debug("channel:"+chn_name+" units:"+units)
