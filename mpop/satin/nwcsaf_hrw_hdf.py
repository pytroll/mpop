"""Loader for MSG, nwcsaf high resolution hdf5 format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
from numpy import array as np_array
from numpy import empty as np_empty
from numpy import append as np_append
from numpy import dtype as np_dtype
from numpy import append as np_append
from numpy import where as np_where
from numpy import in1d as np_in1d
from numpy import logical_and as np_logical_and
from glob import glob
from mpop.projector import get_area_def
import datetime 

from copy import deepcopy

try:
    import h5py
except ImportError:
    print "... module h5py needs to be installed"
    quit()

from mipp.xrit.MSG import _Calibrator

import logging
LOG = logging.getLogger(__name__)
#from mpop.utils import debug_on
#debug_on()

GP_IDs = { 321: '08',  # Meteosat 8 
           322: '09',  # Meteosat 9 
           323: '10',  # Meteosat 10
           324: '11' } # Meteosat 11

dict_channel = {'CHANNEL00':'HRV',   'CHANNEL01':'VIS006','CHANNEL02':'VIS008','CHANNEL03':'IR_016','CHANNEL04':'IR_039','CHANNEL05':'WV_062',\
                'CHANNEL06':'WV_073','CHANNEL07':'IR_087','CHANNEL08':'IR_097','CHANNEL09':'IR_108','CHANNEL10':'IR_120','CHANNEL11':'IR_134'}


# class definition of a high resolution wind data
class HRW_class:

    def __init__(self):
        # see http://docs.scipy.org/doc/numpy/reference/generated/numpy.dtype.html
        self.date            = None # datetime of the observation  
        self.detailed        = None # False-> basic, True -> detailed 
        self.channel         = np_array([], dtype='|S6') 
        self.wind_id         = np_array([], dtype=int)
        self.prev_wind_id    = np_array([], dtype=int)
        self.segment_X       = np_array([], dtype='f') 
        self.segment_Y       = np_array([], dtype='f') 
        self.t_corr_method   = np_array([], dtype=int) 
        self.lon             = np_array([], dtype='f')     #  6.3862 [longitude in degree E]
        self.lat             = np_array([], dtype='f')     # 46.8823 [latitude in degree N]
        self.dlon            = np_array([], dtype='f')     # -0.011  [longitude in degree E]
        self.dlat            = np_array([], dtype='f')     #  0.01   [latitude in degree N]
        self.pressure        = np_array([], dtype='f')     # 64200.0 [p in Pa]
        self.wind_speed      = np_array([], dtype='f')     #   3.1   [v in m/s]
        self.wind_direction  = np_array([], dtype='f')     #  313.0  [v_dir in deg]
        self.temperature     = np_array([], dtype='f')     #  272.4  [T in K]
        self.conf_nwp        = np_array([], dtype='f')
        self.conf_no_nwp     = np_array([], dtype='f')
        self.t_type          = np_array([], dtype=int)
        self.t_level_method  = np_array([], dtype=int)
        self.t_winds         = np_array([], dtype=int)
        self.t_corr_test     = np_array([], dtype=int)
        self.applied_QI      = np_array([], dtype=int)
        self.NWP_wind_levels = np_array([], dtype=int)
        self.num_prev_winds  = np_array([], dtype=int)
        self.orographic_index= np_array([], dtype=int)
        self.cloud_type      = np_array([], dtype=int)
        self.wind_channel    = np_array([], dtype=int)
        self.correlation     = np_array([], dtype=int)
        self.pressure_error  = np_array([], dtype='f')

    # ---------------- add two data sets e.g. time steps ---------------------
    def __add__(self, HRW_class2):

        HRW_new = HRW_class()

        HRW_new.date            = self.date      # !!! does not make sense !!! 
        HRW_new.detailed        = self.detailed  # !!! does not make sense !!!
        HRW_new.channel         = np_append(self.channel,         HRW_class2.channel)
        HRW_new.wind_id         = np_append(self.wind_id,         HRW_class2.wind_id)
        HRW_new.prev_wind_id    = np_append(self.prev_wind_id,    HRW_class2.prev_wind_id)
        HRW_new.segment_X       = np_append(self.segment_X,       HRW_class2.segment_X)
        HRW_new.segment_Y       = np_append(self.segment_Y,       HRW_class2.segment_Y)
        HRW_new.t_corr_method   = np_append(self.t_corr_method,   HRW_class2.t_corr_method)
        HRW_new.lon             = np_append(self.lon,             HRW_class2.lon)
        HRW_new.lat             = np_append(self.lat,             HRW_class2.lat)
        HRW_new.dlon            = np_append(self.dlon,            HRW_class2.dlon)
        HRW_new.dlat            = np_append(self.dlat,            HRW_class2.dlat)
        HRW_new.pressure        = np_append(self.pressure,        HRW_class2.pressure)
        HRW_new.wind_speed      = np_append(self.wind_speed,      HRW_class2.wind_speed)
        HRW_new.wind_direction  = np_append(self.wind_direction,  HRW_class2.wind_direction)       
        HRW_new.temperature     = np_append(self.temperature,     HRW_class2.temperature)
        HRW_new.conf_nwp        = np_append(self.conf_nwp,        HRW_class2.conf_nwp)
        HRW_new.conf_no_nwp     = np_append(self.conf_no_nwp,     HRW_class2.conf_no_nwp)
        HRW_new.t_type          = np_append(self.t_type,          HRW_class2.t_type)
        HRW_new.t_level_method  = np_append(self.t_level_method,  HRW_class2.t_level_method)
        HRW_new.t_winds         = np_append(self.t_winds,         HRW_class2.t_winds)
        HRW_new.t_corr_test     = np_append(self.t_corr_test,     HRW_class2.t_corr_test)
        HRW_new.applied_QI      = np_append(self.applied_QI,      HRW_class2.applied_QI)
        HRW_new.NWP_wind_levels = np_append(self.NWP_wind_levels, HRW_class2.NWP_wind_levels)
        HRW_new.num_prev_winds  = np_append(self.num_prev_winds,  HRW_class2.num_prev_winds)
        HRW_new.orographic_index= np_append(self.orographic_index,HRW_class2.orographic_index)
        HRW_new.cloud_type      = np_append(self.cloud_type,      HRW_class2.cloud_type)
        HRW_new.wind_channel    = np_append(self.wind_channel,    HRW_class2.wind_channel)
        HRW_new.correlation     = np_append(self.correlation,     HRW_class2.correlation)
        HRW_new.pressure_error  = np_append(self.pressure_error,  HRW_class2.pressure_error)

        return HRW_new

    # ---------------- filter for certain criterias  ---------------------
    def filter(self, **kwargs):

        # if empty than return self (already empty)
        if self.channel.size == 0:
            return self

        HRW_new = deepcopy(self)

        for key_filter in ['min_correlation', 'min_conf_nwp', 'min_conf_no_nwp', 'cloud_type', 'level']:
            if key_filter in kwargs.keys():
                
                # if argument given is None or all keyword than skip this filter 
                if kwargs[key_filter] == None or kwargs[key_filter] == 'all' or kwargs[key_filter] == 'ALL' or kwargs[key_filter] == 'A':
                    continue

                n1 = str(HRW_new.channel.size)

                if key_filter == 'min_correlation':
                    inds = np_where(HRW_new.correlation > kwargs[key_filter])
                elif key_filter == 'min_conf_nwp':
                    inds = np_where(HRW_new.conf_nwp    > kwargs[key_filter])
                elif key_filter == 'min_conf_no_nwp':
                    inds = np_where(HRW_new.conf_no_nwp > kwargs[key_filter])
                elif key_filter == 'cloud_type':
                    mask = np_in1d(HRW_new.cloud_type, kwargs[key_filter]) 
                    inds = np_where(mask)[0]
                elif key_filter == 'level':
                    if kwargs[key_filter] == 'H': # high level: < 440hPa like in the ISCCP
                        inds = np_where(HRW_new.pressure < 44000 ) 
                    elif kwargs[key_filter] == 'M': # mid level: 440hPa ... 680hPa like in the ISCCP
                        inds = np_where( np_logical_and(44000 < HRW_new.pressure, HRW_new.pressure < 68000) ) 
                    elif kwargs[key_filter] == 'L': # low level: > 680hPa like in the ISCCP
                        inds = np_where(68000 < HRW_new.pressure)

                HRW_new.subset(inds)
                print "    filter for "+key_filter+" = ", kwargs[key_filter],' ('+n1+'->'+str(HRW_new.channel.size)+')'

        return HRW_new

    # ---------------- reduce HRW_dataset to the given indices inds ---------------------
    def subset(self, inds):

        self.channel          = self.channel         [inds]    
        self.wind_id          = self.wind_id         [inds]    
        self.prev_wind_id     = self.prev_wind_id    [inds]    
        self.segment_X        = self.segment_X       [inds]   
        self.segment_Y        = self.segment_Y       [inds]   
        self.t_corr_method    = self.t_corr_method   [inds]   
        self.lon              = self.lon             [inds]   
        self.lat              = self.lat             [inds]   
        self.dlon             = self.dlon            [inds]  
        self.dlat             = self.dlat            [inds]   
        self.pressure         = self.pressure        [inds]   
        self.wind_speed       = self.wind_speed      [inds]   
        self.wind_direction   = self.wind_direction  [inds]  
        self.temperature      = self.temperature     [inds]   
        self.conf_nwp         = self.conf_nwp        [inds]   
        self.conf_no_nwp      = self.conf_no_nwp     [inds]   
        self.t_type           = self.t_type          [inds]  
        self.t_level_method   = self.t_level_method  [inds]  
        self.t_winds          = self.t_winds         [inds] 
        self.t_corr_test      = self.t_corr_test     [inds]  
        self.applied_QI       = self.applied_QI      [inds]  
        self.NWP_wind_levels  = self.NWP_wind_levels [inds] 
        self.num_prev_winds   = self.num_prev_winds  [inds]
        self.orographic_index = self.orographic_index[inds]
        self.cloud_type       = self.cloud_type      [inds]
        self.wind_channel     = self.wind_channel    [inds]
        self.correlation      = self.correlation     [inds]
        self.pressure_error   = self.pressure_error  [inds]

        return self


def load(satscene, calibrate=True, area_extent=None, read_basic_or_detailed='both', **kwargs):
    """Load MSG SEVIRI High Resolution Wind (HRW) data from hdf5 format.
    """

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"orbit": satscene.orbit,
    "satname": satscene.satname,
    "number": satscene.number,
    "instrument": satscene.instrument_name,
    "satellite": satscene.fullname
    }

    LOG.info("assume seviri-level5")
    print "... assume seviri-level5"

    satscene.add_to_history("hdf5 data read by mpop/nwcsaf_hrw_hdf.py")

    # end of scan time 4 min after start 
    end_time = satscene.time_slot + datetime.timedelta(minutes=4)

    # area !!! satscene.area

    filename = os.path.join( satscene.time_slot.strftime(conf.get("seviri-level5", "dir", raw=True)),
                             satscene.time_slot.strftime(conf.get("seviri-level5", "filename", raw=True)) % values )

    # define classes before we search for files (in order to return empty class if no file is found)
    HRW_basic             = HRW_class()
    HRW_basic.detailed    = False 
    HRW_basic.date        = satscene.time_slot
    HRW_detailed          = HRW_class()
    HRW_detailed.detailed = True
    HRW_detailed.date     = satscene.time_slot

    print "... search for file: ", filename
    filenames=glob(str(filename))

    if len(filenames) != 0:

        if len(filenames) > 1:
            print "*** Warning, more than 1 datafile found: ", filenames 

        filename = filenames[0]
        print("... read data from %s" % str(filename))

        # create an instant of the HRW_class
        m_per_s_to_knots = 1.944

        ## limit channels to read 
        #hrw_channels=['HRV']
        # limit basic or detailed or both
        #read_basic_or_detailed='detailed'
        #read_basic_or_detailed='basic'


        with h5py.File(filename,'r') as hf:

            #print hf.attrs.keys()
            #print hf.attrs.values()

            region_name = hf.attrs['REGION_NAME'].replace("_", "")
            print "... read HRW data for region ", region_name
            LOG.info("... read HRW data for region "+region_name)
            sat_ID = GP_IDs[int(hf.attrs["GP_SC_ID"])]
            print "... derived from Meteosat ", sat_ID
            LOG.info("... derived from Meteosat "+sat_ID)

            # print('List of arrays in this file: \n', hf.keys()), len(hf.keys())

            if len(hf.keys()) == 0:
                print "*** Warning, empty file ", filename
                print ""
            else:
                for key in hf.keys():

                    if key[4:9] == "BASIC":
                        if 'read_basic_or_detailed' in locals():
                            if read_basic_or_detailed.lower() == "detailed":
                                continue
                        HRW_data = HRW_basic   # shallow copy 
                    elif key[4:12] == "DETAILED":
                        if 'read_basic_or_detailed' in locals():
                            if read_basic_or_detailed.lower() == "basic":
                                continue
                        HRW_data = HRW_detailed # shallow copy 

                    hrw_chn = dict_channel[key[len(key)-9:]]

                    if 'hrw_channels' in locals():
                        if hrw_channels != None:
                            if hrw_chn not in hrw_channels:
                                print "... "+hrw_chn+" is not in hrw_channels", hrw_channels 
                                print "    skip reading this channel" 
                                continue 

                    # read all  data 
                    channel = hf.get(key)
                    # print '... read wind vectors of channel ', channel.name, hrw_chn
                    # print  "  i    lon        lat      speed[kn] dir   pressure"
                    #for i in range(channel.len()):
                    #    print '%3d %10.7f %10.7f %7.2f %7.1f %8.1f' % (channel[i]['wind_id'], channel[i]['lon'], channel[i]['lat'], \
                    #                                                   channel[i]['wind_speed']*m_per_s_to_knots, \
                    #                                                   channel[i]['wind_direction'], channel[i]['pressure'])
                    # create string array with channel names 
                    channel_chararray = np_empty(channel.len(), dtype='|S6')
                    channel_chararray[:] = hrw_chn

                    HRW_data.channel          = np_append(HRW_data.channel         , channel_chararray              )
                    HRW_data.wind_id          = np_append(HRW_data.wind_id         , channel[:]['wind_id']          )    
                    HRW_data.prev_wind_id     = np_append(HRW_data.prev_wind_id    , channel[:]['prev_wind_id']     )    
                    HRW_data.segment_X        = np_append(HRW_data.segment_X       , channel[:]['segment_X']        )   
                    HRW_data.segment_Y        = np_append(HRW_data.segment_Y       , channel[:]['segment_Y']        )   
                    HRW_data.t_corr_method    = np_append(HRW_data.t_corr_method   , channel[:]['t_corr_method']    )   
                    HRW_data.lon              = np_append(HRW_data.lon             , channel[:]['lon']              )   
                    HRW_data.lat              = np_append(HRW_data.lat             , channel[:]['lat']              )   
                    HRW_data.dlon             = np_append(HRW_data.dlon            , channel[:]['dlon']             )  
                    HRW_data.dlat             = np_append(HRW_data.dlat            , channel[:]['dlat']             )   
                    HRW_data.pressure         = np_append(HRW_data.pressure        , channel[:]['pressure']         )   
                    HRW_data.wind_speed       = np_append(HRW_data.wind_speed      , channel[:]['wind_speed']       )   
                    HRW_data.wind_direction   = np_append(HRW_data.wind_direction  , channel[:]['wind_direction']   )  
                    HRW_data.temperature      = np_append(HRW_data.temperature     , channel[:]['temperature']      )   
                    HRW_data.conf_nwp         = np_append(HRW_data.conf_nwp        , channel[:]['conf_nwp']         )   
                    HRW_data.conf_no_nwp      = np_append(HRW_data.conf_no_nwp     , channel[:]['conf_no_nwp']      )   
                    HRW_data.t_type           = np_append(HRW_data.t_type          , channel[:]['t_type']           )  
                    HRW_data.t_level_method   = np_append(HRW_data.t_level_method  , channel[:]['t_level_method']   )  
                    HRW_data.t_winds          = np_append(HRW_data.t_winds         , channel[:]['t_winds']          ) 
                    HRW_data.t_corr_test      = np_append(HRW_data.t_corr_test     , channel[:]['t_corr_test']      )   
                    HRW_data.applied_QI       = np_append(HRW_data.applied_QI      , channel[:]['applied_QI']       )  
                    HRW_data.NWP_wind_levels  = np_append(HRW_data.NWP_wind_levels , channel[:]['NWP_wind_levels']  ) 
                    HRW_data.num_prev_winds   = np_append(HRW_data.num_prev_winds  , channel[:]['num_prev_winds']   )
                    HRW_data.orographic_index = np_append(HRW_data.orographic_index, channel[:]['orographic_index'] )
                    HRW_data.cloud_type       = np_append(HRW_data.cloud_type      , channel[:]['cloud_type']       )
                    HRW_data.wind_channel     = np_append(HRW_data.wind_channel    , channel[:]['wind_channel']     )
                    HRW_data.correlation      = np_append(HRW_data.correlation     , channel[:]['correlation']      )
                    HRW_data.pressure_error   = np_append(HRW_data.pressure_error  , channel[:]['pressure_error']   )

                # sort according to wind_id
                inds = HRW_data.wind_id.argsort()
                HRW_data.subset(inds) # changes HRW_data itself

                # sorting without conversion to numpy arrays 
                #[e for (wid,pwid) in sorted(zip(HRW_data.wind_id,HRW_data.prev_wind_id))]

    else:
        print "*** Error, no file found"
        print ""
        sat_ID = "no file"
        # but we continue the program in order to add an empty channel below 


    ## filter data according to the given optional arguments 
    #n1 = str(HRW_data.channel.size)
    #HRW_data = HRW_data.filter(**kwargs)   
    #print "    apply filters "+' ('+n1+'->'+str(HRW_data.channel.size)+')'

    chn_name="HRW"
    satscene[chn_name].HRW_basic    = HRW_basic.filter(**kwargs)     # returns new object (deepcopied and filtered)
    satscene[chn_name].HRW_detailed = HRW_detailed.filter(**kwargs)  # returns new object (deepcopied and filtered)
    satscene[chn_name].info['units'] = 'm/s'
    satscene[chn_name].info['satname'] = 'meteosat'
    satscene[chn_name].info['satnumber'] = sat_ID
    satscene[chn_name].info['instrument_name'] = 'seviri'
    satscene[chn_name].info['time'] = satscene.time_slot
    satscene[chn_name].info['is_calibrated'] = True
