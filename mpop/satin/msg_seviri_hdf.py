"""Loader for MSG, netcdf format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
import numpy.ma as ma
from numpy import array as np_array
from numpy import nan as np_nan
from glob import glob
from mpop.projector import get_area_def
import datetime 

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

SatelliteIds = { '08': 321,  # Meteosat 8 
                  '8': 321,  # Meteosat 8 
                 '09': 322,  # Meteosat 9 
                  '9': 322,  # Meteosat 9 
                 '10': 323,  # Meteosat 10
                 '11': 324 } # Meteosat 11

channel_numbers = {"VIS006": 1,
                   "VIS008": 2,
                   "IR_016": 3,
                   "IR_039": 4,
                   "WV_062": 5,
                   "WV_073": 6,
                   "IR_087": 7,
                   "IR_097": 8,
                   "IR_108": 9,
                   "IR_120": 10,
                   "IR_134": 11,
                   "HRV": 12}

dict_channel= {'VIS006':'Channel 01','VIS008':'Channel 02','IR_016':'Channel 03','IR_039':'Channel 04','WV_062':'Channel 05','WV_073':'Channel 06',\
               'IR_087':'Channel 07','IR_097':'Channel 08','IR_108':'Channel 09','IR_120':'Channel 10','IR_134':'Channel 11','HRV':'Channel 12'}


def load(satscene, calibrate=True, area_extent=None, **kwargs):
    """Load MSG SEVIRI data from hdf5 format.
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

    LOG.info("assume seviri-level4")
    print "... assume seviri-level4"

    satscene.add_to_history("hdf5 data read by mpop/msg_seviri_hdf.py")


    if "reader_level" in kwargs.keys():
        reader_level = kwargs["reader_level"]
    else:
        reader_level = "seviri-level4"

    if "RSS" in kwargs.keys():
        if kwargs["RSS"]:
            dt_end =  4
        else:
            dt_end = 12
    else:
        from my_msg_module import check_RSS
        RSS = check_RSS(satscene.sat_nr(), satscene.time_slot)
        if RSS == None:
            print "*** Error in mpop/satin/msg_seviri_hdf.py"
            print "    satellite MSG", satscene.sat_nr() ," is not active yet"
            quit()
        else:
            if RSS:
                dt_end =  4
            else:
                dt_end = 12

    print "... hdf file name is specified by observation end time"
    print "    assume ", dt_end, " min between start and end time of observation"

    # end of scan time 4 min after start 
    end_time = satscene.time_slot + datetime.timedelta(minutes=dt_end)

    filename = os.path.join( end_time.strftime(conf.get(reader_level, "dir", raw=True)),
                             end_time.strftime(conf.get(reader_level, "filename", raw=True)) % values )
    
    print "... search for file: ", filename
    filenames=glob(str(filename))
    if len(filenames) == 0:
        print "*** Error, no file found"
        return # just return without exit the program 
    elif len(filenames) > 1:
        print "*** Warning, more than 1 datafile found: ", filenames 
    filename = filenames[0]
    print("... read data from %s" % str(filename))

    # read data from hdf5 file 
    data_folder='U-MARF/MSG/Level1.5/'

    # Load data from hdf file
    with h5py.File(filename,'r') as hf:

        subset_info=hf.get(data_folder+'METADATA/SUBSET')
        for i in range(subset_info.len()):
            #print subset_info[i]['EntryName'], subset_info[i]['Value']
            if subset_info[i]['EntryName'] == "VIS_IRSouthLineSelectedRectangle":
                VIS_IRSouthLine = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "VIS_IRNorthLineSelectedRectangle":
                VIS_IRNorthLine = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "VIS_IREastColumnSelectedRectangle":
                VIS_IREastColumn = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "VIS_IRWestColumnSelectedRectangle":
                VIS_IRWestColumn = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "HRVLowerNorthLineSelectedRectangle":
                HRVLowerNorthLine = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "HRVLowerSouthLineSelectedRectangle":
                HRVLowerSouthLine = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "HRVLowerEastColumnSelectedRectangle":
                HRVLowerEastColumn = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "HRVLowerWestColumnSelectedRectangle":
                HRVLowerWestColumn = int(subset_info[i]['Value'])
            if subset_info[i]['EntryName'] == "HRVUpperSouthLineSelectedRectangle":
                HRVUpperSouthLine = int(subset_info[i]['Value'])  # 0
            if subset_info[i]['EntryName'] == "HRVUpperNorthLineSelectedRectangle":
                HRVUpperNorthLine = int(subset_info[i]['Value'])  # 0
            if subset_info[i]['EntryName'] == "HRVUpperEastColumnSelectedRectangle":
                HRVUpperEastColumn = int(subset_info[i]['Value']) # 0
            if subset_info[i]['EntryName'] == "HRVUpperWestColumnSelectedRectangle":
                HRVUpperWestColumn = int(subset_info[i]['Value']) # 0

        sat_status=hf.get(data_folder+'METADATA/HEADER/SatelliteStatus/SatelliteStatus_DESCR')
        for i in range(subset_info.len()):
            if sat_status[i]['EntryName']=="SatelliteDefinition-NominalLongitude":
                sat_lon = sat_status[i]['Value']
                break

        #print 'VIS_IRSouthLine', VIS_IRSouthLine
        #print 'VIS_IRNorthLine', VIS_IRNorthLine
        #print 'VIS_IREastColumn', VIS_IREastColumn
        #print 'VIS_IRWestColumn', VIS_IRWestColumn
        #print 'sat_longitude', sat_lon, type(sat_lon), 'GEOS<'+'{:+06.1f}'.format(sat_lon)+'>' 

        if 1 == 0:
            # works only if all pixels are on the disk 
            from msg_pixcoord2area import msg_pixcoord2area
            print "VIS_IRNorthLine, VIS_IRWestColumn, VIS_IRSouthLine, VIS_IREastColumn: ", VIS_IRNorthLine, VIS_IRWestColumn, VIS_IRSouthLine, VIS_IREastColumn
            area_def = msg_pixcoord2area ( VIS_IRNorthLine, VIS_IRWestColumn, VIS_IRSouthLine, VIS_IREastColumn, "vis", sat_lon )
        else:
            # works also for pixels outside of the disk 
            pname = 'GEOS<'+'{:+06.1f}'.format(sat_lon)+'>'  # "GEOS<+009.5>"
            proj = {'proj': 'geos', 'a': '6378169.0', 'b': '6356583.8', 'h': '35785831.0', 'lon_0': str(sat_lon)}
            aex=(-5570248.4773392612, -5567248.074173444, 5567248.074173444, 5570248.4773392612)

            # define full disk projection 
            from pyresample.geometry import AreaDefinition
            full_disk_def = AreaDefinition('full_disk',
                                           'full_disk',
                                           pname,
                                           proj,
                                           3712,
                                           3712,
                                           aex )

            # define name and calculate area for sub-demain 
            area_name= 'MSG_'+'{:04d}'.format(VIS_IRNorthLine)+'_'+'{:04d}'.format(VIS_IRWestColumn)+'_'+'{:04d}'.format(VIS_IRSouthLine)+'_'+'{:04d}'.format(VIS_IREastColumn)
            aex = full_disk_def.get_area_extent_for_subset(3712-VIS_IRSouthLine,3712-VIS_IRWestColumn,3712-VIS_IRNorthLine,3712-VIS_IREastColumn)

            area_def = AreaDefinition(area_name,
                                      area_name,
                                      pname,
                                      proj,
                                      (VIS_IRWestColumn-VIS_IREastColumn)+1,
                                      (VIS_IRNorthLine-VIS_IRSouthLine)+1,
                                      aex )

        #print area_def
        #print "REGION:", area_def.area_id, "{"
        #print "\tNAME:\t", area_def.name
        #print "\tPCS_ID:\t", area_def.proj_id
        #print ("\tPCS_DEF:\tproj="+area_def.proj_dict['proj']+", lon_0=" + area_def.proj_dict['lon_0'] + ", a="+area_def.proj_dict['a']+", b="+area_def.proj_dict['b']+", h="+area_def.proj_dict['h'])
        #print "\tXSIZE:\t", area_def.x_size
        #print "\tYSIZE:\t", area_def.y_size
        #print "\tAREA_EXTENT:\t", area_def.area_extent
        #print "};"

        # copy area to satscene 
        satscene.area = area_def

        # write information used by mipp.xrit.MSG._Calibrator in a fake header file
        hdr = dict()

        # satellite ID number 
        hdr["SatelliteDefinition"] = dict()
        hdr["SatelliteDefinition"]["SatelliteId"] = SatelliteIds[str(satscene.sat_nr())]
        
        # processing 
        hdr["Level 1_5 ImageProduction"] = dict()
        hdr["Level 1_5 ImageProduction"]["PlannedChanProcessing"] = np_array([2,2,2,2,2,2,2,2,2,2,2,2], int)
        
        # calibration factors  
        Level15ImageCalibration = hf.get(data_folder+'METADATA/HEADER/RadiometricProcessing/Level15ImageCalibration_ARRAY')
        hdr["Level1_5ImageCalibration"] = dict()

        for chn_name in channel_numbers.keys():
            chn_nb = channel_numbers[chn_name]-1
            hdr["Level1_5ImageCalibration"][chn_nb] = dict()
            #print chn_name, chn_nb, Level15ImageCalibration[chn_nb]['Cal_Slope'], Level15ImageCalibration[chn_nb]['Cal_Offset']
            hdr["Level1_5ImageCalibration"][chn_nb]['Cal_Slope']  = Level15ImageCalibration[chn_nb]['Cal_Slope']
            hdr["Level1_5ImageCalibration"][chn_nb]['Cal_Offset'] = Level15ImageCalibration[chn_nb]['Cal_Offset']

        # loop over channels to load 
        for chn_name in satscene.channels_to_load:

            dataset_name = data_folder+'DATA/'+dict_channel[chn_name]+'/IMAGE_DATA'
            if dataset_name in hf:
                data_tmp = hf.get(data_folder+'DATA/'+dict_channel[chn_name]+'/IMAGE_DATA')

                LOG.info('hdr["SatelliteDefinition"]["SatelliteId"]: '+str(hdr["SatelliteDefinition"]["SatelliteId"]))
                #LOG.info('hdr["Level 1_5 ImageProduction"]["PlannedChanProcessing"]', hdr["Level 1_5 ImageProduction"]["PlannedChanProcessing"])
                chn_nb = channel_numbers[chn_name]-1
                LOG.info('hdr["Level1_5ImageCalibration"][chn_nb]["Cal_Slope"]:  '+str(hdr["Level1_5ImageCalibration"][chn_nb]["Cal_Slope"]))
                LOG.info('hdr["Level1_5ImageCalibration"][chn_nb]["Cal_Offset"]: '+str(hdr["Level1_5ImageCalibration"][chn_nb]["Cal_Offset"]))

                if calibrate:
                    #Calibrator = _Calibrator(hdr, chn_name)
                    bits_per_pixel = 10   ### !!! I have no idea if this is correct !!!
                    Calibrator = _Calibrator(hdr, chn_name, bits_per_pixel) ## changed call in mipp/xrit/MSG.py
                    data, calibration_unit = Calibrator (data_tmp, calibrate=1)
                else:
                    data = data_tmp
                    calibration_unit = "counts"

                LOG.info(chn_name+ " min/max: "+str(data.min())+","+str(data.max())+" "+calibration_unit )

                satscene[chn_name] = ma.asarray(data)

                satscene[chn_name].info['units'] = calibration_unit
                satscene[chn_name].info['satname'] = satscene.satname
                satscene[chn_name].info['satnumber'] = satscene.number
                satscene[chn_name].info['instrument_name'] = satscene.instrument_name
                satscene[chn_name].info['time'] = satscene.time_slot
                satscene[chn_name].info['is_calibrated'] = True

            else: 
                print "*** Warning, no data for channel "+ chn_name+ " in file "+ filename
                data = np_nan
                calibration_unit = ""
                LOG.info("*** Warning, no data for channel "+ chn_name+" in file "+filename)
                # do not append the channel chn_name

