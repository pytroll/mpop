import matplotlib as mpl   # this HAS TO BE the very first lines (before any other matplotlib functions are imported) 
mpl.use('Agg')             # this HAS TO BE the very first lines (before any other matplotlib functions are imported) 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import rcParams
from PIL import Image as PIL_Image
from TRTimage import fig2data, fig2img
from pylab import text as pylab_text
from numpy import sin, cos, radians, where, nonzero, transpose, arange, append, meshgrid, mgrid, empty, isnan, nan, percentile
from numpy import sum as np_sum
from scipy.interpolate import griddata  # interp2d
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

def prepare_figure(obj_area):

    # create new figure 
    #fig = Figure()  # old version, does not work for the stream plot 

    ## Turn interactive plotting off
    #plt.ioff()
    fig = plt.figure()   # needs a DISPLAY environment variable (simulated here with mpl.use('Agg'))

    # define size of image 
    nx = obj_area.x_size
    ny = obj_area.y_size
    # canvas figure 
    canvas = FigureCanvas(fig)
    # get dots per inch of the screen
    DPI = fig.get_dpi()
    # print "DPI", DPI
    fig.set_size_inches(nx/float(DPI),ny/float(DPI))
    # set fonts to bold
    plt.rc('font', weight='bold')
    # get axis object 
    ax = fig.add_subplot(111, aspect='equal')
    ## eliminates margins totally 
    fig.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0, wspace=0, hspace=0)
    # set limits of the axis
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # set transparent backgroud 
    fig.patch.set_alpha(0.0)        # transparent outside of diagram  
    ax.set_axis_bgcolor([1,0,0,0])  # transparent color inside diagram 

    return fig, ax 


def HRWimage( HRW_data, obj_area, hrw_channels=None, min_correlation=None, cloud_type=None, style='barbs', \
                        barb_length=None, color_mode='channel', pivot='middle', legend=True, legend_loc=3):

    """Create a PIL image from high resolution wind data
       required input:
         HRW data [HRW_class]:      HRW_class instant containing the data, see mpop/satin/nwcsaf_hrw_hdf.py
         obj_area [area_class]:     instant of area class, returned by area_def
       optional input
       hrw_channels [string array]: giving the channels that are used derive the HRW vectors
                                    e.g. hrw_channels['HRV','WV_073']
       min_correlation [int]:       minimum correlation of tracking, if below arrow is not shown
       cloud_type [int array]:      cloud types of the wind vectors, e.g. cloud_type=[8,10,11]
       style [string]:              different styles of plotting
                                    style='barbs' or style='5min_displacement' or style='15min_displacement'
       color_mode [string]:         choose color of the wind symbols, possible choises:
                                    color_mode='channel'      -> one color per SEVIRI channel used to derive HRW
                                    color_mode='pressure'     -> colorcoded cloud top pressure
                                    color_mode='temperature'  -> colorcoded cloud top temperature
                                    color_mode='cloud_type'   -> NWC-SAF cloud types
                                    color_mode='correlation'  80 ... 100
                                    color_mode='conf_nwp'      70 ... 100
                                    color_mode='conf_no_nwp'   70 ... 100
       pivot [string]:              position of the barb, e.g. pivot='middle' == center of barb at origin
       legend [True or False] :     show legend or not 
       legend_loc [string or int]:  location of the legend 
                                    upper right    1
                                    upper left     2
                                    lower left     3
                                    lower right    4
                                    right          5
                                    center left    6
                                    center right   7
                                    lower center   8
                                    upper center   9
                                    center         10
                                    best
    """

    #if min_correlation != None:
    #    print "    filter for min_correlation = ", min_correlation
    #    inds = where(HRW_data.correlation > min_correlation)
    #    HRW_data = HRW_data.subset(inds)
     
    print "... create HRWimage, color_mode = ", color_mode

    # get a empty figure with transparent background, no axis and no margins outside the diagram
    fig, ax = prepare_figure(obj_area)

    # define arrow properties 
    head_width  = 0.006 * min(obj_area.x_size,obj_area.x_size)
    head_length = 2 * head_width
    m_per_s_to_knots = 1.944

    #barb_length = 0.008 * min(obj_area.x_size,obj_area.x_size)
    
    if barb_length == None:
        n_winds = len(HRW_data.wind_id)
        if n_winds < 300:
            barb_length = 5.68
        elif n_winds < 500:
            barb_length = 5.43
        elif n_winds < 700:
            barb_length = 5.18
        elif n_winds < 900:
            barb_length = 4.68
        else:          
            barb_length = 4.00
    print "barb_length", barb_length

    if color_mode == 'channel':
        classes = ('HRV',          'VIS008 ', 'WV_062 ',   'WV_073 ',   'IR_120 ')
        colors   = ['mediumorchid', 'red',     'limegreen', 'darkgreen', 'darkturquoise']
    elif color_mode == 'pressure':
        classes = ['<200hPa',  '200-300hPa','300-400hPa','400-500hPa','500-600hPa','600-700hPa','700-800hPa', '800-900hPa','>900hPa']
        colors   = ['darksalmon', 'red'     ,'darkorange','yellow'    ,'lime'      ,'seagreen',  'deepskyblue','blue',      'mediumorchid']
        classes = tuple(['CTP '+cl for cl in classes])
    elif color_mode == 'cloud_type' or color_mode == 'cloudtype':
        classes=['non-processed','cloud free land', 'cloud free sea', 'land snow', 'sea ice',\
                 'very low cum.', 'very low', 'low cum.', 'low', 'med cum.', 'med', 'high cum.', 'high', 'very high cum.', 'very high', \
                 'sem. thin', 'sem. med.', 'sem. thick', 'sem. above', 'broken', 'undefined']

        colors = empty( (len(classes),3), dtype=int )
        colors[ 0,:] = [100, 100, 100]
        colors[ 1,:] = [  0, 120,   0]
        colors[ 2,:] = [  0,   0,   0]
        colors[ 3,:] = [250, 190, 250]
        colors[ 4,:] = [220, 160, 220]
        colors[ 5,:] = [255, 150,   0]
        colors[ 6,:] = [255, 100,   0]
        colors[ 7,:] = [255, 220,   0]
        colors[ 8,:] = [255, 180,   0]
        colors[ 9,:] = [255, 255, 140]
        colors[10,:] = [240, 240,   0]
        colors[11,:] = [250, 240, 200]
        colors[12,:] = [215, 215, 150]
        colors[13,:] = [255, 255, 255]
        colors[14,:] = [230, 230, 230]
        colors[15,:] = [  0,  80, 215]
        colors[16,:] = [  0, 180, 230]
        colors[17,:] = [  0, 240, 240]
        colors[18,:] = [ 90, 200, 160]
        colors[19,:] = [200,   0, 200]
        colors[20,:] = [ 95,  60,  30]
        colors = colors/255.
    elif color_mode in ['correlation','conf_nwp','conf_no_nwp']:
        classes  = ['<70',    '<75',     '<80',    '<85',   '<90',  '<95',   '>95' ]
        colors   = ['indigo', 'darkred', 'red','darkorange','gold', 'lime', 'green']
        classes = tuple([color_mode+' '+cl for cl in classes])
    else:
          print "*** Error in HRW_streamplot (mpop/imageo/HRWimage.py)"
          print "    unknown color_mode"
          quit()

    for wid in range(len(HRW_data.wind_id)):

        if color_mode == 'channel':

            if HRW_data.channel[wid].find('HRV') != -1: # HRV
                barbcolor = colors[0]
            elif HRW_data.channel[wid].find('VIS008') != -1: #  0.8 micro m
                barbcolor = colors[1]
            elif HRW_data.channel[wid].find('WV_062') != -1: #  6.2 micro m
                barbcolor = colors[2]
            elif HRW_data.channel[wid].find('WV_073') != -1: #  7.3 micro m
                barbcolor = colors[3]
            elif HRW_data.channel[wid].find('IR_120') != -1: # 12.0 micro m
                barbcolor = colors[4]

        elif color_mode == 'pressure':

            if HRW_data.pressure[wid] < 20000:
                barbcolor = colors[0]
            elif HRW_data.pressure[wid] < 30000:
                barbcolor = colors[1]
            elif HRW_data.pressure[wid] < 40000:
                barbcolor = colors[2]
            elif HRW_data.pressure[wid] < 50000:
                barbcolor = colors[3] 
            elif HRW_data.pressure[wid] < 60000:
                barbcolor = colors[4]
            elif HRW_data.pressure[wid] < 70000:
                barbcolor = colors[5]
            elif HRW_data.pressure[wid] < 80000:
                barbcolor = colors[6]
            elif HRW_data.pressure[wid] < 90000:
                barbcolor = colors[7]
            else:
                barbcolor = colors[8]

        elif color_mode == 'cloud_type' or color_mode == 'cloudtype':

            barbcolor = list(colors[HRW_data.cloud_type[wid], :])

        elif color_mode in ['correlation','conf_nwp','conf_no_nwp']:
            if color_mode == 'correlation':
                cdata = HRW_data.correlation
            elif color_mode == 'conf_nwp':
                cdata = HRW_data.conf_nwp
            elif color_mode == 'conf_no_nwp':
                cdata = HRW_data.conf_no_nwp

            if cdata[wid] < 70:
                barbcolor = colors[0]
            elif cdata[wid] < 75:
                barbcolor = colors[1]
            elif cdata[wid] < 80:
                barbcolor = colors[2]
            elif cdata[wid] < 85:
                barbcolor = colors[3] 
            elif cdata[wid] < 90:
                barbcolor = colors[4]
            elif cdata[wid] < 95:
                barbcolor = colors[5]
            else:
                barbcolor = colors[6]
        else:
              print "*** Error in HRW_streamplot (mpop/imageo/HRWimage.py)"
              print "    unknown color_mode"
              quit()
      
        x0, y0 = obj_area.get_xy_from_lonlat( HRW_data.lon[wid], HRW_data.lat[wid], outside_error=False) #, return_int=True

        u = HRW_data.wind_speed[wid] * -1 * sin(radians(HRW_data.wind_direction[wid])) 
        v = HRW_data.wind_speed[wid] * -1 * cos(radians(HRW_data.wind_direction[wid]))

        #print '%6s %3d %10.7f %10.7f %7.2f %7.1f %8.1f %10s' % (HRW_data.channel[wid], HRW_data.wind_id[wid], \
        #                                                        HRW_data.lon[wid], HRW_data.lat[wid], \
        #                                                        HRW_data.wind_speed[wid]*m_per_s_to_knots, \
        #                                                        HRW_data.wind_direction[wid], HRW_data.pressure[wid], barbcolor)


        if style == 'barbs':
            u = HRW_data.wind_speed[wid] * -1 * sin(radians(HRW_data.wind_direction[wid])) * m_per_s_to_knots
            v = HRW_data.wind_speed[wid] * -1 * cos(radians(HRW_data.wind_direction[wid])) * m_per_s_to_knots
            ax.barbs(x0, obj_area.y_size - y0, u * m_per_s_to_knots, v * m_per_s_to_knots, length = barb_length, pivot='middle', barbcolor=barbcolor)

        elif style == '5min_displacement' or style == '15min_displacement':
            if style == '5min_displacement':
                t_in_s =  5*60
            elif style == '15min_displacement':
                t_in_s = 15*60
            dx = u * t_in_s / obj_area.pixel_size_x
            dy = v * t_in_s / obj_area.pixel_size_y
            ax.arrow(x0, y0, dx, dy, head_width = head_width, head_length = head_length, fc=barbcolor, ec=barbcolor)

    if legend:

        rcParams['legend.handlelength'] = 0
        rcParams['legend.numpoints'] = 1

        # create blank rectangle
        rec = Rectangle((0, 0), 0, 0, fc="w", fill=False, edgecolor='none', linewidth=0)

        ##  *fontsize*: [size in points | 'xx-small' | 'x-small' | 'small' |
        ##              'medium' | 'large' | 'x-large' | 'xx-large']


        alpha=1.0
        bbox={'facecolor':'white', 'alpha':alpha, 'pad':10}

        print "... add legend: color is a function of ",  color_mode
        
        recs = empty( len(classes), dtype=object)
        recs[:] = rec 

        #if color_mode == 'pressure':
        #    recs = [rec, rec, rec, rec, rec, rec, rec, rec, rec]
        #if color_mode == 'channel':
        #    recs = [rec, rec, rec, rec, rec]
        #if color_mode in ['correlation','conf_nwp','conf_no_nwp']:
        #    recs = [rec, rec, rec, rec, rec, rec, rec]

        size=12
        if color_mode=='cloud_type':
            size=10

        leg = ax.legend(recs, classes, loc=legend_loc, prop={'size':size})

        for color,text in zip(colors,leg.get_texts()):
            text.set_color(color)


    return fig2img ( fig )


# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

def HRWstreamplot( u2d, v2d, obj_area, interpol_method, color_mode='speed', density=(3,3), linewidth_mode="scaled", linewidth_max=2.5, \
                   min_correlation=None, legend=True, legend_loc=3, vmax=None, colorbar=True, fontcolor='w'):

    """ Create a streamplot image in PIL format of a 2d wind field
           color_mode [string]:     choose color of the stream lines,
                                    color_mode='speed' -> wind speed
                                    color_mode='u'     -> u-wind component
                                    color_mode='v'     -> v-wind component
           density [2 int tuple]    density of stream lines, default density = (4,4)
           linewidth_mode [string]  "scaled" to color_mode data
                                    "const" always linewith_max
    """


    ## get a empty figure with transparent background, no axis and no margins outside the diagram
    fig, ax = prepare_figure(obj_area)
    #print dir(ax)

    # check if there is there have been enough observation (or if data is not a number array)
    if isnan(np_sum(u2d)):
        print "... there are not enough observations"
        ax.text(0.95, 0.01, 'currently not enough observations',
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, color='red', fontsize=15)
    else:
        print "there is enough data, interpolation method: ", interpol_method

        # create grid for the wind data
        [nx, ny] = u2d.shape  # for ccs4 this is (640, 710)
        Y, X = mgrid[nx-1:-1:-1, 0:ny] # watch out for Y->nx and X->ny
        #print "X.shape ", Y.shape
        #print Y[:,0]

        print "   calculate color data ", color_mode
        if color_mode == 'speed':
            from numpy import sqrt
            cdata = sqrt(u2d*u2d + v2d*v2d)
        elif color_mode == 'u':
            cdata = u2d
        elif color_mode == 'v':
            cdata = v2d
        else:
            print "*** Error in HRW_streamplot (mpop/imageo/HRWimage.py)"
            print "    unknown color_mode"
            quit()

        print "   calculate linewidth ", linewidth_mode
        if linewidth_mode == "const":
            linewidth = linewidth_max
        elif linewidth_mode == "scaled":
            if vmax != None:
                linewidth = 1 + linewidth_max*(cdata) / vmax
            else:
                linewidth = 1 + linewidth_max*(cdata) / cdata.max()
        else:
            print "*** Error in HRW_streamplot (mpop/imageo/HRWimage.py)"
            print "    unknown linewidth_mode"
            quit()

        print "... data_max =", cdata.max() ,", vmax=", vmax

        if vmax != None:
            norm = Normalize(vmin=0, vmax=vmax)
        else:
            norm = Normalize(vmin=0, vmax=cdata.max())

        #optional arguments of streamplot
        #           density=1, linewidth=None, color=None,
        #           cmap=None, norm=None, arrowsize=1, arrowstyle='-|>',
        #           minlength=0.1, transform=None, zorder=1, start_points=None, INTEGRATOR='RK4',density=(10,10)
        plt.streamplot(X, Y, u2d, v2d, color=cdata, linewidth=linewidth, cmap=plt.cm.rainbow, density=density, norm=norm) 

        if colorbar:
            colorbar_ax = fig.add_axes([0.9, 0.1, 0.05, 0.8])
            cbar = plt.colorbar(cax=colorbar_ax)
            plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=fontcolor)
            xlabel = cbar.ax.set_xlabel('m/s', weight='bold') #get the title property handler
            plt.setp(xlabel, color=fontcolor) 

        #plt.savefig("test_streamplot.png")

        # add information about interpolation method
        ax.text(0.95, 0.01, interpol_method,
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, color='green', fontsize=15)

    return fig2img ( fig )


# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------


def fill_with_closest_pixel(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """

    from numpy import isnan 
    from scipy.ndimage import distance_transform_edt

    if invalid is None: invalid = isnan(data)

    ind = distance_transform_edt(invalid, return_distances=False, return_indices=True)

    return data[tuple(ind)]

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
def HRWscatterplot( HRW_data, title='', hrw_channels=None, min_correlation=None, cloud_type=None, color_mode='direction'):

    ## get a empty figure with transparent background, no axis and no margins outside the diagram
    # fig = plt.figure()
    import pylab 
    fig = pylab.figure()
    ax = plt.subplot(111)
    ax.set_yscale("log", nonposx='clip')
    plt.scatter(HRW_data.wind_speed, HRW_data.pressure/100, s=5, c=HRW_data.wind_direction, alpha=0.5, edgecolor='none')
    pylab.title(title)
    pylab.ylim([1000,100])
    plt.yticks([1000,900,800,700,600,500,400,300,200,100], ['1000','900','800','700','600','500','400','300','200','100'], rotation='horizontal')

    p = percentile(HRW_data.wind_speed, 95)
    vmax = (round(p/10)+1)*10
    print "... vmax:", vmax 

    plt.plot([0,vmax], [680,680], color='g')
    plt.plot([0,vmax], [440,440], color='b')

    pylab.xlim([0,vmax])
    ax.set_xlabel('HRW [m/s]')
    ax.set_ylabel('p [hPa]')

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('wind direction')

    return fig2img ( fig )

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

def HRW_2dfield( HRW_data, obj_area, interpol_method=None, hrw_channels=None, min_correlation=None, level=''):

    print "... calculate 2d wind field (HRW_2dfield)"

    if min_correlation != None:
        print "    filter for min_correlation = ", min_correlation
        inds = where(HRW_data.correlation > min_correlation)
        HRW_data.subset(inds)

    xx, yy = obj_area.get_xy_from_lonlat( HRW_data.lon, HRW_data.lat, outside_error=False, return_int=False) #, return_int=True

    yy = obj_area.y_size - yy  

    uu = - HRW_data.wind_speed * sin(radians(HRW_data.wind_direction))
    vv = - HRW_data.wind_speed * cos(radians(HRW_data.wind_direction))

    # get rid of all vectors outside of the field 
    index = nonzero(xx)
    xx = xx[index]
    yy = yy[index]
    uu = uu[index]
    vv = vv[index]

    points = transpose(append([xx], [yy], axis=0))
    #print type(uu), uu.shape
    #print type(points), points.shape
    #print points[0], yy[0], xx[0]
    #print uu[0]

    nx = obj_area.x_size
    ny = obj_area.y_size

    x2 = arange(nx)
    y2 = (ny-1) - arange(ny)

    grid_x, grid_y = meshgrid(x2, y2)
    
    if interpol_method == None:
        # we need at least 2 winds to interpolate 
        if uu.size < 4:  
            print "*** Warning, not wnough wind data available, n_winds = ", uu.size
            fake = empty(grid_x.shape)
            fake[:,:] = nan
            HRW_data.interpol_method = None
            return fake, fake

        elif uu.size < 50:
            interpol_method = "RBF"

        else:
            interpol_method = "linear + nearest"
     
    #interpol_method = "nearest"
    #interpol_method = "cubic + nearest" # might cause unrealistic overshoots
    #interpol_method = "kriging"
    #interpol_method = "..."

    print "... min windspeed (org data): ", HRW_data.wind_speed.min()
    print "... max windspeed (org data): ", HRW_data.wind_speed.max()

    for i_iteration in [0,1]:

        if interpol_method == "nearest":

            print '... fill with nearest neighbour'
            # griddata, see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html
            grid_u1x = griddata(points, uu, (grid_x, grid_y), method='nearest')
            grid_v1x = griddata(points, vv, (grid_x, grid_y), method='nearest')

        elif interpol_method == "RBF":

            print '... inter- and extrapolation using radial basis functions'
            # https://www.youtube.com/watch?v=_cJLVhdj0j4
            print "... start Rbf"
            from scipy.interpolate import Rbf
            # rbfu = Rbf(xx, yy, uu, epsilon=0.1) #
            rbfu = Rbf(xx, yy, uu, epsilon=0.2)
            grid_u1x = rbfu(grid_x, grid_y)
            rbfv = Rbf(xx, yy, vv, epsilon=0.1) #
            grid_v1x = rbfv(grid_x, grid_y)
            print "... finish Rbf"
            # !very! slow for a large number of observations 

        elif interpol_method == "linear + nearest" or interpol_method == "cubic + nearest":

            if interpol_method == "linear + nearest":
                print '... calculate linear interpolation'
                # griddata, see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html
                grid_u1 = griddata(points, uu, (grid_x, grid_y), method='linear')
                grid_v1 = griddata(points, vv, (grid_x, grid_y), method='linear')
            elif interpol_method == "cubic + nearest":
                # smoother, but can cause unrealistic overshoots
                print '... calculate cubic interpolation'
                grid_u1 = griddata(points, uu, (grid_x, grid_y), method='cubic')
                grid_v1 = griddata(points, vv, (grid_x, grid_y), method='cubic')
            else:
                print "*** Error in mpop/imageo/HRWimage.py"
                print "    unknown interpolation method: ", interpol_method
                quit()

            if 1==1:
                # use faster function to extrapolate with closest neighbour
                print "... fill outside area with closest value"
                grid_u1x = fill_with_closest_pixel(grid_u1, invalid=None) 
                grid_v1x = fill_with_closest_pixel(grid_v1, invalid=None) 
            else:
                # use griddata to extrapolate with closest neighbour
                points2 = transpose(append([grid_x.flatten()], [grid_y.flatten()], axis=0))
                print type(grid_x.flatten()), grid_x.flatten().shape
                print type(points2), points2.shape
                mask = ~isnan(grid_v1.flatten())
                inds = where(mask)[0]
                grid_u1x = griddata(points2[inds], grid_u1.flatten()[inds], (grid_x, grid_y), method='nearest')
                grid_v1x = griddata(points2[inds], grid_v1.flatten()[inds], (grid_x, grid_y), method='nearest')

            if 1==0:
                # add othermost points as additional data
                y_add = [0,    0, ny-1, ny-1]
                x_add = [0, nx-1,    0, nx-1]
                for (i,j) in zip(x_add,y_add):
                    uu = append(uu, grid_u0[i,j])
                    vv = append(vv, grid_v0[i,j])
                xx = append(xx, x_add)
                yy = append(yy, y_add)
                points = transpose(append([yy], [xx], axis=0))

                print 'calc extent1'
                grid_u1e = griddata(points, uu, (grid_x, grid_y), method='linear')
                grid_v1e = griddata(points, vv, (grid_x, grid_y), method='linear')

        else:
            print "*** Error in mpop/imageo/HRWimage.py"
            print "    unknown interpol_method", interpol_method
            quit()

        ##http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
        ##http://stackoverflow.com/questions/3526514/problem-with-2d-interpolation-in-scipy-non-rectangular-grid
        #print "SmoothBivariateSpline:"
        #from scipy.interpolate import SmoothBivariateSpline
        #fitu = SmoothBivariateSpline( xx, yy, uu, s=1000) # , kx=3, ky=3, s = smooth * z2sum m 
        #from numpy import empty 
        #grid_u_SBP = empty(grid_x.shape)
        #for i in range(0,nx-1):       # starting upper right going down
        #    for j in range(0,ny-1):   # starting lower right going right
        #        #print i,j
        #        grid_u_SBP[j,i] = fitu(j,i)

        #grid_u_SBP = np.array([k.predict([x,y]) for x,y in zip(np.ravel(grid_x), np.ravel(grid_y))])
        #grid_u_SBP = grid_u_SBP.reshape(grid_x.shape)

        ##print x2
        ##print y2
        #grid_u_SBP = fitu(x2,y2)
        ##print "grid_u_SBP.shape", grid_u_SBP.shape
        ###print grid_u_SBP
        #print "End SmoothBivariateSpline:"

        #print "bisplrep:"
        #from scipy import interpolate
        #tck = interpolate.bisplrep(xx, yy, uu)
        #grid_u_BSR = interpolate.bisplev(grid_x[:,0], grid_y[0,:], tck)
        #print grid_u_BSR.shape
        #print "bisplrep"
        #print "grid_v1x.shape", grid_v1x.shape

        extent=(0,nx,0,ny)
        origin='lower'
        origin='upper'
        origin=None 

        # show different stages of 2d inter- and extra-polation 
        if 1==0:
            print 'make matplotlib.pyplot'
            import matplotlib.pyplot as plt
            vmin=-10
            vmax=10
            fig = plt.figure()
            plt.subplot(221)
            plt.title('u '+interpol_method)
            plt.plot(points[:,0], ny-1-points[:,1], 'k.', ms=1)
            plt.imshow(grid_u1x, vmin=vmin, vmax=vmax) #, extent=extent
            #plt.colorbar()
            plt.subplot(222)
            plt.title('v '+interpol_method)
            plt.plot(points[:,0], ny-1-points[:,1], 'k.', ms=1)
            plt.imshow(grid_v1x, origin=origin, vmin=vmin, vmax=vmax) #, extent=extent
            #plt.colorbar()

            # standard calculation for comparison 
            print '... calculate linear interpolation'
            grid_u1 = griddata(points, uu, (grid_x, grid_y), method='linear')
            grid_v1 = griddata(points, vv, (grid_x, grid_y), method='linear')
            grid_u1xx = fill_with_closest_pixel(grid_u1, invalid=None) 
            grid_v1xx = fill_with_closest_pixel(grid_v1, invalid=None) 

            plt.subplot(223)
            plt.title('U Linear+Nearest')
            plt.plot(points[:,0], ny-1-points[:,1], 'k.', ms=1)
            plt.imshow(grid_u1xx, origin=origin, vmin=vmin, vmax=vmax) #, extent=extent
            #plt.colorbar()
            plt.subplot(224)
            plt.title('V Linear+Nearest') 
            plt.plot(points[:,0], ny-1-points[:,1], 'k.', ms=1)
            #plt.title('Cubic')
            plt.imshow(grid_v1xx, origin=origin, vmin=vmin, vmax=vmax) #, extent=extent
            #plt.colorbar()
            plt.gcf().set_size_inches(6, 6)
            #plt.show()  # does not work with AGG
            tmpfile="test_hrw"+level+".png"
            fig.savefig(tmpfile)
            print "display "+tmpfile+" &"


        if grid_u1x.min() < -150 or grid_v1x.min() < -150 or grid_u1x.max() > 150 or grid_v1x.max() > 150:
            print "*** Warning, numerical instability detected, interpolation method: ", interpol_method
            print "    min u windspeed (u 2dimensional): ", grid_u1x.min()
            print "    min v windspeed (v 2dimensional): ", grid_v1x.min()
            print "    max u windspeed (u 2dimensional): ", grid_u1x.max()
            print "    max v windspeed (v 2dimensional): ", grid_v1x.max()
            interpol_method = "glinear + nearest"
            print "... try another interpolation method: ", interpol_method
        else:
            # (hopefully) numerical stable interpolation, exit the interpolation loop
            break 

    HRW_data.interpol_method = interpol_method

    return grid_u1x, grid_v1x
