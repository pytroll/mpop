from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
#import matplotlib as mpl
#mpl.use('Agg')
#from pylab import figure
from pylab import rand
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from uuid import uuid4
import subprocess
from PIL import Image as PIL_Image
import numpy

def fig2data ( fig ):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )

    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring ( fig.canvas.tostring_argb(), dtype=numpy.uint8 )
    buf.shape = ( w, h, 4 )

    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = numpy.roll ( buf, 3, axis = 2 )
    return buf

def fig2img ( fig ):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data ( fig )
    w, h, d = buf.shape
    return PIL_Image.frombytes( "RGBA", ( w ,h ), buf.tostring( ) )


def TRTimage( TRTcell_IDs, TRTcells, obj_area, minRank=8, alpha_max=1.0, plot_vel=True):

    # define size of image 
    nx = obj_area.x_size
    ny = obj_area.y_size

    # create new figure 
    fig = Figure()
    # canvas figure 
    canvas = FigureCanvas(fig)
    # get dots per inch of the screen
    DPI = fig.get_dpi()
    # print "DPI", DPI
    fig.set_size_inches(nx/float(DPI),ny/float(DPI))
    # get axis object 
    ax = fig.add_subplot(111, aspect='equal')
    ## eliminates margins totally 
    fig.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0, wspace=0, hspace=0)
    #plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0) # does only work with x11 display 
    # set limits of the axis
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # set transparent backgroud 
    fig.patch.set_alpha(0.0)        # transparent outside of diagram  
    ax.set_axis_bgcolor([1,0,0,0])  # transparent color inside diagram 

    # define arrow properties 
    head_width  = 0.006 * min(obj_area.x_size,obj_area.x_size)
    head_length = 2 * head_width

    pixel_size_x_km = 0.001 * obj_area.pixel_size_x
    pixel_size_y_km = 0.001 * obj_area.pixel_size_y

    for cell in TRTcell_IDs:

        if TRTcells[cell].RANKr > minRank:

            (x0,y0) = obj_area.get_xy_from_lonlat(TRTcells[cell].lon, TRTcells[cell].lat, outside_error=False, return_int=False)
            y0 = (obj_area.y_size-1)-y0
            # print (x0,y0)

            vx = TRTcells[cell].vel_x
            vy = TRTcells[cell].vel_y
   
            # !!!scaling of width and height is not correct, that is on map projection, but not on the ground!!!
            e = Ellipse( xy     =  (x0, y0),                           \
                         width  =  2*TRTcells[cell].ell_S / pixel_size_x_km, \
                         height =  2*TRTcells[cell].ell_L / pixel_size_y_km, \
                         angle  = -TRTcells[cell].angle )
            
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            
            if TRTcells[cell].RANKr <= 12:
                cell_color="white"
                alpha = (alpha_max-0.2) / 12. * TRTcells[cell].RANKr
            elif TRTcells[cell].RANKr <= 15:
                cell_color="white"
                alpha = alpha_max
            elif TRTcells[cell].RANKr <= 25:
                cell_color="green"
                alpha = alpha_max
            elif TRTcells[cell].RANKr <= 35:
                cell_color="yellow"
                alpha = alpha_max
            else:
                cell_color="red"
                alpha = alpha_max
            # print "cell ID: %s, cell rank: %2d, cell_color:%7s, alpha = %4.1f" % (cell, TRTcells[cell].RANKr, cell_color, alpha)
            e.set_alpha(alpha)       # transparency: 0.0 transparent, 1 total visible  
            e.set_facecolor(cell_color)  # "white" or [1,1,1]

            if plot_vel:
                ax.arrow(x0, y0, vx, vy, head_width = head_width, head_length = head_length, fc=cell_color, ec=cell_color)

    if 1==1:
        # print " !!! convert fig to image by function fig2img !!!"
        ### this would avoid saving into a file, but it fills the transparent areas with "white"
        PIL_image = fig2img ( fig )  
    else: 
        tmp_file = '/tmp/TRT_'+str(uuid4())+'.png'
        # print tmp_file
        plt.savefig(tmp_file, dpi=DPI, transparent=True) #, bbox_inches='tight'
        # subprocess.call("display "+tmp_file+" &", shell=True) 
        PIL_image = PIL_Image.open(tmp_file)
        subprocess.call("rm "+tmp_file+" &", shell=True) 

    return PIL_image
