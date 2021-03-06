# -*- coding: utf-8 -*-
"""Class holding a set of TwoDSpectra



    Class Details
    -------------

"""
#import h5py
#import matplotlib.pyplot as plt  
import numpy

from ..core.time import TimeAxis
from ..core.valueaxis import ValueAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction
from .twod2 import TwoDSpectrum

from ..core.managers import energy_units
from .. import COMPLEX
from .. import REAL

from ..core.saveable import Saveable


class TwoDSpectrumContainer(Saveable):
    """Class holding a set of TwoDSpectra
    

    Parameters
    ----------
    
    t2axis: TimeAxis
       object holding waiting times at which spectra are calculated
       
    keep_pathways: bool
       if set True, the container will keep all types of Liouville pathways
       stored separately
       
    keep_stypes: bool
       if se t True, the container will keep rephasing and non-rephasing 
       spectra stored separately
       
       
    """
    
    def __init__(self, t2axis=None, keep_pathways=False, keep_stypes=True):
        
        self.t2axis = t2axis
        self.keep_pathways = keep_pathways
        self.keep_stypes = keep_stypes
        
        self.axis = None
        
        self.itype = None
        self.index = 0
        self.tags = []
        
        if self.keep_pathways:
            raise Exception("Container keeping pathways not available yet")
            
        self.spectra = {}
        self._which = None
        
        
    def use_indexing_type(self, itype):
        """Sets the type of indices used to identify spectra
        
        Parameters
        ----------
        
        itype : string, ValueAxis, TimeAxis, FrequencyAxis
            Type of indexig. If itype is a string, it should have values
            either 'integer' or 'string' in which case the specra will be 
            stored by integer index or by string (as in dictionary). If 
            itype is a ValueAxis (TimeAxis, FrequencyAxis), spectra will
            be indexed by the values in the axis object.
        
        """
        
        if isinstance(itype, str):
            if itype == "integer":
                self.itype = "integer"
            elif itype == "string":
                self.itype = "string"
            else:
                raise Exception("Unknown indexing type")
        elif isinstance(itype, ValueAxis):
            if isinstance(itype, TimeAxis):
                self.itype = "TimeAxis"
                self.axis = itype
                # This axis must FFT in the "standard" way 
                # -- it must of the "complete" type 
                self.axis.atype = "complete" 
            elif isinstance(itype, FrequencyAxis):
                self.itype = "FrequencyAxis"
                self.axis = itype
            else:
                self.itype = "ValueAxis"
                self.axis = itype
        else:
            raise Exception("Unknown indexing type")
        

    def set_spectrum(self, spect, tag=None):
        """Stores spectrum with a tag (time, index, etc.)
        
        Stores the spectrum according to present indexing scheme
        
        Parameters
        ----------
        
        spect : TwoDSpectrum
            Object holding the spectrum; when not tag is specified for a
            spectrum which as its t2 time set, the tag is set to t2 time.
            
        tag : {int, string, ValuesAxis, TimeAxis, FrequencyAxis}
            Tag which will be used for retrieval of the spectrum from 
            the container.
        
        """
        
        
        if self.itype == "integer":
            self.spectra[str(self.index)] = spect
            self.index += 1
            return self.index
        
        elif self.itype in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            
            if tag is None:
                # we will read the spectrum intrinsic t2 time and set it as tag
                if spect.t2 >= 0.0:
                    tag = spect.t2
                    
            if tag is not None:
                if tag in self.axis.data:
                    self.spectra[tag] = spect
                    self.tags.append(tag)
                    self.index += 1
                else:
                    raise Exception("Tag not compatible with the ValueAxis")
            else:
                raise Exception("No tag specified, and spectrum"
                                +" does not have t2 time set"
                                +" - cannot store spectrum")
            return self.index
        
        elif self.itype == "string":
            if tag is not None:
                stag = str(tag)
                self.spectra[stag] = spect
                self.tags.append(stag)
                self.index += 1
            else:
                raise Exception("No tag specified - cannot store spectrum")
            return self.index

        else:
            
            raise Exception("Unknown type of indexing")    


    def _lousy_equal(self, x1, x2, dx, frac=0.25):
        """Equals up to fraction of dx
        
        This function returns True if x1 is closer to x2 than `frac` of 
        a specified interval. In addition it saves the value of x2 to which
        x1 is equal in the attribute _which of the present class.
        
        
        """
        if abs(x1-x2) < dx*frac: 
            self._which = x2
            return True
        
        self._which = None
        return False


    def get_spectrum_by_index(self, indx):
        """Returns spectrum by integet index
        
        The integer index is assigned to all spectra in the order they were
        saved to the container. They can be retrieved in this order
        
        Parameters
        ----------
        
        indx : int
            Index of the spectrum to be retrieved
            
        """
        
        if self.itype == "integer":
            
            return self.get_spectrum(indx)
        
        else:

            return self.spectra[self.tags[indx]]


    def get_spectrum(self, tag):
        """Returns spectrum corresponing to time t2
        
        Checks if the time t2 is present in the t2axis
        
        Parameters
        ----------
        
        t2 : float
            Waiting time for which spectrum should be returned
            
            
        """        
        if self.itype in ["integer", "string"]:
            
            return self.spectra[str(tag)]

        elif self.itype in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            if any(self._lousy_equal(tag, li, self.axis.step) 
               for li in self.axis.data):

                return self.spectra[self._which]     
            else:
                raise Exception("Tag not compatible with the ValueAxis")
        
        else:
            
            raise Exception("Unknown type of indexing")


    def length(self):
        """Returns the length of the container
        
        
        """
        return len(self.spectra.keys())


    def get_spectra(self, start=None, end=None):
        """Returns a list of the calculated spectra
        
        Returns all spectra or an interval of spectra when `start` and `end`
        are specified
        
        Parameters
        ----------
        
        start : int
            Index of the first spectrum to be returned
            
        end : int
            Index of the last spectrum to be returned

        """
        
        ven = [value for (key, value) in sorted(self.spectra.items())]
        
        if (start is None) and (end is None): 
            return ven
        else:
            ven2 = []
            vkeys = [key for (key, value) in sorted(self.spectra.items())]
            for k in vkeys:
                if k >= start and k <= end:
                    ven2.append(self.spectra[k])
            return ven2

        
    def get_PumpProbeSpectrumContainer(self, skip=0):
        """Converts this container into PumpProbeSpectrumContainer
        
        """
        
        from .pumpprobe import PumpProbeSpectrumContainer
        
        k = 0
        ppc = []
        for sp in self.get_spectra():
            if k == 0:
                pp = sp.get_PumpProbeSpectrum()
                ppc.append(pp)
            k += 1
            if k > skip:
                k = 0
            
        length = len(ppc)
        start = ppc[0].get_t2()
        step = ppc[1].get_t2()-start

        naxis = TimeAxis(start,length,step)        
        ppcont = PumpProbeSpectrumContainer(t2axis=naxis)

        for sp in ppc:
            ppcont.set_spectrum(sp)
            
        return ppcont

    
    def get_point_evolution(self, x, y, times):
        """Tracks an evolution of a single point on the 2D spectrum
        
        
        Parameters
        ----------
        
        x : float
            x coordinate in the 2D spectrum (usually omega_1 axis)

        y : float
            y coordinate in the 2D spectrum (usually omega_3 axis)
            
        times : ValueAxis
            Times (usually waiting t_2 times) in which spectra are taken
            
        """
        
        vals = numpy.zeros(times.length, dtype=COMPLEX)
        k = 0
        for t2 in times.data:
            
            sp = self.get_spectrum(t2)
            vals[k] = sp.get_value_at(x,y)
            k +=1
            
        return vals

    
    def fft(self, ffttype="complex-positive", window=None, dtype="total", tag=None):
        """Fourier transform in t2 time
        
        This method performs FFT on the container data determined by the
        value of the `dtype` argument. The new container is created and 
        the storage resolution of its components is set `off`. This means
        that the container and its spectra have no idea about what data they
        store. Even when plotting the spectra, one has to set plotting of
        the `total` spectrum.
        
        Parameters
        ----------
        
        ffttype : string
            Specifies the type Fourier transform we perform
            
        window : DFunction
            Windowing function for the data. Default is None
        
        """
        
        legacy = False
        
        if self.itype not in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            raise Exception("FFT cannot be performed for"+
                            " this type of indexing")

        # even when no window function is supplied, we create one with
        # all elements equal to one
        if window is None:
            winfce = DFunction(self.axis, 
                               numpy.ones(self.axis.length, dtype=REAL))
        else:
            winfce = window
                
        # put all data into one array
        tags = self.axis.data
        Nos = self.length()

        if len(tags) == Nos:
            tag1 = self.axis.data[0]
            sp1 = self.get_spectrum(tag1)
            
            if legacy:
                N1, N2 = sp1.data.shape
                data = numpy.zeros((N1, N2, Nos), dtype=sp1.data.dtype)
            else:
                sp1.set_data_flag(dtype)
                N1, N2 = sp1.d__data.shape
                data = numpy.zeros((N1, N2, Nos), dtype=sp1.d__data.dtype)
            for k_n in range(Nos):
                tag = self.axis.data[k_n]
                
                if legacy:
                    spect = self.get_spectrum(tag)
                    data[:,:,k_n] = spect.data
                else:
                    spect = self.get_spectrum(tag)
                    spect.set_data_flag(dtype)
                    data[:,:,k_n] = spect.d__data
            
        else:
            raise Exception("Number of spectra not consistent"+
                            " with ValueAxis object")
            
        #
        # FFT of the axis
        #
        if isinstance(self.axis, TimeAxis):
            new_axis = self.axis.get_FrequencyAxis()
            
        elif isinstance(self.axis, FrequencyAxis):
            new_axis = self.axis.get_TimeAxis()
            
        else: 
            # this must be ValueAxis

            ftaxdata = (2.0*numpy.pi)*numpy.fft.fftfreq(self.axis.length,
                                                        self.axis.step)
            ftaxdata = numpy.fft.fftshift(ftaxdata)
            
            start = ftaxdata[0]
            length = len(ftaxdata)
            step = ftaxdata[1]-ftaxdata[0]

            new_axis = ValueAxis(start, length, step)            

        
        #
        # FFT of the data
        #
        
        # window function
        ftdata = numpy.zeros(data.shape, dtype=data.dtype)
        for i_n in range(data.shape[0]):
            for j_n in range(data.shape[1]):
                ftdata[i_n,j_n,:] = data[i_n,j_n,:]*winfce.data
        ftdata = numpy.fft.fft(ftdata, axis=2)
        ftdata = numpy.fft.fftshift(ftdata, axes=2)
        
        # save it to a new container
        new_container = TwoDSpectrumContainer()
        new_container.use_indexing_type(new_axis)

        for k_n in range(Nos):
            tag = new_axis.data[k_n]
            spect = TwoDSpectrum()
            spect.set_axis_1(sp1.xaxis)
            spect.set_axis_3(sp1.yaxis)
            
            if legacy:
                spect.set_data(ftdata[:, :, k_n])
            else:
                spect.set_resolution("off")
                spect._add_data(ftdata[:, :, k_n],dtype="total")

            new_container.set_spectrum(spect, tag=tag)
        
        
        return new_container

          
    def _create_root_group(self, start, name):
        return start.create_group(name)


    def _save_axis(self, rt, name, ax):
        axdir = rt.create_group(name)
        axdir.attrs.create("start",ax.start)
        axdir.attrs.create("length",ax.length)
        axdir.attrs.create("step",ax.step)


    def _load_axis(self, rt, name):
        axdir = rt[name]
        start = axdir.attrs["start"]
        length = axdir.attrs["length"]
        step = axdir.attrs["step"]
        return TimeAxis(start, length, step) 

        
#    # FIXME: this through Savable
#    def save(self, filename):
#        """Saves the whole object into file
#        
#        
#        """
#        with energy_units("int"):
#            with h5py.File(filename,"w") as f:
#                self._save_axis(f,"t2axis",self.t2axis)
#                rt = self._create_root_group(f, "spectra")            
#                for sp in self.get_spectra():
#                    t2 = sp.get_t2
#                    rgname = "spectrum_"+str(t2)
#                    srt = sp._create_root_group(rt,rgname)
#                    sp._save_attributes(srt)
#                    sp._save_data(srt)
#                    sp._save_axis(srt,"xaxis",sp.xaxis,)
#                    sp._save_axis(srt,"yaxis",sp.yaxis)
#            
#      
#    # FIXME: this through Savable    
#    def load(self, filename):
#        """Loads the whole object from a file
#        
#        
#        """
#        with energy_units("int"):
#            with h5py.File(filename,"r") as f:
#                self.t2axis = self._load_axis(f, "t2axis")
#                rt = f["spectra"]
#                for key in rt.keys():
#                    sp = TwoDSpectrum()
#                    srt = rt[key]
#                    sp._load_attributes(srt)
#                    sp._load_data(srt)
#                    sp.xaxis = sp._load_axis(srt,"xaxis")
#                    sp.yaxis = sp._load_axis(srt,"yaxis")
#                    
#                    self.set_spectrum(sp)
#

    def trimall_to(self, window=None):
        """Trims all spectra in the container

        Parameters
        ----------
        
        window: list of floats
            Window, specified by four float number, to which all spectra
            in the container should be trimmed
            
        """
        if window is not None:
            axes = window
            for s in self.get_spectra():
                s.trim_to(window=axes)

    # FIXME: this is still legacy version        
    def amax(self, spart="real"):
        """Returns maximum amplitude of the spectra in the container
        
        """
        mxs = []
        for s in self.get_spectra():       
            spect2D = numpy.real(s.reph2D) + numpy.real(s.nonr2D)   
            mx = numpy.amax(spect2D)
            mxs.append(mx)
        return numpy.amax(numpy.array(mxs))


    # Print iterations progress
    def _printProgressBar(self, iteration, total, 
                          prefix = '', suffix = '', 
                          decimals = 1, length = 100,
                          fill='*'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            
        Based on: 
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
        """
#                          fill = '█'):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()
    
             
    def make_movie(self, filename, window=None,
                   stype="total", spart="real", 
                   cmap=None, Npos_contours=10,
                   frate=20, dpi=100, start=None, end=None,
                   show_states=None, progressbar=False, vmax=None):
        
        import matplotlib.pyplot as plt
        import matplotlib.animation as manimation
        
        FFMpegWriter = manimation.writers["ffmpeg"]
        metadata = dict(title="Test Movie", artist='Matplotlib',
                comment='Movie support!')
        writer = FFMpegWriter(fps=frate, metadata=metadata)
        
        fig = plt.figure() 
        
        spctr = self.get_spectra()
        l = len(spctr)
        last_t2 = spctr[l-1].get_t2()
        first_t2 = spctr[0].get_t2()
        
        if vmax is None:
            mx = self.amax()
        else:
            mx = vmax
        
        if start is None:
            start = first_t2
        if end is None:
            end = last_t2
        
                
        with writer.saving(fig, filename, dpi):  
            k = 0
            # Initial call to print 0% progress
            sp2write = self.get_spectra(start=start, end=end)
            l = len(sp2write)
            if progressbar:
                self._printProgressBar(0, l, prefix = 'Progress:',
                                       suffix = 'Complete', length = 50)
            for sp in self.get_spectra(start=start, end=end):
                sp.plot(fig=fig, window=window, cmap=cmap, vmax=mx, 
                        Npos_contours=Npos_contours,
                        stype=stype,spart=spart,
                        show_states=show_states,
                        label="T="+str(sp.get_t2())+"fs")
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1
#                if k == 20:
#                    return
