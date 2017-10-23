# -*- coding: utf-8 -*-
from .saveable import Saveable

class Plotable:
    """Plotable objects can be plotted
    
    """
 
    def __init__(self):
        self.init_plotconf()
    
    def init_plotconf(self):
        """Sets default settings for the plotting.
        
        This method needs to be implemented by subclasses. Its only purpose is to set the plotconf attribute of the Plotable class
        """
        self.plotconf=PlotConfiguration()
    
    def _get_axes(self):
        """Returns one or two axis for plotting
        
        """
        # here we construct the axes
        xaxis = 0
        yaxis = 0
        if self.plotconf.dim == 2:
            return xaxis
        elif self.plotconf.dim == 3:
            return (xaxis, yaxis)
        
        
        
    def _get_values(self):
        """Returns values for plotting
        
        """
        # here we construct the values
        yvalues = 0
        zvalues = 0
        if self.plotconf.dim == 2:
            return yvalues
        elif self.plotconf.dim == 3:
            return zvalues
            
        
    def plot(self, show=False, conf=None, save_conf=True):
        """Plots the Plotable object """
        if conf is not None:
            if not save_conf:
                backup = self.plotconf
            self.plotconf = conf
        if self.plotconf is None:
            self.init_plotconf()
            
        axes = self._get_axes()
        vals = self._get_values()
        
        import matplotlib.pyplot as plt
        if self.plotconf.dim == 2:
             plt.plot(axes[0], vals)
        elif self.plotconf.dim == 3:
             plt.contourf(axes[0], axes[1], vals)
             
        if show:
            # configure the plot here
            self._configure_plot()
            plt.show()
        if not save_conf:
            self.plotconf = backup
            
    
    # can be made a function
    def _configure_plot(self):
        """Configures current figure """
        ...
        self._plot_configured = True
        
            
    def save_fig(self, file, conf=None):
        """Saves the figure into a file"""
        if not self._plot_configured:
            self._configure_plot()
        ...
    
        
    def get_PlotConfiguration(self):
        """Returns PlotConfiguration object
        """
        if self.plotconf is None:
            self._init_plotable()
        return self.plotconf



class PlotConfiguration(Saveable):
    """Configuration for plotting instances of Plotable class
    """

    def __init__(self):
        self.title = ""
        self.subtitle = ""
        self.colormap = ...
        self.texts = {"text1":[0.2, 0.8, "Sanf-serif", 10]}
        self.text_font = "Times"
        self.text_size = 15
        self.xlabel = ["x-axis", ["Times", 20]]
        self.ylabel = ["y-axis", ["Times", 20]]
        
        
        
        
        
    def set_axis_attribute(self, x_attr_name=None, y_attr_name=None):
        """Sets names of the object attributes to be used as plotting axes
        """
        self.x_attr_name = x_attr_name
        self.y_attr_name = y_attr_name
        
        
    def set_values_attribute(self, val_attr_name=None):
        """Sets the name of the object attribute to be used as values for plotting
        """
        self.val_attr_name = val_attr_name
        
    def set_data_dimension(self, dim=1):
        """Sets data dimensionality
        """
        self.dim = dim
        
    
        
    def __str__(self):
        """Show all options for configuration
        
        """
        pass
        

