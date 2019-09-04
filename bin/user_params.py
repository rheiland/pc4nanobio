 
# This file is auto-generated from a Python script that parses a PhysiCell configuration (.xml) file.
#
# Edit at your own risk.
#
import os
from ipywidgets import Label,Text,Checkbox,Button,HBox,VBox,FloatText,IntText,BoundedIntText,BoundedFloatText,Layout,Box
    
class UserTab(object):

    def __init__(self):
        
        micron_units = Label('micron')   # use "option m" (Mac, for micro symbol)

        constWidth = '180px'
        tab_height = '500px'
        stepsize = 10

        #style = {'description_width': '250px'}
        style = {'description_width': '25%'}
        layout = {'width': '400px'}

        name_button_layout={'width':'25%'}
        widget_layout = {'width': '15%'}
        units_button_layout ={'width':'15%'}
        desc_button_layout={'width':'45%'}

        param_name1 = Button(description='tumor_radius', disabled=True, layout=name_button_layout)
        param_name1.style.button_color = 'lightgreen'

        self.tumor_radius = FloatText(
          value=750,
          step=10,
          style=style, layout=widget_layout)

        param_name2 = Button(description='enable_pH_effects', disabled=True, layout=name_button_layout)
        param_name2.style.button_color = 'tan'

        self.enable_pH_effects = Checkbox(
          value=False, disabled=True,
          style=style, layout=widget_layout)

        units_button1 = Button(description='micron', disabled=True, layout=units_button_layout) 
        units_button1.style.button_color = 'lightgreen'
        units_button2 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button2.style.button_color = 'tan'

        desc_button1 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button1.style.button_color = 'lightgreen'
        desc_button2 = Button(description='(disabled in this version)', disabled=True, layout=desc_button_layout) 
        desc_button2.style.button_color = 'tan'

        row1 = [param_name1, self.tumor_radius, units_button1, desc_button1] 
        row2 = [param_name2, self.enable_pH_effects, units_button2, desc_button2] 

        box_layout = Layout(display='flex', flex_flow='row', align_items='stretch', width='100%')
        box1 = Box(children=row1, layout=box_layout)
        box2 = Box(children=row2, layout=box_layout)

        self.tab = VBox([
          box1,
          box2,
        ])

    # Populate the GUI widgets with values from the XML
    def fill_gui(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        self.tumor_radius.value = float(uep.find('.//tumor_radius').text)
        self.enable_pH_effects.value = ('true' == (uep.find('.//enable_pH_effects').text.lower()) )


    # Read values from the GUI widgets to enable editing XML
    def fill_xml(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        uep.find('.//tumor_radius').text = str(self.tumor_radius.value)
        uep.find('.//enable_pH_effects').text = str(self.enable_pH_effects.value)
