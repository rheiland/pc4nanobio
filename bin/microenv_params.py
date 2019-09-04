 
# This file is auto-generated from a Python script that parses a PhysiCell configuration (.xml) file.
#
# Edit at your own risk.
#
import os
from ipywidgets import Label,Text,Checkbox,Button,HBox,VBox,FloatText,IntText,BoundedIntText,BoundedFloatText,Layout,Box
    
class MicroenvTab(object):

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


        menv_var1 = Button(description='oxygen (mmHg)', disabled=True, layout=name_button_layout)
        menv_var1.style.button_color = 'tan'

        param_name1 = Button(description='diffusion_coefficient', disabled=True, layout=name_button_layout)

        self.oxygen_diffusion_coefficient = FloatText(value=100000.000000,
          step=10000,style=style, layout=widget_layout)

        param_name2 = Button(description='decay_rate', disabled=True, layout=name_button_layout)

        self.oxygen_decay_rate = FloatText(value=.1,
          step=0.01,style=style, layout=widget_layout)

        menv_var2 = Button(description='glucose', disabled=True, layout=name_button_layout)
        menv_var2.style.button_color = 'lightgreen'

        param_name3 = Button(description='diffusion_coefficient', disabled=True, layout=name_button_layout)

        self.glucose_diffusion_coefficient = FloatText(value=30000,
          step=1000,style=style, layout=widget_layout)

        param_name4 = Button(description='decay_rate', disabled=True, layout=name_button_layout)

        self.glucose_decay_rate = FloatText(value=0.00003,
          step=1e-06,style=style, layout=widget_layout)

        menv_var3 = Button(description='ECM', disabled=True, layout=name_button_layout)
        menv_var3.style.button_color = 'tan'

        param_name5 = Button(description='diffusion_coefficient', disabled=True, layout=name_button_layout)

        self.ECM_diffusion_coefficient = FloatText(value=0,
          step=0.01,style=style, layout=widget_layout)

        param_name6 = Button(description='decay_rate', disabled=True, layout=name_button_layout)

        self.ECM_decay_rate = FloatText(value=0.0,
          step=0.01,style=style, layout=widget_layout)

        menv_var4 = Button(description='NP1', disabled=True, layout=name_button_layout)
        menv_var4.style.button_color = 'lightgreen'

        param_name7 = Button(description='diffusion_coefficient', disabled=True, layout=name_button_layout)

        self.NP1_diffusion_coefficient = FloatText(value=0,
          step=0.01,style=style, layout=widget_layout)

        param_name8 = Button(description='decay_rate', disabled=True, layout=name_button_layout)

        self.NP1_decay_rate = FloatText(value=0.0,
          step=0.01,style=style, layout=widget_layout)

        menv_var5 = Button(description='pH', disabled=True, layout=name_button_layout)
        menv_var5.style.button_color = 'tan'

        param_name9 = Button(description='diffusion_coefficient', disabled=True, layout=name_button_layout)

        self.pH_diffusion_coefficient = FloatText(value=0,
          step=0.01,style=style, layout=widget_layout)

        param_name10 = Button(description='decay_rate', disabled=True, layout=name_button_layout)

        self.pH_decay_rate = FloatText(value=0.0,
          step=0.01,style=style, layout=widget_layout)
        self.calculate_gradient = Checkbox(description='calculate_gradients', disabled=False, layout=desc_button_layout)
        self.track_internal = Checkbox(description='track_in_agents', disabled=False, layout=desc_button_layout)


         #  ------- micronenv info
        menv_units_button1 = Button(description='micron^2/min', disabled=True, layout=units_button_layout) 
        menv_units_button2 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        menv_units_button3 = Button(description='micron^2/min', disabled=True, layout=units_button_layout) 
        menv_units_button4 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        menv_units_button5 = Button(description='micron^2/min', disabled=True, layout=units_button_layout) 
        menv_units_button6 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        menv_units_button7 = Button(description='micron^2/min', disabled=True, layout=units_button_layout) 
        menv_units_button8 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        menv_units_button9 = Button(description='micron^2/min', disabled=True, layout=units_button_layout) 
        menv_units_button10 = Button(description='1/min', disabled=True, layout=units_button_layout) 




        row_oxygen = [menv_var1,  ] 
        row1 = [param_name1, self.oxygen_diffusion_coefficient, menv_units_button1]
        row2 = [param_name2, self.oxygen_decay_rate, menv_units_button2]
        row_glucose = [menv_var2,  ] 
        row3 = [param_name3, self.glucose_diffusion_coefficient, menv_units_button3]
        row4 = [param_name4, self.glucose_decay_rate, menv_units_button4]
        row_ECM = [menv_var3,  ] 
        row5 = [param_name5, self.ECM_diffusion_coefficient, menv_units_button5]
        row6 = [param_name6, self.ECM_decay_rate, menv_units_button6]
        row_NP1 = [menv_var4,  ] 
        row7 = [param_name7, self.NP1_diffusion_coefficient, menv_units_button7]
        row8 = [param_name8, self.NP1_decay_rate, menv_units_button8]
        row_pH = [menv_var5,  ] 
        row9 = [param_name9, self.pH_diffusion_coefficient, menv_units_button9]
        row10 = [param_name10, self.pH_decay_rate, menv_units_button10]
        row11 = [self.calculate_gradient,]
        row12 = [self.track_internal,]


        box_layout = Layout(display='flex', flex_flow='row', align_items='stretch', width='100%')
        box_oxygen = Box(children=row_oxygen, layout=box_layout)
        box1 = Box(children=row1, layout=box_layout)
        box2 = Box(children=row2, layout=box_layout)
        box_glucose = Box(children=row_glucose, layout=box_layout)
        box3 = Box(children=row3, layout=box_layout)
        box4 = Box(children=row4, layout=box_layout)
        box_ECM = Box(children=row_ECM, layout=box_layout)
        box5 = Box(children=row5, layout=box_layout)
        box6 = Box(children=row6, layout=box_layout)
        box_NP1 = Box(children=row_NP1, layout=box_layout)
        box7 = Box(children=row7, layout=box_layout)
        box8 = Box(children=row8, layout=box_layout)
        box_pH = Box(children=row_pH, layout=box_layout)
        box9 = Box(children=row9, layout=box_layout)
        box10 = Box(children=row10, layout=box_layout)
        box11 = Box(children=row11, layout=box_layout)
        box12 = Box(children=row12, layout=box_layout)

        self.tab = VBox([
          box_oxygen,
          box1,
          box2,
          box_glucose,
          box3,
          box4,
          box_ECM,
          box5,
          box6,
          box_NP1,
          box7,
          box8,
          box_pH,
          box9,
          box10,
          box11,
          box12,
        ])

    # Populate the GUI widgets with values from the XML
    def fill_gui(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point

        self.oxygen_diffusion_coefficient.value = float(vp[0].find('.//diffusion_coefficient').text)
        self.oxygen_decay_rate.value = float(vp[0].find('.//decay_rate').text)

        self.glucose_diffusion_coefficient.value = float(vp[1].find('.//diffusion_coefficient').text)
        self.glucose_decay_rate.value = float(vp[1].find('.//decay_rate').text)

        self.ECM_diffusion_coefficient.value = float(vp[2].find('.//diffusion_coefficient').text)
        self.ECM_decay_rate.value = float(vp[2].find('.//decay_rate').text)

        self.NP1_diffusion_coefficient.value = float(vp[3].find('.//diffusion_coefficient').text)
        self.NP1_decay_rate.value = float(vp[3].find('.//decay_rate').text)

        self.pH_diffusion_coefficient.value = float(vp[4].find('.//diffusion_coefficient').text)
        self.pH_decay_rate.value = float(vp[4].find('.//decay_rate').text)

        if uep.find('.//options//calculate_gradients').text.lower() == 'true':
          self.calculate_gradient.value = True
        else:
          self.calculate_gradient.value = False
        if uep.find('.//options//track_internalized_substrates_in_each_agent').text.lower() == 'true':
          self.track_internal.value = True
        else:
          self.track_internal.value = False
        


    # Read values from the GUI widgets to enable editing XML
    def fill_xml(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp[0].find('.//diffusion_coefficient').text = str(self.oxygen_diffusion_coefficient.value)
        vp[0].find('.//decay_rate').text = str(self.oxygen_decay_rate.value)
        vp[1].find('.//diffusion_coefficient').text = str(self.glucose_diffusion_coefficient.value)
        vp[1].find('.//decay_rate').text = str(self.glucose_decay_rate.value)
        vp[2].find('.//diffusion_coefficient').text = str(self.ECM_diffusion_coefficient.value)
        vp[2].find('.//decay_rate').text = str(self.ECM_decay_rate.value)
        vp[3].find('.//diffusion_coefficient').text = str(self.NP1_diffusion_coefficient.value)
        vp[3].find('.//decay_rate').text = str(self.NP1_decay_rate.value)
        vp[4].find('.//diffusion_coefficient').text = str(self.pH_diffusion_coefficient.value)
        vp[4].find('.//decay_rate').text = str(self.pH_decay_rate.value)

        uep.find('.//options//calculate_gradients').text = str(self.calculate_gradient.value)
        uep.find('.//options//track_internalized_substrates_in_each_agent').text = str(self.track_internal.value)
