# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:33:18 2016

"""
import os
# open source module
import sys
import ctypes
from PyQt4.uic import loadUiType
from PyQt4 import QtCore, QtGui, Qt
import PyQt4.Qwt5 as Qwt
#from matplotlib import colors
from matplotlib.patches import Ellipse
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import numpy as np
import time
Ui_MainWindow, QMainWindow = loadUiType('Experiment_GUI.ui')

# personal module
os.chdir('Z:/Users/Maxime Joos/Python_wd/')
from PyAPT import APTMotor
from Experiment import Experiment
import hamamatsu_camera as hc
from pyduino import Arduino
from MyFitFunctions import Satlmfit
from SupportFunctions import ConvertImage, GetColorTable

class Thread(QtCore.QThread):
    def __init__(self, fun):
        super(Thread, self).__init__()
        self.fun = fun
    def run(self):
        self.fun()

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        # Create the main window
        super(Main, self).__init__()
        self.setupUi(self)
        
        self.Exp = Experiment()# lunch the experiment
        self.list_connected_APT()
        if self.Exp.hcam.n_cameras == 0: # If no Hammamatsu camera connected
            self.checkBox_Camera.setChecked(False) # Uncheck the Qcheckbox.
        
        if self.Exp.SPD.number_devices == 1: # If Aurea TDC is connected ...
            self.checkBox_TDC.setChecked(True) # Check the Qcheckbox.
            
        self.Import_parameters()
        self.Save_path = self.lineEdit_Save_path.text()
        
        #Link each button of the graphical interface to a specific function
        #Polarisation Tab
        self.lineEdit_hwp_start.editingFinished.connect(self.Set_hwp_start)
        self.lineEdit_hwp_stop.editingFinished.connect(self.Set_hwp_stop)
        self.lineEdit_hwp_step.editingFinished.connect(self.Set_hwp_step)
        self.lineEdit_qwp_start.editingFinished.connect(self.Set_qwp_start)
        self.lineEdit_qwp_stop.editingFinished.connect(self.Set_qwp_stop)
        self.lineEdit_qwp_step.editingFinished.connect(self.Set_qwp_step)
        self.lineEdit_hwp_out_start.editingFinished.connect(self.Set_hwp_out_start)
        self.lineEdit_hwp_out_stop.editingFinished.connect(self.Set_hwp_out_stop)
        self.lineEdit_hwp_out_step.editingFinished.connect(self.Set_hwp_out_step)
        self.lineEdit_x0.editingFinished.connect(self.Set_x0)
        self.lineEdit_y0.editingFinished.connect(self.Set_y0)
        self.lineEdit_width.editingFinished.connect(self.Set_width)
        self.lineEdit_height.editingFinished.connect(self.Set_height)
        self.lineEdit_PG_x0.editingFinished.connect(self.Set_PG_x0)
        self.lineEdit_PG_y0.editingFinished.connect(self.Set_PG_y0)
        self.lineEdit_PG_x1.editingFinished.connect(self.Set_PG_x1)
        self.lineEdit_PG_y1.editingFinished.connect(self.Set_PG_y1)
        self.lineEdit_PG_x0_V.editingFinished.connect(self.Set_PG_V_x0)
        self.lineEdit_PG_y0_V.editingFinished.connect(self.Set_PG_V_y0)
        self.lineEdit_PG_x1_V.editingFinished.connect(self.Set_PG_V_x1)
        self.lineEdit_PG_y1_V.editingFinished.connect(self.Set_PG_V_y1)
        self.lineEdit_Exp_time.editingFinished.connect(self.Set_Exp_time)
        self.lineEdit_Exp_time_PG.editingFinished.connect(self.Set_Exp_time_PG)
        self.lineEdit_Exp_time_PG_V.editingFinished.connect(self.Set_Exp_time_PG_V)
        self.lineEdit_hcam_max.editingFinished.connect(self.Set_hcam_max)
        self.lineEdit_hcam_min.editingFinished.connect(self.Set_hcam_min)
        self.lineEdit_PG45_max.editingFinished.connect(self.Set_PG45_max)
        self.lineEdit_PG45_min.editingFinished.connect(self.Set_PG45_min)
        self.lineEdit_PGvert_max.editingFinished.connect(self.Set_PGvert_max)
        self.lineEdit_PGvert_min.editingFinished.connect(self.Set_PGvert_min)
        
        self.pushButton_Snapshot.clicked.connect(self.Snapshot_raw)
        self.pushButton_Dark_substracted.clicked.connect(self.Snapshot_Dark)
        self.pushButton_Stokes.clicked.connect(self.Measure_Stokes)
        self.pushButton_polarisation_map.clicked.connect(self.Polarisation_map)
        self.pushButton_Stop_map.clicked.connect(self.Stop_Polarisation_map)
        self.pushButton_Save.clicked.connect(self.Save)
        self.pushButton_Start_Csignal.clicked.connect(self.Start_Csignal)
        self.pushButton_Stop_Csignal.clicked.connect(self.Stop_Csignal)
        self.pushButton_hwp_fextremum.clicked.connect(self.hwp_fextremum)
        self.pushButton_qwp_fextremum.clicked.connect(self.qwp_fextremum)
        self.pushButton_calibrate_wp.clicked.connect(self.calibrate_wp)
        self.pushButton_Psat.clicked.connect(self.Psat_measurement)
        self.pushButton_Fit_sat.clicked.connect(self.Fit_sat)
        self.pushButton_start_apd.clicked.connect(self.start_APD)
        self.pushButton_stop_apd.clicked.connect(self.stop_APD)
        
        self.comboBox_select_camera.currentIndexChanged.connect(self.Set_CurrentCam)
        self.comboBox_resolution.currentIndexChanged.connect(self.Set_Plot_resolution)
        self.CurrentCam = self.comboBox_select_camera.currentText() # String name of the current camera
        self.Camera_1 = self.comboBox_select_camera.itemText(0) # Hamamatsu camera
        self.Camera_2 = self.comboBox_select_camera.itemText(1) # PointGrey camera
        self.Camera_3 = self.comboBox_select_camera.itemText(2) # PointGrey camera
        self.Camera_all = self.comboBox_select_camera.itemText(3) # PointGrey camera
        self.checkBox_normalisation.clicked.connect(self.Switch_normalisation)
#        self.checkBox_Dark.clicked.connect(self.switch_dark)
        # Devices Tab
        self.checkBox_Shutter.clicked.connect(self.Switch_Shutter)
        self.checkBox_Shutter_out.clicked.connect(self.Switch_Shutter_out)
        self.checkBox_Shutter_microscope.clicked.connect(self.Switch_Shutter_microscope)
        self.checkBox_Arduino_Pin_6.clicked.connect(self.Switch_Pin_6)
        self.checkBox_Arduino.clicked.connect(self.Switch_Arduino)
        self.checkBox_Camera.clicked.connect(self.Switch_camera)
        self.checkBox_Camera_PG.clicked.connect(self.Switch_camera_PG)
        self.checkBox_APT_hwp.clicked.connect(self.Switch_APT_hwp)
        self.checkBox_APT_qwp.clicked.connect(self.Switch_APT_qwp)
        self.checkBox_APT_hwp_out.clicked.connect(self.Switch_APT_hwp_out)
        self.checkBox_APT_hwp_laser.clicked.connect(self.Switch_APT_hwp_laser)
        self.checkBox_Laser.clicked.connect(self.Switch_Laser)
        self.checkBox_TDC.clicked.connect(self.switch_SPD)
        self.listWidget_APT.itemClicked.connect(self.change_APT)
        self.pushButton_Home.clicked.connect(self.Home)
        self.pushButton_mAbs.clicked.connect(self.MoveAbs)
        self.pushButton_mRel.clicked.connect(self.MoveRel)
        self.pushButton_getPos.clicked.connect(self.GetPos)
        
        # Parameters tab
        self.lineEdit_Arduino_port.editingFinished.connect(self.Set_Arduino_port)
        self.lineEdit_Save_path.editingFinished.connect(self.Set_Save_path)
        self.lineEdit_Channel.editingFinished.connect(self.Set_Oscillo_channel)
        self.lineEdit_Acquisition_time.editingFinished.connect(self.Set_Oscillo_Acquisition_time)
        self.lineEdit_APD_index.editingFinished.connect(self.Set_APD_index)
        self.lineEdit_SPD_sample_time.editingFinished.connect(self.Set_SPD_sample_time)
        
        # Acquisition thread for polarisation map
        self.acquisition_thread = Thread(self.Acquire_map)
        self.acquiring_map = False
        self.acquisition_thread.finished.connect(self.callback)
        
        # Acquisition thread for the ROI Camera signal
        self.acquisition_thread_Csignal = Thread(self.Acquire_Csignal)
        self.acquiring_Csignal = False
        self.acquisition_thread_Csignal.finished.connect(self.callback_Csignal)
        self.N_Csignal = 100
        self.Csignal_hcam = np.zeros((self.N_Csignal), dtype = np.float32)
        self.Csignal_PGcam = np.zeros((self.N_Csignal), dtype = np.float32)
        self.Csignal_PGcam_V = np.zeros((self.N_Csignal), dtype = np.float32)
        self.C_x_axis = np.array(range(self.N_Csignal))
        
        # Acquisition thread for calibration of the waveplate
        self.acquisition_thread_calibration = Thread(self.Acquire_calibration)
        self.acquiring_calibration = False
        self.acquisition_thread_calibration.finished.connect(self.callback_calibration)
        
        # Acquisition thread for Stokes measurement
        self.acquisition_thread_Stokes = Thread(self.Acquire_Stokes)
        self.acquiring_Stokes = False
        self.acquisition_thread_Stokes.finished.connect(self.callback_Stokes)
        
        # Acquisition thread for Psat measurement
        self.acquisition_thread_Psat = Thread(self.Acquire_Psat)
        self.acquiring_Psat = False
        self.acquisition_thread_Psat.finished.connect(self.callback_Psat)
        
        # Acquisition thread for APD count rate ploting
        self.N_APD = 100
        self.APD1_CountRate = np.zeros(self.N_APD)
        self.APD2_CountRate = np.zeros(self.N_APD)
        self.acquisition_thread_APD = Thread(self.Acquire_APD)
        self.acquiring_APD = False
        self.acquisition_thread_APD.finished.connect(self.callback_APD)
        
        #setup the hamamatsu Camera signal plot
        self.curve_hcam = Qwt.QwtPlotCurve()
        self.curve_hcam.attach(self.qwtPlot_hcam)
        #setup the PG Camera signal plot
        self.curve_PGcam = Qwt.QwtPlotCurve()
        self.curve_PGcam.attach(self.qwtPlot_PGcam)
        #setup the PG Camera signal plot
        self.curve_PGcam_V = Qwt.QwtPlotCurve()
        self.curve_PGcam_V.attach(self.qwtPlot_PGcam_V)
        #setup the saturation curve 
        self.curve_sat = Qwt.QwtPlotCurve()
        self.curve_sat.attach(self.qwtPlot_hcam)
        #setup the saturation Fit curve 
        self.curve_Sat_fit = Qwt.QwtPlotCurve()
        self.curve_Sat_fit.attach(self.qwtPlot_hcam)
#        self.curve.setPen(Qt.QPen(Qt.Qt.red, 3))
        
        #setup the polarisation resolved measurement plot
        self.curve_pola = Qwt.QwtPlotCurve()
        self.curve_pola.attach(self.qwtPlot_pola)
        self.curve_pola_fit = Qwt.QwtPlotCurve()
        self.curve_pola_fit.attach(self.qwtPlot_pola)
        
        fig_map = Figure()
        self.canvas_map = FigureCanvas(fig_map)
        self.mplvl_map.addWidget(self.canvas_map)
        self.canvas_map.draw()
        self.toolbar_map = NavigationToolbar(self.canvas_map, self.mplwindow_map, coordinates=True)
        self.mplvl_map.addWidget(self.toolbar_map)

        fig_map_Sd = Figure()
        self.canvas_map_Sd = FigureCanvas(fig_map_Sd)
        self.mplvl_map_Sd.addWidget(self.canvas_map_Sd)
        self.canvas_map_Sd.draw()
        self.toolbar_map_Sd = NavigationToolbar(self.canvas_map_Sd, self.mplwindow_map_Sd, coordinates=True)
        self.mplvl_map_Sd.addWidget(self.toolbar_map_Sd)

        fig_map_St = Figure()
        self.canvas_map_St = FigureCanvas(fig_map_St)
        self.mplvl_map_St.addWidget(self.canvas_map_St)
        self.canvas_map_St.draw()
        self.toolbar_map_St = NavigationToolbar(self.canvas_map_St, self.mplwindow_map_St, coordinates=True)
        self.mplvl_map_St.addWidget(self.toolbar_map_St)
        
    def Acquire_APD(self):
        CountRate_APD1,  CountRate_APD2 = self.Exp.getCountRateData(self.APD_plot_sample_time, 1, False)[0] / self.APD_plot_sample_time # [Counts/s]
        self.APD1_CountRate = np.roll(self.APD1_CountRate,-1) # Shift all points to the left
        self.APD2_CountRate = np.roll(self.APD2_CountRate,-1) # Shift all points to the left
        self.APD1_CountRate[self.N_APD-1] = CountRate_APD1 # Update data
        self.APD2_CountRate[self.N_APD-1] = CountRate_APD2 # Update data
        
    def Acquire_calibration(self):
        self.pos_hwp_out = range(self.hwp_out_start_Pos, self.hwp_out_stop_Pos, self.hwp_out_step)
        Oscillo_name = self.comboBox_select_source_calibration.currentText()
        dev = getattr(self.Exp, Oscillo_name)
        Ch = int(self.comboBox_Channel.currentText())
#        Angle = int(self.lineEdit_Analyser_angle.text())
        motor_name = self.comboBox_APT.currentText()
        RM = getattr(self.Exp, motor_name)
        dark = self.checkBox_Dark.isChecked()
        shutter_pin = int(self.lineEdit_shutter_pin.text())
        self.Exp.calibrate_wp(self.pos_hwp_out, dev, Ch, RM, dark, shutter_pin)
        
    def Acquire_Csignal(self): # acquire the signal in the thread
        time.sleep(0.05)
        if self.CurrentCam == self.Camera_1 or self.CurrentCam == self.Camera_all: # if at least Hamamatsu camera is selected ...
            pic = self.Exp.GrabHamamatsuPic(self.Exp.hcam)
            S_hcam = np.sum(pic)
            self.Csignal_hcam = np.roll(self.Csignal_hcam,-1) # Roll all elements to the left
            self.Csignal_hcam[-1::] = S_hcam # Insert latest signal measurement at the end of array
            if not self.normalisation: self.lineEdit_Csignal_hcam.setText("{0:.{1}e}".format(S_hcam, 2))
            else:
                S_hcam_norm = (S_hcam - self.S_hcam_min)/(self.S_hcam_max - self.S_hcam_min) # normalised signal
                self.lineEdit_Csignal_hcam.setText(str(np.round(S_hcam_norm,2)))
            
        if self.CurrentCam == self.Camera_2 or self.CurrentCam == self.Camera_all: # if at least PG camera is selected ...
            pic = self.Exp.PGcam.camera.GrabNumPyImage('gray')
            pic = pic[self.Exp.PGcam_roi[1]:self.Exp.PGcam_roi[3], self.Exp.PGcam_roi[0]:self.Exp.PGcam_roi[2]]
            S_PGcam = np.sum(pic)
            self.Csignal_PGcam = np.roll(self.Csignal_PGcam,-1) # Roll all elements to the left
            self.Csignal_PGcam[-1::] = S_PGcam # Insert latest signal measurement at the end of array
#            self.lineEdit_Csignal_PGcam.setText("{0:.{1}e}".format(S_PGcam, 2))
            if not self.normalisation: self.lineEdit_Csignal_PGcam.setText("{0:.{1}e}".format(S_PGcam, 2))
            else:
                S_PGcam_norm = (S_PGcam - self.S_PG45_min)/(self.S_PG45_max - self.S_PG45_min) # normalised signal
                self.lineEdit_Csignal_PGcam.setText(str(np.round(S_PGcam_norm,2)))
            
        if self.CurrentCam == self.Camera_3 or self.CurrentCam == self.Camera_all: # if at least PG camera is selected ...
            pic = self.Exp.PGcam_V.camera.GrabNumPyImage('gray')
            pic = pic[self.Exp.PGcam_V_roi[1]:self.Exp.PGcam_V_roi[3], self.Exp.PGcam_V_roi[0]:self.Exp.PGcam_V_roi[2]]
            S_PGcam = np.sum(pic)
            self.Csignal_PGcam_V = np.roll(self.Csignal_PGcam_V,-1) # Roll all elements to the left
            self.Csignal_PGcam_V[-1::] = S_PGcam # Insert latest signal measurement at the end of array
#            self.lineEdit_Csignal_PGcam_V.setText("{0:.{1}e}".format(S_PGcam, 2))
            if not self.normalisation: self.lineEdit_Csignal_PGcam_V.setText("{0:.{1}e}".format(S_PGcam, 2))
            else:
                S_PGcam_norm = (S_PGcam - self.S_PGvert_min)/(self.S_PGvert_max - self.S_PGvert_min) # normalised signal
                self.lineEdit_Csignal_PGcam_V.setText(str(np.round(S_PGcam_norm,2)))
        
    def Acquire_Degree_of_pola(self):
        time.sleep(0.05)
        self.Exp.Shutter('close')
        self.Exp.SPD.SetVisibleModuleDetectionMode(self.APD_index, 0) # Continuous mode
        Dark_cnt = self.Exp.getCountRateData(self.Exp.SPD_sample_time, 1)
        for i in range(len(self.pos_hwp_out)):
            self.Exp.RMHWP_out.mAbs(self.pos_hwp_out[i])
            self.Exp.wait_motor(self.Exp.RMHWP_out)
            Fluo_cnt = self.Exp.getCountRateData(self.Exp.SPD_sample_time, 1)
            self.S_degree_of_pola[i] = Fluo_cnt[0][self.APD_index]
        self.S_degree_of_pola = self.S_degree_of_pola - Dark_cnt[0][self.APD_index]
        
    def Acquire_map(self): # acquire the signal in the thread
        time.sleep(0.05)
        self.Exp.RMQWP.mAbs(self.pos_qwp[self.cntr])
        d, t = self.Exp.Scan_polarisation_1axis(self.pos_hwp, self.Exp.RMHWP)
        self.Sd = np.column_stack((self.Sd, d))
        self.St = np.column_stack((self.St, t))
        if self.cntr == 0:
            self.Sd = self.Sd[:, 1::] # Remove first dummy raw
            self.St = self.St[:, 1::] # Remove first dummy raw
        self.Sn = self.Sd/self.St
        
    def Acquire_Psat(self):
        self.hwp_laser_start_Pos = int(self.from_lineEdit('lineEdit_hwp_laser_start'))
        self.hwp_laser_stop_Pos = int(self.from_lineEdit('lineEdit_hwp_laser_stop'))
        self.hwp_laser_step = int(self.from_lineEdit('lineEdit_hwp_laser_step'))
#        Pos = range(self.hwp_laser_start_Pos, self.hwp_laser_stop_Pos, self.hwp_laser_step)
        Pos = np.logspace(np.log10(self.hwp_laser_start_Pos), np.log10(self.hwp_laser_stop_Pos), 20)
        Oscillo_name = self.comboBox_select_source_saturation.currentText()
        dev = getattr(self.Exp, Oscillo_name) # desired oscilloscope
        Ch = int(self.comboBox_Input_S.currentText()) # desired channel
        Sout_Source = self.comboBox_Output_S.currentText()
        shutter_pin = int(self.lineEdit_shutter_pin.text())
        N_samples = int(self.from_lineEdit('lineEdit_Saturation_N_samples'))
        N_bright_points = int(self.from_lineEdit('lineEdit_Saturation_bright_points'))
        self.Sin_Psat, self.Sout_Psat = self.Exp.Psat(Pos, dev, Ch, Sout_Source, shutter_pin, N_samples, N_bright_points)
        
    def Acquire_Stokes(self):
        shutter_pin = int(self.lineEdit_shutter_pin.text())
        motor_name = self.comboBox_APT.currentText()
        RM = getattr(self.Exp, motor_name)
        N_samples = int(self.from_lineEdit('lineEdit_Saturation_N_samples'))
        N_bright_points = int(self.from_lineEdit('lineEdit_Saturation_bright_points'))
        self.Stokes_raw = self.Exp.getStokes(self.pos_hwp_out, self.retardance, self.phi_0, self.comboBox_select_source.currentText(), shutter_pin, RM, N_samples, N_bright_points)
        
    def addAPT(self, l):
        try:
            self.APT_dict[l] = getattr(self.Exp, l)
            self.listWidget_APT.addItem(l)
        except: pass

    def addmpl_map(self, fig_map): #Plot
        self.canvas_map = FigureCanvas(fig_map)
        self.mplvl_map.addWidget(self.canvas_map)
        self.canvas_map.draw()

    def addmpl_map_Sd(self, fig_map): #Plot
        self.canvas_map_Sd = FigureCanvas(fig_map)
        self.mplvl_map_Sd.addWidget(self.canvas_map_Sd)
        self.canvas_map_Sd.draw()
        
    def addmpl_map_St(self, fig_map): #Plot
        self.canvas_map_St = FigureCanvas(fig_map)
        self.mplvl_map_St.addWidget(self.canvas_map_St)
        self.canvas_map_St.draw()

    def addmpl_image(self): #Plot
        self.canvas_image = FigureCanvas(self.fig_image)
        self.mplvl_image.addWidget(self.canvas_image)
        self.canvas_image.draw()
        self.toolbar_image = NavigationToolbar(self.canvas_image, self.mplwindow_image, coordinates=True)
        self.mplvl_image.addWidget(self.toolbar_image)

    def BuildColorbar(self, ScaleWidgetName, Min, Max):
        ScaleWidget = getattr(self, ScaleWidgetName)
        # Define colors of the colorbar
        CM = Qwt.QwtLinearColorMap(Qt.Qt.blue, Qt.Qt.red)
        CM.addColorStop(0.333, Qt.Qt.cyan)
        CM.addColorStop(0.667, Qt.Qt.yellow)
        ScaleWidget.setColorMap(Qwt.QwtDoubleInterval(Min, Max), CM)        
        # graduation of the scale
        Scale_E = Qwt.QwtLinearScaleEngine()
        
        ScaleWidget.setScaleDiv(Scale_E.transformation(), Scale_E.divideScale(Min-(1/(Max-Min+1))*0.5, Max+(1/(Max-Min+1))*0.5, 10, 0))
        #Draw the colorbar
        ScaleWidget.setColorBarEnabled(True)
        
    def calibrate_wp(self):
        self.acquiring_calibration = True
        self.acquisition_thread_calibration.start()    
        
    def callback(self): # plot in the callback
        self.cntr+=1
        self.Plot_Sn()
        self.Plot_Sd()
        self.Plot_St()
        if self.cntr >= len(self.pos_qwp): # if all positions have been measured
            self.acquiring_map = False
#            if self.checkBox_Save_and_close.isChecked():
#                self.Save()
#                self.Exp.Laser('off')
#                self.Exp.Shutter('open')
        if self.acquiring_map:
            self.acquisition_thread.start() # restart the polarisation map acquisition.
            
    def callback_APD(self):
        if self.acquiring_APD:
            self.Plot_curve_hcam( np.arange(0, self.N_APD), self.APD1_CountRate)
            self.Plot_curve_PGcam( np.arange(0, self.N_APD), self.APD2_CountRate)
            self.acquisition_thread_APD.start() # Restart acquisition thread
            
    def callback_calibration(self):
        self.lineEdit_retardance.setText(str(np.round(self.Exp.retardance, 4)))
        self.retardance = self.Exp.retardance
        self.phi_0 = self.Exp.FastAxis_origin/4.0
        self.lineEdit_phi_0.setText(str(np.round(self.Exp.FastAxis_origin/4.0, 4)))
        self.plot_pola(self.pos_hwp_out, self.Exp.S_cal, self.Exp.Fit_cal)
        self.acquiring_calibration = False
            
    def callback_Csignal(self):
        if self.CurrentCam == self.Camera_1 or self.CurrentCam == self.Camera_all: # if at least Hamamatsu camera is selected ...
            self.Plot_curve_hcam(self.C_x_axis, self.Csignal_hcam)
        if self.CurrentCam == self.Camera_2 or self.CurrentCam == self.Camera_all: # if at least PointGrey camera is selected ...
            self.Plot_curve_PGcam(self.C_x_axis, self.Csignal_PGcam)
        if self.CurrentCam == self.Camera_3 or self.CurrentCam == self.Camera_all: # if at least PointGrey camera is selected ...
            self.Plot_curve_PGcam_V(self.C_x_axis, self.Csignal_PGcam_V)
        
        if self.acquiring_Csignal:
            self.acquisition_thread_Csignal.start()

    def callback_Psat(self):
        self.Plot_Psat(self.Sin_Psat, self.Sout_Psat)
        
    def callback_Stokes(self):
        self.Stokes_norm = self.Stokes_raw / self.Stokes_raw[0]
        self.Degree_of_pol = np.sqrt(self.Stokes_raw[1]**2 + (self.Stokes_raw[2])**2 + (self.Stokes_raw[3])**2)/self.Stokes_raw[0]
        self.Degree_of_lin_pol = np.sqrt(self.Stokes_norm[1]**2 + (self.Stokes_norm[2])**2)
        
        self.lineEdit_Degree_of_pol.setText(str(np.round(self.Degree_of_pol, 2)))
        self.lineEdit_Degree_of_lin_pol.setText(str(np.round(self.Degree_of_lin_pol, 2)))
        self.lineEdit_Stokes_0.setText(str(np.round(self.Stokes_raw[0], 2)))
        self.lineEdit_Stokes_1.setText(str(np.round(self.Stokes_raw[1], 2)))
        self.lineEdit_Stokes_2.setText(str(np.round(self.Stokes_raw[2], 2)))
        self.lineEdit_Stokes_3.setText(str(np.round(self.Stokes_raw[3], 2)))
#        self.lineEdit_Stokes_0_norm.setText(str(np.round(main.Stokes_norm[0], 2)))
        self.lineEdit_Stokes_1_norm.setText(str(np.round(self.Stokes_norm[1], 2)))
        self.lineEdit_Stokes_2_norm.setText(str(np.round(self.Stokes_norm[2], 2)))
        self.lineEdit_Stokes_3_norm.setText(str(np.round(self.Stokes_norm[3], 2)))
        self.lineEdit_Stokes_1_std.setText(str(np.round(self.Exp.Std_Stokes[1], 3)))
        self.lineEdit_Stokes_2_std.setText(str(np.round(self.Exp.Std_Stokes[2], 3)))
        self.lineEdit_Stokes_3_std.setText(str(np.round(self.Exp.Std_Stokes[3], 3)))
#        self.Plot_Stokes()
        self.plot_pola(self.pos_hwp_out, self.Exp.I_t, self.Exp.Fit_Stokes)
        
    def change_APT(self, item): # Change the current APT motor
        APTname = item.text() # read the selected list element 
        self.CurrentAPT = self.Exp.connected_apt_dict[APTname] # Search in the dictionary the corresponding instance
        
    def closeCamera(self):
        self.Exp.hcam.shutdown()
        self.Exp.hcam.terminate()
        
    def closeEvent(self, event):
        self.Exp.CloseExp()
        
    def Connect_camera(self):
        self.Exp.hcam = hc.HamamatsuCamera()
        
    def Fit_sat(self):
        self.Fit, self.par, self.par_guess, self.par_Std = Satlmfit(self.Sin_Psat, self.Sout_Psat)
        self.lineEdit_Psat.setText(str(np.round(self.par[1], 3)))
        self.Plot_fit_Psat()
    
    def from_lineEdit(self, lineEditName):
        """
        lineEdit : String. QLineEdit name in the Qt .ui file.
        """
        return getattr(self, lineEditName).text() # return the string text of desired QLineEdit widget
        
    def GetPos(self):
        pos = self.CurrentAPT.getPos()
        self.lineEdit_getPos.setText(str(np.round(pos, 5)))
    
    def Home(self):
        self.CurrentAPT.go_home()
    
    def hwp_fextremum(self):
        self.Set_hwp_in_red_start()
        self.Set_hwp_in_red_stop()
        self.Set_hwp_in_red_step()
        Pos = range(self.hwp_in_red_start_Pos, self.hwp_in_red_stop_Pos, self.hwp_in_red_step)
        self.S = np.zeros(len(Pos))
        i = 0
        for p in Pos:
            self.Exp.RMHWP_in_red.mAbs(p, True)
            if self.CurrentCam == self.Camera_1 : # if Hamamatsu camera is selected ...
                self.S[i] = np.sum(self.Exp.getImageHamamatsu(self.Exp.hcam))
            
            elif self.CurrentCam == self.Camera_2: # if  PG camera is selected ...
                pic = self.Exp.getImagePointGrey(self.Exp.PGcam, self.Exp.ExpTime_PGcam, self.Exp.PGcam_roi)
                self.S[i] = np.sum(pic)
            
            elif self.CurrentCam == self.Camera_3: # if PG camera is selected ...
                pic = self.Exp.getImagePointGrey(self.Exp.PGcam_V, self.Exp.ExpTime_PGcam_V, self.Exp.PGcam_V_roi)
                self.S[i] = np.sum(pic)
            
            elif self.CurrentCam == self.Camera_all:
                self.ShowMessage('Select a unique Camera !')
                return
            i+=1
        p = np.polyfit(Pos, self.S, 2)
        extremum_pos = -p[1]/(2*p[0])
        Smax = p[2] - p[1]**2/(4*p[0])
        self.Exp.RMHWP_in_red.mAbs(extremum_pos)
        self.Plot_curve_hcam(Pos, self.S)
        self.lineEdit_Csignal_hcam.setText(str(int(extremum_pos)))
        self.lineEdit_Csignal_PGcam_V.setText("{0:.{1}e}".format(Smax, 2))
        self.lineEdit_hwp_in_red_start.setText(str(int(extremum_pos - 13)))
        self.lineEdit_hwp_in_red_stop.setText(str(int(extremum_pos + 17)))
        
    def qwp_fextremum(self):
        self.qwp_in_red_start_Pos = self.from_lineEdit('lineEdit_qwp_in_red_start')
        self.qwp_in_red_stop_Pos = self.from_lineEdit('lineEdit_qwp_in_red_stop')
        self.qwp_in_red_step_Pos = self.from_lineEdit('lineEdit_qwp_in_red_step')
        Pos = range(self.qwp_in_red_start_Pos, self.qwp_in_red_stop_Pos, self.qwp_in_red_step_Pos)
        self.S = np.zeros(len(Pos))
        i = 0
        for p in Pos:
            self.Exp.RMQWP_in_red.mAbs(p, True)
            if self.CurrentCam == self.Camera_1 : # if Hamamatsu camera is selected ...
                self.S[i] = np.sum(self.Exp.getImageHamamatsu(self.Exp.hcam))
            
            if self.CurrentCam == self.Camera_2: # if  PG camera is selected ...
                pic = self.Exp.getImagePointGrey(self.Exp.PGcam, self.Exp.ExpTime_PGcam, self.Exp.PGcam_roi)
                self.S[i] = np.sum(pic)
            
            if self.CurrentCam == self.Camera_3: # if PG camera is selected ...
                pic = self.Exp.getImagePointGrey(self.Exp.PGcam_V, self.Exp.ExpTime_PGcam_V, self.Exp.PGcam_V_roi)
                self.S[i] = np.sum(pic)
            
            elif self.CurrentCam == self.Camera_all:
                self.ShowMessage('Select a unique Camera !')
                return
            i+=1
        p = np.polyfit(Pos, self.S, 2)
        extremum_pos = -p[1]/(2*p[0])
        Smax = p[2] - p[1]**2/(4*p[0])
        self.Exp.RMQWP_in_red.mAbs(extremum_pos)
        self.Plot_curve_PGcam(Pos, self.S)
        self.lineEdit_Csignal_PGcam.setText(str(int(extremum_pos)))
        self.lineEdit_Csignal_PGcam_V.setText("{0:.{1}e}".format(Smax, 2))
        self.lineEdit_hwp_laser_start.setText(str(int(extremum_pos - 15)))
        self.lineEdit_hwp_laser_stop.setText(str(int(extremum_pos + 17)))
        
    def Import_parameters(self):
        # Polarisation tab
        self.Set_hwp_start()
        self.Set_hwp_stop()
        self.Set_hwp_step()
        self.Set_qwp_start()
        self.Set_qwp_stop()
        self.Set_qwp_step()
        self.Set_hwp_out_start()
        self.Set_hwp_out_stop()
        self.Set_hwp_out_step()
        # Hamamatsu ROI
        self.Set_x0()
        self.Set_y0()
        self.Set_width()
        self.Set_height()
        self.Set_Exp_time()
        # PG ROI
        self.Set_PG_x0()
        self.Set_PG_y0()
        self.Set_PG_x1()
        self.Set_PG_y1()
        self.Set_Exp_time_PG()
        # PG ROI
        self.Set_PG_V_x0()
        self.Set_PG_V_y0()
        self.Set_PG_V_x1()
        self.Set_PG_V_y1()
        self.Set_Exp_time_PG_V()
        # Parameters Tab
        self.lineEdit_SPD_sample_time.setText(str(self.Exp.SPD_sample_time))
        
        self.Set_Plot_resolution()
        self.normalisation = self.checkBox_normalisation.isChecked()
        
    def list_connected_APT(self):
        for key in self.Exp.connected_apt_dict: # For each connected apt controller ...
            self.listWidget_APT.addItem(key) #add APT label to the list widget of the GUI.
            self.comboBox_APT.addItem(key)
            
    def Measure_Stokes(self):
        if self.acquiring_Csignal: #If signal acquisition is going on ...
            self.Stop_Csignal() # Stop it.
#        self.Exp.Shutter('open')
        self.pos_hwp_out = range(self.hwp_out_start_Pos, self.hwp_out_stop_Pos, self.hwp_out_step)
        self.Set_retardance()
        self.Set_phi_0()
        self.acquisition_thread_Stokes.start()
        
    def MoveAbs(self):
        self.CurrentAPT.mAbs(float(self.lineEdit_mAbs.text()))

    def MoveRel(self):
        self.CurrentAPT.mRel(float(self.lineEdit_mRel.text()))
    
    def plot_pola(self, x, y, fit):
        self.refresh_qwtPlot_pola()
        #Plot the measurement points
        self.curve_pola.setData(x, y)
        self.curve_pola.setStyle(Qwt.QwtPlotCurve.Dots)
        self.curve_pola.setPen(Qt.QPen(Qt.Qt.blue, 3))
        #Plot the fit
        self.curve_pola_fit.setData(x, fit)
        self.curve_pola_fit.setPen(Qt.QPen(Qt.Qt.red, 2))
        # Set axis labels
        self.qwtPlot_pola.setAxisTitle(0, u'Transmitted Signal (V)')
        self.qwtPlot_pola.setAxisTitle(2, u'WP position (°)')
        #Plot
        self.qwtPlot_pola.replot()

#    def Plot_camera(self, im):
#        if self.plot_factor != 1:
#            im = self.reshape_image(np.array(im), self.plot_factor)        
#
#        self.h.set_array(im)
#        if self.CurrentCam == self.Camera_1:
##            self.h.set_array(im)
#            self.h.set_extent([self.Exp.x0, self.Exp.x0 + self.Exp.width , self.Exp.y0 + self.Exp.height, self.Exp.y0])
##            a = image.imshow(im, interpolation='None', extent = [self.Exp.x0, self.Exp.x0 + self.Exp.width , self.Exp.y0 + self.Exp.height, self.Exp.y0])
#        elif self.CurrentCam == self.Camera_2:
#            self.h.set_extent([self.Exp.PGcam_roi[0], self.Exp.PGcam_roi[2] , self.Exp.PGcam_roi[3], self.Exp.PGcam_roi[1]])
##            a = image.imshow(im, interpolation='None', extent = [self.Exp.PG_x0, self.Exp.PG_x1 , self.Exp.PG_y1, self.Exp.PG_y0])
##            a = image.imshow(im, interpolation='None', extent = [self.Exp.PGcam_roi[0], self.Exp.PGcam_roi[2] , self.Exp.PGcam_roi[3], self.Exp.PGcam_roi[1]])
#        elif self.CurrentCam == self.Camera_3:
#            self.h.set_extent([self.Exp.PGcam_V_roi[0], self.Exp.PGcam_V_roi[2] , self.Exp.PGcam_V_roi[3], self.Exp.PGcam_V_roi[1]])
#
#        self.h.autoscale()
#        self.rmmpl_image()
##        self.addmpl_image(fig_image)
#        self.addmpl_image()

    def Plot_camera(self, im, LabelWidgetName): # plot image using setPixmap instance of a QLabel widget.
        LabelWidget = getattr(self, LabelWidgetName)
        if self.plot_factor != 1:
            im = self.reshape_image(np.array(im), self.plot_factor)

        Imax = np.max(im)
        Imin = np.min(im)
        Im_width = im.shape[1] # Number of columns
        Im_height = im.shape[0] # Nummber of rows
        
        if im.dtype.name != 'uint8': im = ConvertImage(im) # If image is not of uint8 type, convert it to uint8.
        
        im = np.array(im-np.min(im), copy = True)
        
        ColorTable = GetColorTable(np.max(im)-np.min(im)+1)
        QI = QtGui.QImage(im.data, Im_width, Im_height, Im_width, QtGui.QImage.Format_Indexed8)
        QI.setColorTable(ColorTable)
        LabelWidget.setPixmap(QtGui.QPixmap.fromImage(QI))
        self.BuildColorbar('ScaleWidget_ColorBar_Camera', Imin, Imax)

#        if self.CurrentCam == self.Camera_1:
#            im = ConvertImage(im)
#            ColorTable = GetColorTable(Imax, Imin)
#            QI = QtGui.QImage(im.data, Im_width, Im_height, QtGui.QImage.Format_Indexed8)
#            QI.setColorTable(ColorTable)
#            self.label_50.setPixmap(QtGui.QPixmap.fromImage(QI))
#            self.BuildColorbar('ScaleWidget', Imin, Imax)
#
#        elif self.CurrentCam == self.Camera_2:
#            ColorTable = GetColorTable(Imax, Imin)
#            QI = QtGui.QImage(im.data, Im_width, Im_height, QtGui.QImage.Format_Indexed8)
#            QI.setColorTable(ColorTable)
#            self.label_50.setPixmap(QtGui.QPixmap.fromImage(QI))
#            self.BuildColorbar('ScaleWidget', Imin, Imax)
#
#        elif self.CurrentCam == self.Camera_3:

#        self.h.autoscale()
#        self.rmmpl_image()
##        self.addmpl_image(fig_image)
#        self.addmpl_image()

    def Plot_curve_hcam(self, x, y):
        self.curve_hcam.setData(x, y)
        self.curve_hcam.setPen(Qt.QPen(Qt.Qt.red, 3))
        self.qwtPlot_hcam.replot()

    def Plot_curve_PGcam(self, x, y):
        self.curve_PGcam.setData(x, y)
        self.curve_PGcam.setPen(Qt.QPen(Qt.Qt.red, 3))
        self.qwtPlot_PGcam.replot()
        
    def Plot_curve_PGcam_V(self, x, y):
        self.curve_PGcam_V.setData(x, y)
        self.curve_PGcam_V.setPen(Qt.QPen(Qt.Qt.red, 3))
        self.qwtPlot_PGcam_V.replot()
    
    def Plot_fit_Psat(self):
        self.curve_Sat_fit.setData(self.Sin_Psat, self.Fit)
#        self.curve_Sat_fit.setStyle(Qwt.QwtPlotCurve.Lines )
        self.curve_Sat_fit.setPen(Qt.QPen(Qt.Qt.red, 3))
        self.qwtPlot_hcam.replot()
        
    def Plot_ellipse(self, S1, S2, S3):
        Dp = np.sqrt(S1**2 + S2**2 + S3**2) # Degree of polarisation
        if S1 > 0 : Psi = np.arctan(S2/S1) / 2
        elif S1 < 0: Psi = np.arctan(S2/S1) / 2 + np.pi/2
        elif S1 == 0: Psi = np.pi/4 - np.sign(S2) * np.pi/4 
        Chi = np.arctan(S3/np.sqrt(S1**2 + S2**2)) / 2

        grand_axe = Dp*np.cos(Chi)
        petit_axe = Dp*np.sin(Chi)

        ell = Ellipse(xy = np.zeros(2), width = petit_axe, height = grand_axe, angle = Psi*180/np.pi - 90)
        fig_ellipse = Figure()
        ax_ellipse = fig_ellipse.add_subplot(111, aspect = 'equal')
        ax_ellipse.add_artist(ell)
        ell.set_facecolor('w')

        ax_ellipse.set_xlim(-0.5, 0.5)
        ax_ellipse.set_ylim(-0.5, 0.5)
        self.rmmpl_map()
        self.addmpl_map(fig_ellipse)
        
    def Plot_Psat(self, x, y):
        self.curve_sat.setData(x, y)
        self.curve_sat.setStyle(Qwt.QwtPlotCurve.Dots)
        self.curve_sat.setPen(Qt.QPen(Qt.Qt.blue, 3))
        self.qwtPlot_hcam.setAxisTitle(0, u'Pout (counts/s)')
        self.qwtPlot_hcam.setAxisTitle(2, u'Pin (V)')
        self.qwtPlot_hcam.replot()
        
    def Plot_Sn(self):
        fig_map = Figure()
        pola_map = fig_map.add_subplot(111)
        axes_limit = [self.qwp_start_Pos - self.qwp_step/2., self.qwp_start_Pos + (np.shape(self.Sn)[1]-1)*self.qwp_step + self.qwp_step/2., self.hwp_stop_Pos - self.hwp_step/2., self.hwp_start_Pos - self.hwp_step/2.]
        im_map = pola_map.imshow(self.Sn, interpolation='None', extent = axes_limit)
        fig_map.colorbar(im_map) # plot colorbar
        pola_map.set_title('Normalised Signal')
        pola_map.set_xlabel(u'QWP position (°)', fontsize = 10)
        pola_map.set_ylabel(u'HWP position (°)', fontsize = 10)
        self.rmmpl_map()
        self.addmpl_map(fig_map)
    
    def Plot_Sd(self):
        fig_map = Figure()
        pola_map = fig_map.add_subplot(111)
        axes_limit = [self.qwp_start_Pos - self.qwp_step/2., self.qwp_start_Pos + (np.shape(self.Sd)[1]-1)*self.qwp_step + self.qwp_step/2., self.hwp_stop_Pos - self.hwp_step/2., self.hwp_start_Pos - self.hwp_step/2.]
        im_map = pola_map.imshow(self.Sd, interpolation='None', extent = axes_limit)
        fig_map.colorbar(im_map) # plot colorbar
        pola_map.set_title('Diffusion Signal')
        pola_map.set_xlabel(u'QWP position (°)', fontsize = 10)
        pola_map.set_ylabel(u'HWP position (°)', fontsize = 10)
        self.rmmpl_map_Sd()
        self.addmpl_map_Sd(fig_map)
    
    def Plot_St(self):
        fig_map = Figure()
        pola_map = fig_map.add_subplot(111)
        axes_limit = [self.qwp_start_Pos - self.qwp_step/2., self.qwp_start_Pos + (np.shape(self.St)[1]-1)*self.qwp_step + self.qwp_step/2., self.hwp_stop_Pos - self.hwp_step/2., self.hwp_start_Pos - self.hwp_step/2.]
        im_map = pola_map.imshow(self.St, interpolation='None', extent = axes_limit)
        fig_map.colorbar(im_map) # plot colorbar
        pola_map.set_title('Transmited Signal')
        pola_map.set_xlabel(u'QWP position (°)', fontsize = 10)
        pola_map.set_ylabel(u'HWP position (°)', fontsize = 10)
        self.rmmpl_map_St()
        self.addmpl_map_St(fig_map)
    
    def Polarisation_map(self):
        self.Exp.RMQWP.mAbs(self.qwp_start_Pos)
        self.Exp.RMHWP.mAbs(self.hwp_start_Pos)
        self.pos_qwp = range(self.qwp_start_Pos, self.qwp_stop_Pos, self.qwp_step)
        self.pos_hwp = range(self.hwp_start_Pos, self.hwp_stop_Pos, self.hwp_step)

        self.Sd = np.zeros((len(self.pos_hwp)), dtype = np.float32)
        self.St = np.zeros((len(self.pos_hwp)), dtype = np.float32)
        
        self.cntr=0
        self.acquiring_map = True
        self.acquisition_thread.start()   
        
    def Psat_measurement(self):
        self.refresh_qwtPlot_hcam()
        self.acquiring_Psat = True
        self.acquisition_thread_Psat.start()
        
    def refresh_qwtPlot_hcam(self):
        self.qwtPlot_hcam.detachItems()
        #setup the hamamatsu Camera signal plot
        self.curve_hcam = Qwt.QwtPlotCurve()
        self.curve_hcam.attach(self.qwtPlot_hcam)
        #setup the saturation curve 
        self.curve_sat = Qwt.QwtPlotCurve()
        self.curve_sat.attach(self.qwtPlot_hcam)
        #setup the saturation Fit curve 
        self.curve_Sat_fit = Qwt.QwtPlotCurve()
        self.curve_Sat_fit.attach(self.qwtPlot_hcam)
        
    def refresh_qwtPlot_pola(self):
        self.qwtPlot_pola.detachItems()

        self.curve_pola = Qwt.QwtPlotCurve()
        self.curve_pola.attach(self.qwtPlot_pola)
        self.curve_pola_fit = Qwt.QwtPlotCurve()
        self.curve_pola_fit.attach(self.qwtPlot_pola)
        
    def reshape_image(self, im, plot_factor): # reduce image by a factor plot_factor
        [height, width]= np.shape(im)
        reshaped = np.reshape(im, (height/plot_factor, plot_factor, width/plot_factor, plot_factor))
        pic = np.array(np.sum(np.sum(reshaped, axis=1), axis=2)/plot_factor**2, dtype = np.int32)
        return pic
        
    def rmmpl_map(self): # Clear figure
        self.mplvl_map.removeWidget(self.canvas_map)
        self.canvas_map.close()
        
    def rmmpl_map_Sd(self): # Clear figure
        self.mplvl_map_Sd.removeWidget(self.canvas_map_Sd)
        self.canvas_map_Sd.close()
        
    def rmmpl_map_St(self): # Clear figure
        self.mplvl_map_St.removeWidget(self.canvas_map_St)
        self.canvas_map_St.close()

    def rmmpl_image(self): # Clear figure
        self.mplvl_image.removeWidget(self.canvas_image)
        self.canvas_image.close()
        self.mplvl_image.removeWidget(self.toolbar_image)
        self.toolbar_image.close()
        
    def Save(self):
        import datetime
        today = datetime.date.today()  # get today's date as a datetime type
        todaystr = today.isoformat()   # get string representation: YYYY-MM-DD
        today_folder = os.path.join(self.Save_path, todaystr)
        if not os.path.exists(today_folder): #If today's folder does not exists
            os.makedirs(today_folder) # Create it.
        Time = time.strftime("%Hh%M") # Time at the format HH:MM
        Time_folder = os.path.join(today_folder, Time)
        os.makedirs(Time_folder) # Create time folder
        
        try: np.save(os.path.join(Time_folder,'retardance'), self.retardance)
        except: pass
        try: np.save(os.path.join(Time_folder,'phi_0'), self.phi_0)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes_I_t'), self.Exp.I_t)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes_pos'), self.pos_hwp_out)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes'), self.Stokes_raw)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes_norm'), self.Stokes_norm)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes_std'), self.Exp.Std_Stokes)
        except: pass
        try: np.save(os.path.join(Time_folder,'Stokes_Fit'), self.Exp.Fit_Stokes)
        except: pass
        try: np.save(os.path.join(Time_folder,'Deg_of_pol'), self.Degree_of_pol)
        except: pass
        try: np.save(os.path.join(Time_folder,'Deg_of_lin_pol'), self.Degree_of_lin_pol)
        except: pass
        try: np.save(os.path.join(Time_folder,'Sn'), self.Sn)
        except: pass
        try: np.save(os.path.join(Time_folder,'Sd'), self.Sd)
        except: pass
        try: np.save(os.path.join(Time_folder,'St'), self.St)
        except: pass
        try: np.save(os.path.join(Time_folder,'ROI'), np.array([self.Exp.x0, self.Exp.y0, self.Exp.width, self.Exp.height]))
        except: pass
        try: # Save the polarisation plot as a png file
            pixmap = QtGui.QPixmap.grabWidget(self.qwtPlot_pola) # create a Pixmap
            pixmap.save(os.path.join(Time_folder,'Plot_pola_' + Time + '.png'), "PNG") # Save pixmap
        except: pass
        try: np.save(os.path.join(Time_folder,'Sat_parameters'), self.par)
        except: pass
        try: np.save(os.path.join(Time_folder,'Sat_Sout'), self.Sout_Psat)
        except: pass
        try: np.save(os.path.join(Time_folder,'Sat_Sin'), self.Sin_Psat)
        except: pass
        try: # Save the Saturation plot as a png file
            pixmap = QtGui.QPixmap.grabWidget(self.qwtPlot_hcam) # create a Pixmap
            pixmap.save(os.path.join(Time_folder,'Plot_Psat_' + Time + '.png'), "PNG") # Save pixmap
        except: pass

    def Set_APD_index(self):
        self.APD_index = int(self.lineEdit_APD_index.text()) # index 0 (1) for APD 1 (2)
    def Set_CurrentCam(self):
        self.CurrentCam = self.comboBox_select_camera.currentText()
    def Set_SPD_sample_time(self):
        self.Exp.SPD_sample_time = float(self.lineEdit_SPD_sample_time.text())
    def Set_Arduino_port(self):
        self.Arduino_port = self.lineEdit_Arduino_port.text()
    def Set_hwp_start(self):
        self.hwp_start_Pos = int(self.lineEdit_hwp_start.text())
    def Set_hwp_stop(self):
        self.hwp_stop_Pos = int(self.lineEdit_hwp_stop.text()) 
    def Set_hwp_step(self):
        self.hwp_step= int(self.lineEdit_hwp_step.text())
    def Set_hwp_in_red_start(self):
        self.hwp_in_red_start_Pos = int(self.lineEdit_hwp_in_red_start.text())
    def Set_hwp_in_red_stop(self):
        self.hwp_in_red_stop_Pos = int(self.lineEdit_hwp_in_red_stop.text()) 
    def Set_hwp_in_red_step(self):
        self.hwp_in_red_step = int(self.lineEdit_hwp_in_red_step.text())
    def Set_qwp_start(self):
        self.qwp_start_Pos = int(self.lineEdit_qwp_start.text())
    def Set_qwp_stop(self):
        self.qwp_stop_Pos = int(self.lineEdit_qwp_stop.text())     
    def Set_qwp_step(self):
        self.qwp_step= int(self.lineEdit_qwp_step.text())
    def Set_hwp_out_start(self):
        self.hwp_out_start_Pos = int(self.lineEdit_hwp_out_start.text())
    def Set_hwp_out_stop(self):
        self.hwp_out_stop_Pos = int(self.lineEdit_hwp_out_stop.text()) 
    def Set_hwp_out_step(self):
        self.hwp_out_step= int(self.lineEdit_hwp_out_step.text())
    def Set_x0(self):
        self.Exp.x0 = int(self.lineEdit_x0.text())
    def Set_y0(self):
        self.Exp.y0 = int(self.lineEdit_y0.text())
    def Set_width(self):
        self.Exp.width = int(self.lineEdit_width.text())
    def Set_height(self):
        self.Exp.height = int(self.lineEdit_height.text())
    def Set_phi_0(self):
        self.phi_0 = float(self.lineEdit_phi_0.text())
    def Set_PG_x0(self):
#        self.Exp.PG_x0 = int(self.lineEdit_PG_x0.text())
        self.Exp.PGcam_roi[0] = int(self.lineEdit_PG_x0.text())
    def Set_PG_y0(self):
#        self.Exp.PG_y0 = int(self.lineEdit_PG_y0.text())
        self.Exp.PGcam_roi[1] = int(self.lineEdit_PG_y0.text())
    def Set_PG_x1(self):
#        self.Exp.PG_x1 = int(self.lineEdit_PG_x1.text())
        self.Exp.PGcam_roi[2] = int(self.lineEdit_PG_x1.text())
    def Set_PG_y1(self):
#        self.Exp.PG_y1 = int(self.lineEdit_PG_y1.text())
        self.Exp.PGcam_roi[3] = int(self.lineEdit_PG_y1.text())
    def Set_PG_V_x0(self):
#        self.Exp.PG_x0 = int(self.lineEdit_PG_x0.text())
        self.Exp.PGcam_V_roi[0] = int(self.lineEdit_PG_x0_V.text())
    def Set_PG_V_y0(self):
#        self.Exp.PG_y0 = int(self.lineEdit_PG_y0.text())
        self.Exp.PGcam_V_roi[1] = int(self.lineEdit_PG_y0_V.text())
    def Set_PG_V_x1(self):
#        self.Exp.PG_x1 = int(self.lineEdit_PG_x1.text())
        self.Exp.PGcam_V_roi[2] = int(self.lineEdit_PG_x1_V.text())
    def Set_PG_V_y1(self):
#        self.Exp.PG_y1 = int(self.lineEdit_PG_y1.text())
        self.Exp.PGcam_V_roi[3] = int(self.lineEdit_PG_y1_V.text())
    def Set_Exp_time(self):
        self.Exp.ExpTime = float(self.lineEdit_Exp_time.text())
    def Set_Exp_time_PG(self):
        self.Exp.ExpTime_PGcam = float(self.lineEdit_Exp_time_PG.text())
    def Set_Exp_time_PG_V(self):
        self.Exp.ExpTime_PGcam_V = float(self.lineEdit_Exp_time_PG_V.text())
    def Set_Save_path(self):
        self.Save_path = self.lineEdit_Save_path.text()
    def Set_Oscillo_Acquisition_time(self):
        self.Exp.Acquisition_time = int(self.lineEdit_Acquisition_time.text())
    def Set_Oscillo_channel(self):
        self.Exp.trans_chan = int(self.lineEdit_Channel.text())
    def Set_retardance(self):
        self.retardance= float(self.lineEdit_retardance.text())
    def Set_Plot_resolution(self):
        self.plot_factor = int(self.comboBox_resolution.currentText())
        
    def ShowMessage(self, Str):
        ctypes.windll.user32.MessageBoxA(0, Str, "Warning", 0)
        
    def Show_ROI_par(self):
        self.lineEdit_x0.setText(str(self.Exp.x0))
        self.lineEdit_y0.setText(str(self.Exp.y0))
        self.lineEdit_width.setText(str(self.Exp.width))
        self.lineEdit_height.setText(str(self.Exp.height))

    def Snapshot_raw(self):
        if self.CurrentCam == self.Camera_1: # if Hamamatsu camera is selected ...
            if self.Exp.hcam.hamamatsu_open: # If one camera is opened ...
                self.Exp.hcam_compatible_ROI([[self.Exp.x0, self.Exp.y0], [self.Exp.x0 + self.Exp.width, self.Exp.y0 + self.Exp.height]])
                self.Show_ROI_par()
#                self.Snapshot_pic = self.Exp.getImageHamamatsu(self.Exp.hcam)
                Snapshot_pic = self.Exp.getImageHamamatsu(self.Exp.hcam)
            else:
                self.ShowMessage('Camera not opened !') # else show warning.
                return
        elif self.CurrentCam == self.Camera_2: # if PointGrey camera is selected ...
            if self.Exp.PGcam.PointGrey_open: # If one PG camera is opened ...
#                self.Snapshot_pic = self.Exp.getImagePointGrey(self.Exp.PGcam, self.Exp.ExpTime_PGcam, self.Exp.PGcam_roi)
                Snapshot_pic = self.Exp.getImagePointGrey(self.Exp.PGcam, self.Exp.ExpTime_PGcam, self.Exp.PGcam_roi)
            else:
                self.ShowMessage('Camera not opened !') # else show warning.
                return
        elif self.CurrentCam == self.Camera_3: # if PointGrey camera is selected ...
            if self.Exp.PGcam_V.PointGrey_open: # If one PG camera is opened ...
#                self.Snapshot_pic = self.Exp.getImagePointGrey(self.Exp.PGcam_V, self.Exp.ExpTime_PGcam_V, self.Exp.PGcam_V_roi)
                Snapshot_pic = self.Exp.getImagePointGrey(self.Exp.PGcam_V, self.Exp.ExpTime_PGcam_V, self.Exp.PGcam_V_roi)
            else:
                self.ShowMessage('Camera not opened !') # else show warning.
                return
        elif self.CurrentCam == self.Camera_all: # if all camera are selected ...
            self.ShowMessage('Select a unique camera!') # Show warning.
            return
#        self.Plot_camera(self.Snapshot_pic)
        self.Plot_camera(Snapshot_pic, 'label_Plot_Camera')
    
    def Snapshot_Dark(self):
        if self.CurrentCam == self.Camera_1: # if Hamamatsu camera is selected ...
            self.Exp.hcam_compatible_ROI([[self.Exp.x0, self.Exp.y0], [self.Exp.x0 + self.Exp.width, self.Exp.y0 + self.Exp.height]])
            self.Show_ROI_par()
            Snapshot_pic = self.Exp.GetDarkSubstractedImage()
        elif self.CurrentCam == self.Camera_2: # if PointGrey camera is selected ...
            Snapshot_pic = self.Exp.GetDarkSubstractedImage_PG()
        self.Plot_camera(Snapshot_pic, 'label_Plot_Camera')
        
    def start_APD(self):
        self.acquiring_APD = True
        self.APD_plot_sample_time = float(self.from_lineEdit('lineEdit_plot_sample_time'))
        self.Exp.SPD.SetVisibleModuleDetectionMode(0,0) # Continuous mode for APD 1
        self.Exp.SPD.SetVisibleModuleDetectionMode(1,0) # Continuous mode for APD 2
        self.refresh_qwtPlot_hcam()
        self.acquisition_thread_APD.start()
        
    def Start_Csignal(self):
        
        if (self.CurrentCam == self.Camera_1 and self.Exp.hcam.hamamatsu_open) or self.CurrentCam == self.Camera_all: # if at least the Hamamatsu camera is selected ...
            self.Exp.Set_hamamatsu_properties(self.Exp.hcam)
            self.Exp.hcam.setSubArrayMode() # turn on subarray mode
            self.refresh_qwtPlot_hcam()
            self.Exp.hcam.startAcquisition()
            
        if (self.CurrentCam == self.Camera_2 and self.Exp.PGcam.PointGrey_open) or self.CurrentCam == self.Camera_all: # if at least the PointGrey camera is selected ...
            self.Exp.PGcam.camera.SetPropertyValue('shutter', self.Exp.ExpTime_PGcam) # [ms] 
#            self.refresh_qwtPlot_PGcam()
            self.Exp.PGcam.camera.StartCapture()
            
        if (self.CurrentCam == self.Camera_3 and self.Exp.PGcam_V.PointGrey_open) or self.CurrentCam == self.Camera_all: # if at least the PointGrey camera is selected ...
            self.Exp.PGcam_V.camera.SetPropertyValue('shutter', self.Exp.ExpTime_PGcam_V) # [ms] 
#            self.refresh_qwtPlot_PGcam()
            self.Exp.PGcam_V.camera.StartCapture()

        self.acquiring_Csignal = True
        self.acquisition_thread_Csignal.start()
    
    def stop_APD(self):
        self.acquiring_APD = False
        
    def Stop_Csignal(self):
        self.acquiring_Csignal = False
        time.sleep(0.2)
        if self.CurrentCam == self.Camera_1 or self.CurrentCam == self.Camera_all: # if at least Hamamatsu camera is selected ...
            self.Exp.hcam.stopAcquisition()
        if self.CurrentCam == self.Camera_2 or self.CurrentCam == self.Camera_all: # if at least PG camera camera is selected ...
            self.Exp.PGcam.camera.StopCapture()
        if self.CurrentCam == self.Camera_3 or self.CurrentCam == self.Camera_all: # if at least PG camera camera is selected ...
            self.Exp.PGcam_V.camera.StopCapture()
#        self.acquiring_Csignal = False
        
    def Stop_Polarisation_map(self):
        self.acquiring_map = False
        
    def Switch_APT_hwp(self):
        if self.checkBox_APT_hwp.isChecked():
            self.Exp.RMHWP = APTMotor(self.Exp.APT_hwp_ID, HWTYPE=31)
        elif not self.checkBox_APT_hwp.isChecked():
            self.Exp.RMHWP.cleanUpAPT()
    
    def Switch_APT_hwp_laser(self):
        if self.checkBox_APT_hwp_laser.isChecked():
            self.Exp.RMHWP_in = APTMotor(self.APT_hwp_in_ID, HWTYPE=31)
        elif not self.checkBox_APT_hwp_laser.isChecked():
            self.Exp.RMHWP_in.cleanUpAPT()
            
    def Switch_APT_hwp_out(self):
        if self.checkBox_APT_hwp_out.isChecked():
            self.Exp.RMHWP_out = APTMotor(self.Exp.APT_hwp_out_ID, HWTYPE=31)
        elif not self.checkBox_APT_hwp_out.isChecked():
            self.Exp.RMHWP_out.cleanUpAPT()
            
    def Switch_APT_qwp(self):
        if self.checkBox_APT_qwp.isChecked():
            self.Exp.RMQWP = APTMotor(self.Exp.APT_qwp_ID, HWTYPE=31)
        elif not self.checkBox_APT_qwp.isChecked():
            self.Exp.RMQWP.cleanUpAPT()
            
    def Switch_Arduino(self):
        if self.checkBox_Arduino.isChecked():
            self.Exp.arduino = Arduino(self.Arduino_port)
        else:
            self.Exp.arduino.close()
        
    def Switch_camera(self):
        if self.checkBox_Camera.isChecked():
            self.Connect_camera()
        elif not self.checkBox_Camera.isChecked():
            # Close camera connection
            self.Exp.hcam.shutdown()
            self.Exp.hcam.terminate()
            
    def Switch_camera_PG(self):
        if self.checkBox_Camera_PG.isChecked():
            self.Connect_camera()
        elif not self.checkBox_Camera.isChecked():
            # Close camera connection
            self.Exp.hcam.shutdown()
            self.Exp.hcam.terminate()
    
    def switch_SPD(self):
        if self.checkBox_TDC.isChecked():
            dev_list, dev_numb = self.Exp.SPD.ListATdevices()
            if dev_numb == 0: self.ShowMessage('No Aurea devices detected.')
            else : self.ShowMessage('Aurea devices connected : {l}'.format(l = dev_list))
            if dev_numb == 1: self.Exp.SPD.OpenATdevice(0) # Open device by default
        if (not self.checkBox_TDC.isChecked()) and self.Exp.SPD.AT_device_opened:
            self.Exp.SPD.CloseATdevice() 
    
    def Set_hcam_max(self):
        self.S_hcam_max = float(self.lineEdit_hcam_max.text())
    
    def Set_hcam_min(self):
        self.S_hcam_min = float(self.lineEdit_hcam_min.text())

    def Switch_Laser(self):
        if self.checkBox_Laser.isChecked():
            self.Exp.Laser('on')
        elif not self.checkBox_Laser.isChecked():
            self.Exp.Laser('off')
    
    def Switch_Pin_6(self):
        if self.checkBox_Arduino_Pin_6.isChecked():
            self.Exp.arduino.digital_write(self.Exp.Pin_6, 1)
        else:
            self.Exp.arduino.digital_write(self.Exp.Pin_6, 0)

    def Set_PG45_max(self):
        self.S_PG45_max = float(self.lineEdit_PG45_max.text())
    
    def Set_PG45_min(self):
        self.S_PG45_min = float(self.lineEdit_PG45_min.text())

    def Set_PGvert_max(self):
        self.S_PGvert_max = float(self.lineEdit_PGvert_max.text())
    
    def Set_PGvert_min(self):
        self.S_PGvert_min = float(self.lineEdit_PGvert_min.text())
            
    def Switch_normalisation(self):
        self.normalisation = self.checkBox_normalisation.isChecked()
       
    def Switch_Shutter(self):
        if self.checkBox_Shutter.isChecked():
            self.Exp.Shutter('close')
        else:
            self.Exp.Shutter('open')
            
    def Switch_Shutter_microscope(self):
        if self.checkBox_Shutter_microscope.isChecked():
            self.Exp.Shutter_(self.Exp.Pin_5, 'close')
        else:
            self.Exp.Shutter_(self.Exp.Pin_5, 'open')

    def Switch_Shutter_out(self):
        if self.checkBox_Shutter_out.isChecked():
            self.Exp.Shutter_out('close')
        else:
            self.Exp.Shutter_out('open')
        
if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)
    win_plot = QtGui.QMainWindow()

    main = Main()
    main.show()
    sys.exit(app.exec_())
    