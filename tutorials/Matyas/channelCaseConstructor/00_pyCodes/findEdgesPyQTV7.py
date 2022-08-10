#!/usr/bin/python

#FILE DESCRIPTION=======================================================
# Python script to automatically find edges in an input image - with GUI
#
# CHANGELOG
# 2018-11-26
# - divided SELF.COMPUTE function into several steps to get cleaner
#   code structure
# - finished the code input/output functions (these are not ideal, but
#   they work
# - implemented the hovering value display when a slider is moved
# - TO DO:
#   * try to identify which parameter was changed and update only the
#     needed data (to save the computational time)
#   * go through the "verticality" and "horizontality", clean it up and
#     correct all the bugs
#   * implement correction for rotated images - this might be tough and
#     it will be dependent on some hardcoded threshold to distinguish
#     between "horozintal" and "vertical" lines
# 2018-11-24
# - converted the code from tkinter to pyQt (it seems easier to work with)
# - first finished GUI layout
# Notes:
# - you need to create gui elements in a same structure as your data
#   in order to have a clear and easily understandable link between them
# -> rewrite the gui elements to reflect the data structure
# - due to the slider properties (only integers), I had to rewrite a few
#   things as percentual ratios
# 2018-11-20
# - implement automatic gui generation - you have identified the main
#   parameters - draw gui on paper and construct it
# - look into update of gui images after recomputing - it seems that you
#   can recompute results from gui, but the image outputs are not updated
# - would it be usefull to have all the outputs in a canvas?
# - implemented a data structure to store all the important gui parameters
#   (this class has a number of subclasses - could it be simplified?)
# 2018-11-14
# - including cropping of data at edges (to not to overflow from channel)
# - started working on GUI
# 2018-11-13
# - made the program outputs consistent with Tomas data
# 2018-11-12
# - initial version
# - I need to go through the matrix indexing and clean it up and make
#   it consistent (at the moment, it works only thanks to the fact, that
#   do domain is basically a squara in the middle of the image)
# - it would be worthwhile NOT TO FLIP the corners - while it looks nice
#   my origin of the coordinate system is in (nRows,0), which is not
#   ideal

#HELP===================================================================
#
# 1. input parameters description
# -- IO data
# inDir, outDir     ... input/output directory
# imNm              ... input image name
# outFlNm           ... output file name (loadable in Tomases functions)
#
# -- process parameters
# wallThickness     ... media wall thickness, in (m), without coating
# objectSize        ... total size of the channel,
#                     * = 2 x wallThickness + channelWidth, in (m)
#                     * !!! ALWAYS ASSUMES SQUARE CHANNELS !!!
# wallRatioInImg    ... how much of the wall is visible in the image
#                     * if wallRatioInImg = 1.0, I always assume that 
#                       the outer edges of the channel are wall without
#                       coating (this is used to move the center of the
#                       coordinate system properly)
#
# -- image thresholding and closing
# bwThresh          ... light intensity to create binary image
#                     * in [0,255] val < bwThresh -> 0, default = 127
# kernSizeR         ... ratio of kernel for closing w.r.t. image size
#                     * default is 5%
#
# -- canny edge detection and hough transform
# cannyMin, cannyMax... minimum and maximum threshold values for
#                       hysteresis in cv2 canny edge detection
#                     * default = 100, 200
# minLenR           ... ratio of the minimal line length w.r.t. image
#                       size (for Hough lines detection)
#                     * default = 5%
#
# -- custom edge detection function
# outliersR         ... ratio to identify outliers w.r.t. mean arc grad
#                       default = 110%
# edgeBufferR       ... cut-off interpolation points from channel center
#                     * the points with distance from one of the edges
#                       smaller than edgeBufferR*wallThickness will be 
#                       removed from the interpolation points file
#                       (outFileNm)
#                     * default = 10%
# ptFreq            ... save every ptFreq point
#                     * default = 5
#
# 2. program usage
# a. select folder from which the processed image will be loaded
# b. select image to be loaded
# c. select folder in which to store the data
# d. select fileName of the file in which the data will be stored
# Note: the IO of the folders/files is rather cumbersome, but in the
#       future, I would like to extend the program to work with sets
#       of images (e.g. load all the images from a given folder...)
#       and to possibly write several files into given output folder
# Note: however, the IO behavior might be simplified
# e. press the "Load image" button (or Alt+l), the computation will be
#    performed automatically
# f. adjust the image processing parameters to get the wanted result
#    (or to see the effects of changes of image processing parameters
#    on the result)
# g. export the results into given file


#LICENSE================================================================
#  findEdgesPyQT.py
#
#  Copyright 2018 Martin Isoz & Tomas Hlavaty
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

#########DO NOT EDIT####################################################
# -- gui imports
from PyQt4.QtGui import *
from PyQt4.QtCore import *

# -- plotting imports
import matplotlib
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# -- image processing imports
import cv2
import numpy as np
import math

# -- other imports
import sys

#IMPORT BLOCK-CUSTOM====================================================

#CUSTOM FUNCTIONS=======================================================

#-----------------------------------------------------------------------
# Data class
#-----------------------------------------------------------------------
class guiDataClass:
    """ Python class to store the guiData in """
  
    def __init__(self,data=None):
        
        if not data:
           data = self.setDefaultData() 
        
        # -- input/output data
        self.IO = ioDataClass(data[0])
        self.TC = tcDataClass(data[1])
        self.CH = chDataClass(data[2])
        self.EE = eeDataClass(data[3])
        self.PP = ppDataClass(data[4])
        
        self.printClassContent()
        
    # Note: the data in TC,CH and EE store values for sliders
    #       therefore, I need [value,minValue,maxValue] for each
        
        
    def setDefaultData(self):
        # -- I/O data
        inDir   = '01_testImages/'
        imNm    = 'testCF3V2.png'
        
        outDir  = inDir
        outFlNm = 'edgSAutoV1.py'
        
        # -- process parameters
        wallThickness = 216.0e-6                                    #(m)
        objectSize    = 1250e-6+2*wallThickness                     #complete object size (m)
        wallRatioInImg= 1.0                                         #how much of the wall is in the image
        
        # -- image thresholding and closing
        bwThresh      = [127,50,200]
        kernSizeR     = [5,1,10]
        
        # -- canny edge detection and houdh lines transform
        cannyMin      = [100,50,150]
        cannyMax      = [200,cannyMin[0],250]
        minLenR       = [5,1,10]
        
        # -- custom edge detection function
        outliersR     = [130,50,150]
        edgeBufferR   = [10,5,100]
        ptFreq        = [5,1,10]  
        
        # -- data gathering and guiData class initialization
        inData = []

        ioData = [inDir,outDir,imNm,outFlNm]
        inData.append(ioData)
        
        tcData =[bwThresh,kernSizeR]
        inData.append(tcData)
        
        chData = [cannyMin,cannyMax,minLenR]
        inData.append(chData)
        
        eeData = [outliersR,edgeBufferR,ptFreq]
        inData.append(eeData)
        
        ppData = [wallThickness,objectSize,wallRatioInImg]
        inData.append(ppData)
        return inData
        
    def printClassContent(self): 
        outStr = []       
        for attr in self.__dict__:
            for val in self.__dict__[attr].__dict__.items():
                outStr.append(str(val))
        return outStr
        

# -- auxiliary data classes
class ioDataClass:
    def __init__(self,data):
        self.inFolder        = data[0]
        self.outFolder       = data[1]
        self.inFile          = data[2]
        self.outFile         = data[3]
        
class tcDataClass:
    def __init__(self,data):
        self.bwThresh        = data[0]
        self.kernSizeR       = data[1]
        
class chDataClass:
    def __init__(self,data):
        self.cannyMin        = data[0]
        self.cannyMax        = data[1]
        self.minLenR         = data[2]
        
class eeDataClass:
    def __init__(self,data):
        self.outliersR       = data[0]
        self.edgeBufferR     = data[1]
        self.ptFreq          = data[2]
        
class ppDataClass:
    def __init__(self,data):
        self.wallThickness   = data[0]
        self.objectSize      = data[1]
        self.wallRatioInImg  = data[2]
        
class compDataClass:
    def __init__(self):
        # -- images
        self.initIm     = None
        self.closing    = None
        self.edges      = None
        # -- computation data
        self.edges      = None
        self.minHCor    = None
        self.maxHCor    = None
        self.minVCor    = None
        self.maxVCor    = None
        self.origLst    = None
        self.pixelSize  = None
        # -- auxiliary hardcoded
        self.lThickness = 10
        # Note: this is just to have all the variables that will be
        #       stored in the compData at one place
        
class plot2dCanvas(FigureCanvas):
    def __init__(self,canvLabel=None):
        if not canvLabel:
            canvLabel = 'unsetLabel'
        self.idLabel = canvLabel
        self.data = []
        self.fig = Figure()
        self.fig.suptitle('To be updated')
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(
            self,
            QSizePolicy.Expanding,
            QSizePolicy.Expanding,
        )
        FigureCanvas.updateGeometry(self)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel('')
        self.ax.set_ylabel('')
        self.toolbar = NavigationToolbar(self, self)
        self.toolbar.setOrientation(Qt.Vertical)
        
    def showImInPlot(self,im,imName):
        self.ax.imshow(im)
        self.fig.suptitle(imName)
        self.draw()
        
    def showFinalPlot(self,xyData,chanDims=None):
        # Inputs
        # xyData ... list of lists with data to plot
        # chanDims = [minX maxX minY maxY]
        #        ... original channel dimensions to get the scales
        if not chanDims:
            chanDims = [
                min(xyData[0]),max(xyData[0]),
                min(xyData[1]),max(xyData[1]),
            ]
        minX,maxX,minY,maxY = chanDims
        self.ax.clear()
        self.ax.plot(xyData[0],xyData[1],linewidth=3,color='red')
        self.ax.set_xlim(0,maxX)
        self.ax.set_ylim(0,maxY)
        self.ax.set_aspect('equal')
        
        rect1 = patches.Rectangle((0,0),minX,minY,edgecolor='none',facecolor='blue',alpha=0.5)
        rect2 = patches.Rectangle((minX,0),maxX,minY,edgecolor='none',facecolor='blue',alpha=0.5)
        rect3 = patches.Rectangle((0,minY),minX,maxY,edgecolor='none',facecolor='blue',alpha=0.5)
        
        self.ax.add_patch(rect1)
        self.ax.add_patch(rect2)
        self.ax.add_patch(rect3)
        
        self.fig.suptitle('Identified channel cross-section')
        self.ax.set_xlabel('x [microns]')
        self.ax.set_ylabel('y [microns]')
        self.draw()
        
class hoverDialog(QMessageBox):
    def __init__(self):
        super(hoverDialog,self).__init__()
        self.setModal(False)

        
    def showEvent(self,event):
        geom = self.frameGeometry()
        showPos = QCursor.pos()
        showPos -= QPoint(-int(self.size().width()*0.6),int(self.size().height()*0.8))
        geom.moveCenter(showPos)
        self.setGeometry(geom)
        super(hoverDialog,self).showEvent(event)
        
class initGUI(QMainWindow):
    def __init__(self):
        super(initGUI, self).__init__()
        self.centralWidget = QWidget()
        
        # -- initialize the data
        self.guiData = guiDataClass()
        self.compData= compDataClass()
        
        # -- prepare the fields to save the data in
        self.ioGUI = []
        self.tcGUI = []
        self.chGUI = []
        self.eeGUI = []
        self.ppGUI = []
        # Note: this is a mirror of guiDataClass
        
        grid = QGridLayout()
        #~ grid.setSpacing(10)
        self.centralWidget.setLayout(grid)
         
        # -- input/output data
        ioNames = [
            'Input folder:',
            'Image titles:',
            'Output folder:',
            'Output file:',
        ]
        ioButtonsNames = [
            'Set folder',
            'Set image',
            'Set folder',
            'Set file',
        ]
        ioOrder = ['inFolder','inFile','outFolder','outFile']
        
        nLines = len(ioNames)
        
        ioBoxLst = []
        ioBoxLbls = []
        ioButtonsLst = []
        for ind in range(nLines):
            ioBoxLbls.append(QLabel(ioNames[ind]))
            ioBoxLst.append(QLineEdit(self.guiData.__dict__['IO'].__dict__[ioOrder[ind]]))
            ioBoxLst[-1].setObjectName(ioOrder[ind]+'LineEdit')
            ioButtonsLst.append(QPushButton(ioButtonsNames[ind]))
            ioButtonsLst[-1].setObjectName(ioOrder[ind]+'Button')
            grid.addWidget(ioBoxLbls[-1],ind,0,1,1)
            grid.addWidget(ioBoxLst[-1],ind,1,1,4)
            grid.addWidget(ioButtonsLst[-1],ind,5,1,1)
        
            if ind % 2 != 0:
                ioButtonsLst[-1].clicked.connect(self.setFile)
            else:
                ioButtonsLst[-1].clicked.connect(self.setDir)
            # Note: this is not particularly nice... (but it seems
            #       sufficient for my cause)
            
            ioBoxLst[-1].editingFinished.connect(self.harvestData)
            
        self.ioGUI = [ioBoxLbls,ioBoxLst,ioButtonsLst,ioOrder]
        
        grid.addWidget(QHLine(),nLines,0,1,5)
        nLines += 1
        
        # -- setting placeholder for images
        nRowsFigs = 5
        nColsFigs = 5
        
        tabs        = QTabWidget()
        tabInit     = self.genTabCanvas('tabInit')
        tabClosed   = self.genTabCanvas('tabClosed')
        tabEdges    = self.genTabCanvas('tabEdges')
        tabFinal    = self.genTabCanvas('tabFinal')
        
        tabs.addTab(tabInit,'Initial ima&ge')
        tabs.addTab(tabClosed,'&Closed image')
        tabs.addTab(tabEdges,'Identified &edges')
        tabs.addTab(tabFinal,'&Result')
        tabs.setUpdatesEnabled(True)
        
        grid.addWidget(tabs,nLines,0,nRowsFigs,nColsFigs)
        
        self.tabs = tabs
        
        # -- prepare log output text field
        logOut = QPlainTextEdit()
        logOut.setReadOnly(True)
        logOut.setPlainText("Program successfully started\n")
        logOut.appendPlainText("Initial parameter values")
        for line in self.guiData.printClassContent():
            logOut.appendPlainText(line)
        logOut.appendPlainText("\n")
        
        grid.addWidget(logOut,nLines,nColsFigs,nRowsFigs,2)
        
        self.logOut = logOut
        
        nLines += nRowsFigs
        
        grid.addWidget(QHLine(),nLines,0,1,5)
        
        nLines += 1
        
        # -- setting textBoxes for process parameters
        ppNames = [
            'Wall thickness',
            'Object size',
            'Wall ratio in image',
        ]
        ppOrder = ['wallThickness','objectSize','wallRatioInImg']
        ppLst = []
        ppLbls= []
        nPPs  = len(ppNames)
        for ind in range(nPPs):
            ppLbls.append(QLabel(ppNames[ind]))
            ppLst.append(QLineEdit(repr(self.guiData.__dict__['PP'].__dict__[ppOrder[ind]])))
            grid.addWidget(ppLbls[-1],nLines,2*ind,1,1)
            grid.addWidget(ppLst[-1],nLines,2*ind+1,1,1)
            
            ppLst[-1].editingFinished.connect(self.harvestData)
            
        self.ppGUI = [ppLbls,ppLst,ppOrder]
            
        nLines += 1
        
        # Note: keep in on one line
        
        grid.addWidget(QHLine(),nLines,0,1,4)
        nLines += 1
        
        # -- setting sliders
        
        slNames = [
            'bwThresh',
            'kernSizeR',
            'cannyMin ',
            'cannyMax ',
            'minLenR',
            'outliersR',
            'edgeBufferR',
            'ptFreq',
        ]
        slLbls = []
        slLst  = []
        nSldrs = len(slNames)
        for ind in range(0,nSldrs,2):
            slLbls.append(QLabel(slNames[ind]))
            slLst.append(QSlider(Qt.Horizontal))
            grid.addWidget(slLbls[-1],nLines+ind,0,1,1)
            grid.addWidget(slLst[-1],nLines+ind,1,1,1)
            slLbls.append(QLabel(slNames[ind+1]))
            slLst.append(QSlider(Qt.Horizontal))
            grid.addWidget(slLbls[-1],nLines+ind,2,1,1)
            grid.addWidget(slLst[-1],nLines+ind,3,1,1)
            
            map(slLst[-2].sliderReleased.connect, [self.harvestData,self.closeHoverDialog])
            map(slLst[-1].sliderReleased.connect, [self.harvestData,self.closeHoverDialog])
            slLst[-2].sliderMoved.connect(self.showHoverDialog)
            slLst[-1].sliderMoved.connect(self.showHoverDialog)
            
            
        tcOrder = ['bwThresh','kernSizeR']
        chOrder = ['cannyMin','cannyMax','minLenR']
        eeOrder = ['outliersR','edgeBufferR','ptFreq']
            
        self.tcGUI = [slLbls[0:2],slLst[0:2],tcOrder]
        self.chGUI = [slLbls[2:5],slLst[2:5],chOrder]
        self.eeGUI = [slLbls[5::],slLst[5::],eeOrder]
        # Note: rather nasty
        
        self.setSliderVals(tcOrder,'tcGUI','TC')
        self.setSliderVals(chOrder,'chGUI','CH')
        self.setSliderVals(eeOrder,'eeGUI','EE')
        
        nLines += nSldrs
        
        # -- set the load image button
        loadImButton = QPushButton('&Load image')
        loadImButton.clicked.connect(self.loadImage)
        grid.addWidget(loadImButton,nLines-4,5,1,1)
        
        # -- set the export button
        exportButton = QPushButton('Export &data')
        exportButton.clicked.connect(self.exportEdgeFile)
        grid.addWidget(exportButton,nLines-2,5,1,1)
                                
        self.setCentralWidget(self.centralWidget)
        self.setWindowTitle('Testing GUI')
        self.show()
        
    def setFile(self,objectName):
        fileStr   = QFileDialog.getOpenFileName(self,'Set file')        #get the file name
        fileStr   = fileStr.split('/')[-1]                              #get only file name
        self.QFileDialogUpdateIOGuiAndData(fileStr)
        # Note: this will let you select any file type!
        
    def setDir(self):
        dirPath   = QFileDialog.getExistingDirectory(self,'Select directory') + '/'
        self.QFileDialogUpdateIOGuiAndData(dirPath)
        # Note: due to my coding style, I neeed '/' at the end of dir path
        
    def QFileDialogUpdateIOGuiAndData(self,wQStr):
        wStr      = str(wQStr)                                          #convert QStr to str
        sender    = self.sender()                                       #get the signal sender
        fieldName = str(sender.objectName().split('Button')[0])         #get the concerned field
        self.guiData.IO.__dict__[fieldName] = wStr                      #update guiData
        boxLst    = self.ioGUI[1]                                       #[label,box,button]
        for box in boxLst:
            if box.objectName().split('LineEdit')[0] == fieldName:
                box.setText(wStr)
                break
            self.logOut.appendPlainText("!!!UNABLE TO PROCESS\n")
        
    def genTabCanvas(self,canvLabel):
        tab = plot2dCanvas(canvLabel)
        return tab
        
    def setSliderVals(self,slOrder,slClass,dataClass):
        for ind in range(len(slOrder)):
            sl = self.__dict__[slClass][1][ind]
            sl.setMinimum(self.guiData.__dict__[dataClass].__dict__[slOrder[ind]][1])
            sl.setMaximum(self.guiData.__dict__[dataClass].__dict__[slOrder[ind]][2])
            sl.setValue(self.guiData.__dict__[dataClass].__dict__[slOrder[ind]][0])
            
    def harvestData(self):
        # -- get data from textFields
        
        # (a) io textFields
        self.harvestDataText('ioGUI','IO')
        
        # (b) pp textFields
        self.harvestDataText('ppGUI','PP')
            
        # -- get data from sliders
        # (a) tc sliders            
        self.harvestDataSlider('tcGUI','TC')
        
        # (b) ch sliders
        self.harvestDataSlider('chGUI','CH')
        
        # (c) ee sliders
        self.harvestDataSlider('eeGUI','EE')
        
        # -- print the new data into log
        self.logOut.appendPlainText('Updated parameter values')
        for line in self.guiData.printClassContent():
            self.logOut.appendPlainText(line)
        self.logOut.appendPlainText('\n')
        
        if isinstance(self.compData.initIm,np.ndarray):
            self.compute()
        
    def harvestDataText(self,leGroup,dataClass):
        leOrder = self.__dict__[leGroup][-1]
        le      = self.__dict__[leGroup][1]
        for ind in range(len(le)):
            auxText = le[ind].text()
            try:
                self.guiData.__dict__[dataClass].__dict__[leOrder[ind]] = float(auxText)
            except:
                self.guiData.__dict__[dataClass].__dict__[leOrder[ind]] = str(auxText)
        # Note: I can have floats in my textFields
        
    def harvestDataSlider(self,slGroup,dataClass):
        slOrder = self.__dict__[slGroup][-1]
        sl      = self.__dict__[slGroup][1]
        for ind in range(len(sl)):
            self.guiData.__dict__[dataClass].__dict__[slOrder[ind]][0] = sl[ind].value()
    
    def showHoverDialog(self):
        sender     = self.sender()
        cValue     = sender.value()
        minValue   = sender.minimum()
        maxValue   = sender.maximum()
        self.hDiag = hoverDialog()
        self.hDiag.setText('min: %d < %d < %d :max'%(minValue,cValue,maxValue))
        self.hDiag.setWindowFlags(Qt.FramelessWindowHint)
        self.hDiag.setStandardButtons(QMessageBox.Ok)
        self.hDiag.removeButton(self.hDiag.button(QMessageBox.Ok))
        self.hDiag.show()
        self.hDiag.raise_()
        # Note: removing the OK button is kind of a hack. Ok button is
        #       created by default, but it does not have a pointer???
        #       so I create it manually to have a pointer and then
        #       delete it.
        
    def closeHoverDialog(self):
        self.hDiag.accept()
    
    def loadImage(self):
        inFolder = self.guiData.IO.inFolder
        inFile   = self.guiData.IO.inFile
        self.compData.initIm = cv2.imread(inFolder+inFile,0)
        
        # -- get and display the input image
        initIm = self.compData.initIm
        self.showImage('tabInit',initIm,'Initial image')
        
        self.compute()
        
    def compute(self):
        """ this function is the shit"""
                
        # -- get and display the closed image
        self.computeTC()
        
        # -- get and display the found lines
        self.computeCH()
        
        # -- find the edges and draw the result
        self.computeEE()
        
        
        self.logOut.appendPlainText("Finished\n-----------\n")
        
    def computeTC(self):
        """ apply thresholding and closing of the image """
        
        # -- get the input image (I do not need to redraw it)
        img = self.compData.initIm
        
        # -- get the paramaters
        bwThresh    = self.guiData.TC.bwThresh[0]
        kernSizeR   = self.guiData.TC.kernSizeR[0]*0.01
        
        # -- get and display the closed image
        # -- threshold the image -- black and white conversion
        ret,img = cv2.threshold(img,bwThresh,255,cv2.THRESH_BINARY)
        
        # -- close the holes in the image
        kSize  = int(round(np.average([img.shape])*kernSizeR))
        kernel = np.ones((kSize,kSize),np.uint8)
        self.compData.closing= cv2.morphologyEx(img, cv2.MORPH_CLOSE, kernel)
        self.showImage('tabClosed',self.compData.closing,'Closed image')
            
    def computeCH(self):
        """ perform canny edge detection and hough transform """
        
        # -- get the parameters
        wallThickness   = self.guiData.PP.wallThickness
        objectSize      = self.guiData.PP.objectSize
        wallRatioInImg  = self.guiData.PP.wallRatioInImg
        
        cannyMin        = self.guiData.CH.cannyMin[0]
        cannyMax        = self.guiData.CH.cannyMax[0]
        minLenR         = self.guiData.CH.minLenR[0]*0.01
        
        lThickness      = self.compData.lThickness
        
        # -- find edges via canny edge detection
        edges  = cv2.Canny(self.compData.closing,cannyMin,cannyMax)
        # Note: what do these parameters do?
        
        # -- look for the lines in the image
        minLLen= int(round(np.average([self.compData.closing.shape])*minLenR))                       #5% of average image resolution
        lines  = cv2.HoughLinesP(edges, 1, np.pi / 180, minLLen, None, 50, 10)
        self.compData.lines = lines
        
        nLines = (len(lines[0]))
        
        self.logOut.appendPlainText("Found %d lines\n"%nLines)
        
        # -- find image edges (outer lines)
        # a) sort the lines in group (horizontal/vertical)
        horVec = np.array([1.0,0.0])
        verVec = np.array([0.0,1.0])
        threshC= 0.05                                                           #colinearity threshold
        
        verLines = []
        horLines = []
        dLinesC  = 0
        for lInd in range(0, nLines):
            cLine  = lines[0][lInd]
            curVec = np.array([float(cLine[2]-cLine[0]),float(cLine[3]-cLine[1])])#direction vector of current line
            curVec/= np.linalg.norm(curVec)
            if abs(horVec.dot(curVec)) > 1.0-threshC:
                horLines.append(cLine)
            elif abs(verVec.dot(curVec)) > 1.0-threshC:
                verLines.append(cLine)
            else:
                self.logOut.appendPlainText("Cannot decide on line verticality - discarding")
                dLinesC += 1
                
        # -- test output
        self.logOut.appendPlainText("\n")
        self.logOut.appendPlainText("Found     %03d horizontal lines"%len(horLines))
        self.logOut.appendPlainText("Found     %03d vertical   lines"%len(verLines))
        self.logOut.appendPlainText("Discarded %03d            lines"%dLinesC)
        self.logOut.appendPlainText("\n")
                
        # b) group the lines with similar defining coordinate
        thresD  = int(round(np.average([self.compData.closing.shape])*0.01))                      #1% of average image resolution
        self.logOut.appendPlainText("Grouping lines with defining coordinate within %03d pixels"%thresD)
        meanHCor= [int(round(np.average([cLine[1],cLine[3]]))) for cLine in horLines]
        meanVCor= [int(round(np.average([cLine[0],cLine[2]]))) for cLine in verLines]
        
        minHCor,maxHCor = min(meanHCor),max(meanHCor)
        minVCor,maxVCor = min(meanVCor),max(meanVCor)
        
        # ## intermezzo ## -- get the wall width in pixels and so on
        objectSizePixels= np.average([maxHCor-minHCor,maxVCor-minVCor])         #average object size in pixels
        pixelSize       = objectSize/objectSizePixels
        
        wallThicknessPixels = int(round(wallThickness/pixelSize))               #wall size in pixels
        
        moveByPixels    = wallThicknessPixels*(wallRatioInImg-0.5)            #where is the wall center
        # ## \intermezzo##
        
        topHorLst = []
        botHorLst = []
        lefVerLst = []
        rigVerLst = []
        if lines is not None:
            for lInd in range(len(verLines)):
                if abs(meanVCor[lInd] - minVCor) < thresD:
                    lefVerLst.append(verLines[lInd])
                if abs(meanVCor[lInd] - maxVCor) < thresD:
                    rigVerLst.append(verLines[lInd])
            for lInd in range(len(horLines)):
                if abs(meanHCor[lInd] - minHCor) < thresD:
                    botHorLst.append(horLines[lInd])
                if abs(meanHCor[lInd] - maxHCor) < thresD:
                    topHorLst.append(horLines[lInd])
                    
        self.logOut.appendPlainText("\n")
        self.logOut.appendPlainText("Left   edge consits of %03d lines"%len(lefVerLst))
        self.logOut.appendPlainText("Right  edge consits of %03d lines"%len(rigVerLst))
        self.logOut.appendPlainText("Top    edge consits of %03d lines"%len(botHorLst))
        self.logOut.appendPlainText("Bottom edge consits of %03d lines"%len(topHorLst))
        
        # c) create average edge lines
        leftEdge = [
            int(round(np.average([cLine[0]+moveByPixels for cLine in lefVerLst]))),
            minVCor,
            int(round(np.average([cLine[2]+moveByPixels for cLine in lefVerLst]))),
            maxHCor
        ]
        rightEdge = [
            int(round(np.average([cLine[0]-moveByPixels for cLine in rigVerLst]))),
            minVCor,
            int(round(np.average([cLine[2]-moveByPixels for cLine in rigVerLst]))),
            maxHCor
        ]
        
        botEdge = [
            minHCor,
            int(round(np.average([cLine[1]+moveByPixels for cLine in botHorLst]))),
            maxHCor,
            int(round(np.average([cLine[3]+moveByPixels for cLine in botHorLst]))),
        ]
        topEdge = [
            minHCor,
            int(round(np.average([cLine[1]-moveByPixels for cLine in topHorLst]))),
            maxHCor,
            int(round(np.average([cLine[3]-moveByPixels for cLine in topHorLst]))),
        ]
        
        # d) locate the "centers of coordinate system" in each corner
        # Note: and clean up the confusion with top-bottom and left-right
        topLeft = [
            int(round(np.average([cLine[0]+moveByPixels for cLine in lefVerLst]))),
            int(round(np.average([cLine[1]+moveByPixels for cLine in botHorLst])))
        ]
        
        topRight = [
            int(round(np.average([cLine[0]-moveByPixels for cLine in rigVerLst]))),
            int(round(np.average([cLine[1]+moveByPixels for cLine in botHorLst])))
        ]
        
        botLeft = [
            int(round(np.average([cLine[0]+moveByPixels for cLine in lefVerLst]))),
            int(round(np.average([cLine[1]-moveByPixels for cLine in topHorLst])))
        ]
        
        botRight = [
            int(round(np.average([cLine[0]-moveByPixels for cLine in rigVerLst]))),
            int(round(np.average([cLine[1]-moveByPixels for cLine in topHorLst])))
        ]
        
        
        # e) prepare the image corner coordinate "after division and fliping"
        nRows,nCols = self.compData.closing.shape                                                 #matrix dimensions
        splitH = int(np.average([minVCor,maxVCor]))
        splitV = int(np.average([minHCor,maxHCor]))
        
        topLeft = [nRows-splitV-topLeft[1],topLeft[0]]                          #!!! row vs column
        topRight= [nRows-splitV-topRight[1],nCols-topRight[0]]
        botLeft = [botLeft[1]-splitV,botLeft[0]]
        botRight= [botRight[1]-splitV,nCols-botRight[0]] 
        
        origLst = [botLeft,botRight,topLeft,topRight]                   #list of coordinate systems origins
        # Note: must be in the same order as it is in the splitAndFlip function
        #       (to be moved there)
        
        # z) test image output
        drawImageWithLines(edges,[lefVerLst,rigVerLst,botHorLst,topHorLst],lThickness,False)
        drawImageWithLines(edges,[[leftEdge],[rightEdge],[botEdge],[topEdge]],int(lThickness*0.5),False)
        
        self.showImage('tabEdges',edges,'Identified edges')
        
        # -- save the important data into compData class
        
        self.compData.edges     = edges                                 #image with identified edges drawn in
        self.compData.origLst   = origLst                               #list of coordinate systems origins
        self.compData.minHCor   = minHCor                               #coordinate edges
        self.compData.maxHCor   = maxHCor
        self.compData.minVCor   = minVCor
        self.compData.maxVCor   = maxVCor
        self.compData.pixelSize = pixelSize                             #object size in pixels
        
    
    def computeEE(self):
        """ do edges position estimation """
        
        # -- get the parameters (user defined)
        wallThickness   = self.guiData.PP.wallThickness
        objectSize      = self.guiData.PP.objectSize
        wallRatioInImg  = self.guiData.PP.wallRatioInImg
        
        outliersR       = self.guiData.EE.outliersR[0]*0.01
        edgeBuffer      = self.guiData.EE.edgeBufferR[0]*0.01*wallThickness
        ptFreq          = self.guiData.EE.ptFreq[0]
        
        # -- get the parameters (computed during program run)
        edges           = self.compData.edges
        minHCor         = self.compData.minHCor
        maxHCor         = self.compData.maxHCor
        minVCor         = self.compData.minVCor
        maxVCor         = self.compData.maxVCor
        origLst         = self.compData.origLst
        pixelSize       = self.compData.pixelSize
        
        # -- get the parameters (auxiliary line thickness)
        lThickness      = self.compData.lThickness
        
        # -- split the image to four parts (I want to get average wall position)
        # Note: in the future, if I will process more images at once, I can get
        #       average positions over several cuts. Or I can save several
        #       channel cross-sections and create curved edges ALONG the channel
        
        splitImLst = splitAndFlip(edges,[minHCor,maxHCor],[minVCor,maxVCor])
        
        [botLeftCorner,botRightCorner,topLeftCorner,topRightCorner] = splitImLst
        
        # -- crop the image and get the edge coordinates
        edgePos = []
        minLen  = 1e5                                                   #ugly trick
        for iInd in range(len(splitImLst)):
            cOrig = origLst[iInd]
            splitImLst[iInd] = splitImLst[iInd][:cOrig[0],cOrig[1]:]
            nRows,nCols = splitImLst[iInd].shape
            cEdgePos = []
            for iRow in range(nRows-lThickness*2):
                for jCol in range(lThickness*2,nCols):
                    #~ print splitImLst[iInd][iRow,jCol]
                    if splitImLst[iInd][iRow,jCol] != 0:
                        cEdgePos.append([jCol*pixelSize,(nRows-iRow)*pixelSize])
            edgePos.append(cEdgePos)
            minLen = min(minLen,len(cEdgePos))
            # Note: skip image edges - I do not want the coordinate system lines
        
        # -- get mean edge position for all the four corner
        # Note: this is not at all exact
        meanEdgePos = []
        arcLenLst = []
        for lInd in range(minLen):
            meanEdgePos.append([
                -np.average([edgePos[j][lInd][0] for j in range(len(edgePos[0][lInd]))]),
                np.average([edgePos[j][lInd][1] for j in range(len(edgePos[0][lInd]))]),
                #~ -edgePos[0][lInd][0],
                #~ edgePos[0][lInd][1],
                0.0
            ])
            arcLenLst.append(np.linalg.norm(np.array(meanEdgePos[lInd])))
        meanEdgePos = np.array(meanEdgePos)
        arcLenLst = np.array(arcLenLst)
        # Note: arcLenLst variable is used to identify "bad data" - I do not 
        #       want any jumps in my data
        # Note: the minus in x-coordinate is added to make the program outputs
        #       consistent with Tomas format of edgS.py
        
        for i in range(500):
            mask =  findOutliers(arcLenLst,outliersR)                   #identify outliers
            arcLenLst = arcLenLst[mask]                                 #reduce the data
            meanEdgePos = meanEdgePos[mask,:]
        
        edgeBufferPixels = edgeBuffer/pixelSize                         #edgeBuffer size in pixels
        
        chanDims = [
            min(-meanEdgePos[:,0]),max(-meanEdgePos[:,0]),              #min/max X !!!SIGNS!!!
            min(meanEdgePos[:,1]),max(meanEdgePos[:,1])                 #min/max Y
        ]
        chanDims = [val*1e6 for val in chanDims]                        #convert to microns
        
        arcLenLst   = arcLenLst[edgeBufferPixels:-edgeBufferPixels+1:ptFreq]
        meanEdgePos = meanEdgePos[edgeBufferPixels:-edgeBufferPixels+1:ptFreq,:]
        # Note: the results visualization is not prepared for inclusion of
        #       edgeBuffers
        
        self.logOut.appendPlainText("\n")
        self.logOut.appendPlainText("Found %04d usable points for edge interpolation"%len(mask))
        self.logOut.appendPlainText("Using %04d of them (every %02d-th point)"%(len(arcLenLst),ptFreq))
        self.logOut.appendPlainText("\n")
        
        # -- get the position of arc central point and top and bottom arcs
        indVerPoint2 = arcLenLst.argmin()
        verPoint2    = meanEdgePos[indVerPoint2,:]
        inPB01       = meanEdgePos[indVerPoint2+1:,:]
        inPB12       = meanEdgePos[:indVerPoint2,:]
        # Note: the naming is consistent with Tomas scripts
        
        self.compData.verPoint2 = verPoint2
        self.compData.inPB01    = inPB01
        self.compData.inPB12    = inPB12
            
        [botLeftCorner,botRightCorner,topLeftCorner,topRightCorner] = splitImLst
        
        
        plotPos = [[],[]]
        for lInd in range(len(meanEdgePos)):
            plotPos[0].append(-meanEdgePos[lInd][0]*1e6)
            plotPos[1].append(meanEdgePos[lInd][1]*1e6)
            
        self.showFinalPlot('tabFinal',plotPos,chanDims)
        
    def showImage(self,tabId,im,imName):
        for link in self.tabs.children()[0].children():
            if isinstance(link,plot2dCanvas):
                if link.idLabel == tabId:
                    link.showImInPlot(im,imName)
        # Note: this is probably shit programming
        
    def showFinalPlot(self,tabId,xydata,chanDims=None):
        for link in self.tabs.children()[0].children():
            if isinstance(link,plot2dCanvas):
                if link.idLabel == tabId:
                    link.showFinalPlot(xydata,chanDims)
        # Note: this is not at all generic
                    
    def exportEdgeFile(self):
        try:
            # -- write out the results
            outFolder = self.guiData.IO.outFolder
            outFile   = self.guiData.IO.outFile
            outFile   = open(outFolder+outFile,'w')
            
            verPoint2 = self.compData.verPoint2
            inPB01    = self.compData.inPB01
            inPB12    = self.compData.inPB12
            
            outFile.write('verPoint2=[%10.5e,%10.5e,%10.5e]\n'%tuple(verPoint2))
            
            outFile.write('inPB01 = [\n')
            for pos in np.flipud(inPB01):
                outFile.write('\t[%10.5e,%10.5e,%10.5e],\n'%tuple(pos))
            outFile.write(']\n')
            
            outFile.write('inPB12 = [\n')
            for pos in np.flipud(inPB12):
                outFile.write('\t[%10.5e,%10.5e,%10.5e],\n'%tuple(pos))
            outFile.write(']\n')
            
            outFile.close()
            
            self.logOut.appendPlainText("\n")
            self.logOut.appendPlainText("Data export succesfull")
            self.logOut.appendPlainText("output file: %s"%self.guiData.IO.outFile)
            self.logOut.appendPlainText("folder     : %s"%outFolder)
            self.logOut.appendPlainText("\n")
            # Note: the upside-down flips are to be consistent with Tomas
        except:
            noDataErr = QMessageBox()
            noDataErr.setIcon(QMessageBox.Critical)
            noDataErr.setWindowTitle('Error - no data to export')
            noDataErr.setText('There is no data to export yet')
            noDataErr.setInformativeText('First, load an image and compute the output')
            noDataErr.setStandardButtons(QMessageBox.Ok)
            noDataErr.exec_()
        
class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)
        
class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)

#-----------------------------------------------------------------------
# Processing functions
#-----------------------------------------------------------------------
def drawLineLst(bGround,lineLst,lThickness,lClr):
    """ draws a list of lines into background image """
    
    # Note: lClr is a single scalar that works in a mysterious way
    
    for cLine in lineLst:
        iPt   = (cLine[0],cLine[1])
        ePt   = (cLine[2],cLine[3])
        cv2.line(bGround, iPt, ePt, (lClr,0,0), lThickness,4,0)
        
    return bGround
    
def drawImageWithLines(bGround,lineLstLst,lThickness,showPlot):
    """ draws given lists of lines over image and shows the plot """
    
    # -- auxiliary calculations
    nGroups = len(lineLstLst)                                           #number of different line groups
    minClr,maxClr = 50,255
    clrStep = int(math.floor((maxClr-minClr)/nGroups))
    
    # -- cycle over the lineLsts
    cClr    = minClr
    for lineLst in lineLstLst:
        drawLineLst(bGround,lineLst,lThickness,cClr)
        cClr+=clrStep
    
    if showPlot:    
        # -- show the plots
        plt.imshow(bGround)
        plt.title('Image with edges'), plt.xticks([]), plt.yticks([])
        
        plt.show()
        
def splitAndFlip(img,horLimits,verLimits):
    """ splits an image in the middle of given limits and flips
        sub-images to be "bottom-left" """
        
    splitH = int(np.average(verLimits))
    splitV = int(np.average(horLimits))
    
    # Note: this is due to the consistency with the rest of the code
    #       vertical limits were obtained for vertical lines. i.e., they
    #       are basically horizontal limits
    
    topLeftCorner = img[:splitH,:splitV]
    botLeftCorner = img[splitH:,:splitV]
    
    topRightCorner = img[:splitH,splitV:]
    botRightCorner = img[splitH:,splitV:]
    
    # b) rotate everything to be "bottom left"
    topLeftCorner  = np.flipud(topLeftCorner)
    botRightCorner = np.fliplr(botRightCorner)
    topRightCorner = np.flipud(np.fliplr(topRightCorner))
    
    return [botLeftCorner,botRightCorner,topLeftCorner,topRightCorner]
    
def dirtySmooth(y, boxPts):
    """ a quick and dirty way to smooth data based on a moving average
        box """
        
    box     = np.ones(boxPts)/boxPts
    ySmooth = np.convolve(y, box, mode='same')
    
    return ySmooth
    
def findOutliers(y,thresO):
    """ a quick and dirty way to find outliers in data """
    
    y = np.array(y)
    
    lefDer = np.r_[0, np.abs(y[1:] - y[:-1])]
    rigDer = np.r_[np.abs(y[1:] - y[:-1]), 0]
    
    thresO *= np.average([np.average(vec) for vec in [lefDer,rigDer]])
        
    mask = (lefDer < thresO) | (rigDer < thresO)                        #data to be kept
    
    return mask
    
#-----------------------------------------------------------------------
# GUI functions
#-----------------------------------------------------------------------
    
    
def setGuiVals(guiData):
    """ function to set values in gui based on the guiData structure"""
    
    # Note: to be implemented - is it needed?
    
#-----------------------------------------------------------------------
# Actual computation function
#-----------------------------------------------------------------------
if __name__ == '__main__':
    
    winSizeX,winSizeY = 1200,800
    
    # -- put in gui elements
    # a) create fields to specify folders and fileNames
    formFields = 'Source folder', 'Image', 'Out folder', 'Edge file'
    
    slNames = 'kernel size ratio (%)','Hough - minLen ratio (%)'
    
    app = QApplication([])
    window = initGUI()
    
    window.resize(winSizeX,winSizeY)
        
    window.show()
    app.exit(app.exec_())

#=======================================================================
#							EDITABLE
#=======================================================================
