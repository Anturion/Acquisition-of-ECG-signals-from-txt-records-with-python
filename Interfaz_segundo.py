# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:11:08 2019

@author: Alejandro
"""

# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import (QMainWindow, QApplication, QLabel, QLineEdit, QDialog,
       QPushButton, QSizePolicy,QVBoxLayout, QWidget)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import os
global usuario, clave
import Filtros_EEG as filtros

#Se crea la funcion principal donde se ingresa el nombre del archiuvo y el canal
class MainWindow(QMainWindow):
    
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("nombre del archivo de señales EEG")
        self.setFixedSize(400, 400)
        self.Label = QLabel(self)
        self.Label.setText('Ingrese nombre del archivo de señales EEG')
        self.Label.resize(700,20)
        self.Label.move(100, 50)
        self.archivo=QLineEdit(self)
        self.archivo.resize(250,20)
        self.archivo.move(100,100)
        self.Label1 = QLabel(self)
        self.Label1.setText('Ingrese Canal')
        self.Label1.move(100, 150)
        self.canal =QLineEdit(self)
        self.canal.move(100, 200)
        self.button = QPushButton('Aceptar', self)
        self.button.move(150,250)
        self.button.clicked.connect(self.ingresar)
        
   
#Se crea un archivo csv con el nombre de la señal y la clave
    def ingresar(self):
        self.archivo1=self.archivo.text()
        self.canal1=self.canal.text()
        self.nuevoregistro=pd.DataFrame({"senal":[self.archivo1],"canal":[self.canal1]})
        self.nuevoregistro.to_csv("Senales.csv")
        funcion=Procesamiento(self)  # Dirige a la clase Procesamiento
        funcion.show() 
        

    def onclick_ECG2(self):
        ECG=PlotECG(self)    # Dirige a la clase PlotECG
        ECG.show()           # Muestra la grafica de la señal filtrada sin épocas
    
    def onclick_EEG2(self):
        EEG=PlotEEG(self)    # Dirige a la clase PlotEEG
        EEG.show()        # Muestra la grafica del periodograma de la señal original
    
    def onclick_EMG2(self):
        EMG=PlotEMG(self)    # Dirige a la clase PlotEMG
        EMG.show()          # Muestra la gráfica del hisograma de la señal filtrada con épocas
        
    def onclick_senal_respiratoria2(self):
        senal=Plotsenal_respiratoria(self) # Dirige a la clase Plotsenal_respiratoria
        senal.show()            # Muestra la grafica de la señal sin épocas filtradas
 
    def onclick_grafica_periodograma_resultante(self):
        senal=Plot_periodograma_resultante(self) # Dirige a la clase Plot_periodograma_resultante
        senal.show() #MUestra la gráfica del periodograma de la señal resultante

#Se define la ventana en donde estarán los botones para graficar los elementos dados
class Procesamiento(QMainWindow): 
    
     def __init__(self, *args, **kwargs):
        super(Procesamiento, self).__init__(*args, **kwargs)
        self.setWindowTitle("Graficos")
        self.setFixedSize(680,300 )
        self.button3 = QPushButton('Graficar señal filtrada mediante umbrales', self)
        self.button3.move(200,10)
        self.button3.resize(250,30)
        self.button4 = QPushButton('Graficar periodograma señal original', self)
        self.button4.move(200,40)
        self.button4.resize(250,30)
        self.button5 = QPushButton('Graficar histograma de señal con filtros', self)
        self.button5.move(200,70)
        self.button5.resize(250,30)
        self.button6 = QPushButton('Eliminación ruido con pasa altas y bajas', self)
        self.button6.move(200,100)
        self.button6.resize(250,30)
        self.button7 = QPushButton('Graficar periodograma señal filtrada', self)
        self.button7.move(200,130)
        self.button7.resize(250,30)
        self.button3.clicked.connect(self.Graficar_canal_filtrado)
        self.button4.clicked.connect(self.Graficar_periodograma)
        self.button5.clicked.connect(self.Graficar_histograma)
        self.button6.clicked.connect(self.Graficar_canal)
        self.button7.clicked.connect(self.Graficar_periodograma_resultante)
#Funcion que retrona a la clase principal donde despues se va hacia la funcion para graficar la señal filtrada sin épocas con ruido   
     def Graficar_canal_filtrado(self):
        MainWindow.onclick_ECG2(self)
#Funcion que retrona a la clase principal donde despues se va hacia la funcion para graficar el periodograma de la señal original
     def Graficar_periodograma(self):
         MainWindow.onclick_EEG2(self)
#Funcion que retrona a la clase principal donde despues se va hacia la funcion para graficar el histograma de la señal filtrada conépocas
     def Graficar_histograma(self):
         MainWindow.onclick_EMG2(self)
#Funcion que retrona a la clase principal donde despues se va hacia la funcion para graficar la señal filtrada con épocas
     def Graficar_canal(self):
         MainWindow.onclick_senal_respiratoria2(self)
#Funcion que retrona a la clase principal donde despues se va hacia la funcion para graficar el periodograma de la señal filtrada sin épocas
     def Graficar_periodograma_resultante(self):
         MainWindow.onclick_grafica_periodograma_resultante(self)
    
        
        
class PlotECG(QMainWindow):
    
    
     def __init__(self,  parent=None):
        super(PlotECG, self).__init__(parent)
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 =filtros.matrices_archivos(canal,archivo)
        self.setWindowTitle("Gráfica filtrada")
        self.setFixedSize(800,600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = Figure((140.0, 6.0), dpi=70, facecolor="#F6F4F2")
        self.Label=QLabel(self)
        self.Label.setText('Umbral Superior:  '+str(cota_superior)+" mV")
        self.Label.move(50,450)
        self.Label.resize(250,30)
        self.Label1=QLabel(self)
        self.Label1.setText('Umbral Inferior:  '+str(cota_inferior)+" mV")
        self.Label1.move(500,450)
        self.Label1.resize(250,30)
        
        
        m = PlotCanvasECG(self)
        m.move(0,0)

#Se configura la ventana donde se va a graficar la señal EEG
class PlotEEG(QMainWindow):
    
    
     def __init__(self,  parent=None):
        super(PlotEEG, self).__init__(parent)
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 =filtros.matrices_archivos(canal,archivo)
        self.setWindowTitle("Periodograma señal original")
        self.setFixedSize(800,600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = Figure((130.0, 6.0), dpi=70, facecolor="#F6F4F2")
        self.Label=QLabel(self)
        self.Label.setText('Frecuencia de corte filtro pasa-bajas:  '+str(f_pasa_bajas)+"Hz")
        self.Label.move(50,450)
        self.Label.resize(250,30)
        self.Label1=QLabel(self)
        self.Label1.setText('Frecuencia de corte filtro pasa-altas:  '+str(f_pasa_altas)+"Hz")
        self.Label1.move(500,450)
        self.Label1.resize(250,30)
        
        m = PlotCanvasEEG(self)
        m.move(0,0)

#Se configura la ventana donde se va a graficar la señal EMG        
class PlotEMG(QMainWindow):
    
    def __init__(self,  parent=None):
        super(PlotEMG, self).__init__(parent)
        self.setWindowTitle("Histograma")
        self.setFixedSize(800,600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = Figure((130.0, 6.0), dpi=70, facecolor="#F6F4F2")
        
        m = PlotCanvasEMG(self)
        m.move(0,0)
 
#Se configura la ventana donde se va a graficar la señal Señal_respiratoria     
class Plotsenal_respiratoria(QMainWindow):
    
     def __init__(self,  parent=None):
        super(Plotsenal_respiratoria, self).__init__(parent)
        self.setWindowTitle("Eliminación ruido ambiental")
        self.setFixedSize(800,600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = Figure((130.0, 6.0), dpi=70, facecolor="#F6F4F2")
        
        m = PlotCanvasSenal_respiratoria(self)
        m.move(0,0)
 
class Plot_periodograma_resultante(QMainWindow):
    
    
     def __init__(self,  parent=None):
        super(Plot_periodograma_resultante, self).__init__(parent)
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 =filtros.matrices_archivos(canal,archivo)
        self.setWindowTitle("Periodograma señal filtrada")
        self.setFixedSize(800,600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = Figure((140.0, 6.0), dpi=70, facecolor="#F6F4F2")
        
        
        m = PlotCanvasECG_periodograma_resultante(self)
        m.move(0,0)

    
#Esta clase realiza el plot sobre la ventana ya configurada      
class PlotCanvasECG(FigureCanvas):
 
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()
    
#Se define la función que genera la grafica
    def plot(self):
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 =filtros.matrices_archivos(canal,archivo)
        ax = self.figure.add_subplot(111)
        ax.plot(tiempo,senal_resultante)
        ax.set_xlabel('Segundos')
        ax.set_ylabel('mV')
        ax.grid()
        ax.set_title("Señal:"+ archivo + "  Canal:" +str(canal))
        self.draw()   #Se muestra el ECG sobre la ventana
 
#Esta clase realiza el plot sobre la ventana ya configurada         
class PlotCanvasEEG(FigureCanvas):
 
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()
    
#Se define la función que genera la grafica
    def plot(self):
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1=filtros.matrices_archivos(canal,archivo)
        ax = self.figure.add_subplot(111)
        ax.plot(f,pot)
        ax.set_xlabel('Frecuencia')
        ax.set_ylabel('Potencia')
        ax.grid()
        ax.set_title("Señal:"+ archivo + "  Canal:" +str(canal))
        self.draw()   #Se muestra el ECG sobre la ventana

#Esta clase realiza el plot sobre la ventana ya configurada 
class PlotCanvasEMG(FigureCanvas):
 
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()

#Se define la función que genera la grafica
    def plot(self):
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo,tiempo1=filtros.matrices_archivos(canal,archivo)
        ax = self.figure.add_subplot(111)
        ax.hist(senal_resultante)
        ax.grid()
        ax.set_title("Señal:"+ archivo + "  Canal:" +str(canal))
        self.draw()   #Se muestra el ECG sobre la ventanatana
        

#Esta clase realiza el plot sobre la ventana ya configurada 
class PlotCanvasSenal_respiratoria(FigureCanvas):
    
     def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()

#Se define la función que genera la grafica final
     def plot(self):
          
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo,tiempo1=filtros.matrices_archivos(canal,archivo)
        ax = self.figure.add_subplot(111)
        ax.plot(tiempo1,senal_filtrada_pasabajas1)
        ax.set_xlabel('Segundos')
        ax.set_ylabel('mV')
        ax.grid()
        ax.set_title("Señal:"+ archivo + "  Canal:" +str(canal))
        self.draw()   #Se muestra la grafica sobre la ventana

class PlotCanvasECG_periodograma_resultante(FigureCanvas):
 
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()
    
#Se define la función que genera la grafica final
    def plot(self):
        datos=pd.read_csv('Senales.csv')
        canal=int(datos['canal'][0])
        archivo=datos['senal'][0]
        senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 =filtros.matrices_archivos(canal,archivo)
        ax = self.figure.add_subplot(111)
        ax.plot(f2,pot2)
        ax.set_xlabel('Hz')
        ax.set_ylabel('mV')
        ax.grid()
        ax.set_title("Señal:"+ archivo + "  Canal:" +str(canal))
        self.draw()   #Se muestra sobre la ventana
        
if __name__ == "__main__":  
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()