# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 08:47:39 2019

@author: Alejandro
"""
#Se importan las librerias necesarias
import numpy as np
import scipy.signal as signal

#Función como parte de los filtros a aplicar
def fkernel(m, f, w):
    m = np.arange(-m/2, (m/2)+1)
    b = np.zeros((m.shape[0]))
    b[m==0] = 2*np.pi*f # No division by zero
    b[m!=0] = np.sin(2*np.pi*f*m[m!=0]) / m[m!=0] # Sinc
    b = b * w # Windowing
    b = b / np.sum(b) # Normalization to unity gain at DC
    return b

#Función como parte de los filtros a aplicar dada por el docente
def firws(m, f , w , t = None):
    """
    Designs windowed sinc type I linear phase FIR filter.
    Parameters:        
        m: filter order.
        f: cutoff frequency/ies (-6 dB;pi rad / sample).
        w: vector of length m + 1 defining window. 
        t: 'high' for highpass, 'stop' for bandstop filter. {default low-/bandpass}
    Returns:
        b: numpy.ndarray
            filter coefficients 
    """
    f = np.squeeze(f)
    f = f / 2; 
    w = np.squeeze(w)
    if (f.ndim == 0): #low pass
        b = fkernel(m, f, w)
    else:
        b = fkernel(m, f[0], w) #band
    
    if (f.ndim == 0) and (t == 'high'):
        b = fspecinv(b)
    elif (f.size == 2):
        b = b + fspecinv(fkernel(m, f[1], w)) #reject
        if t == None or (t != 'stop'):
            b = fspecinv(b) #bandpass        
    return b

## Spectral inversion
def fspecinv(b):
    b = -b
    b[int((b.shape[0]-1)/2)] = b[int((b.shape[0]-1)/2)]+1
    return b
#%%
def mfreqz(b,a,order,nyq_rate = 1):
    
    """
    Plot the impulse response of the filter in the frequency domain

    Parameters:
        
        b: numerator values of the transfer function (coefficients of the filter)
        a: denominator values of the transfer function (coefficients of the filter)
        
        order: order of the filter 
                
        nyq_rate = nyquist frequency
    """
    
    w,h = signal.freqz(b,a);
    h_dB = 20 * np.log10 (abs(h));
    
  
#%% Diseño de los filtros
    
def filter_design(srate, locutoff = 0, hicutoff = 0, revfilt = 0):
    #Constants
    TRANSWIDTHRATIO = 0.25;
    fNyquist = srate/2;  
    
    #The prototipical filter is the low-pass, we design a low pass and transform it
    if hicutoff == 0: #Convert highpass to inverted lowpass
        hicutoff = locutoff
        locutoff = 0
        revfilt = 1 #invert the logic for low-pass to high-pass and for
                    #band-pass to notch
    if locutoff > 0 and hicutoff > 0:
        edgeArray = np.array([locutoff , hicutoff])
    else:
        edgeArray = np.array([hicutoff]);
    
    #Not negative frequencies and not frequencies above Nyquist
    if np.any(edgeArray<0) or np.any(edgeArray >= fNyquist):
        print('Cutoff frequency out of range')
        return False  
    
    # Max stop-band width
    maxBWArray = edgeArray.copy() # Band-/highpass
    if revfilt == 0: # Band-/lowpass
        maxBWArray[-1] = fNyquist - edgeArray[-1];
    elif len(edgeArray) == 2: # Bandstop
        maxBWArray = np.diff(edgeArray) / 2;
    maxDf = np.min(maxBWArray);
    
    # Default filter order heuristic
    if revfilt == 1: # Highpass and bandstop
        df = np.min([np.max([maxDf * TRANSWIDTHRATIO, 2]) , maxDf]);
    else: # Lowpass and bandpass
        df = np.min([np.max([edgeArray[0] * TRANSWIDTHRATIO, 2]) , maxDf]);
    
    print(df)
    
    filtorder = 3.3 / (df / srate); # Hamming window
    filtorder = np.ceil(filtorder / 2) * 2; # Filter order must be even.
    
    # Passband edge to cutoff (transition band center; -6 dB)
    dfArray = [[df, [-df, df]] , [-df, [df, -df]]];
    cutoffArray = edgeArray + np.array(dfArray[revfilt][len(edgeArray) - 1]) / 2;
    print('pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): '+str(cutoffArray)+' Hz\n');
    # Window
    winArray = signal.hamming(int(filtorder) + 1);
    # Filter coefficients
    if revfilt == 1:
        filterTypeArray = ['high', 'stop'];
        b = firws(filtorder, cutoffArray / fNyquist, winArray, filterTypeArray[len(edgeArray) - 1]);
    else:
        b = firws(filtorder, cutoffArray / fNyquist, winArray);

    return filtorder, b;    


#design

#Se define función para procesar los datos dados en el archivo de las señales
#como parámetros se ingresa el canal que se quiere visualizar y el nombre del archivo entre comillas
def matrices_archivos(n,x):
    
    data = open(x+".txt", 'r') #Se lee el archivo de texto
    fila=data.readlines() #Se leen las filas del archivo
    EEG_C=np.zeros(len(fila)) #Se crea un vector donde se van a guardar los datos de las columnas
    for i in np.arange(6,len(fila)): #Se crea un ciclo for que recorre las filas del archivo
        dato=fila[i].find(' ')
        for j in (np.arange(n)): #Este ciclo for extrae los datos de la columna deseada
            dato1=fila[i].find(' ',dato+j)
            dato2=fila[i].find(',',dato1)
            dato=dato1
        EEG_C[i]=float(fila[i][dato1:dato2])
    EEG_C=EEG_C[6:len(EEG_C)]#se eliminan los ceros iniciales que quedan
    fs = 250 #Se define frecuencia de muestreo
        
    tiempo1=np.arange(0,len(EEG_C)/fs,1/fs) #se crea el vector de tiempo con la frecuencia de muestreo
    
    f,pot=senal_frecuencia=signal.welch(EEG_C,fs,'hanning', fs*2, fs) #Se extrae el periodograma de la señal
    
    Bajas_frecuencias=f[0:20] #Se identifican las frecuencias bajas con el periodograma
    Altas_frecuencias=f[20:len(f)] #Se identifian las frecuencias altas 
    potencia_bajas=pot[0:20] #Se identifican las potencias correspondientes a las frecuencias bajas
    potencia_altas=pot[20:len(pot)] #Se identifican las potencias correspondientes a las frecuencias altas
    
    for i in np.arange(len(potencia_bajas)): #Se crea una rutina para identificar la frecuencia baja con la mayor potencia
        if potencia_bajas[i]==np.max(potencia_bajas):
            f_pasa_altas=Bajas_frecuencias[i]
    
    for i in np.arange(len(potencia_altas)): #Se crea una rutina para identificar la freciencia alta con la mayor potencia
        if potencia_altas[i]==np.max(potencia_altas):
            f_pasa_bajas=Altas_frecuencias[i-20]
            
    order, lowpass = filter_design(fs, locutoff = 0, hicutoff = f_pasa_bajas, revfilt = 0); #Se define el filtro pasabajas
    #plot
    mfreqz(lowpass,1,order, fs/2);
    
    order, highpass = filter_design(fs, locutoff = f_pasa_altas, hicutoff = 0, revfilt = 1); #Se define el filtro pasa altas
    #plot
    mfreqz(highpass,1,order, fs/2);
    
    order, bandpass = filter_design(fs, locutoff = 4, hicutoff = 50, revfilt = 0);#Se define el filtro pasa banda
    #plot
    mfreqz(bandpass,1,order, fs/2);
    
    order, notch = filter_design(fs, locutoff = 50, hicutoff = 70, revfilt = 1);#Se define el rechaza banda
    
    mfreqz(notch,1,order, fs/2);
    
    senal_filtrada_pasabajas = signal.filtfilt(highpass, 1, EEG_C); #Se filtra con pasa altas la señal
    senal_filtrada_pasabajas = signal.filtfilt(lowpass, 1, senal_filtrada_pasabajas);#Se filtra con pasa bajas la señal
    senal_filtrada_pasabajas1=signal.filtfilt(lowpass, 1, senal_filtrada_pasabajas);
    
    histograma=np.histogram(senal_filtrada_pasabajas) #Se halla el histograma de la señal filtrada
    hist_muestras=np.zeros(0) #Se define un vector para guardar datos
    
    for i in np.arange(len(histograma[0])): #Se crea una rutina para guardar los valores de las voltajes en donde halla menos de 2000 muestras
        if histograma[0][i] < 2000:
            hist_muestras=np.append(hist_muestras,histograma[1][i])
            
    negativos=np.zeros(0) #Se crea matriz de ceros
    positivos=np.zeros(0) #Se crea matriz de ceros
    
    for i in hist_muestras: #Se crea rutina para guardar los valores de los voltajes negativos y positivos de las muestras filtradas del histograma
        if i<0:
            negativos=np.append(negativos,i)
        else:
            positivos=np.append(positivos,i)
    
    if positivos.any()==[0]: #Si no se encuentran valores donde halla menos de 2000 muestras se toma el mayor valor de voltaje como umbral superior
       positivos=np.max(senal_filtrada_pasabajas)
      
    if negativos.any()==[0]: #Si no se encuentran valores donde halla menos de 2000 muestras se toma el mayor valor de voltaje como umbral inferior
       negativos=np.min(senal_filtrada_pasabajas)
        
    cota_inferior=np.max(negativos) #Se extra la cota inferior para realizar el filtrado por umbrales
    cota_superior=np.min(positivos) #Se extrae la cota superior para realizar el filtrado mediante umbrales
    
    f2,pot2=senal_frecuencia=signal.welch(senal_filtrada_pasabajas,fs,'hanning', fs*2, fs) #Se halla el periodograma de la señal filtrada sin sacar las épocas
    
    multiplos=np.zeros(0) #Se crea una rutina para hallar los multiplos por 500 de la longitud de la señal filtrada
    for i in np.arange(len(EEG_C)):
        resto=i % 500
        if resto==0:
            multiplos=np.append(multiplos,i)
            
    valor=int(np.max(multiplos)) #Se halla el mallor multiplo de 500 de la señal filtrada         
    senal_filtrada_pasabajas=senal_filtrada_pasabajas[0:valor]#Se genera un vector de la señal filtrada hasta el valor del maximo multiplo
    
    
    muestras=np.split(senal_filtrada_pasabajas,500,axis=0) #Se generan las épocas cada 50o muestras
    promedio=np.mean(senal_filtrada_pasabajas)
    senal_resultante=np.zeros(0) #Se genera un vector vacio
    
    #plt.hist(senal_filtrada_pasabajas)
    #plt.show()
    
    for i in muestras: #Se genera una rutina para extraer las épocas y concatenar las de interes y generar la señal filtrada definitiva
        bandera=0
        for j in i:
            if j > cota_superior or j < cota_inferior:
                bandera=1
        if bandera==0:
            senal_resultante=np.append(senal_resultante,i)
    tiempo=np.arange(0,len(senal_resultante)/fs,1/fs)     
    
    return senal_filtrada_pasabajas1, senal_resultante, histograma,f,pot,f2,pot2,cota_inferior,cota_superior,f_pasa_bajas, f_pasa_altas, tiempo, tiempo1 #Retorna todas las variables calculadas durante la definición de la función



