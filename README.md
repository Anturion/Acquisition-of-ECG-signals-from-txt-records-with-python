# Acquisition-of-ECG-signals-from-txt-records-with-python
-The project takes the data from different ECG channels during, before and after surgery, exports them and filters them to obtain a signal with relevant information  

-El proyecto toma los datos provenientes de diferentes canales de ECG durante, antes y después de una cirugía, los exporta y los filtra para obtener una señal con información de relevancia

You must run the file Interfaz_segundo.py once in the PyQt interface you can put the name of the file with a format similar to that of the .txt files found in the confirmation example (P1_RAWEEG_2018-11-15_DuranteCirugia_42min.txt) Then click on Enterokay.

The program delivers the signal filtered by epochs and applies digital filters to eliminate noise from the electrical network and only leaves the relevant information from the ECG, it is also capable of calculating the upper and lower bounds for each channel and eliminating those portions to extract only the signal of interest

Debes correr el archivo Interfaz_segundo.py una vez en la interfaz de PyQt puedes poner el nombre del archivo con un formato parecido al de los archivos .txt que se encuentran en el ejemplo de confirmación (P1_RAWEEG_2018-11-15_DuranteCirugia_42min.txt) Posteriormente click en ingresar.

El programa entrega la señal filtrada por épocas y aplica filtros digitales para eliminar el ruido de la red eléctrica y solo deja la información relevante del ECG, también es capaz de calcular las cotas superriores e inferiores para cáda canal y eliminar esas porciones para lograr extraer solo la señal de interes
