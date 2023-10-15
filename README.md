# Silpakorn-Forearm-Vital-Signs
The contents of this repository describes the structure of our proposed method (ica+)for remote pulse signal measurement in non-facial regions. 
Notably, we have observed that the subtle color variations in non-facial region images exhibit reduced intensity compared to facial regions due to the skin's thickness and lower blood flow.
To address this, our method aims to enhance the accuracy and reliability of pulse signal detection by carefully smoothing and detrending the noise without deforming the signal.
We achieve this by employing Independent Component Analysis (ICA) to separate the most pulse-like component. 
Furthermore, we utilize similarity measurements based on the auto-correlation function to emphasize the most cyclic component. 
The results indicated that our proposed method effectively handles noise signals , particularly under fluorescent lighting conditions, yielding reliable results when compared to established iPPG techniques.
