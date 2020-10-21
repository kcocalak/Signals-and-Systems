#Nurettin Kağan Çocalak
#150160529

import numpy as np
import wave
import sys
import math
import warnings
warnings.simplefilter("ignore", DeprecationWarning)

def hpf(samples,fc,delta,gain):                     #high pass filter function definition with samples, cut off, delta and gain inputs
    alpha = delta/(delta + 1 /(2*np.pi*fc))*gain    #calculation of alpha with gain for filters
    resultFirstOrder = samples.copy()               #resultFirstOrder arrayz
    resultSecondOrder = samples.copy()              #resultSecondOrder array
    resulthpf = samples.copy()                      #resulthpf array

    resultFirstOrder[0] = alpha * samples[0]            #samples multiply with alpha for the first order filter
    resultSecondOrder[0] = alpha * resultFirstOrder[0]  #created first order filter multiply with alpha for the creation of second order filter
    resulthpf[0] = 0

    for i in range(len(samples) - 1):
        resultFirstOrder[i+1] = resultFirstOrder[i]*(1-alpha) + alpha*samples[i+1]              #calculation for first order low pass filter for len(samples) inputs
        resultSecondOrder[i+1] = resultSecondOrder[i]*(1-alpha) + alpha*resultFirstOrder[i+1]   #calculation for second order low pass filter thanks to first order
        resulthpf[i+1] = samples[i+1]-resultSecondOrder[i+1]                                    #calculation of high pass filter with subtruct low pass filter from samples

    return resulthpf

wavFileNameFirst ='Africa.wav'
wavFileNameSecond = 'WinnerTakesAll.wav'
outputFileNameWithoutGain = 'HPF.wav'

#read wav file
obj = wave.open(wavFileNameFirst,'rb')                  #Africa.wav is opened and readed
amplitudeWidth = obj.getsampwidth()                     #sample width to in bytes
frameRate = obj.getframerate()                          #sampling frequency
nTimesFrames = obj.getnframes()                         #number of audio frames
readFrames = obj.readframes(nTimesFrames)               #reads and returns at most n frames of audio, as a bytes object
samples = np.fromstring(readFrames, np.int16)           #create array with frames
obj.close()                                             #close the stream if it was opened by wave module

obj2 = wave.open(wavFileNameSecond,'rb')                #WinnerTakesAll.wav is opened and readed
amplitudeWidth2 = obj2.getsampwidth()                   #sample width to in bytes
frameRate2 = obj2.getframerate()                        #sampling frequency
nTimesFrames2 = obj2.getnframes()                       #number of audio frames
readFrames2 = obj2.readframes(nTimesFrames2)            #reads and returns at most n frames of audio, as a bytes object
samples2 = np.fromstring(readFrames2, np.int16)         #create array with frames
obj2.close()                                            #close the stream if it was opened by wave module

fc = 2000.0
delta = 1.0/44100.0
#calculate low pass filter
gain = 1
hpFiltered = hpf(samples,fc,delta,gain).astype(samples.dtype)           #hpf function is called for Africa.wav file
hpFiltered2 = hpf(samples2,fc,delta,gain).astype(samples2.dtype)        #hpf function is called for WinnerTakesAll.wav file

#write file with outputname
fileNameFirst = wavFileNameFirst.split(".",1)                           #write for outputfilename file name split for just Africa
outputFileNameNoGain = fileNameFirst[0] + outputFileNameWithoutGain     #and combine with HPF.wav
outputName = wave.open(outputFileNameNoGain,'w')                        #filtered file result open and write with new name which is created with outputFileNameNoGain
outputName.setparams((1,amplitudeWidth,frameRate,nTimesFrames,obj.getcomptype(),obj.getcompname()))  #accepts parameter tuple
outputName.writeframes(hpFiltered.tobytes('C'))                         #write audio frames and make sure they are correct
outputName.close()                                                      #close the stream if it was opened by wave module

fileNameSecond = wavFileNameSecond.split(".",1)                          #write for outputfilename file name split for just WinnerTakesAll
outputFileNameNoGain = fileNameSecond[0] + outputFileNameWithoutGain     #and combine with HPF.wav
outputName = wave.open(outputFileNameNoGain,'w')                         #filtered file result open and write with new name which is created with outputFileNameNoGain
outputName.setparams((1,amplitudeWidth2,frameRate2,nTimesFrames2,obj2.getcomptype(),obj2.getcompname()))    #accepts parameter tuple
outputName.writeframes(hpFiltered2.tobytes('C'))                         #write audio frames and make sure they are correct
outputName.close()                                                       #close the stream if it was opened by wave module
