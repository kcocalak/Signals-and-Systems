#Nurettin Kağan Çocalak
#150160529

import numpy as np
import wave
import sys
import math
import warnings
warnings.simplefilter("ignore", DeprecationWarning)

def lpf(samples,fc,delta,gain):                             #low pass filter function definition with samples, cut off, delta and gain inputs
    alpha = delta/(delta + 1 /(2*np.pi*fc))*gain            #calculation of alpha with gain for filters
    resultFirstOrder = samples.copy()                       #resultFirstOrder arrayz
    resultSecondOrder = samples.copy()                      #resultSecondOrder array

    resultFirstOrder[0] = alpha * samples[0]                #samples multiply with alpha for the first order filter
    resultSecondOrder[0] = alpha * resultFirstOrder[0]      #created first order filter multiply with alpha for the creation of second order filter

    for i in range(len(samples) - 1):

        resultFirstOrder[i+1] = resultFirstOrder[i]*(1-alpha) + alpha*samples[i+1]              #calculation for first order low pass filter for len(samples) inputs
        resultSecondOrder[i+1] = resultSecondOrder[i]*(1-alpha) + alpha*resultFirstOrder[i+1]   #calculation for second order low pass filter thanks to first order

    return resultSecondOrder

wavFileNameFirst ='Africa.wav'
wavFileNameSecond = 'WinnerTakesAll.wav'
outputFileNameWithGain = '5dBGainLPF.wav'
outputFileNameWithoutGain = '0dBGainLPF.wav'

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
lpFilteredNoGain = lpf(samples,fc,delta,gain).astype(samples.dtype)             #lpf function is called for Africa.wav file with no gain
gain = pow(10,0.25)
lpFiltered = lpf(samples,fc,delta,gain).astype(samples.dtype)                   #lpf function is called for Africa.wav file with 5dB gain
gain = 1
lpFilteredNoGain2 = lpf(samples2,fc,delta,gain).astype(samples2.dtype)          #lpf function is called for WinnerTakesAll.wav file with no gain
gain = pow(10,0.25)
lpFiltered2 = lpf(samples2,fc,delta,gain).astype(samples2.dtype)                #lpf function is called for WinnerTakesAll.wav file with 5dB gain

#write file with outputname for Africa0dBGainLPF
fileNameFirst = wavFileNameFirst.split(".",1)                           #write for outputfilename file name split for just Africa
outputFileNameNoGain = fileNameFirst[0] + outputFileNameWithoutGain     #and combine with 0dBGainLPF
outputName = wave.open(outputFileNameNoGain,'w')                        #filtered file result open and write with new name which is created with outputFileNameWithoutGain
outputName.setparams((1,amplitudeWidth,frameRate,nTimesFrames,obj.getcomptype(),obj.getcompname()))  #accepts parameter tuple
outputName.writeframes(lpFiltered.tobytes('C'))                         #write audio frames and make sure they are correct
outputName.close()                                                      #close the stream if it was opened by wave module

#write file with outputname for Africa5dBGainLPF
outputFileName = fileNameFirst[0] + outputFileNameWithGain              #combine with 5dBGainLPF
outputName = wave.open(outputFileName,'w')                              #filtered file result open and write with new name which is created with outputFileNameWithGain
outputName.setparams((1,amplitudeWidth,frameRate,nTimesFrames,obj.getcomptype(),obj.getcompname()))     #accepts parameter tuple
outputName.writeframes(lpFiltered.tobytes('C'))                         #write audio frames and make sure they are correct
outputName.close()                                                      #close the stream if it was opened by wave module

#write file with outputname for WinnerTakesAll0dBGainLPF
fileNameSecond = wavFileNameSecond.split(".",1)                          #write for outputfilename file name split for just WinnerTakesAll
outputFileNameNoGain = fileNameSecond[0] + outputFileNameWithoutGain     #and combine with HPF.wav
outputName = wave.open(outputFileNameNoGain,'w')                         #filtered file result open and write with new name which is created with outputFileNameWithoutGain
outputName.setparams((1,amplitudeWidth2,frameRate2,nTimesFrames2,obj2.getcomptype(),obj2.getcompname()))    #accepts parameter tuple
outputName.writeframes(lpFiltered2.tobytes('C'))                         #write audio frames and make sure they are correct
outputName.close()                                                       #close the stream if it was opened by wave module

#write file with outputname for WinnerTakesAll5dBGainLPF
outputFileName = fileNameSecond[0] + outputFileNameWithGain             #combine with 5dBGainLPF
outputName = wave.open(outputFileName,'w')                              #filtered file result open and write with new name which is created with outputFileNameWithGain
outputName.setparams((1,amplitudeWidth2,frameRate2,nTimesFrames2,obj.getcomptype(),obj.getcompname()))     #accepts parameter tuple
outputName.writeframes(lpFiltered2.tobytes('C'))                        #write audio frames and make sure they are correct
outputName.close()                                                      #close the stream if it was opened by wave module
