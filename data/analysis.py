# -*- coding: cp1251 -*-
# получить среднее значение осцилляции

import os, sys
import getopt
import string
import time
import msvcrt
import re
import random  ????

from numpy import *
from numpy.fft import *
from pylab import *

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt


#Нижняя и верхняя границы номера в названий входных файлов 
clipB = 1
clipT = 15


#Функция поиска комплексов (каждый комплекс есть пара [точка минимума, точка максимума])
def FindRawShapes (data):	
	dx = 10
	do = 0.60
	
	complexDict = {}
	complexes = []
	
	minx = -1
	maxx = -1
	for i in range(1, len(data) - 1):
		j = len(data) - i - 1
		if data[j] >= data[j + 1] and data[j] > data[j - 1]: 
			maxx = j
			minx = -1
		if data[j] <= data[j - 1] and data[j] < data[j + 1]: minx = j
		if data[j] < data[j + 1] and maxx != -1: descentingLength =  descentingLength + 1
		else: descentingLength = 0
		if descentingLength >= dx and minx != -1 and maxx != -1:
			h = -1
			for hh in complexDict.keys(): 
				if min(data[maxx] - data[minx], hh) / max(data[maxx] - data[minx], hh) > do: h = hh 
			if h == -1: complexDict[data[maxx] - data[minx]] = array([minx,maxx], ndmin = 2)
			else: complexDict[h] = insert(complexDict[h], 0, array([minx, maxx], ndmin = 2), axis = 0)
			descentingLength = 0
			minx = -1
			maxx = -1
	for l in complexDict.values():
		if len(l) > len(complexes): complexes = l
	complexes = delete(complexes, 0, 0)
	return complexes
	
	
#Начало


pressArray = loadtxt("pressArray.txt") #Массив со средним давлением для каждого из файлов clip
allIndices = [] #Список массивов комплексов
allArr1 = [] #Список входных массивов

#Парсим M_25
file = open("M_25.txt", "r")
text = file.read()
m = re.findall('^[\d,]+\t[\d,]+\t([\d,]+)\t[\d,]+\n', text, flags = re.MULTILINE)
referenceAverage = map(lambda x : float(x.replace(',', '.')), m)

#Считываем входные массивы, находим комплексы и выводим все графики
NameList = map (lambda x: 'clip_' + '%2.2d'%x + '.txt', xrange (clipB, clipT))
for i in xrange(len(NameList)):
	arr1 = loadtxt (NameList[i])
	allArr1.append(arr1)
	indices = FindRawShapes (arr1)
	allIndices.append(indices)
	subplot(len(NameList), 1, i)
	plot(arr1, 'b-')
	for index in indices:
		plot(index, map(lambda x: arr1[x], index), 'ro')
show()






avgArrayList = [] #Список из средних комплексов
avgNameList = [] #Список имен средних комплексов

#Находим средние комплексы

for j in range(len(allIndices)):

    
	OutFileName = "average_" + str(j + 1) + ".txt"
    
	indices = allIndices[j]
	arr1 = allArr1[j]

    # Берем последние 8 комплексов и находим среднюю ширину комплекса
    
	summ_width = 0
	if(len(indices) > 9): indices = indices[len(indices) - 9:]
	for i in xrange(len(indices) - 1):
		summ_width = summ_width + indices[i + 1, 0] - indices[i, 0]
	averageWidth = 	summ_width / (len(indices) - 1)

	# Проводим корреляцию
	badComplexes = []
	ref = arr1[indices[0, 0] : indices[0, 0] + averageWidth]
	for i in xrange (len(indices) - 1):
		start = indices[i, 0] - 20
		corr = correlate(ref, arr1[start: start + averageWidth + 20], mode = "same")
		delta = int (averageWidth / 2 - argmax(corr))
		if abs(delta) > 10: badComplexes.append(i)
		indices[i, 0] = indices[i, 0] + delta
		indices[i, 1] = indices[i, 1] + delta	
	
    #Фильтруем нескореллированные комплексы	
	k = 0
	for i in badComplexes:
		indices = np.delete(indices, i - k, 0)
		k = k + 1
		
	if len(indices) > 4:


		# Создаем двумерный массив из комплексов, который представляет собой набор точек, а не просто пары [точка минимума, точка максимума]
		shapeArr = empty ([averageWidth, len(indices) - 1])
		for i in xrange (len(indices) - 1):
			#print 'i = ', i
			shapeArr[: ,i] = arr1[indices[i, 0] : indices[i, 0] + averageWidth]
			
		# Получаем среднее
		averageArr = empty(averageWidth)
		for i in xrange(averageWidth):
			averageArr[i] = average (shapeArr[i])
		savetxt (OutFileName, averageArr, "%f")
		
		
		# Рисуем
		subplot(4, 4, j) 
		plt.title("average_" + str(j + 1) + ".txt with p = " + str(pressArray[j]))
		for i in xrange(len(indices) - 1):
			plot(xrange(int(averageWidth)), shapeArr[:,i], 'b-')
		plot(xrange(int(averageWidth)), averageArr, 'ro')
		
		#Нормируем среднее
		averageArr = map(lambda x : x - averageArr[0], averageArr)
		averageArr = map(lambda x : x / max(averageArr) * 1000, averageArr)
	
		avgArrayList.append (averageArr)
		avgNameList.append("average_" + str(j + 1) + ".txt with p = " + str(pressArray[j]))
show()




#Выводим отнормированные средние и референсный комплекс (его тоже нормируем)

for j in range(len(avgArrayList)):
	subplot(4, 4, j) 
	plt.title(avgNameList[j])
	plot(avgArrayList[j], 'r-')
	
referenceAverage = map(lambda x : x - referenceAverage[0], referenceAverage)
referenceAverage = map(lambda x : x / max(referenceAverage) * 1000, referenceAverage)
subplot(4, 4, len(avgArrayList))
plt.title("referenceAverage")
plot(referenceAverage, 'r-')

show()

#Все средние делаем одной длины
minLength = min(map(lambda x : len(x), avgArrayList))
avgArrayList = map(lambda x : x[:minLength], avgArrayList)


#Отображаем среднее в 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X = [range(minLength)] * len(avgArrayList)
Y = map(lambda x : [x] ,range(len(avgArrayList)))
Z = avgArrayList

ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
plt.show()






#Находим спектры
AllFFTList = []
for avg in avgArrayList:
	avgFFT = fft(avg);
	avgFFTAmp = map (lambda x: abs(x), avgFFT)
	avgFFTAmp = avgFFTAmp[:12]
	semilogy (avgFFTAmp, 'b')
	avgFFTAmpLog = map(lambda x: math.log(x, 10), avgFFTAmp)
	AllFFTList.append(avgFFTAmpLog)

# avgFFT = fft(avgArrayList[1]);
# avgFFTAmp = map (lambda x: abs(x), avgFFT)
# semilogy (avgFFTAmp[:12], 'b-')

# avgFFT = fft(map(lambda x: x * 10, avgArrayList[1]));
# avgFFTAmp = map (lambda x: abs(x), avgFFT)
# semilogy (avgFFTAmp[:12], 'g-')

referenceFFT = fft(referenceAverage)
referenceFFTAmp = map (lambda x: abs(x), referenceFFT)[:12]
semilogy(referenceFFTAmp, 'ro')
	
show()

#Отобразить спектры в 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.zaxis.set_scale('log')

X = [range(12)] * len(AllFFTList)
Y = map(lambda x : [x] ,range(len(AllFFTList)))
Z = AllFFTList

ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
plt.show()

sys.exit()


	
	
#print averArray
arr2 = ifft(fft1)
fft2 = fft(arr2)

file = open (OutFileName, "wt")
#file.write ("PeriphAmp\tPeriphPhase\tAortalAmp\tAortalPhase\tFilterTime\tFilterAmp\tFilterPhase\tFrq\n");
for i in xrange (len(fft1)):
	file.write ("%f\t%f\t%f\t%f\t%f\t%f\n" % (arr1[i], abs(fft1[i]), angle(fft1[i]), real(arr2[i]), abs(fft2[i]), angle(fft2[i])))
file.close()

sys.exit()
	
file = open (InFileName, "rt")
for i in file.readlines()[4:]:
	j = i.replace(',','.').split()
	if len(j) > 2:
		sig1.append(float(j[2]))
		sig2.append(float(j[3]))		
file.close()

# for i in xrange (lenArray - len(sig1)):
	# sig1.append(0.0)
	# sig2.append(0.0)

frqSlice = 128.0 / len(sig1)	
	
signal1 = array(sig1)
signal2 = array(sig2)
fourier1 = fft(signal1)
fourier2 = fft(signal2)

for a, b in zip (fourier1, fourier2):
	filterFreq.append (b / a)

filterTime = ifft (array(filterFreq))





