# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 07:32:30 2017

@author: Huda Ibrhm
"""



import csv
import numpy

#read file training
file1 = open('dataTraining.csv')
csv_reader1= csv.reader(file1)

#inisiasi list training
n = []
for row in csv_reader1:
    n.append(row)

njml = len(n)

A1 = []
A2 = []
for b1 in n:
    A1.append(b1[1])
    A2.append(b1[2])
print(A1)
A1setosa = []
A1versi = []
A2setosa = []
A2versi = []

for x in range(1,njml):
    if x<=40:
        A1setosa.append(float(A1[x]))
        A2setosa.append(float(A2[x]))
    else:
        A1versi.append(float(A1[x]))
        A2versi.append(float(A2[x]))        


file1.close()

#data testing
file2 = open('dataTesting.csv')
csv_reader2 = csv.reader(file2)

nTest = []
for row in csv_reader2:
    nTest.append(row)

njmlTest = len(nTest)

A1Test = []
A2Test = []
for b2 in nTest:
    A1Test.append(b2[1])
    A2Test.append(b2[2])


A1setosaTest = []
A1versiTest = []
A2setosaTest = []
A2versiTest = []
for x in range(0,njmlTest):
    if x<10:
        A1setosaTest.append(float(A1Test[x]))
        A2setosaTest.append(float(A2Test[x]))                
    else:
        A1versiTest.append(float(A1Test[x]))
        A2versiTest.append(float(A2Test[x]))
        
njmlTestSetosa = len(A1setosaTest);
njmlTestVersi = len(A2versiTest);

print(njmlTestVersi)

file2.close()


#Bayes

"""Mencari rata rata"""
mean01 = numpy.mean(A1setosa)
mean02 = numpy.mean(A2setosa)
mean11 = numpy.mean(A1versi)
mean12 = numpy.mean(A2versi)

"""Mencari variansi"""
var01 = numpy.var(A1setosa)
var02 = numpy.var(A2setosa)
var11 = numpy.var(A1versi)
var12 = numpy.var(A2versi)

""" Mencari 1/sqrt(2*phi*data training)"""
akar10 = 1/numpy.sqrt(var01*2*22/7)
akar20 = 1/numpy.sqrt(var02*2*22/7)
akar11 = 1/numpy.sqrt(var11*2*22/7)
akar21 = 1/numpy.sqrt(var12*2*22/7)

#menghitung exp data training 

expA10 = [0]*njmlTestSetosa
expA20 = [0]*njmlTestSetosa
expA11 = [0]*njmlTestVersi
expA21 = [0]*njmlTestVersi

for x in range(0,njmlTestSetosa):
    expA10[x]=numpy.exp((-1*((A1setosaTest[x]-mean01)**2))/(2*var01))
    expA20[x]=numpy.exp((-1*((A2setosaTest[x]-mean02)**2))/(2*var02))
    expA11[x]=numpy.exp((-1*((A1versiTest[x]-mean11)**2))/(2*var11))
    expA21[x]=numpy.exp((-1*((A2versiTest[x]-mean12)**2))/(2*var12))
print("nilaiexp")
print(expA21)
#kelas lainnya

expA10v2 = [0]*njmlTestSetosa
expA20v2 = [0]*njmlTestSetosa
expA11v2 = [0]*njmlTestVersi
expA21v2 = [0]*njmlTestVersi

for x in range(0,njmlTestSetosa):
    expA10v2[x]=numpy.exp((-1*((A1setosaTest[x]-mean11)**2))/(2*var11))
    expA20v2[x]=numpy.exp((-1*((A2setosaTest[x]-mean12)**2))/(2*var12))
    expA11v2[x]=numpy.exp((-1*((A1versiTest[x]-mean01)**2))/(2*var01))
    expA21v2[x]=numpy.exp((-1*((A2versiTest[x]-mean02)**2))/(2*var02))
print("a")
print(A2versiTest[1])
#menghitung hasil peluang kelasnya    
pel_hasil0 = [0]*njmlTestSetosa
pel_hasil1 = [0]*njmlTestSetosa
pel_hasil0v2 = [0]*njmlTestVersi
pel_hasil1v2 = [0]*njmlTestVersi

for x in range(0,njmlTestSetosa):
    pel_hasil0[x]=akar10*expA10[x]*akar20*expA20[x]
    pel_hasil1[x]=akar11*expA11[x]*akar21*expA21[x]
    pel_hasil0v2[x]=akar11*expA10v2[x]*akar21*expA20v2[x]
    pel_hasil1v2[x]=akar10*expA11v2[x]*akar20*expA21v2[x]

#mencetak nilai2 training
print("Nilai data Training :")
print("A1 Bunga Setosa :",A1setosa)
print()
print("A2 Bunga Setosa :",A2setosa)
print()
print("A1 Bunga Versicolor :",A1versi)
print()
print("A2 Bunga Versicolor :",A2versi)

#mencetak nilai rata2 data training
print("Rata2 data training :")
print(mean01)
print(mean02)
print(mean11)
print(mean12)
print()

#mencetak nilai variansi data training
print("Varians data training :")
print(var01)
print(var02)
print(var11)
print(var12)
print()

#mencetak nilai2 data testing
print("Nilai data Testing :")
print("A1 Bunga Setosa :",A1setosaTest)
print()
print("A2 Bunga Setosa :",A2setosaTest)
print()
print("A1 Bunga Versicolor :",A1versiTest)
print()
print("A2 Bunga Versicolor :",A2versiTest)

#mengoutputkan hasil p dari perhitungan        
print("Hasil P setosa = ",pel_hasil0);
print("Hasil P setosa v2 =",pel_hasil0v2);    
print("Hasil P versicolor = ",pel_hasil1);
print("Hasil P versicolor v2 = ",pel_hasil1v2); 

print("hasil Perbandingan Setosa && SetosaV2 =")
for x in range(0,njmlTestSetosa):

    if pel_hasil0[x] > pel_hasil0v2[x] :
        print("TRUE")
    else:
        print("FALSE")
        
print()
print("hasil Perbandingan Versicolor && VersicolorV2 =")
for x in range(0,njmlTestSetosa):
    if pel_hasil1[x] > pel_hasil1v2[x] :
        print("TRUE")
    else:
        print("FALSE")

#mengukur Akurasi
akurasi = 0
for x in range(0,njmlTestSetosa):
    if pel_hasil0[x] > pel_hasil0v2[x]:
        akurasi += 1 
    
    if pel_hasil1[x] > pel_hasil1v2[x]:
        akurasi += 1


print()
print("Akurasi yang didapat adalah =",akurasi*10/2, " persen")












