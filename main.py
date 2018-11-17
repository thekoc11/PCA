import csv
import matplotlib.pyplot as plt
from numpy import transpose
from numpy import diag
from numpy import dot
from numpy.linalg import inv
from numpy import array
from numpy.linalg import eig
from numpy import mean
from numpy import std
import numpy as np


def quickSort(alist, c):
    quickSortHelper(alist, c, 0, len(alist) - 1)


def quickSortHelper(alist, c, first, last):
    if first < last:
        split = partition(alist, c, first, last)
        quickSortHelper(alist, c, first, split - 1)
        quickSortHelper(alist, c, split + 1, last)


def partition(alist, c, first, last):
    pivot = alist[first]
    leftmark = first + 1
    rightmark = last

    done = False

    while not done:
        while leftmark <= rightmark and alist[leftmark] >= pivot:
            leftmark = leftmark + 1

        while alist[rightmark] <= pivot and rightmark >= leftmark:
            rightmark = rightmark - 1

        if rightmark < leftmark:
            done = True
        else:
            temp = alist[leftmark]
            alist[leftmark] = alist[rightmark]
            alist[rightmark] = temp

            ctemp = c[leftmark]
            c[leftmark] = c[rightmark]
            c[rightmark] = ctemp
    temp = alist[first]
    alist[first] = alist[rightmark]
    alist[rightmark] = temp

    ctemp = c[first]
    c[first] = c[rightmark]
    c[rightmark] = ctemp

    return rightmark


data = []
rlen = 0
timedata = []
trlen = 0

try:
    with open('gene.csv', 'rb') as csvfile:
        datareader = csv.reader(csvfile)
        for row in datareader:
            data.append(row)
            rlen = len(row)

    for i in range(len(data)):
        for j in range(rlen):
            data[i][j] = float(data[i][j])

    print(len(data), rlen)

    with open('Meta.csv', 'rb') as csvfile:
        datareader = csv.reader(csvfile)
        for row in datareader:
            timedata.append(row)
            trlen = len(row)

    for i in range(len(timedata)):
        for j in range(trlen):
            timedata[i][j] = float(timedata[i][j])


    print(len(timedata), trlen)


    X = array(data)
    Y = array(timedata)
    Z = []

    xbar = mean(X, axis=0, dtype=np.float32)
    sigma = std(X, axis=0, dtype=np.float32)

    for row in X:
        row = row - xbar
        for i in range(len(row)):
            row[i] = row[i]/sigma[i]
        Z.append(row)

    Z = array(Z)
    Zt = transpose(Z)

    M = Zt.dot(Z)

    d, P = eig(M)
    D = diag(d)
    Pinv = inv(P)
    max = d[0]

    # test = array([31,26,20,17,44,77,55,93])
    c = np.arange(len(d))
    # print(d)
    # print(c)

    quickSort(d, c)

    # print(d)
    # print(c)
    Pstar = []

    for i in range(len(P)):
        Pstar.append(P[c[i]])

    Dp = diag(d)
    Pstar = array(Pstar)

    Zstar = Z.dot(Pstar)

    # print (Zstar)

    pc1 = Zstar[0]
    pc2 = Y

    with plt.style.context('seaborn-whitegrid'):
        plt.scatter(pc1, pc2)
        plt.tight_layout()
        plt.show()
except:
    print("Wrong Format!")
