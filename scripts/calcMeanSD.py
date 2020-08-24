#! /usr/bin/python

def calcMeanSD(data):
    ''' calculate standard deviation of a list of values
    Inputs are mean of list, length of list, and actual list'''
    mean = sum(data) / float(len(data))
    n = len(data)
    
    sum_diff_square_mean = 0
    for i in range(len(data)):
        sum_diff_square_mean += (data[i] - mean)**2
    variance = sum_diff_square_mean / float(n - 1.0)
    sd = variance**(0.5)
    return (mean, sd)
