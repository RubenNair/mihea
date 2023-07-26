# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:13:24 2023

@author: ruben
"""

import matplotlib.pyplot as plt

# F1
x = [40, 80, 120, 160]
ys = [[4216, 15142, 34456, 60076], #[[19932, 61612, 137139, 217356],   # 0.25 
     [14063, 53117, 124221, 211174],    # 0.5  
     [29644, 110814, 261951, 450035]]    # 0.75

plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
plt.loglog(x, ys[1], 'ro-', label='f=0.5', color='green')
plt.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue')
plt.xticks(x)
plt.title('F1 OneMaxSphere')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()

# F2
plt.figure()
x = [40, 80, 120, 160]
ys = [[6627, 24968, 49312, 87558],   # 0.25 
     [24174, 83250, 177455, 322287],    # 0.5  
     [50256, 181585, 401100, 833628]]    # 0.75

plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
plt.loglog(x, ys[1], 'ro-', label='f=0.5', color='green')
plt.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue')
plt.xticks(x)
plt.title('F2 OneMaxEllipse')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()

# F3
plt.figure()
x = [40, 80, 120, 160]
ys = [[39333, 90051, 138947, 228831],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [33597, 82496, 142046, 246402],    # 0.5  
     [41442, 136503, 288474, 529023]]    # 0.75

plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
plt.loglog(x, ys[1], 'ro-', label='f=0.5', color='green')
plt.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue')
plt.xticks(x)
plt.title('F3 TrapSphere')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()


# F4
plt.figure()
x = [40, 80, 120, 160]
ys = [[38827, 91340, 175779, 260535],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [33967, 94766, 203402, 345474],    # 0.5  
     [59661, 184399, 425302, 926610]]    # 0.75

plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
plt.loglog(x, ys[1], 'ro-', label='f=0.5', color='green')
plt.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue')
plt.xticks(x)
plt.title('F4 TrapEllipse')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()

# F5
plt.figure()
x = [20, 40, 60]
ys = [[20114, -1, -1],   # a = 1.15
     [21818, 53865, 100150],    # a = 2.75  
     [21530, 67267, 108328]]    # a = 3

plt.loglog(x, ys[0], 'ro-', label='a=1.15', color='red')
plt.loglog(x, ys[1], 'ro-', label='a=2.75', color='green')
plt.loglog(x, ys[2], 'ro-', label='a=3.00', color='blue')
plt.xticks(x)
plt.title('F5 Trap Block Ellipse')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()