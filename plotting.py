# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:13:24 2023

@author: ruben
"""

import matplotlib.pyplot as plt

# F1
x = [40, 80, 120, 160]
ys = [[5688, 17758, 37789, 64212], #[[19932, 61612, 137139, 217356],   # 0.25 
     [13695, 51318, 121494, 206816],    # 0.5  
     [27927, 107064, 254429, 434079]]    # 0.75

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
ys = [[7806, 27268, 55326, 95748],   # 0.25 
     [25246, 89258, 192247, 339028],    # 0.5  
     [53418, 190459, 413796, 730680]]    # 0.75

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
ys = [[52090, 110611, 203312, 344011],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [47760, 107410, 174776, 286718],    # 0.5  
     [44755, 139756, 296978, 518617]]    # 0.75

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
ys = [[51521, 117306, 279272, 370100],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [41714, 109080, 241091, 395885],    # 0.5  
     [61871, 199232, 450671, 746933]]    # 0.75

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