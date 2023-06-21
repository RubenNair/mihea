# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:13:24 2023

@author: ruben
"""

import matplotlib.pyplot as plt

# F1
x = [40, 80, 120, 160]
ys = [[4236, 14849, 33709, 58528], #[[19932, 61612, 137139, 217356],   # 0.25 
     [14058, 51816, 114000, 184206],    # 0.5  
     [29527, 101776, 225688, 355719]]    # 0.75

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
ys = [[6383, 24019, 47861, 83815],   # 0.25 
     [23373, 82483, 171398, 319308],    # 0.5  
     [46893, 175098, 385119, 684601]]    # 0.75

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
ys = [[-1, -1, -1, -1],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [1154892, 2498129, 518911, 889026],    # 0.5  
     [39761, 130092, 275686, 469221]]    # 0.75

#plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
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
ys = [[-1, -1, -1, -1],   # 0.25  DIDN'T HAVE ANY VTR HITS
     [34013, 87202, 195221, 324518],    # 0.5  
     [55594, 175926, 413803, 709015]]    # 0.75

#plt.loglog(x, ys[0], 'ro-', label='f=0.25', color='red')
plt.loglog(x, ys[1], 'ro-', label='f=0.5', color='green')
plt.loglog(x, ys[2], 'ro-', label='f=0.75', color='blue')
plt.xticks(x)
plt.title('F4 TrapEllipse')
plt.xlabel('Problem Size')
plt.ylabel('Number of Evaluations')
plt.legend()
plt.show()