import pandas as pd
filename = 'naca633618.txt'
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv(filename,
                 sep="\s+",
                 skiprows=1,
                 usecols=[0,1,2,3,4,5,6,7,8,9,10,11],
                 names=['Alpha' ,'CL', 'CDi', 'CDv', 'CD','CY','Cl', 'Cm','Cn','Cni','QInf','XCP'])

print(df['Alpha'])

Alpha_AF = np.array(df['Alpha'])
CL_AF = np.array(df['CL'])
CD = np.array(df['CDi'])
CM = np.array(df['Cm'])
# CL_CD = CL/CD
# print(CL_CD)
# print(CL_CD.max())
# ind = np.where(CL_CD == CL_CD.max())
# print(CL[ind])
# print(Alpha[ind])
# print(CL, CD)
# print(CL_AF.max())
# ind = np.where(CL == CL.max())
# print(CL[ind])
# print(Alpha[ind])
# plt.plot(Alpha_AF, CL_AF)
# plt.show()

