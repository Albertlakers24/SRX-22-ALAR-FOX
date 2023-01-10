import pandas as pd
filenameAF = 'naca633618.txt'
filenameWing = 'naca633618_wing.txt'
filenameAF_xfoil = 'naca633618_XFOIL.txt'
import numpy as np
import matplotlib.pyplot as plt


df_AF_XFLR = pd.read_csv(filenameAF,
                 sep="\s+",
                 skiprows=1,
                 usecols=[0,1,2,3,4,5,6,7,8,9,10,11],
                 names=['Alpha' ,'CL', 'CD', 'CDp', 'Cm','Top Xtr','Bot Xtr', 'Cpmin','Chinge','Cni','QInf','XCP'])

df_wing = pd.read_csv(filenameWing,
                 sep="\s+",
                 skiprows=1,
                 usecols=[0,1,2,3,4,5,6,7,8,9,10,11, 12],
                 names=['Alpha' ,'beta','CL', 'CDi', 'CDv', 'CD','CY','Cl', 'Cm','Cn','Cni','QInf','XCP'])

df_AF_Xfoil = pd.read_csv(filenameAF_xfoil,
                 sep="\s+",
                 skiprows=1,
                 usecols=[0,1,2,3,4,5,6],
                 names=['Alpha' ,'CL', 'CD', 'CDp', 'CM','Top Xtr','Bot Xtr'])


Alpha_AF_xflr = np.array(df_AF_XFLR['Alpha'])
CL_AF_xflr = np.array(df_AF_XFLR['CL'])
CD_AF_xflr = np.array(df_AF_XFLR['CD'])

Alpha_Wing = np.array(df_wing['Alpha'])
CL_Wing = np.array(df_wing['CL'])

Alpha_AF_xfoil = np.array(df_AF_Xfoil['Alpha'])
CL_AF_xfoil = np.array(df_AF_Xfoil['CL'])

# CL_CD = CL_AF_xflr/CD_AF_xflr
# print(CL_CD)
# print(CL_CD.max())
# ind = np.where(CL_CD == CL_CD.max())
# print(CL_AF_xflr[ind])
# print(Alpha_AF_xflr[ind])

# print(CL_AF, CD)
# print(CL_AF.max())
# ind = np.where(CL_AF == CL_AF.max())
# print(CL_AF[ind])
# print(Alpha_AF[ind])
# plt.plot(Alpha_AF, CL_AF)
# plt.show()

