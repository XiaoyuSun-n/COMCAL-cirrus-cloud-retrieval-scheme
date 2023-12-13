# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 14:12:13 2022

@author: Lenovo
"""

import netCDF4 as nc
import numpy as np
import os,time
import matplotlib.pyplot as plt
import colorcet as cc
import scipy.io as scio

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def WCT_calculation(a,z_,var):
    
    z=z_
    P=var
    
    WCT=np.zeros(z.shape[0])
    
    for i in range(z.shape[0]):
        zb=z[i]-a/2
        zt=z[i]+a/2
        zin=z[np.where((z>=zb)&(z<=zt))]
        Pin=P[np.where((z>=zb)&(z<=zt))]
        
        dp=np.zeros(zin.shape[0])
        h=np.zeros(zin.shape[0])
        b=z[i]
        # b=14000*np.ones(np.shape(b)[0])
        for j in range (zin.shape[0]):
            if ( zin[j]>=(b-a/2)) & ( zin[j]<b):
                h[j]=1
            elif ( zin[j] >b) & ( zin[j] <=(b+a/2)) :
                h[j]=-1
            else:
                h[j]=0
            #print(h[j])
            dp[j]=Pin[j]*h[j]
            WCT[i]=WCT[i]+dp[j]
        WCT[i]=WCT[i]
    return WCT

def findzbandzt (WCT,z):
 
    startheightmax=10
    # startheightmin=5
    startheightmin=10
    
    endheight=np.where(z<18500)[0][-1]
    zmax=np.where(WCT==np.max(WCT[startheightmax:endheight]))[0]
    zmin=np.where(WCT==np.min(WCT[startheightmin:endheight]))[0]
    # zt=find_nearest(WCT[zmax[0]:endheight],Ws)+zmax[0]
    # zb=find_nearest(WCT[startheight:zmin[0]],-1*Ws)+startheight
    WCT_rangeup=np.sort(WCT[zmax[0]:endheight])
    WCT_rangedown=np.sort(WCT[startheightmin:zmin[0]])
    
    Wsmax=np.max(WCT[startheightmax:endheight])
    Wsmin=np.min(WCT[startheightmax:endheight])
    
    Ws1=np.max(WCT[startheightmax:endheight])-Wsmax/10.0
    #Ws1=2.5
        # Ws2=WCT[zmin]+10.0
    Ws2=np.min(WCT[startheightmax:endheight])+Wsmin/-5.0
    
    try:
        zt_temp=find_closet_ele(WCT_rangeup, Ws1, 4)
        zb_temp=find_closet_ele(WCT_rangedown,Ws2, 4)
        ztidx_temp=np.empty(4,dtype=int)
        zbidx_temp=np.empty(4,dtype=int)
        for i in range(4):
            if zt_temp[i] !=np.inf:
                ztidx_temp[i]=np.where(WCT==zt_temp[i])[0][0]
            else:
                ztidx_temp[i]=999
            if zb_temp[i] !=np.inf:
                zbidx_temp[i]=np.where(WCT==zb_temp[i])[0][0]
            else:
                zbidx_temp[i]=999
        zt=np.max(ztidx_temp) 
        # if np.sum(zbidx_temp==1)>=0:
        #     zb=np.sort(zbidx_temp)[-2]
        # else:
        #     zb=np.min(zbidx_temp)
        zb=np.min(zbidx_temp)
    except:
        zt=999
        zb=999
    if zt>=endheight:
        zt=999
    if zb<=startheightmax:
        zb=999
    if zb>=zt:
        zb=999
        zt=999 
    return zb,zt

# WCT=WCT[0,:]
# z=z_[:]

# startheightmax=1
# startheightmin=0
# endheight=np.where(z<18000)[0][-1]
# zmax=np.where(WCT==np.max(WCT[startheightmax:endheight]))[0]
# zmin=np.where(WCT==np.min(WCT[startheightmin:endheight]))[0]
# # zt=find_nearest(WCT[zmax[0]:endheight],Ws)+zmax[0]
# # zb=find_nearest(WCT[startheight:zmin[0]],-1*Ws)+startheight
# WCT_rangeup=np.sort(WCT[zmax[0]:endheight])
# WCT_rangedown=np.sort(WCT[startheightmin:zmin[0]])

# Wsmax=np.max(WCT[startheightmax:endheight])
# Wsmin=np.min(WCT[startheightmax:endheight])

# Ws1=np.max(WCT[startheightmax:endheight])-Wsmax/6.0
    
#     # Ws2=WCT[zmin]+10.0
# Ws2=np.min(WCT[startheightmin:endheight])+Wsmin/-2.0
# try:
#     zt_temp=find_closet_ele(WCT_rangeup, Ws1, 4)
#     zb_temp=find_closet_ele(WCT_rangedown,Ws2, 4)
#     ztidx_temp=np.empty(4,dtype=int)
#     zbidx_temp=np.empty(4,dtype=int)
#     for i in range(4):
#         if zt_temp[i] !=np.inf:
#             ztidx_temp[i]=np.where(WCT==zt_temp[i])[0][0]
#         else:
#             ztidx_temp[i]=999
#         if zb_temp[i] !=np.inf:
#             zbidx_temp[i]=np.where(WCT==zb_temp[i])[0][0]
#         else:
#             zbidx_temp[i]=999
#     zt=np.max(ztidx_temp) 
#     # if np.sum(zbidx_temp==1)>=0:
#     #     zb=np.sort(zbidx_temp)[-2]
#     # else:
#     #     zb=np.min(zbidx_temp)
#     zb=np.min(zbidx_temp)
# except:
#     zt=999
#     zb=999
# if zt>=endheight:
#     zt=999
# if zb<=startheightmax:
#     zb=999
# if zb>=zt:
#     zb=999
#     zt=999



def find_closet_ele(lst, target, k):
    closet_ele_lst = []
    closest_index = find_closest_index(lst, target)   # 距离target最近的元素的位置
    closet_ele_lst.append(lst[closest_index])

    ele_count = 1
    left_index = closest_index - 1
    right_index = closest_index + 1
    # 借鉴归并排序思路,向两侧滑动遍历
    while ele_count < k:
        left_ele = get_ele(lst, left_index)
        right_ele = get_ele(lst, right_index)
        # 哪边元素距离target更近,哪边就走一步
        if target - left_ele <= right_ele - target:
            closet_ele_lst.append(left_ele)
            left_index -= 1
        else:
            closet_ele_lst.append(right_ele)
            right_index += 1
        ele_count += 1

    closet_ele_lst.sort()
    return closet_ele_lst


def get_ele(lst, index):
    # 索引超出范围的,返回无穷大
    if index < 0 or index >= len(lst):
        return float("inf")
    return lst[index]

def find_closest_index(lst, target):
    """
    寻找距离target最近的位置
    :param lst:
    :param target:
    :return:
    """
    eg_index = find_eg_index(lst, 0, len(lst)-1, target)
    if eg_index == len(lst):
        return eg_index - 1
    elif eg_index == 0:
        return 0
    else:
        if target - lst[eg_index-1] <= lst[eg_index] - target:
            return eg_index-1
        else:
            return eg_index

def find_eg_index(lst, start, end, target):
    """
    找到第一个大于等于target的元素的位置, 如果lst中最大的元素小于target
    则返回len(lst)
    :param lst:
    :param start:
    :param end:
    :param target:
    :return:
    """
    if start > end:
        return end + 1
    middle = (start + end)//2
    middle_value = lst[middle]
    if middle_value == target:
        return middle
    elif middle_value < target:
        return find_eg_index(lst, middle+1, end, target)
    else:
        return find_eg_index(lst, start, middle-1, target)
    

fn = 'D:\\Lidar\\palau\\ComCAL\\auswertung\\netcdf\\2022\\ComCAL220805_2.nc'
data = nc.Dataset(fn, 'r')
time = data.variables['Time'][:].data
Height = data.variables['Height'][:].data

timesize=time.shape[1]
Hindex=np.where((Height>=10000) & (Height<=19000))
#Hindex=np.where((Height>=10000))

hourtime=np.empty(timesize,dtype=int)
minutetime=np.empty(timesize,dtype=float)
for i in range(timesize):
    hourtime[i]=int(time[12:13,i])*10+int(time[13:14,i])
    minutetime[i]=int(time[12:13,i])*10+int(time[13:14,i])+(int(time[15:16,i])*10+int(time[16:17,i]))/60

#tindex=np.where((hourtime>=11)&(hourtime<=21) )
tindex=np.where((hourtime>=0)&(hourtime<=24) )[0]
# tindex=tindex[:-4]
# tindex=np.where((hourtime<=24))

hourtime=hourtime[tindex]
minutetime=minutetime[tindex]
BSR532tot=data.variables['BSR532tot'][:].data[tindex,:][:,Hindex[0]]
BetaAer532tot=data.variables['BetaAer532tot'][:].data[tindex,:][:,Hindex[0]]
BetaAer355tot=data.variables['BetaAer355tot'][:].data[tindex,:][:,Hindex[0]]
BetaAer532S=data.variables['BetaAer532S'][:].data[tindex,:][:,Hindex[0]]
LR532arr=data.variables['LR532arr'][:].data[tindex,:][:,Hindex[0]]
AlphaAer532S=data.variables['AlphaAer532S'][:].data[tindex,:][:,Hindex[0]]
AlphaAer532P=data.variables['AlphaAer532P'][:].data[tindex,:][:,Hindex[0]]
BetaAer532toterr=data.variables['BetaAer532toterr'][:].data[tindex,:][:,Hindex[0]]

filterBSRidx=0
for i in range(np.shape(BSR532tot)[0]):
    filterBSR=np.max(BSR532tot[i,:])
    if filterBSR<=1000:
        filterBSRidx=np.append(filterBSRidx,i)


filterBSRidx=filterBSRidx[1:]
time=minutetime[filterBSRidx]
BSR532tot=BSR532tot[filterBSRidx,:]
BetaAer532tot=BetaAer532tot[filterBSRidx]
BetaAer355tot=BetaAer355tot[filterBSRidx]
BetaAer532S=BetaAer532S[filterBSRidx]
LR532arr=LR532arr[filterBSRidx]
AlphaAer532S=AlphaAer532S[filterBSRidx]
AlphaAer532P=AlphaAer532P[filterBSRidx]
Alpha532tot=BetaAer532tot*LR532arr
BetaAer532toterr=BetaAer532toterr[filterBSRidx]

z_=Height[Hindex[0]]
heightsize=z_.shape[0]
timesize=filterBSRidx.shape[0]

print(BSR532tot.shape)
print(z_.shape)

cmap = plt.cm.get_cmap('cet_coolwarm')  
plt.figure()
timeforplot=np.arange(0,BSR532tot.shape[0],1)
plt.contourf(time,z_[:],BSR532tot[:,:].T,cmap=cmap,levels=50,vmin=1,vmax=40)
plt.colorbar()


# plt.figure()
# plt.plot(BetaAer532S[35,:]/(BetaAer532tot[35,:]-BetaAer532S[35,:]))

# plt.figure()
# plt.plot(BetaAer355tot[35,:]/(BetaAer532tot[35,:]))

a=1000

WCT=np.zeros((timesize,heightsize))
for i in range(timesize):
    WCT[i,:]=WCT_calculation(a,z_,BSR532tot[i,:])





#first layer
zbori=np.empty(timesize,dtype=int)
ztori=np.empty(timesize,dtype=int)
for i in range (timesize):
    (zbori[i],ztori[i])=findzbandzt (WCT[i,:],z_)


zb=zbori[(zbori!=999)&(ztori!=999)]
zt=ztori[(zbori!=999)&(ztori!=999)]
time=time[(zbori!=999)&(ztori!=999)]
Alpha532tot=Alpha532tot[(zbori!=999)&(ztori!=999)]
BSR532tot2=BSR532tot[(zbori!=999)&(ztori!=999)]
BetaAer532toterr=BetaAer532toterr[(zbori!=999)&(ztori!=999)]
BetaAer532tot=BetaAer532tot[(zbori!=999)&(ztori!=999)]
cldbase=z_[zb]
cldtop=z_[zt]
actsize=np.shape(zb)[0]

fig, ax = plt.subplots(figsize=(3, 3))
plt.plot(z_/1000,BSR532tot[10,:])
plt.scatter(z_[zbori[10]]/1000,BSR532tot[10,zbori[10]])
plt.scatter(z_[ztori[10]]/1000,BSR532tot[10,ztori[10]])
plt.xlabel('Height (km)')

plt.ylabel('BSR@532nm')
# fig.savefig('D://Lidar//plot//BSR221206BSR.png', dpi=300, bbox_inches='tight')
fig, ax = plt.subplots(figsize=(3, 3))
plt.scatter(z_[:]/1000,WCT[10,:])
plt.scatter(z_[zbori[10]]/1000,WCT[10,zbori[10]])
plt.scatter(z_[ztori[10]]/1000,WCT[10,ztori[10]])
plt.xlabel('Height (km)')
plt.ylabel('WCT value')
# fig.savefig('D://Lidar//plot//WCT221206BSR.png', dpi=300, bbox_inches='tight')


# plt.figure()
# Ws=150
# plt.plot(z_,WCT[108,:])
# plt.scatter(z_[zbori[108]],-Ws+50)
# plt.scatter(z_[ztori[108]],Ws)
# plt.xlabel('Height (m)')
# plt.ylabel('WCT')
# # plt.savefig('D://Lidar//plot//2022-08//BSR0805WCT.png',dpi=300,bbox_inches='tight')

#calculate tau
taucld532=np.zeros(actsize,dtype=float)
Betaerr532=[]
BetaAercld532=[]
for i in range(actsize):
    wherecld=z_[zb[i]:zt[i]]
    clddepth=np.shape(wherecld)[0]
    for j in range(clddepth-1):
        temptau=Alpha532tot[i,(zb[i]+j)]*(wherecld[j+1]-wherecld[j])
        # print(temptau)
        taucld532[i]=taucld532[i]+temptau
    Betaerr532.append(BetaAer532toterr[i,zb[i]:zt[i]])
    BetaAercld532.append(BetaAer532tot[i,zb[i]:zt[i]])
    
        # print(taucld532[i])

for i in range(actsize):
    errsize=len(Betaerr532[i])
    # xxline=np.ones(errsize)*taucld532[i]
    # plt.scatter(xxline,Betaerr532[i],color='black')
    # plt.scatter(xxline,BetaAercld532[i],color='red')
    
    plt.scatter(taucld532[i],np.std(BetaAercld532[i]),color='blue')
      
plt.figure()
time1=np.arange(0,actsize,1)
plt.pcolormesh(time[:],z_[:]/1000,BSR532tot2[:,:].T/1.15,cmap=cmap,vmin=1.0)
cb=plt.colorbar()
# plt.scatter(time,z_[zb[:]],s=10,color='blue')
# plt.scatter(time,z_[zt[:]],s=10,color='orange')
plt.ylim(10,19)
plt.xlabel('Time UTC-hour')
plt.ylabel('Height (m)')
cb.set_label('Backscatter ratio @532 nm')
# plt.savefig('D://Lidar//plot//BSR1206.png',dpi=300,bbox_inches='tight')


plt.figure()
plt.contourf(time1,z_,Alpha532tot.T,levels=50,cmap=cc.cm.CET_L20)


plt.figure()
plt.scatter(time,taucld532)
plt.xlabel('Time UTC-hour')
# plt.savefig('D://Lidar//plot//2022-08//OD0805',dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(z_[zb])
plt.plot(z_[zt[:]])

cldbase=z_[zb]
cldtop=z_[zt]
# savefile
# fn2='D:/Lidar/palau/cloudlayer/'
# np.savetxt(fn2+'cloudbase'+fn[51:57]+'_3.txt',cldbase[:-1])
# np.savetxt(fn2+'cloudtop'+fn[51:57]+'_3.txt',cldtop[:-1])
# np.savetxt(fn2+'cloudtau532'+fn[51:57]+'_3.txt',taucld532[:-1])
# np.savetxt(fn2+'time'+fn[51:57]+'_3.txt',time[:-1])

# fn2='D:/Lidar/palau/cloudlayer/'
# np.savetxt(fn2+'cloudbase'+fn[51:57]+'.txt',cldbase[:])
# np.savetxt(fn2+'cloudtop'+fn[51:57]+'.txt',cldtop[:])
# np.savetxt(fn2+'cloudtau532'+fn[51:57]+'.txt',taucld532[:])
# np.savetxt(fn2+'time'+fn[51:57]+'.txt',time[:])


#second cloud # # # #
# tsdidx=0
# cldly2=67
# WCTsd=WCT[tsdidx:,cldly2:]
# Alpha532totsd=Alpha532tot[tsdidx:,cldly2:]
# BSR532totsd=BSR532tot[tsdidx:,cldly2:]
# timesd=time[tsdidx:]
# z=z_[cldly2:]
# tsdsize=np.shape(WCTsd)[0]
# zbori=np.empty(tsdsize,dtype=int)
# ztori=np.empty(tsdsize,dtype=int)
# for i in range (tsdsize):
#     (zbori[i],ztori[i])=findzbandzt (WCTsd[i,:],z)


# zb=zbori[(zbori!=999)&(ztori!=999)]
# zt=ztori[(zbori!=999)&(ztori!=999)]
# timesd=timesd[(zbori!=999)&(ztori!=999)]
# Alpha532totsd=Alpha532totsd[(zbori!=999)&(ztori!=999)]
# BSR532totsd=BSR532totsd[(zbori!=999)&(ztori!=999)]
# cldbase=z[zb]
# cldtop=z[zt]

# #calculate tau
# actsize=np.shape(zb)[0]
# taucld532=np.zeros(actsize,dtype=float)

# for i in range(actsize):
#     wherecld=z[zb[i]:zt[i]]
#     clddepth=np.shape(wherecld)[0]
#     for j in range(clddepth-1):
#         temptau=Alpha532totsd[i,(zb[i]+j)]*(wherecld[j+1]-wherecld[j])
#         # print(temptau)
#         taucld532[i]=taucld532[i]+temptau
#         print(taucld532[i])
        
# # # # # savefile2
# # fn2='D:/Lidar/palau/cloudlayer/'
# # np.savetxt(fn2+'cloudbase'+fn[51:57]+'_2.txt',cldbase)
# # np.savetxt(fn2+'cloudtop'+fn[51:57]+'_2.txt',cldtop)
# # np.savetxt(fn2+'cloudtau532'+fn[51:57]+'_2.txt',taucld532)
# # np.savetxt(fn2+'time'+fn[51:57]+'_2.txt',timesd)


# # # # # plt.figure()
# # # # # plt.plot(WCTsd[130,:])

# plt.figure()
# time1=np.arange(0,actsize,1)
# plt.contourf(timesd,z,BSR532totsd.T,levels=50,cmap=cc.cm.CET_L20)
# plt.colorbar()
# plt.scatter(timesd,z[zb],s=10,color='red')
# plt.scatter(timesd,z[zt],s=10,color='darkblue')
# plt.xlabel('Time UTC-hour')
# plt.ylabel('Height (m)')
# # # plt.savefig('D://Lidar//plot//BSR1213.png',dpi=300,bbox_inches='tight')

# plt.figure()
# plt.scatter(timesd,taucld532)
# plt.ylabel('Cloud OD')
# plt.xlabel('Time UTC-hour')
# plt.savefig('D://Lidar//plot//2022-08//OD0826.png',dpi=300,bbox_inches='tight')


