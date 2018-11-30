## This is program to project sky coordinate to a plane [2018 Nov 01]
## Introduced the shift of coordinate in the midpoint calculation, required due to the default return of arctan/arccos [2018 Nov 21] 
## pbera.phy@gmail.com
###########################################
import numpy as np
import pylab as pl

DEC = -1.0
RA = 2.0

if (abs(DEC)>0.5*np.pi):
	print '** abs(DEC) should be less than 0.5*np.pi !!'

delta_theta = 0.1
delta_phi = 0.1

## 1st Source :: S1
S1_theta = (0.5*np.pi - DEC) - 0.5*delta_theta
S1_phi = RA - 0.5*delta_phi

## 2nd Source :: S2
S2_theta = (0.5*np.pi - DEC) + 0.5*delta_theta
S2_phi = RA + 0.5*delta_phi

## Assumed distance between S1 and S2 on the projected plane
D_S1S2 = 2.0

## Certesian coordinate of the midpoint between the sorces S1 and S2 on the unit sphere
Xa = 0.5*(np.sin(S1_theta)*np.cos(S1_phi)+np.sin(S2_theta)*np.cos(S2_phi))
Ya = 0.5*(np.sin(S1_theta)*np.sin(S1_phi)+np.sin(S2_theta)*np.sin(S2_phi))
#Za = 0.5*(S1_theta/abs(S1_theta)*np.cos(S1_theta) + S2_theta/abs(S2_theta)*np.cos(S2_theta))

Za = 0.5*(np.cos(S1_theta) + np.cos(S2_theta))

## Angular coordinate of the mid point bteween S1 & S2, The projected plane will be tangent to the sphere of a specified radius
Sa_theta = np.arccos(Za/(Xa*Xa+Ya*Ya+Za*Za)**0.5)
Sa_phi = np.arctan(Ya/Xa)

if (Za<-0.0):
	Sa_theta = Sa_theta-np.pi

if (Ya>0.0 and Sa_phi<0.0): ## Shift mid point to the second quardent if required
	Sa_phi = np.pi + Sa_phi

if (Ya<0.0 and Sa_phi>0.0):
	Sa_phi = np.pi + Sa_phi

print 'mid point', Sa_theta, Sa_phi

## Distance of the sources S1 and S2
r0 = D_S1S2/((np.sin(S1_theta)*np.cos(S1_phi)-np.sin(S2_theta)*np.cos(S2_phi))**2 + (np.sin(S1_theta)*np.sin(S1_phi)-np.sin(S2_theta)*np.sin(S2_phi))**2 + (np.cos(S1_theta) - np.cos(S2_theta))**2)**0.5

## Subtended angle by the sources S1 and S2
gamma = np.arccos(np.sin(S1_theta)*np.cos(S1_phi)*np.sin(S2_theta)*np.cos(S2_phi) + np.sin(S1_theta)*np.sin(S1_phi)*np.sin(S2_theta)*np.sin(S2_phi) + np.cos(S1_theta)*np.cos(S2_theta))

## Projected distance of the plane
ra = r0*np.cos(0.5*gamma)
#ra = r0*(np.sin(S1_theta)*np.cos(S1_phi)*np.sin(Sa_theta)*np.cos(Sa_phi) + np.sin(S1_theta)*np.sin(S1_phi)*np.sin(Sa_theta)*np.sin(Sa_phi) + np.cos(S1_theta)*np.cos(Sa_theta))

def distance_rr(theta, phi):
	## A function to find out the distance of the projected point on the plane for a specified direction (theta_i, phi_i); outpur is distance 'r_i'
	rr = ra/(np.sin(Sa_theta)*np.cos(Sa_phi)*np.sin(theta)*np.cos(phi) + np.sin(Sa_theta)*np.sin(Sa_phi)*np.sin(theta)*np.sin(phi) + np.cos(Sa_theta)*np.cos(theta))
	return rr


## Vector components of the X-axis
XXi = r0*(np.sin(S1_theta)*np.cos(S1_phi)-np.sin(S2_theta)*np.cos(S2_phi))/D_S1S2
XXj = r0*(np.sin(S1_theta)*np.sin(S1_phi)-np.sin(S2_theta)*np.sin(S2_phi))/D_S1S2
XXk = r0*(np.cos(S1_theta) - np.cos(S2_theta))/D_S1S2
'''
## Alternative construction of X-axis
epsilon = 0.09
Sd1_theta = (1.0-epsilon)*Sa_theta + epsilon*S1_theta
Sd1_phi = (1.0-epsilon)*Sa_phi + epsilon*S1_phi
rrSd1 = distance_rr(Sd1_theta, Sd1_phi)

XXi = (rrSd1*np.sin(Sd1_theta)*np.cos(Sd1_phi)-ra*np.sin(Sa_theta)*np.cos(Sa_phi))
XXj = (rrSd1*np.sin(Sd1_theta)*np.sin(Sd1_phi)-ra*np.sin(Sa_theta)*np.sin(Sa_phi))
XXk = (rrSd1*np.cos(Sd1_theta) - ra*np.cos(Sa_theta))

arg_XX = (XXi*XXi + XXj*XXj + XXk*XXk)**0.5
XXi  = XXi/arg_XX
XXj = XXj/arg_XX
XXk = XXk/arg_XX
'''

## Vecotr components of Z-axis (in the vertical direction to the projected plane)
ZZi = np.sin(Sa_theta)*np.cos(Sa_phi)
ZZj = np.sin(Sa_theta)*np.sin(Sa_phi)
ZZk = np.cos(Sa_theta)

## Vector components of Y-axis, i.e. the other axis perpendicular to X-axis and in the projection plane
YYi = ZZj*XXk - ZZk*XXj
YYj = ZZk*XXi - ZZi*XXk
YYk = ZZi*XXj - ZZj*XXi

## Check the Vector comenents 
print 'XX', XXi, XXj, XXk
print 'XX.XX',XXi*XXi + XXj*XXj + XXk*XXk, 'XX.YY', XXi*YYi + XXj*YYj + XXk*YYk, 'YY.YY', YYi*YYi + YYj*YYj + YYk*YYk

def PlaneCo_from_SkyCo(theta, phi):
	## Function to find out the coordinate on the projected plane for a given sky coordinate
	## input : sky coordinate = (theta_i, phi_i)
	## output : projected coordinate = (XX_i, YY_i)
	#rr = ra/(np.sin(Sa_theta)*np.cos(Sa_phi)*np.sin(theta)*np.cos(phi) + np.sin(Sa_theta)*np.sin(Sa_phi)*np.sin(theta)*np.sin(phi) + np.cos(Sa_theta)*np.cos(theta))
	#rr = ra*(np.sin(Sa_theta)*np.cos(Sa_phi)*np.sin(Sa_theta)*np.cos(Sa_phi) + np.sin(Sa_theta)*np.sin(Sa_phi)*np.sin(Sa_theta)*np.sin(Sa_phi) + np.cos(Sa_theta)*np.cos(Sa_theta))/(np.sin(Sa_theta)*np.cos(Sa_phi)*np.sin(theta)*np.cos(phi) + np.sin(Sa_theta)*np.sin(Sa_phi)*np.sin(theta)*np.sin(phi) + np.cos(Sa_theta)*np.cos(theta))
	rr = distance_rr(theta, phi)
	
	OPi = rr*np.sin(theta)*np.cos(phi) - ra*np.sin(Sa_theta)*np.cos(Sa_phi)
	OPj = rr*np.sin(theta)*np.sin(phi) - ra*np.sin(Sa_theta)*np.sin(Sa_phi)
	OPk = rr*np.cos(theta) - ra*np.cos(Sa_theta)
	
	print 'rr = ',rr, 'OP',OPi, OPj, OPk
	
	return OPi*XXi+OPj*XXj+OPk*XXk, OPi*YYi+OPj*YYj+OPk*YYk
	

print 'GAMMA : ',gamma, 'r0 : ', r0, 'ra : ', ra

print PlaneCo_from_SkyCo(S1_theta, S1_phi)
print PlaneCo_from_SkyCo(S2_theta, S2_phi)
