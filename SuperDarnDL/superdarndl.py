
from davitpy.pydarn.sdio.sdDataTypes import sdDataPtr,gridData
from davitpy.pydarn.sdio.sdDataRead import * #May want to check for collisions from sdDataRead
#import davitpy.gme.ind as gmeind
from geospacepy import omnireader
import apexpy
import igrfpy
import datetime
import aacgmv2
import numpy as np

def _sdGridLoad(sTime,eTime,deltaT,hemi,fileType,src,fileName,custType,noCache,estd=False,TWidth=None):
	#Internal function to get SuperDARN grid data for desired times, and adjust time resolution if necessary
	
	# Read in SuperDARN grid data
	myPtr  = sdDataOpen(sTime,hemi=hemi,eTime=eTime,fileType=fileType,src=src,
						fileName=fileName,noCache=noCache)
	#LMK remove custType kwarg because it's not in sdDataRead.sdDataOpen anymore as of DavitPy 0.5
	sdList = sdDataReadAll(myPtr)
	
	#This code has been seen to create problems with new versions of DavitPy
	#as the ptr attribute is now 'private'
	#if not myPtr.ptr.closed:
	#    myPtr.ptr.close()

	dTnative = (sdList[0].eTime - sdList[0].sTime).seconds/60

	#LMK added new optional argument TWidth to allow adjacent steps to possibily overlap in time
	#(e.g. each entry has 4 minutes of data but each adjacent is only 2 minutes apart)
	if TWidth is None:
		TWidth = deltaT #Default behaviour

	if (deltaT > dTnative): #will need to concatentate multiple records
		#Create new list using longer time steps
		sdListNew = []

		nn = 0
		iT = sTime
		while iT <= eTime:
			sdI = gridData()
			sdI.sTime = iT
			sdI.eTime = (iT + dt.timedelta(minutes=TWidth))
			sdI.vector.mlat      = []
			sdI.vector.mlon      = []
			sdI.vector.velmedian = []
			sdI.vector.velsd     = []
			sdI.vector.kvect     = []
			sdI.vector.stid      = []

			# Concatenate entries within larger time window
			while (nn < len(sdList)) and (sdList[nn].eTime <= (iT + dt.timedelta(minutes=TWidth))):
				sdI.vector.mlat.extend(sdList[nn].vector.mlat)
				sdI.vector.mlon.extend(sdList[nn].vector.mlon)
				sdI.vector.velmedian.extend(sdList[nn].vector.velmedian)
				sdI.vector.velsd.extend(sdList[nn].vector.velsd)
				sdI.vector.kvect.extend(sdList[nn].vector.kvect)
				sdI.vector.stid.extend(sdList[nn].vector.stid)
				nn += 1

			sdListNew.append(sdI)
			iT += dt.timedelta(minutes=deltaT)

		sdList = sdListNew

	if estd: # Need to recalulated errors based on regional std dev values
		err_min = 50.
		for sdI in sdList:

			lat = np.asarray(sdI.vector.mlat)
			lon = np.asarray(sdI.vector.mlon)
			vel = np.asarray(sdI.vector.velmedian)
			kaz = np.asarray(sdI.vector.kvect)
			
			rid_arr = np.asarray(sdI.vector.stid)
			rid_set = np.unique(rid_arr)
			for rad in rid_set:
				q = np.where(rid_arr == rad)
				q = q[0]

				for i in range(len(q)):
					# Find neighbooring data
					qq = np.where(np.logical_and(abs(lat[q] - lat[q[i]]) <= 1,
												 abs(lon[q] - lon[q[i]]) <= 7.5))
					qq = qq[0]

					velqq = vel[q[qq]]
					kazqq = kaz[q[qq]]

					# Not enough to calc variance
					if (len(qq) < 5): 
						sdI.vector.velsd[q[i]] *= 2

					else:
						var0 = np.std(velqq) #Simple std dev
						
						# Try magnitude variance taking into account variance in directions
						(mvel,mkaz) = _merge_one(velqq,kazqq)
						if mvel is not None:
							var1 = np.std(velqq/np.cos(np.radians(kazqq-mkaz)))
							sdI.vector.velsd[q[i]] = min((var0,var1))
						else: 
							sdI.vector.velsd[q[i]] = var0

					#err_min < err < 1000
					sdI.vector.velsd[q[i]] = max((min((sdI.vector.velsd[q[i]],1000)),err_min))

	return sdList

def readDay(hemiSN,dt,
			deltaT=2,fileType='grdex',src=None,fileName=None,
			custType='grdex',noCache=False):
	"""Read one day and hemisphere of SuperDARN data
	
	(All of these other options are passed through just_sam to DavitPy)
	
	Parameters
	----------
        hemiSN: N or S
	dt : TYPE
	    Description
	deltaT : int, optional
	    Description
	fileType : str, optional
	    Description
	src : None, optional
	    Description
	fileName : None, optional
	    Description
	custType : str, optional
	    Description
	noCache : bool, optional
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""
	#import davitpy.models.aacgm as aacgm

	#def _sdGridLoad(sTime,eTime,deltaT,hemi,fileType,src,fileName,custType,noCache,estd=False,TWidth=None):
	
	#SAM uses different hemisphere naming convention
	if hemiSN is 'N':
		hemi = 'north' 
	elif hemiSN is 'S':
		hemi = 'south'

	#Load in data in 2 minute chunks
	sdList = _sdGridLoad(dt,
									dt+datetime.timedelta(days=1),
									deltaT,
									hemi,
									fileType,
									src,
									fileName,
									custType,
									noCache,
									estd=False,TWidth=None)

	#Split out data into numpy arrays
	#Each sdI is a griddata object from pydarn.sdio.sdDataTypes
	startdts,enddts = [],[]
	mlats,mlons,mltlons,vels,verrs,azms,rids = [],[],[],[],[],[],[]
	eloss,eerrs = [],[] #Electric fields will be compute inplace
	bigrfs = []

	for sdI in sdList:
		if sdI.vector.mlat is None or sdI.vector.mlon is None:
			print sdI.vector
			continue

		startdts.append(sdI.sTime)
		enddts.append(sdI.eTime)
		mlat = np.asarray(sdI.vector.mlat)
		mlats.append(mlat)
		#Do in-band MLT conversion
		mlon = np.asarray(sdI.vector.mlon)
		y,M,d = sdI.sTime.year,sdI.sTime.month,sdI.sTime.day
		h,m,s = sdI.sTime.hour,sdI.sTime.minute,sdI.sTime.second
		# mltDef = aacgm.mltFromYmdhms(y,M,d,h,m,s,0.0) * 15. 
		# mltDef is the rotation that needs to be applied, 
		# and lon is the AACGM longitude.
		# use modulo so new longitude is between 0 & 360
		# mlt_lon = np.mod((mlon + mltDef), 360.)
		
                mlt_lon=aacgmv2.convert_mlt(mlon, sdI.sTime, False) ## sishen added 
                mlons.append(mlon) #Actual mlon
		mltlons.append(mlt_lon) #MLT in degrees
		#Convert to electric field
		vel = np.asarray(sdI.vector.velmedian)
		verr = np.asarray(sdI.vector.velsd)
		vels.append(vel)
		verrs.append(verr)

		# Where do superDARN returns come from?
		approx_alt = 300. 
		
		#Could use either Apex-Python or AACGMv2 package
		#for magnetic to geocentric conversion
		gdlat,gdlon = aacgmv2.convert(mlat,mlon,approx_alt,date=sdI.sTime,a2g=True)
		#gdlat,gdlon = ac.apex2geo(mlat,mlon,approx_alt)
		Be,Bn,Bu = igrfpy.getmainfield(sdI.sTime,
										gdlat,gdlon,
										np.ones_like(gdlat)*approx_alt,
										geocentric=False,silent=True)
		Bu = np.abs(np.array(Bu)) # It returns lists

		# Convert SuperDARN Vlos to E-field
		# (Should be using IGRF...just using constant instead)
		# Will change after compare with SAM
		#Bconst = 0.5e-4
		
		bigrfs.append(Bu*1.0e-9)
		eloss.append(vel*Bu*1.0e-9)
		eerrs.append(verr*Bu*1.0e-9)

		azms.append(np.asarray(sdI.vector.kvect)) #degrees
		rids.append(np.asarray(sdI.vector.stid)) #station id num

	return (startdts,enddts,mlats,mlons,mltlons,
		vels,verrs,eloss,eerrs,azms,rids,bigrfs)



## just test
if __name__ == "__main__":
        dt=datetime.datetime(2017,2,1,0,0,0)
        (startdts, enddts, mlats, mlons, mltlons, vels, verrs, eloss, eerrs, azms, rids, bigrfs) = readDay("N", dt)
        print(len(startdts))
        print("0:"+str(startdts[0])+" -> "+str(enddts[0]))
        mid=len(startdts)//2
        print("mid:"+str(startdts[mid])+" -> "+str(enddts[mid]))
        print("end:"+str(startdts[-1])+" -> "+str(enddts[-1]))
        np.set_printoptions(threshold=50)
        len_each_mlat=np.zeros(len(mlats))
        len_each_vel=np.zeros(len(vels))
        len_each_elos=np.zeros(len(eloss))
        for i in range(len(mlats)):
            len_each_mlat[i]=len(mlats[i])
        for i in range(len(vels)):
            len_each_vel[i]=len(vels[i])
        for i in range(len(eloss)):
            len_each_elos[i]=len(eloss[i])
        print("mlat len:")
        print(len_each_mlat)
        print("vel len:")
        print(len_each_vel)
        print("elos len:")
        print(len_each_elos)


