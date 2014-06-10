"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module computes sunrise and sunset times, start and end of
twilight, and the length of the day at any date and latitude. It is
based on code written by Paul Schlyter available at
http://kortis.to/radix/python/code/Sun.py

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-01 18:13 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from numpy import (pi, deg2rad, rad2deg, sin, cos, tan, arcsin, arccos,
    arctan, arctan2, floor, sqrt, array)

import calendar

class Sun:
    
    def __init__(self):
        # Some conversion factors between radians and degrees
        self.INV360 = 1.0 / 360.0

    
    def daysSince2000Jan0(self, y, m, d):
        """A macro to compute the number of days elapsed since
        2000 Jan 0.0 (which is equal to 1999 Dec 31, 0h UT)

        """
        return (367*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530)
    

    # The trigonometric functions in degrees
    def sind(self, x):
        """Returns the sine in degrees."""
        return sin(deg2rad(x))


    def cosd(self, x):
        """Returns the cosine in degrees."""
        return cos(deg2rad(x))


    def tand(self, x):
        """Returns the tangent in degrees."""
        return tan(deg2rad(x))


    def atand(self, x):
        """Returns the arctangent in degrees."""
        return rad2deg(arctan(x))

    
    def asind(self, x):
        """Returns the arcsine in degrees."""
        return rad2deg(arcsin(x))


    def acosd(self, x):
        """Returns the arccosine in degrees"""
        return rad2deg(arccos(x))


    def atan2d(self, y, x):
        """Returns the arctangent2 in degrees."""
        return rad2deg(arctan2(y, x))

    # Some lambda functions for setting longitude limits
    lon_n = lambda self, x, n: x + (x <= (n - 360)) * 360 - (x >= n) * 360
    lon180 = lambda self, x: self.lon_n(x, 180)
    lon360 = lambda self, x: self.lon_n(x, 360)


    # Following are some macros around the "workhorse" function __daylen__ 
    # They mainly fill in the desired values for the reference altitude    
    # below the horizon, and also selects whether this altitude should     
    # refer to the Sun's center or its upper limb.                         
    def dayLength(self, year, month, day, lon, lat):
        """This macro computes the length of the day, from sunrise to
        sunset. Sunrise/set is considered to occur when the Sun's upper
        limb is 35 arc minutes below the horizon (this accounts for the
        refraction of the Earth's atmosphere).
        
        """
        return self.__daylen__(year, month, day, lon, lat, -35.0/60.0, 1)

    
    def dayCivilTwilightLength(self, year, month, day, lon, lat):
        """This macro computes the length of the day, including civil
        twilight. Civil twilight starts/ends when the Sun's center is 6
        degrees below the horizon.
        
        """
        return self.__daylen__(year, month, day, lon, lat, -6.0, 0)


    def dayNauticalTwilightLength(self, year, month, day, lon, lat):
        """This macro computes the length of the day, incl. nautical
        twilight. Nautical twilight starts/ends when the Sun's center is
        12 degrees below the horizon.
        
        """
        return self.__daylen__(year, month, day, lon, lat, -12.0, 0)


    def dayAstronomicalTwilightLength(self, year, month, day, lon, lat):
        """This macro computes the length of the day, incl. astronomical
        twilight. Astronomical twilight starts/ends when the Sun's
        center is 18 degrees below the horizon.
        
        """
        return self.__daylen__(year, month, day, lon, lat, -18.0, 0)

    
    def sunRiseSet(self, year, month, day, lon, lat):
        """This macro computes times for sunrise/sunset. Sunrise/set is
        considered to occur when the Sun's upper limb is 35 arc minutes
        below the horizon (this accounts for the refraction of the
        Earth's atmosphere).

        """
        return self.__sunriset__(year, month, day, lon, lat, -35.0/60.0, 1)

    
    def civilTwilight(self, year, month, day, lon, lat):
        """This macro computes the start and end times of civil
        twilight. Civil twilight starts/ends when the Sun's center is 6
        degrees below the horizon.

        """
        return self.__sunriset__(year, month, day, lon, lat, -6.0, 0)

    
    def nauticalTwilight(self, year, month, day, lon, lat):
        """This macro computes the start and end times of nautical
        twilight. Nautical twilight starts/ends when the Sun's center is
        12 degrees below the horizon.

        """
        return self.__sunriset__(year, month, day, lon, lat, -12.0, 0)

    
    def astronomicalTwilight(self, year, month, day, lon, lat):
        """This macro computes the start and end times of astronomical
        twilight. Astronomical twilight starts/ends when the Sun's
        center is 18 degrees below the horizon.

        """
        return self.__sunriset__(year, month, day, lon, lat, -18.0, 0)

    
    # The "workhorse" function for sun rise/set times
    def __sunriset__(self, year, month, day, lon, lat, altit, upper_limb):
        """Note: year, month, date = calendar date, 1801-2099 only.

        Eastern longitude positive, western longitude negative, northern
        latitude positive, southern latitude negative.
        
        The longitude value IS critical in this function!

        altit = the altitude which the Sun should cross. Set to -35/60
        degrees for rise/set, -6 degrees for civil, -12 degrees for
        nautical and -18 degrees for astronomical twilight.

        upper_limb: non-zero -> upper limb, zero -> center.  Set to
        non-zero (e.g. 1) when computing rise/set times, and to zero
        when computing start/end of twilight.

        *rise = where to store the rise time.

        *set  = where to store the set time.  Both times are relative to
        the specified altitude, and thus this function can be used to
        compute various twilight times, as well as rise/set times.
        
        RETURNS
            0 = sun rises/sets this day, times stored at *trise and
            *tset.

            +1 = sun above the specified 'horizon' 24 hours. *trise set
            to time when the sun is at south, minus 12 hours while *tset
            is set to the south time plus 12 hours. 'Day' length = 24
            hours

            -1 = sun is below the specified 'horizon' 24 hours 'Day'
            length = 0 hours, *trise and *tset are both set to the time
            when the sun is at south.
        
        """
        # Compute d of 12h local mean solar time
        d = self.daysSince2000Jan0(year, month, day) + 0.5 - (lon / 360.0)
        
        # Compute local sidereal time of this moment 
        sidtime = self.revolution(self.GMST0(d) + 180.0 + lon)
        
        # Compute Sun's RA + Decl at this moment 
        res = self.sunRADec(d)
        sRA = res[0]
        sdec = res[1]
        sr = res[2]
        
        # Compute time when Sun is at south - in hours UT 
        tsouth = 12.0 - self.rev180(sidtime - sRA) / 15.0;
        
        # Compute the Sun's apparent radius, degrees 
        sradius = 0.2666 / sr;
        
        # Do correction to upper limb, if necessary 
        if upper_limb:
            altit = altit - sradius
        
        # Compute the diurnal arc that the Sun traverses to reach 
        # the specified altitude altit: 
        cost = ((self.sind(altit) - self.sind(lat) * self.sind(sdec)) /
            (self.cosd(lat) * self.cosd(sdec)))

        # Makes sure that when sun is below altit (cost >= 1.0),
        # when sun is above altit (cost <= -1.0) and otherwise, the length
        # is properly calculated.
        t = ((cost >= 1.0) * 0 + (cost <= -1.0) * 12 +
            ((cost > -1.0) & ( cost < 1.0)) * self.acosd(cost) / 15.0)
        rc = ((cost >= 1.0) * (-1) + (cost <= -1.0) * 1 +
            ((cost > -1.0) & ( cost < 1.0)) * 0)
        # THIS IS THE OLD CODE:
        #if cost >= 1.0:
        #    rc = -1 
        #    t = 0.0         # Sun always below altit
        #elif cost <= -1.0:
        #    rc = +1 
        #    t = 12.0;       # Sun always above altit
        #else:
        #    t = self.acosd(cost)/15.0   # The diurnal arc, hours
        
        # Store rise and set times - in hours UT 
        return (tsouth-t, tsouth+t)


    def __daylen__(self, year, month, day, lon, lat, altit, upper_limb):
        """Note: year,month,date = calendar date, 1801-2099 only.

        Eastern longitude positive, western longitude negative, northern
        latitude positive, southern latitude negative. The longitude
        value is not critical. Set it to the correct longitude if you're
        picky, otherwise set to, say, 0.0.

        The latitude however IS critical - be sure to get it correct.

        altit = the altitude which the Sun should cross. Set to -35/60
        degrees for rise/set, -6 degrees for civil, -12 degrees for
        nautical and -18 degrees for astronomical twilight.

        upper_limb: non-zero -> upper limb, zero -> center. Set to
        non-zero (e.g. 1) when computing day length and to zero when
        computing day+twilight length.

        """
        # Compute d of 12h local mean solar time 
        d = self.daysSince2000Jan0(year, month, day) + 0.5 - (lon / 360.0)
                 
        # Compute obliquity of ecliptic (inclination of Earth's axis) 
        obl_ecl = 23.4393 - 3.563E-7 * d
        
        # Compute Sun's position 
        res = self.sunpos(d)
        slon = res[0]
        sr = res[1]
        
        # Compute sine and cosine of Sun's declination 
        sin_sdecl = self.sind(obl_ecl) * self.sind(slon)
        cos_sdecl = sqrt(1.0 - sin_sdecl * sin_sdecl)
        
        # Compute the Sun's apparent radius, degrees 
        sradius = 0.2666 / sr
        
        # Do correction to upper limb, if necessary 
        if upper_limb:
            altit = altit - sradius
        
        cost = ((self.sind(altit) - self.sind(lat) * sin_sdecl) /
            (self.cosd(lat) * cos_sdecl))

        # Makes sure that when sun is below altit (cost >= 1.0),
        # when sun is above altit (cost <= -1.0) and otherwise, the length
        # is properly calculated.
        mask = (cost > -1.0) & ( cost < 1.0)
        mask_below = (cost >= 1.0)
        mask_above = (cost <= -1.0)
        t = (mask_below * 0 + mask_above * 24 +
            mask * (2.0/15.0) * self.acosd(mask * cost))
        # THIS IS THE OLD CODE:
        #if cost >= 1.0:
        #    t = 0.0         # Sun always below altit
        #elif cost <= -1.0:
        #    t = 24.0        # Sun always above altit
        #else:
        #    t = (2.0/15.0) * self.acosd(cost);     # The diurnal arc, hours
        
        return t


    def sunpos(self, d):
        """
        Computes the Sun's ecliptic longitude and distance 
        at an instant given in d, number of days since     
        2000 Jan 0.0.  The Sun's ecliptic latitude is not  
        computed, since it's always very near 0.           
        """

        # Compute mean elements 
        M = self.revolution(356.0470 + 0.9856002585 * d)
        w = 282.9404 + 4.70935E-5 * d
        e = 0.016709 - 1.151E-9 * d
        
        # Compute true longitude and radius vector 
        E = M + rad2deg(e * self.sind(M) * (1.0 + e * self.cosd(M)))
        x = self.cosd(E) - e
        y = sqrt(1.0 - e*e) * self.sind(E)
        r = sqrt(x*x + y*y)                #Solar distance 
        v = self.atan2d(y, x)                   # True anomaly 
        lon = v + w                             # True solar longitude
        lon = self.lon360(lon)                  # Make it 0..360 degrees
            
        return (lon, r)
    

    def sunRADec(self, d):
        """Returns the angle of the Sun (RA) the declination (dec)
        and the distance of the Sun (r) for a given day d.
        
        """
        # Compute Sun's ecliptical coordinates
        res = self.sunpos(d)
        lon = res[0]        # True solar longitude
        r = res[1]          # Solar distance
        
        # Compute ecliptic rectangular coordinates (z=0) 
        x = r * self.cosd(lon)
        y = r * self.sind(lon)
        
        # Compute obliquity of ecliptic (inclination of Earth's axis) 
        obl_ecl = 23.4393 - 3.563E-7 * d
        
        # Convert to equatorial rectangular coordinates - x is unchanged 
        z = y * self.sind(obl_ecl)
        y = y * self.cosd(obl_ecl)

        # Convert to spherical coordinates 
        RA = self.atan2d(y, x)
        dec = self.atan2d(z, sqrt(x*x + y*y))

        return (RA, dec, r)


    def revolution(self, x):
        """This function reduces any angle to within the first
        revolution by subtracting or adding even multiples of 360.0
        until the result is >= 0.0 and < 360.0.
        
        Reduce angle to within 0..360 degrees.
        
        """
        return (x - 360.0 * floor(x * self.INV360))


    def rev180(self, x):
        """Reduce angle to within +180..+180 degrees.""" 
        return (x - 360.0 * floor(x * self.INV360 + 0.5))


    def GMST0(self, d):
        """This function computes GMST0, the Greenwich Mean Sidereal
        Time at 0h UT (i.e. the sidereal time at the Greenwhich meridian
        at 0h UT).  GMST is then the sidereal time at Greenwich at any
        time of the day.  I've generalized GMST0 as well, and define it
        as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at
        other times than 0h UT as well.  While this sounds somewhat
        contradictory, it is very practical:  instead of computing GMST
        like:

            GMST = (GMST0) + UT * (366.2422/365.2422),
        
        where (GMST0) is the GMST last time UT was 0 hours, one simply
        computes:
        
            GMST = GMST0 + UT,
        
        where GMST0 is the GMST "at 0h UT" but at the current moment!
        Defined in this way, GMST0 will increase with about 4 min a day.
        It also happens that GMST0 (in degrees, 1 hr = 15 degr) is equal
        to the Sun's mean longitude plus/minus 180 degrees! (if we
        neglect aberration, which amounts to 20 seconds of arc or 1.33
        seconds of time)
        
        """
        # Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  
        # L = M + w, as defined in sunpos().  Since I'm too lazy to 
        # add these numbers, I'll let the C compiler do it for me.  
        # Any decent C compiler will add the constants at compile   
        # time, imposing no runtime or code overhead.               
                                                
        sidtim0 = self.revolution((180.0 + 356.0470 + 282.9404) +
            (0.9856002585 + 4.70935E-5) * d)
        return sidtim0;


    def solar_altitude(self, latitude, year, month, day):
        """Compute the altitude of the sun. No atmospherical refraction
        taken in account. Altitude of the southern hemisphere are given
        relative to true north. Altitude of the northern hemisphere are
        given relative to true south. Declination is between 23.5N and
        23.5S depending on the period of the year. Source of formula
        for altitude is PhysicalGeography.net available at
        http://www.physicalgeography.net/fundamentals/6h.html
        
        """
        # Compute declination
        N = self.daysSince2000Jan0(year, month, day)
        res =  self.sunRADec(N)
        declination = res[1]
        sr = res[2]

        # Compute the altitude
        altitude = 90.0 - latitude  + declination

        # In the tropical and  in extreme latitude, values over 90 may occurs.
        if altitude > 90:
            altitude = 90 - (altitude-90)

        if altitude < 0:
            altitude = 0

        return altitude


    def get_max_solar_flux(self, latitude, year, month, day):
        """Compute the maximal solar flux to reach the ground for this
        date and latitude. Originaly comes from Environment Canada
        weather forecast model. Information was of the public domain
        before release by Environment Canada Output is in W M-2.

        """

        (fEot, fR0r, tDeclsc) = self.equation_of_time(year, month, day,
            latitude)
        fSF = (tDeclsc[0] + tDeclsc[1]) * fR0r

        # In the case of a negative declinaison, solar flux is null
        fCoeff = ((fSF < 0) * 0 +
            (fSF >= 0) * (-1.56e-12 * fSF**4 + 5.972e-9 * fSF**3 -
            8.364e-6 * fSF**2 + 5.183e-3 * fSF - 0.435))

        fSFT = fSF * fCoeff
        fSFT = (fSFT < 0) * 0 + (fSFT >= 0) * fSFT
        return fSFT


    def equation_of_time(self, year, month, day, latitude):
        """Subroutine computing the part of the equation of time needed
        in the computing of the theoritical solar flux correction
        originating of the CMC GEM model.

        PARAMETERS
            int nTime : cTime for the correction of the time.

        RETURNS
            double dEot: Correction for the equation of time 
            double dR0r: Corrected solar constant for the equation of time
            tuple tDeclsc: Declinaison
        
        """
        # Julian date 
        nJulianDate = self.Julian(year, month, day)
        # Check if it is a leap year
        if(calendar.isleap(year)):
            fDivide = 366.0
        else:
            fDivide = 365.0
        # Correction for "equation of time"
        fA = nJulianDate / fDivide * 2 * pi
        fR0r = self.__Solcons(fA) * 0.1367e4
        fRdecl = 0.412 * cos((nJulianDate + 10.0) * 2.0 * pi / fDivide - pi)
        fDeclsc1 = self.sind(latitude) * sin(fRdecl)
        fDeclsc2 = self.cosd(latitude) * cos(fRdecl)
        tDeclsc = (fDeclsc1, fDeclsc2)
        # in minutes
        fEot = (0.002733 - 7.343 * sin(fA)+ .5519 * cos(fA) -
            9.47 * sin(2.0 * fA) - 3.02 * cos(2.0*fA) - 0.3289 * sin(3. * fA) -
            0.07581 * cos(3.0 * fA) -0.1935 * sin(4.0 * fA) -
            0.1245 * cos(4.0 * fA))
        # Express in fraction of hour
        fEot = fEot / 60.0
        # Express in radians
        fEot = fEot * 15 * pi / 180.0

        return (fEot, fR0r, tDeclsc)


    def __Solcons(self, dAlf):
        """Statement function that calculates the variation of the solar
        constant as a function of the julian day. (dAlf, in radians).
        
        PARAMETERS
            double dAlf : Solar constant to correct the excentricity.

        RETURNS
            double dVar : Variation of the solar constant
        
        """
        
        dVar = 1.0 / (1.0 - 9.464e-4 * sin(dAlf) - 0.01671 * cos(dAlf) -
            1.489e-4 * cos(2.0 * dAlf) - 2.917e-5 * sin(3.0 * dAlf) -
            3.438e-4 * cos(4.0 * dAlf))**2
        return dVar


    def Julian(self, year, month, day):
        """Return julian day."""
        if calendar.isleap(year): # Bissextil year, 366 days
            lMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335,
                366]
        else: # Normal year, 365 days
            lMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334,
                365]
        nJulian = lMonth[month - 1] + day
        return nJulian


if __name__ == "__main__":
    k = Sun()
    lat = array([46.2, 22.5, 62.15])
    lon = array([25.42, 20.4, 187.0])
    print k.get_max_solar_flux(lat, 2004, 01, 30)
    print k.sunRiseSet(2002, 3, 22, lon, lat)
