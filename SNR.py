#!/usr/bin/env python
import math
from __future__ import print_function # Python 2 & 3 compatibility
from __future__ import division # Python 2 & 3 compatibility

   
# Sky surface brightness, Gunn g & r bands
sky_mu = {'g':22.0, 'r':21.0} # AB mag/arcsec^2

# SBIG STD-8300M properties
dark_current = 0.04 # e/pixel/second
read_noise = 10.0 # e/pixel/read

# Number of photons/second/pixel at the top of atmosphere 
# that corresponds to 0 AB mag/arcsec^2, Gunn g & r bands.
gamma0 = {'g':1.79e9, 'r':1.16e9}

# Fractional overall system efficiency, Gunn g & r bands
efficiency = {'g':0.34, 'r':0.35}


def snr(sci_mu, band, total_exp_time, sub_exp_time=600, binning=1, N=1, round_up=True):
    
    # Science count rates, uses equations 8 & 9 from Dragonfly paper to convert AB mag/arcsec^2 to e/pixel/s
    rate_sci = gamma0[band] * N * efficiency[band] * 10**(-0.4 * sci_mu) # e/pixel/second
    
    # Sky count rates, uses equations 8 & 9 from Dragonfly paper to convert AB mag/arcsec^2 to e/pixel/s
    rate_sky =  gamma0[band] * N * efficiency[band] * 10**(-0.4 * sky_mu[band]) # e/pixel/second
    
    # Number of sub-exposures
    number_subs = int(math.ceil(total_exp_time/sub_exp_time))
    if round_up and (total_exp_time != number_subs * sub_exp_time):
        total_exp_time = number_subs * sub_exp_time
        print('Rounding up total exposure time to next integer multiple of sub-exposure time:', total_exp_time)
    
    # Noise sources
    signal = rate_sci * total_exp_time # e/pixel
    sky_counts = rate_sky * total_exp_time # e/pixel
    dark_counts = dark_current * total_exp_time # e/pixel
    total_read_noise = math.sqrt(number_subs) * read_noise  # e/pixel
    
    noise = math.sqrt(signal + sky_counts + dark_counts + total_read_noise**2) # e/pixel
    
    s = binning * signal/noise # Binning increases signal to noise by the binning factor
    
    #print(s, signal, sky_counts, dark_counts, total_read_noise**2)
        
    return s


def etc(sci_mu, band, snr_target, sub_exp_time=600, binning=1, N=1):
    
    # Convert target SNR per binned pixel to SNR per unbinned pixel
    snr_target /= binning
    
    # Science count rates, uses equations 8 & 9 from Dragonfly paper to convert AB mag/arcsec^2 to e/pixel/s
    rate_sci = gamma0[band] * N * efficiency[band] * 10**(-0.4 * sci_mu) # e/pixel/second
    
    # Sky count rates, uses equations 8 & 9 from Dragonfly paper to convert AB mag/arcsec^2 to e/pixel/s
    rate_sky =  gamma0[band] * N * efficiency[band] * 10**(-0.4 * sky_mu[band]) # e/pixel/second
    
    # If required total exposure time is much greater than the length of a sub-exposure then
    # all noise sources (including read noise) are proportional to t^0.5 and we can use a 
    # simplified expression to estimate total exposure time.
    total_exp_time = snr_target**2 * (rate_sci + rate_sky + dark_current + read_noise**2/sub_exp_time) / rate_sci**2
    
    # The simplified expression underestimates read noise due to fractional number of sub-exposure,
    # the effect will be neglible unless the total exposure time is very short but we can fix it anyway...
    # First round up to the next integer number of sub-exposures:
    number_subs = int(math.ceil(total_exp_time / sub_exp_time))
    # If the SNR has dropped below the target value as a result of the extra read noise add another sub
    # Note: calling snr() here is horribly inefficient as it recalculates a bunch of stuff but I don't care.
    while snr(sci_mu, band, number_subs*sub_exp_time, sub_exp_time, binning, N, round_up=False) < snr_target:
        print("Adding a sub-exposure to overcome read noise!")
        number_subs += 1
    
    return number_subs*sub_exp_time, number_subs


def limit(band, total_exp_time, snr_target, sub_exp_time=600, binning=1, N=1, round_up=True):
    
    # Convert target SNR per binned pixel to SNR per unbinned pixel
    snr_target /= binning
    
    # Sky count rates, uses equations 8 & 9 from Dragonfly paper to convert AB mag/arcsec^2 to e/pixel/s
    rate_sky =  gamma0[band] * N * efficiency[band] * 10**(-0.4 * sky_mu[band]) # e/pixel/second
    
    # Number of sub-exposures
    number_subs = int(math.ceil(total_exp_time/sub_exp_time))
    if round_up and (total_exp_time != number_subs * sub_exp_time):
        total_exp_time = number_subs * sub_exp_time
        print('Rounding up total exposure time to next integer multiple of sub-exposure time:', total_exp_time)
        
    # Noise sources
    sky_counts = rate_sky * total_exp_time # e/pixel
    dark_counts = dark_current * total_exp_time # e/pixel
    total_read_noise = math.sqrt(number_subs) * read_noise  # e/pixel
    
    # Calculate science count rate for target signal to noise ratio
    a = total_exp_time**2
    b = -snr_target**2 * total_exp_time
    c = -snr_target**2 * (sky_counts + dark_counts + total_read_noise**2)
    rate_sci = (-b + math.sqrt(b**2 - 4*a*c))/(2*a) # e/pixel/second
    
    # Convert science count rate to surface brightness with inverse of equations 8 & 9 from Dragonfly paper
    sci_mu = -2.5 * math.log10(rate_sci / (gamma0[band] * N * efficiency[band])) # AB mag/arcsec^2
    
    return sci_mu

if __name__ == '__main__':
    snr_g, snr_r = snr(3.0, 'g', 150000), snr(30.0, 'r', 150000)
    print('g\'-band S/N =', snr_g)
    print('r\'-band S/N =', snr_r)




