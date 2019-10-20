# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 20:21:05 2019

@author: davidd
"""
###########################@author: davidd###############################
'''Binomial pricing Model '''
'''opt: 1=call, -1=put
   exer= option type (1=European, 2=American)
   div=dividend
   S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   sigma=volatility
   T=maturity
   nstep=number of steps'''

def BinOptCRR_Euro_Ameri(opt, exer, S, K, r, div, T, sigma, nstep):
    import numpy as np
    if S > 0.0 and K > 0.0 and T > 0.0 and sigma > 0.0:
        delta = T / nstep
        act = np.exp(r * delta)
        actd = np.exp((r - div) * delta)
        u = np.exp(sigma * np.sqrt(delta))
        d = 1 / u
        p = (actd - d) / (u - d)
        pstar = 1 - p
        
        
        #bb=np.empty([3,2], dtype=int)
        a=np.empty([nstep+1], dtype=float)
        
        i=0
        while i<a.size:
            a[i]=np.maximum(opt*(S*(u**i) * (d**(nstep-i)) -K), 0)
            #print(a[i], end=' ')
            i+=1
        
        '''a=np.arange(6)
        for x in np.nditer(a):
            print(x, end='')           
            
        i=0
        while i<a.size:
            print(a[i])
            i+=1'''
        
        #j=nstep-1
        for j in range(nstep,-1,-1):
            i=0
            while i<j:
                a[i]= ((p*a[i+1]) + (pstar * a[i]))/act
                #print("i:j", i,":",j, a[i])
                if exer == 2:
                    a[i] = np.maximum(a[i], opt * (S * (u ** i) * (d ** (j - i)) - K))
                i+=1
                
        if opt ==1:
            tyo="Call"
        else:
            tyo="Put"
            
        if exer ==1:
            echop="European"
        else:
            echop="American"

        print("My " + echop + " " + tyo, " option value is : " + "%.6f"%(a[0]))
        

#BinOptCRR_Euro_Ameri(opt, exer, S, K, r, div, T, sigma, nstep):            
#BinOptCRR_Euro_Ameri(1,1,50,50,0.05,0,0.5,0.3,100)=4.806953
BinOptCRR_Euro_Ameri(1,2,0.79,0.7950,-0.04,0,0.75,0.04,3)    
###########################@author: davidd###############################
  
'''Black and Scholes Option Valuation '''
'''opt: 1=call, -1=put
   S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   sigma=volatility
   T=maturity'''

def BSOpValuation(opt, S, K, r, q, T, sigma):
    import numpy as np
    from scipy import stats
    
    act = np.exp(-r * T)
    actd = np.exp(-q * T)
    S=float(S)
    
    if S > 0 and K > 0 and T > 0 and sigma > 0:
        Done = (np.log(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
        Dtwo = (np.log(S / K) + (r - q - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
        
        Ndone = stats.norm.cdf(opt * Done, 0.0, 1.0)
        Ndtwo = stats.norm.cdf(opt * Dtwo, 0.0, 1.0)
        
        return opt * (S * actd * Ndone - K * act * Ndtwo)

    elif S > 0 and K > 0 and sigma > 0 and T == 0:
        return np.maximum(0, opt * (S - K))
    
    else:
        print("Warning : Share price(S) and Exercice price(X) and time to maturity(tyr) and volatility(sigma) must be positive")
    
'''value = (S0 * stats.norm.cdf(d1, 0.0, 1.0)
- K * exp(-r * T) * stats.norm.cdf(d2, 0.0, 1.0))'''

'''BSOpValuation(1;100;95;8%;3%;0.5;20%)'''

#BSOpValuation(1,50,50,0.05,0,0.5,30/100)=4.817438
###########################@author: davidd###############################
'''S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   sigma=volatility
   T=maturity'''

'''NDONE Function: returns the Black - Scholes d1 value '''
def BSDOne(S, K, r, q, T, sigma):
    import numpy as np
    return (np.log(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))


'''NDTWO Function : return the Black - Scholes d2 value'''
def BSDTwo(S, K, r, q, T, sigma):
    import numpy as np
    return (np.log(S / K) + (r - q - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

'''NDASHONE Function '''
def BSNdashOne(S, K, r, q, T, sigma):
    import numpy as np
    return (np.exp(-(BSDOne(S, K, r, q, T, sigma)) ** 2 / 2) / np.sqrt(2 * np.pi))
###########################@author: davidd###############################

'''Black and Scholes Greeks Valuation '''
''' greek: delta(1), gamma(2), rho(3), theta(4) or vega(5)
opt=1 for call, opt=-1 for put; q=dividend yield; r= interest rate; T=maturity
S=spot or current price; K=Strike, sigma=volatility '''

def BSValGreeks(greek, opt, S, K, r, q, T, sigma):
    from scipy import stats
    import numpy as np
    act = np.exp(-q * T)
    
    c = BSOpValuation(opt, S, K, r, q, T, sigma)
    c1 = stats.norm.cdf(opt * BSDOne(S, K, r, q, T, sigma))
    c1d = BSNdashOne(S, K, r, q, T, sigma)
    c2 = stats.norm.cdf(opt * BSDTwo(S, K, r, q, T, sigma))
    d = opt * act * c1
    g = c1d * act / (S * sigma * np.sqrt(T))
    
    if greek == 1:
        return d
    elif greek == 2:
        return g
    elif greek == 3: 
        return opt * K * T * np.exp(-r * T) * c2
    elif greek == 4:
        return r * c - (r - q) * S * d - 0.5 * (sigma * S) ** 2 * g
    elif greek == 5:
        return S * np.sqrt(T) * c1d * act

#BSValGreeks(1,1,49,50,0.05,0,0.3846,20/100)=0.5216016
        
###########################@author: davidd###############################
'''Monte Carlo Quasi Random valuation option '''
'''S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   sigma=volatility
   T=maturity
   nsim: number of simulations'''

def MCValQuasiRandom(S, K, r, q, sigma, T, nsim):
    import numpy as np
    from scipy import stats
    from scipy.stats import norm
    import random as rdm
    
    drift = (r - q - 0.5 * sigma ** 2) * T
    volatility = sigma * np.sqrt(T)
    summ = 0
    
    i=0
    for i in range(nsim+1):
        a = norm.ppf(rdm.random()) #inverse of normal distribution of numbers between 0 and 1 excluded
        b = S * np.exp(drift + volatility * a)
        summ = summ + np.maximum((b - K), 0)
    
    return np.exp(-r * T) * summ / nsim

MCValQuasiRandom(50,50,0.05,0,0.3,0.5,1000)
#MCOptionValueQuasiRandom(50,50,0.05,0,0.3,0.5,1000)=4 et 5
#MCOptionValueQuasiRandom(1;100;0.08;0.03;20%;0.5;95;1000)=9 or 10
###########################@author: davidd###############################

'''Monte Carlo Greeks Valuation '''
'''opt: 1=call, 2=put
   greek: 1=delta, 3=rho, 5=vega
   price = price of the option whose we are looking for volatility  
   S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   sigma=volatility
   T=maturity
   nsim: number of simulations'''
#code of sign using in Monte Carlo greek valuation
def Sign(x):
    if x>0:
        return 1
    elif x<0:
        return -1
    else:
        return 0    
#true code        
def MCValGreek(greek, opt, S, K, r, q, sigma, T, nsim):
    import numpy as np
    from scipy import stats
    from scipy.stats import norm
    import random as rdm

    drift = (r - q - 0.5 * sigma ** 2) * T
    volatility = sigma * np.sqrt(T)
    summ = 0
    ggreek=0
    
    i=0
    for i in range(nsim+1):
        a = norm.ppf(rdm.random()) #inverse of normal distribution of numbers between 0 and 1 excluded
        b = S * np.exp(drift + volatility * a)
        summ = summ + np.maximum((b - K), 0)
        
        if (greek == 1 and Sign(opt * (b - K)) == 1):
            ggreek = ggreek + b
            #print(Sign(opt * (b - K))
        elif (greek == 3 and Sign(opt * (b - K)) >= 0):
            ggreek = ggreek + 1
        elif (greek == 5 and Sign(opt * (b - K)) >= 0):
            ggreek = ggreek + b * (np.log(b / S) - drift)

    if greek == 1:
        return np.exp(-r * T) * (ggreek / S) / nsim
    elif greek == 3: 
        return np.exp(-r * T) * K * T * ggreek / nsim
    elif greek == 5:
        return np.exp(-r * T) * (ggreek / sigma) / nsim
    
#MCValQuasiRandom(1,100,95,0.08,0.03,0.2,0.5,10000)
#MCValGreek(greek, opt, S, K, r, q, sigma, T, nsim):
#QMCOptionGreek135(1;1;49;0.05;0;20%;0.3846;50;50000)

#MCValGreek(1,1,49,50,0.05,0.0,0.2,0.3846,10000)=0.51 and 0.53
###########################@author: davidd###############################

'''Implied Volatility'''
'''opt: 1=call, 2=put
   price = price of the option whose we are looking for volatility  
   S=spot price, current price
   K=strike
   r=interest rate
   q=dividend
   T=maturity'''

def ImpliedVol(opt, price, S, K, r, q, T):
    sigma = 0.01
    for i in range(300):
        priceEsti = BSOpValuation(opt, S, K, r, q, T, sigma)
        if priceEsti >= price:
            print("Implied Volatility is : " + "%4.2f"%(sigma*100), "%")
            break
        sigma = sigma + 0.01

    if i > 300:
        print("Volatility over 300%")

price=9.77936 #I suppose option price equals to 9.77936
#ImpliedVol(1,price,100,95,0.08,0.03,0.5)=21%
###########################@author: davidd###############################

''' Bond Valuation '''
'''T = maturity ie T=10 means 10 years 
   coupon = coupon
   faceValue= face value'''

def BondValue(T, coupon, faceValue):
    import numpy as np
    intervalle = 1
    total = 0
    
    #aa=np.empty([nstep+1, 1], dtype=float)
    #T=10
    matzeroprice=np.empty([T], dtype=float)
    matcoupon=np.empty([T], dtype=float)
    matpvbond=np.empty([T], dtype=float)
    mattime=np.empty([T], dtype=float)

    maturity=np.array([1,2,3,4,5,6,7,8,9,10], dtype=int)
    #Zero_price=np.exp(-rt)
    zero_price=np.array([0.941,0.885,0.83,0.777,0.726,0.677,0.63,0.585,0.542,0.502], dtype=float)
    
    zero_yield=np.divide(-np.log(zero_price),maturity)
    
    for i in range(1,T+1,1):
        if i==1:
            mattime[i-1]=i
        if i>1:
            mattime[i-1]=mattime[i-2] + intervalle
            
        matzeroprice[i-1]=faceValue * np.exp(-mattime[i-1]*zero_yield[i-1])
        
        if i < T:
            matcoupon[i-1]=coupon
        else:
            matcoupon[i-1]=1+ coupon
        
        matpvbond[i-1] = matzeroprice[i-1] * matcoupon[i-1]
        
        total=total+matpvbond[i-1]
        
        
    print("Bond value is : ", total)

#BondValue(10,0.05,1)=0.85675

###########################@author: davidd###############################
''' VaR Valuation '''
#dist (1=Loss distribution, 2=Profit and Loss distribution)
'''obs=number of observation
   alpha=confidence level'''

def ValAtRisk(obs, alpha):
    import numpy as np
    import random as rdm
    
    dist=int(input("Please enter a type of distribution, 1=Loss, 2=PnL : "))
    i=0
    obser=np.empty([obs], dtype=float)
    
    for i in range(obs):
        obser[i] = rdm.uniform(-3000000,1000000) #random numbers between -300000 and 100000 excluded
        #print(obser[i])
    
    if dist == 2:
        alpha = 1 - alpha
        
    #print(obser.size)

    varr = np.percentile(obser, alpha)
    
    print("VaR is : " + "%6.6f"%(varr))
    

ValAtRisk(1000,0.99)    
###########################@author: davidd###############################
'''CVaR valuation '''
'''obs=number of observation
   alpha=confidence level'''

def CVaRisk(obs, alpha):
    import numpy as np
    import random as rdm
    
    dist=int(input("Please enter a type of distribution, 1=Loss, 2=PnL : "))
    i=0
    obser=np.empty([obs], dtype=float)
    
    for i in range(obs):
        obser[i] = rdm.uniform(-3000000,1000000) #random numbers between -300000 and 100000 excluded
        #print(obser[i])
    
    if dist == 2:
        alpha = 1 - alpha
        
    #print(alpha)
    
    varr = np.percentile(obser, alpha)
        
    count = 0
    average = 0
    
    for number in obser:
        if number < varr:
            count = count + 1
            average = average + number
    
    cvarr = average / count

    print("CVaR is : " + "%6.6f"%(cvarr))

CVaRisk(1000, 0.99)    
###########################@author: davidd###############################
'''Sources
John Hull 10th Edition Options, Futures, and Other Derivatives 2017 (french version)
Mary Jackson and Mike Staunton : Advanced Modelling in Finance using Excel-VBA -2001
Patrice Poncet and Roland Portait : Finance de Marché 4th Edition '''
