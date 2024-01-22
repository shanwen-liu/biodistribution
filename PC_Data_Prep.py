#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Wade P Winton"
__version__ = "0.1.0"
__license__ = "MIT"

import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
import csv
import os
import time

def mono(x,a,b):
    """Returns results of an exponential in the form A*e^(b*x)"""
    return a*np.exp(b*x)
def bi(x,a,b,c,d):                          #define exponential functions for curve fitting and decay extrapolation
    """Returns results of bioexponential A*e^(b*x) + C*e^(d*x)"""
    return a*np.exp(b*x)+c*np.exp(d*x)
def tri(x,a,b,c,d,e,f):
    """Returns results of triexponential A*e^(b*x) + C*e^(d*x) + E*e^(f*x)"""
    return a*np.exp(b*x)+c*np.exp(d*x)+e*np.exp(f*x)

def fileimport(descriptor):
    """Accepts descriptor of .csv file to be imported. Returns 2d list of contents"""
    end = False
    while(not(end)):
        path = input("Please enter the file path for the .csv of the "+descriptor+": ")
        if(os.path.exists(path)):
            with open(path, 'r') as datafile:                           #Open file
                datareader = csv.reader(datafile, delimiter=',')
                data = []
                for row in datareader:
                    row = list(filter(None,row))    #Filters empty entries
                    if(row!=[]):                    #Filters empty rows
                        data.append(row)                                 #Read in values from .csv
                end = True
        else:
            print("\nPlease enter a valid file path.\n")

    return data

def fileexport(descriptor, data):
    """Accepts descriptor of .csv file and 2d list to be exported"""
    path = input("\nPlease specify a file path to save the "+ descriptor +" .csv to (or \"N/A\" if you don't want to save): ")
    if (path == "N/A"):                             #Allows user to input N/A if they don't want to save the data
        return
    else:
        with open(path,"w",newline='') as datafile: 
            datawriter = csv.writer(datafile)       #After all organs analyzed (above while loop exited), save results to new csv
            datawriter.writerows(data)

def isotopeentry():
    """Returns decay coefficient & gamma counter efficiency for user inputted isotope OR allows user to input these values"""
    end = False
    while (end == False):
        isotope = input("For the purposes of extrapolating radioactive decay, select an isotope (18F, 11C, 68Ga, or halflife in hours): ")
        if (isotope == "18F"):
            decaycoeff = np.log(2)/1.8295
            efficiency = 0.48
            end = True
        elif (isotope == "11C"):
            decaycoeff = np.log(2)/0.3398                                   #calculates lambda value in radioactive decay equation 
            efficiency = 0.48
            end = True
        elif (isotope == "68Ga"):
            decaycoeff = np.log(2)/1.12715
            efficiency = 0.48
            end = True
        elif (isotope.isnumeric()):
            decaycoeff = np.log(2)/float(isotope)
            efficiency = float(input("Please enter the efficiency for your selected isotope from the gamma counter's manual (if unknown and isotope is a common PET emitter, enter 0.48):"))
        else:
            print("\nPlease enter a valid input\n")
    return decaycoeff, efficiency

def biodistribution():
    """Function for processing gamma counter and weight data into %ID/g & (optionally) humanized %ID/organ"""
    decaycoeff, efficiency = isotopeentry()         #User inputs isotope information
    injinfo = fileimport("injection information")   #Import injection information
    maleidg = [["Timepoint", "Animal ID"]]
    femaleidg = [["Timepoint", "Animal ID"]]        #initialize 2d lists "idg" to append results to
    mlist = ["",""]                                 #initialize placeholder lists - filled with each animal's results, 
    flist = ["",""]                                     #then appended to the final "idg list"
    timeid = ""
    for a in range(1,len(injinfo)):     #a = each animal in injection information file
        mlist[1] = injinfo[a][0]        #fill in animal ID
        flist[1] = injinfo[a][0]
        if (timeid != injinfo[a][1]):   #if this animal is from a new timepoint:
            timeid = injinfo[a][1]      #update timeid to new value
            weightdata = fileimport("weight data for the "+injinfo[a][1]+" hr timepoint")   #import new weight file for this timepoint
            weightindex = 0             #reinitialize iterator for weightdata looking for each organ
            mlist[0] = timeid           #fill in timepoint into header of each animal's entry     
            flist[0] = timeid
        bq = float(injinfo[a][4])*37000 #calculates dose in Bq from injection info                         
        injtime = time.strptime(injinfo[a][5],"%m/%d/%y %H:%M") #converts str injection time into time struct (WARNING: possible mismatch between yy and YYYY. If it's a problem, check your data file and adjust!)
        countdata = fileimport("gamma counter data for animal #"+injinfo[a][0]) #imports gamma counter file for each animal
        gammaindex = 3                  #reinitialize iterator for each organ in countdata
        while(gammaindex < len(countdata)-1):
            if (str.endswith(weightdata[weightindex][0]," g") and str.endswith(weightdata[weightindex-1][0]," g")   #checks if 4 trailing g's in a row (2 tube organ)
                and str.endswith(weightdata[weightindex-2][0]," g") and str.endswith(weightdata[weightindex-3][0]," g")):

                mass1 = float(str.strip(str.removesuffix(weightdata[weightindex][0]," g")))
                mass2 = float(str.strip(str.removesuffix(weightdata[weightindex-1][0]," g")))
                mass3 = float(str.strip(str.removesuffix(weightdata[weightindex-2][0]," g")))
                mass4 = float(str.strip(str.removesuffix(weightdata[weightindex-3][0]," g")))
                mass = (mass1-mass2)+(mass3-mass4)      #calculates total organ mass

                cttime = time.strptime(countdata[gammaindex][1],"%m/%d/%Y %H:%M:%S")        #converts count time for decay correction
                timediff = (time.mktime(cttime)-time.mktime(injtime))/3600                  #finds time difference in hours for decay correction
                counts = float(countdata[gammaindex][6])+float(countdata[gammaindex+1][6])  #finds total CPM (remember, 2 tubes = 2 gamma counter entries)
                idg = mono(timediff,counts/60,decaycoeff)/efficiency/bq/mass                #finds %ID/g (converts CPM -> CPS -> Bq with detector efficiency, decay corrects, /injected dose in Bq/mass = %ID/g)

                if (injinfo[a][2]=="M"):                                #if it's a male, append to the male lists
                    if (weightdata[weightindex-4][0] in maleidg[0]):    #and if organ label already present in maleidg
                        mlist.append(idg)                               #append only calculated %ID/g value
                    else:
                        maleidg[0].append(weightdata[weightindex-4][0])     #if organ label isn't present, append organ label into row
                        mlist.append(idg)                                   #and append calculated %ID/g
                elif (injinfo[a][2]=="F"):
                    if (weightdata[weightindex-4][0] in femaleidg[0]):      #repeat with female list if the animal is a female
                        flist.append(idg)    
                    else:
                        femaleidg[0].append(weightdata[weightindex-4][0])
                        flist.append(idg)
                weightindex+=1                  #move on to next weight index
                gammaindex+=2                   #advance gamma index by 2 to account for 2 tubes
            elif (str.endswith(weightdata[weightindex][0]," g") and str.endswith(weightdata[weightindex-1][0]," g")
                  and not str.endswith(weightdata[weightindex-2][0]," g") and not str.endswith(weightdata[weightindex+1][0]," g")):
                    #if two trailing g's in a row without g's above or below (otherwise 2 tube organs mess things up)
                mass1 = float(str.strip(str.removesuffix(weightdata[weightindex][0]," g")))
                mass2 = float(str.strip(str.removesuffix(weightdata[weightindex-1][0]," g")))
                mass = mass1-mass2                                           

                cttime = time.strptime(countdata[gammaindex][1],"%m/%d/%Y %H:%M:%S")    #see description of each step above
                timediff = (time.mktime(cttime)-time.mktime(injtime))/3600
                counts = float(countdata[gammaindex][6])
                idg = mono(timediff,counts/60,decaycoeff)/efficiency/bq/mass

                if (injinfo[a][2]=="M"):
                    if (weightdata[weightindex-2][0] in maleidg[0]):
                        mlist.append(idg)    
                    else:
                        maleidg[0].append(weightdata[weightindex-2][0])
                        mlist.append(idg)
                elif (injinfo[a][2]=="F"):
                    if (weightdata[weightindex-2][0] in femaleidg[0]):
                        flist.append(idg)    
                    else:
                        femaleidg[0].append(weightdata[weightindex-2][0])
                        flist.append(idg)
                weightindex+=1
                gammaindex+=1
            else:
                weightindex+=1

        if (injinfo[a][2]=="M"):    #append this animal's data to the corresponding results list
            maleidg.append(mlist)
            mlist = [mlist[0],""]     
        elif (injinfo[a][2]=="F"):
            femaleidg.append(flist)
            flist = [flist[0],""] 
    
    maleidg = np.array(maleidg).T.tolist()      #transpose the lists for readability
    femaleidg = np.array(femaleidg).T.tolist()
    
    fileexport("Male %ID/g data",maleidg)       #export the lists to .csv files
    fileexport("Female %ID/g data",femaleidg)

    if(input("Would you like to output humanized data for dosimetry (See instructions)?  (y/n): ") == "y"):
        icrpmale = {"BRAI":1450,"EYES":15,"HEAR":840,"LUNG":1200,"LIVE":1800,
                    "PANC":140,"SPLE":150,"ADRE":14,"KIDN":310,"ADIP":18200,
                    "STOM":150,"COST":250,"SINT":650,"COSI":350,"CEAC":150,
                    "COCE":150,"LINT":150,"COLI":75,"TEST":35,"MUSC":29000,
                    "BONE":10450,"BLOO":5600,"BW":73000}
        icrpfemale = {"BRAI":1300,"EYES":15,"HEAR":840,"LUNG":950,"LIVE":1400,
                        "PANC":120,"SPLE":130,"ADRE":13,"KIDN":275,"ADIP":22500,    #dictionary of ICRP organ weights
                        "STOM":140,"COST":230,"SINT":600,"COSI":280,"CEAC":145,
                        "COCE":160,"LINT":145,"COLI":80,"OVAR":11,"UTER":80,
                        "MUSC":17500,"BONE":7760,"BLOO":4100,"BW":60000}
        for a in range(1,len(maleidg[0])):                              #for each animal                        
            for b in range(2, len(maleidg)):                            #for each organ (this being [b][a] was less than intuitive to me, hence transposing the lists)
                maleidg[b][a] = mono(float(maleidg[0][a]), float(maleidg[b][a])*
                                    float(injinfo[injinfo.index(next(animal for animal in injinfo if animal[0] == maleidg[1][a]))][3])/
                                    icrpmale["BW"]*icrpmale[maleidg[b][0]],np.negative(decaycoeff)) 
                #lots of lists, but basically, takes %ID/g, multiplies by animal bodyweight, divides by ICRP bodyweight, and multiplies by ICRP organ weight.
                #does all of this in a call to mono() to reapply radioactive decay
        for a in range(1,len(femaleidg[0])):    
            for b in range(2, len(femaleidg)):                      #same calculation for female list
                femaleidg[b][a] = mono(float(femaleidg[0][a]), float(femaleidg[b][a])*
                                    float(injinfo[injinfo.index(next(animal for animal in injinfo if animal[0] == femaleidg[1][a]))][3])/
                                    icrpfemale["BW"]*icrpfemale[femaleidg[b][0]],np.negative(decaycoeff))
        maleidg.pop(1)
        femaleidg.pop(1)    #removes Animal ID
        mhumanized = []     #initializes lists for averaged, humanized %ID/organ data
        fhumanized = []

        for a in range(0,len(maleidg)):
            mhumanized.append([maleidg[a][0]])      #initializes header/organ information
            maleidg[a].append("")                   #appends an empty row so next loop can cycle through every entry
        for a in range(0,len(femaleidg)):
            fhumanized.append([femaleidg[a][0]])
            femaleidg[a].append("")

        #next section averages humanized %ID/organ results of animals in each timepoint
        firstindex = 1              #index of first animal in a timepoint to be averaged
        lastindex = 0               #index of last animal
        timeid = maleidg[0][1]      #initializes value of first timepoint for comparison
        for a in range(1,len(maleidg[0])):      #looks at each timepoint label in the header
            if (maleidg[0][a] != timeid):       #if this timepoint isn't the same as the previous one (i.e. you just passed the last animal at that timepoint)
                mhumanized[0].append(timeid)    #append the timepoint to the header of the results tab
                lastindex = a                   #update the index of the last animal to be averaged
                for b in range(1,len(maleidg)): #for each organ
                    mhumanized[b].append(np.average(maleidg[b][firstindex:lastindex]))  #append the average organ values for this range of animals at this timepoint
                timeid = maleidg[0][a]          #update timeid
                firstindex = a                  #update first index of the first animal of the next timepoint
        firstindex = 1
        lastindex = 0
        timeid = femaleidg[0][1]
        for a in range(1,len(femaleidg[0])):    #repeat the averaging process for the female list
            if (femaleidg[0][a] != timeid):
                fhumanized[0].append(timeid)
                lastindex = a
                for b in range(1,len(femaleidg)):
                    fhumanized[b].append(np.average(femaleidg[b][firstindex:lastindex]))
                timeid = femaleidg[0][a]
                firstindex = a
        fileexport("male humanized %ID/organ data",mhumanized)      #export the humanized %ID/organ data to .csv
        fileexport("female humanized %ID/organ data",fhumanized)
        
def dosimetry():
    """Inputs the humanized %ID/organ data, calculates integral 0 -> infinity based on either fitted curve or trapezoidal method"""
    decaycoeff, filler = isotopeentry()
    remainder = 1/decaycoeff
    data = fileimport("Humanized %ID/organ data")
    bladderflag = True
    counter =0
    
    results = []                                                        #Initialize 2d list "results"
    results.append(["Organ","A","b","C","d","E","f","Integral"])        #Header with column labels
    for i in range(1,len(data)):
        results.append([data[i][0],"","","","","","",""])               #Organ labels + blank spaces for coefficients and integral
    organid = 1
    while (bladderflag):
        for i in range(0,len(data)):
            for j in range(1,len(data[i])):
                data[i][j]=float(data[i][j])                #Typecast all values str --> float

        if (len(data[0])-1<4):
            print("\nBecause your data contains fewer than 4 timepoints, you can only fit 1 exponential term.\n")
            maxterms = 1
        elif (len(data[0])-1<6):
            print("\nBecause your data contains fewer than 6 timepoints, you can only fit up to 2 exponential terms.\n")
            maxterms = 2
        else:
            maxterms = 3            #Inform user of limited number of fittable coefficients (inbuilt limitation of "optimize")
    
        fit_t=np.linspace(0,(data[0][-1]),100)      #Create linearly spaced range of time from t=0 to final timepoint for graphical display of fitted curve
        pretime = np.linspace(0,data[0][1],100)     #For Trap method: create linearly spaced ranges of time from t=0 to first timepoint            
        posttime = np.linspace(data[0][-1],2*data[0][-1],100)       #and last timepoint to 2*last timepoint for graphical representation of extrapolated integrals

        end1 = False
        while (not(end1)):      #cycles through each organ allowing user to fit/integrate each organ's time activity curve (TAC)
            fit_data = []
            mode = input("\nSelect a mode for integrating the "+data[organid][0]+" data (\"trap\" or a number of exponential terms up to "+str(maxterms)+"): ")
            if (mode == "1"):
                fit, filler = optimize.curve_fit(mono,data[0][1:],data[organid][1:])    #Fits monoexponential curve to time and organ data
                plt.plot(data[0][1:],data[organid][1:],"o")                             #scatter plot of time vs organ data
                for i in range (0,len(fit_t)):
                    fit_data.append(mono(fit_t[i],fit[0],fit[1]))                       #fills fit_data with corresponding values predicted by fitted curve
                plt.plot(fit_t,fit_data)
                plt.title(data[organid][0])
                plt.xlabel("Time (h)")
                plt.ylabel("("+data[organid][0]+" MBq/Injected MBq)")
                plt.show()

                if(input("\nWould you like to reanalyze this organ in a different mode? (y/n): ")=="n"):    #if user is satisfied (i.e. curve converges),
                    integral = integrate.quad(mono,0,np.inf,(fit[0],fit[1]),1)                              #computes definite integral of TAC from 0 --> inf
                    results[organid]=[results[organid][0],fit[0],fit[1],"","","","",integral[0]]            #stores coefficients and integral to "results"
                    remainder-=integral[0]
                    organid+=1                                                                              #increments organid so the next organ will be analyzed

            elif (mode == "2" and int(mode) <= maxterms):
                fit, filler = optimize.curve_fit(bi,data[0][1:],data[organid][1:])                #repeats process above for bi and tri exponentials as desired
                plt.plot(data[0][1:],data[organid][1:],"o")
                for i in range (0,len(fit_t)):
                    fit_data.append(bi(fit_t[i],fit[0],fit[1],fit[2],fit[3]))
                plt.plot(fit_t,fit_data)
                plt.title(data[organid][0])
                plt.xlabel("Time (h)")
                plt.ylabel("("+data[organid][0]+" MBq/Injected MBq)")
                plt.show()

                if(input("\nWould you like to reanalyze this organ in a different mode? (y/n): ")=="n"):
                    integral = integrate.quad(bi,0,np.inf,(fit[0],fit[1],fit[2],fit[3]),1)
                    if (len(integral)>3):                                                       #checks if "integral" contains divergence warning message (bi and tri
                        x_int = optimize.root_scalar(bi,(fit[0],fit[1],fit[2],fit[3]),bracket=[data[0][-1],data[0][-1]*2])             #occassionally diverge to -inf)
                        integral = integrate.quad(bi,0,x_int.root,(fit[0],fit[1],fit[2],fit[3]))            #if 0->inf diverges, uses root_scalar to find x-int,
                    results[organid]=[results[organid][0],fit[0],fit[1],fit[2],fit[3],"","",integral[0]]    #and recalculates integral from 0->x-int
                    remainder-=integral[0]
                    organid+=1
            elif (mode == "3" and int(mode) <= maxterms):
                fit, filler = optimize.curve_fit(tri,data[0][1:],data[organid][1:])
                plt.plot(data[0][1:],data[organid][1:],"o")
                for i in range (0,len(fit_t)):
                    fit_data.append(tri(fit_t[i],fit[0],fit[1],fit[2],fit[3],fit[4],fit[5]))
                plt.plot(fit_t,fit_data)
                plt.title(data[organid][0])
                plt.xlabel("Time (h)")
                plt.ylabel("("+data[organid][0]+" MBq/Injected MBq)")
                plt.show()

                if(input("\nWould you like to reanalyze this organ in a different mode? (y/n): ")=="n"):
                    integral = integrate.quad(tri,0,np.inf,(fit[0],fit[1],fit[2],fit[3],fit[4],fit[5]),1)
                    if (len(integral)>3):
                        x_int = optimize.root_scalar(bi,(fit[0],fit[1],fit[2],fit[3]),bracket=[data[0][-1],data[0][-1]*2])
                        integral = integrate.quad(tri,0,x_int.root,(fit[0],fit[1],fit[2],fit[3],fit[4],fit[5]))
                    results[organid]=[results[organid][0],fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],integral[0]]
                    remainder-=integral[0]
                    organid+=1
            elif (mode == "trap"):      #trap mode = trapezoidal method, used for data where no good TAC can be calculated (e.g. poor fits, divergent integrals, etc.)
                predata = []
                postdata = []                   
                for i in reversed(range(0,len(pretime))):               #Use decay function to back-calculate initial activity at time = 0, and all values in between 
                    predata.append(mono(pretime[i],data[organid][1],decaycoeff))    #This produces a smooth exponential decay curve between time = 0 and first timepoint
                for i in range(0,len(posttime)):                                                    #repeat process of generating smooth exponential decay curve between 
                    postdata.append(mono(data[0][-1]-posttime[i],data[organid][-1],decaycoeff))    #final timepoint and (final timepoint) * 2
            
                integral = integrate.trapezoid(data[organid][1:],data[0][1:])                 #Compute interpolative integral using trapezoidal method
                preint = integrate.quad(mono,0,data[0][1],(mono(0,data[organid][1],decaycoeff),decaycoeff))     #Compute integral of extrapolated data between 0 and first timepoint
                postint = integrate.quad(mono,data[0][-1],np.inf,(data[organid][-1],np.negative(decaycoeff)))    #and extrapolated data between final timepoint and inf
                integral += preint[0]+postint[0]                        #sum interpolatative integral with extrapolated integrals.
    
                plt.plot(data[0][1:],data[organid][1:])
                plt.plot(pretime,predata)
                plt.plot(posttime,postdata)
                plt.fill_between(data[0][1:],data[organid][1:])                #plot the graphical approximations of interpolated and extrapolated integrals
                plt.fill_between(pretime,predata)
                plt.fill_between(posttime,postdata)
                plt.title(data[organid][0])
                plt.xlabel("Time (h)")
                plt.ylabel("("+data[organid][0]+" MBq/Injected MBq)")
                plt.show()

                if(input("\nWould you like to reanalyze this organ in a different mode? (y/n): ")=="n"):    #if user is satisfied
                    results[organid][-1]=integral               #save integral in last slot in "results" row corresponding to the organ
                    remainder-=integral
                    organid+=1
            else:
                print("\nPlease input a valid mode\n")

            if(organid > len(data)-1):          #After a given organ is analyzed and user is satisfied (the only way organid is incremented), check if final organ
                end1 = True
        if(counter == 0 and input("Do you have bladder imaging data to fit? (y/n): ")=="y"):
            bladderflag = True
            results.append(["Bladder","","","","","","",""])
            bladderdata = fileimport("bladder %ID/organ data")
            data[0]=bladderdata[0]
            data.append(bladderdata[1])
            counter+=1
        else:
            bladderflag = False
    results.append(["Remainder","","","","","","",remainder])
    fileexport("coefficient and integral data",results)

def main():
    end = False
    while(not(end)):
        mode = input("Select a mode (\"bio\" for biodistribution, \"dose\" for dosimetry preprocessing, or \"end\"): ")
        if(mode == "bio"):
            biodistribution()
        elif(mode == "dose"):
            dosimetry()
        elif(mode == "end"):
            end = not(end)
        else:
            print("\nPlease enter a valid mode.\n")


if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
