import sys
sys.path.insert(0, '/home/jaimedgp/J_Graphics_P/') # some other scripts from other project are going to be used

folderPath = '/home/jaimedgp/Dropbox/VerbTeX/Astronomy/Pleiades/' #folder path 

# absolute visuak magnitudes given by the teacher
MVgiven = [-5.8, -4.1, -1.1, -0.7, 2.0, 2.6, 3.4, 4.4, 5.1, 5.9,
                                            7.3, 9.0, 11.8, 16.0]
# color index (B-V) given by the teacher
BVgiven = [-0.35, -0.31, -0.16, 0.00, 0.13, 0.27, 0.42, 0.58,
                               0.70, 0.89, 1.18, 1.45, 1.63, 1.80]

def ObtainData():
	
	# Use another script to open and obtain the data from a TXT
	# OpenScript in: https://github.com/Jaimedgp/J_Graphics_P/OpenScript.py
	from OpenScript import Open_file_TXT

	# Stock the data in 2 dictionaris with titles and values
	# table = {title:values} ; index = {0:title}
	table, index = Open_file_TXT(folderPath+'Data/OBSERVATIONS.TXT')

	OurData = {'Object-ID':table['Object-ID'],'$m_B$':table['B'],
											  '$m_V$':table['V']} # in this experiment only need B, V and the star ID data
	Ourindex = {0:'Object-ID', 1:'$m_B$', 2:'$m_V$'}

	for i in range(len(OurData['$m_B$'])):
		 OurData['$m_B$'][i] = float(OurData['$m_B$'][i])
	for i in range(len(OurData['$m_V$'])):
		 OurData['$m_V$'][i] = float(OurData['$m_V$'][i])

	# Use another script to save athe data in LaTex table format
	# SaveScript in: https://github.com/Jaimedgp/J_Graphics_P/SaveScript.py
	from SaveScript import saveLaTex

	saveLaTex(OurData, Ourindex, [0,1,2], folderPath+'TexTables/data.tex') #save the data in Latex table format

	return OurData, Ourindex

def Plot(data, index, axes, name, types):
	import matplotlib.pyplot as plt

	fig = plt.figure(1)
	
	# diferents axes scale depends on if are 
	# the theorical or the experimental vales
	if types == "A":
		plt.ylim(17, 0)
		marker="b*"
	elif types == "B":
		plt.ylim(17, -8)
		marker="ro"

	plt.plot(data[index[axes[0]]], data[index[axes[1]]], marker)
	plt.xlabel(index[axes[0]], fontsize=20)
	plt.ylabel(index[axes[1]], fontsize=20)
	plt.xlim(-0.4, 1.8)
	plt.show()

	fig.set_size_inches(9,6)
	fig.savefig(folderPath+"Figures/"+name+".png")

def Plot2(data, index, axes, name):
	import matplotlib.pyplot as plt
	from scipy import stats
	import numpy as np

	##############################################################
	## THEORICAL ABSOLUTE VISUAL MAGNITUDE M_V VS THEORICAL B-V ##
	##############################################################

	MVth = data[index[axes[1]]]
	BVth = data[index[axes[0]]]
	BV = []
	MV = []

	# Take only the center values
	i = 0
	while i < len(BVth):
		if BVth[i] < 1.6 and BVth[i] > 0.2:
			BV.append(BVth[i])
			MV.append(MVth[i])
		i = i+1

	# Calculate the linear regression
	slope, intercept, r_value, p_value, error = stats.linregress(BV, MV)

	BVmean = np.mean(BV)
	sx2 = ((BV-BVmean)**2).sum()
	std_intercept = error * np.sqrt(1./len(BV) + BVmean*BVmean/sx2)

	slopes = {"M":[slope, error]}
	intercepts = {"M":[intercept, std_intercept]}

	# Plot the results
	fig, ax1 = plt.subplots()
	ax1.plot(BV, MV, 'ro')

	#Represent the regression
	x = np.linspace(min(BV), max(BV), 1000)
	plt.plot(x, slope*x+intercept, 'r')
	
	ax1.set_xlabel(index[axes[0]], fontsize=20)
	ax1.set_ylim(17, 0)
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_ylabel(index[axes[1]], color='r', fontsize=20)
	ax1.set_xlim(0.2, 1.6)
	ax1.tick_params('y', colors='r')

	##################################################################
	# EXPERIMENTAL APPARENT VISUAL MAGNITUDE m_V VS EXPERIMENTAL B-V #
	##################################################################

	mVexp = data[index[axes[3]]]
	BVexp = data[index[axes[2]]]
	BVs = []
	mV = []

	# Take only the center values
	i = 0
	while i < len(BVexp):
		if BVexp[i] < 1.6 and BVexp[i] > 0.2:
			if mVexp[i] > 8:
				BVs.append(BVexp[i])
				mV.append(mVexp[i])
		i = i+1

	# Calculate the linear regression
	slope, intercept, r_value, p_value, error = stats.linregress(BVs, mV)

	
	BVmean = np.mean(BVs)
	sx2 = ((BVs-BVmean)**2).sum()
	std_intercept = error * np.sqrt(1./len(BVs) + BVmean*BVmean/sx2)

	slopes["m"]= [slope, error]
	intercepts["m"]= [intercept, std_intercept]

	ax2 = ax1.twinx()
	ax2.plot(BVs, mV, 'b*')

	#Represent the regression
	x = np.linspace(min(BVs), max(BVs), 1000)
	plt.plot(x, slope*x+intercept, 'b')

	ax2.set_ylabel(index[axes[3]], color='b', fontsize=20)
	ax2.set_ylim(25, 0)
	ax2.tick_params('y', colors='b')

	fig.tight_layout()
	plt.show()

	fig.set_size_inches(9,6)
	fig.savefig(folderPath+"Figures/"+name+".png") #save graph

	return slopes, intercepts

Data, index = ObtainData()

BlessV = []

for i in range(len(Data['$m_B$'])):
	# the color index is obtained by rest B-V
	BlessV.append(Data['$m_B$'][i]-Data['$m_V$'][i])

# add the experimental color index to data
index[3] = '$m_B-m_V$'
Data['$m_B-m_V$'] = BlessV

# add theorical color index for a given group of stars
index[4] = '$M_V$'
Data['$M_V$'] = MVgiven
index[5] = '$[m(B)-m(V)]_{teo}$'
Data['$[m(B)-m(V)]_{teo}$'] = BVgiven

# represent apparent visual magnitude and experimental B-V
#Plot(Data, index, [3, 2], 'expAparent', "A")

# represent absolute visual magnitude and theorical B-V
#Plot(Data, index, [5, 4], 'theoabsolute', "B")

#calculate and plot both apparent and absolute magnitude
slopes, intercepts = Plot2(Data, index, [5, 4, 3, 2], 'both')

print ("Slopes: \n \t m = %g +- %g" %(slopes["m"][0], slopes["m"][1]) 
	               + "\n \t M = %g +- %g \n" %(slopes["M"][0], slopes["M"][1]))
print ("Intercept: \n \t m = %g +- %g" %(intercepts["m"][0], intercepts["m"][1]) 
	               + "\n \t M = %g +- %g \n" %(intercepts["M"][0], intercepts["M"][1]))

########################################
#   CALCULATE THE DISTANCE MODULUS    ##
########################################

# Use another script to calculate the error from the expression
# WidgetsScript in: https://github.com/Jaimedgp/J_Graphics_P/WidgetsScript.py
from WidgetsScript import ErrorsCalculator

#distance modulus
mM = intercepts["m"][0] - intercepts["M"][0]

STD = ErrorsCalculator("m M",[intercepts["m"][0], intercepts["M"][0]],
                           [intercepts["m"][1], intercepts["M"][1]], "m - M")

print "M-m = %g +- %g \n" %(mM, STD)

########################################
#     CALCUALTE THE DISTANCE IN PC    ##
########################################

D = 10**((mM/5)+1)
STD_D = ErrorsCalculator("mM",[mM], [STD], "10**((mM/5)+1)")

print "D = %g +- %g \n" %(D, STD_D)
