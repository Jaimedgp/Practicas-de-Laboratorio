archive = ["A", "B", "C", "D", "Esp"]


def mergeTables():

	frA = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/A.csv", "r")
	frB = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/B.csv", "r")
	frC = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/C.csv", "r")
	frD = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/D.csv", "r")
	frEsp = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/Esp.csv", "r")

	fw = open("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/All.csv", "w")

	for i in range(0, 15):
		A = frA.readline()
		Esp = frEsp.readline()
		if A != '':
			if Esp != '':
				line = A.split("\n")[0]+","+frB.readline().split("\n")[0]+","+frC.readline().split("\n")[0]+","+frD.readline().split("\n")[0]+","+Esp
			else:
				line = A.split("\n")[0]+","+frB.readline().split("\n")[0]+","+frC.readline().split("\n")[0]+","+frD.readline().split("\n")[0]+", \n"
			fw.write(line)

	frA.close()
	frB.close()
	frC.close()
	frD.close()
	frEsp.close()
	fw.close()


def Calculations():
	from OpenScript import Open_file_CSV
	from SaveScript import saveLaTex

	table, index = Open_file_CSV("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/Data/All.csv")
	#saveLaTex(table, index, index.keys(), "/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/TableTex/All(1).tex")

	import matplotlib.pyplot as plt
	import numpy as np
	from scipy import stats
	from GraphPlot import nodosChebyshev

	color = ["r", "b", "g", 'c', "k"]
	slopes = [0,0,0,0,0]
	errors = [0,0,0,0,0]
	fig = plt.figure(1)

	for i in range(0,4):
		if i != 2:
			x = np.linspace(min(table[index[0]]), max(table[index[0]]), 1000)
			slopes[i], intercept, r_value, p_value, errors[i] = stats.linregress(
	                                                         table[index[3*i]], table[index[(3*i)+1]])

			plt.plot(table[index[3*i]], table[index[(3*i)+1]], color[i]+"o", label="$"+archive[i]+"$")
			plt.plot(x, slopes[i]*x+intercept, color[i])
			print archive[i]+": "+ "%g , %g" %(slopes[i], errors[i])
		else:
			x = np.linspace(min(table[index[3*i]][:-1]), max(table[index[(3*i)]][:-1]), 1000)

			slopes[i], intercept, r_value, p_value, errors[i] = stats.linregress(
	                                                         table[index[3*i]][:-1], table[index[(3*i)+1]][:-1])

			print archive[i]+": "+ "%g , %g" %(slopes[i], errors[i])
			plt.plot(x, slopes[i]*x+intercept, color[i])
			plt.plot(table[index[3*i]][:-1], table[index[(3*i)+1]][:-1], color[i]+"o", label="$"+archive[i]+"$")


	x = np.linspace(min(table[index[12]][:-4]), max(table[index[12]][:-4]), 1000)
	slopes[4], intercept, r_value, p_value, errors[4] = stats.linregress(
	                                                         table[index[12]][:-4], table[index[13]][:-4])
	print archive[4]+": "+ "%g , %g" %(slopes[4], errors[4])
	plt.plot(x, slopes[4]*x+intercept, color[4])
	plt.plot(table[index[12]][:-4], table[index[13]][:-4], color[4]+"o", label="$"+archive[4]+"$")

	slope, error = np.mean(slopes), np.std(slopes)
	slopes.append(slope)
	errors.append(error)

	plt.axhline(0, color='black').set_dashes([1,1,1,1])
	plt.xlabel("$Time_{Julian Date}$")
	plt.ylabel("$Longitude$")
	plt.legend()
	plt.show()
	fig.savefig("/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/Figure.png") # Save plot figure


	saveLaTex({"$Dots$":archive, "$Slopes/^o days^{-1}$":slopes, "$\sigma$":errors}, {0:"$Dots$", 1:"$Slopes/^o days^{-1}$", 2:"$\sigma$"}, [0,1,2], "/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/All.tex")

	synodic = [0,0,0,0,0]
	synodicer = [0,0,0,0,0]
	sideralP = [0,0,0,0,0]
	for i in range(len(synodic)):
		synodic[i] = 360/slopes[i]
		synodicer[i] = 360*errors[i]/(slopes[i]**2)

		sideralP[i] = (synodic[i] * 365.25) / (synodic[i] + 365.25)

	print np.mean(synodic[:-1]), np.std(synodic[:-1]) 
	print np.mean(synodic), np.std(synodic) 
	print np.mean(sideralP), np.std(sideralP)
	print np.mean(sideralP[:-1]), np.std(sideralP[:-1]) 


	#saveLaTex({"$Sunspot$":archive, "$Synodic/days^{-1}$":synodic, "$\sigma_{sync}/days^{-1}$":synodicer}, {0:"$Sunspot$", 1:"$Synodic/days^{-1}$", 2:"$\sigma_{sync}/days^{-1}$"}, [0,1,2], "/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/synodic.tex")
	saveLaTex({"$Sunspot$":archive, "$P_{sideral}/days$":sideralP}, {0:"$Sunspot$", 1:"$P_{sideral}/days$"}, [0,1], "/home/jaimedgp/Dropbox/VerbTeX/Astronomy/P1/sideral.tex")
Calculations()