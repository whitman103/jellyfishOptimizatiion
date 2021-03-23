import matplotlib.pyplot as plt

inFile=open("results.txt",'r')

trueParams=inFile.readline()
trueData=trueParams.split(" ")
trueData.pop()

paramHist=[[] for i in range(len(trueData))]
for line in inFile:
	data=line.split(" ")
	data.pop()
	for index,element in enumerate(data):
		paramHist[index].append(float(element))

fitnessValues=paramHist.pop(0)
for index,element in enumerate(paramHist):
	plt.scatter(element,[index for i in range(len(element))])
plt.show()