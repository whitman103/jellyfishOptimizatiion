import matplotlib.pyplot as plt

inFile=open("results.txt",'r')

trueParams=inFile.readline()
trueData=trueParams.split(" ")
print(trueData)

paramHist=[]
for line in inFile:
	data=line.split(" ")
	data.pop()
	paramHist.append(float(data[0]))

plt.hist(paramHist,bins=20)
plt.show()