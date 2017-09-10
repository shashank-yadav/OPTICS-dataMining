import matplotlib.pyplot as plt
import numpy as np

with open('result.txt') as f:
    array = [float(line) for line in f]

x = np.arange(len(array))
plt.bar(x , array,  align='center', alpha=0.5, color='purple' )

plt.plot( array )
plt.title("Reachability Plot")
plt.xlabel("Data Point")
plt.ylabel("Reachability Distance")
plt.savefig('plot.png')

# plt.show()
