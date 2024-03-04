import matplotlib.pyplot as plt
import numpy as np

a = [20, 40, 60, 80, 100]
c = [12.609, 11.55, 11.187, 10.945, 10.833]

plt.plot(a,c)
plt.xlabel('#cores')
plt.ylabel('Runtime')
plt.title('Runtime vs number of cores')
plt.show()
