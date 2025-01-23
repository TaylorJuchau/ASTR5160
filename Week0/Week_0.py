#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
import matplotlib.pyplot as plt
data = np.arange(1,11) #create array of numbers from 1-10
#print(data) #verify list is what it should be

#create txt file called example_text, which will be rewritten every time this is run
with open("example_text.txt", "w") as file:
    for i in range(0,len(data)):
        file.write(f"{data[i]}\t {data[i]} \n")
#read file contents to check
with open("example_text.txt", 'r') as file:
        txt_file = file.read()
print(txt_file)


# In[38]:


data = np.loadtxt('example_text.txt')
plt.plot(data[:,0],data[:,1], label = 'boring data as a line')
plt.scatter(data[:,0],data[:,1], marker = 'x', color = 'y', label = 'boring data as yellow crosses')
plt.xlabel('first column')
plt.ylabel('second column')
plt.title('plotting this rather boring dataset')
plt.legend()


# In[39]:


def plot_the_boring_data():
    data = np.loadtxt('example_text.txt')
    plt.plot(data[:,0],data[:,1], label = 'boring data as a line')
    plt.scatter(data[:,0],data[:,1], marker = 'x', color = 'y', label = 'boring data as yellow crosses')
    plt.xlabel('first column')
    plt.ylabel('second column')
    plt.title('plotting this rather boring dataset')
    plt.legend()
    return None
if __name__ == "__main__":
    plot_the_boring_data()
    


# In[ ]:





# In[ ]:




