import numpy as np
from scipy import optimize
from scipy.linalg import sqrtm

p = 10 #Number of nodes in the graph

n = 100 #Number of datapoints

mu = 1

lamb = 0.1

np.random.seed(2018)

theta = np.random.rand(p,p)

X = np.random.multivariate_normal([0]*p, np.identity(p), (n))

X = np.where(X >= 0,1,-1)

def cond(X,theta,r):
    """ Return the Conditional Probability of x_r given x/r for one sample of data. 
        Arguments: 
        X - One sample of data, 
        theta - The Corresponding theta Vector from the full theta matrix
        r = The node whose neighbourhood is to be estimated 
    """
    x_r = X[r]

    x_nr = np.array([X[i] for i in range(len(X)) if i!=r]) # leave out 'r'th column when writing to x_nr

    temp = np.exp(2 * x_r * (theta.dot(x_nr)))      #Evaluating the numerator in the conditional probability
    return(temp/(temp+1))                           #Return Conditional probability of x_r given x_nr

def reg_func(X_nr,theta_nr, param=None):
    """ Return the Normalized Regularization term
        Arguments:
        X_nr     : One sample of data with the 'r'th column removed
        theta_nr : The corresponding theta vector with the rth element removed
        param    : If param=1 do normalized l1 regularization else do trace lasso
    """
    if(param==1):
        l1_norm = 0
        for i in range(len(theta_nr)):
            l1_norm += np.abs(theta_nr[i]) * np.linalg.norm(X_nr[:,i])  #Return normalized l1 norm
        return l1_norm
    else:
        global mu
        M = X_nr.dot(np.diag(theta_nr))
        S = sqrtm(M.dot(np.transpose(M)) + mu * np.identity(n))
        Sinv = np.linalg.inv(S)
        return 0.5 * theta_nr.dot(np.diag(np.diagonal(np.transpose(X_nr).dot(Sinv).dot(X_nr)))).dot(theta_nr) + 0.5 * np.trace(S)

def obj_func(theta,r):
    objVal = 0

    for i in range(n):
        objVal -= cond(X[i,:],theta,r)
    
    X_nr = np.delete(X,r,1)

    #objVal += lamb*reg_func(X_nr,theta,1) #lasso  
    objVal += lamb*reg_func(X_nr,theta)  #Trace lasso

    return objVal

""""Perform Signed Edged recovery for every edge in the graph"""
for i in range(p):
    theta_arg = []
    """ The following loop constructs the theta vector for this particular node
        and basically eliminates duplicates.
    """
    for j in range(p):
        if(i > j):
            theta_arg.append(theta[j,i])
        if(i < j):
            theta_arg.append(theta[i,j])

    theta_arg = np.array(theta_arg) 

    #Run the Optimizer on the Objective function, method can be changed
    theta_opt = optimize.minimize(obj_func, theta_arg,args=(i,),method='CG')

    #Write the elements of minimized theta to their right place in the overall theta matrix
    #This step is mainly to avoid confusion between the nodes
    for j in range(p-1):
        if(i <= j):
            theta[i,j+1] = theta_opt.x[j]
        if(i > j):
            theta[j,i] = theta_opt.x[j]

#print the p choose 2 unique theta values of the graph
#Set a margin around 0 so that elements within the margin are made 0
#For the rest, print 1 if positive, -1 if negative
for i in range(p):
    temp = np.array([theta[i,j] for j in range(p) if i < j])

    temp = np.where(np.abs(temp) < 1e-8, 0,temp)
    temp = np.where(temp > 0, 1, temp)
    temp = np.where(temp < 0, -1, temp)
    print(temp)
