# GF_direct.py
#   This program is made to calculate the Time Independent
# Green's Function of two interacting hard core bosons on a lattice
# through the exact diagonalization of the Hamiltonian.
# The on-site energy, tunnelling parameter, and interaction
# parameter can all be changed or set to a random distribution.
# The tunnelling range and interaction range can also be modified,
# in addition to the crystal size.

#First date of completion: 12 Feb 2015 - JTC
#13 Feb 2015 - JTC - Switched GF calculation from S*I*S^t to only using the specific rows/columns necessary

from __future__ import division
import numpy as np
import sys
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

np.set_printoptions(threshold=10000,linewidth=2000,precision=4,suppress=False)


E = 5.0 #NOTE: This is the on-site energy of a SINGLE particle
t = 1.0
d = 2.0

N = 50

E_random = True
T_random = False
D_random = False
rand_seed = 52

T_range = N #NOTE: the form is t / dist^3
D_range = N #NOTE: the form is d / dist^3

#These are for the Green's function G(z), z = omega + i*eta
omega_start = -8.0
omega_stop = 8.0 #NOTE: np.arange is used, so the endpoint could be exceeded a bit
omega_step = 0.025
eta = 0.001
# omega : [0.1,5.5] step 0.025

zeroErrorTolerance = 1E-7 #Means results are single precision or better

####################################################################
####################################################################
# Functions
####################################################################
####################################################################

####################################################################
# Calculate the Green's Function <mn|G(z)|ij>
####################################################################

# <mn|G(z)|ij> = sum_a <mn|a><a|ij>/(z-E_a) = S*(I*1/(z-E_a))*S^t
# where S is the matrix of eigenvectors where the columns are <mn|a>
# and I is the identity matrix, thus (I*1/(z-E_a)) is a diagonal matrix
# in the eigenvector basis with values 1/(z-E_a)
# BUT, I can only calculate the relevant rows, instead of the whole matrices.
# Thus, I can do S_(mn)a * E_a(z) * S_(ij)a^* (Einstein notation) for each
# desired pair of (m,n) and (i,j) instead of S_(mn)a * E_ab(z) * S_(ij)b^*
# for all pairs.
# Returns GF_2DArray[nKet,z]. ie. first index corresponds to the matrix element, and the second to the value of z

def GF_Calc (omega_start, omega_stop, omega_step, eta, leftKets, rightKets, eval_Array_srtd, evec_2DArray_srtd):
    
    global basis_to_n
    
    # Make the z-value array; NOTE: the interval is [start, stop)
    z_Array = np.arange(omega_start, omega_stop + omega_step, omega_step) + 1j*eta
    num_z = z_Array.size
#    print num_z

    #Get the relvant indices
    leftKetsSize = len(leftKets)
    #    print leftKetsSize
    if leftKetsSize != len(rightKets):
        sys.exit("ERROR in GF_Calc: Number of left and right kets do not match.")
    
    leftIndices = np.zeros(leftKetsSize)
    rightIndices = np.zeros(leftKetsSize)
    
    for n in range(0, leftKetsSize):
        leftIndices[n] = basis_to_n[leftKets[n]]
        rightIndices[n] = basis_to_n[rightKets[n]]
    
    #    print n_to_basis
    #    print leftKets
    #    print leftIndices
    #    print rightKets
    #    print rightIndices
    
    # GF_2DArray[nKet,z]
    GF_2DArray = np.zeros((leftKetsSize,num_z),dtype=np.complex128)

    for i in range(0, num_z):
        
        #1/(z-E_a) array
        E_Array = 1/(z_Array[i] - eval_Array_srtd)
            
        for nKet in range(0, leftKetsSize):
            
            #tmp_a = E_a * S_(ij)a^* (element by element multiplication, no sum)
            tmp = E_Array*(evec_2DArray_srtd[rightIndices[nKet],:].conj())
            
            #GF = S_(mn)a * tmp_a = S_(mn)a * E_a * S_(ij)a^* (Einstein Notation)
            GF_2DArray[nKet,i] = np.dot(evec_2DArray_srtd[leftIndices[nKet],:], tmp)
        
    return GF_2DArray, z_Array


####################################################################
####################################################################
# Program
####################################################################
####################################################################

ti = time.time()
print "--------------------------------------------------------------"
print "GF_direct.py"
print "------------"

print " "
print "System Parameters"
print "-----------------"
print "E (On-Site Energy/Disorder Width): %.3F" % (E)
print "t (Tunnelling Parameter): %.3F" % (t)
print "d (Interaction Parameter): %.3F" % (d)
print " "

print "Number of Lattice Sites: %d" % (N)
print " "

print "On-site Energy Randomized: %r" % (E_random)
print "Tunnelling Parameter Randomized: %r" % (T_random)
print "Interaction Parameter Randomized: %r" % (D_random)
print "Random Seed: %s" % (rand_seed)
print " "

print "Tunnelling Range: %d lattice sites" % (T_range)
print "Interaction Range: %d lattice sites" % (D_range)
print " "

print "Green's Function Parameters"
print "---------------------------"
print "z = omega + i*eta"
print "Omega Start: %.3F" % (omega_start)
print "Omega Stop: %.3F" % (omega_stop)
print "Omega Step: %.3F" % (omega_step)
print "Eta: %.3G" % (eta)
print "Number of z-values: %d" % (np.ceil( ((omega_stop + omega_step) - omega_start)/omega_step)) #See np.arange documentation

print " "
print "Calculations Beginning:"
print "-----------------------"



####################################################################
# Set up the basis/index conversion arrays
####################################################################

if N < 2:
    sys.exit("Error: Can't have fewer lattice sites than hard core bosons.")
#This converts the tuple list to a numpy array
#n_to_basis_tmp = []
#
#for i in range(0,N):
#    for j in range(i+1,N):
#        n_to_basis_tmp.append((i,j))
#
#n_to_basis = np.array(n_to_basis_tmp)
#
#print n_to_basis_tmp
#print n_to_basis

#This uses the tuple list directly
n_to_basis = []

for i in range(0,N):
    for j in range(i+1,N):
        n_to_basis.append((i,j))

basis_size = len(n_to_basis)
#print basis_size

#print n_to_basis

#Set up a dictionary for the reverse look up (could this be too slow?)
basis_to_n = {}

i = 0
for ket in n_to_basis:
    basis_to_n[ket] = i
    i += 1

#print basis_to_n
#print basis_to_n[(1,2)]
#print basis_to_n[(1,4)] #This line will error out with a key error (if N =< 4), which is the desired behaviour

print "Basis built, time: %G s" % (time.time() - ti)
t2 = time.time()

####################################################################
# Set up the basis/index conversion arrays
####################################################################


# Set up the On Site Energy array (diagonal in Hamiltonian)
####################################################################
if E_random == False:
    E_onSite_Array = np.ones(basis_size) * 2.0 * E

elif E_random == True:
    np.random.seed(rand_seed)
    
    e_array = np.random.uniform(-E/2.,E/2.,N) #Set a random value for each site

    E_onSite_Array = np.zeros(basis_size)
    for i in range(0,basis_size):
        m, n = n_to_basis[i]
        E_onSite_Array[i] = e_array[m] + e_array[n] #<m,n|E|m,n> = E_m + E_n

#    print n_to_basis
#    print e_array
#    print E_onSite_Array

else:
#    sys.exit("Error: Random on-site energy not yet implemented.")
    sys.exit("Error: Invalid value for 'E_random'.")

#print E_onSite_Array

# Set up the Ineteraction Energy array (diagonal in Hamiltonian)
####################################################################
if D_random == False:
    D_interaction_Array = np.zeros(basis_size)

    for i in range(0,basis_size):
        dist = np.abs(n_to_basis[i][1] - n_to_basis[i][0])

        if dist <= D_range:
            D_interaction_Array[i] += d / (dist**3)
#            print "Yes, dist: ", dist
#
#        else:
#            print "No, dist: ", dist

else:
    sys.exit("Error: Random interaction energy not yet implemented.")

#print D_interaction_Array


# Set up the Tunnelling Probability Matrix
####################################################################
if T_random == False:

# The below looks at element <mn|T|ab> and the following logic is in the
# if statements:
# if either a or b equals m or n, then <mn|T|ab> = <an|T|ab> = <ma|T|ab> = <bn|T|ab> = <mb|T|ab> = t / dist**3 (ie. one particle hops)
# if neither a or b equals m or n, then <mn|T|ab> = 0 (ie. both particles would hop)
# if a and b equal m or n, then <mn|T|ab> = <ab|T|ab> = <ba|T|ab> = 0

    T_tunnelling_2DArray = np.zeros((basis_size,basis_size))

    for i in range(0, basis_size):
        for j in range(i+1, basis_size):
            m, n = n_to_basis[i]
            a, b = n_to_basis[j]

#            print m, n
#            print a, b
#            print "____"

            if (a == m) and (b != n):
                dist = b - n

            elif (a == n) and (b != m):
                dist = b - m

            elif (b == m) and (a != n):
                dist = a - n

            elif (b == n) and (a != m):
                dist = a - m

            else:
                dist = 0

            dist = np.abs(dist)

            #put in 0.5 below instead of 0.0 in case there is a numerical error and dist is slightly bigger than 0; it would never be 0.5 larger, however. Ideally, I would have "if dist > 0:"
            if dist > 0.5:
                T_tunnelling_2DArray[i,j] += t / (dist**3)

    T_tunnelling_2DArray += T_tunnelling_2DArray.T

#    print T_tunnelling_2DArray
#
#    T2 = np.copy(T_tunnelling_2DArray)
#    print T2.T - T_tunnelling_2DArray



else:
    sys.exit("Error: Random tunnelling probability not yet implemented.")



# Set up the Hamiltonian Matrix
####################################################################
H_2DArray = np.diag(E_onSite_Array,0)
#print H_2DArray

H_2DArray += np.diag(D_interaction_Array,0)
#print H_2DArray

H_2DArray += T_tunnelling_2DArray
#print H_2DArray

print "Hamiltonian built, time: %G s" % (time.time() - t2)
t2 = time.time()

####################################################################
# Diagonalize the Hamiltonian
####################################################################
# eval_Array is an array of Eigenvalues, NOT sorted
# evec_2DArray is the corresponding "array" of eigenvectors, hence a 2D array with elements <ij|n>, where |n> is an eigenvector
eval_Array, evec_2DArray = np.linalg.eig(H_2DArray)

print "Hamiltonian diagonalized, time: %G s" % (time.time() - t2)
t2 = time.time()

eval_Array_srtd = np.sort(eval_Array)

#sort the eigenvectors, so that they are in order of increasing energy
evec_2DArray_srtd = np.copy(evec_2DArray[:,eval_Array.argsort()])

#print eval_Array
#print evec_2DArray
#print eval_Array_srtd
#print evec_2DArray_srtd

#Check normality
normCheck = np.dot(evec_2DArray_srtd.conj().T,evec_2DArray_srtd)
normCheck2 = np.abs(normCheck - np.identity(basis_size))
if np.any(normCheck2 > zeroErrorTolerance) :
    print "ERROR: Normality problem"
    print normCheck
    print "Largest value: ", np.max(normCheck2)
    sys.exit("ERROR: Normality problem, check the above matrix")
else:
    print "Normality of wavefunctions confirmed, all elements of |S^t*S - I| < %G" % (zeroErrorTolerance)

#Show range of eigenvalues
eval_min = np.min(eval_Array_srtd)
eval_max = np.max(eval_Array_srtd)

print "Eigenvalue range: %G to %G" % (eval_min, eval_max)


#####################################################################
## Calculate the Green's Function <mn|G(z)|ij>
#####################################################################
#
## <mn|G(z)|ij> = sum_a <mn|a><a|ij>/(z-E_a) = S*(I*1/(z-E_a))*S^t
## where S is the matrix of eigenvectors where the columns are <mn|a>
## and I is the identity matrix, thus (I*1/(z-E_a)) is a diagonal matrix
## in the eigenvector basis with values 1/(z-E_a)
#
## Make the z-value array; NOTE: the interval is [start, stop)
#z_Array = np.arange(omega_start, omega_stop + omega_step, omega_step) + 1j*eta
#num_z = z_Array.size
#
## Make the 1/(z-E_a) matrix for each value of z first; thus a 3D Array
#GF_3DArray = np.zeros((basis_size,basis_size,num_z),dtype=np.complex128)
#
#for i in range(0, num_z):
#    #1/(z-E_a) matrix
#    GF_3DArray[:,:,i] = np.diag(1/(z_Array[i] - eval_Array_srtd))
#
#    #1/(z-E_a) matrix * S^t
#    GF_3DArray[:,:,i] = np.dot(GF_3DArray[:,:,i], evec_2DArray_srtd.conj().T)
#
#    #S *(1/(z-E_a) matrix * S^t)
#    GF_3DArray[:,:,i] = np.dot(evec_2DArray_srtd, GF_3DArray[:,:,i])
#
##    print GF_3DArray[:,:,i]
#
## This was to check that the arithmetic is correct, confirmed by the below code and with a calculator.
##    den = (z_Array[i].real - eval_Array_srtd).astype(np.complex128)
##    den += 1j*z_Array[i].imag
##    print 1/den
##    print np.diag(GF_3DArray[:,:,i]) - 1/den
##    print GF_3DArray[:,:,i]
#
##print z_Array
##print num_z
##print GF_3DArray
##print eval_Array_srtd
#
#leftKets = [(0,1), (2,3), (1,3)]
#rightKets = [(2,3), (0,1), (1,3)]
#
#GF_2DArray = GF_Calc(omega_start, omega_stop, omega_step, eta, leftKets, rightKets, eval_Array_srtd, evec_2DArray_srtd)
#
##Compare the result with the old method
#
##Get the relvant indices
#
#leftKetsSize = len(leftKets)
#
#leftIndices = np.zeros(leftKetsSize,dtype=np.int64)
#rightIndices = np.zeros(leftKetsSize,dtype=np.int64)
#
#for n in range(0, leftKetsSize):
#    leftIndices[n] = basis_to_n[leftKets[n]]
#    rightIndices[n] = basis_to_n[rightKets[n]]
#
#diff = GF_2DArray - GF_3DArray[leftIndices,rightIndices,:]
##print diff
#zeroErrorTolerance2 = 1E-13
#if np.any(np.abs(diff) > zeroErrorTolerance2) :
#    print "ERROR: GF's Different"
#    print np.abs(diff)
#    print "Largest value: ", np.max(np.abs(diff))
#    sys.exit("ERROR: GF's Different, check the above matrix")
#else:
#    print "GF's are the same; all differences < %G" % (zeroErrorTolerance2)
#
##print GF_2DArray
##print GF_3DArray[leftIndices,rightIndices,:]
##print GF_2DArray.shape
##print GF_3DArray.shape
##print GF_3DArray[leftIndices,rightIndices,:].shape
#
#####################################################################
##I now have the Green's Function, but to what do I compare the result?
#####################################################################
#
#print "Green's Function built, time: %G s" % (time.time() - t2)
#t2 = time.time()








#print "Program completed, time: %G s" % (time.time() - ti)
#print "All times are the Wall Time."
#print "--------------------------------------------------------------"

####################################################################
#Plotting
####################################################################



figNum = 0

# Plot |<m,n|G(z)|i,j>| as a function of omega

leftKets = [(0,1)]
rightKets = [(0,1)]

m,n = leftKets[0]
i,j = rightKets[0]

leftIndex = basis_to_n[(m,n)]
rightIndex = basis_to_n[(i,j)]

GF_2DArray, z_Array = GF_Calc(omega_start, omega_stop, omega_step, eta, leftKets, rightKets, eval_Array_srtd, evec_2DArray_srtd)

#print GF_2DArray.shape

GF_squared = GF_2DArray[0].conj()*GF_2DArray[0]

#print leftIndex
#print rightIndex
#print n_to_basis
#print basis_to_n

#GF_squared = GF_3DArray[leftIndex,rightIndex,:].conj()*GF_3DArray[leftIndex,rightIndex,:]
#print GF_squared

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
line, = ax.plot(z_Array.real, np.sqrt(GF_squared.real))

custom_yLabel = r'$|<{0},{1}|G(z)|{2},{3}>|$'.format(m+1,n+1,i+1,j+1)

#plt.ylim(0,50)
plt.xlabel(r'$\omega + {0}i$'.format(eta), fontsize=16)
plt.ylabel(custom_yLabel, fontsize=16)
plt.title("Green's Function", fontsize=18)



####################################################################
#Calculate the localization length
####################################################################

#Calculate |<m,m+a|G|n,n+a>| ~ exp(-|n-m|/xi_2) as |m-n| -> inf

#Make a 2D plot with x-axis (m-n), then make a 3D plot with an additional axis with a

a = 1

n = int(np.round(N/2))
m_start = 0
m_stop = N-1

m_Array = np.arange(m_start, m_stop-a+1, 1) #interval [start, stop)

#print m_Array

rightIndexLocLen = basis_to_n[(n,n+a)]

leftIndexLocLen_Array = np.zeros(m_Array.size)
leftKetsLocLen_Array = []
rightKetsLocLen_Array = []

for i in range(0,m_Array.size):
    leftIndexLocLen_Array[i] = basis_to_n[(m_Array[i],m_Array[i]+a)]

    leftKetsLocLen_Array.append((m_Array[i],m_Array[i]+a))
    rightKetsLocLen_Array.append((n,n+a))


#print rightIndexLocLen
#print leftIndexLocLen_Array
#print n_to_basis
#print leftKetsLocLen_Array
#print rightKetsLocLen_Array

GF_2DArray, z_Array = GF_Calc(omega_start, omega_stop, omega_step, eta, leftKetsLocLen_Array, rightKetsLocLen_Array, eval_Array_srtd, evec_2DArray_srtd)

#print GF_2DArray.shape

GF_squared_2DArray = GF_2DArray.conj()*GF_2DArray

#GF_partial_2DArray = np.zeros((m_Array.size, num_z)) # (m,z) indexing
#
#for i in range(0,m_Array.size):
#    GF_partial_2DArray[i,:] = (GF_3DArray[leftIndexLocLen_Array[i],rightIndexLocLen,:].conj()*GF_3DArray[leftIndexLocLen_Array[i],rightIndexLocLen,:]).real

print "Relevant Green's Functions built, time: %G s" % (time.time() - t2)
t2 = time.time()

figNum += 1
fig = plt.figure(figNum)
ax = fig.gca(projection='3d')
Y = m_Array.copy()
X = (z_Array.real).copy()

#print Y

X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, np.sqrt(GF_squared_2DArray.real), rstride=100, cstride=1, cmap=cm.hsv, linewidth=0, antialiased=False)

custom_zLabel = r'$|<m,m+{0}|G(z)|{1},{1}+{0}>|$'.format(a,n+1)

#plt.ylim(0,50)
ax.set_ylabel(r'$m$', fontsize=16)
ax.set_xlabel(r'$\omega + {0}i$'.format(eta), fontsize=16)
#ax.set_zlabel(custom_zLabel, fontsize=16)
ax.set_title(custom_zLabel, fontsize=18)

skip = 10
skip_Array = [i for i in range(0,m_Array.size,skip)]
#print skip_Array
ax.set_yticks(m_Array[skip_Array])
ax.set_yticklabels(m_Array[skip_Array])

#ax.set_zlim(0, 25)


fig.colorbar(surf, shrink=0.5, aspect=5)







#print "Green's Function built, time: %G s" % (time.time() - t2)
#t2 = time.time()


print "Program completed, time: %G s" % (time.time() - ti)
print "All times are the Wall Time."
print "--------------------------------------------------------------"


print "Now for plots."

plt.show()





























