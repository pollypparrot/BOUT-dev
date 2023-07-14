# Holly Chandler
# Masterclass Research Project
# Code to support maths required in the X-Z laplacian for various functions

#import all neccesary libraries and allow use of unicode symbols
import sympy
from sympy import symbols, expand, factor, latex
from sympy import *
init_printing(use_unicode=True)
#import library allowing ext and latex to be on the same line:
from IPython.display import display, Math


#set it so you can have pretty outputs for multiple different functions in one cell
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = 'all'


#initialise all the symbols used for each coordinate system
alpha, gamma, beta, x, y, z, r, theta, phi, R, A, B, D, h, I, Bpol  = symbols('alpha gamma beta x y z r theta phi R A B D h_theta I B_pol')


#order for all lists of components follows:

#0 cartesian (xyz),
#1 cyllindrical (r,theta,z),
#2 spherical (r,theta,phi),
#3 toroidal (r,theta,phi),
#4 Clebsch Field Alligned Coordinates (x,y,z)

#initialsise coordinate system order:
C = ["cartesian (xyz)","cyllindrical (r,theta,z)","spherical (r,theta,phi)","toroidal (r,theta,phi)","Clebsch Field Alligned Coordinates (x,y,z)" ]

#Initialise all the known Jacobians
J = [1, r, r**(2)*sin(theta),-r*(R+r*cos(theta)),h/Bpol]

#initialise the main diagonal metric tensors in order g11,g22,g33 where 1,2,3 refer to the coordinate system variables in order as below
g = [[1,1,1],[1,1/(r**2),1],[1,1/(r**2),1],[1,1/(r**2),1/(R+r*cos(theta))**2],[(R*Bpol)**2, 1/(h)**2,B**2/(R*Bpol)**2]]

#all the coordinate system variables.
coords = [[x,y,z],[r,theta,z],[r,theta,phi],[r,theta,phi],[x,y,z]]


#create a function to independently calculate b
def calculate_b_1(A,D,f,J,gxx,gzz,z_coord,x_coord):
    dx = diff(J*gxx*diff(f,x_coord),x_coord)
    dz = diff(J*gzz*diff(f,z_coord),z_coord)
    b = A/J*(dx+dz)+D*f
    return b

def calculate_b_2(A,D,f,J,gxx,gzz,z_coord,x_coord,gxz):
    b = A/J*(diff(J*gxx*diff(f, x_coord),x_coord)+(diff(J*gzz*diff(f, z_coord),z_coord))+diff(J*gxz*diff(f,z),x)+diff(J*gxz*diff(f,z),z))+D*f


def convertPowers(function):
    #first replace all the ** with commas as this is what splits the item to the power
    function = function.replace("**",",")

    #add space to start and end in case of digit being a power-involved digit
    function = " " + function + " "

    #count the number of powers in the string to loop through them
    numberOfPowers = function.count(",")

    #inistialse the start of the search to be the beginning of the string
    whereToSearch = 0

    #loop through each power situation
    for item in range (0,numberOfPowers):

        #find the next instance of a comma (demonstrating a power)
        powerIndex = function.find(',',whereToSearch)

        #change where to start the search start point to just after the location of the comma just found
        # ensure different instances are checked each time
        # +5 as +1 to search after the current comma and +4 to account for the extra indexs' from inserting pow( before comma 
        whereToSearch = powerIndex + 5

        #first start checking where the bracket should go after/before the comma i.e. where the power ends.
        numberIndexEnd = powerIndex + 1

        #START WITH EXPONENT
        #start by checking if the exponent is a function i.e. in brackets.
        if function[numberIndexEnd]=="(":
            #now increase the index so we are looking at brackets past the first open one
            #if there are brackets within the exponent bracket then these need to be considered- numberClosed counts them
            finished = False
            numberIndexEnd +=1
            numberClosed = 0
            #while loop to find end of function
            while not finished:
                #look for the next closed and open bracket after the initial open one
                closedBracketIndexBack = function.find(")",numberIndexEnd)
                openBracketIndexBack = function.find("(",numberIndexEnd)

                if (openBracketIndexBack > closedBracketIndexBack)or(openBracketIndexBack == -1):
                    #then the open bracket is after the closed for the exponent
                    #then there are no open brackets after the last point

                    #now check if there are any seconary opened brackets
                    #if numberClosed == 0 then there arent any to be considered- end loop
                    if numberClosed == 0:
                        function = str(function[ :closedBracketIndexBack]+")"+function[closedBracketIndexBack: ])
                        finished = True
                    #if there are opened brackets then we close one off numberClosed -1
                    #we change where we look for the next closed bracket to be one past the current closed bracket
                    else:
                        numberClosed = numberClosed -1
                        numberIndexEnd = closedBracketIndexBack +1
                #otherwise there is an open bracket before the closed for the exponent
                #we consider this by looking for closed brackets after and noting we need more than one closed bracket now
                else:
                    numberClosed +=1
                    numberIndexEnd = openBracketIndexBack+1

        #or if there is no bracket
        else:
            #if it is not bracketed then it is just a number so we check each charater until we don't reach a digit.
            #when we reach a non-digit we know its the end of the exponent and add the bracket
            end = False
            while not end:
                if function[numberIndexEnd].isdigit():
                    numberIndexEnd +=1
                else:
                    function = str(function[ :numberIndexEnd]+")"+function[numberIndexEnd: ])
                    end = True

        #NOW DO BASE of power:

        numberIndexStart =  powerIndex - 1
        #if the base is bracketed then we search for the open bracket and insert pow( there

        if function[numberIndexStart]==")":
            found = False
            index = 1
            numberClosed = 0
            #while loop that looks at each character indiviually and continues till brackets are closed
            while not found:
                #identify the next character
                nextIndex = numberIndexStart-index
                nextCharacter = function[nextIndex]

                #if there are no more brackets to closed and we found the final closed bracket- insert pow(
                if (nextCharacter == "(") and (numberClosed ==0):

                    #if there is a sin or cos etc in front of the brackets we want to put the pow before this
                    #awkwardSymbold[] is where you can input such symbols
                    awkwardSymbols = ["cos","sin","tan","cosh", "sinh", "tanh", "log", "cot"]
                    #as they can be 3 or 4 length units we have both options
                    find3Letter = function[nextIndex-3:nextIndex]
                    find4Letter = function[nextIndex-4:nextIndex]
                    #check if it meets each item in the list
                    #if it does we change the nextIndex to that position so that pow is inserted there
                    for item in awkwardSymbols:
                        if (str(item) == str(find3Letter)):
                            nextIndex = nextIndex -3
                        if item ==find4Letter:
                            nextIndex = nextIndex -4
                    found = True
                    function = str(function[ :nextIndex]+"pow("+function[nextIndex: ])
                #otherwise add one to index - look at next character before this one
                else:
                    index +=1
                    #if the current character was a close brackets- we need to find another open in order to insert pow(
                    #add one to numberClosed
                    if nextCharacter == ")":
                        numberClosed +=1
                    #otherwise if it was a closed bracket and numberClosed !=1 then we remove one from number closed
                    #i.e we are closer to closing the surrounding bracket of the base
                    if nextCharacter =="(":
                        numberClosed -=1

        #if it is a number base we loop until we find a non number. Insert pow there.
        else:
            end = False
            while not end:
                if function[numberIndexStart].isdigit():
                    numberIndexStart =numberIndexStart -1
                else:
                    function = str(function[ :numberIndexStart]+"pow("+function[numberIndexStart: ])
                    end = True
    print(function)

#example
convertPowers("23**2+(236-(x+y))**(4+x)+x*cos(x)**(2+5*(x+2)+6(x-2))")


#function finds the symbols of the coordinate system and outputs the result of calculating the forward laplacian

def calculateForwardLaplacian(f,coordSystemChoice,xCoordIndex,zCoordIndex):
    #choose x coordinate unit based on coordinate system and which combination
    x_coord = coords[coordSystemChoice][xCoordIndex]
    z_coord = coords[coordSystemChoice][zCoordIndex]

    #assign x and z the units found above to ensure that the function inputted is changed to have those units
    x= x_coord
    z= z_coord

    #calculate result of laplacian by calling calculae_b function
    b = calculate_b_1(A,B,f,J[coordSystemChoice],g[coordSystemChoice][xCoordIndex],g[coordSystemChoice][zCoordIndex],z_coord,x_coord)

    #print out all the values used
    print("Coordinate system choice is", C[coordSystemChoice])
    display(Math(f'Jacobian = {latex(J[coordSystemChoice])}'))
    display(Math(f'X-Coordinate = {latex(x)}'))
    display(Math(f'Z-Coordinate = {latex(z)}'))
    display(Math(f'f = {latex(f)}'))
    display(Math(f'b = {latex(b)}'))
    print ()
    print("sympy result")
    print(b)
    print()
    print("c++ conversion")
    print(convertPowers(str(b)))
    print()

#example
calculateForwardLaplacian(sin(x*y)*y*z,2,1,0)


def loopAllCoordinatesAndSystems(f):
    #for each coordinate system all combinations of units are used i.e 01, 12, 20.
    #coordinate systems are shown as above
    #to change the function that is being inputted into the laplacian, you change f which is seen within the for loops

    #for loop to go through each coordinate system
    for i in range(0,4):
        #for loop for x to go through each basis vector
        for j in range (0,3):
            #for loop to ensure that each combination of basis vectors are seen i.e 01,12,21,02,10,20
            for k in range(0,2):

                #set the coordinate system and coord choice through for loop
                coordSystemChoice = i
                coordChoiceX = j

                #% is modulus which is the remainder when (j+1) is divided by 3. Y goes 1,2,0.
                coordChoiceZ = (j+k+1)%3

                #calculate result:
                calculateForwardLaplacian(f,coordSystemChoice,coordChoiceX,coordChoiceZ)

#example
loopAllCoordinatesAndSystems(sin(x)*x*y*1/z)



def clebschCoordinatesForwardLaplacian(f):
    #set the coordinate system and coord choice through for loop
    coordSystemChoice = 4

    #choose x and z coordinate to e x and z.
    x_coord = coords[coordSystemChoice][0]
    z_coord = coords[coordSystemChoice][2]

    #assign x and z the units found above to ensure that the function inputted is changed to have those units
    x= x_coord
    z= z_coord

    #NEW STUFF
    #write Jacobian, metric tensors and functions for B and R
    #start with B functions as then, when the Jacobian and metric tensors are applied, the functions are subbed in.

    eps = 1
    r = x
    theta = y
    Rxy = r
    Btor = Rxy
    Bpol = r
    Bxy = sqrt(Bpol**2 + Btor**2)
    h = r
    JClebsch = h/Bpol
    clebschMetricTensor = [(Rxy*Bpol)**2, 1/(h)**2,Bxy**2/(Rxy*Bpol)**2]

    A = 1
    D = 1

    #calculate result of laplacian by calling calculae_b function
    b = calculate_b_1(A,D,f,JClebsch,clebschMetricTensor[0],clebschMetricTensor[2],z_coord,x_coord)

    #print out all the values used
    print("Coordinate system choice is Clebsch")
    print("Facts about the system")
    display(Math(f'X-Coordinate = {latex(x)}'))
    display(Math(f'Z-Coordinate = {latex(z)}'))
    display(Math(f'Bpol = {latex(Bpol)}'))
    display(Math(f'Btor = {latex(Btor)}'))
    display(Math(f'Bxy = {latex(Bxy)}'))
    display(Math(f'gxx = {latex(clebschMetricTensor[0])}'))
    display(Math(f'gzz = {latex(clebschMetricTensor[2])}'))
    print("Inputted function")
    display(Math(f'f = {latex(f)}'))
    print("Result")
    display(Math(f'b = {latex(b)}'))
    print("sympy result")
    print(b)
    print()
    print("c++ conversion")
    print(convertPowers(str(b)))
    print()

#example
clebschCoordinatesForwardLaplacian(x)
