#!/usr/bin/python

# This is part of the Vertigo suite.
# Niko Suenderhauf
# Chemnitz University of Technology
# niko@etit.tu-chemnitz.de


from optparse import OptionParser
import sys
import random
from math import *


# =================================================================
def checkOptions(options):
    """Make sure the options entered by the user make sense."""

    if options.outliers<0:
        print("Number of outliers (--outliers) must be >=0.")
        return False
    
    if options.groupsize<0:
        print("Groupsize (--groupsize) must be >=0.")
        return False

    if options.filename == "" or options.filename==None:
        print("Dataset to read (--in) must be given.")
        return False

    if options.information == "":
        options.information=None
    
    
    if options.information:
        if (options.information.count(",") != 0) and (options.information.count(",") != 5) and  (options.information.count(",") != 20):
            print("Information matrix must be given in full upper-triangular form. E.g. --information=42,0,0,42,0,42 or as a single value that is used for all diagonal entries, e.g. --information=42.")
            return False
    


    return True
    
    
# =================================================================
def readDataset(filename, vertexStr='VERTEX_SE2', edgeStr='EDGE_SE2'):

    # read the complete file
    f = open(filename,'r')
    lines=f.readlines()


    # determine whether this is a 3D or 2D dataset
    mode=None
    
    for i in range(len(lines)):
        if lines[i].startswith("VERTEX_SE2"):
            mode=2
            break
        elif lines[i].startswith("VERTEX_SE3"):
            mode=3
            break
        
                                       
    if mode == 2:
        vertexStr='VERTEX_SE2'
        edgeStr='EDGE_SE2'               
    elif mode == 3:
        vertexStr='VERTEX_SE3:QUAT'
        edgeStr='EDGE_SE3:QUAT'
       
    # build a dictionary of vertices and edges
    v=[]
    e=[]
    
    for line in lines:
        if line.startswith(vertexStr):
            v.append(line)

        elif line.startswith(edgeStr):
            e.append(line)

    f.close()

    return (v,e, mode)


# ==================================================================
def euler_to_quat(yaw,  pitch,  roll):
   sy = sin(yaw*0.5);
   cy = cos(yaw*0.5);
   sp = sin(pitch*0.5);
   cp = cos(pitch*0.5);
   sr = sin(roll*0.5);
   cr = cos(roll*0.5);
   w = cr*cp*cy + sr*sp*sy;
   x = sr*cp*cy - cr*sp*sy;
   y = cr*sp*cy + sr*cp*sy;
   z = cr*cp*sy - sr*sp*cy;

   return (w,x,y,z)



# ==================================================================
def writeDataset(filename, vertices, edges, mode, outliers=0, groupSize=1, doLocal=0, informationMatrix="42,0,0,42,0,42", perfectMatch=False):
  
    # first write out all pose vertices (no need to change them)
    f = open(filename, 'w')
    for n in vertices:
        f.write(n)


    # edges and vertices are called differently, depending on 2D or 3D mode ...
    if mode == 2:
        edgeStr='EDGE_SE2'

        # check entries for information matrix for additional loop closure constraints (outliers)
        if not informationMatrix:
            print("Determining information matrix automatically...")                                      
        else:        
            if informationMatrix.count(",")==0:
                # if there is only a single value, convert it into full upper triangular form with that value on the diagonal
                try:                
                    diagEntry=float(informationMatrix)
                except:
                    print("! Invalid value for information matrix. If you give only a single value, it must be a number, e.g. --information=42")
                    return False
                informationMatrix = "%f,0,0,%f,0,%f" % (diagEntry,diagEntry,diagEntry)
        
            elif informationMatrix.count(",")!=5:
                print("! Invalid number of entries in information matrix. Full upper triangular form has to be given, e.g. --information=42,0,0,42,0,42.")
                return False
            
    elif mode == 3:
        edgeStr='EDGE_SE3'

        # check entries for information matrix for additional loop closure constraints (outliers)        
        if not informationMatrix:
            print("Determining information matrix automatically")                                      
        else:        
            if informationMatrix.count(",")==0:
                # if there is only a single value, convert it into full upper triangular form with that value on the diagonal
                try:                
                    diagEntry=float(informationMatrix)
                except:
                    print("! Invalid value for information matrix. If you give only a single value, it must be a number, e.g. --information=42")
                    return False
                informationMatrix = "%f,0,0,0,0,0,%f,0,0,0,0,%f,0,0,0,%f,0,0,%f,0,%f" % (diagEntry,diagEntry,diagEntry,diagEntry, diagEntry, diagEntry)
            elif informationMatrix.count(",")!=20:
                print("! Invalid number of entries in information matrix (--information). Full upper triangular form has to be given.")
                return False
    else:
        print("! Invalid mode. It must be either 2 or 3 but was", mode)
        return False
    

        
    # now for every edge we need to write out a switchable edge
    # therefore we need a switch vertex and an associated prior edge


    # how many poses are there?
    poseCount = len(vertices)

    # for every edge, create a new switch node and its associated prior and write the new switchable edges
    for oldStr in edges:

      (a,b) = oldStr.split()[1:3]
      if int(a) != int(b)-1: 
        isOdometryEdge = False
      else:
        isOdometryEdge = True

      # auto determine information matrix for additional outliers from the first loop closure edge we find in the dataset
      if not isOdometryEdge and informationMatrix==None:
          elem = oldStr.split()
          if mode==2:
              informationMatrix = ' '.join(elem[-6:])
          else:
              informationMatrix = ' '.join(elem[-21:])
          print(informationMatrix)
      
      f.write(oldStr)


    # now create the desired number of additional outlier edges
    for _ in range(outliers):        

        elem = oldStr.split()

        # determine random indices for the two vertices that are connected by an outlier edge
        v1=0
        v2=0
        while v1==v2:
            v1=random.randint(0,poseCount-1-groupSize)
            if doLocal<1:
              v2=random.randint(0,poseCount-1-groupSize)
            else: 
              v2=random.randint(v1,min(poseCount-1-groupSize, v1+20)) 

            if v1>v2:
              tmp=v1
              v1=v2
              v2=tmp                       
            if v2==v1+1:
              v2=v1+2
        
        # determine coordinates of the loop closure constraint
        if mode == 2:
            x1=random.gauss(0,0.3)
            x2=random.gauss(0,0.3)
            x3=random.gauss(0,10*pi/180.0)

        if mode == 3:
            x1=random.gauss(0,0.3)
            x2=random.gauss(0,0.3)
            x3=random.gauss(0,0.3)
            
            
            sigma = 10.0*pi/180.0
            roll = random.gauss(0,sigma)
            pitch = random.gauss(0,sigma)
            yaw = random.gauss(0,sigma)            
            (q0, q1, q2, q3) = euler_to_quat(yaw, pitch, roll)

        if perfectMatch:
          x1=x2=x3=0
          q0=1
          q1=q2=q3=0

        
        for _ in range(groupSize):

            info_str = informationMatrix.replace(",", " ")

            n=[v1, v2, x1, x2, x3]
            if mode == 3:
                n.extend([q0, q1, q2, q3])                    
                s = ' '.join([edgeStr+":QUAT"] + [str(x) for x in n]) + " " + info_str
            else:
                s = ' '.join([edgeStr] + [str(x) for x in n]) + " " + info_str
            f.write(s+'\n')

            v1=v1+1
            v2=v2+1    

    f.close()

    return True

# ==================================================================    
# ==================================================================
# ==================================================================

if __name__ == "__main__":


    # let's start by preparing to parse the command line options
    parser = OptionParser()
    
    # string or numeral options
    parser.add_option("-i", "--in", help = "Path to the original dataset file (in g2o format).", dest="filename")
    parser.add_option("-o", "--out", help = "Results will be written into this file.", default="new.g2o")
    parser.add_option("-n", "--outliers", help = "Spoil the dataset with this many outliers. Default = 100.", default=100, type="int")
    parser.add_option("-g", "--groupsize", help = "Use this groupsize. Default = 1.", default=1, type="int")
    parser.add_option("--information", help = "Set the information matrix for the additional false positive loop closure constraints. Provide either a single value e.g. --information=42 that will be used for all diagonal entries. Or provide the full upper triangular matrix using values separated by commas, but no spaces: --information=42,0,0,42,0,42 etc.")  #, default="42.7,0,0,42.7,0,42.7")
    parser.add_option("--seed", help = "Random seed. If >0 it will be used to initialize the random number generator to create repeatable random false positive loop closures.", default=None, type="int")

    # boolean options
    parser.add_option("-l", "--local", help = "Create only local false positive loop closure constraints.", action="store_true", default=False)
    parser.add_option("-p", "--perfectMatch", help = "Loop closures match perfectly, i.e. the transformation between both poses is (0,0,0).", action="store_true", default=False)

    # parse the command line options
    (options, args) = parser.parse_args()

    # check if the options are valid and sound       
    if checkOptions(options):

        # initialize the random number generator
        random.seed(options.seed)

        # read the original dataset
        print("Reading original dataset", options.filename, "...")        
        (vertices, edges, mode) = readDataset(options.filename)

        # build and save the modified dataset with additional false positive loop closures
        print("Writing modified dataset", options.out, "...")
        if writeDataset(options.out, vertices, edges, mode,
                     options.outliers,
                     options.groupsize,
                     options.local,
                     options.information,
                     options.perfectMatch):
            print("Done.")


    # command line options were not ok
    else: 
        print("Please use --help to see all available command line parameters.")


    
