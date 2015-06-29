#***************************************************************************
#*                                                                         *
#*   Copyright (c) 2010 Dan Falck <ddfalck@gmail.com>                      *
#*   derived from Yorik van Havre's <yorik@gmx.fr>  importDXF.py           *
#*   script that is part of the Draft plugin for FreeCAD                   *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License (GPL)            *
#*   as published by the Free Software Foundation; either version 2 of     *
#*   the License, or (at your option) any later version.                   *
#*   for detail see the LICENCE text file.                                 *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU Library General Public License for more details.                  *
#*                                                                         *
#*   You should have received a copy of the GNU Library General Public     *
#*   License along with this program; if not, write to the Free Software   *
#*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
#*   USA                                                                   *
#*                                                                         *
#***************************************************************************

'''
This script uses a DXF-parsing library created by Stani,
Kitsu and Migius for Blender

It is also based off the Heeks DXF importer
'''

import sys
sys.path.append('./lib')
import math
import numpy
from dxfReader import readDXF
from random import randint

def drawLine(line):
    if (len(line.points) > 1):
        v1=(line.points[0][0],line.points[0][1],line.points[0][2])
        v2=(line.points[1][0],line.points[1][1],line.points[1][2])

        x0=line.points[0][0];y0=line.points[0][1];z0=line.points[0][2]
        x1=line.points[1][0];y1=line.points[1][1];z1=line.points[1][2]

        if not equals(v1,v2):
            try: return ("cad.line3d("+str(x0)+","+str(y0)+","+str(z0)+","+str(x1)+","+str(y1)+","+str(z1)+")\n")
            except: warn('line',line)
    return None

def drawArc(arc):
    cen=(arc.loc[0],arc.loc[1],arc.loc[2])
    firstangle=(arc.start_angle/180)*math.pi
    lastangle=(arc.end_angle/180)*math.pi
    rad = arc.radius

    center = (str(arc.loc[0])+ ", " +str(arc.loc[1])+ ", " + str(arc.loc[2]))
    radius =  str(rad)
    angle1 = str(firstangle)
    angle2 = str(lastangle)
    dir_vec = "0,0,1"
    try: return ("cad.arc(" + center + ", " + radius + ", " + angle1 + ", " + angle2 + ", " + dir_vec + ")\n")
    except: warn('arc',arc)
    return None

#***************************************************************************
# functions copied from fcvec for polyline bulges
#***************************************************************************
def precision():
    return 6

def isNull(vector):
    '''isNull(vector): Tests if a vector is nul vector'''
    p = precision()
    return (round(vector[0],p)==0.0 and round(vector[1],p)==0.0 and round(vector[2],p)==0.0)

def equals(u,v):
    "returns True if vectors differ by less than precision (from ParamGet), elementwise "
    #typecheck ([(u,Vector), (v,Vector)], "equals")
    return isNull(numpy.subtract(u,v))

def isColinear(vlist):
    '''isColinear(list_of_vectors): checks if vectors in given list are colinear'''
    #typecheck ([(vlist,list)], "isColinear");
    if len(vlist) < 3: return True
    #first = vlist[1].sub(vlist[0])
    first = numpy.subtract(vlist[1],vlist[0])
    for i in range(2,len(vlist)):
        #if angle(vlist[i].sub(vlist[0]),first) != 0:
        if angle(numpy.subtract(vlist[i],vlist[0]),first) != 0:
            return False
    return True

def angle(u,v=(1,0,0),normal=(0,0,1)):
    '''angle(Vector,[Vector],[Vector]) - returns the angle in radians between the two vectors.
    If only one is given, angle is between the vector and the horizontal East direction.
    If a third vector is given, it is the normal used to determine the sign of the angle.
    '''
    #typecheck ([(u,Vector), (v,Vector)], "angle")
    #ll = u.Length*v.Length
    ll = numpy.linalg.norm(u)*numpy.linalg.norm(v)
    if ll==0: return 0
    #dp=u.dot(v)/ll
    dp=numpy.dot(u,v)/ll
    if (dp < -1): dp = -1 # roundoff errors can push dp out of the ...
    elif (dp > 1): dp = 1 # ...geometrically meaningful interval [-1,1]
    ang = math.acos(dp)
    #normal1 = u.cross(v)
    normal1 = numpy.cross(u,v)
    #coeff = normal.dot(normal1)
    coeff = numpy.dot(normal,normal1)
    if coeff >= 0:
        return ang
    else:
        return -ang

#*******************************************************************
#   polyline related functions
#*******************************************************************

def calc_center(v1,bulge,v2):
    '''
    calculates center of arc- this one works
    ''' 
    chord = numpy.subtract(v2,v1)
    chord_length= numpy.linalg.norm(chord)
    sagitta = (bulge*chord_length)/2.0
    inc_angle = numpy.arctan(bulge)*4.0
    radius = (chord_length/2.0)/numpy.sin(inc_angle/2.0)
    if bulge >= 0:
        perp = (numpy.cross(chord,(0,0,-1)))
        #perp = (numpy.cross(chord,(-1,-1,-1)))
    else:
        perp = (numpy.cross(chord,(0,0,1)))
        #perp = (numpy.cross(chord,(1,1,1)))
        #sagitta = sagitta*8
        radius = -radius
    chord_mid_pt = numpy.add(numpy.multiply(chord,(.5,.5,.5)),v1)
    unit_vec = perp/ numpy.linalg.norm(perp)
    arc_center = numpy.add(numpy.multiply((radius-sagitta),unit_vec),chord_mid_pt)
    return arc_center

def distanceFormula(x1,y1,x2,y2):
    a = (x2-x1)*(x2-x1)
    b = (y2-y1)*(y2-y1)
    return math.sqrt(a+b)

def addDegrees(base,mod):
    # this function expects a 360 degree number
    # base and mod must be between 0-360
    v = base+mod
    if v > 360:
        v = 360-v
    elif v < 0:
        v = 360+v
    return abs(v)

def newPointFromDistanceAndAngle(x,y,ang,distance):
    r = []
    r.append(x+(distance*math.cos(ang*math.pi/180)))
    r.append(y+(distance*math.sin(ang*math.pi/180)))
    return r

def drawPolyline(polyline):
    # returns a list of points for a polyline with line segments replacing arcs

    global minX
    global maxX
    global minY
    global maxY

    if (len(polyline.points) > 1):
        collector = []
        for p in range(len(polyline.points)-1):
            p1 = polyline.points[p]
            p2 = polyline.points[p+1]
            v1 = (p1[0],p1[1],p1[2])
            v2 = (p2[0],p2[1],p2[2])

            if p == 0:
                minX = p1[0]
                maxX = p1[0]
                minY = p1[1]
                maxY = p1[1]
            else:
                if p1[0] < minX:
                    minX = p1[0]
                elif p1[0] > maxX:
                    maxX = p1[0]
                if p1[1] < minY:
                    minY = p1[1]
                elif p1[1] > maxY:
                    maxY = p1[1]

            if not equals(v1,v2):
                
                if polyline.points[p].bulge:
                    # get the center point for the bulge
                    cv = calc_center(v1,polyline.points[p].bulge,v2)

                    if isColinear([v1,cv,v2]):
                        # the bulge points are on the same line, stupid editor
			#print 'LINE IS COLINEAR'
                        try: 
			    collector.append([p1[0],p1[1]])
                            collector.append([p2[0],p2[1]])
                        except: 
                            warn('polyline',polyline)
                    else:
                        start=(str(v1[0])+","+str(v1[1])+","+str(v1[2]))
                        center=(str(cv[0])+","+str(cv[1])+","+str(cv[2]))
                        end = (str(v2[0])+","+str(v2[1])+","+str(v2[2]))

                        r = distanceFormula(p1[0],p1[1],cv[0],cv[1])

                        # first add the start point
                        # this is not needed as the first point is at startAng
                        #collector.append([p1[0],p1[1]])

                        #
                        #   2   |   1
                        #       |
                        # ------|-------
                        #       |
                        #   3   |   4
                        #
                        # first find start point quadrant relative to cv
                        # and end point quadrant relative to cv

                        # start point
                        if p1[0] > cv[0] and p1[1] > cv[1]:
                            startPointQuad = 1
                        elif p1[0] < cv[0] and p1[1] > cv[1]:
                            startPointQuad = 2
                        elif p1[0] < cv[0] and p1[1] < cv[1]:
                            startPointQuad = 3
                        elif p1[0] > cv[0] and p1[1] < cv[1]:
                            startPointQuad = 4

                        # end point
                        if p2[0] > cv[0] and p2[1] > cv[1]:
                            endPointQuad = 1
                        elif p2[0] < cv[0] and p2[1] > cv[1]:
                            endPointQuad = 2
                        elif p2[0] < cv[0] and p2[1] < cv[1]:
                            endPointQuad = 3
                        elif p2[0] > cv[0] and p2[1] < cv[1]:
                            endPointQuad = 4

                        # start angle from cv to p1
                        startSlope = (cv[1]-p1[1]) / (cv[0]-p1[0])
                        startAng = 180*math.atan(startSlope)/math.pi
                        if polyline.points[p].bulge >= 0:
                            # positive bulge
                            if startPointQuad == 2:
                                startAng = 180+startAng
                            elif startPointQuad == 3:
                                startAng = 180+startAng
                            elif startPointQuad == 4:
                                startAng = 360+startAng
                        else:
                            # negative bulge
                            if startPointQuad == 2:
                                startAng = 180+startAng
                            elif startPointQuad == 3:
                                startAng = 180+startAng
                            elif startPointQuad == 4:
                                startAng = 360+startAng
                        #print "startAng", startAng

                        # end angle from cv to p2
                        endSlope = (cv[1]-p2[1]) / (cv[0]-p2[0])
                        endAng = 180*math.atan(endSlope)/math.pi
                        if polyline.points[p].bulge < 0:
                            # negative bulge
                            if endPointQuad == 2:
                                endAng = 180+endAng;
                            elif endPointQuad == 3:
                                endAng = 180+endAng;
                            elif endPointQuad == 4:
                                endAng = 360+endAng;
                        else:
                            # positive bulge
                            if endPointQuad == 2:
                                endAng = 180+endAng;
                            elif endPointQuad == 3:
                                endAng = 180+endAng;
                            elif endPointQuad == 4:
                                endAng = 360+endAng;
                        #print "endAng", endAng

                        if polyline.points[p].bulge < 0:
                            # this is a negative bulge so it will be an arc that goes from p1 to p2 except
                            # it will be bulging toward the cv point and not away from it like normal
                            arcTotalDeg = 360-addDegrees(endAng,-startAng)
                        else:
                            arcTotalDeg = addDegrees(endAng,-startAng)

                        # now we need to create the line segments in the arc
                        numSegments = 40
                        degreeStep = arcTotalDeg/numSegments

                        #print "\nPOINT",p
			#print "p1",p1
			#print "p2",p2
			#print "cv",cv
			#print 'startPointQuad',startPointQuad
			#print 'endPointQuad',endPointQuad
			#print 'bulge',polyline.points[p].bulge
			#print 'startAng',startAng
			#print 'endAng',endAng
			#print 'arcTotalDeg',arcTotalDeg

                        # now loop through each degreeStep
                        for c in range(1, numSegments+1):
                            # for a positive bulge the start point is always a lower number of degrees
                            if polyline.points[p].bulge < 0:
                                # so for a negative bulge we need to subtract degreeStep
                                pt = newPointFromDistanceAndAngle(cv[0],cv[1],addDegrees(startAng,-(degreeStep*c)),r)
                            else:
                                # and for a positive bulge we add degreeStep
                                pt = newPointFromDistanceAndAngle(cv[0],cv[1],addDegrees(startAng,(degreeStep*c)),r)
                                # reverse direction
                            collector.append([pt[0],pt[1]])

                else:

                    try: 
                        # this is a normal line segment, add it
			collector.append([p1[0],p1[1]])
                        collector.append([p2[0],p2[1]])
                    except: 
                        warn('polyline',polyline)
                
            else:
                return "no way!\n"

        #string = "".join(collector)
        return collector
    else:
        return "No\n"
#*******************************************************************
#   end of polyline related functions
#*******************************************************************
def warn(type,dxfobject):
    print "dxf: couldn't import type "+type+":", dxfobject.layer

# process file
def process(filename):
    # process a file
    global drawing 
    global minX
    global maxX
    global minY
    global maxY
    drawing = readDXF(filename)
    mcOut = ''

    lines = drawing.entities.get_type("line")    
    for line in lines:
        shape = drawLine(line)
        if shape:
            print "\nFOUND UNSUPPORTED LINE OBJECT\n"

    arcs = drawing.entities.get_type("arc")    
    for arc in arcs:
        shape = drawArc(arc)
        if shape:
            print "\nFOUND UNSUPPORTED ARC OBJECT\n"

    polylines = drawing.entities.get_type("polyline")
    polylines.extend(drawing.entities.get_type("lwpolyline"))
    i = 0
    for polyline in polylines:
        shape = drawPolyline(polyline)
        if shape:

            mcOut += "\n// "+polyline.layer+"\n"
            mcOut += "var polygon"+str(i)+" = {type:'polygon',points:["
            for oooo in shape:
		mcOut += '['+str(oooo[0])+','+str(oooo[1])+'],';
            mcOut += "]};"
            mcOut += "\nmc.cut('centerOnPath',polygon"+str(i)+", 4, [0,0]);\n"
            i += 1

    mcOut += '\nmc.get();\n'

    totalX = maxX-minX
    totalY = maxY-minY
    mcOut = 'var tool = {units:\'mm\',diameter:6.35,passDepth:4,step:1,rapid:2000,plunge:100,cut:600,zClearance:5,returnHome:true};\n\nvar mc = new Millcrum(tool);\n\nmc.surface('+str(totalX*2)+','+str(totalY*2)+');\n' + mcOut

    if len(sys.argv) > 2:
        # output to file
        print "PROCESSING DXF FILE "+sys.argv[1]+" TO "+sys.argv[2]+"\n"
        f = open(sys.argv[2], 'w')
        f.write(mcOut+"\n")
        f.close()
    else:
        # output to STDOUT
        print "PROCESSING DXF FILE "+sys.argv[1]+" TO STDOUT\n"
        print mcOut+"\n"

if len(sys.argv) > 1:
    process(sys.argv[1])
else:
    print 'USAGE: python dxfToMillcrum.py input.dxf output.millcrum\n'
    print 'If an output file is not specified the Millcrum code will be output to STDOUT\n'
