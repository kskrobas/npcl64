import OpenGL.GL as gl
import OpenGL.GLU as glu
import OpenGL.GLUT as glut
import numpy as np

from numpy_minmax import minmax
from order_of_magnitude import order_of_magnitude as om
import pygame.freetype  # Import the freetype module.

import pixDistance



    
#..............................................................................

def orderOfMagn(nmb):
    return np.floor(np.log10(nmb))
    
#..............................................................................
class StAxisPrms:
    def __init__(self,xo,yo,zo,W,L,H):
        self.xo=xo
        self.yo=yo
        self.zo=zo
        self.W=W
        self.L=L
        self.H=H
        
                
    def getPrms(self):
        return self.xo,self.yo,self.zo,self.W,self.L,self.H
#..............................................................................
    
    
def getScreenPos(coord3d,viewport,modelmat,projmat):
    
    if viewport == None:
        viewport = gl.glGetIntegerv (gl.GL_VIEWPORT);						# //get actual viewport
        modelmat = gl.glGetDoublev (gl.GL_MODELVIEW_MATRIX);				# //get actual model view matrix
        projmat  = gl.glGetDoublev (gl.GL_PROJECTION_MATRIX);			# //get actual projiection matrix    
                
    tx,ty,tz= glu.gluProject(coord3d[0],coord3d[1],coord3d[2], modelmat, projmat, viewport)    
    return np.array([tx,ty,tz])
    
    
    
def drawText(x, y, text):                                                
    textSurface = font.render(text, True, (25,25,22,25)).convert_alpha()
    textData = pygame.image.tostring(textSurface, "RGBA", True)
    gl.glWindowPos2d(x, y)
    gl.glDrawPixels(textSurface.get_width(), textSurface.get_height(), gl.GL_RGBA, gl.GL_UNSIGNED_BYTE, textData)
    
    
def pdist2(a,b):
    return np.sum((a-b)**2)
        
#------------------------------------------------------------------------------
def drawLines():
    gl.glBegin(gl.GL_LINE_STRIP)
    gl.glVertex3f(0, 0, 0)
    gl.glVertex3f(0.5,0,0)
    gl.glVertex3f(0.25,0.75,0)
    gl.glVertex3f(0,0,0)    
    gl.glEnd()
    
#------------------------------------------------------------------------------
def readXYZ(fileName):
    fin=open(fileName,'r')
    
    nOfa=fin.readline()
    fin.readline()
    
    atomPos=np.ndarray((int(nOfa),3),float)
    atomNames=np.ndarray((int(nOfa),3),int)
    
    for i,apos in enumerate(fin):
        fline=apos.split()
        
        atomPos[i,0]=float(fline[1])
        atomPos[i,1]=float(fline[2])
        atomPos[i,2]=float(fline[3])
        
        if fline[0]=='Zn':
            atomNames[i,:]=np.array([1,1,0])
        else:
            atomNames[i,:]=np.array([0,1,1])
            
    fin.close()
    
    return atomPos,atomNames
        
#------------------------------------------------------------------------------    
def drawSphere():
              
    Quadric = glu.gluNewQuadric()    
        
    gl.glPushMatrix()
        
    gl.glTranslate(-0.25,-0.25,0)        
    gl.glColor(1,0,0)
    glu.gluSphere(Quadric, 0.2, 25, 25)
    
    
    gl.glTranslate(.5,0,0)        
    gl.glColor(1,1,0)
    glu.gluSphere(Quadric, 0.2, 25, 25)
    

    gl.glTranslate(0,0.5,0)        
    gl.glColor(0,1,0)
    glu.gluSphere(Quadric, 0.2, 25, 25)
    
    gl.glTranslate(-0.5,0,0)        
    gl.glColor(1,0,1)
    glu.gluSphere(Quadric, 0.2, 25, 25)
    
    gl.glPopMatrix()
    
#------------------------------------------------------------------------------
def drawAtoms(apos,anames):
    Quadric = glu.gluNewQuadric()  
    
    
    for aitr,axyz in enumerate(apos):
        gl.glPushMatrix()
        
        gl.glTranslatef(axyz[0],axyz[1],axyz[2])        
        gl.glColor(anames[aitr,0],anames[aitr,1],anames[aitr,2])
        glu.gluSphere(Quadric,1,10,10)
        
        gl.glPopMatrix()
  
#------------------------------------------------------------------------------
def createListOfAtoms(apos,anames,listId):
  
    gl.glNewList(listId,gl.GL_COMPILE)
    Quadric = glu.gluNewQuadric()  
    
    for aitr,axyz in enumerate(apos):
        gl.glPushMatrix()
        
        gl.glTranslatef(axyz[0],axyz[1],axyz[2])        
        gl.glColor(anames[aitr,0],anames[aitr,1],anames[aitr,2])
        glu.gluSphere(Quadric,1,10,10)
        
        gl.glPopMatrix()
        
    gl.glEndList()    
#------------------------------------------------------------------------------
def drawAxes(axisPrm):
    xo,yo,zo,w,l,h=axisPrm.getPrms()  
    
                
    xf,xw=xo,xo+w
    
    xnsp=getScreenPos(np.array([xo,yo,-zo]),None,None,None)
    xpsp=getScreenPos(np.array([xw,yo,-zo]),None,None,None)
    
    xPixDist=pdist2(xnsp,xpsp)**(0.5)    
    minDistX,nOfStepsX=pixDistance.pixelDistance(xPixDist,w)
    
    if nOfStepsX>0:
        xf=np.floor(xo/minDistX)*minDistX
        xw=np.ceil((xw)/minDistX)*minDistX
        
        
        
    yf,yl=yo,yo+l
    
    ynsp=getScreenPos(np.array([xo,yo,-zo]),None,None,None)
    ypsp=getScreenPos(np.array([xo,yl,-zo]),None,None,None)
    
    yPixDist=pdist2(ynsp,ypsp)**(0.5)    
    minDistY,nOfStepsY=pixDistance.pixelDistance(yPixDist,l)
    
    if nOfStepsY>0:
        yf=np.floor(yo/minDistY)*minDistY
        yl=np.ceil((yl)/minDistY)*minDistY
                
     
    
    gl.glPushMatrix()  
    gl.glLineWidth(5)    
    
    
    gl.glBegin(gl.GL_LINES)
    
    
    gl.glColor(1,0,0)
    gl.glVertex3f(xf,yf,-zo)
    gl.glVertex3f(xw,yf,-zo)
    
    gl.glColor(0,1,0)        
    gl.glVertex3f(xf,yf,-zo)
    gl.glVertex3f(xf,yl,-zo)
    
    gl.glColor(0,0,1)   
    gl.glVertex3f(xf,yf,-zo)
    gl.glVertex3f(xf,yf,-zo-h)
    
    gl.glEnd()
    
    gl.glPopMatrix()
   
    '''
    xorou=np.round(xo,0)
    
    xyo=getScreenPos(np.array([xorou,yo,-zo]),None,None,None)
    xoc=getScreenPos(np.array([0,yo,-zo]),None,None,None)
    xyw=getScreenPos(np.array([xorou+w,yo,-zo]),None,None,None)
    sxyo=f'{xorou:5.1f}'
    sxyw=f'{xorou+w:5.1f}'
    #print(sxyo)
    
    pygame.font.SysFont('arial', 14)
    if pdist2(xyo,xyw)>10000:
        drawText(xyo[0]-50,xyo[1]-50,sxyo)
        drawText(xyw[0]-50,xyw[1]-50,sxyw)
        
    drawText(xoc[0]-10,xoc[1]-50,"0")
    '''
    
        
    
    
    
#------------------------------------------------------------------------------
def drawBox(w,l,h,xo,yo,zo):
    
          
    
    gl.glPushMatrix()    
    gl.glPushAttrib (gl.GL_LINE_BIT);
    gl.glLineStipple (3, 0xAAAA);
    gl.glColor(0.2,0.2,0.2)
    gl.glLineWidth(2)
        
    gl.glBegin(gl.GL_LINE_LOOP)
        
    gl.glVertex3f(xo,yo,zo)
    gl.glVertex3f(xo+w,yo,zo)
    gl.glVertex3f(xo+w,yo+l,zo)
    gl.glVertex3f(xo,yo+l,zo)
        
    gl.glEnd()
    
    
    gl.glBegin(gl.GL_LINE_LOOP)
        
    gl.glVertex3f(xo,yo,zo+h)
    gl.glVertex3f(xo+w,yo,zo+h)
    gl.glVertex3f(xo+w,yo+l,zo+h)
    gl.glVertex3f(xo,yo+l,zo+h)
        
    gl.glEnd()
        
    
    gl.glBegin(gl.GL_LINES)
    
    gl.glVertex3f(xo,yo,zo)
    gl.glVertex3f(xo,yo,zo+h)
        
    gl.glVertex3f(xo+w,yo,zo)
    gl.glVertex3f(xo+w,yo,zo+h)    
    
    gl.glVertex3f(xo+w,yo+l,zo)
    gl.glVertex3f(xo+w,yo+l,zo+h)
    
    gl.glVertex3f(xo,yo+l,zo)
    gl.glVertex3f(xo,yo+l,zo+h)
    
    gl.glEnd()
    
    gl.glPopMatrix()
    gl.glPopAttrib ();
        
    
#------------------------------------------------------------------------------
def drawPoint(xo,yo,zo):
    gl.glColor3f(0,0,0)
    gl.glPointSize(10)
    
    gl.glBegin(gl.GL_POINTS)    
    
    gl.glVertex3f(xo,yo,zo)
    
    gl.glEnd()
    
    
    

    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------    

apos,anames=readXYZ('/home/fizyk/python/pygnome/znsDiam.xyz')

center=np.average(apos,axis=0)
bpos=apos-center
bX,bY,bZ=minmax(bpos[:,0]),minmax(bpos[:,1]),minmax(bpos[:,2])
boxW=bX[1]-bX[0]
boxL=bY[1]-bY[0]
boxH=bZ[1]-bZ[0]
bsize=np.max(np.array([boxW,boxL,boxH]))

axisPrms=StAxisPrms(bX[0],bY[0], bZ[0], boxW,boxL,boxH)


boxW*=1.5
boxL*=1.5
boxH*=1.5
bsize*=1.5
xo,yo,zo=-0.5*bsize,-0.5*bsize,-0.5*bsize

angleX=-1
angleY=-1

mLastPos=(0,0)
mcurrPos=(0,0)
rotX=0
rotY=0
rotZ=0
drx=dry=0.25

red = (255,0,0)


import pygame
#from Mesh3D import *
from pygame.locals import *
done = False
pygame.init()
screen_width = 800
screen_height = 800
screen = pygame.display.set_mode((screen_width, screen_height), DOUBLEBUF | OPENGL)

gl.glEnable(gl.GL_BLEND)
gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
pygame.display.set_caption("Opengl.gl w Pythonie")
font = pygame.font.SysFont('arial', 32)
gl.glMatrixMode(gl.GL_PROJECTION)
#gl.gluPerspective(45, screen_width/screen_height, 0.1, 1000)
gl.glOrtho(-bsize,bsize,-bsize,bsize,-bsize,bsize)

transZ=-2
gl.glMatrixMode(gl.GL_MODELVIEW)
gl.glTranslatef(0,0,transZ)

gl.glEnable (gl.GL_LINE_STIPPLE)
gl.glEnable(gl.GL_DEPTH_TEST)
gl.glEnable(gl.GL_LIGHTING)
gl.glEnable(gl.GL_LIGHT0)

#glLightfv(GL_LIGHT0, GL_AMBIENT, (1,1,1, 1))
#glLightfv(GL_LIGHT0, GL_DIFFUSE, (bsize,bsize,bsize, 0.5))
#glLightfv(GL_LIGHT0, GL_SPECULAR, (0, 1, 0, 1))

gl.glEnable(gl.GL_COLOR_MATERIAL)
gl.glColorMaterial(gl.GL_FRONT,gl.GL_AMBIENT_AND_DIFFUSE)
gl.glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, (0, 1, 0, 1))
gl.glShadeModel(gl.GL_SMOOTH)

gl.glLight(gl.GL_LIGHT0, gl.GL_POSITION, (bsize,bsize,bsize,1))
gl.glLineWidth(2)
 


gList=gl.glGenLists(1)
createListOfAtoms(bpos, anames, gList)
mLastPos=None
#------------------------------------------------------------------------------


while not done:
            
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    gl.glClearColor(0.8,0.8,0.8,0)
    
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
        if event.type == pygame.VIDEORESIZE:
            gl.glMatrixMode(gl.GL_PROJECTION)
            #gluPerspective(45, 1, 0.1, 1000)
            gl.glOrtho(-0.5*boxW,0.5*boxW,-0.5*boxL,.5*boxL,-0.5*boxH,0.5*boxH)

            gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glTranslatef(0,0,transZ)
            gl.glRotatef(rotY , 1.0, 0.0, 0.0);
            gl.glRotatef(rotX , 0.0, 1.0, 0.0);
            
            
        keys=pygame.key.get_pressed()
        if keys[pygame.K_UP]:
            gl.glScalef(0.9,0.9,0.9)
        if keys[pygame.K_DOWN]:
            gl.glScalef(1.11,1.11,1.11)
            
        if keys[pygame.K_LEFT]:
            angleX+=0.25
            gl.glLight(gl.GL_LIGHT0, gl.GL_POSITION, (angleX,angleY,1, 1))
            
        if keys[pygame.K_RIGHT]:
            angleX-=0.25
            gl.glLight(gl.GL_LIGHT0, gl.GL_POSITION, (angleX,angleY,1, 1))
        
        
        if event.type == MOUSEWHEEL:
           if event.y>0:
               bsize*=1.125
           else:
               bsize/=1.125
               
           
           gl.glMatrixMode(gl.GL_PROJECTION)
           gl.glLoadIdentity()
           gl.glOrtho(-bsize,bsize,-bsize,bsize,-bsize,bsize)

           gl.glMatrixMode(gl.GL_MODELVIEW)
           gl.glLoadIdentity()
           gl.glTranslatef(0,0,transZ)
           gl.glRotatef(rotY , 1.0, 0.0, 0.0);
           gl.glRotatef(rotX , 0.0, 1.0, 0.0);           
      


        if event.type == MOUSEBUTTONUP:
            mLastPos=None

        if pygame.mouse.get_pressed()[0]:                                    
            if mLastPos == None:
                mLastPos=pygame.mouse.get_pos()
            
            if event.type == MOUSEMOTION:                                              
            #print(mcurrPos,len(mcurrPos))
                mcurrPos=pygame.mouse.get_pos()
                mdiffPosX=float(mcurrPos[0])-float(mLastPos[0])
                mdiffPosY=float(mcurrPos[1])-float(mLastPos[1])
                mLastPos= pygame.mouse.get_pos()
                
                rotX+=drx*mdiffPosX
                rotY+=dry*mdiffPosY
                
                gl.glLoadIdentity()
                gl.glTranslatef(0,0,transZ)
                gl.glRotatef(rotY , 1.0, 0.0, 0.0);
                gl.glRotatef(rotX , 0.0, 1.0, 0.0);
                
                


    gl.glCallList(gList)
    drawBox(boxW,boxL,boxH,xo, yo, zo)
    drawAxes(axisPrms)
    #drawPoint(0,0,zo)
    #drawText(screen_width*0.5,screen_height*0.5, "(0,0,0)")
    #drawAxes(boxW,boxL,boxH,xo,yo,zo)
        
    pygame.display.flip()
    
pygame.quit()
